#include "ElectronRepulsionGeom0100ContrRecHKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_hkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hkxx,
                                            const size_t idx_gkxx,
                                            const size_t idx_geom_01_gkxx,
                                            const size_t idx_geom_01_glxx,
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
            /// Set up components of auxilary buffer : GKSS

            const auto gk_off = idx_gkxx + i * dcomps + j;

            auto g_xxxx_xxxxxxx = cbuffer.data(gk_off + 0 * ccomps * dcomps);

            auto g_xxxx_xxxxxxy = cbuffer.data(gk_off + 1 * ccomps * dcomps);

            auto g_xxxx_xxxxxxz = cbuffer.data(gk_off + 2 * ccomps * dcomps);

            auto g_xxxx_xxxxxyy = cbuffer.data(gk_off + 3 * ccomps * dcomps);

            auto g_xxxx_xxxxxyz = cbuffer.data(gk_off + 4 * ccomps * dcomps);

            auto g_xxxx_xxxxxzz = cbuffer.data(gk_off + 5 * ccomps * dcomps);

            auto g_xxxx_xxxxyyy = cbuffer.data(gk_off + 6 * ccomps * dcomps);

            auto g_xxxx_xxxxyyz = cbuffer.data(gk_off + 7 * ccomps * dcomps);

            auto g_xxxx_xxxxyzz = cbuffer.data(gk_off + 8 * ccomps * dcomps);

            auto g_xxxx_xxxxzzz = cbuffer.data(gk_off + 9 * ccomps * dcomps);

            auto g_xxxx_xxxyyyy = cbuffer.data(gk_off + 10 * ccomps * dcomps);

            auto g_xxxx_xxxyyyz = cbuffer.data(gk_off + 11 * ccomps * dcomps);

            auto g_xxxx_xxxyyzz = cbuffer.data(gk_off + 12 * ccomps * dcomps);

            auto g_xxxx_xxxyzzz = cbuffer.data(gk_off + 13 * ccomps * dcomps);

            auto g_xxxx_xxxzzzz = cbuffer.data(gk_off + 14 * ccomps * dcomps);

            auto g_xxxx_xxyyyyy = cbuffer.data(gk_off + 15 * ccomps * dcomps);

            auto g_xxxx_xxyyyyz = cbuffer.data(gk_off + 16 * ccomps * dcomps);

            auto g_xxxx_xxyyyzz = cbuffer.data(gk_off + 17 * ccomps * dcomps);

            auto g_xxxx_xxyyzzz = cbuffer.data(gk_off + 18 * ccomps * dcomps);

            auto g_xxxx_xxyzzzz = cbuffer.data(gk_off + 19 * ccomps * dcomps);

            auto g_xxxx_xxzzzzz = cbuffer.data(gk_off + 20 * ccomps * dcomps);

            auto g_xxxx_xyyyyyy = cbuffer.data(gk_off + 21 * ccomps * dcomps);

            auto g_xxxx_xyyyyyz = cbuffer.data(gk_off + 22 * ccomps * dcomps);

            auto g_xxxx_xyyyyzz = cbuffer.data(gk_off + 23 * ccomps * dcomps);

            auto g_xxxx_xyyyzzz = cbuffer.data(gk_off + 24 * ccomps * dcomps);

            auto g_xxxx_xyyzzzz = cbuffer.data(gk_off + 25 * ccomps * dcomps);

            auto g_xxxx_xyzzzzz = cbuffer.data(gk_off + 26 * ccomps * dcomps);

            auto g_xxxx_xzzzzzz = cbuffer.data(gk_off + 27 * ccomps * dcomps);

            auto g_xxxx_yyyyyyy = cbuffer.data(gk_off + 28 * ccomps * dcomps);

            auto g_xxxx_yyyyyyz = cbuffer.data(gk_off + 29 * ccomps * dcomps);

            auto g_xxxx_yyyyyzz = cbuffer.data(gk_off + 30 * ccomps * dcomps);

            auto g_xxxx_yyyyzzz = cbuffer.data(gk_off + 31 * ccomps * dcomps);

            auto g_xxxx_yyyzzzz = cbuffer.data(gk_off + 32 * ccomps * dcomps);

            auto g_xxxx_yyzzzzz = cbuffer.data(gk_off + 33 * ccomps * dcomps);

            auto g_xxxx_yzzzzzz = cbuffer.data(gk_off + 34 * ccomps * dcomps);

            auto g_xxxx_zzzzzzz = cbuffer.data(gk_off + 35 * ccomps * dcomps);

            auto g_xxxy_xxxxxxx = cbuffer.data(gk_off + 36 * ccomps * dcomps);

            auto g_xxxy_xxxxxxy = cbuffer.data(gk_off + 37 * ccomps * dcomps);

            auto g_xxxy_xxxxxxz = cbuffer.data(gk_off + 38 * ccomps * dcomps);

            auto g_xxxy_xxxxxyy = cbuffer.data(gk_off + 39 * ccomps * dcomps);

            auto g_xxxy_xxxxxyz = cbuffer.data(gk_off + 40 * ccomps * dcomps);

            auto g_xxxy_xxxxxzz = cbuffer.data(gk_off + 41 * ccomps * dcomps);

            auto g_xxxy_xxxxyyy = cbuffer.data(gk_off + 42 * ccomps * dcomps);

            auto g_xxxy_xxxxyyz = cbuffer.data(gk_off + 43 * ccomps * dcomps);

            auto g_xxxy_xxxxyzz = cbuffer.data(gk_off + 44 * ccomps * dcomps);

            auto g_xxxy_xxxxzzz = cbuffer.data(gk_off + 45 * ccomps * dcomps);

            auto g_xxxy_xxxyyyy = cbuffer.data(gk_off + 46 * ccomps * dcomps);

            auto g_xxxy_xxxyyyz = cbuffer.data(gk_off + 47 * ccomps * dcomps);

            auto g_xxxy_xxxyyzz = cbuffer.data(gk_off + 48 * ccomps * dcomps);

            auto g_xxxy_xxxyzzz = cbuffer.data(gk_off + 49 * ccomps * dcomps);

            auto g_xxxy_xxxzzzz = cbuffer.data(gk_off + 50 * ccomps * dcomps);

            auto g_xxxy_xxyyyyy = cbuffer.data(gk_off + 51 * ccomps * dcomps);

            auto g_xxxy_xxyyyyz = cbuffer.data(gk_off + 52 * ccomps * dcomps);

            auto g_xxxy_xxyyyzz = cbuffer.data(gk_off + 53 * ccomps * dcomps);

            auto g_xxxy_xxyyzzz = cbuffer.data(gk_off + 54 * ccomps * dcomps);

            auto g_xxxy_xxyzzzz = cbuffer.data(gk_off + 55 * ccomps * dcomps);

            auto g_xxxy_xxzzzzz = cbuffer.data(gk_off + 56 * ccomps * dcomps);

            auto g_xxxy_xyyyyyy = cbuffer.data(gk_off + 57 * ccomps * dcomps);

            auto g_xxxy_xyyyyyz = cbuffer.data(gk_off + 58 * ccomps * dcomps);

            auto g_xxxy_xyyyyzz = cbuffer.data(gk_off + 59 * ccomps * dcomps);

            auto g_xxxy_xyyyzzz = cbuffer.data(gk_off + 60 * ccomps * dcomps);

            auto g_xxxy_xyyzzzz = cbuffer.data(gk_off + 61 * ccomps * dcomps);

            auto g_xxxy_xyzzzzz = cbuffer.data(gk_off + 62 * ccomps * dcomps);

            auto g_xxxy_xzzzzzz = cbuffer.data(gk_off + 63 * ccomps * dcomps);

            auto g_xxxy_yyyyyyy = cbuffer.data(gk_off + 64 * ccomps * dcomps);

            auto g_xxxy_yyyyyyz = cbuffer.data(gk_off + 65 * ccomps * dcomps);

            auto g_xxxy_yyyyyzz = cbuffer.data(gk_off + 66 * ccomps * dcomps);

            auto g_xxxy_yyyyzzz = cbuffer.data(gk_off + 67 * ccomps * dcomps);

            auto g_xxxy_yyyzzzz = cbuffer.data(gk_off + 68 * ccomps * dcomps);

            auto g_xxxy_yyzzzzz = cbuffer.data(gk_off + 69 * ccomps * dcomps);

            auto g_xxxy_yzzzzzz = cbuffer.data(gk_off + 70 * ccomps * dcomps);

            auto g_xxxy_zzzzzzz = cbuffer.data(gk_off + 71 * ccomps * dcomps);

            auto g_xxxz_xxxxxxx = cbuffer.data(gk_off + 72 * ccomps * dcomps);

            auto g_xxxz_xxxxxxy = cbuffer.data(gk_off + 73 * ccomps * dcomps);

            auto g_xxxz_xxxxxxz = cbuffer.data(gk_off + 74 * ccomps * dcomps);

            auto g_xxxz_xxxxxyy = cbuffer.data(gk_off + 75 * ccomps * dcomps);

            auto g_xxxz_xxxxxyz = cbuffer.data(gk_off + 76 * ccomps * dcomps);

            auto g_xxxz_xxxxxzz = cbuffer.data(gk_off + 77 * ccomps * dcomps);

            auto g_xxxz_xxxxyyy = cbuffer.data(gk_off + 78 * ccomps * dcomps);

            auto g_xxxz_xxxxyyz = cbuffer.data(gk_off + 79 * ccomps * dcomps);

            auto g_xxxz_xxxxyzz = cbuffer.data(gk_off + 80 * ccomps * dcomps);

            auto g_xxxz_xxxxzzz = cbuffer.data(gk_off + 81 * ccomps * dcomps);

            auto g_xxxz_xxxyyyy = cbuffer.data(gk_off + 82 * ccomps * dcomps);

            auto g_xxxz_xxxyyyz = cbuffer.data(gk_off + 83 * ccomps * dcomps);

            auto g_xxxz_xxxyyzz = cbuffer.data(gk_off + 84 * ccomps * dcomps);

            auto g_xxxz_xxxyzzz = cbuffer.data(gk_off + 85 * ccomps * dcomps);

            auto g_xxxz_xxxzzzz = cbuffer.data(gk_off + 86 * ccomps * dcomps);

            auto g_xxxz_xxyyyyy = cbuffer.data(gk_off + 87 * ccomps * dcomps);

            auto g_xxxz_xxyyyyz = cbuffer.data(gk_off + 88 * ccomps * dcomps);

            auto g_xxxz_xxyyyzz = cbuffer.data(gk_off + 89 * ccomps * dcomps);

            auto g_xxxz_xxyyzzz = cbuffer.data(gk_off + 90 * ccomps * dcomps);

            auto g_xxxz_xxyzzzz = cbuffer.data(gk_off + 91 * ccomps * dcomps);

            auto g_xxxz_xxzzzzz = cbuffer.data(gk_off + 92 * ccomps * dcomps);

            auto g_xxxz_xyyyyyy = cbuffer.data(gk_off + 93 * ccomps * dcomps);

            auto g_xxxz_xyyyyyz = cbuffer.data(gk_off + 94 * ccomps * dcomps);

            auto g_xxxz_xyyyyzz = cbuffer.data(gk_off + 95 * ccomps * dcomps);

            auto g_xxxz_xyyyzzz = cbuffer.data(gk_off + 96 * ccomps * dcomps);

            auto g_xxxz_xyyzzzz = cbuffer.data(gk_off + 97 * ccomps * dcomps);

            auto g_xxxz_xyzzzzz = cbuffer.data(gk_off + 98 * ccomps * dcomps);

            auto g_xxxz_xzzzzzz = cbuffer.data(gk_off + 99 * ccomps * dcomps);

            auto g_xxxz_yyyyyyy = cbuffer.data(gk_off + 100 * ccomps * dcomps);

            auto g_xxxz_yyyyyyz = cbuffer.data(gk_off + 101 * ccomps * dcomps);

            auto g_xxxz_yyyyyzz = cbuffer.data(gk_off + 102 * ccomps * dcomps);

            auto g_xxxz_yyyyzzz = cbuffer.data(gk_off + 103 * ccomps * dcomps);

            auto g_xxxz_yyyzzzz = cbuffer.data(gk_off + 104 * ccomps * dcomps);

            auto g_xxxz_yyzzzzz = cbuffer.data(gk_off + 105 * ccomps * dcomps);

            auto g_xxxz_yzzzzzz = cbuffer.data(gk_off + 106 * ccomps * dcomps);

            auto g_xxxz_zzzzzzz = cbuffer.data(gk_off + 107 * ccomps * dcomps);

            auto g_xxyy_xxxxxxx = cbuffer.data(gk_off + 108 * ccomps * dcomps);

            auto g_xxyy_xxxxxxy = cbuffer.data(gk_off + 109 * ccomps * dcomps);

            auto g_xxyy_xxxxxxz = cbuffer.data(gk_off + 110 * ccomps * dcomps);

            auto g_xxyy_xxxxxyy = cbuffer.data(gk_off + 111 * ccomps * dcomps);

            auto g_xxyy_xxxxxyz = cbuffer.data(gk_off + 112 * ccomps * dcomps);

            auto g_xxyy_xxxxxzz = cbuffer.data(gk_off + 113 * ccomps * dcomps);

            auto g_xxyy_xxxxyyy = cbuffer.data(gk_off + 114 * ccomps * dcomps);

            auto g_xxyy_xxxxyyz = cbuffer.data(gk_off + 115 * ccomps * dcomps);

            auto g_xxyy_xxxxyzz = cbuffer.data(gk_off + 116 * ccomps * dcomps);

            auto g_xxyy_xxxxzzz = cbuffer.data(gk_off + 117 * ccomps * dcomps);

            auto g_xxyy_xxxyyyy = cbuffer.data(gk_off + 118 * ccomps * dcomps);

            auto g_xxyy_xxxyyyz = cbuffer.data(gk_off + 119 * ccomps * dcomps);

            auto g_xxyy_xxxyyzz = cbuffer.data(gk_off + 120 * ccomps * dcomps);

            auto g_xxyy_xxxyzzz = cbuffer.data(gk_off + 121 * ccomps * dcomps);

            auto g_xxyy_xxxzzzz = cbuffer.data(gk_off + 122 * ccomps * dcomps);

            auto g_xxyy_xxyyyyy = cbuffer.data(gk_off + 123 * ccomps * dcomps);

            auto g_xxyy_xxyyyyz = cbuffer.data(gk_off + 124 * ccomps * dcomps);

            auto g_xxyy_xxyyyzz = cbuffer.data(gk_off + 125 * ccomps * dcomps);

            auto g_xxyy_xxyyzzz = cbuffer.data(gk_off + 126 * ccomps * dcomps);

            auto g_xxyy_xxyzzzz = cbuffer.data(gk_off + 127 * ccomps * dcomps);

            auto g_xxyy_xxzzzzz = cbuffer.data(gk_off + 128 * ccomps * dcomps);

            auto g_xxyy_xyyyyyy = cbuffer.data(gk_off + 129 * ccomps * dcomps);

            auto g_xxyy_xyyyyyz = cbuffer.data(gk_off + 130 * ccomps * dcomps);

            auto g_xxyy_xyyyyzz = cbuffer.data(gk_off + 131 * ccomps * dcomps);

            auto g_xxyy_xyyyzzz = cbuffer.data(gk_off + 132 * ccomps * dcomps);

            auto g_xxyy_xyyzzzz = cbuffer.data(gk_off + 133 * ccomps * dcomps);

            auto g_xxyy_xyzzzzz = cbuffer.data(gk_off + 134 * ccomps * dcomps);

            auto g_xxyy_xzzzzzz = cbuffer.data(gk_off + 135 * ccomps * dcomps);

            auto g_xxyy_yyyyyyy = cbuffer.data(gk_off + 136 * ccomps * dcomps);

            auto g_xxyy_yyyyyyz = cbuffer.data(gk_off + 137 * ccomps * dcomps);

            auto g_xxyy_yyyyyzz = cbuffer.data(gk_off + 138 * ccomps * dcomps);

            auto g_xxyy_yyyyzzz = cbuffer.data(gk_off + 139 * ccomps * dcomps);

            auto g_xxyy_yyyzzzz = cbuffer.data(gk_off + 140 * ccomps * dcomps);

            auto g_xxyy_yyzzzzz = cbuffer.data(gk_off + 141 * ccomps * dcomps);

            auto g_xxyy_yzzzzzz = cbuffer.data(gk_off + 142 * ccomps * dcomps);

            auto g_xxyy_zzzzzzz = cbuffer.data(gk_off + 143 * ccomps * dcomps);

            auto g_xxyz_xxxxxxx = cbuffer.data(gk_off + 144 * ccomps * dcomps);

            auto g_xxyz_xxxxxxy = cbuffer.data(gk_off + 145 * ccomps * dcomps);

            auto g_xxyz_xxxxxxz = cbuffer.data(gk_off + 146 * ccomps * dcomps);

            auto g_xxyz_xxxxxyy = cbuffer.data(gk_off + 147 * ccomps * dcomps);

            auto g_xxyz_xxxxxyz = cbuffer.data(gk_off + 148 * ccomps * dcomps);

            auto g_xxyz_xxxxxzz = cbuffer.data(gk_off + 149 * ccomps * dcomps);

            auto g_xxyz_xxxxyyy = cbuffer.data(gk_off + 150 * ccomps * dcomps);

            auto g_xxyz_xxxxyyz = cbuffer.data(gk_off + 151 * ccomps * dcomps);

            auto g_xxyz_xxxxyzz = cbuffer.data(gk_off + 152 * ccomps * dcomps);

            auto g_xxyz_xxxxzzz = cbuffer.data(gk_off + 153 * ccomps * dcomps);

            auto g_xxyz_xxxyyyy = cbuffer.data(gk_off + 154 * ccomps * dcomps);

            auto g_xxyz_xxxyyyz = cbuffer.data(gk_off + 155 * ccomps * dcomps);

            auto g_xxyz_xxxyyzz = cbuffer.data(gk_off + 156 * ccomps * dcomps);

            auto g_xxyz_xxxyzzz = cbuffer.data(gk_off + 157 * ccomps * dcomps);

            auto g_xxyz_xxxzzzz = cbuffer.data(gk_off + 158 * ccomps * dcomps);

            auto g_xxyz_xxyyyyy = cbuffer.data(gk_off + 159 * ccomps * dcomps);

            auto g_xxyz_xxyyyyz = cbuffer.data(gk_off + 160 * ccomps * dcomps);

            auto g_xxyz_xxyyyzz = cbuffer.data(gk_off + 161 * ccomps * dcomps);

            auto g_xxyz_xxyyzzz = cbuffer.data(gk_off + 162 * ccomps * dcomps);

            auto g_xxyz_xxyzzzz = cbuffer.data(gk_off + 163 * ccomps * dcomps);

            auto g_xxyz_xxzzzzz = cbuffer.data(gk_off + 164 * ccomps * dcomps);

            auto g_xxyz_xyyyyyy = cbuffer.data(gk_off + 165 * ccomps * dcomps);

            auto g_xxyz_xyyyyyz = cbuffer.data(gk_off + 166 * ccomps * dcomps);

            auto g_xxyz_xyyyyzz = cbuffer.data(gk_off + 167 * ccomps * dcomps);

            auto g_xxyz_xyyyzzz = cbuffer.data(gk_off + 168 * ccomps * dcomps);

            auto g_xxyz_xyyzzzz = cbuffer.data(gk_off + 169 * ccomps * dcomps);

            auto g_xxyz_xyzzzzz = cbuffer.data(gk_off + 170 * ccomps * dcomps);

            auto g_xxyz_xzzzzzz = cbuffer.data(gk_off + 171 * ccomps * dcomps);

            auto g_xxyz_yyyyyyy = cbuffer.data(gk_off + 172 * ccomps * dcomps);

            auto g_xxyz_yyyyyyz = cbuffer.data(gk_off + 173 * ccomps * dcomps);

            auto g_xxyz_yyyyyzz = cbuffer.data(gk_off + 174 * ccomps * dcomps);

            auto g_xxyz_yyyyzzz = cbuffer.data(gk_off + 175 * ccomps * dcomps);

            auto g_xxyz_yyyzzzz = cbuffer.data(gk_off + 176 * ccomps * dcomps);

            auto g_xxyz_yyzzzzz = cbuffer.data(gk_off + 177 * ccomps * dcomps);

            auto g_xxyz_yzzzzzz = cbuffer.data(gk_off + 178 * ccomps * dcomps);

            auto g_xxyz_zzzzzzz = cbuffer.data(gk_off + 179 * ccomps * dcomps);

            auto g_xxzz_xxxxxxx = cbuffer.data(gk_off + 180 * ccomps * dcomps);

            auto g_xxzz_xxxxxxy = cbuffer.data(gk_off + 181 * ccomps * dcomps);

            auto g_xxzz_xxxxxxz = cbuffer.data(gk_off + 182 * ccomps * dcomps);

            auto g_xxzz_xxxxxyy = cbuffer.data(gk_off + 183 * ccomps * dcomps);

            auto g_xxzz_xxxxxyz = cbuffer.data(gk_off + 184 * ccomps * dcomps);

            auto g_xxzz_xxxxxzz = cbuffer.data(gk_off + 185 * ccomps * dcomps);

            auto g_xxzz_xxxxyyy = cbuffer.data(gk_off + 186 * ccomps * dcomps);

            auto g_xxzz_xxxxyyz = cbuffer.data(gk_off + 187 * ccomps * dcomps);

            auto g_xxzz_xxxxyzz = cbuffer.data(gk_off + 188 * ccomps * dcomps);

            auto g_xxzz_xxxxzzz = cbuffer.data(gk_off + 189 * ccomps * dcomps);

            auto g_xxzz_xxxyyyy = cbuffer.data(gk_off + 190 * ccomps * dcomps);

            auto g_xxzz_xxxyyyz = cbuffer.data(gk_off + 191 * ccomps * dcomps);

            auto g_xxzz_xxxyyzz = cbuffer.data(gk_off + 192 * ccomps * dcomps);

            auto g_xxzz_xxxyzzz = cbuffer.data(gk_off + 193 * ccomps * dcomps);

            auto g_xxzz_xxxzzzz = cbuffer.data(gk_off + 194 * ccomps * dcomps);

            auto g_xxzz_xxyyyyy = cbuffer.data(gk_off + 195 * ccomps * dcomps);

            auto g_xxzz_xxyyyyz = cbuffer.data(gk_off + 196 * ccomps * dcomps);

            auto g_xxzz_xxyyyzz = cbuffer.data(gk_off + 197 * ccomps * dcomps);

            auto g_xxzz_xxyyzzz = cbuffer.data(gk_off + 198 * ccomps * dcomps);

            auto g_xxzz_xxyzzzz = cbuffer.data(gk_off + 199 * ccomps * dcomps);

            auto g_xxzz_xxzzzzz = cbuffer.data(gk_off + 200 * ccomps * dcomps);

            auto g_xxzz_xyyyyyy = cbuffer.data(gk_off + 201 * ccomps * dcomps);

            auto g_xxzz_xyyyyyz = cbuffer.data(gk_off + 202 * ccomps * dcomps);

            auto g_xxzz_xyyyyzz = cbuffer.data(gk_off + 203 * ccomps * dcomps);

            auto g_xxzz_xyyyzzz = cbuffer.data(gk_off + 204 * ccomps * dcomps);

            auto g_xxzz_xyyzzzz = cbuffer.data(gk_off + 205 * ccomps * dcomps);

            auto g_xxzz_xyzzzzz = cbuffer.data(gk_off + 206 * ccomps * dcomps);

            auto g_xxzz_xzzzzzz = cbuffer.data(gk_off + 207 * ccomps * dcomps);

            auto g_xxzz_yyyyyyy = cbuffer.data(gk_off + 208 * ccomps * dcomps);

            auto g_xxzz_yyyyyyz = cbuffer.data(gk_off + 209 * ccomps * dcomps);

            auto g_xxzz_yyyyyzz = cbuffer.data(gk_off + 210 * ccomps * dcomps);

            auto g_xxzz_yyyyzzz = cbuffer.data(gk_off + 211 * ccomps * dcomps);

            auto g_xxzz_yyyzzzz = cbuffer.data(gk_off + 212 * ccomps * dcomps);

            auto g_xxzz_yyzzzzz = cbuffer.data(gk_off + 213 * ccomps * dcomps);

            auto g_xxzz_yzzzzzz = cbuffer.data(gk_off + 214 * ccomps * dcomps);

            auto g_xxzz_zzzzzzz = cbuffer.data(gk_off + 215 * ccomps * dcomps);

            auto g_xyyy_xxxxxxx = cbuffer.data(gk_off + 216 * ccomps * dcomps);

            auto g_xyyy_xxxxxxy = cbuffer.data(gk_off + 217 * ccomps * dcomps);

            auto g_xyyy_xxxxxxz = cbuffer.data(gk_off + 218 * ccomps * dcomps);

            auto g_xyyy_xxxxxyy = cbuffer.data(gk_off + 219 * ccomps * dcomps);

            auto g_xyyy_xxxxxyz = cbuffer.data(gk_off + 220 * ccomps * dcomps);

            auto g_xyyy_xxxxxzz = cbuffer.data(gk_off + 221 * ccomps * dcomps);

            auto g_xyyy_xxxxyyy = cbuffer.data(gk_off + 222 * ccomps * dcomps);

            auto g_xyyy_xxxxyyz = cbuffer.data(gk_off + 223 * ccomps * dcomps);

            auto g_xyyy_xxxxyzz = cbuffer.data(gk_off + 224 * ccomps * dcomps);

            auto g_xyyy_xxxxzzz = cbuffer.data(gk_off + 225 * ccomps * dcomps);

            auto g_xyyy_xxxyyyy = cbuffer.data(gk_off + 226 * ccomps * dcomps);

            auto g_xyyy_xxxyyyz = cbuffer.data(gk_off + 227 * ccomps * dcomps);

            auto g_xyyy_xxxyyzz = cbuffer.data(gk_off + 228 * ccomps * dcomps);

            auto g_xyyy_xxxyzzz = cbuffer.data(gk_off + 229 * ccomps * dcomps);

            auto g_xyyy_xxxzzzz = cbuffer.data(gk_off + 230 * ccomps * dcomps);

            auto g_xyyy_xxyyyyy = cbuffer.data(gk_off + 231 * ccomps * dcomps);

            auto g_xyyy_xxyyyyz = cbuffer.data(gk_off + 232 * ccomps * dcomps);

            auto g_xyyy_xxyyyzz = cbuffer.data(gk_off + 233 * ccomps * dcomps);

            auto g_xyyy_xxyyzzz = cbuffer.data(gk_off + 234 * ccomps * dcomps);

            auto g_xyyy_xxyzzzz = cbuffer.data(gk_off + 235 * ccomps * dcomps);

            auto g_xyyy_xxzzzzz = cbuffer.data(gk_off + 236 * ccomps * dcomps);

            auto g_xyyy_xyyyyyy = cbuffer.data(gk_off + 237 * ccomps * dcomps);

            auto g_xyyy_xyyyyyz = cbuffer.data(gk_off + 238 * ccomps * dcomps);

            auto g_xyyy_xyyyyzz = cbuffer.data(gk_off + 239 * ccomps * dcomps);

            auto g_xyyy_xyyyzzz = cbuffer.data(gk_off + 240 * ccomps * dcomps);

            auto g_xyyy_xyyzzzz = cbuffer.data(gk_off + 241 * ccomps * dcomps);

            auto g_xyyy_xyzzzzz = cbuffer.data(gk_off + 242 * ccomps * dcomps);

            auto g_xyyy_xzzzzzz = cbuffer.data(gk_off + 243 * ccomps * dcomps);

            auto g_xyyy_yyyyyyy = cbuffer.data(gk_off + 244 * ccomps * dcomps);

            auto g_xyyy_yyyyyyz = cbuffer.data(gk_off + 245 * ccomps * dcomps);

            auto g_xyyy_yyyyyzz = cbuffer.data(gk_off + 246 * ccomps * dcomps);

            auto g_xyyy_yyyyzzz = cbuffer.data(gk_off + 247 * ccomps * dcomps);

            auto g_xyyy_yyyzzzz = cbuffer.data(gk_off + 248 * ccomps * dcomps);

            auto g_xyyy_yyzzzzz = cbuffer.data(gk_off + 249 * ccomps * dcomps);

            auto g_xyyy_yzzzzzz = cbuffer.data(gk_off + 250 * ccomps * dcomps);

            auto g_xyyy_zzzzzzz = cbuffer.data(gk_off + 251 * ccomps * dcomps);

            auto g_xyyz_xxxxxxx = cbuffer.data(gk_off + 252 * ccomps * dcomps);

            auto g_xyyz_xxxxxxy = cbuffer.data(gk_off + 253 * ccomps * dcomps);

            auto g_xyyz_xxxxxxz = cbuffer.data(gk_off + 254 * ccomps * dcomps);

            auto g_xyyz_xxxxxyy = cbuffer.data(gk_off + 255 * ccomps * dcomps);

            auto g_xyyz_xxxxxyz = cbuffer.data(gk_off + 256 * ccomps * dcomps);

            auto g_xyyz_xxxxxzz = cbuffer.data(gk_off + 257 * ccomps * dcomps);

            auto g_xyyz_xxxxyyy = cbuffer.data(gk_off + 258 * ccomps * dcomps);

            auto g_xyyz_xxxxyyz = cbuffer.data(gk_off + 259 * ccomps * dcomps);

            auto g_xyyz_xxxxyzz = cbuffer.data(gk_off + 260 * ccomps * dcomps);

            auto g_xyyz_xxxxzzz = cbuffer.data(gk_off + 261 * ccomps * dcomps);

            auto g_xyyz_xxxyyyy = cbuffer.data(gk_off + 262 * ccomps * dcomps);

            auto g_xyyz_xxxyyyz = cbuffer.data(gk_off + 263 * ccomps * dcomps);

            auto g_xyyz_xxxyyzz = cbuffer.data(gk_off + 264 * ccomps * dcomps);

            auto g_xyyz_xxxyzzz = cbuffer.data(gk_off + 265 * ccomps * dcomps);

            auto g_xyyz_xxxzzzz = cbuffer.data(gk_off + 266 * ccomps * dcomps);

            auto g_xyyz_xxyyyyy = cbuffer.data(gk_off + 267 * ccomps * dcomps);

            auto g_xyyz_xxyyyyz = cbuffer.data(gk_off + 268 * ccomps * dcomps);

            auto g_xyyz_xxyyyzz = cbuffer.data(gk_off + 269 * ccomps * dcomps);

            auto g_xyyz_xxyyzzz = cbuffer.data(gk_off + 270 * ccomps * dcomps);

            auto g_xyyz_xxyzzzz = cbuffer.data(gk_off + 271 * ccomps * dcomps);

            auto g_xyyz_xxzzzzz = cbuffer.data(gk_off + 272 * ccomps * dcomps);

            auto g_xyyz_xyyyyyy = cbuffer.data(gk_off + 273 * ccomps * dcomps);

            auto g_xyyz_xyyyyyz = cbuffer.data(gk_off + 274 * ccomps * dcomps);

            auto g_xyyz_xyyyyzz = cbuffer.data(gk_off + 275 * ccomps * dcomps);

            auto g_xyyz_xyyyzzz = cbuffer.data(gk_off + 276 * ccomps * dcomps);

            auto g_xyyz_xyyzzzz = cbuffer.data(gk_off + 277 * ccomps * dcomps);

            auto g_xyyz_xyzzzzz = cbuffer.data(gk_off + 278 * ccomps * dcomps);

            auto g_xyyz_xzzzzzz = cbuffer.data(gk_off + 279 * ccomps * dcomps);

            auto g_xyyz_yyyyyyy = cbuffer.data(gk_off + 280 * ccomps * dcomps);

            auto g_xyyz_yyyyyyz = cbuffer.data(gk_off + 281 * ccomps * dcomps);

            auto g_xyyz_yyyyyzz = cbuffer.data(gk_off + 282 * ccomps * dcomps);

            auto g_xyyz_yyyyzzz = cbuffer.data(gk_off + 283 * ccomps * dcomps);

            auto g_xyyz_yyyzzzz = cbuffer.data(gk_off + 284 * ccomps * dcomps);

            auto g_xyyz_yyzzzzz = cbuffer.data(gk_off + 285 * ccomps * dcomps);

            auto g_xyyz_yzzzzzz = cbuffer.data(gk_off + 286 * ccomps * dcomps);

            auto g_xyyz_zzzzzzz = cbuffer.data(gk_off + 287 * ccomps * dcomps);

            auto g_xyzz_xxxxxxx = cbuffer.data(gk_off + 288 * ccomps * dcomps);

            auto g_xyzz_xxxxxxy = cbuffer.data(gk_off + 289 * ccomps * dcomps);

            auto g_xyzz_xxxxxxz = cbuffer.data(gk_off + 290 * ccomps * dcomps);

            auto g_xyzz_xxxxxyy = cbuffer.data(gk_off + 291 * ccomps * dcomps);

            auto g_xyzz_xxxxxyz = cbuffer.data(gk_off + 292 * ccomps * dcomps);

            auto g_xyzz_xxxxxzz = cbuffer.data(gk_off + 293 * ccomps * dcomps);

            auto g_xyzz_xxxxyyy = cbuffer.data(gk_off + 294 * ccomps * dcomps);

            auto g_xyzz_xxxxyyz = cbuffer.data(gk_off + 295 * ccomps * dcomps);

            auto g_xyzz_xxxxyzz = cbuffer.data(gk_off + 296 * ccomps * dcomps);

            auto g_xyzz_xxxxzzz = cbuffer.data(gk_off + 297 * ccomps * dcomps);

            auto g_xyzz_xxxyyyy = cbuffer.data(gk_off + 298 * ccomps * dcomps);

            auto g_xyzz_xxxyyyz = cbuffer.data(gk_off + 299 * ccomps * dcomps);

            auto g_xyzz_xxxyyzz = cbuffer.data(gk_off + 300 * ccomps * dcomps);

            auto g_xyzz_xxxyzzz = cbuffer.data(gk_off + 301 * ccomps * dcomps);

            auto g_xyzz_xxxzzzz = cbuffer.data(gk_off + 302 * ccomps * dcomps);

            auto g_xyzz_xxyyyyy = cbuffer.data(gk_off + 303 * ccomps * dcomps);

            auto g_xyzz_xxyyyyz = cbuffer.data(gk_off + 304 * ccomps * dcomps);

            auto g_xyzz_xxyyyzz = cbuffer.data(gk_off + 305 * ccomps * dcomps);

            auto g_xyzz_xxyyzzz = cbuffer.data(gk_off + 306 * ccomps * dcomps);

            auto g_xyzz_xxyzzzz = cbuffer.data(gk_off + 307 * ccomps * dcomps);

            auto g_xyzz_xxzzzzz = cbuffer.data(gk_off + 308 * ccomps * dcomps);

            auto g_xyzz_xyyyyyy = cbuffer.data(gk_off + 309 * ccomps * dcomps);

            auto g_xyzz_xyyyyyz = cbuffer.data(gk_off + 310 * ccomps * dcomps);

            auto g_xyzz_xyyyyzz = cbuffer.data(gk_off + 311 * ccomps * dcomps);

            auto g_xyzz_xyyyzzz = cbuffer.data(gk_off + 312 * ccomps * dcomps);

            auto g_xyzz_xyyzzzz = cbuffer.data(gk_off + 313 * ccomps * dcomps);

            auto g_xyzz_xyzzzzz = cbuffer.data(gk_off + 314 * ccomps * dcomps);

            auto g_xyzz_xzzzzzz = cbuffer.data(gk_off + 315 * ccomps * dcomps);

            auto g_xyzz_yyyyyyy = cbuffer.data(gk_off + 316 * ccomps * dcomps);

            auto g_xyzz_yyyyyyz = cbuffer.data(gk_off + 317 * ccomps * dcomps);

            auto g_xyzz_yyyyyzz = cbuffer.data(gk_off + 318 * ccomps * dcomps);

            auto g_xyzz_yyyyzzz = cbuffer.data(gk_off + 319 * ccomps * dcomps);

            auto g_xyzz_yyyzzzz = cbuffer.data(gk_off + 320 * ccomps * dcomps);

            auto g_xyzz_yyzzzzz = cbuffer.data(gk_off + 321 * ccomps * dcomps);

            auto g_xyzz_yzzzzzz = cbuffer.data(gk_off + 322 * ccomps * dcomps);

            auto g_xyzz_zzzzzzz = cbuffer.data(gk_off + 323 * ccomps * dcomps);

            auto g_xzzz_xxxxxxx = cbuffer.data(gk_off + 324 * ccomps * dcomps);

            auto g_xzzz_xxxxxxy = cbuffer.data(gk_off + 325 * ccomps * dcomps);

            auto g_xzzz_xxxxxxz = cbuffer.data(gk_off + 326 * ccomps * dcomps);

            auto g_xzzz_xxxxxyy = cbuffer.data(gk_off + 327 * ccomps * dcomps);

            auto g_xzzz_xxxxxyz = cbuffer.data(gk_off + 328 * ccomps * dcomps);

            auto g_xzzz_xxxxxzz = cbuffer.data(gk_off + 329 * ccomps * dcomps);

            auto g_xzzz_xxxxyyy = cbuffer.data(gk_off + 330 * ccomps * dcomps);

            auto g_xzzz_xxxxyyz = cbuffer.data(gk_off + 331 * ccomps * dcomps);

            auto g_xzzz_xxxxyzz = cbuffer.data(gk_off + 332 * ccomps * dcomps);

            auto g_xzzz_xxxxzzz = cbuffer.data(gk_off + 333 * ccomps * dcomps);

            auto g_xzzz_xxxyyyy = cbuffer.data(gk_off + 334 * ccomps * dcomps);

            auto g_xzzz_xxxyyyz = cbuffer.data(gk_off + 335 * ccomps * dcomps);

            auto g_xzzz_xxxyyzz = cbuffer.data(gk_off + 336 * ccomps * dcomps);

            auto g_xzzz_xxxyzzz = cbuffer.data(gk_off + 337 * ccomps * dcomps);

            auto g_xzzz_xxxzzzz = cbuffer.data(gk_off + 338 * ccomps * dcomps);

            auto g_xzzz_xxyyyyy = cbuffer.data(gk_off + 339 * ccomps * dcomps);

            auto g_xzzz_xxyyyyz = cbuffer.data(gk_off + 340 * ccomps * dcomps);

            auto g_xzzz_xxyyyzz = cbuffer.data(gk_off + 341 * ccomps * dcomps);

            auto g_xzzz_xxyyzzz = cbuffer.data(gk_off + 342 * ccomps * dcomps);

            auto g_xzzz_xxyzzzz = cbuffer.data(gk_off + 343 * ccomps * dcomps);

            auto g_xzzz_xxzzzzz = cbuffer.data(gk_off + 344 * ccomps * dcomps);

            auto g_xzzz_xyyyyyy = cbuffer.data(gk_off + 345 * ccomps * dcomps);

            auto g_xzzz_xyyyyyz = cbuffer.data(gk_off + 346 * ccomps * dcomps);

            auto g_xzzz_xyyyyzz = cbuffer.data(gk_off + 347 * ccomps * dcomps);

            auto g_xzzz_xyyyzzz = cbuffer.data(gk_off + 348 * ccomps * dcomps);

            auto g_xzzz_xyyzzzz = cbuffer.data(gk_off + 349 * ccomps * dcomps);

            auto g_xzzz_xyzzzzz = cbuffer.data(gk_off + 350 * ccomps * dcomps);

            auto g_xzzz_xzzzzzz = cbuffer.data(gk_off + 351 * ccomps * dcomps);

            auto g_xzzz_yyyyyyy = cbuffer.data(gk_off + 352 * ccomps * dcomps);

            auto g_xzzz_yyyyyyz = cbuffer.data(gk_off + 353 * ccomps * dcomps);

            auto g_xzzz_yyyyyzz = cbuffer.data(gk_off + 354 * ccomps * dcomps);

            auto g_xzzz_yyyyzzz = cbuffer.data(gk_off + 355 * ccomps * dcomps);

            auto g_xzzz_yyyzzzz = cbuffer.data(gk_off + 356 * ccomps * dcomps);

            auto g_xzzz_yyzzzzz = cbuffer.data(gk_off + 357 * ccomps * dcomps);

            auto g_xzzz_yzzzzzz = cbuffer.data(gk_off + 358 * ccomps * dcomps);

            auto g_xzzz_zzzzzzz = cbuffer.data(gk_off + 359 * ccomps * dcomps);

            auto g_yyyy_xxxxxxx = cbuffer.data(gk_off + 360 * ccomps * dcomps);

            auto g_yyyy_xxxxxxy = cbuffer.data(gk_off + 361 * ccomps * dcomps);

            auto g_yyyy_xxxxxxz = cbuffer.data(gk_off + 362 * ccomps * dcomps);

            auto g_yyyy_xxxxxyy = cbuffer.data(gk_off + 363 * ccomps * dcomps);

            auto g_yyyy_xxxxxyz = cbuffer.data(gk_off + 364 * ccomps * dcomps);

            auto g_yyyy_xxxxxzz = cbuffer.data(gk_off + 365 * ccomps * dcomps);

            auto g_yyyy_xxxxyyy = cbuffer.data(gk_off + 366 * ccomps * dcomps);

            auto g_yyyy_xxxxyyz = cbuffer.data(gk_off + 367 * ccomps * dcomps);

            auto g_yyyy_xxxxyzz = cbuffer.data(gk_off + 368 * ccomps * dcomps);

            auto g_yyyy_xxxxzzz = cbuffer.data(gk_off + 369 * ccomps * dcomps);

            auto g_yyyy_xxxyyyy = cbuffer.data(gk_off + 370 * ccomps * dcomps);

            auto g_yyyy_xxxyyyz = cbuffer.data(gk_off + 371 * ccomps * dcomps);

            auto g_yyyy_xxxyyzz = cbuffer.data(gk_off + 372 * ccomps * dcomps);

            auto g_yyyy_xxxyzzz = cbuffer.data(gk_off + 373 * ccomps * dcomps);

            auto g_yyyy_xxxzzzz = cbuffer.data(gk_off + 374 * ccomps * dcomps);

            auto g_yyyy_xxyyyyy = cbuffer.data(gk_off + 375 * ccomps * dcomps);

            auto g_yyyy_xxyyyyz = cbuffer.data(gk_off + 376 * ccomps * dcomps);

            auto g_yyyy_xxyyyzz = cbuffer.data(gk_off + 377 * ccomps * dcomps);

            auto g_yyyy_xxyyzzz = cbuffer.data(gk_off + 378 * ccomps * dcomps);

            auto g_yyyy_xxyzzzz = cbuffer.data(gk_off + 379 * ccomps * dcomps);

            auto g_yyyy_xxzzzzz = cbuffer.data(gk_off + 380 * ccomps * dcomps);

            auto g_yyyy_xyyyyyy = cbuffer.data(gk_off + 381 * ccomps * dcomps);

            auto g_yyyy_xyyyyyz = cbuffer.data(gk_off + 382 * ccomps * dcomps);

            auto g_yyyy_xyyyyzz = cbuffer.data(gk_off + 383 * ccomps * dcomps);

            auto g_yyyy_xyyyzzz = cbuffer.data(gk_off + 384 * ccomps * dcomps);

            auto g_yyyy_xyyzzzz = cbuffer.data(gk_off + 385 * ccomps * dcomps);

            auto g_yyyy_xyzzzzz = cbuffer.data(gk_off + 386 * ccomps * dcomps);

            auto g_yyyy_xzzzzzz = cbuffer.data(gk_off + 387 * ccomps * dcomps);

            auto g_yyyy_yyyyyyy = cbuffer.data(gk_off + 388 * ccomps * dcomps);

            auto g_yyyy_yyyyyyz = cbuffer.data(gk_off + 389 * ccomps * dcomps);

            auto g_yyyy_yyyyyzz = cbuffer.data(gk_off + 390 * ccomps * dcomps);

            auto g_yyyy_yyyyzzz = cbuffer.data(gk_off + 391 * ccomps * dcomps);

            auto g_yyyy_yyyzzzz = cbuffer.data(gk_off + 392 * ccomps * dcomps);

            auto g_yyyy_yyzzzzz = cbuffer.data(gk_off + 393 * ccomps * dcomps);

            auto g_yyyy_yzzzzzz = cbuffer.data(gk_off + 394 * ccomps * dcomps);

            auto g_yyyy_zzzzzzz = cbuffer.data(gk_off + 395 * ccomps * dcomps);

            auto g_yyyz_xxxxxxx = cbuffer.data(gk_off + 396 * ccomps * dcomps);

            auto g_yyyz_xxxxxxy = cbuffer.data(gk_off + 397 * ccomps * dcomps);

            auto g_yyyz_xxxxxxz = cbuffer.data(gk_off + 398 * ccomps * dcomps);

            auto g_yyyz_xxxxxyy = cbuffer.data(gk_off + 399 * ccomps * dcomps);

            auto g_yyyz_xxxxxyz = cbuffer.data(gk_off + 400 * ccomps * dcomps);

            auto g_yyyz_xxxxxzz = cbuffer.data(gk_off + 401 * ccomps * dcomps);

            auto g_yyyz_xxxxyyy = cbuffer.data(gk_off + 402 * ccomps * dcomps);

            auto g_yyyz_xxxxyyz = cbuffer.data(gk_off + 403 * ccomps * dcomps);

            auto g_yyyz_xxxxyzz = cbuffer.data(gk_off + 404 * ccomps * dcomps);

            auto g_yyyz_xxxxzzz = cbuffer.data(gk_off + 405 * ccomps * dcomps);

            auto g_yyyz_xxxyyyy = cbuffer.data(gk_off + 406 * ccomps * dcomps);

            auto g_yyyz_xxxyyyz = cbuffer.data(gk_off + 407 * ccomps * dcomps);

            auto g_yyyz_xxxyyzz = cbuffer.data(gk_off + 408 * ccomps * dcomps);

            auto g_yyyz_xxxyzzz = cbuffer.data(gk_off + 409 * ccomps * dcomps);

            auto g_yyyz_xxxzzzz = cbuffer.data(gk_off + 410 * ccomps * dcomps);

            auto g_yyyz_xxyyyyy = cbuffer.data(gk_off + 411 * ccomps * dcomps);

            auto g_yyyz_xxyyyyz = cbuffer.data(gk_off + 412 * ccomps * dcomps);

            auto g_yyyz_xxyyyzz = cbuffer.data(gk_off + 413 * ccomps * dcomps);

            auto g_yyyz_xxyyzzz = cbuffer.data(gk_off + 414 * ccomps * dcomps);

            auto g_yyyz_xxyzzzz = cbuffer.data(gk_off + 415 * ccomps * dcomps);

            auto g_yyyz_xxzzzzz = cbuffer.data(gk_off + 416 * ccomps * dcomps);

            auto g_yyyz_xyyyyyy = cbuffer.data(gk_off + 417 * ccomps * dcomps);

            auto g_yyyz_xyyyyyz = cbuffer.data(gk_off + 418 * ccomps * dcomps);

            auto g_yyyz_xyyyyzz = cbuffer.data(gk_off + 419 * ccomps * dcomps);

            auto g_yyyz_xyyyzzz = cbuffer.data(gk_off + 420 * ccomps * dcomps);

            auto g_yyyz_xyyzzzz = cbuffer.data(gk_off + 421 * ccomps * dcomps);

            auto g_yyyz_xyzzzzz = cbuffer.data(gk_off + 422 * ccomps * dcomps);

            auto g_yyyz_xzzzzzz = cbuffer.data(gk_off + 423 * ccomps * dcomps);

            auto g_yyyz_yyyyyyy = cbuffer.data(gk_off + 424 * ccomps * dcomps);

            auto g_yyyz_yyyyyyz = cbuffer.data(gk_off + 425 * ccomps * dcomps);

            auto g_yyyz_yyyyyzz = cbuffer.data(gk_off + 426 * ccomps * dcomps);

            auto g_yyyz_yyyyzzz = cbuffer.data(gk_off + 427 * ccomps * dcomps);

            auto g_yyyz_yyyzzzz = cbuffer.data(gk_off + 428 * ccomps * dcomps);

            auto g_yyyz_yyzzzzz = cbuffer.data(gk_off + 429 * ccomps * dcomps);

            auto g_yyyz_yzzzzzz = cbuffer.data(gk_off + 430 * ccomps * dcomps);

            auto g_yyyz_zzzzzzz = cbuffer.data(gk_off + 431 * ccomps * dcomps);

            auto g_yyzz_xxxxxxx = cbuffer.data(gk_off + 432 * ccomps * dcomps);

            auto g_yyzz_xxxxxxy = cbuffer.data(gk_off + 433 * ccomps * dcomps);

            auto g_yyzz_xxxxxxz = cbuffer.data(gk_off + 434 * ccomps * dcomps);

            auto g_yyzz_xxxxxyy = cbuffer.data(gk_off + 435 * ccomps * dcomps);

            auto g_yyzz_xxxxxyz = cbuffer.data(gk_off + 436 * ccomps * dcomps);

            auto g_yyzz_xxxxxzz = cbuffer.data(gk_off + 437 * ccomps * dcomps);

            auto g_yyzz_xxxxyyy = cbuffer.data(gk_off + 438 * ccomps * dcomps);

            auto g_yyzz_xxxxyyz = cbuffer.data(gk_off + 439 * ccomps * dcomps);

            auto g_yyzz_xxxxyzz = cbuffer.data(gk_off + 440 * ccomps * dcomps);

            auto g_yyzz_xxxxzzz = cbuffer.data(gk_off + 441 * ccomps * dcomps);

            auto g_yyzz_xxxyyyy = cbuffer.data(gk_off + 442 * ccomps * dcomps);

            auto g_yyzz_xxxyyyz = cbuffer.data(gk_off + 443 * ccomps * dcomps);

            auto g_yyzz_xxxyyzz = cbuffer.data(gk_off + 444 * ccomps * dcomps);

            auto g_yyzz_xxxyzzz = cbuffer.data(gk_off + 445 * ccomps * dcomps);

            auto g_yyzz_xxxzzzz = cbuffer.data(gk_off + 446 * ccomps * dcomps);

            auto g_yyzz_xxyyyyy = cbuffer.data(gk_off + 447 * ccomps * dcomps);

            auto g_yyzz_xxyyyyz = cbuffer.data(gk_off + 448 * ccomps * dcomps);

            auto g_yyzz_xxyyyzz = cbuffer.data(gk_off + 449 * ccomps * dcomps);

            auto g_yyzz_xxyyzzz = cbuffer.data(gk_off + 450 * ccomps * dcomps);

            auto g_yyzz_xxyzzzz = cbuffer.data(gk_off + 451 * ccomps * dcomps);

            auto g_yyzz_xxzzzzz = cbuffer.data(gk_off + 452 * ccomps * dcomps);

            auto g_yyzz_xyyyyyy = cbuffer.data(gk_off + 453 * ccomps * dcomps);

            auto g_yyzz_xyyyyyz = cbuffer.data(gk_off + 454 * ccomps * dcomps);

            auto g_yyzz_xyyyyzz = cbuffer.data(gk_off + 455 * ccomps * dcomps);

            auto g_yyzz_xyyyzzz = cbuffer.data(gk_off + 456 * ccomps * dcomps);

            auto g_yyzz_xyyzzzz = cbuffer.data(gk_off + 457 * ccomps * dcomps);

            auto g_yyzz_xyzzzzz = cbuffer.data(gk_off + 458 * ccomps * dcomps);

            auto g_yyzz_xzzzzzz = cbuffer.data(gk_off + 459 * ccomps * dcomps);

            auto g_yyzz_yyyyyyy = cbuffer.data(gk_off + 460 * ccomps * dcomps);

            auto g_yyzz_yyyyyyz = cbuffer.data(gk_off + 461 * ccomps * dcomps);

            auto g_yyzz_yyyyyzz = cbuffer.data(gk_off + 462 * ccomps * dcomps);

            auto g_yyzz_yyyyzzz = cbuffer.data(gk_off + 463 * ccomps * dcomps);

            auto g_yyzz_yyyzzzz = cbuffer.data(gk_off + 464 * ccomps * dcomps);

            auto g_yyzz_yyzzzzz = cbuffer.data(gk_off + 465 * ccomps * dcomps);

            auto g_yyzz_yzzzzzz = cbuffer.data(gk_off + 466 * ccomps * dcomps);

            auto g_yyzz_zzzzzzz = cbuffer.data(gk_off + 467 * ccomps * dcomps);

            auto g_yzzz_xxxxxxx = cbuffer.data(gk_off + 468 * ccomps * dcomps);

            auto g_yzzz_xxxxxxy = cbuffer.data(gk_off + 469 * ccomps * dcomps);

            auto g_yzzz_xxxxxxz = cbuffer.data(gk_off + 470 * ccomps * dcomps);

            auto g_yzzz_xxxxxyy = cbuffer.data(gk_off + 471 * ccomps * dcomps);

            auto g_yzzz_xxxxxyz = cbuffer.data(gk_off + 472 * ccomps * dcomps);

            auto g_yzzz_xxxxxzz = cbuffer.data(gk_off + 473 * ccomps * dcomps);

            auto g_yzzz_xxxxyyy = cbuffer.data(gk_off + 474 * ccomps * dcomps);

            auto g_yzzz_xxxxyyz = cbuffer.data(gk_off + 475 * ccomps * dcomps);

            auto g_yzzz_xxxxyzz = cbuffer.data(gk_off + 476 * ccomps * dcomps);

            auto g_yzzz_xxxxzzz = cbuffer.data(gk_off + 477 * ccomps * dcomps);

            auto g_yzzz_xxxyyyy = cbuffer.data(gk_off + 478 * ccomps * dcomps);

            auto g_yzzz_xxxyyyz = cbuffer.data(gk_off + 479 * ccomps * dcomps);

            auto g_yzzz_xxxyyzz = cbuffer.data(gk_off + 480 * ccomps * dcomps);

            auto g_yzzz_xxxyzzz = cbuffer.data(gk_off + 481 * ccomps * dcomps);

            auto g_yzzz_xxxzzzz = cbuffer.data(gk_off + 482 * ccomps * dcomps);

            auto g_yzzz_xxyyyyy = cbuffer.data(gk_off + 483 * ccomps * dcomps);

            auto g_yzzz_xxyyyyz = cbuffer.data(gk_off + 484 * ccomps * dcomps);

            auto g_yzzz_xxyyyzz = cbuffer.data(gk_off + 485 * ccomps * dcomps);

            auto g_yzzz_xxyyzzz = cbuffer.data(gk_off + 486 * ccomps * dcomps);

            auto g_yzzz_xxyzzzz = cbuffer.data(gk_off + 487 * ccomps * dcomps);

            auto g_yzzz_xxzzzzz = cbuffer.data(gk_off + 488 * ccomps * dcomps);

            auto g_yzzz_xyyyyyy = cbuffer.data(gk_off + 489 * ccomps * dcomps);

            auto g_yzzz_xyyyyyz = cbuffer.data(gk_off + 490 * ccomps * dcomps);

            auto g_yzzz_xyyyyzz = cbuffer.data(gk_off + 491 * ccomps * dcomps);

            auto g_yzzz_xyyyzzz = cbuffer.data(gk_off + 492 * ccomps * dcomps);

            auto g_yzzz_xyyzzzz = cbuffer.data(gk_off + 493 * ccomps * dcomps);

            auto g_yzzz_xyzzzzz = cbuffer.data(gk_off + 494 * ccomps * dcomps);

            auto g_yzzz_xzzzzzz = cbuffer.data(gk_off + 495 * ccomps * dcomps);

            auto g_yzzz_yyyyyyy = cbuffer.data(gk_off + 496 * ccomps * dcomps);

            auto g_yzzz_yyyyyyz = cbuffer.data(gk_off + 497 * ccomps * dcomps);

            auto g_yzzz_yyyyyzz = cbuffer.data(gk_off + 498 * ccomps * dcomps);

            auto g_yzzz_yyyyzzz = cbuffer.data(gk_off + 499 * ccomps * dcomps);

            auto g_yzzz_yyyzzzz = cbuffer.data(gk_off + 500 * ccomps * dcomps);

            auto g_yzzz_yyzzzzz = cbuffer.data(gk_off + 501 * ccomps * dcomps);

            auto g_yzzz_yzzzzzz = cbuffer.data(gk_off + 502 * ccomps * dcomps);

            auto g_yzzz_zzzzzzz = cbuffer.data(gk_off + 503 * ccomps * dcomps);

            auto g_zzzz_xxxxxxx = cbuffer.data(gk_off + 504 * ccomps * dcomps);

            auto g_zzzz_xxxxxxy = cbuffer.data(gk_off + 505 * ccomps * dcomps);

            auto g_zzzz_xxxxxxz = cbuffer.data(gk_off + 506 * ccomps * dcomps);

            auto g_zzzz_xxxxxyy = cbuffer.data(gk_off + 507 * ccomps * dcomps);

            auto g_zzzz_xxxxxyz = cbuffer.data(gk_off + 508 * ccomps * dcomps);

            auto g_zzzz_xxxxxzz = cbuffer.data(gk_off + 509 * ccomps * dcomps);

            auto g_zzzz_xxxxyyy = cbuffer.data(gk_off + 510 * ccomps * dcomps);

            auto g_zzzz_xxxxyyz = cbuffer.data(gk_off + 511 * ccomps * dcomps);

            auto g_zzzz_xxxxyzz = cbuffer.data(gk_off + 512 * ccomps * dcomps);

            auto g_zzzz_xxxxzzz = cbuffer.data(gk_off + 513 * ccomps * dcomps);

            auto g_zzzz_xxxyyyy = cbuffer.data(gk_off + 514 * ccomps * dcomps);

            auto g_zzzz_xxxyyyz = cbuffer.data(gk_off + 515 * ccomps * dcomps);

            auto g_zzzz_xxxyyzz = cbuffer.data(gk_off + 516 * ccomps * dcomps);

            auto g_zzzz_xxxyzzz = cbuffer.data(gk_off + 517 * ccomps * dcomps);

            auto g_zzzz_xxxzzzz = cbuffer.data(gk_off + 518 * ccomps * dcomps);

            auto g_zzzz_xxyyyyy = cbuffer.data(gk_off + 519 * ccomps * dcomps);

            auto g_zzzz_xxyyyyz = cbuffer.data(gk_off + 520 * ccomps * dcomps);

            auto g_zzzz_xxyyyzz = cbuffer.data(gk_off + 521 * ccomps * dcomps);

            auto g_zzzz_xxyyzzz = cbuffer.data(gk_off + 522 * ccomps * dcomps);

            auto g_zzzz_xxyzzzz = cbuffer.data(gk_off + 523 * ccomps * dcomps);

            auto g_zzzz_xxzzzzz = cbuffer.data(gk_off + 524 * ccomps * dcomps);

            auto g_zzzz_xyyyyyy = cbuffer.data(gk_off + 525 * ccomps * dcomps);

            auto g_zzzz_xyyyyyz = cbuffer.data(gk_off + 526 * ccomps * dcomps);

            auto g_zzzz_xyyyyzz = cbuffer.data(gk_off + 527 * ccomps * dcomps);

            auto g_zzzz_xyyyzzz = cbuffer.data(gk_off + 528 * ccomps * dcomps);

            auto g_zzzz_xyyzzzz = cbuffer.data(gk_off + 529 * ccomps * dcomps);

            auto g_zzzz_xyzzzzz = cbuffer.data(gk_off + 530 * ccomps * dcomps);

            auto g_zzzz_xzzzzzz = cbuffer.data(gk_off + 531 * ccomps * dcomps);

            auto g_zzzz_yyyyyyy = cbuffer.data(gk_off + 532 * ccomps * dcomps);

            auto g_zzzz_yyyyyyz = cbuffer.data(gk_off + 533 * ccomps * dcomps);

            auto g_zzzz_yyyyyzz = cbuffer.data(gk_off + 534 * ccomps * dcomps);

            auto g_zzzz_yyyyzzz = cbuffer.data(gk_off + 535 * ccomps * dcomps);

            auto g_zzzz_yyyzzzz = cbuffer.data(gk_off + 536 * ccomps * dcomps);

            auto g_zzzz_yyzzzzz = cbuffer.data(gk_off + 537 * ccomps * dcomps);

            auto g_zzzz_yzzzzzz = cbuffer.data(gk_off + 538 * ccomps * dcomps);

            auto g_zzzz_zzzzzzz = cbuffer.data(gk_off + 539 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GKSS

            const auto gk_geom_01_off = idx_geom_01_gkxx + i * dcomps + j;

            auto g_0_x_xxxx_xxxxxxx = cbuffer.data(gk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxy = cbuffer.data(gk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxz = cbuffer.data(gk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxyy = cbuffer.data(gk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxyz = cbuffer.data(gk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxzz = cbuffer.data(gk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyyy = cbuffer.data(gk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyyz = cbuffer.data(gk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyzz = cbuffer.data(gk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxzzz = cbuffer.data(gk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyyy = cbuffer.data(gk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyyz = cbuffer.data(gk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyzz = cbuffer.data(gk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyzzz = cbuffer.data(gk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxzzzz = cbuffer.data(gk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyyy = cbuffer.data(gk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyyz = cbuffer.data(gk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyzz = cbuffer.data(gk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyzzz = cbuffer.data(gk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyzzzz = cbuffer.data(gk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxx_xxzzzzz = cbuffer.data(gk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyyy = cbuffer.data(gk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyyz = cbuffer.data(gk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyzz = cbuffer.data(gk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyzzz = cbuffer.data(gk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyzzzz = cbuffer.data(gk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxx_xyzzzzz = cbuffer.data(gk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxx_xzzzzzz = cbuffer.data(gk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyyy = cbuffer.data(gk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyyz = cbuffer.data(gk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyzz = cbuffer.data(gk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyzzz = cbuffer.data(gk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyzzzz = cbuffer.data(gk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxx_yyzzzzz = cbuffer.data(gk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxx_yzzzzzz = cbuffer.data(gk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxx_zzzzzzz = cbuffer.data(gk_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxx = cbuffer.data(gk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxy = cbuffer.data(gk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxz = cbuffer.data(gk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxyy = cbuffer.data(gk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxyz = cbuffer.data(gk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxzz = cbuffer.data(gk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyyy = cbuffer.data(gk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyyz = cbuffer.data(gk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyzz = cbuffer.data(gk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxzzz = cbuffer.data(gk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyyy = cbuffer.data(gk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyyz = cbuffer.data(gk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyzz = cbuffer.data(gk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyzzz = cbuffer.data(gk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxzzzz = cbuffer.data(gk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyyy = cbuffer.data(gk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyyz = cbuffer.data(gk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyzz = cbuffer.data(gk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyzzz = cbuffer.data(gk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyzzzz = cbuffer.data(gk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxy_xxzzzzz = cbuffer.data(gk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyyy = cbuffer.data(gk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyyz = cbuffer.data(gk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyzz = cbuffer.data(gk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyzzz = cbuffer.data(gk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyzzzz = cbuffer.data(gk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxy_xyzzzzz = cbuffer.data(gk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxy_xzzzzzz = cbuffer.data(gk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyyy = cbuffer.data(gk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyyz = cbuffer.data(gk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyzz = cbuffer.data(gk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyzzz = cbuffer.data(gk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyzzzz = cbuffer.data(gk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxy_yyzzzzz = cbuffer.data(gk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxy_yzzzzzz = cbuffer.data(gk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxy_zzzzzzz = cbuffer.data(gk_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxx = cbuffer.data(gk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxy = cbuffer.data(gk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxz = cbuffer.data(gk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxyy = cbuffer.data(gk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxyz = cbuffer.data(gk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxzz = cbuffer.data(gk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyyy = cbuffer.data(gk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyyz = cbuffer.data(gk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyzz = cbuffer.data(gk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxzzz = cbuffer.data(gk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyyy = cbuffer.data(gk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyyz = cbuffer.data(gk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyzz = cbuffer.data(gk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyzzz = cbuffer.data(gk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxzzzz = cbuffer.data(gk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyyy = cbuffer.data(gk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyyz = cbuffer.data(gk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyzz = cbuffer.data(gk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyzzz = cbuffer.data(gk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyzzzz = cbuffer.data(gk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxz_xxzzzzz = cbuffer.data(gk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyyy = cbuffer.data(gk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyyz = cbuffer.data(gk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyzz = cbuffer.data(gk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyzzz = cbuffer.data(gk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyzzzz = cbuffer.data(gk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxz_xyzzzzz = cbuffer.data(gk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxz_xzzzzzz = cbuffer.data(gk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyyy = cbuffer.data(gk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyyz = cbuffer.data(gk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyzz = cbuffer.data(gk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyzzz = cbuffer.data(gk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyzzzz = cbuffer.data(gk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxxz_yyzzzzz = cbuffer.data(gk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxz_yzzzzzz = cbuffer.data(gk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxz_zzzzzzz = cbuffer.data(gk_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_xyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_xyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_xyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_xyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_xyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_xzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_xzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_xzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_xzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_xzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_xzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_yyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_yyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_yyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_yyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_yyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_yyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_yyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_yyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_yyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_yyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_yyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_yyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_x_yyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_x_yyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_x_yyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_x_yyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_x_yyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_x_yyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_x_yzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_x_yzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_x_yzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_x_yzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_x_yzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_x_yzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_x_zzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_x_zzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_x_zzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_x_zzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_x_zzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_x_zzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxx = cbuffer.data(gk_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxy = cbuffer.data(gk_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxz = cbuffer.data(gk_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxyy = cbuffer.data(gk_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxyz = cbuffer.data(gk_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxzz = cbuffer.data(gk_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyyy = cbuffer.data(gk_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyyz = cbuffer.data(gk_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyzz = cbuffer.data(gk_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxzzz = cbuffer.data(gk_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyyy = cbuffer.data(gk_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyyz = cbuffer.data(gk_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyzz = cbuffer.data(gk_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyzzz = cbuffer.data(gk_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxzzzz = cbuffer.data(gk_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyyy = cbuffer.data(gk_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyyz = cbuffer.data(gk_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyzz = cbuffer.data(gk_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyzzz = cbuffer.data(gk_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyzzzz = cbuffer.data(gk_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxxx_xxzzzzz = cbuffer.data(gk_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyyy = cbuffer.data(gk_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyyz = cbuffer.data(gk_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyzz = cbuffer.data(gk_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyzzz = cbuffer.data(gk_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyzzzz = cbuffer.data(gk_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxxx_xyzzzzz = cbuffer.data(gk_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxxx_xzzzzzz = cbuffer.data(gk_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyyy = cbuffer.data(gk_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyyz = cbuffer.data(gk_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyzz = cbuffer.data(gk_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyzzz = cbuffer.data(gk_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyzzzz = cbuffer.data(gk_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxxx_yyzzzzz = cbuffer.data(gk_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxxx_yzzzzzz = cbuffer.data(gk_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxxx_zzzzzzz = cbuffer.data(gk_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxx = cbuffer.data(gk_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxy = cbuffer.data(gk_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxz = cbuffer.data(gk_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxyy = cbuffer.data(gk_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxyz = cbuffer.data(gk_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxzz = cbuffer.data(gk_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyyy = cbuffer.data(gk_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyyz = cbuffer.data(gk_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyzz = cbuffer.data(gk_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxzzz = cbuffer.data(gk_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyyy = cbuffer.data(gk_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyyz = cbuffer.data(gk_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyzz = cbuffer.data(gk_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyzzz = cbuffer.data(gk_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxzzzz = cbuffer.data(gk_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyyy = cbuffer.data(gk_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyyz = cbuffer.data(gk_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyzz = cbuffer.data(gk_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyzzz = cbuffer.data(gk_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyzzzz = cbuffer.data(gk_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xxxy_xxzzzzz = cbuffer.data(gk_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyyy = cbuffer.data(gk_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyyz = cbuffer.data(gk_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyzz = cbuffer.data(gk_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyzzz = cbuffer.data(gk_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyzzzz = cbuffer.data(gk_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xxxy_xyzzzzz = cbuffer.data(gk_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xxxy_xzzzzzz = cbuffer.data(gk_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyyy = cbuffer.data(gk_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyyz = cbuffer.data(gk_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyzz = cbuffer.data(gk_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyzzz = cbuffer.data(gk_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyzzzz = cbuffer.data(gk_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xxxy_yyzzzzz = cbuffer.data(gk_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xxxy_yzzzzzz = cbuffer.data(gk_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xxxy_zzzzzzz = cbuffer.data(gk_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxx = cbuffer.data(gk_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxy = cbuffer.data(gk_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxz = cbuffer.data(gk_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxyy = cbuffer.data(gk_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxyz = cbuffer.data(gk_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxzz = cbuffer.data(gk_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyyy = cbuffer.data(gk_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyyz = cbuffer.data(gk_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyzz = cbuffer.data(gk_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxzzz = cbuffer.data(gk_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyyy = cbuffer.data(gk_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyyz = cbuffer.data(gk_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyzz = cbuffer.data(gk_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyzzz = cbuffer.data(gk_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxzzzz = cbuffer.data(gk_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyyy = cbuffer.data(gk_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyyz = cbuffer.data(gk_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyzz = cbuffer.data(gk_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyzzz = cbuffer.data(gk_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyzzzz = cbuffer.data(gk_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xxxz_xxzzzzz = cbuffer.data(gk_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyyy = cbuffer.data(gk_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyyz = cbuffer.data(gk_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyzz = cbuffer.data(gk_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyzzz = cbuffer.data(gk_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyzzzz = cbuffer.data(gk_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xxxz_xyzzzzz = cbuffer.data(gk_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xxxz_xzzzzzz = cbuffer.data(gk_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyyy = cbuffer.data(gk_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyyz = cbuffer.data(gk_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyzz = cbuffer.data(gk_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyzzz = cbuffer.data(gk_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyzzzz = cbuffer.data(gk_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_xxxz_yyzzzzz = cbuffer.data(gk_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xxxz_yzzzzzz = cbuffer.data(gk_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xxxz_zzzzzzz = cbuffer.data(gk_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xxyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xxyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xxyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xxyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xxyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xxyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xxyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xxyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xxyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xxyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xxyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xxyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_xxzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_xxzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_xxzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_xxzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_xxzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_xxzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_xyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_xyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_xyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_xyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_xyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_xyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_xyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_xyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_xyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_xyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_xyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_xyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_xyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_xyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_xyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_xyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_xyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_xyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_xzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_xzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_xzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_xzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_xzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_xzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_y_yyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_y_yyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_y_yyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_y_yyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_y_yyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_y_yyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_y_yyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_y_yyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_y_yyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_y_yyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_y_yyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_y_yyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_y_yyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_y_yyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_y_yyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_y_yyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_y_yyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_y_yyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_y_yzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_y_yzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_y_yzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_y_yzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_y_yzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_y_yzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_y_zzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_y_zzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_y_zzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_y_zzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_y_zzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_y_zzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxx = cbuffer.data(gk_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxy = cbuffer.data(gk_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxz = cbuffer.data(gk_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxyy = cbuffer.data(gk_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxyz = cbuffer.data(gk_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxzz = cbuffer.data(gk_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyyy = cbuffer.data(gk_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyyz = cbuffer.data(gk_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyzz = cbuffer.data(gk_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxzzz = cbuffer.data(gk_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyyy = cbuffer.data(gk_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyyz = cbuffer.data(gk_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyzz = cbuffer.data(gk_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyzzz = cbuffer.data(gk_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxzzzz = cbuffer.data(gk_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyyy = cbuffer.data(gk_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyyz = cbuffer.data(gk_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyzz = cbuffer.data(gk_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyzzz = cbuffer.data(gk_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyzzzz = cbuffer.data(gk_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xxxx_xxzzzzz = cbuffer.data(gk_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyyy = cbuffer.data(gk_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyyz = cbuffer.data(gk_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyzz = cbuffer.data(gk_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyzzz = cbuffer.data(gk_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyzzzz = cbuffer.data(gk_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xxxx_xyzzzzz = cbuffer.data(gk_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xxxx_xzzzzzz = cbuffer.data(gk_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyyy = cbuffer.data(gk_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyyz = cbuffer.data(gk_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyzz = cbuffer.data(gk_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyzzz = cbuffer.data(gk_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyzzzz = cbuffer.data(gk_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_z_xxxx_yyzzzzz = cbuffer.data(gk_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xxxx_yzzzzzz = cbuffer.data(gk_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xxxx_zzzzzzz = cbuffer.data(gk_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxx = cbuffer.data(gk_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxy = cbuffer.data(gk_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxz = cbuffer.data(gk_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxyy = cbuffer.data(gk_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxyz = cbuffer.data(gk_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxzz = cbuffer.data(gk_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyyy = cbuffer.data(gk_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyyz = cbuffer.data(gk_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyzz = cbuffer.data(gk_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxzzz = cbuffer.data(gk_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyyy = cbuffer.data(gk_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyyz = cbuffer.data(gk_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyzz = cbuffer.data(gk_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyzzz = cbuffer.data(gk_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxzzzz = cbuffer.data(gk_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyyy = cbuffer.data(gk_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyyz = cbuffer.data(gk_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyzz = cbuffer.data(gk_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyzzz = cbuffer.data(gk_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyzzzz = cbuffer.data(gk_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_xxxy_xxzzzzz = cbuffer.data(gk_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyyy = cbuffer.data(gk_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyyz = cbuffer.data(gk_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyzz = cbuffer.data(gk_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyzzz = cbuffer.data(gk_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyzzzz = cbuffer.data(gk_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_xxxy_xyzzzzz = cbuffer.data(gk_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_xxxy_xzzzzzz = cbuffer.data(gk_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyyy = cbuffer.data(gk_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyyz = cbuffer.data(gk_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyzz = cbuffer.data(gk_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyzzz = cbuffer.data(gk_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyzzzz = cbuffer.data(gk_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_xxxy_yyzzzzz = cbuffer.data(gk_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_xxxy_yzzzzzz = cbuffer.data(gk_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_xxxy_zzzzzzz = cbuffer.data(gk_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_xxxz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_xxxz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_xxxz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_z_xxxz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_xxxz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_xxxz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_xxyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_xxyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_z_xxyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_xxyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_xxyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_xxyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_xxyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_xxyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_xxyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_xxyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_xxyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_xxyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1259 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_z_xxzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_z_xxzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_z_xxzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_z_xxzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_z_xxzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_z_xxzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_z_xyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_z_xyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_z_xyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_z_xyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_z_xyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_z_xyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 1331 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1334 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1349 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1350 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1351 * ccomps * dcomps);

            auto g_0_z_xyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1352 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1353 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1354 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1355 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1356 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1357 * ccomps * dcomps);

            auto g_0_z_xyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1358 * ccomps * dcomps);

            auto g_0_z_xyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1359 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1360 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1361 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1362 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1363 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1364 * ccomps * dcomps);

            auto g_0_z_xyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1365 * ccomps * dcomps);

            auto g_0_z_xyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1366 * ccomps * dcomps);

            auto g_0_z_xyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1367 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1368 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1369 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1370 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1371 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1372 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1373 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1374 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1375 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1376 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1377 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1378 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1379 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1380 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1381 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1382 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1383 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1384 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1385 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1386 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1387 * ccomps * dcomps);

            auto g_0_z_xyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1388 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1389 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1390 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1391 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1392 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1393 * ccomps * dcomps);

            auto g_0_z_xyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1394 * ccomps * dcomps);

            auto g_0_z_xyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1395 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1396 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1397 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1398 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1399 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1400 * ccomps * dcomps);

            auto g_0_z_xyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1401 * ccomps * dcomps);

            auto g_0_z_xyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1402 * ccomps * dcomps);

            auto g_0_z_xyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1403 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1404 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1405 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1406 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1407 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1408 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1409 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1410 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1411 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1412 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1413 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1414 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1415 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1416 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1417 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1418 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1419 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1420 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1421 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1422 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1423 * ccomps * dcomps);

            auto g_0_z_xzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1424 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1425 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1426 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1427 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1428 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1429 * ccomps * dcomps);

            auto g_0_z_xzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1430 * ccomps * dcomps);

            auto g_0_z_xzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1431 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1432 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1433 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1434 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1435 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1436 * ccomps * dcomps);

            auto g_0_z_xzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1437 * ccomps * dcomps);

            auto g_0_z_xzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1438 * ccomps * dcomps);

            auto g_0_z_xzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1439 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxx = cbuffer.data(gk_geom_01_off + 1440 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxy = cbuffer.data(gk_geom_01_off + 1441 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxz = cbuffer.data(gk_geom_01_off + 1442 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxyy = cbuffer.data(gk_geom_01_off + 1443 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxyz = cbuffer.data(gk_geom_01_off + 1444 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxzz = cbuffer.data(gk_geom_01_off + 1445 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyyy = cbuffer.data(gk_geom_01_off + 1446 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyyz = cbuffer.data(gk_geom_01_off + 1447 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyzz = cbuffer.data(gk_geom_01_off + 1448 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxzzz = cbuffer.data(gk_geom_01_off + 1449 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyyy = cbuffer.data(gk_geom_01_off + 1450 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyyz = cbuffer.data(gk_geom_01_off + 1451 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyzz = cbuffer.data(gk_geom_01_off + 1452 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyzzz = cbuffer.data(gk_geom_01_off + 1453 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxzzzz = cbuffer.data(gk_geom_01_off + 1454 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyyy = cbuffer.data(gk_geom_01_off + 1455 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyyz = cbuffer.data(gk_geom_01_off + 1456 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyzz = cbuffer.data(gk_geom_01_off + 1457 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyzzz = cbuffer.data(gk_geom_01_off + 1458 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyzzzz = cbuffer.data(gk_geom_01_off + 1459 * ccomps * dcomps);

            auto g_0_z_yyyy_xxzzzzz = cbuffer.data(gk_geom_01_off + 1460 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyyy = cbuffer.data(gk_geom_01_off + 1461 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyyz = cbuffer.data(gk_geom_01_off + 1462 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyzz = cbuffer.data(gk_geom_01_off + 1463 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyzzz = cbuffer.data(gk_geom_01_off + 1464 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyzzzz = cbuffer.data(gk_geom_01_off + 1465 * ccomps * dcomps);

            auto g_0_z_yyyy_xyzzzzz = cbuffer.data(gk_geom_01_off + 1466 * ccomps * dcomps);

            auto g_0_z_yyyy_xzzzzzz = cbuffer.data(gk_geom_01_off + 1467 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyyy = cbuffer.data(gk_geom_01_off + 1468 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyyz = cbuffer.data(gk_geom_01_off + 1469 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyzz = cbuffer.data(gk_geom_01_off + 1470 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyzzz = cbuffer.data(gk_geom_01_off + 1471 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyzzzz = cbuffer.data(gk_geom_01_off + 1472 * ccomps * dcomps);

            auto g_0_z_yyyy_yyzzzzz = cbuffer.data(gk_geom_01_off + 1473 * ccomps * dcomps);

            auto g_0_z_yyyy_yzzzzzz = cbuffer.data(gk_geom_01_off + 1474 * ccomps * dcomps);

            auto g_0_z_yyyy_zzzzzzz = cbuffer.data(gk_geom_01_off + 1475 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1476 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1477 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1478 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1479 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1480 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1481 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1482 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1483 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1484 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1485 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1486 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1487 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1488 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1489 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1490 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1491 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1492 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1493 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1494 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1495 * ccomps * dcomps);

            auto g_0_z_yyyz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1496 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1497 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1498 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1499 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1500 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1501 * ccomps * dcomps);

            auto g_0_z_yyyz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1502 * ccomps * dcomps);

            auto g_0_z_yyyz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1503 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1504 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1505 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1506 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1507 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1508 * ccomps * dcomps);

            auto g_0_z_yyyz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1509 * ccomps * dcomps);

            auto g_0_z_yyyz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1510 * ccomps * dcomps);

            auto g_0_z_yyyz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1511 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1512 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1513 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1514 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1515 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1516 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1517 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1518 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1519 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1520 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1521 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1522 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1523 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1524 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1525 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1526 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1527 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1528 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1529 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1530 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1531 * ccomps * dcomps);

            auto g_0_z_yyzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1532 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1533 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1534 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1535 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1536 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1537 * ccomps * dcomps);

            auto g_0_z_yyzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1538 * ccomps * dcomps);

            auto g_0_z_yyzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1539 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1540 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1541 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1542 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1543 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1544 * ccomps * dcomps);

            auto g_0_z_yyzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1545 * ccomps * dcomps);

            auto g_0_z_yyzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1546 * ccomps * dcomps);

            auto g_0_z_yyzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1547 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1548 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1549 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1550 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1551 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1552 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1553 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1554 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1555 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1556 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1557 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1558 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1559 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1560 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1561 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1562 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1563 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1564 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1565 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1566 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1567 * ccomps * dcomps);

            auto g_0_z_yzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1568 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1569 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1570 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1571 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1572 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1573 * ccomps * dcomps);

            auto g_0_z_yzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1574 * ccomps * dcomps);

            auto g_0_z_yzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1575 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1576 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1577 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1578 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1579 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1580 * ccomps * dcomps);

            auto g_0_z_yzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1581 * ccomps * dcomps);

            auto g_0_z_yzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1582 * ccomps * dcomps);

            auto g_0_z_yzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1583 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxx = cbuffer.data(gk_geom_01_off + 1584 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxy = cbuffer.data(gk_geom_01_off + 1585 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxz = cbuffer.data(gk_geom_01_off + 1586 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxyy = cbuffer.data(gk_geom_01_off + 1587 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxyz = cbuffer.data(gk_geom_01_off + 1588 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxzz = cbuffer.data(gk_geom_01_off + 1589 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyyy = cbuffer.data(gk_geom_01_off + 1590 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyyz = cbuffer.data(gk_geom_01_off + 1591 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyzz = cbuffer.data(gk_geom_01_off + 1592 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxzzz = cbuffer.data(gk_geom_01_off + 1593 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyyy = cbuffer.data(gk_geom_01_off + 1594 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyyz = cbuffer.data(gk_geom_01_off + 1595 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyzz = cbuffer.data(gk_geom_01_off + 1596 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyzzz = cbuffer.data(gk_geom_01_off + 1597 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxzzzz = cbuffer.data(gk_geom_01_off + 1598 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyyy = cbuffer.data(gk_geom_01_off + 1599 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyyz = cbuffer.data(gk_geom_01_off + 1600 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyzz = cbuffer.data(gk_geom_01_off + 1601 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyzzz = cbuffer.data(gk_geom_01_off + 1602 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyzzzz = cbuffer.data(gk_geom_01_off + 1603 * ccomps * dcomps);

            auto g_0_z_zzzz_xxzzzzz = cbuffer.data(gk_geom_01_off + 1604 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyyy = cbuffer.data(gk_geom_01_off + 1605 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyyz = cbuffer.data(gk_geom_01_off + 1606 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyzz = cbuffer.data(gk_geom_01_off + 1607 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyzzz = cbuffer.data(gk_geom_01_off + 1608 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyzzzz = cbuffer.data(gk_geom_01_off + 1609 * ccomps * dcomps);

            auto g_0_z_zzzz_xyzzzzz = cbuffer.data(gk_geom_01_off + 1610 * ccomps * dcomps);

            auto g_0_z_zzzz_xzzzzzz = cbuffer.data(gk_geom_01_off + 1611 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyyy = cbuffer.data(gk_geom_01_off + 1612 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyyz = cbuffer.data(gk_geom_01_off + 1613 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyzz = cbuffer.data(gk_geom_01_off + 1614 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyzzz = cbuffer.data(gk_geom_01_off + 1615 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyzzzz = cbuffer.data(gk_geom_01_off + 1616 * ccomps * dcomps);

            auto g_0_z_zzzz_yyzzzzz = cbuffer.data(gk_geom_01_off + 1617 * ccomps * dcomps);

            auto g_0_z_zzzz_yzzzzzz = cbuffer.data(gk_geom_01_off + 1618 * ccomps * dcomps);

            auto g_0_z_zzzz_zzzzzzz = cbuffer.data(gk_geom_01_off + 1619 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GLSS

            const auto gl_geom_01_off = idx_geom_01_glxx + i * dcomps + j;

            auto g_0_x_xxxx_xxxxxxxx = cbuffer.data(gl_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxxz = cbuffer.data(gl_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxyy = cbuffer.data(gl_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxyz = cbuffer.data(gl_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxxzz = cbuffer.data(gl_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxyyy = cbuffer.data(gl_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxyyz = cbuffer.data(gl_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxyzz = cbuffer.data(gl_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxzzz = cbuffer.data(gl_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyyyy = cbuffer.data(gl_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyyyz = cbuffer.data(gl_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyyzz = cbuffer.data(gl_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyzzz = cbuffer.data(gl_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxzzzz = cbuffer.data(gl_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyyyy = cbuffer.data(gl_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyyyz = cbuffer.data(gl_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyyzz = cbuffer.data(gl_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyzzz = cbuffer.data(gl_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyzzzz = cbuffer.data(gl_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxzzzzz = cbuffer.data(gl_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyyyy = cbuffer.data(gl_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyyyz = cbuffer.data(gl_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyyzz = cbuffer.data(gl_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyzzz = cbuffer.data(gl_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyzzzz = cbuffer.data(gl_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyzzzzz = cbuffer.data(gl_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxx_xxzzzzzz = cbuffer.data(gl_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyyyy = cbuffer.data(gl_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyyyz = cbuffer.data(gl_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyyzz = cbuffer.data(gl_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyzzz = cbuffer.data(gl_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyzzzz = cbuffer.data(gl_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyzzzzz = cbuffer.data(gl_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxx_xyzzzzzz = cbuffer.data(gl_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxx_xzzzzzzz = cbuffer.data(gl_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyyyy = cbuffer.data(gl_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyyyz = cbuffer.data(gl_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyyzz = cbuffer.data(gl_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyzzz = cbuffer.data(gl_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyzzzz = cbuffer.data(gl_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyzzzzz = cbuffer.data(gl_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxx_yyzzzzzz = cbuffer.data(gl_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxx_yzzzzzzz = cbuffer.data(gl_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxx_zzzzzzzz = cbuffer.data(gl_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xxyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xxyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xxyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xxzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xxzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xxzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xxzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xxzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xxzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_xyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_xyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_xyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_xyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_xyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_xyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_xyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_xyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_xyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_xyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_xyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_xyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_xzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_xzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_xzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_xzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_xzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_xzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_x_yyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_x_yyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_x_yyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_x_yyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_x_yyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_x_yyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_x_yyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_x_yyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_x_yyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_x_yyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_x_yyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_x_yyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_x_yyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_x_yyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_x_yyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_x_yyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_x_yyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_x_yyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_x_yzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_x_yzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_x_yzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_x_yzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_x_yzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_x_yzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_x_zzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_x_zzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_x_zzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_x_zzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_x_zzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_x_zzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxxx = cbuffer.data(gl_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxxy = cbuffer.data(gl_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxxz = cbuffer.data(gl_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxyy = cbuffer.data(gl_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxyz = cbuffer.data(gl_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxxzz = cbuffer.data(gl_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxyyy = cbuffer.data(gl_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxyyz = cbuffer.data(gl_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxyzz = cbuffer.data(gl_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxzzz = cbuffer.data(gl_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyyyy = cbuffer.data(gl_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyyyz = cbuffer.data(gl_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyyzz = cbuffer.data(gl_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyzzz = cbuffer.data(gl_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxzzzz = cbuffer.data(gl_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyyyy = cbuffer.data(gl_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyyyz = cbuffer.data(gl_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyyzz = cbuffer.data(gl_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyzzz = cbuffer.data(gl_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyzzzz = cbuffer.data(gl_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxzzzzz = cbuffer.data(gl_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyyyy = cbuffer.data(gl_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyyyz = cbuffer.data(gl_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyyzz = cbuffer.data(gl_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyzzz = cbuffer.data(gl_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyzzzz = cbuffer.data(gl_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyzzzzz = cbuffer.data(gl_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xxxx_xxzzzzzz = cbuffer.data(gl_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyyyy = cbuffer.data(gl_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyyyz = cbuffer.data(gl_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyyzz = cbuffer.data(gl_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyzzz = cbuffer.data(gl_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyzzzz = cbuffer.data(gl_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyzzzzz = cbuffer.data(gl_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xxxx_xyzzzzzz = cbuffer.data(gl_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xxxx_xzzzzzzz = cbuffer.data(gl_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyyyy = cbuffer.data(gl_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyyyz = cbuffer.data(gl_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyyzz = cbuffer.data(gl_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyzzz = cbuffer.data(gl_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyzzzz = cbuffer.data(gl_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyzzzzz = cbuffer.data(gl_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xxxx_yyzzzzzz = cbuffer.data(gl_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xxxx_yzzzzzzz = cbuffer.data(gl_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xxxx_zzzzzzzz = cbuffer.data(gl_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_xxxy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_xxxy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_xxxy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_xxxy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_xxxy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_xxxy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_xxxz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_xxxz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_xxxz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_xxxz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_xxxz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_xxxz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_xxyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_xxyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_xxyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_xxyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_xxyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_xxyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_xxyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_xxyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_xxyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_xxyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_xxyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_xxyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_y_xxzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_y_xxzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_y_xxzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_y_xxzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_y_xxzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_y_xxzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_y_xyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_y_xyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_y_xyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_y_xyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_y_xyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_y_xyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_y_xyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_y_xyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_y_xyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_y_xyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_y_xyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_y_xyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_y_xyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_y_xyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_y_xyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_y_xyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_y_xyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_y_xyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_y_xzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_y_xzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_y_xzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_y_xzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_y_xzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_y_xzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_y_yyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_y_yyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_y_yyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_y_yyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_y_yyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_y_yyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_y_yyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_y_yyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_y_yyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_y_yyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_y_yyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_y_yyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_y_yyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_y_yyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_y_yyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_y_yyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_y_yyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_y_yyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1259 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_y_yzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_y_yzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_y_yzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_y_yzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_y_yzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_y_yzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1331 * ccomps * dcomps);

            auto g_0_y_zzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1334 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_y_zzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_y_zzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_y_zzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_y_zzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_y_zzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1349 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1350 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1351 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1352 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1353 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1354 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1355 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1356 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1357 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1358 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1359 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1360 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1361 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1362 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1363 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1364 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1365 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1366 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1367 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1368 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1369 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1370 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1371 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1372 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1373 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1374 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1375 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1376 * ccomps * dcomps);

            auto g_0_z_xxxx_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1377 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1378 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1379 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1380 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1381 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1382 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1383 * ccomps * dcomps);

            auto g_0_z_xxxx_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1384 * ccomps * dcomps);

            auto g_0_z_xxxx_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1385 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1386 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1387 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1388 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1389 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1390 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1391 * ccomps * dcomps);

            auto g_0_z_xxxx_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1392 * ccomps * dcomps);

            auto g_0_z_xxxx_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1393 * ccomps * dcomps);

            auto g_0_z_xxxx_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1394 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1395 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1396 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1397 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1398 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1399 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1400 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1401 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1402 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1403 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1404 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1405 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1406 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1407 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1408 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1409 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1410 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1411 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1412 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1413 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1414 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1415 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1416 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1417 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1418 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1419 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1420 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1421 * ccomps * dcomps);

            auto g_0_z_xxxy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1422 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1423 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1424 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1425 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1426 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1427 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1428 * ccomps * dcomps);

            auto g_0_z_xxxy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1429 * ccomps * dcomps);

            auto g_0_z_xxxy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1430 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1431 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1432 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1433 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1434 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1435 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1436 * ccomps * dcomps);

            auto g_0_z_xxxy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1437 * ccomps * dcomps);

            auto g_0_z_xxxy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1438 * ccomps * dcomps);

            auto g_0_z_xxxy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1439 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1440 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1441 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1442 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1443 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1444 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1445 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1446 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1447 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1448 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1449 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1450 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1451 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1452 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1453 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1454 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1455 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1456 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1457 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1458 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1459 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1460 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1461 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1462 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1463 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1464 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1465 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1466 * ccomps * dcomps);

            auto g_0_z_xxxz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1467 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1468 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1469 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1470 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1471 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1472 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1473 * ccomps * dcomps);

            auto g_0_z_xxxz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1474 * ccomps * dcomps);

            auto g_0_z_xxxz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1475 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1476 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1477 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1478 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1479 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1480 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1481 * ccomps * dcomps);

            auto g_0_z_xxxz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1482 * ccomps * dcomps);

            auto g_0_z_xxxz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1483 * ccomps * dcomps);

            auto g_0_z_xxxz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1484 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1485 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1486 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1487 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1488 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1489 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1490 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1491 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1492 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1493 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1494 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1495 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1496 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1497 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1498 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1499 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1500 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1501 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1502 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1503 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1504 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1505 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1506 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1507 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1508 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1509 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1510 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1511 * ccomps * dcomps);

            auto g_0_z_xxyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1512 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1513 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1514 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1515 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1516 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1517 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1518 * ccomps * dcomps);

            auto g_0_z_xxyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1519 * ccomps * dcomps);

            auto g_0_z_xxyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1520 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1521 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1522 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1523 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1524 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1525 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1526 * ccomps * dcomps);

            auto g_0_z_xxyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1527 * ccomps * dcomps);

            auto g_0_z_xxyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1528 * ccomps * dcomps);

            auto g_0_z_xxyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1529 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1530 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1531 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1532 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1533 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1534 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1535 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1536 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1537 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1538 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1539 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1540 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1541 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1542 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1543 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1544 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1545 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1546 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1547 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1548 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1549 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1550 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1551 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1552 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1553 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1554 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1555 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1556 * ccomps * dcomps);

            auto g_0_z_xxyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1557 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1558 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1559 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1560 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1561 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1562 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1563 * ccomps * dcomps);

            auto g_0_z_xxyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1564 * ccomps * dcomps);

            auto g_0_z_xxyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1565 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1566 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1567 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1568 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1569 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1570 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1571 * ccomps * dcomps);

            auto g_0_z_xxyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1572 * ccomps * dcomps);

            auto g_0_z_xxyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1573 * ccomps * dcomps);

            auto g_0_z_xxyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1574 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1575 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1576 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1577 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1578 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1579 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1580 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1581 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1582 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1583 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1584 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1585 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1586 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1587 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1588 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1589 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1590 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1591 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1592 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1593 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1594 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1595 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1596 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1597 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1598 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1599 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1600 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1601 * ccomps * dcomps);

            auto g_0_z_xxzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1602 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1603 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1604 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1605 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1606 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1607 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1608 * ccomps * dcomps);

            auto g_0_z_xxzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1609 * ccomps * dcomps);

            auto g_0_z_xxzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1610 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1611 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1612 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1613 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1614 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1615 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1616 * ccomps * dcomps);

            auto g_0_z_xxzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1617 * ccomps * dcomps);

            auto g_0_z_xxzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1618 * ccomps * dcomps);

            auto g_0_z_xxzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1619 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1620 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1621 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1622 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1623 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1624 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1625 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1626 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1627 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1628 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1629 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1630 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1631 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1632 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1633 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1634 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1635 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1636 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1637 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1638 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1639 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1640 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1641 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1642 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1643 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1644 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1645 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1646 * ccomps * dcomps);

            auto g_0_z_xyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1647 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1648 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1649 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1650 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1651 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1652 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1653 * ccomps * dcomps);

            auto g_0_z_xyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1654 * ccomps * dcomps);

            auto g_0_z_xyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1655 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1656 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1657 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1658 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1659 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1660 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1661 * ccomps * dcomps);

            auto g_0_z_xyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1662 * ccomps * dcomps);

            auto g_0_z_xyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1663 * ccomps * dcomps);

            auto g_0_z_xyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1664 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1665 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1666 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1667 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1668 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1669 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1670 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1671 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1672 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1673 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1674 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1675 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1676 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1677 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1678 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1679 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1680 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1681 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1682 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1683 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1684 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1685 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1686 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1687 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1688 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1689 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1690 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1691 * ccomps * dcomps);

            auto g_0_z_xyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1692 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1693 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1694 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1695 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1696 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1697 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1698 * ccomps * dcomps);

            auto g_0_z_xyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1699 * ccomps * dcomps);

            auto g_0_z_xyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1700 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1701 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1702 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1703 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1704 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1705 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1706 * ccomps * dcomps);

            auto g_0_z_xyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1707 * ccomps * dcomps);

            auto g_0_z_xyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1708 * ccomps * dcomps);

            auto g_0_z_xyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1709 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1710 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1711 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1712 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1713 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1714 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1715 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1716 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1717 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1718 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1719 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1720 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1721 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1722 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1723 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1724 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1725 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1726 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1727 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1728 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1729 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1730 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1731 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1732 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1733 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1734 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1735 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1736 * ccomps * dcomps);

            auto g_0_z_xyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1737 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1738 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1739 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1740 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1741 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1742 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1743 * ccomps * dcomps);

            auto g_0_z_xyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1744 * ccomps * dcomps);

            auto g_0_z_xyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1745 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1746 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1747 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1748 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1749 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1750 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1751 * ccomps * dcomps);

            auto g_0_z_xyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1752 * ccomps * dcomps);

            auto g_0_z_xyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1753 * ccomps * dcomps);

            auto g_0_z_xyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1754 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1755 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1756 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1757 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1758 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1759 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1760 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1761 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1762 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1763 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1764 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1765 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1766 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1767 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1768 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1769 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1770 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1771 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1772 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1773 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1774 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1775 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1776 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1777 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1778 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1779 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1780 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1781 * ccomps * dcomps);

            auto g_0_z_xzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1782 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1783 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1784 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1785 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1786 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1787 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1788 * ccomps * dcomps);

            auto g_0_z_xzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1789 * ccomps * dcomps);

            auto g_0_z_xzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1790 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1791 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1792 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1793 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1794 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1795 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1796 * ccomps * dcomps);

            auto g_0_z_xzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1797 * ccomps * dcomps);

            auto g_0_z_xzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1798 * ccomps * dcomps);

            auto g_0_z_xzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1799 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1800 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1801 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1802 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1803 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1804 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1805 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1806 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1807 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1808 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1809 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1810 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1811 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1812 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1813 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1814 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1815 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1816 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1817 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1818 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1819 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1820 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1821 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1822 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1823 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1824 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1825 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1826 * ccomps * dcomps);

            auto g_0_z_yyyy_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1827 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1828 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1829 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1830 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1831 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1832 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1833 * ccomps * dcomps);

            auto g_0_z_yyyy_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1834 * ccomps * dcomps);

            auto g_0_z_yyyy_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1835 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1836 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1837 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1838 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1839 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1840 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1841 * ccomps * dcomps);

            auto g_0_z_yyyy_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1842 * ccomps * dcomps);

            auto g_0_z_yyyy_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1843 * ccomps * dcomps);

            auto g_0_z_yyyy_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1844 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1845 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1846 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1847 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1848 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1849 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1850 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1851 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1852 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1853 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1854 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1855 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1856 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1857 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1858 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1859 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1860 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1861 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1862 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1863 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1864 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1865 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1866 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1867 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1868 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1869 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1870 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1871 * ccomps * dcomps);

            auto g_0_z_yyyz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1872 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1873 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1874 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1875 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1876 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1877 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1878 * ccomps * dcomps);

            auto g_0_z_yyyz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1879 * ccomps * dcomps);

            auto g_0_z_yyyz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1880 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1881 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1882 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1883 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1884 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1885 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1886 * ccomps * dcomps);

            auto g_0_z_yyyz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1887 * ccomps * dcomps);

            auto g_0_z_yyyz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1888 * ccomps * dcomps);

            auto g_0_z_yyyz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1889 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1890 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1891 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1892 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1893 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1894 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1895 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1896 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1897 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1898 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1899 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1900 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1901 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1902 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1903 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1904 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1905 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1906 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1907 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1908 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1909 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1910 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1911 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1912 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1913 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1914 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1915 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1916 * ccomps * dcomps);

            auto g_0_z_yyzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1917 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1918 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1919 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1920 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1921 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1922 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1923 * ccomps * dcomps);

            auto g_0_z_yyzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1924 * ccomps * dcomps);

            auto g_0_z_yyzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1925 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1926 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1927 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1928 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1929 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1930 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1931 * ccomps * dcomps);

            auto g_0_z_yyzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1932 * ccomps * dcomps);

            auto g_0_z_yyzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1933 * ccomps * dcomps);

            auto g_0_z_yyzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1934 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1935 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1936 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1937 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1938 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1939 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1940 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1941 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1942 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1943 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1944 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1945 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1946 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1947 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1948 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1949 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1950 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1951 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1952 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1953 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1954 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 1955 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 1956 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 1957 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 1958 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 1959 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 1960 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 1961 * ccomps * dcomps);

            auto g_0_z_yzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 1962 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 1963 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 1964 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 1965 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 1966 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 1967 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 1968 * ccomps * dcomps);

            auto g_0_z_yzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 1969 * ccomps * dcomps);

            auto g_0_z_yzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 1970 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 1971 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 1972 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 1973 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 1974 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 1975 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 1976 * ccomps * dcomps);

            auto g_0_z_yzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 1977 * ccomps * dcomps);

            auto g_0_z_yzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 1978 * ccomps * dcomps);

            auto g_0_z_yzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 1979 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxxx = cbuffer.data(gl_geom_01_off + 1980 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxxy = cbuffer.data(gl_geom_01_off + 1981 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxxz = cbuffer.data(gl_geom_01_off + 1982 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxyy = cbuffer.data(gl_geom_01_off + 1983 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxyz = cbuffer.data(gl_geom_01_off + 1984 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxxzz = cbuffer.data(gl_geom_01_off + 1985 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxyyy = cbuffer.data(gl_geom_01_off + 1986 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxyyz = cbuffer.data(gl_geom_01_off + 1987 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxyzz = cbuffer.data(gl_geom_01_off + 1988 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxzzz = cbuffer.data(gl_geom_01_off + 1989 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyyyy = cbuffer.data(gl_geom_01_off + 1990 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyyyz = cbuffer.data(gl_geom_01_off + 1991 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyyzz = cbuffer.data(gl_geom_01_off + 1992 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyzzz = cbuffer.data(gl_geom_01_off + 1993 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxzzzz = cbuffer.data(gl_geom_01_off + 1994 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyyyy = cbuffer.data(gl_geom_01_off + 1995 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyyyz = cbuffer.data(gl_geom_01_off + 1996 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyyzz = cbuffer.data(gl_geom_01_off + 1997 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyzzz = cbuffer.data(gl_geom_01_off + 1998 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyzzzz = cbuffer.data(gl_geom_01_off + 1999 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxzzzzz = cbuffer.data(gl_geom_01_off + 2000 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyyyy = cbuffer.data(gl_geom_01_off + 2001 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyyyz = cbuffer.data(gl_geom_01_off + 2002 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyyzz = cbuffer.data(gl_geom_01_off + 2003 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyzzz = cbuffer.data(gl_geom_01_off + 2004 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyzzzz = cbuffer.data(gl_geom_01_off + 2005 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyzzzzz = cbuffer.data(gl_geom_01_off + 2006 * ccomps * dcomps);

            auto g_0_z_zzzz_xxzzzzzz = cbuffer.data(gl_geom_01_off + 2007 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyyyy = cbuffer.data(gl_geom_01_off + 2008 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyyyz = cbuffer.data(gl_geom_01_off + 2009 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyyzz = cbuffer.data(gl_geom_01_off + 2010 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyzzz = cbuffer.data(gl_geom_01_off + 2011 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyzzzz = cbuffer.data(gl_geom_01_off + 2012 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyzzzzz = cbuffer.data(gl_geom_01_off + 2013 * ccomps * dcomps);

            auto g_0_z_zzzz_xyzzzzzz = cbuffer.data(gl_geom_01_off + 2014 * ccomps * dcomps);

            auto g_0_z_zzzz_xzzzzzzz = cbuffer.data(gl_geom_01_off + 2015 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyyyy = cbuffer.data(gl_geom_01_off + 2016 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyyyz = cbuffer.data(gl_geom_01_off + 2017 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyyzz = cbuffer.data(gl_geom_01_off + 2018 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyzzz = cbuffer.data(gl_geom_01_off + 2019 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyzzzz = cbuffer.data(gl_geom_01_off + 2020 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyzzzzz = cbuffer.data(gl_geom_01_off + 2021 * ccomps * dcomps);

            auto g_0_z_zzzz_yyzzzzzz = cbuffer.data(gl_geom_01_off + 2022 * ccomps * dcomps);

            auto g_0_z_zzzz_yzzzzzzz = cbuffer.data(gl_geom_01_off + 2023 * ccomps * dcomps);

            auto g_0_z_zzzz_zzzzzzzz = cbuffer.data(gl_geom_01_off + 2024 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hkxx

            const auto hk_geom_01_off = idx_geom_01_hkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxx_xxxxxxx = cbuffer.data(hk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxxxy = cbuffer.data(hk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxxxz = cbuffer.data(hk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxxyy = cbuffer.data(hk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxxyz = cbuffer.data(hk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxxzz = cbuffer.data(hk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxyyy = cbuffer.data(hk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxyyz = cbuffer.data(hk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxyzz = cbuffer.data(hk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxzzz = cbuffer.data(hk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyyyy = cbuffer.data(hk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyyyz = cbuffer.data(hk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyyzz = cbuffer.data(hk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyzzz = cbuffer.data(hk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxzzzz = cbuffer.data(hk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyyyy = cbuffer.data(hk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyyyz = cbuffer.data(hk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyyzz = cbuffer.data(hk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyzzz = cbuffer.data(hk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyzzzz = cbuffer.data(hk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxzzzzz = cbuffer.data(hk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyyyy = cbuffer.data(hk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyyyz = cbuffer.data(hk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyyzz = cbuffer.data(hk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyzzz = cbuffer.data(hk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyzzzz = cbuffer.data(hk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyzzzzz = cbuffer.data(hk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxx_xzzzzzz = cbuffer.data(hk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyyyy = cbuffer.data(hk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyyyz = cbuffer.data(hk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyyzz = cbuffer.data(hk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyzzz = cbuffer.data(hk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyzzzz = cbuffer.data(hk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyzzzzz = cbuffer.data(hk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxx_yzzzzzz = cbuffer.data(hk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxx_zzzzzzz = cbuffer.data(hk_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxxxx, g_0_x_xxxx_xxxxxxxx, g_0_x_xxxx_xxxxxxxy, g_0_x_xxxx_xxxxxxxz, g_0_x_xxxx_xxxxxxy, g_0_x_xxxx_xxxxxxyy, g_0_x_xxxx_xxxxxxyz, g_0_x_xxxx_xxxxxxz, g_0_x_xxxx_xxxxxxzz, g_0_x_xxxx_xxxxxyy, g_0_x_xxxx_xxxxxyyy, g_0_x_xxxx_xxxxxyyz, g_0_x_xxxx_xxxxxyz, g_0_x_xxxx_xxxxxyzz, g_0_x_xxxx_xxxxxzz, g_0_x_xxxx_xxxxxzzz, g_0_x_xxxx_xxxxyyy, g_0_x_xxxx_xxxxyyyy, g_0_x_xxxx_xxxxyyyz, g_0_x_xxxx_xxxxyyz, g_0_x_xxxx_xxxxyyzz, g_0_x_xxxx_xxxxyzz, g_0_x_xxxx_xxxxyzzz, g_0_x_xxxx_xxxxzzz, g_0_x_xxxx_xxxxzzzz, g_0_x_xxxx_xxxyyyy, g_0_x_xxxx_xxxyyyyy, g_0_x_xxxx_xxxyyyyz, g_0_x_xxxx_xxxyyyz, g_0_x_xxxx_xxxyyyzz, g_0_x_xxxx_xxxyyzz, g_0_x_xxxx_xxxyyzzz, g_0_x_xxxx_xxxyzzz, g_0_x_xxxx_xxxyzzzz, g_0_x_xxxx_xxxzzzz, g_0_x_xxxx_xxxzzzzz, g_0_x_xxxx_xxyyyyy, g_0_x_xxxx_xxyyyyyy, g_0_x_xxxx_xxyyyyyz, g_0_x_xxxx_xxyyyyz, g_0_x_xxxx_xxyyyyzz, g_0_x_xxxx_xxyyyzz, g_0_x_xxxx_xxyyyzzz, g_0_x_xxxx_xxyyzzz, g_0_x_xxxx_xxyyzzzz, g_0_x_xxxx_xxyzzzz, g_0_x_xxxx_xxyzzzzz, g_0_x_xxxx_xxzzzzz, g_0_x_xxxx_xxzzzzzz, g_0_x_xxxx_xyyyyyy, g_0_x_xxxx_xyyyyyyy, g_0_x_xxxx_xyyyyyyz, g_0_x_xxxx_xyyyyyz, g_0_x_xxxx_xyyyyyzz, g_0_x_xxxx_xyyyyzz, g_0_x_xxxx_xyyyyzzz, g_0_x_xxxx_xyyyzzz, g_0_x_xxxx_xyyyzzzz, g_0_x_xxxx_xyyzzzz, g_0_x_xxxx_xyyzzzzz, g_0_x_xxxx_xyzzzzz, g_0_x_xxxx_xyzzzzzz, g_0_x_xxxx_xzzzzzz, g_0_x_xxxx_xzzzzzzz, g_0_x_xxxx_yyyyyyy, g_0_x_xxxx_yyyyyyz, g_0_x_xxxx_yyyyyzz, g_0_x_xxxx_yyyyzzz, g_0_x_xxxx_yyyzzzz, g_0_x_xxxx_yyzzzzz, g_0_x_xxxx_yzzzzzz, g_0_x_xxxx_zzzzzzz, g_0_x_xxxxx_xxxxxxx, g_0_x_xxxxx_xxxxxxy, g_0_x_xxxxx_xxxxxxz, g_0_x_xxxxx_xxxxxyy, g_0_x_xxxxx_xxxxxyz, g_0_x_xxxxx_xxxxxzz, g_0_x_xxxxx_xxxxyyy, g_0_x_xxxxx_xxxxyyz, g_0_x_xxxxx_xxxxyzz, g_0_x_xxxxx_xxxxzzz, g_0_x_xxxxx_xxxyyyy, g_0_x_xxxxx_xxxyyyz, g_0_x_xxxxx_xxxyyzz, g_0_x_xxxxx_xxxyzzz, g_0_x_xxxxx_xxxzzzz, g_0_x_xxxxx_xxyyyyy, g_0_x_xxxxx_xxyyyyz, g_0_x_xxxxx_xxyyyzz, g_0_x_xxxxx_xxyyzzz, g_0_x_xxxxx_xxyzzzz, g_0_x_xxxxx_xxzzzzz, g_0_x_xxxxx_xyyyyyy, g_0_x_xxxxx_xyyyyyz, g_0_x_xxxxx_xyyyyzz, g_0_x_xxxxx_xyyyzzz, g_0_x_xxxxx_xyyzzzz, g_0_x_xxxxx_xyzzzzz, g_0_x_xxxxx_xzzzzzz, g_0_x_xxxxx_yyyyyyy, g_0_x_xxxxx_yyyyyyz, g_0_x_xxxxx_yyyyyzz, g_0_x_xxxxx_yyyyzzz, g_0_x_xxxxx_yyyzzzz, g_0_x_xxxxx_yyzzzzz, g_0_x_xxxxx_yzzzzzz, g_0_x_xxxxx_zzzzzzz, g_xxxx_xxxxxxx, g_xxxx_xxxxxxy, g_xxxx_xxxxxxz, g_xxxx_xxxxxyy, g_xxxx_xxxxxyz, g_xxxx_xxxxxzz, g_xxxx_xxxxyyy, g_xxxx_xxxxyyz, g_xxxx_xxxxyzz, g_xxxx_xxxxzzz, g_xxxx_xxxyyyy, g_xxxx_xxxyyyz, g_xxxx_xxxyyzz, g_xxxx_xxxyzzz, g_xxxx_xxxzzzz, g_xxxx_xxyyyyy, g_xxxx_xxyyyyz, g_xxxx_xxyyyzz, g_xxxx_xxyyzzz, g_xxxx_xxyzzzz, g_xxxx_xxzzzzz, g_xxxx_xyyyyyy, g_xxxx_xyyyyyz, g_xxxx_xyyyyzz, g_xxxx_xyyyzzz, g_xxxx_xyyzzzz, g_xxxx_xyzzzzz, g_xxxx_xzzzzzz, g_xxxx_yyyyyyy, g_xxxx_yyyyyyz, g_xxxx_yyyyyzz, g_xxxx_yyyyzzz, g_xxxx_yyyzzzz, g_xxxx_yyzzzzz, g_xxxx_yzzzzzz, g_xxxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxx_xxxxxxx[k] = g_xxxx_xxxxxxx[k] - g_0_x_xxxx_xxxxxxx[k] * ab_x + g_0_x_xxxx_xxxxxxxx[k];

                g_0_x_xxxxx_xxxxxxy[k] = g_xxxx_xxxxxxy[k] - g_0_x_xxxx_xxxxxxy[k] * ab_x + g_0_x_xxxx_xxxxxxxy[k];

                g_0_x_xxxxx_xxxxxxz[k] = g_xxxx_xxxxxxz[k] - g_0_x_xxxx_xxxxxxz[k] * ab_x + g_0_x_xxxx_xxxxxxxz[k];

                g_0_x_xxxxx_xxxxxyy[k] = g_xxxx_xxxxxyy[k] - g_0_x_xxxx_xxxxxyy[k] * ab_x + g_0_x_xxxx_xxxxxxyy[k];

                g_0_x_xxxxx_xxxxxyz[k] = g_xxxx_xxxxxyz[k] - g_0_x_xxxx_xxxxxyz[k] * ab_x + g_0_x_xxxx_xxxxxxyz[k];

                g_0_x_xxxxx_xxxxxzz[k] = g_xxxx_xxxxxzz[k] - g_0_x_xxxx_xxxxxzz[k] * ab_x + g_0_x_xxxx_xxxxxxzz[k];

                g_0_x_xxxxx_xxxxyyy[k] = g_xxxx_xxxxyyy[k] - g_0_x_xxxx_xxxxyyy[k] * ab_x + g_0_x_xxxx_xxxxxyyy[k];

                g_0_x_xxxxx_xxxxyyz[k] = g_xxxx_xxxxyyz[k] - g_0_x_xxxx_xxxxyyz[k] * ab_x + g_0_x_xxxx_xxxxxyyz[k];

                g_0_x_xxxxx_xxxxyzz[k] = g_xxxx_xxxxyzz[k] - g_0_x_xxxx_xxxxyzz[k] * ab_x + g_0_x_xxxx_xxxxxyzz[k];

                g_0_x_xxxxx_xxxxzzz[k] = g_xxxx_xxxxzzz[k] - g_0_x_xxxx_xxxxzzz[k] * ab_x + g_0_x_xxxx_xxxxxzzz[k];

                g_0_x_xxxxx_xxxyyyy[k] = g_xxxx_xxxyyyy[k] - g_0_x_xxxx_xxxyyyy[k] * ab_x + g_0_x_xxxx_xxxxyyyy[k];

                g_0_x_xxxxx_xxxyyyz[k] = g_xxxx_xxxyyyz[k] - g_0_x_xxxx_xxxyyyz[k] * ab_x + g_0_x_xxxx_xxxxyyyz[k];

                g_0_x_xxxxx_xxxyyzz[k] = g_xxxx_xxxyyzz[k] - g_0_x_xxxx_xxxyyzz[k] * ab_x + g_0_x_xxxx_xxxxyyzz[k];

                g_0_x_xxxxx_xxxyzzz[k] = g_xxxx_xxxyzzz[k] - g_0_x_xxxx_xxxyzzz[k] * ab_x + g_0_x_xxxx_xxxxyzzz[k];

                g_0_x_xxxxx_xxxzzzz[k] = g_xxxx_xxxzzzz[k] - g_0_x_xxxx_xxxzzzz[k] * ab_x + g_0_x_xxxx_xxxxzzzz[k];

                g_0_x_xxxxx_xxyyyyy[k] = g_xxxx_xxyyyyy[k] - g_0_x_xxxx_xxyyyyy[k] * ab_x + g_0_x_xxxx_xxxyyyyy[k];

                g_0_x_xxxxx_xxyyyyz[k] = g_xxxx_xxyyyyz[k] - g_0_x_xxxx_xxyyyyz[k] * ab_x + g_0_x_xxxx_xxxyyyyz[k];

                g_0_x_xxxxx_xxyyyzz[k] = g_xxxx_xxyyyzz[k] - g_0_x_xxxx_xxyyyzz[k] * ab_x + g_0_x_xxxx_xxxyyyzz[k];

                g_0_x_xxxxx_xxyyzzz[k] = g_xxxx_xxyyzzz[k] - g_0_x_xxxx_xxyyzzz[k] * ab_x + g_0_x_xxxx_xxxyyzzz[k];

                g_0_x_xxxxx_xxyzzzz[k] = g_xxxx_xxyzzzz[k] - g_0_x_xxxx_xxyzzzz[k] * ab_x + g_0_x_xxxx_xxxyzzzz[k];

                g_0_x_xxxxx_xxzzzzz[k] = g_xxxx_xxzzzzz[k] - g_0_x_xxxx_xxzzzzz[k] * ab_x + g_0_x_xxxx_xxxzzzzz[k];

                g_0_x_xxxxx_xyyyyyy[k] = g_xxxx_xyyyyyy[k] - g_0_x_xxxx_xyyyyyy[k] * ab_x + g_0_x_xxxx_xxyyyyyy[k];

                g_0_x_xxxxx_xyyyyyz[k] = g_xxxx_xyyyyyz[k] - g_0_x_xxxx_xyyyyyz[k] * ab_x + g_0_x_xxxx_xxyyyyyz[k];

                g_0_x_xxxxx_xyyyyzz[k] = g_xxxx_xyyyyzz[k] - g_0_x_xxxx_xyyyyzz[k] * ab_x + g_0_x_xxxx_xxyyyyzz[k];

                g_0_x_xxxxx_xyyyzzz[k] = g_xxxx_xyyyzzz[k] - g_0_x_xxxx_xyyyzzz[k] * ab_x + g_0_x_xxxx_xxyyyzzz[k];

                g_0_x_xxxxx_xyyzzzz[k] = g_xxxx_xyyzzzz[k] - g_0_x_xxxx_xyyzzzz[k] * ab_x + g_0_x_xxxx_xxyyzzzz[k];

                g_0_x_xxxxx_xyzzzzz[k] = g_xxxx_xyzzzzz[k] - g_0_x_xxxx_xyzzzzz[k] * ab_x + g_0_x_xxxx_xxyzzzzz[k];

                g_0_x_xxxxx_xzzzzzz[k] = g_xxxx_xzzzzzz[k] - g_0_x_xxxx_xzzzzzz[k] * ab_x + g_0_x_xxxx_xxzzzzzz[k];

                g_0_x_xxxxx_yyyyyyy[k] = g_xxxx_yyyyyyy[k] - g_0_x_xxxx_yyyyyyy[k] * ab_x + g_0_x_xxxx_xyyyyyyy[k];

                g_0_x_xxxxx_yyyyyyz[k] = g_xxxx_yyyyyyz[k] - g_0_x_xxxx_yyyyyyz[k] * ab_x + g_0_x_xxxx_xyyyyyyz[k];

                g_0_x_xxxxx_yyyyyzz[k] = g_xxxx_yyyyyzz[k] - g_0_x_xxxx_yyyyyzz[k] * ab_x + g_0_x_xxxx_xyyyyyzz[k];

                g_0_x_xxxxx_yyyyzzz[k] = g_xxxx_yyyyzzz[k] - g_0_x_xxxx_yyyyzzz[k] * ab_x + g_0_x_xxxx_xyyyyzzz[k];

                g_0_x_xxxxx_yyyzzzz[k] = g_xxxx_yyyzzzz[k] - g_0_x_xxxx_yyyzzzz[k] * ab_x + g_0_x_xxxx_xyyyzzzz[k];

                g_0_x_xxxxx_yyzzzzz[k] = g_xxxx_yyzzzzz[k] - g_0_x_xxxx_yyzzzzz[k] * ab_x + g_0_x_xxxx_xyyzzzzz[k];

                g_0_x_xxxxx_yzzzzzz[k] = g_xxxx_yzzzzzz[k] - g_0_x_xxxx_yzzzzzz[k] * ab_x + g_0_x_xxxx_xyzzzzzz[k];

                g_0_x_xxxxx_zzzzzzz[k] = g_xxxx_zzzzzzz[k] - g_0_x_xxxx_zzzzzzz[k] * ab_x + g_0_x_xxxx_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxy_xxxxxxx = cbuffer.data(hk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxxxy = cbuffer.data(hk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxxxz = cbuffer.data(hk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxxyy = cbuffer.data(hk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxxyz = cbuffer.data(hk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxxzz = cbuffer.data(hk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxyyy = cbuffer.data(hk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxyyz = cbuffer.data(hk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxyzz = cbuffer.data(hk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxzzz = cbuffer.data(hk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyyyy = cbuffer.data(hk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyyyz = cbuffer.data(hk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyyzz = cbuffer.data(hk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyzzz = cbuffer.data(hk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxzzzz = cbuffer.data(hk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyyyy = cbuffer.data(hk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyyyz = cbuffer.data(hk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyyzz = cbuffer.data(hk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyzzz = cbuffer.data(hk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyzzzz = cbuffer.data(hk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxzzzzz = cbuffer.data(hk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyyyy = cbuffer.data(hk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyyyz = cbuffer.data(hk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyyzz = cbuffer.data(hk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyzzz = cbuffer.data(hk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyzzzz = cbuffer.data(hk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyzzzzz = cbuffer.data(hk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxxy_xzzzzzz = cbuffer.data(hk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyyyy = cbuffer.data(hk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyyyz = cbuffer.data(hk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyyzz = cbuffer.data(hk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyzzz = cbuffer.data(hk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyzzzz = cbuffer.data(hk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyzzzzz = cbuffer.data(hk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxxy_yzzzzzz = cbuffer.data(hk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxxy_zzzzzzz = cbuffer.data(hk_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxxxx, g_0_x_xxxx_xxxxxxxy, g_0_x_xxxx_xxxxxxy, g_0_x_xxxx_xxxxxxyy, g_0_x_xxxx_xxxxxxyz, g_0_x_xxxx_xxxxxxz, g_0_x_xxxx_xxxxxyy, g_0_x_xxxx_xxxxxyyy, g_0_x_xxxx_xxxxxyyz, g_0_x_xxxx_xxxxxyz, g_0_x_xxxx_xxxxxyzz, g_0_x_xxxx_xxxxxzz, g_0_x_xxxx_xxxxyyy, g_0_x_xxxx_xxxxyyyy, g_0_x_xxxx_xxxxyyyz, g_0_x_xxxx_xxxxyyz, g_0_x_xxxx_xxxxyyzz, g_0_x_xxxx_xxxxyzz, g_0_x_xxxx_xxxxyzzz, g_0_x_xxxx_xxxxzzz, g_0_x_xxxx_xxxyyyy, g_0_x_xxxx_xxxyyyyy, g_0_x_xxxx_xxxyyyyz, g_0_x_xxxx_xxxyyyz, g_0_x_xxxx_xxxyyyzz, g_0_x_xxxx_xxxyyzz, g_0_x_xxxx_xxxyyzzz, g_0_x_xxxx_xxxyzzz, g_0_x_xxxx_xxxyzzzz, g_0_x_xxxx_xxxzzzz, g_0_x_xxxx_xxyyyyy, g_0_x_xxxx_xxyyyyyy, g_0_x_xxxx_xxyyyyyz, g_0_x_xxxx_xxyyyyz, g_0_x_xxxx_xxyyyyzz, g_0_x_xxxx_xxyyyzz, g_0_x_xxxx_xxyyyzzz, g_0_x_xxxx_xxyyzzz, g_0_x_xxxx_xxyyzzzz, g_0_x_xxxx_xxyzzzz, g_0_x_xxxx_xxyzzzzz, g_0_x_xxxx_xxzzzzz, g_0_x_xxxx_xyyyyyy, g_0_x_xxxx_xyyyyyyy, g_0_x_xxxx_xyyyyyyz, g_0_x_xxxx_xyyyyyz, g_0_x_xxxx_xyyyyyzz, g_0_x_xxxx_xyyyyzz, g_0_x_xxxx_xyyyyzzz, g_0_x_xxxx_xyyyzzz, g_0_x_xxxx_xyyyzzzz, g_0_x_xxxx_xyyzzzz, g_0_x_xxxx_xyyzzzzz, g_0_x_xxxx_xyzzzzz, g_0_x_xxxx_xyzzzzzz, g_0_x_xxxx_xzzzzzz, g_0_x_xxxx_yyyyyyy, g_0_x_xxxx_yyyyyyyy, g_0_x_xxxx_yyyyyyyz, g_0_x_xxxx_yyyyyyz, g_0_x_xxxx_yyyyyyzz, g_0_x_xxxx_yyyyyzz, g_0_x_xxxx_yyyyyzzz, g_0_x_xxxx_yyyyzzz, g_0_x_xxxx_yyyyzzzz, g_0_x_xxxx_yyyzzzz, g_0_x_xxxx_yyyzzzzz, g_0_x_xxxx_yyzzzzz, g_0_x_xxxx_yyzzzzzz, g_0_x_xxxx_yzzzzzz, g_0_x_xxxx_yzzzzzzz, g_0_x_xxxx_zzzzzzz, g_0_x_xxxxy_xxxxxxx, g_0_x_xxxxy_xxxxxxy, g_0_x_xxxxy_xxxxxxz, g_0_x_xxxxy_xxxxxyy, g_0_x_xxxxy_xxxxxyz, g_0_x_xxxxy_xxxxxzz, g_0_x_xxxxy_xxxxyyy, g_0_x_xxxxy_xxxxyyz, g_0_x_xxxxy_xxxxyzz, g_0_x_xxxxy_xxxxzzz, g_0_x_xxxxy_xxxyyyy, g_0_x_xxxxy_xxxyyyz, g_0_x_xxxxy_xxxyyzz, g_0_x_xxxxy_xxxyzzz, g_0_x_xxxxy_xxxzzzz, g_0_x_xxxxy_xxyyyyy, g_0_x_xxxxy_xxyyyyz, g_0_x_xxxxy_xxyyyzz, g_0_x_xxxxy_xxyyzzz, g_0_x_xxxxy_xxyzzzz, g_0_x_xxxxy_xxzzzzz, g_0_x_xxxxy_xyyyyyy, g_0_x_xxxxy_xyyyyyz, g_0_x_xxxxy_xyyyyzz, g_0_x_xxxxy_xyyyzzz, g_0_x_xxxxy_xyyzzzz, g_0_x_xxxxy_xyzzzzz, g_0_x_xxxxy_xzzzzzz, g_0_x_xxxxy_yyyyyyy, g_0_x_xxxxy_yyyyyyz, g_0_x_xxxxy_yyyyyzz, g_0_x_xxxxy_yyyyzzz, g_0_x_xxxxy_yyyzzzz, g_0_x_xxxxy_yyzzzzz, g_0_x_xxxxy_yzzzzzz, g_0_x_xxxxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxy_xxxxxxx[k] = -g_0_x_xxxx_xxxxxxx[k] * ab_y + g_0_x_xxxx_xxxxxxxy[k];

                g_0_x_xxxxy_xxxxxxy[k] = -g_0_x_xxxx_xxxxxxy[k] * ab_y + g_0_x_xxxx_xxxxxxyy[k];

                g_0_x_xxxxy_xxxxxxz[k] = -g_0_x_xxxx_xxxxxxz[k] * ab_y + g_0_x_xxxx_xxxxxxyz[k];

                g_0_x_xxxxy_xxxxxyy[k] = -g_0_x_xxxx_xxxxxyy[k] * ab_y + g_0_x_xxxx_xxxxxyyy[k];

                g_0_x_xxxxy_xxxxxyz[k] = -g_0_x_xxxx_xxxxxyz[k] * ab_y + g_0_x_xxxx_xxxxxyyz[k];

                g_0_x_xxxxy_xxxxxzz[k] = -g_0_x_xxxx_xxxxxzz[k] * ab_y + g_0_x_xxxx_xxxxxyzz[k];

                g_0_x_xxxxy_xxxxyyy[k] = -g_0_x_xxxx_xxxxyyy[k] * ab_y + g_0_x_xxxx_xxxxyyyy[k];

                g_0_x_xxxxy_xxxxyyz[k] = -g_0_x_xxxx_xxxxyyz[k] * ab_y + g_0_x_xxxx_xxxxyyyz[k];

                g_0_x_xxxxy_xxxxyzz[k] = -g_0_x_xxxx_xxxxyzz[k] * ab_y + g_0_x_xxxx_xxxxyyzz[k];

                g_0_x_xxxxy_xxxxzzz[k] = -g_0_x_xxxx_xxxxzzz[k] * ab_y + g_0_x_xxxx_xxxxyzzz[k];

                g_0_x_xxxxy_xxxyyyy[k] = -g_0_x_xxxx_xxxyyyy[k] * ab_y + g_0_x_xxxx_xxxyyyyy[k];

                g_0_x_xxxxy_xxxyyyz[k] = -g_0_x_xxxx_xxxyyyz[k] * ab_y + g_0_x_xxxx_xxxyyyyz[k];

                g_0_x_xxxxy_xxxyyzz[k] = -g_0_x_xxxx_xxxyyzz[k] * ab_y + g_0_x_xxxx_xxxyyyzz[k];

                g_0_x_xxxxy_xxxyzzz[k] = -g_0_x_xxxx_xxxyzzz[k] * ab_y + g_0_x_xxxx_xxxyyzzz[k];

                g_0_x_xxxxy_xxxzzzz[k] = -g_0_x_xxxx_xxxzzzz[k] * ab_y + g_0_x_xxxx_xxxyzzzz[k];

                g_0_x_xxxxy_xxyyyyy[k] = -g_0_x_xxxx_xxyyyyy[k] * ab_y + g_0_x_xxxx_xxyyyyyy[k];

                g_0_x_xxxxy_xxyyyyz[k] = -g_0_x_xxxx_xxyyyyz[k] * ab_y + g_0_x_xxxx_xxyyyyyz[k];

                g_0_x_xxxxy_xxyyyzz[k] = -g_0_x_xxxx_xxyyyzz[k] * ab_y + g_0_x_xxxx_xxyyyyzz[k];

                g_0_x_xxxxy_xxyyzzz[k] = -g_0_x_xxxx_xxyyzzz[k] * ab_y + g_0_x_xxxx_xxyyyzzz[k];

                g_0_x_xxxxy_xxyzzzz[k] = -g_0_x_xxxx_xxyzzzz[k] * ab_y + g_0_x_xxxx_xxyyzzzz[k];

                g_0_x_xxxxy_xxzzzzz[k] = -g_0_x_xxxx_xxzzzzz[k] * ab_y + g_0_x_xxxx_xxyzzzzz[k];

                g_0_x_xxxxy_xyyyyyy[k] = -g_0_x_xxxx_xyyyyyy[k] * ab_y + g_0_x_xxxx_xyyyyyyy[k];

                g_0_x_xxxxy_xyyyyyz[k] = -g_0_x_xxxx_xyyyyyz[k] * ab_y + g_0_x_xxxx_xyyyyyyz[k];

                g_0_x_xxxxy_xyyyyzz[k] = -g_0_x_xxxx_xyyyyzz[k] * ab_y + g_0_x_xxxx_xyyyyyzz[k];

                g_0_x_xxxxy_xyyyzzz[k] = -g_0_x_xxxx_xyyyzzz[k] * ab_y + g_0_x_xxxx_xyyyyzzz[k];

                g_0_x_xxxxy_xyyzzzz[k] = -g_0_x_xxxx_xyyzzzz[k] * ab_y + g_0_x_xxxx_xyyyzzzz[k];

                g_0_x_xxxxy_xyzzzzz[k] = -g_0_x_xxxx_xyzzzzz[k] * ab_y + g_0_x_xxxx_xyyzzzzz[k];

                g_0_x_xxxxy_xzzzzzz[k] = -g_0_x_xxxx_xzzzzzz[k] * ab_y + g_0_x_xxxx_xyzzzzzz[k];

                g_0_x_xxxxy_yyyyyyy[k] = -g_0_x_xxxx_yyyyyyy[k] * ab_y + g_0_x_xxxx_yyyyyyyy[k];

                g_0_x_xxxxy_yyyyyyz[k] = -g_0_x_xxxx_yyyyyyz[k] * ab_y + g_0_x_xxxx_yyyyyyyz[k];

                g_0_x_xxxxy_yyyyyzz[k] = -g_0_x_xxxx_yyyyyzz[k] * ab_y + g_0_x_xxxx_yyyyyyzz[k];

                g_0_x_xxxxy_yyyyzzz[k] = -g_0_x_xxxx_yyyyzzz[k] * ab_y + g_0_x_xxxx_yyyyyzzz[k];

                g_0_x_xxxxy_yyyzzzz[k] = -g_0_x_xxxx_yyyzzzz[k] * ab_y + g_0_x_xxxx_yyyyzzzz[k];

                g_0_x_xxxxy_yyzzzzz[k] = -g_0_x_xxxx_yyzzzzz[k] * ab_y + g_0_x_xxxx_yyyzzzzz[k];

                g_0_x_xxxxy_yzzzzzz[k] = -g_0_x_xxxx_yzzzzzz[k] * ab_y + g_0_x_xxxx_yyzzzzzz[k];

                g_0_x_xxxxy_zzzzzzz[k] = -g_0_x_xxxx_zzzzzzz[k] * ab_y + g_0_x_xxxx_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxz_xxxxxxx = cbuffer.data(hk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxxxy = cbuffer.data(hk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxxxz = cbuffer.data(hk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxxyy = cbuffer.data(hk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxxyz = cbuffer.data(hk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxxzz = cbuffer.data(hk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxyyy = cbuffer.data(hk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxyyz = cbuffer.data(hk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxyzz = cbuffer.data(hk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxzzz = cbuffer.data(hk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyyyy = cbuffer.data(hk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyyyz = cbuffer.data(hk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyyzz = cbuffer.data(hk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyzzz = cbuffer.data(hk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxzzzz = cbuffer.data(hk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyyyy = cbuffer.data(hk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyyyz = cbuffer.data(hk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyyzz = cbuffer.data(hk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyzzz = cbuffer.data(hk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyzzzz = cbuffer.data(hk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxzzzzz = cbuffer.data(hk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyyyy = cbuffer.data(hk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyyyz = cbuffer.data(hk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyyzz = cbuffer.data(hk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyzzz = cbuffer.data(hk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyzzzz = cbuffer.data(hk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyzzzzz = cbuffer.data(hk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxxz_xzzzzzz = cbuffer.data(hk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyyyy = cbuffer.data(hk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyyyz = cbuffer.data(hk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyyzz = cbuffer.data(hk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyzzz = cbuffer.data(hk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyzzzz = cbuffer.data(hk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyzzzzz = cbuffer.data(hk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxxz_yzzzzzz = cbuffer.data(hk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxxz_zzzzzzz = cbuffer.data(hk_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxxxx, g_0_x_xxxx_xxxxxxxz, g_0_x_xxxx_xxxxxxy, g_0_x_xxxx_xxxxxxyz, g_0_x_xxxx_xxxxxxz, g_0_x_xxxx_xxxxxxzz, g_0_x_xxxx_xxxxxyy, g_0_x_xxxx_xxxxxyyz, g_0_x_xxxx_xxxxxyz, g_0_x_xxxx_xxxxxyzz, g_0_x_xxxx_xxxxxzz, g_0_x_xxxx_xxxxxzzz, g_0_x_xxxx_xxxxyyy, g_0_x_xxxx_xxxxyyyz, g_0_x_xxxx_xxxxyyz, g_0_x_xxxx_xxxxyyzz, g_0_x_xxxx_xxxxyzz, g_0_x_xxxx_xxxxyzzz, g_0_x_xxxx_xxxxzzz, g_0_x_xxxx_xxxxzzzz, g_0_x_xxxx_xxxyyyy, g_0_x_xxxx_xxxyyyyz, g_0_x_xxxx_xxxyyyz, g_0_x_xxxx_xxxyyyzz, g_0_x_xxxx_xxxyyzz, g_0_x_xxxx_xxxyyzzz, g_0_x_xxxx_xxxyzzz, g_0_x_xxxx_xxxyzzzz, g_0_x_xxxx_xxxzzzz, g_0_x_xxxx_xxxzzzzz, g_0_x_xxxx_xxyyyyy, g_0_x_xxxx_xxyyyyyz, g_0_x_xxxx_xxyyyyz, g_0_x_xxxx_xxyyyyzz, g_0_x_xxxx_xxyyyzz, g_0_x_xxxx_xxyyyzzz, g_0_x_xxxx_xxyyzzz, g_0_x_xxxx_xxyyzzzz, g_0_x_xxxx_xxyzzzz, g_0_x_xxxx_xxyzzzzz, g_0_x_xxxx_xxzzzzz, g_0_x_xxxx_xxzzzzzz, g_0_x_xxxx_xyyyyyy, g_0_x_xxxx_xyyyyyyz, g_0_x_xxxx_xyyyyyz, g_0_x_xxxx_xyyyyyzz, g_0_x_xxxx_xyyyyzz, g_0_x_xxxx_xyyyyzzz, g_0_x_xxxx_xyyyzzz, g_0_x_xxxx_xyyyzzzz, g_0_x_xxxx_xyyzzzz, g_0_x_xxxx_xyyzzzzz, g_0_x_xxxx_xyzzzzz, g_0_x_xxxx_xyzzzzzz, g_0_x_xxxx_xzzzzzz, g_0_x_xxxx_xzzzzzzz, g_0_x_xxxx_yyyyyyy, g_0_x_xxxx_yyyyyyyz, g_0_x_xxxx_yyyyyyz, g_0_x_xxxx_yyyyyyzz, g_0_x_xxxx_yyyyyzz, g_0_x_xxxx_yyyyyzzz, g_0_x_xxxx_yyyyzzz, g_0_x_xxxx_yyyyzzzz, g_0_x_xxxx_yyyzzzz, g_0_x_xxxx_yyyzzzzz, g_0_x_xxxx_yyzzzzz, g_0_x_xxxx_yyzzzzzz, g_0_x_xxxx_yzzzzzz, g_0_x_xxxx_yzzzzzzz, g_0_x_xxxx_zzzzzzz, g_0_x_xxxx_zzzzzzzz, g_0_x_xxxxz_xxxxxxx, g_0_x_xxxxz_xxxxxxy, g_0_x_xxxxz_xxxxxxz, g_0_x_xxxxz_xxxxxyy, g_0_x_xxxxz_xxxxxyz, g_0_x_xxxxz_xxxxxzz, g_0_x_xxxxz_xxxxyyy, g_0_x_xxxxz_xxxxyyz, g_0_x_xxxxz_xxxxyzz, g_0_x_xxxxz_xxxxzzz, g_0_x_xxxxz_xxxyyyy, g_0_x_xxxxz_xxxyyyz, g_0_x_xxxxz_xxxyyzz, g_0_x_xxxxz_xxxyzzz, g_0_x_xxxxz_xxxzzzz, g_0_x_xxxxz_xxyyyyy, g_0_x_xxxxz_xxyyyyz, g_0_x_xxxxz_xxyyyzz, g_0_x_xxxxz_xxyyzzz, g_0_x_xxxxz_xxyzzzz, g_0_x_xxxxz_xxzzzzz, g_0_x_xxxxz_xyyyyyy, g_0_x_xxxxz_xyyyyyz, g_0_x_xxxxz_xyyyyzz, g_0_x_xxxxz_xyyyzzz, g_0_x_xxxxz_xyyzzzz, g_0_x_xxxxz_xyzzzzz, g_0_x_xxxxz_xzzzzzz, g_0_x_xxxxz_yyyyyyy, g_0_x_xxxxz_yyyyyyz, g_0_x_xxxxz_yyyyyzz, g_0_x_xxxxz_yyyyzzz, g_0_x_xxxxz_yyyzzzz, g_0_x_xxxxz_yyzzzzz, g_0_x_xxxxz_yzzzzzz, g_0_x_xxxxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxz_xxxxxxx[k] = -g_0_x_xxxx_xxxxxxx[k] * ab_z + g_0_x_xxxx_xxxxxxxz[k];

                g_0_x_xxxxz_xxxxxxy[k] = -g_0_x_xxxx_xxxxxxy[k] * ab_z + g_0_x_xxxx_xxxxxxyz[k];

                g_0_x_xxxxz_xxxxxxz[k] = -g_0_x_xxxx_xxxxxxz[k] * ab_z + g_0_x_xxxx_xxxxxxzz[k];

                g_0_x_xxxxz_xxxxxyy[k] = -g_0_x_xxxx_xxxxxyy[k] * ab_z + g_0_x_xxxx_xxxxxyyz[k];

                g_0_x_xxxxz_xxxxxyz[k] = -g_0_x_xxxx_xxxxxyz[k] * ab_z + g_0_x_xxxx_xxxxxyzz[k];

                g_0_x_xxxxz_xxxxxzz[k] = -g_0_x_xxxx_xxxxxzz[k] * ab_z + g_0_x_xxxx_xxxxxzzz[k];

                g_0_x_xxxxz_xxxxyyy[k] = -g_0_x_xxxx_xxxxyyy[k] * ab_z + g_0_x_xxxx_xxxxyyyz[k];

                g_0_x_xxxxz_xxxxyyz[k] = -g_0_x_xxxx_xxxxyyz[k] * ab_z + g_0_x_xxxx_xxxxyyzz[k];

                g_0_x_xxxxz_xxxxyzz[k] = -g_0_x_xxxx_xxxxyzz[k] * ab_z + g_0_x_xxxx_xxxxyzzz[k];

                g_0_x_xxxxz_xxxxzzz[k] = -g_0_x_xxxx_xxxxzzz[k] * ab_z + g_0_x_xxxx_xxxxzzzz[k];

                g_0_x_xxxxz_xxxyyyy[k] = -g_0_x_xxxx_xxxyyyy[k] * ab_z + g_0_x_xxxx_xxxyyyyz[k];

                g_0_x_xxxxz_xxxyyyz[k] = -g_0_x_xxxx_xxxyyyz[k] * ab_z + g_0_x_xxxx_xxxyyyzz[k];

                g_0_x_xxxxz_xxxyyzz[k] = -g_0_x_xxxx_xxxyyzz[k] * ab_z + g_0_x_xxxx_xxxyyzzz[k];

                g_0_x_xxxxz_xxxyzzz[k] = -g_0_x_xxxx_xxxyzzz[k] * ab_z + g_0_x_xxxx_xxxyzzzz[k];

                g_0_x_xxxxz_xxxzzzz[k] = -g_0_x_xxxx_xxxzzzz[k] * ab_z + g_0_x_xxxx_xxxzzzzz[k];

                g_0_x_xxxxz_xxyyyyy[k] = -g_0_x_xxxx_xxyyyyy[k] * ab_z + g_0_x_xxxx_xxyyyyyz[k];

                g_0_x_xxxxz_xxyyyyz[k] = -g_0_x_xxxx_xxyyyyz[k] * ab_z + g_0_x_xxxx_xxyyyyzz[k];

                g_0_x_xxxxz_xxyyyzz[k] = -g_0_x_xxxx_xxyyyzz[k] * ab_z + g_0_x_xxxx_xxyyyzzz[k];

                g_0_x_xxxxz_xxyyzzz[k] = -g_0_x_xxxx_xxyyzzz[k] * ab_z + g_0_x_xxxx_xxyyzzzz[k];

                g_0_x_xxxxz_xxyzzzz[k] = -g_0_x_xxxx_xxyzzzz[k] * ab_z + g_0_x_xxxx_xxyzzzzz[k];

                g_0_x_xxxxz_xxzzzzz[k] = -g_0_x_xxxx_xxzzzzz[k] * ab_z + g_0_x_xxxx_xxzzzzzz[k];

                g_0_x_xxxxz_xyyyyyy[k] = -g_0_x_xxxx_xyyyyyy[k] * ab_z + g_0_x_xxxx_xyyyyyyz[k];

                g_0_x_xxxxz_xyyyyyz[k] = -g_0_x_xxxx_xyyyyyz[k] * ab_z + g_0_x_xxxx_xyyyyyzz[k];

                g_0_x_xxxxz_xyyyyzz[k] = -g_0_x_xxxx_xyyyyzz[k] * ab_z + g_0_x_xxxx_xyyyyzzz[k];

                g_0_x_xxxxz_xyyyzzz[k] = -g_0_x_xxxx_xyyyzzz[k] * ab_z + g_0_x_xxxx_xyyyzzzz[k];

                g_0_x_xxxxz_xyyzzzz[k] = -g_0_x_xxxx_xyyzzzz[k] * ab_z + g_0_x_xxxx_xyyzzzzz[k];

                g_0_x_xxxxz_xyzzzzz[k] = -g_0_x_xxxx_xyzzzzz[k] * ab_z + g_0_x_xxxx_xyzzzzzz[k];

                g_0_x_xxxxz_xzzzzzz[k] = -g_0_x_xxxx_xzzzzzz[k] * ab_z + g_0_x_xxxx_xzzzzzzz[k];

                g_0_x_xxxxz_yyyyyyy[k] = -g_0_x_xxxx_yyyyyyy[k] * ab_z + g_0_x_xxxx_yyyyyyyz[k];

                g_0_x_xxxxz_yyyyyyz[k] = -g_0_x_xxxx_yyyyyyz[k] * ab_z + g_0_x_xxxx_yyyyyyzz[k];

                g_0_x_xxxxz_yyyyyzz[k] = -g_0_x_xxxx_yyyyyzz[k] * ab_z + g_0_x_xxxx_yyyyyzzz[k];

                g_0_x_xxxxz_yyyyzzz[k] = -g_0_x_xxxx_yyyyzzz[k] * ab_z + g_0_x_xxxx_yyyyzzzz[k];

                g_0_x_xxxxz_yyyzzzz[k] = -g_0_x_xxxx_yyyzzzz[k] * ab_z + g_0_x_xxxx_yyyzzzzz[k];

                g_0_x_xxxxz_yyzzzzz[k] = -g_0_x_xxxx_yyzzzzz[k] * ab_z + g_0_x_xxxx_yyzzzzzz[k];

                g_0_x_xxxxz_yzzzzzz[k] = -g_0_x_xxxx_yzzzzzz[k] * ab_z + g_0_x_xxxx_yzzzzzzz[k];

                g_0_x_xxxxz_zzzzzzz[k] = -g_0_x_xxxx_zzzzzzz[k] * ab_z + g_0_x_xxxx_zzzzzzzz[k];
            }

            /// Set up 108-144 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxxyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxxyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxxyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxy_xxxxxxx, g_0_x_xxxy_xxxxxxxy, g_0_x_xxxy_xxxxxxy, g_0_x_xxxy_xxxxxxyy, g_0_x_xxxy_xxxxxxyz, g_0_x_xxxy_xxxxxxz, g_0_x_xxxy_xxxxxyy, g_0_x_xxxy_xxxxxyyy, g_0_x_xxxy_xxxxxyyz, g_0_x_xxxy_xxxxxyz, g_0_x_xxxy_xxxxxyzz, g_0_x_xxxy_xxxxxzz, g_0_x_xxxy_xxxxyyy, g_0_x_xxxy_xxxxyyyy, g_0_x_xxxy_xxxxyyyz, g_0_x_xxxy_xxxxyyz, g_0_x_xxxy_xxxxyyzz, g_0_x_xxxy_xxxxyzz, g_0_x_xxxy_xxxxyzzz, g_0_x_xxxy_xxxxzzz, g_0_x_xxxy_xxxyyyy, g_0_x_xxxy_xxxyyyyy, g_0_x_xxxy_xxxyyyyz, g_0_x_xxxy_xxxyyyz, g_0_x_xxxy_xxxyyyzz, g_0_x_xxxy_xxxyyzz, g_0_x_xxxy_xxxyyzzz, g_0_x_xxxy_xxxyzzz, g_0_x_xxxy_xxxyzzzz, g_0_x_xxxy_xxxzzzz, g_0_x_xxxy_xxyyyyy, g_0_x_xxxy_xxyyyyyy, g_0_x_xxxy_xxyyyyyz, g_0_x_xxxy_xxyyyyz, g_0_x_xxxy_xxyyyyzz, g_0_x_xxxy_xxyyyzz, g_0_x_xxxy_xxyyyzzz, g_0_x_xxxy_xxyyzzz, g_0_x_xxxy_xxyyzzzz, g_0_x_xxxy_xxyzzzz, g_0_x_xxxy_xxyzzzzz, g_0_x_xxxy_xxzzzzz, g_0_x_xxxy_xyyyyyy, g_0_x_xxxy_xyyyyyyy, g_0_x_xxxy_xyyyyyyz, g_0_x_xxxy_xyyyyyz, g_0_x_xxxy_xyyyyyzz, g_0_x_xxxy_xyyyyzz, g_0_x_xxxy_xyyyyzzz, g_0_x_xxxy_xyyyzzz, g_0_x_xxxy_xyyyzzzz, g_0_x_xxxy_xyyzzzz, g_0_x_xxxy_xyyzzzzz, g_0_x_xxxy_xyzzzzz, g_0_x_xxxy_xyzzzzzz, g_0_x_xxxy_xzzzzzz, g_0_x_xxxy_yyyyyyy, g_0_x_xxxy_yyyyyyyy, g_0_x_xxxy_yyyyyyyz, g_0_x_xxxy_yyyyyyz, g_0_x_xxxy_yyyyyyzz, g_0_x_xxxy_yyyyyzz, g_0_x_xxxy_yyyyyzzz, g_0_x_xxxy_yyyyzzz, g_0_x_xxxy_yyyyzzzz, g_0_x_xxxy_yyyzzzz, g_0_x_xxxy_yyyzzzzz, g_0_x_xxxy_yyzzzzz, g_0_x_xxxy_yyzzzzzz, g_0_x_xxxy_yzzzzzz, g_0_x_xxxy_yzzzzzzz, g_0_x_xxxy_zzzzzzz, g_0_x_xxxyy_xxxxxxx, g_0_x_xxxyy_xxxxxxy, g_0_x_xxxyy_xxxxxxz, g_0_x_xxxyy_xxxxxyy, g_0_x_xxxyy_xxxxxyz, g_0_x_xxxyy_xxxxxzz, g_0_x_xxxyy_xxxxyyy, g_0_x_xxxyy_xxxxyyz, g_0_x_xxxyy_xxxxyzz, g_0_x_xxxyy_xxxxzzz, g_0_x_xxxyy_xxxyyyy, g_0_x_xxxyy_xxxyyyz, g_0_x_xxxyy_xxxyyzz, g_0_x_xxxyy_xxxyzzz, g_0_x_xxxyy_xxxzzzz, g_0_x_xxxyy_xxyyyyy, g_0_x_xxxyy_xxyyyyz, g_0_x_xxxyy_xxyyyzz, g_0_x_xxxyy_xxyyzzz, g_0_x_xxxyy_xxyzzzz, g_0_x_xxxyy_xxzzzzz, g_0_x_xxxyy_xyyyyyy, g_0_x_xxxyy_xyyyyyz, g_0_x_xxxyy_xyyyyzz, g_0_x_xxxyy_xyyyzzz, g_0_x_xxxyy_xyyzzzz, g_0_x_xxxyy_xyzzzzz, g_0_x_xxxyy_xzzzzzz, g_0_x_xxxyy_yyyyyyy, g_0_x_xxxyy_yyyyyyz, g_0_x_xxxyy_yyyyyzz, g_0_x_xxxyy_yyyyzzz, g_0_x_xxxyy_yyyzzzz, g_0_x_xxxyy_yyzzzzz, g_0_x_xxxyy_yzzzzzz, g_0_x_xxxyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyy_xxxxxxx[k] = -g_0_x_xxxy_xxxxxxx[k] * ab_y + g_0_x_xxxy_xxxxxxxy[k];

                g_0_x_xxxyy_xxxxxxy[k] = -g_0_x_xxxy_xxxxxxy[k] * ab_y + g_0_x_xxxy_xxxxxxyy[k];

                g_0_x_xxxyy_xxxxxxz[k] = -g_0_x_xxxy_xxxxxxz[k] * ab_y + g_0_x_xxxy_xxxxxxyz[k];

                g_0_x_xxxyy_xxxxxyy[k] = -g_0_x_xxxy_xxxxxyy[k] * ab_y + g_0_x_xxxy_xxxxxyyy[k];

                g_0_x_xxxyy_xxxxxyz[k] = -g_0_x_xxxy_xxxxxyz[k] * ab_y + g_0_x_xxxy_xxxxxyyz[k];

                g_0_x_xxxyy_xxxxxzz[k] = -g_0_x_xxxy_xxxxxzz[k] * ab_y + g_0_x_xxxy_xxxxxyzz[k];

                g_0_x_xxxyy_xxxxyyy[k] = -g_0_x_xxxy_xxxxyyy[k] * ab_y + g_0_x_xxxy_xxxxyyyy[k];

                g_0_x_xxxyy_xxxxyyz[k] = -g_0_x_xxxy_xxxxyyz[k] * ab_y + g_0_x_xxxy_xxxxyyyz[k];

                g_0_x_xxxyy_xxxxyzz[k] = -g_0_x_xxxy_xxxxyzz[k] * ab_y + g_0_x_xxxy_xxxxyyzz[k];

                g_0_x_xxxyy_xxxxzzz[k] = -g_0_x_xxxy_xxxxzzz[k] * ab_y + g_0_x_xxxy_xxxxyzzz[k];

                g_0_x_xxxyy_xxxyyyy[k] = -g_0_x_xxxy_xxxyyyy[k] * ab_y + g_0_x_xxxy_xxxyyyyy[k];

                g_0_x_xxxyy_xxxyyyz[k] = -g_0_x_xxxy_xxxyyyz[k] * ab_y + g_0_x_xxxy_xxxyyyyz[k];

                g_0_x_xxxyy_xxxyyzz[k] = -g_0_x_xxxy_xxxyyzz[k] * ab_y + g_0_x_xxxy_xxxyyyzz[k];

                g_0_x_xxxyy_xxxyzzz[k] = -g_0_x_xxxy_xxxyzzz[k] * ab_y + g_0_x_xxxy_xxxyyzzz[k];

                g_0_x_xxxyy_xxxzzzz[k] = -g_0_x_xxxy_xxxzzzz[k] * ab_y + g_0_x_xxxy_xxxyzzzz[k];

                g_0_x_xxxyy_xxyyyyy[k] = -g_0_x_xxxy_xxyyyyy[k] * ab_y + g_0_x_xxxy_xxyyyyyy[k];

                g_0_x_xxxyy_xxyyyyz[k] = -g_0_x_xxxy_xxyyyyz[k] * ab_y + g_0_x_xxxy_xxyyyyyz[k];

                g_0_x_xxxyy_xxyyyzz[k] = -g_0_x_xxxy_xxyyyzz[k] * ab_y + g_0_x_xxxy_xxyyyyzz[k];

                g_0_x_xxxyy_xxyyzzz[k] = -g_0_x_xxxy_xxyyzzz[k] * ab_y + g_0_x_xxxy_xxyyyzzz[k];

                g_0_x_xxxyy_xxyzzzz[k] = -g_0_x_xxxy_xxyzzzz[k] * ab_y + g_0_x_xxxy_xxyyzzzz[k];

                g_0_x_xxxyy_xxzzzzz[k] = -g_0_x_xxxy_xxzzzzz[k] * ab_y + g_0_x_xxxy_xxyzzzzz[k];

                g_0_x_xxxyy_xyyyyyy[k] = -g_0_x_xxxy_xyyyyyy[k] * ab_y + g_0_x_xxxy_xyyyyyyy[k];

                g_0_x_xxxyy_xyyyyyz[k] = -g_0_x_xxxy_xyyyyyz[k] * ab_y + g_0_x_xxxy_xyyyyyyz[k];

                g_0_x_xxxyy_xyyyyzz[k] = -g_0_x_xxxy_xyyyyzz[k] * ab_y + g_0_x_xxxy_xyyyyyzz[k];

                g_0_x_xxxyy_xyyyzzz[k] = -g_0_x_xxxy_xyyyzzz[k] * ab_y + g_0_x_xxxy_xyyyyzzz[k];

                g_0_x_xxxyy_xyyzzzz[k] = -g_0_x_xxxy_xyyzzzz[k] * ab_y + g_0_x_xxxy_xyyyzzzz[k];

                g_0_x_xxxyy_xyzzzzz[k] = -g_0_x_xxxy_xyzzzzz[k] * ab_y + g_0_x_xxxy_xyyzzzzz[k];

                g_0_x_xxxyy_xzzzzzz[k] = -g_0_x_xxxy_xzzzzzz[k] * ab_y + g_0_x_xxxy_xyzzzzzz[k];

                g_0_x_xxxyy_yyyyyyy[k] = -g_0_x_xxxy_yyyyyyy[k] * ab_y + g_0_x_xxxy_yyyyyyyy[k];

                g_0_x_xxxyy_yyyyyyz[k] = -g_0_x_xxxy_yyyyyyz[k] * ab_y + g_0_x_xxxy_yyyyyyyz[k];

                g_0_x_xxxyy_yyyyyzz[k] = -g_0_x_xxxy_yyyyyzz[k] * ab_y + g_0_x_xxxy_yyyyyyzz[k];

                g_0_x_xxxyy_yyyyzzz[k] = -g_0_x_xxxy_yyyyzzz[k] * ab_y + g_0_x_xxxy_yyyyyzzz[k];

                g_0_x_xxxyy_yyyzzzz[k] = -g_0_x_xxxy_yyyzzzz[k] * ab_y + g_0_x_xxxy_yyyyzzzz[k];

                g_0_x_xxxyy_yyzzzzz[k] = -g_0_x_xxxy_yyzzzzz[k] * ab_y + g_0_x_xxxy_yyyzzzzz[k];

                g_0_x_xxxyy_yzzzzzz[k] = -g_0_x_xxxy_yzzzzzz[k] * ab_y + g_0_x_xxxy_yyzzzzzz[k];

                g_0_x_xxxyy_zzzzzzz[k] = -g_0_x_xxxy_zzzzzzz[k] * ab_y + g_0_x_xxxy_yzzzzzzz[k];
            }

            /// Set up 144-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxxyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxxyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxxyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyz_xxxxxxx, g_0_x_xxxyz_xxxxxxy, g_0_x_xxxyz_xxxxxxz, g_0_x_xxxyz_xxxxxyy, g_0_x_xxxyz_xxxxxyz, g_0_x_xxxyz_xxxxxzz, g_0_x_xxxyz_xxxxyyy, g_0_x_xxxyz_xxxxyyz, g_0_x_xxxyz_xxxxyzz, g_0_x_xxxyz_xxxxzzz, g_0_x_xxxyz_xxxyyyy, g_0_x_xxxyz_xxxyyyz, g_0_x_xxxyz_xxxyyzz, g_0_x_xxxyz_xxxyzzz, g_0_x_xxxyz_xxxzzzz, g_0_x_xxxyz_xxyyyyy, g_0_x_xxxyz_xxyyyyz, g_0_x_xxxyz_xxyyyzz, g_0_x_xxxyz_xxyyzzz, g_0_x_xxxyz_xxyzzzz, g_0_x_xxxyz_xxzzzzz, g_0_x_xxxyz_xyyyyyy, g_0_x_xxxyz_xyyyyyz, g_0_x_xxxyz_xyyyyzz, g_0_x_xxxyz_xyyyzzz, g_0_x_xxxyz_xyyzzzz, g_0_x_xxxyz_xyzzzzz, g_0_x_xxxyz_xzzzzzz, g_0_x_xxxyz_yyyyyyy, g_0_x_xxxyz_yyyyyyz, g_0_x_xxxyz_yyyyyzz, g_0_x_xxxyz_yyyyzzz, g_0_x_xxxyz_yyyzzzz, g_0_x_xxxyz_yyzzzzz, g_0_x_xxxyz_yzzzzzz, g_0_x_xxxyz_zzzzzzz, g_0_x_xxxz_xxxxxxx, g_0_x_xxxz_xxxxxxxy, g_0_x_xxxz_xxxxxxy, g_0_x_xxxz_xxxxxxyy, g_0_x_xxxz_xxxxxxyz, g_0_x_xxxz_xxxxxxz, g_0_x_xxxz_xxxxxyy, g_0_x_xxxz_xxxxxyyy, g_0_x_xxxz_xxxxxyyz, g_0_x_xxxz_xxxxxyz, g_0_x_xxxz_xxxxxyzz, g_0_x_xxxz_xxxxxzz, g_0_x_xxxz_xxxxyyy, g_0_x_xxxz_xxxxyyyy, g_0_x_xxxz_xxxxyyyz, g_0_x_xxxz_xxxxyyz, g_0_x_xxxz_xxxxyyzz, g_0_x_xxxz_xxxxyzz, g_0_x_xxxz_xxxxyzzz, g_0_x_xxxz_xxxxzzz, g_0_x_xxxz_xxxyyyy, g_0_x_xxxz_xxxyyyyy, g_0_x_xxxz_xxxyyyyz, g_0_x_xxxz_xxxyyyz, g_0_x_xxxz_xxxyyyzz, g_0_x_xxxz_xxxyyzz, g_0_x_xxxz_xxxyyzzz, g_0_x_xxxz_xxxyzzz, g_0_x_xxxz_xxxyzzzz, g_0_x_xxxz_xxxzzzz, g_0_x_xxxz_xxyyyyy, g_0_x_xxxz_xxyyyyyy, g_0_x_xxxz_xxyyyyyz, g_0_x_xxxz_xxyyyyz, g_0_x_xxxz_xxyyyyzz, g_0_x_xxxz_xxyyyzz, g_0_x_xxxz_xxyyyzzz, g_0_x_xxxz_xxyyzzz, g_0_x_xxxz_xxyyzzzz, g_0_x_xxxz_xxyzzzz, g_0_x_xxxz_xxyzzzzz, g_0_x_xxxz_xxzzzzz, g_0_x_xxxz_xyyyyyy, g_0_x_xxxz_xyyyyyyy, g_0_x_xxxz_xyyyyyyz, g_0_x_xxxz_xyyyyyz, g_0_x_xxxz_xyyyyyzz, g_0_x_xxxz_xyyyyzz, g_0_x_xxxz_xyyyyzzz, g_0_x_xxxz_xyyyzzz, g_0_x_xxxz_xyyyzzzz, g_0_x_xxxz_xyyzzzz, g_0_x_xxxz_xyyzzzzz, g_0_x_xxxz_xyzzzzz, g_0_x_xxxz_xyzzzzzz, g_0_x_xxxz_xzzzzzz, g_0_x_xxxz_yyyyyyy, g_0_x_xxxz_yyyyyyyy, g_0_x_xxxz_yyyyyyyz, g_0_x_xxxz_yyyyyyz, g_0_x_xxxz_yyyyyyzz, g_0_x_xxxz_yyyyyzz, g_0_x_xxxz_yyyyyzzz, g_0_x_xxxz_yyyyzzz, g_0_x_xxxz_yyyyzzzz, g_0_x_xxxz_yyyzzzz, g_0_x_xxxz_yyyzzzzz, g_0_x_xxxz_yyzzzzz, g_0_x_xxxz_yyzzzzzz, g_0_x_xxxz_yzzzzzz, g_0_x_xxxz_yzzzzzzz, g_0_x_xxxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyz_xxxxxxx[k] = -g_0_x_xxxz_xxxxxxx[k] * ab_y + g_0_x_xxxz_xxxxxxxy[k];

                g_0_x_xxxyz_xxxxxxy[k] = -g_0_x_xxxz_xxxxxxy[k] * ab_y + g_0_x_xxxz_xxxxxxyy[k];

                g_0_x_xxxyz_xxxxxxz[k] = -g_0_x_xxxz_xxxxxxz[k] * ab_y + g_0_x_xxxz_xxxxxxyz[k];

                g_0_x_xxxyz_xxxxxyy[k] = -g_0_x_xxxz_xxxxxyy[k] * ab_y + g_0_x_xxxz_xxxxxyyy[k];

                g_0_x_xxxyz_xxxxxyz[k] = -g_0_x_xxxz_xxxxxyz[k] * ab_y + g_0_x_xxxz_xxxxxyyz[k];

                g_0_x_xxxyz_xxxxxzz[k] = -g_0_x_xxxz_xxxxxzz[k] * ab_y + g_0_x_xxxz_xxxxxyzz[k];

                g_0_x_xxxyz_xxxxyyy[k] = -g_0_x_xxxz_xxxxyyy[k] * ab_y + g_0_x_xxxz_xxxxyyyy[k];

                g_0_x_xxxyz_xxxxyyz[k] = -g_0_x_xxxz_xxxxyyz[k] * ab_y + g_0_x_xxxz_xxxxyyyz[k];

                g_0_x_xxxyz_xxxxyzz[k] = -g_0_x_xxxz_xxxxyzz[k] * ab_y + g_0_x_xxxz_xxxxyyzz[k];

                g_0_x_xxxyz_xxxxzzz[k] = -g_0_x_xxxz_xxxxzzz[k] * ab_y + g_0_x_xxxz_xxxxyzzz[k];

                g_0_x_xxxyz_xxxyyyy[k] = -g_0_x_xxxz_xxxyyyy[k] * ab_y + g_0_x_xxxz_xxxyyyyy[k];

                g_0_x_xxxyz_xxxyyyz[k] = -g_0_x_xxxz_xxxyyyz[k] * ab_y + g_0_x_xxxz_xxxyyyyz[k];

                g_0_x_xxxyz_xxxyyzz[k] = -g_0_x_xxxz_xxxyyzz[k] * ab_y + g_0_x_xxxz_xxxyyyzz[k];

                g_0_x_xxxyz_xxxyzzz[k] = -g_0_x_xxxz_xxxyzzz[k] * ab_y + g_0_x_xxxz_xxxyyzzz[k];

                g_0_x_xxxyz_xxxzzzz[k] = -g_0_x_xxxz_xxxzzzz[k] * ab_y + g_0_x_xxxz_xxxyzzzz[k];

                g_0_x_xxxyz_xxyyyyy[k] = -g_0_x_xxxz_xxyyyyy[k] * ab_y + g_0_x_xxxz_xxyyyyyy[k];

                g_0_x_xxxyz_xxyyyyz[k] = -g_0_x_xxxz_xxyyyyz[k] * ab_y + g_0_x_xxxz_xxyyyyyz[k];

                g_0_x_xxxyz_xxyyyzz[k] = -g_0_x_xxxz_xxyyyzz[k] * ab_y + g_0_x_xxxz_xxyyyyzz[k];

                g_0_x_xxxyz_xxyyzzz[k] = -g_0_x_xxxz_xxyyzzz[k] * ab_y + g_0_x_xxxz_xxyyyzzz[k];

                g_0_x_xxxyz_xxyzzzz[k] = -g_0_x_xxxz_xxyzzzz[k] * ab_y + g_0_x_xxxz_xxyyzzzz[k];

                g_0_x_xxxyz_xxzzzzz[k] = -g_0_x_xxxz_xxzzzzz[k] * ab_y + g_0_x_xxxz_xxyzzzzz[k];

                g_0_x_xxxyz_xyyyyyy[k] = -g_0_x_xxxz_xyyyyyy[k] * ab_y + g_0_x_xxxz_xyyyyyyy[k];

                g_0_x_xxxyz_xyyyyyz[k] = -g_0_x_xxxz_xyyyyyz[k] * ab_y + g_0_x_xxxz_xyyyyyyz[k];

                g_0_x_xxxyz_xyyyyzz[k] = -g_0_x_xxxz_xyyyyzz[k] * ab_y + g_0_x_xxxz_xyyyyyzz[k];

                g_0_x_xxxyz_xyyyzzz[k] = -g_0_x_xxxz_xyyyzzz[k] * ab_y + g_0_x_xxxz_xyyyyzzz[k];

                g_0_x_xxxyz_xyyzzzz[k] = -g_0_x_xxxz_xyyzzzz[k] * ab_y + g_0_x_xxxz_xyyyzzzz[k];

                g_0_x_xxxyz_xyzzzzz[k] = -g_0_x_xxxz_xyzzzzz[k] * ab_y + g_0_x_xxxz_xyyzzzzz[k];

                g_0_x_xxxyz_xzzzzzz[k] = -g_0_x_xxxz_xzzzzzz[k] * ab_y + g_0_x_xxxz_xyzzzzzz[k];

                g_0_x_xxxyz_yyyyyyy[k] = -g_0_x_xxxz_yyyyyyy[k] * ab_y + g_0_x_xxxz_yyyyyyyy[k];

                g_0_x_xxxyz_yyyyyyz[k] = -g_0_x_xxxz_yyyyyyz[k] * ab_y + g_0_x_xxxz_yyyyyyyz[k];

                g_0_x_xxxyz_yyyyyzz[k] = -g_0_x_xxxz_yyyyyzz[k] * ab_y + g_0_x_xxxz_yyyyyyzz[k];

                g_0_x_xxxyz_yyyyzzz[k] = -g_0_x_xxxz_yyyyzzz[k] * ab_y + g_0_x_xxxz_yyyyyzzz[k];

                g_0_x_xxxyz_yyyzzzz[k] = -g_0_x_xxxz_yyyzzzz[k] * ab_y + g_0_x_xxxz_yyyyzzzz[k];

                g_0_x_xxxyz_yyzzzzz[k] = -g_0_x_xxxz_yyzzzzz[k] * ab_y + g_0_x_xxxz_yyyzzzzz[k];

                g_0_x_xxxyz_yzzzzzz[k] = -g_0_x_xxxz_yzzzzzz[k] * ab_y + g_0_x_xxxz_yyzzzzzz[k];

                g_0_x_xxxyz_zzzzzzz[k] = -g_0_x_xxxz_zzzzzzz[k] * ab_y + g_0_x_xxxz_yzzzzzzz[k];
            }

            /// Set up 180-216 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxxzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxxzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxxzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxz_xxxxxxx, g_0_x_xxxz_xxxxxxxz, g_0_x_xxxz_xxxxxxy, g_0_x_xxxz_xxxxxxyz, g_0_x_xxxz_xxxxxxz, g_0_x_xxxz_xxxxxxzz, g_0_x_xxxz_xxxxxyy, g_0_x_xxxz_xxxxxyyz, g_0_x_xxxz_xxxxxyz, g_0_x_xxxz_xxxxxyzz, g_0_x_xxxz_xxxxxzz, g_0_x_xxxz_xxxxxzzz, g_0_x_xxxz_xxxxyyy, g_0_x_xxxz_xxxxyyyz, g_0_x_xxxz_xxxxyyz, g_0_x_xxxz_xxxxyyzz, g_0_x_xxxz_xxxxyzz, g_0_x_xxxz_xxxxyzzz, g_0_x_xxxz_xxxxzzz, g_0_x_xxxz_xxxxzzzz, g_0_x_xxxz_xxxyyyy, g_0_x_xxxz_xxxyyyyz, g_0_x_xxxz_xxxyyyz, g_0_x_xxxz_xxxyyyzz, g_0_x_xxxz_xxxyyzz, g_0_x_xxxz_xxxyyzzz, g_0_x_xxxz_xxxyzzz, g_0_x_xxxz_xxxyzzzz, g_0_x_xxxz_xxxzzzz, g_0_x_xxxz_xxxzzzzz, g_0_x_xxxz_xxyyyyy, g_0_x_xxxz_xxyyyyyz, g_0_x_xxxz_xxyyyyz, g_0_x_xxxz_xxyyyyzz, g_0_x_xxxz_xxyyyzz, g_0_x_xxxz_xxyyyzzz, g_0_x_xxxz_xxyyzzz, g_0_x_xxxz_xxyyzzzz, g_0_x_xxxz_xxyzzzz, g_0_x_xxxz_xxyzzzzz, g_0_x_xxxz_xxzzzzz, g_0_x_xxxz_xxzzzzzz, g_0_x_xxxz_xyyyyyy, g_0_x_xxxz_xyyyyyyz, g_0_x_xxxz_xyyyyyz, g_0_x_xxxz_xyyyyyzz, g_0_x_xxxz_xyyyyzz, g_0_x_xxxz_xyyyyzzz, g_0_x_xxxz_xyyyzzz, g_0_x_xxxz_xyyyzzzz, g_0_x_xxxz_xyyzzzz, g_0_x_xxxz_xyyzzzzz, g_0_x_xxxz_xyzzzzz, g_0_x_xxxz_xyzzzzzz, g_0_x_xxxz_xzzzzzz, g_0_x_xxxz_xzzzzzzz, g_0_x_xxxz_yyyyyyy, g_0_x_xxxz_yyyyyyyz, g_0_x_xxxz_yyyyyyz, g_0_x_xxxz_yyyyyyzz, g_0_x_xxxz_yyyyyzz, g_0_x_xxxz_yyyyyzzz, g_0_x_xxxz_yyyyzzz, g_0_x_xxxz_yyyyzzzz, g_0_x_xxxz_yyyzzzz, g_0_x_xxxz_yyyzzzzz, g_0_x_xxxz_yyzzzzz, g_0_x_xxxz_yyzzzzzz, g_0_x_xxxz_yzzzzzz, g_0_x_xxxz_yzzzzzzz, g_0_x_xxxz_zzzzzzz, g_0_x_xxxz_zzzzzzzz, g_0_x_xxxzz_xxxxxxx, g_0_x_xxxzz_xxxxxxy, g_0_x_xxxzz_xxxxxxz, g_0_x_xxxzz_xxxxxyy, g_0_x_xxxzz_xxxxxyz, g_0_x_xxxzz_xxxxxzz, g_0_x_xxxzz_xxxxyyy, g_0_x_xxxzz_xxxxyyz, g_0_x_xxxzz_xxxxyzz, g_0_x_xxxzz_xxxxzzz, g_0_x_xxxzz_xxxyyyy, g_0_x_xxxzz_xxxyyyz, g_0_x_xxxzz_xxxyyzz, g_0_x_xxxzz_xxxyzzz, g_0_x_xxxzz_xxxzzzz, g_0_x_xxxzz_xxyyyyy, g_0_x_xxxzz_xxyyyyz, g_0_x_xxxzz_xxyyyzz, g_0_x_xxxzz_xxyyzzz, g_0_x_xxxzz_xxyzzzz, g_0_x_xxxzz_xxzzzzz, g_0_x_xxxzz_xyyyyyy, g_0_x_xxxzz_xyyyyyz, g_0_x_xxxzz_xyyyyzz, g_0_x_xxxzz_xyyyzzz, g_0_x_xxxzz_xyyzzzz, g_0_x_xxxzz_xyzzzzz, g_0_x_xxxzz_xzzzzzz, g_0_x_xxxzz_yyyyyyy, g_0_x_xxxzz_yyyyyyz, g_0_x_xxxzz_yyyyyzz, g_0_x_xxxzz_yyyyzzz, g_0_x_xxxzz_yyyzzzz, g_0_x_xxxzz_yyzzzzz, g_0_x_xxxzz_yzzzzzz, g_0_x_xxxzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzz_xxxxxxx[k] = -g_0_x_xxxz_xxxxxxx[k] * ab_z + g_0_x_xxxz_xxxxxxxz[k];

                g_0_x_xxxzz_xxxxxxy[k] = -g_0_x_xxxz_xxxxxxy[k] * ab_z + g_0_x_xxxz_xxxxxxyz[k];

                g_0_x_xxxzz_xxxxxxz[k] = -g_0_x_xxxz_xxxxxxz[k] * ab_z + g_0_x_xxxz_xxxxxxzz[k];

                g_0_x_xxxzz_xxxxxyy[k] = -g_0_x_xxxz_xxxxxyy[k] * ab_z + g_0_x_xxxz_xxxxxyyz[k];

                g_0_x_xxxzz_xxxxxyz[k] = -g_0_x_xxxz_xxxxxyz[k] * ab_z + g_0_x_xxxz_xxxxxyzz[k];

                g_0_x_xxxzz_xxxxxzz[k] = -g_0_x_xxxz_xxxxxzz[k] * ab_z + g_0_x_xxxz_xxxxxzzz[k];

                g_0_x_xxxzz_xxxxyyy[k] = -g_0_x_xxxz_xxxxyyy[k] * ab_z + g_0_x_xxxz_xxxxyyyz[k];

                g_0_x_xxxzz_xxxxyyz[k] = -g_0_x_xxxz_xxxxyyz[k] * ab_z + g_0_x_xxxz_xxxxyyzz[k];

                g_0_x_xxxzz_xxxxyzz[k] = -g_0_x_xxxz_xxxxyzz[k] * ab_z + g_0_x_xxxz_xxxxyzzz[k];

                g_0_x_xxxzz_xxxxzzz[k] = -g_0_x_xxxz_xxxxzzz[k] * ab_z + g_0_x_xxxz_xxxxzzzz[k];

                g_0_x_xxxzz_xxxyyyy[k] = -g_0_x_xxxz_xxxyyyy[k] * ab_z + g_0_x_xxxz_xxxyyyyz[k];

                g_0_x_xxxzz_xxxyyyz[k] = -g_0_x_xxxz_xxxyyyz[k] * ab_z + g_0_x_xxxz_xxxyyyzz[k];

                g_0_x_xxxzz_xxxyyzz[k] = -g_0_x_xxxz_xxxyyzz[k] * ab_z + g_0_x_xxxz_xxxyyzzz[k];

                g_0_x_xxxzz_xxxyzzz[k] = -g_0_x_xxxz_xxxyzzz[k] * ab_z + g_0_x_xxxz_xxxyzzzz[k];

                g_0_x_xxxzz_xxxzzzz[k] = -g_0_x_xxxz_xxxzzzz[k] * ab_z + g_0_x_xxxz_xxxzzzzz[k];

                g_0_x_xxxzz_xxyyyyy[k] = -g_0_x_xxxz_xxyyyyy[k] * ab_z + g_0_x_xxxz_xxyyyyyz[k];

                g_0_x_xxxzz_xxyyyyz[k] = -g_0_x_xxxz_xxyyyyz[k] * ab_z + g_0_x_xxxz_xxyyyyzz[k];

                g_0_x_xxxzz_xxyyyzz[k] = -g_0_x_xxxz_xxyyyzz[k] * ab_z + g_0_x_xxxz_xxyyyzzz[k];

                g_0_x_xxxzz_xxyyzzz[k] = -g_0_x_xxxz_xxyyzzz[k] * ab_z + g_0_x_xxxz_xxyyzzzz[k];

                g_0_x_xxxzz_xxyzzzz[k] = -g_0_x_xxxz_xxyzzzz[k] * ab_z + g_0_x_xxxz_xxyzzzzz[k];

                g_0_x_xxxzz_xxzzzzz[k] = -g_0_x_xxxz_xxzzzzz[k] * ab_z + g_0_x_xxxz_xxzzzzzz[k];

                g_0_x_xxxzz_xyyyyyy[k] = -g_0_x_xxxz_xyyyyyy[k] * ab_z + g_0_x_xxxz_xyyyyyyz[k];

                g_0_x_xxxzz_xyyyyyz[k] = -g_0_x_xxxz_xyyyyyz[k] * ab_z + g_0_x_xxxz_xyyyyyzz[k];

                g_0_x_xxxzz_xyyyyzz[k] = -g_0_x_xxxz_xyyyyzz[k] * ab_z + g_0_x_xxxz_xyyyyzzz[k];

                g_0_x_xxxzz_xyyyzzz[k] = -g_0_x_xxxz_xyyyzzz[k] * ab_z + g_0_x_xxxz_xyyyzzzz[k];

                g_0_x_xxxzz_xyyzzzz[k] = -g_0_x_xxxz_xyyzzzz[k] * ab_z + g_0_x_xxxz_xyyzzzzz[k];

                g_0_x_xxxzz_xyzzzzz[k] = -g_0_x_xxxz_xyzzzzz[k] * ab_z + g_0_x_xxxz_xyzzzzzz[k];

                g_0_x_xxxzz_xzzzzzz[k] = -g_0_x_xxxz_xzzzzzz[k] * ab_z + g_0_x_xxxz_xzzzzzzz[k];

                g_0_x_xxxzz_yyyyyyy[k] = -g_0_x_xxxz_yyyyyyy[k] * ab_z + g_0_x_xxxz_yyyyyyyz[k];

                g_0_x_xxxzz_yyyyyyz[k] = -g_0_x_xxxz_yyyyyyz[k] * ab_z + g_0_x_xxxz_yyyyyyzz[k];

                g_0_x_xxxzz_yyyyyzz[k] = -g_0_x_xxxz_yyyyyzz[k] * ab_z + g_0_x_xxxz_yyyyyzzz[k];

                g_0_x_xxxzz_yyyyzzz[k] = -g_0_x_xxxz_yyyyzzz[k] * ab_z + g_0_x_xxxz_yyyyzzzz[k];

                g_0_x_xxxzz_yyyzzzz[k] = -g_0_x_xxxz_yyyzzzz[k] * ab_z + g_0_x_xxxz_yyyzzzzz[k];

                g_0_x_xxxzz_yyzzzzz[k] = -g_0_x_xxxz_yyzzzzz[k] * ab_z + g_0_x_xxxz_yyzzzzzz[k];

                g_0_x_xxxzz_yzzzzzz[k] = -g_0_x_xxxz_yzzzzzz[k] * ab_z + g_0_x_xxxz_yzzzzzzz[k];

                g_0_x_xxxzz_zzzzzzz[k] = -g_0_x_xxxz_zzzzzzz[k] * ab_z + g_0_x_xxxz_zzzzzzzz[k];
            }

            /// Set up 216-252 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xxyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xxyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xxyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyy_xxxxxxx, g_0_x_xxyy_xxxxxxxy, g_0_x_xxyy_xxxxxxy, g_0_x_xxyy_xxxxxxyy, g_0_x_xxyy_xxxxxxyz, g_0_x_xxyy_xxxxxxz, g_0_x_xxyy_xxxxxyy, g_0_x_xxyy_xxxxxyyy, g_0_x_xxyy_xxxxxyyz, g_0_x_xxyy_xxxxxyz, g_0_x_xxyy_xxxxxyzz, g_0_x_xxyy_xxxxxzz, g_0_x_xxyy_xxxxyyy, g_0_x_xxyy_xxxxyyyy, g_0_x_xxyy_xxxxyyyz, g_0_x_xxyy_xxxxyyz, g_0_x_xxyy_xxxxyyzz, g_0_x_xxyy_xxxxyzz, g_0_x_xxyy_xxxxyzzz, g_0_x_xxyy_xxxxzzz, g_0_x_xxyy_xxxyyyy, g_0_x_xxyy_xxxyyyyy, g_0_x_xxyy_xxxyyyyz, g_0_x_xxyy_xxxyyyz, g_0_x_xxyy_xxxyyyzz, g_0_x_xxyy_xxxyyzz, g_0_x_xxyy_xxxyyzzz, g_0_x_xxyy_xxxyzzz, g_0_x_xxyy_xxxyzzzz, g_0_x_xxyy_xxxzzzz, g_0_x_xxyy_xxyyyyy, g_0_x_xxyy_xxyyyyyy, g_0_x_xxyy_xxyyyyyz, g_0_x_xxyy_xxyyyyz, g_0_x_xxyy_xxyyyyzz, g_0_x_xxyy_xxyyyzz, g_0_x_xxyy_xxyyyzzz, g_0_x_xxyy_xxyyzzz, g_0_x_xxyy_xxyyzzzz, g_0_x_xxyy_xxyzzzz, g_0_x_xxyy_xxyzzzzz, g_0_x_xxyy_xxzzzzz, g_0_x_xxyy_xyyyyyy, g_0_x_xxyy_xyyyyyyy, g_0_x_xxyy_xyyyyyyz, g_0_x_xxyy_xyyyyyz, g_0_x_xxyy_xyyyyyzz, g_0_x_xxyy_xyyyyzz, g_0_x_xxyy_xyyyyzzz, g_0_x_xxyy_xyyyzzz, g_0_x_xxyy_xyyyzzzz, g_0_x_xxyy_xyyzzzz, g_0_x_xxyy_xyyzzzzz, g_0_x_xxyy_xyzzzzz, g_0_x_xxyy_xyzzzzzz, g_0_x_xxyy_xzzzzzz, g_0_x_xxyy_yyyyyyy, g_0_x_xxyy_yyyyyyyy, g_0_x_xxyy_yyyyyyyz, g_0_x_xxyy_yyyyyyz, g_0_x_xxyy_yyyyyyzz, g_0_x_xxyy_yyyyyzz, g_0_x_xxyy_yyyyyzzz, g_0_x_xxyy_yyyyzzz, g_0_x_xxyy_yyyyzzzz, g_0_x_xxyy_yyyzzzz, g_0_x_xxyy_yyyzzzzz, g_0_x_xxyy_yyzzzzz, g_0_x_xxyy_yyzzzzzz, g_0_x_xxyy_yzzzzzz, g_0_x_xxyy_yzzzzzzz, g_0_x_xxyy_zzzzzzz, g_0_x_xxyyy_xxxxxxx, g_0_x_xxyyy_xxxxxxy, g_0_x_xxyyy_xxxxxxz, g_0_x_xxyyy_xxxxxyy, g_0_x_xxyyy_xxxxxyz, g_0_x_xxyyy_xxxxxzz, g_0_x_xxyyy_xxxxyyy, g_0_x_xxyyy_xxxxyyz, g_0_x_xxyyy_xxxxyzz, g_0_x_xxyyy_xxxxzzz, g_0_x_xxyyy_xxxyyyy, g_0_x_xxyyy_xxxyyyz, g_0_x_xxyyy_xxxyyzz, g_0_x_xxyyy_xxxyzzz, g_0_x_xxyyy_xxxzzzz, g_0_x_xxyyy_xxyyyyy, g_0_x_xxyyy_xxyyyyz, g_0_x_xxyyy_xxyyyzz, g_0_x_xxyyy_xxyyzzz, g_0_x_xxyyy_xxyzzzz, g_0_x_xxyyy_xxzzzzz, g_0_x_xxyyy_xyyyyyy, g_0_x_xxyyy_xyyyyyz, g_0_x_xxyyy_xyyyyzz, g_0_x_xxyyy_xyyyzzz, g_0_x_xxyyy_xyyzzzz, g_0_x_xxyyy_xyzzzzz, g_0_x_xxyyy_xzzzzzz, g_0_x_xxyyy_yyyyyyy, g_0_x_xxyyy_yyyyyyz, g_0_x_xxyyy_yyyyyzz, g_0_x_xxyyy_yyyyzzz, g_0_x_xxyyy_yyyzzzz, g_0_x_xxyyy_yyzzzzz, g_0_x_xxyyy_yzzzzzz, g_0_x_xxyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyy_xxxxxxx[k] = -g_0_x_xxyy_xxxxxxx[k] * ab_y + g_0_x_xxyy_xxxxxxxy[k];

                g_0_x_xxyyy_xxxxxxy[k] = -g_0_x_xxyy_xxxxxxy[k] * ab_y + g_0_x_xxyy_xxxxxxyy[k];

                g_0_x_xxyyy_xxxxxxz[k] = -g_0_x_xxyy_xxxxxxz[k] * ab_y + g_0_x_xxyy_xxxxxxyz[k];

                g_0_x_xxyyy_xxxxxyy[k] = -g_0_x_xxyy_xxxxxyy[k] * ab_y + g_0_x_xxyy_xxxxxyyy[k];

                g_0_x_xxyyy_xxxxxyz[k] = -g_0_x_xxyy_xxxxxyz[k] * ab_y + g_0_x_xxyy_xxxxxyyz[k];

                g_0_x_xxyyy_xxxxxzz[k] = -g_0_x_xxyy_xxxxxzz[k] * ab_y + g_0_x_xxyy_xxxxxyzz[k];

                g_0_x_xxyyy_xxxxyyy[k] = -g_0_x_xxyy_xxxxyyy[k] * ab_y + g_0_x_xxyy_xxxxyyyy[k];

                g_0_x_xxyyy_xxxxyyz[k] = -g_0_x_xxyy_xxxxyyz[k] * ab_y + g_0_x_xxyy_xxxxyyyz[k];

                g_0_x_xxyyy_xxxxyzz[k] = -g_0_x_xxyy_xxxxyzz[k] * ab_y + g_0_x_xxyy_xxxxyyzz[k];

                g_0_x_xxyyy_xxxxzzz[k] = -g_0_x_xxyy_xxxxzzz[k] * ab_y + g_0_x_xxyy_xxxxyzzz[k];

                g_0_x_xxyyy_xxxyyyy[k] = -g_0_x_xxyy_xxxyyyy[k] * ab_y + g_0_x_xxyy_xxxyyyyy[k];

                g_0_x_xxyyy_xxxyyyz[k] = -g_0_x_xxyy_xxxyyyz[k] * ab_y + g_0_x_xxyy_xxxyyyyz[k];

                g_0_x_xxyyy_xxxyyzz[k] = -g_0_x_xxyy_xxxyyzz[k] * ab_y + g_0_x_xxyy_xxxyyyzz[k];

                g_0_x_xxyyy_xxxyzzz[k] = -g_0_x_xxyy_xxxyzzz[k] * ab_y + g_0_x_xxyy_xxxyyzzz[k];

                g_0_x_xxyyy_xxxzzzz[k] = -g_0_x_xxyy_xxxzzzz[k] * ab_y + g_0_x_xxyy_xxxyzzzz[k];

                g_0_x_xxyyy_xxyyyyy[k] = -g_0_x_xxyy_xxyyyyy[k] * ab_y + g_0_x_xxyy_xxyyyyyy[k];

                g_0_x_xxyyy_xxyyyyz[k] = -g_0_x_xxyy_xxyyyyz[k] * ab_y + g_0_x_xxyy_xxyyyyyz[k];

                g_0_x_xxyyy_xxyyyzz[k] = -g_0_x_xxyy_xxyyyzz[k] * ab_y + g_0_x_xxyy_xxyyyyzz[k];

                g_0_x_xxyyy_xxyyzzz[k] = -g_0_x_xxyy_xxyyzzz[k] * ab_y + g_0_x_xxyy_xxyyyzzz[k];

                g_0_x_xxyyy_xxyzzzz[k] = -g_0_x_xxyy_xxyzzzz[k] * ab_y + g_0_x_xxyy_xxyyzzzz[k];

                g_0_x_xxyyy_xxzzzzz[k] = -g_0_x_xxyy_xxzzzzz[k] * ab_y + g_0_x_xxyy_xxyzzzzz[k];

                g_0_x_xxyyy_xyyyyyy[k] = -g_0_x_xxyy_xyyyyyy[k] * ab_y + g_0_x_xxyy_xyyyyyyy[k];

                g_0_x_xxyyy_xyyyyyz[k] = -g_0_x_xxyy_xyyyyyz[k] * ab_y + g_0_x_xxyy_xyyyyyyz[k];

                g_0_x_xxyyy_xyyyyzz[k] = -g_0_x_xxyy_xyyyyzz[k] * ab_y + g_0_x_xxyy_xyyyyyzz[k];

                g_0_x_xxyyy_xyyyzzz[k] = -g_0_x_xxyy_xyyyzzz[k] * ab_y + g_0_x_xxyy_xyyyyzzz[k];

                g_0_x_xxyyy_xyyzzzz[k] = -g_0_x_xxyy_xyyzzzz[k] * ab_y + g_0_x_xxyy_xyyyzzzz[k];

                g_0_x_xxyyy_xyzzzzz[k] = -g_0_x_xxyy_xyzzzzz[k] * ab_y + g_0_x_xxyy_xyyzzzzz[k];

                g_0_x_xxyyy_xzzzzzz[k] = -g_0_x_xxyy_xzzzzzz[k] * ab_y + g_0_x_xxyy_xyzzzzzz[k];

                g_0_x_xxyyy_yyyyyyy[k] = -g_0_x_xxyy_yyyyyyy[k] * ab_y + g_0_x_xxyy_yyyyyyyy[k];

                g_0_x_xxyyy_yyyyyyz[k] = -g_0_x_xxyy_yyyyyyz[k] * ab_y + g_0_x_xxyy_yyyyyyyz[k];

                g_0_x_xxyyy_yyyyyzz[k] = -g_0_x_xxyy_yyyyyzz[k] * ab_y + g_0_x_xxyy_yyyyyyzz[k];

                g_0_x_xxyyy_yyyyzzz[k] = -g_0_x_xxyy_yyyyzzz[k] * ab_y + g_0_x_xxyy_yyyyyzzz[k];

                g_0_x_xxyyy_yyyzzzz[k] = -g_0_x_xxyy_yyyzzzz[k] * ab_y + g_0_x_xxyy_yyyyzzzz[k];

                g_0_x_xxyyy_yyzzzzz[k] = -g_0_x_xxyy_yyzzzzz[k] * ab_y + g_0_x_xxyy_yyyzzzzz[k];

                g_0_x_xxyyy_yzzzzzz[k] = -g_0_x_xxyy_yzzzzzz[k] * ab_y + g_0_x_xxyy_yyzzzzzz[k];

                g_0_x_xxyyy_zzzzzzz[k] = -g_0_x_xxyy_zzzzzzz[k] * ab_y + g_0_x_xxyy_yzzzzzzz[k];
            }

            /// Set up 252-288 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xxyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xxyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xxyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyz_xxxxxxx, g_0_x_xxyyz_xxxxxxy, g_0_x_xxyyz_xxxxxxz, g_0_x_xxyyz_xxxxxyy, g_0_x_xxyyz_xxxxxyz, g_0_x_xxyyz_xxxxxzz, g_0_x_xxyyz_xxxxyyy, g_0_x_xxyyz_xxxxyyz, g_0_x_xxyyz_xxxxyzz, g_0_x_xxyyz_xxxxzzz, g_0_x_xxyyz_xxxyyyy, g_0_x_xxyyz_xxxyyyz, g_0_x_xxyyz_xxxyyzz, g_0_x_xxyyz_xxxyzzz, g_0_x_xxyyz_xxxzzzz, g_0_x_xxyyz_xxyyyyy, g_0_x_xxyyz_xxyyyyz, g_0_x_xxyyz_xxyyyzz, g_0_x_xxyyz_xxyyzzz, g_0_x_xxyyz_xxyzzzz, g_0_x_xxyyz_xxzzzzz, g_0_x_xxyyz_xyyyyyy, g_0_x_xxyyz_xyyyyyz, g_0_x_xxyyz_xyyyyzz, g_0_x_xxyyz_xyyyzzz, g_0_x_xxyyz_xyyzzzz, g_0_x_xxyyz_xyzzzzz, g_0_x_xxyyz_xzzzzzz, g_0_x_xxyyz_yyyyyyy, g_0_x_xxyyz_yyyyyyz, g_0_x_xxyyz_yyyyyzz, g_0_x_xxyyz_yyyyzzz, g_0_x_xxyyz_yyyzzzz, g_0_x_xxyyz_yyzzzzz, g_0_x_xxyyz_yzzzzzz, g_0_x_xxyyz_zzzzzzz, g_0_x_xxyz_xxxxxxx, g_0_x_xxyz_xxxxxxxy, g_0_x_xxyz_xxxxxxy, g_0_x_xxyz_xxxxxxyy, g_0_x_xxyz_xxxxxxyz, g_0_x_xxyz_xxxxxxz, g_0_x_xxyz_xxxxxyy, g_0_x_xxyz_xxxxxyyy, g_0_x_xxyz_xxxxxyyz, g_0_x_xxyz_xxxxxyz, g_0_x_xxyz_xxxxxyzz, g_0_x_xxyz_xxxxxzz, g_0_x_xxyz_xxxxyyy, g_0_x_xxyz_xxxxyyyy, g_0_x_xxyz_xxxxyyyz, g_0_x_xxyz_xxxxyyz, g_0_x_xxyz_xxxxyyzz, g_0_x_xxyz_xxxxyzz, g_0_x_xxyz_xxxxyzzz, g_0_x_xxyz_xxxxzzz, g_0_x_xxyz_xxxyyyy, g_0_x_xxyz_xxxyyyyy, g_0_x_xxyz_xxxyyyyz, g_0_x_xxyz_xxxyyyz, g_0_x_xxyz_xxxyyyzz, g_0_x_xxyz_xxxyyzz, g_0_x_xxyz_xxxyyzzz, g_0_x_xxyz_xxxyzzz, g_0_x_xxyz_xxxyzzzz, g_0_x_xxyz_xxxzzzz, g_0_x_xxyz_xxyyyyy, g_0_x_xxyz_xxyyyyyy, g_0_x_xxyz_xxyyyyyz, g_0_x_xxyz_xxyyyyz, g_0_x_xxyz_xxyyyyzz, g_0_x_xxyz_xxyyyzz, g_0_x_xxyz_xxyyyzzz, g_0_x_xxyz_xxyyzzz, g_0_x_xxyz_xxyyzzzz, g_0_x_xxyz_xxyzzzz, g_0_x_xxyz_xxyzzzzz, g_0_x_xxyz_xxzzzzz, g_0_x_xxyz_xyyyyyy, g_0_x_xxyz_xyyyyyyy, g_0_x_xxyz_xyyyyyyz, g_0_x_xxyz_xyyyyyz, g_0_x_xxyz_xyyyyyzz, g_0_x_xxyz_xyyyyzz, g_0_x_xxyz_xyyyyzzz, g_0_x_xxyz_xyyyzzz, g_0_x_xxyz_xyyyzzzz, g_0_x_xxyz_xyyzzzz, g_0_x_xxyz_xyyzzzzz, g_0_x_xxyz_xyzzzzz, g_0_x_xxyz_xyzzzzzz, g_0_x_xxyz_xzzzzzz, g_0_x_xxyz_yyyyyyy, g_0_x_xxyz_yyyyyyyy, g_0_x_xxyz_yyyyyyyz, g_0_x_xxyz_yyyyyyz, g_0_x_xxyz_yyyyyyzz, g_0_x_xxyz_yyyyyzz, g_0_x_xxyz_yyyyyzzz, g_0_x_xxyz_yyyyzzz, g_0_x_xxyz_yyyyzzzz, g_0_x_xxyz_yyyzzzz, g_0_x_xxyz_yyyzzzzz, g_0_x_xxyz_yyzzzzz, g_0_x_xxyz_yyzzzzzz, g_0_x_xxyz_yzzzzzz, g_0_x_xxyz_yzzzzzzz, g_0_x_xxyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyz_xxxxxxx[k] = -g_0_x_xxyz_xxxxxxx[k] * ab_y + g_0_x_xxyz_xxxxxxxy[k];

                g_0_x_xxyyz_xxxxxxy[k] = -g_0_x_xxyz_xxxxxxy[k] * ab_y + g_0_x_xxyz_xxxxxxyy[k];

                g_0_x_xxyyz_xxxxxxz[k] = -g_0_x_xxyz_xxxxxxz[k] * ab_y + g_0_x_xxyz_xxxxxxyz[k];

                g_0_x_xxyyz_xxxxxyy[k] = -g_0_x_xxyz_xxxxxyy[k] * ab_y + g_0_x_xxyz_xxxxxyyy[k];

                g_0_x_xxyyz_xxxxxyz[k] = -g_0_x_xxyz_xxxxxyz[k] * ab_y + g_0_x_xxyz_xxxxxyyz[k];

                g_0_x_xxyyz_xxxxxzz[k] = -g_0_x_xxyz_xxxxxzz[k] * ab_y + g_0_x_xxyz_xxxxxyzz[k];

                g_0_x_xxyyz_xxxxyyy[k] = -g_0_x_xxyz_xxxxyyy[k] * ab_y + g_0_x_xxyz_xxxxyyyy[k];

                g_0_x_xxyyz_xxxxyyz[k] = -g_0_x_xxyz_xxxxyyz[k] * ab_y + g_0_x_xxyz_xxxxyyyz[k];

                g_0_x_xxyyz_xxxxyzz[k] = -g_0_x_xxyz_xxxxyzz[k] * ab_y + g_0_x_xxyz_xxxxyyzz[k];

                g_0_x_xxyyz_xxxxzzz[k] = -g_0_x_xxyz_xxxxzzz[k] * ab_y + g_0_x_xxyz_xxxxyzzz[k];

                g_0_x_xxyyz_xxxyyyy[k] = -g_0_x_xxyz_xxxyyyy[k] * ab_y + g_0_x_xxyz_xxxyyyyy[k];

                g_0_x_xxyyz_xxxyyyz[k] = -g_0_x_xxyz_xxxyyyz[k] * ab_y + g_0_x_xxyz_xxxyyyyz[k];

                g_0_x_xxyyz_xxxyyzz[k] = -g_0_x_xxyz_xxxyyzz[k] * ab_y + g_0_x_xxyz_xxxyyyzz[k];

                g_0_x_xxyyz_xxxyzzz[k] = -g_0_x_xxyz_xxxyzzz[k] * ab_y + g_0_x_xxyz_xxxyyzzz[k];

                g_0_x_xxyyz_xxxzzzz[k] = -g_0_x_xxyz_xxxzzzz[k] * ab_y + g_0_x_xxyz_xxxyzzzz[k];

                g_0_x_xxyyz_xxyyyyy[k] = -g_0_x_xxyz_xxyyyyy[k] * ab_y + g_0_x_xxyz_xxyyyyyy[k];

                g_0_x_xxyyz_xxyyyyz[k] = -g_0_x_xxyz_xxyyyyz[k] * ab_y + g_0_x_xxyz_xxyyyyyz[k];

                g_0_x_xxyyz_xxyyyzz[k] = -g_0_x_xxyz_xxyyyzz[k] * ab_y + g_0_x_xxyz_xxyyyyzz[k];

                g_0_x_xxyyz_xxyyzzz[k] = -g_0_x_xxyz_xxyyzzz[k] * ab_y + g_0_x_xxyz_xxyyyzzz[k];

                g_0_x_xxyyz_xxyzzzz[k] = -g_0_x_xxyz_xxyzzzz[k] * ab_y + g_0_x_xxyz_xxyyzzzz[k];

                g_0_x_xxyyz_xxzzzzz[k] = -g_0_x_xxyz_xxzzzzz[k] * ab_y + g_0_x_xxyz_xxyzzzzz[k];

                g_0_x_xxyyz_xyyyyyy[k] = -g_0_x_xxyz_xyyyyyy[k] * ab_y + g_0_x_xxyz_xyyyyyyy[k];

                g_0_x_xxyyz_xyyyyyz[k] = -g_0_x_xxyz_xyyyyyz[k] * ab_y + g_0_x_xxyz_xyyyyyyz[k];

                g_0_x_xxyyz_xyyyyzz[k] = -g_0_x_xxyz_xyyyyzz[k] * ab_y + g_0_x_xxyz_xyyyyyzz[k];

                g_0_x_xxyyz_xyyyzzz[k] = -g_0_x_xxyz_xyyyzzz[k] * ab_y + g_0_x_xxyz_xyyyyzzz[k];

                g_0_x_xxyyz_xyyzzzz[k] = -g_0_x_xxyz_xyyzzzz[k] * ab_y + g_0_x_xxyz_xyyyzzzz[k];

                g_0_x_xxyyz_xyzzzzz[k] = -g_0_x_xxyz_xyzzzzz[k] * ab_y + g_0_x_xxyz_xyyzzzzz[k];

                g_0_x_xxyyz_xzzzzzz[k] = -g_0_x_xxyz_xzzzzzz[k] * ab_y + g_0_x_xxyz_xyzzzzzz[k];

                g_0_x_xxyyz_yyyyyyy[k] = -g_0_x_xxyz_yyyyyyy[k] * ab_y + g_0_x_xxyz_yyyyyyyy[k];

                g_0_x_xxyyz_yyyyyyz[k] = -g_0_x_xxyz_yyyyyyz[k] * ab_y + g_0_x_xxyz_yyyyyyyz[k];

                g_0_x_xxyyz_yyyyyzz[k] = -g_0_x_xxyz_yyyyyzz[k] * ab_y + g_0_x_xxyz_yyyyyyzz[k];

                g_0_x_xxyyz_yyyyzzz[k] = -g_0_x_xxyz_yyyyzzz[k] * ab_y + g_0_x_xxyz_yyyyyzzz[k];

                g_0_x_xxyyz_yyyzzzz[k] = -g_0_x_xxyz_yyyzzzz[k] * ab_y + g_0_x_xxyz_yyyyzzzz[k];

                g_0_x_xxyyz_yyzzzzz[k] = -g_0_x_xxyz_yyzzzzz[k] * ab_y + g_0_x_xxyz_yyyzzzzz[k];

                g_0_x_xxyyz_yzzzzzz[k] = -g_0_x_xxyz_yzzzzzz[k] * ab_y + g_0_x_xxyz_yyzzzzzz[k];

                g_0_x_xxyyz_zzzzzzz[k] = -g_0_x_xxyz_zzzzzzz[k] * ab_y + g_0_x_xxyz_yzzzzzzz[k];
            }

            /// Set up 288-324 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_xxyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_xxyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_xxyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzz_xxxxxxx, g_0_x_xxyzz_xxxxxxy, g_0_x_xxyzz_xxxxxxz, g_0_x_xxyzz_xxxxxyy, g_0_x_xxyzz_xxxxxyz, g_0_x_xxyzz_xxxxxzz, g_0_x_xxyzz_xxxxyyy, g_0_x_xxyzz_xxxxyyz, g_0_x_xxyzz_xxxxyzz, g_0_x_xxyzz_xxxxzzz, g_0_x_xxyzz_xxxyyyy, g_0_x_xxyzz_xxxyyyz, g_0_x_xxyzz_xxxyyzz, g_0_x_xxyzz_xxxyzzz, g_0_x_xxyzz_xxxzzzz, g_0_x_xxyzz_xxyyyyy, g_0_x_xxyzz_xxyyyyz, g_0_x_xxyzz_xxyyyzz, g_0_x_xxyzz_xxyyzzz, g_0_x_xxyzz_xxyzzzz, g_0_x_xxyzz_xxzzzzz, g_0_x_xxyzz_xyyyyyy, g_0_x_xxyzz_xyyyyyz, g_0_x_xxyzz_xyyyyzz, g_0_x_xxyzz_xyyyzzz, g_0_x_xxyzz_xyyzzzz, g_0_x_xxyzz_xyzzzzz, g_0_x_xxyzz_xzzzzzz, g_0_x_xxyzz_yyyyyyy, g_0_x_xxyzz_yyyyyyz, g_0_x_xxyzz_yyyyyzz, g_0_x_xxyzz_yyyyzzz, g_0_x_xxyzz_yyyzzzz, g_0_x_xxyzz_yyzzzzz, g_0_x_xxyzz_yzzzzzz, g_0_x_xxyzz_zzzzzzz, g_0_x_xxzz_xxxxxxx, g_0_x_xxzz_xxxxxxxy, g_0_x_xxzz_xxxxxxy, g_0_x_xxzz_xxxxxxyy, g_0_x_xxzz_xxxxxxyz, g_0_x_xxzz_xxxxxxz, g_0_x_xxzz_xxxxxyy, g_0_x_xxzz_xxxxxyyy, g_0_x_xxzz_xxxxxyyz, g_0_x_xxzz_xxxxxyz, g_0_x_xxzz_xxxxxyzz, g_0_x_xxzz_xxxxxzz, g_0_x_xxzz_xxxxyyy, g_0_x_xxzz_xxxxyyyy, g_0_x_xxzz_xxxxyyyz, g_0_x_xxzz_xxxxyyz, g_0_x_xxzz_xxxxyyzz, g_0_x_xxzz_xxxxyzz, g_0_x_xxzz_xxxxyzzz, g_0_x_xxzz_xxxxzzz, g_0_x_xxzz_xxxyyyy, g_0_x_xxzz_xxxyyyyy, g_0_x_xxzz_xxxyyyyz, g_0_x_xxzz_xxxyyyz, g_0_x_xxzz_xxxyyyzz, g_0_x_xxzz_xxxyyzz, g_0_x_xxzz_xxxyyzzz, g_0_x_xxzz_xxxyzzz, g_0_x_xxzz_xxxyzzzz, g_0_x_xxzz_xxxzzzz, g_0_x_xxzz_xxyyyyy, g_0_x_xxzz_xxyyyyyy, g_0_x_xxzz_xxyyyyyz, g_0_x_xxzz_xxyyyyz, g_0_x_xxzz_xxyyyyzz, g_0_x_xxzz_xxyyyzz, g_0_x_xxzz_xxyyyzzz, g_0_x_xxzz_xxyyzzz, g_0_x_xxzz_xxyyzzzz, g_0_x_xxzz_xxyzzzz, g_0_x_xxzz_xxyzzzzz, g_0_x_xxzz_xxzzzzz, g_0_x_xxzz_xyyyyyy, g_0_x_xxzz_xyyyyyyy, g_0_x_xxzz_xyyyyyyz, g_0_x_xxzz_xyyyyyz, g_0_x_xxzz_xyyyyyzz, g_0_x_xxzz_xyyyyzz, g_0_x_xxzz_xyyyyzzz, g_0_x_xxzz_xyyyzzz, g_0_x_xxzz_xyyyzzzz, g_0_x_xxzz_xyyzzzz, g_0_x_xxzz_xyyzzzzz, g_0_x_xxzz_xyzzzzz, g_0_x_xxzz_xyzzzzzz, g_0_x_xxzz_xzzzzzz, g_0_x_xxzz_yyyyyyy, g_0_x_xxzz_yyyyyyyy, g_0_x_xxzz_yyyyyyyz, g_0_x_xxzz_yyyyyyz, g_0_x_xxzz_yyyyyyzz, g_0_x_xxzz_yyyyyzz, g_0_x_xxzz_yyyyyzzz, g_0_x_xxzz_yyyyzzz, g_0_x_xxzz_yyyyzzzz, g_0_x_xxzz_yyyzzzz, g_0_x_xxzz_yyyzzzzz, g_0_x_xxzz_yyzzzzz, g_0_x_xxzz_yyzzzzzz, g_0_x_xxzz_yzzzzzz, g_0_x_xxzz_yzzzzzzz, g_0_x_xxzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzz_xxxxxxx[k] = -g_0_x_xxzz_xxxxxxx[k] * ab_y + g_0_x_xxzz_xxxxxxxy[k];

                g_0_x_xxyzz_xxxxxxy[k] = -g_0_x_xxzz_xxxxxxy[k] * ab_y + g_0_x_xxzz_xxxxxxyy[k];

                g_0_x_xxyzz_xxxxxxz[k] = -g_0_x_xxzz_xxxxxxz[k] * ab_y + g_0_x_xxzz_xxxxxxyz[k];

                g_0_x_xxyzz_xxxxxyy[k] = -g_0_x_xxzz_xxxxxyy[k] * ab_y + g_0_x_xxzz_xxxxxyyy[k];

                g_0_x_xxyzz_xxxxxyz[k] = -g_0_x_xxzz_xxxxxyz[k] * ab_y + g_0_x_xxzz_xxxxxyyz[k];

                g_0_x_xxyzz_xxxxxzz[k] = -g_0_x_xxzz_xxxxxzz[k] * ab_y + g_0_x_xxzz_xxxxxyzz[k];

                g_0_x_xxyzz_xxxxyyy[k] = -g_0_x_xxzz_xxxxyyy[k] * ab_y + g_0_x_xxzz_xxxxyyyy[k];

                g_0_x_xxyzz_xxxxyyz[k] = -g_0_x_xxzz_xxxxyyz[k] * ab_y + g_0_x_xxzz_xxxxyyyz[k];

                g_0_x_xxyzz_xxxxyzz[k] = -g_0_x_xxzz_xxxxyzz[k] * ab_y + g_0_x_xxzz_xxxxyyzz[k];

                g_0_x_xxyzz_xxxxzzz[k] = -g_0_x_xxzz_xxxxzzz[k] * ab_y + g_0_x_xxzz_xxxxyzzz[k];

                g_0_x_xxyzz_xxxyyyy[k] = -g_0_x_xxzz_xxxyyyy[k] * ab_y + g_0_x_xxzz_xxxyyyyy[k];

                g_0_x_xxyzz_xxxyyyz[k] = -g_0_x_xxzz_xxxyyyz[k] * ab_y + g_0_x_xxzz_xxxyyyyz[k];

                g_0_x_xxyzz_xxxyyzz[k] = -g_0_x_xxzz_xxxyyzz[k] * ab_y + g_0_x_xxzz_xxxyyyzz[k];

                g_0_x_xxyzz_xxxyzzz[k] = -g_0_x_xxzz_xxxyzzz[k] * ab_y + g_0_x_xxzz_xxxyyzzz[k];

                g_0_x_xxyzz_xxxzzzz[k] = -g_0_x_xxzz_xxxzzzz[k] * ab_y + g_0_x_xxzz_xxxyzzzz[k];

                g_0_x_xxyzz_xxyyyyy[k] = -g_0_x_xxzz_xxyyyyy[k] * ab_y + g_0_x_xxzz_xxyyyyyy[k];

                g_0_x_xxyzz_xxyyyyz[k] = -g_0_x_xxzz_xxyyyyz[k] * ab_y + g_0_x_xxzz_xxyyyyyz[k];

                g_0_x_xxyzz_xxyyyzz[k] = -g_0_x_xxzz_xxyyyzz[k] * ab_y + g_0_x_xxzz_xxyyyyzz[k];

                g_0_x_xxyzz_xxyyzzz[k] = -g_0_x_xxzz_xxyyzzz[k] * ab_y + g_0_x_xxzz_xxyyyzzz[k];

                g_0_x_xxyzz_xxyzzzz[k] = -g_0_x_xxzz_xxyzzzz[k] * ab_y + g_0_x_xxzz_xxyyzzzz[k];

                g_0_x_xxyzz_xxzzzzz[k] = -g_0_x_xxzz_xxzzzzz[k] * ab_y + g_0_x_xxzz_xxyzzzzz[k];

                g_0_x_xxyzz_xyyyyyy[k] = -g_0_x_xxzz_xyyyyyy[k] * ab_y + g_0_x_xxzz_xyyyyyyy[k];

                g_0_x_xxyzz_xyyyyyz[k] = -g_0_x_xxzz_xyyyyyz[k] * ab_y + g_0_x_xxzz_xyyyyyyz[k];

                g_0_x_xxyzz_xyyyyzz[k] = -g_0_x_xxzz_xyyyyzz[k] * ab_y + g_0_x_xxzz_xyyyyyzz[k];

                g_0_x_xxyzz_xyyyzzz[k] = -g_0_x_xxzz_xyyyzzz[k] * ab_y + g_0_x_xxzz_xyyyyzzz[k];

                g_0_x_xxyzz_xyyzzzz[k] = -g_0_x_xxzz_xyyzzzz[k] * ab_y + g_0_x_xxzz_xyyyzzzz[k];

                g_0_x_xxyzz_xyzzzzz[k] = -g_0_x_xxzz_xyzzzzz[k] * ab_y + g_0_x_xxzz_xyyzzzzz[k];

                g_0_x_xxyzz_xzzzzzz[k] = -g_0_x_xxzz_xzzzzzz[k] * ab_y + g_0_x_xxzz_xyzzzzzz[k];

                g_0_x_xxyzz_yyyyyyy[k] = -g_0_x_xxzz_yyyyyyy[k] * ab_y + g_0_x_xxzz_yyyyyyyy[k];

                g_0_x_xxyzz_yyyyyyz[k] = -g_0_x_xxzz_yyyyyyz[k] * ab_y + g_0_x_xxzz_yyyyyyyz[k];

                g_0_x_xxyzz_yyyyyzz[k] = -g_0_x_xxzz_yyyyyzz[k] * ab_y + g_0_x_xxzz_yyyyyyzz[k];

                g_0_x_xxyzz_yyyyzzz[k] = -g_0_x_xxzz_yyyyzzz[k] * ab_y + g_0_x_xxzz_yyyyyzzz[k];

                g_0_x_xxyzz_yyyzzzz[k] = -g_0_x_xxzz_yyyzzzz[k] * ab_y + g_0_x_xxzz_yyyyzzzz[k];

                g_0_x_xxyzz_yyzzzzz[k] = -g_0_x_xxzz_yyzzzzz[k] * ab_y + g_0_x_xxzz_yyyzzzzz[k];

                g_0_x_xxyzz_yzzzzzz[k] = -g_0_x_xxzz_yzzzzzz[k] * ab_y + g_0_x_xxzz_yyzzzzzz[k];

                g_0_x_xxyzz_zzzzzzz[k] = -g_0_x_xxzz_zzzzzzz[k] * ab_y + g_0_x_xxzz_yzzzzzzz[k];
            }

            /// Set up 324-360 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_xxzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_xxzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_xxzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzz_xxxxxxx, g_0_x_xxzz_xxxxxxxz, g_0_x_xxzz_xxxxxxy, g_0_x_xxzz_xxxxxxyz, g_0_x_xxzz_xxxxxxz, g_0_x_xxzz_xxxxxxzz, g_0_x_xxzz_xxxxxyy, g_0_x_xxzz_xxxxxyyz, g_0_x_xxzz_xxxxxyz, g_0_x_xxzz_xxxxxyzz, g_0_x_xxzz_xxxxxzz, g_0_x_xxzz_xxxxxzzz, g_0_x_xxzz_xxxxyyy, g_0_x_xxzz_xxxxyyyz, g_0_x_xxzz_xxxxyyz, g_0_x_xxzz_xxxxyyzz, g_0_x_xxzz_xxxxyzz, g_0_x_xxzz_xxxxyzzz, g_0_x_xxzz_xxxxzzz, g_0_x_xxzz_xxxxzzzz, g_0_x_xxzz_xxxyyyy, g_0_x_xxzz_xxxyyyyz, g_0_x_xxzz_xxxyyyz, g_0_x_xxzz_xxxyyyzz, g_0_x_xxzz_xxxyyzz, g_0_x_xxzz_xxxyyzzz, g_0_x_xxzz_xxxyzzz, g_0_x_xxzz_xxxyzzzz, g_0_x_xxzz_xxxzzzz, g_0_x_xxzz_xxxzzzzz, g_0_x_xxzz_xxyyyyy, g_0_x_xxzz_xxyyyyyz, g_0_x_xxzz_xxyyyyz, g_0_x_xxzz_xxyyyyzz, g_0_x_xxzz_xxyyyzz, g_0_x_xxzz_xxyyyzzz, g_0_x_xxzz_xxyyzzz, g_0_x_xxzz_xxyyzzzz, g_0_x_xxzz_xxyzzzz, g_0_x_xxzz_xxyzzzzz, g_0_x_xxzz_xxzzzzz, g_0_x_xxzz_xxzzzzzz, g_0_x_xxzz_xyyyyyy, g_0_x_xxzz_xyyyyyyz, g_0_x_xxzz_xyyyyyz, g_0_x_xxzz_xyyyyyzz, g_0_x_xxzz_xyyyyzz, g_0_x_xxzz_xyyyyzzz, g_0_x_xxzz_xyyyzzz, g_0_x_xxzz_xyyyzzzz, g_0_x_xxzz_xyyzzzz, g_0_x_xxzz_xyyzzzzz, g_0_x_xxzz_xyzzzzz, g_0_x_xxzz_xyzzzzzz, g_0_x_xxzz_xzzzzzz, g_0_x_xxzz_xzzzzzzz, g_0_x_xxzz_yyyyyyy, g_0_x_xxzz_yyyyyyyz, g_0_x_xxzz_yyyyyyz, g_0_x_xxzz_yyyyyyzz, g_0_x_xxzz_yyyyyzz, g_0_x_xxzz_yyyyyzzz, g_0_x_xxzz_yyyyzzz, g_0_x_xxzz_yyyyzzzz, g_0_x_xxzz_yyyzzzz, g_0_x_xxzz_yyyzzzzz, g_0_x_xxzz_yyzzzzz, g_0_x_xxzz_yyzzzzzz, g_0_x_xxzz_yzzzzzz, g_0_x_xxzz_yzzzzzzz, g_0_x_xxzz_zzzzzzz, g_0_x_xxzz_zzzzzzzz, g_0_x_xxzzz_xxxxxxx, g_0_x_xxzzz_xxxxxxy, g_0_x_xxzzz_xxxxxxz, g_0_x_xxzzz_xxxxxyy, g_0_x_xxzzz_xxxxxyz, g_0_x_xxzzz_xxxxxzz, g_0_x_xxzzz_xxxxyyy, g_0_x_xxzzz_xxxxyyz, g_0_x_xxzzz_xxxxyzz, g_0_x_xxzzz_xxxxzzz, g_0_x_xxzzz_xxxyyyy, g_0_x_xxzzz_xxxyyyz, g_0_x_xxzzz_xxxyyzz, g_0_x_xxzzz_xxxyzzz, g_0_x_xxzzz_xxxzzzz, g_0_x_xxzzz_xxyyyyy, g_0_x_xxzzz_xxyyyyz, g_0_x_xxzzz_xxyyyzz, g_0_x_xxzzz_xxyyzzz, g_0_x_xxzzz_xxyzzzz, g_0_x_xxzzz_xxzzzzz, g_0_x_xxzzz_xyyyyyy, g_0_x_xxzzz_xyyyyyz, g_0_x_xxzzz_xyyyyzz, g_0_x_xxzzz_xyyyzzz, g_0_x_xxzzz_xyyzzzz, g_0_x_xxzzz_xyzzzzz, g_0_x_xxzzz_xzzzzzz, g_0_x_xxzzz_yyyyyyy, g_0_x_xxzzz_yyyyyyz, g_0_x_xxzzz_yyyyyzz, g_0_x_xxzzz_yyyyzzz, g_0_x_xxzzz_yyyzzzz, g_0_x_xxzzz_yyzzzzz, g_0_x_xxzzz_yzzzzzz, g_0_x_xxzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzz_xxxxxxx[k] = -g_0_x_xxzz_xxxxxxx[k] * ab_z + g_0_x_xxzz_xxxxxxxz[k];

                g_0_x_xxzzz_xxxxxxy[k] = -g_0_x_xxzz_xxxxxxy[k] * ab_z + g_0_x_xxzz_xxxxxxyz[k];

                g_0_x_xxzzz_xxxxxxz[k] = -g_0_x_xxzz_xxxxxxz[k] * ab_z + g_0_x_xxzz_xxxxxxzz[k];

                g_0_x_xxzzz_xxxxxyy[k] = -g_0_x_xxzz_xxxxxyy[k] * ab_z + g_0_x_xxzz_xxxxxyyz[k];

                g_0_x_xxzzz_xxxxxyz[k] = -g_0_x_xxzz_xxxxxyz[k] * ab_z + g_0_x_xxzz_xxxxxyzz[k];

                g_0_x_xxzzz_xxxxxzz[k] = -g_0_x_xxzz_xxxxxzz[k] * ab_z + g_0_x_xxzz_xxxxxzzz[k];

                g_0_x_xxzzz_xxxxyyy[k] = -g_0_x_xxzz_xxxxyyy[k] * ab_z + g_0_x_xxzz_xxxxyyyz[k];

                g_0_x_xxzzz_xxxxyyz[k] = -g_0_x_xxzz_xxxxyyz[k] * ab_z + g_0_x_xxzz_xxxxyyzz[k];

                g_0_x_xxzzz_xxxxyzz[k] = -g_0_x_xxzz_xxxxyzz[k] * ab_z + g_0_x_xxzz_xxxxyzzz[k];

                g_0_x_xxzzz_xxxxzzz[k] = -g_0_x_xxzz_xxxxzzz[k] * ab_z + g_0_x_xxzz_xxxxzzzz[k];

                g_0_x_xxzzz_xxxyyyy[k] = -g_0_x_xxzz_xxxyyyy[k] * ab_z + g_0_x_xxzz_xxxyyyyz[k];

                g_0_x_xxzzz_xxxyyyz[k] = -g_0_x_xxzz_xxxyyyz[k] * ab_z + g_0_x_xxzz_xxxyyyzz[k];

                g_0_x_xxzzz_xxxyyzz[k] = -g_0_x_xxzz_xxxyyzz[k] * ab_z + g_0_x_xxzz_xxxyyzzz[k];

                g_0_x_xxzzz_xxxyzzz[k] = -g_0_x_xxzz_xxxyzzz[k] * ab_z + g_0_x_xxzz_xxxyzzzz[k];

                g_0_x_xxzzz_xxxzzzz[k] = -g_0_x_xxzz_xxxzzzz[k] * ab_z + g_0_x_xxzz_xxxzzzzz[k];

                g_0_x_xxzzz_xxyyyyy[k] = -g_0_x_xxzz_xxyyyyy[k] * ab_z + g_0_x_xxzz_xxyyyyyz[k];

                g_0_x_xxzzz_xxyyyyz[k] = -g_0_x_xxzz_xxyyyyz[k] * ab_z + g_0_x_xxzz_xxyyyyzz[k];

                g_0_x_xxzzz_xxyyyzz[k] = -g_0_x_xxzz_xxyyyzz[k] * ab_z + g_0_x_xxzz_xxyyyzzz[k];

                g_0_x_xxzzz_xxyyzzz[k] = -g_0_x_xxzz_xxyyzzz[k] * ab_z + g_0_x_xxzz_xxyyzzzz[k];

                g_0_x_xxzzz_xxyzzzz[k] = -g_0_x_xxzz_xxyzzzz[k] * ab_z + g_0_x_xxzz_xxyzzzzz[k];

                g_0_x_xxzzz_xxzzzzz[k] = -g_0_x_xxzz_xxzzzzz[k] * ab_z + g_0_x_xxzz_xxzzzzzz[k];

                g_0_x_xxzzz_xyyyyyy[k] = -g_0_x_xxzz_xyyyyyy[k] * ab_z + g_0_x_xxzz_xyyyyyyz[k];

                g_0_x_xxzzz_xyyyyyz[k] = -g_0_x_xxzz_xyyyyyz[k] * ab_z + g_0_x_xxzz_xyyyyyzz[k];

                g_0_x_xxzzz_xyyyyzz[k] = -g_0_x_xxzz_xyyyyzz[k] * ab_z + g_0_x_xxzz_xyyyyzzz[k];

                g_0_x_xxzzz_xyyyzzz[k] = -g_0_x_xxzz_xyyyzzz[k] * ab_z + g_0_x_xxzz_xyyyzzzz[k];

                g_0_x_xxzzz_xyyzzzz[k] = -g_0_x_xxzz_xyyzzzz[k] * ab_z + g_0_x_xxzz_xyyzzzzz[k];

                g_0_x_xxzzz_xyzzzzz[k] = -g_0_x_xxzz_xyzzzzz[k] * ab_z + g_0_x_xxzz_xyzzzzzz[k];

                g_0_x_xxzzz_xzzzzzz[k] = -g_0_x_xxzz_xzzzzzz[k] * ab_z + g_0_x_xxzz_xzzzzzzz[k];

                g_0_x_xxzzz_yyyyyyy[k] = -g_0_x_xxzz_yyyyyyy[k] * ab_z + g_0_x_xxzz_yyyyyyyz[k];

                g_0_x_xxzzz_yyyyyyz[k] = -g_0_x_xxzz_yyyyyyz[k] * ab_z + g_0_x_xxzz_yyyyyyzz[k];

                g_0_x_xxzzz_yyyyyzz[k] = -g_0_x_xxzz_yyyyyzz[k] * ab_z + g_0_x_xxzz_yyyyyzzz[k];

                g_0_x_xxzzz_yyyyzzz[k] = -g_0_x_xxzz_yyyyzzz[k] * ab_z + g_0_x_xxzz_yyyyzzzz[k];

                g_0_x_xxzzz_yyyzzzz[k] = -g_0_x_xxzz_yyyzzzz[k] * ab_z + g_0_x_xxzz_yyyzzzzz[k];

                g_0_x_xxzzz_yyzzzzz[k] = -g_0_x_xxzz_yyzzzzz[k] * ab_z + g_0_x_xxzz_yyzzzzzz[k];

                g_0_x_xxzzz_yzzzzzz[k] = -g_0_x_xxzz_yzzzzzz[k] * ab_z + g_0_x_xxzz_yzzzzzzz[k];

                g_0_x_xxzzz_zzzzzzz[k] = -g_0_x_xxzz_zzzzzzz[k] * ab_z + g_0_x_xxzz_zzzzzzzz[k];
            }

            /// Set up 360-396 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_xyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_xyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_xyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyy_xxxxxxx, g_0_x_xyyy_xxxxxxxy, g_0_x_xyyy_xxxxxxy, g_0_x_xyyy_xxxxxxyy, g_0_x_xyyy_xxxxxxyz, g_0_x_xyyy_xxxxxxz, g_0_x_xyyy_xxxxxyy, g_0_x_xyyy_xxxxxyyy, g_0_x_xyyy_xxxxxyyz, g_0_x_xyyy_xxxxxyz, g_0_x_xyyy_xxxxxyzz, g_0_x_xyyy_xxxxxzz, g_0_x_xyyy_xxxxyyy, g_0_x_xyyy_xxxxyyyy, g_0_x_xyyy_xxxxyyyz, g_0_x_xyyy_xxxxyyz, g_0_x_xyyy_xxxxyyzz, g_0_x_xyyy_xxxxyzz, g_0_x_xyyy_xxxxyzzz, g_0_x_xyyy_xxxxzzz, g_0_x_xyyy_xxxyyyy, g_0_x_xyyy_xxxyyyyy, g_0_x_xyyy_xxxyyyyz, g_0_x_xyyy_xxxyyyz, g_0_x_xyyy_xxxyyyzz, g_0_x_xyyy_xxxyyzz, g_0_x_xyyy_xxxyyzzz, g_0_x_xyyy_xxxyzzz, g_0_x_xyyy_xxxyzzzz, g_0_x_xyyy_xxxzzzz, g_0_x_xyyy_xxyyyyy, g_0_x_xyyy_xxyyyyyy, g_0_x_xyyy_xxyyyyyz, g_0_x_xyyy_xxyyyyz, g_0_x_xyyy_xxyyyyzz, g_0_x_xyyy_xxyyyzz, g_0_x_xyyy_xxyyyzzz, g_0_x_xyyy_xxyyzzz, g_0_x_xyyy_xxyyzzzz, g_0_x_xyyy_xxyzzzz, g_0_x_xyyy_xxyzzzzz, g_0_x_xyyy_xxzzzzz, g_0_x_xyyy_xyyyyyy, g_0_x_xyyy_xyyyyyyy, g_0_x_xyyy_xyyyyyyz, g_0_x_xyyy_xyyyyyz, g_0_x_xyyy_xyyyyyzz, g_0_x_xyyy_xyyyyzz, g_0_x_xyyy_xyyyyzzz, g_0_x_xyyy_xyyyzzz, g_0_x_xyyy_xyyyzzzz, g_0_x_xyyy_xyyzzzz, g_0_x_xyyy_xyyzzzzz, g_0_x_xyyy_xyzzzzz, g_0_x_xyyy_xyzzzzzz, g_0_x_xyyy_xzzzzzz, g_0_x_xyyy_yyyyyyy, g_0_x_xyyy_yyyyyyyy, g_0_x_xyyy_yyyyyyyz, g_0_x_xyyy_yyyyyyz, g_0_x_xyyy_yyyyyyzz, g_0_x_xyyy_yyyyyzz, g_0_x_xyyy_yyyyyzzz, g_0_x_xyyy_yyyyzzz, g_0_x_xyyy_yyyyzzzz, g_0_x_xyyy_yyyzzzz, g_0_x_xyyy_yyyzzzzz, g_0_x_xyyy_yyzzzzz, g_0_x_xyyy_yyzzzzzz, g_0_x_xyyy_yzzzzzz, g_0_x_xyyy_yzzzzzzz, g_0_x_xyyy_zzzzzzz, g_0_x_xyyyy_xxxxxxx, g_0_x_xyyyy_xxxxxxy, g_0_x_xyyyy_xxxxxxz, g_0_x_xyyyy_xxxxxyy, g_0_x_xyyyy_xxxxxyz, g_0_x_xyyyy_xxxxxzz, g_0_x_xyyyy_xxxxyyy, g_0_x_xyyyy_xxxxyyz, g_0_x_xyyyy_xxxxyzz, g_0_x_xyyyy_xxxxzzz, g_0_x_xyyyy_xxxyyyy, g_0_x_xyyyy_xxxyyyz, g_0_x_xyyyy_xxxyyzz, g_0_x_xyyyy_xxxyzzz, g_0_x_xyyyy_xxxzzzz, g_0_x_xyyyy_xxyyyyy, g_0_x_xyyyy_xxyyyyz, g_0_x_xyyyy_xxyyyzz, g_0_x_xyyyy_xxyyzzz, g_0_x_xyyyy_xxyzzzz, g_0_x_xyyyy_xxzzzzz, g_0_x_xyyyy_xyyyyyy, g_0_x_xyyyy_xyyyyyz, g_0_x_xyyyy_xyyyyzz, g_0_x_xyyyy_xyyyzzz, g_0_x_xyyyy_xyyzzzz, g_0_x_xyyyy_xyzzzzz, g_0_x_xyyyy_xzzzzzz, g_0_x_xyyyy_yyyyyyy, g_0_x_xyyyy_yyyyyyz, g_0_x_xyyyy_yyyyyzz, g_0_x_xyyyy_yyyyzzz, g_0_x_xyyyy_yyyzzzz, g_0_x_xyyyy_yyzzzzz, g_0_x_xyyyy_yzzzzzz, g_0_x_xyyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyy_xxxxxxx[k] = -g_0_x_xyyy_xxxxxxx[k] * ab_y + g_0_x_xyyy_xxxxxxxy[k];

                g_0_x_xyyyy_xxxxxxy[k] = -g_0_x_xyyy_xxxxxxy[k] * ab_y + g_0_x_xyyy_xxxxxxyy[k];

                g_0_x_xyyyy_xxxxxxz[k] = -g_0_x_xyyy_xxxxxxz[k] * ab_y + g_0_x_xyyy_xxxxxxyz[k];

                g_0_x_xyyyy_xxxxxyy[k] = -g_0_x_xyyy_xxxxxyy[k] * ab_y + g_0_x_xyyy_xxxxxyyy[k];

                g_0_x_xyyyy_xxxxxyz[k] = -g_0_x_xyyy_xxxxxyz[k] * ab_y + g_0_x_xyyy_xxxxxyyz[k];

                g_0_x_xyyyy_xxxxxzz[k] = -g_0_x_xyyy_xxxxxzz[k] * ab_y + g_0_x_xyyy_xxxxxyzz[k];

                g_0_x_xyyyy_xxxxyyy[k] = -g_0_x_xyyy_xxxxyyy[k] * ab_y + g_0_x_xyyy_xxxxyyyy[k];

                g_0_x_xyyyy_xxxxyyz[k] = -g_0_x_xyyy_xxxxyyz[k] * ab_y + g_0_x_xyyy_xxxxyyyz[k];

                g_0_x_xyyyy_xxxxyzz[k] = -g_0_x_xyyy_xxxxyzz[k] * ab_y + g_0_x_xyyy_xxxxyyzz[k];

                g_0_x_xyyyy_xxxxzzz[k] = -g_0_x_xyyy_xxxxzzz[k] * ab_y + g_0_x_xyyy_xxxxyzzz[k];

                g_0_x_xyyyy_xxxyyyy[k] = -g_0_x_xyyy_xxxyyyy[k] * ab_y + g_0_x_xyyy_xxxyyyyy[k];

                g_0_x_xyyyy_xxxyyyz[k] = -g_0_x_xyyy_xxxyyyz[k] * ab_y + g_0_x_xyyy_xxxyyyyz[k];

                g_0_x_xyyyy_xxxyyzz[k] = -g_0_x_xyyy_xxxyyzz[k] * ab_y + g_0_x_xyyy_xxxyyyzz[k];

                g_0_x_xyyyy_xxxyzzz[k] = -g_0_x_xyyy_xxxyzzz[k] * ab_y + g_0_x_xyyy_xxxyyzzz[k];

                g_0_x_xyyyy_xxxzzzz[k] = -g_0_x_xyyy_xxxzzzz[k] * ab_y + g_0_x_xyyy_xxxyzzzz[k];

                g_0_x_xyyyy_xxyyyyy[k] = -g_0_x_xyyy_xxyyyyy[k] * ab_y + g_0_x_xyyy_xxyyyyyy[k];

                g_0_x_xyyyy_xxyyyyz[k] = -g_0_x_xyyy_xxyyyyz[k] * ab_y + g_0_x_xyyy_xxyyyyyz[k];

                g_0_x_xyyyy_xxyyyzz[k] = -g_0_x_xyyy_xxyyyzz[k] * ab_y + g_0_x_xyyy_xxyyyyzz[k];

                g_0_x_xyyyy_xxyyzzz[k] = -g_0_x_xyyy_xxyyzzz[k] * ab_y + g_0_x_xyyy_xxyyyzzz[k];

                g_0_x_xyyyy_xxyzzzz[k] = -g_0_x_xyyy_xxyzzzz[k] * ab_y + g_0_x_xyyy_xxyyzzzz[k];

                g_0_x_xyyyy_xxzzzzz[k] = -g_0_x_xyyy_xxzzzzz[k] * ab_y + g_0_x_xyyy_xxyzzzzz[k];

                g_0_x_xyyyy_xyyyyyy[k] = -g_0_x_xyyy_xyyyyyy[k] * ab_y + g_0_x_xyyy_xyyyyyyy[k];

                g_0_x_xyyyy_xyyyyyz[k] = -g_0_x_xyyy_xyyyyyz[k] * ab_y + g_0_x_xyyy_xyyyyyyz[k];

                g_0_x_xyyyy_xyyyyzz[k] = -g_0_x_xyyy_xyyyyzz[k] * ab_y + g_0_x_xyyy_xyyyyyzz[k];

                g_0_x_xyyyy_xyyyzzz[k] = -g_0_x_xyyy_xyyyzzz[k] * ab_y + g_0_x_xyyy_xyyyyzzz[k];

                g_0_x_xyyyy_xyyzzzz[k] = -g_0_x_xyyy_xyyzzzz[k] * ab_y + g_0_x_xyyy_xyyyzzzz[k];

                g_0_x_xyyyy_xyzzzzz[k] = -g_0_x_xyyy_xyzzzzz[k] * ab_y + g_0_x_xyyy_xyyzzzzz[k];

                g_0_x_xyyyy_xzzzzzz[k] = -g_0_x_xyyy_xzzzzzz[k] * ab_y + g_0_x_xyyy_xyzzzzzz[k];

                g_0_x_xyyyy_yyyyyyy[k] = -g_0_x_xyyy_yyyyyyy[k] * ab_y + g_0_x_xyyy_yyyyyyyy[k];

                g_0_x_xyyyy_yyyyyyz[k] = -g_0_x_xyyy_yyyyyyz[k] * ab_y + g_0_x_xyyy_yyyyyyyz[k];

                g_0_x_xyyyy_yyyyyzz[k] = -g_0_x_xyyy_yyyyyzz[k] * ab_y + g_0_x_xyyy_yyyyyyzz[k];

                g_0_x_xyyyy_yyyyzzz[k] = -g_0_x_xyyy_yyyyzzz[k] * ab_y + g_0_x_xyyy_yyyyyzzz[k];

                g_0_x_xyyyy_yyyzzzz[k] = -g_0_x_xyyy_yyyzzzz[k] * ab_y + g_0_x_xyyy_yyyyzzzz[k];

                g_0_x_xyyyy_yyzzzzz[k] = -g_0_x_xyyy_yyzzzzz[k] * ab_y + g_0_x_xyyy_yyyzzzzz[k];

                g_0_x_xyyyy_yzzzzzz[k] = -g_0_x_xyyy_yzzzzzz[k] * ab_y + g_0_x_xyyy_yyzzzzzz[k];

                g_0_x_xyyyy_zzzzzzz[k] = -g_0_x_xyyy_zzzzzzz[k] * ab_y + g_0_x_xyyy_yzzzzzzz[k];
            }

            /// Set up 396-432 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_xyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_xyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_xyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyz_xxxxxxx, g_0_x_xyyyz_xxxxxxy, g_0_x_xyyyz_xxxxxxz, g_0_x_xyyyz_xxxxxyy, g_0_x_xyyyz_xxxxxyz, g_0_x_xyyyz_xxxxxzz, g_0_x_xyyyz_xxxxyyy, g_0_x_xyyyz_xxxxyyz, g_0_x_xyyyz_xxxxyzz, g_0_x_xyyyz_xxxxzzz, g_0_x_xyyyz_xxxyyyy, g_0_x_xyyyz_xxxyyyz, g_0_x_xyyyz_xxxyyzz, g_0_x_xyyyz_xxxyzzz, g_0_x_xyyyz_xxxzzzz, g_0_x_xyyyz_xxyyyyy, g_0_x_xyyyz_xxyyyyz, g_0_x_xyyyz_xxyyyzz, g_0_x_xyyyz_xxyyzzz, g_0_x_xyyyz_xxyzzzz, g_0_x_xyyyz_xxzzzzz, g_0_x_xyyyz_xyyyyyy, g_0_x_xyyyz_xyyyyyz, g_0_x_xyyyz_xyyyyzz, g_0_x_xyyyz_xyyyzzz, g_0_x_xyyyz_xyyzzzz, g_0_x_xyyyz_xyzzzzz, g_0_x_xyyyz_xzzzzzz, g_0_x_xyyyz_yyyyyyy, g_0_x_xyyyz_yyyyyyz, g_0_x_xyyyz_yyyyyzz, g_0_x_xyyyz_yyyyzzz, g_0_x_xyyyz_yyyzzzz, g_0_x_xyyyz_yyzzzzz, g_0_x_xyyyz_yzzzzzz, g_0_x_xyyyz_zzzzzzz, g_0_x_xyyz_xxxxxxx, g_0_x_xyyz_xxxxxxxy, g_0_x_xyyz_xxxxxxy, g_0_x_xyyz_xxxxxxyy, g_0_x_xyyz_xxxxxxyz, g_0_x_xyyz_xxxxxxz, g_0_x_xyyz_xxxxxyy, g_0_x_xyyz_xxxxxyyy, g_0_x_xyyz_xxxxxyyz, g_0_x_xyyz_xxxxxyz, g_0_x_xyyz_xxxxxyzz, g_0_x_xyyz_xxxxxzz, g_0_x_xyyz_xxxxyyy, g_0_x_xyyz_xxxxyyyy, g_0_x_xyyz_xxxxyyyz, g_0_x_xyyz_xxxxyyz, g_0_x_xyyz_xxxxyyzz, g_0_x_xyyz_xxxxyzz, g_0_x_xyyz_xxxxyzzz, g_0_x_xyyz_xxxxzzz, g_0_x_xyyz_xxxyyyy, g_0_x_xyyz_xxxyyyyy, g_0_x_xyyz_xxxyyyyz, g_0_x_xyyz_xxxyyyz, g_0_x_xyyz_xxxyyyzz, g_0_x_xyyz_xxxyyzz, g_0_x_xyyz_xxxyyzzz, g_0_x_xyyz_xxxyzzz, g_0_x_xyyz_xxxyzzzz, g_0_x_xyyz_xxxzzzz, g_0_x_xyyz_xxyyyyy, g_0_x_xyyz_xxyyyyyy, g_0_x_xyyz_xxyyyyyz, g_0_x_xyyz_xxyyyyz, g_0_x_xyyz_xxyyyyzz, g_0_x_xyyz_xxyyyzz, g_0_x_xyyz_xxyyyzzz, g_0_x_xyyz_xxyyzzz, g_0_x_xyyz_xxyyzzzz, g_0_x_xyyz_xxyzzzz, g_0_x_xyyz_xxyzzzzz, g_0_x_xyyz_xxzzzzz, g_0_x_xyyz_xyyyyyy, g_0_x_xyyz_xyyyyyyy, g_0_x_xyyz_xyyyyyyz, g_0_x_xyyz_xyyyyyz, g_0_x_xyyz_xyyyyyzz, g_0_x_xyyz_xyyyyzz, g_0_x_xyyz_xyyyyzzz, g_0_x_xyyz_xyyyzzz, g_0_x_xyyz_xyyyzzzz, g_0_x_xyyz_xyyzzzz, g_0_x_xyyz_xyyzzzzz, g_0_x_xyyz_xyzzzzz, g_0_x_xyyz_xyzzzzzz, g_0_x_xyyz_xzzzzzz, g_0_x_xyyz_yyyyyyy, g_0_x_xyyz_yyyyyyyy, g_0_x_xyyz_yyyyyyyz, g_0_x_xyyz_yyyyyyz, g_0_x_xyyz_yyyyyyzz, g_0_x_xyyz_yyyyyzz, g_0_x_xyyz_yyyyyzzz, g_0_x_xyyz_yyyyzzz, g_0_x_xyyz_yyyyzzzz, g_0_x_xyyz_yyyzzzz, g_0_x_xyyz_yyyzzzzz, g_0_x_xyyz_yyzzzzz, g_0_x_xyyz_yyzzzzzz, g_0_x_xyyz_yzzzzzz, g_0_x_xyyz_yzzzzzzz, g_0_x_xyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyz_xxxxxxx[k] = -g_0_x_xyyz_xxxxxxx[k] * ab_y + g_0_x_xyyz_xxxxxxxy[k];

                g_0_x_xyyyz_xxxxxxy[k] = -g_0_x_xyyz_xxxxxxy[k] * ab_y + g_0_x_xyyz_xxxxxxyy[k];

                g_0_x_xyyyz_xxxxxxz[k] = -g_0_x_xyyz_xxxxxxz[k] * ab_y + g_0_x_xyyz_xxxxxxyz[k];

                g_0_x_xyyyz_xxxxxyy[k] = -g_0_x_xyyz_xxxxxyy[k] * ab_y + g_0_x_xyyz_xxxxxyyy[k];

                g_0_x_xyyyz_xxxxxyz[k] = -g_0_x_xyyz_xxxxxyz[k] * ab_y + g_0_x_xyyz_xxxxxyyz[k];

                g_0_x_xyyyz_xxxxxzz[k] = -g_0_x_xyyz_xxxxxzz[k] * ab_y + g_0_x_xyyz_xxxxxyzz[k];

                g_0_x_xyyyz_xxxxyyy[k] = -g_0_x_xyyz_xxxxyyy[k] * ab_y + g_0_x_xyyz_xxxxyyyy[k];

                g_0_x_xyyyz_xxxxyyz[k] = -g_0_x_xyyz_xxxxyyz[k] * ab_y + g_0_x_xyyz_xxxxyyyz[k];

                g_0_x_xyyyz_xxxxyzz[k] = -g_0_x_xyyz_xxxxyzz[k] * ab_y + g_0_x_xyyz_xxxxyyzz[k];

                g_0_x_xyyyz_xxxxzzz[k] = -g_0_x_xyyz_xxxxzzz[k] * ab_y + g_0_x_xyyz_xxxxyzzz[k];

                g_0_x_xyyyz_xxxyyyy[k] = -g_0_x_xyyz_xxxyyyy[k] * ab_y + g_0_x_xyyz_xxxyyyyy[k];

                g_0_x_xyyyz_xxxyyyz[k] = -g_0_x_xyyz_xxxyyyz[k] * ab_y + g_0_x_xyyz_xxxyyyyz[k];

                g_0_x_xyyyz_xxxyyzz[k] = -g_0_x_xyyz_xxxyyzz[k] * ab_y + g_0_x_xyyz_xxxyyyzz[k];

                g_0_x_xyyyz_xxxyzzz[k] = -g_0_x_xyyz_xxxyzzz[k] * ab_y + g_0_x_xyyz_xxxyyzzz[k];

                g_0_x_xyyyz_xxxzzzz[k] = -g_0_x_xyyz_xxxzzzz[k] * ab_y + g_0_x_xyyz_xxxyzzzz[k];

                g_0_x_xyyyz_xxyyyyy[k] = -g_0_x_xyyz_xxyyyyy[k] * ab_y + g_0_x_xyyz_xxyyyyyy[k];

                g_0_x_xyyyz_xxyyyyz[k] = -g_0_x_xyyz_xxyyyyz[k] * ab_y + g_0_x_xyyz_xxyyyyyz[k];

                g_0_x_xyyyz_xxyyyzz[k] = -g_0_x_xyyz_xxyyyzz[k] * ab_y + g_0_x_xyyz_xxyyyyzz[k];

                g_0_x_xyyyz_xxyyzzz[k] = -g_0_x_xyyz_xxyyzzz[k] * ab_y + g_0_x_xyyz_xxyyyzzz[k];

                g_0_x_xyyyz_xxyzzzz[k] = -g_0_x_xyyz_xxyzzzz[k] * ab_y + g_0_x_xyyz_xxyyzzzz[k];

                g_0_x_xyyyz_xxzzzzz[k] = -g_0_x_xyyz_xxzzzzz[k] * ab_y + g_0_x_xyyz_xxyzzzzz[k];

                g_0_x_xyyyz_xyyyyyy[k] = -g_0_x_xyyz_xyyyyyy[k] * ab_y + g_0_x_xyyz_xyyyyyyy[k];

                g_0_x_xyyyz_xyyyyyz[k] = -g_0_x_xyyz_xyyyyyz[k] * ab_y + g_0_x_xyyz_xyyyyyyz[k];

                g_0_x_xyyyz_xyyyyzz[k] = -g_0_x_xyyz_xyyyyzz[k] * ab_y + g_0_x_xyyz_xyyyyyzz[k];

                g_0_x_xyyyz_xyyyzzz[k] = -g_0_x_xyyz_xyyyzzz[k] * ab_y + g_0_x_xyyz_xyyyyzzz[k];

                g_0_x_xyyyz_xyyzzzz[k] = -g_0_x_xyyz_xyyzzzz[k] * ab_y + g_0_x_xyyz_xyyyzzzz[k];

                g_0_x_xyyyz_xyzzzzz[k] = -g_0_x_xyyz_xyzzzzz[k] * ab_y + g_0_x_xyyz_xyyzzzzz[k];

                g_0_x_xyyyz_xzzzzzz[k] = -g_0_x_xyyz_xzzzzzz[k] * ab_y + g_0_x_xyyz_xyzzzzzz[k];

                g_0_x_xyyyz_yyyyyyy[k] = -g_0_x_xyyz_yyyyyyy[k] * ab_y + g_0_x_xyyz_yyyyyyyy[k];

                g_0_x_xyyyz_yyyyyyz[k] = -g_0_x_xyyz_yyyyyyz[k] * ab_y + g_0_x_xyyz_yyyyyyyz[k];

                g_0_x_xyyyz_yyyyyzz[k] = -g_0_x_xyyz_yyyyyzz[k] * ab_y + g_0_x_xyyz_yyyyyyzz[k];

                g_0_x_xyyyz_yyyyzzz[k] = -g_0_x_xyyz_yyyyzzz[k] * ab_y + g_0_x_xyyz_yyyyyzzz[k];

                g_0_x_xyyyz_yyyzzzz[k] = -g_0_x_xyyz_yyyzzzz[k] * ab_y + g_0_x_xyyz_yyyyzzzz[k];

                g_0_x_xyyyz_yyzzzzz[k] = -g_0_x_xyyz_yyzzzzz[k] * ab_y + g_0_x_xyyz_yyyzzzzz[k];

                g_0_x_xyyyz_yzzzzzz[k] = -g_0_x_xyyz_yzzzzzz[k] * ab_y + g_0_x_xyyz_yyzzzzzz[k];

                g_0_x_xyyyz_zzzzzzz[k] = -g_0_x_xyyz_zzzzzzz[k] * ab_y + g_0_x_xyyz_yzzzzzzz[k];
            }

            /// Set up 432-468 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_x_xyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_x_xyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_x_xyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzz_xxxxxxx, g_0_x_xyyzz_xxxxxxy, g_0_x_xyyzz_xxxxxxz, g_0_x_xyyzz_xxxxxyy, g_0_x_xyyzz_xxxxxyz, g_0_x_xyyzz_xxxxxzz, g_0_x_xyyzz_xxxxyyy, g_0_x_xyyzz_xxxxyyz, g_0_x_xyyzz_xxxxyzz, g_0_x_xyyzz_xxxxzzz, g_0_x_xyyzz_xxxyyyy, g_0_x_xyyzz_xxxyyyz, g_0_x_xyyzz_xxxyyzz, g_0_x_xyyzz_xxxyzzz, g_0_x_xyyzz_xxxzzzz, g_0_x_xyyzz_xxyyyyy, g_0_x_xyyzz_xxyyyyz, g_0_x_xyyzz_xxyyyzz, g_0_x_xyyzz_xxyyzzz, g_0_x_xyyzz_xxyzzzz, g_0_x_xyyzz_xxzzzzz, g_0_x_xyyzz_xyyyyyy, g_0_x_xyyzz_xyyyyyz, g_0_x_xyyzz_xyyyyzz, g_0_x_xyyzz_xyyyzzz, g_0_x_xyyzz_xyyzzzz, g_0_x_xyyzz_xyzzzzz, g_0_x_xyyzz_xzzzzzz, g_0_x_xyyzz_yyyyyyy, g_0_x_xyyzz_yyyyyyz, g_0_x_xyyzz_yyyyyzz, g_0_x_xyyzz_yyyyzzz, g_0_x_xyyzz_yyyzzzz, g_0_x_xyyzz_yyzzzzz, g_0_x_xyyzz_yzzzzzz, g_0_x_xyyzz_zzzzzzz, g_0_x_xyzz_xxxxxxx, g_0_x_xyzz_xxxxxxxy, g_0_x_xyzz_xxxxxxy, g_0_x_xyzz_xxxxxxyy, g_0_x_xyzz_xxxxxxyz, g_0_x_xyzz_xxxxxxz, g_0_x_xyzz_xxxxxyy, g_0_x_xyzz_xxxxxyyy, g_0_x_xyzz_xxxxxyyz, g_0_x_xyzz_xxxxxyz, g_0_x_xyzz_xxxxxyzz, g_0_x_xyzz_xxxxxzz, g_0_x_xyzz_xxxxyyy, g_0_x_xyzz_xxxxyyyy, g_0_x_xyzz_xxxxyyyz, g_0_x_xyzz_xxxxyyz, g_0_x_xyzz_xxxxyyzz, g_0_x_xyzz_xxxxyzz, g_0_x_xyzz_xxxxyzzz, g_0_x_xyzz_xxxxzzz, g_0_x_xyzz_xxxyyyy, g_0_x_xyzz_xxxyyyyy, g_0_x_xyzz_xxxyyyyz, g_0_x_xyzz_xxxyyyz, g_0_x_xyzz_xxxyyyzz, g_0_x_xyzz_xxxyyzz, g_0_x_xyzz_xxxyyzzz, g_0_x_xyzz_xxxyzzz, g_0_x_xyzz_xxxyzzzz, g_0_x_xyzz_xxxzzzz, g_0_x_xyzz_xxyyyyy, g_0_x_xyzz_xxyyyyyy, g_0_x_xyzz_xxyyyyyz, g_0_x_xyzz_xxyyyyz, g_0_x_xyzz_xxyyyyzz, g_0_x_xyzz_xxyyyzz, g_0_x_xyzz_xxyyyzzz, g_0_x_xyzz_xxyyzzz, g_0_x_xyzz_xxyyzzzz, g_0_x_xyzz_xxyzzzz, g_0_x_xyzz_xxyzzzzz, g_0_x_xyzz_xxzzzzz, g_0_x_xyzz_xyyyyyy, g_0_x_xyzz_xyyyyyyy, g_0_x_xyzz_xyyyyyyz, g_0_x_xyzz_xyyyyyz, g_0_x_xyzz_xyyyyyzz, g_0_x_xyzz_xyyyyzz, g_0_x_xyzz_xyyyyzzz, g_0_x_xyzz_xyyyzzz, g_0_x_xyzz_xyyyzzzz, g_0_x_xyzz_xyyzzzz, g_0_x_xyzz_xyyzzzzz, g_0_x_xyzz_xyzzzzz, g_0_x_xyzz_xyzzzzzz, g_0_x_xyzz_xzzzzzz, g_0_x_xyzz_yyyyyyy, g_0_x_xyzz_yyyyyyyy, g_0_x_xyzz_yyyyyyyz, g_0_x_xyzz_yyyyyyz, g_0_x_xyzz_yyyyyyzz, g_0_x_xyzz_yyyyyzz, g_0_x_xyzz_yyyyyzzz, g_0_x_xyzz_yyyyzzz, g_0_x_xyzz_yyyyzzzz, g_0_x_xyzz_yyyzzzz, g_0_x_xyzz_yyyzzzzz, g_0_x_xyzz_yyzzzzz, g_0_x_xyzz_yyzzzzzz, g_0_x_xyzz_yzzzzzz, g_0_x_xyzz_yzzzzzzz, g_0_x_xyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzz_xxxxxxx[k] = -g_0_x_xyzz_xxxxxxx[k] * ab_y + g_0_x_xyzz_xxxxxxxy[k];

                g_0_x_xyyzz_xxxxxxy[k] = -g_0_x_xyzz_xxxxxxy[k] * ab_y + g_0_x_xyzz_xxxxxxyy[k];

                g_0_x_xyyzz_xxxxxxz[k] = -g_0_x_xyzz_xxxxxxz[k] * ab_y + g_0_x_xyzz_xxxxxxyz[k];

                g_0_x_xyyzz_xxxxxyy[k] = -g_0_x_xyzz_xxxxxyy[k] * ab_y + g_0_x_xyzz_xxxxxyyy[k];

                g_0_x_xyyzz_xxxxxyz[k] = -g_0_x_xyzz_xxxxxyz[k] * ab_y + g_0_x_xyzz_xxxxxyyz[k];

                g_0_x_xyyzz_xxxxxzz[k] = -g_0_x_xyzz_xxxxxzz[k] * ab_y + g_0_x_xyzz_xxxxxyzz[k];

                g_0_x_xyyzz_xxxxyyy[k] = -g_0_x_xyzz_xxxxyyy[k] * ab_y + g_0_x_xyzz_xxxxyyyy[k];

                g_0_x_xyyzz_xxxxyyz[k] = -g_0_x_xyzz_xxxxyyz[k] * ab_y + g_0_x_xyzz_xxxxyyyz[k];

                g_0_x_xyyzz_xxxxyzz[k] = -g_0_x_xyzz_xxxxyzz[k] * ab_y + g_0_x_xyzz_xxxxyyzz[k];

                g_0_x_xyyzz_xxxxzzz[k] = -g_0_x_xyzz_xxxxzzz[k] * ab_y + g_0_x_xyzz_xxxxyzzz[k];

                g_0_x_xyyzz_xxxyyyy[k] = -g_0_x_xyzz_xxxyyyy[k] * ab_y + g_0_x_xyzz_xxxyyyyy[k];

                g_0_x_xyyzz_xxxyyyz[k] = -g_0_x_xyzz_xxxyyyz[k] * ab_y + g_0_x_xyzz_xxxyyyyz[k];

                g_0_x_xyyzz_xxxyyzz[k] = -g_0_x_xyzz_xxxyyzz[k] * ab_y + g_0_x_xyzz_xxxyyyzz[k];

                g_0_x_xyyzz_xxxyzzz[k] = -g_0_x_xyzz_xxxyzzz[k] * ab_y + g_0_x_xyzz_xxxyyzzz[k];

                g_0_x_xyyzz_xxxzzzz[k] = -g_0_x_xyzz_xxxzzzz[k] * ab_y + g_0_x_xyzz_xxxyzzzz[k];

                g_0_x_xyyzz_xxyyyyy[k] = -g_0_x_xyzz_xxyyyyy[k] * ab_y + g_0_x_xyzz_xxyyyyyy[k];

                g_0_x_xyyzz_xxyyyyz[k] = -g_0_x_xyzz_xxyyyyz[k] * ab_y + g_0_x_xyzz_xxyyyyyz[k];

                g_0_x_xyyzz_xxyyyzz[k] = -g_0_x_xyzz_xxyyyzz[k] * ab_y + g_0_x_xyzz_xxyyyyzz[k];

                g_0_x_xyyzz_xxyyzzz[k] = -g_0_x_xyzz_xxyyzzz[k] * ab_y + g_0_x_xyzz_xxyyyzzz[k];

                g_0_x_xyyzz_xxyzzzz[k] = -g_0_x_xyzz_xxyzzzz[k] * ab_y + g_0_x_xyzz_xxyyzzzz[k];

                g_0_x_xyyzz_xxzzzzz[k] = -g_0_x_xyzz_xxzzzzz[k] * ab_y + g_0_x_xyzz_xxyzzzzz[k];

                g_0_x_xyyzz_xyyyyyy[k] = -g_0_x_xyzz_xyyyyyy[k] * ab_y + g_0_x_xyzz_xyyyyyyy[k];

                g_0_x_xyyzz_xyyyyyz[k] = -g_0_x_xyzz_xyyyyyz[k] * ab_y + g_0_x_xyzz_xyyyyyyz[k];

                g_0_x_xyyzz_xyyyyzz[k] = -g_0_x_xyzz_xyyyyzz[k] * ab_y + g_0_x_xyzz_xyyyyyzz[k];

                g_0_x_xyyzz_xyyyzzz[k] = -g_0_x_xyzz_xyyyzzz[k] * ab_y + g_0_x_xyzz_xyyyyzzz[k];

                g_0_x_xyyzz_xyyzzzz[k] = -g_0_x_xyzz_xyyzzzz[k] * ab_y + g_0_x_xyzz_xyyyzzzz[k];

                g_0_x_xyyzz_xyzzzzz[k] = -g_0_x_xyzz_xyzzzzz[k] * ab_y + g_0_x_xyzz_xyyzzzzz[k];

                g_0_x_xyyzz_xzzzzzz[k] = -g_0_x_xyzz_xzzzzzz[k] * ab_y + g_0_x_xyzz_xyzzzzzz[k];

                g_0_x_xyyzz_yyyyyyy[k] = -g_0_x_xyzz_yyyyyyy[k] * ab_y + g_0_x_xyzz_yyyyyyyy[k];

                g_0_x_xyyzz_yyyyyyz[k] = -g_0_x_xyzz_yyyyyyz[k] * ab_y + g_0_x_xyzz_yyyyyyyz[k];

                g_0_x_xyyzz_yyyyyzz[k] = -g_0_x_xyzz_yyyyyzz[k] * ab_y + g_0_x_xyzz_yyyyyyzz[k];

                g_0_x_xyyzz_yyyyzzz[k] = -g_0_x_xyzz_yyyyzzz[k] * ab_y + g_0_x_xyzz_yyyyyzzz[k];

                g_0_x_xyyzz_yyyzzzz[k] = -g_0_x_xyzz_yyyzzzz[k] * ab_y + g_0_x_xyzz_yyyyzzzz[k];

                g_0_x_xyyzz_yyzzzzz[k] = -g_0_x_xyzz_yyzzzzz[k] * ab_y + g_0_x_xyzz_yyyzzzzz[k];

                g_0_x_xyyzz_yzzzzzz[k] = -g_0_x_xyzz_yzzzzzz[k] * ab_y + g_0_x_xyzz_yyzzzzzz[k];

                g_0_x_xyyzz_zzzzzzz[k] = -g_0_x_xyzz_zzzzzzz[k] * ab_y + g_0_x_xyzz_yzzzzzzz[k];
            }

            /// Set up 468-504 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_x_xyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_x_xyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_x_xyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzz_xxxxxxx, g_0_x_xyzzz_xxxxxxy, g_0_x_xyzzz_xxxxxxz, g_0_x_xyzzz_xxxxxyy, g_0_x_xyzzz_xxxxxyz, g_0_x_xyzzz_xxxxxzz, g_0_x_xyzzz_xxxxyyy, g_0_x_xyzzz_xxxxyyz, g_0_x_xyzzz_xxxxyzz, g_0_x_xyzzz_xxxxzzz, g_0_x_xyzzz_xxxyyyy, g_0_x_xyzzz_xxxyyyz, g_0_x_xyzzz_xxxyyzz, g_0_x_xyzzz_xxxyzzz, g_0_x_xyzzz_xxxzzzz, g_0_x_xyzzz_xxyyyyy, g_0_x_xyzzz_xxyyyyz, g_0_x_xyzzz_xxyyyzz, g_0_x_xyzzz_xxyyzzz, g_0_x_xyzzz_xxyzzzz, g_0_x_xyzzz_xxzzzzz, g_0_x_xyzzz_xyyyyyy, g_0_x_xyzzz_xyyyyyz, g_0_x_xyzzz_xyyyyzz, g_0_x_xyzzz_xyyyzzz, g_0_x_xyzzz_xyyzzzz, g_0_x_xyzzz_xyzzzzz, g_0_x_xyzzz_xzzzzzz, g_0_x_xyzzz_yyyyyyy, g_0_x_xyzzz_yyyyyyz, g_0_x_xyzzz_yyyyyzz, g_0_x_xyzzz_yyyyzzz, g_0_x_xyzzz_yyyzzzz, g_0_x_xyzzz_yyzzzzz, g_0_x_xyzzz_yzzzzzz, g_0_x_xyzzz_zzzzzzz, g_0_x_xzzz_xxxxxxx, g_0_x_xzzz_xxxxxxxy, g_0_x_xzzz_xxxxxxy, g_0_x_xzzz_xxxxxxyy, g_0_x_xzzz_xxxxxxyz, g_0_x_xzzz_xxxxxxz, g_0_x_xzzz_xxxxxyy, g_0_x_xzzz_xxxxxyyy, g_0_x_xzzz_xxxxxyyz, g_0_x_xzzz_xxxxxyz, g_0_x_xzzz_xxxxxyzz, g_0_x_xzzz_xxxxxzz, g_0_x_xzzz_xxxxyyy, g_0_x_xzzz_xxxxyyyy, g_0_x_xzzz_xxxxyyyz, g_0_x_xzzz_xxxxyyz, g_0_x_xzzz_xxxxyyzz, g_0_x_xzzz_xxxxyzz, g_0_x_xzzz_xxxxyzzz, g_0_x_xzzz_xxxxzzz, g_0_x_xzzz_xxxyyyy, g_0_x_xzzz_xxxyyyyy, g_0_x_xzzz_xxxyyyyz, g_0_x_xzzz_xxxyyyz, g_0_x_xzzz_xxxyyyzz, g_0_x_xzzz_xxxyyzz, g_0_x_xzzz_xxxyyzzz, g_0_x_xzzz_xxxyzzz, g_0_x_xzzz_xxxyzzzz, g_0_x_xzzz_xxxzzzz, g_0_x_xzzz_xxyyyyy, g_0_x_xzzz_xxyyyyyy, g_0_x_xzzz_xxyyyyyz, g_0_x_xzzz_xxyyyyz, g_0_x_xzzz_xxyyyyzz, g_0_x_xzzz_xxyyyzz, g_0_x_xzzz_xxyyyzzz, g_0_x_xzzz_xxyyzzz, g_0_x_xzzz_xxyyzzzz, g_0_x_xzzz_xxyzzzz, g_0_x_xzzz_xxyzzzzz, g_0_x_xzzz_xxzzzzz, g_0_x_xzzz_xyyyyyy, g_0_x_xzzz_xyyyyyyy, g_0_x_xzzz_xyyyyyyz, g_0_x_xzzz_xyyyyyz, g_0_x_xzzz_xyyyyyzz, g_0_x_xzzz_xyyyyzz, g_0_x_xzzz_xyyyyzzz, g_0_x_xzzz_xyyyzzz, g_0_x_xzzz_xyyyzzzz, g_0_x_xzzz_xyyzzzz, g_0_x_xzzz_xyyzzzzz, g_0_x_xzzz_xyzzzzz, g_0_x_xzzz_xyzzzzzz, g_0_x_xzzz_xzzzzzz, g_0_x_xzzz_yyyyyyy, g_0_x_xzzz_yyyyyyyy, g_0_x_xzzz_yyyyyyyz, g_0_x_xzzz_yyyyyyz, g_0_x_xzzz_yyyyyyzz, g_0_x_xzzz_yyyyyzz, g_0_x_xzzz_yyyyyzzz, g_0_x_xzzz_yyyyzzz, g_0_x_xzzz_yyyyzzzz, g_0_x_xzzz_yyyzzzz, g_0_x_xzzz_yyyzzzzz, g_0_x_xzzz_yyzzzzz, g_0_x_xzzz_yyzzzzzz, g_0_x_xzzz_yzzzzzz, g_0_x_xzzz_yzzzzzzz, g_0_x_xzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzz_xxxxxxx[k] = -g_0_x_xzzz_xxxxxxx[k] * ab_y + g_0_x_xzzz_xxxxxxxy[k];

                g_0_x_xyzzz_xxxxxxy[k] = -g_0_x_xzzz_xxxxxxy[k] * ab_y + g_0_x_xzzz_xxxxxxyy[k];

                g_0_x_xyzzz_xxxxxxz[k] = -g_0_x_xzzz_xxxxxxz[k] * ab_y + g_0_x_xzzz_xxxxxxyz[k];

                g_0_x_xyzzz_xxxxxyy[k] = -g_0_x_xzzz_xxxxxyy[k] * ab_y + g_0_x_xzzz_xxxxxyyy[k];

                g_0_x_xyzzz_xxxxxyz[k] = -g_0_x_xzzz_xxxxxyz[k] * ab_y + g_0_x_xzzz_xxxxxyyz[k];

                g_0_x_xyzzz_xxxxxzz[k] = -g_0_x_xzzz_xxxxxzz[k] * ab_y + g_0_x_xzzz_xxxxxyzz[k];

                g_0_x_xyzzz_xxxxyyy[k] = -g_0_x_xzzz_xxxxyyy[k] * ab_y + g_0_x_xzzz_xxxxyyyy[k];

                g_0_x_xyzzz_xxxxyyz[k] = -g_0_x_xzzz_xxxxyyz[k] * ab_y + g_0_x_xzzz_xxxxyyyz[k];

                g_0_x_xyzzz_xxxxyzz[k] = -g_0_x_xzzz_xxxxyzz[k] * ab_y + g_0_x_xzzz_xxxxyyzz[k];

                g_0_x_xyzzz_xxxxzzz[k] = -g_0_x_xzzz_xxxxzzz[k] * ab_y + g_0_x_xzzz_xxxxyzzz[k];

                g_0_x_xyzzz_xxxyyyy[k] = -g_0_x_xzzz_xxxyyyy[k] * ab_y + g_0_x_xzzz_xxxyyyyy[k];

                g_0_x_xyzzz_xxxyyyz[k] = -g_0_x_xzzz_xxxyyyz[k] * ab_y + g_0_x_xzzz_xxxyyyyz[k];

                g_0_x_xyzzz_xxxyyzz[k] = -g_0_x_xzzz_xxxyyzz[k] * ab_y + g_0_x_xzzz_xxxyyyzz[k];

                g_0_x_xyzzz_xxxyzzz[k] = -g_0_x_xzzz_xxxyzzz[k] * ab_y + g_0_x_xzzz_xxxyyzzz[k];

                g_0_x_xyzzz_xxxzzzz[k] = -g_0_x_xzzz_xxxzzzz[k] * ab_y + g_0_x_xzzz_xxxyzzzz[k];

                g_0_x_xyzzz_xxyyyyy[k] = -g_0_x_xzzz_xxyyyyy[k] * ab_y + g_0_x_xzzz_xxyyyyyy[k];

                g_0_x_xyzzz_xxyyyyz[k] = -g_0_x_xzzz_xxyyyyz[k] * ab_y + g_0_x_xzzz_xxyyyyyz[k];

                g_0_x_xyzzz_xxyyyzz[k] = -g_0_x_xzzz_xxyyyzz[k] * ab_y + g_0_x_xzzz_xxyyyyzz[k];

                g_0_x_xyzzz_xxyyzzz[k] = -g_0_x_xzzz_xxyyzzz[k] * ab_y + g_0_x_xzzz_xxyyyzzz[k];

                g_0_x_xyzzz_xxyzzzz[k] = -g_0_x_xzzz_xxyzzzz[k] * ab_y + g_0_x_xzzz_xxyyzzzz[k];

                g_0_x_xyzzz_xxzzzzz[k] = -g_0_x_xzzz_xxzzzzz[k] * ab_y + g_0_x_xzzz_xxyzzzzz[k];

                g_0_x_xyzzz_xyyyyyy[k] = -g_0_x_xzzz_xyyyyyy[k] * ab_y + g_0_x_xzzz_xyyyyyyy[k];

                g_0_x_xyzzz_xyyyyyz[k] = -g_0_x_xzzz_xyyyyyz[k] * ab_y + g_0_x_xzzz_xyyyyyyz[k];

                g_0_x_xyzzz_xyyyyzz[k] = -g_0_x_xzzz_xyyyyzz[k] * ab_y + g_0_x_xzzz_xyyyyyzz[k];

                g_0_x_xyzzz_xyyyzzz[k] = -g_0_x_xzzz_xyyyzzz[k] * ab_y + g_0_x_xzzz_xyyyyzzz[k];

                g_0_x_xyzzz_xyyzzzz[k] = -g_0_x_xzzz_xyyzzzz[k] * ab_y + g_0_x_xzzz_xyyyzzzz[k];

                g_0_x_xyzzz_xyzzzzz[k] = -g_0_x_xzzz_xyzzzzz[k] * ab_y + g_0_x_xzzz_xyyzzzzz[k];

                g_0_x_xyzzz_xzzzzzz[k] = -g_0_x_xzzz_xzzzzzz[k] * ab_y + g_0_x_xzzz_xyzzzzzz[k];

                g_0_x_xyzzz_yyyyyyy[k] = -g_0_x_xzzz_yyyyyyy[k] * ab_y + g_0_x_xzzz_yyyyyyyy[k];

                g_0_x_xyzzz_yyyyyyz[k] = -g_0_x_xzzz_yyyyyyz[k] * ab_y + g_0_x_xzzz_yyyyyyyz[k];

                g_0_x_xyzzz_yyyyyzz[k] = -g_0_x_xzzz_yyyyyzz[k] * ab_y + g_0_x_xzzz_yyyyyyzz[k];

                g_0_x_xyzzz_yyyyzzz[k] = -g_0_x_xzzz_yyyyzzz[k] * ab_y + g_0_x_xzzz_yyyyyzzz[k];

                g_0_x_xyzzz_yyyzzzz[k] = -g_0_x_xzzz_yyyzzzz[k] * ab_y + g_0_x_xzzz_yyyyzzzz[k];

                g_0_x_xyzzz_yyzzzzz[k] = -g_0_x_xzzz_yyzzzzz[k] * ab_y + g_0_x_xzzz_yyyzzzzz[k];

                g_0_x_xyzzz_yzzzzzz[k] = -g_0_x_xzzz_yzzzzzz[k] * ab_y + g_0_x_xzzz_yyzzzzzz[k];

                g_0_x_xyzzz_zzzzzzz[k] = -g_0_x_xzzz_zzzzzzz[k] * ab_y + g_0_x_xzzz_yzzzzzzz[k];
            }

            /// Set up 504-540 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_x_xzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_x_xzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_x_xzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzz_xxxxxxx, g_0_x_xzzz_xxxxxxxz, g_0_x_xzzz_xxxxxxy, g_0_x_xzzz_xxxxxxyz, g_0_x_xzzz_xxxxxxz, g_0_x_xzzz_xxxxxxzz, g_0_x_xzzz_xxxxxyy, g_0_x_xzzz_xxxxxyyz, g_0_x_xzzz_xxxxxyz, g_0_x_xzzz_xxxxxyzz, g_0_x_xzzz_xxxxxzz, g_0_x_xzzz_xxxxxzzz, g_0_x_xzzz_xxxxyyy, g_0_x_xzzz_xxxxyyyz, g_0_x_xzzz_xxxxyyz, g_0_x_xzzz_xxxxyyzz, g_0_x_xzzz_xxxxyzz, g_0_x_xzzz_xxxxyzzz, g_0_x_xzzz_xxxxzzz, g_0_x_xzzz_xxxxzzzz, g_0_x_xzzz_xxxyyyy, g_0_x_xzzz_xxxyyyyz, g_0_x_xzzz_xxxyyyz, g_0_x_xzzz_xxxyyyzz, g_0_x_xzzz_xxxyyzz, g_0_x_xzzz_xxxyyzzz, g_0_x_xzzz_xxxyzzz, g_0_x_xzzz_xxxyzzzz, g_0_x_xzzz_xxxzzzz, g_0_x_xzzz_xxxzzzzz, g_0_x_xzzz_xxyyyyy, g_0_x_xzzz_xxyyyyyz, g_0_x_xzzz_xxyyyyz, g_0_x_xzzz_xxyyyyzz, g_0_x_xzzz_xxyyyzz, g_0_x_xzzz_xxyyyzzz, g_0_x_xzzz_xxyyzzz, g_0_x_xzzz_xxyyzzzz, g_0_x_xzzz_xxyzzzz, g_0_x_xzzz_xxyzzzzz, g_0_x_xzzz_xxzzzzz, g_0_x_xzzz_xxzzzzzz, g_0_x_xzzz_xyyyyyy, g_0_x_xzzz_xyyyyyyz, g_0_x_xzzz_xyyyyyz, g_0_x_xzzz_xyyyyyzz, g_0_x_xzzz_xyyyyzz, g_0_x_xzzz_xyyyyzzz, g_0_x_xzzz_xyyyzzz, g_0_x_xzzz_xyyyzzzz, g_0_x_xzzz_xyyzzzz, g_0_x_xzzz_xyyzzzzz, g_0_x_xzzz_xyzzzzz, g_0_x_xzzz_xyzzzzzz, g_0_x_xzzz_xzzzzzz, g_0_x_xzzz_xzzzzzzz, g_0_x_xzzz_yyyyyyy, g_0_x_xzzz_yyyyyyyz, g_0_x_xzzz_yyyyyyz, g_0_x_xzzz_yyyyyyzz, g_0_x_xzzz_yyyyyzz, g_0_x_xzzz_yyyyyzzz, g_0_x_xzzz_yyyyzzz, g_0_x_xzzz_yyyyzzzz, g_0_x_xzzz_yyyzzzz, g_0_x_xzzz_yyyzzzzz, g_0_x_xzzz_yyzzzzz, g_0_x_xzzz_yyzzzzzz, g_0_x_xzzz_yzzzzzz, g_0_x_xzzz_yzzzzzzz, g_0_x_xzzz_zzzzzzz, g_0_x_xzzz_zzzzzzzz, g_0_x_xzzzz_xxxxxxx, g_0_x_xzzzz_xxxxxxy, g_0_x_xzzzz_xxxxxxz, g_0_x_xzzzz_xxxxxyy, g_0_x_xzzzz_xxxxxyz, g_0_x_xzzzz_xxxxxzz, g_0_x_xzzzz_xxxxyyy, g_0_x_xzzzz_xxxxyyz, g_0_x_xzzzz_xxxxyzz, g_0_x_xzzzz_xxxxzzz, g_0_x_xzzzz_xxxyyyy, g_0_x_xzzzz_xxxyyyz, g_0_x_xzzzz_xxxyyzz, g_0_x_xzzzz_xxxyzzz, g_0_x_xzzzz_xxxzzzz, g_0_x_xzzzz_xxyyyyy, g_0_x_xzzzz_xxyyyyz, g_0_x_xzzzz_xxyyyzz, g_0_x_xzzzz_xxyyzzz, g_0_x_xzzzz_xxyzzzz, g_0_x_xzzzz_xxzzzzz, g_0_x_xzzzz_xyyyyyy, g_0_x_xzzzz_xyyyyyz, g_0_x_xzzzz_xyyyyzz, g_0_x_xzzzz_xyyyzzz, g_0_x_xzzzz_xyyzzzz, g_0_x_xzzzz_xyzzzzz, g_0_x_xzzzz_xzzzzzz, g_0_x_xzzzz_yyyyyyy, g_0_x_xzzzz_yyyyyyz, g_0_x_xzzzz_yyyyyzz, g_0_x_xzzzz_yyyyzzz, g_0_x_xzzzz_yyyzzzz, g_0_x_xzzzz_yyzzzzz, g_0_x_xzzzz_yzzzzzz, g_0_x_xzzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzz_xxxxxxx[k] = -g_0_x_xzzz_xxxxxxx[k] * ab_z + g_0_x_xzzz_xxxxxxxz[k];

                g_0_x_xzzzz_xxxxxxy[k] = -g_0_x_xzzz_xxxxxxy[k] * ab_z + g_0_x_xzzz_xxxxxxyz[k];

                g_0_x_xzzzz_xxxxxxz[k] = -g_0_x_xzzz_xxxxxxz[k] * ab_z + g_0_x_xzzz_xxxxxxzz[k];

                g_0_x_xzzzz_xxxxxyy[k] = -g_0_x_xzzz_xxxxxyy[k] * ab_z + g_0_x_xzzz_xxxxxyyz[k];

                g_0_x_xzzzz_xxxxxyz[k] = -g_0_x_xzzz_xxxxxyz[k] * ab_z + g_0_x_xzzz_xxxxxyzz[k];

                g_0_x_xzzzz_xxxxxzz[k] = -g_0_x_xzzz_xxxxxzz[k] * ab_z + g_0_x_xzzz_xxxxxzzz[k];

                g_0_x_xzzzz_xxxxyyy[k] = -g_0_x_xzzz_xxxxyyy[k] * ab_z + g_0_x_xzzz_xxxxyyyz[k];

                g_0_x_xzzzz_xxxxyyz[k] = -g_0_x_xzzz_xxxxyyz[k] * ab_z + g_0_x_xzzz_xxxxyyzz[k];

                g_0_x_xzzzz_xxxxyzz[k] = -g_0_x_xzzz_xxxxyzz[k] * ab_z + g_0_x_xzzz_xxxxyzzz[k];

                g_0_x_xzzzz_xxxxzzz[k] = -g_0_x_xzzz_xxxxzzz[k] * ab_z + g_0_x_xzzz_xxxxzzzz[k];

                g_0_x_xzzzz_xxxyyyy[k] = -g_0_x_xzzz_xxxyyyy[k] * ab_z + g_0_x_xzzz_xxxyyyyz[k];

                g_0_x_xzzzz_xxxyyyz[k] = -g_0_x_xzzz_xxxyyyz[k] * ab_z + g_0_x_xzzz_xxxyyyzz[k];

                g_0_x_xzzzz_xxxyyzz[k] = -g_0_x_xzzz_xxxyyzz[k] * ab_z + g_0_x_xzzz_xxxyyzzz[k];

                g_0_x_xzzzz_xxxyzzz[k] = -g_0_x_xzzz_xxxyzzz[k] * ab_z + g_0_x_xzzz_xxxyzzzz[k];

                g_0_x_xzzzz_xxxzzzz[k] = -g_0_x_xzzz_xxxzzzz[k] * ab_z + g_0_x_xzzz_xxxzzzzz[k];

                g_0_x_xzzzz_xxyyyyy[k] = -g_0_x_xzzz_xxyyyyy[k] * ab_z + g_0_x_xzzz_xxyyyyyz[k];

                g_0_x_xzzzz_xxyyyyz[k] = -g_0_x_xzzz_xxyyyyz[k] * ab_z + g_0_x_xzzz_xxyyyyzz[k];

                g_0_x_xzzzz_xxyyyzz[k] = -g_0_x_xzzz_xxyyyzz[k] * ab_z + g_0_x_xzzz_xxyyyzzz[k];

                g_0_x_xzzzz_xxyyzzz[k] = -g_0_x_xzzz_xxyyzzz[k] * ab_z + g_0_x_xzzz_xxyyzzzz[k];

                g_0_x_xzzzz_xxyzzzz[k] = -g_0_x_xzzz_xxyzzzz[k] * ab_z + g_0_x_xzzz_xxyzzzzz[k];

                g_0_x_xzzzz_xxzzzzz[k] = -g_0_x_xzzz_xxzzzzz[k] * ab_z + g_0_x_xzzz_xxzzzzzz[k];

                g_0_x_xzzzz_xyyyyyy[k] = -g_0_x_xzzz_xyyyyyy[k] * ab_z + g_0_x_xzzz_xyyyyyyz[k];

                g_0_x_xzzzz_xyyyyyz[k] = -g_0_x_xzzz_xyyyyyz[k] * ab_z + g_0_x_xzzz_xyyyyyzz[k];

                g_0_x_xzzzz_xyyyyzz[k] = -g_0_x_xzzz_xyyyyzz[k] * ab_z + g_0_x_xzzz_xyyyyzzz[k];

                g_0_x_xzzzz_xyyyzzz[k] = -g_0_x_xzzz_xyyyzzz[k] * ab_z + g_0_x_xzzz_xyyyzzzz[k];

                g_0_x_xzzzz_xyyzzzz[k] = -g_0_x_xzzz_xyyzzzz[k] * ab_z + g_0_x_xzzz_xyyzzzzz[k];

                g_0_x_xzzzz_xyzzzzz[k] = -g_0_x_xzzz_xyzzzzz[k] * ab_z + g_0_x_xzzz_xyzzzzzz[k];

                g_0_x_xzzzz_xzzzzzz[k] = -g_0_x_xzzz_xzzzzzz[k] * ab_z + g_0_x_xzzz_xzzzzzzz[k];

                g_0_x_xzzzz_yyyyyyy[k] = -g_0_x_xzzz_yyyyyyy[k] * ab_z + g_0_x_xzzz_yyyyyyyz[k];

                g_0_x_xzzzz_yyyyyyz[k] = -g_0_x_xzzz_yyyyyyz[k] * ab_z + g_0_x_xzzz_yyyyyyzz[k];

                g_0_x_xzzzz_yyyyyzz[k] = -g_0_x_xzzz_yyyyyzz[k] * ab_z + g_0_x_xzzz_yyyyyzzz[k];

                g_0_x_xzzzz_yyyyzzz[k] = -g_0_x_xzzz_yyyyzzz[k] * ab_z + g_0_x_xzzz_yyyyzzzz[k];

                g_0_x_xzzzz_yyyzzzz[k] = -g_0_x_xzzz_yyyzzzz[k] * ab_z + g_0_x_xzzz_yyyzzzzz[k];

                g_0_x_xzzzz_yyzzzzz[k] = -g_0_x_xzzz_yyzzzzz[k] * ab_z + g_0_x_xzzz_yyzzzzzz[k];

                g_0_x_xzzzz_yzzzzzz[k] = -g_0_x_xzzz_yzzzzzz[k] * ab_z + g_0_x_xzzz_yzzzzzzz[k];

                g_0_x_xzzzz_zzzzzzz[k] = -g_0_x_xzzz_zzzzzzz[k] * ab_z + g_0_x_xzzz_zzzzzzzz[k];
            }

            /// Set up 540-576 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_x_yyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_x_yyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_x_yyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_xxxxxxx, g_0_x_yyyy_xxxxxxxy, g_0_x_yyyy_xxxxxxy, g_0_x_yyyy_xxxxxxyy, g_0_x_yyyy_xxxxxxyz, g_0_x_yyyy_xxxxxxz, g_0_x_yyyy_xxxxxyy, g_0_x_yyyy_xxxxxyyy, g_0_x_yyyy_xxxxxyyz, g_0_x_yyyy_xxxxxyz, g_0_x_yyyy_xxxxxyzz, g_0_x_yyyy_xxxxxzz, g_0_x_yyyy_xxxxyyy, g_0_x_yyyy_xxxxyyyy, g_0_x_yyyy_xxxxyyyz, g_0_x_yyyy_xxxxyyz, g_0_x_yyyy_xxxxyyzz, g_0_x_yyyy_xxxxyzz, g_0_x_yyyy_xxxxyzzz, g_0_x_yyyy_xxxxzzz, g_0_x_yyyy_xxxyyyy, g_0_x_yyyy_xxxyyyyy, g_0_x_yyyy_xxxyyyyz, g_0_x_yyyy_xxxyyyz, g_0_x_yyyy_xxxyyyzz, g_0_x_yyyy_xxxyyzz, g_0_x_yyyy_xxxyyzzz, g_0_x_yyyy_xxxyzzz, g_0_x_yyyy_xxxyzzzz, g_0_x_yyyy_xxxzzzz, g_0_x_yyyy_xxyyyyy, g_0_x_yyyy_xxyyyyyy, g_0_x_yyyy_xxyyyyyz, g_0_x_yyyy_xxyyyyz, g_0_x_yyyy_xxyyyyzz, g_0_x_yyyy_xxyyyzz, g_0_x_yyyy_xxyyyzzz, g_0_x_yyyy_xxyyzzz, g_0_x_yyyy_xxyyzzzz, g_0_x_yyyy_xxyzzzz, g_0_x_yyyy_xxyzzzzz, g_0_x_yyyy_xxzzzzz, g_0_x_yyyy_xyyyyyy, g_0_x_yyyy_xyyyyyyy, g_0_x_yyyy_xyyyyyyz, g_0_x_yyyy_xyyyyyz, g_0_x_yyyy_xyyyyyzz, g_0_x_yyyy_xyyyyzz, g_0_x_yyyy_xyyyyzzz, g_0_x_yyyy_xyyyzzz, g_0_x_yyyy_xyyyzzzz, g_0_x_yyyy_xyyzzzz, g_0_x_yyyy_xyyzzzzz, g_0_x_yyyy_xyzzzzz, g_0_x_yyyy_xyzzzzzz, g_0_x_yyyy_xzzzzzz, g_0_x_yyyy_yyyyyyy, g_0_x_yyyy_yyyyyyyy, g_0_x_yyyy_yyyyyyyz, g_0_x_yyyy_yyyyyyz, g_0_x_yyyy_yyyyyyzz, g_0_x_yyyy_yyyyyzz, g_0_x_yyyy_yyyyyzzz, g_0_x_yyyy_yyyyzzz, g_0_x_yyyy_yyyyzzzz, g_0_x_yyyy_yyyzzzz, g_0_x_yyyy_yyyzzzzz, g_0_x_yyyy_yyzzzzz, g_0_x_yyyy_yyzzzzzz, g_0_x_yyyy_yzzzzzz, g_0_x_yyyy_yzzzzzzz, g_0_x_yyyy_zzzzzzz, g_0_x_yyyyy_xxxxxxx, g_0_x_yyyyy_xxxxxxy, g_0_x_yyyyy_xxxxxxz, g_0_x_yyyyy_xxxxxyy, g_0_x_yyyyy_xxxxxyz, g_0_x_yyyyy_xxxxxzz, g_0_x_yyyyy_xxxxyyy, g_0_x_yyyyy_xxxxyyz, g_0_x_yyyyy_xxxxyzz, g_0_x_yyyyy_xxxxzzz, g_0_x_yyyyy_xxxyyyy, g_0_x_yyyyy_xxxyyyz, g_0_x_yyyyy_xxxyyzz, g_0_x_yyyyy_xxxyzzz, g_0_x_yyyyy_xxxzzzz, g_0_x_yyyyy_xxyyyyy, g_0_x_yyyyy_xxyyyyz, g_0_x_yyyyy_xxyyyzz, g_0_x_yyyyy_xxyyzzz, g_0_x_yyyyy_xxyzzzz, g_0_x_yyyyy_xxzzzzz, g_0_x_yyyyy_xyyyyyy, g_0_x_yyyyy_xyyyyyz, g_0_x_yyyyy_xyyyyzz, g_0_x_yyyyy_xyyyzzz, g_0_x_yyyyy_xyyzzzz, g_0_x_yyyyy_xyzzzzz, g_0_x_yyyyy_xzzzzzz, g_0_x_yyyyy_yyyyyyy, g_0_x_yyyyy_yyyyyyz, g_0_x_yyyyy_yyyyyzz, g_0_x_yyyyy_yyyyzzz, g_0_x_yyyyy_yyyzzzz, g_0_x_yyyyy_yyzzzzz, g_0_x_yyyyy_yzzzzzz, g_0_x_yyyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyy_xxxxxxx[k] = -g_0_x_yyyy_xxxxxxx[k] * ab_y + g_0_x_yyyy_xxxxxxxy[k];

                g_0_x_yyyyy_xxxxxxy[k] = -g_0_x_yyyy_xxxxxxy[k] * ab_y + g_0_x_yyyy_xxxxxxyy[k];

                g_0_x_yyyyy_xxxxxxz[k] = -g_0_x_yyyy_xxxxxxz[k] * ab_y + g_0_x_yyyy_xxxxxxyz[k];

                g_0_x_yyyyy_xxxxxyy[k] = -g_0_x_yyyy_xxxxxyy[k] * ab_y + g_0_x_yyyy_xxxxxyyy[k];

                g_0_x_yyyyy_xxxxxyz[k] = -g_0_x_yyyy_xxxxxyz[k] * ab_y + g_0_x_yyyy_xxxxxyyz[k];

                g_0_x_yyyyy_xxxxxzz[k] = -g_0_x_yyyy_xxxxxzz[k] * ab_y + g_0_x_yyyy_xxxxxyzz[k];

                g_0_x_yyyyy_xxxxyyy[k] = -g_0_x_yyyy_xxxxyyy[k] * ab_y + g_0_x_yyyy_xxxxyyyy[k];

                g_0_x_yyyyy_xxxxyyz[k] = -g_0_x_yyyy_xxxxyyz[k] * ab_y + g_0_x_yyyy_xxxxyyyz[k];

                g_0_x_yyyyy_xxxxyzz[k] = -g_0_x_yyyy_xxxxyzz[k] * ab_y + g_0_x_yyyy_xxxxyyzz[k];

                g_0_x_yyyyy_xxxxzzz[k] = -g_0_x_yyyy_xxxxzzz[k] * ab_y + g_0_x_yyyy_xxxxyzzz[k];

                g_0_x_yyyyy_xxxyyyy[k] = -g_0_x_yyyy_xxxyyyy[k] * ab_y + g_0_x_yyyy_xxxyyyyy[k];

                g_0_x_yyyyy_xxxyyyz[k] = -g_0_x_yyyy_xxxyyyz[k] * ab_y + g_0_x_yyyy_xxxyyyyz[k];

                g_0_x_yyyyy_xxxyyzz[k] = -g_0_x_yyyy_xxxyyzz[k] * ab_y + g_0_x_yyyy_xxxyyyzz[k];

                g_0_x_yyyyy_xxxyzzz[k] = -g_0_x_yyyy_xxxyzzz[k] * ab_y + g_0_x_yyyy_xxxyyzzz[k];

                g_0_x_yyyyy_xxxzzzz[k] = -g_0_x_yyyy_xxxzzzz[k] * ab_y + g_0_x_yyyy_xxxyzzzz[k];

                g_0_x_yyyyy_xxyyyyy[k] = -g_0_x_yyyy_xxyyyyy[k] * ab_y + g_0_x_yyyy_xxyyyyyy[k];

                g_0_x_yyyyy_xxyyyyz[k] = -g_0_x_yyyy_xxyyyyz[k] * ab_y + g_0_x_yyyy_xxyyyyyz[k];

                g_0_x_yyyyy_xxyyyzz[k] = -g_0_x_yyyy_xxyyyzz[k] * ab_y + g_0_x_yyyy_xxyyyyzz[k];

                g_0_x_yyyyy_xxyyzzz[k] = -g_0_x_yyyy_xxyyzzz[k] * ab_y + g_0_x_yyyy_xxyyyzzz[k];

                g_0_x_yyyyy_xxyzzzz[k] = -g_0_x_yyyy_xxyzzzz[k] * ab_y + g_0_x_yyyy_xxyyzzzz[k];

                g_0_x_yyyyy_xxzzzzz[k] = -g_0_x_yyyy_xxzzzzz[k] * ab_y + g_0_x_yyyy_xxyzzzzz[k];

                g_0_x_yyyyy_xyyyyyy[k] = -g_0_x_yyyy_xyyyyyy[k] * ab_y + g_0_x_yyyy_xyyyyyyy[k];

                g_0_x_yyyyy_xyyyyyz[k] = -g_0_x_yyyy_xyyyyyz[k] * ab_y + g_0_x_yyyy_xyyyyyyz[k];

                g_0_x_yyyyy_xyyyyzz[k] = -g_0_x_yyyy_xyyyyzz[k] * ab_y + g_0_x_yyyy_xyyyyyzz[k];

                g_0_x_yyyyy_xyyyzzz[k] = -g_0_x_yyyy_xyyyzzz[k] * ab_y + g_0_x_yyyy_xyyyyzzz[k];

                g_0_x_yyyyy_xyyzzzz[k] = -g_0_x_yyyy_xyyzzzz[k] * ab_y + g_0_x_yyyy_xyyyzzzz[k];

                g_0_x_yyyyy_xyzzzzz[k] = -g_0_x_yyyy_xyzzzzz[k] * ab_y + g_0_x_yyyy_xyyzzzzz[k];

                g_0_x_yyyyy_xzzzzzz[k] = -g_0_x_yyyy_xzzzzzz[k] * ab_y + g_0_x_yyyy_xyzzzzzz[k];

                g_0_x_yyyyy_yyyyyyy[k] = -g_0_x_yyyy_yyyyyyy[k] * ab_y + g_0_x_yyyy_yyyyyyyy[k];

                g_0_x_yyyyy_yyyyyyz[k] = -g_0_x_yyyy_yyyyyyz[k] * ab_y + g_0_x_yyyy_yyyyyyyz[k];

                g_0_x_yyyyy_yyyyyzz[k] = -g_0_x_yyyy_yyyyyzz[k] * ab_y + g_0_x_yyyy_yyyyyyzz[k];

                g_0_x_yyyyy_yyyyzzz[k] = -g_0_x_yyyy_yyyyzzz[k] * ab_y + g_0_x_yyyy_yyyyyzzz[k];

                g_0_x_yyyyy_yyyzzzz[k] = -g_0_x_yyyy_yyyzzzz[k] * ab_y + g_0_x_yyyy_yyyyzzzz[k];

                g_0_x_yyyyy_yyzzzzz[k] = -g_0_x_yyyy_yyzzzzz[k] * ab_y + g_0_x_yyyy_yyyzzzzz[k];

                g_0_x_yyyyy_yzzzzzz[k] = -g_0_x_yyyy_yzzzzzz[k] * ab_y + g_0_x_yyyy_yyzzzzzz[k];

                g_0_x_yyyyy_zzzzzzz[k] = -g_0_x_yyyy_zzzzzzz[k] * ab_y + g_0_x_yyyy_yzzzzzzz[k];
            }

            /// Set up 576-612 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_x_yyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_x_yyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_x_yyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyz_xxxxxxx, g_0_x_yyyyz_xxxxxxy, g_0_x_yyyyz_xxxxxxz, g_0_x_yyyyz_xxxxxyy, g_0_x_yyyyz_xxxxxyz, g_0_x_yyyyz_xxxxxzz, g_0_x_yyyyz_xxxxyyy, g_0_x_yyyyz_xxxxyyz, g_0_x_yyyyz_xxxxyzz, g_0_x_yyyyz_xxxxzzz, g_0_x_yyyyz_xxxyyyy, g_0_x_yyyyz_xxxyyyz, g_0_x_yyyyz_xxxyyzz, g_0_x_yyyyz_xxxyzzz, g_0_x_yyyyz_xxxzzzz, g_0_x_yyyyz_xxyyyyy, g_0_x_yyyyz_xxyyyyz, g_0_x_yyyyz_xxyyyzz, g_0_x_yyyyz_xxyyzzz, g_0_x_yyyyz_xxyzzzz, g_0_x_yyyyz_xxzzzzz, g_0_x_yyyyz_xyyyyyy, g_0_x_yyyyz_xyyyyyz, g_0_x_yyyyz_xyyyyzz, g_0_x_yyyyz_xyyyzzz, g_0_x_yyyyz_xyyzzzz, g_0_x_yyyyz_xyzzzzz, g_0_x_yyyyz_xzzzzzz, g_0_x_yyyyz_yyyyyyy, g_0_x_yyyyz_yyyyyyz, g_0_x_yyyyz_yyyyyzz, g_0_x_yyyyz_yyyyzzz, g_0_x_yyyyz_yyyzzzz, g_0_x_yyyyz_yyzzzzz, g_0_x_yyyyz_yzzzzzz, g_0_x_yyyyz_zzzzzzz, g_0_x_yyyz_xxxxxxx, g_0_x_yyyz_xxxxxxxy, g_0_x_yyyz_xxxxxxy, g_0_x_yyyz_xxxxxxyy, g_0_x_yyyz_xxxxxxyz, g_0_x_yyyz_xxxxxxz, g_0_x_yyyz_xxxxxyy, g_0_x_yyyz_xxxxxyyy, g_0_x_yyyz_xxxxxyyz, g_0_x_yyyz_xxxxxyz, g_0_x_yyyz_xxxxxyzz, g_0_x_yyyz_xxxxxzz, g_0_x_yyyz_xxxxyyy, g_0_x_yyyz_xxxxyyyy, g_0_x_yyyz_xxxxyyyz, g_0_x_yyyz_xxxxyyz, g_0_x_yyyz_xxxxyyzz, g_0_x_yyyz_xxxxyzz, g_0_x_yyyz_xxxxyzzz, g_0_x_yyyz_xxxxzzz, g_0_x_yyyz_xxxyyyy, g_0_x_yyyz_xxxyyyyy, g_0_x_yyyz_xxxyyyyz, g_0_x_yyyz_xxxyyyz, g_0_x_yyyz_xxxyyyzz, g_0_x_yyyz_xxxyyzz, g_0_x_yyyz_xxxyyzzz, g_0_x_yyyz_xxxyzzz, g_0_x_yyyz_xxxyzzzz, g_0_x_yyyz_xxxzzzz, g_0_x_yyyz_xxyyyyy, g_0_x_yyyz_xxyyyyyy, g_0_x_yyyz_xxyyyyyz, g_0_x_yyyz_xxyyyyz, g_0_x_yyyz_xxyyyyzz, g_0_x_yyyz_xxyyyzz, g_0_x_yyyz_xxyyyzzz, g_0_x_yyyz_xxyyzzz, g_0_x_yyyz_xxyyzzzz, g_0_x_yyyz_xxyzzzz, g_0_x_yyyz_xxyzzzzz, g_0_x_yyyz_xxzzzzz, g_0_x_yyyz_xyyyyyy, g_0_x_yyyz_xyyyyyyy, g_0_x_yyyz_xyyyyyyz, g_0_x_yyyz_xyyyyyz, g_0_x_yyyz_xyyyyyzz, g_0_x_yyyz_xyyyyzz, g_0_x_yyyz_xyyyyzzz, g_0_x_yyyz_xyyyzzz, g_0_x_yyyz_xyyyzzzz, g_0_x_yyyz_xyyzzzz, g_0_x_yyyz_xyyzzzzz, g_0_x_yyyz_xyzzzzz, g_0_x_yyyz_xyzzzzzz, g_0_x_yyyz_xzzzzzz, g_0_x_yyyz_yyyyyyy, g_0_x_yyyz_yyyyyyyy, g_0_x_yyyz_yyyyyyyz, g_0_x_yyyz_yyyyyyz, g_0_x_yyyz_yyyyyyzz, g_0_x_yyyz_yyyyyzz, g_0_x_yyyz_yyyyyzzz, g_0_x_yyyz_yyyyzzz, g_0_x_yyyz_yyyyzzzz, g_0_x_yyyz_yyyzzzz, g_0_x_yyyz_yyyzzzzz, g_0_x_yyyz_yyzzzzz, g_0_x_yyyz_yyzzzzzz, g_0_x_yyyz_yzzzzzz, g_0_x_yyyz_yzzzzzzz, g_0_x_yyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyz_xxxxxxx[k] = -g_0_x_yyyz_xxxxxxx[k] * ab_y + g_0_x_yyyz_xxxxxxxy[k];

                g_0_x_yyyyz_xxxxxxy[k] = -g_0_x_yyyz_xxxxxxy[k] * ab_y + g_0_x_yyyz_xxxxxxyy[k];

                g_0_x_yyyyz_xxxxxxz[k] = -g_0_x_yyyz_xxxxxxz[k] * ab_y + g_0_x_yyyz_xxxxxxyz[k];

                g_0_x_yyyyz_xxxxxyy[k] = -g_0_x_yyyz_xxxxxyy[k] * ab_y + g_0_x_yyyz_xxxxxyyy[k];

                g_0_x_yyyyz_xxxxxyz[k] = -g_0_x_yyyz_xxxxxyz[k] * ab_y + g_0_x_yyyz_xxxxxyyz[k];

                g_0_x_yyyyz_xxxxxzz[k] = -g_0_x_yyyz_xxxxxzz[k] * ab_y + g_0_x_yyyz_xxxxxyzz[k];

                g_0_x_yyyyz_xxxxyyy[k] = -g_0_x_yyyz_xxxxyyy[k] * ab_y + g_0_x_yyyz_xxxxyyyy[k];

                g_0_x_yyyyz_xxxxyyz[k] = -g_0_x_yyyz_xxxxyyz[k] * ab_y + g_0_x_yyyz_xxxxyyyz[k];

                g_0_x_yyyyz_xxxxyzz[k] = -g_0_x_yyyz_xxxxyzz[k] * ab_y + g_0_x_yyyz_xxxxyyzz[k];

                g_0_x_yyyyz_xxxxzzz[k] = -g_0_x_yyyz_xxxxzzz[k] * ab_y + g_0_x_yyyz_xxxxyzzz[k];

                g_0_x_yyyyz_xxxyyyy[k] = -g_0_x_yyyz_xxxyyyy[k] * ab_y + g_0_x_yyyz_xxxyyyyy[k];

                g_0_x_yyyyz_xxxyyyz[k] = -g_0_x_yyyz_xxxyyyz[k] * ab_y + g_0_x_yyyz_xxxyyyyz[k];

                g_0_x_yyyyz_xxxyyzz[k] = -g_0_x_yyyz_xxxyyzz[k] * ab_y + g_0_x_yyyz_xxxyyyzz[k];

                g_0_x_yyyyz_xxxyzzz[k] = -g_0_x_yyyz_xxxyzzz[k] * ab_y + g_0_x_yyyz_xxxyyzzz[k];

                g_0_x_yyyyz_xxxzzzz[k] = -g_0_x_yyyz_xxxzzzz[k] * ab_y + g_0_x_yyyz_xxxyzzzz[k];

                g_0_x_yyyyz_xxyyyyy[k] = -g_0_x_yyyz_xxyyyyy[k] * ab_y + g_0_x_yyyz_xxyyyyyy[k];

                g_0_x_yyyyz_xxyyyyz[k] = -g_0_x_yyyz_xxyyyyz[k] * ab_y + g_0_x_yyyz_xxyyyyyz[k];

                g_0_x_yyyyz_xxyyyzz[k] = -g_0_x_yyyz_xxyyyzz[k] * ab_y + g_0_x_yyyz_xxyyyyzz[k];

                g_0_x_yyyyz_xxyyzzz[k] = -g_0_x_yyyz_xxyyzzz[k] * ab_y + g_0_x_yyyz_xxyyyzzz[k];

                g_0_x_yyyyz_xxyzzzz[k] = -g_0_x_yyyz_xxyzzzz[k] * ab_y + g_0_x_yyyz_xxyyzzzz[k];

                g_0_x_yyyyz_xxzzzzz[k] = -g_0_x_yyyz_xxzzzzz[k] * ab_y + g_0_x_yyyz_xxyzzzzz[k];

                g_0_x_yyyyz_xyyyyyy[k] = -g_0_x_yyyz_xyyyyyy[k] * ab_y + g_0_x_yyyz_xyyyyyyy[k];

                g_0_x_yyyyz_xyyyyyz[k] = -g_0_x_yyyz_xyyyyyz[k] * ab_y + g_0_x_yyyz_xyyyyyyz[k];

                g_0_x_yyyyz_xyyyyzz[k] = -g_0_x_yyyz_xyyyyzz[k] * ab_y + g_0_x_yyyz_xyyyyyzz[k];

                g_0_x_yyyyz_xyyyzzz[k] = -g_0_x_yyyz_xyyyzzz[k] * ab_y + g_0_x_yyyz_xyyyyzzz[k];

                g_0_x_yyyyz_xyyzzzz[k] = -g_0_x_yyyz_xyyzzzz[k] * ab_y + g_0_x_yyyz_xyyyzzzz[k];

                g_0_x_yyyyz_xyzzzzz[k] = -g_0_x_yyyz_xyzzzzz[k] * ab_y + g_0_x_yyyz_xyyzzzzz[k];

                g_0_x_yyyyz_xzzzzzz[k] = -g_0_x_yyyz_xzzzzzz[k] * ab_y + g_0_x_yyyz_xyzzzzzz[k];

                g_0_x_yyyyz_yyyyyyy[k] = -g_0_x_yyyz_yyyyyyy[k] * ab_y + g_0_x_yyyz_yyyyyyyy[k];

                g_0_x_yyyyz_yyyyyyz[k] = -g_0_x_yyyz_yyyyyyz[k] * ab_y + g_0_x_yyyz_yyyyyyyz[k];

                g_0_x_yyyyz_yyyyyzz[k] = -g_0_x_yyyz_yyyyyzz[k] * ab_y + g_0_x_yyyz_yyyyyyzz[k];

                g_0_x_yyyyz_yyyyzzz[k] = -g_0_x_yyyz_yyyyzzz[k] * ab_y + g_0_x_yyyz_yyyyyzzz[k];

                g_0_x_yyyyz_yyyzzzz[k] = -g_0_x_yyyz_yyyzzzz[k] * ab_y + g_0_x_yyyz_yyyyzzzz[k];

                g_0_x_yyyyz_yyzzzzz[k] = -g_0_x_yyyz_yyzzzzz[k] * ab_y + g_0_x_yyyz_yyyzzzzz[k];

                g_0_x_yyyyz_yzzzzzz[k] = -g_0_x_yyyz_yzzzzzz[k] * ab_y + g_0_x_yyyz_yyzzzzzz[k];

                g_0_x_yyyyz_zzzzzzz[k] = -g_0_x_yyyz_zzzzzzz[k] * ab_y + g_0_x_yyyz_yzzzzzzz[k];
            }

            /// Set up 612-648 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_x_yyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_x_yyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_x_yyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzz_xxxxxxx, g_0_x_yyyzz_xxxxxxy, g_0_x_yyyzz_xxxxxxz, g_0_x_yyyzz_xxxxxyy, g_0_x_yyyzz_xxxxxyz, g_0_x_yyyzz_xxxxxzz, g_0_x_yyyzz_xxxxyyy, g_0_x_yyyzz_xxxxyyz, g_0_x_yyyzz_xxxxyzz, g_0_x_yyyzz_xxxxzzz, g_0_x_yyyzz_xxxyyyy, g_0_x_yyyzz_xxxyyyz, g_0_x_yyyzz_xxxyyzz, g_0_x_yyyzz_xxxyzzz, g_0_x_yyyzz_xxxzzzz, g_0_x_yyyzz_xxyyyyy, g_0_x_yyyzz_xxyyyyz, g_0_x_yyyzz_xxyyyzz, g_0_x_yyyzz_xxyyzzz, g_0_x_yyyzz_xxyzzzz, g_0_x_yyyzz_xxzzzzz, g_0_x_yyyzz_xyyyyyy, g_0_x_yyyzz_xyyyyyz, g_0_x_yyyzz_xyyyyzz, g_0_x_yyyzz_xyyyzzz, g_0_x_yyyzz_xyyzzzz, g_0_x_yyyzz_xyzzzzz, g_0_x_yyyzz_xzzzzzz, g_0_x_yyyzz_yyyyyyy, g_0_x_yyyzz_yyyyyyz, g_0_x_yyyzz_yyyyyzz, g_0_x_yyyzz_yyyyzzz, g_0_x_yyyzz_yyyzzzz, g_0_x_yyyzz_yyzzzzz, g_0_x_yyyzz_yzzzzzz, g_0_x_yyyzz_zzzzzzz, g_0_x_yyzz_xxxxxxx, g_0_x_yyzz_xxxxxxxy, g_0_x_yyzz_xxxxxxy, g_0_x_yyzz_xxxxxxyy, g_0_x_yyzz_xxxxxxyz, g_0_x_yyzz_xxxxxxz, g_0_x_yyzz_xxxxxyy, g_0_x_yyzz_xxxxxyyy, g_0_x_yyzz_xxxxxyyz, g_0_x_yyzz_xxxxxyz, g_0_x_yyzz_xxxxxyzz, g_0_x_yyzz_xxxxxzz, g_0_x_yyzz_xxxxyyy, g_0_x_yyzz_xxxxyyyy, g_0_x_yyzz_xxxxyyyz, g_0_x_yyzz_xxxxyyz, g_0_x_yyzz_xxxxyyzz, g_0_x_yyzz_xxxxyzz, g_0_x_yyzz_xxxxyzzz, g_0_x_yyzz_xxxxzzz, g_0_x_yyzz_xxxyyyy, g_0_x_yyzz_xxxyyyyy, g_0_x_yyzz_xxxyyyyz, g_0_x_yyzz_xxxyyyz, g_0_x_yyzz_xxxyyyzz, g_0_x_yyzz_xxxyyzz, g_0_x_yyzz_xxxyyzzz, g_0_x_yyzz_xxxyzzz, g_0_x_yyzz_xxxyzzzz, g_0_x_yyzz_xxxzzzz, g_0_x_yyzz_xxyyyyy, g_0_x_yyzz_xxyyyyyy, g_0_x_yyzz_xxyyyyyz, g_0_x_yyzz_xxyyyyz, g_0_x_yyzz_xxyyyyzz, g_0_x_yyzz_xxyyyzz, g_0_x_yyzz_xxyyyzzz, g_0_x_yyzz_xxyyzzz, g_0_x_yyzz_xxyyzzzz, g_0_x_yyzz_xxyzzzz, g_0_x_yyzz_xxyzzzzz, g_0_x_yyzz_xxzzzzz, g_0_x_yyzz_xyyyyyy, g_0_x_yyzz_xyyyyyyy, g_0_x_yyzz_xyyyyyyz, g_0_x_yyzz_xyyyyyz, g_0_x_yyzz_xyyyyyzz, g_0_x_yyzz_xyyyyzz, g_0_x_yyzz_xyyyyzzz, g_0_x_yyzz_xyyyzzz, g_0_x_yyzz_xyyyzzzz, g_0_x_yyzz_xyyzzzz, g_0_x_yyzz_xyyzzzzz, g_0_x_yyzz_xyzzzzz, g_0_x_yyzz_xyzzzzzz, g_0_x_yyzz_xzzzzzz, g_0_x_yyzz_yyyyyyy, g_0_x_yyzz_yyyyyyyy, g_0_x_yyzz_yyyyyyyz, g_0_x_yyzz_yyyyyyz, g_0_x_yyzz_yyyyyyzz, g_0_x_yyzz_yyyyyzz, g_0_x_yyzz_yyyyyzzz, g_0_x_yyzz_yyyyzzz, g_0_x_yyzz_yyyyzzzz, g_0_x_yyzz_yyyzzzz, g_0_x_yyzz_yyyzzzzz, g_0_x_yyzz_yyzzzzz, g_0_x_yyzz_yyzzzzzz, g_0_x_yyzz_yzzzzzz, g_0_x_yyzz_yzzzzzzz, g_0_x_yyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzz_xxxxxxx[k] = -g_0_x_yyzz_xxxxxxx[k] * ab_y + g_0_x_yyzz_xxxxxxxy[k];

                g_0_x_yyyzz_xxxxxxy[k] = -g_0_x_yyzz_xxxxxxy[k] * ab_y + g_0_x_yyzz_xxxxxxyy[k];

                g_0_x_yyyzz_xxxxxxz[k] = -g_0_x_yyzz_xxxxxxz[k] * ab_y + g_0_x_yyzz_xxxxxxyz[k];

                g_0_x_yyyzz_xxxxxyy[k] = -g_0_x_yyzz_xxxxxyy[k] * ab_y + g_0_x_yyzz_xxxxxyyy[k];

                g_0_x_yyyzz_xxxxxyz[k] = -g_0_x_yyzz_xxxxxyz[k] * ab_y + g_0_x_yyzz_xxxxxyyz[k];

                g_0_x_yyyzz_xxxxxzz[k] = -g_0_x_yyzz_xxxxxzz[k] * ab_y + g_0_x_yyzz_xxxxxyzz[k];

                g_0_x_yyyzz_xxxxyyy[k] = -g_0_x_yyzz_xxxxyyy[k] * ab_y + g_0_x_yyzz_xxxxyyyy[k];

                g_0_x_yyyzz_xxxxyyz[k] = -g_0_x_yyzz_xxxxyyz[k] * ab_y + g_0_x_yyzz_xxxxyyyz[k];

                g_0_x_yyyzz_xxxxyzz[k] = -g_0_x_yyzz_xxxxyzz[k] * ab_y + g_0_x_yyzz_xxxxyyzz[k];

                g_0_x_yyyzz_xxxxzzz[k] = -g_0_x_yyzz_xxxxzzz[k] * ab_y + g_0_x_yyzz_xxxxyzzz[k];

                g_0_x_yyyzz_xxxyyyy[k] = -g_0_x_yyzz_xxxyyyy[k] * ab_y + g_0_x_yyzz_xxxyyyyy[k];

                g_0_x_yyyzz_xxxyyyz[k] = -g_0_x_yyzz_xxxyyyz[k] * ab_y + g_0_x_yyzz_xxxyyyyz[k];

                g_0_x_yyyzz_xxxyyzz[k] = -g_0_x_yyzz_xxxyyzz[k] * ab_y + g_0_x_yyzz_xxxyyyzz[k];

                g_0_x_yyyzz_xxxyzzz[k] = -g_0_x_yyzz_xxxyzzz[k] * ab_y + g_0_x_yyzz_xxxyyzzz[k];

                g_0_x_yyyzz_xxxzzzz[k] = -g_0_x_yyzz_xxxzzzz[k] * ab_y + g_0_x_yyzz_xxxyzzzz[k];

                g_0_x_yyyzz_xxyyyyy[k] = -g_0_x_yyzz_xxyyyyy[k] * ab_y + g_0_x_yyzz_xxyyyyyy[k];

                g_0_x_yyyzz_xxyyyyz[k] = -g_0_x_yyzz_xxyyyyz[k] * ab_y + g_0_x_yyzz_xxyyyyyz[k];

                g_0_x_yyyzz_xxyyyzz[k] = -g_0_x_yyzz_xxyyyzz[k] * ab_y + g_0_x_yyzz_xxyyyyzz[k];

                g_0_x_yyyzz_xxyyzzz[k] = -g_0_x_yyzz_xxyyzzz[k] * ab_y + g_0_x_yyzz_xxyyyzzz[k];

                g_0_x_yyyzz_xxyzzzz[k] = -g_0_x_yyzz_xxyzzzz[k] * ab_y + g_0_x_yyzz_xxyyzzzz[k];

                g_0_x_yyyzz_xxzzzzz[k] = -g_0_x_yyzz_xxzzzzz[k] * ab_y + g_0_x_yyzz_xxyzzzzz[k];

                g_0_x_yyyzz_xyyyyyy[k] = -g_0_x_yyzz_xyyyyyy[k] * ab_y + g_0_x_yyzz_xyyyyyyy[k];

                g_0_x_yyyzz_xyyyyyz[k] = -g_0_x_yyzz_xyyyyyz[k] * ab_y + g_0_x_yyzz_xyyyyyyz[k];

                g_0_x_yyyzz_xyyyyzz[k] = -g_0_x_yyzz_xyyyyzz[k] * ab_y + g_0_x_yyzz_xyyyyyzz[k];

                g_0_x_yyyzz_xyyyzzz[k] = -g_0_x_yyzz_xyyyzzz[k] * ab_y + g_0_x_yyzz_xyyyyzzz[k];

                g_0_x_yyyzz_xyyzzzz[k] = -g_0_x_yyzz_xyyzzzz[k] * ab_y + g_0_x_yyzz_xyyyzzzz[k];

                g_0_x_yyyzz_xyzzzzz[k] = -g_0_x_yyzz_xyzzzzz[k] * ab_y + g_0_x_yyzz_xyyzzzzz[k];

                g_0_x_yyyzz_xzzzzzz[k] = -g_0_x_yyzz_xzzzzzz[k] * ab_y + g_0_x_yyzz_xyzzzzzz[k];

                g_0_x_yyyzz_yyyyyyy[k] = -g_0_x_yyzz_yyyyyyy[k] * ab_y + g_0_x_yyzz_yyyyyyyy[k];

                g_0_x_yyyzz_yyyyyyz[k] = -g_0_x_yyzz_yyyyyyz[k] * ab_y + g_0_x_yyzz_yyyyyyyz[k];

                g_0_x_yyyzz_yyyyyzz[k] = -g_0_x_yyzz_yyyyyzz[k] * ab_y + g_0_x_yyzz_yyyyyyzz[k];

                g_0_x_yyyzz_yyyyzzz[k] = -g_0_x_yyzz_yyyyzzz[k] * ab_y + g_0_x_yyzz_yyyyyzzz[k];

                g_0_x_yyyzz_yyyzzzz[k] = -g_0_x_yyzz_yyyzzzz[k] * ab_y + g_0_x_yyzz_yyyyzzzz[k];

                g_0_x_yyyzz_yyzzzzz[k] = -g_0_x_yyzz_yyzzzzz[k] * ab_y + g_0_x_yyzz_yyyzzzzz[k];

                g_0_x_yyyzz_yzzzzzz[k] = -g_0_x_yyzz_yzzzzzz[k] * ab_y + g_0_x_yyzz_yyzzzzzz[k];

                g_0_x_yyyzz_zzzzzzz[k] = -g_0_x_yyzz_zzzzzzz[k] * ab_y + g_0_x_yyzz_yzzzzzzz[k];
            }

            /// Set up 648-684 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_x_yyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_x_yyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_x_yyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzz_xxxxxxx, g_0_x_yyzzz_xxxxxxy, g_0_x_yyzzz_xxxxxxz, g_0_x_yyzzz_xxxxxyy, g_0_x_yyzzz_xxxxxyz, g_0_x_yyzzz_xxxxxzz, g_0_x_yyzzz_xxxxyyy, g_0_x_yyzzz_xxxxyyz, g_0_x_yyzzz_xxxxyzz, g_0_x_yyzzz_xxxxzzz, g_0_x_yyzzz_xxxyyyy, g_0_x_yyzzz_xxxyyyz, g_0_x_yyzzz_xxxyyzz, g_0_x_yyzzz_xxxyzzz, g_0_x_yyzzz_xxxzzzz, g_0_x_yyzzz_xxyyyyy, g_0_x_yyzzz_xxyyyyz, g_0_x_yyzzz_xxyyyzz, g_0_x_yyzzz_xxyyzzz, g_0_x_yyzzz_xxyzzzz, g_0_x_yyzzz_xxzzzzz, g_0_x_yyzzz_xyyyyyy, g_0_x_yyzzz_xyyyyyz, g_0_x_yyzzz_xyyyyzz, g_0_x_yyzzz_xyyyzzz, g_0_x_yyzzz_xyyzzzz, g_0_x_yyzzz_xyzzzzz, g_0_x_yyzzz_xzzzzzz, g_0_x_yyzzz_yyyyyyy, g_0_x_yyzzz_yyyyyyz, g_0_x_yyzzz_yyyyyzz, g_0_x_yyzzz_yyyyzzz, g_0_x_yyzzz_yyyzzzz, g_0_x_yyzzz_yyzzzzz, g_0_x_yyzzz_yzzzzzz, g_0_x_yyzzz_zzzzzzz, g_0_x_yzzz_xxxxxxx, g_0_x_yzzz_xxxxxxxy, g_0_x_yzzz_xxxxxxy, g_0_x_yzzz_xxxxxxyy, g_0_x_yzzz_xxxxxxyz, g_0_x_yzzz_xxxxxxz, g_0_x_yzzz_xxxxxyy, g_0_x_yzzz_xxxxxyyy, g_0_x_yzzz_xxxxxyyz, g_0_x_yzzz_xxxxxyz, g_0_x_yzzz_xxxxxyzz, g_0_x_yzzz_xxxxxzz, g_0_x_yzzz_xxxxyyy, g_0_x_yzzz_xxxxyyyy, g_0_x_yzzz_xxxxyyyz, g_0_x_yzzz_xxxxyyz, g_0_x_yzzz_xxxxyyzz, g_0_x_yzzz_xxxxyzz, g_0_x_yzzz_xxxxyzzz, g_0_x_yzzz_xxxxzzz, g_0_x_yzzz_xxxyyyy, g_0_x_yzzz_xxxyyyyy, g_0_x_yzzz_xxxyyyyz, g_0_x_yzzz_xxxyyyz, g_0_x_yzzz_xxxyyyzz, g_0_x_yzzz_xxxyyzz, g_0_x_yzzz_xxxyyzzz, g_0_x_yzzz_xxxyzzz, g_0_x_yzzz_xxxyzzzz, g_0_x_yzzz_xxxzzzz, g_0_x_yzzz_xxyyyyy, g_0_x_yzzz_xxyyyyyy, g_0_x_yzzz_xxyyyyyz, g_0_x_yzzz_xxyyyyz, g_0_x_yzzz_xxyyyyzz, g_0_x_yzzz_xxyyyzz, g_0_x_yzzz_xxyyyzzz, g_0_x_yzzz_xxyyzzz, g_0_x_yzzz_xxyyzzzz, g_0_x_yzzz_xxyzzzz, g_0_x_yzzz_xxyzzzzz, g_0_x_yzzz_xxzzzzz, g_0_x_yzzz_xyyyyyy, g_0_x_yzzz_xyyyyyyy, g_0_x_yzzz_xyyyyyyz, g_0_x_yzzz_xyyyyyz, g_0_x_yzzz_xyyyyyzz, g_0_x_yzzz_xyyyyzz, g_0_x_yzzz_xyyyyzzz, g_0_x_yzzz_xyyyzzz, g_0_x_yzzz_xyyyzzzz, g_0_x_yzzz_xyyzzzz, g_0_x_yzzz_xyyzzzzz, g_0_x_yzzz_xyzzzzz, g_0_x_yzzz_xyzzzzzz, g_0_x_yzzz_xzzzzzz, g_0_x_yzzz_yyyyyyy, g_0_x_yzzz_yyyyyyyy, g_0_x_yzzz_yyyyyyyz, g_0_x_yzzz_yyyyyyz, g_0_x_yzzz_yyyyyyzz, g_0_x_yzzz_yyyyyzz, g_0_x_yzzz_yyyyyzzz, g_0_x_yzzz_yyyyzzz, g_0_x_yzzz_yyyyzzzz, g_0_x_yzzz_yyyzzzz, g_0_x_yzzz_yyyzzzzz, g_0_x_yzzz_yyzzzzz, g_0_x_yzzz_yyzzzzzz, g_0_x_yzzz_yzzzzzz, g_0_x_yzzz_yzzzzzzz, g_0_x_yzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzz_xxxxxxx[k] = -g_0_x_yzzz_xxxxxxx[k] * ab_y + g_0_x_yzzz_xxxxxxxy[k];

                g_0_x_yyzzz_xxxxxxy[k] = -g_0_x_yzzz_xxxxxxy[k] * ab_y + g_0_x_yzzz_xxxxxxyy[k];

                g_0_x_yyzzz_xxxxxxz[k] = -g_0_x_yzzz_xxxxxxz[k] * ab_y + g_0_x_yzzz_xxxxxxyz[k];

                g_0_x_yyzzz_xxxxxyy[k] = -g_0_x_yzzz_xxxxxyy[k] * ab_y + g_0_x_yzzz_xxxxxyyy[k];

                g_0_x_yyzzz_xxxxxyz[k] = -g_0_x_yzzz_xxxxxyz[k] * ab_y + g_0_x_yzzz_xxxxxyyz[k];

                g_0_x_yyzzz_xxxxxzz[k] = -g_0_x_yzzz_xxxxxzz[k] * ab_y + g_0_x_yzzz_xxxxxyzz[k];

                g_0_x_yyzzz_xxxxyyy[k] = -g_0_x_yzzz_xxxxyyy[k] * ab_y + g_0_x_yzzz_xxxxyyyy[k];

                g_0_x_yyzzz_xxxxyyz[k] = -g_0_x_yzzz_xxxxyyz[k] * ab_y + g_0_x_yzzz_xxxxyyyz[k];

                g_0_x_yyzzz_xxxxyzz[k] = -g_0_x_yzzz_xxxxyzz[k] * ab_y + g_0_x_yzzz_xxxxyyzz[k];

                g_0_x_yyzzz_xxxxzzz[k] = -g_0_x_yzzz_xxxxzzz[k] * ab_y + g_0_x_yzzz_xxxxyzzz[k];

                g_0_x_yyzzz_xxxyyyy[k] = -g_0_x_yzzz_xxxyyyy[k] * ab_y + g_0_x_yzzz_xxxyyyyy[k];

                g_0_x_yyzzz_xxxyyyz[k] = -g_0_x_yzzz_xxxyyyz[k] * ab_y + g_0_x_yzzz_xxxyyyyz[k];

                g_0_x_yyzzz_xxxyyzz[k] = -g_0_x_yzzz_xxxyyzz[k] * ab_y + g_0_x_yzzz_xxxyyyzz[k];

                g_0_x_yyzzz_xxxyzzz[k] = -g_0_x_yzzz_xxxyzzz[k] * ab_y + g_0_x_yzzz_xxxyyzzz[k];

                g_0_x_yyzzz_xxxzzzz[k] = -g_0_x_yzzz_xxxzzzz[k] * ab_y + g_0_x_yzzz_xxxyzzzz[k];

                g_0_x_yyzzz_xxyyyyy[k] = -g_0_x_yzzz_xxyyyyy[k] * ab_y + g_0_x_yzzz_xxyyyyyy[k];

                g_0_x_yyzzz_xxyyyyz[k] = -g_0_x_yzzz_xxyyyyz[k] * ab_y + g_0_x_yzzz_xxyyyyyz[k];

                g_0_x_yyzzz_xxyyyzz[k] = -g_0_x_yzzz_xxyyyzz[k] * ab_y + g_0_x_yzzz_xxyyyyzz[k];

                g_0_x_yyzzz_xxyyzzz[k] = -g_0_x_yzzz_xxyyzzz[k] * ab_y + g_0_x_yzzz_xxyyyzzz[k];

                g_0_x_yyzzz_xxyzzzz[k] = -g_0_x_yzzz_xxyzzzz[k] * ab_y + g_0_x_yzzz_xxyyzzzz[k];

                g_0_x_yyzzz_xxzzzzz[k] = -g_0_x_yzzz_xxzzzzz[k] * ab_y + g_0_x_yzzz_xxyzzzzz[k];

                g_0_x_yyzzz_xyyyyyy[k] = -g_0_x_yzzz_xyyyyyy[k] * ab_y + g_0_x_yzzz_xyyyyyyy[k];

                g_0_x_yyzzz_xyyyyyz[k] = -g_0_x_yzzz_xyyyyyz[k] * ab_y + g_0_x_yzzz_xyyyyyyz[k];

                g_0_x_yyzzz_xyyyyzz[k] = -g_0_x_yzzz_xyyyyzz[k] * ab_y + g_0_x_yzzz_xyyyyyzz[k];

                g_0_x_yyzzz_xyyyzzz[k] = -g_0_x_yzzz_xyyyzzz[k] * ab_y + g_0_x_yzzz_xyyyyzzz[k];

                g_0_x_yyzzz_xyyzzzz[k] = -g_0_x_yzzz_xyyzzzz[k] * ab_y + g_0_x_yzzz_xyyyzzzz[k];

                g_0_x_yyzzz_xyzzzzz[k] = -g_0_x_yzzz_xyzzzzz[k] * ab_y + g_0_x_yzzz_xyyzzzzz[k];

                g_0_x_yyzzz_xzzzzzz[k] = -g_0_x_yzzz_xzzzzzz[k] * ab_y + g_0_x_yzzz_xyzzzzzz[k];

                g_0_x_yyzzz_yyyyyyy[k] = -g_0_x_yzzz_yyyyyyy[k] * ab_y + g_0_x_yzzz_yyyyyyyy[k];

                g_0_x_yyzzz_yyyyyyz[k] = -g_0_x_yzzz_yyyyyyz[k] * ab_y + g_0_x_yzzz_yyyyyyyz[k];

                g_0_x_yyzzz_yyyyyzz[k] = -g_0_x_yzzz_yyyyyzz[k] * ab_y + g_0_x_yzzz_yyyyyyzz[k];

                g_0_x_yyzzz_yyyyzzz[k] = -g_0_x_yzzz_yyyyzzz[k] * ab_y + g_0_x_yzzz_yyyyyzzz[k];

                g_0_x_yyzzz_yyyzzzz[k] = -g_0_x_yzzz_yyyzzzz[k] * ab_y + g_0_x_yzzz_yyyyzzzz[k];

                g_0_x_yyzzz_yyzzzzz[k] = -g_0_x_yzzz_yyzzzzz[k] * ab_y + g_0_x_yzzz_yyyzzzzz[k];

                g_0_x_yyzzz_yzzzzzz[k] = -g_0_x_yzzz_yzzzzzz[k] * ab_y + g_0_x_yzzz_yyzzzzzz[k];

                g_0_x_yyzzz_zzzzzzz[k] = -g_0_x_yzzz_zzzzzzz[k] * ab_y + g_0_x_yzzz_yzzzzzzz[k];
            }

            /// Set up 684-720 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_x_yzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_x_yzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_x_yzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzz_xxxxxxx, g_0_x_yzzzz_xxxxxxy, g_0_x_yzzzz_xxxxxxz, g_0_x_yzzzz_xxxxxyy, g_0_x_yzzzz_xxxxxyz, g_0_x_yzzzz_xxxxxzz, g_0_x_yzzzz_xxxxyyy, g_0_x_yzzzz_xxxxyyz, g_0_x_yzzzz_xxxxyzz, g_0_x_yzzzz_xxxxzzz, g_0_x_yzzzz_xxxyyyy, g_0_x_yzzzz_xxxyyyz, g_0_x_yzzzz_xxxyyzz, g_0_x_yzzzz_xxxyzzz, g_0_x_yzzzz_xxxzzzz, g_0_x_yzzzz_xxyyyyy, g_0_x_yzzzz_xxyyyyz, g_0_x_yzzzz_xxyyyzz, g_0_x_yzzzz_xxyyzzz, g_0_x_yzzzz_xxyzzzz, g_0_x_yzzzz_xxzzzzz, g_0_x_yzzzz_xyyyyyy, g_0_x_yzzzz_xyyyyyz, g_0_x_yzzzz_xyyyyzz, g_0_x_yzzzz_xyyyzzz, g_0_x_yzzzz_xyyzzzz, g_0_x_yzzzz_xyzzzzz, g_0_x_yzzzz_xzzzzzz, g_0_x_yzzzz_yyyyyyy, g_0_x_yzzzz_yyyyyyz, g_0_x_yzzzz_yyyyyzz, g_0_x_yzzzz_yyyyzzz, g_0_x_yzzzz_yyyzzzz, g_0_x_yzzzz_yyzzzzz, g_0_x_yzzzz_yzzzzzz, g_0_x_yzzzz_zzzzzzz, g_0_x_zzzz_xxxxxxx, g_0_x_zzzz_xxxxxxxy, g_0_x_zzzz_xxxxxxy, g_0_x_zzzz_xxxxxxyy, g_0_x_zzzz_xxxxxxyz, g_0_x_zzzz_xxxxxxz, g_0_x_zzzz_xxxxxyy, g_0_x_zzzz_xxxxxyyy, g_0_x_zzzz_xxxxxyyz, g_0_x_zzzz_xxxxxyz, g_0_x_zzzz_xxxxxyzz, g_0_x_zzzz_xxxxxzz, g_0_x_zzzz_xxxxyyy, g_0_x_zzzz_xxxxyyyy, g_0_x_zzzz_xxxxyyyz, g_0_x_zzzz_xxxxyyz, g_0_x_zzzz_xxxxyyzz, g_0_x_zzzz_xxxxyzz, g_0_x_zzzz_xxxxyzzz, g_0_x_zzzz_xxxxzzz, g_0_x_zzzz_xxxyyyy, g_0_x_zzzz_xxxyyyyy, g_0_x_zzzz_xxxyyyyz, g_0_x_zzzz_xxxyyyz, g_0_x_zzzz_xxxyyyzz, g_0_x_zzzz_xxxyyzz, g_0_x_zzzz_xxxyyzzz, g_0_x_zzzz_xxxyzzz, g_0_x_zzzz_xxxyzzzz, g_0_x_zzzz_xxxzzzz, g_0_x_zzzz_xxyyyyy, g_0_x_zzzz_xxyyyyyy, g_0_x_zzzz_xxyyyyyz, g_0_x_zzzz_xxyyyyz, g_0_x_zzzz_xxyyyyzz, g_0_x_zzzz_xxyyyzz, g_0_x_zzzz_xxyyyzzz, g_0_x_zzzz_xxyyzzz, g_0_x_zzzz_xxyyzzzz, g_0_x_zzzz_xxyzzzz, g_0_x_zzzz_xxyzzzzz, g_0_x_zzzz_xxzzzzz, g_0_x_zzzz_xyyyyyy, g_0_x_zzzz_xyyyyyyy, g_0_x_zzzz_xyyyyyyz, g_0_x_zzzz_xyyyyyz, g_0_x_zzzz_xyyyyyzz, g_0_x_zzzz_xyyyyzz, g_0_x_zzzz_xyyyyzzz, g_0_x_zzzz_xyyyzzz, g_0_x_zzzz_xyyyzzzz, g_0_x_zzzz_xyyzzzz, g_0_x_zzzz_xyyzzzzz, g_0_x_zzzz_xyzzzzz, g_0_x_zzzz_xyzzzzzz, g_0_x_zzzz_xzzzzzz, g_0_x_zzzz_yyyyyyy, g_0_x_zzzz_yyyyyyyy, g_0_x_zzzz_yyyyyyyz, g_0_x_zzzz_yyyyyyz, g_0_x_zzzz_yyyyyyzz, g_0_x_zzzz_yyyyyzz, g_0_x_zzzz_yyyyyzzz, g_0_x_zzzz_yyyyzzz, g_0_x_zzzz_yyyyzzzz, g_0_x_zzzz_yyyzzzz, g_0_x_zzzz_yyyzzzzz, g_0_x_zzzz_yyzzzzz, g_0_x_zzzz_yyzzzzzz, g_0_x_zzzz_yzzzzzz, g_0_x_zzzz_yzzzzzzz, g_0_x_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzz_xxxxxxx[k] = -g_0_x_zzzz_xxxxxxx[k] * ab_y + g_0_x_zzzz_xxxxxxxy[k];

                g_0_x_yzzzz_xxxxxxy[k] = -g_0_x_zzzz_xxxxxxy[k] * ab_y + g_0_x_zzzz_xxxxxxyy[k];

                g_0_x_yzzzz_xxxxxxz[k] = -g_0_x_zzzz_xxxxxxz[k] * ab_y + g_0_x_zzzz_xxxxxxyz[k];

                g_0_x_yzzzz_xxxxxyy[k] = -g_0_x_zzzz_xxxxxyy[k] * ab_y + g_0_x_zzzz_xxxxxyyy[k];

                g_0_x_yzzzz_xxxxxyz[k] = -g_0_x_zzzz_xxxxxyz[k] * ab_y + g_0_x_zzzz_xxxxxyyz[k];

                g_0_x_yzzzz_xxxxxzz[k] = -g_0_x_zzzz_xxxxxzz[k] * ab_y + g_0_x_zzzz_xxxxxyzz[k];

                g_0_x_yzzzz_xxxxyyy[k] = -g_0_x_zzzz_xxxxyyy[k] * ab_y + g_0_x_zzzz_xxxxyyyy[k];

                g_0_x_yzzzz_xxxxyyz[k] = -g_0_x_zzzz_xxxxyyz[k] * ab_y + g_0_x_zzzz_xxxxyyyz[k];

                g_0_x_yzzzz_xxxxyzz[k] = -g_0_x_zzzz_xxxxyzz[k] * ab_y + g_0_x_zzzz_xxxxyyzz[k];

                g_0_x_yzzzz_xxxxzzz[k] = -g_0_x_zzzz_xxxxzzz[k] * ab_y + g_0_x_zzzz_xxxxyzzz[k];

                g_0_x_yzzzz_xxxyyyy[k] = -g_0_x_zzzz_xxxyyyy[k] * ab_y + g_0_x_zzzz_xxxyyyyy[k];

                g_0_x_yzzzz_xxxyyyz[k] = -g_0_x_zzzz_xxxyyyz[k] * ab_y + g_0_x_zzzz_xxxyyyyz[k];

                g_0_x_yzzzz_xxxyyzz[k] = -g_0_x_zzzz_xxxyyzz[k] * ab_y + g_0_x_zzzz_xxxyyyzz[k];

                g_0_x_yzzzz_xxxyzzz[k] = -g_0_x_zzzz_xxxyzzz[k] * ab_y + g_0_x_zzzz_xxxyyzzz[k];

                g_0_x_yzzzz_xxxzzzz[k] = -g_0_x_zzzz_xxxzzzz[k] * ab_y + g_0_x_zzzz_xxxyzzzz[k];

                g_0_x_yzzzz_xxyyyyy[k] = -g_0_x_zzzz_xxyyyyy[k] * ab_y + g_0_x_zzzz_xxyyyyyy[k];

                g_0_x_yzzzz_xxyyyyz[k] = -g_0_x_zzzz_xxyyyyz[k] * ab_y + g_0_x_zzzz_xxyyyyyz[k];

                g_0_x_yzzzz_xxyyyzz[k] = -g_0_x_zzzz_xxyyyzz[k] * ab_y + g_0_x_zzzz_xxyyyyzz[k];

                g_0_x_yzzzz_xxyyzzz[k] = -g_0_x_zzzz_xxyyzzz[k] * ab_y + g_0_x_zzzz_xxyyyzzz[k];

                g_0_x_yzzzz_xxyzzzz[k] = -g_0_x_zzzz_xxyzzzz[k] * ab_y + g_0_x_zzzz_xxyyzzzz[k];

                g_0_x_yzzzz_xxzzzzz[k] = -g_0_x_zzzz_xxzzzzz[k] * ab_y + g_0_x_zzzz_xxyzzzzz[k];

                g_0_x_yzzzz_xyyyyyy[k] = -g_0_x_zzzz_xyyyyyy[k] * ab_y + g_0_x_zzzz_xyyyyyyy[k];

                g_0_x_yzzzz_xyyyyyz[k] = -g_0_x_zzzz_xyyyyyz[k] * ab_y + g_0_x_zzzz_xyyyyyyz[k];

                g_0_x_yzzzz_xyyyyzz[k] = -g_0_x_zzzz_xyyyyzz[k] * ab_y + g_0_x_zzzz_xyyyyyzz[k];

                g_0_x_yzzzz_xyyyzzz[k] = -g_0_x_zzzz_xyyyzzz[k] * ab_y + g_0_x_zzzz_xyyyyzzz[k];

                g_0_x_yzzzz_xyyzzzz[k] = -g_0_x_zzzz_xyyzzzz[k] * ab_y + g_0_x_zzzz_xyyyzzzz[k];

                g_0_x_yzzzz_xyzzzzz[k] = -g_0_x_zzzz_xyzzzzz[k] * ab_y + g_0_x_zzzz_xyyzzzzz[k];

                g_0_x_yzzzz_xzzzzzz[k] = -g_0_x_zzzz_xzzzzzz[k] * ab_y + g_0_x_zzzz_xyzzzzzz[k];

                g_0_x_yzzzz_yyyyyyy[k] = -g_0_x_zzzz_yyyyyyy[k] * ab_y + g_0_x_zzzz_yyyyyyyy[k];

                g_0_x_yzzzz_yyyyyyz[k] = -g_0_x_zzzz_yyyyyyz[k] * ab_y + g_0_x_zzzz_yyyyyyyz[k];

                g_0_x_yzzzz_yyyyyzz[k] = -g_0_x_zzzz_yyyyyzz[k] * ab_y + g_0_x_zzzz_yyyyyyzz[k];

                g_0_x_yzzzz_yyyyzzz[k] = -g_0_x_zzzz_yyyyzzz[k] * ab_y + g_0_x_zzzz_yyyyyzzz[k];

                g_0_x_yzzzz_yyyzzzz[k] = -g_0_x_zzzz_yyyzzzz[k] * ab_y + g_0_x_zzzz_yyyyzzzz[k];

                g_0_x_yzzzz_yyzzzzz[k] = -g_0_x_zzzz_yyzzzzz[k] * ab_y + g_0_x_zzzz_yyyzzzzz[k];

                g_0_x_yzzzz_yzzzzzz[k] = -g_0_x_zzzz_yzzzzzz[k] * ab_y + g_0_x_zzzz_yyzzzzzz[k];

                g_0_x_yzzzz_zzzzzzz[k] = -g_0_x_zzzz_zzzzzzz[k] * ab_y + g_0_x_zzzz_yzzzzzzz[k];
            }

            /// Set up 720-756 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_x_zzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_x_zzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_x_zzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_xxxxxxx, g_0_x_zzzz_xxxxxxxz, g_0_x_zzzz_xxxxxxy, g_0_x_zzzz_xxxxxxyz, g_0_x_zzzz_xxxxxxz, g_0_x_zzzz_xxxxxxzz, g_0_x_zzzz_xxxxxyy, g_0_x_zzzz_xxxxxyyz, g_0_x_zzzz_xxxxxyz, g_0_x_zzzz_xxxxxyzz, g_0_x_zzzz_xxxxxzz, g_0_x_zzzz_xxxxxzzz, g_0_x_zzzz_xxxxyyy, g_0_x_zzzz_xxxxyyyz, g_0_x_zzzz_xxxxyyz, g_0_x_zzzz_xxxxyyzz, g_0_x_zzzz_xxxxyzz, g_0_x_zzzz_xxxxyzzz, g_0_x_zzzz_xxxxzzz, g_0_x_zzzz_xxxxzzzz, g_0_x_zzzz_xxxyyyy, g_0_x_zzzz_xxxyyyyz, g_0_x_zzzz_xxxyyyz, g_0_x_zzzz_xxxyyyzz, g_0_x_zzzz_xxxyyzz, g_0_x_zzzz_xxxyyzzz, g_0_x_zzzz_xxxyzzz, g_0_x_zzzz_xxxyzzzz, g_0_x_zzzz_xxxzzzz, g_0_x_zzzz_xxxzzzzz, g_0_x_zzzz_xxyyyyy, g_0_x_zzzz_xxyyyyyz, g_0_x_zzzz_xxyyyyz, g_0_x_zzzz_xxyyyyzz, g_0_x_zzzz_xxyyyzz, g_0_x_zzzz_xxyyyzzz, g_0_x_zzzz_xxyyzzz, g_0_x_zzzz_xxyyzzzz, g_0_x_zzzz_xxyzzzz, g_0_x_zzzz_xxyzzzzz, g_0_x_zzzz_xxzzzzz, g_0_x_zzzz_xxzzzzzz, g_0_x_zzzz_xyyyyyy, g_0_x_zzzz_xyyyyyyz, g_0_x_zzzz_xyyyyyz, g_0_x_zzzz_xyyyyyzz, g_0_x_zzzz_xyyyyzz, g_0_x_zzzz_xyyyyzzz, g_0_x_zzzz_xyyyzzz, g_0_x_zzzz_xyyyzzzz, g_0_x_zzzz_xyyzzzz, g_0_x_zzzz_xyyzzzzz, g_0_x_zzzz_xyzzzzz, g_0_x_zzzz_xyzzzzzz, g_0_x_zzzz_xzzzzzz, g_0_x_zzzz_xzzzzzzz, g_0_x_zzzz_yyyyyyy, g_0_x_zzzz_yyyyyyyz, g_0_x_zzzz_yyyyyyz, g_0_x_zzzz_yyyyyyzz, g_0_x_zzzz_yyyyyzz, g_0_x_zzzz_yyyyyzzz, g_0_x_zzzz_yyyyzzz, g_0_x_zzzz_yyyyzzzz, g_0_x_zzzz_yyyzzzz, g_0_x_zzzz_yyyzzzzz, g_0_x_zzzz_yyzzzzz, g_0_x_zzzz_yyzzzzzz, g_0_x_zzzz_yzzzzzz, g_0_x_zzzz_yzzzzzzz, g_0_x_zzzz_zzzzzzz, g_0_x_zzzz_zzzzzzzz, g_0_x_zzzzz_xxxxxxx, g_0_x_zzzzz_xxxxxxy, g_0_x_zzzzz_xxxxxxz, g_0_x_zzzzz_xxxxxyy, g_0_x_zzzzz_xxxxxyz, g_0_x_zzzzz_xxxxxzz, g_0_x_zzzzz_xxxxyyy, g_0_x_zzzzz_xxxxyyz, g_0_x_zzzzz_xxxxyzz, g_0_x_zzzzz_xxxxzzz, g_0_x_zzzzz_xxxyyyy, g_0_x_zzzzz_xxxyyyz, g_0_x_zzzzz_xxxyyzz, g_0_x_zzzzz_xxxyzzz, g_0_x_zzzzz_xxxzzzz, g_0_x_zzzzz_xxyyyyy, g_0_x_zzzzz_xxyyyyz, g_0_x_zzzzz_xxyyyzz, g_0_x_zzzzz_xxyyzzz, g_0_x_zzzzz_xxyzzzz, g_0_x_zzzzz_xxzzzzz, g_0_x_zzzzz_xyyyyyy, g_0_x_zzzzz_xyyyyyz, g_0_x_zzzzz_xyyyyzz, g_0_x_zzzzz_xyyyzzz, g_0_x_zzzzz_xyyzzzz, g_0_x_zzzzz_xyzzzzz, g_0_x_zzzzz_xzzzzzz, g_0_x_zzzzz_yyyyyyy, g_0_x_zzzzz_yyyyyyz, g_0_x_zzzzz_yyyyyzz, g_0_x_zzzzz_yyyyzzz, g_0_x_zzzzz_yyyzzzz, g_0_x_zzzzz_yyzzzzz, g_0_x_zzzzz_yzzzzzz, g_0_x_zzzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzz_xxxxxxx[k] = -g_0_x_zzzz_xxxxxxx[k] * ab_z + g_0_x_zzzz_xxxxxxxz[k];

                g_0_x_zzzzz_xxxxxxy[k] = -g_0_x_zzzz_xxxxxxy[k] * ab_z + g_0_x_zzzz_xxxxxxyz[k];

                g_0_x_zzzzz_xxxxxxz[k] = -g_0_x_zzzz_xxxxxxz[k] * ab_z + g_0_x_zzzz_xxxxxxzz[k];

                g_0_x_zzzzz_xxxxxyy[k] = -g_0_x_zzzz_xxxxxyy[k] * ab_z + g_0_x_zzzz_xxxxxyyz[k];

                g_0_x_zzzzz_xxxxxyz[k] = -g_0_x_zzzz_xxxxxyz[k] * ab_z + g_0_x_zzzz_xxxxxyzz[k];

                g_0_x_zzzzz_xxxxxzz[k] = -g_0_x_zzzz_xxxxxzz[k] * ab_z + g_0_x_zzzz_xxxxxzzz[k];

                g_0_x_zzzzz_xxxxyyy[k] = -g_0_x_zzzz_xxxxyyy[k] * ab_z + g_0_x_zzzz_xxxxyyyz[k];

                g_0_x_zzzzz_xxxxyyz[k] = -g_0_x_zzzz_xxxxyyz[k] * ab_z + g_0_x_zzzz_xxxxyyzz[k];

                g_0_x_zzzzz_xxxxyzz[k] = -g_0_x_zzzz_xxxxyzz[k] * ab_z + g_0_x_zzzz_xxxxyzzz[k];

                g_0_x_zzzzz_xxxxzzz[k] = -g_0_x_zzzz_xxxxzzz[k] * ab_z + g_0_x_zzzz_xxxxzzzz[k];

                g_0_x_zzzzz_xxxyyyy[k] = -g_0_x_zzzz_xxxyyyy[k] * ab_z + g_0_x_zzzz_xxxyyyyz[k];

                g_0_x_zzzzz_xxxyyyz[k] = -g_0_x_zzzz_xxxyyyz[k] * ab_z + g_0_x_zzzz_xxxyyyzz[k];

                g_0_x_zzzzz_xxxyyzz[k] = -g_0_x_zzzz_xxxyyzz[k] * ab_z + g_0_x_zzzz_xxxyyzzz[k];

                g_0_x_zzzzz_xxxyzzz[k] = -g_0_x_zzzz_xxxyzzz[k] * ab_z + g_0_x_zzzz_xxxyzzzz[k];

                g_0_x_zzzzz_xxxzzzz[k] = -g_0_x_zzzz_xxxzzzz[k] * ab_z + g_0_x_zzzz_xxxzzzzz[k];

                g_0_x_zzzzz_xxyyyyy[k] = -g_0_x_zzzz_xxyyyyy[k] * ab_z + g_0_x_zzzz_xxyyyyyz[k];

                g_0_x_zzzzz_xxyyyyz[k] = -g_0_x_zzzz_xxyyyyz[k] * ab_z + g_0_x_zzzz_xxyyyyzz[k];

                g_0_x_zzzzz_xxyyyzz[k] = -g_0_x_zzzz_xxyyyzz[k] * ab_z + g_0_x_zzzz_xxyyyzzz[k];

                g_0_x_zzzzz_xxyyzzz[k] = -g_0_x_zzzz_xxyyzzz[k] * ab_z + g_0_x_zzzz_xxyyzzzz[k];

                g_0_x_zzzzz_xxyzzzz[k] = -g_0_x_zzzz_xxyzzzz[k] * ab_z + g_0_x_zzzz_xxyzzzzz[k];

                g_0_x_zzzzz_xxzzzzz[k] = -g_0_x_zzzz_xxzzzzz[k] * ab_z + g_0_x_zzzz_xxzzzzzz[k];

                g_0_x_zzzzz_xyyyyyy[k] = -g_0_x_zzzz_xyyyyyy[k] * ab_z + g_0_x_zzzz_xyyyyyyz[k];

                g_0_x_zzzzz_xyyyyyz[k] = -g_0_x_zzzz_xyyyyyz[k] * ab_z + g_0_x_zzzz_xyyyyyzz[k];

                g_0_x_zzzzz_xyyyyzz[k] = -g_0_x_zzzz_xyyyyzz[k] * ab_z + g_0_x_zzzz_xyyyyzzz[k];

                g_0_x_zzzzz_xyyyzzz[k] = -g_0_x_zzzz_xyyyzzz[k] * ab_z + g_0_x_zzzz_xyyyzzzz[k];

                g_0_x_zzzzz_xyyzzzz[k] = -g_0_x_zzzz_xyyzzzz[k] * ab_z + g_0_x_zzzz_xyyzzzzz[k];

                g_0_x_zzzzz_xyzzzzz[k] = -g_0_x_zzzz_xyzzzzz[k] * ab_z + g_0_x_zzzz_xyzzzzzz[k];

                g_0_x_zzzzz_xzzzzzz[k] = -g_0_x_zzzz_xzzzzzz[k] * ab_z + g_0_x_zzzz_xzzzzzzz[k];

                g_0_x_zzzzz_yyyyyyy[k] = -g_0_x_zzzz_yyyyyyy[k] * ab_z + g_0_x_zzzz_yyyyyyyz[k];

                g_0_x_zzzzz_yyyyyyz[k] = -g_0_x_zzzz_yyyyyyz[k] * ab_z + g_0_x_zzzz_yyyyyyzz[k];

                g_0_x_zzzzz_yyyyyzz[k] = -g_0_x_zzzz_yyyyyzz[k] * ab_z + g_0_x_zzzz_yyyyyzzz[k];

                g_0_x_zzzzz_yyyyzzz[k] = -g_0_x_zzzz_yyyyzzz[k] * ab_z + g_0_x_zzzz_yyyyzzzz[k];

                g_0_x_zzzzz_yyyzzzz[k] = -g_0_x_zzzz_yyyzzzz[k] * ab_z + g_0_x_zzzz_yyyzzzzz[k];

                g_0_x_zzzzz_yyzzzzz[k] = -g_0_x_zzzz_yyzzzzz[k] * ab_z + g_0_x_zzzz_yyzzzzzz[k];

                g_0_x_zzzzz_yzzzzzz[k] = -g_0_x_zzzz_yzzzzzz[k] * ab_z + g_0_x_zzzz_yzzzzzzz[k];

                g_0_x_zzzzz_zzzzzzz[k] = -g_0_x_zzzz_zzzzzzz[k] * ab_z + g_0_x_zzzz_zzzzzzzz[k];
            }

            /// Set up 756-792 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxx_xxxxxxx = cbuffer.data(hk_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxxxy = cbuffer.data(hk_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxxxz = cbuffer.data(hk_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxxyy = cbuffer.data(hk_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxxyz = cbuffer.data(hk_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxxzz = cbuffer.data(hk_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxyyy = cbuffer.data(hk_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxyyz = cbuffer.data(hk_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxyzz = cbuffer.data(hk_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxzzz = cbuffer.data(hk_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyyyy = cbuffer.data(hk_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyyyz = cbuffer.data(hk_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyyzz = cbuffer.data(hk_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyzzz = cbuffer.data(hk_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxzzzz = cbuffer.data(hk_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyyyy = cbuffer.data(hk_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyyyz = cbuffer.data(hk_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyyzz = cbuffer.data(hk_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyzzz = cbuffer.data(hk_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyzzzz = cbuffer.data(hk_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxzzzzz = cbuffer.data(hk_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyyyy = cbuffer.data(hk_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyyyz = cbuffer.data(hk_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyyzz = cbuffer.data(hk_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyzzz = cbuffer.data(hk_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyzzzz = cbuffer.data(hk_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyzzzzz = cbuffer.data(hk_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_xxxxx_xzzzzzz = cbuffer.data(hk_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyyyy = cbuffer.data(hk_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyyyz = cbuffer.data(hk_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyyzz = cbuffer.data(hk_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyzzz = cbuffer.data(hk_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyzzzz = cbuffer.data(hk_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyzzzzz = cbuffer.data(hk_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_xxxxx_yzzzzzz = cbuffer.data(hk_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_xxxxx_zzzzzzz = cbuffer.data(hk_geom_01_off + 791 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_xxxxxxx, g_0_y_xxxx_xxxxxxxx, g_0_y_xxxx_xxxxxxxy, g_0_y_xxxx_xxxxxxxz, g_0_y_xxxx_xxxxxxy, g_0_y_xxxx_xxxxxxyy, g_0_y_xxxx_xxxxxxyz, g_0_y_xxxx_xxxxxxz, g_0_y_xxxx_xxxxxxzz, g_0_y_xxxx_xxxxxyy, g_0_y_xxxx_xxxxxyyy, g_0_y_xxxx_xxxxxyyz, g_0_y_xxxx_xxxxxyz, g_0_y_xxxx_xxxxxyzz, g_0_y_xxxx_xxxxxzz, g_0_y_xxxx_xxxxxzzz, g_0_y_xxxx_xxxxyyy, g_0_y_xxxx_xxxxyyyy, g_0_y_xxxx_xxxxyyyz, g_0_y_xxxx_xxxxyyz, g_0_y_xxxx_xxxxyyzz, g_0_y_xxxx_xxxxyzz, g_0_y_xxxx_xxxxyzzz, g_0_y_xxxx_xxxxzzz, g_0_y_xxxx_xxxxzzzz, g_0_y_xxxx_xxxyyyy, g_0_y_xxxx_xxxyyyyy, g_0_y_xxxx_xxxyyyyz, g_0_y_xxxx_xxxyyyz, g_0_y_xxxx_xxxyyyzz, g_0_y_xxxx_xxxyyzz, g_0_y_xxxx_xxxyyzzz, g_0_y_xxxx_xxxyzzz, g_0_y_xxxx_xxxyzzzz, g_0_y_xxxx_xxxzzzz, g_0_y_xxxx_xxxzzzzz, g_0_y_xxxx_xxyyyyy, g_0_y_xxxx_xxyyyyyy, g_0_y_xxxx_xxyyyyyz, g_0_y_xxxx_xxyyyyz, g_0_y_xxxx_xxyyyyzz, g_0_y_xxxx_xxyyyzz, g_0_y_xxxx_xxyyyzzz, g_0_y_xxxx_xxyyzzz, g_0_y_xxxx_xxyyzzzz, g_0_y_xxxx_xxyzzzz, g_0_y_xxxx_xxyzzzzz, g_0_y_xxxx_xxzzzzz, g_0_y_xxxx_xxzzzzzz, g_0_y_xxxx_xyyyyyy, g_0_y_xxxx_xyyyyyyy, g_0_y_xxxx_xyyyyyyz, g_0_y_xxxx_xyyyyyz, g_0_y_xxxx_xyyyyyzz, g_0_y_xxxx_xyyyyzz, g_0_y_xxxx_xyyyyzzz, g_0_y_xxxx_xyyyzzz, g_0_y_xxxx_xyyyzzzz, g_0_y_xxxx_xyyzzzz, g_0_y_xxxx_xyyzzzzz, g_0_y_xxxx_xyzzzzz, g_0_y_xxxx_xyzzzzzz, g_0_y_xxxx_xzzzzzz, g_0_y_xxxx_xzzzzzzz, g_0_y_xxxx_yyyyyyy, g_0_y_xxxx_yyyyyyz, g_0_y_xxxx_yyyyyzz, g_0_y_xxxx_yyyyzzz, g_0_y_xxxx_yyyzzzz, g_0_y_xxxx_yyzzzzz, g_0_y_xxxx_yzzzzzz, g_0_y_xxxx_zzzzzzz, g_0_y_xxxxx_xxxxxxx, g_0_y_xxxxx_xxxxxxy, g_0_y_xxxxx_xxxxxxz, g_0_y_xxxxx_xxxxxyy, g_0_y_xxxxx_xxxxxyz, g_0_y_xxxxx_xxxxxzz, g_0_y_xxxxx_xxxxyyy, g_0_y_xxxxx_xxxxyyz, g_0_y_xxxxx_xxxxyzz, g_0_y_xxxxx_xxxxzzz, g_0_y_xxxxx_xxxyyyy, g_0_y_xxxxx_xxxyyyz, g_0_y_xxxxx_xxxyyzz, g_0_y_xxxxx_xxxyzzz, g_0_y_xxxxx_xxxzzzz, g_0_y_xxxxx_xxyyyyy, g_0_y_xxxxx_xxyyyyz, g_0_y_xxxxx_xxyyyzz, g_0_y_xxxxx_xxyyzzz, g_0_y_xxxxx_xxyzzzz, g_0_y_xxxxx_xxzzzzz, g_0_y_xxxxx_xyyyyyy, g_0_y_xxxxx_xyyyyyz, g_0_y_xxxxx_xyyyyzz, g_0_y_xxxxx_xyyyzzz, g_0_y_xxxxx_xyyzzzz, g_0_y_xxxxx_xyzzzzz, g_0_y_xxxxx_xzzzzzz, g_0_y_xxxxx_yyyyyyy, g_0_y_xxxxx_yyyyyyz, g_0_y_xxxxx_yyyyyzz, g_0_y_xxxxx_yyyyzzz, g_0_y_xxxxx_yyyzzzz, g_0_y_xxxxx_yyzzzzz, g_0_y_xxxxx_yzzzzzz, g_0_y_xxxxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxx_xxxxxxx[k] = -g_0_y_xxxx_xxxxxxx[k] * ab_x + g_0_y_xxxx_xxxxxxxx[k];

                g_0_y_xxxxx_xxxxxxy[k] = -g_0_y_xxxx_xxxxxxy[k] * ab_x + g_0_y_xxxx_xxxxxxxy[k];

                g_0_y_xxxxx_xxxxxxz[k] = -g_0_y_xxxx_xxxxxxz[k] * ab_x + g_0_y_xxxx_xxxxxxxz[k];

                g_0_y_xxxxx_xxxxxyy[k] = -g_0_y_xxxx_xxxxxyy[k] * ab_x + g_0_y_xxxx_xxxxxxyy[k];

                g_0_y_xxxxx_xxxxxyz[k] = -g_0_y_xxxx_xxxxxyz[k] * ab_x + g_0_y_xxxx_xxxxxxyz[k];

                g_0_y_xxxxx_xxxxxzz[k] = -g_0_y_xxxx_xxxxxzz[k] * ab_x + g_0_y_xxxx_xxxxxxzz[k];

                g_0_y_xxxxx_xxxxyyy[k] = -g_0_y_xxxx_xxxxyyy[k] * ab_x + g_0_y_xxxx_xxxxxyyy[k];

                g_0_y_xxxxx_xxxxyyz[k] = -g_0_y_xxxx_xxxxyyz[k] * ab_x + g_0_y_xxxx_xxxxxyyz[k];

                g_0_y_xxxxx_xxxxyzz[k] = -g_0_y_xxxx_xxxxyzz[k] * ab_x + g_0_y_xxxx_xxxxxyzz[k];

                g_0_y_xxxxx_xxxxzzz[k] = -g_0_y_xxxx_xxxxzzz[k] * ab_x + g_0_y_xxxx_xxxxxzzz[k];

                g_0_y_xxxxx_xxxyyyy[k] = -g_0_y_xxxx_xxxyyyy[k] * ab_x + g_0_y_xxxx_xxxxyyyy[k];

                g_0_y_xxxxx_xxxyyyz[k] = -g_0_y_xxxx_xxxyyyz[k] * ab_x + g_0_y_xxxx_xxxxyyyz[k];

                g_0_y_xxxxx_xxxyyzz[k] = -g_0_y_xxxx_xxxyyzz[k] * ab_x + g_0_y_xxxx_xxxxyyzz[k];

                g_0_y_xxxxx_xxxyzzz[k] = -g_0_y_xxxx_xxxyzzz[k] * ab_x + g_0_y_xxxx_xxxxyzzz[k];

                g_0_y_xxxxx_xxxzzzz[k] = -g_0_y_xxxx_xxxzzzz[k] * ab_x + g_0_y_xxxx_xxxxzzzz[k];

                g_0_y_xxxxx_xxyyyyy[k] = -g_0_y_xxxx_xxyyyyy[k] * ab_x + g_0_y_xxxx_xxxyyyyy[k];

                g_0_y_xxxxx_xxyyyyz[k] = -g_0_y_xxxx_xxyyyyz[k] * ab_x + g_0_y_xxxx_xxxyyyyz[k];

                g_0_y_xxxxx_xxyyyzz[k] = -g_0_y_xxxx_xxyyyzz[k] * ab_x + g_0_y_xxxx_xxxyyyzz[k];

                g_0_y_xxxxx_xxyyzzz[k] = -g_0_y_xxxx_xxyyzzz[k] * ab_x + g_0_y_xxxx_xxxyyzzz[k];

                g_0_y_xxxxx_xxyzzzz[k] = -g_0_y_xxxx_xxyzzzz[k] * ab_x + g_0_y_xxxx_xxxyzzzz[k];

                g_0_y_xxxxx_xxzzzzz[k] = -g_0_y_xxxx_xxzzzzz[k] * ab_x + g_0_y_xxxx_xxxzzzzz[k];

                g_0_y_xxxxx_xyyyyyy[k] = -g_0_y_xxxx_xyyyyyy[k] * ab_x + g_0_y_xxxx_xxyyyyyy[k];

                g_0_y_xxxxx_xyyyyyz[k] = -g_0_y_xxxx_xyyyyyz[k] * ab_x + g_0_y_xxxx_xxyyyyyz[k];

                g_0_y_xxxxx_xyyyyzz[k] = -g_0_y_xxxx_xyyyyzz[k] * ab_x + g_0_y_xxxx_xxyyyyzz[k];

                g_0_y_xxxxx_xyyyzzz[k] = -g_0_y_xxxx_xyyyzzz[k] * ab_x + g_0_y_xxxx_xxyyyzzz[k];

                g_0_y_xxxxx_xyyzzzz[k] = -g_0_y_xxxx_xyyzzzz[k] * ab_x + g_0_y_xxxx_xxyyzzzz[k];

                g_0_y_xxxxx_xyzzzzz[k] = -g_0_y_xxxx_xyzzzzz[k] * ab_x + g_0_y_xxxx_xxyzzzzz[k];

                g_0_y_xxxxx_xzzzzzz[k] = -g_0_y_xxxx_xzzzzzz[k] * ab_x + g_0_y_xxxx_xxzzzzzz[k];

                g_0_y_xxxxx_yyyyyyy[k] = -g_0_y_xxxx_yyyyyyy[k] * ab_x + g_0_y_xxxx_xyyyyyyy[k];

                g_0_y_xxxxx_yyyyyyz[k] = -g_0_y_xxxx_yyyyyyz[k] * ab_x + g_0_y_xxxx_xyyyyyyz[k];

                g_0_y_xxxxx_yyyyyzz[k] = -g_0_y_xxxx_yyyyyzz[k] * ab_x + g_0_y_xxxx_xyyyyyzz[k];

                g_0_y_xxxxx_yyyyzzz[k] = -g_0_y_xxxx_yyyyzzz[k] * ab_x + g_0_y_xxxx_xyyyyzzz[k];

                g_0_y_xxxxx_yyyzzzz[k] = -g_0_y_xxxx_yyyzzzz[k] * ab_x + g_0_y_xxxx_xyyyzzzz[k];

                g_0_y_xxxxx_yyzzzzz[k] = -g_0_y_xxxx_yyzzzzz[k] * ab_x + g_0_y_xxxx_xyyzzzzz[k];

                g_0_y_xxxxx_yzzzzzz[k] = -g_0_y_xxxx_yzzzzzz[k] * ab_x + g_0_y_xxxx_xyzzzzzz[k];

                g_0_y_xxxxx_zzzzzzz[k] = -g_0_y_xxxx_zzzzzzz[k] * ab_x + g_0_y_xxxx_xzzzzzzz[k];
            }

            /// Set up 792-828 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxy_xxxxxxx = cbuffer.data(hk_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxxxy = cbuffer.data(hk_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxxxz = cbuffer.data(hk_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxxyy = cbuffer.data(hk_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxxyz = cbuffer.data(hk_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxxzz = cbuffer.data(hk_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxyyy = cbuffer.data(hk_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxyyz = cbuffer.data(hk_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxyzz = cbuffer.data(hk_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxzzz = cbuffer.data(hk_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyyyy = cbuffer.data(hk_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyyyz = cbuffer.data(hk_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyyzz = cbuffer.data(hk_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyzzz = cbuffer.data(hk_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxzzzz = cbuffer.data(hk_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyyyy = cbuffer.data(hk_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyyyz = cbuffer.data(hk_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyyzz = cbuffer.data(hk_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyzzz = cbuffer.data(hk_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyzzzz = cbuffer.data(hk_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxzzzzz = cbuffer.data(hk_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyyyy = cbuffer.data(hk_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyyyz = cbuffer.data(hk_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyyzz = cbuffer.data(hk_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyzzz = cbuffer.data(hk_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyzzzz = cbuffer.data(hk_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyzzzzz = cbuffer.data(hk_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_xxxxy_xzzzzzz = cbuffer.data(hk_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyyyy = cbuffer.data(hk_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyyyz = cbuffer.data(hk_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyyzz = cbuffer.data(hk_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyzzz = cbuffer.data(hk_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyzzzz = cbuffer.data(hk_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyzzzzz = cbuffer.data(hk_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_xxxxy_yzzzzzz = cbuffer.data(hk_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_xxxxy_zzzzzzz = cbuffer.data(hk_geom_01_off + 827 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxy_xxxxxxx, g_0_y_xxxxy_xxxxxxy, g_0_y_xxxxy_xxxxxxz, g_0_y_xxxxy_xxxxxyy, g_0_y_xxxxy_xxxxxyz, g_0_y_xxxxy_xxxxxzz, g_0_y_xxxxy_xxxxyyy, g_0_y_xxxxy_xxxxyyz, g_0_y_xxxxy_xxxxyzz, g_0_y_xxxxy_xxxxzzz, g_0_y_xxxxy_xxxyyyy, g_0_y_xxxxy_xxxyyyz, g_0_y_xxxxy_xxxyyzz, g_0_y_xxxxy_xxxyzzz, g_0_y_xxxxy_xxxzzzz, g_0_y_xxxxy_xxyyyyy, g_0_y_xxxxy_xxyyyyz, g_0_y_xxxxy_xxyyyzz, g_0_y_xxxxy_xxyyzzz, g_0_y_xxxxy_xxyzzzz, g_0_y_xxxxy_xxzzzzz, g_0_y_xxxxy_xyyyyyy, g_0_y_xxxxy_xyyyyyz, g_0_y_xxxxy_xyyyyzz, g_0_y_xxxxy_xyyyzzz, g_0_y_xxxxy_xyyzzzz, g_0_y_xxxxy_xyzzzzz, g_0_y_xxxxy_xzzzzzz, g_0_y_xxxxy_yyyyyyy, g_0_y_xxxxy_yyyyyyz, g_0_y_xxxxy_yyyyyzz, g_0_y_xxxxy_yyyyzzz, g_0_y_xxxxy_yyyzzzz, g_0_y_xxxxy_yyzzzzz, g_0_y_xxxxy_yzzzzzz, g_0_y_xxxxy_zzzzzzz, g_0_y_xxxy_xxxxxxx, g_0_y_xxxy_xxxxxxxx, g_0_y_xxxy_xxxxxxxy, g_0_y_xxxy_xxxxxxxz, g_0_y_xxxy_xxxxxxy, g_0_y_xxxy_xxxxxxyy, g_0_y_xxxy_xxxxxxyz, g_0_y_xxxy_xxxxxxz, g_0_y_xxxy_xxxxxxzz, g_0_y_xxxy_xxxxxyy, g_0_y_xxxy_xxxxxyyy, g_0_y_xxxy_xxxxxyyz, g_0_y_xxxy_xxxxxyz, g_0_y_xxxy_xxxxxyzz, g_0_y_xxxy_xxxxxzz, g_0_y_xxxy_xxxxxzzz, g_0_y_xxxy_xxxxyyy, g_0_y_xxxy_xxxxyyyy, g_0_y_xxxy_xxxxyyyz, g_0_y_xxxy_xxxxyyz, g_0_y_xxxy_xxxxyyzz, g_0_y_xxxy_xxxxyzz, g_0_y_xxxy_xxxxyzzz, g_0_y_xxxy_xxxxzzz, g_0_y_xxxy_xxxxzzzz, g_0_y_xxxy_xxxyyyy, g_0_y_xxxy_xxxyyyyy, g_0_y_xxxy_xxxyyyyz, g_0_y_xxxy_xxxyyyz, g_0_y_xxxy_xxxyyyzz, g_0_y_xxxy_xxxyyzz, g_0_y_xxxy_xxxyyzzz, g_0_y_xxxy_xxxyzzz, g_0_y_xxxy_xxxyzzzz, g_0_y_xxxy_xxxzzzz, g_0_y_xxxy_xxxzzzzz, g_0_y_xxxy_xxyyyyy, g_0_y_xxxy_xxyyyyyy, g_0_y_xxxy_xxyyyyyz, g_0_y_xxxy_xxyyyyz, g_0_y_xxxy_xxyyyyzz, g_0_y_xxxy_xxyyyzz, g_0_y_xxxy_xxyyyzzz, g_0_y_xxxy_xxyyzzz, g_0_y_xxxy_xxyyzzzz, g_0_y_xxxy_xxyzzzz, g_0_y_xxxy_xxyzzzzz, g_0_y_xxxy_xxzzzzz, g_0_y_xxxy_xxzzzzzz, g_0_y_xxxy_xyyyyyy, g_0_y_xxxy_xyyyyyyy, g_0_y_xxxy_xyyyyyyz, g_0_y_xxxy_xyyyyyz, g_0_y_xxxy_xyyyyyzz, g_0_y_xxxy_xyyyyzz, g_0_y_xxxy_xyyyyzzz, g_0_y_xxxy_xyyyzzz, g_0_y_xxxy_xyyyzzzz, g_0_y_xxxy_xyyzzzz, g_0_y_xxxy_xyyzzzzz, g_0_y_xxxy_xyzzzzz, g_0_y_xxxy_xyzzzzzz, g_0_y_xxxy_xzzzzzz, g_0_y_xxxy_xzzzzzzz, g_0_y_xxxy_yyyyyyy, g_0_y_xxxy_yyyyyyz, g_0_y_xxxy_yyyyyzz, g_0_y_xxxy_yyyyzzz, g_0_y_xxxy_yyyzzzz, g_0_y_xxxy_yyzzzzz, g_0_y_xxxy_yzzzzzz, g_0_y_xxxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxy_xxxxxxx[k] = -g_0_y_xxxy_xxxxxxx[k] * ab_x + g_0_y_xxxy_xxxxxxxx[k];

                g_0_y_xxxxy_xxxxxxy[k] = -g_0_y_xxxy_xxxxxxy[k] * ab_x + g_0_y_xxxy_xxxxxxxy[k];

                g_0_y_xxxxy_xxxxxxz[k] = -g_0_y_xxxy_xxxxxxz[k] * ab_x + g_0_y_xxxy_xxxxxxxz[k];

                g_0_y_xxxxy_xxxxxyy[k] = -g_0_y_xxxy_xxxxxyy[k] * ab_x + g_0_y_xxxy_xxxxxxyy[k];

                g_0_y_xxxxy_xxxxxyz[k] = -g_0_y_xxxy_xxxxxyz[k] * ab_x + g_0_y_xxxy_xxxxxxyz[k];

                g_0_y_xxxxy_xxxxxzz[k] = -g_0_y_xxxy_xxxxxzz[k] * ab_x + g_0_y_xxxy_xxxxxxzz[k];

                g_0_y_xxxxy_xxxxyyy[k] = -g_0_y_xxxy_xxxxyyy[k] * ab_x + g_0_y_xxxy_xxxxxyyy[k];

                g_0_y_xxxxy_xxxxyyz[k] = -g_0_y_xxxy_xxxxyyz[k] * ab_x + g_0_y_xxxy_xxxxxyyz[k];

                g_0_y_xxxxy_xxxxyzz[k] = -g_0_y_xxxy_xxxxyzz[k] * ab_x + g_0_y_xxxy_xxxxxyzz[k];

                g_0_y_xxxxy_xxxxzzz[k] = -g_0_y_xxxy_xxxxzzz[k] * ab_x + g_0_y_xxxy_xxxxxzzz[k];

                g_0_y_xxxxy_xxxyyyy[k] = -g_0_y_xxxy_xxxyyyy[k] * ab_x + g_0_y_xxxy_xxxxyyyy[k];

                g_0_y_xxxxy_xxxyyyz[k] = -g_0_y_xxxy_xxxyyyz[k] * ab_x + g_0_y_xxxy_xxxxyyyz[k];

                g_0_y_xxxxy_xxxyyzz[k] = -g_0_y_xxxy_xxxyyzz[k] * ab_x + g_0_y_xxxy_xxxxyyzz[k];

                g_0_y_xxxxy_xxxyzzz[k] = -g_0_y_xxxy_xxxyzzz[k] * ab_x + g_0_y_xxxy_xxxxyzzz[k];

                g_0_y_xxxxy_xxxzzzz[k] = -g_0_y_xxxy_xxxzzzz[k] * ab_x + g_0_y_xxxy_xxxxzzzz[k];

                g_0_y_xxxxy_xxyyyyy[k] = -g_0_y_xxxy_xxyyyyy[k] * ab_x + g_0_y_xxxy_xxxyyyyy[k];

                g_0_y_xxxxy_xxyyyyz[k] = -g_0_y_xxxy_xxyyyyz[k] * ab_x + g_0_y_xxxy_xxxyyyyz[k];

                g_0_y_xxxxy_xxyyyzz[k] = -g_0_y_xxxy_xxyyyzz[k] * ab_x + g_0_y_xxxy_xxxyyyzz[k];

                g_0_y_xxxxy_xxyyzzz[k] = -g_0_y_xxxy_xxyyzzz[k] * ab_x + g_0_y_xxxy_xxxyyzzz[k];

                g_0_y_xxxxy_xxyzzzz[k] = -g_0_y_xxxy_xxyzzzz[k] * ab_x + g_0_y_xxxy_xxxyzzzz[k];

                g_0_y_xxxxy_xxzzzzz[k] = -g_0_y_xxxy_xxzzzzz[k] * ab_x + g_0_y_xxxy_xxxzzzzz[k];

                g_0_y_xxxxy_xyyyyyy[k] = -g_0_y_xxxy_xyyyyyy[k] * ab_x + g_0_y_xxxy_xxyyyyyy[k];

                g_0_y_xxxxy_xyyyyyz[k] = -g_0_y_xxxy_xyyyyyz[k] * ab_x + g_0_y_xxxy_xxyyyyyz[k];

                g_0_y_xxxxy_xyyyyzz[k] = -g_0_y_xxxy_xyyyyzz[k] * ab_x + g_0_y_xxxy_xxyyyyzz[k];

                g_0_y_xxxxy_xyyyzzz[k] = -g_0_y_xxxy_xyyyzzz[k] * ab_x + g_0_y_xxxy_xxyyyzzz[k];

                g_0_y_xxxxy_xyyzzzz[k] = -g_0_y_xxxy_xyyzzzz[k] * ab_x + g_0_y_xxxy_xxyyzzzz[k];

                g_0_y_xxxxy_xyzzzzz[k] = -g_0_y_xxxy_xyzzzzz[k] * ab_x + g_0_y_xxxy_xxyzzzzz[k];

                g_0_y_xxxxy_xzzzzzz[k] = -g_0_y_xxxy_xzzzzzz[k] * ab_x + g_0_y_xxxy_xxzzzzzz[k];

                g_0_y_xxxxy_yyyyyyy[k] = -g_0_y_xxxy_yyyyyyy[k] * ab_x + g_0_y_xxxy_xyyyyyyy[k];

                g_0_y_xxxxy_yyyyyyz[k] = -g_0_y_xxxy_yyyyyyz[k] * ab_x + g_0_y_xxxy_xyyyyyyz[k];

                g_0_y_xxxxy_yyyyyzz[k] = -g_0_y_xxxy_yyyyyzz[k] * ab_x + g_0_y_xxxy_xyyyyyzz[k];

                g_0_y_xxxxy_yyyyzzz[k] = -g_0_y_xxxy_yyyyzzz[k] * ab_x + g_0_y_xxxy_xyyyyzzz[k];

                g_0_y_xxxxy_yyyzzzz[k] = -g_0_y_xxxy_yyyzzzz[k] * ab_x + g_0_y_xxxy_xyyyzzzz[k];

                g_0_y_xxxxy_yyzzzzz[k] = -g_0_y_xxxy_yyzzzzz[k] * ab_x + g_0_y_xxxy_xyyzzzzz[k];

                g_0_y_xxxxy_yzzzzzz[k] = -g_0_y_xxxy_yzzzzzz[k] * ab_x + g_0_y_xxxy_xyzzzzzz[k];

                g_0_y_xxxxy_zzzzzzz[k] = -g_0_y_xxxy_zzzzzzz[k] * ab_x + g_0_y_xxxy_xzzzzzzz[k];
            }

            /// Set up 828-864 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxz_xxxxxxx = cbuffer.data(hk_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxxxy = cbuffer.data(hk_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxxxz = cbuffer.data(hk_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxxyy = cbuffer.data(hk_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxxyz = cbuffer.data(hk_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxxzz = cbuffer.data(hk_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxyyy = cbuffer.data(hk_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxyyz = cbuffer.data(hk_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxyzz = cbuffer.data(hk_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxzzz = cbuffer.data(hk_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyyyy = cbuffer.data(hk_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyyyz = cbuffer.data(hk_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyyzz = cbuffer.data(hk_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyzzz = cbuffer.data(hk_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxzzzz = cbuffer.data(hk_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyyyy = cbuffer.data(hk_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyyyz = cbuffer.data(hk_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyyzz = cbuffer.data(hk_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyzzz = cbuffer.data(hk_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyzzzz = cbuffer.data(hk_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxzzzzz = cbuffer.data(hk_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyyyy = cbuffer.data(hk_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyyyz = cbuffer.data(hk_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyyzz = cbuffer.data(hk_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyzzz = cbuffer.data(hk_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyzzzz = cbuffer.data(hk_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyzzzzz = cbuffer.data(hk_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_xxxxz_xzzzzzz = cbuffer.data(hk_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyyyy = cbuffer.data(hk_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyyyz = cbuffer.data(hk_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyyzz = cbuffer.data(hk_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyzzz = cbuffer.data(hk_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyzzzz = cbuffer.data(hk_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyzzzzz = cbuffer.data(hk_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_xxxxz_yzzzzzz = cbuffer.data(hk_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_xxxxz_zzzzzzz = cbuffer.data(hk_geom_01_off + 863 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxz_xxxxxxx, g_0_y_xxxxz_xxxxxxy, g_0_y_xxxxz_xxxxxxz, g_0_y_xxxxz_xxxxxyy, g_0_y_xxxxz_xxxxxyz, g_0_y_xxxxz_xxxxxzz, g_0_y_xxxxz_xxxxyyy, g_0_y_xxxxz_xxxxyyz, g_0_y_xxxxz_xxxxyzz, g_0_y_xxxxz_xxxxzzz, g_0_y_xxxxz_xxxyyyy, g_0_y_xxxxz_xxxyyyz, g_0_y_xxxxz_xxxyyzz, g_0_y_xxxxz_xxxyzzz, g_0_y_xxxxz_xxxzzzz, g_0_y_xxxxz_xxyyyyy, g_0_y_xxxxz_xxyyyyz, g_0_y_xxxxz_xxyyyzz, g_0_y_xxxxz_xxyyzzz, g_0_y_xxxxz_xxyzzzz, g_0_y_xxxxz_xxzzzzz, g_0_y_xxxxz_xyyyyyy, g_0_y_xxxxz_xyyyyyz, g_0_y_xxxxz_xyyyyzz, g_0_y_xxxxz_xyyyzzz, g_0_y_xxxxz_xyyzzzz, g_0_y_xxxxz_xyzzzzz, g_0_y_xxxxz_xzzzzzz, g_0_y_xxxxz_yyyyyyy, g_0_y_xxxxz_yyyyyyz, g_0_y_xxxxz_yyyyyzz, g_0_y_xxxxz_yyyyzzz, g_0_y_xxxxz_yyyzzzz, g_0_y_xxxxz_yyzzzzz, g_0_y_xxxxz_yzzzzzz, g_0_y_xxxxz_zzzzzzz, g_0_y_xxxz_xxxxxxx, g_0_y_xxxz_xxxxxxxx, g_0_y_xxxz_xxxxxxxy, g_0_y_xxxz_xxxxxxxz, g_0_y_xxxz_xxxxxxy, g_0_y_xxxz_xxxxxxyy, g_0_y_xxxz_xxxxxxyz, g_0_y_xxxz_xxxxxxz, g_0_y_xxxz_xxxxxxzz, g_0_y_xxxz_xxxxxyy, g_0_y_xxxz_xxxxxyyy, g_0_y_xxxz_xxxxxyyz, g_0_y_xxxz_xxxxxyz, g_0_y_xxxz_xxxxxyzz, g_0_y_xxxz_xxxxxzz, g_0_y_xxxz_xxxxxzzz, g_0_y_xxxz_xxxxyyy, g_0_y_xxxz_xxxxyyyy, g_0_y_xxxz_xxxxyyyz, g_0_y_xxxz_xxxxyyz, g_0_y_xxxz_xxxxyyzz, g_0_y_xxxz_xxxxyzz, g_0_y_xxxz_xxxxyzzz, g_0_y_xxxz_xxxxzzz, g_0_y_xxxz_xxxxzzzz, g_0_y_xxxz_xxxyyyy, g_0_y_xxxz_xxxyyyyy, g_0_y_xxxz_xxxyyyyz, g_0_y_xxxz_xxxyyyz, g_0_y_xxxz_xxxyyyzz, g_0_y_xxxz_xxxyyzz, g_0_y_xxxz_xxxyyzzz, g_0_y_xxxz_xxxyzzz, g_0_y_xxxz_xxxyzzzz, g_0_y_xxxz_xxxzzzz, g_0_y_xxxz_xxxzzzzz, g_0_y_xxxz_xxyyyyy, g_0_y_xxxz_xxyyyyyy, g_0_y_xxxz_xxyyyyyz, g_0_y_xxxz_xxyyyyz, g_0_y_xxxz_xxyyyyzz, g_0_y_xxxz_xxyyyzz, g_0_y_xxxz_xxyyyzzz, g_0_y_xxxz_xxyyzzz, g_0_y_xxxz_xxyyzzzz, g_0_y_xxxz_xxyzzzz, g_0_y_xxxz_xxyzzzzz, g_0_y_xxxz_xxzzzzz, g_0_y_xxxz_xxzzzzzz, g_0_y_xxxz_xyyyyyy, g_0_y_xxxz_xyyyyyyy, g_0_y_xxxz_xyyyyyyz, g_0_y_xxxz_xyyyyyz, g_0_y_xxxz_xyyyyyzz, g_0_y_xxxz_xyyyyzz, g_0_y_xxxz_xyyyyzzz, g_0_y_xxxz_xyyyzzz, g_0_y_xxxz_xyyyzzzz, g_0_y_xxxz_xyyzzzz, g_0_y_xxxz_xyyzzzzz, g_0_y_xxxz_xyzzzzz, g_0_y_xxxz_xyzzzzzz, g_0_y_xxxz_xzzzzzz, g_0_y_xxxz_xzzzzzzz, g_0_y_xxxz_yyyyyyy, g_0_y_xxxz_yyyyyyz, g_0_y_xxxz_yyyyyzz, g_0_y_xxxz_yyyyzzz, g_0_y_xxxz_yyyzzzz, g_0_y_xxxz_yyzzzzz, g_0_y_xxxz_yzzzzzz, g_0_y_xxxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxz_xxxxxxx[k] = -g_0_y_xxxz_xxxxxxx[k] * ab_x + g_0_y_xxxz_xxxxxxxx[k];

                g_0_y_xxxxz_xxxxxxy[k] = -g_0_y_xxxz_xxxxxxy[k] * ab_x + g_0_y_xxxz_xxxxxxxy[k];

                g_0_y_xxxxz_xxxxxxz[k] = -g_0_y_xxxz_xxxxxxz[k] * ab_x + g_0_y_xxxz_xxxxxxxz[k];

                g_0_y_xxxxz_xxxxxyy[k] = -g_0_y_xxxz_xxxxxyy[k] * ab_x + g_0_y_xxxz_xxxxxxyy[k];

                g_0_y_xxxxz_xxxxxyz[k] = -g_0_y_xxxz_xxxxxyz[k] * ab_x + g_0_y_xxxz_xxxxxxyz[k];

                g_0_y_xxxxz_xxxxxzz[k] = -g_0_y_xxxz_xxxxxzz[k] * ab_x + g_0_y_xxxz_xxxxxxzz[k];

                g_0_y_xxxxz_xxxxyyy[k] = -g_0_y_xxxz_xxxxyyy[k] * ab_x + g_0_y_xxxz_xxxxxyyy[k];

                g_0_y_xxxxz_xxxxyyz[k] = -g_0_y_xxxz_xxxxyyz[k] * ab_x + g_0_y_xxxz_xxxxxyyz[k];

                g_0_y_xxxxz_xxxxyzz[k] = -g_0_y_xxxz_xxxxyzz[k] * ab_x + g_0_y_xxxz_xxxxxyzz[k];

                g_0_y_xxxxz_xxxxzzz[k] = -g_0_y_xxxz_xxxxzzz[k] * ab_x + g_0_y_xxxz_xxxxxzzz[k];

                g_0_y_xxxxz_xxxyyyy[k] = -g_0_y_xxxz_xxxyyyy[k] * ab_x + g_0_y_xxxz_xxxxyyyy[k];

                g_0_y_xxxxz_xxxyyyz[k] = -g_0_y_xxxz_xxxyyyz[k] * ab_x + g_0_y_xxxz_xxxxyyyz[k];

                g_0_y_xxxxz_xxxyyzz[k] = -g_0_y_xxxz_xxxyyzz[k] * ab_x + g_0_y_xxxz_xxxxyyzz[k];

                g_0_y_xxxxz_xxxyzzz[k] = -g_0_y_xxxz_xxxyzzz[k] * ab_x + g_0_y_xxxz_xxxxyzzz[k];

                g_0_y_xxxxz_xxxzzzz[k] = -g_0_y_xxxz_xxxzzzz[k] * ab_x + g_0_y_xxxz_xxxxzzzz[k];

                g_0_y_xxxxz_xxyyyyy[k] = -g_0_y_xxxz_xxyyyyy[k] * ab_x + g_0_y_xxxz_xxxyyyyy[k];

                g_0_y_xxxxz_xxyyyyz[k] = -g_0_y_xxxz_xxyyyyz[k] * ab_x + g_0_y_xxxz_xxxyyyyz[k];

                g_0_y_xxxxz_xxyyyzz[k] = -g_0_y_xxxz_xxyyyzz[k] * ab_x + g_0_y_xxxz_xxxyyyzz[k];

                g_0_y_xxxxz_xxyyzzz[k] = -g_0_y_xxxz_xxyyzzz[k] * ab_x + g_0_y_xxxz_xxxyyzzz[k];

                g_0_y_xxxxz_xxyzzzz[k] = -g_0_y_xxxz_xxyzzzz[k] * ab_x + g_0_y_xxxz_xxxyzzzz[k];

                g_0_y_xxxxz_xxzzzzz[k] = -g_0_y_xxxz_xxzzzzz[k] * ab_x + g_0_y_xxxz_xxxzzzzz[k];

                g_0_y_xxxxz_xyyyyyy[k] = -g_0_y_xxxz_xyyyyyy[k] * ab_x + g_0_y_xxxz_xxyyyyyy[k];

                g_0_y_xxxxz_xyyyyyz[k] = -g_0_y_xxxz_xyyyyyz[k] * ab_x + g_0_y_xxxz_xxyyyyyz[k];

                g_0_y_xxxxz_xyyyyzz[k] = -g_0_y_xxxz_xyyyyzz[k] * ab_x + g_0_y_xxxz_xxyyyyzz[k];

                g_0_y_xxxxz_xyyyzzz[k] = -g_0_y_xxxz_xyyyzzz[k] * ab_x + g_0_y_xxxz_xxyyyzzz[k];

                g_0_y_xxxxz_xyyzzzz[k] = -g_0_y_xxxz_xyyzzzz[k] * ab_x + g_0_y_xxxz_xxyyzzzz[k];

                g_0_y_xxxxz_xyzzzzz[k] = -g_0_y_xxxz_xyzzzzz[k] * ab_x + g_0_y_xxxz_xxyzzzzz[k];

                g_0_y_xxxxz_xzzzzzz[k] = -g_0_y_xxxz_xzzzzzz[k] * ab_x + g_0_y_xxxz_xxzzzzzz[k];

                g_0_y_xxxxz_yyyyyyy[k] = -g_0_y_xxxz_yyyyyyy[k] * ab_x + g_0_y_xxxz_xyyyyyyy[k];

                g_0_y_xxxxz_yyyyyyz[k] = -g_0_y_xxxz_yyyyyyz[k] * ab_x + g_0_y_xxxz_xyyyyyyz[k];

                g_0_y_xxxxz_yyyyyzz[k] = -g_0_y_xxxz_yyyyyzz[k] * ab_x + g_0_y_xxxz_xyyyyyzz[k];

                g_0_y_xxxxz_yyyyzzz[k] = -g_0_y_xxxz_yyyyzzz[k] * ab_x + g_0_y_xxxz_xyyyyzzz[k];

                g_0_y_xxxxz_yyyzzzz[k] = -g_0_y_xxxz_yyyzzzz[k] * ab_x + g_0_y_xxxz_xyyyzzzz[k];

                g_0_y_xxxxz_yyzzzzz[k] = -g_0_y_xxxz_yyzzzzz[k] * ab_x + g_0_y_xxxz_xyyzzzzz[k];

                g_0_y_xxxxz_yzzzzzz[k] = -g_0_y_xxxz_yzzzzzz[k] * ab_x + g_0_y_xxxz_xyzzzzzz[k];

                g_0_y_xxxxz_zzzzzzz[k] = -g_0_y_xxxz_zzzzzzz[k] * ab_x + g_0_y_xxxz_xzzzzzzz[k];
            }

            /// Set up 864-900 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_xxxyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_xxxyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_xxxyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyy_xxxxxxx, g_0_y_xxxyy_xxxxxxy, g_0_y_xxxyy_xxxxxxz, g_0_y_xxxyy_xxxxxyy, g_0_y_xxxyy_xxxxxyz, g_0_y_xxxyy_xxxxxzz, g_0_y_xxxyy_xxxxyyy, g_0_y_xxxyy_xxxxyyz, g_0_y_xxxyy_xxxxyzz, g_0_y_xxxyy_xxxxzzz, g_0_y_xxxyy_xxxyyyy, g_0_y_xxxyy_xxxyyyz, g_0_y_xxxyy_xxxyyzz, g_0_y_xxxyy_xxxyzzz, g_0_y_xxxyy_xxxzzzz, g_0_y_xxxyy_xxyyyyy, g_0_y_xxxyy_xxyyyyz, g_0_y_xxxyy_xxyyyzz, g_0_y_xxxyy_xxyyzzz, g_0_y_xxxyy_xxyzzzz, g_0_y_xxxyy_xxzzzzz, g_0_y_xxxyy_xyyyyyy, g_0_y_xxxyy_xyyyyyz, g_0_y_xxxyy_xyyyyzz, g_0_y_xxxyy_xyyyzzz, g_0_y_xxxyy_xyyzzzz, g_0_y_xxxyy_xyzzzzz, g_0_y_xxxyy_xzzzzzz, g_0_y_xxxyy_yyyyyyy, g_0_y_xxxyy_yyyyyyz, g_0_y_xxxyy_yyyyyzz, g_0_y_xxxyy_yyyyzzz, g_0_y_xxxyy_yyyzzzz, g_0_y_xxxyy_yyzzzzz, g_0_y_xxxyy_yzzzzzz, g_0_y_xxxyy_zzzzzzz, g_0_y_xxyy_xxxxxxx, g_0_y_xxyy_xxxxxxxx, g_0_y_xxyy_xxxxxxxy, g_0_y_xxyy_xxxxxxxz, g_0_y_xxyy_xxxxxxy, g_0_y_xxyy_xxxxxxyy, g_0_y_xxyy_xxxxxxyz, g_0_y_xxyy_xxxxxxz, g_0_y_xxyy_xxxxxxzz, g_0_y_xxyy_xxxxxyy, g_0_y_xxyy_xxxxxyyy, g_0_y_xxyy_xxxxxyyz, g_0_y_xxyy_xxxxxyz, g_0_y_xxyy_xxxxxyzz, g_0_y_xxyy_xxxxxzz, g_0_y_xxyy_xxxxxzzz, g_0_y_xxyy_xxxxyyy, g_0_y_xxyy_xxxxyyyy, g_0_y_xxyy_xxxxyyyz, g_0_y_xxyy_xxxxyyz, g_0_y_xxyy_xxxxyyzz, g_0_y_xxyy_xxxxyzz, g_0_y_xxyy_xxxxyzzz, g_0_y_xxyy_xxxxzzz, g_0_y_xxyy_xxxxzzzz, g_0_y_xxyy_xxxyyyy, g_0_y_xxyy_xxxyyyyy, g_0_y_xxyy_xxxyyyyz, g_0_y_xxyy_xxxyyyz, g_0_y_xxyy_xxxyyyzz, g_0_y_xxyy_xxxyyzz, g_0_y_xxyy_xxxyyzzz, g_0_y_xxyy_xxxyzzz, g_0_y_xxyy_xxxyzzzz, g_0_y_xxyy_xxxzzzz, g_0_y_xxyy_xxxzzzzz, g_0_y_xxyy_xxyyyyy, g_0_y_xxyy_xxyyyyyy, g_0_y_xxyy_xxyyyyyz, g_0_y_xxyy_xxyyyyz, g_0_y_xxyy_xxyyyyzz, g_0_y_xxyy_xxyyyzz, g_0_y_xxyy_xxyyyzzz, g_0_y_xxyy_xxyyzzz, g_0_y_xxyy_xxyyzzzz, g_0_y_xxyy_xxyzzzz, g_0_y_xxyy_xxyzzzzz, g_0_y_xxyy_xxzzzzz, g_0_y_xxyy_xxzzzzzz, g_0_y_xxyy_xyyyyyy, g_0_y_xxyy_xyyyyyyy, g_0_y_xxyy_xyyyyyyz, g_0_y_xxyy_xyyyyyz, g_0_y_xxyy_xyyyyyzz, g_0_y_xxyy_xyyyyzz, g_0_y_xxyy_xyyyyzzz, g_0_y_xxyy_xyyyzzz, g_0_y_xxyy_xyyyzzzz, g_0_y_xxyy_xyyzzzz, g_0_y_xxyy_xyyzzzzz, g_0_y_xxyy_xyzzzzz, g_0_y_xxyy_xyzzzzzz, g_0_y_xxyy_xzzzzzz, g_0_y_xxyy_xzzzzzzz, g_0_y_xxyy_yyyyyyy, g_0_y_xxyy_yyyyyyz, g_0_y_xxyy_yyyyyzz, g_0_y_xxyy_yyyyzzz, g_0_y_xxyy_yyyzzzz, g_0_y_xxyy_yyzzzzz, g_0_y_xxyy_yzzzzzz, g_0_y_xxyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyy_xxxxxxx[k] = -g_0_y_xxyy_xxxxxxx[k] * ab_x + g_0_y_xxyy_xxxxxxxx[k];

                g_0_y_xxxyy_xxxxxxy[k] = -g_0_y_xxyy_xxxxxxy[k] * ab_x + g_0_y_xxyy_xxxxxxxy[k];

                g_0_y_xxxyy_xxxxxxz[k] = -g_0_y_xxyy_xxxxxxz[k] * ab_x + g_0_y_xxyy_xxxxxxxz[k];

                g_0_y_xxxyy_xxxxxyy[k] = -g_0_y_xxyy_xxxxxyy[k] * ab_x + g_0_y_xxyy_xxxxxxyy[k];

                g_0_y_xxxyy_xxxxxyz[k] = -g_0_y_xxyy_xxxxxyz[k] * ab_x + g_0_y_xxyy_xxxxxxyz[k];

                g_0_y_xxxyy_xxxxxzz[k] = -g_0_y_xxyy_xxxxxzz[k] * ab_x + g_0_y_xxyy_xxxxxxzz[k];

                g_0_y_xxxyy_xxxxyyy[k] = -g_0_y_xxyy_xxxxyyy[k] * ab_x + g_0_y_xxyy_xxxxxyyy[k];

                g_0_y_xxxyy_xxxxyyz[k] = -g_0_y_xxyy_xxxxyyz[k] * ab_x + g_0_y_xxyy_xxxxxyyz[k];

                g_0_y_xxxyy_xxxxyzz[k] = -g_0_y_xxyy_xxxxyzz[k] * ab_x + g_0_y_xxyy_xxxxxyzz[k];

                g_0_y_xxxyy_xxxxzzz[k] = -g_0_y_xxyy_xxxxzzz[k] * ab_x + g_0_y_xxyy_xxxxxzzz[k];

                g_0_y_xxxyy_xxxyyyy[k] = -g_0_y_xxyy_xxxyyyy[k] * ab_x + g_0_y_xxyy_xxxxyyyy[k];

                g_0_y_xxxyy_xxxyyyz[k] = -g_0_y_xxyy_xxxyyyz[k] * ab_x + g_0_y_xxyy_xxxxyyyz[k];

                g_0_y_xxxyy_xxxyyzz[k] = -g_0_y_xxyy_xxxyyzz[k] * ab_x + g_0_y_xxyy_xxxxyyzz[k];

                g_0_y_xxxyy_xxxyzzz[k] = -g_0_y_xxyy_xxxyzzz[k] * ab_x + g_0_y_xxyy_xxxxyzzz[k];

                g_0_y_xxxyy_xxxzzzz[k] = -g_0_y_xxyy_xxxzzzz[k] * ab_x + g_0_y_xxyy_xxxxzzzz[k];

                g_0_y_xxxyy_xxyyyyy[k] = -g_0_y_xxyy_xxyyyyy[k] * ab_x + g_0_y_xxyy_xxxyyyyy[k];

                g_0_y_xxxyy_xxyyyyz[k] = -g_0_y_xxyy_xxyyyyz[k] * ab_x + g_0_y_xxyy_xxxyyyyz[k];

                g_0_y_xxxyy_xxyyyzz[k] = -g_0_y_xxyy_xxyyyzz[k] * ab_x + g_0_y_xxyy_xxxyyyzz[k];

                g_0_y_xxxyy_xxyyzzz[k] = -g_0_y_xxyy_xxyyzzz[k] * ab_x + g_0_y_xxyy_xxxyyzzz[k];

                g_0_y_xxxyy_xxyzzzz[k] = -g_0_y_xxyy_xxyzzzz[k] * ab_x + g_0_y_xxyy_xxxyzzzz[k];

                g_0_y_xxxyy_xxzzzzz[k] = -g_0_y_xxyy_xxzzzzz[k] * ab_x + g_0_y_xxyy_xxxzzzzz[k];

                g_0_y_xxxyy_xyyyyyy[k] = -g_0_y_xxyy_xyyyyyy[k] * ab_x + g_0_y_xxyy_xxyyyyyy[k];

                g_0_y_xxxyy_xyyyyyz[k] = -g_0_y_xxyy_xyyyyyz[k] * ab_x + g_0_y_xxyy_xxyyyyyz[k];

                g_0_y_xxxyy_xyyyyzz[k] = -g_0_y_xxyy_xyyyyzz[k] * ab_x + g_0_y_xxyy_xxyyyyzz[k];

                g_0_y_xxxyy_xyyyzzz[k] = -g_0_y_xxyy_xyyyzzz[k] * ab_x + g_0_y_xxyy_xxyyyzzz[k];

                g_0_y_xxxyy_xyyzzzz[k] = -g_0_y_xxyy_xyyzzzz[k] * ab_x + g_0_y_xxyy_xxyyzzzz[k];

                g_0_y_xxxyy_xyzzzzz[k] = -g_0_y_xxyy_xyzzzzz[k] * ab_x + g_0_y_xxyy_xxyzzzzz[k];

                g_0_y_xxxyy_xzzzzzz[k] = -g_0_y_xxyy_xzzzzzz[k] * ab_x + g_0_y_xxyy_xxzzzzzz[k];

                g_0_y_xxxyy_yyyyyyy[k] = -g_0_y_xxyy_yyyyyyy[k] * ab_x + g_0_y_xxyy_xyyyyyyy[k];

                g_0_y_xxxyy_yyyyyyz[k] = -g_0_y_xxyy_yyyyyyz[k] * ab_x + g_0_y_xxyy_xyyyyyyz[k];

                g_0_y_xxxyy_yyyyyzz[k] = -g_0_y_xxyy_yyyyyzz[k] * ab_x + g_0_y_xxyy_xyyyyyzz[k];

                g_0_y_xxxyy_yyyyzzz[k] = -g_0_y_xxyy_yyyyzzz[k] * ab_x + g_0_y_xxyy_xyyyyzzz[k];

                g_0_y_xxxyy_yyyzzzz[k] = -g_0_y_xxyy_yyyzzzz[k] * ab_x + g_0_y_xxyy_xyyyzzzz[k];

                g_0_y_xxxyy_yyzzzzz[k] = -g_0_y_xxyy_yyzzzzz[k] * ab_x + g_0_y_xxyy_xyyzzzzz[k];

                g_0_y_xxxyy_yzzzzzz[k] = -g_0_y_xxyy_yzzzzzz[k] * ab_x + g_0_y_xxyy_xyzzzzzz[k];

                g_0_y_xxxyy_zzzzzzz[k] = -g_0_y_xxyy_zzzzzzz[k] * ab_x + g_0_y_xxyy_xzzzzzzz[k];
            }

            /// Set up 900-936 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_y_xxxyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_y_xxxyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_y_xxxyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 935 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyz_xxxxxxx, g_0_y_xxxyz_xxxxxxy, g_0_y_xxxyz_xxxxxxz, g_0_y_xxxyz_xxxxxyy, g_0_y_xxxyz_xxxxxyz, g_0_y_xxxyz_xxxxxzz, g_0_y_xxxyz_xxxxyyy, g_0_y_xxxyz_xxxxyyz, g_0_y_xxxyz_xxxxyzz, g_0_y_xxxyz_xxxxzzz, g_0_y_xxxyz_xxxyyyy, g_0_y_xxxyz_xxxyyyz, g_0_y_xxxyz_xxxyyzz, g_0_y_xxxyz_xxxyzzz, g_0_y_xxxyz_xxxzzzz, g_0_y_xxxyz_xxyyyyy, g_0_y_xxxyz_xxyyyyz, g_0_y_xxxyz_xxyyyzz, g_0_y_xxxyz_xxyyzzz, g_0_y_xxxyz_xxyzzzz, g_0_y_xxxyz_xxzzzzz, g_0_y_xxxyz_xyyyyyy, g_0_y_xxxyz_xyyyyyz, g_0_y_xxxyz_xyyyyzz, g_0_y_xxxyz_xyyyzzz, g_0_y_xxxyz_xyyzzzz, g_0_y_xxxyz_xyzzzzz, g_0_y_xxxyz_xzzzzzz, g_0_y_xxxyz_yyyyyyy, g_0_y_xxxyz_yyyyyyz, g_0_y_xxxyz_yyyyyzz, g_0_y_xxxyz_yyyyzzz, g_0_y_xxxyz_yyyzzzz, g_0_y_xxxyz_yyzzzzz, g_0_y_xxxyz_yzzzzzz, g_0_y_xxxyz_zzzzzzz, g_0_y_xxyz_xxxxxxx, g_0_y_xxyz_xxxxxxxx, g_0_y_xxyz_xxxxxxxy, g_0_y_xxyz_xxxxxxxz, g_0_y_xxyz_xxxxxxy, g_0_y_xxyz_xxxxxxyy, g_0_y_xxyz_xxxxxxyz, g_0_y_xxyz_xxxxxxz, g_0_y_xxyz_xxxxxxzz, g_0_y_xxyz_xxxxxyy, g_0_y_xxyz_xxxxxyyy, g_0_y_xxyz_xxxxxyyz, g_0_y_xxyz_xxxxxyz, g_0_y_xxyz_xxxxxyzz, g_0_y_xxyz_xxxxxzz, g_0_y_xxyz_xxxxxzzz, g_0_y_xxyz_xxxxyyy, g_0_y_xxyz_xxxxyyyy, g_0_y_xxyz_xxxxyyyz, g_0_y_xxyz_xxxxyyz, g_0_y_xxyz_xxxxyyzz, g_0_y_xxyz_xxxxyzz, g_0_y_xxyz_xxxxyzzz, g_0_y_xxyz_xxxxzzz, g_0_y_xxyz_xxxxzzzz, g_0_y_xxyz_xxxyyyy, g_0_y_xxyz_xxxyyyyy, g_0_y_xxyz_xxxyyyyz, g_0_y_xxyz_xxxyyyz, g_0_y_xxyz_xxxyyyzz, g_0_y_xxyz_xxxyyzz, g_0_y_xxyz_xxxyyzzz, g_0_y_xxyz_xxxyzzz, g_0_y_xxyz_xxxyzzzz, g_0_y_xxyz_xxxzzzz, g_0_y_xxyz_xxxzzzzz, g_0_y_xxyz_xxyyyyy, g_0_y_xxyz_xxyyyyyy, g_0_y_xxyz_xxyyyyyz, g_0_y_xxyz_xxyyyyz, g_0_y_xxyz_xxyyyyzz, g_0_y_xxyz_xxyyyzz, g_0_y_xxyz_xxyyyzzz, g_0_y_xxyz_xxyyzzz, g_0_y_xxyz_xxyyzzzz, g_0_y_xxyz_xxyzzzz, g_0_y_xxyz_xxyzzzzz, g_0_y_xxyz_xxzzzzz, g_0_y_xxyz_xxzzzzzz, g_0_y_xxyz_xyyyyyy, g_0_y_xxyz_xyyyyyyy, g_0_y_xxyz_xyyyyyyz, g_0_y_xxyz_xyyyyyz, g_0_y_xxyz_xyyyyyzz, g_0_y_xxyz_xyyyyzz, g_0_y_xxyz_xyyyyzzz, g_0_y_xxyz_xyyyzzz, g_0_y_xxyz_xyyyzzzz, g_0_y_xxyz_xyyzzzz, g_0_y_xxyz_xyyzzzzz, g_0_y_xxyz_xyzzzzz, g_0_y_xxyz_xyzzzzzz, g_0_y_xxyz_xzzzzzz, g_0_y_xxyz_xzzzzzzz, g_0_y_xxyz_yyyyyyy, g_0_y_xxyz_yyyyyyz, g_0_y_xxyz_yyyyyzz, g_0_y_xxyz_yyyyzzz, g_0_y_xxyz_yyyzzzz, g_0_y_xxyz_yyzzzzz, g_0_y_xxyz_yzzzzzz, g_0_y_xxyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyz_xxxxxxx[k] = -g_0_y_xxyz_xxxxxxx[k] * ab_x + g_0_y_xxyz_xxxxxxxx[k];

                g_0_y_xxxyz_xxxxxxy[k] = -g_0_y_xxyz_xxxxxxy[k] * ab_x + g_0_y_xxyz_xxxxxxxy[k];

                g_0_y_xxxyz_xxxxxxz[k] = -g_0_y_xxyz_xxxxxxz[k] * ab_x + g_0_y_xxyz_xxxxxxxz[k];

                g_0_y_xxxyz_xxxxxyy[k] = -g_0_y_xxyz_xxxxxyy[k] * ab_x + g_0_y_xxyz_xxxxxxyy[k];

                g_0_y_xxxyz_xxxxxyz[k] = -g_0_y_xxyz_xxxxxyz[k] * ab_x + g_0_y_xxyz_xxxxxxyz[k];

                g_0_y_xxxyz_xxxxxzz[k] = -g_0_y_xxyz_xxxxxzz[k] * ab_x + g_0_y_xxyz_xxxxxxzz[k];

                g_0_y_xxxyz_xxxxyyy[k] = -g_0_y_xxyz_xxxxyyy[k] * ab_x + g_0_y_xxyz_xxxxxyyy[k];

                g_0_y_xxxyz_xxxxyyz[k] = -g_0_y_xxyz_xxxxyyz[k] * ab_x + g_0_y_xxyz_xxxxxyyz[k];

                g_0_y_xxxyz_xxxxyzz[k] = -g_0_y_xxyz_xxxxyzz[k] * ab_x + g_0_y_xxyz_xxxxxyzz[k];

                g_0_y_xxxyz_xxxxzzz[k] = -g_0_y_xxyz_xxxxzzz[k] * ab_x + g_0_y_xxyz_xxxxxzzz[k];

                g_0_y_xxxyz_xxxyyyy[k] = -g_0_y_xxyz_xxxyyyy[k] * ab_x + g_0_y_xxyz_xxxxyyyy[k];

                g_0_y_xxxyz_xxxyyyz[k] = -g_0_y_xxyz_xxxyyyz[k] * ab_x + g_0_y_xxyz_xxxxyyyz[k];

                g_0_y_xxxyz_xxxyyzz[k] = -g_0_y_xxyz_xxxyyzz[k] * ab_x + g_0_y_xxyz_xxxxyyzz[k];

                g_0_y_xxxyz_xxxyzzz[k] = -g_0_y_xxyz_xxxyzzz[k] * ab_x + g_0_y_xxyz_xxxxyzzz[k];

                g_0_y_xxxyz_xxxzzzz[k] = -g_0_y_xxyz_xxxzzzz[k] * ab_x + g_0_y_xxyz_xxxxzzzz[k];

                g_0_y_xxxyz_xxyyyyy[k] = -g_0_y_xxyz_xxyyyyy[k] * ab_x + g_0_y_xxyz_xxxyyyyy[k];

                g_0_y_xxxyz_xxyyyyz[k] = -g_0_y_xxyz_xxyyyyz[k] * ab_x + g_0_y_xxyz_xxxyyyyz[k];

                g_0_y_xxxyz_xxyyyzz[k] = -g_0_y_xxyz_xxyyyzz[k] * ab_x + g_0_y_xxyz_xxxyyyzz[k];

                g_0_y_xxxyz_xxyyzzz[k] = -g_0_y_xxyz_xxyyzzz[k] * ab_x + g_0_y_xxyz_xxxyyzzz[k];

                g_0_y_xxxyz_xxyzzzz[k] = -g_0_y_xxyz_xxyzzzz[k] * ab_x + g_0_y_xxyz_xxxyzzzz[k];

                g_0_y_xxxyz_xxzzzzz[k] = -g_0_y_xxyz_xxzzzzz[k] * ab_x + g_0_y_xxyz_xxxzzzzz[k];

                g_0_y_xxxyz_xyyyyyy[k] = -g_0_y_xxyz_xyyyyyy[k] * ab_x + g_0_y_xxyz_xxyyyyyy[k];

                g_0_y_xxxyz_xyyyyyz[k] = -g_0_y_xxyz_xyyyyyz[k] * ab_x + g_0_y_xxyz_xxyyyyyz[k];

                g_0_y_xxxyz_xyyyyzz[k] = -g_0_y_xxyz_xyyyyzz[k] * ab_x + g_0_y_xxyz_xxyyyyzz[k];

                g_0_y_xxxyz_xyyyzzz[k] = -g_0_y_xxyz_xyyyzzz[k] * ab_x + g_0_y_xxyz_xxyyyzzz[k];

                g_0_y_xxxyz_xyyzzzz[k] = -g_0_y_xxyz_xyyzzzz[k] * ab_x + g_0_y_xxyz_xxyyzzzz[k];

                g_0_y_xxxyz_xyzzzzz[k] = -g_0_y_xxyz_xyzzzzz[k] * ab_x + g_0_y_xxyz_xxyzzzzz[k];

                g_0_y_xxxyz_xzzzzzz[k] = -g_0_y_xxyz_xzzzzzz[k] * ab_x + g_0_y_xxyz_xxzzzzzz[k];

                g_0_y_xxxyz_yyyyyyy[k] = -g_0_y_xxyz_yyyyyyy[k] * ab_x + g_0_y_xxyz_xyyyyyyy[k];

                g_0_y_xxxyz_yyyyyyz[k] = -g_0_y_xxyz_yyyyyyz[k] * ab_x + g_0_y_xxyz_xyyyyyyz[k];

                g_0_y_xxxyz_yyyyyzz[k] = -g_0_y_xxyz_yyyyyzz[k] * ab_x + g_0_y_xxyz_xyyyyyzz[k];

                g_0_y_xxxyz_yyyyzzz[k] = -g_0_y_xxyz_yyyyzzz[k] * ab_x + g_0_y_xxyz_xyyyyzzz[k];

                g_0_y_xxxyz_yyyzzzz[k] = -g_0_y_xxyz_yyyzzzz[k] * ab_x + g_0_y_xxyz_xyyyzzzz[k];

                g_0_y_xxxyz_yyzzzzz[k] = -g_0_y_xxyz_yyzzzzz[k] * ab_x + g_0_y_xxyz_xyyzzzzz[k];

                g_0_y_xxxyz_yzzzzzz[k] = -g_0_y_xxyz_yzzzzzz[k] * ab_x + g_0_y_xxyz_xyzzzzzz[k];

                g_0_y_xxxyz_zzzzzzz[k] = -g_0_y_xxyz_zzzzzzz[k] * ab_x + g_0_y_xxyz_xzzzzzzz[k];
            }

            /// Set up 936-972 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_y_xxxzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_y_xxxzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_y_xxxzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 971 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzz_xxxxxxx, g_0_y_xxxzz_xxxxxxy, g_0_y_xxxzz_xxxxxxz, g_0_y_xxxzz_xxxxxyy, g_0_y_xxxzz_xxxxxyz, g_0_y_xxxzz_xxxxxzz, g_0_y_xxxzz_xxxxyyy, g_0_y_xxxzz_xxxxyyz, g_0_y_xxxzz_xxxxyzz, g_0_y_xxxzz_xxxxzzz, g_0_y_xxxzz_xxxyyyy, g_0_y_xxxzz_xxxyyyz, g_0_y_xxxzz_xxxyyzz, g_0_y_xxxzz_xxxyzzz, g_0_y_xxxzz_xxxzzzz, g_0_y_xxxzz_xxyyyyy, g_0_y_xxxzz_xxyyyyz, g_0_y_xxxzz_xxyyyzz, g_0_y_xxxzz_xxyyzzz, g_0_y_xxxzz_xxyzzzz, g_0_y_xxxzz_xxzzzzz, g_0_y_xxxzz_xyyyyyy, g_0_y_xxxzz_xyyyyyz, g_0_y_xxxzz_xyyyyzz, g_0_y_xxxzz_xyyyzzz, g_0_y_xxxzz_xyyzzzz, g_0_y_xxxzz_xyzzzzz, g_0_y_xxxzz_xzzzzzz, g_0_y_xxxzz_yyyyyyy, g_0_y_xxxzz_yyyyyyz, g_0_y_xxxzz_yyyyyzz, g_0_y_xxxzz_yyyyzzz, g_0_y_xxxzz_yyyzzzz, g_0_y_xxxzz_yyzzzzz, g_0_y_xxxzz_yzzzzzz, g_0_y_xxxzz_zzzzzzz, g_0_y_xxzz_xxxxxxx, g_0_y_xxzz_xxxxxxxx, g_0_y_xxzz_xxxxxxxy, g_0_y_xxzz_xxxxxxxz, g_0_y_xxzz_xxxxxxy, g_0_y_xxzz_xxxxxxyy, g_0_y_xxzz_xxxxxxyz, g_0_y_xxzz_xxxxxxz, g_0_y_xxzz_xxxxxxzz, g_0_y_xxzz_xxxxxyy, g_0_y_xxzz_xxxxxyyy, g_0_y_xxzz_xxxxxyyz, g_0_y_xxzz_xxxxxyz, g_0_y_xxzz_xxxxxyzz, g_0_y_xxzz_xxxxxzz, g_0_y_xxzz_xxxxxzzz, g_0_y_xxzz_xxxxyyy, g_0_y_xxzz_xxxxyyyy, g_0_y_xxzz_xxxxyyyz, g_0_y_xxzz_xxxxyyz, g_0_y_xxzz_xxxxyyzz, g_0_y_xxzz_xxxxyzz, g_0_y_xxzz_xxxxyzzz, g_0_y_xxzz_xxxxzzz, g_0_y_xxzz_xxxxzzzz, g_0_y_xxzz_xxxyyyy, g_0_y_xxzz_xxxyyyyy, g_0_y_xxzz_xxxyyyyz, g_0_y_xxzz_xxxyyyz, g_0_y_xxzz_xxxyyyzz, g_0_y_xxzz_xxxyyzz, g_0_y_xxzz_xxxyyzzz, g_0_y_xxzz_xxxyzzz, g_0_y_xxzz_xxxyzzzz, g_0_y_xxzz_xxxzzzz, g_0_y_xxzz_xxxzzzzz, g_0_y_xxzz_xxyyyyy, g_0_y_xxzz_xxyyyyyy, g_0_y_xxzz_xxyyyyyz, g_0_y_xxzz_xxyyyyz, g_0_y_xxzz_xxyyyyzz, g_0_y_xxzz_xxyyyzz, g_0_y_xxzz_xxyyyzzz, g_0_y_xxzz_xxyyzzz, g_0_y_xxzz_xxyyzzzz, g_0_y_xxzz_xxyzzzz, g_0_y_xxzz_xxyzzzzz, g_0_y_xxzz_xxzzzzz, g_0_y_xxzz_xxzzzzzz, g_0_y_xxzz_xyyyyyy, g_0_y_xxzz_xyyyyyyy, g_0_y_xxzz_xyyyyyyz, g_0_y_xxzz_xyyyyyz, g_0_y_xxzz_xyyyyyzz, g_0_y_xxzz_xyyyyzz, g_0_y_xxzz_xyyyyzzz, g_0_y_xxzz_xyyyzzz, g_0_y_xxzz_xyyyzzzz, g_0_y_xxzz_xyyzzzz, g_0_y_xxzz_xyyzzzzz, g_0_y_xxzz_xyzzzzz, g_0_y_xxzz_xyzzzzzz, g_0_y_xxzz_xzzzzzz, g_0_y_xxzz_xzzzzzzz, g_0_y_xxzz_yyyyyyy, g_0_y_xxzz_yyyyyyz, g_0_y_xxzz_yyyyyzz, g_0_y_xxzz_yyyyzzz, g_0_y_xxzz_yyyzzzz, g_0_y_xxzz_yyzzzzz, g_0_y_xxzz_yzzzzzz, g_0_y_xxzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzz_xxxxxxx[k] = -g_0_y_xxzz_xxxxxxx[k] * ab_x + g_0_y_xxzz_xxxxxxxx[k];

                g_0_y_xxxzz_xxxxxxy[k] = -g_0_y_xxzz_xxxxxxy[k] * ab_x + g_0_y_xxzz_xxxxxxxy[k];

                g_0_y_xxxzz_xxxxxxz[k] = -g_0_y_xxzz_xxxxxxz[k] * ab_x + g_0_y_xxzz_xxxxxxxz[k];

                g_0_y_xxxzz_xxxxxyy[k] = -g_0_y_xxzz_xxxxxyy[k] * ab_x + g_0_y_xxzz_xxxxxxyy[k];

                g_0_y_xxxzz_xxxxxyz[k] = -g_0_y_xxzz_xxxxxyz[k] * ab_x + g_0_y_xxzz_xxxxxxyz[k];

                g_0_y_xxxzz_xxxxxzz[k] = -g_0_y_xxzz_xxxxxzz[k] * ab_x + g_0_y_xxzz_xxxxxxzz[k];

                g_0_y_xxxzz_xxxxyyy[k] = -g_0_y_xxzz_xxxxyyy[k] * ab_x + g_0_y_xxzz_xxxxxyyy[k];

                g_0_y_xxxzz_xxxxyyz[k] = -g_0_y_xxzz_xxxxyyz[k] * ab_x + g_0_y_xxzz_xxxxxyyz[k];

                g_0_y_xxxzz_xxxxyzz[k] = -g_0_y_xxzz_xxxxyzz[k] * ab_x + g_0_y_xxzz_xxxxxyzz[k];

                g_0_y_xxxzz_xxxxzzz[k] = -g_0_y_xxzz_xxxxzzz[k] * ab_x + g_0_y_xxzz_xxxxxzzz[k];

                g_0_y_xxxzz_xxxyyyy[k] = -g_0_y_xxzz_xxxyyyy[k] * ab_x + g_0_y_xxzz_xxxxyyyy[k];

                g_0_y_xxxzz_xxxyyyz[k] = -g_0_y_xxzz_xxxyyyz[k] * ab_x + g_0_y_xxzz_xxxxyyyz[k];

                g_0_y_xxxzz_xxxyyzz[k] = -g_0_y_xxzz_xxxyyzz[k] * ab_x + g_0_y_xxzz_xxxxyyzz[k];

                g_0_y_xxxzz_xxxyzzz[k] = -g_0_y_xxzz_xxxyzzz[k] * ab_x + g_0_y_xxzz_xxxxyzzz[k];

                g_0_y_xxxzz_xxxzzzz[k] = -g_0_y_xxzz_xxxzzzz[k] * ab_x + g_0_y_xxzz_xxxxzzzz[k];

                g_0_y_xxxzz_xxyyyyy[k] = -g_0_y_xxzz_xxyyyyy[k] * ab_x + g_0_y_xxzz_xxxyyyyy[k];

                g_0_y_xxxzz_xxyyyyz[k] = -g_0_y_xxzz_xxyyyyz[k] * ab_x + g_0_y_xxzz_xxxyyyyz[k];

                g_0_y_xxxzz_xxyyyzz[k] = -g_0_y_xxzz_xxyyyzz[k] * ab_x + g_0_y_xxzz_xxxyyyzz[k];

                g_0_y_xxxzz_xxyyzzz[k] = -g_0_y_xxzz_xxyyzzz[k] * ab_x + g_0_y_xxzz_xxxyyzzz[k];

                g_0_y_xxxzz_xxyzzzz[k] = -g_0_y_xxzz_xxyzzzz[k] * ab_x + g_0_y_xxzz_xxxyzzzz[k];

                g_0_y_xxxzz_xxzzzzz[k] = -g_0_y_xxzz_xxzzzzz[k] * ab_x + g_0_y_xxzz_xxxzzzzz[k];

                g_0_y_xxxzz_xyyyyyy[k] = -g_0_y_xxzz_xyyyyyy[k] * ab_x + g_0_y_xxzz_xxyyyyyy[k];

                g_0_y_xxxzz_xyyyyyz[k] = -g_0_y_xxzz_xyyyyyz[k] * ab_x + g_0_y_xxzz_xxyyyyyz[k];

                g_0_y_xxxzz_xyyyyzz[k] = -g_0_y_xxzz_xyyyyzz[k] * ab_x + g_0_y_xxzz_xxyyyyzz[k];

                g_0_y_xxxzz_xyyyzzz[k] = -g_0_y_xxzz_xyyyzzz[k] * ab_x + g_0_y_xxzz_xxyyyzzz[k];

                g_0_y_xxxzz_xyyzzzz[k] = -g_0_y_xxzz_xyyzzzz[k] * ab_x + g_0_y_xxzz_xxyyzzzz[k];

                g_0_y_xxxzz_xyzzzzz[k] = -g_0_y_xxzz_xyzzzzz[k] * ab_x + g_0_y_xxzz_xxyzzzzz[k];

                g_0_y_xxxzz_xzzzzzz[k] = -g_0_y_xxzz_xzzzzzz[k] * ab_x + g_0_y_xxzz_xxzzzzzz[k];

                g_0_y_xxxzz_yyyyyyy[k] = -g_0_y_xxzz_yyyyyyy[k] * ab_x + g_0_y_xxzz_xyyyyyyy[k];

                g_0_y_xxxzz_yyyyyyz[k] = -g_0_y_xxzz_yyyyyyz[k] * ab_x + g_0_y_xxzz_xyyyyyyz[k];

                g_0_y_xxxzz_yyyyyzz[k] = -g_0_y_xxzz_yyyyyzz[k] * ab_x + g_0_y_xxzz_xyyyyyzz[k];

                g_0_y_xxxzz_yyyyzzz[k] = -g_0_y_xxzz_yyyyzzz[k] * ab_x + g_0_y_xxzz_xyyyyzzz[k];

                g_0_y_xxxzz_yyyzzzz[k] = -g_0_y_xxzz_yyyzzzz[k] * ab_x + g_0_y_xxzz_xyyyzzzz[k];

                g_0_y_xxxzz_yyzzzzz[k] = -g_0_y_xxzz_yyzzzzz[k] * ab_x + g_0_y_xxzz_xyyzzzzz[k];

                g_0_y_xxxzz_yzzzzzz[k] = -g_0_y_xxzz_yzzzzzz[k] * ab_x + g_0_y_xxzz_xyzzzzzz[k];

                g_0_y_xxxzz_zzzzzzz[k] = -g_0_y_xxzz_zzzzzzz[k] * ab_x + g_0_y_xxzz_xzzzzzzz[k];
            }

            /// Set up 972-1008 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_y_xxyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_y_xxyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_y_xxyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyy_xxxxxxx, g_0_y_xxyyy_xxxxxxy, g_0_y_xxyyy_xxxxxxz, g_0_y_xxyyy_xxxxxyy, g_0_y_xxyyy_xxxxxyz, g_0_y_xxyyy_xxxxxzz, g_0_y_xxyyy_xxxxyyy, g_0_y_xxyyy_xxxxyyz, g_0_y_xxyyy_xxxxyzz, g_0_y_xxyyy_xxxxzzz, g_0_y_xxyyy_xxxyyyy, g_0_y_xxyyy_xxxyyyz, g_0_y_xxyyy_xxxyyzz, g_0_y_xxyyy_xxxyzzz, g_0_y_xxyyy_xxxzzzz, g_0_y_xxyyy_xxyyyyy, g_0_y_xxyyy_xxyyyyz, g_0_y_xxyyy_xxyyyzz, g_0_y_xxyyy_xxyyzzz, g_0_y_xxyyy_xxyzzzz, g_0_y_xxyyy_xxzzzzz, g_0_y_xxyyy_xyyyyyy, g_0_y_xxyyy_xyyyyyz, g_0_y_xxyyy_xyyyyzz, g_0_y_xxyyy_xyyyzzz, g_0_y_xxyyy_xyyzzzz, g_0_y_xxyyy_xyzzzzz, g_0_y_xxyyy_xzzzzzz, g_0_y_xxyyy_yyyyyyy, g_0_y_xxyyy_yyyyyyz, g_0_y_xxyyy_yyyyyzz, g_0_y_xxyyy_yyyyzzz, g_0_y_xxyyy_yyyzzzz, g_0_y_xxyyy_yyzzzzz, g_0_y_xxyyy_yzzzzzz, g_0_y_xxyyy_zzzzzzz, g_0_y_xyyy_xxxxxxx, g_0_y_xyyy_xxxxxxxx, g_0_y_xyyy_xxxxxxxy, g_0_y_xyyy_xxxxxxxz, g_0_y_xyyy_xxxxxxy, g_0_y_xyyy_xxxxxxyy, g_0_y_xyyy_xxxxxxyz, g_0_y_xyyy_xxxxxxz, g_0_y_xyyy_xxxxxxzz, g_0_y_xyyy_xxxxxyy, g_0_y_xyyy_xxxxxyyy, g_0_y_xyyy_xxxxxyyz, g_0_y_xyyy_xxxxxyz, g_0_y_xyyy_xxxxxyzz, g_0_y_xyyy_xxxxxzz, g_0_y_xyyy_xxxxxzzz, g_0_y_xyyy_xxxxyyy, g_0_y_xyyy_xxxxyyyy, g_0_y_xyyy_xxxxyyyz, g_0_y_xyyy_xxxxyyz, g_0_y_xyyy_xxxxyyzz, g_0_y_xyyy_xxxxyzz, g_0_y_xyyy_xxxxyzzz, g_0_y_xyyy_xxxxzzz, g_0_y_xyyy_xxxxzzzz, g_0_y_xyyy_xxxyyyy, g_0_y_xyyy_xxxyyyyy, g_0_y_xyyy_xxxyyyyz, g_0_y_xyyy_xxxyyyz, g_0_y_xyyy_xxxyyyzz, g_0_y_xyyy_xxxyyzz, g_0_y_xyyy_xxxyyzzz, g_0_y_xyyy_xxxyzzz, g_0_y_xyyy_xxxyzzzz, g_0_y_xyyy_xxxzzzz, g_0_y_xyyy_xxxzzzzz, g_0_y_xyyy_xxyyyyy, g_0_y_xyyy_xxyyyyyy, g_0_y_xyyy_xxyyyyyz, g_0_y_xyyy_xxyyyyz, g_0_y_xyyy_xxyyyyzz, g_0_y_xyyy_xxyyyzz, g_0_y_xyyy_xxyyyzzz, g_0_y_xyyy_xxyyzzz, g_0_y_xyyy_xxyyzzzz, g_0_y_xyyy_xxyzzzz, g_0_y_xyyy_xxyzzzzz, g_0_y_xyyy_xxzzzzz, g_0_y_xyyy_xxzzzzzz, g_0_y_xyyy_xyyyyyy, g_0_y_xyyy_xyyyyyyy, g_0_y_xyyy_xyyyyyyz, g_0_y_xyyy_xyyyyyz, g_0_y_xyyy_xyyyyyzz, g_0_y_xyyy_xyyyyzz, g_0_y_xyyy_xyyyyzzz, g_0_y_xyyy_xyyyzzz, g_0_y_xyyy_xyyyzzzz, g_0_y_xyyy_xyyzzzz, g_0_y_xyyy_xyyzzzzz, g_0_y_xyyy_xyzzzzz, g_0_y_xyyy_xyzzzzzz, g_0_y_xyyy_xzzzzzz, g_0_y_xyyy_xzzzzzzz, g_0_y_xyyy_yyyyyyy, g_0_y_xyyy_yyyyyyz, g_0_y_xyyy_yyyyyzz, g_0_y_xyyy_yyyyzzz, g_0_y_xyyy_yyyzzzz, g_0_y_xyyy_yyzzzzz, g_0_y_xyyy_yzzzzzz, g_0_y_xyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyy_xxxxxxx[k] = -g_0_y_xyyy_xxxxxxx[k] * ab_x + g_0_y_xyyy_xxxxxxxx[k];

                g_0_y_xxyyy_xxxxxxy[k] = -g_0_y_xyyy_xxxxxxy[k] * ab_x + g_0_y_xyyy_xxxxxxxy[k];

                g_0_y_xxyyy_xxxxxxz[k] = -g_0_y_xyyy_xxxxxxz[k] * ab_x + g_0_y_xyyy_xxxxxxxz[k];

                g_0_y_xxyyy_xxxxxyy[k] = -g_0_y_xyyy_xxxxxyy[k] * ab_x + g_0_y_xyyy_xxxxxxyy[k];

                g_0_y_xxyyy_xxxxxyz[k] = -g_0_y_xyyy_xxxxxyz[k] * ab_x + g_0_y_xyyy_xxxxxxyz[k];

                g_0_y_xxyyy_xxxxxzz[k] = -g_0_y_xyyy_xxxxxzz[k] * ab_x + g_0_y_xyyy_xxxxxxzz[k];

                g_0_y_xxyyy_xxxxyyy[k] = -g_0_y_xyyy_xxxxyyy[k] * ab_x + g_0_y_xyyy_xxxxxyyy[k];

                g_0_y_xxyyy_xxxxyyz[k] = -g_0_y_xyyy_xxxxyyz[k] * ab_x + g_0_y_xyyy_xxxxxyyz[k];

                g_0_y_xxyyy_xxxxyzz[k] = -g_0_y_xyyy_xxxxyzz[k] * ab_x + g_0_y_xyyy_xxxxxyzz[k];

                g_0_y_xxyyy_xxxxzzz[k] = -g_0_y_xyyy_xxxxzzz[k] * ab_x + g_0_y_xyyy_xxxxxzzz[k];

                g_0_y_xxyyy_xxxyyyy[k] = -g_0_y_xyyy_xxxyyyy[k] * ab_x + g_0_y_xyyy_xxxxyyyy[k];

                g_0_y_xxyyy_xxxyyyz[k] = -g_0_y_xyyy_xxxyyyz[k] * ab_x + g_0_y_xyyy_xxxxyyyz[k];

                g_0_y_xxyyy_xxxyyzz[k] = -g_0_y_xyyy_xxxyyzz[k] * ab_x + g_0_y_xyyy_xxxxyyzz[k];

                g_0_y_xxyyy_xxxyzzz[k] = -g_0_y_xyyy_xxxyzzz[k] * ab_x + g_0_y_xyyy_xxxxyzzz[k];

                g_0_y_xxyyy_xxxzzzz[k] = -g_0_y_xyyy_xxxzzzz[k] * ab_x + g_0_y_xyyy_xxxxzzzz[k];

                g_0_y_xxyyy_xxyyyyy[k] = -g_0_y_xyyy_xxyyyyy[k] * ab_x + g_0_y_xyyy_xxxyyyyy[k];

                g_0_y_xxyyy_xxyyyyz[k] = -g_0_y_xyyy_xxyyyyz[k] * ab_x + g_0_y_xyyy_xxxyyyyz[k];

                g_0_y_xxyyy_xxyyyzz[k] = -g_0_y_xyyy_xxyyyzz[k] * ab_x + g_0_y_xyyy_xxxyyyzz[k];

                g_0_y_xxyyy_xxyyzzz[k] = -g_0_y_xyyy_xxyyzzz[k] * ab_x + g_0_y_xyyy_xxxyyzzz[k];

                g_0_y_xxyyy_xxyzzzz[k] = -g_0_y_xyyy_xxyzzzz[k] * ab_x + g_0_y_xyyy_xxxyzzzz[k];

                g_0_y_xxyyy_xxzzzzz[k] = -g_0_y_xyyy_xxzzzzz[k] * ab_x + g_0_y_xyyy_xxxzzzzz[k];

                g_0_y_xxyyy_xyyyyyy[k] = -g_0_y_xyyy_xyyyyyy[k] * ab_x + g_0_y_xyyy_xxyyyyyy[k];

                g_0_y_xxyyy_xyyyyyz[k] = -g_0_y_xyyy_xyyyyyz[k] * ab_x + g_0_y_xyyy_xxyyyyyz[k];

                g_0_y_xxyyy_xyyyyzz[k] = -g_0_y_xyyy_xyyyyzz[k] * ab_x + g_0_y_xyyy_xxyyyyzz[k];

                g_0_y_xxyyy_xyyyzzz[k] = -g_0_y_xyyy_xyyyzzz[k] * ab_x + g_0_y_xyyy_xxyyyzzz[k];

                g_0_y_xxyyy_xyyzzzz[k] = -g_0_y_xyyy_xyyzzzz[k] * ab_x + g_0_y_xyyy_xxyyzzzz[k];

                g_0_y_xxyyy_xyzzzzz[k] = -g_0_y_xyyy_xyzzzzz[k] * ab_x + g_0_y_xyyy_xxyzzzzz[k];

                g_0_y_xxyyy_xzzzzzz[k] = -g_0_y_xyyy_xzzzzzz[k] * ab_x + g_0_y_xyyy_xxzzzzzz[k];

                g_0_y_xxyyy_yyyyyyy[k] = -g_0_y_xyyy_yyyyyyy[k] * ab_x + g_0_y_xyyy_xyyyyyyy[k];

                g_0_y_xxyyy_yyyyyyz[k] = -g_0_y_xyyy_yyyyyyz[k] * ab_x + g_0_y_xyyy_xyyyyyyz[k];

                g_0_y_xxyyy_yyyyyzz[k] = -g_0_y_xyyy_yyyyyzz[k] * ab_x + g_0_y_xyyy_xyyyyyzz[k];

                g_0_y_xxyyy_yyyyzzz[k] = -g_0_y_xyyy_yyyyzzz[k] * ab_x + g_0_y_xyyy_xyyyyzzz[k];

                g_0_y_xxyyy_yyyzzzz[k] = -g_0_y_xyyy_yyyzzzz[k] * ab_x + g_0_y_xyyy_xyyyzzzz[k];

                g_0_y_xxyyy_yyzzzzz[k] = -g_0_y_xyyy_yyzzzzz[k] * ab_x + g_0_y_xyyy_xyyzzzzz[k];

                g_0_y_xxyyy_yzzzzzz[k] = -g_0_y_xyyy_yzzzzzz[k] * ab_x + g_0_y_xyyy_xyzzzzzz[k];

                g_0_y_xxyyy_zzzzzzz[k] = -g_0_y_xyyy_zzzzzzz[k] * ab_x + g_0_y_xyyy_xzzzzzzz[k];
            }

            /// Set up 1008-1044 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_y_xxyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_y_xxyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_y_xxyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1043 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyz_xxxxxxx, g_0_y_xxyyz_xxxxxxy, g_0_y_xxyyz_xxxxxxz, g_0_y_xxyyz_xxxxxyy, g_0_y_xxyyz_xxxxxyz, g_0_y_xxyyz_xxxxxzz, g_0_y_xxyyz_xxxxyyy, g_0_y_xxyyz_xxxxyyz, g_0_y_xxyyz_xxxxyzz, g_0_y_xxyyz_xxxxzzz, g_0_y_xxyyz_xxxyyyy, g_0_y_xxyyz_xxxyyyz, g_0_y_xxyyz_xxxyyzz, g_0_y_xxyyz_xxxyzzz, g_0_y_xxyyz_xxxzzzz, g_0_y_xxyyz_xxyyyyy, g_0_y_xxyyz_xxyyyyz, g_0_y_xxyyz_xxyyyzz, g_0_y_xxyyz_xxyyzzz, g_0_y_xxyyz_xxyzzzz, g_0_y_xxyyz_xxzzzzz, g_0_y_xxyyz_xyyyyyy, g_0_y_xxyyz_xyyyyyz, g_0_y_xxyyz_xyyyyzz, g_0_y_xxyyz_xyyyzzz, g_0_y_xxyyz_xyyzzzz, g_0_y_xxyyz_xyzzzzz, g_0_y_xxyyz_xzzzzzz, g_0_y_xxyyz_yyyyyyy, g_0_y_xxyyz_yyyyyyz, g_0_y_xxyyz_yyyyyzz, g_0_y_xxyyz_yyyyzzz, g_0_y_xxyyz_yyyzzzz, g_0_y_xxyyz_yyzzzzz, g_0_y_xxyyz_yzzzzzz, g_0_y_xxyyz_zzzzzzz, g_0_y_xyyz_xxxxxxx, g_0_y_xyyz_xxxxxxxx, g_0_y_xyyz_xxxxxxxy, g_0_y_xyyz_xxxxxxxz, g_0_y_xyyz_xxxxxxy, g_0_y_xyyz_xxxxxxyy, g_0_y_xyyz_xxxxxxyz, g_0_y_xyyz_xxxxxxz, g_0_y_xyyz_xxxxxxzz, g_0_y_xyyz_xxxxxyy, g_0_y_xyyz_xxxxxyyy, g_0_y_xyyz_xxxxxyyz, g_0_y_xyyz_xxxxxyz, g_0_y_xyyz_xxxxxyzz, g_0_y_xyyz_xxxxxzz, g_0_y_xyyz_xxxxxzzz, g_0_y_xyyz_xxxxyyy, g_0_y_xyyz_xxxxyyyy, g_0_y_xyyz_xxxxyyyz, g_0_y_xyyz_xxxxyyz, g_0_y_xyyz_xxxxyyzz, g_0_y_xyyz_xxxxyzz, g_0_y_xyyz_xxxxyzzz, g_0_y_xyyz_xxxxzzz, g_0_y_xyyz_xxxxzzzz, g_0_y_xyyz_xxxyyyy, g_0_y_xyyz_xxxyyyyy, g_0_y_xyyz_xxxyyyyz, g_0_y_xyyz_xxxyyyz, g_0_y_xyyz_xxxyyyzz, g_0_y_xyyz_xxxyyzz, g_0_y_xyyz_xxxyyzzz, g_0_y_xyyz_xxxyzzz, g_0_y_xyyz_xxxyzzzz, g_0_y_xyyz_xxxzzzz, g_0_y_xyyz_xxxzzzzz, g_0_y_xyyz_xxyyyyy, g_0_y_xyyz_xxyyyyyy, g_0_y_xyyz_xxyyyyyz, g_0_y_xyyz_xxyyyyz, g_0_y_xyyz_xxyyyyzz, g_0_y_xyyz_xxyyyzz, g_0_y_xyyz_xxyyyzzz, g_0_y_xyyz_xxyyzzz, g_0_y_xyyz_xxyyzzzz, g_0_y_xyyz_xxyzzzz, g_0_y_xyyz_xxyzzzzz, g_0_y_xyyz_xxzzzzz, g_0_y_xyyz_xxzzzzzz, g_0_y_xyyz_xyyyyyy, g_0_y_xyyz_xyyyyyyy, g_0_y_xyyz_xyyyyyyz, g_0_y_xyyz_xyyyyyz, g_0_y_xyyz_xyyyyyzz, g_0_y_xyyz_xyyyyzz, g_0_y_xyyz_xyyyyzzz, g_0_y_xyyz_xyyyzzz, g_0_y_xyyz_xyyyzzzz, g_0_y_xyyz_xyyzzzz, g_0_y_xyyz_xyyzzzzz, g_0_y_xyyz_xyzzzzz, g_0_y_xyyz_xyzzzzzz, g_0_y_xyyz_xzzzzzz, g_0_y_xyyz_xzzzzzzz, g_0_y_xyyz_yyyyyyy, g_0_y_xyyz_yyyyyyz, g_0_y_xyyz_yyyyyzz, g_0_y_xyyz_yyyyzzz, g_0_y_xyyz_yyyzzzz, g_0_y_xyyz_yyzzzzz, g_0_y_xyyz_yzzzzzz, g_0_y_xyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyz_xxxxxxx[k] = -g_0_y_xyyz_xxxxxxx[k] * ab_x + g_0_y_xyyz_xxxxxxxx[k];

                g_0_y_xxyyz_xxxxxxy[k] = -g_0_y_xyyz_xxxxxxy[k] * ab_x + g_0_y_xyyz_xxxxxxxy[k];

                g_0_y_xxyyz_xxxxxxz[k] = -g_0_y_xyyz_xxxxxxz[k] * ab_x + g_0_y_xyyz_xxxxxxxz[k];

                g_0_y_xxyyz_xxxxxyy[k] = -g_0_y_xyyz_xxxxxyy[k] * ab_x + g_0_y_xyyz_xxxxxxyy[k];

                g_0_y_xxyyz_xxxxxyz[k] = -g_0_y_xyyz_xxxxxyz[k] * ab_x + g_0_y_xyyz_xxxxxxyz[k];

                g_0_y_xxyyz_xxxxxzz[k] = -g_0_y_xyyz_xxxxxzz[k] * ab_x + g_0_y_xyyz_xxxxxxzz[k];

                g_0_y_xxyyz_xxxxyyy[k] = -g_0_y_xyyz_xxxxyyy[k] * ab_x + g_0_y_xyyz_xxxxxyyy[k];

                g_0_y_xxyyz_xxxxyyz[k] = -g_0_y_xyyz_xxxxyyz[k] * ab_x + g_0_y_xyyz_xxxxxyyz[k];

                g_0_y_xxyyz_xxxxyzz[k] = -g_0_y_xyyz_xxxxyzz[k] * ab_x + g_0_y_xyyz_xxxxxyzz[k];

                g_0_y_xxyyz_xxxxzzz[k] = -g_0_y_xyyz_xxxxzzz[k] * ab_x + g_0_y_xyyz_xxxxxzzz[k];

                g_0_y_xxyyz_xxxyyyy[k] = -g_0_y_xyyz_xxxyyyy[k] * ab_x + g_0_y_xyyz_xxxxyyyy[k];

                g_0_y_xxyyz_xxxyyyz[k] = -g_0_y_xyyz_xxxyyyz[k] * ab_x + g_0_y_xyyz_xxxxyyyz[k];

                g_0_y_xxyyz_xxxyyzz[k] = -g_0_y_xyyz_xxxyyzz[k] * ab_x + g_0_y_xyyz_xxxxyyzz[k];

                g_0_y_xxyyz_xxxyzzz[k] = -g_0_y_xyyz_xxxyzzz[k] * ab_x + g_0_y_xyyz_xxxxyzzz[k];

                g_0_y_xxyyz_xxxzzzz[k] = -g_0_y_xyyz_xxxzzzz[k] * ab_x + g_0_y_xyyz_xxxxzzzz[k];

                g_0_y_xxyyz_xxyyyyy[k] = -g_0_y_xyyz_xxyyyyy[k] * ab_x + g_0_y_xyyz_xxxyyyyy[k];

                g_0_y_xxyyz_xxyyyyz[k] = -g_0_y_xyyz_xxyyyyz[k] * ab_x + g_0_y_xyyz_xxxyyyyz[k];

                g_0_y_xxyyz_xxyyyzz[k] = -g_0_y_xyyz_xxyyyzz[k] * ab_x + g_0_y_xyyz_xxxyyyzz[k];

                g_0_y_xxyyz_xxyyzzz[k] = -g_0_y_xyyz_xxyyzzz[k] * ab_x + g_0_y_xyyz_xxxyyzzz[k];

                g_0_y_xxyyz_xxyzzzz[k] = -g_0_y_xyyz_xxyzzzz[k] * ab_x + g_0_y_xyyz_xxxyzzzz[k];

                g_0_y_xxyyz_xxzzzzz[k] = -g_0_y_xyyz_xxzzzzz[k] * ab_x + g_0_y_xyyz_xxxzzzzz[k];

                g_0_y_xxyyz_xyyyyyy[k] = -g_0_y_xyyz_xyyyyyy[k] * ab_x + g_0_y_xyyz_xxyyyyyy[k];

                g_0_y_xxyyz_xyyyyyz[k] = -g_0_y_xyyz_xyyyyyz[k] * ab_x + g_0_y_xyyz_xxyyyyyz[k];

                g_0_y_xxyyz_xyyyyzz[k] = -g_0_y_xyyz_xyyyyzz[k] * ab_x + g_0_y_xyyz_xxyyyyzz[k];

                g_0_y_xxyyz_xyyyzzz[k] = -g_0_y_xyyz_xyyyzzz[k] * ab_x + g_0_y_xyyz_xxyyyzzz[k];

                g_0_y_xxyyz_xyyzzzz[k] = -g_0_y_xyyz_xyyzzzz[k] * ab_x + g_0_y_xyyz_xxyyzzzz[k];

                g_0_y_xxyyz_xyzzzzz[k] = -g_0_y_xyyz_xyzzzzz[k] * ab_x + g_0_y_xyyz_xxyzzzzz[k];

                g_0_y_xxyyz_xzzzzzz[k] = -g_0_y_xyyz_xzzzzzz[k] * ab_x + g_0_y_xyyz_xxzzzzzz[k];

                g_0_y_xxyyz_yyyyyyy[k] = -g_0_y_xyyz_yyyyyyy[k] * ab_x + g_0_y_xyyz_xyyyyyyy[k];

                g_0_y_xxyyz_yyyyyyz[k] = -g_0_y_xyyz_yyyyyyz[k] * ab_x + g_0_y_xyyz_xyyyyyyz[k];

                g_0_y_xxyyz_yyyyyzz[k] = -g_0_y_xyyz_yyyyyzz[k] * ab_x + g_0_y_xyyz_xyyyyyzz[k];

                g_0_y_xxyyz_yyyyzzz[k] = -g_0_y_xyyz_yyyyzzz[k] * ab_x + g_0_y_xyyz_xyyyyzzz[k];

                g_0_y_xxyyz_yyyzzzz[k] = -g_0_y_xyyz_yyyzzzz[k] * ab_x + g_0_y_xyyz_xyyyzzzz[k];

                g_0_y_xxyyz_yyzzzzz[k] = -g_0_y_xyyz_yyzzzzz[k] * ab_x + g_0_y_xyyz_xyyzzzzz[k];

                g_0_y_xxyyz_yzzzzzz[k] = -g_0_y_xyyz_yzzzzzz[k] * ab_x + g_0_y_xyyz_xyzzzzzz[k];

                g_0_y_xxyyz_zzzzzzz[k] = -g_0_y_xyyz_zzzzzzz[k] * ab_x + g_0_y_xyyz_xzzzzzzz[k];
            }

            /// Set up 1044-1080 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_y_xxyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_y_xxyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_y_xxyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzz_xxxxxxx, g_0_y_xxyzz_xxxxxxy, g_0_y_xxyzz_xxxxxxz, g_0_y_xxyzz_xxxxxyy, g_0_y_xxyzz_xxxxxyz, g_0_y_xxyzz_xxxxxzz, g_0_y_xxyzz_xxxxyyy, g_0_y_xxyzz_xxxxyyz, g_0_y_xxyzz_xxxxyzz, g_0_y_xxyzz_xxxxzzz, g_0_y_xxyzz_xxxyyyy, g_0_y_xxyzz_xxxyyyz, g_0_y_xxyzz_xxxyyzz, g_0_y_xxyzz_xxxyzzz, g_0_y_xxyzz_xxxzzzz, g_0_y_xxyzz_xxyyyyy, g_0_y_xxyzz_xxyyyyz, g_0_y_xxyzz_xxyyyzz, g_0_y_xxyzz_xxyyzzz, g_0_y_xxyzz_xxyzzzz, g_0_y_xxyzz_xxzzzzz, g_0_y_xxyzz_xyyyyyy, g_0_y_xxyzz_xyyyyyz, g_0_y_xxyzz_xyyyyzz, g_0_y_xxyzz_xyyyzzz, g_0_y_xxyzz_xyyzzzz, g_0_y_xxyzz_xyzzzzz, g_0_y_xxyzz_xzzzzzz, g_0_y_xxyzz_yyyyyyy, g_0_y_xxyzz_yyyyyyz, g_0_y_xxyzz_yyyyyzz, g_0_y_xxyzz_yyyyzzz, g_0_y_xxyzz_yyyzzzz, g_0_y_xxyzz_yyzzzzz, g_0_y_xxyzz_yzzzzzz, g_0_y_xxyzz_zzzzzzz, g_0_y_xyzz_xxxxxxx, g_0_y_xyzz_xxxxxxxx, g_0_y_xyzz_xxxxxxxy, g_0_y_xyzz_xxxxxxxz, g_0_y_xyzz_xxxxxxy, g_0_y_xyzz_xxxxxxyy, g_0_y_xyzz_xxxxxxyz, g_0_y_xyzz_xxxxxxz, g_0_y_xyzz_xxxxxxzz, g_0_y_xyzz_xxxxxyy, g_0_y_xyzz_xxxxxyyy, g_0_y_xyzz_xxxxxyyz, g_0_y_xyzz_xxxxxyz, g_0_y_xyzz_xxxxxyzz, g_0_y_xyzz_xxxxxzz, g_0_y_xyzz_xxxxxzzz, g_0_y_xyzz_xxxxyyy, g_0_y_xyzz_xxxxyyyy, g_0_y_xyzz_xxxxyyyz, g_0_y_xyzz_xxxxyyz, g_0_y_xyzz_xxxxyyzz, g_0_y_xyzz_xxxxyzz, g_0_y_xyzz_xxxxyzzz, g_0_y_xyzz_xxxxzzz, g_0_y_xyzz_xxxxzzzz, g_0_y_xyzz_xxxyyyy, g_0_y_xyzz_xxxyyyyy, g_0_y_xyzz_xxxyyyyz, g_0_y_xyzz_xxxyyyz, g_0_y_xyzz_xxxyyyzz, g_0_y_xyzz_xxxyyzz, g_0_y_xyzz_xxxyyzzz, g_0_y_xyzz_xxxyzzz, g_0_y_xyzz_xxxyzzzz, g_0_y_xyzz_xxxzzzz, g_0_y_xyzz_xxxzzzzz, g_0_y_xyzz_xxyyyyy, g_0_y_xyzz_xxyyyyyy, g_0_y_xyzz_xxyyyyyz, g_0_y_xyzz_xxyyyyz, g_0_y_xyzz_xxyyyyzz, g_0_y_xyzz_xxyyyzz, g_0_y_xyzz_xxyyyzzz, g_0_y_xyzz_xxyyzzz, g_0_y_xyzz_xxyyzzzz, g_0_y_xyzz_xxyzzzz, g_0_y_xyzz_xxyzzzzz, g_0_y_xyzz_xxzzzzz, g_0_y_xyzz_xxzzzzzz, g_0_y_xyzz_xyyyyyy, g_0_y_xyzz_xyyyyyyy, g_0_y_xyzz_xyyyyyyz, g_0_y_xyzz_xyyyyyz, g_0_y_xyzz_xyyyyyzz, g_0_y_xyzz_xyyyyzz, g_0_y_xyzz_xyyyyzzz, g_0_y_xyzz_xyyyzzz, g_0_y_xyzz_xyyyzzzz, g_0_y_xyzz_xyyzzzz, g_0_y_xyzz_xyyzzzzz, g_0_y_xyzz_xyzzzzz, g_0_y_xyzz_xyzzzzzz, g_0_y_xyzz_xzzzzzz, g_0_y_xyzz_xzzzzzzz, g_0_y_xyzz_yyyyyyy, g_0_y_xyzz_yyyyyyz, g_0_y_xyzz_yyyyyzz, g_0_y_xyzz_yyyyzzz, g_0_y_xyzz_yyyzzzz, g_0_y_xyzz_yyzzzzz, g_0_y_xyzz_yzzzzzz, g_0_y_xyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzz_xxxxxxx[k] = -g_0_y_xyzz_xxxxxxx[k] * ab_x + g_0_y_xyzz_xxxxxxxx[k];

                g_0_y_xxyzz_xxxxxxy[k] = -g_0_y_xyzz_xxxxxxy[k] * ab_x + g_0_y_xyzz_xxxxxxxy[k];

                g_0_y_xxyzz_xxxxxxz[k] = -g_0_y_xyzz_xxxxxxz[k] * ab_x + g_0_y_xyzz_xxxxxxxz[k];

                g_0_y_xxyzz_xxxxxyy[k] = -g_0_y_xyzz_xxxxxyy[k] * ab_x + g_0_y_xyzz_xxxxxxyy[k];

                g_0_y_xxyzz_xxxxxyz[k] = -g_0_y_xyzz_xxxxxyz[k] * ab_x + g_0_y_xyzz_xxxxxxyz[k];

                g_0_y_xxyzz_xxxxxzz[k] = -g_0_y_xyzz_xxxxxzz[k] * ab_x + g_0_y_xyzz_xxxxxxzz[k];

                g_0_y_xxyzz_xxxxyyy[k] = -g_0_y_xyzz_xxxxyyy[k] * ab_x + g_0_y_xyzz_xxxxxyyy[k];

                g_0_y_xxyzz_xxxxyyz[k] = -g_0_y_xyzz_xxxxyyz[k] * ab_x + g_0_y_xyzz_xxxxxyyz[k];

                g_0_y_xxyzz_xxxxyzz[k] = -g_0_y_xyzz_xxxxyzz[k] * ab_x + g_0_y_xyzz_xxxxxyzz[k];

                g_0_y_xxyzz_xxxxzzz[k] = -g_0_y_xyzz_xxxxzzz[k] * ab_x + g_0_y_xyzz_xxxxxzzz[k];

                g_0_y_xxyzz_xxxyyyy[k] = -g_0_y_xyzz_xxxyyyy[k] * ab_x + g_0_y_xyzz_xxxxyyyy[k];

                g_0_y_xxyzz_xxxyyyz[k] = -g_0_y_xyzz_xxxyyyz[k] * ab_x + g_0_y_xyzz_xxxxyyyz[k];

                g_0_y_xxyzz_xxxyyzz[k] = -g_0_y_xyzz_xxxyyzz[k] * ab_x + g_0_y_xyzz_xxxxyyzz[k];

                g_0_y_xxyzz_xxxyzzz[k] = -g_0_y_xyzz_xxxyzzz[k] * ab_x + g_0_y_xyzz_xxxxyzzz[k];

                g_0_y_xxyzz_xxxzzzz[k] = -g_0_y_xyzz_xxxzzzz[k] * ab_x + g_0_y_xyzz_xxxxzzzz[k];

                g_0_y_xxyzz_xxyyyyy[k] = -g_0_y_xyzz_xxyyyyy[k] * ab_x + g_0_y_xyzz_xxxyyyyy[k];

                g_0_y_xxyzz_xxyyyyz[k] = -g_0_y_xyzz_xxyyyyz[k] * ab_x + g_0_y_xyzz_xxxyyyyz[k];

                g_0_y_xxyzz_xxyyyzz[k] = -g_0_y_xyzz_xxyyyzz[k] * ab_x + g_0_y_xyzz_xxxyyyzz[k];

                g_0_y_xxyzz_xxyyzzz[k] = -g_0_y_xyzz_xxyyzzz[k] * ab_x + g_0_y_xyzz_xxxyyzzz[k];

                g_0_y_xxyzz_xxyzzzz[k] = -g_0_y_xyzz_xxyzzzz[k] * ab_x + g_0_y_xyzz_xxxyzzzz[k];

                g_0_y_xxyzz_xxzzzzz[k] = -g_0_y_xyzz_xxzzzzz[k] * ab_x + g_0_y_xyzz_xxxzzzzz[k];

                g_0_y_xxyzz_xyyyyyy[k] = -g_0_y_xyzz_xyyyyyy[k] * ab_x + g_0_y_xyzz_xxyyyyyy[k];

                g_0_y_xxyzz_xyyyyyz[k] = -g_0_y_xyzz_xyyyyyz[k] * ab_x + g_0_y_xyzz_xxyyyyyz[k];

                g_0_y_xxyzz_xyyyyzz[k] = -g_0_y_xyzz_xyyyyzz[k] * ab_x + g_0_y_xyzz_xxyyyyzz[k];

                g_0_y_xxyzz_xyyyzzz[k] = -g_0_y_xyzz_xyyyzzz[k] * ab_x + g_0_y_xyzz_xxyyyzzz[k];

                g_0_y_xxyzz_xyyzzzz[k] = -g_0_y_xyzz_xyyzzzz[k] * ab_x + g_0_y_xyzz_xxyyzzzz[k];

                g_0_y_xxyzz_xyzzzzz[k] = -g_0_y_xyzz_xyzzzzz[k] * ab_x + g_0_y_xyzz_xxyzzzzz[k];

                g_0_y_xxyzz_xzzzzzz[k] = -g_0_y_xyzz_xzzzzzz[k] * ab_x + g_0_y_xyzz_xxzzzzzz[k];

                g_0_y_xxyzz_yyyyyyy[k] = -g_0_y_xyzz_yyyyyyy[k] * ab_x + g_0_y_xyzz_xyyyyyyy[k];

                g_0_y_xxyzz_yyyyyyz[k] = -g_0_y_xyzz_yyyyyyz[k] * ab_x + g_0_y_xyzz_xyyyyyyz[k];

                g_0_y_xxyzz_yyyyyzz[k] = -g_0_y_xyzz_yyyyyzz[k] * ab_x + g_0_y_xyzz_xyyyyyzz[k];

                g_0_y_xxyzz_yyyyzzz[k] = -g_0_y_xyzz_yyyyzzz[k] * ab_x + g_0_y_xyzz_xyyyyzzz[k];

                g_0_y_xxyzz_yyyzzzz[k] = -g_0_y_xyzz_yyyzzzz[k] * ab_x + g_0_y_xyzz_xyyyzzzz[k];

                g_0_y_xxyzz_yyzzzzz[k] = -g_0_y_xyzz_yyzzzzz[k] * ab_x + g_0_y_xyzz_xyyzzzzz[k];

                g_0_y_xxyzz_yzzzzzz[k] = -g_0_y_xyzz_yzzzzzz[k] * ab_x + g_0_y_xyzz_xyzzzzzz[k];

                g_0_y_xxyzz_zzzzzzz[k] = -g_0_y_xyzz_zzzzzzz[k] * ab_x + g_0_y_xyzz_xzzzzzzz[k];
            }

            /// Set up 1080-1116 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_y_xxzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_y_xxzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_y_xxzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1115 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzz_xxxxxxx, g_0_y_xxzzz_xxxxxxy, g_0_y_xxzzz_xxxxxxz, g_0_y_xxzzz_xxxxxyy, g_0_y_xxzzz_xxxxxyz, g_0_y_xxzzz_xxxxxzz, g_0_y_xxzzz_xxxxyyy, g_0_y_xxzzz_xxxxyyz, g_0_y_xxzzz_xxxxyzz, g_0_y_xxzzz_xxxxzzz, g_0_y_xxzzz_xxxyyyy, g_0_y_xxzzz_xxxyyyz, g_0_y_xxzzz_xxxyyzz, g_0_y_xxzzz_xxxyzzz, g_0_y_xxzzz_xxxzzzz, g_0_y_xxzzz_xxyyyyy, g_0_y_xxzzz_xxyyyyz, g_0_y_xxzzz_xxyyyzz, g_0_y_xxzzz_xxyyzzz, g_0_y_xxzzz_xxyzzzz, g_0_y_xxzzz_xxzzzzz, g_0_y_xxzzz_xyyyyyy, g_0_y_xxzzz_xyyyyyz, g_0_y_xxzzz_xyyyyzz, g_0_y_xxzzz_xyyyzzz, g_0_y_xxzzz_xyyzzzz, g_0_y_xxzzz_xyzzzzz, g_0_y_xxzzz_xzzzzzz, g_0_y_xxzzz_yyyyyyy, g_0_y_xxzzz_yyyyyyz, g_0_y_xxzzz_yyyyyzz, g_0_y_xxzzz_yyyyzzz, g_0_y_xxzzz_yyyzzzz, g_0_y_xxzzz_yyzzzzz, g_0_y_xxzzz_yzzzzzz, g_0_y_xxzzz_zzzzzzz, g_0_y_xzzz_xxxxxxx, g_0_y_xzzz_xxxxxxxx, g_0_y_xzzz_xxxxxxxy, g_0_y_xzzz_xxxxxxxz, g_0_y_xzzz_xxxxxxy, g_0_y_xzzz_xxxxxxyy, g_0_y_xzzz_xxxxxxyz, g_0_y_xzzz_xxxxxxz, g_0_y_xzzz_xxxxxxzz, g_0_y_xzzz_xxxxxyy, g_0_y_xzzz_xxxxxyyy, g_0_y_xzzz_xxxxxyyz, g_0_y_xzzz_xxxxxyz, g_0_y_xzzz_xxxxxyzz, g_0_y_xzzz_xxxxxzz, g_0_y_xzzz_xxxxxzzz, g_0_y_xzzz_xxxxyyy, g_0_y_xzzz_xxxxyyyy, g_0_y_xzzz_xxxxyyyz, g_0_y_xzzz_xxxxyyz, g_0_y_xzzz_xxxxyyzz, g_0_y_xzzz_xxxxyzz, g_0_y_xzzz_xxxxyzzz, g_0_y_xzzz_xxxxzzz, g_0_y_xzzz_xxxxzzzz, g_0_y_xzzz_xxxyyyy, g_0_y_xzzz_xxxyyyyy, g_0_y_xzzz_xxxyyyyz, g_0_y_xzzz_xxxyyyz, g_0_y_xzzz_xxxyyyzz, g_0_y_xzzz_xxxyyzz, g_0_y_xzzz_xxxyyzzz, g_0_y_xzzz_xxxyzzz, g_0_y_xzzz_xxxyzzzz, g_0_y_xzzz_xxxzzzz, g_0_y_xzzz_xxxzzzzz, g_0_y_xzzz_xxyyyyy, g_0_y_xzzz_xxyyyyyy, g_0_y_xzzz_xxyyyyyz, g_0_y_xzzz_xxyyyyz, g_0_y_xzzz_xxyyyyzz, g_0_y_xzzz_xxyyyzz, g_0_y_xzzz_xxyyyzzz, g_0_y_xzzz_xxyyzzz, g_0_y_xzzz_xxyyzzzz, g_0_y_xzzz_xxyzzzz, g_0_y_xzzz_xxyzzzzz, g_0_y_xzzz_xxzzzzz, g_0_y_xzzz_xxzzzzzz, g_0_y_xzzz_xyyyyyy, g_0_y_xzzz_xyyyyyyy, g_0_y_xzzz_xyyyyyyz, g_0_y_xzzz_xyyyyyz, g_0_y_xzzz_xyyyyyzz, g_0_y_xzzz_xyyyyzz, g_0_y_xzzz_xyyyyzzz, g_0_y_xzzz_xyyyzzz, g_0_y_xzzz_xyyyzzzz, g_0_y_xzzz_xyyzzzz, g_0_y_xzzz_xyyzzzzz, g_0_y_xzzz_xyzzzzz, g_0_y_xzzz_xyzzzzzz, g_0_y_xzzz_xzzzzzz, g_0_y_xzzz_xzzzzzzz, g_0_y_xzzz_yyyyyyy, g_0_y_xzzz_yyyyyyz, g_0_y_xzzz_yyyyyzz, g_0_y_xzzz_yyyyzzz, g_0_y_xzzz_yyyzzzz, g_0_y_xzzz_yyzzzzz, g_0_y_xzzz_yzzzzzz, g_0_y_xzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzz_xxxxxxx[k] = -g_0_y_xzzz_xxxxxxx[k] * ab_x + g_0_y_xzzz_xxxxxxxx[k];

                g_0_y_xxzzz_xxxxxxy[k] = -g_0_y_xzzz_xxxxxxy[k] * ab_x + g_0_y_xzzz_xxxxxxxy[k];

                g_0_y_xxzzz_xxxxxxz[k] = -g_0_y_xzzz_xxxxxxz[k] * ab_x + g_0_y_xzzz_xxxxxxxz[k];

                g_0_y_xxzzz_xxxxxyy[k] = -g_0_y_xzzz_xxxxxyy[k] * ab_x + g_0_y_xzzz_xxxxxxyy[k];

                g_0_y_xxzzz_xxxxxyz[k] = -g_0_y_xzzz_xxxxxyz[k] * ab_x + g_0_y_xzzz_xxxxxxyz[k];

                g_0_y_xxzzz_xxxxxzz[k] = -g_0_y_xzzz_xxxxxzz[k] * ab_x + g_0_y_xzzz_xxxxxxzz[k];

                g_0_y_xxzzz_xxxxyyy[k] = -g_0_y_xzzz_xxxxyyy[k] * ab_x + g_0_y_xzzz_xxxxxyyy[k];

                g_0_y_xxzzz_xxxxyyz[k] = -g_0_y_xzzz_xxxxyyz[k] * ab_x + g_0_y_xzzz_xxxxxyyz[k];

                g_0_y_xxzzz_xxxxyzz[k] = -g_0_y_xzzz_xxxxyzz[k] * ab_x + g_0_y_xzzz_xxxxxyzz[k];

                g_0_y_xxzzz_xxxxzzz[k] = -g_0_y_xzzz_xxxxzzz[k] * ab_x + g_0_y_xzzz_xxxxxzzz[k];

                g_0_y_xxzzz_xxxyyyy[k] = -g_0_y_xzzz_xxxyyyy[k] * ab_x + g_0_y_xzzz_xxxxyyyy[k];

                g_0_y_xxzzz_xxxyyyz[k] = -g_0_y_xzzz_xxxyyyz[k] * ab_x + g_0_y_xzzz_xxxxyyyz[k];

                g_0_y_xxzzz_xxxyyzz[k] = -g_0_y_xzzz_xxxyyzz[k] * ab_x + g_0_y_xzzz_xxxxyyzz[k];

                g_0_y_xxzzz_xxxyzzz[k] = -g_0_y_xzzz_xxxyzzz[k] * ab_x + g_0_y_xzzz_xxxxyzzz[k];

                g_0_y_xxzzz_xxxzzzz[k] = -g_0_y_xzzz_xxxzzzz[k] * ab_x + g_0_y_xzzz_xxxxzzzz[k];

                g_0_y_xxzzz_xxyyyyy[k] = -g_0_y_xzzz_xxyyyyy[k] * ab_x + g_0_y_xzzz_xxxyyyyy[k];

                g_0_y_xxzzz_xxyyyyz[k] = -g_0_y_xzzz_xxyyyyz[k] * ab_x + g_0_y_xzzz_xxxyyyyz[k];

                g_0_y_xxzzz_xxyyyzz[k] = -g_0_y_xzzz_xxyyyzz[k] * ab_x + g_0_y_xzzz_xxxyyyzz[k];

                g_0_y_xxzzz_xxyyzzz[k] = -g_0_y_xzzz_xxyyzzz[k] * ab_x + g_0_y_xzzz_xxxyyzzz[k];

                g_0_y_xxzzz_xxyzzzz[k] = -g_0_y_xzzz_xxyzzzz[k] * ab_x + g_0_y_xzzz_xxxyzzzz[k];

                g_0_y_xxzzz_xxzzzzz[k] = -g_0_y_xzzz_xxzzzzz[k] * ab_x + g_0_y_xzzz_xxxzzzzz[k];

                g_0_y_xxzzz_xyyyyyy[k] = -g_0_y_xzzz_xyyyyyy[k] * ab_x + g_0_y_xzzz_xxyyyyyy[k];

                g_0_y_xxzzz_xyyyyyz[k] = -g_0_y_xzzz_xyyyyyz[k] * ab_x + g_0_y_xzzz_xxyyyyyz[k];

                g_0_y_xxzzz_xyyyyzz[k] = -g_0_y_xzzz_xyyyyzz[k] * ab_x + g_0_y_xzzz_xxyyyyzz[k];

                g_0_y_xxzzz_xyyyzzz[k] = -g_0_y_xzzz_xyyyzzz[k] * ab_x + g_0_y_xzzz_xxyyyzzz[k];

                g_0_y_xxzzz_xyyzzzz[k] = -g_0_y_xzzz_xyyzzzz[k] * ab_x + g_0_y_xzzz_xxyyzzzz[k];

                g_0_y_xxzzz_xyzzzzz[k] = -g_0_y_xzzz_xyzzzzz[k] * ab_x + g_0_y_xzzz_xxyzzzzz[k];

                g_0_y_xxzzz_xzzzzzz[k] = -g_0_y_xzzz_xzzzzzz[k] * ab_x + g_0_y_xzzz_xxzzzzzz[k];

                g_0_y_xxzzz_yyyyyyy[k] = -g_0_y_xzzz_yyyyyyy[k] * ab_x + g_0_y_xzzz_xyyyyyyy[k];

                g_0_y_xxzzz_yyyyyyz[k] = -g_0_y_xzzz_yyyyyyz[k] * ab_x + g_0_y_xzzz_xyyyyyyz[k];

                g_0_y_xxzzz_yyyyyzz[k] = -g_0_y_xzzz_yyyyyzz[k] * ab_x + g_0_y_xzzz_xyyyyyzz[k];

                g_0_y_xxzzz_yyyyzzz[k] = -g_0_y_xzzz_yyyyzzz[k] * ab_x + g_0_y_xzzz_xyyyyzzz[k];

                g_0_y_xxzzz_yyyzzzz[k] = -g_0_y_xzzz_yyyzzzz[k] * ab_x + g_0_y_xzzz_xyyyzzzz[k];

                g_0_y_xxzzz_yyzzzzz[k] = -g_0_y_xzzz_yyzzzzz[k] * ab_x + g_0_y_xzzz_xyyzzzzz[k];

                g_0_y_xxzzz_yzzzzzz[k] = -g_0_y_xzzz_yzzzzzz[k] * ab_x + g_0_y_xzzz_xyzzzzzz[k];

                g_0_y_xxzzz_zzzzzzz[k] = -g_0_y_xzzz_zzzzzzz[k] * ab_x + g_0_y_xzzz_xzzzzzzz[k];
            }

            /// Set up 1116-1152 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_y_xyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_y_xyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_y_xyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1151 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyy_xxxxxxx, g_0_y_xyyyy_xxxxxxy, g_0_y_xyyyy_xxxxxxz, g_0_y_xyyyy_xxxxxyy, g_0_y_xyyyy_xxxxxyz, g_0_y_xyyyy_xxxxxzz, g_0_y_xyyyy_xxxxyyy, g_0_y_xyyyy_xxxxyyz, g_0_y_xyyyy_xxxxyzz, g_0_y_xyyyy_xxxxzzz, g_0_y_xyyyy_xxxyyyy, g_0_y_xyyyy_xxxyyyz, g_0_y_xyyyy_xxxyyzz, g_0_y_xyyyy_xxxyzzz, g_0_y_xyyyy_xxxzzzz, g_0_y_xyyyy_xxyyyyy, g_0_y_xyyyy_xxyyyyz, g_0_y_xyyyy_xxyyyzz, g_0_y_xyyyy_xxyyzzz, g_0_y_xyyyy_xxyzzzz, g_0_y_xyyyy_xxzzzzz, g_0_y_xyyyy_xyyyyyy, g_0_y_xyyyy_xyyyyyz, g_0_y_xyyyy_xyyyyzz, g_0_y_xyyyy_xyyyzzz, g_0_y_xyyyy_xyyzzzz, g_0_y_xyyyy_xyzzzzz, g_0_y_xyyyy_xzzzzzz, g_0_y_xyyyy_yyyyyyy, g_0_y_xyyyy_yyyyyyz, g_0_y_xyyyy_yyyyyzz, g_0_y_xyyyy_yyyyzzz, g_0_y_xyyyy_yyyzzzz, g_0_y_xyyyy_yyzzzzz, g_0_y_xyyyy_yzzzzzz, g_0_y_xyyyy_zzzzzzz, g_0_y_yyyy_xxxxxxx, g_0_y_yyyy_xxxxxxxx, g_0_y_yyyy_xxxxxxxy, g_0_y_yyyy_xxxxxxxz, g_0_y_yyyy_xxxxxxy, g_0_y_yyyy_xxxxxxyy, g_0_y_yyyy_xxxxxxyz, g_0_y_yyyy_xxxxxxz, g_0_y_yyyy_xxxxxxzz, g_0_y_yyyy_xxxxxyy, g_0_y_yyyy_xxxxxyyy, g_0_y_yyyy_xxxxxyyz, g_0_y_yyyy_xxxxxyz, g_0_y_yyyy_xxxxxyzz, g_0_y_yyyy_xxxxxzz, g_0_y_yyyy_xxxxxzzz, g_0_y_yyyy_xxxxyyy, g_0_y_yyyy_xxxxyyyy, g_0_y_yyyy_xxxxyyyz, g_0_y_yyyy_xxxxyyz, g_0_y_yyyy_xxxxyyzz, g_0_y_yyyy_xxxxyzz, g_0_y_yyyy_xxxxyzzz, g_0_y_yyyy_xxxxzzz, g_0_y_yyyy_xxxxzzzz, g_0_y_yyyy_xxxyyyy, g_0_y_yyyy_xxxyyyyy, g_0_y_yyyy_xxxyyyyz, g_0_y_yyyy_xxxyyyz, g_0_y_yyyy_xxxyyyzz, g_0_y_yyyy_xxxyyzz, g_0_y_yyyy_xxxyyzzz, g_0_y_yyyy_xxxyzzz, g_0_y_yyyy_xxxyzzzz, g_0_y_yyyy_xxxzzzz, g_0_y_yyyy_xxxzzzzz, g_0_y_yyyy_xxyyyyy, g_0_y_yyyy_xxyyyyyy, g_0_y_yyyy_xxyyyyyz, g_0_y_yyyy_xxyyyyz, g_0_y_yyyy_xxyyyyzz, g_0_y_yyyy_xxyyyzz, g_0_y_yyyy_xxyyyzzz, g_0_y_yyyy_xxyyzzz, g_0_y_yyyy_xxyyzzzz, g_0_y_yyyy_xxyzzzz, g_0_y_yyyy_xxyzzzzz, g_0_y_yyyy_xxzzzzz, g_0_y_yyyy_xxzzzzzz, g_0_y_yyyy_xyyyyyy, g_0_y_yyyy_xyyyyyyy, g_0_y_yyyy_xyyyyyyz, g_0_y_yyyy_xyyyyyz, g_0_y_yyyy_xyyyyyzz, g_0_y_yyyy_xyyyyzz, g_0_y_yyyy_xyyyyzzz, g_0_y_yyyy_xyyyzzz, g_0_y_yyyy_xyyyzzzz, g_0_y_yyyy_xyyzzzz, g_0_y_yyyy_xyyzzzzz, g_0_y_yyyy_xyzzzzz, g_0_y_yyyy_xyzzzzzz, g_0_y_yyyy_xzzzzzz, g_0_y_yyyy_xzzzzzzz, g_0_y_yyyy_yyyyyyy, g_0_y_yyyy_yyyyyyz, g_0_y_yyyy_yyyyyzz, g_0_y_yyyy_yyyyzzz, g_0_y_yyyy_yyyzzzz, g_0_y_yyyy_yyzzzzz, g_0_y_yyyy_yzzzzzz, g_0_y_yyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyy_xxxxxxx[k] = -g_0_y_yyyy_xxxxxxx[k] * ab_x + g_0_y_yyyy_xxxxxxxx[k];

                g_0_y_xyyyy_xxxxxxy[k] = -g_0_y_yyyy_xxxxxxy[k] * ab_x + g_0_y_yyyy_xxxxxxxy[k];

                g_0_y_xyyyy_xxxxxxz[k] = -g_0_y_yyyy_xxxxxxz[k] * ab_x + g_0_y_yyyy_xxxxxxxz[k];

                g_0_y_xyyyy_xxxxxyy[k] = -g_0_y_yyyy_xxxxxyy[k] * ab_x + g_0_y_yyyy_xxxxxxyy[k];

                g_0_y_xyyyy_xxxxxyz[k] = -g_0_y_yyyy_xxxxxyz[k] * ab_x + g_0_y_yyyy_xxxxxxyz[k];

                g_0_y_xyyyy_xxxxxzz[k] = -g_0_y_yyyy_xxxxxzz[k] * ab_x + g_0_y_yyyy_xxxxxxzz[k];

                g_0_y_xyyyy_xxxxyyy[k] = -g_0_y_yyyy_xxxxyyy[k] * ab_x + g_0_y_yyyy_xxxxxyyy[k];

                g_0_y_xyyyy_xxxxyyz[k] = -g_0_y_yyyy_xxxxyyz[k] * ab_x + g_0_y_yyyy_xxxxxyyz[k];

                g_0_y_xyyyy_xxxxyzz[k] = -g_0_y_yyyy_xxxxyzz[k] * ab_x + g_0_y_yyyy_xxxxxyzz[k];

                g_0_y_xyyyy_xxxxzzz[k] = -g_0_y_yyyy_xxxxzzz[k] * ab_x + g_0_y_yyyy_xxxxxzzz[k];

                g_0_y_xyyyy_xxxyyyy[k] = -g_0_y_yyyy_xxxyyyy[k] * ab_x + g_0_y_yyyy_xxxxyyyy[k];

                g_0_y_xyyyy_xxxyyyz[k] = -g_0_y_yyyy_xxxyyyz[k] * ab_x + g_0_y_yyyy_xxxxyyyz[k];

                g_0_y_xyyyy_xxxyyzz[k] = -g_0_y_yyyy_xxxyyzz[k] * ab_x + g_0_y_yyyy_xxxxyyzz[k];

                g_0_y_xyyyy_xxxyzzz[k] = -g_0_y_yyyy_xxxyzzz[k] * ab_x + g_0_y_yyyy_xxxxyzzz[k];

                g_0_y_xyyyy_xxxzzzz[k] = -g_0_y_yyyy_xxxzzzz[k] * ab_x + g_0_y_yyyy_xxxxzzzz[k];

                g_0_y_xyyyy_xxyyyyy[k] = -g_0_y_yyyy_xxyyyyy[k] * ab_x + g_0_y_yyyy_xxxyyyyy[k];

                g_0_y_xyyyy_xxyyyyz[k] = -g_0_y_yyyy_xxyyyyz[k] * ab_x + g_0_y_yyyy_xxxyyyyz[k];

                g_0_y_xyyyy_xxyyyzz[k] = -g_0_y_yyyy_xxyyyzz[k] * ab_x + g_0_y_yyyy_xxxyyyzz[k];

                g_0_y_xyyyy_xxyyzzz[k] = -g_0_y_yyyy_xxyyzzz[k] * ab_x + g_0_y_yyyy_xxxyyzzz[k];

                g_0_y_xyyyy_xxyzzzz[k] = -g_0_y_yyyy_xxyzzzz[k] * ab_x + g_0_y_yyyy_xxxyzzzz[k];

                g_0_y_xyyyy_xxzzzzz[k] = -g_0_y_yyyy_xxzzzzz[k] * ab_x + g_0_y_yyyy_xxxzzzzz[k];

                g_0_y_xyyyy_xyyyyyy[k] = -g_0_y_yyyy_xyyyyyy[k] * ab_x + g_0_y_yyyy_xxyyyyyy[k];

                g_0_y_xyyyy_xyyyyyz[k] = -g_0_y_yyyy_xyyyyyz[k] * ab_x + g_0_y_yyyy_xxyyyyyz[k];

                g_0_y_xyyyy_xyyyyzz[k] = -g_0_y_yyyy_xyyyyzz[k] * ab_x + g_0_y_yyyy_xxyyyyzz[k];

                g_0_y_xyyyy_xyyyzzz[k] = -g_0_y_yyyy_xyyyzzz[k] * ab_x + g_0_y_yyyy_xxyyyzzz[k];

                g_0_y_xyyyy_xyyzzzz[k] = -g_0_y_yyyy_xyyzzzz[k] * ab_x + g_0_y_yyyy_xxyyzzzz[k];

                g_0_y_xyyyy_xyzzzzz[k] = -g_0_y_yyyy_xyzzzzz[k] * ab_x + g_0_y_yyyy_xxyzzzzz[k];

                g_0_y_xyyyy_xzzzzzz[k] = -g_0_y_yyyy_xzzzzzz[k] * ab_x + g_0_y_yyyy_xxzzzzzz[k];

                g_0_y_xyyyy_yyyyyyy[k] = -g_0_y_yyyy_yyyyyyy[k] * ab_x + g_0_y_yyyy_xyyyyyyy[k];

                g_0_y_xyyyy_yyyyyyz[k] = -g_0_y_yyyy_yyyyyyz[k] * ab_x + g_0_y_yyyy_xyyyyyyz[k];

                g_0_y_xyyyy_yyyyyzz[k] = -g_0_y_yyyy_yyyyyzz[k] * ab_x + g_0_y_yyyy_xyyyyyzz[k];

                g_0_y_xyyyy_yyyyzzz[k] = -g_0_y_yyyy_yyyyzzz[k] * ab_x + g_0_y_yyyy_xyyyyzzz[k];

                g_0_y_xyyyy_yyyzzzz[k] = -g_0_y_yyyy_yyyzzzz[k] * ab_x + g_0_y_yyyy_xyyyzzzz[k];

                g_0_y_xyyyy_yyzzzzz[k] = -g_0_y_yyyy_yyzzzzz[k] * ab_x + g_0_y_yyyy_xyyzzzzz[k];

                g_0_y_xyyyy_yzzzzzz[k] = -g_0_y_yyyy_yzzzzzz[k] * ab_x + g_0_y_yyyy_xyzzzzzz[k];

                g_0_y_xyyyy_zzzzzzz[k] = -g_0_y_yyyy_zzzzzzz[k] * ab_x + g_0_y_yyyy_xzzzzzzz[k];
            }

            /// Set up 1152-1188 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_y_xyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_y_xyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_y_xyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1187 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyz_xxxxxxx, g_0_y_xyyyz_xxxxxxy, g_0_y_xyyyz_xxxxxxz, g_0_y_xyyyz_xxxxxyy, g_0_y_xyyyz_xxxxxyz, g_0_y_xyyyz_xxxxxzz, g_0_y_xyyyz_xxxxyyy, g_0_y_xyyyz_xxxxyyz, g_0_y_xyyyz_xxxxyzz, g_0_y_xyyyz_xxxxzzz, g_0_y_xyyyz_xxxyyyy, g_0_y_xyyyz_xxxyyyz, g_0_y_xyyyz_xxxyyzz, g_0_y_xyyyz_xxxyzzz, g_0_y_xyyyz_xxxzzzz, g_0_y_xyyyz_xxyyyyy, g_0_y_xyyyz_xxyyyyz, g_0_y_xyyyz_xxyyyzz, g_0_y_xyyyz_xxyyzzz, g_0_y_xyyyz_xxyzzzz, g_0_y_xyyyz_xxzzzzz, g_0_y_xyyyz_xyyyyyy, g_0_y_xyyyz_xyyyyyz, g_0_y_xyyyz_xyyyyzz, g_0_y_xyyyz_xyyyzzz, g_0_y_xyyyz_xyyzzzz, g_0_y_xyyyz_xyzzzzz, g_0_y_xyyyz_xzzzzzz, g_0_y_xyyyz_yyyyyyy, g_0_y_xyyyz_yyyyyyz, g_0_y_xyyyz_yyyyyzz, g_0_y_xyyyz_yyyyzzz, g_0_y_xyyyz_yyyzzzz, g_0_y_xyyyz_yyzzzzz, g_0_y_xyyyz_yzzzzzz, g_0_y_xyyyz_zzzzzzz, g_0_y_yyyz_xxxxxxx, g_0_y_yyyz_xxxxxxxx, g_0_y_yyyz_xxxxxxxy, g_0_y_yyyz_xxxxxxxz, g_0_y_yyyz_xxxxxxy, g_0_y_yyyz_xxxxxxyy, g_0_y_yyyz_xxxxxxyz, g_0_y_yyyz_xxxxxxz, g_0_y_yyyz_xxxxxxzz, g_0_y_yyyz_xxxxxyy, g_0_y_yyyz_xxxxxyyy, g_0_y_yyyz_xxxxxyyz, g_0_y_yyyz_xxxxxyz, g_0_y_yyyz_xxxxxyzz, g_0_y_yyyz_xxxxxzz, g_0_y_yyyz_xxxxxzzz, g_0_y_yyyz_xxxxyyy, g_0_y_yyyz_xxxxyyyy, g_0_y_yyyz_xxxxyyyz, g_0_y_yyyz_xxxxyyz, g_0_y_yyyz_xxxxyyzz, g_0_y_yyyz_xxxxyzz, g_0_y_yyyz_xxxxyzzz, g_0_y_yyyz_xxxxzzz, g_0_y_yyyz_xxxxzzzz, g_0_y_yyyz_xxxyyyy, g_0_y_yyyz_xxxyyyyy, g_0_y_yyyz_xxxyyyyz, g_0_y_yyyz_xxxyyyz, g_0_y_yyyz_xxxyyyzz, g_0_y_yyyz_xxxyyzz, g_0_y_yyyz_xxxyyzzz, g_0_y_yyyz_xxxyzzz, g_0_y_yyyz_xxxyzzzz, g_0_y_yyyz_xxxzzzz, g_0_y_yyyz_xxxzzzzz, g_0_y_yyyz_xxyyyyy, g_0_y_yyyz_xxyyyyyy, g_0_y_yyyz_xxyyyyyz, g_0_y_yyyz_xxyyyyz, g_0_y_yyyz_xxyyyyzz, g_0_y_yyyz_xxyyyzz, g_0_y_yyyz_xxyyyzzz, g_0_y_yyyz_xxyyzzz, g_0_y_yyyz_xxyyzzzz, g_0_y_yyyz_xxyzzzz, g_0_y_yyyz_xxyzzzzz, g_0_y_yyyz_xxzzzzz, g_0_y_yyyz_xxzzzzzz, g_0_y_yyyz_xyyyyyy, g_0_y_yyyz_xyyyyyyy, g_0_y_yyyz_xyyyyyyz, g_0_y_yyyz_xyyyyyz, g_0_y_yyyz_xyyyyyzz, g_0_y_yyyz_xyyyyzz, g_0_y_yyyz_xyyyyzzz, g_0_y_yyyz_xyyyzzz, g_0_y_yyyz_xyyyzzzz, g_0_y_yyyz_xyyzzzz, g_0_y_yyyz_xyyzzzzz, g_0_y_yyyz_xyzzzzz, g_0_y_yyyz_xyzzzzzz, g_0_y_yyyz_xzzzzzz, g_0_y_yyyz_xzzzzzzz, g_0_y_yyyz_yyyyyyy, g_0_y_yyyz_yyyyyyz, g_0_y_yyyz_yyyyyzz, g_0_y_yyyz_yyyyzzz, g_0_y_yyyz_yyyzzzz, g_0_y_yyyz_yyzzzzz, g_0_y_yyyz_yzzzzzz, g_0_y_yyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyz_xxxxxxx[k] = -g_0_y_yyyz_xxxxxxx[k] * ab_x + g_0_y_yyyz_xxxxxxxx[k];

                g_0_y_xyyyz_xxxxxxy[k] = -g_0_y_yyyz_xxxxxxy[k] * ab_x + g_0_y_yyyz_xxxxxxxy[k];

                g_0_y_xyyyz_xxxxxxz[k] = -g_0_y_yyyz_xxxxxxz[k] * ab_x + g_0_y_yyyz_xxxxxxxz[k];

                g_0_y_xyyyz_xxxxxyy[k] = -g_0_y_yyyz_xxxxxyy[k] * ab_x + g_0_y_yyyz_xxxxxxyy[k];

                g_0_y_xyyyz_xxxxxyz[k] = -g_0_y_yyyz_xxxxxyz[k] * ab_x + g_0_y_yyyz_xxxxxxyz[k];

                g_0_y_xyyyz_xxxxxzz[k] = -g_0_y_yyyz_xxxxxzz[k] * ab_x + g_0_y_yyyz_xxxxxxzz[k];

                g_0_y_xyyyz_xxxxyyy[k] = -g_0_y_yyyz_xxxxyyy[k] * ab_x + g_0_y_yyyz_xxxxxyyy[k];

                g_0_y_xyyyz_xxxxyyz[k] = -g_0_y_yyyz_xxxxyyz[k] * ab_x + g_0_y_yyyz_xxxxxyyz[k];

                g_0_y_xyyyz_xxxxyzz[k] = -g_0_y_yyyz_xxxxyzz[k] * ab_x + g_0_y_yyyz_xxxxxyzz[k];

                g_0_y_xyyyz_xxxxzzz[k] = -g_0_y_yyyz_xxxxzzz[k] * ab_x + g_0_y_yyyz_xxxxxzzz[k];

                g_0_y_xyyyz_xxxyyyy[k] = -g_0_y_yyyz_xxxyyyy[k] * ab_x + g_0_y_yyyz_xxxxyyyy[k];

                g_0_y_xyyyz_xxxyyyz[k] = -g_0_y_yyyz_xxxyyyz[k] * ab_x + g_0_y_yyyz_xxxxyyyz[k];

                g_0_y_xyyyz_xxxyyzz[k] = -g_0_y_yyyz_xxxyyzz[k] * ab_x + g_0_y_yyyz_xxxxyyzz[k];

                g_0_y_xyyyz_xxxyzzz[k] = -g_0_y_yyyz_xxxyzzz[k] * ab_x + g_0_y_yyyz_xxxxyzzz[k];

                g_0_y_xyyyz_xxxzzzz[k] = -g_0_y_yyyz_xxxzzzz[k] * ab_x + g_0_y_yyyz_xxxxzzzz[k];

                g_0_y_xyyyz_xxyyyyy[k] = -g_0_y_yyyz_xxyyyyy[k] * ab_x + g_0_y_yyyz_xxxyyyyy[k];

                g_0_y_xyyyz_xxyyyyz[k] = -g_0_y_yyyz_xxyyyyz[k] * ab_x + g_0_y_yyyz_xxxyyyyz[k];

                g_0_y_xyyyz_xxyyyzz[k] = -g_0_y_yyyz_xxyyyzz[k] * ab_x + g_0_y_yyyz_xxxyyyzz[k];

                g_0_y_xyyyz_xxyyzzz[k] = -g_0_y_yyyz_xxyyzzz[k] * ab_x + g_0_y_yyyz_xxxyyzzz[k];

                g_0_y_xyyyz_xxyzzzz[k] = -g_0_y_yyyz_xxyzzzz[k] * ab_x + g_0_y_yyyz_xxxyzzzz[k];

                g_0_y_xyyyz_xxzzzzz[k] = -g_0_y_yyyz_xxzzzzz[k] * ab_x + g_0_y_yyyz_xxxzzzzz[k];

                g_0_y_xyyyz_xyyyyyy[k] = -g_0_y_yyyz_xyyyyyy[k] * ab_x + g_0_y_yyyz_xxyyyyyy[k];

                g_0_y_xyyyz_xyyyyyz[k] = -g_0_y_yyyz_xyyyyyz[k] * ab_x + g_0_y_yyyz_xxyyyyyz[k];

                g_0_y_xyyyz_xyyyyzz[k] = -g_0_y_yyyz_xyyyyzz[k] * ab_x + g_0_y_yyyz_xxyyyyzz[k];

                g_0_y_xyyyz_xyyyzzz[k] = -g_0_y_yyyz_xyyyzzz[k] * ab_x + g_0_y_yyyz_xxyyyzzz[k];

                g_0_y_xyyyz_xyyzzzz[k] = -g_0_y_yyyz_xyyzzzz[k] * ab_x + g_0_y_yyyz_xxyyzzzz[k];

                g_0_y_xyyyz_xyzzzzz[k] = -g_0_y_yyyz_xyzzzzz[k] * ab_x + g_0_y_yyyz_xxyzzzzz[k];

                g_0_y_xyyyz_xzzzzzz[k] = -g_0_y_yyyz_xzzzzzz[k] * ab_x + g_0_y_yyyz_xxzzzzzz[k];

                g_0_y_xyyyz_yyyyyyy[k] = -g_0_y_yyyz_yyyyyyy[k] * ab_x + g_0_y_yyyz_xyyyyyyy[k];

                g_0_y_xyyyz_yyyyyyz[k] = -g_0_y_yyyz_yyyyyyz[k] * ab_x + g_0_y_yyyz_xyyyyyyz[k];

                g_0_y_xyyyz_yyyyyzz[k] = -g_0_y_yyyz_yyyyyzz[k] * ab_x + g_0_y_yyyz_xyyyyyzz[k];

                g_0_y_xyyyz_yyyyzzz[k] = -g_0_y_yyyz_yyyyzzz[k] * ab_x + g_0_y_yyyz_xyyyyzzz[k];

                g_0_y_xyyyz_yyyzzzz[k] = -g_0_y_yyyz_yyyzzzz[k] * ab_x + g_0_y_yyyz_xyyyzzzz[k];

                g_0_y_xyyyz_yyzzzzz[k] = -g_0_y_yyyz_yyzzzzz[k] * ab_x + g_0_y_yyyz_xyyzzzzz[k];

                g_0_y_xyyyz_yzzzzzz[k] = -g_0_y_yyyz_yzzzzzz[k] * ab_x + g_0_y_yyyz_xyzzzzzz[k];

                g_0_y_xyyyz_zzzzzzz[k] = -g_0_y_yyyz_zzzzzzz[k] * ab_x + g_0_y_yyyz_xzzzzzzz[k];
            }

            /// Set up 1188-1224 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_y_xyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_y_xyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_y_xyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1223 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzz_xxxxxxx, g_0_y_xyyzz_xxxxxxy, g_0_y_xyyzz_xxxxxxz, g_0_y_xyyzz_xxxxxyy, g_0_y_xyyzz_xxxxxyz, g_0_y_xyyzz_xxxxxzz, g_0_y_xyyzz_xxxxyyy, g_0_y_xyyzz_xxxxyyz, g_0_y_xyyzz_xxxxyzz, g_0_y_xyyzz_xxxxzzz, g_0_y_xyyzz_xxxyyyy, g_0_y_xyyzz_xxxyyyz, g_0_y_xyyzz_xxxyyzz, g_0_y_xyyzz_xxxyzzz, g_0_y_xyyzz_xxxzzzz, g_0_y_xyyzz_xxyyyyy, g_0_y_xyyzz_xxyyyyz, g_0_y_xyyzz_xxyyyzz, g_0_y_xyyzz_xxyyzzz, g_0_y_xyyzz_xxyzzzz, g_0_y_xyyzz_xxzzzzz, g_0_y_xyyzz_xyyyyyy, g_0_y_xyyzz_xyyyyyz, g_0_y_xyyzz_xyyyyzz, g_0_y_xyyzz_xyyyzzz, g_0_y_xyyzz_xyyzzzz, g_0_y_xyyzz_xyzzzzz, g_0_y_xyyzz_xzzzzzz, g_0_y_xyyzz_yyyyyyy, g_0_y_xyyzz_yyyyyyz, g_0_y_xyyzz_yyyyyzz, g_0_y_xyyzz_yyyyzzz, g_0_y_xyyzz_yyyzzzz, g_0_y_xyyzz_yyzzzzz, g_0_y_xyyzz_yzzzzzz, g_0_y_xyyzz_zzzzzzz, g_0_y_yyzz_xxxxxxx, g_0_y_yyzz_xxxxxxxx, g_0_y_yyzz_xxxxxxxy, g_0_y_yyzz_xxxxxxxz, g_0_y_yyzz_xxxxxxy, g_0_y_yyzz_xxxxxxyy, g_0_y_yyzz_xxxxxxyz, g_0_y_yyzz_xxxxxxz, g_0_y_yyzz_xxxxxxzz, g_0_y_yyzz_xxxxxyy, g_0_y_yyzz_xxxxxyyy, g_0_y_yyzz_xxxxxyyz, g_0_y_yyzz_xxxxxyz, g_0_y_yyzz_xxxxxyzz, g_0_y_yyzz_xxxxxzz, g_0_y_yyzz_xxxxxzzz, g_0_y_yyzz_xxxxyyy, g_0_y_yyzz_xxxxyyyy, g_0_y_yyzz_xxxxyyyz, g_0_y_yyzz_xxxxyyz, g_0_y_yyzz_xxxxyyzz, g_0_y_yyzz_xxxxyzz, g_0_y_yyzz_xxxxyzzz, g_0_y_yyzz_xxxxzzz, g_0_y_yyzz_xxxxzzzz, g_0_y_yyzz_xxxyyyy, g_0_y_yyzz_xxxyyyyy, g_0_y_yyzz_xxxyyyyz, g_0_y_yyzz_xxxyyyz, g_0_y_yyzz_xxxyyyzz, g_0_y_yyzz_xxxyyzz, g_0_y_yyzz_xxxyyzzz, g_0_y_yyzz_xxxyzzz, g_0_y_yyzz_xxxyzzzz, g_0_y_yyzz_xxxzzzz, g_0_y_yyzz_xxxzzzzz, g_0_y_yyzz_xxyyyyy, g_0_y_yyzz_xxyyyyyy, g_0_y_yyzz_xxyyyyyz, g_0_y_yyzz_xxyyyyz, g_0_y_yyzz_xxyyyyzz, g_0_y_yyzz_xxyyyzz, g_0_y_yyzz_xxyyyzzz, g_0_y_yyzz_xxyyzzz, g_0_y_yyzz_xxyyzzzz, g_0_y_yyzz_xxyzzzz, g_0_y_yyzz_xxyzzzzz, g_0_y_yyzz_xxzzzzz, g_0_y_yyzz_xxzzzzzz, g_0_y_yyzz_xyyyyyy, g_0_y_yyzz_xyyyyyyy, g_0_y_yyzz_xyyyyyyz, g_0_y_yyzz_xyyyyyz, g_0_y_yyzz_xyyyyyzz, g_0_y_yyzz_xyyyyzz, g_0_y_yyzz_xyyyyzzz, g_0_y_yyzz_xyyyzzz, g_0_y_yyzz_xyyyzzzz, g_0_y_yyzz_xyyzzzz, g_0_y_yyzz_xyyzzzzz, g_0_y_yyzz_xyzzzzz, g_0_y_yyzz_xyzzzzzz, g_0_y_yyzz_xzzzzzz, g_0_y_yyzz_xzzzzzzz, g_0_y_yyzz_yyyyyyy, g_0_y_yyzz_yyyyyyz, g_0_y_yyzz_yyyyyzz, g_0_y_yyzz_yyyyzzz, g_0_y_yyzz_yyyzzzz, g_0_y_yyzz_yyzzzzz, g_0_y_yyzz_yzzzzzz, g_0_y_yyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzz_xxxxxxx[k] = -g_0_y_yyzz_xxxxxxx[k] * ab_x + g_0_y_yyzz_xxxxxxxx[k];

                g_0_y_xyyzz_xxxxxxy[k] = -g_0_y_yyzz_xxxxxxy[k] * ab_x + g_0_y_yyzz_xxxxxxxy[k];

                g_0_y_xyyzz_xxxxxxz[k] = -g_0_y_yyzz_xxxxxxz[k] * ab_x + g_0_y_yyzz_xxxxxxxz[k];

                g_0_y_xyyzz_xxxxxyy[k] = -g_0_y_yyzz_xxxxxyy[k] * ab_x + g_0_y_yyzz_xxxxxxyy[k];

                g_0_y_xyyzz_xxxxxyz[k] = -g_0_y_yyzz_xxxxxyz[k] * ab_x + g_0_y_yyzz_xxxxxxyz[k];

                g_0_y_xyyzz_xxxxxzz[k] = -g_0_y_yyzz_xxxxxzz[k] * ab_x + g_0_y_yyzz_xxxxxxzz[k];

                g_0_y_xyyzz_xxxxyyy[k] = -g_0_y_yyzz_xxxxyyy[k] * ab_x + g_0_y_yyzz_xxxxxyyy[k];

                g_0_y_xyyzz_xxxxyyz[k] = -g_0_y_yyzz_xxxxyyz[k] * ab_x + g_0_y_yyzz_xxxxxyyz[k];

                g_0_y_xyyzz_xxxxyzz[k] = -g_0_y_yyzz_xxxxyzz[k] * ab_x + g_0_y_yyzz_xxxxxyzz[k];

                g_0_y_xyyzz_xxxxzzz[k] = -g_0_y_yyzz_xxxxzzz[k] * ab_x + g_0_y_yyzz_xxxxxzzz[k];

                g_0_y_xyyzz_xxxyyyy[k] = -g_0_y_yyzz_xxxyyyy[k] * ab_x + g_0_y_yyzz_xxxxyyyy[k];

                g_0_y_xyyzz_xxxyyyz[k] = -g_0_y_yyzz_xxxyyyz[k] * ab_x + g_0_y_yyzz_xxxxyyyz[k];

                g_0_y_xyyzz_xxxyyzz[k] = -g_0_y_yyzz_xxxyyzz[k] * ab_x + g_0_y_yyzz_xxxxyyzz[k];

                g_0_y_xyyzz_xxxyzzz[k] = -g_0_y_yyzz_xxxyzzz[k] * ab_x + g_0_y_yyzz_xxxxyzzz[k];

                g_0_y_xyyzz_xxxzzzz[k] = -g_0_y_yyzz_xxxzzzz[k] * ab_x + g_0_y_yyzz_xxxxzzzz[k];

                g_0_y_xyyzz_xxyyyyy[k] = -g_0_y_yyzz_xxyyyyy[k] * ab_x + g_0_y_yyzz_xxxyyyyy[k];

                g_0_y_xyyzz_xxyyyyz[k] = -g_0_y_yyzz_xxyyyyz[k] * ab_x + g_0_y_yyzz_xxxyyyyz[k];

                g_0_y_xyyzz_xxyyyzz[k] = -g_0_y_yyzz_xxyyyzz[k] * ab_x + g_0_y_yyzz_xxxyyyzz[k];

                g_0_y_xyyzz_xxyyzzz[k] = -g_0_y_yyzz_xxyyzzz[k] * ab_x + g_0_y_yyzz_xxxyyzzz[k];

                g_0_y_xyyzz_xxyzzzz[k] = -g_0_y_yyzz_xxyzzzz[k] * ab_x + g_0_y_yyzz_xxxyzzzz[k];

                g_0_y_xyyzz_xxzzzzz[k] = -g_0_y_yyzz_xxzzzzz[k] * ab_x + g_0_y_yyzz_xxxzzzzz[k];

                g_0_y_xyyzz_xyyyyyy[k] = -g_0_y_yyzz_xyyyyyy[k] * ab_x + g_0_y_yyzz_xxyyyyyy[k];

                g_0_y_xyyzz_xyyyyyz[k] = -g_0_y_yyzz_xyyyyyz[k] * ab_x + g_0_y_yyzz_xxyyyyyz[k];

                g_0_y_xyyzz_xyyyyzz[k] = -g_0_y_yyzz_xyyyyzz[k] * ab_x + g_0_y_yyzz_xxyyyyzz[k];

                g_0_y_xyyzz_xyyyzzz[k] = -g_0_y_yyzz_xyyyzzz[k] * ab_x + g_0_y_yyzz_xxyyyzzz[k];

                g_0_y_xyyzz_xyyzzzz[k] = -g_0_y_yyzz_xyyzzzz[k] * ab_x + g_0_y_yyzz_xxyyzzzz[k];

                g_0_y_xyyzz_xyzzzzz[k] = -g_0_y_yyzz_xyzzzzz[k] * ab_x + g_0_y_yyzz_xxyzzzzz[k];

                g_0_y_xyyzz_xzzzzzz[k] = -g_0_y_yyzz_xzzzzzz[k] * ab_x + g_0_y_yyzz_xxzzzzzz[k];

                g_0_y_xyyzz_yyyyyyy[k] = -g_0_y_yyzz_yyyyyyy[k] * ab_x + g_0_y_yyzz_xyyyyyyy[k];

                g_0_y_xyyzz_yyyyyyz[k] = -g_0_y_yyzz_yyyyyyz[k] * ab_x + g_0_y_yyzz_xyyyyyyz[k];

                g_0_y_xyyzz_yyyyyzz[k] = -g_0_y_yyzz_yyyyyzz[k] * ab_x + g_0_y_yyzz_xyyyyyzz[k];

                g_0_y_xyyzz_yyyyzzz[k] = -g_0_y_yyzz_yyyyzzz[k] * ab_x + g_0_y_yyzz_xyyyyzzz[k];

                g_0_y_xyyzz_yyyzzzz[k] = -g_0_y_yyzz_yyyzzzz[k] * ab_x + g_0_y_yyzz_xyyyzzzz[k];

                g_0_y_xyyzz_yyzzzzz[k] = -g_0_y_yyzz_yyzzzzz[k] * ab_x + g_0_y_yyzz_xyyzzzzz[k];

                g_0_y_xyyzz_yzzzzzz[k] = -g_0_y_yyzz_yzzzzzz[k] * ab_x + g_0_y_yyzz_xyzzzzzz[k];

                g_0_y_xyyzz_zzzzzzz[k] = -g_0_y_yyzz_zzzzzzz[k] * ab_x + g_0_y_yyzz_xzzzzzzz[k];
            }

            /// Set up 1224-1260 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_y_xyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_y_xyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_y_xyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzz_xxxxxxx, g_0_y_xyzzz_xxxxxxy, g_0_y_xyzzz_xxxxxxz, g_0_y_xyzzz_xxxxxyy, g_0_y_xyzzz_xxxxxyz, g_0_y_xyzzz_xxxxxzz, g_0_y_xyzzz_xxxxyyy, g_0_y_xyzzz_xxxxyyz, g_0_y_xyzzz_xxxxyzz, g_0_y_xyzzz_xxxxzzz, g_0_y_xyzzz_xxxyyyy, g_0_y_xyzzz_xxxyyyz, g_0_y_xyzzz_xxxyyzz, g_0_y_xyzzz_xxxyzzz, g_0_y_xyzzz_xxxzzzz, g_0_y_xyzzz_xxyyyyy, g_0_y_xyzzz_xxyyyyz, g_0_y_xyzzz_xxyyyzz, g_0_y_xyzzz_xxyyzzz, g_0_y_xyzzz_xxyzzzz, g_0_y_xyzzz_xxzzzzz, g_0_y_xyzzz_xyyyyyy, g_0_y_xyzzz_xyyyyyz, g_0_y_xyzzz_xyyyyzz, g_0_y_xyzzz_xyyyzzz, g_0_y_xyzzz_xyyzzzz, g_0_y_xyzzz_xyzzzzz, g_0_y_xyzzz_xzzzzzz, g_0_y_xyzzz_yyyyyyy, g_0_y_xyzzz_yyyyyyz, g_0_y_xyzzz_yyyyyzz, g_0_y_xyzzz_yyyyzzz, g_0_y_xyzzz_yyyzzzz, g_0_y_xyzzz_yyzzzzz, g_0_y_xyzzz_yzzzzzz, g_0_y_xyzzz_zzzzzzz, g_0_y_yzzz_xxxxxxx, g_0_y_yzzz_xxxxxxxx, g_0_y_yzzz_xxxxxxxy, g_0_y_yzzz_xxxxxxxz, g_0_y_yzzz_xxxxxxy, g_0_y_yzzz_xxxxxxyy, g_0_y_yzzz_xxxxxxyz, g_0_y_yzzz_xxxxxxz, g_0_y_yzzz_xxxxxxzz, g_0_y_yzzz_xxxxxyy, g_0_y_yzzz_xxxxxyyy, g_0_y_yzzz_xxxxxyyz, g_0_y_yzzz_xxxxxyz, g_0_y_yzzz_xxxxxyzz, g_0_y_yzzz_xxxxxzz, g_0_y_yzzz_xxxxxzzz, g_0_y_yzzz_xxxxyyy, g_0_y_yzzz_xxxxyyyy, g_0_y_yzzz_xxxxyyyz, g_0_y_yzzz_xxxxyyz, g_0_y_yzzz_xxxxyyzz, g_0_y_yzzz_xxxxyzz, g_0_y_yzzz_xxxxyzzz, g_0_y_yzzz_xxxxzzz, g_0_y_yzzz_xxxxzzzz, g_0_y_yzzz_xxxyyyy, g_0_y_yzzz_xxxyyyyy, g_0_y_yzzz_xxxyyyyz, g_0_y_yzzz_xxxyyyz, g_0_y_yzzz_xxxyyyzz, g_0_y_yzzz_xxxyyzz, g_0_y_yzzz_xxxyyzzz, g_0_y_yzzz_xxxyzzz, g_0_y_yzzz_xxxyzzzz, g_0_y_yzzz_xxxzzzz, g_0_y_yzzz_xxxzzzzz, g_0_y_yzzz_xxyyyyy, g_0_y_yzzz_xxyyyyyy, g_0_y_yzzz_xxyyyyyz, g_0_y_yzzz_xxyyyyz, g_0_y_yzzz_xxyyyyzz, g_0_y_yzzz_xxyyyzz, g_0_y_yzzz_xxyyyzzz, g_0_y_yzzz_xxyyzzz, g_0_y_yzzz_xxyyzzzz, g_0_y_yzzz_xxyzzzz, g_0_y_yzzz_xxyzzzzz, g_0_y_yzzz_xxzzzzz, g_0_y_yzzz_xxzzzzzz, g_0_y_yzzz_xyyyyyy, g_0_y_yzzz_xyyyyyyy, g_0_y_yzzz_xyyyyyyz, g_0_y_yzzz_xyyyyyz, g_0_y_yzzz_xyyyyyzz, g_0_y_yzzz_xyyyyzz, g_0_y_yzzz_xyyyyzzz, g_0_y_yzzz_xyyyzzz, g_0_y_yzzz_xyyyzzzz, g_0_y_yzzz_xyyzzzz, g_0_y_yzzz_xyyzzzzz, g_0_y_yzzz_xyzzzzz, g_0_y_yzzz_xyzzzzzz, g_0_y_yzzz_xzzzzzz, g_0_y_yzzz_xzzzzzzz, g_0_y_yzzz_yyyyyyy, g_0_y_yzzz_yyyyyyz, g_0_y_yzzz_yyyyyzz, g_0_y_yzzz_yyyyzzz, g_0_y_yzzz_yyyzzzz, g_0_y_yzzz_yyzzzzz, g_0_y_yzzz_yzzzzzz, g_0_y_yzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzz_xxxxxxx[k] = -g_0_y_yzzz_xxxxxxx[k] * ab_x + g_0_y_yzzz_xxxxxxxx[k];

                g_0_y_xyzzz_xxxxxxy[k] = -g_0_y_yzzz_xxxxxxy[k] * ab_x + g_0_y_yzzz_xxxxxxxy[k];

                g_0_y_xyzzz_xxxxxxz[k] = -g_0_y_yzzz_xxxxxxz[k] * ab_x + g_0_y_yzzz_xxxxxxxz[k];

                g_0_y_xyzzz_xxxxxyy[k] = -g_0_y_yzzz_xxxxxyy[k] * ab_x + g_0_y_yzzz_xxxxxxyy[k];

                g_0_y_xyzzz_xxxxxyz[k] = -g_0_y_yzzz_xxxxxyz[k] * ab_x + g_0_y_yzzz_xxxxxxyz[k];

                g_0_y_xyzzz_xxxxxzz[k] = -g_0_y_yzzz_xxxxxzz[k] * ab_x + g_0_y_yzzz_xxxxxxzz[k];

                g_0_y_xyzzz_xxxxyyy[k] = -g_0_y_yzzz_xxxxyyy[k] * ab_x + g_0_y_yzzz_xxxxxyyy[k];

                g_0_y_xyzzz_xxxxyyz[k] = -g_0_y_yzzz_xxxxyyz[k] * ab_x + g_0_y_yzzz_xxxxxyyz[k];

                g_0_y_xyzzz_xxxxyzz[k] = -g_0_y_yzzz_xxxxyzz[k] * ab_x + g_0_y_yzzz_xxxxxyzz[k];

                g_0_y_xyzzz_xxxxzzz[k] = -g_0_y_yzzz_xxxxzzz[k] * ab_x + g_0_y_yzzz_xxxxxzzz[k];

                g_0_y_xyzzz_xxxyyyy[k] = -g_0_y_yzzz_xxxyyyy[k] * ab_x + g_0_y_yzzz_xxxxyyyy[k];

                g_0_y_xyzzz_xxxyyyz[k] = -g_0_y_yzzz_xxxyyyz[k] * ab_x + g_0_y_yzzz_xxxxyyyz[k];

                g_0_y_xyzzz_xxxyyzz[k] = -g_0_y_yzzz_xxxyyzz[k] * ab_x + g_0_y_yzzz_xxxxyyzz[k];

                g_0_y_xyzzz_xxxyzzz[k] = -g_0_y_yzzz_xxxyzzz[k] * ab_x + g_0_y_yzzz_xxxxyzzz[k];

                g_0_y_xyzzz_xxxzzzz[k] = -g_0_y_yzzz_xxxzzzz[k] * ab_x + g_0_y_yzzz_xxxxzzzz[k];

                g_0_y_xyzzz_xxyyyyy[k] = -g_0_y_yzzz_xxyyyyy[k] * ab_x + g_0_y_yzzz_xxxyyyyy[k];

                g_0_y_xyzzz_xxyyyyz[k] = -g_0_y_yzzz_xxyyyyz[k] * ab_x + g_0_y_yzzz_xxxyyyyz[k];

                g_0_y_xyzzz_xxyyyzz[k] = -g_0_y_yzzz_xxyyyzz[k] * ab_x + g_0_y_yzzz_xxxyyyzz[k];

                g_0_y_xyzzz_xxyyzzz[k] = -g_0_y_yzzz_xxyyzzz[k] * ab_x + g_0_y_yzzz_xxxyyzzz[k];

                g_0_y_xyzzz_xxyzzzz[k] = -g_0_y_yzzz_xxyzzzz[k] * ab_x + g_0_y_yzzz_xxxyzzzz[k];

                g_0_y_xyzzz_xxzzzzz[k] = -g_0_y_yzzz_xxzzzzz[k] * ab_x + g_0_y_yzzz_xxxzzzzz[k];

                g_0_y_xyzzz_xyyyyyy[k] = -g_0_y_yzzz_xyyyyyy[k] * ab_x + g_0_y_yzzz_xxyyyyyy[k];

                g_0_y_xyzzz_xyyyyyz[k] = -g_0_y_yzzz_xyyyyyz[k] * ab_x + g_0_y_yzzz_xxyyyyyz[k];

                g_0_y_xyzzz_xyyyyzz[k] = -g_0_y_yzzz_xyyyyzz[k] * ab_x + g_0_y_yzzz_xxyyyyzz[k];

                g_0_y_xyzzz_xyyyzzz[k] = -g_0_y_yzzz_xyyyzzz[k] * ab_x + g_0_y_yzzz_xxyyyzzz[k];

                g_0_y_xyzzz_xyyzzzz[k] = -g_0_y_yzzz_xyyzzzz[k] * ab_x + g_0_y_yzzz_xxyyzzzz[k];

                g_0_y_xyzzz_xyzzzzz[k] = -g_0_y_yzzz_xyzzzzz[k] * ab_x + g_0_y_yzzz_xxyzzzzz[k];

                g_0_y_xyzzz_xzzzzzz[k] = -g_0_y_yzzz_xzzzzzz[k] * ab_x + g_0_y_yzzz_xxzzzzzz[k];

                g_0_y_xyzzz_yyyyyyy[k] = -g_0_y_yzzz_yyyyyyy[k] * ab_x + g_0_y_yzzz_xyyyyyyy[k];

                g_0_y_xyzzz_yyyyyyz[k] = -g_0_y_yzzz_yyyyyyz[k] * ab_x + g_0_y_yzzz_xyyyyyyz[k];

                g_0_y_xyzzz_yyyyyzz[k] = -g_0_y_yzzz_yyyyyzz[k] * ab_x + g_0_y_yzzz_xyyyyyzz[k];

                g_0_y_xyzzz_yyyyzzz[k] = -g_0_y_yzzz_yyyyzzz[k] * ab_x + g_0_y_yzzz_xyyyyzzz[k];

                g_0_y_xyzzz_yyyzzzz[k] = -g_0_y_yzzz_yyyzzzz[k] * ab_x + g_0_y_yzzz_xyyyzzzz[k];

                g_0_y_xyzzz_yyzzzzz[k] = -g_0_y_yzzz_yyzzzzz[k] * ab_x + g_0_y_yzzz_xyyzzzzz[k];

                g_0_y_xyzzz_yzzzzzz[k] = -g_0_y_yzzz_yzzzzzz[k] * ab_x + g_0_y_yzzz_xyzzzzzz[k];

                g_0_y_xyzzz_zzzzzzz[k] = -g_0_y_yzzz_zzzzzzz[k] * ab_x + g_0_y_yzzz_xzzzzzzz[k];
            }

            /// Set up 1260-1296 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_y_xzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_y_xzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_y_xzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1295 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzz_xxxxxxx, g_0_y_xzzzz_xxxxxxy, g_0_y_xzzzz_xxxxxxz, g_0_y_xzzzz_xxxxxyy, g_0_y_xzzzz_xxxxxyz, g_0_y_xzzzz_xxxxxzz, g_0_y_xzzzz_xxxxyyy, g_0_y_xzzzz_xxxxyyz, g_0_y_xzzzz_xxxxyzz, g_0_y_xzzzz_xxxxzzz, g_0_y_xzzzz_xxxyyyy, g_0_y_xzzzz_xxxyyyz, g_0_y_xzzzz_xxxyyzz, g_0_y_xzzzz_xxxyzzz, g_0_y_xzzzz_xxxzzzz, g_0_y_xzzzz_xxyyyyy, g_0_y_xzzzz_xxyyyyz, g_0_y_xzzzz_xxyyyzz, g_0_y_xzzzz_xxyyzzz, g_0_y_xzzzz_xxyzzzz, g_0_y_xzzzz_xxzzzzz, g_0_y_xzzzz_xyyyyyy, g_0_y_xzzzz_xyyyyyz, g_0_y_xzzzz_xyyyyzz, g_0_y_xzzzz_xyyyzzz, g_0_y_xzzzz_xyyzzzz, g_0_y_xzzzz_xyzzzzz, g_0_y_xzzzz_xzzzzzz, g_0_y_xzzzz_yyyyyyy, g_0_y_xzzzz_yyyyyyz, g_0_y_xzzzz_yyyyyzz, g_0_y_xzzzz_yyyyzzz, g_0_y_xzzzz_yyyzzzz, g_0_y_xzzzz_yyzzzzz, g_0_y_xzzzz_yzzzzzz, g_0_y_xzzzz_zzzzzzz, g_0_y_zzzz_xxxxxxx, g_0_y_zzzz_xxxxxxxx, g_0_y_zzzz_xxxxxxxy, g_0_y_zzzz_xxxxxxxz, g_0_y_zzzz_xxxxxxy, g_0_y_zzzz_xxxxxxyy, g_0_y_zzzz_xxxxxxyz, g_0_y_zzzz_xxxxxxz, g_0_y_zzzz_xxxxxxzz, g_0_y_zzzz_xxxxxyy, g_0_y_zzzz_xxxxxyyy, g_0_y_zzzz_xxxxxyyz, g_0_y_zzzz_xxxxxyz, g_0_y_zzzz_xxxxxyzz, g_0_y_zzzz_xxxxxzz, g_0_y_zzzz_xxxxxzzz, g_0_y_zzzz_xxxxyyy, g_0_y_zzzz_xxxxyyyy, g_0_y_zzzz_xxxxyyyz, g_0_y_zzzz_xxxxyyz, g_0_y_zzzz_xxxxyyzz, g_0_y_zzzz_xxxxyzz, g_0_y_zzzz_xxxxyzzz, g_0_y_zzzz_xxxxzzz, g_0_y_zzzz_xxxxzzzz, g_0_y_zzzz_xxxyyyy, g_0_y_zzzz_xxxyyyyy, g_0_y_zzzz_xxxyyyyz, g_0_y_zzzz_xxxyyyz, g_0_y_zzzz_xxxyyyzz, g_0_y_zzzz_xxxyyzz, g_0_y_zzzz_xxxyyzzz, g_0_y_zzzz_xxxyzzz, g_0_y_zzzz_xxxyzzzz, g_0_y_zzzz_xxxzzzz, g_0_y_zzzz_xxxzzzzz, g_0_y_zzzz_xxyyyyy, g_0_y_zzzz_xxyyyyyy, g_0_y_zzzz_xxyyyyyz, g_0_y_zzzz_xxyyyyz, g_0_y_zzzz_xxyyyyzz, g_0_y_zzzz_xxyyyzz, g_0_y_zzzz_xxyyyzzz, g_0_y_zzzz_xxyyzzz, g_0_y_zzzz_xxyyzzzz, g_0_y_zzzz_xxyzzzz, g_0_y_zzzz_xxyzzzzz, g_0_y_zzzz_xxzzzzz, g_0_y_zzzz_xxzzzzzz, g_0_y_zzzz_xyyyyyy, g_0_y_zzzz_xyyyyyyy, g_0_y_zzzz_xyyyyyyz, g_0_y_zzzz_xyyyyyz, g_0_y_zzzz_xyyyyyzz, g_0_y_zzzz_xyyyyzz, g_0_y_zzzz_xyyyyzzz, g_0_y_zzzz_xyyyzzz, g_0_y_zzzz_xyyyzzzz, g_0_y_zzzz_xyyzzzz, g_0_y_zzzz_xyyzzzzz, g_0_y_zzzz_xyzzzzz, g_0_y_zzzz_xyzzzzzz, g_0_y_zzzz_xzzzzzz, g_0_y_zzzz_xzzzzzzz, g_0_y_zzzz_yyyyyyy, g_0_y_zzzz_yyyyyyz, g_0_y_zzzz_yyyyyzz, g_0_y_zzzz_yyyyzzz, g_0_y_zzzz_yyyzzzz, g_0_y_zzzz_yyzzzzz, g_0_y_zzzz_yzzzzzz, g_0_y_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzz_xxxxxxx[k] = -g_0_y_zzzz_xxxxxxx[k] * ab_x + g_0_y_zzzz_xxxxxxxx[k];

                g_0_y_xzzzz_xxxxxxy[k] = -g_0_y_zzzz_xxxxxxy[k] * ab_x + g_0_y_zzzz_xxxxxxxy[k];

                g_0_y_xzzzz_xxxxxxz[k] = -g_0_y_zzzz_xxxxxxz[k] * ab_x + g_0_y_zzzz_xxxxxxxz[k];

                g_0_y_xzzzz_xxxxxyy[k] = -g_0_y_zzzz_xxxxxyy[k] * ab_x + g_0_y_zzzz_xxxxxxyy[k];

                g_0_y_xzzzz_xxxxxyz[k] = -g_0_y_zzzz_xxxxxyz[k] * ab_x + g_0_y_zzzz_xxxxxxyz[k];

                g_0_y_xzzzz_xxxxxzz[k] = -g_0_y_zzzz_xxxxxzz[k] * ab_x + g_0_y_zzzz_xxxxxxzz[k];

                g_0_y_xzzzz_xxxxyyy[k] = -g_0_y_zzzz_xxxxyyy[k] * ab_x + g_0_y_zzzz_xxxxxyyy[k];

                g_0_y_xzzzz_xxxxyyz[k] = -g_0_y_zzzz_xxxxyyz[k] * ab_x + g_0_y_zzzz_xxxxxyyz[k];

                g_0_y_xzzzz_xxxxyzz[k] = -g_0_y_zzzz_xxxxyzz[k] * ab_x + g_0_y_zzzz_xxxxxyzz[k];

                g_0_y_xzzzz_xxxxzzz[k] = -g_0_y_zzzz_xxxxzzz[k] * ab_x + g_0_y_zzzz_xxxxxzzz[k];

                g_0_y_xzzzz_xxxyyyy[k] = -g_0_y_zzzz_xxxyyyy[k] * ab_x + g_0_y_zzzz_xxxxyyyy[k];

                g_0_y_xzzzz_xxxyyyz[k] = -g_0_y_zzzz_xxxyyyz[k] * ab_x + g_0_y_zzzz_xxxxyyyz[k];

                g_0_y_xzzzz_xxxyyzz[k] = -g_0_y_zzzz_xxxyyzz[k] * ab_x + g_0_y_zzzz_xxxxyyzz[k];

                g_0_y_xzzzz_xxxyzzz[k] = -g_0_y_zzzz_xxxyzzz[k] * ab_x + g_0_y_zzzz_xxxxyzzz[k];

                g_0_y_xzzzz_xxxzzzz[k] = -g_0_y_zzzz_xxxzzzz[k] * ab_x + g_0_y_zzzz_xxxxzzzz[k];

                g_0_y_xzzzz_xxyyyyy[k] = -g_0_y_zzzz_xxyyyyy[k] * ab_x + g_0_y_zzzz_xxxyyyyy[k];

                g_0_y_xzzzz_xxyyyyz[k] = -g_0_y_zzzz_xxyyyyz[k] * ab_x + g_0_y_zzzz_xxxyyyyz[k];

                g_0_y_xzzzz_xxyyyzz[k] = -g_0_y_zzzz_xxyyyzz[k] * ab_x + g_0_y_zzzz_xxxyyyzz[k];

                g_0_y_xzzzz_xxyyzzz[k] = -g_0_y_zzzz_xxyyzzz[k] * ab_x + g_0_y_zzzz_xxxyyzzz[k];

                g_0_y_xzzzz_xxyzzzz[k] = -g_0_y_zzzz_xxyzzzz[k] * ab_x + g_0_y_zzzz_xxxyzzzz[k];

                g_0_y_xzzzz_xxzzzzz[k] = -g_0_y_zzzz_xxzzzzz[k] * ab_x + g_0_y_zzzz_xxxzzzzz[k];

                g_0_y_xzzzz_xyyyyyy[k] = -g_0_y_zzzz_xyyyyyy[k] * ab_x + g_0_y_zzzz_xxyyyyyy[k];

                g_0_y_xzzzz_xyyyyyz[k] = -g_0_y_zzzz_xyyyyyz[k] * ab_x + g_0_y_zzzz_xxyyyyyz[k];

                g_0_y_xzzzz_xyyyyzz[k] = -g_0_y_zzzz_xyyyyzz[k] * ab_x + g_0_y_zzzz_xxyyyyzz[k];

                g_0_y_xzzzz_xyyyzzz[k] = -g_0_y_zzzz_xyyyzzz[k] * ab_x + g_0_y_zzzz_xxyyyzzz[k];

                g_0_y_xzzzz_xyyzzzz[k] = -g_0_y_zzzz_xyyzzzz[k] * ab_x + g_0_y_zzzz_xxyyzzzz[k];

                g_0_y_xzzzz_xyzzzzz[k] = -g_0_y_zzzz_xyzzzzz[k] * ab_x + g_0_y_zzzz_xxyzzzzz[k];

                g_0_y_xzzzz_xzzzzzz[k] = -g_0_y_zzzz_xzzzzzz[k] * ab_x + g_0_y_zzzz_xxzzzzzz[k];

                g_0_y_xzzzz_yyyyyyy[k] = -g_0_y_zzzz_yyyyyyy[k] * ab_x + g_0_y_zzzz_xyyyyyyy[k];

                g_0_y_xzzzz_yyyyyyz[k] = -g_0_y_zzzz_yyyyyyz[k] * ab_x + g_0_y_zzzz_xyyyyyyz[k];

                g_0_y_xzzzz_yyyyyzz[k] = -g_0_y_zzzz_yyyyyzz[k] * ab_x + g_0_y_zzzz_xyyyyyzz[k];

                g_0_y_xzzzz_yyyyzzz[k] = -g_0_y_zzzz_yyyyzzz[k] * ab_x + g_0_y_zzzz_xyyyyzzz[k];

                g_0_y_xzzzz_yyyzzzz[k] = -g_0_y_zzzz_yyyzzzz[k] * ab_x + g_0_y_zzzz_xyyyzzzz[k];

                g_0_y_xzzzz_yyzzzzz[k] = -g_0_y_zzzz_yyzzzzz[k] * ab_x + g_0_y_zzzz_xyyzzzzz[k];

                g_0_y_xzzzz_yzzzzzz[k] = -g_0_y_zzzz_yzzzzzz[k] * ab_x + g_0_y_zzzz_xyzzzzzz[k];

                g_0_y_xzzzz_zzzzzzz[k] = -g_0_y_zzzz_zzzzzzz[k] * ab_x + g_0_y_zzzz_xzzzzzzz[k];
            }

            /// Set up 1296-1332 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_y_yyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_y_yyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_y_yyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1331 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xxxxxxx, g_0_y_yyyy_xxxxxxxy, g_0_y_yyyy_xxxxxxy, g_0_y_yyyy_xxxxxxyy, g_0_y_yyyy_xxxxxxyz, g_0_y_yyyy_xxxxxxz, g_0_y_yyyy_xxxxxyy, g_0_y_yyyy_xxxxxyyy, g_0_y_yyyy_xxxxxyyz, g_0_y_yyyy_xxxxxyz, g_0_y_yyyy_xxxxxyzz, g_0_y_yyyy_xxxxxzz, g_0_y_yyyy_xxxxyyy, g_0_y_yyyy_xxxxyyyy, g_0_y_yyyy_xxxxyyyz, g_0_y_yyyy_xxxxyyz, g_0_y_yyyy_xxxxyyzz, g_0_y_yyyy_xxxxyzz, g_0_y_yyyy_xxxxyzzz, g_0_y_yyyy_xxxxzzz, g_0_y_yyyy_xxxyyyy, g_0_y_yyyy_xxxyyyyy, g_0_y_yyyy_xxxyyyyz, g_0_y_yyyy_xxxyyyz, g_0_y_yyyy_xxxyyyzz, g_0_y_yyyy_xxxyyzz, g_0_y_yyyy_xxxyyzzz, g_0_y_yyyy_xxxyzzz, g_0_y_yyyy_xxxyzzzz, g_0_y_yyyy_xxxzzzz, g_0_y_yyyy_xxyyyyy, g_0_y_yyyy_xxyyyyyy, g_0_y_yyyy_xxyyyyyz, g_0_y_yyyy_xxyyyyz, g_0_y_yyyy_xxyyyyzz, g_0_y_yyyy_xxyyyzz, g_0_y_yyyy_xxyyyzzz, g_0_y_yyyy_xxyyzzz, g_0_y_yyyy_xxyyzzzz, g_0_y_yyyy_xxyzzzz, g_0_y_yyyy_xxyzzzzz, g_0_y_yyyy_xxzzzzz, g_0_y_yyyy_xyyyyyy, g_0_y_yyyy_xyyyyyyy, g_0_y_yyyy_xyyyyyyz, g_0_y_yyyy_xyyyyyz, g_0_y_yyyy_xyyyyyzz, g_0_y_yyyy_xyyyyzz, g_0_y_yyyy_xyyyyzzz, g_0_y_yyyy_xyyyzzz, g_0_y_yyyy_xyyyzzzz, g_0_y_yyyy_xyyzzzz, g_0_y_yyyy_xyyzzzzz, g_0_y_yyyy_xyzzzzz, g_0_y_yyyy_xyzzzzzz, g_0_y_yyyy_xzzzzzz, g_0_y_yyyy_yyyyyyy, g_0_y_yyyy_yyyyyyyy, g_0_y_yyyy_yyyyyyyz, g_0_y_yyyy_yyyyyyz, g_0_y_yyyy_yyyyyyzz, g_0_y_yyyy_yyyyyzz, g_0_y_yyyy_yyyyyzzz, g_0_y_yyyy_yyyyzzz, g_0_y_yyyy_yyyyzzzz, g_0_y_yyyy_yyyzzzz, g_0_y_yyyy_yyyzzzzz, g_0_y_yyyy_yyzzzzz, g_0_y_yyyy_yyzzzzzz, g_0_y_yyyy_yzzzzzz, g_0_y_yyyy_yzzzzzzz, g_0_y_yyyy_zzzzzzz, g_0_y_yyyyy_xxxxxxx, g_0_y_yyyyy_xxxxxxy, g_0_y_yyyyy_xxxxxxz, g_0_y_yyyyy_xxxxxyy, g_0_y_yyyyy_xxxxxyz, g_0_y_yyyyy_xxxxxzz, g_0_y_yyyyy_xxxxyyy, g_0_y_yyyyy_xxxxyyz, g_0_y_yyyyy_xxxxyzz, g_0_y_yyyyy_xxxxzzz, g_0_y_yyyyy_xxxyyyy, g_0_y_yyyyy_xxxyyyz, g_0_y_yyyyy_xxxyyzz, g_0_y_yyyyy_xxxyzzz, g_0_y_yyyyy_xxxzzzz, g_0_y_yyyyy_xxyyyyy, g_0_y_yyyyy_xxyyyyz, g_0_y_yyyyy_xxyyyzz, g_0_y_yyyyy_xxyyzzz, g_0_y_yyyyy_xxyzzzz, g_0_y_yyyyy_xxzzzzz, g_0_y_yyyyy_xyyyyyy, g_0_y_yyyyy_xyyyyyz, g_0_y_yyyyy_xyyyyzz, g_0_y_yyyyy_xyyyzzz, g_0_y_yyyyy_xyyzzzz, g_0_y_yyyyy_xyzzzzz, g_0_y_yyyyy_xzzzzzz, g_0_y_yyyyy_yyyyyyy, g_0_y_yyyyy_yyyyyyz, g_0_y_yyyyy_yyyyyzz, g_0_y_yyyyy_yyyyzzz, g_0_y_yyyyy_yyyzzzz, g_0_y_yyyyy_yyzzzzz, g_0_y_yyyyy_yzzzzzz, g_0_y_yyyyy_zzzzzzz, g_yyyy_xxxxxxx, g_yyyy_xxxxxxy, g_yyyy_xxxxxxz, g_yyyy_xxxxxyy, g_yyyy_xxxxxyz, g_yyyy_xxxxxzz, g_yyyy_xxxxyyy, g_yyyy_xxxxyyz, g_yyyy_xxxxyzz, g_yyyy_xxxxzzz, g_yyyy_xxxyyyy, g_yyyy_xxxyyyz, g_yyyy_xxxyyzz, g_yyyy_xxxyzzz, g_yyyy_xxxzzzz, g_yyyy_xxyyyyy, g_yyyy_xxyyyyz, g_yyyy_xxyyyzz, g_yyyy_xxyyzzz, g_yyyy_xxyzzzz, g_yyyy_xxzzzzz, g_yyyy_xyyyyyy, g_yyyy_xyyyyyz, g_yyyy_xyyyyzz, g_yyyy_xyyyzzz, g_yyyy_xyyzzzz, g_yyyy_xyzzzzz, g_yyyy_xzzzzzz, g_yyyy_yyyyyyy, g_yyyy_yyyyyyz, g_yyyy_yyyyyzz, g_yyyy_yyyyzzz, g_yyyy_yyyzzzz, g_yyyy_yyzzzzz, g_yyyy_yzzzzzz, g_yyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyy_xxxxxxx[k] = g_yyyy_xxxxxxx[k] - g_0_y_yyyy_xxxxxxx[k] * ab_y + g_0_y_yyyy_xxxxxxxy[k];

                g_0_y_yyyyy_xxxxxxy[k] = g_yyyy_xxxxxxy[k] - g_0_y_yyyy_xxxxxxy[k] * ab_y + g_0_y_yyyy_xxxxxxyy[k];

                g_0_y_yyyyy_xxxxxxz[k] = g_yyyy_xxxxxxz[k] - g_0_y_yyyy_xxxxxxz[k] * ab_y + g_0_y_yyyy_xxxxxxyz[k];

                g_0_y_yyyyy_xxxxxyy[k] = g_yyyy_xxxxxyy[k] - g_0_y_yyyy_xxxxxyy[k] * ab_y + g_0_y_yyyy_xxxxxyyy[k];

                g_0_y_yyyyy_xxxxxyz[k] = g_yyyy_xxxxxyz[k] - g_0_y_yyyy_xxxxxyz[k] * ab_y + g_0_y_yyyy_xxxxxyyz[k];

                g_0_y_yyyyy_xxxxxzz[k] = g_yyyy_xxxxxzz[k] - g_0_y_yyyy_xxxxxzz[k] * ab_y + g_0_y_yyyy_xxxxxyzz[k];

                g_0_y_yyyyy_xxxxyyy[k] = g_yyyy_xxxxyyy[k] - g_0_y_yyyy_xxxxyyy[k] * ab_y + g_0_y_yyyy_xxxxyyyy[k];

                g_0_y_yyyyy_xxxxyyz[k] = g_yyyy_xxxxyyz[k] - g_0_y_yyyy_xxxxyyz[k] * ab_y + g_0_y_yyyy_xxxxyyyz[k];

                g_0_y_yyyyy_xxxxyzz[k] = g_yyyy_xxxxyzz[k] - g_0_y_yyyy_xxxxyzz[k] * ab_y + g_0_y_yyyy_xxxxyyzz[k];

                g_0_y_yyyyy_xxxxzzz[k] = g_yyyy_xxxxzzz[k] - g_0_y_yyyy_xxxxzzz[k] * ab_y + g_0_y_yyyy_xxxxyzzz[k];

                g_0_y_yyyyy_xxxyyyy[k] = g_yyyy_xxxyyyy[k] - g_0_y_yyyy_xxxyyyy[k] * ab_y + g_0_y_yyyy_xxxyyyyy[k];

                g_0_y_yyyyy_xxxyyyz[k] = g_yyyy_xxxyyyz[k] - g_0_y_yyyy_xxxyyyz[k] * ab_y + g_0_y_yyyy_xxxyyyyz[k];

                g_0_y_yyyyy_xxxyyzz[k] = g_yyyy_xxxyyzz[k] - g_0_y_yyyy_xxxyyzz[k] * ab_y + g_0_y_yyyy_xxxyyyzz[k];

                g_0_y_yyyyy_xxxyzzz[k] = g_yyyy_xxxyzzz[k] - g_0_y_yyyy_xxxyzzz[k] * ab_y + g_0_y_yyyy_xxxyyzzz[k];

                g_0_y_yyyyy_xxxzzzz[k] = g_yyyy_xxxzzzz[k] - g_0_y_yyyy_xxxzzzz[k] * ab_y + g_0_y_yyyy_xxxyzzzz[k];

                g_0_y_yyyyy_xxyyyyy[k] = g_yyyy_xxyyyyy[k] - g_0_y_yyyy_xxyyyyy[k] * ab_y + g_0_y_yyyy_xxyyyyyy[k];

                g_0_y_yyyyy_xxyyyyz[k] = g_yyyy_xxyyyyz[k] - g_0_y_yyyy_xxyyyyz[k] * ab_y + g_0_y_yyyy_xxyyyyyz[k];

                g_0_y_yyyyy_xxyyyzz[k] = g_yyyy_xxyyyzz[k] - g_0_y_yyyy_xxyyyzz[k] * ab_y + g_0_y_yyyy_xxyyyyzz[k];

                g_0_y_yyyyy_xxyyzzz[k] = g_yyyy_xxyyzzz[k] - g_0_y_yyyy_xxyyzzz[k] * ab_y + g_0_y_yyyy_xxyyyzzz[k];

                g_0_y_yyyyy_xxyzzzz[k] = g_yyyy_xxyzzzz[k] - g_0_y_yyyy_xxyzzzz[k] * ab_y + g_0_y_yyyy_xxyyzzzz[k];

                g_0_y_yyyyy_xxzzzzz[k] = g_yyyy_xxzzzzz[k] - g_0_y_yyyy_xxzzzzz[k] * ab_y + g_0_y_yyyy_xxyzzzzz[k];

                g_0_y_yyyyy_xyyyyyy[k] = g_yyyy_xyyyyyy[k] - g_0_y_yyyy_xyyyyyy[k] * ab_y + g_0_y_yyyy_xyyyyyyy[k];

                g_0_y_yyyyy_xyyyyyz[k] = g_yyyy_xyyyyyz[k] - g_0_y_yyyy_xyyyyyz[k] * ab_y + g_0_y_yyyy_xyyyyyyz[k];

                g_0_y_yyyyy_xyyyyzz[k] = g_yyyy_xyyyyzz[k] - g_0_y_yyyy_xyyyyzz[k] * ab_y + g_0_y_yyyy_xyyyyyzz[k];

                g_0_y_yyyyy_xyyyzzz[k] = g_yyyy_xyyyzzz[k] - g_0_y_yyyy_xyyyzzz[k] * ab_y + g_0_y_yyyy_xyyyyzzz[k];

                g_0_y_yyyyy_xyyzzzz[k] = g_yyyy_xyyzzzz[k] - g_0_y_yyyy_xyyzzzz[k] * ab_y + g_0_y_yyyy_xyyyzzzz[k];

                g_0_y_yyyyy_xyzzzzz[k] = g_yyyy_xyzzzzz[k] - g_0_y_yyyy_xyzzzzz[k] * ab_y + g_0_y_yyyy_xyyzzzzz[k];

                g_0_y_yyyyy_xzzzzzz[k] = g_yyyy_xzzzzzz[k] - g_0_y_yyyy_xzzzzzz[k] * ab_y + g_0_y_yyyy_xyzzzzzz[k];

                g_0_y_yyyyy_yyyyyyy[k] = g_yyyy_yyyyyyy[k] - g_0_y_yyyy_yyyyyyy[k] * ab_y + g_0_y_yyyy_yyyyyyyy[k];

                g_0_y_yyyyy_yyyyyyz[k] = g_yyyy_yyyyyyz[k] - g_0_y_yyyy_yyyyyyz[k] * ab_y + g_0_y_yyyy_yyyyyyyz[k];

                g_0_y_yyyyy_yyyyyzz[k] = g_yyyy_yyyyyzz[k] - g_0_y_yyyy_yyyyyzz[k] * ab_y + g_0_y_yyyy_yyyyyyzz[k];

                g_0_y_yyyyy_yyyyzzz[k] = g_yyyy_yyyyzzz[k] - g_0_y_yyyy_yyyyzzz[k] * ab_y + g_0_y_yyyy_yyyyyzzz[k];

                g_0_y_yyyyy_yyyzzzz[k] = g_yyyy_yyyzzzz[k] - g_0_y_yyyy_yyyzzzz[k] * ab_y + g_0_y_yyyy_yyyyzzzz[k];

                g_0_y_yyyyy_yyzzzzz[k] = g_yyyy_yyzzzzz[k] - g_0_y_yyyy_yyzzzzz[k] * ab_y + g_0_y_yyyy_yyyzzzzz[k];

                g_0_y_yyyyy_yzzzzzz[k] = g_yyyy_yzzzzzz[k] - g_0_y_yyyy_yzzzzzz[k] * ab_y + g_0_y_yyyy_yyzzzzzz[k];

                g_0_y_yyyyy_zzzzzzz[k] = g_yyyy_zzzzzzz[k] - g_0_y_yyyy_zzzzzzz[k] * ab_y + g_0_y_yyyy_yzzzzzzz[k];
            }

            /// Set up 1332-1368 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1334 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1349 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1350 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1351 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1352 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1353 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1354 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1355 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1356 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1357 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1358 * ccomps * dcomps);

            auto g_0_y_yyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1359 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1360 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1361 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1362 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1363 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1364 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1365 * ccomps * dcomps);

            auto g_0_y_yyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1366 * ccomps * dcomps);

            auto g_0_y_yyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1367 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xxxxxxx, g_0_y_yyyy_xxxxxxxz, g_0_y_yyyy_xxxxxxy, g_0_y_yyyy_xxxxxxyz, g_0_y_yyyy_xxxxxxz, g_0_y_yyyy_xxxxxxzz, g_0_y_yyyy_xxxxxyy, g_0_y_yyyy_xxxxxyyz, g_0_y_yyyy_xxxxxyz, g_0_y_yyyy_xxxxxyzz, g_0_y_yyyy_xxxxxzz, g_0_y_yyyy_xxxxxzzz, g_0_y_yyyy_xxxxyyy, g_0_y_yyyy_xxxxyyyz, g_0_y_yyyy_xxxxyyz, g_0_y_yyyy_xxxxyyzz, g_0_y_yyyy_xxxxyzz, g_0_y_yyyy_xxxxyzzz, g_0_y_yyyy_xxxxzzz, g_0_y_yyyy_xxxxzzzz, g_0_y_yyyy_xxxyyyy, g_0_y_yyyy_xxxyyyyz, g_0_y_yyyy_xxxyyyz, g_0_y_yyyy_xxxyyyzz, g_0_y_yyyy_xxxyyzz, g_0_y_yyyy_xxxyyzzz, g_0_y_yyyy_xxxyzzz, g_0_y_yyyy_xxxyzzzz, g_0_y_yyyy_xxxzzzz, g_0_y_yyyy_xxxzzzzz, g_0_y_yyyy_xxyyyyy, g_0_y_yyyy_xxyyyyyz, g_0_y_yyyy_xxyyyyz, g_0_y_yyyy_xxyyyyzz, g_0_y_yyyy_xxyyyzz, g_0_y_yyyy_xxyyyzzz, g_0_y_yyyy_xxyyzzz, g_0_y_yyyy_xxyyzzzz, g_0_y_yyyy_xxyzzzz, g_0_y_yyyy_xxyzzzzz, g_0_y_yyyy_xxzzzzz, g_0_y_yyyy_xxzzzzzz, g_0_y_yyyy_xyyyyyy, g_0_y_yyyy_xyyyyyyz, g_0_y_yyyy_xyyyyyz, g_0_y_yyyy_xyyyyyzz, g_0_y_yyyy_xyyyyzz, g_0_y_yyyy_xyyyyzzz, g_0_y_yyyy_xyyyzzz, g_0_y_yyyy_xyyyzzzz, g_0_y_yyyy_xyyzzzz, g_0_y_yyyy_xyyzzzzz, g_0_y_yyyy_xyzzzzz, g_0_y_yyyy_xyzzzzzz, g_0_y_yyyy_xzzzzzz, g_0_y_yyyy_xzzzzzzz, g_0_y_yyyy_yyyyyyy, g_0_y_yyyy_yyyyyyyz, g_0_y_yyyy_yyyyyyz, g_0_y_yyyy_yyyyyyzz, g_0_y_yyyy_yyyyyzz, g_0_y_yyyy_yyyyyzzz, g_0_y_yyyy_yyyyzzz, g_0_y_yyyy_yyyyzzzz, g_0_y_yyyy_yyyzzzz, g_0_y_yyyy_yyyzzzzz, g_0_y_yyyy_yyzzzzz, g_0_y_yyyy_yyzzzzzz, g_0_y_yyyy_yzzzzzz, g_0_y_yyyy_yzzzzzzz, g_0_y_yyyy_zzzzzzz, g_0_y_yyyy_zzzzzzzz, g_0_y_yyyyz_xxxxxxx, g_0_y_yyyyz_xxxxxxy, g_0_y_yyyyz_xxxxxxz, g_0_y_yyyyz_xxxxxyy, g_0_y_yyyyz_xxxxxyz, g_0_y_yyyyz_xxxxxzz, g_0_y_yyyyz_xxxxyyy, g_0_y_yyyyz_xxxxyyz, g_0_y_yyyyz_xxxxyzz, g_0_y_yyyyz_xxxxzzz, g_0_y_yyyyz_xxxyyyy, g_0_y_yyyyz_xxxyyyz, g_0_y_yyyyz_xxxyyzz, g_0_y_yyyyz_xxxyzzz, g_0_y_yyyyz_xxxzzzz, g_0_y_yyyyz_xxyyyyy, g_0_y_yyyyz_xxyyyyz, g_0_y_yyyyz_xxyyyzz, g_0_y_yyyyz_xxyyzzz, g_0_y_yyyyz_xxyzzzz, g_0_y_yyyyz_xxzzzzz, g_0_y_yyyyz_xyyyyyy, g_0_y_yyyyz_xyyyyyz, g_0_y_yyyyz_xyyyyzz, g_0_y_yyyyz_xyyyzzz, g_0_y_yyyyz_xyyzzzz, g_0_y_yyyyz_xyzzzzz, g_0_y_yyyyz_xzzzzzz, g_0_y_yyyyz_yyyyyyy, g_0_y_yyyyz_yyyyyyz, g_0_y_yyyyz_yyyyyzz, g_0_y_yyyyz_yyyyzzz, g_0_y_yyyyz_yyyzzzz, g_0_y_yyyyz_yyzzzzz, g_0_y_yyyyz_yzzzzzz, g_0_y_yyyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyz_xxxxxxx[k] = -g_0_y_yyyy_xxxxxxx[k] * ab_z + g_0_y_yyyy_xxxxxxxz[k];

                g_0_y_yyyyz_xxxxxxy[k] = -g_0_y_yyyy_xxxxxxy[k] * ab_z + g_0_y_yyyy_xxxxxxyz[k];

                g_0_y_yyyyz_xxxxxxz[k] = -g_0_y_yyyy_xxxxxxz[k] * ab_z + g_0_y_yyyy_xxxxxxzz[k];

                g_0_y_yyyyz_xxxxxyy[k] = -g_0_y_yyyy_xxxxxyy[k] * ab_z + g_0_y_yyyy_xxxxxyyz[k];

                g_0_y_yyyyz_xxxxxyz[k] = -g_0_y_yyyy_xxxxxyz[k] * ab_z + g_0_y_yyyy_xxxxxyzz[k];

                g_0_y_yyyyz_xxxxxzz[k] = -g_0_y_yyyy_xxxxxzz[k] * ab_z + g_0_y_yyyy_xxxxxzzz[k];

                g_0_y_yyyyz_xxxxyyy[k] = -g_0_y_yyyy_xxxxyyy[k] * ab_z + g_0_y_yyyy_xxxxyyyz[k];

                g_0_y_yyyyz_xxxxyyz[k] = -g_0_y_yyyy_xxxxyyz[k] * ab_z + g_0_y_yyyy_xxxxyyzz[k];

                g_0_y_yyyyz_xxxxyzz[k] = -g_0_y_yyyy_xxxxyzz[k] * ab_z + g_0_y_yyyy_xxxxyzzz[k];

                g_0_y_yyyyz_xxxxzzz[k] = -g_0_y_yyyy_xxxxzzz[k] * ab_z + g_0_y_yyyy_xxxxzzzz[k];

                g_0_y_yyyyz_xxxyyyy[k] = -g_0_y_yyyy_xxxyyyy[k] * ab_z + g_0_y_yyyy_xxxyyyyz[k];

                g_0_y_yyyyz_xxxyyyz[k] = -g_0_y_yyyy_xxxyyyz[k] * ab_z + g_0_y_yyyy_xxxyyyzz[k];

                g_0_y_yyyyz_xxxyyzz[k] = -g_0_y_yyyy_xxxyyzz[k] * ab_z + g_0_y_yyyy_xxxyyzzz[k];

                g_0_y_yyyyz_xxxyzzz[k] = -g_0_y_yyyy_xxxyzzz[k] * ab_z + g_0_y_yyyy_xxxyzzzz[k];

                g_0_y_yyyyz_xxxzzzz[k] = -g_0_y_yyyy_xxxzzzz[k] * ab_z + g_0_y_yyyy_xxxzzzzz[k];

                g_0_y_yyyyz_xxyyyyy[k] = -g_0_y_yyyy_xxyyyyy[k] * ab_z + g_0_y_yyyy_xxyyyyyz[k];

                g_0_y_yyyyz_xxyyyyz[k] = -g_0_y_yyyy_xxyyyyz[k] * ab_z + g_0_y_yyyy_xxyyyyzz[k];

                g_0_y_yyyyz_xxyyyzz[k] = -g_0_y_yyyy_xxyyyzz[k] * ab_z + g_0_y_yyyy_xxyyyzzz[k];

                g_0_y_yyyyz_xxyyzzz[k] = -g_0_y_yyyy_xxyyzzz[k] * ab_z + g_0_y_yyyy_xxyyzzzz[k];

                g_0_y_yyyyz_xxyzzzz[k] = -g_0_y_yyyy_xxyzzzz[k] * ab_z + g_0_y_yyyy_xxyzzzzz[k];

                g_0_y_yyyyz_xxzzzzz[k] = -g_0_y_yyyy_xxzzzzz[k] * ab_z + g_0_y_yyyy_xxzzzzzz[k];

                g_0_y_yyyyz_xyyyyyy[k] = -g_0_y_yyyy_xyyyyyy[k] * ab_z + g_0_y_yyyy_xyyyyyyz[k];

                g_0_y_yyyyz_xyyyyyz[k] = -g_0_y_yyyy_xyyyyyz[k] * ab_z + g_0_y_yyyy_xyyyyyzz[k];

                g_0_y_yyyyz_xyyyyzz[k] = -g_0_y_yyyy_xyyyyzz[k] * ab_z + g_0_y_yyyy_xyyyyzzz[k];

                g_0_y_yyyyz_xyyyzzz[k] = -g_0_y_yyyy_xyyyzzz[k] * ab_z + g_0_y_yyyy_xyyyzzzz[k];

                g_0_y_yyyyz_xyyzzzz[k] = -g_0_y_yyyy_xyyzzzz[k] * ab_z + g_0_y_yyyy_xyyzzzzz[k];

                g_0_y_yyyyz_xyzzzzz[k] = -g_0_y_yyyy_xyzzzzz[k] * ab_z + g_0_y_yyyy_xyzzzzzz[k];

                g_0_y_yyyyz_xzzzzzz[k] = -g_0_y_yyyy_xzzzzzz[k] * ab_z + g_0_y_yyyy_xzzzzzzz[k];

                g_0_y_yyyyz_yyyyyyy[k] = -g_0_y_yyyy_yyyyyyy[k] * ab_z + g_0_y_yyyy_yyyyyyyz[k];

                g_0_y_yyyyz_yyyyyyz[k] = -g_0_y_yyyy_yyyyyyz[k] * ab_z + g_0_y_yyyy_yyyyyyzz[k];

                g_0_y_yyyyz_yyyyyzz[k] = -g_0_y_yyyy_yyyyyzz[k] * ab_z + g_0_y_yyyy_yyyyyzzz[k];

                g_0_y_yyyyz_yyyyzzz[k] = -g_0_y_yyyy_yyyyzzz[k] * ab_z + g_0_y_yyyy_yyyyzzzz[k];

                g_0_y_yyyyz_yyyzzzz[k] = -g_0_y_yyyy_yyyzzzz[k] * ab_z + g_0_y_yyyy_yyyzzzzz[k];

                g_0_y_yyyyz_yyzzzzz[k] = -g_0_y_yyyy_yyzzzzz[k] * ab_z + g_0_y_yyyy_yyzzzzzz[k];

                g_0_y_yyyyz_yzzzzzz[k] = -g_0_y_yyyy_yzzzzzz[k] * ab_z + g_0_y_yyyy_yzzzzzzz[k];

                g_0_y_yyyyz_zzzzzzz[k] = -g_0_y_yyyy_zzzzzzz[k] * ab_z + g_0_y_yyyy_zzzzzzzz[k];
            }

            /// Set up 1368-1404 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1368 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1369 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1370 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1371 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1372 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1373 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1374 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1375 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1376 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1377 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1378 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1379 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1380 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1381 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1382 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1383 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1384 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1385 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1386 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1387 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1388 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1389 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1390 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1391 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1392 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1393 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1394 * ccomps * dcomps);

            auto g_0_y_yyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1395 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1396 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1397 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1398 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1399 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1400 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1401 * ccomps * dcomps);

            auto g_0_y_yyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1402 * ccomps * dcomps);

            auto g_0_y_yyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1403 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyz_xxxxxxx, g_0_y_yyyz_xxxxxxxz, g_0_y_yyyz_xxxxxxy, g_0_y_yyyz_xxxxxxyz, g_0_y_yyyz_xxxxxxz, g_0_y_yyyz_xxxxxxzz, g_0_y_yyyz_xxxxxyy, g_0_y_yyyz_xxxxxyyz, g_0_y_yyyz_xxxxxyz, g_0_y_yyyz_xxxxxyzz, g_0_y_yyyz_xxxxxzz, g_0_y_yyyz_xxxxxzzz, g_0_y_yyyz_xxxxyyy, g_0_y_yyyz_xxxxyyyz, g_0_y_yyyz_xxxxyyz, g_0_y_yyyz_xxxxyyzz, g_0_y_yyyz_xxxxyzz, g_0_y_yyyz_xxxxyzzz, g_0_y_yyyz_xxxxzzz, g_0_y_yyyz_xxxxzzzz, g_0_y_yyyz_xxxyyyy, g_0_y_yyyz_xxxyyyyz, g_0_y_yyyz_xxxyyyz, g_0_y_yyyz_xxxyyyzz, g_0_y_yyyz_xxxyyzz, g_0_y_yyyz_xxxyyzzz, g_0_y_yyyz_xxxyzzz, g_0_y_yyyz_xxxyzzzz, g_0_y_yyyz_xxxzzzz, g_0_y_yyyz_xxxzzzzz, g_0_y_yyyz_xxyyyyy, g_0_y_yyyz_xxyyyyyz, g_0_y_yyyz_xxyyyyz, g_0_y_yyyz_xxyyyyzz, g_0_y_yyyz_xxyyyzz, g_0_y_yyyz_xxyyyzzz, g_0_y_yyyz_xxyyzzz, g_0_y_yyyz_xxyyzzzz, g_0_y_yyyz_xxyzzzz, g_0_y_yyyz_xxyzzzzz, g_0_y_yyyz_xxzzzzz, g_0_y_yyyz_xxzzzzzz, g_0_y_yyyz_xyyyyyy, g_0_y_yyyz_xyyyyyyz, g_0_y_yyyz_xyyyyyz, g_0_y_yyyz_xyyyyyzz, g_0_y_yyyz_xyyyyzz, g_0_y_yyyz_xyyyyzzz, g_0_y_yyyz_xyyyzzz, g_0_y_yyyz_xyyyzzzz, g_0_y_yyyz_xyyzzzz, g_0_y_yyyz_xyyzzzzz, g_0_y_yyyz_xyzzzzz, g_0_y_yyyz_xyzzzzzz, g_0_y_yyyz_xzzzzzz, g_0_y_yyyz_xzzzzzzz, g_0_y_yyyz_yyyyyyy, g_0_y_yyyz_yyyyyyyz, g_0_y_yyyz_yyyyyyz, g_0_y_yyyz_yyyyyyzz, g_0_y_yyyz_yyyyyzz, g_0_y_yyyz_yyyyyzzz, g_0_y_yyyz_yyyyzzz, g_0_y_yyyz_yyyyzzzz, g_0_y_yyyz_yyyzzzz, g_0_y_yyyz_yyyzzzzz, g_0_y_yyyz_yyzzzzz, g_0_y_yyyz_yyzzzzzz, g_0_y_yyyz_yzzzzzz, g_0_y_yyyz_yzzzzzzz, g_0_y_yyyz_zzzzzzz, g_0_y_yyyz_zzzzzzzz, g_0_y_yyyzz_xxxxxxx, g_0_y_yyyzz_xxxxxxy, g_0_y_yyyzz_xxxxxxz, g_0_y_yyyzz_xxxxxyy, g_0_y_yyyzz_xxxxxyz, g_0_y_yyyzz_xxxxxzz, g_0_y_yyyzz_xxxxyyy, g_0_y_yyyzz_xxxxyyz, g_0_y_yyyzz_xxxxyzz, g_0_y_yyyzz_xxxxzzz, g_0_y_yyyzz_xxxyyyy, g_0_y_yyyzz_xxxyyyz, g_0_y_yyyzz_xxxyyzz, g_0_y_yyyzz_xxxyzzz, g_0_y_yyyzz_xxxzzzz, g_0_y_yyyzz_xxyyyyy, g_0_y_yyyzz_xxyyyyz, g_0_y_yyyzz_xxyyyzz, g_0_y_yyyzz_xxyyzzz, g_0_y_yyyzz_xxyzzzz, g_0_y_yyyzz_xxzzzzz, g_0_y_yyyzz_xyyyyyy, g_0_y_yyyzz_xyyyyyz, g_0_y_yyyzz_xyyyyzz, g_0_y_yyyzz_xyyyzzz, g_0_y_yyyzz_xyyzzzz, g_0_y_yyyzz_xyzzzzz, g_0_y_yyyzz_xzzzzzz, g_0_y_yyyzz_yyyyyyy, g_0_y_yyyzz_yyyyyyz, g_0_y_yyyzz_yyyyyzz, g_0_y_yyyzz_yyyyzzz, g_0_y_yyyzz_yyyzzzz, g_0_y_yyyzz_yyzzzzz, g_0_y_yyyzz_yzzzzzz, g_0_y_yyyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzz_xxxxxxx[k] = -g_0_y_yyyz_xxxxxxx[k] * ab_z + g_0_y_yyyz_xxxxxxxz[k];

                g_0_y_yyyzz_xxxxxxy[k] = -g_0_y_yyyz_xxxxxxy[k] * ab_z + g_0_y_yyyz_xxxxxxyz[k];

                g_0_y_yyyzz_xxxxxxz[k] = -g_0_y_yyyz_xxxxxxz[k] * ab_z + g_0_y_yyyz_xxxxxxzz[k];

                g_0_y_yyyzz_xxxxxyy[k] = -g_0_y_yyyz_xxxxxyy[k] * ab_z + g_0_y_yyyz_xxxxxyyz[k];

                g_0_y_yyyzz_xxxxxyz[k] = -g_0_y_yyyz_xxxxxyz[k] * ab_z + g_0_y_yyyz_xxxxxyzz[k];

                g_0_y_yyyzz_xxxxxzz[k] = -g_0_y_yyyz_xxxxxzz[k] * ab_z + g_0_y_yyyz_xxxxxzzz[k];

                g_0_y_yyyzz_xxxxyyy[k] = -g_0_y_yyyz_xxxxyyy[k] * ab_z + g_0_y_yyyz_xxxxyyyz[k];

                g_0_y_yyyzz_xxxxyyz[k] = -g_0_y_yyyz_xxxxyyz[k] * ab_z + g_0_y_yyyz_xxxxyyzz[k];

                g_0_y_yyyzz_xxxxyzz[k] = -g_0_y_yyyz_xxxxyzz[k] * ab_z + g_0_y_yyyz_xxxxyzzz[k];

                g_0_y_yyyzz_xxxxzzz[k] = -g_0_y_yyyz_xxxxzzz[k] * ab_z + g_0_y_yyyz_xxxxzzzz[k];

                g_0_y_yyyzz_xxxyyyy[k] = -g_0_y_yyyz_xxxyyyy[k] * ab_z + g_0_y_yyyz_xxxyyyyz[k];

                g_0_y_yyyzz_xxxyyyz[k] = -g_0_y_yyyz_xxxyyyz[k] * ab_z + g_0_y_yyyz_xxxyyyzz[k];

                g_0_y_yyyzz_xxxyyzz[k] = -g_0_y_yyyz_xxxyyzz[k] * ab_z + g_0_y_yyyz_xxxyyzzz[k];

                g_0_y_yyyzz_xxxyzzz[k] = -g_0_y_yyyz_xxxyzzz[k] * ab_z + g_0_y_yyyz_xxxyzzzz[k];

                g_0_y_yyyzz_xxxzzzz[k] = -g_0_y_yyyz_xxxzzzz[k] * ab_z + g_0_y_yyyz_xxxzzzzz[k];

                g_0_y_yyyzz_xxyyyyy[k] = -g_0_y_yyyz_xxyyyyy[k] * ab_z + g_0_y_yyyz_xxyyyyyz[k];

                g_0_y_yyyzz_xxyyyyz[k] = -g_0_y_yyyz_xxyyyyz[k] * ab_z + g_0_y_yyyz_xxyyyyzz[k];

                g_0_y_yyyzz_xxyyyzz[k] = -g_0_y_yyyz_xxyyyzz[k] * ab_z + g_0_y_yyyz_xxyyyzzz[k];

                g_0_y_yyyzz_xxyyzzz[k] = -g_0_y_yyyz_xxyyzzz[k] * ab_z + g_0_y_yyyz_xxyyzzzz[k];

                g_0_y_yyyzz_xxyzzzz[k] = -g_0_y_yyyz_xxyzzzz[k] * ab_z + g_0_y_yyyz_xxyzzzzz[k];

                g_0_y_yyyzz_xxzzzzz[k] = -g_0_y_yyyz_xxzzzzz[k] * ab_z + g_0_y_yyyz_xxzzzzzz[k];

                g_0_y_yyyzz_xyyyyyy[k] = -g_0_y_yyyz_xyyyyyy[k] * ab_z + g_0_y_yyyz_xyyyyyyz[k];

                g_0_y_yyyzz_xyyyyyz[k] = -g_0_y_yyyz_xyyyyyz[k] * ab_z + g_0_y_yyyz_xyyyyyzz[k];

                g_0_y_yyyzz_xyyyyzz[k] = -g_0_y_yyyz_xyyyyzz[k] * ab_z + g_0_y_yyyz_xyyyyzzz[k];

                g_0_y_yyyzz_xyyyzzz[k] = -g_0_y_yyyz_xyyyzzz[k] * ab_z + g_0_y_yyyz_xyyyzzzz[k];

                g_0_y_yyyzz_xyyzzzz[k] = -g_0_y_yyyz_xyyzzzz[k] * ab_z + g_0_y_yyyz_xyyzzzzz[k];

                g_0_y_yyyzz_xyzzzzz[k] = -g_0_y_yyyz_xyzzzzz[k] * ab_z + g_0_y_yyyz_xyzzzzzz[k];

                g_0_y_yyyzz_xzzzzzz[k] = -g_0_y_yyyz_xzzzzzz[k] * ab_z + g_0_y_yyyz_xzzzzzzz[k];

                g_0_y_yyyzz_yyyyyyy[k] = -g_0_y_yyyz_yyyyyyy[k] * ab_z + g_0_y_yyyz_yyyyyyyz[k];

                g_0_y_yyyzz_yyyyyyz[k] = -g_0_y_yyyz_yyyyyyz[k] * ab_z + g_0_y_yyyz_yyyyyyzz[k];

                g_0_y_yyyzz_yyyyyzz[k] = -g_0_y_yyyz_yyyyyzz[k] * ab_z + g_0_y_yyyz_yyyyyzzz[k];

                g_0_y_yyyzz_yyyyzzz[k] = -g_0_y_yyyz_yyyyzzz[k] * ab_z + g_0_y_yyyz_yyyyzzzz[k];

                g_0_y_yyyzz_yyyzzzz[k] = -g_0_y_yyyz_yyyzzzz[k] * ab_z + g_0_y_yyyz_yyyzzzzz[k];

                g_0_y_yyyzz_yyzzzzz[k] = -g_0_y_yyyz_yyzzzzz[k] * ab_z + g_0_y_yyyz_yyzzzzzz[k];

                g_0_y_yyyzz_yzzzzzz[k] = -g_0_y_yyyz_yzzzzzz[k] * ab_z + g_0_y_yyyz_yzzzzzzz[k];

                g_0_y_yyyzz_zzzzzzz[k] = -g_0_y_yyyz_zzzzzzz[k] * ab_z + g_0_y_yyyz_zzzzzzzz[k];
            }

            /// Set up 1404-1440 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1404 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1405 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1406 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1407 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1408 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1409 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1410 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1411 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1412 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1413 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1414 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1415 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1416 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1417 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1418 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1419 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1420 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1421 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1422 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1423 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1424 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1425 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1426 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1427 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1428 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1429 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1430 * ccomps * dcomps);

            auto g_0_y_yyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1431 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1432 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1433 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1434 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1435 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1436 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1437 * ccomps * dcomps);

            auto g_0_y_yyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1438 * ccomps * dcomps);

            auto g_0_y_yyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1439 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzz_xxxxxxx, g_0_y_yyzz_xxxxxxxz, g_0_y_yyzz_xxxxxxy, g_0_y_yyzz_xxxxxxyz, g_0_y_yyzz_xxxxxxz, g_0_y_yyzz_xxxxxxzz, g_0_y_yyzz_xxxxxyy, g_0_y_yyzz_xxxxxyyz, g_0_y_yyzz_xxxxxyz, g_0_y_yyzz_xxxxxyzz, g_0_y_yyzz_xxxxxzz, g_0_y_yyzz_xxxxxzzz, g_0_y_yyzz_xxxxyyy, g_0_y_yyzz_xxxxyyyz, g_0_y_yyzz_xxxxyyz, g_0_y_yyzz_xxxxyyzz, g_0_y_yyzz_xxxxyzz, g_0_y_yyzz_xxxxyzzz, g_0_y_yyzz_xxxxzzz, g_0_y_yyzz_xxxxzzzz, g_0_y_yyzz_xxxyyyy, g_0_y_yyzz_xxxyyyyz, g_0_y_yyzz_xxxyyyz, g_0_y_yyzz_xxxyyyzz, g_0_y_yyzz_xxxyyzz, g_0_y_yyzz_xxxyyzzz, g_0_y_yyzz_xxxyzzz, g_0_y_yyzz_xxxyzzzz, g_0_y_yyzz_xxxzzzz, g_0_y_yyzz_xxxzzzzz, g_0_y_yyzz_xxyyyyy, g_0_y_yyzz_xxyyyyyz, g_0_y_yyzz_xxyyyyz, g_0_y_yyzz_xxyyyyzz, g_0_y_yyzz_xxyyyzz, g_0_y_yyzz_xxyyyzzz, g_0_y_yyzz_xxyyzzz, g_0_y_yyzz_xxyyzzzz, g_0_y_yyzz_xxyzzzz, g_0_y_yyzz_xxyzzzzz, g_0_y_yyzz_xxzzzzz, g_0_y_yyzz_xxzzzzzz, g_0_y_yyzz_xyyyyyy, g_0_y_yyzz_xyyyyyyz, g_0_y_yyzz_xyyyyyz, g_0_y_yyzz_xyyyyyzz, g_0_y_yyzz_xyyyyzz, g_0_y_yyzz_xyyyyzzz, g_0_y_yyzz_xyyyzzz, g_0_y_yyzz_xyyyzzzz, g_0_y_yyzz_xyyzzzz, g_0_y_yyzz_xyyzzzzz, g_0_y_yyzz_xyzzzzz, g_0_y_yyzz_xyzzzzzz, g_0_y_yyzz_xzzzzzz, g_0_y_yyzz_xzzzzzzz, g_0_y_yyzz_yyyyyyy, g_0_y_yyzz_yyyyyyyz, g_0_y_yyzz_yyyyyyz, g_0_y_yyzz_yyyyyyzz, g_0_y_yyzz_yyyyyzz, g_0_y_yyzz_yyyyyzzz, g_0_y_yyzz_yyyyzzz, g_0_y_yyzz_yyyyzzzz, g_0_y_yyzz_yyyzzzz, g_0_y_yyzz_yyyzzzzz, g_0_y_yyzz_yyzzzzz, g_0_y_yyzz_yyzzzzzz, g_0_y_yyzz_yzzzzzz, g_0_y_yyzz_yzzzzzzz, g_0_y_yyzz_zzzzzzz, g_0_y_yyzz_zzzzzzzz, g_0_y_yyzzz_xxxxxxx, g_0_y_yyzzz_xxxxxxy, g_0_y_yyzzz_xxxxxxz, g_0_y_yyzzz_xxxxxyy, g_0_y_yyzzz_xxxxxyz, g_0_y_yyzzz_xxxxxzz, g_0_y_yyzzz_xxxxyyy, g_0_y_yyzzz_xxxxyyz, g_0_y_yyzzz_xxxxyzz, g_0_y_yyzzz_xxxxzzz, g_0_y_yyzzz_xxxyyyy, g_0_y_yyzzz_xxxyyyz, g_0_y_yyzzz_xxxyyzz, g_0_y_yyzzz_xxxyzzz, g_0_y_yyzzz_xxxzzzz, g_0_y_yyzzz_xxyyyyy, g_0_y_yyzzz_xxyyyyz, g_0_y_yyzzz_xxyyyzz, g_0_y_yyzzz_xxyyzzz, g_0_y_yyzzz_xxyzzzz, g_0_y_yyzzz_xxzzzzz, g_0_y_yyzzz_xyyyyyy, g_0_y_yyzzz_xyyyyyz, g_0_y_yyzzz_xyyyyzz, g_0_y_yyzzz_xyyyzzz, g_0_y_yyzzz_xyyzzzz, g_0_y_yyzzz_xyzzzzz, g_0_y_yyzzz_xzzzzzz, g_0_y_yyzzz_yyyyyyy, g_0_y_yyzzz_yyyyyyz, g_0_y_yyzzz_yyyyyzz, g_0_y_yyzzz_yyyyzzz, g_0_y_yyzzz_yyyzzzz, g_0_y_yyzzz_yyzzzzz, g_0_y_yyzzz_yzzzzzz, g_0_y_yyzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzz_xxxxxxx[k] = -g_0_y_yyzz_xxxxxxx[k] * ab_z + g_0_y_yyzz_xxxxxxxz[k];

                g_0_y_yyzzz_xxxxxxy[k] = -g_0_y_yyzz_xxxxxxy[k] * ab_z + g_0_y_yyzz_xxxxxxyz[k];

                g_0_y_yyzzz_xxxxxxz[k] = -g_0_y_yyzz_xxxxxxz[k] * ab_z + g_0_y_yyzz_xxxxxxzz[k];

                g_0_y_yyzzz_xxxxxyy[k] = -g_0_y_yyzz_xxxxxyy[k] * ab_z + g_0_y_yyzz_xxxxxyyz[k];

                g_0_y_yyzzz_xxxxxyz[k] = -g_0_y_yyzz_xxxxxyz[k] * ab_z + g_0_y_yyzz_xxxxxyzz[k];

                g_0_y_yyzzz_xxxxxzz[k] = -g_0_y_yyzz_xxxxxzz[k] * ab_z + g_0_y_yyzz_xxxxxzzz[k];

                g_0_y_yyzzz_xxxxyyy[k] = -g_0_y_yyzz_xxxxyyy[k] * ab_z + g_0_y_yyzz_xxxxyyyz[k];

                g_0_y_yyzzz_xxxxyyz[k] = -g_0_y_yyzz_xxxxyyz[k] * ab_z + g_0_y_yyzz_xxxxyyzz[k];

                g_0_y_yyzzz_xxxxyzz[k] = -g_0_y_yyzz_xxxxyzz[k] * ab_z + g_0_y_yyzz_xxxxyzzz[k];

                g_0_y_yyzzz_xxxxzzz[k] = -g_0_y_yyzz_xxxxzzz[k] * ab_z + g_0_y_yyzz_xxxxzzzz[k];

                g_0_y_yyzzz_xxxyyyy[k] = -g_0_y_yyzz_xxxyyyy[k] * ab_z + g_0_y_yyzz_xxxyyyyz[k];

                g_0_y_yyzzz_xxxyyyz[k] = -g_0_y_yyzz_xxxyyyz[k] * ab_z + g_0_y_yyzz_xxxyyyzz[k];

                g_0_y_yyzzz_xxxyyzz[k] = -g_0_y_yyzz_xxxyyzz[k] * ab_z + g_0_y_yyzz_xxxyyzzz[k];

                g_0_y_yyzzz_xxxyzzz[k] = -g_0_y_yyzz_xxxyzzz[k] * ab_z + g_0_y_yyzz_xxxyzzzz[k];

                g_0_y_yyzzz_xxxzzzz[k] = -g_0_y_yyzz_xxxzzzz[k] * ab_z + g_0_y_yyzz_xxxzzzzz[k];

                g_0_y_yyzzz_xxyyyyy[k] = -g_0_y_yyzz_xxyyyyy[k] * ab_z + g_0_y_yyzz_xxyyyyyz[k];

                g_0_y_yyzzz_xxyyyyz[k] = -g_0_y_yyzz_xxyyyyz[k] * ab_z + g_0_y_yyzz_xxyyyyzz[k];

                g_0_y_yyzzz_xxyyyzz[k] = -g_0_y_yyzz_xxyyyzz[k] * ab_z + g_0_y_yyzz_xxyyyzzz[k];

                g_0_y_yyzzz_xxyyzzz[k] = -g_0_y_yyzz_xxyyzzz[k] * ab_z + g_0_y_yyzz_xxyyzzzz[k];

                g_0_y_yyzzz_xxyzzzz[k] = -g_0_y_yyzz_xxyzzzz[k] * ab_z + g_0_y_yyzz_xxyzzzzz[k];

                g_0_y_yyzzz_xxzzzzz[k] = -g_0_y_yyzz_xxzzzzz[k] * ab_z + g_0_y_yyzz_xxzzzzzz[k];

                g_0_y_yyzzz_xyyyyyy[k] = -g_0_y_yyzz_xyyyyyy[k] * ab_z + g_0_y_yyzz_xyyyyyyz[k];

                g_0_y_yyzzz_xyyyyyz[k] = -g_0_y_yyzz_xyyyyyz[k] * ab_z + g_0_y_yyzz_xyyyyyzz[k];

                g_0_y_yyzzz_xyyyyzz[k] = -g_0_y_yyzz_xyyyyzz[k] * ab_z + g_0_y_yyzz_xyyyyzzz[k];

                g_0_y_yyzzz_xyyyzzz[k] = -g_0_y_yyzz_xyyyzzz[k] * ab_z + g_0_y_yyzz_xyyyzzzz[k];

                g_0_y_yyzzz_xyyzzzz[k] = -g_0_y_yyzz_xyyzzzz[k] * ab_z + g_0_y_yyzz_xyyzzzzz[k];

                g_0_y_yyzzz_xyzzzzz[k] = -g_0_y_yyzz_xyzzzzz[k] * ab_z + g_0_y_yyzz_xyzzzzzz[k];

                g_0_y_yyzzz_xzzzzzz[k] = -g_0_y_yyzz_xzzzzzz[k] * ab_z + g_0_y_yyzz_xzzzzzzz[k];

                g_0_y_yyzzz_yyyyyyy[k] = -g_0_y_yyzz_yyyyyyy[k] * ab_z + g_0_y_yyzz_yyyyyyyz[k];

                g_0_y_yyzzz_yyyyyyz[k] = -g_0_y_yyzz_yyyyyyz[k] * ab_z + g_0_y_yyzz_yyyyyyzz[k];

                g_0_y_yyzzz_yyyyyzz[k] = -g_0_y_yyzz_yyyyyzz[k] * ab_z + g_0_y_yyzz_yyyyyzzz[k];

                g_0_y_yyzzz_yyyyzzz[k] = -g_0_y_yyzz_yyyyzzz[k] * ab_z + g_0_y_yyzz_yyyyzzzz[k];

                g_0_y_yyzzz_yyyzzzz[k] = -g_0_y_yyzz_yyyzzzz[k] * ab_z + g_0_y_yyzz_yyyzzzzz[k];

                g_0_y_yyzzz_yyzzzzz[k] = -g_0_y_yyzz_yyzzzzz[k] * ab_z + g_0_y_yyzz_yyzzzzzz[k];

                g_0_y_yyzzz_yzzzzzz[k] = -g_0_y_yyzz_yzzzzzz[k] * ab_z + g_0_y_yyzz_yzzzzzzz[k];

                g_0_y_yyzzz_zzzzzzz[k] = -g_0_y_yyzz_zzzzzzz[k] * ab_z + g_0_y_yyzz_zzzzzzzz[k];
            }

            /// Set up 1440-1476 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1440 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1441 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1442 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1443 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1444 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1445 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1446 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1447 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1448 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1449 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1450 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1451 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1452 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1453 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1454 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1455 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1456 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1457 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1458 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1459 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1460 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1461 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1462 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1463 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1464 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1465 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1466 * ccomps * dcomps);

            auto g_0_y_yzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1467 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1468 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1469 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1470 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1471 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1472 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1473 * ccomps * dcomps);

            auto g_0_y_yzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1474 * ccomps * dcomps);

            auto g_0_y_yzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1475 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzz_xxxxxxx, g_0_y_yzzz_xxxxxxxz, g_0_y_yzzz_xxxxxxy, g_0_y_yzzz_xxxxxxyz, g_0_y_yzzz_xxxxxxz, g_0_y_yzzz_xxxxxxzz, g_0_y_yzzz_xxxxxyy, g_0_y_yzzz_xxxxxyyz, g_0_y_yzzz_xxxxxyz, g_0_y_yzzz_xxxxxyzz, g_0_y_yzzz_xxxxxzz, g_0_y_yzzz_xxxxxzzz, g_0_y_yzzz_xxxxyyy, g_0_y_yzzz_xxxxyyyz, g_0_y_yzzz_xxxxyyz, g_0_y_yzzz_xxxxyyzz, g_0_y_yzzz_xxxxyzz, g_0_y_yzzz_xxxxyzzz, g_0_y_yzzz_xxxxzzz, g_0_y_yzzz_xxxxzzzz, g_0_y_yzzz_xxxyyyy, g_0_y_yzzz_xxxyyyyz, g_0_y_yzzz_xxxyyyz, g_0_y_yzzz_xxxyyyzz, g_0_y_yzzz_xxxyyzz, g_0_y_yzzz_xxxyyzzz, g_0_y_yzzz_xxxyzzz, g_0_y_yzzz_xxxyzzzz, g_0_y_yzzz_xxxzzzz, g_0_y_yzzz_xxxzzzzz, g_0_y_yzzz_xxyyyyy, g_0_y_yzzz_xxyyyyyz, g_0_y_yzzz_xxyyyyz, g_0_y_yzzz_xxyyyyzz, g_0_y_yzzz_xxyyyzz, g_0_y_yzzz_xxyyyzzz, g_0_y_yzzz_xxyyzzz, g_0_y_yzzz_xxyyzzzz, g_0_y_yzzz_xxyzzzz, g_0_y_yzzz_xxyzzzzz, g_0_y_yzzz_xxzzzzz, g_0_y_yzzz_xxzzzzzz, g_0_y_yzzz_xyyyyyy, g_0_y_yzzz_xyyyyyyz, g_0_y_yzzz_xyyyyyz, g_0_y_yzzz_xyyyyyzz, g_0_y_yzzz_xyyyyzz, g_0_y_yzzz_xyyyyzzz, g_0_y_yzzz_xyyyzzz, g_0_y_yzzz_xyyyzzzz, g_0_y_yzzz_xyyzzzz, g_0_y_yzzz_xyyzzzzz, g_0_y_yzzz_xyzzzzz, g_0_y_yzzz_xyzzzzzz, g_0_y_yzzz_xzzzzzz, g_0_y_yzzz_xzzzzzzz, g_0_y_yzzz_yyyyyyy, g_0_y_yzzz_yyyyyyyz, g_0_y_yzzz_yyyyyyz, g_0_y_yzzz_yyyyyyzz, g_0_y_yzzz_yyyyyzz, g_0_y_yzzz_yyyyyzzz, g_0_y_yzzz_yyyyzzz, g_0_y_yzzz_yyyyzzzz, g_0_y_yzzz_yyyzzzz, g_0_y_yzzz_yyyzzzzz, g_0_y_yzzz_yyzzzzz, g_0_y_yzzz_yyzzzzzz, g_0_y_yzzz_yzzzzzz, g_0_y_yzzz_yzzzzzzz, g_0_y_yzzz_zzzzzzz, g_0_y_yzzz_zzzzzzzz, g_0_y_yzzzz_xxxxxxx, g_0_y_yzzzz_xxxxxxy, g_0_y_yzzzz_xxxxxxz, g_0_y_yzzzz_xxxxxyy, g_0_y_yzzzz_xxxxxyz, g_0_y_yzzzz_xxxxxzz, g_0_y_yzzzz_xxxxyyy, g_0_y_yzzzz_xxxxyyz, g_0_y_yzzzz_xxxxyzz, g_0_y_yzzzz_xxxxzzz, g_0_y_yzzzz_xxxyyyy, g_0_y_yzzzz_xxxyyyz, g_0_y_yzzzz_xxxyyzz, g_0_y_yzzzz_xxxyzzz, g_0_y_yzzzz_xxxzzzz, g_0_y_yzzzz_xxyyyyy, g_0_y_yzzzz_xxyyyyz, g_0_y_yzzzz_xxyyyzz, g_0_y_yzzzz_xxyyzzz, g_0_y_yzzzz_xxyzzzz, g_0_y_yzzzz_xxzzzzz, g_0_y_yzzzz_xyyyyyy, g_0_y_yzzzz_xyyyyyz, g_0_y_yzzzz_xyyyyzz, g_0_y_yzzzz_xyyyzzz, g_0_y_yzzzz_xyyzzzz, g_0_y_yzzzz_xyzzzzz, g_0_y_yzzzz_xzzzzzz, g_0_y_yzzzz_yyyyyyy, g_0_y_yzzzz_yyyyyyz, g_0_y_yzzzz_yyyyyzz, g_0_y_yzzzz_yyyyzzz, g_0_y_yzzzz_yyyzzzz, g_0_y_yzzzz_yyzzzzz, g_0_y_yzzzz_yzzzzzz, g_0_y_yzzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzz_xxxxxxx[k] = -g_0_y_yzzz_xxxxxxx[k] * ab_z + g_0_y_yzzz_xxxxxxxz[k];

                g_0_y_yzzzz_xxxxxxy[k] = -g_0_y_yzzz_xxxxxxy[k] * ab_z + g_0_y_yzzz_xxxxxxyz[k];

                g_0_y_yzzzz_xxxxxxz[k] = -g_0_y_yzzz_xxxxxxz[k] * ab_z + g_0_y_yzzz_xxxxxxzz[k];

                g_0_y_yzzzz_xxxxxyy[k] = -g_0_y_yzzz_xxxxxyy[k] * ab_z + g_0_y_yzzz_xxxxxyyz[k];

                g_0_y_yzzzz_xxxxxyz[k] = -g_0_y_yzzz_xxxxxyz[k] * ab_z + g_0_y_yzzz_xxxxxyzz[k];

                g_0_y_yzzzz_xxxxxzz[k] = -g_0_y_yzzz_xxxxxzz[k] * ab_z + g_0_y_yzzz_xxxxxzzz[k];

                g_0_y_yzzzz_xxxxyyy[k] = -g_0_y_yzzz_xxxxyyy[k] * ab_z + g_0_y_yzzz_xxxxyyyz[k];

                g_0_y_yzzzz_xxxxyyz[k] = -g_0_y_yzzz_xxxxyyz[k] * ab_z + g_0_y_yzzz_xxxxyyzz[k];

                g_0_y_yzzzz_xxxxyzz[k] = -g_0_y_yzzz_xxxxyzz[k] * ab_z + g_0_y_yzzz_xxxxyzzz[k];

                g_0_y_yzzzz_xxxxzzz[k] = -g_0_y_yzzz_xxxxzzz[k] * ab_z + g_0_y_yzzz_xxxxzzzz[k];

                g_0_y_yzzzz_xxxyyyy[k] = -g_0_y_yzzz_xxxyyyy[k] * ab_z + g_0_y_yzzz_xxxyyyyz[k];

                g_0_y_yzzzz_xxxyyyz[k] = -g_0_y_yzzz_xxxyyyz[k] * ab_z + g_0_y_yzzz_xxxyyyzz[k];

                g_0_y_yzzzz_xxxyyzz[k] = -g_0_y_yzzz_xxxyyzz[k] * ab_z + g_0_y_yzzz_xxxyyzzz[k];

                g_0_y_yzzzz_xxxyzzz[k] = -g_0_y_yzzz_xxxyzzz[k] * ab_z + g_0_y_yzzz_xxxyzzzz[k];

                g_0_y_yzzzz_xxxzzzz[k] = -g_0_y_yzzz_xxxzzzz[k] * ab_z + g_0_y_yzzz_xxxzzzzz[k];

                g_0_y_yzzzz_xxyyyyy[k] = -g_0_y_yzzz_xxyyyyy[k] * ab_z + g_0_y_yzzz_xxyyyyyz[k];

                g_0_y_yzzzz_xxyyyyz[k] = -g_0_y_yzzz_xxyyyyz[k] * ab_z + g_0_y_yzzz_xxyyyyzz[k];

                g_0_y_yzzzz_xxyyyzz[k] = -g_0_y_yzzz_xxyyyzz[k] * ab_z + g_0_y_yzzz_xxyyyzzz[k];

                g_0_y_yzzzz_xxyyzzz[k] = -g_0_y_yzzz_xxyyzzz[k] * ab_z + g_0_y_yzzz_xxyyzzzz[k];

                g_0_y_yzzzz_xxyzzzz[k] = -g_0_y_yzzz_xxyzzzz[k] * ab_z + g_0_y_yzzz_xxyzzzzz[k];

                g_0_y_yzzzz_xxzzzzz[k] = -g_0_y_yzzz_xxzzzzz[k] * ab_z + g_0_y_yzzz_xxzzzzzz[k];

                g_0_y_yzzzz_xyyyyyy[k] = -g_0_y_yzzz_xyyyyyy[k] * ab_z + g_0_y_yzzz_xyyyyyyz[k];

                g_0_y_yzzzz_xyyyyyz[k] = -g_0_y_yzzz_xyyyyyz[k] * ab_z + g_0_y_yzzz_xyyyyyzz[k];

                g_0_y_yzzzz_xyyyyzz[k] = -g_0_y_yzzz_xyyyyzz[k] * ab_z + g_0_y_yzzz_xyyyyzzz[k];

                g_0_y_yzzzz_xyyyzzz[k] = -g_0_y_yzzz_xyyyzzz[k] * ab_z + g_0_y_yzzz_xyyyzzzz[k];

                g_0_y_yzzzz_xyyzzzz[k] = -g_0_y_yzzz_xyyzzzz[k] * ab_z + g_0_y_yzzz_xyyzzzzz[k];

                g_0_y_yzzzz_xyzzzzz[k] = -g_0_y_yzzz_xyzzzzz[k] * ab_z + g_0_y_yzzz_xyzzzzzz[k];

                g_0_y_yzzzz_xzzzzzz[k] = -g_0_y_yzzz_xzzzzzz[k] * ab_z + g_0_y_yzzz_xzzzzzzz[k];

                g_0_y_yzzzz_yyyyyyy[k] = -g_0_y_yzzz_yyyyyyy[k] * ab_z + g_0_y_yzzz_yyyyyyyz[k];

                g_0_y_yzzzz_yyyyyyz[k] = -g_0_y_yzzz_yyyyyyz[k] * ab_z + g_0_y_yzzz_yyyyyyzz[k];

                g_0_y_yzzzz_yyyyyzz[k] = -g_0_y_yzzz_yyyyyzz[k] * ab_z + g_0_y_yzzz_yyyyyzzz[k];

                g_0_y_yzzzz_yyyyzzz[k] = -g_0_y_yzzz_yyyyzzz[k] * ab_z + g_0_y_yzzz_yyyyzzzz[k];

                g_0_y_yzzzz_yyyzzzz[k] = -g_0_y_yzzz_yyyzzzz[k] * ab_z + g_0_y_yzzz_yyyzzzzz[k];

                g_0_y_yzzzz_yyzzzzz[k] = -g_0_y_yzzz_yyzzzzz[k] * ab_z + g_0_y_yzzz_yyzzzzzz[k];

                g_0_y_yzzzz_yzzzzzz[k] = -g_0_y_yzzz_yzzzzzz[k] * ab_z + g_0_y_yzzz_yzzzzzzz[k];

                g_0_y_yzzzz_zzzzzzz[k] = -g_0_y_yzzz_zzzzzzz[k] * ab_z + g_0_y_yzzz_zzzzzzzz[k];
            }

            /// Set up 1476-1512 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1476 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1477 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1478 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1479 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1480 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1481 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1482 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1483 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1484 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1485 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1486 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1487 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1488 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1489 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1490 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1491 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1492 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1493 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1494 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1495 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1496 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1497 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1498 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1499 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1500 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1501 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1502 * ccomps * dcomps);

            auto g_0_y_zzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1503 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1504 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1505 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1506 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1507 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1508 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1509 * ccomps * dcomps);

            auto g_0_y_zzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1510 * ccomps * dcomps);

            auto g_0_y_zzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_xxxxxxx, g_0_y_zzzz_xxxxxxxz, g_0_y_zzzz_xxxxxxy, g_0_y_zzzz_xxxxxxyz, g_0_y_zzzz_xxxxxxz, g_0_y_zzzz_xxxxxxzz, g_0_y_zzzz_xxxxxyy, g_0_y_zzzz_xxxxxyyz, g_0_y_zzzz_xxxxxyz, g_0_y_zzzz_xxxxxyzz, g_0_y_zzzz_xxxxxzz, g_0_y_zzzz_xxxxxzzz, g_0_y_zzzz_xxxxyyy, g_0_y_zzzz_xxxxyyyz, g_0_y_zzzz_xxxxyyz, g_0_y_zzzz_xxxxyyzz, g_0_y_zzzz_xxxxyzz, g_0_y_zzzz_xxxxyzzz, g_0_y_zzzz_xxxxzzz, g_0_y_zzzz_xxxxzzzz, g_0_y_zzzz_xxxyyyy, g_0_y_zzzz_xxxyyyyz, g_0_y_zzzz_xxxyyyz, g_0_y_zzzz_xxxyyyzz, g_0_y_zzzz_xxxyyzz, g_0_y_zzzz_xxxyyzzz, g_0_y_zzzz_xxxyzzz, g_0_y_zzzz_xxxyzzzz, g_0_y_zzzz_xxxzzzz, g_0_y_zzzz_xxxzzzzz, g_0_y_zzzz_xxyyyyy, g_0_y_zzzz_xxyyyyyz, g_0_y_zzzz_xxyyyyz, g_0_y_zzzz_xxyyyyzz, g_0_y_zzzz_xxyyyzz, g_0_y_zzzz_xxyyyzzz, g_0_y_zzzz_xxyyzzz, g_0_y_zzzz_xxyyzzzz, g_0_y_zzzz_xxyzzzz, g_0_y_zzzz_xxyzzzzz, g_0_y_zzzz_xxzzzzz, g_0_y_zzzz_xxzzzzzz, g_0_y_zzzz_xyyyyyy, g_0_y_zzzz_xyyyyyyz, g_0_y_zzzz_xyyyyyz, g_0_y_zzzz_xyyyyyzz, g_0_y_zzzz_xyyyyzz, g_0_y_zzzz_xyyyyzzz, g_0_y_zzzz_xyyyzzz, g_0_y_zzzz_xyyyzzzz, g_0_y_zzzz_xyyzzzz, g_0_y_zzzz_xyyzzzzz, g_0_y_zzzz_xyzzzzz, g_0_y_zzzz_xyzzzzzz, g_0_y_zzzz_xzzzzzz, g_0_y_zzzz_xzzzzzzz, g_0_y_zzzz_yyyyyyy, g_0_y_zzzz_yyyyyyyz, g_0_y_zzzz_yyyyyyz, g_0_y_zzzz_yyyyyyzz, g_0_y_zzzz_yyyyyzz, g_0_y_zzzz_yyyyyzzz, g_0_y_zzzz_yyyyzzz, g_0_y_zzzz_yyyyzzzz, g_0_y_zzzz_yyyzzzz, g_0_y_zzzz_yyyzzzzz, g_0_y_zzzz_yyzzzzz, g_0_y_zzzz_yyzzzzzz, g_0_y_zzzz_yzzzzzz, g_0_y_zzzz_yzzzzzzz, g_0_y_zzzz_zzzzzzz, g_0_y_zzzz_zzzzzzzz, g_0_y_zzzzz_xxxxxxx, g_0_y_zzzzz_xxxxxxy, g_0_y_zzzzz_xxxxxxz, g_0_y_zzzzz_xxxxxyy, g_0_y_zzzzz_xxxxxyz, g_0_y_zzzzz_xxxxxzz, g_0_y_zzzzz_xxxxyyy, g_0_y_zzzzz_xxxxyyz, g_0_y_zzzzz_xxxxyzz, g_0_y_zzzzz_xxxxzzz, g_0_y_zzzzz_xxxyyyy, g_0_y_zzzzz_xxxyyyz, g_0_y_zzzzz_xxxyyzz, g_0_y_zzzzz_xxxyzzz, g_0_y_zzzzz_xxxzzzz, g_0_y_zzzzz_xxyyyyy, g_0_y_zzzzz_xxyyyyz, g_0_y_zzzzz_xxyyyzz, g_0_y_zzzzz_xxyyzzz, g_0_y_zzzzz_xxyzzzz, g_0_y_zzzzz_xxzzzzz, g_0_y_zzzzz_xyyyyyy, g_0_y_zzzzz_xyyyyyz, g_0_y_zzzzz_xyyyyzz, g_0_y_zzzzz_xyyyzzz, g_0_y_zzzzz_xyyzzzz, g_0_y_zzzzz_xyzzzzz, g_0_y_zzzzz_xzzzzzz, g_0_y_zzzzz_yyyyyyy, g_0_y_zzzzz_yyyyyyz, g_0_y_zzzzz_yyyyyzz, g_0_y_zzzzz_yyyyzzz, g_0_y_zzzzz_yyyzzzz, g_0_y_zzzzz_yyzzzzz, g_0_y_zzzzz_yzzzzzz, g_0_y_zzzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzz_xxxxxxx[k] = -g_0_y_zzzz_xxxxxxx[k] * ab_z + g_0_y_zzzz_xxxxxxxz[k];

                g_0_y_zzzzz_xxxxxxy[k] = -g_0_y_zzzz_xxxxxxy[k] * ab_z + g_0_y_zzzz_xxxxxxyz[k];

                g_0_y_zzzzz_xxxxxxz[k] = -g_0_y_zzzz_xxxxxxz[k] * ab_z + g_0_y_zzzz_xxxxxxzz[k];

                g_0_y_zzzzz_xxxxxyy[k] = -g_0_y_zzzz_xxxxxyy[k] * ab_z + g_0_y_zzzz_xxxxxyyz[k];

                g_0_y_zzzzz_xxxxxyz[k] = -g_0_y_zzzz_xxxxxyz[k] * ab_z + g_0_y_zzzz_xxxxxyzz[k];

                g_0_y_zzzzz_xxxxxzz[k] = -g_0_y_zzzz_xxxxxzz[k] * ab_z + g_0_y_zzzz_xxxxxzzz[k];

                g_0_y_zzzzz_xxxxyyy[k] = -g_0_y_zzzz_xxxxyyy[k] * ab_z + g_0_y_zzzz_xxxxyyyz[k];

                g_0_y_zzzzz_xxxxyyz[k] = -g_0_y_zzzz_xxxxyyz[k] * ab_z + g_0_y_zzzz_xxxxyyzz[k];

                g_0_y_zzzzz_xxxxyzz[k] = -g_0_y_zzzz_xxxxyzz[k] * ab_z + g_0_y_zzzz_xxxxyzzz[k];

                g_0_y_zzzzz_xxxxzzz[k] = -g_0_y_zzzz_xxxxzzz[k] * ab_z + g_0_y_zzzz_xxxxzzzz[k];

                g_0_y_zzzzz_xxxyyyy[k] = -g_0_y_zzzz_xxxyyyy[k] * ab_z + g_0_y_zzzz_xxxyyyyz[k];

                g_0_y_zzzzz_xxxyyyz[k] = -g_0_y_zzzz_xxxyyyz[k] * ab_z + g_0_y_zzzz_xxxyyyzz[k];

                g_0_y_zzzzz_xxxyyzz[k] = -g_0_y_zzzz_xxxyyzz[k] * ab_z + g_0_y_zzzz_xxxyyzzz[k];

                g_0_y_zzzzz_xxxyzzz[k] = -g_0_y_zzzz_xxxyzzz[k] * ab_z + g_0_y_zzzz_xxxyzzzz[k];

                g_0_y_zzzzz_xxxzzzz[k] = -g_0_y_zzzz_xxxzzzz[k] * ab_z + g_0_y_zzzz_xxxzzzzz[k];

                g_0_y_zzzzz_xxyyyyy[k] = -g_0_y_zzzz_xxyyyyy[k] * ab_z + g_0_y_zzzz_xxyyyyyz[k];

                g_0_y_zzzzz_xxyyyyz[k] = -g_0_y_zzzz_xxyyyyz[k] * ab_z + g_0_y_zzzz_xxyyyyzz[k];

                g_0_y_zzzzz_xxyyyzz[k] = -g_0_y_zzzz_xxyyyzz[k] * ab_z + g_0_y_zzzz_xxyyyzzz[k];

                g_0_y_zzzzz_xxyyzzz[k] = -g_0_y_zzzz_xxyyzzz[k] * ab_z + g_0_y_zzzz_xxyyzzzz[k];

                g_0_y_zzzzz_xxyzzzz[k] = -g_0_y_zzzz_xxyzzzz[k] * ab_z + g_0_y_zzzz_xxyzzzzz[k];

                g_0_y_zzzzz_xxzzzzz[k] = -g_0_y_zzzz_xxzzzzz[k] * ab_z + g_0_y_zzzz_xxzzzzzz[k];

                g_0_y_zzzzz_xyyyyyy[k] = -g_0_y_zzzz_xyyyyyy[k] * ab_z + g_0_y_zzzz_xyyyyyyz[k];

                g_0_y_zzzzz_xyyyyyz[k] = -g_0_y_zzzz_xyyyyyz[k] * ab_z + g_0_y_zzzz_xyyyyyzz[k];

                g_0_y_zzzzz_xyyyyzz[k] = -g_0_y_zzzz_xyyyyzz[k] * ab_z + g_0_y_zzzz_xyyyyzzz[k];

                g_0_y_zzzzz_xyyyzzz[k] = -g_0_y_zzzz_xyyyzzz[k] * ab_z + g_0_y_zzzz_xyyyzzzz[k];

                g_0_y_zzzzz_xyyzzzz[k] = -g_0_y_zzzz_xyyzzzz[k] * ab_z + g_0_y_zzzz_xyyzzzzz[k];

                g_0_y_zzzzz_xyzzzzz[k] = -g_0_y_zzzz_xyzzzzz[k] * ab_z + g_0_y_zzzz_xyzzzzzz[k];

                g_0_y_zzzzz_xzzzzzz[k] = -g_0_y_zzzz_xzzzzzz[k] * ab_z + g_0_y_zzzz_xzzzzzzz[k];

                g_0_y_zzzzz_yyyyyyy[k] = -g_0_y_zzzz_yyyyyyy[k] * ab_z + g_0_y_zzzz_yyyyyyyz[k];

                g_0_y_zzzzz_yyyyyyz[k] = -g_0_y_zzzz_yyyyyyz[k] * ab_z + g_0_y_zzzz_yyyyyyzz[k];

                g_0_y_zzzzz_yyyyyzz[k] = -g_0_y_zzzz_yyyyyzz[k] * ab_z + g_0_y_zzzz_yyyyyzzz[k];

                g_0_y_zzzzz_yyyyzzz[k] = -g_0_y_zzzz_yyyyzzz[k] * ab_z + g_0_y_zzzz_yyyyzzzz[k];

                g_0_y_zzzzz_yyyzzzz[k] = -g_0_y_zzzz_yyyzzzz[k] * ab_z + g_0_y_zzzz_yyyzzzzz[k];

                g_0_y_zzzzz_yyzzzzz[k] = -g_0_y_zzzz_yyzzzzz[k] * ab_z + g_0_y_zzzz_yyzzzzzz[k];

                g_0_y_zzzzz_yzzzzzz[k] = -g_0_y_zzzz_yzzzzzz[k] * ab_z + g_0_y_zzzz_yzzzzzzz[k];

                g_0_y_zzzzz_zzzzzzz[k] = -g_0_y_zzzz_zzzzzzz[k] * ab_z + g_0_y_zzzz_zzzzzzzz[k];
            }

            /// Set up 1512-1548 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxx_xxxxxxx = cbuffer.data(hk_geom_01_off + 1512 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxxxy = cbuffer.data(hk_geom_01_off + 1513 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxxxz = cbuffer.data(hk_geom_01_off + 1514 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxxyy = cbuffer.data(hk_geom_01_off + 1515 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxxyz = cbuffer.data(hk_geom_01_off + 1516 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxxzz = cbuffer.data(hk_geom_01_off + 1517 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxyyy = cbuffer.data(hk_geom_01_off + 1518 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxyyz = cbuffer.data(hk_geom_01_off + 1519 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxyzz = cbuffer.data(hk_geom_01_off + 1520 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxzzz = cbuffer.data(hk_geom_01_off + 1521 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyyyy = cbuffer.data(hk_geom_01_off + 1522 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyyyz = cbuffer.data(hk_geom_01_off + 1523 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyyzz = cbuffer.data(hk_geom_01_off + 1524 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyzzz = cbuffer.data(hk_geom_01_off + 1525 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxzzzz = cbuffer.data(hk_geom_01_off + 1526 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyyyy = cbuffer.data(hk_geom_01_off + 1527 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyyyz = cbuffer.data(hk_geom_01_off + 1528 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyyzz = cbuffer.data(hk_geom_01_off + 1529 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyzzz = cbuffer.data(hk_geom_01_off + 1530 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyzzzz = cbuffer.data(hk_geom_01_off + 1531 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxzzzzz = cbuffer.data(hk_geom_01_off + 1532 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyyyy = cbuffer.data(hk_geom_01_off + 1533 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyyyz = cbuffer.data(hk_geom_01_off + 1534 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyyzz = cbuffer.data(hk_geom_01_off + 1535 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyzzz = cbuffer.data(hk_geom_01_off + 1536 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyzzzz = cbuffer.data(hk_geom_01_off + 1537 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyzzzzz = cbuffer.data(hk_geom_01_off + 1538 * ccomps * dcomps);

            auto g_0_z_xxxxx_xzzzzzz = cbuffer.data(hk_geom_01_off + 1539 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyyyy = cbuffer.data(hk_geom_01_off + 1540 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyyyz = cbuffer.data(hk_geom_01_off + 1541 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyyzz = cbuffer.data(hk_geom_01_off + 1542 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyzzz = cbuffer.data(hk_geom_01_off + 1543 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyzzzz = cbuffer.data(hk_geom_01_off + 1544 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyzzzzz = cbuffer.data(hk_geom_01_off + 1545 * ccomps * dcomps);

            auto g_0_z_xxxxx_yzzzzzz = cbuffer.data(hk_geom_01_off + 1546 * ccomps * dcomps);

            auto g_0_z_xxxxx_zzzzzzz = cbuffer.data(hk_geom_01_off + 1547 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_xxxxxxx, g_0_z_xxxx_xxxxxxxx, g_0_z_xxxx_xxxxxxxy, g_0_z_xxxx_xxxxxxxz, g_0_z_xxxx_xxxxxxy, g_0_z_xxxx_xxxxxxyy, g_0_z_xxxx_xxxxxxyz, g_0_z_xxxx_xxxxxxz, g_0_z_xxxx_xxxxxxzz, g_0_z_xxxx_xxxxxyy, g_0_z_xxxx_xxxxxyyy, g_0_z_xxxx_xxxxxyyz, g_0_z_xxxx_xxxxxyz, g_0_z_xxxx_xxxxxyzz, g_0_z_xxxx_xxxxxzz, g_0_z_xxxx_xxxxxzzz, g_0_z_xxxx_xxxxyyy, g_0_z_xxxx_xxxxyyyy, g_0_z_xxxx_xxxxyyyz, g_0_z_xxxx_xxxxyyz, g_0_z_xxxx_xxxxyyzz, g_0_z_xxxx_xxxxyzz, g_0_z_xxxx_xxxxyzzz, g_0_z_xxxx_xxxxzzz, g_0_z_xxxx_xxxxzzzz, g_0_z_xxxx_xxxyyyy, g_0_z_xxxx_xxxyyyyy, g_0_z_xxxx_xxxyyyyz, g_0_z_xxxx_xxxyyyz, g_0_z_xxxx_xxxyyyzz, g_0_z_xxxx_xxxyyzz, g_0_z_xxxx_xxxyyzzz, g_0_z_xxxx_xxxyzzz, g_0_z_xxxx_xxxyzzzz, g_0_z_xxxx_xxxzzzz, g_0_z_xxxx_xxxzzzzz, g_0_z_xxxx_xxyyyyy, g_0_z_xxxx_xxyyyyyy, g_0_z_xxxx_xxyyyyyz, g_0_z_xxxx_xxyyyyz, g_0_z_xxxx_xxyyyyzz, g_0_z_xxxx_xxyyyzz, g_0_z_xxxx_xxyyyzzz, g_0_z_xxxx_xxyyzzz, g_0_z_xxxx_xxyyzzzz, g_0_z_xxxx_xxyzzzz, g_0_z_xxxx_xxyzzzzz, g_0_z_xxxx_xxzzzzz, g_0_z_xxxx_xxzzzzzz, g_0_z_xxxx_xyyyyyy, g_0_z_xxxx_xyyyyyyy, g_0_z_xxxx_xyyyyyyz, g_0_z_xxxx_xyyyyyz, g_0_z_xxxx_xyyyyyzz, g_0_z_xxxx_xyyyyzz, g_0_z_xxxx_xyyyyzzz, g_0_z_xxxx_xyyyzzz, g_0_z_xxxx_xyyyzzzz, g_0_z_xxxx_xyyzzzz, g_0_z_xxxx_xyyzzzzz, g_0_z_xxxx_xyzzzzz, g_0_z_xxxx_xyzzzzzz, g_0_z_xxxx_xzzzzzz, g_0_z_xxxx_xzzzzzzz, g_0_z_xxxx_yyyyyyy, g_0_z_xxxx_yyyyyyz, g_0_z_xxxx_yyyyyzz, g_0_z_xxxx_yyyyzzz, g_0_z_xxxx_yyyzzzz, g_0_z_xxxx_yyzzzzz, g_0_z_xxxx_yzzzzzz, g_0_z_xxxx_zzzzzzz, g_0_z_xxxxx_xxxxxxx, g_0_z_xxxxx_xxxxxxy, g_0_z_xxxxx_xxxxxxz, g_0_z_xxxxx_xxxxxyy, g_0_z_xxxxx_xxxxxyz, g_0_z_xxxxx_xxxxxzz, g_0_z_xxxxx_xxxxyyy, g_0_z_xxxxx_xxxxyyz, g_0_z_xxxxx_xxxxyzz, g_0_z_xxxxx_xxxxzzz, g_0_z_xxxxx_xxxyyyy, g_0_z_xxxxx_xxxyyyz, g_0_z_xxxxx_xxxyyzz, g_0_z_xxxxx_xxxyzzz, g_0_z_xxxxx_xxxzzzz, g_0_z_xxxxx_xxyyyyy, g_0_z_xxxxx_xxyyyyz, g_0_z_xxxxx_xxyyyzz, g_0_z_xxxxx_xxyyzzz, g_0_z_xxxxx_xxyzzzz, g_0_z_xxxxx_xxzzzzz, g_0_z_xxxxx_xyyyyyy, g_0_z_xxxxx_xyyyyyz, g_0_z_xxxxx_xyyyyzz, g_0_z_xxxxx_xyyyzzz, g_0_z_xxxxx_xyyzzzz, g_0_z_xxxxx_xyzzzzz, g_0_z_xxxxx_xzzzzzz, g_0_z_xxxxx_yyyyyyy, g_0_z_xxxxx_yyyyyyz, g_0_z_xxxxx_yyyyyzz, g_0_z_xxxxx_yyyyzzz, g_0_z_xxxxx_yyyzzzz, g_0_z_xxxxx_yyzzzzz, g_0_z_xxxxx_yzzzzzz, g_0_z_xxxxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxx_xxxxxxx[k] = -g_0_z_xxxx_xxxxxxx[k] * ab_x + g_0_z_xxxx_xxxxxxxx[k];

                g_0_z_xxxxx_xxxxxxy[k] = -g_0_z_xxxx_xxxxxxy[k] * ab_x + g_0_z_xxxx_xxxxxxxy[k];

                g_0_z_xxxxx_xxxxxxz[k] = -g_0_z_xxxx_xxxxxxz[k] * ab_x + g_0_z_xxxx_xxxxxxxz[k];

                g_0_z_xxxxx_xxxxxyy[k] = -g_0_z_xxxx_xxxxxyy[k] * ab_x + g_0_z_xxxx_xxxxxxyy[k];

                g_0_z_xxxxx_xxxxxyz[k] = -g_0_z_xxxx_xxxxxyz[k] * ab_x + g_0_z_xxxx_xxxxxxyz[k];

                g_0_z_xxxxx_xxxxxzz[k] = -g_0_z_xxxx_xxxxxzz[k] * ab_x + g_0_z_xxxx_xxxxxxzz[k];

                g_0_z_xxxxx_xxxxyyy[k] = -g_0_z_xxxx_xxxxyyy[k] * ab_x + g_0_z_xxxx_xxxxxyyy[k];

                g_0_z_xxxxx_xxxxyyz[k] = -g_0_z_xxxx_xxxxyyz[k] * ab_x + g_0_z_xxxx_xxxxxyyz[k];

                g_0_z_xxxxx_xxxxyzz[k] = -g_0_z_xxxx_xxxxyzz[k] * ab_x + g_0_z_xxxx_xxxxxyzz[k];

                g_0_z_xxxxx_xxxxzzz[k] = -g_0_z_xxxx_xxxxzzz[k] * ab_x + g_0_z_xxxx_xxxxxzzz[k];

                g_0_z_xxxxx_xxxyyyy[k] = -g_0_z_xxxx_xxxyyyy[k] * ab_x + g_0_z_xxxx_xxxxyyyy[k];

                g_0_z_xxxxx_xxxyyyz[k] = -g_0_z_xxxx_xxxyyyz[k] * ab_x + g_0_z_xxxx_xxxxyyyz[k];

                g_0_z_xxxxx_xxxyyzz[k] = -g_0_z_xxxx_xxxyyzz[k] * ab_x + g_0_z_xxxx_xxxxyyzz[k];

                g_0_z_xxxxx_xxxyzzz[k] = -g_0_z_xxxx_xxxyzzz[k] * ab_x + g_0_z_xxxx_xxxxyzzz[k];

                g_0_z_xxxxx_xxxzzzz[k] = -g_0_z_xxxx_xxxzzzz[k] * ab_x + g_0_z_xxxx_xxxxzzzz[k];

                g_0_z_xxxxx_xxyyyyy[k] = -g_0_z_xxxx_xxyyyyy[k] * ab_x + g_0_z_xxxx_xxxyyyyy[k];

                g_0_z_xxxxx_xxyyyyz[k] = -g_0_z_xxxx_xxyyyyz[k] * ab_x + g_0_z_xxxx_xxxyyyyz[k];

                g_0_z_xxxxx_xxyyyzz[k] = -g_0_z_xxxx_xxyyyzz[k] * ab_x + g_0_z_xxxx_xxxyyyzz[k];

                g_0_z_xxxxx_xxyyzzz[k] = -g_0_z_xxxx_xxyyzzz[k] * ab_x + g_0_z_xxxx_xxxyyzzz[k];

                g_0_z_xxxxx_xxyzzzz[k] = -g_0_z_xxxx_xxyzzzz[k] * ab_x + g_0_z_xxxx_xxxyzzzz[k];

                g_0_z_xxxxx_xxzzzzz[k] = -g_0_z_xxxx_xxzzzzz[k] * ab_x + g_0_z_xxxx_xxxzzzzz[k];

                g_0_z_xxxxx_xyyyyyy[k] = -g_0_z_xxxx_xyyyyyy[k] * ab_x + g_0_z_xxxx_xxyyyyyy[k];

                g_0_z_xxxxx_xyyyyyz[k] = -g_0_z_xxxx_xyyyyyz[k] * ab_x + g_0_z_xxxx_xxyyyyyz[k];

                g_0_z_xxxxx_xyyyyzz[k] = -g_0_z_xxxx_xyyyyzz[k] * ab_x + g_0_z_xxxx_xxyyyyzz[k];

                g_0_z_xxxxx_xyyyzzz[k] = -g_0_z_xxxx_xyyyzzz[k] * ab_x + g_0_z_xxxx_xxyyyzzz[k];

                g_0_z_xxxxx_xyyzzzz[k] = -g_0_z_xxxx_xyyzzzz[k] * ab_x + g_0_z_xxxx_xxyyzzzz[k];

                g_0_z_xxxxx_xyzzzzz[k] = -g_0_z_xxxx_xyzzzzz[k] * ab_x + g_0_z_xxxx_xxyzzzzz[k];

                g_0_z_xxxxx_xzzzzzz[k] = -g_0_z_xxxx_xzzzzzz[k] * ab_x + g_0_z_xxxx_xxzzzzzz[k];

                g_0_z_xxxxx_yyyyyyy[k] = -g_0_z_xxxx_yyyyyyy[k] * ab_x + g_0_z_xxxx_xyyyyyyy[k];

                g_0_z_xxxxx_yyyyyyz[k] = -g_0_z_xxxx_yyyyyyz[k] * ab_x + g_0_z_xxxx_xyyyyyyz[k];

                g_0_z_xxxxx_yyyyyzz[k] = -g_0_z_xxxx_yyyyyzz[k] * ab_x + g_0_z_xxxx_xyyyyyzz[k];

                g_0_z_xxxxx_yyyyzzz[k] = -g_0_z_xxxx_yyyyzzz[k] * ab_x + g_0_z_xxxx_xyyyyzzz[k];

                g_0_z_xxxxx_yyyzzzz[k] = -g_0_z_xxxx_yyyzzzz[k] * ab_x + g_0_z_xxxx_xyyyzzzz[k];

                g_0_z_xxxxx_yyzzzzz[k] = -g_0_z_xxxx_yyzzzzz[k] * ab_x + g_0_z_xxxx_xyyzzzzz[k];

                g_0_z_xxxxx_yzzzzzz[k] = -g_0_z_xxxx_yzzzzzz[k] * ab_x + g_0_z_xxxx_xyzzzzzz[k];

                g_0_z_xxxxx_zzzzzzz[k] = -g_0_z_xxxx_zzzzzzz[k] * ab_x + g_0_z_xxxx_xzzzzzzz[k];
            }

            /// Set up 1548-1584 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1548 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1549 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1550 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1551 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1552 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1553 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1554 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1555 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1556 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1557 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1558 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1559 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1560 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1561 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1562 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1563 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1564 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1565 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1566 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1567 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1568 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1569 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1570 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1571 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1572 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1573 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1574 * ccomps * dcomps);

            auto g_0_z_xxxxy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1575 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1576 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1577 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1578 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1579 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1580 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1581 * ccomps * dcomps);

            auto g_0_z_xxxxy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1582 * ccomps * dcomps);

            auto g_0_z_xxxxy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1583 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxy_xxxxxxx, g_0_z_xxxxy_xxxxxxy, g_0_z_xxxxy_xxxxxxz, g_0_z_xxxxy_xxxxxyy, g_0_z_xxxxy_xxxxxyz, g_0_z_xxxxy_xxxxxzz, g_0_z_xxxxy_xxxxyyy, g_0_z_xxxxy_xxxxyyz, g_0_z_xxxxy_xxxxyzz, g_0_z_xxxxy_xxxxzzz, g_0_z_xxxxy_xxxyyyy, g_0_z_xxxxy_xxxyyyz, g_0_z_xxxxy_xxxyyzz, g_0_z_xxxxy_xxxyzzz, g_0_z_xxxxy_xxxzzzz, g_0_z_xxxxy_xxyyyyy, g_0_z_xxxxy_xxyyyyz, g_0_z_xxxxy_xxyyyzz, g_0_z_xxxxy_xxyyzzz, g_0_z_xxxxy_xxyzzzz, g_0_z_xxxxy_xxzzzzz, g_0_z_xxxxy_xyyyyyy, g_0_z_xxxxy_xyyyyyz, g_0_z_xxxxy_xyyyyzz, g_0_z_xxxxy_xyyyzzz, g_0_z_xxxxy_xyyzzzz, g_0_z_xxxxy_xyzzzzz, g_0_z_xxxxy_xzzzzzz, g_0_z_xxxxy_yyyyyyy, g_0_z_xxxxy_yyyyyyz, g_0_z_xxxxy_yyyyyzz, g_0_z_xxxxy_yyyyzzz, g_0_z_xxxxy_yyyzzzz, g_0_z_xxxxy_yyzzzzz, g_0_z_xxxxy_yzzzzzz, g_0_z_xxxxy_zzzzzzz, g_0_z_xxxy_xxxxxxx, g_0_z_xxxy_xxxxxxxx, g_0_z_xxxy_xxxxxxxy, g_0_z_xxxy_xxxxxxxz, g_0_z_xxxy_xxxxxxy, g_0_z_xxxy_xxxxxxyy, g_0_z_xxxy_xxxxxxyz, g_0_z_xxxy_xxxxxxz, g_0_z_xxxy_xxxxxxzz, g_0_z_xxxy_xxxxxyy, g_0_z_xxxy_xxxxxyyy, g_0_z_xxxy_xxxxxyyz, g_0_z_xxxy_xxxxxyz, g_0_z_xxxy_xxxxxyzz, g_0_z_xxxy_xxxxxzz, g_0_z_xxxy_xxxxxzzz, g_0_z_xxxy_xxxxyyy, g_0_z_xxxy_xxxxyyyy, g_0_z_xxxy_xxxxyyyz, g_0_z_xxxy_xxxxyyz, g_0_z_xxxy_xxxxyyzz, g_0_z_xxxy_xxxxyzz, g_0_z_xxxy_xxxxyzzz, g_0_z_xxxy_xxxxzzz, g_0_z_xxxy_xxxxzzzz, g_0_z_xxxy_xxxyyyy, g_0_z_xxxy_xxxyyyyy, g_0_z_xxxy_xxxyyyyz, g_0_z_xxxy_xxxyyyz, g_0_z_xxxy_xxxyyyzz, g_0_z_xxxy_xxxyyzz, g_0_z_xxxy_xxxyyzzz, g_0_z_xxxy_xxxyzzz, g_0_z_xxxy_xxxyzzzz, g_0_z_xxxy_xxxzzzz, g_0_z_xxxy_xxxzzzzz, g_0_z_xxxy_xxyyyyy, g_0_z_xxxy_xxyyyyyy, g_0_z_xxxy_xxyyyyyz, g_0_z_xxxy_xxyyyyz, g_0_z_xxxy_xxyyyyzz, g_0_z_xxxy_xxyyyzz, g_0_z_xxxy_xxyyyzzz, g_0_z_xxxy_xxyyzzz, g_0_z_xxxy_xxyyzzzz, g_0_z_xxxy_xxyzzzz, g_0_z_xxxy_xxyzzzzz, g_0_z_xxxy_xxzzzzz, g_0_z_xxxy_xxzzzzzz, g_0_z_xxxy_xyyyyyy, g_0_z_xxxy_xyyyyyyy, g_0_z_xxxy_xyyyyyyz, g_0_z_xxxy_xyyyyyz, g_0_z_xxxy_xyyyyyzz, g_0_z_xxxy_xyyyyzz, g_0_z_xxxy_xyyyyzzz, g_0_z_xxxy_xyyyzzz, g_0_z_xxxy_xyyyzzzz, g_0_z_xxxy_xyyzzzz, g_0_z_xxxy_xyyzzzzz, g_0_z_xxxy_xyzzzzz, g_0_z_xxxy_xyzzzzzz, g_0_z_xxxy_xzzzzzz, g_0_z_xxxy_xzzzzzzz, g_0_z_xxxy_yyyyyyy, g_0_z_xxxy_yyyyyyz, g_0_z_xxxy_yyyyyzz, g_0_z_xxxy_yyyyzzz, g_0_z_xxxy_yyyzzzz, g_0_z_xxxy_yyzzzzz, g_0_z_xxxy_yzzzzzz, g_0_z_xxxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxy_xxxxxxx[k] = -g_0_z_xxxy_xxxxxxx[k] * ab_x + g_0_z_xxxy_xxxxxxxx[k];

                g_0_z_xxxxy_xxxxxxy[k] = -g_0_z_xxxy_xxxxxxy[k] * ab_x + g_0_z_xxxy_xxxxxxxy[k];

                g_0_z_xxxxy_xxxxxxz[k] = -g_0_z_xxxy_xxxxxxz[k] * ab_x + g_0_z_xxxy_xxxxxxxz[k];

                g_0_z_xxxxy_xxxxxyy[k] = -g_0_z_xxxy_xxxxxyy[k] * ab_x + g_0_z_xxxy_xxxxxxyy[k];

                g_0_z_xxxxy_xxxxxyz[k] = -g_0_z_xxxy_xxxxxyz[k] * ab_x + g_0_z_xxxy_xxxxxxyz[k];

                g_0_z_xxxxy_xxxxxzz[k] = -g_0_z_xxxy_xxxxxzz[k] * ab_x + g_0_z_xxxy_xxxxxxzz[k];

                g_0_z_xxxxy_xxxxyyy[k] = -g_0_z_xxxy_xxxxyyy[k] * ab_x + g_0_z_xxxy_xxxxxyyy[k];

                g_0_z_xxxxy_xxxxyyz[k] = -g_0_z_xxxy_xxxxyyz[k] * ab_x + g_0_z_xxxy_xxxxxyyz[k];

                g_0_z_xxxxy_xxxxyzz[k] = -g_0_z_xxxy_xxxxyzz[k] * ab_x + g_0_z_xxxy_xxxxxyzz[k];

                g_0_z_xxxxy_xxxxzzz[k] = -g_0_z_xxxy_xxxxzzz[k] * ab_x + g_0_z_xxxy_xxxxxzzz[k];

                g_0_z_xxxxy_xxxyyyy[k] = -g_0_z_xxxy_xxxyyyy[k] * ab_x + g_0_z_xxxy_xxxxyyyy[k];

                g_0_z_xxxxy_xxxyyyz[k] = -g_0_z_xxxy_xxxyyyz[k] * ab_x + g_0_z_xxxy_xxxxyyyz[k];

                g_0_z_xxxxy_xxxyyzz[k] = -g_0_z_xxxy_xxxyyzz[k] * ab_x + g_0_z_xxxy_xxxxyyzz[k];

                g_0_z_xxxxy_xxxyzzz[k] = -g_0_z_xxxy_xxxyzzz[k] * ab_x + g_0_z_xxxy_xxxxyzzz[k];

                g_0_z_xxxxy_xxxzzzz[k] = -g_0_z_xxxy_xxxzzzz[k] * ab_x + g_0_z_xxxy_xxxxzzzz[k];

                g_0_z_xxxxy_xxyyyyy[k] = -g_0_z_xxxy_xxyyyyy[k] * ab_x + g_0_z_xxxy_xxxyyyyy[k];

                g_0_z_xxxxy_xxyyyyz[k] = -g_0_z_xxxy_xxyyyyz[k] * ab_x + g_0_z_xxxy_xxxyyyyz[k];

                g_0_z_xxxxy_xxyyyzz[k] = -g_0_z_xxxy_xxyyyzz[k] * ab_x + g_0_z_xxxy_xxxyyyzz[k];

                g_0_z_xxxxy_xxyyzzz[k] = -g_0_z_xxxy_xxyyzzz[k] * ab_x + g_0_z_xxxy_xxxyyzzz[k];

                g_0_z_xxxxy_xxyzzzz[k] = -g_0_z_xxxy_xxyzzzz[k] * ab_x + g_0_z_xxxy_xxxyzzzz[k];

                g_0_z_xxxxy_xxzzzzz[k] = -g_0_z_xxxy_xxzzzzz[k] * ab_x + g_0_z_xxxy_xxxzzzzz[k];

                g_0_z_xxxxy_xyyyyyy[k] = -g_0_z_xxxy_xyyyyyy[k] * ab_x + g_0_z_xxxy_xxyyyyyy[k];

                g_0_z_xxxxy_xyyyyyz[k] = -g_0_z_xxxy_xyyyyyz[k] * ab_x + g_0_z_xxxy_xxyyyyyz[k];

                g_0_z_xxxxy_xyyyyzz[k] = -g_0_z_xxxy_xyyyyzz[k] * ab_x + g_0_z_xxxy_xxyyyyzz[k];

                g_0_z_xxxxy_xyyyzzz[k] = -g_0_z_xxxy_xyyyzzz[k] * ab_x + g_0_z_xxxy_xxyyyzzz[k];

                g_0_z_xxxxy_xyyzzzz[k] = -g_0_z_xxxy_xyyzzzz[k] * ab_x + g_0_z_xxxy_xxyyzzzz[k];

                g_0_z_xxxxy_xyzzzzz[k] = -g_0_z_xxxy_xyzzzzz[k] * ab_x + g_0_z_xxxy_xxyzzzzz[k];

                g_0_z_xxxxy_xzzzzzz[k] = -g_0_z_xxxy_xzzzzzz[k] * ab_x + g_0_z_xxxy_xxzzzzzz[k];

                g_0_z_xxxxy_yyyyyyy[k] = -g_0_z_xxxy_yyyyyyy[k] * ab_x + g_0_z_xxxy_xyyyyyyy[k];

                g_0_z_xxxxy_yyyyyyz[k] = -g_0_z_xxxy_yyyyyyz[k] * ab_x + g_0_z_xxxy_xyyyyyyz[k];

                g_0_z_xxxxy_yyyyyzz[k] = -g_0_z_xxxy_yyyyyzz[k] * ab_x + g_0_z_xxxy_xyyyyyzz[k];

                g_0_z_xxxxy_yyyyzzz[k] = -g_0_z_xxxy_yyyyzzz[k] * ab_x + g_0_z_xxxy_xyyyyzzz[k];

                g_0_z_xxxxy_yyyzzzz[k] = -g_0_z_xxxy_yyyzzzz[k] * ab_x + g_0_z_xxxy_xyyyzzzz[k];

                g_0_z_xxxxy_yyzzzzz[k] = -g_0_z_xxxy_yyzzzzz[k] * ab_x + g_0_z_xxxy_xyyzzzzz[k];

                g_0_z_xxxxy_yzzzzzz[k] = -g_0_z_xxxy_yzzzzzz[k] * ab_x + g_0_z_xxxy_xyzzzzzz[k];

                g_0_z_xxxxy_zzzzzzz[k] = -g_0_z_xxxy_zzzzzzz[k] * ab_x + g_0_z_xxxy_xzzzzzzz[k];
            }

            /// Set up 1584-1620 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1584 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1585 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1586 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1587 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1588 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1589 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1590 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1591 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1592 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1593 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1594 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1595 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1596 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1597 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1598 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1599 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1600 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1601 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1602 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1603 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1604 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1605 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1606 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1607 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1608 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1609 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1610 * ccomps * dcomps);

            auto g_0_z_xxxxz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1611 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1612 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1613 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1614 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1615 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1616 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1617 * ccomps * dcomps);

            auto g_0_z_xxxxz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1618 * ccomps * dcomps);

            auto g_0_z_xxxxz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1619 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxz_xxxxxxx, g_0_z_xxxxz_xxxxxxy, g_0_z_xxxxz_xxxxxxz, g_0_z_xxxxz_xxxxxyy, g_0_z_xxxxz_xxxxxyz, g_0_z_xxxxz_xxxxxzz, g_0_z_xxxxz_xxxxyyy, g_0_z_xxxxz_xxxxyyz, g_0_z_xxxxz_xxxxyzz, g_0_z_xxxxz_xxxxzzz, g_0_z_xxxxz_xxxyyyy, g_0_z_xxxxz_xxxyyyz, g_0_z_xxxxz_xxxyyzz, g_0_z_xxxxz_xxxyzzz, g_0_z_xxxxz_xxxzzzz, g_0_z_xxxxz_xxyyyyy, g_0_z_xxxxz_xxyyyyz, g_0_z_xxxxz_xxyyyzz, g_0_z_xxxxz_xxyyzzz, g_0_z_xxxxz_xxyzzzz, g_0_z_xxxxz_xxzzzzz, g_0_z_xxxxz_xyyyyyy, g_0_z_xxxxz_xyyyyyz, g_0_z_xxxxz_xyyyyzz, g_0_z_xxxxz_xyyyzzz, g_0_z_xxxxz_xyyzzzz, g_0_z_xxxxz_xyzzzzz, g_0_z_xxxxz_xzzzzzz, g_0_z_xxxxz_yyyyyyy, g_0_z_xxxxz_yyyyyyz, g_0_z_xxxxz_yyyyyzz, g_0_z_xxxxz_yyyyzzz, g_0_z_xxxxz_yyyzzzz, g_0_z_xxxxz_yyzzzzz, g_0_z_xxxxz_yzzzzzz, g_0_z_xxxxz_zzzzzzz, g_0_z_xxxz_xxxxxxx, g_0_z_xxxz_xxxxxxxx, g_0_z_xxxz_xxxxxxxy, g_0_z_xxxz_xxxxxxxz, g_0_z_xxxz_xxxxxxy, g_0_z_xxxz_xxxxxxyy, g_0_z_xxxz_xxxxxxyz, g_0_z_xxxz_xxxxxxz, g_0_z_xxxz_xxxxxxzz, g_0_z_xxxz_xxxxxyy, g_0_z_xxxz_xxxxxyyy, g_0_z_xxxz_xxxxxyyz, g_0_z_xxxz_xxxxxyz, g_0_z_xxxz_xxxxxyzz, g_0_z_xxxz_xxxxxzz, g_0_z_xxxz_xxxxxzzz, g_0_z_xxxz_xxxxyyy, g_0_z_xxxz_xxxxyyyy, g_0_z_xxxz_xxxxyyyz, g_0_z_xxxz_xxxxyyz, g_0_z_xxxz_xxxxyyzz, g_0_z_xxxz_xxxxyzz, g_0_z_xxxz_xxxxyzzz, g_0_z_xxxz_xxxxzzz, g_0_z_xxxz_xxxxzzzz, g_0_z_xxxz_xxxyyyy, g_0_z_xxxz_xxxyyyyy, g_0_z_xxxz_xxxyyyyz, g_0_z_xxxz_xxxyyyz, g_0_z_xxxz_xxxyyyzz, g_0_z_xxxz_xxxyyzz, g_0_z_xxxz_xxxyyzzz, g_0_z_xxxz_xxxyzzz, g_0_z_xxxz_xxxyzzzz, g_0_z_xxxz_xxxzzzz, g_0_z_xxxz_xxxzzzzz, g_0_z_xxxz_xxyyyyy, g_0_z_xxxz_xxyyyyyy, g_0_z_xxxz_xxyyyyyz, g_0_z_xxxz_xxyyyyz, g_0_z_xxxz_xxyyyyzz, g_0_z_xxxz_xxyyyzz, g_0_z_xxxz_xxyyyzzz, g_0_z_xxxz_xxyyzzz, g_0_z_xxxz_xxyyzzzz, g_0_z_xxxz_xxyzzzz, g_0_z_xxxz_xxyzzzzz, g_0_z_xxxz_xxzzzzz, g_0_z_xxxz_xxzzzzzz, g_0_z_xxxz_xyyyyyy, g_0_z_xxxz_xyyyyyyy, g_0_z_xxxz_xyyyyyyz, g_0_z_xxxz_xyyyyyz, g_0_z_xxxz_xyyyyyzz, g_0_z_xxxz_xyyyyzz, g_0_z_xxxz_xyyyyzzz, g_0_z_xxxz_xyyyzzz, g_0_z_xxxz_xyyyzzzz, g_0_z_xxxz_xyyzzzz, g_0_z_xxxz_xyyzzzzz, g_0_z_xxxz_xyzzzzz, g_0_z_xxxz_xyzzzzzz, g_0_z_xxxz_xzzzzzz, g_0_z_xxxz_xzzzzzzz, g_0_z_xxxz_yyyyyyy, g_0_z_xxxz_yyyyyyz, g_0_z_xxxz_yyyyyzz, g_0_z_xxxz_yyyyzzz, g_0_z_xxxz_yyyzzzz, g_0_z_xxxz_yyzzzzz, g_0_z_xxxz_yzzzzzz, g_0_z_xxxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxz_xxxxxxx[k] = -g_0_z_xxxz_xxxxxxx[k] * ab_x + g_0_z_xxxz_xxxxxxxx[k];

                g_0_z_xxxxz_xxxxxxy[k] = -g_0_z_xxxz_xxxxxxy[k] * ab_x + g_0_z_xxxz_xxxxxxxy[k];

                g_0_z_xxxxz_xxxxxxz[k] = -g_0_z_xxxz_xxxxxxz[k] * ab_x + g_0_z_xxxz_xxxxxxxz[k];

                g_0_z_xxxxz_xxxxxyy[k] = -g_0_z_xxxz_xxxxxyy[k] * ab_x + g_0_z_xxxz_xxxxxxyy[k];

                g_0_z_xxxxz_xxxxxyz[k] = -g_0_z_xxxz_xxxxxyz[k] * ab_x + g_0_z_xxxz_xxxxxxyz[k];

                g_0_z_xxxxz_xxxxxzz[k] = -g_0_z_xxxz_xxxxxzz[k] * ab_x + g_0_z_xxxz_xxxxxxzz[k];

                g_0_z_xxxxz_xxxxyyy[k] = -g_0_z_xxxz_xxxxyyy[k] * ab_x + g_0_z_xxxz_xxxxxyyy[k];

                g_0_z_xxxxz_xxxxyyz[k] = -g_0_z_xxxz_xxxxyyz[k] * ab_x + g_0_z_xxxz_xxxxxyyz[k];

                g_0_z_xxxxz_xxxxyzz[k] = -g_0_z_xxxz_xxxxyzz[k] * ab_x + g_0_z_xxxz_xxxxxyzz[k];

                g_0_z_xxxxz_xxxxzzz[k] = -g_0_z_xxxz_xxxxzzz[k] * ab_x + g_0_z_xxxz_xxxxxzzz[k];

                g_0_z_xxxxz_xxxyyyy[k] = -g_0_z_xxxz_xxxyyyy[k] * ab_x + g_0_z_xxxz_xxxxyyyy[k];

                g_0_z_xxxxz_xxxyyyz[k] = -g_0_z_xxxz_xxxyyyz[k] * ab_x + g_0_z_xxxz_xxxxyyyz[k];

                g_0_z_xxxxz_xxxyyzz[k] = -g_0_z_xxxz_xxxyyzz[k] * ab_x + g_0_z_xxxz_xxxxyyzz[k];

                g_0_z_xxxxz_xxxyzzz[k] = -g_0_z_xxxz_xxxyzzz[k] * ab_x + g_0_z_xxxz_xxxxyzzz[k];

                g_0_z_xxxxz_xxxzzzz[k] = -g_0_z_xxxz_xxxzzzz[k] * ab_x + g_0_z_xxxz_xxxxzzzz[k];

                g_0_z_xxxxz_xxyyyyy[k] = -g_0_z_xxxz_xxyyyyy[k] * ab_x + g_0_z_xxxz_xxxyyyyy[k];

                g_0_z_xxxxz_xxyyyyz[k] = -g_0_z_xxxz_xxyyyyz[k] * ab_x + g_0_z_xxxz_xxxyyyyz[k];

                g_0_z_xxxxz_xxyyyzz[k] = -g_0_z_xxxz_xxyyyzz[k] * ab_x + g_0_z_xxxz_xxxyyyzz[k];

                g_0_z_xxxxz_xxyyzzz[k] = -g_0_z_xxxz_xxyyzzz[k] * ab_x + g_0_z_xxxz_xxxyyzzz[k];

                g_0_z_xxxxz_xxyzzzz[k] = -g_0_z_xxxz_xxyzzzz[k] * ab_x + g_0_z_xxxz_xxxyzzzz[k];

                g_0_z_xxxxz_xxzzzzz[k] = -g_0_z_xxxz_xxzzzzz[k] * ab_x + g_0_z_xxxz_xxxzzzzz[k];

                g_0_z_xxxxz_xyyyyyy[k] = -g_0_z_xxxz_xyyyyyy[k] * ab_x + g_0_z_xxxz_xxyyyyyy[k];

                g_0_z_xxxxz_xyyyyyz[k] = -g_0_z_xxxz_xyyyyyz[k] * ab_x + g_0_z_xxxz_xxyyyyyz[k];

                g_0_z_xxxxz_xyyyyzz[k] = -g_0_z_xxxz_xyyyyzz[k] * ab_x + g_0_z_xxxz_xxyyyyzz[k];

                g_0_z_xxxxz_xyyyzzz[k] = -g_0_z_xxxz_xyyyzzz[k] * ab_x + g_0_z_xxxz_xxyyyzzz[k];

                g_0_z_xxxxz_xyyzzzz[k] = -g_0_z_xxxz_xyyzzzz[k] * ab_x + g_0_z_xxxz_xxyyzzzz[k];

                g_0_z_xxxxz_xyzzzzz[k] = -g_0_z_xxxz_xyzzzzz[k] * ab_x + g_0_z_xxxz_xxyzzzzz[k];

                g_0_z_xxxxz_xzzzzzz[k] = -g_0_z_xxxz_xzzzzzz[k] * ab_x + g_0_z_xxxz_xxzzzzzz[k];

                g_0_z_xxxxz_yyyyyyy[k] = -g_0_z_xxxz_yyyyyyy[k] * ab_x + g_0_z_xxxz_xyyyyyyy[k];

                g_0_z_xxxxz_yyyyyyz[k] = -g_0_z_xxxz_yyyyyyz[k] * ab_x + g_0_z_xxxz_xyyyyyyz[k];

                g_0_z_xxxxz_yyyyyzz[k] = -g_0_z_xxxz_yyyyyzz[k] * ab_x + g_0_z_xxxz_xyyyyyzz[k];

                g_0_z_xxxxz_yyyyzzz[k] = -g_0_z_xxxz_yyyyzzz[k] * ab_x + g_0_z_xxxz_xyyyyzzz[k];

                g_0_z_xxxxz_yyyzzzz[k] = -g_0_z_xxxz_yyyzzzz[k] * ab_x + g_0_z_xxxz_xyyyzzzz[k];

                g_0_z_xxxxz_yyzzzzz[k] = -g_0_z_xxxz_yyzzzzz[k] * ab_x + g_0_z_xxxz_xyyzzzzz[k];

                g_0_z_xxxxz_yzzzzzz[k] = -g_0_z_xxxz_yzzzzzz[k] * ab_x + g_0_z_xxxz_xyzzzzzz[k];

                g_0_z_xxxxz_zzzzzzz[k] = -g_0_z_xxxz_zzzzzzz[k] * ab_x + g_0_z_xxxz_xzzzzzzz[k];
            }

            /// Set up 1620-1656 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1620 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1621 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1622 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1623 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1624 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1625 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1626 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1627 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1628 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1629 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1630 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1631 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1632 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1633 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1634 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1635 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1636 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1637 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1638 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1639 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1640 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1641 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1642 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1643 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1644 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1645 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1646 * ccomps * dcomps);

            auto g_0_z_xxxyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1647 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1648 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1649 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1650 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1651 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1652 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1653 * ccomps * dcomps);

            auto g_0_z_xxxyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1654 * ccomps * dcomps);

            auto g_0_z_xxxyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1655 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyy_xxxxxxx, g_0_z_xxxyy_xxxxxxy, g_0_z_xxxyy_xxxxxxz, g_0_z_xxxyy_xxxxxyy, g_0_z_xxxyy_xxxxxyz, g_0_z_xxxyy_xxxxxzz, g_0_z_xxxyy_xxxxyyy, g_0_z_xxxyy_xxxxyyz, g_0_z_xxxyy_xxxxyzz, g_0_z_xxxyy_xxxxzzz, g_0_z_xxxyy_xxxyyyy, g_0_z_xxxyy_xxxyyyz, g_0_z_xxxyy_xxxyyzz, g_0_z_xxxyy_xxxyzzz, g_0_z_xxxyy_xxxzzzz, g_0_z_xxxyy_xxyyyyy, g_0_z_xxxyy_xxyyyyz, g_0_z_xxxyy_xxyyyzz, g_0_z_xxxyy_xxyyzzz, g_0_z_xxxyy_xxyzzzz, g_0_z_xxxyy_xxzzzzz, g_0_z_xxxyy_xyyyyyy, g_0_z_xxxyy_xyyyyyz, g_0_z_xxxyy_xyyyyzz, g_0_z_xxxyy_xyyyzzz, g_0_z_xxxyy_xyyzzzz, g_0_z_xxxyy_xyzzzzz, g_0_z_xxxyy_xzzzzzz, g_0_z_xxxyy_yyyyyyy, g_0_z_xxxyy_yyyyyyz, g_0_z_xxxyy_yyyyyzz, g_0_z_xxxyy_yyyyzzz, g_0_z_xxxyy_yyyzzzz, g_0_z_xxxyy_yyzzzzz, g_0_z_xxxyy_yzzzzzz, g_0_z_xxxyy_zzzzzzz, g_0_z_xxyy_xxxxxxx, g_0_z_xxyy_xxxxxxxx, g_0_z_xxyy_xxxxxxxy, g_0_z_xxyy_xxxxxxxz, g_0_z_xxyy_xxxxxxy, g_0_z_xxyy_xxxxxxyy, g_0_z_xxyy_xxxxxxyz, g_0_z_xxyy_xxxxxxz, g_0_z_xxyy_xxxxxxzz, g_0_z_xxyy_xxxxxyy, g_0_z_xxyy_xxxxxyyy, g_0_z_xxyy_xxxxxyyz, g_0_z_xxyy_xxxxxyz, g_0_z_xxyy_xxxxxyzz, g_0_z_xxyy_xxxxxzz, g_0_z_xxyy_xxxxxzzz, g_0_z_xxyy_xxxxyyy, g_0_z_xxyy_xxxxyyyy, g_0_z_xxyy_xxxxyyyz, g_0_z_xxyy_xxxxyyz, g_0_z_xxyy_xxxxyyzz, g_0_z_xxyy_xxxxyzz, g_0_z_xxyy_xxxxyzzz, g_0_z_xxyy_xxxxzzz, g_0_z_xxyy_xxxxzzzz, g_0_z_xxyy_xxxyyyy, g_0_z_xxyy_xxxyyyyy, g_0_z_xxyy_xxxyyyyz, g_0_z_xxyy_xxxyyyz, g_0_z_xxyy_xxxyyyzz, g_0_z_xxyy_xxxyyzz, g_0_z_xxyy_xxxyyzzz, g_0_z_xxyy_xxxyzzz, g_0_z_xxyy_xxxyzzzz, g_0_z_xxyy_xxxzzzz, g_0_z_xxyy_xxxzzzzz, g_0_z_xxyy_xxyyyyy, g_0_z_xxyy_xxyyyyyy, g_0_z_xxyy_xxyyyyyz, g_0_z_xxyy_xxyyyyz, g_0_z_xxyy_xxyyyyzz, g_0_z_xxyy_xxyyyzz, g_0_z_xxyy_xxyyyzzz, g_0_z_xxyy_xxyyzzz, g_0_z_xxyy_xxyyzzzz, g_0_z_xxyy_xxyzzzz, g_0_z_xxyy_xxyzzzzz, g_0_z_xxyy_xxzzzzz, g_0_z_xxyy_xxzzzzzz, g_0_z_xxyy_xyyyyyy, g_0_z_xxyy_xyyyyyyy, g_0_z_xxyy_xyyyyyyz, g_0_z_xxyy_xyyyyyz, g_0_z_xxyy_xyyyyyzz, g_0_z_xxyy_xyyyyzz, g_0_z_xxyy_xyyyyzzz, g_0_z_xxyy_xyyyzzz, g_0_z_xxyy_xyyyzzzz, g_0_z_xxyy_xyyzzzz, g_0_z_xxyy_xyyzzzzz, g_0_z_xxyy_xyzzzzz, g_0_z_xxyy_xyzzzzzz, g_0_z_xxyy_xzzzzzz, g_0_z_xxyy_xzzzzzzz, g_0_z_xxyy_yyyyyyy, g_0_z_xxyy_yyyyyyz, g_0_z_xxyy_yyyyyzz, g_0_z_xxyy_yyyyzzz, g_0_z_xxyy_yyyzzzz, g_0_z_xxyy_yyzzzzz, g_0_z_xxyy_yzzzzzz, g_0_z_xxyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyy_xxxxxxx[k] = -g_0_z_xxyy_xxxxxxx[k] * ab_x + g_0_z_xxyy_xxxxxxxx[k];

                g_0_z_xxxyy_xxxxxxy[k] = -g_0_z_xxyy_xxxxxxy[k] * ab_x + g_0_z_xxyy_xxxxxxxy[k];

                g_0_z_xxxyy_xxxxxxz[k] = -g_0_z_xxyy_xxxxxxz[k] * ab_x + g_0_z_xxyy_xxxxxxxz[k];

                g_0_z_xxxyy_xxxxxyy[k] = -g_0_z_xxyy_xxxxxyy[k] * ab_x + g_0_z_xxyy_xxxxxxyy[k];

                g_0_z_xxxyy_xxxxxyz[k] = -g_0_z_xxyy_xxxxxyz[k] * ab_x + g_0_z_xxyy_xxxxxxyz[k];

                g_0_z_xxxyy_xxxxxzz[k] = -g_0_z_xxyy_xxxxxzz[k] * ab_x + g_0_z_xxyy_xxxxxxzz[k];

                g_0_z_xxxyy_xxxxyyy[k] = -g_0_z_xxyy_xxxxyyy[k] * ab_x + g_0_z_xxyy_xxxxxyyy[k];

                g_0_z_xxxyy_xxxxyyz[k] = -g_0_z_xxyy_xxxxyyz[k] * ab_x + g_0_z_xxyy_xxxxxyyz[k];

                g_0_z_xxxyy_xxxxyzz[k] = -g_0_z_xxyy_xxxxyzz[k] * ab_x + g_0_z_xxyy_xxxxxyzz[k];

                g_0_z_xxxyy_xxxxzzz[k] = -g_0_z_xxyy_xxxxzzz[k] * ab_x + g_0_z_xxyy_xxxxxzzz[k];

                g_0_z_xxxyy_xxxyyyy[k] = -g_0_z_xxyy_xxxyyyy[k] * ab_x + g_0_z_xxyy_xxxxyyyy[k];

                g_0_z_xxxyy_xxxyyyz[k] = -g_0_z_xxyy_xxxyyyz[k] * ab_x + g_0_z_xxyy_xxxxyyyz[k];

                g_0_z_xxxyy_xxxyyzz[k] = -g_0_z_xxyy_xxxyyzz[k] * ab_x + g_0_z_xxyy_xxxxyyzz[k];

                g_0_z_xxxyy_xxxyzzz[k] = -g_0_z_xxyy_xxxyzzz[k] * ab_x + g_0_z_xxyy_xxxxyzzz[k];

                g_0_z_xxxyy_xxxzzzz[k] = -g_0_z_xxyy_xxxzzzz[k] * ab_x + g_0_z_xxyy_xxxxzzzz[k];

                g_0_z_xxxyy_xxyyyyy[k] = -g_0_z_xxyy_xxyyyyy[k] * ab_x + g_0_z_xxyy_xxxyyyyy[k];

                g_0_z_xxxyy_xxyyyyz[k] = -g_0_z_xxyy_xxyyyyz[k] * ab_x + g_0_z_xxyy_xxxyyyyz[k];

                g_0_z_xxxyy_xxyyyzz[k] = -g_0_z_xxyy_xxyyyzz[k] * ab_x + g_0_z_xxyy_xxxyyyzz[k];

                g_0_z_xxxyy_xxyyzzz[k] = -g_0_z_xxyy_xxyyzzz[k] * ab_x + g_0_z_xxyy_xxxyyzzz[k];

                g_0_z_xxxyy_xxyzzzz[k] = -g_0_z_xxyy_xxyzzzz[k] * ab_x + g_0_z_xxyy_xxxyzzzz[k];

                g_0_z_xxxyy_xxzzzzz[k] = -g_0_z_xxyy_xxzzzzz[k] * ab_x + g_0_z_xxyy_xxxzzzzz[k];

                g_0_z_xxxyy_xyyyyyy[k] = -g_0_z_xxyy_xyyyyyy[k] * ab_x + g_0_z_xxyy_xxyyyyyy[k];

                g_0_z_xxxyy_xyyyyyz[k] = -g_0_z_xxyy_xyyyyyz[k] * ab_x + g_0_z_xxyy_xxyyyyyz[k];

                g_0_z_xxxyy_xyyyyzz[k] = -g_0_z_xxyy_xyyyyzz[k] * ab_x + g_0_z_xxyy_xxyyyyzz[k];

                g_0_z_xxxyy_xyyyzzz[k] = -g_0_z_xxyy_xyyyzzz[k] * ab_x + g_0_z_xxyy_xxyyyzzz[k];

                g_0_z_xxxyy_xyyzzzz[k] = -g_0_z_xxyy_xyyzzzz[k] * ab_x + g_0_z_xxyy_xxyyzzzz[k];

                g_0_z_xxxyy_xyzzzzz[k] = -g_0_z_xxyy_xyzzzzz[k] * ab_x + g_0_z_xxyy_xxyzzzzz[k];

                g_0_z_xxxyy_xzzzzzz[k] = -g_0_z_xxyy_xzzzzzz[k] * ab_x + g_0_z_xxyy_xxzzzzzz[k];

                g_0_z_xxxyy_yyyyyyy[k] = -g_0_z_xxyy_yyyyyyy[k] * ab_x + g_0_z_xxyy_xyyyyyyy[k];

                g_0_z_xxxyy_yyyyyyz[k] = -g_0_z_xxyy_yyyyyyz[k] * ab_x + g_0_z_xxyy_xyyyyyyz[k];

                g_0_z_xxxyy_yyyyyzz[k] = -g_0_z_xxyy_yyyyyzz[k] * ab_x + g_0_z_xxyy_xyyyyyzz[k];

                g_0_z_xxxyy_yyyyzzz[k] = -g_0_z_xxyy_yyyyzzz[k] * ab_x + g_0_z_xxyy_xyyyyzzz[k];

                g_0_z_xxxyy_yyyzzzz[k] = -g_0_z_xxyy_yyyzzzz[k] * ab_x + g_0_z_xxyy_xyyyzzzz[k];

                g_0_z_xxxyy_yyzzzzz[k] = -g_0_z_xxyy_yyzzzzz[k] * ab_x + g_0_z_xxyy_xyyzzzzz[k];

                g_0_z_xxxyy_yzzzzzz[k] = -g_0_z_xxyy_yzzzzzz[k] * ab_x + g_0_z_xxyy_xyzzzzzz[k];

                g_0_z_xxxyy_zzzzzzz[k] = -g_0_z_xxyy_zzzzzzz[k] * ab_x + g_0_z_xxyy_xzzzzzzz[k];
            }

            /// Set up 1656-1692 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1656 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1657 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1658 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1659 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1660 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1661 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1662 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1663 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1664 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1665 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1666 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1667 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1668 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1669 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1670 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1671 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1672 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1673 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1674 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1675 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1676 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1677 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1678 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1679 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1680 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1681 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1682 * ccomps * dcomps);

            auto g_0_z_xxxyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1683 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1684 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1685 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1686 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1687 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1688 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1689 * ccomps * dcomps);

            auto g_0_z_xxxyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1690 * ccomps * dcomps);

            auto g_0_z_xxxyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1691 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyz_xxxxxxx, g_0_z_xxxyz_xxxxxxy, g_0_z_xxxyz_xxxxxxz, g_0_z_xxxyz_xxxxxyy, g_0_z_xxxyz_xxxxxyz, g_0_z_xxxyz_xxxxxzz, g_0_z_xxxyz_xxxxyyy, g_0_z_xxxyz_xxxxyyz, g_0_z_xxxyz_xxxxyzz, g_0_z_xxxyz_xxxxzzz, g_0_z_xxxyz_xxxyyyy, g_0_z_xxxyz_xxxyyyz, g_0_z_xxxyz_xxxyyzz, g_0_z_xxxyz_xxxyzzz, g_0_z_xxxyz_xxxzzzz, g_0_z_xxxyz_xxyyyyy, g_0_z_xxxyz_xxyyyyz, g_0_z_xxxyz_xxyyyzz, g_0_z_xxxyz_xxyyzzz, g_0_z_xxxyz_xxyzzzz, g_0_z_xxxyz_xxzzzzz, g_0_z_xxxyz_xyyyyyy, g_0_z_xxxyz_xyyyyyz, g_0_z_xxxyz_xyyyyzz, g_0_z_xxxyz_xyyyzzz, g_0_z_xxxyz_xyyzzzz, g_0_z_xxxyz_xyzzzzz, g_0_z_xxxyz_xzzzzzz, g_0_z_xxxyz_yyyyyyy, g_0_z_xxxyz_yyyyyyz, g_0_z_xxxyz_yyyyyzz, g_0_z_xxxyz_yyyyzzz, g_0_z_xxxyz_yyyzzzz, g_0_z_xxxyz_yyzzzzz, g_0_z_xxxyz_yzzzzzz, g_0_z_xxxyz_zzzzzzz, g_0_z_xxyz_xxxxxxx, g_0_z_xxyz_xxxxxxxx, g_0_z_xxyz_xxxxxxxy, g_0_z_xxyz_xxxxxxxz, g_0_z_xxyz_xxxxxxy, g_0_z_xxyz_xxxxxxyy, g_0_z_xxyz_xxxxxxyz, g_0_z_xxyz_xxxxxxz, g_0_z_xxyz_xxxxxxzz, g_0_z_xxyz_xxxxxyy, g_0_z_xxyz_xxxxxyyy, g_0_z_xxyz_xxxxxyyz, g_0_z_xxyz_xxxxxyz, g_0_z_xxyz_xxxxxyzz, g_0_z_xxyz_xxxxxzz, g_0_z_xxyz_xxxxxzzz, g_0_z_xxyz_xxxxyyy, g_0_z_xxyz_xxxxyyyy, g_0_z_xxyz_xxxxyyyz, g_0_z_xxyz_xxxxyyz, g_0_z_xxyz_xxxxyyzz, g_0_z_xxyz_xxxxyzz, g_0_z_xxyz_xxxxyzzz, g_0_z_xxyz_xxxxzzz, g_0_z_xxyz_xxxxzzzz, g_0_z_xxyz_xxxyyyy, g_0_z_xxyz_xxxyyyyy, g_0_z_xxyz_xxxyyyyz, g_0_z_xxyz_xxxyyyz, g_0_z_xxyz_xxxyyyzz, g_0_z_xxyz_xxxyyzz, g_0_z_xxyz_xxxyyzzz, g_0_z_xxyz_xxxyzzz, g_0_z_xxyz_xxxyzzzz, g_0_z_xxyz_xxxzzzz, g_0_z_xxyz_xxxzzzzz, g_0_z_xxyz_xxyyyyy, g_0_z_xxyz_xxyyyyyy, g_0_z_xxyz_xxyyyyyz, g_0_z_xxyz_xxyyyyz, g_0_z_xxyz_xxyyyyzz, g_0_z_xxyz_xxyyyzz, g_0_z_xxyz_xxyyyzzz, g_0_z_xxyz_xxyyzzz, g_0_z_xxyz_xxyyzzzz, g_0_z_xxyz_xxyzzzz, g_0_z_xxyz_xxyzzzzz, g_0_z_xxyz_xxzzzzz, g_0_z_xxyz_xxzzzzzz, g_0_z_xxyz_xyyyyyy, g_0_z_xxyz_xyyyyyyy, g_0_z_xxyz_xyyyyyyz, g_0_z_xxyz_xyyyyyz, g_0_z_xxyz_xyyyyyzz, g_0_z_xxyz_xyyyyzz, g_0_z_xxyz_xyyyyzzz, g_0_z_xxyz_xyyyzzz, g_0_z_xxyz_xyyyzzzz, g_0_z_xxyz_xyyzzzz, g_0_z_xxyz_xyyzzzzz, g_0_z_xxyz_xyzzzzz, g_0_z_xxyz_xyzzzzzz, g_0_z_xxyz_xzzzzzz, g_0_z_xxyz_xzzzzzzz, g_0_z_xxyz_yyyyyyy, g_0_z_xxyz_yyyyyyz, g_0_z_xxyz_yyyyyzz, g_0_z_xxyz_yyyyzzz, g_0_z_xxyz_yyyzzzz, g_0_z_xxyz_yyzzzzz, g_0_z_xxyz_yzzzzzz, g_0_z_xxyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyz_xxxxxxx[k] = -g_0_z_xxyz_xxxxxxx[k] * ab_x + g_0_z_xxyz_xxxxxxxx[k];

                g_0_z_xxxyz_xxxxxxy[k] = -g_0_z_xxyz_xxxxxxy[k] * ab_x + g_0_z_xxyz_xxxxxxxy[k];

                g_0_z_xxxyz_xxxxxxz[k] = -g_0_z_xxyz_xxxxxxz[k] * ab_x + g_0_z_xxyz_xxxxxxxz[k];

                g_0_z_xxxyz_xxxxxyy[k] = -g_0_z_xxyz_xxxxxyy[k] * ab_x + g_0_z_xxyz_xxxxxxyy[k];

                g_0_z_xxxyz_xxxxxyz[k] = -g_0_z_xxyz_xxxxxyz[k] * ab_x + g_0_z_xxyz_xxxxxxyz[k];

                g_0_z_xxxyz_xxxxxzz[k] = -g_0_z_xxyz_xxxxxzz[k] * ab_x + g_0_z_xxyz_xxxxxxzz[k];

                g_0_z_xxxyz_xxxxyyy[k] = -g_0_z_xxyz_xxxxyyy[k] * ab_x + g_0_z_xxyz_xxxxxyyy[k];

                g_0_z_xxxyz_xxxxyyz[k] = -g_0_z_xxyz_xxxxyyz[k] * ab_x + g_0_z_xxyz_xxxxxyyz[k];

                g_0_z_xxxyz_xxxxyzz[k] = -g_0_z_xxyz_xxxxyzz[k] * ab_x + g_0_z_xxyz_xxxxxyzz[k];

                g_0_z_xxxyz_xxxxzzz[k] = -g_0_z_xxyz_xxxxzzz[k] * ab_x + g_0_z_xxyz_xxxxxzzz[k];

                g_0_z_xxxyz_xxxyyyy[k] = -g_0_z_xxyz_xxxyyyy[k] * ab_x + g_0_z_xxyz_xxxxyyyy[k];

                g_0_z_xxxyz_xxxyyyz[k] = -g_0_z_xxyz_xxxyyyz[k] * ab_x + g_0_z_xxyz_xxxxyyyz[k];

                g_0_z_xxxyz_xxxyyzz[k] = -g_0_z_xxyz_xxxyyzz[k] * ab_x + g_0_z_xxyz_xxxxyyzz[k];

                g_0_z_xxxyz_xxxyzzz[k] = -g_0_z_xxyz_xxxyzzz[k] * ab_x + g_0_z_xxyz_xxxxyzzz[k];

                g_0_z_xxxyz_xxxzzzz[k] = -g_0_z_xxyz_xxxzzzz[k] * ab_x + g_0_z_xxyz_xxxxzzzz[k];

                g_0_z_xxxyz_xxyyyyy[k] = -g_0_z_xxyz_xxyyyyy[k] * ab_x + g_0_z_xxyz_xxxyyyyy[k];

                g_0_z_xxxyz_xxyyyyz[k] = -g_0_z_xxyz_xxyyyyz[k] * ab_x + g_0_z_xxyz_xxxyyyyz[k];

                g_0_z_xxxyz_xxyyyzz[k] = -g_0_z_xxyz_xxyyyzz[k] * ab_x + g_0_z_xxyz_xxxyyyzz[k];

                g_0_z_xxxyz_xxyyzzz[k] = -g_0_z_xxyz_xxyyzzz[k] * ab_x + g_0_z_xxyz_xxxyyzzz[k];

                g_0_z_xxxyz_xxyzzzz[k] = -g_0_z_xxyz_xxyzzzz[k] * ab_x + g_0_z_xxyz_xxxyzzzz[k];

                g_0_z_xxxyz_xxzzzzz[k] = -g_0_z_xxyz_xxzzzzz[k] * ab_x + g_0_z_xxyz_xxxzzzzz[k];

                g_0_z_xxxyz_xyyyyyy[k] = -g_0_z_xxyz_xyyyyyy[k] * ab_x + g_0_z_xxyz_xxyyyyyy[k];

                g_0_z_xxxyz_xyyyyyz[k] = -g_0_z_xxyz_xyyyyyz[k] * ab_x + g_0_z_xxyz_xxyyyyyz[k];

                g_0_z_xxxyz_xyyyyzz[k] = -g_0_z_xxyz_xyyyyzz[k] * ab_x + g_0_z_xxyz_xxyyyyzz[k];

                g_0_z_xxxyz_xyyyzzz[k] = -g_0_z_xxyz_xyyyzzz[k] * ab_x + g_0_z_xxyz_xxyyyzzz[k];

                g_0_z_xxxyz_xyyzzzz[k] = -g_0_z_xxyz_xyyzzzz[k] * ab_x + g_0_z_xxyz_xxyyzzzz[k];

                g_0_z_xxxyz_xyzzzzz[k] = -g_0_z_xxyz_xyzzzzz[k] * ab_x + g_0_z_xxyz_xxyzzzzz[k];

                g_0_z_xxxyz_xzzzzzz[k] = -g_0_z_xxyz_xzzzzzz[k] * ab_x + g_0_z_xxyz_xxzzzzzz[k];

                g_0_z_xxxyz_yyyyyyy[k] = -g_0_z_xxyz_yyyyyyy[k] * ab_x + g_0_z_xxyz_xyyyyyyy[k];

                g_0_z_xxxyz_yyyyyyz[k] = -g_0_z_xxyz_yyyyyyz[k] * ab_x + g_0_z_xxyz_xyyyyyyz[k];

                g_0_z_xxxyz_yyyyyzz[k] = -g_0_z_xxyz_yyyyyzz[k] * ab_x + g_0_z_xxyz_xyyyyyzz[k];

                g_0_z_xxxyz_yyyyzzz[k] = -g_0_z_xxyz_yyyyzzz[k] * ab_x + g_0_z_xxyz_xyyyyzzz[k];

                g_0_z_xxxyz_yyyzzzz[k] = -g_0_z_xxyz_yyyzzzz[k] * ab_x + g_0_z_xxyz_xyyyzzzz[k];

                g_0_z_xxxyz_yyzzzzz[k] = -g_0_z_xxyz_yyzzzzz[k] * ab_x + g_0_z_xxyz_xyyzzzzz[k];

                g_0_z_xxxyz_yzzzzzz[k] = -g_0_z_xxyz_yzzzzzz[k] * ab_x + g_0_z_xxyz_xyzzzzzz[k];

                g_0_z_xxxyz_zzzzzzz[k] = -g_0_z_xxyz_zzzzzzz[k] * ab_x + g_0_z_xxyz_xzzzzzzz[k];
            }

            /// Set up 1692-1728 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1692 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1693 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1694 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1695 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1696 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1697 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1698 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1699 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1700 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1701 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1702 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1703 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1704 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1705 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1706 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1707 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1708 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1709 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1710 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1711 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1712 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1713 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1714 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1715 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1716 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1717 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1718 * ccomps * dcomps);

            auto g_0_z_xxxzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1719 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1720 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1721 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1722 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1723 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1724 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1725 * ccomps * dcomps);

            auto g_0_z_xxxzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1726 * ccomps * dcomps);

            auto g_0_z_xxxzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1727 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzz_xxxxxxx, g_0_z_xxxzz_xxxxxxy, g_0_z_xxxzz_xxxxxxz, g_0_z_xxxzz_xxxxxyy, g_0_z_xxxzz_xxxxxyz, g_0_z_xxxzz_xxxxxzz, g_0_z_xxxzz_xxxxyyy, g_0_z_xxxzz_xxxxyyz, g_0_z_xxxzz_xxxxyzz, g_0_z_xxxzz_xxxxzzz, g_0_z_xxxzz_xxxyyyy, g_0_z_xxxzz_xxxyyyz, g_0_z_xxxzz_xxxyyzz, g_0_z_xxxzz_xxxyzzz, g_0_z_xxxzz_xxxzzzz, g_0_z_xxxzz_xxyyyyy, g_0_z_xxxzz_xxyyyyz, g_0_z_xxxzz_xxyyyzz, g_0_z_xxxzz_xxyyzzz, g_0_z_xxxzz_xxyzzzz, g_0_z_xxxzz_xxzzzzz, g_0_z_xxxzz_xyyyyyy, g_0_z_xxxzz_xyyyyyz, g_0_z_xxxzz_xyyyyzz, g_0_z_xxxzz_xyyyzzz, g_0_z_xxxzz_xyyzzzz, g_0_z_xxxzz_xyzzzzz, g_0_z_xxxzz_xzzzzzz, g_0_z_xxxzz_yyyyyyy, g_0_z_xxxzz_yyyyyyz, g_0_z_xxxzz_yyyyyzz, g_0_z_xxxzz_yyyyzzz, g_0_z_xxxzz_yyyzzzz, g_0_z_xxxzz_yyzzzzz, g_0_z_xxxzz_yzzzzzz, g_0_z_xxxzz_zzzzzzz, g_0_z_xxzz_xxxxxxx, g_0_z_xxzz_xxxxxxxx, g_0_z_xxzz_xxxxxxxy, g_0_z_xxzz_xxxxxxxz, g_0_z_xxzz_xxxxxxy, g_0_z_xxzz_xxxxxxyy, g_0_z_xxzz_xxxxxxyz, g_0_z_xxzz_xxxxxxz, g_0_z_xxzz_xxxxxxzz, g_0_z_xxzz_xxxxxyy, g_0_z_xxzz_xxxxxyyy, g_0_z_xxzz_xxxxxyyz, g_0_z_xxzz_xxxxxyz, g_0_z_xxzz_xxxxxyzz, g_0_z_xxzz_xxxxxzz, g_0_z_xxzz_xxxxxzzz, g_0_z_xxzz_xxxxyyy, g_0_z_xxzz_xxxxyyyy, g_0_z_xxzz_xxxxyyyz, g_0_z_xxzz_xxxxyyz, g_0_z_xxzz_xxxxyyzz, g_0_z_xxzz_xxxxyzz, g_0_z_xxzz_xxxxyzzz, g_0_z_xxzz_xxxxzzz, g_0_z_xxzz_xxxxzzzz, g_0_z_xxzz_xxxyyyy, g_0_z_xxzz_xxxyyyyy, g_0_z_xxzz_xxxyyyyz, g_0_z_xxzz_xxxyyyz, g_0_z_xxzz_xxxyyyzz, g_0_z_xxzz_xxxyyzz, g_0_z_xxzz_xxxyyzzz, g_0_z_xxzz_xxxyzzz, g_0_z_xxzz_xxxyzzzz, g_0_z_xxzz_xxxzzzz, g_0_z_xxzz_xxxzzzzz, g_0_z_xxzz_xxyyyyy, g_0_z_xxzz_xxyyyyyy, g_0_z_xxzz_xxyyyyyz, g_0_z_xxzz_xxyyyyz, g_0_z_xxzz_xxyyyyzz, g_0_z_xxzz_xxyyyzz, g_0_z_xxzz_xxyyyzzz, g_0_z_xxzz_xxyyzzz, g_0_z_xxzz_xxyyzzzz, g_0_z_xxzz_xxyzzzz, g_0_z_xxzz_xxyzzzzz, g_0_z_xxzz_xxzzzzz, g_0_z_xxzz_xxzzzzzz, g_0_z_xxzz_xyyyyyy, g_0_z_xxzz_xyyyyyyy, g_0_z_xxzz_xyyyyyyz, g_0_z_xxzz_xyyyyyz, g_0_z_xxzz_xyyyyyzz, g_0_z_xxzz_xyyyyzz, g_0_z_xxzz_xyyyyzzz, g_0_z_xxzz_xyyyzzz, g_0_z_xxzz_xyyyzzzz, g_0_z_xxzz_xyyzzzz, g_0_z_xxzz_xyyzzzzz, g_0_z_xxzz_xyzzzzz, g_0_z_xxzz_xyzzzzzz, g_0_z_xxzz_xzzzzzz, g_0_z_xxzz_xzzzzzzz, g_0_z_xxzz_yyyyyyy, g_0_z_xxzz_yyyyyyz, g_0_z_xxzz_yyyyyzz, g_0_z_xxzz_yyyyzzz, g_0_z_xxzz_yyyzzzz, g_0_z_xxzz_yyzzzzz, g_0_z_xxzz_yzzzzzz, g_0_z_xxzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzz_xxxxxxx[k] = -g_0_z_xxzz_xxxxxxx[k] * ab_x + g_0_z_xxzz_xxxxxxxx[k];

                g_0_z_xxxzz_xxxxxxy[k] = -g_0_z_xxzz_xxxxxxy[k] * ab_x + g_0_z_xxzz_xxxxxxxy[k];

                g_0_z_xxxzz_xxxxxxz[k] = -g_0_z_xxzz_xxxxxxz[k] * ab_x + g_0_z_xxzz_xxxxxxxz[k];

                g_0_z_xxxzz_xxxxxyy[k] = -g_0_z_xxzz_xxxxxyy[k] * ab_x + g_0_z_xxzz_xxxxxxyy[k];

                g_0_z_xxxzz_xxxxxyz[k] = -g_0_z_xxzz_xxxxxyz[k] * ab_x + g_0_z_xxzz_xxxxxxyz[k];

                g_0_z_xxxzz_xxxxxzz[k] = -g_0_z_xxzz_xxxxxzz[k] * ab_x + g_0_z_xxzz_xxxxxxzz[k];

                g_0_z_xxxzz_xxxxyyy[k] = -g_0_z_xxzz_xxxxyyy[k] * ab_x + g_0_z_xxzz_xxxxxyyy[k];

                g_0_z_xxxzz_xxxxyyz[k] = -g_0_z_xxzz_xxxxyyz[k] * ab_x + g_0_z_xxzz_xxxxxyyz[k];

                g_0_z_xxxzz_xxxxyzz[k] = -g_0_z_xxzz_xxxxyzz[k] * ab_x + g_0_z_xxzz_xxxxxyzz[k];

                g_0_z_xxxzz_xxxxzzz[k] = -g_0_z_xxzz_xxxxzzz[k] * ab_x + g_0_z_xxzz_xxxxxzzz[k];

                g_0_z_xxxzz_xxxyyyy[k] = -g_0_z_xxzz_xxxyyyy[k] * ab_x + g_0_z_xxzz_xxxxyyyy[k];

                g_0_z_xxxzz_xxxyyyz[k] = -g_0_z_xxzz_xxxyyyz[k] * ab_x + g_0_z_xxzz_xxxxyyyz[k];

                g_0_z_xxxzz_xxxyyzz[k] = -g_0_z_xxzz_xxxyyzz[k] * ab_x + g_0_z_xxzz_xxxxyyzz[k];

                g_0_z_xxxzz_xxxyzzz[k] = -g_0_z_xxzz_xxxyzzz[k] * ab_x + g_0_z_xxzz_xxxxyzzz[k];

                g_0_z_xxxzz_xxxzzzz[k] = -g_0_z_xxzz_xxxzzzz[k] * ab_x + g_0_z_xxzz_xxxxzzzz[k];

                g_0_z_xxxzz_xxyyyyy[k] = -g_0_z_xxzz_xxyyyyy[k] * ab_x + g_0_z_xxzz_xxxyyyyy[k];

                g_0_z_xxxzz_xxyyyyz[k] = -g_0_z_xxzz_xxyyyyz[k] * ab_x + g_0_z_xxzz_xxxyyyyz[k];

                g_0_z_xxxzz_xxyyyzz[k] = -g_0_z_xxzz_xxyyyzz[k] * ab_x + g_0_z_xxzz_xxxyyyzz[k];

                g_0_z_xxxzz_xxyyzzz[k] = -g_0_z_xxzz_xxyyzzz[k] * ab_x + g_0_z_xxzz_xxxyyzzz[k];

                g_0_z_xxxzz_xxyzzzz[k] = -g_0_z_xxzz_xxyzzzz[k] * ab_x + g_0_z_xxzz_xxxyzzzz[k];

                g_0_z_xxxzz_xxzzzzz[k] = -g_0_z_xxzz_xxzzzzz[k] * ab_x + g_0_z_xxzz_xxxzzzzz[k];

                g_0_z_xxxzz_xyyyyyy[k] = -g_0_z_xxzz_xyyyyyy[k] * ab_x + g_0_z_xxzz_xxyyyyyy[k];

                g_0_z_xxxzz_xyyyyyz[k] = -g_0_z_xxzz_xyyyyyz[k] * ab_x + g_0_z_xxzz_xxyyyyyz[k];

                g_0_z_xxxzz_xyyyyzz[k] = -g_0_z_xxzz_xyyyyzz[k] * ab_x + g_0_z_xxzz_xxyyyyzz[k];

                g_0_z_xxxzz_xyyyzzz[k] = -g_0_z_xxzz_xyyyzzz[k] * ab_x + g_0_z_xxzz_xxyyyzzz[k];

                g_0_z_xxxzz_xyyzzzz[k] = -g_0_z_xxzz_xyyzzzz[k] * ab_x + g_0_z_xxzz_xxyyzzzz[k];

                g_0_z_xxxzz_xyzzzzz[k] = -g_0_z_xxzz_xyzzzzz[k] * ab_x + g_0_z_xxzz_xxyzzzzz[k];

                g_0_z_xxxzz_xzzzzzz[k] = -g_0_z_xxzz_xzzzzzz[k] * ab_x + g_0_z_xxzz_xxzzzzzz[k];

                g_0_z_xxxzz_yyyyyyy[k] = -g_0_z_xxzz_yyyyyyy[k] * ab_x + g_0_z_xxzz_xyyyyyyy[k];

                g_0_z_xxxzz_yyyyyyz[k] = -g_0_z_xxzz_yyyyyyz[k] * ab_x + g_0_z_xxzz_xyyyyyyz[k];

                g_0_z_xxxzz_yyyyyzz[k] = -g_0_z_xxzz_yyyyyzz[k] * ab_x + g_0_z_xxzz_xyyyyyzz[k];

                g_0_z_xxxzz_yyyyzzz[k] = -g_0_z_xxzz_yyyyzzz[k] * ab_x + g_0_z_xxzz_xyyyyzzz[k];

                g_0_z_xxxzz_yyyzzzz[k] = -g_0_z_xxzz_yyyzzzz[k] * ab_x + g_0_z_xxzz_xyyyzzzz[k];

                g_0_z_xxxzz_yyzzzzz[k] = -g_0_z_xxzz_yyzzzzz[k] * ab_x + g_0_z_xxzz_xyyzzzzz[k];

                g_0_z_xxxzz_yzzzzzz[k] = -g_0_z_xxzz_yzzzzzz[k] * ab_x + g_0_z_xxzz_xyzzzzzz[k];

                g_0_z_xxxzz_zzzzzzz[k] = -g_0_z_xxzz_zzzzzzz[k] * ab_x + g_0_z_xxzz_xzzzzzzz[k];
            }

            /// Set up 1728-1764 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1728 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1729 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1730 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1731 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1732 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1733 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1734 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1735 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1736 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1737 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1738 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1739 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1740 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1741 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1742 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1743 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1744 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1745 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1746 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1747 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1748 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1749 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1750 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1751 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1752 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1753 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1754 * ccomps * dcomps);

            auto g_0_z_xxyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1755 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1756 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1757 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1758 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1759 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1760 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1761 * ccomps * dcomps);

            auto g_0_z_xxyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1762 * ccomps * dcomps);

            auto g_0_z_xxyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1763 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyy_xxxxxxx, g_0_z_xxyyy_xxxxxxy, g_0_z_xxyyy_xxxxxxz, g_0_z_xxyyy_xxxxxyy, g_0_z_xxyyy_xxxxxyz, g_0_z_xxyyy_xxxxxzz, g_0_z_xxyyy_xxxxyyy, g_0_z_xxyyy_xxxxyyz, g_0_z_xxyyy_xxxxyzz, g_0_z_xxyyy_xxxxzzz, g_0_z_xxyyy_xxxyyyy, g_0_z_xxyyy_xxxyyyz, g_0_z_xxyyy_xxxyyzz, g_0_z_xxyyy_xxxyzzz, g_0_z_xxyyy_xxxzzzz, g_0_z_xxyyy_xxyyyyy, g_0_z_xxyyy_xxyyyyz, g_0_z_xxyyy_xxyyyzz, g_0_z_xxyyy_xxyyzzz, g_0_z_xxyyy_xxyzzzz, g_0_z_xxyyy_xxzzzzz, g_0_z_xxyyy_xyyyyyy, g_0_z_xxyyy_xyyyyyz, g_0_z_xxyyy_xyyyyzz, g_0_z_xxyyy_xyyyzzz, g_0_z_xxyyy_xyyzzzz, g_0_z_xxyyy_xyzzzzz, g_0_z_xxyyy_xzzzzzz, g_0_z_xxyyy_yyyyyyy, g_0_z_xxyyy_yyyyyyz, g_0_z_xxyyy_yyyyyzz, g_0_z_xxyyy_yyyyzzz, g_0_z_xxyyy_yyyzzzz, g_0_z_xxyyy_yyzzzzz, g_0_z_xxyyy_yzzzzzz, g_0_z_xxyyy_zzzzzzz, g_0_z_xyyy_xxxxxxx, g_0_z_xyyy_xxxxxxxx, g_0_z_xyyy_xxxxxxxy, g_0_z_xyyy_xxxxxxxz, g_0_z_xyyy_xxxxxxy, g_0_z_xyyy_xxxxxxyy, g_0_z_xyyy_xxxxxxyz, g_0_z_xyyy_xxxxxxz, g_0_z_xyyy_xxxxxxzz, g_0_z_xyyy_xxxxxyy, g_0_z_xyyy_xxxxxyyy, g_0_z_xyyy_xxxxxyyz, g_0_z_xyyy_xxxxxyz, g_0_z_xyyy_xxxxxyzz, g_0_z_xyyy_xxxxxzz, g_0_z_xyyy_xxxxxzzz, g_0_z_xyyy_xxxxyyy, g_0_z_xyyy_xxxxyyyy, g_0_z_xyyy_xxxxyyyz, g_0_z_xyyy_xxxxyyz, g_0_z_xyyy_xxxxyyzz, g_0_z_xyyy_xxxxyzz, g_0_z_xyyy_xxxxyzzz, g_0_z_xyyy_xxxxzzz, g_0_z_xyyy_xxxxzzzz, g_0_z_xyyy_xxxyyyy, g_0_z_xyyy_xxxyyyyy, g_0_z_xyyy_xxxyyyyz, g_0_z_xyyy_xxxyyyz, g_0_z_xyyy_xxxyyyzz, g_0_z_xyyy_xxxyyzz, g_0_z_xyyy_xxxyyzzz, g_0_z_xyyy_xxxyzzz, g_0_z_xyyy_xxxyzzzz, g_0_z_xyyy_xxxzzzz, g_0_z_xyyy_xxxzzzzz, g_0_z_xyyy_xxyyyyy, g_0_z_xyyy_xxyyyyyy, g_0_z_xyyy_xxyyyyyz, g_0_z_xyyy_xxyyyyz, g_0_z_xyyy_xxyyyyzz, g_0_z_xyyy_xxyyyzz, g_0_z_xyyy_xxyyyzzz, g_0_z_xyyy_xxyyzzz, g_0_z_xyyy_xxyyzzzz, g_0_z_xyyy_xxyzzzz, g_0_z_xyyy_xxyzzzzz, g_0_z_xyyy_xxzzzzz, g_0_z_xyyy_xxzzzzzz, g_0_z_xyyy_xyyyyyy, g_0_z_xyyy_xyyyyyyy, g_0_z_xyyy_xyyyyyyz, g_0_z_xyyy_xyyyyyz, g_0_z_xyyy_xyyyyyzz, g_0_z_xyyy_xyyyyzz, g_0_z_xyyy_xyyyyzzz, g_0_z_xyyy_xyyyzzz, g_0_z_xyyy_xyyyzzzz, g_0_z_xyyy_xyyzzzz, g_0_z_xyyy_xyyzzzzz, g_0_z_xyyy_xyzzzzz, g_0_z_xyyy_xyzzzzzz, g_0_z_xyyy_xzzzzzz, g_0_z_xyyy_xzzzzzzz, g_0_z_xyyy_yyyyyyy, g_0_z_xyyy_yyyyyyz, g_0_z_xyyy_yyyyyzz, g_0_z_xyyy_yyyyzzz, g_0_z_xyyy_yyyzzzz, g_0_z_xyyy_yyzzzzz, g_0_z_xyyy_yzzzzzz, g_0_z_xyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyy_xxxxxxx[k] = -g_0_z_xyyy_xxxxxxx[k] * ab_x + g_0_z_xyyy_xxxxxxxx[k];

                g_0_z_xxyyy_xxxxxxy[k] = -g_0_z_xyyy_xxxxxxy[k] * ab_x + g_0_z_xyyy_xxxxxxxy[k];

                g_0_z_xxyyy_xxxxxxz[k] = -g_0_z_xyyy_xxxxxxz[k] * ab_x + g_0_z_xyyy_xxxxxxxz[k];

                g_0_z_xxyyy_xxxxxyy[k] = -g_0_z_xyyy_xxxxxyy[k] * ab_x + g_0_z_xyyy_xxxxxxyy[k];

                g_0_z_xxyyy_xxxxxyz[k] = -g_0_z_xyyy_xxxxxyz[k] * ab_x + g_0_z_xyyy_xxxxxxyz[k];

                g_0_z_xxyyy_xxxxxzz[k] = -g_0_z_xyyy_xxxxxzz[k] * ab_x + g_0_z_xyyy_xxxxxxzz[k];

                g_0_z_xxyyy_xxxxyyy[k] = -g_0_z_xyyy_xxxxyyy[k] * ab_x + g_0_z_xyyy_xxxxxyyy[k];

                g_0_z_xxyyy_xxxxyyz[k] = -g_0_z_xyyy_xxxxyyz[k] * ab_x + g_0_z_xyyy_xxxxxyyz[k];

                g_0_z_xxyyy_xxxxyzz[k] = -g_0_z_xyyy_xxxxyzz[k] * ab_x + g_0_z_xyyy_xxxxxyzz[k];

                g_0_z_xxyyy_xxxxzzz[k] = -g_0_z_xyyy_xxxxzzz[k] * ab_x + g_0_z_xyyy_xxxxxzzz[k];

                g_0_z_xxyyy_xxxyyyy[k] = -g_0_z_xyyy_xxxyyyy[k] * ab_x + g_0_z_xyyy_xxxxyyyy[k];

                g_0_z_xxyyy_xxxyyyz[k] = -g_0_z_xyyy_xxxyyyz[k] * ab_x + g_0_z_xyyy_xxxxyyyz[k];

                g_0_z_xxyyy_xxxyyzz[k] = -g_0_z_xyyy_xxxyyzz[k] * ab_x + g_0_z_xyyy_xxxxyyzz[k];

                g_0_z_xxyyy_xxxyzzz[k] = -g_0_z_xyyy_xxxyzzz[k] * ab_x + g_0_z_xyyy_xxxxyzzz[k];

                g_0_z_xxyyy_xxxzzzz[k] = -g_0_z_xyyy_xxxzzzz[k] * ab_x + g_0_z_xyyy_xxxxzzzz[k];

                g_0_z_xxyyy_xxyyyyy[k] = -g_0_z_xyyy_xxyyyyy[k] * ab_x + g_0_z_xyyy_xxxyyyyy[k];

                g_0_z_xxyyy_xxyyyyz[k] = -g_0_z_xyyy_xxyyyyz[k] * ab_x + g_0_z_xyyy_xxxyyyyz[k];

                g_0_z_xxyyy_xxyyyzz[k] = -g_0_z_xyyy_xxyyyzz[k] * ab_x + g_0_z_xyyy_xxxyyyzz[k];

                g_0_z_xxyyy_xxyyzzz[k] = -g_0_z_xyyy_xxyyzzz[k] * ab_x + g_0_z_xyyy_xxxyyzzz[k];

                g_0_z_xxyyy_xxyzzzz[k] = -g_0_z_xyyy_xxyzzzz[k] * ab_x + g_0_z_xyyy_xxxyzzzz[k];

                g_0_z_xxyyy_xxzzzzz[k] = -g_0_z_xyyy_xxzzzzz[k] * ab_x + g_0_z_xyyy_xxxzzzzz[k];

                g_0_z_xxyyy_xyyyyyy[k] = -g_0_z_xyyy_xyyyyyy[k] * ab_x + g_0_z_xyyy_xxyyyyyy[k];

                g_0_z_xxyyy_xyyyyyz[k] = -g_0_z_xyyy_xyyyyyz[k] * ab_x + g_0_z_xyyy_xxyyyyyz[k];

                g_0_z_xxyyy_xyyyyzz[k] = -g_0_z_xyyy_xyyyyzz[k] * ab_x + g_0_z_xyyy_xxyyyyzz[k];

                g_0_z_xxyyy_xyyyzzz[k] = -g_0_z_xyyy_xyyyzzz[k] * ab_x + g_0_z_xyyy_xxyyyzzz[k];

                g_0_z_xxyyy_xyyzzzz[k] = -g_0_z_xyyy_xyyzzzz[k] * ab_x + g_0_z_xyyy_xxyyzzzz[k];

                g_0_z_xxyyy_xyzzzzz[k] = -g_0_z_xyyy_xyzzzzz[k] * ab_x + g_0_z_xyyy_xxyzzzzz[k];

                g_0_z_xxyyy_xzzzzzz[k] = -g_0_z_xyyy_xzzzzzz[k] * ab_x + g_0_z_xyyy_xxzzzzzz[k];

                g_0_z_xxyyy_yyyyyyy[k] = -g_0_z_xyyy_yyyyyyy[k] * ab_x + g_0_z_xyyy_xyyyyyyy[k];

                g_0_z_xxyyy_yyyyyyz[k] = -g_0_z_xyyy_yyyyyyz[k] * ab_x + g_0_z_xyyy_xyyyyyyz[k];

                g_0_z_xxyyy_yyyyyzz[k] = -g_0_z_xyyy_yyyyyzz[k] * ab_x + g_0_z_xyyy_xyyyyyzz[k];

                g_0_z_xxyyy_yyyyzzz[k] = -g_0_z_xyyy_yyyyzzz[k] * ab_x + g_0_z_xyyy_xyyyyzzz[k];

                g_0_z_xxyyy_yyyzzzz[k] = -g_0_z_xyyy_yyyzzzz[k] * ab_x + g_0_z_xyyy_xyyyzzzz[k];

                g_0_z_xxyyy_yyzzzzz[k] = -g_0_z_xyyy_yyzzzzz[k] * ab_x + g_0_z_xyyy_xyyzzzzz[k];

                g_0_z_xxyyy_yzzzzzz[k] = -g_0_z_xyyy_yzzzzzz[k] * ab_x + g_0_z_xyyy_xyzzzzzz[k];

                g_0_z_xxyyy_zzzzzzz[k] = -g_0_z_xyyy_zzzzzzz[k] * ab_x + g_0_z_xyyy_xzzzzzzz[k];
            }

            /// Set up 1764-1800 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1764 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1765 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1766 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1767 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1768 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1769 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1770 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1771 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1772 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1773 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1774 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1775 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1776 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1777 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1778 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1779 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1780 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1781 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1782 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1783 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1784 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1785 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1786 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1787 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1788 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1789 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1790 * ccomps * dcomps);

            auto g_0_z_xxyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1791 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1792 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1793 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1794 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1795 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1796 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1797 * ccomps * dcomps);

            auto g_0_z_xxyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1798 * ccomps * dcomps);

            auto g_0_z_xxyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1799 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyz_xxxxxxx, g_0_z_xxyyz_xxxxxxy, g_0_z_xxyyz_xxxxxxz, g_0_z_xxyyz_xxxxxyy, g_0_z_xxyyz_xxxxxyz, g_0_z_xxyyz_xxxxxzz, g_0_z_xxyyz_xxxxyyy, g_0_z_xxyyz_xxxxyyz, g_0_z_xxyyz_xxxxyzz, g_0_z_xxyyz_xxxxzzz, g_0_z_xxyyz_xxxyyyy, g_0_z_xxyyz_xxxyyyz, g_0_z_xxyyz_xxxyyzz, g_0_z_xxyyz_xxxyzzz, g_0_z_xxyyz_xxxzzzz, g_0_z_xxyyz_xxyyyyy, g_0_z_xxyyz_xxyyyyz, g_0_z_xxyyz_xxyyyzz, g_0_z_xxyyz_xxyyzzz, g_0_z_xxyyz_xxyzzzz, g_0_z_xxyyz_xxzzzzz, g_0_z_xxyyz_xyyyyyy, g_0_z_xxyyz_xyyyyyz, g_0_z_xxyyz_xyyyyzz, g_0_z_xxyyz_xyyyzzz, g_0_z_xxyyz_xyyzzzz, g_0_z_xxyyz_xyzzzzz, g_0_z_xxyyz_xzzzzzz, g_0_z_xxyyz_yyyyyyy, g_0_z_xxyyz_yyyyyyz, g_0_z_xxyyz_yyyyyzz, g_0_z_xxyyz_yyyyzzz, g_0_z_xxyyz_yyyzzzz, g_0_z_xxyyz_yyzzzzz, g_0_z_xxyyz_yzzzzzz, g_0_z_xxyyz_zzzzzzz, g_0_z_xyyz_xxxxxxx, g_0_z_xyyz_xxxxxxxx, g_0_z_xyyz_xxxxxxxy, g_0_z_xyyz_xxxxxxxz, g_0_z_xyyz_xxxxxxy, g_0_z_xyyz_xxxxxxyy, g_0_z_xyyz_xxxxxxyz, g_0_z_xyyz_xxxxxxz, g_0_z_xyyz_xxxxxxzz, g_0_z_xyyz_xxxxxyy, g_0_z_xyyz_xxxxxyyy, g_0_z_xyyz_xxxxxyyz, g_0_z_xyyz_xxxxxyz, g_0_z_xyyz_xxxxxyzz, g_0_z_xyyz_xxxxxzz, g_0_z_xyyz_xxxxxzzz, g_0_z_xyyz_xxxxyyy, g_0_z_xyyz_xxxxyyyy, g_0_z_xyyz_xxxxyyyz, g_0_z_xyyz_xxxxyyz, g_0_z_xyyz_xxxxyyzz, g_0_z_xyyz_xxxxyzz, g_0_z_xyyz_xxxxyzzz, g_0_z_xyyz_xxxxzzz, g_0_z_xyyz_xxxxzzzz, g_0_z_xyyz_xxxyyyy, g_0_z_xyyz_xxxyyyyy, g_0_z_xyyz_xxxyyyyz, g_0_z_xyyz_xxxyyyz, g_0_z_xyyz_xxxyyyzz, g_0_z_xyyz_xxxyyzz, g_0_z_xyyz_xxxyyzzz, g_0_z_xyyz_xxxyzzz, g_0_z_xyyz_xxxyzzzz, g_0_z_xyyz_xxxzzzz, g_0_z_xyyz_xxxzzzzz, g_0_z_xyyz_xxyyyyy, g_0_z_xyyz_xxyyyyyy, g_0_z_xyyz_xxyyyyyz, g_0_z_xyyz_xxyyyyz, g_0_z_xyyz_xxyyyyzz, g_0_z_xyyz_xxyyyzz, g_0_z_xyyz_xxyyyzzz, g_0_z_xyyz_xxyyzzz, g_0_z_xyyz_xxyyzzzz, g_0_z_xyyz_xxyzzzz, g_0_z_xyyz_xxyzzzzz, g_0_z_xyyz_xxzzzzz, g_0_z_xyyz_xxzzzzzz, g_0_z_xyyz_xyyyyyy, g_0_z_xyyz_xyyyyyyy, g_0_z_xyyz_xyyyyyyz, g_0_z_xyyz_xyyyyyz, g_0_z_xyyz_xyyyyyzz, g_0_z_xyyz_xyyyyzz, g_0_z_xyyz_xyyyyzzz, g_0_z_xyyz_xyyyzzz, g_0_z_xyyz_xyyyzzzz, g_0_z_xyyz_xyyzzzz, g_0_z_xyyz_xyyzzzzz, g_0_z_xyyz_xyzzzzz, g_0_z_xyyz_xyzzzzzz, g_0_z_xyyz_xzzzzzz, g_0_z_xyyz_xzzzzzzz, g_0_z_xyyz_yyyyyyy, g_0_z_xyyz_yyyyyyz, g_0_z_xyyz_yyyyyzz, g_0_z_xyyz_yyyyzzz, g_0_z_xyyz_yyyzzzz, g_0_z_xyyz_yyzzzzz, g_0_z_xyyz_yzzzzzz, g_0_z_xyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyz_xxxxxxx[k] = -g_0_z_xyyz_xxxxxxx[k] * ab_x + g_0_z_xyyz_xxxxxxxx[k];

                g_0_z_xxyyz_xxxxxxy[k] = -g_0_z_xyyz_xxxxxxy[k] * ab_x + g_0_z_xyyz_xxxxxxxy[k];

                g_0_z_xxyyz_xxxxxxz[k] = -g_0_z_xyyz_xxxxxxz[k] * ab_x + g_0_z_xyyz_xxxxxxxz[k];

                g_0_z_xxyyz_xxxxxyy[k] = -g_0_z_xyyz_xxxxxyy[k] * ab_x + g_0_z_xyyz_xxxxxxyy[k];

                g_0_z_xxyyz_xxxxxyz[k] = -g_0_z_xyyz_xxxxxyz[k] * ab_x + g_0_z_xyyz_xxxxxxyz[k];

                g_0_z_xxyyz_xxxxxzz[k] = -g_0_z_xyyz_xxxxxzz[k] * ab_x + g_0_z_xyyz_xxxxxxzz[k];

                g_0_z_xxyyz_xxxxyyy[k] = -g_0_z_xyyz_xxxxyyy[k] * ab_x + g_0_z_xyyz_xxxxxyyy[k];

                g_0_z_xxyyz_xxxxyyz[k] = -g_0_z_xyyz_xxxxyyz[k] * ab_x + g_0_z_xyyz_xxxxxyyz[k];

                g_0_z_xxyyz_xxxxyzz[k] = -g_0_z_xyyz_xxxxyzz[k] * ab_x + g_0_z_xyyz_xxxxxyzz[k];

                g_0_z_xxyyz_xxxxzzz[k] = -g_0_z_xyyz_xxxxzzz[k] * ab_x + g_0_z_xyyz_xxxxxzzz[k];

                g_0_z_xxyyz_xxxyyyy[k] = -g_0_z_xyyz_xxxyyyy[k] * ab_x + g_0_z_xyyz_xxxxyyyy[k];

                g_0_z_xxyyz_xxxyyyz[k] = -g_0_z_xyyz_xxxyyyz[k] * ab_x + g_0_z_xyyz_xxxxyyyz[k];

                g_0_z_xxyyz_xxxyyzz[k] = -g_0_z_xyyz_xxxyyzz[k] * ab_x + g_0_z_xyyz_xxxxyyzz[k];

                g_0_z_xxyyz_xxxyzzz[k] = -g_0_z_xyyz_xxxyzzz[k] * ab_x + g_0_z_xyyz_xxxxyzzz[k];

                g_0_z_xxyyz_xxxzzzz[k] = -g_0_z_xyyz_xxxzzzz[k] * ab_x + g_0_z_xyyz_xxxxzzzz[k];

                g_0_z_xxyyz_xxyyyyy[k] = -g_0_z_xyyz_xxyyyyy[k] * ab_x + g_0_z_xyyz_xxxyyyyy[k];

                g_0_z_xxyyz_xxyyyyz[k] = -g_0_z_xyyz_xxyyyyz[k] * ab_x + g_0_z_xyyz_xxxyyyyz[k];

                g_0_z_xxyyz_xxyyyzz[k] = -g_0_z_xyyz_xxyyyzz[k] * ab_x + g_0_z_xyyz_xxxyyyzz[k];

                g_0_z_xxyyz_xxyyzzz[k] = -g_0_z_xyyz_xxyyzzz[k] * ab_x + g_0_z_xyyz_xxxyyzzz[k];

                g_0_z_xxyyz_xxyzzzz[k] = -g_0_z_xyyz_xxyzzzz[k] * ab_x + g_0_z_xyyz_xxxyzzzz[k];

                g_0_z_xxyyz_xxzzzzz[k] = -g_0_z_xyyz_xxzzzzz[k] * ab_x + g_0_z_xyyz_xxxzzzzz[k];

                g_0_z_xxyyz_xyyyyyy[k] = -g_0_z_xyyz_xyyyyyy[k] * ab_x + g_0_z_xyyz_xxyyyyyy[k];

                g_0_z_xxyyz_xyyyyyz[k] = -g_0_z_xyyz_xyyyyyz[k] * ab_x + g_0_z_xyyz_xxyyyyyz[k];

                g_0_z_xxyyz_xyyyyzz[k] = -g_0_z_xyyz_xyyyyzz[k] * ab_x + g_0_z_xyyz_xxyyyyzz[k];

                g_0_z_xxyyz_xyyyzzz[k] = -g_0_z_xyyz_xyyyzzz[k] * ab_x + g_0_z_xyyz_xxyyyzzz[k];

                g_0_z_xxyyz_xyyzzzz[k] = -g_0_z_xyyz_xyyzzzz[k] * ab_x + g_0_z_xyyz_xxyyzzzz[k];

                g_0_z_xxyyz_xyzzzzz[k] = -g_0_z_xyyz_xyzzzzz[k] * ab_x + g_0_z_xyyz_xxyzzzzz[k];

                g_0_z_xxyyz_xzzzzzz[k] = -g_0_z_xyyz_xzzzzzz[k] * ab_x + g_0_z_xyyz_xxzzzzzz[k];

                g_0_z_xxyyz_yyyyyyy[k] = -g_0_z_xyyz_yyyyyyy[k] * ab_x + g_0_z_xyyz_xyyyyyyy[k];

                g_0_z_xxyyz_yyyyyyz[k] = -g_0_z_xyyz_yyyyyyz[k] * ab_x + g_0_z_xyyz_xyyyyyyz[k];

                g_0_z_xxyyz_yyyyyzz[k] = -g_0_z_xyyz_yyyyyzz[k] * ab_x + g_0_z_xyyz_xyyyyyzz[k];

                g_0_z_xxyyz_yyyyzzz[k] = -g_0_z_xyyz_yyyyzzz[k] * ab_x + g_0_z_xyyz_xyyyyzzz[k];

                g_0_z_xxyyz_yyyzzzz[k] = -g_0_z_xyyz_yyyzzzz[k] * ab_x + g_0_z_xyyz_xyyyzzzz[k];

                g_0_z_xxyyz_yyzzzzz[k] = -g_0_z_xyyz_yyzzzzz[k] * ab_x + g_0_z_xyyz_xyyzzzzz[k];

                g_0_z_xxyyz_yzzzzzz[k] = -g_0_z_xyyz_yzzzzzz[k] * ab_x + g_0_z_xyyz_xyzzzzzz[k];

                g_0_z_xxyyz_zzzzzzz[k] = -g_0_z_xyyz_zzzzzzz[k] * ab_x + g_0_z_xyyz_xzzzzzzz[k];
            }

            /// Set up 1800-1836 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1800 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1801 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1802 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1803 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1804 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1805 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1806 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1807 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1808 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1809 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1810 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1811 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1812 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1813 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1814 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1815 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1816 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1817 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1818 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1819 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1820 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1821 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1822 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1823 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1824 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1825 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1826 * ccomps * dcomps);

            auto g_0_z_xxyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1827 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1828 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1829 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1830 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1831 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1832 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1833 * ccomps * dcomps);

            auto g_0_z_xxyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1834 * ccomps * dcomps);

            auto g_0_z_xxyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1835 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzz_xxxxxxx, g_0_z_xxyzz_xxxxxxy, g_0_z_xxyzz_xxxxxxz, g_0_z_xxyzz_xxxxxyy, g_0_z_xxyzz_xxxxxyz, g_0_z_xxyzz_xxxxxzz, g_0_z_xxyzz_xxxxyyy, g_0_z_xxyzz_xxxxyyz, g_0_z_xxyzz_xxxxyzz, g_0_z_xxyzz_xxxxzzz, g_0_z_xxyzz_xxxyyyy, g_0_z_xxyzz_xxxyyyz, g_0_z_xxyzz_xxxyyzz, g_0_z_xxyzz_xxxyzzz, g_0_z_xxyzz_xxxzzzz, g_0_z_xxyzz_xxyyyyy, g_0_z_xxyzz_xxyyyyz, g_0_z_xxyzz_xxyyyzz, g_0_z_xxyzz_xxyyzzz, g_0_z_xxyzz_xxyzzzz, g_0_z_xxyzz_xxzzzzz, g_0_z_xxyzz_xyyyyyy, g_0_z_xxyzz_xyyyyyz, g_0_z_xxyzz_xyyyyzz, g_0_z_xxyzz_xyyyzzz, g_0_z_xxyzz_xyyzzzz, g_0_z_xxyzz_xyzzzzz, g_0_z_xxyzz_xzzzzzz, g_0_z_xxyzz_yyyyyyy, g_0_z_xxyzz_yyyyyyz, g_0_z_xxyzz_yyyyyzz, g_0_z_xxyzz_yyyyzzz, g_0_z_xxyzz_yyyzzzz, g_0_z_xxyzz_yyzzzzz, g_0_z_xxyzz_yzzzzzz, g_0_z_xxyzz_zzzzzzz, g_0_z_xyzz_xxxxxxx, g_0_z_xyzz_xxxxxxxx, g_0_z_xyzz_xxxxxxxy, g_0_z_xyzz_xxxxxxxz, g_0_z_xyzz_xxxxxxy, g_0_z_xyzz_xxxxxxyy, g_0_z_xyzz_xxxxxxyz, g_0_z_xyzz_xxxxxxz, g_0_z_xyzz_xxxxxxzz, g_0_z_xyzz_xxxxxyy, g_0_z_xyzz_xxxxxyyy, g_0_z_xyzz_xxxxxyyz, g_0_z_xyzz_xxxxxyz, g_0_z_xyzz_xxxxxyzz, g_0_z_xyzz_xxxxxzz, g_0_z_xyzz_xxxxxzzz, g_0_z_xyzz_xxxxyyy, g_0_z_xyzz_xxxxyyyy, g_0_z_xyzz_xxxxyyyz, g_0_z_xyzz_xxxxyyz, g_0_z_xyzz_xxxxyyzz, g_0_z_xyzz_xxxxyzz, g_0_z_xyzz_xxxxyzzz, g_0_z_xyzz_xxxxzzz, g_0_z_xyzz_xxxxzzzz, g_0_z_xyzz_xxxyyyy, g_0_z_xyzz_xxxyyyyy, g_0_z_xyzz_xxxyyyyz, g_0_z_xyzz_xxxyyyz, g_0_z_xyzz_xxxyyyzz, g_0_z_xyzz_xxxyyzz, g_0_z_xyzz_xxxyyzzz, g_0_z_xyzz_xxxyzzz, g_0_z_xyzz_xxxyzzzz, g_0_z_xyzz_xxxzzzz, g_0_z_xyzz_xxxzzzzz, g_0_z_xyzz_xxyyyyy, g_0_z_xyzz_xxyyyyyy, g_0_z_xyzz_xxyyyyyz, g_0_z_xyzz_xxyyyyz, g_0_z_xyzz_xxyyyyzz, g_0_z_xyzz_xxyyyzz, g_0_z_xyzz_xxyyyzzz, g_0_z_xyzz_xxyyzzz, g_0_z_xyzz_xxyyzzzz, g_0_z_xyzz_xxyzzzz, g_0_z_xyzz_xxyzzzzz, g_0_z_xyzz_xxzzzzz, g_0_z_xyzz_xxzzzzzz, g_0_z_xyzz_xyyyyyy, g_0_z_xyzz_xyyyyyyy, g_0_z_xyzz_xyyyyyyz, g_0_z_xyzz_xyyyyyz, g_0_z_xyzz_xyyyyyzz, g_0_z_xyzz_xyyyyzz, g_0_z_xyzz_xyyyyzzz, g_0_z_xyzz_xyyyzzz, g_0_z_xyzz_xyyyzzzz, g_0_z_xyzz_xyyzzzz, g_0_z_xyzz_xyyzzzzz, g_0_z_xyzz_xyzzzzz, g_0_z_xyzz_xyzzzzzz, g_0_z_xyzz_xzzzzzz, g_0_z_xyzz_xzzzzzzz, g_0_z_xyzz_yyyyyyy, g_0_z_xyzz_yyyyyyz, g_0_z_xyzz_yyyyyzz, g_0_z_xyzz_yyyyzzz, g_0_z_xyzz_yyyzzzz, g_0_z_xyzz_yyzzzzz, g_0_z_xyzz_yzzzzzz, g_0_z_xyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzz_xxxxxxx[k] = -g_0_z_xyzz_xxxxxxx[k] * ab_x + g_0_z_xyzz_xxxxxxxx[k];

                g_0_z_xxyzz_xxxxxxy[k] = -g_0_z_xyzz_xxxxxxy[k] * ab_x + g_0_z_xyzz_xxxxxxxy[k];

                g_0_z_xxyzz_xxxxxxz[k] = -g_0_z_xyzz_xxxxxxz[k] * ab_x + g_0_z_xyzz_xxxxxxxz[k];

                g_0_z_xxyzz_xxxxxyy[k] = -g_0_z_xyzz_xxxxxyy[k] * ab_x + g_0_z_xyzz_xxxxxxyy[k];

                g_0_z_xxyzz_xxxxxyz[k] = -g_0_z_xyzz_xxxxxyz[k] * ab_x + g_0_z_xyzz_xxxxxxyz[k];

                g_0_z_xxyzz_xxxxxzz[k] = -g_0_z_xyzz_xxxxxzz[k] * ab_x + g_0_z_xyzz_xxxxxxzz[k];

                g_0_z_xxyzz_xxxxyyy[k] = -g_0_z_xyzz_xxxxyyy[k] * ab_x + g_0_z_xyzz_xxxxxyyy[k];

                g_0_z_xxyzz_xxxxyyz[k] = -g_0_z_xyzz_xxxxyyz[k] * ab_x + g_0_z_xyzz_xxxxxyyz[k];

                g_0_z_xxyzz_xxxxyzz[k] = -g_0_z_xyzz_xxxxyzz[k] * ab_x + g_0_z_xyzz_xxxxxyzz[k];

                g_0_z_xxyzz_xxxxzzz[k] = -g_0_z_xyzz_xxxxzzz[k] * ab_x + g_0_z_xyzz_xxxxxzzz[k];

                g_0_z_xxyzz_xxxyyyy[k] = -g_0_z_xyzz_xxxyyyy[k] * ab_x + g_0_z_xyzz_xxxxyyyy[k];

                g_0_z_xxyzz_xxxyyyz[k] = -g_0_z_xyzz_xxxyyyz[k] * ab_x + g_0_z_xyzz_xxxxyyyz[k];

                g_0_z_xxyzz_xxxyyzz[k] = -g_0_z_xyzz_xxxyyzz[k] * ab_x + g_0_z_xyzz_xxxxyyzz[k];

                g_0_z_xxyzz_xxxyzzz[k] = -g_0_z_xyzz_xxxyzzz[k] * ab_x + g_0_z_xyzz_xxxxyzzz[k];

                g_0_z_xxyzz_xxxzzzz[k] = -g_0_z_xyzz_xxxzzzz[k] * ab_x + g_0_z_xyzz_xxxxzzzz[k];

                g_0_z_xxyzz_xxyyyyy[k] = -g_0_z_xyzz_xxyyyyy[k] * ab_x + g_0_z_xyzz_xxxyyyyy[k];

                g_0_z_xxyzz_xxyyyyz[k] = -g_0_z_xyzz_xxyyyyz[k] * ab_x + g_0_z_xyzz_xxxyyyyz[k];

                g_0_z_xxyzz_xxyyyzz[k] = -g_0_z_xyzz_xxyyyzz[k] * ab_x + g_0_z_xyzz_xxxyyyzz[k];

                g_0_z_xxyzz_xxyyzzz[k] = -g_0_z_xyzz_xxyyzzz[k] * ab_x + g_0_z_xyzz_xxxyyzzz[k];

                g_0_z_xxyzz_xxyzzzz[k] = -g_0_z_xyzz_xxyzzzz[k] * ab_x + g_0_z_xyzz_xxxyzzzz[k];

                g_0_z_xxyzz_xxzzzzz[k] = -g_0_z_xyzz_xxzzzzz[k] * ab_x + g_0_z_xyzz_xxxzzzzz[k];

                g_0_z_xxyzz_xyyyyyy[k] = -g_0_z_xyzz_xyyyyyy[k] * ab_x + g_0_z_xyzz_xxyyyyyy[k];

                g_0_z_xxyzz_xyyyyyz[k] = -g_0_z_xyzz_xyyyyyz[k] * ab_x + g_0_z_xyzz_xxyyyyyz[k];

                g_0_z_xxyzz_xyyyyzz[k] = -g_0_z_xyzz_xyyyyzz[k] * ab_x + g_0_z_xyzz_xxyyyyzz[k];

                g_0_z_xxyzz_xyyyzzz[k] = -g_0_z_xyzz_xyyyzzz[k] * ab_x + g_0_z_xyzz_xxyyyzzz[k];

                g_0_z_xxyzz_xyyzzzz[k] = -g_0_z_xyzz_xyyzzzz[k] * ab_x + g_0_z_xyzz_xxyyzzzz[k];

                g_0_z_xxyzz_xyzzzzz[k] = -g_0_z_xyzz_xyzzzzz[k] * ab_x + g_0_z_xyzz_xxyzzzzz[k];

                g_0_z_xxyzz_xzzzzzz[k] = -g_0_z_xyzz_xzzzzzz[k] * ab_x + g_0_z_xyzz_xxzzzzzz[k];

                g_0_z_xxyzz_yyyyyyy[k] = -g_0_z_xyzz_yyyyyyy[k] * ab_x + g_0_z_xyzz_xyyyyyyy[k];

                g_0_z_xxyzz_yyyyyyz[k] = -g_0_z_xyzz_yyyyyyz[k] * ab_x + g_0_z_xyzz_xyyyyyyz[k];

                g_0_z_xxyzz_yyyyyzz[k] = -g_0_z_xyzz_yyyyyzz[k] * ab_x + g_0_z_xyzz_xyyyyyzz[k];

                g_0_z_xxyzz_yyyyzzz[k] = -g_0_z_xyzz_yyyyzzz[k] * ab_x + g_0_z_xyzz_xyyyyzzz[k];

                g_0_z_xxyzz_yyyzzzz[k] = -g_0_z_xyzz_yyyzzzz[k] * ab_x + g_0_z_xyzz_xyyyzzzz[k];

                g_0_z_xxyzz_yyzzzzz[k] = -g_0_z_xyzz_yyzzzzz[k] * ab_x + g_0_z_xyzz_xyyzzzzz[k];

                g_0_z_xxyzz_yzzzzzz[k] = -g_0_z_xyzz_yzzzzzz[k] * ab_x + g_0_z_xyzz_xyzzzzzz[k];

                g_0_z_xxyzz_zzzzzzz[k] = -g_0_z_xyzz_zzzzzzz[k] * ab_x + g_0_z_xyzz_xzzzzzzz[k];
            }

            /// Set up 1836-1872 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1836 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1837 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1838 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1839 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1840 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1841 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1842 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1843 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1844 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1845 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1846 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1847 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1848 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1849 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1850 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1851 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1852 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1853 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1854 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1855 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1856 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1857 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1858 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1859 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1860 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1861 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1862 * ccomps * dcomps);

            auto g_0_z_xxzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1863 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1864 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1865 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1866 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1867 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1868 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1869 * ccomps * dcomps);

            auto g_0_z_xxzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1870 * ccomps * dcomps);

            auto g_0_z_xxzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1871 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzz_xxxxxxx, g_0_z_xxzzz_xxxxxxy, g_0_z_xxzzz_xxxxxxz, g_0_z_xxzzz_xxxxxyy, g_0_z_xxzzz_xxxxxyz, g_0_z_xxzzz_xxxxxzz, g_0_z_xxzzz_xxxxyyy, g_0_z_xxzzz_xxxxyyz, g_0_z_xxzzz_xxxxyzz, g_0_z_xxzzz_xxxxzzz, g_0_z_xxzzz_xxxyyyy, g_0_z_xxzzz_xxxyyyz, g_0_z_xxzzz_xxxyyzz, g_0_z_xxzzz_xxxyzzz, g_0_z_xxzzz_xxxzzzz, g_0_z_xxzzz_xxyyyyy, g_0_z_xxzzz_xxyyyyz, g_0_z_xxzzz_xxyyyzz, g_0_z_xxzzz_xxyyzzz, g_0_z_xxzzz_xxyzzzz, g_0_z_xxzzz_xxzzzzz, g_0_z_xxzzz_xyyyyyy, g_0_z_xxzzz_xyyyyyz, g_0_z_xxzzz_xyyyyzz, g_0_z_xxzzz_xyyyzzz, g_0_z_xxzzz_xyyzzzz, g_0_z_xxzzz_xyzzzzz, g_0_z_xxzzz_xzzzzzz, g_0_z_xxzzz_yyyyyyy, g_0_z_xxzzz_yyyyyyz, g_0_z_xxzzz_yyyyyzz, g_0_z_xxzzz_yyyyzzz, g_0_z_xxzzz_yyyzzzz, g_0_z_xxzzz_yyzzzzz, g_0_z_xxzzz_yzzzzzz, g_0_z_xxzzz_zzzzzzz, g_0_z_xzzz_xxxxxxx, g_0_z_xzzz_xxxxxxxx, g_0_z_xzzz_xxxxxxxy, g_0_z_xzzz_xxxxxxxz, g_0_z_xzzz_xxxxxxy, g_0_z_xzzz_xxxxxxyy, g_0_z_xzzz_xxxxxxyz, g_0_z_xzzz_xxxxxxz, g_0_z_xzzz_xxxxxxzz, g_0_z_xzzz_xxxxxyy, g_0_z_xzzz_xxxxxyyy, g_0_z_xzzz_xxxxxyyz, g_0_z_xzzz_xxxxxyz, g_0_z_xzzz_xxxxxyzz, g_0_z_xzzz_xxxxxzz, g_0_z_xzzz_xxxxxzzz, g_0_z_xzzz_xxxxyyy, g_0_z_xzzz_xxxxyyyy, g_0_z_xzzz_xxxxyyyz, g_0_z_xzzz_xxxxyyz, g_0_z_xzzz_xxxxyyzz, g_0_z_xzzz_xxxxyzz, g_0_z_xzzz_xxxxyzzz, g_0_z_xzzz_xxxxzzz, g_0_z_xzzz_xxxxzzzz, g_0_z_xzzz_xxxyyyy, g_0_z_xzzz_xxxyyyyy, g_0_z_xzzz_xxxyyyyz, g_0_z_xzzz_xxxyyyz, g_0_z_xzzz_xxxyyyzz, g_0_z_xzzz_xxxyyzz, g_0_z_xzzz_xxxyyzzz, g_0_z_xzzz_xxxyzzz, g_0_z_xzzz_xxxyzzzz, g_0_z_xzzz_xxxzzzz, g_0_z_xzzz_xxxzzzzz, g_0_z_xzzz_xxyyyyy, g_0_z_xzzz_xxyyyyyy, g_0_z_xzzz_xxyyyyyz, g_0_z_xzzz_xxyyyyz, g_0_z_xzzz_xxyyyyzz, g_0_z_xzzz_xxyyyzz, g_0_z_xzzz_xxyyyzzz, g_0_z_xzzz_xxyyzzz, g_0_z_xzzz_xxyyzzzz, g_0_z_xzzz_xxyzzzz, g_0_z_xzzz_xxyzzzzz, g_0_z_xzzz_xxzzzzz, g_0_z_xzzz_xxzzzzzz, g_0_z_xzzz_xyyyyyy, g_0_z_xzzz_xyyyyyyy, g_0_z_xzzz_xyyyyyyz, g_0_z_xzzz_xyyyyyz, g_0_z_xzzz_xyyyyyzz, g_0_z_xzzz_xyyyyzz, g_0_z_xzzz_xyyyyzzz, g_0_z_xzzz_xyyyzzz, g_0_z_xzzz_xyyyzzzz, g_0_z_xzzz_xyyzzzz, g_0_z_xzzz_xyyzzzzz, g_0_z_xzzz_xyzzzzz, g_0_z_xzzz_xyzzzzzz, g_0_z_xzzz_xzzzzzz, g_0_z_xzzz_xzzzzzzz, g_0_z_xzzz_yyyyyyy, g_0_z_xzzz_yyyyyyz, g_0_z_xzzz_yyyyyzz, g_0_z_xzzz_yyyyzzz, g_0_z_xzzz_yyyzzzz, g_0_z_xzzz_yyzzzzz, g_0_z_xzzz_yzzzzzz, g_0_z_xzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzz_xxxxxxx[k] = -g_0_z_xzzz_xxxxxxx[k] * ab_x + g_0_z_xzzz_xxxxxxxx[k];

                g_0_z_xxzzz_xxxxxxy[k] = -g_0_z_xzzz_xxxxxxy[k] * ab_x + g_0_z_xzzz_xxxxxxxy[k];

                g_0_z_xxzzz_xxxxxxz[k] = -g_0_z_xzzz_xxxxxxz[k] * ab_x + g_0_z_xzzz_xxxxxxxz[k];

                g_0_z_xxzzz_xxxxxyy[k] = -g_0_z_xzzz_xxxxxyy[k] * ab_x + g_0_z_xzzz_xxxxxxyy[k];

                g_0_z_xxzzz_xxxxxyz[k] = -g_0_z_xzzz_xxxxxyz[k] * ab_x + g_0_z_xzzz_xxxxxxyz[k];

                g_0_z_xxzzz_xxxxxzz[k] = -g_0_z_xzzz_xxxxxzz[k] * ab_x + g_0_z_xzzz_xxxxxxzz[k];

                g_0_z_xxzzz_xxxxyyy[k] = -g_0_z_xzzz_xxxxyyy[k] * ab_x + g_0_z_xzzz_xxxxxyyy[k];

                g_0_z_xxzzz_xxxxyyz[k] = -g_0_z_xzzz_xxxxyyz[k] * ab_x + g_0_z_xzzz_xxxxxyyz[k];

                g_0_z_xxzzz_xxxxyzz[k] = -g_0_z_xzzz_xxxxyzz[k] * ab_x + g_0_z_xzzz_xxxxxyzz[k];

                g_0_z_xxzzz_xxxxzzz[k] = -g_0_z_xzzz_xxxxzzz[k] * ab_x + g_0_z_xzzz_xxxxxzzz[k];

                g_0_z_xxzzz_xxxyyyy[k] = -g_0_z_xzzz_xxxyyyy[k] * ab_x + g_0_z_xzzz_xxxxyyyy[k];

                g_0_z_xxzzz_xxxyyyz[k] = -g_0_z_xzzz_xxxyyyz[k] * ab_x + g_0_z_xzzz_xxxxyyyz[k];

                g_0_z_xxzzz_xxxyyzz[k] = -g_0_z_xzzz_xxxyyzz[k] * ab_x + g_0_z_xzzz_xxxxyyzz[k];

                g_0_z_xxzzz_xxxyzzz[k] = -g_0_z_xzzz_xxxyzzz[k] * ab_x + g_0_z_xzzz_xxxxyzzz[k];

                g_0_z_xxzzz_xxxzzzz[k] = -g_0_z_xzzz_xxxzzzz[k] * ab_x + g_0_z_xzzz_xxxxzzzz[k];

                g_0_z_xxzzz_xxyyyyy[k] = -g_0_z_xzzz_xxyyyyy[k] * ab_x + g_0_z_xzzz_xxxyyyyy[k];

                g_0_z_xxzzz_xxyyyyz[k] = -g_0_z_xzzz_xxyyyyz[k] * ab_x + g_0_z_xzzz_xxxyyyyz[k];

                g_0_z_xxzzz_xxyyyzz[k] = -g_0_z_xzzz_xxyyyzz[k] * ab_x + g_0_z_xzzz_xxxyyyzz[k];

                g_0_z_xxzzz_xxyyzzz[k] = -g_0_z_xzzz_xxyyzzz[k] * ab_x + g_0_z_xzzz_xxxyyzzz[k];

                g_0_z_xxzzz_xxyzzzz[k] = -g_0_z_xzzz_xxyzzzz[k] * ab_x + g_0_z_xzzz_xxxyzzzz[k];

                g_0_z_xxzzz_xxzzzzz[k] = -g_0_z_xzzz_xxzzzzz[k] * ab_x + g_0_z_xzzz_xxxzzzzz[k];

                g_0_z_xxzzz_xyyyyyy[k] = -g_0_z_xzzz_xyyyyyy[k] * ab_x + g_0_z_xzzz_xxyyyyyy[k];

                g_0_z_xxzzz_xyyyyyz[k] = -g_0_z_xzzz_xyyyyyz[k] * ab_x + g_0_z_xzzz_xxyyyyyz[k];

                g_0_z_xxzzz_xyyyyzz[k] = -g_0_z_xzzz_xyyyyzz[k] * ab_x + g_0_z_xzzz_xxyyyyzz[k];

                g_0_z_xxzzz_xyyyzzz[k] = -g_0_z_xzzz_xyyyzzz[k] * ab_x + g_0_z_xzzz_xxyyyzzz[k];

                g_0_z_xxzzz_xyyzzzz[k] = -g_0_z_xzzz_xyyzzzz[k] * ab_x + g_0_z_xzzz_xxyyzzzz[k];

                g_0_z_xxzzz_xyzzzzz[k] = -g_0_z_xzzz_xyzzzzz[k] * ab_x + g_0_z_xzzz_xxyzzzzz[k];

                g_0_z_xxzzz_xzzzzzz[k] = -g_0_z_xzzz_xzzzzzz[k] * ab_x + g_0_z_xzzz_xxzzzzzz[k];

                g_0_z_xxzzz_yyyyyyy[k] = -g_0_z_xzzz_yyyyyyy[k] * ab_x + g_0_z_xzzz_xyyyyyyy[k];

                g_0_z_xxzzz_yyyyyyz[k] = -g_0_z_xzzz_yyyyyyz[k] * ab_x + g_0_z_xzzz_xyyyyyyz[k];

                g_0_z_xxzzz_yyyyyzz[k] = -g_0_z_xzzz_yyyyyzz[k] * ab_x + g_0_z_xzzz_xyyyyyzz[k];

                g_0_z_xxzzz_yyyyzzz[k] = -g_0_z_xzzz_yyyyzzz[k] * ab_x + g_0_z_xzzz_xyyyyzzz[k];

                g_0_z_xxzzz_yyyzzzz[k] = -g_0_z_xzzz_yyyzzzz[k] * ab_x + g_0_z_xzzz_xyyyzzzz[k];

                g_0_z_xxzzz_yyzzzzz[k] = -g_0_z_xzzz_yyzzzzz[k] * ab_x + g_0_z_xzzz_xyyzzzzz[k];

                g_0_z_xxzzz_yzzzzzz[k] = -g_0_z_xzzz_yzzzzzz[k] * ab_x + g_0_z_xzzz_xyzzzzzz[k];

                g_0_z_xxzzz_zzzzzzz[k] = -g_0_z_xzzz_zzzzzzz[k] * ab_x + g_0_z_xzzz_xzzzzzzz[k];
            }

            /// Set up 1872-1908 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 1872 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 1873 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 1874 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 1875 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 1876 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 1877 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 1878 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 1879 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 1880 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 1881 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 1882 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 1883 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 1884 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 1885 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 1886 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 1887 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 1888 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 1889 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 1890 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 1891 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 1892 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 1893 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 1894 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 1895 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 1896 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 1897 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 1898 * ccomps * dcomps);

            auto g_0_z_xyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 1899 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 1900 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 1901 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 1902 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 1903 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 1904 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 1905 * ccomps * dcomps);

            auto g_0_z_xyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 1906 * ccomps * dcomps);

            auto g_0_z_xyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 1907 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyy_xxxxxxx, g_0_z_xyyyy_xxxxxxy, g_0_z_xyyyy_xxxxxxz, g_0_z_xyyyy_xxxxxyy, g_0_z_xyyyy_xxxxxyz, g_0_z_xyyyy_xxxxxzz, g_0_z_xyyyy_xxxxyyy, g_0_z_xyyyy_xxxxyyz, g_0_z_xyyyy_xxxxyzz, g_0_z_xyyyy_xxxxzzz, g_0_z_xyyyy_xxxyyyy, g_0_z_xyyyy_xxxyyyz, g_0_z_xyyyy_xxxyyzz, g_0_z_xyyyy_xxxyzzz, g_0_z_xyyyy_xxxzzzz, g_0_z_xyyyy_xxyyyyy, g_0_z_xyyyy_xxyyyyz, g_0_z_xyyyy_xxyyyzz, g_0_z_xyyyy_xxyyzzz, g_0_z_xyyyy_xxyzzzz, g_0_z_xyyyy_xxzzzzz, g_0_z_xyyyy_xyyyyyy, g_0_z_xyyyy_xyyyyyz, g_0_z_xyyyy_xyyyyzz, g_0_z_xyyyy_xyyyzzz, g_0_z_xyyyy_xyyzzzz, g_0_z_xyyyy_xyzzzzz, g_0_z_xyyyy_xzzzzzz, g_0_z_xyyyy_yyyyyyy, g_0_z_xyyyy_yyyyyyz, g_0_z_xyyyy_yyyyyzz, g_0_z_xyyyy_yyyyzzz, g_0_z_xyyyy_yyyzzzz, g_0_z_xyyyy_yyzzzzz, g_0_z_xyyyy_yzzzzzz, g_0_z_xyyyy_zzzzzzz, g_0_z_yyyy_xxxxxxx, g_0_z_yyyy_xxxxxxxx, g_0_z_yyyy_xxxxxxxy, g_0_z_yyyy_xxxxxxxz, g_0_z_yyyy_xxxxxxy, g_0_z_yyyy_xxxxxxyy, g_0_z_yyyy_xxxxxxyz, g_0_z_yyyy_xxxxxxz, g_0_z_yyyy_xxxxxxzz, g_0_z_yyyy_xxxxxyy, g_0_z_yyyy_xxxxxyyy, g_0_z_yyyy_xxxxxyyz, g_0_z_yyyy_xxxxxyz, g_0_z_yyyy_xxxxxyzz, g_0_z_yyyy_xxxxxzz, g_0_z_yyyy_xxxxxzzz, g_0_z_yyyy_xxxxyyy, g_0_z_yyyy_xxxxyyyy, g_0_z_yyyy_xxxxyyyz, g_0_z_yyyy_xxxxyyz, g_0_z_yyyy_xxxxyyzz, g_0_z_yyyy_xxxxyzz, g_0_z_yyyy_xxxxyzzz, g_0_z_yyyy_xxxxzzz, g_0_z_yyyy_xxxxzzzz, g_0_z_yyyy_xxxyyyy, g_0_z_yyyy_xxxyyyyy, g_0_z_yyyy_xxxyyyyz, g_0_z_yyyy_xxxyyyz, g_0_z_yyyy_xxxyyyzz, g_0_z_yyyy_xxxyyzz, g_0_z_yyyy_xxxyyzzz, g_0_z_yyyy_xxxyzzz, g_0_z_yyyy_xxxyzzzz, g_0_z_yyyy_xxxzzzz, g_0_z_yyyy_xxxzzzzz, g_0_z_yyyy_xxyyyyy, g_0_z_yyyy_xxyyyyyy, g_0_z_yyyy_xxyyyyyz, g_0_z_yyyy_xxyyyyz, g_0_z_yyyy_xxyyyyzz, g_0_z_yyyy_xxyyyzz, g_0_z_yyyy_xxyyyzzz, g_0_z_yyyy_xxyyzzz, g_0_z_yyyy_xxyyzzzz, g_0_z_yyyy_xxyzzzz, g_0_z_yyyy_xxyzzzzz, g_0_z_yyyy_xxzzzzz, g_0_z_yyyy_xxzzzzzz, g_0_z_yyyy_xyyyyyy, g_0_z_yyyy_xyyyyyyy, g_0_z_yyyy_xyyyyyyz, g_0_z_yyyy_xyyyyyz, g_0_z_yyyy_xyyyyyzz, g_0_z_yyyy_xyyyyzz, g_0_z_yyyy_xyyyyzzz, g_0_z_yyyy_xyyyzzz, g_0_z_yyyy_xyyyzzzz, g_0_z_yyyy_xyyzzzz, g_0_z_yyyy_xyyzzzzz, g_0_z_yyyy_xyzzzzz, g_0_z_yyyy_xyzzzzzz, g_0_z_yyyy_xzzzzzz, g_0_z_yyyy_xzzzzzzz, g_0_z_yyyy_yyyyyyy, g_0_z_yyyy_yyyyyyz, g_0_z_yyyy_yyyyyzz, g_0_z_yyyy_yyyyzzz, g_0_z_yyyy_yyyzzzz, g_0_z_yyyy_yyzzzzz, g_0_z_yyyy_yzzzzzz, g_0_z_yyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyy_xxxxxxx[k] = -g_0_z_yyyy_xxxxxxx[k] * ab_x + g_0_z_yyyy_xxxxxxxx[k];

                g_0_z_xyyyy_xxxxxxy[k] = -g_0_z_yyyy_xxxxxxy[k] * ab_x + g_0_z_yyyy_xxxxxxxy[k];

                g_0_z_xyyyy_xxxxxxz[k] = -g_0_z_yyyy_xxxxxxz[k] * ab_x + g_0_z_yyyy_xxxxxxxz[k];

                g_0_z_xyyyy_xxxxxyy[k] = -g_0_z_yyyy_xxxxxyy[k] * ab_x + g_0_z_yyyy_xxxxxxyy[k];

                g_0_z_xyyyy_xxxxxyz[k] = -g_0_z_yyyy_xxxxxyz[k] * ab_x + g_0_z_yyyy_xxxxxxyz[k];

                g_0_z_xyyyy_xxxxxzz[k] = -g_0_z_yyyy_xxxxxzz[k] * ab_x + g_0_z_yyyy_xxxxxxzz[k];

                g_0_z_xyyyy_xxxxyyy[k] = -g_0_z_yyyy_xxxxyyy[k] * ab_x + g_0_z_yyyy_xxxxxyyy[k];

                g_0_z_xyyyy_xxxxyyz[k] = -g_0_z_yyyy_xxxxyyz[k] * ab_x + g_0_z_yyyy_xxxxxyyz[k];

                g_0_z_xyyyy_xxxxyzz[k] = -g_0_z_yyyy_xxxxyzz[k] * ab_x + g_0_z_yyyy_xxxxxyzz[k];

                g_0_z_xyyyy_xxxxzzz[k] = -g_0_z_yyyy_xxxxzzz[k] * ab_x + g_0_z_yyyy_xxxxxzzz[k];

                g_0_z_xyyyy_xxxyyyy[k] = -g_0_z_yyyy_xxxyyyy[k] * ab_x + g_0_z_yyyy_xxxxyyyy[k];

                g_0_z_xyyyy_xxxyyyz[k] = -g_0_z_yyyy_xxxyyyz[k] * ab_x + g_0_z_yyyy_xxxxyyyz[k];

                g_0_z_xyyyy_xxxyyzz[k] = -g_0_z_yyyy_xxxyyzz[k] * ab_x + g_0_z_yyyy_xxxxyyzz[k];

                g_0_z_xyyyy_xxxyzzz[k] = -g_0_z_yyyy_xxxyzzz[k] * ab_x + g_0_z_yyyy_xxxxyzzz[k];

                g_0_z_xyyyy_xxxzzzz[k] = -g_0_z_yyyy_xxxzzzz[k] * ab_x + g_0_z_yyyy_xxxxzzzz[k];

                g_0_z_xyyyy_xxyyyyy[k] = -g_0_z_yyyy_xxyyyyy[k] * ab_x + g_0_z_yyyy_xxxyyyyy[k];

                g_0_z_xyyyy_xxyyyyz[k] = -g_0_z_yyyy_xxyyyyz[k] * ab_x + g_0_z_yyyy_xxxyyyyz[k];

                g_0_z_xyyyy_xxyyyzz[k] = -g_0_z_yyyy_xxyyyzz[k] * ab_x + g_0_z_yyyy_xxxyyyzz[k];

                g_0_z_xyyyy_xxyyzzz[k] = -g_0_z_yyyy_xxyyzzz[k] * ab_x + g_0_z_yyyy_xxxyyzzz[k];

                g_0_z_xyyyy_xxyzzzz[k] = -g_0_z_yyyy_xxyzzzz[k] * ab_x + g_0_z_yyyy_xxxyzzzz[k];

                g_0_z_xyyyy_xxzzzzz[k] = -g_0_z_yyyy_xxzzzzz[k] * ab_x + g_0_z_yyyy_xxxzzzzz[k];

                g_0_z_xyyyy_xyyyyyy[k] = -g_0_z_yyyy_xyyyyyy[k] * ab_x + g_0_z_yyyy_xxyyyyyy[k];

                g_0_z_xyyyy_xyyyyyz[k] = -g_0_z_yyyy_xyyyyyz[k] * ab_x + g_0_z_yyyy_xxyyyyyz[k];

                g_0_z_xyyyy_xyyyyzz[k] = -g_0_z_yyyy_xyyyyzz[k] * ab_x + g_0_z_yyyy_xxyyyyzz[k];

                g_0_z_xyyyy_xyyyzzz[k] = -g_0_z_yyyy_xyyyzzz[k] * ab_x + g_0_z_yyyy_xxyyyzzz[k];

                g_0_z_xyyyy_xyyzzzz[k] = -g_0_z_yyyy_xyyzzzz[k] * ab_x + g_0_z_yyyy_xxyyzzzz[k];

                g_0_z_xyyyy_xyzzzzz[k] = -g_0_z_yyyy_xyzzzzz[k] * ab_x + g_0_z_yyyy_xxyzzzzz[k];

                g_0_z_xyyyy_xzzzzzz[k] = -g_0_z_yyyy_xzzzzzz[k] * ab_x + g_0_z_yyyy_xxzzzzzz[k];

                g_0_z_xyyyy_yyyyyyy[k] = -g_0_z_yyyy_yyyyyyy[k] * ab_x + g_0_z_yyyy_xyyyyyyy[k];

                g_0_z_xyyyy_yyyyyyz[k] = -g_0_z_yyyy_yyyyyyz[k] * ab_x + g_0_z_yyyy_xyyyyyyz[k];

                g_0_z_xyyyy_yyyyyzz[k] = -g_0_z_yyyy_yyyyyzz[k] * ab_x + g_0_z_yyyy_xyyyyyzz[k];

                g_0_z_xyyyy_yyyyzzz[k] = -g_0_z_yyyy_yyyyzzz[k] * ab_x + g_0_z_yyyy_xyyyyzzz[k];

                g_0_z_xyyyy_yyyzzzz[k] = -g_0_z_yyyy_yyyzzzz[k] * ab_x + g_0_z_yyyy_xyyyzzzz[k];

                g_0_z_xyyyy_yyzzzzz[k] = -g_0_z_yyyy_yyzzzzz[k] * ab_x + g_0_z_yyyy_xyyzzzzz[k];

                g_0_z_xyyyy_yzzzzzz[k] = -g_0_z_yyyy_yzzzzzz[k] * ab_x + g_0_z_yyyy_xyzzzzzz[k];

                g_0_z_xyyyy_zzzzzzz[k] = -g_0_z_yyyy_zzzzzzz[k] * ab_x + g_0_z_yyyy_xzzzzzzz[k];
            }

            /// Set up 1908-1944 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1908 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1909 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1910 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1911 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1912 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1913 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1914 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1915 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1916 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1917 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1918 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1919 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1920 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1921 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1922 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1923 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1924 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1925 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1926 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1927 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1928 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1929 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1930 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1931 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1932 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1933 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1934 * ccomps * dcomps);

            auto g_0_z_xyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1935 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1936 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1937 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1938 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1939 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1940 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1941 * ccomps * dcomps);

            auto g_0_z_xyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1942 * ccomps * dcomps);

            auto g_0_z_xyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1943 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyz_xxxxxxx, g_0_z_xyyyz_xxxxxxy, g_0_z_xyyyz_xxxxxxz, g_0_z_xyyyz_xxxxxyy, g_0_z_xyyyz_xxxxxyz, g_0_z_xyyyz_xxxxxzz, g_0_z_xyyyz_xxxxyyy, g_0_z_xyyyz_xxxxyyz, g_0_z_xyyyz_xxxxyzz, g_0_z_xyyyz_xxxxzzz, g_0_z_xyyyz_xxxyyyy, g_0_z_xyyyz_xxxyyyz, g_0_z_xyyyz_xxxyyzz, g_0_z_xyyyz_xxxyzzz, g_0_z_xyyyz_xxxzzzz, g_0_z_xyyyz_xxyyyyy, g_0_z_xyyyz_xxyyyyz, g_0_z_xyyyz_xxyyyzz, g_0_z_xyyyz_xxyyzzz, g_0_z_xyyyz_xxyzzzz, g_0_z_xyyyz_xxzzzzz, g_0_z_xyyyz_xyyyyyy, g_0_z_xyyyz_xyyyyyz, g_0_z_xyyyz_xyyyyzz, g_0_z_xyyyz_xyyyzzz, g_0_z_xyyyz_xyyzzzz, g_0_z_xyyyz_xyzzzzz, g_0_z_xyyyz_xzzzzzz, g_0_z_xyyyz_yyyyyyy, g_0_z_xyyyz_yyyyyyz, g_0_z_xyyyz_yyyyyzz, g_0_z_xyyyz_yyyyzzz, g_0_z_xyyyz_yyyzzzz, g_0_z_xyyyz_yyzzzzz, g_0_z_xyyyz_yzzzzzz, g_0_z_xyyyz_zzzzzzz, g_0_z_yyyz_xxxxxxx, g_0_z_yyyz_xxxxxxxx, g_0_z_yyyz_xxxxxxxy, g_0_z_yyyz_xxxxxxxz, g_0_z_yyyz_xxxxxxy, g_0_z_yyyz_xxxxxxyy, g_0_z_yyyz_xxxxxxyz, g_0_z_yyyz_xxxxxxz, g_0_z_yyyz_xxxxxxzz, g_0_z_yyyz_xxxxxyy, g_0_z_yyyz_xxxxxyyy, g_0_z_yyyz_xxxxxyyz, g_0_z_yyyz_xxxxxyz, g_0_z_yyyz_xxxxxyzz, g_0_z_yyyz_xxxxxzz, g_0_z_yyyz_xxxxxzzz, g_0_z_yyyz_xxxxyyy, g_0_z_yyyz_xxxxyyyy, g_0_z_yyyz_xxxxyyyz, g_0_z_yyyz_xxxxyyz, g_0_z_yyyz_xxxxyyzz, g_0_z_yyyz_xxxxyzz, g_0_z_yyyz_xxxxyzzz, g_0_z_yyyz_xxxxzzz, g_0_z_yyyz_xxxxzzzz, g_0_z_yyyz_xxxyyyy, g_0_z_yyyz_xxxyyyyy, g_0_z_yyyz_xxxyyyyz, g_0_z_yyyz_xxxyyyz, g_0_z_yyyz_xxxyyyzz, g_0_z_yyyz_xxxyyzz, g_0_z_yyyz_xxxyyzzz, g_0_z_yyyz_xxxyzzz, g_0_z_yyyz_xxxyzzzz, g_0_z_yyyz_xxxzzzz, g_0_z_yyyz_xxxzzzzz, g_0_z_yyyz_xxyyyyy, g_0_z_yyyz_xxyyyyyy, g_0_z_yyyz_xxyyyyyz, g_0_z_yyyz_xxyyyyz, g_0_z_yyyz_xxyyyyzz, g_0_z_yyyz_xxyyyzz, g_0_z_yyyz_xxyyyzzz, g_0_z_yyyz_xxyyzzz, g_0_z_yyyz_xxyyzzzz, g_0_z_yyyz_xxyzzzz, g_0_z_yyyz_xxyzzzzz, g_0_z_yyyz_xxzzzzz, g_0_z_yyyz_xxzzzzzz, g_0_z_yyyz_xyyyyyy, g_0_z_yyyz_xyyyyyyy, g_0_z_yyyz_xyyyyyyz, g_0_z_yyyz_xyyyyyz, g_0_z_yyyz_xyyyyyzz, g_0_z_yyyz_xyyyyzz, g_0_z_yyyz_xyyyyzzz, g_0_z_yyyz_xyyyzzz, g_0_z_yyyz_xyyyzzzz, g_0_z_yyyz_xyyzzzz, g_0_z_yyyz_xyyzzzzz, g_0_z_yyyz_xyzzzzz, g_0_z_yyyz_xyzzzzzz, g_0_z_yyyz_xzzzzzz, g_0_z_yyyz_xzzzzzzz, g_0_z_yyyz_yyyyyyy, g_0_z_yyyz_yyyyyyz, g_0_z_yyyz_yyyyyzz, g_0_z_yyyz_yyyyzzz, g_0_z_yyyz_yyyzzzz, g_0_z_yyyz_yyzzzzz, g_0_z_yyyz_yzzzzzz, g_0_z_yyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyz_xxxxxxx[k] = -g_0_z_yyyz_xxxxxxx[k] * ab_x + g_0_z_yyyz_xxxxxxxx[k];

                g_0_z_xyyyz_xxxxxxy[k] = -g_0_z_yyyz_xxxxxxy[k] * ab_x + g_0_z_yyyz_xxxxxxxy[k];

                g_0_z_xyyyz_xxxxxxz[k] = -g_0_z_yyyz_xxxxxxz[k] * ab_x + g_0_z_yyyz_xxxxxxxz[k];

                g_0_z_xyyyz_xxxxxyy[k] = -g_0_z_yyyz_xxxxxyy[k] * ab_x + g_0_z_yyyz_xxxxxxyy[k];

                g_0_z_xyyyz_xxxxxyz[k] = -g_0_z_yyyz_xxxxxyz[k] * ab_x + g_0_z_yyyz_xxxxxxyz[k];

                g_0_z_xyyyz_xxxxxzz[k] = -g_0_z_yyyz_xxxxxzz[k] * ab_x + g_0_z_yyyz_xxxxxxzz[k];

                g_0_z_xyyyz_xxxxyyy[k] = -g_0_z_yyyz_xxxxyyy[k] * ab_x + g_0_z_yyyz_xxxxxyyy[k];

                g_0_z_xyyyz_xxxxyyz[k] = -g_0_z_yyyz_xxxxyyz[k] * ab_x + g_0_z_yyyz_xxxxxyyz[k];

                g_0_z_xyyyz_xxxxyzz[k] = -g_0_z_yyyz_xxxxyzz[k] * ab_x + g_0_z_yyyz_xxxxxyzz[k];

                g_0_z_xyyyz_xxxxzzz[k] = -g_0_z_yyyz_xxxxzzz[k] * ab_x + g_0_z_yyyz_xxxxxzzz[k];

                g_0_z_xyyyz_xxxyyyy[k] = -g_0_z_yyyz_xxxyyyy[k] * ab_x + g_0_z_yyyz_xxxxyyyy[k];

                g_0_z_xyyyz_xxxyyyz[k] = -g_0_z_yyyz_xxxyyyz[k] * ab_x + g_0_z_yyyz_xxxxyyyz[k];

                g_0_z_xyyyz_xxxyyzz[k] = -g_0_z_yyyz_xxxyyzz[k] * ab_x + g_0_z_yyyz_xxxxyyzz[k];

                g_0_z_xyyyz_xxxyzzz[k] = -g_0_z_yyyz_xxxyzzz[k] * ab_x + g_0_z_yyyz_xxxxyzzz[k];

                g_0_z_xyyyz_xxxzzzz[k] = -g_0_z_yyyz_xxxzzzz[k] * ab_x + g_0_z_yyyz_xxxxzzzz[k];

                g_0_z_xyyyz_xxyyyyy[k] = -g_0_z_yyyz_xxyyyyy[k] * ab_x + g_0_z_yyyz_xxxyyyyy[k];

                g_0_z_xyyyz_xxyyyyz[k] = -g_0_z_yyyz_xxyyyyz[k] * ab_x + g_0_z_yyyz_xxxyyyyz[k];

                g_0_z_xyyyz_xxyyyzz[k] = -g_0_z_yyyz_xxyyyzz[k] * ab_x + g_0_z_yyyz_xxxyyyzz[k];

                g_0_z_xyyyz_xxyyzzz[k] = -g_0_z_yyyz_xxyyzzz[k] * ab_x + g_0_z_yyyz_xxxyyzzz[k];

                g_0_z_xyyyz_xxyzzzz[k] = -g_0_z_yyyz_xxyzzzz[k] * ab_x + g_0_z_yyyz_xxxyzzzz[k];

                g_0_z_xyyyz_xxzzzzz[k] = -g_0_z_yyyz_xxzzzzz[k] * ab_x + g_0_z_yyyz_xxxzzzzz[k];

                g_0_z_xyyyz_xyyyyyy[k] = -g_0_z_yyyz_xyyyyyy[k] * ab_x + g_0_z_yyyz_xxyyyyyy[k];

                g_0_z_xyyyz_xyyyyyz[k] = -g_0_z_yyyz_xyyyyyz[k] * ab_x + g_0_z_yyyz_xxyyyyyz[k];

                g_0_z_xyyyz_xyyyyzz[k] = -g_0_z_yyyz_xyyyyzz[k] * ab_x + g_0_z_yyyz_xxyyyyzz[k];

                g_0_z_xyyyz_xyyyzzz[k] = -g_0_z_yyyz_xyyyzzz[k] * ab_x + g_0_z_yyyz_xxyyyzzz[k];

                g_0_z_xyyyz_xyyzzzz[k] = -g_0_z_yyyz_xyyzzzz[k] * ab_x + g_0_z_yyyz_xxyyzzzz[k];

                g_0_z_xyyyz_xyzzzzz[k] = -g_0_z_yyyz_xyzzzzz[k] * ab_x + g_0_z_yyyz_xxyzzzzz[k];

                g_0_z_xyyyz_xzzzzzz[k] = -g_0_z_yyyz_xzzzzzz[k] * ab_x + g_0_z_yyyz_xxzzzzzz[k];

                g_0_z_xyyyz_yyyyyyy[k] = -g_0_z_yyyz_yyyyyyy[k] * ab_x + g_0_z_yyyz_xyyyyyyy[k];

                g_0_z_xyyyz_yyyyyyz[k] = -g_0_z_yyyz_yyyyyyz[k] * ab_x + g_0_z_yyyz_xyyyyyyz[k];

                g_0_z_xyyyz_yyyyyzz[k] = -g_0_z_yyyz_yyyyyzz[k] * ab_x + g_0_z_yyyz_xyyyyyzz[k];

                g_0_z_xyyyz_yyyyzzz[k] = -g_0_z_yyyz_yyyyzzz[k] * ab_x + g_0_z_yyyz_xyyyyzzz[k];

                g_0_z_xyyyz_yyyzzzz[k] = -g_0_z_yyyz_yyyzzzz[k] * ab_x + g_0_z_yyyz_xyyyzzzz[k];

                g_0_z_xyyyz_yyzzzzz[k] = -g_0_z_yyyz_yyzzzzz[k] * ab_x + g_0_z_yyyz_xyyzzzzz[k];

                g_0_z_xyyyz_yzzzzzz[k] = -g_0_z_yyyz_yzzzzzz[k] * ab_x + g_0_z_yyyz_xyzzzzzz[k];

                g_0_z_xyyyz_zzzzzzz[k] = -g_0_z_yyyz_zzzzzzz[k] * ab_x + g_0_z_yyyz_xzzzzzzz[k];
            }

            /// Set up 1944-1980 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1944 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1945 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1946 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1947 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1948 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1949 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1950 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1951 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1952 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1953 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1954 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1955 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1956 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1957 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1958 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1959 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1960 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1961 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1962 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1963 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 1964 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 1965 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 1966 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 1967 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 1968 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 1969 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 1970 * ccomps * dcomps);

            auto g_0_z_xyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 1971 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 1972 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 1973 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 1974 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 1975 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 1976 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 1977 * ccomps * dcomps);

            auto g_0_z_xyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 1978 * ccomps * dcomps);

            auto g_0_z_xyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 1979 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzz_xxxxxxx, g_0_z_xyyzz_xxxxxxy, g_0_z_xyyzz_xxxxxxz, g_0_z_xyyzz_xxxxxyy, g_0_z_xyyzz_xxxxxyz, g_0_z_xyyzz_xxxxxzz, g_0_z_xyyzz_xxxxyyy, g_0_z_xyyzz_xxxxyyz, g_0_z_xyyzz_xxxxyzz, g_0_z_xyyzz_xxxxzzz, g_0_z_xyyzz_xxxyyyy, g_0_z_xyyzz_xxxyyyz, g_0_z_xyyzz_xxxyyzz, g_0_z_xyyzz_xxxyzzz, g_0_z_xyyzz_xxxzzzz, g_0_z_xyyzz_xxyyyyy, g_0_z_xyyzz_xxyyyyz, g_0_z_xyyzz_xxyyyzz, g_0_z_xyyzz_xxyyzzz, g_0_z_xyyzz_xxyzzzz, g_0_z_xyyzz_xxzzzzz, g_0_z_xyyzz_xyyyyyy, g_0_z_xyyzz_xyyyyyz, g_0_z_xyyzz_xyyyyzz, g_0_z_xyyzz_xyyyzzz, g_0_z_xyyzz_xyyzzzz, g_0_z_xyyzz_xyzzzzz, g_0_z_xyyzz_xzzzzzz, g_0_z_xyyzz_yyyyyyy, g_0_z_xyyzz_yyyyyyz, g_0_z_xyyzz_yyyyyzz, g_0_z_xyyzz_yyyyzzz, g_0_z_xyyzz_yyyzzzz, g_0_z_xyyzz_yyzzzzz, g_0_z_xyyzz_yzzzzzz, g_0_z_xyyzz_zzzzzzz, g_0_z_yyzz_xxxxxxx, g_0_z_yyzz_xxxxxxxx, g_0_z_yyzz_xxxxxxxy, g_0_z_yyzz_xxxxxxxz, g_0_z_yyzz_xxxxxxy, g_0_z_yyzz_xxxxxxyy, g_0_z_yyzz_xxxxxxyz, g_0_z_yyzz_xxxxxxz, g_0_z_yyzz_xxxxxxzz, g_0_z_yyzz_xxxxxyy, g_0_z_yyzz_xxxxxyyy, g_0_z_yyzz_xxxxxyyz, g_0_z_yyzz_xxxxxyz, g_0_z_yyzz_xxxxxyzz, g_0_z_yyzz_xxxxxzz, g_0_z_yyzz_xxxxxzzz, g_0_z_yyzz_xxxxyyy, g_0_z_yyzz_xxxxyyyy, g_0_z_yyzz_xxxxyyyz, g_0_z_yyzz_xxxxyyz, g_0_z_yyzz_xxxxyyzz, g_0_z_yyzz_xxxxyzz, g_0_z_yyzz_xxxxyzzz, g_0_z_yyzz_xxxxzzz, g_0_z_yyzz_xxxxzzzz, g_0_z_yyzz_xxxyyyy, g_0_z_yyzz_xxxyyyyy, g_0_z_yyzz_xxxyyyyz, g_0_z_yyzz_xxxyyyz, g_0_z_yyzz_xxxyyyzz, g_0_z_yyzz_xxxyyzz, g_0_z_yyzz_xxxyyzzz, g_0_z_yyzz_xxxyzzz, g_0_z_yyzz_xxxyzzzz, g_0_z_yyzz_xxxzzzz, g_0_z_yyzz_xxxzzzzz, g_0_z_yyzz_xxyyyyy, g_0_z_yyzz_xxyyyyyy, g_0_z_yyzz_xxyyyyyz, g_0_z_yyzz_xxyyyyz, g_0_z_yyzz_xxyyyyzz, g_0_z_yyzz_xxyyyzz, g_0_z_yyzz_xxyyyzzz, g_0_z_yyzz_xxyyzzz, g_0_z_yyzz_xxyyzzzz, g_0_z_yyzz_xxyzzzz, g_0_z_yyzz_xxyzzzzz, g_0_z_yyzz_xxzzzzz, g_0_z_yyzz_xxzzzzzz, g_0_z_yyzz_xyyyyyy, g_0_z_yyzz_xyyyyyyy, g_0_z_yyzz_xyyyyyyz, g_0_z_yyzz_xyyyyyz, g_0_z_yyzz_xyyyyyzz, g_0_z_yyzz_xyyyyzz, g_0_z_yyzz_xyyyyzzz, g_0_z_yyzz_xyyyzzz, g_0_z_yyzz_xyyyzzzz, g_0_z_yyzz_xyyzzzz, g_0_z_yyzz_xyyzzzzz, g_0_z_yyzz_xyzzzzz, g_0_z_yyzz_xyzzzzzz, g_0_z_yyzz_xzzzzzz, g_0_z_yyzz_xzzzzzzz, g_0_z_yyzz_yyyyyyy, g_0_z_yyzz_yyyyyyz, g_0_z_yyzz_yyyyyzz, g_0_z_yyzz_yyyyzzz, g_0_z_yyzz_yyyzzzz, g_0_z_yyzz_yyzzzzz, g_0_z_yyzz_yzzzzzz, g_0_z_yyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzz_xxxxxxx[k] = -g_0_z_yyzz_xxxxxxx[k] * ab_x + g_0_z_yyzz_xxxxxxxx[k];

                g_0_z_xyyzz_xxxxxxy[k] = -g_0_z_yyzz_xxxxxxy[k] * ab_x + g_0_z_yyzz_xxxxxxxy[k];

                g_0_z_xyyzz_xxxxxxz[k] = -g_0_z_yyzz_xxxxxxz[k] * ab_x + g_0_z_yyzz_xxxxxxxz[k];

                g_0_z_xyyzz_xxxxxyy[k] = -g_0_z_yyzz_xxxxxyy[k] * ab_x + g_0_z_yyzz_xxxxxxyy[k];

                g_0_z_xyyzz_xxxxxyz[k] = -g_0_z_yyzz_xxxxxyz[k] * ab_x + g_0_z_yyzz_xxxxxxyz[k];

                g_0_z_xyyzz_xxxxxzz[k] = -g_0_z_yyzz_xxxxxzz[k] * ab_x + g_0_z_yyzz_xxxxxxzz[k];

                g_0_z_xyyzz_xxxxyyy[k] = -g_0_z_yyzz_xxxxyyy[k] * ab_x + g_0_z_yyzz_xxxxxyyy[k];

                g_0_z_xyyzz_xxxxyyz[k] = -g_0_z_yyzz_xxxxyyz[k] * ab_x + g_0_z_yyzz_xxxxxyyz[k];

                g_0_z_xyyzz_xxxxyzz[k] = -g_0_z_yyzz_xxxxyzz[k] * ab_x + g_0_z_yyzz_xxxxxyzz[k];

                g_0_z_xyyzz_xxxxzzz[k] = -g_0_z_yyzz_xxxxzzz[k] * ab_x + g_0_z_yyzz_xxxxxzzz[k];

                g_0_z_xyyzz_xxxyyyy[k] = -g_0_z_yyzz_xxxyyyy[k] * ab_x + g_0_z_yyzz_xxxxyyyy[k];

                g_0_z_xyyzz_xxxyyyz[k] = -g_0_z_yyzz_xxxyyyz[k] * ab_x + g_0_z_yyzz_xxxxyyyz[k];

                g_0_z_xyyzz_xxxyyzz[k] = -g_0_z_yyzz_xxxyyzz[k] * ab_x + g_0_z_yyzz_xxxxyyzz[k];

                g_0_z_xyyzz_xxxyzzz[k] = -g_0_z_yyzz_xxxyzzz[k] * ab_x + g_0_z_yyzz_xxxxyzzz[k];

                g_0_z_xyyzz_xxxzzzz[k] = -g_0_z_yyzz_xxxzzzz[k] * ab_x + g_0_z_yyzz_xxxxzzzz[k];

                g_0_z_xyyzz_xxyyyyy[k] = -g_0_z_yyzz_xxyyyyy[k] * ab_x + g_0_z_yyzz_xxxyyyyy[k];

                g_0_z_xyyzz_xxyyyyz[k] = -g_0_z_yyzz_xxyyyyz[k] * ab_x + g_0_z_yyzz_xxxyyyyz[k];

                g_0_z_xyyzz_xxyyyzz[k] = -g_0_z_yyzz_xxyyyzz[k] * ab_x + g_0_z_yyzz_xxxyyyzz[k];

                g_0_z_xyyzz_xxyyzzz[k] = -g_0_z_yyzz_xxyyzzz[k] * ab_x + g_0_z_yyzz_xxxyyzzz[k];

                g_0_z_xyyzz_xxyzzzz[k] = -g_0_z_yyzz_xxyzzzz[k] * ab_x + g_0_z_yyzz_xxxyzzzz[k];

                g_0_z_xyyzz_xxzzzzz[k] = -g_0_z_yyzz_xxzzzzz[k] * ab_x + g_0_z_yyzz_xxxzzzzz[k];

                g_0_z_xyyzz_xyyyyyy[k] = -g_0_z_yyzz_xyyyyyy[k] * ab_x + g_0_z_yyzz_xxyyyyyy[k];

                g_0_z_xyyzz_xyyyyyz[k] = -g_0_z_yyzz_xyyyyyz[k] * ab_x + g_0_z_yyzz_xxyyyyyz[k];

                g_0_z_xyyzz_xyyyyzz[k] = -g_0_z_yyzz_xyyyyzz[k] * ab_x + g_0_z_yyzz_xxyyyyzz[k];

                g_0_z_xyyzz_xyyyzzz[k] = -g_0_z_yyzz_xyyyzzz[k] * ab_x + g_0_z_yyzz_xxyyyzzz[k];

                g_0_z_xyyzz_xyyzzzz[k] = -g_0_z_yyzz_xyyzzzz[k] * ab_x + g_0_z_yyzz_xxyyzzzz[k];

                g_0_z_xyyzz_xyzzzzz[k] = -g_0_z_yyzz_xyzzzzz[k] * ab_x + g_0_z_yyzz_xxyzzzzz[k];

                g_0_z_xyyzz_xzzzzzz[k] = -g_0_z_yyzz_xzzzzzz[k] * ab_x + g_0_z_yyzz_xxzzzzzz[k];

                g_0_z_xyyzz_yyyyyyy[k] = -g_0_z_yyzz_yyyyyyy[k] * ab_x + g_0_z_yyzz_xyyyyyyy[k];

                g_0_z_xyyzz_yyyyyyz[k] = -g_0_z_yyzz_yyyyyyz[k] * ab_x + g_0_z_yyzz_xyyyyyyz[k];

                g_0_z_xyyzz_yyyyyzz[k] = -g_0_z_yyzz_yyyyyzz[k] * ab_x + g_0_z_yyzz_xyyyyyzz[k];

                g_0_z_xyyzz_yyyyzzz[k] = -g_0_z_yyzz_yyyyzzz[k] * ab_x + g_0_z_yyzz_xyyyyzzz[k];

                g_0_z_xyyzz_yyyzzzz[k] = -g_0_z_yyzz_yyyzzzz[k] * ab_x + g_0_z_yyzz_xyyyzzzz[k];

                g_0_z_xyyzz_yyzzzzz[k] = -g_0_z_yyzz_yyzzzzz[k] * ab_x + g_0_z_yyzz_xyyzzzzz[k];

                g_0_z_xyyzz_yzzzzzz[k] = -g_0_z_yyzz_yzzzzzz[k] * ab_x + g_0_z_yyzz_xyzzzzzz[k];

                g_0_z_xyyzz_zzzzzzz[k] = -g_0_z_yyzz_zzzzzzz[k] * ab_x + g_0_z_yyzz_xzzzzzzz[k];
            }

            /// Set up 1980-2016 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 1980 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 1981 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 1982 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 1983 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 1984 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 1985 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 1986 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 1987 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 1988 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 1989 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 1990 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 1991 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 1992 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 1993 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 1994 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 1995 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 1996 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 1997 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 1998 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 1999 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2000 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2001 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2002 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2003 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2004 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2005 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2006 * ccomps * dcomps);

            auto g_0_z_xyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2007 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2008 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2009 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2010 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2011 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2012 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2013 * ccomps * dcomps);

            auto g_0_z_xyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2014 * ccomps * dcomps);

            auto g_0_z_xyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2015 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzz_xxxxxxx, g_0_z_xyzzz_xxxxxxy, g_0_z_xyzzz_xxxxxxz, g_0_z_xyzzz_xxxxxyy, g_0_z_xyzzz_xxxxxyz, g_0_z_xyzzz_xxxxxzz, g_0_z_xyzzz_xxxxyyy, g_0_z_xyzzz_xxxxyyz, g_0_z_xyzzz_xxxxyzz, g_0_z_xyzzz_xxxxzzz, g_0_z_xyzzz_xxxyyyy, g_0_z_xyzzz_xxxyyyz, g_0_z_xyzzz_xxxyyzz, g_0_z_xyzzz_xxxyzzz, g_0_z_xyzzz_xxxzzzz, g_0_z_xyzzz_xxyyyyy, g_0_z_xyzzz_xxyyyyz, g_0_z_xyzzz_xxyyyzz, g_0_z_xyzzz_xxyyzzz, g_0_z_xyzzz_xxyzzzz, g_0_z_xyzzz_xxzzzzz, g_0_z_xyzzz_xyyyyyy, g_0_z_xyzzz_xyyyyyz, g_0_z_xyzzz_xyyyyzz, g_0_z_xyzzz_xyyyzzz, g_0_z_xyzzz_xyyzzzz, g_0_z_xyzzz_xyzzzzz, g_0_z_xyzzz_xzzzzzz, g_0_z_xyzzz_yyyyyyy, g_0_z_xyzzz_yyyyyyz, g_0_z_xyzzz_yyyyyzz, g_0_z_xyzzz_yyyyzzz, g_0_z_xyzzz_yyyzzzz, g_0_z_xyzzz_yyzzzzz, g_0_z_xyzzz_yzzzzzz, g_0_z_xyzzz_zzzzzzz, g_0_z_yzzz_xxxxxxx, g_0_z_yzzz_xxxxxxxx, g_0_z_yzzz_xxxxxxxy, g_0_z_yzzz_xxxxxxxz, g_0_z_yzzz_xxxxxxy, g_0_z_yzzz_xxxxxxyy, g_0_z_yzzz_xxxxxxyz, g_0_z_yzzz_xxxxxxz, g_0_z_yzzz_xxxxxxzz, g_0_z_yzzz_xxxxxyy, g_0_z_yzzz_xxxxxyyy, g_0_z_yzzz_xxxxxyyz, g_0_z_yzzz_xxxxxyz, g_0_z_yzzz_xxxxxyzz, g_0_z_yzzz_xxxxxzz, g_0_z_yzzz_xxxxxzzz, g_0_z_yzzz_xxxxyyy, g_0_z_yzzz_xxxxyyyy, g_0_z_yzzz_xxxxyyyz, g_0_z_yzzz_xxxxyyz, g_0_z_yzzz_xxxxyyzz, g_0_z_yzzz_xxxxyzz, g_0_z_yzzz_xxxxyzzz, g_0_z_yzzz_xxxxzzz, g_0_z_yzzz_xxxxzzzz, g_0_z_yzzz_xxxyyyy, g_0_z_yzzz_xxxyyyyy, g_0_z_yzzz_xxxyyyyz, g_0_z_yzzz_xxxyyyz, g_0_z_yzzz_xxxyyyzz, g_0_z_yzzz_xxxyyzz, g_0_z_yzzz_xxxyyzzz, g_0_z_yzzz_xxxyzzz, g_0_z_yzzz_xxxyzzzz, g_0_z_yzzz_xxxzzzz, g_0_z_yzzz_xxxzzzzz, g_0_z_yzzz_xxyyyyy, g_0_z_yzzz_xxyyyyyy, g_0_z_yzzz_xxyyyyyz, g_0_z_yzzz_xxyyyyz, g_0_z_yzzz_xxyyyyzz, g_0_z_yzzz_xxyyyzz, g_0_z_yzzz_xxyyyzzz, g_0_z_yzzz_xxyyzzz, g_0_z_yzzz_xxyyzzzz, g_0_z_yzzz_xxyzzzz, g_0_z_yzzz_xxyzzzzz, g_0_z_yzzz_xxzzzzz, g_0_z_yzzz_xxzzzzzz, g_0_z_yzzz_xyyyyyy, g_0_z_yzzz_xyyyyyyy, g_0_z_yzzz_xyyyyyyz, g_0_z_yzzz_xyyyyyz, g_0_z_yzzz_xyyyyyzz, g_0_z_yzzz_xyyyyzz, g_0_z_yzzz_xyyyyzzz, g_0_z_yzzz_xyyyzzz, g_0_z_yzzz_xyyyzzzz, g_0_z_yzzz_xyyzzzz, g_0_z_yzzz_xyyzzzzz, g_0_z_yzzz_xyzzzzz, g_0_z_yzzz_xyzzzzzz, g_0_z_yzzz_xzzzzzz, g_0_z_yzzz_xzzzzzzz, g_0_z_yzzz_yyyyyyy, g_0_z_yzzz_yyyyyyz, g_0_z_yzzz_yyyyyzz, g_0_z_yzzz_yyyyzzz, g_0_z_yzzz_yyyzzzz, g_0_z_yzzz_yyzzzzz, g_0_z_yzzz_yzzzzzz, g_0_z_yzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzz_xxxxxxx[k] = -g_0_z_yzzz_xxxxxxx[k] * ab_x + g_0_z_yzzz_xxxxxxxx[k];

                g_0_z_xyzzz_xxxxxxy[k] = -g_0_z_yzzz_xxxxxxy[k] * ab_x + g_0_z_yzzz_xxxxxxxy[k];

                g_0_z_xyzzz_xxxxxxz[k] = -g_0_z_yzzz_xxxxxxz[k] * ab_x + g_0_z_yzzz_xxxxxxxz[k];

                g_0_z_xyzzz_xxxxxyy[k] = -g_0_z_yzzz_xxxxxyy[k] * ab_x + g_0_z_yzzz_xxxxxxyy[k];

                g_0_z_xyzzz_xxxxxyz[k] = -g_0_z_yzzz_xxxxxyz[k] * ab_x + g_0_z_yzzz_xxxxxxyz[k];

                g_0_z_xyzzz_xxxxxzz[k] = -g_0_z_yzzz_xxxxxzz[k] * ab_x + g_0_z_yzzz_xxxxxxzz[k];

                g_0_z_xyzzz_xxxxyyy[k] = -g_0_z_yzzz_xxxxyyy[k] * ab_x + g_0_z_yzzz_xxxxxyyy[k];

                g_0_z_xyzzz_xxxxyyz[k] = -g_0_z_yzzz_xxxxyyz[k] * ab_x + g_0_z_yzzz_xxxxxyyz[k];

                g_0_z_xyzzz_xxxxyzz[k] = -g_0_z_yzzz_xxxxyzz[k] * ab_x + g_0_z_yzzz_xxxxxyzz[k];

                g_0_z_xyzzz_xxxxzzz[k] = -g_0_z_yzzz_xxxxzzz[k] * ab_x + g_0_z_yzzz_xxxxxzzz[k];

                g_0_z_xyzzz_xxxyyyy[k] = -g_0_z_yzzz_xxxyyyy[k] * ab_x + g_0_z_yzzz_xxxxyyyy[k];

                g_0_z_xyzzz_xxxyyyz[k] = -g_0_z_yzzz_xxxyyyz[k] * ab_x + g_0_z_yzzz_xxxxyyyz[k];

                g_0_z_xyzzz_xxxyyzz[k] = -g_0_z_yzzz_xxxyyzz[k] * ab_x + g_0_z_yzzz_xxxxyyzz[k];

                g_0_z_xyzzz_xxxyzzz[k] = -g_0_z_yzzz_xxxyzzz[k] * ab_x + g_0_z_yzzz_xxxxyzzz[k];

                g_0_z_xyzzz_xxxzzzz[k] = -g_0_z_yzzz_xxxzzzz[k] * ab_x + g_0_z_yzzz_xxxxzzzz[k];

                g_0_z_xyzzz_xxyyyyy[k] = -g_0_z_yzzz_xxyyyyy[k] * ab_x + g_0_z_yzzz_xxxyyyyy[k];

                g_0_z_xyzzz_xxyyyyz[k] = -g_0_z_yzzz_xxyyyyz[k] * ab_x + g_0_z_yzzz_xxxyyyyz[k];

                g_0_z_xyzzz_xxyyyzz[k] = -g_0_z_yzzz_xxyyyzz[k] * ab_x + g_0_z_yzzz_xxxyyyzz[k];

                g_0_z_xyzzz_xxyyzzz[k] = -g_0_z_yzzz_xxyyzzz[k] * ab_x + g_0_z_yzzz_xxxyyzzz[k];

                g_0_z_xyzzz_xxyzzzz[k] = -g_0_z_yzzz_xxyzzzz[k] * ab_x + g_0_z_yzzz_xxxyzzzz[k];

                g_0_z_xyzzz_xxzzzzz[k] = -g_0_z_yzzz_xxzzzzz[k] * ab_x + g_0_z_yzzz_xxxzzzzz[k];

                g_0_z_xyzzz_xyyyyyy[k] = -g_0_z_yzzz_xyyyyyy[k] * ab_x + g_0_z_yzzz_xxyyyyyy[k];

                g_0_z_xyzzz_xyyyyyz[k] = -g_0_z_yzzz_xyyyyyz[k] * ab_x + g_0_z_yzzz_xxyyyyyz[k];

                g_0_z_xyzzz_xyyyyzz[k] = -g_0_z_yzzz_xyyyyzz[k] * ab_x + g_0_z_yzzz_xxyyyyzz[k];

                g_0_z_xyzzz_xyyyzzz[k] = -g_0_z_yzzz_xyyyzzz[k] * ab_x + g_0_z_yzzz_xxyyyzzz[k];

                g_0_z_xyzzz_xyyzzzz[k] = -g_0_z_yzzz_xyyzzzz[k] * ab_x + g_0_z_yzzz_xxyyzzzz[k];

                g_0_z_xyzzz_xyzzzzz[k] = -g_0_z_yzzz_xyzzzzz[k] * ab_x + g_0_z_yzzz_xxyzzzzz[k];

                g_0_z_xyzzz_xzzzzzz[k] = -g_0_z_yzzz_xzzzzzz[k] * ab_x + g_0_z_yzzz_xxzzzzzz[k];

                g_0_z_xyzzz_yyyyyyy[k] = -g_0_z_yzzz_yyyyyyy[k] * ab_x + g_0_z_yzzz_xyyyyyyy[k];

                g_0_z_xyzzz_yyyyyyz[k] = -g_0_z_yzzz_yyyyyyz[k] * ab_x + g_0_z_yzzz_xyyyyyyz[k];

                g_0_z_xyzzz_yyyyyzz[k] = -g_0_z_yzzz_yyyyyzz[k] * ab_x + g_0_z_yzzz_xyyyyyzz[k];

                g_0_z_xyzzz_yyyyzzz[k] = -g_0_z_yzzz_yyyyzzz[k] * ab_x + g_0_z_yzzz_xyyyyzzz[k];

                g_0_z_xyzzz_yyyzzzz[k] = -g_0_z_yzzz_yyyzzzz[k] * ab_x + g_0_z_yzzz_xyyyzzzz[k];

                g_0_z_xyzzz_yyzzzzz[k] = -g_0_z_yzzz_yyzzzzz[k] * ab_x + g_0_z_yzzz_xyyzzzzz[k];

                g_0_z_xyzzz_yzzzzzz[k] = -g_0_z_yzzz_yzzzzzz[k] * ab_x + g_0_z_yzzz_xyzzzzzz[k];

                g_0_z_xyzzz_zzzzzzz[k] = -g_0_z_yzzz_zzzzzzz[k] * ab_x + g_0_z_yzzz_xzzzzzzz[k];
            }

            /// Set up 2016-2052 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2016 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2017 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2018 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2019 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2020 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2021 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2022 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2023 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2024 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2025 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2026 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2027 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2028 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2029 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2030 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2031 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2032 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2033 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2034 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2035 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2036 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2037 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2038 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2039 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2040 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2041 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2042 * ccomps * dcomps);

            auto g_0_z_xzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2043 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2044 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2045 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2046 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2047 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2048 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2049 * ccomps * dcomps);

            auto g_0_z_xzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2050 * ccomps * dcomps);

            auto g_0_z_xzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2051 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzz_xxxxxxx, g_0_z_xzzzz_xxxxxxy, g_0_z_xzzzz_xxxxxxz, g_0_z_xzzzz_xxxxxyy, g_0_z_xzzzz_xxxxxyz, g_0_z_xzzzz_xxxxxzz, g_0_z_xzzzz_xxxxyyy, g_0_z_xzzzz_xxxxyyz, g_0_z_xzzzz_xxxxyzz, g_0_z_xzzzz_xxxxzzz, g_0_z_xzzzz_xxxyyyy, g_0_z_xzzzz_xxxyyyz, g_0_z_xzzzz_xxxyyzz, g_0_z_xzzzz_xxxyzzz, g_0_z_xzzzz_xxxzzzz, g_0_z_xzzzz_xxyyyyy, g_0_z_xzzzz_xxyyyyz, g_0_z_xzzzz_xxyyyzz, g_0_z_xzzzz_xxyyzzz, g_0_z_xzzzz_xxyzzzz, g_0_z_xzzzz_xxzzzzz, g_0_z_xzzzz_xyyyyyy, g_0_z_xzzzz_xyyyyyz, g_0_z_xzzzz_xyyyyzz, g_0_z_xzzzz_xyyyzzz, g_0_z_xzzzz_xyyzzzz, g_0_z_xzzzz_xyzzzzz, g_0_z_xzzzz_xzzzzzz, g_0_z_xzzzz_yyyyyyy, g_0_z_xzzzz_yyyyyyz, g_0_z_xzzzz_yyyyyzz, g_0_z_xzzzz_yyyyzzz, g_0_z_xzzzz_yyyzzzz, g_0_z_xzzzz_yyzzzzz, g_0_z_xzzzz_yzzzzzz, g_0_z_xzzzz_zzzzzzz, g_0_z_zzzz_xxxxxxx, g_0_z_zzzz_xxxxxxxx, g_0_z_zzzz_xxxxxxxy, g_0_z_zzzz_xxxxxxxz, g_0_z_zzzz_xxxxxxy, g_0_z_zzzz_xxxxxxyy, g_0_z_zzzz_xxxxxxyz, g_0_z_zzzz_xxxxxxz, g_0_z_zzzz_xxxxxxzz, g_0_z_zzzz_xxxxxyy, g_0_z_zzzz_xxxxxyyy, g_0_z_zzzz_xxxxxyyz, g_0_z_zzzz_xxxxxyz, g_0_z_zzzz_xxxxxyzz, g_0_z_zzzz_xxxxxzz, g_0_z_zzzz_xxxxxzzz, g_0_z_zzzz_xxxxyyy, g_0_z_zzzz_xxxxyyyy, g_0_z_zzzz_xxxxyyyz, g_0_z_zzzz_xxxxyyz, g_0_z_zzzz_xxxxyyzz, g_0_z_zzzz_xxxxyzz, g_0_z_zzzz_xxxxyzzz, g_0_z_zzzz_xxxxzzz, g_0_z_zzzz_xxxxzzzz, g_0_z_zzzz_xxxyyyy, g_0_z_zzzz_xxxyyyyy, g_0_z_zzzz_xxxyyyyz, g_0_z_zzzz_xxxyyyz, g_0_z_zzzz_xxxyyyzz, g_0_z_zzzz_xxxyyzz, g_0_z_zzzz_xxxyyzzz, g_0_z_zzzz_xxxyzzz, g_0_z_zzzz_xxxyzzzz, g_0_z_zzzz_xxxzzzz, g_0_z_zzzz_xxxzzzzz, g_0_z_zzzz_xxyyyyy, g_0_z_zzzz_xxyyyyyy, g_0_z_zzzz_xxyyyyyz, g_0_z_zzzz_xxyyyyz, g_0_z_zzzz_xxyyyyzz, g_0_z_zzzz_xxyyyzz, g_0_z_zzzz_xxyyyzzz, g_0_z_zzzz_xxyyzzz, g_0_z_zzzz_xxyyzzzz, g_0_z_zzzz_xxyzzzz, g_0_z_zzzz_xxyzzzzz, g_0_z_zzzz_xxzzzzz, g_0_z_zzzz_xxzzzzzz, g_0_z_zzzz_xyyyyyy, g_0_z_zzzz_xyyyyyyy, g_0_z_zzzz_xyyyyyyz, g_0_z_zzzz_xyyyyyz, g_0_z_zzzz_xyyyyyzz, g_0_z_zzzz_xyyyyzz, g_0_z_zzzz_xyyyyzzz, g_0_z_zzzz_xyyyzzz, g_0_z_zzzz_xyyyzzzz, g_0_z_zzzz_xyyzzzz, g_0_z_zzzz_xyyzzzzz, g_0_z_zzzz_xyzzzzz, g_0_z_zzzz_xyzzzzzz, g_0_z_zzzz_xzzzzzz, g_0_z_zzzz_xzzzzzzz, g_0_z_zzzz_yyyyyyy, g_0_z_zzzz_yyyyyyz, g_0_z_zzzz_yyyyyzz, g_0_z_zzzz_yyyyzzz, g_0_z_zzzz_yyyzzzz, g_0_z_zzzz_yyzzzzz, g_0_z_zzzz_yzzzzzz, g_0_z_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzz_xxxxxxx[k] = -g_0_z_zzzz_xxxxxxx[k] * ab_x + g_0_z_zzzz_xxxxxxxx[k];

                g_0_z_xzzzz_xxxxxxy[k] = -g_0_z_zzzz_xxxxxxy[k] * ab_x + g_0_z_zzzz_xxxxxxxy[k];

                g_0_z_xzzzz_xxxxxxz[k] = -g_0_z_zzzz_xxxxxxz[k] * ab_x + g_0_z_zzzz_xxxxxxxz[k];

                g_0_z_xzzzz_xxxxxyy[k] = -g_0_z_zzzz_xxxxxyy[k] * ab_x + g_0_z_zzzz_xxxxxxyy[k];

                g_0_z_xzzzz_xxxxxyz[k] = -g_0_z_zzzz_xxxxxyz[k] * ab_x + g_0_z_zzzz_xxxxxxyz[k];

                g_0_z_xzzzz_xxxxxzz[k] = -g_0_z_zzzz_xxxxxzz[k] * ab_x + g_0_z_zzzz_xxxxxxzz[k];

                g_0_z_xzzzz_xxxxyyy[k] = -g_0_z_zzzz_xxxxyyy[k] * ab_x + g_0_z_zzzz_xxxxxyyy[k];

                g_0_z_xzzzz_xxxxyyz[k] = -g_0_z_zzzz_xxxxyyz[k] * ab_x + g_0_z_zzzz_xxxxxyyz[k];

                g_0_z_xzzzz_xxxxyzz[k] = -g_0_z_zzzz_xxxxyzz[k] * ab_x + g_0_z_zzzz_xxxxxyzz[k];

                g_0_z_xzzzz_xxxxzzz[k] = -g_0_z_zzzz_xxxxzzz[k] * ab_x + g_0_z_zzzz_xxxxxzzz[k];

                g_0_z_xzzzz_xxxyyyy[k] = -g_0_z_zzzz_xxxyyyy[k] * ab_x + g_0_z_zzzz_xxxxyyyy[k];

                g_0_z_xzzzz_xxxyyyz[k] = -g_0_z_zzzz_xxxyyyz[k] * ab_x + g_0_z_zzzz_xxxxyyyz[k];

                g_0_z_xzzzz_xxxyyzz[k] = -g_0_z_zzzz_xxxyyzz[k] * ab_x + g_0_z_zzzz_xxxxyyzz[k];

                g_0_z_xzzzz_xxxyzzz[k] = -g_0_z_zzzz_xxxyzzz[k] * ab_x + g_0_z_zzzz_xxxxyzzz[k];

                g_0_z_xzzzz_xxxzzzz[k] = -g_0_z_zzzz_xxxzzzz[k] * ab_x + g_0_z_zzzz_xxxxzzzz[k];

                g_0_z_xzzzz_xxyyyyy[k] = -g_0_z_zzzz_xxyyyyy[k] * ab_x + g_0_z_zzzz_xxxyyyyy[k];

                g_0_z_xzzzz_xxyyyyz[k] = -g_0_z_zzzz_xxyyyyz[k] * ab_x + g_0_z_zzzz_xxxyyyyz[k];

                g_0_z_xzzzz_xxyyyzz[k] = -g_0_z_zzzz_xxyyyzz[k] * ab_x + g_0_z_zzzz_xxxyyyzz[k];

                g_0_z_xzzzz_xxyyzzz[k] = -g_0_z_zzzz_xxyyzzz[k] * ab_x + g_0_z_zzzz_xxxyyzzz[k];

                g_0_z_xzzzz_xxyzzzz[k] = -g_0_z_zzzz_xxyzzzz[k] * ab_x + g_0_z_zzzz_xxxyzzzz[k];

                g_0_z_xzzzz_xxzzzzz[k] = -g_0_z_zzzz_xxzzzzz[k] * ab_x + g_0_z_zzzz_xxxzzzzz[k];

                g_0_z_xzzzz_xyyyyyy[k] = -g_0_z_zzzz_xyyyyyy[k] * ab_x + g_0_z_zzzz_xxyyyyyy[k];

                g_0_z_xzzzz_xyyyyyz[k] = -g_0_z_zzzz_xyyyyyz[k] * ab_x + g_0_z_zzzz_xxyyyyyz[k];

                g_0_z_xzzzz_xyyyyzz[k] = -g_0_z_zzzz_xyyyyzz[k] * ab_x + g_0_z_zzzz_xxyyyyzz[k];

                g_0_z_xzzzz_xyyyzzz[k] = -g_0_z_zzzz_xyyyzzz[k] * ab_x + g_0_z_zzzz_xxyyyzzz[k];

                g_0_z_xzzzz_xyyzzzz[k] = -g_0_z_zzzz_xyyzzzz[k] * ab_x + g_0_z_zzzz_xxyyzzzz[k];

                g_0_z_xzzzz_xyzzzzz[k] = -g_0_z_zzzz_xyzzzzz[k] * ab_x + g_0_z_zzzz_xxyzzzzz[k];

                g_0_z_xzzzz_xzzzzzz[k] = -g_0_z_zzzz_xzzzzzz[k] * ab_x + g_0_z_zzzz_xxzzzzzz[k];

                g_0_z_xzzzz_yyyyyyy[k] = -g_0_z_zzzz_yyyyyyy[k] * ab_x + g_0_z_zzzz_xyyyyyyy[k];

                g_0_z_xzzzz_yyyyyyz[k] = -g_0_z_zzzz_yyyyyyz[k] * ab_x + g_0_z_zzzz_xyyyyyyz[k];

                g_0_z_xzzzz_yyyyyzz[k] = -g_0_z_zzzz_yyyyyzz[k] * ab_x + g_0_z_zzzz_xyyyyyzz[k];

                g_0_z_xzzzz_yyyyzzz[k] = -g_0_z_zzzz_yyyyzzz[k] * ab_x + g_0_z_zzzz_xyyyyzzz[k];

                g_0_z_xzzzz_yyyzzzz[k] = -g_0_z_zzzz_yyyzzzz[k] * ab_x + g_0_z_zzzz_xyyyzzzz[k];

                g_0_z_xzzzz_yyzzzzz[k] = -g_0_z_zzzz_yyzzzzz[k] * ab_x + g_0_z_zzzz_xyyzzzzz[k];

                g_0_z_xzzzz_yzzzzzz[k] = -g_0_z_zzzz_yzzzzzz[k] * ab_x + g_0_z_zzzz_xyzzzzzz[k];

                g_0_z_xzzzz_zzzzzzz[k] = -g_0_z_zzzz_zzzzzzz[k] * ab_x + g_0_z_zzzz_xzzzzzzz[k];
            }

            /// Set up 2052-2088 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyy_xxxxxxx = cbuffer.data(hk_geom_01_off + 2052 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxxxy = cbuffer.data(hk_geom_01_off + 2053 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxxxz = cbuffer.data(hk_geom_01_off + 2054 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxxyy = cbuffer.data(hk_geom_01_off + 2055 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxxyz = cbuffer.data(hk_geom_01_off + 2056 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxxzz = cbuffer.data(hk_geom_01_off + 2057 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxyyy = cbuffer.data(hk_geom_01_off + 2058 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxyyz = cbuffer.data(hk_geom_01_off + 2059 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxyzz = cbuffer.data(hk_geom_01_off + 2060 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxzzz = cbuffer.data(hk_geom_01_off + 2061 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyyyy = cbuffer.data(hk_geom_01_off + 2062 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyyyz = cbuffer.data(hk_geom_01_off + 2063 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyyzz = cbuffer.data(hk_geom_01_off + 2064 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyzzz = cbuffer.data(hk_geom_01_off + 2065 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxzzzz = cbuffer.data(hk_geom_01_off + 2066 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyyyy = cbuffer.data(hk_geom_01_off + 2067 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyyyz = cbuffer.data(hk_geom_01_off + 2068 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyyzz = cbuffer.data(hk_geom_01_off + 2069 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyzzz = cbuffer.data(hk_geom_01_off + 2070 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyzzzz = cbuffer.data(hk_geom_01_off + 2071 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxzzzzz = cbuffer.data(hk_geom_01_off + 2072 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyyyy = cbuffer.data(hk_geom_01_off + 2073 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyyyz = cbuffer.data(hk_geom_01_off + 2074 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyyzz = cbuffer.data(hk_geom_01_off + 2075 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyzzz = cbuffer.data(hk_geom_01_off + 2076 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyzzzz = cbuffer.data(hk_geom_01_off + 2077 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyzzzzz = cbuffer.data(hk_geom_01_off + 2078 * ccomps * dcomps);

            auto g_0_z_yyyyy_xzzzzzz = cbuffer.data(hk_geom_01_off + 2079 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyyyy = cbuffer.data(hk_geom_01_off + 2080 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyyyz = cbuffer.data(hk_geom_01_off + 2081 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyyzz = cbuffer.data(hk_geom_01_off + 2082 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyzzz = cbuffer.data(hk_geom_01_off + 2083 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyzzzz = cbuffer.data(hk_geom_01_off + 2084 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyzzzzz = cbuffer.data(hk_geom_01_off + 2085 * ccomps * dcomps);

            auto g_0_z_yyyyy_yzzzzzz = cbuffer.data(hk_geom_01_off + 2086 * ccomps * dcomps);

            auto g_0_z_yyyyy_zzzzzzz = cbuffer.data(hk_geom_01_off + 2087 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_xxxxxxx, g_0_z_yyyy_xxxxxxxy, g_0_z_yyyy_xxxxxxy, g_0_z_yyyy_xxxxxxyy, g_0_z_yyyy_xxxxxxyz, g_0_z_yyyy_xxxxxxz, g_0_z_yyyy_xxxxxyy, g_0_z_yyyy_xxxxxyyy, g_0_z_yyyy_xxxxxyyz, g_0_z_yyyy_xxxxxyz, g_0_z_yyyy_xxxxxyzz, g_0_z_yyyy_xxxxxzz, g_0_z_yyyy_xxxxyyy, g_0_z_yyyy_xxxxyyyy, g_0_z_yyyy_xxxxyyyz, g_0_z_yyyy_xxxxyyz, g_0_z_yyyy_xxxxyyzz, g_0_z_yyyy_xxxxyzz, g_0_z_yyyy_xxxxyzzz, g_0_z_yyyy_xxxxzzz, g_0_z_yyyy_xxxyyyy, g_0_z_yyyy_xxxyyyyy, g_0_z_yyyy_xxxyyyyz, g_0_z_yyyy_xxxyyyz, g_0_z_yyyy_xxxyyyzz, g_0_z_yyyy_xxxyyzz, g_0_z_yyyy_xxxyyzzz, g_0_z_yyyy_xxxyzzz, g_0_z_yyyy_xxxyzzzz, g_0_z_yyyy_xxxzzzz, g_0_z_yyyy_xxyyyyy, g_0_z_yyyy_xxyyyyyy, g_0_z_yyyy_xxyyyyyz, g_0_z_yyyy_xxyyyyz, g_0_z_yyyy_xxyyyyzz, g_0_z_yyyy_xxyyyzz, g_0_z_yyyy_xxyyyzzz, g_0_z_yyyy_xxyyzzz, g_0_z_yyyy_xxyyzzzz, g_0_z_yyyy_xxyzzzz, g_0_z_yyyy_xxyzzzzz, g_0_z_yyyy_xxzzzzz, g_0_z_yyyy_xyyyyyy, g_0_z_yyyy_xyyyyyyy, g_0_z_yyyy_xyyyyyyz, g_0_z_yyyy_xyyyyyz, g_0_z_yyyy_xyyyyyzz, g_0_z_yyyy_xyyyyzz, g_0_z_yyyy_xyyyyzzz, g_0_z_yyyy_xyyyzzz, g_0_z_yyyy_xyyyzzzz, g_0_z_yyyy_xyyzzzz, g_0_z_yyyy_xyyzzzzz, g_0_z_yyyy_xyzzzzz, g_0_z_yyyy_xyzzzzzz, g_0_z_yyyy_xzzzzzz, g_0_z_yyyy_yyyyyyy, g_0_z_yyyy_yyyyyyyy, g_0_z_yyyy_yyyyyyyz, g_0_z_yyyy_yyyyyyz, g_0_z_yyyy_yyyyyyzz, g_0_z_yyyy_yyyyyzz, g_0_z_yyyy_yyyyyzzz, g_0_z_yyyy_yyyyzzz, g_0_z_yyyy_yyyyzzzz, g_0_z_yyyy_yyyzzzz, g_0_z_yyyy_yyyzzzzz, g_0_z_yyyy_yyzzzzz, g_0_z_yyyy_yyzzzzzz, g_0_z_yyyy_yzzzzzz, g_0_z_yyyy_yzzzzzzz, g_0_z_yyyy_zzzzzzz, g_0_z_yyyyy_xxxxxxx, g_0_z_yyyyy_xxxxxxy, g_0_z_yyyyy_xxxxxxz, g_0_z_yyyyy_xxxxxyy, g_0_z_yyyyy_xxxxxyz, g_0_z_yyyyy_xxxxxzz, g_0_z_yyyyy_xxxxyyy, g_0_z_yyyyy_xxxxyyz, g_0_z_yyyyy_xxxxyzz, g_0_z_yyyyy_xxxxzzz, g_0_z_yyyyy_xxxyyyy, g_0_z_yyyyy_xxxyyyz, g_0_z_yyyyy_xxxyyzz, g_0_z_yyyyy_xxxyzzz, g_0_z_yyyyy_xxxzzzz, g_0_z_yyyyy_xxyyyyy, g_0_z_yyyyy_xxyyyyz, g_0_z_yyyyy_xxyyyzz, g_0_z_yyyyy_xxyyzzz, g_0_z_yyyyy_xxyzzzz, g_0_z_yyyyy_xxzzzzz, g_0_z_yyyyy_xyyyyyy, g_0_z_yyyyy_xyyyyyz, g_0_z_yyyyy_xyyyyzz, g_0_z_yyyyy_xyyyzzz, g_0_z_yyyyy_xyyzzzz, g_0_z_yyyyy_xyzzzzz, g_0_z_yyyyy_xzzzzzz, g_0_z_yyyyy_yyyyyyy, g_0_z_yyyyy_yyyyyyz, g_0_z_yyyyy_yyyyyzz, g_0_z_yyyyy_yyyyzzz, g_0_z_yyyyy_yyyzzzz, g_0_z_yyyyy_yyzzzzz, g_0_z_yyyyy_yzzzzzz, g_0_z_yyyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyy_xxxxxxx[k] = -g_0_z_yyyy_xxxxxxx[k] * ab_y + g_0_z_yyyy_xxxxxxxy[k];

                g_0_z_yyyyy_xxxxxxy[k] = -g_0_z_yyyy_xxxxxxy[k] * ab_y + g_0_z_yyyy_xxxxxxyy[k];

                g_0_z_yyyyy_xxxxxxz[k] = -g_0_z_yyyy_xxxxxxz[k] * ab_y + g_0_z_yyyy_xxxxxxyz[k];

                g_0_z_yyyyy_xxxxxyy[k] = -g_0_z_yyyy_xxxxxyy[k] * ab_y + g_0_z_yyyy_xxxxxyyy[k];

                g_0_z_yyyyy_xxxxxyz[k] = -g_0_z_yyyy_xxxxxyz[k] * ab_y + g_0_z_yyyy_xxxxxyyz[k];

                g_0_z_yyyyy_xxxxxzz[k] = -g_0_z_yyyy_xxxxxzz[k] * ab_y + g_0_z_yyyy_xxxxxyzz[k];

                g_0_z_yyyyy_xxxxyyy[k] = -g_0_z_yyyy_xxxxyyy[k] * ab_y + g_0_z_yyyy_xxxxyyyy[k];

                g_0_z_yyyyy_xxxxyyz[k] = -g_0_z_yyyy_xxxxyyz[k] * ab_y + g_0_z_yyyy_xxxxyyyz[k];

                g_0_z_yyyyy_xxxxyzz[k] = -g_0_z_yyyy_xxxxyzz[k] * ab_y + g_0_z_yyyy_xxxxyyzz[k];

                g_0_z_yyyyy_xxxxzzz[k] = -g_0_z_yyyy_xxxxzzz[k] * ab_y + g_0_z_yyyy_xxxxyzzz[k];

                g_0_z_yyyyy_xxxyyyy[k] = -g_0_z_yyyy_xxxyyyy[k] * ab_y + g_0_z_yyyy_xxxyyyyy[k];

                g_0_z_yyyyy_xxxyyyz[k] = -g_0_z_yyyy_xxxyyyz[k] * ab_y + g_0_z_yyyy_xxxyyyyz[k];

                g_0_z_yyyyy_xxxyyzz[k] = -g_0_z_yyyy_xxxyyzz[k] * ab_y + g_0_z_yyyy_xxxyyyzz[k];

                g_0_z_yyyyy_xxxyzzz[k] = -g_0_z_yyyy_xxxyzzz[k] * ab_y + g_0_z_yyyy_xxxyyzzz[k];

                g_0_z_yyyyy_xxxzzzz[k] = -g_0_z_yyyy_xxxzzzz[k] * ab_y + g_0_z_yyyy_xxxyzzzz[k];

                g_0_z_yyyyy_xxyyyyy[k] = -g_0_z_yyyy_xxyyyyy[k] * ab_y + g_0_z_yyyy_xxyyyyyy[k];

                g_0_z_yyyyy_xxyyyyz[k] = -g_0_z_yyyy_xxyyyyz[k] * ab_y + g_0_z_yyyy_xxyyyyyz[k];

                g_0_z_yyyyy_xxyyyzz[k] = -g_0_z_yyyy_xxyyyzz[k] * ab_y + g_0_z_yyyy_xxyyyyzz[k];

                g_0_z_yyyyy_xxyyzzz[k] = -g_0_z_yyyy_xxyyzzz[k] * ab_y + g_0_z_yyyy_xxyyyzzz[k];

                g_0_z_yyyyy_xxyzzzz[k] = -g_0_z_yyyy_xxyzzzz[k] * ab_y + g_0_z_yyyy_xxyyzzzz[k];

                g_0_z_yyyyy_xxzzzzz[k] = -g_0_z_yyyy_xxzzzzz[k] * ab_y + g_0_z_yyyy_xxyzzzzz[k];

                g_0_z_yyyyy_xyyyyyy[k] = -g_0_z_yyyy_xyyyyyy[k] * ab_y + g_0_z_yyyy_xyyyyyyy[k];

                g_0_z_yyyyy_xyyyyyz[k] = -g_0_z_yyyy_xyyyyyz[k] * ab_y + g_0_z_yyyy_xyyyyyyz[k];

                g_0_z_yyyyy_xyyyyzz[k] = -g_0_z_yyyy_xyyyyzz[k] * ab_y + g_0_z_yyyy_xyyyyyzz[k];

                g_0_z_yyyyy_xyyyzzz[k] = -g_0_z_yyyy_xyyyzzz[k] * ab_y + g_0_z_yyyy_xyyyyzzz[k];

                g_0_z_yyyyy_xyyzzzz[k] = -g_0_z_yyyy_xyyzzzz[k] * ab_y + g_0_z_yyyy_xyyyzzzz[k];

                g_0_z_yyyyy_xyzzzzz[k] = -g_0_z_yyyy_xyzzzzz[k] * ab_y + g_0_z_yyyy_xyyzzzzz[k];

                g_0_z_yyyyy_xzzzzzz[k] = -g_0_z_yyyy_xzzzzzz[k] * ab_y + g_0_z_yyyy_xyzzzzzz[k];

                g_0_z_yyyyy_yyyyyyy[k] = -g_0_z_yyyy_yyyyyyy[k] * ab_y + g_0_z_yyyy_yyyyyyyy[k];

                g_0_z_yyyyy_yyyyyyz[k] = -g_0_z_yyyy_yyyyyyz[k] * ab_y + g_0_z_yyyy_yyyyyyyz[k];

                g_0_z_yyyyy_yyyyyzz[k] = -g_0_z_yyyy_yyyyyzz[k] * ab_y + g_0_z_yyyy_yyyyyyzz[k];

                g_0_z_yyyyy_yyyyzzz[k] = -g_0_z_yyyy_yyyyzzz[k] * ab_y + g_0_z_yyyy_yyyyyzzz[k];

                g_0_z_yyyyy_yyyzzzz[k] = -g_0_z_yyyy_yyyzzzz[k] * ab_y + g_0_z_yyyy_yyyyzzzz[k];

                g_0_z_yyyyy_yyzzzzz[k] = -g_0_z_yyyy_yyzzzzz[k] * ab_y + g_0_z_yyyy_yyyzzzzz[k];

                g_0_z_yyyyy_yzzzzzz[k] = -g_0_z_yyyy_yzzzzzz[k] * ab_y + g_0_z_yyyy_yyzzzzzz[k];

                g_0_z_yyyyy_zzzzzzz[k] = -g_0_z_yyyy_zzzzzzz[k] * ab_y + g_0_z_yyyy_yzzzzzzz[k];
            }

            /// Set up 2088-2124 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2088 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2089 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2090 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2091 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2092 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2093 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2094 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2095 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2096 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2097 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2098 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2099 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2100 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2101 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2102 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2103 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2104 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2105 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2106 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2107 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2108 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2109 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2110 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2111 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2112 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2113 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2114 * ccomps * dcomps);

            auto g_0_z_yyyyz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2115 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2116 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2117 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2118 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2119 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2120 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2121 * ccomps * dcomps);

            auto g_0_z_yyyyz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2122 * ccomps * dcomps);

            auto g_0_z_yyyyz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2123 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyz_xxxxxxx, g_0_z_yyyyz_xxxxxxy, g_0_z_yyyyz_xxxxxxz, g_0_z_yyyyz_xxxxxyy, g_0_z_yyyyz_xxxxxyz, g_0_z_yyyyz_xxxxxzz, g_0_z_yyyyz_xxxxyyy, g_0_z_yyyyz_xxxxyyz, g_0_z_yyyyz_xxxxyzz, g_0_z_yyyyz_xxxxzzz, g_0_z_yyyyz_xxxyyyy, g_0_z_yyyyz_xxxyyyz, g_0_z_yyyyz_xxxyyzz, g_0_z_yyyyz_xxxyzzz, g_0_z_yyyyz_xxxzzzz, g_0_z_yyyyz_xxyyyyy, g_0_z_yyyyz_xxyyyyz, g_0_z_yyyyz_xxyyyzz, g_0_z_yyyyz_xxyyzzz, g_0_z_yyyyz_xxyzzzz, g_0_z_yyyyz_xxzzzzz, g_0_z_yyyyz_xyyyyyy, g_0_z_yyyyz_xyyyyyz, g_0_z_yyyyz_xyyyyzz, g_0_z_yyyyz_xyyyzzz, g_0_z_yyyyz_xyyzzzz, g_0_z_yyyyz_xyzzzzz, g_0_z_yyyyz_xzzzzzz, g_0_z_yyyyz_yyyyyyy, g_0_z_yyyyz_yyyyyyz, g_0_z_yyyyz_yyyyyzz, g_0_z_yyyyz_yyyyzzz, g_0_z_yyyyz_yyyzzzz, g_0_z_yyyyz_yyzzzzz, g_0_z_yyyyz_yzzzzzz, g_0_z_yyyyz_zzzzzzz, g_0_z_yyyz_xxxxxxx, g_0_z_yyyz_xxxxxxxy, g_0_z_yyyz_xxxxxxy, g_0_z_yyyz_xxxxxxyy, g_0_z_yyyz_xxxxxxyz, g_0_z_yyyz_xxxxxxz, g_0_z_yyyz_xxxxxyy, g_0_z_yyyz_xxxxxyyy, g_0_z_yyyz_xxxxxyyz, g_0_z_yyyz_xxxxxyz, g_0_z_yyyz_xxxxxyzz, g_0_z_yyyz_xxxxxzz, g_0_z_yyyz_xxxxyyy, g_0_z_yyyz_xxxxyyyy, g_0_z_yyyz_xxxxyyyz, g_0_z_yyyz_xxxxyyz, g_0_z_yyyz_xxxxyyzz, g_0_z_yyyz_xxxxyzz, g_0_z_yyyz_xxxxyzzz, g_0_z_yyyz_xxxxzzz, g_0_z_yyyz_xxxyyyy, g_0_z_yyyz_xxxyyyyy, g_0_z_yyyz_xxxyyyyz, g_0_z_yyyz_xxxyyyz, g_0_z_yyyz_xxxyyyzz, g_0_z_yyyz_xxxyyzz, g_0_z_yyyz_xxxyyzzz, g_0_z_yyyz_xxxyzzz, g_0_z_yyyz_xxxyzzzz, g_0_z_yyyz_xxxzzzz, g_0_z_yyyz_xxyyyyy, g_0_z_yyyz_xxyyyyyy, g_0_z_yyyz_xxyyyyyz, g_0_z_yyyz_xxyyyyz, g_0_z_yyyz_xxyyyyzz, g_0_z_yyyz_xxyyyzz, g_0_z_yyyz_xxyyyzzz, g_0_z_yyyz_xxyyzzz, g_0_z_yyyz_xxyyzzzz, g_0_z_yyyz_xxyzzzz, g_0_z_yyyz_xxyzzzzz, g_0_z_yyyz_xxzzzzz, g_0_z_yyyz_xyyyyyy, g_0_z_yyyz_xyyyyyyy, g_0_z_yyyz_xyyyyyyz, g_0_z_yyyz_xyyyyyz, g_0_z_yyyz_xyyyyyzz, g_0_z_yyyz_xyyyyzz, g_0_z_yyyz_xyyyyzzz, g_0_z_yyyz_xyyyzzz, g_0_z_yyyz_xyyyzzzz, g_0_z_yyyz_xyyzzzz, g_0_z_yyyz_xyyzzzzz, g_0_z_yyyz_xyzzzzz, g_0_z_yyyz_xyzzzzzz, g_0_z_yyyz_xzzzzzz, g_0_z_yyyz_yyyyyyy, g_0_z_yyyz_yyyyyyyy, g_0_z_yyyz_yyyyyyyz, g_0_z_yyyz_yyyyyyz, g_0_z_yyyz_yyyyyyzz, g_0_z_yyyz_yyyyyzz, g_0_z_yyyz_yyyyyzzz, g_0_z_yyyz_yyyyzzz, g_0_z_yyyz_yyyyzzzz, g_0_z_yyyz_yyyzzzz, g_0_z_yyyz_yyyzzzzz, g_0_z_yyyz_yyzzzzz, g_0_z_yyyz_yyzzzzzz, g_0_z_yyyz_yzzzzzz, g_0_z_yyyz_yzzzzzzz, g_0_z_yyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyz_xxxxxxx[k] = -g_0_z_yyyz_xxxxxxx[k] * ab_y + g_0_z_yyyz_xxxxxxxy[k];

                g_0_z_yyyyz_xxxxxxy[k] = -g_0_z_yyyz_xxxxxxy[k] * ab_y + g_0_z_yyyz_xxxxxxyy[k];

                g_0_z_yyyyz_xxxxxxz[k] = -g_0_z_yyyz_xxxxxxz[k] * ab_y + g_0_z_yyyz_xxxxxxyz[k];

                g_0_z_yyyyz_xxxxxyy[k] = -g_0_z_yyyz_xxxxxyy[k] * ab_y + g_0_z_yyyz_xxxxxyyy[k];

                g_0_z_yyyyz_xxxxxyz[k] = -g_0_z_yyyz_xxxxxyz[k] * ab_y + g_0_z_yyyz_xxxxxyyz[k];

                g_0_z_yyyyz_xxxxxzz[k] = -g_0_z_yyyz_xxxxxzz[k] * ab_y + g_0_z_yyyz_xxxxxyzz[k];

                g_0_z_yyyyz_xxxxyyy[k] = -g_0_z_yyyz_xxxxyyy[k] * ab_y + g_0_z_yyyz_xxxxyyyy[k];

                g_0_z_yyyyz_xxxxyyz[k] = -g_0_z_yyyz_xxxxyyz[k] * ab_y + g_0_z_yyyz_xxxxyyyz[k];

                g_0_z_yyyyz_xxxxyzz[k] = -g_0_z_yyyz_xxxxyzz[k] * ab_y + g_0_z_yyyz_xxxxyyzz[k];

                g_0_z_yyyyz_xxxxzzz[k] = -g_0_z_yyyz_xxxxzzz[k] * ab_y + g_0_z_yyyz_xxxxyzzz[k];

                g_0_z_yyyyz_xxxyyyy[k] = -g_0_z_yyyz_xxxyyyy[k] * ab_y + g_0_z_yyyz_xxxyyyyy[k];

                g_0_z_yyyyz_xxxyyyz[k] = -g_0_z_yyyz_xxxyyyz[k] * ab_y + g_0_z_yyyz_xxxyyyyz[k];

                g_0_z_yyyyz_xxxyyzz[k] = -g_0_z_yyyz_xxxyyzz[k] * ab_y + g_0_z_yyyz_xxxyyyzz[k];

                g_0_z_yyyyz_xxxyzzz[k] = -g_0_z_yyyz_xxxyzzz[k] * ab_y + g_0_z_yyyz_xxxyyzzz[k];

                g_0_z_yyyyz_xxxzzzz[k] = -g_0_z_yyyz_xxxzzzz[k] * ab_y + g_0_z_yyyz_xxxyzzzz[k];

                g_0_z_yyyyz_xxyyyyy[k] = -g_0_z_yyyz_xxyyyyy[k] * ab_y + g_0_z_yyyz_xxyyyyyy[k];

                g_0_z_yyyyz_xxyyyyz[k] = -g_0_z_yyyz_xxyyyyz[k] * ab_y + g_0_z_yyyz_xxyyyyyz[k];

                g_0_z_yyyyz_xxyyyzz[k] = -g_0_z_yyyz_xxyyyzz[k] * ab_y + g_0_z_yyyz_xxyyyyzz[k];

                g_0_z_yyyyz_xxyyzzz[k] = -g_0_z_yyyz_xxyyzzz[k] * ab_y + g_0_z_yyyz_xxyyyzzz[k];

                g_0_z_yyyyz_xxyzzzz[k] = -g_0_z_yyyz_xxyzzzz[k] * ab_y + g_0_z_yyyz_xxyyzzzz[k];

                g_0_z_yyyyz_xxzzzzz[k] = -g_0_z_yyyz_xxzzzzz[k] * ab_y + g_0_z_yyyz_xxyzzzzz[k];

                g_0_z_yyyyz_xyyyyyy[k] = -g_0_z_yyyz_xyyyyyy[k] * ab_y + g_0_z_yyyz_xyyyyyyy[k];

                g_0_z_yyyyz_xyyyyyz[k] = -g_0_z_yyyz_xyyyyyz[k] * ab_y + g_0_z_yyyz_xyyyyyyz[k];

                g_0_z_yyyyz_xyyyyzz[k] = -g_0_z_yyyz_xyyyyzz[k] * ab_y + g_0_z_yyyz_xyyyyyzz[k];

                g_0_z_yyyyz_xyyyzzz[k] = -g_0_z_yyyz_xyyyzzz[k] * ab_y + g_0_z_yyyz_xyyyyzzz[k];

                g_0_z_yyyyz_xyyzzzz[k] = -g_0_z_yyyz_xyyzzzz[k] * ab_y + g_0_z_yyyz_xyyyzzzz[k];

                g_0_z_yyyyz_xyzzzzz[k] = -g_0_z_yyyz_xyzzzzz[k] * ab_y + g_0_z_yyyz_xyyzzzzz[k];

                g_0_z_yyyyz_xzzzzzz[k] = -g_0_z_yyyz_xzzzzzz[k] * ab_y + g_0_z_yyyz_xyzzzzzz[k];

                g_0_z_yyyyz_yyyyyyy[k] = -g_0_z_yyyz_yyyyyyy[k] * ab_y + g_0_z_yyyz_yyyyyyyy[k];

                g_0_z_yyyyz_yyyyyyz[k] = -g_0_z_yyyz_yyyyyyz[k] * ab_y + g_0_z_yyyz_yyyyyyyz[k];

                g_0_z_yyyyz_yyyyyzz[k] = -g_0_z_yyyz_yyyyyzz[k] * ab_y + g_0_z_yyyz_yyyyyyzz[k];

                g_0_z_yyyyz_yyyyzzz[k] = -g_0_z_yyyz_yyyyzzz[k] * ab_y + g_0_z_yyyz_yyyyyzzz[k];

                g_0_z_yyyyz_yyyzzzz[k] = -g_0_z_yyyz_yyyzzzz[k] * ab_y + g_0_z_yyyz_yyyyzzzz[k];

                g_0_z_yyyyz_yyzzzzz[k] = -g_0_z_yyyz_yyzzzzz[k] * ab_y + g_0_z_yyyz_yyyzzzzz[k];

                g_0_z_yyyyz_yzzzzzz[k] = -g_0_z_yyyz_yzzzzzz[k] * ab_y + g_0_z_yyyz_yyzzzzzz[k];

                g_0_z_yyyyz_zzzzzzz[k] = -g_0_z_yyyz_zzzzzzz[k] * ab_y + g_0_z_yyyz_yzzzzzzz[k];
            }

            /// Set up 2124-2160 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2124 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2125 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2126 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2127 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2128 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2129 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2130 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2131 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2132 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2133 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2134 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2135 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2136 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2137 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2138 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2139 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2140 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2141 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2142 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2143 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2144 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2145 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2146 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2147 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2148 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2149 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2150 * ccomps * dcomps);

            auto g_0_z_yyyzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2151 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2152 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2153 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2154 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2155 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2156 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2157 * ccomps * dcomps);

            auto g_0_z_yyyzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2158 * ccomps * dcomps);

            auto g_0_z_yyyzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2159 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzz_xxxxxxx, g_0_z_yyyzz_xxxxxxy, g_0_z_yyyzz_xxxxxxz, g_0_z_yyyzz_xxxxxyy, g_0_z_yyyzz_xxxxxyz, g_0_z_yyyzz_xxxxxzz, g_0_z_yyyzz_xxxxyyy, g_0_z_yyyzz_xxxxyyz, g_0_z_yyyzz_xxxxyzz, g_0_z_yyyzz_xxxxzzz, g_0_z_yyyzz_xxxyyyy, g_0_z_yyyzz_xxxyyyz, g_0_z_yyyzz_xxxyyzz, g_0_z_yyyzz_xxxyzzz, g_0_z_yyyzz_xxxzzzz, g_0_z_yyyzz_xxyyyyy, g_0_z_yyyzz_xxyyyyz, g_0_z_yyyzz_xxyyyzz, g_0_z_yyyzz_xxyyzzz, g_0_z_yyyzz_xxyzzzz, g_0_z_yyyzz_xxzzzzz, g_0_z_yyyzz_xyyyyyy, g_0_z_yyyzz_xyyyyyz, g_0_z_yyyzz_xyyyyzz, g_0_z_yyyzz_xyyyzzz, g_0_z_yyyzz_xyyzzzz, g_0_z_yyyzz_xyzzzzz, g_0_z_yyyzz_xzzzzzz, g_0_z_yyyzz_yyyyyyy, g_0_z_yyyzz_yyyyyyz, g_0_z_yyyzz_yyyyyzz, g_0_z_yyyzz_yyyyzzz, g_0_z_yyyzz_yyyzzzz, g_0_z_yyyzz_yyzzzzz, g_0_z_yyyzz_yzzzzzz, g_0_z_yyyzz_zzzzzzz, g_0_z_yyzz_xxxxxxx, g_0_z_yyzz_xxxxxxxy, g_0_z_yyzz_xxxxxxy, g_0_z_yyzz_xxxxxxyy, g_0_z_yyzz_xxxxxxyz, g_0_z_yyzz_xxxxxxz, g_0_z_yyzz_xxxxxyy, g_0_z_yyzz_xxxxxyyy, g_0_z_yyzz_xxxxxyyz, g_0_z_yyzz_xxxxxyz, g_0_z_yyzz_xxxxxyzz, g_0_z_yyzz_xxxxxzz, g_0_z_yyzz_xxxxyyy, g_0_z_yyzz_xxxxyyyy, g_0_z_yyzz_xxxxyyyz, g_0_z_yyzz_xxxxyyz, g_0_z_yyzz_xxxxyyzz, g_0_z_yyzz_xxxxyzz, g_0_z_yyzz_xxxxyzzz, g_0_z_yyzz_xxxxzzz, g_0_z_yyzz_xxxyyyy, g_0_z_yyzz_xxxyyyyy, g_0_z_yyzz_xxxyyyyz, g_0_z_yyzz_xxxyyyz, g_0_z_yyzz_xxxyyyzz, g_0_z_yyzz_xxxyyzz, g_0_z_yyzz_xxxyyzzz, g_0_z_yyzz_xxxyzzz, g_0_z_yyzz_xxxyzzzz, g_0_z_yyzz_xxxzzzz, g_0_z_yyzz_xxyyyyy, g_0_z_yyzz_xxyyyyyy, g_0_z_yyzz_xxyyyyyz, g_0_z_yyzz_xxyyyyz, g_0_z_yyzz_xxyyyyzz, g_0_z_yyzz_xxyyyzz, g_0_z_yyzz_xxyyyzzz, g_0_z_yyzz_xxyyzzz, g_0_z_yyzz_xxyyzzzz, g_0_z_yyzz_xxyzzzz, g_0_z_yyzz_xxyzzzzz, g_0_z_yyzz_xxzzzzz, g_0_z_yyzz_xyyyyyy, g_0_z_yyzz_xyyyyyyy, g_0_z_yyzz_xyyyyyyz, g_0_z_yyzz_xyyyyyz, g_0_z_yyzz_xyyyyyzz, g_0_z_yyzz_xyyyyzz, g_0_z_yyzz_xyyyyzzz, g_0_z_yyzz_xyyyzzz, g_0_z_yyzz_xyyyzzzz, g_0_z_yyzz_xyyzzzz, g_0_z_yyzz_xyyzzzzz, g_0_z_yyzz_xyzzzzz, g_0_z_yyzz_xyzzzzzz, g_0_z_yyzz_xzzzzzz, g_0_z_yyzz_yyyyyyy, g_0_z_yyzz_yyyyyyyy, g_0_z_yyzz_yyyyyyyz, g_0_z_yyzz_yyyyyyz, g_0_z_yyzz_yyyyyyzz, g_0_z_yyzz_yyyyyzz, g_0_z_yyzz_yyyyyzzz, g_0_z_yyzz_yyyyzzz, g_0_z_yyzz_yyyyzzzz, g_0_z_yyzz_yyyzzzz, g_0_z_yyzz_yyyzzzzz, g_0_z_yyzz_yyzzzzz, g_0_z_yyzz_yyzzzzzz, g_0_z_yyzz_yzzzzzz, g_0_z_yyzz_yzzzzzzz, g_0_z_yyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzz_xxxxxxx[k] = -g_0_z_yyzz_xxxxxxx[k] * ab_y + g_0_z_yyzz_xxxxxxxy[k];

                g_0_z_yyyzz_xxxxxxy[k] = -g_0_z_yyzz_xxxxxxy[k] * ab_y + g_0_z_yyzz_xxxxxxyy[k];

                g_0_z_yyyzz_xxxxxxz[k] = -g_0_z_yyzz_xxxxxxz[k] * ab_y + g_0_z_yyzz_xxxxxxyz[k];

                g_0_z_yyyzz_xxxxxyy[k] = -g_0_z_yyzz_xxxxxyy[k] * ab_y + g_0_z_yyzz_xxxxxyyy[k];

                g_0_z_yyyzz_xxxxxyz[k] = -g_0_z_yyzz_xxxxxyz[k] * ab_y + g_0_z_yyzz_xxxxxyyz[k];

                g_0_z_yyyzz_xxxxxzz[k] = -g_0_z_yyzz_xxxxxzz[k] * ab_y + g_0_z_yyzz_xxxxxyzz[k];

                g_0_z_yyyzz_xxxxyyy[k] = -g_0_z_yyzz_xxxxyyy[k] * ab_y + g_0_z_yyzz_xxxxyyyy[k];

                g_0_z_yyyzz_xxxxyyz[k] = -g_0_z_yyzz_xxxxyyz[k] * ab_y + g_0_z_yyzz_xxxxyyyz[k];

                g_0_z_yyyzz_xxxxyzz[k] = -g_0_z_yyzz_xxxxyzz[k] * ab_y + g_0_z_yyzz_xxxxyyzz[k];

                g_0_z_yyyzz_xxxxzzz[k] = -g_0_z_yyzz_xxxxzzz[k] * ab_y + g_0_z_yyzz_xxxxyzzz[k];

                g_0_z_yyyzz_xxxyyyy[k] = -g_0_z_yyzz_xxxyyyy[k] * ab_y + g_0_z_yyzz_xxxyyyyy[k];

                g_0_z_yyyzz_xxxyyyz[k] = -g_0_z_yyzz_xxxyyyz[k] * ab_y + g_0_z_yyzz_xxxyyyyz[k];

                g_0_z_yyyzz_xxxyyzz[k] = -g_0_z_yyzz_xxxyyzz[k] * ab_y + g_0_z_yyzz_xxxyyyzz[k];

                g_0_z_yyyzz_xxxyzzz[k] = -g_0_z_yyzz_xxxyzzz[k] * ab_y + g_0_z_yyzz_xxxyyzzz[k];

                g_0_z_yyyzz_xxxzzzz[k] = -g_0_z_yyzz_xxxzzzz[k] * ab_y + g_0_z_yyzz_xxxyzzzz[k];

                g_0_z_yyyzz_xxyyyyy[k] = -g_0_z_yyzz_xxyyyyy[k] * ab_y + g_0_z_yyzz_xxyyyyyy[k];

                g_0_z_yyyzz_xxyyyyz[k] = -g_0_z_yyzz_xxyyyyz[k] * ab_y + g_0_z_yyzz_xxyyyyyz[k];

                g_0_z_yyyzz_xxyyyzz[k] = -g_0_z_yyzz_xxyyyzz[k] * ab_y + g_0_z_yyzz_xxyyyyzz[k];

                g_0_z_yyyzz_xxyyzzz[k] = -g_0_z_yyzz_xxyyzzz[k] * ab_y + g_0_z_yyzz_xxyyyzzz[k];

                g_0_z_yyyzz_xxyzzzz[k] = -g_0_z_yyzz_xxyzzzz[k] * ab_y + g_0_z_yyzz_xxyyzzzz[k];

                g_0_z_yyyzz_xxzzzzz[k] = -g_0_z_yyzz_xxzzzzz[k] * ab_y + g_0_z_yyzz_xxyzzzzz[k];

                g_0_z_yyyzz_xyyyyyy[k] = -g_0_z_yyzz_xyyyyyy[k] * ab_y + g_0_z_yyzz_xyyyyyyy[k];

                g_0_z_yyyzz_xyyyyyz[k] = -g_0_z_yyzz_xyyyyyz[k] * ab_y + g_0_z_yyzz_xyyyyyyz[k];

                g_0_z_yyyzz_xyyyyzz[k] = -g_0_z_yyzz_xyyyyzz[k] * ab_y + g_0_z_yyzz_xyyyyyzz[k];

                g_0_z_yyyzz_xyyyzzz[k] = -g_0_z_yyzz_xyyyzzz[k] * ab_y + g_0_z_yyzz_xyyyyzzz[k];

                g_0_z_yyyzz_xyyzzzz[k] = -g_0_z_yyzz_xyyzzzz[k] * ab_y + g_0_z_yyzz_xyyyzzzz[k];

                g_0_z_yyyzz_xyzzzzz[k] = -g_0_z_yyzz_xyzzzzz[k] * ab_y + g_0_z_yyzz_xyyzzzzz[k];

                g_0_z_yyyzz_xzzzzzz[k] = -g_0_z_yyzz_xzzzzzz[k] * ab_y + g_0_z_yyzz_xyzzzzzz[k];

                g_0_z_yyyzz_yyyyyyy[k] = -g_0_z_yyzz_yyyyyyy[k] * ab_y + g_0_z_yyzz_yyyyyyyy[k];

                g_0_z_yyyzz_yyyyyyz[k] = -g_0_z_yyzz_yyyyyyz[k] * ab_y + g_0_z_yyzz_yyyyyyyz[k];

                g_0_z_yyyzz_yyyyyzz[k] = -g_0_z_yyzz_yyyyyzz[k] * ab_y + g_0_z_yyzz_yyyyyyzz[k];

                g_0_z_yyyzz_yyyyzzz[k] = -g_0_z_yyzz_yyyyzzz[k] * ab_y + g_0_z_yyzz_yyyyyzzz[k];

                g_0_z_yyyzz_yyyzzzz[k] = -g_0_z_yyzz_yyyzzzz[k] * ab_y + g_0_z_yyzz_yyyyzzzz[k];

                g_0_z_yyyzz_yyzzzzz[k] = -g_0_z_yyzz_yyzzzzz[k] * ab_y + g_0_z_yyzz_yyyzzzzz[k];

                g_0_z_yyyzz_yzzzzzz[k] = -g_0_z_yyzz_yzzzzzz[k] * ab_y + g_0_z_yyzz_yyzzzzzz[k];

                g_0_z_yyyzz_zzzzzzz[k] = -g_0_z_yyzz_zzzzzzz[k] * ab_y + g_0_z_yyzz_yzzzzzzz[k];
            }

            /// Set up 2160-2196 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2160 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2161 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2162 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2163 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2164 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2165 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2166 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2167 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2168 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2169 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2170 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2171 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2172 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2173 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2174 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2175 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2176 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2177 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2178 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2179 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2180 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2181 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2182 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2183 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2184 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2185 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2186 * ccomps * dcomps);

            auto g_0_z_yyzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2187 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2188 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2189 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2190 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2191 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2192 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2193 * ccomps * dcomps);

            auto g_0_z_yyzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2194 * ccomps * dcomps);

            auto g_0_z_yyzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzz_xxxxxxx, g_0_z_yyzzz_xxxxxxy, g_0_z_yyzzz_xxxxxxz, g_0_z_yyzzz_xxxxxyy, g_0_z_yyzzz_xxxxxyz, g_0_z_yyzzz_xxxxxzz, g_0_z_yyzzz_xxxxyyy, g_0_z_yyzzz_xxxxyyz, g_0_z_yyzzz_xxxxyzz, g_0_z_yyzzz_xxxxzzz, g_0_z_yyzzz_xxxyyyy, g_0_z_yyzzz_xxxyyyz, g_0_z_yyzzz_xxxyyzz, g_0_z_yyzzz_xxxyzzz, g_0_z_yyzzz_xxxzzzz, g_0_z_yyzzz_xxyyyyy, g_0_z_yyzzz_xxyyyyz, g_0_z_yyzzz_xxyyyzz, g_0_z_yyzzz_xxyyzzz, g_0_z_yyzzz_xxyzzzz, g_0_z_yyzzz_xxzzzzz, g_0_z_yyzzz_xyyyyyy, g_0_z_yyzzz_xyyyyyz, g_0_z_yyzzz_xyyyyzz, g_0_z_yyzzz_xyyyzzz, g_0_z_yyzzz_xyyzzzz, g_0_z_yyzzz_xyzzzzz, g_0_z_yyzzz_xzzzzzz, g_0_z_yyzzz_yyyyyyy, g_0_z_yyzzz_yyyyyyz, g_0_z_yyzzz_yyyyyzz, g_0_z_yyzzz_yyyyzzz, g_0_z_yyzzz_yyyzzzz, g_0_z_yyzzz_yyzzzzz, g_0_z_yyzzz_yzzzzzz, g_0_z_yyzzz_zzzzzzz, g_0_z_yzzz_xxxxxxx, g_0_z_yzzz_xxxxxxxy, g_0_z_yzzz_xxxxxxy, g_0_z_yzzz_xxxxxxyy, g_0_z_yzzz_xxxxxxyz, g_0_z_yzzz_xxxxxxz, g_0_z_yzzz_xxxxxyy, g_0_z_yzzz_xxxxxyyy, g_0_z_yzzz_xxxxxyyz, g_0_z_yzzz_xxxxxyz, g_0_z_yzzz_xxxxxyzz, g_0_z_yzzz_xxxxxzz, g_0_z_yzzz_xxxxyyy, g_0_z_yzzz_xxxxyyyy, g_0_z_yzzz_xxxxyyyz, g_0_z_yzzz_xxxxyyz, g_0_z_yzzz_xxxxyyzz, g_0_z_yzzz_xxxxyzz, g_0_z_yzzz_xxxxyzzz, g_0_z_yzzz_xxxxzzz, g_0_z_yzzz_xxxyyyy, g_0_z_yzzz_xxxyyyyy, g_0_z_yzzz_xxxyyyyz, g_0_z_yzzz_xxxyyyz, g_0_z_yzzz_xxxyyyzz, g_0_z_yzzz_xxxyyzz, g_0_z_yzzz_xxxyyzzz, g_0_z_yzzz_xxxyzzz, g_0_z_yzzz_xxxyzzzz, g_0_z_yzzz_xxxzzzz, g_0_z_yzzz_xxyyyyy, g_0_z_yzzz_xxyyyyyy, g_0_z_yzzz_xxyyyyyz, g_0_z_yzzz_xxyyyyz, g_0_z_yzzz_xxyyyyzz, g_0_z_yzzz_xxyyyzz, g_0_z_yzzz_xxyyyzzz, g_0_z_yzzz_xxyyzzz, g_0_z_yzzz_xxyyzzzz, g_0_z_yzzz_xxyzzzz, g_0_z_yzzz_xxyzzzzz, g_0_z_yzzz_xxzzzzz, g_0_z_yzzz_xyyyyyy, g_0_z_yzzz_xyyyyyyy, g_0_z_yzzz_xyyyyyyz, g_0_z_yzzz_xyyyyyz, g_0_z_yzzz_xyyyyyzz, g_0_z_yzzz_xyyyyzz, g_0_z_yzzz_xyyyyzzz, g_0_z_yzzz_xyyyzzz, g_0_z_yzzz_xyyyzzzz, g_0_z_yzzz_xyyzzzz, g_0_z_yzzz_xyyzzzzz, g_0_z_yzzz_xyzzzzz, g_0_z_yzzz_xyzzzzzz, g_0_z_yzzz_xzzzzzz, g_0_z_yzzz_yyyyyyy, g_0_z_yzzz_yyyyyyyy, g_0_z_yzzz_yyyyyyyz, g_0_z_yzzz_yyyyyyz, g_0_z_yzzz_yyyyyyzz, g_0_z_yzzz_yyyyyzz, g_0_z_yzzz_yyyyyzzz, g_0_z_yzzz_yyyyzzz, g_0_z_yzzz_yyyyzzzz, g_0_z_yzzz_yyyzzzz, g_0_z_yzzz_yyyzzzzz, g_0_z_yzzz_yyzzzzz, g_0_z_yzzz_yyzzzzzz, g_0_z_yzzz_yzzzzzz, g_0_z_yzzz_yzzzzzzz, g_0_z_yzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzz_xxxxxxx[k] = -g_0_z_yzzz_xxxxxxx[k] * ab_y + g_0_z_yzzz_xxxxxxxy[k];

                g_0_z_yyzzz_xxxxxxy[k] = -g_0_z_yzzz_xxxxxxy[k] * ab_y + g_0_z_yzzz_xxxxxxyy[k];

                g_0_z_yyzzz_xxxxxxz[k] = -g_0_z_yzzz_xxxxxxz[k] * ab_y + g_0_z_yzzz_xxxxxxyz[k];

                g_0_z_yyzzz_xxxxxyy[k] = -g_0_z_yzzz_xxxxxyy[k] * ab_y + g_0_z_yzzz_xxxxxyyy[k];

                g_0_z_yyzzz_xxxxxyz[k] = -g_0_z_yzzz_xxxxxyz[k] * ab_y + g_0_z_yzzz_xxxxxyyz[k];

                g_0_z_yyzzz_xxxxxzz[k] = -g_0_z_yzzz_xxxxxzz[k] * ab_y + g_0_z_yzzz_xxxxxyzz[k];

                g_0_z_yyzzz_xxxxyyy[k] = -g_0_z_yzzz_xxxxyyy[k] * ab_y + g_0_z_yzzz_xxxxyyyy[k];

                g_0_z_yyzzz_xxxxyyz[k] = -g_0_z_yzzz_xxxxyyz[k] * ab_y + g_0_z_yzzz_xxxxyyyz[k];

                g_0_z_yyzzz_xxxxyzz[k] = -g_0_z_yzzz_xxxxyzz[k] * ab_y + g_0_z_yzzz_xxxxyyzz[k];

                g_0_z_yyzzz_xxxxzzz[k] = -g_0_z_yzzz_xxxxzzz[k] * ab_y + g_0_z_yzzz_xxxxyzzz[k];

                g_0_z_yyzzz_xxxyyyy[k] = -g_0_z_yzzz_xxxyyyy[k] * ab_y + g_0_z_yzzz_xxxyyyyy[k];

                g_0_z_yyzzz_xxxyyyz[k] = -g_0_z_yzzz_xxxyyyz[k] * ab_y + g_0_z_yzzz_xxxyyyyz[k];

                g_0_z_yyzzz_xxxyyzz[k] = -g_0_z_yzzz_xxxyyzz[k] * ab_y + g_0_z_yzzz_xxxyyyzz[k];

                g_0_z_yyzzz_xxxyzzz[k] = -g_0_z_yzzz_xxxyzzz[k] * ab_y + g_0_z_yzzz_xxxyyzzz[k];

                g_0_z_yyzzz_xxxzzzz[k] = -g_0_z_yzzz_xxxzzzz[k] * ab_y + g_0_z_yzzz_xxxyzzzz[k];

                g_0_z_yyzzz_xxyyyyy[k] = -g_0_z_yzzz_xxyyyyy[k] * ab_y + g_0_z_yzzz_xxyyyyyy[k];

                g_0_z_yyzzz_xxyyyyz[k] = -g_0_z_yzzz_xxyyyyz[k] * ab_y + g_0_z_yzzz_xxyyyyyz[k];

                g_0_z_yyzzz_xxyyyzz[k] = -g_0_z_yzzz_xxyyyzz[k] * ab_y + g_0_z_yzzz_xxyyyyzz[k];

                g_0_z_yyzzz_xxyyzzz[k] = -g_0_z_yzzz_xxyyzzz[k] * ab_y + g_0_z_yzzz_xxyyyzzz[k];

                g_0_z_yyzzz_xxyzzzz[k] = -g_0_z_yzzz_xxyzzzz[k] * ab_y + g_0_z_yzzz_xxyyzzzz[k];

                g_0_z_yyzzz_xxzzzzz[k] = -g_0_z_yzzz_xxzzzzz[k] * ab_y + g_0_z_yzzz_xxyzzzzz[k];

                g_0_z_yyzzz_xyyyyyy[k] = -g_0_z_yzzz_xyyyyyy[k] * ab_y + g_0_z_yzzz_xyyyyyyy[k];

                g_0_z_yyzzz_xyyyyyz[k] = -g_0_z_yzzz_xyyyyyz[k] * ab_y + g_0_z_yzzz_xyyyyyyz[k];

                g_0_z_yyzzz_xyyyyzz[k] = -g_0_z_yzzz_xyyyyzz[k] * ab_y + g_0_z_yzzz_xyyyyyzz[k];

                g_0_z_yyzzz_xyyyzzz[k] = -g_0_z_yzzz_xyyyzzz[k] * ab_y + g_0_z_yzzz_xyyyyzzz[k];

                g_0_z_yyzzz_xyyzzzz[k] = -g_0_z_yzzz_xyyzzzz[k] * ab_y + g_0_z_yzzz_xyyyzzzz[k];

                g_0_z_yyzzz_xyzzzzz[k] = -g_0_z_yzzz_xyzzzzz[k] * ab_y + g_0_z_yzzz_xyyzzzzz[k];

                g_0_z_yyzzz_xzzzzzz[k] = -g_0_z_yzzz_xzzzzzz[k] * ab_y + g_0_z_yzzz_xyzzzzzz[k];

                g_0_z_yyzzz_yyyyyyy[k] = -g_0_z_yzzz_yyyyyyy[k] * ab_y + g_0_z_yzzz_yyyyyyyy[k];

                g_0_z_yyzzz_yyyyyyz[k] = -g_0_z_yzzz_yyyyyyz[k] * ab_y + g_0_z_yzzz_yyyyyyyz[k];

                g_0_z_yyzzz_yyyyyzz[k] = -g_0_z_yzzz_yyyyyzz[k] * ab_y + g_0_z_yzzz_yyyyyyzz[k];

                g_0_z_yyzzz_yyyyzzz[k] = -g_0_z_yzzz_yyyyzzz[k] * ab_y + g_0_z_yzzz_yyyyyzzz[k];

                g_0_z_yyzzz_yyyzzzz[k] = -g_0_z_yzzz_yyyzzzz[k] * ab_y + g_0_z_yzzz_yyyyzzzz[k];

                g_0_z_yyzzz_yyzzzzz[k] = -g_0_z_yzzz_yyzzzzz[k] * ab_y + g_0_z_yzzz_yyyzzzzz[k];

                g_0_z_yyzzz_yzzzzzz[k] = -g_0_z_yzzz_yzzzzzz[k] * ab_y + g_0_z_yzzz_yyzzzzzz[k];

                g_0_z_yyzzz_zzzzzzz[k] = -g_0_z_yzzz_zzzzzzz[k] * ab_y + g_0_z_yzzz_yzzzzzzz[k];
            }

            /// Set up 2196-2232 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2196 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2197 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2198 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2199 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2200 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2201 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2202 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2203 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2204 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2205 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2206 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2207 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2208 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2209 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2210 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2211 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2212 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2213 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2214 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2215 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2216 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2217 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2218 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2219 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2220 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2221 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2222 * ccomps * dcomps);

            auto g_0_z_yzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2223 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2224 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2225 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2226 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2227 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2228 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2229 * ccomps * dcomps);

            auto g_0_z_yzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2230 * ccomps * dcomps);

            auto g_0_z_yzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2231 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzz_xxxxxxx, g_0_z_yzzzz_xxxxxxy, g_0_z_yzzzz_xxxxxxz, g_0_z_yzzzz_xxxxxyy, g_0_z_yzzzz_xxxxxyz, g_0_z_yzzzz_xxxxxzz, g_0_z_yzzzz_xxxxyyy, g_0_z_yzzzz_xxxxyyz, g_0_z_yzzzz_xxxxyzz, g_0_z_yzzzz_xxxxzzz, g_0_z_yzzzz_xxxyyyy, g_0_z_yzzzz_xxxyyyz, g_0_z_yzzzz_xxxyyzz, g_0_z_yzzzz_xxxyzzz, g_0_z_yzzzz_xxxzzzz, g_0_z_yzzzz_xxyyyyy, g_0_z_yzzzz_xxyyyyz, g_0_z_yzzzz_xxyyyzz, g_0_z_yzzzz_xxyyzzz, g_0_z_yzzzz_xxyzzzz, g_0_z_yzzzz_xxzzzzz, g_0_z_yzzzz_xyyyyyy, g_0_z_yzzzz_xyyyyyz, g_0_z_yzzzz_xyyyyzz, g_0_z_yzzzz_xyyyzzz, g_0_z_yzzzz_xyyzzzz, g_0_z_yzzzz_xyzzzzz, g_0_z_yzzzz_xzzzzzz, g_0_z_yzzzz_yyyyyyy, g_0_z_yzzzz_yyyyyyz, g_0_z_yzzzz_yyyyyzz, g_0_z_yzzzz_yyyyzzz, g_0_z_yzzzz_yyyzzzz, g_0_z_yzzzz_yyzzzzz, g_0_z_yzzzz_yzzzzzz, g_0_z_yzzzz_zzzzzzz, g_0_z_zzzz_xxxxxxx, g_0_z_zzzz_xxxxxxxy, g_0_z_zzzz_xxxxxxy, g_0_z_zzzz_xxxxxxyy, g_0_z_zzzz_xxxxxxyz, g_0_z_zzzz_xxxxxxz, g_0_z_zzzz_xxxxxyy, g_0_z_zzzz_xxxxxyyy, g_0_z_zzzz_xxxxxyyz, g_0_z_zzzz_xxxxxyz, g_0_z_zzzz_xxxxxyzz, g_0_z_zzzz_xxxxxzz, g_0_z_zzzz_xxxxyyy, g_0_z_zzzz_xxxxyyyy, g_0_z_zzzz_xxxxyyyz, g_0_z_zzzz_xxxxyyz, g_0_z_zzzz_xxxxyyzz, g_0_z_zzzz_xxxxyzz, g_0_z_zzzz_xxxxyzzz, g_0_z_zzzz_xxxxzzz, g_0_z_zzzz_xxxyyyy, g_0_z_zzzz_xxxyyyyy, g_0_z_zzzz_xxxyyyyz, g_0_z_zzzz_xxxyyyz, g_0_z_zzzz_xxxyyyzz, g_0_z_zzzz_xxxyyzz, g_0_z_zzzz_xxxyyzzz, g_0_z_zzzz_xxxyzzz, g_0_z_zzzz_xxxyzzzz, g_0_z_zzzz_xxxzzzz, g_0_z_zzzz_xxyyyyy, g_0_z_zzzz_xxyyyyyy, g_0_z_zzzz_xxyyyyyz, g_0_z_zzzz_xxyyyyz, g_0_z_zzzz_xxyyyyzz, g_0_z_zzzz_xxyyyzz, g_0_z_zzzz_xxyyyzzz, g_0_z_zzzz_xxyyzzz, g_0_z_zzzz_xxyyzzzz, g_0_z_zzzz_xxyzzzz, g_0_z_zzzz_xxyzzzzz, g_0_z_zzzz_xxzzzzz, g_0_z_zzzz_xyyyyyy, g_0_z_zzzz_xyyyyyyy, g_0_z_zzzz_xyyyyyyz, g_0_z_zzzz_xyyyyyz, g_0_z_zzzz_xyyyyyzz, g_0_z_zzzz_xyyyyzz, g_0_z_zzzz_xyyyyzzz, g_0_z_zzzz_xyyyzzz, g_0_z_zzzz_xyyyzzzz, g_0_z_zzzz_xyyzzzz, g_0_z_zzzz_xyyzzzzz, g_0_z_zzzz_xyzzzzz, g_0_z_zzzz_xyzzzzzz, g_0_z_zzzz_xzzzzzz, g_0_z_zzzz_yyyyyyy, g_0_z_zzzz_yyyyyyyy, g_0_z_zzzz_yyyyyyyz, g_0_z_zzzz_yyyyyyz, g_0_z_zzzz_yyyyyyzz, g_0_z_zzzz_yyyyyzz, g_0_z_zzzz_yyyyyzzz, g_0_z_zzzz_yyyyzzz, g_0_z_zzzz_yyyyzzzz, g_0_z_zzzz_yyyzzzz, g_0_z_zzzz_yyyzzzzz, g_0_z_zzzz_yyzzzzz, g_0_z_zzzz_yyzzzzzz, g_0_z_zzzz_yzzzzzz, g_0_z_zzzz_yzzzzzzz, g_0_z_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzz_xxxxxxx[k] = -g_0_z_zzzz_xxxxxxx[k] * ab_y + g_0_z_zzzz_xxxxxxxy[k];

                g_0_z_yzzzz_xxxxxxy[k] = -g_0_z_zzzz_xxxxxxy[k] * ab_y + g_0_z_zzzz_xxxxxxyy[k];

                g_0_z_yzzzz_xxxxxxz[k] = -g_0_z_zzzz_xxxxxxz[k] * ab_y + g_0_z_zzzz_xxxxxxyz[k];

                g_0_z_yzzzz_xxxxxyy[k] = -g_0_z_zzzz_xxxxxyy[k] * ab_y + g_0_z_zzzz_xxxxxyyy[k];

                g_0_z_yzzzz_xxxxxyz[k] = -g_0_z_zzzz_xxxxxyz[k] * ab_y + g_0_z_zzzz_xxxxxyyz[k];

                g_0_z_yzzzz_xxxxxzz[k] = -g_0_z_zzzz_xxxxxzz[k] * ab_y + g_0_z_zzzz_xxxxxyzz[k];

                g_0_z_yzzzz_xxxxyyy[k] = -g_0_z_zzzz_xxxxyyy[k] * ab_y + g_0_z_zzzz_xxxxyyyy[k];

                g_0_z_yzzzz_xxxxyyz[k] = -g_0_z_zzzz_xxxxyyz[k] * ab_y + g_0_z_zzzz_xxxxyyyz[k];

                g_0_z_yzzzz_xxxxyzz[k] = -g_0_z_zzzz_xxxxyzz[k] * ab_y + g_0_z_zzzz_xxxxyyzz[k];

                g_0_z_yzzzz_xxxxzzz[k] = -g_0_z_zzzz_xxxxzzz[k] * ab_y + g_0_z_zzzz_xxxxyzzz[k];

                g_0_z_yzzzz_xxxyyyy[k] = -g_0_z_zzzz_xxxyyyy[k] * ab_y + g_0_z_zzzz_xxxyyyyy[k];

                g_0_z_yzzzz_xxxyyyz[k] = -g_0_z_zzzz_xxxyyyz[k] * ab_y + g_0_z_zzzz_xxxyyyyz[k];

                g_0_z_yzzzz_xxxyyzz[k] = -g_0_z_zzzz_xxxyyzz[k] * ab_y + g_0_z_zzzz_xxxyyyzz[k];

                g_0_z_yzzzz_xxxyzzz[k] = -g_0_z_zzzz_xxxyzzz[k] * ab_y + g_0_z_zzzz_xxxyyzzz[k];

                g_0_z_yzzzz_xxxzzzz[k] = -g_0_z_zzzz_xxxzzzz[k] * ab_y + g_0_z_zzzz_xxxyzzzz[k];

                g_0_z_yzzzz_xxyyyyy[k] = -g_0_z_zzzz_xxyyyyy[k] * ab_y + g_0_z_zzzz_xxyyyyyy[k];

                g_0_z_yzzzz_xxyyyyz[k] = -g_0_z_zzzz_xxyyyyz[k] * ab_y + g_0_z_zzzz_xxyyyyyz[k];

                g_0_z_yzzzz_xxyyyzz[k] = -g_0_z_zzzz_xxyyyzz[k] * ab_y + g_0_z_zzzz_xxyyyyzz[k];

                g_0_z_yzzzz_xxyyzzz[k] = -g_0_z_zzzz_xxyyzzz[k] * ab_y + g_0_z_zzzz_xxyyyzzz[k];

                g_0_z_yzzzz_xxyzzzz[k] = -g_0_z_zzzz_xxyzzzz[k] * ab_y + g_0_z_zzzz_xxyyzzzz[k];

                g_0_z_yzzzz_xxzzzzz[k] = -g_0_z_zzzz_xxzzzzz[k] * ab_y + g_0_z_zzzz_xxyzzzzz[k];

                g_0_z_yzzzz_xyyyyyy[k] = -g_0_z_zzzz_xyyyyyy[k] * ab_y + g_0_z_zzzz_xyyyyyyy[k];

                g_0_z_yzzzz_xyyyyyz[k] = -g_0_z_zzzz_xyyyyyz[k] * ab_y + g_0_z_zzzz_xyyyyyyz[k];

                g_0_z_yzzzz_xyyyyzz[k] = -g_0_z_zzzz_xyyyyzz[k] * ab_y + g_0_z_zzzz_xyyyyyzz[k];

                g_0_z_yzzzz_xyyyzzz[k] = -g_0_z_zzzz_xyyyzzz[k] * ab_y + g_0_z_zzzz_xyyyyzzz[k];

                g_0_z_yzzzz_xyyzzzz[k] = -g_0_z_zzzz_xyyzzzz[k] * ab_y + g_0_z_zzzz_xyyyzzzz[k];

                g_0_z_yzzzz_xyzzzzz[k] = -g_0_z_zzzz_xyzzzzz[k] * ab_y + g_0_z_zzzz_xyyzzzzz[k];

                g_0_z_yzzzz_xzzzzzz[k] = -g_0_z_zzzz_xzzzzzz[k] * ab_y + g_0_z_zzzz_xyzzzzzz[k];

                g_0_z_yzzzz_yyyyyyy[k] = -g_0_z_zzzz_yyyyyyy[k] * ab_y + g_0_z_zzzz_yyyyyyyy[k];

                g_0_z_yzzzz_yyyyyyz[k] = -g_0_z_zzzz_yyyyyyz[k] * ab_y + g_0_z_zzzz_yyyyyyyz[k];

                g_0_z_yzzzz_yyyyyzz[k] = -g_0_z_zzzz_yyyyyzz[k] * ab_y + g_0_z_zzzz_yyyyyyzz[k];

                g_0_z_yzzzz_yyyyzzz[k] = -g_0_z_zzzz_yyyyzzz[k] * ab_y + g_0_z_zzzz_yyyyyzzz[k];

                g_0_z_yzzzz_yyyzzzz[k] = -g_0_z_zzzz_yyyzzzz[k] * ab_y + g_0_z_zzzz_yyyyzzzz[k];

                g_0_z_yzzzz_yyzzzzz[k] = -g_0_z_zzzz_yyzzzzz[k] * ab_y + g_0_z_zzzz_yyyzzzzz[k];

                g_0_z_yzzzz_yzzzzzz[k] = -g_0_z_zzzz_yzzzzzz[k] * ab_y + g_0_z_zzzz_yyzzzzzz[k];

                g_0_z_yzzzz_zzzzzzz[k] = -g_0_z_zzzz_zzzzzzz[k] * ab_y + g_0_z_zzzz_yzzzzzzz[k];
            }

            /// Set up 2232-2268 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzz_xxxxxxx = cbuffer.data(hk_geom_01_off + 2232 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxxxy = cbuffer.data(hk_geom_01_off + 2233 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxxxz = cbuffer.data(hk_geom_01_off + 2234 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxxyy = cbuffer.data(hk_geom_01_off + 2235 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxxyz = cbuffer.data(hk_geom_01_off + 2236 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxxzz = cbuffer.data(hk_geom_01_off + 2237 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxyyy = cbuffer.data(hk_geom_01_off + 2238 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxyyz = cbuffer.data(hk_geom_01_off + 2239 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxyzz = cbuffer.data(hk_geom_01_off + 2240 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxzzz = cbuffer.data(hk_geom_01_off + 2241 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyyyy = cbuffer.data(hk_geom_01_off + 2242 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyyyz = cbuffer.data(hk_geom_01_off + 2243 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyyzz = cbuffer.data(hk_geom_01_off + 2244 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyzzz = cbuffer.data(hk_geom_01_off + 2245 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxzzzz = cbuffer.data(hk_geom_01_off + 2246 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyyyy = cbuffer.data(hk_geom_01_off + 2247 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyyyz = cbuffer.data(hk_geom_01_off + 2248 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyyzz = cbuffer.data(hk_geom_01_off + 2249 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyzzz = cbuffer.data(hk_geom_01_off + 2250 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyzzzz = cbuffer.data(hk_geom_01_off + 2251 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxzzzzz = cbuffer.data(hk_geom_01_off + 2252 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyyyy = cbuffer.data(hk_geom_01_off + 2253 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyyyz = cbuffer.data(hk_geom_01_off + 2254 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyyzz = cbuffer.data(hk_geom_01_off + 2255 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyzzz = cbuffer.data(hk_geom_01_off + 2256 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyzzzz = cbuffer.data(hk_geom_01_off + 2257 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyzzzzz = cbuffer.data(hk_geom_01_off + 2258 * ccomps * dcomps);

            auto g_0_z_zzzzz_xzzzzzz = cbuffer.data(hk_geom_01_off + 2259 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyyyy = cbuffer.data(hk_geom_01_off + 2260 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyyyz = cbuffer.data(hk_geom_01_off + 2261 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyyzz = cbuffer.data(hk_geom_01_off + 2262 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyzzz = cbuffer.data(hk_geom_01_off + 2263 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyzzzz = cbuffer.data(hk_geom_01_off + 2264 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyzzzzz = cbuffer.data(hk_geom_01_off + 2265 * ccomps * dcomps);

            auto g_0_z_zzzzz_yzzzzzz = cbuffer.data(hk_geom_01_off + 2266 * ccomps * dcomps);

            auto g_0_z_zzzzz_zzzzzzz = cbuffer.data(hk_geom_01_off + 2267 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xxxxxxx, g_0_z_zzzz_xxxxxxxz, g_0_z_zzzz_xxxxxxy, g_0_z_zzzz_xxxxxxyz, g_0_z_zzzz_xxxxxxz, g_0_z_zzzz_xxxxxxzz, g_0_z_zzzz_xxxxxyy, g_0_z_zzzz_xxxxxyyz, g_0_z_zzzz_xxxxxyz, g_0_z_zzzz_xxxxxyzz, g_0_z_zzzz_xxxxxzz, g_0_z_zzzz_xxxxxzzz, g_0_z_zzzz_xxxxyyy, g_0_z_zzzz_xxxxyyyz, g_0_z_zzzz_xxxxyyz, g_0_z_zzzz_xxxxyyzz, g_0_z_zzzz_xxxxyzz, g_0_z_zzzz_xxxxyzzz, g_0_z_zzzz_xxxxzzz, g_0_z_zzzz_xxxxzzzz, g_0_z_zzzz_xxxyyyy, g_0_z_zzzz_xxxyyyyz, g_0_z_zzzz_xxxyyyz, g_0_z_zzzz_xxxyyyzz, g_0_z_zzzz_xxxyyzz, g_0_z_zzzz_xxxyyzzz, g_0_z_zzzz_xxxyzzz, g_0_z_zzzz_xxxyzzzz, g_0_z_zzzz_xxxzzzz, g_0_z_zzzz_xxxzzzzz, g_0_z_zzzz_xxyyyyy, g_0_z_zzzz_xxyyyyyz, g_0_z_zzzz_xxyyyyz, g_0_z_zzzz_xxyyyyzz, g_0_z_zzzz_xxyyyzz, g_0_z_zzzz_xxyyyzzz, g_0_z_zzzz_xxyyzzz, g_0_z_zzzz_xxyyzzzz, g_0_z_zzzz_xxyzzzz, g_0_z_zzzz_xxyzzzzz, g_0_z_zzzz_xxzzzzz, g_0_z_zzzz_xxzzzzzz, g_0_z_zzzz_xyyyyyy, g_0_z_zzzz_xyyyyyyz, g_0_z_zzzz_xyyyyyz, g_0_z_zzzz_xyyyyyzz, g_0_z_zzzz_xyyyyzz, g_0_z_zzzz_xyyyyzzz, g_0_z_zzzz_xyyyzzz, g_0_z_zzzz_xyyyzzzz, g_0_z_zzzz_xyyzzzz, g_0_z_zzzz_xyyzzzzz, g_0_z_zzzz_xyzzzzz, g_0_z_zzzz_xyzzzzzz, g_0_z_zzzz_xzzzzzz, g_0_z_zzzz_xzzzzzzz, g_0_z_zzzz_yyyyyyy, g_0_z_zzzz_yyyyyyyz, g_0_z_zzzz_yyyyyyz, g_0_z_zzzz_yyyyyyzz, g_0_z_zzzz_yyyyyzz, g_0_z_zzzz_yyyyyzzz, g_0_z_zzzz_yyyyzzz, g_0_z_zzzz_yyyyzzzz, g_0_z_zzzz_yyyzzzz, g_0_z_zzzz_yyyzzzzz, g_0_z_zzzz_yyzzzzz, g_0_z_zzzz_yyzzzzzz, g_0_z_zzzz_yzzzzzz, g_0_z_zzzz_yzzzzzzz, g_0_z_zzzz_zzzzzzz, g_0_z_zzzz_zzzzzzzz, g_0_z_zzzzz_xxxxxxx, g_0_z_zzzzz_xxxxxxy, g_0_z_zzzzz_xxxxxxz, g_0_z_zzzzz_xxxxxyy, g_0_z_zzzzz_xxxxxyz, g_0_z_zzzzz_xxxxxzz, g_0_z_zzzzz_xxxxyyy, g_0_z_zzzzz_xxxxyyz, g_0_z_zzzzz_xxxxyzz, g_0_z_zzzzz_xxxxzzz, g_0_z_zzzzz_xxxyyyy, g_0_z_zzzzz_xxxyyyz, g_0_z_zzzzz_xxxyyzz, g_0_z_zzzzz_xxxyzzz, g_0_z_zzzzz_xxxzzzz, g_0_z_zzzzz_xxyyyyy, g_0_z_zzzzz_xxyyyyz, g_0_z_zzzzz_xxyyyzz, g_0_z_zzzzz_xxyyzzz, g_0_z_zzzzz_xxyzzzz, g_0_z_zzzzz_xxzzzzz, g_0_z_zzzzz_xyyyyyy, g_0_z_zzzzz_xyyyyyz, g_0_z_zzzzz_xyyyyzz, g_0_z_zzzzz_xyyyzzz, g_0_z_zzzzz_xyyzzzz, g_0_z_zzzzz_xyzzzzz, g_0_z_zzzzz_xzzzzzz, g_0_z_zzzzz_yyyyyyy, g_0_z_zzzzz_yyyyyyz, g_0_z_zzzzz_yyyyyzz, g_0_z_zzzzz_yyyyzzz, g_0_z_zzzzz_yyyzzzz, g_0_z_zzzzz_yyzzzzz, g_0_z_zzzzz_yzzzzzz, g_0_z_zzzzz_zzzzzzz, g_zzzz_xxxxxxx, g_zzzz_xxxxxxy, g_zzzz_xxxxxxz, g_zzzz_xxxxxyy, g_zzzz_xxxxxyz, g_zzzz_xxxxxzz, g_zzzz_xxxxyyy, g_zzzz_xxxxyyz, g_zzzz_xxxxyzz, g_zzzz_xxxxzzz, g_zzzz_xxxyyyy, g_zzzz_xxxyyyz, g_zzzz_xxxyyzz, g_zzzz_xxxyzzz, g_zzzz_xxxzzzz, g_zzzz_xxyyyyy, g_zzzz_xxyyyyz, g_zzzz_xxyyyzz, g_zzzz_xxyyzzz, g_zzzz_xxyzzzz, g_zzzz_xxzzzzz, g_zzzz_xyyyyyy, g_zzzz_xyyyyyz, g_zzzz_xyyyyzz, g_zzzz_xyyyzzz, g_zzzz_xyyzzzz, g_zzzz_xyzzzzz, g_zzzz_xzzzzzz, g_zzzz_yyyyyyy, g_zzzz_yyyyyyz, g_zzzz_yyyyyzz, g_zzzz_yyyyzzz, g_zzzz_yyyzzzz, g_zzzz_yyzzzzz, g_zzzz_yzzzzzz, g_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzz_xxxxxxx[k] = g_zzzz_xxxxxxx[k] - g_0_z_zzzz_xxxxxxx[k] * ab_z + g_0_z_zzzz_xxxxxxxz[k];

                g_0_z_zzzzz_xxxxxxy[k] = g_zzzz_xxxxxxy[k] - g_0_z_zzzz_xxxxxxy[k] * ab_z + g_0_z_zzzz_xxxxxxyz[k];

                g_0_z_zzzzz_xxxxxxz[k] = g_zzzz_xxxxxxz[k] - g_0_z_zzzz_xxxxxxz[k] * ab_z + g_0_z_zzzz_xxxxxxzz[k];

                g_0_z_zzzzz_xxxxxyy[k] = g_zzzz_xxxxxyy[k] - g_0_z_zzzz_xxxxxyy[k] * ab_z + g_0_z_zzzz_xxxxxyyz[k];

                g_0_z_zzzzz_xxxxxyz[k] = g_zzzz_xxxxxyz[k] - g_0_z_zzzz_xxxxxyz[k] * ab_z + g_0_z_zzzz_xxxxxyzz[k];

                g_0_z_zzzzz_xxxxxzz[k] = g_zzzz_xxxxxzz[k] - g_0_z_zzzz_xxxxxzz[k] * ab_z + g_0_z_zzzz_xxxxxzzz[k];

                g_0_z_zzzzz_xxxxyyy[k] = g_zzzz_xxxxyyy[k] - g_0_z_zzzz_xxxxyyy[k] * ab_z + g_0_z_zzzz_xxxxyyyz[k];

                g_0_z_zzzzz_xxxxyyz[k] = g_zzzz_xxxxyyz[k] - g_0_z_zzzz_xxxxyyz[k] * ab_z + g_0_z_zzzz_xxxxyyzz[k];

                g_0_z_zzzzz_xxxxyzz[k] = g_zzzz_xxxxyzz[k] - g_0_z_zzzz_xxxxyzz[k] * ab_z + g_0_z_zzzz_xxxxyzzz[k];

                g_0_z_zzzzz_xxxxzzz[k] = g_zzzz_xxxxzzz[k] - g_0_z_zzzz_xxxxzzz[k] * ab_z + g_0_z_zzzz_xxxxzzzz[k];

                g_0_z_zzzzz_xxxyyyy[k] = g_zzzz_xxxyyyy[k] - g_0_z_zzzz_xxxyyyy[k] * ab_z + g_0_z_zzzz_xxxyyyyz[k];

                g_0_z_zzzzz_xxxyyyz[k] = g_zzzz_xxxyyyz[k] - g_0_z_zzzz_xxxyyyz[k] * ab_z + g_0_z_zzzz_xxxyyyzz[k];

                g_0_z_zzzzz_xxxyyzz[k] = g_zzzz_xxxyyzz[k] - g_0_z_zzzz_xxxyyzz[k] * ab_z + g_0_z_zzzz_xxxyyzzz[k];

                g_0_z_zzzzz_xxxyzzz[k] = g_zzzz_xxxyzzz[k] - g_0_z_zzzz_xxxyzzz[k] * ab_z + g_0_z_zzzz_xxxyzzzz[k];

                g_0_z_zzzzz_xxxzzzz[k] = g_zzzz_xxxzzzz[k] - g_0_z_zzzz_xxxzzzz[k] * ab_z + g_0_z_zzzz_xxxzzzzz[k];

                g_0_z_zzzzz_xxyyyyy[k] = g_zzzz_xxyyyyy[k] - g_0_z_zzzz_xxyyyyy[k] * ab_z + g_0_z_zzzz_xxyyyyyz[k];

                g_0_z_zzzzz_xxyyyyz[k] = g_zzzz_xxyyyyz[k] - g_0_z_zzzz_xxyyyyz[k] * ab_z + g_0_z_zzzz_xxyyyyzz[k];

                g_0_z_zzzzz_xxyyyzz[k] = g_zzzz_xxyyyzz[k] - g_0_z_zzzz_xxyyyzz[k] * ab_z + g_0_z_zzzz_xxyyyzzz[k];

                g_0_z_zzzzz_xxyyzzz[k] = g_zzzz_xxyyzzz[k] - g_0_z_zzzz_xxyyzzz[k] * ab_z + g_0_z_zzzz_xxyyzzzz[k];

                g_0_z_zzzzz_xxyzzzz[k] = g_zzzz_xxyzzzz[k] - g_0_z_zzzz_xxyzzzz[k] * ab_z + g_0_z_zzzz_xxyzzzzz[k];

                g_0_z_zzzzz_xxzzzzz[k] = g_zzzz_xxzzzzz[k] - g_0_z_zzzz_xxzzzzz[k] * ab_z + g_0_z_zzzz_xxzzzzzz[k];

                g_0_z_zzzzz_xyyyyyy[k] = g_zzzz_xyyyyyy[k] - g_0_z_zzzz_xyyyyyy[k] * ab_z + g_0_z_zzzz_xyyyyyyz[k];

                g_0_z_zzzzz_xyyyyyz[k] = g_zzzz_xyyyyyz[k] - g_0_z_zzzz_xyyyyyz[k] * ab_z + g_0_z_zzzz_xyyyyyzz[k];

                g_0_z_zzzzz_xyyyyzz[k] = g_zzzz_xyyyyzz[k] - g_0_z_zzzz_xyyyyzz[k] * ab_z + g_0_z_zzzz_xyyyyzzz[k];

                g_0_z_zzzzz_xyyyzzz[k] = g_zzzz_xyyyzzz[k] - g_0_z_zzzz_xyyyzzz[k] * ab_z + g_0_z_zzzz_xyyyzzzz[k];

                g_0_z_zzzzz_xyyzzzz[k] = g_zzzz_xyyzzzz[k] - g_0_z_zzzz_xyyzzzz[k] * ab_z + g_0_z_zzzz_xyyzzzzz[k];

                g_0_z_zzzzz_xyzzzzz[k] = g_zzzz_xyzzzzz[k] - g_0_z_zzzz_xyzzzzz[k] * ab_z + g_0_z_zzzz_xyzzzzzz[k];

                g_0_z_zzzzz_xzzzzzz[k] = g_zzzz_xzzzzzz[k] - g_0_z_zzzz_xzzzzzz[k] * ab_z + g_0_z_zzzz_xzzzzzzz[k];

                g_0_z_zzzzz_yyyyyyy[k] = g_zzzz_yyyyyyy[k] - g_0_z_zzzz_yyyyyyy[k] * ab_z + g_0_z_zzzz_yyyyyyyz[k];

                g_0_z_zzzzz_yyyyyyz[k] = g_zzzz_yyyyyyz[k] - g_0_z_zzzz_yyyyyyz[k] * ab_z + g_0_z_zzzz_yyyyyyzz[k];

                g_0_z_zzzzz_yyyyyzz[k] = g_zzzz_yyyyyzz[k] - g_0_z_zzzz_yyyyyzz[k] * ab_z + g_0_z_zzzz_yyyyyzzz[k];

                g_0_z_zzzzz_yyyyzzz[k] = g_zzzz_yyyyzzz[k] - g_0_z_zzzz_yyyyzzz[k] * ab_z + g_0_z_zzzz_yyyyzzzz[k];

                g_0_z_zzzzz_yyyzzzz[k] = g_zzzz_yyyzzzz[k] - g_0_z_zzzz_yyyzzzz[k] * ab_z + g_0_z_zzzz_yyyzzzzz[k];

                g_0_z_zzzzz_yyzzzzz[k] = g_zzzz_yyzzzzz[k] - g_0_z_zzzz_yyzzzzz[k] * ab_z + g_0_z_zzzz_yyzzzzzz[k];

                g_0_z_zzzzz_yzzzzzz[k] = g_zzzz_yzzzzzz[k] - g_0_z_zzzz_yzzzzzz[k] * ab_z + g_0_z_zzzz_yzzzzzzz[k];

                g_0_z_zzzzz_zzzzzzz[k] = g_zzzz_zzzzzzz[k] - g_0_z_zzzz_zzzzzzz[k] * ab_z + g_0_z_zzzz_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

