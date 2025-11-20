#include "ElectronRepulsionGeom0100ContrRecGKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_gkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gkxx,
                                            const size_t idx_fkxx,
                                            const size_t idx_geom_01_fkxx,
                                            const size_t idx_geom_01_flxx,
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
            /// Set up components of auxilary buffer : FKSS

            const auto fk_off = idx_fkxx + i * dcomps + j;

            auto g_xxx_xxxxxxx = cbuffer.data(fk_off + 0 * ccomps * dcomps);

            auto g_xxx_xxxxxxy = cbuffer.data(fk_off + 1 * ccomps * dcomps);

            auto g_xxx_xxxxxxz = cbuffer.data(fk_off + 2 * ccomps * dcomps);

            auto g_xxx_xxxxxyy = cbuffer.data(fk_off + 3 * ccomps * dcomps);

            auto g_xxx_xxxxxyz = cbuffer.data(fk_off + 4 * ccomps * dcomps);

            auto g_xxx_xxxxxzz = cbuffer.data(fk_off + 5 * ccomps * dcomps);

            auto g_xxx_xxxxyyy = cbuffer.data(fk_off + 6 * ccomps * dcomps);

            auto g_xxx_xxxxyyz = cbuffer.data(fk_off + 7 * ccomps * dcomps);

            auto g_xxx_xxxxyzz = cbuffer.data(fk_off + 8 * ccomps * dcomps);

            auto g_xxx_xxxxzzz = cbuffer.data(fk_off + 9 * ccomps * dcomps);

            auto g_xxx_xxxyyyy = cbuffer.data(fk_off + 10 * ccomps * dcomps);

            auto g_xxx_xxxyyyz = cbuffer.data(fk_off + 11 * ccomps * dcomps);

            auto g_xxx_xxxyyzz = cbuffer.data(fk_off + 12 * ccomps * dcomps);

            auto g_xxx_xxxyzzz = cbuffer.data(fk_off + 13 * ccomps * dcomps);

            auto g_xxx_xxxzzzz = cbuffer.data(fk_off + 14 * ccomps * dcomps);

            auto g_xxx_xxyyyyy = cbuffer.data(fk_off + 15 * ccomps * dcomps);

            auto g_xxx_xxyyyyz = cbuffer.data(fk_off + 16 * ccomps * dcomps);

            auto g_xxx_xxyyyzz = cbuffer.data(fk_off + 17 * ccomps * dcomps);

            auto g_xxx_xxyyzzz = cbuffer.data(fk_off + 18 * ccomps * dcomps);

            auto g_xxx_xxyzzzz = cbuffer.data(fk_off + 19 * ccomps * dcomps);

            auto g_xxx_xxzzzzz = cbuffer.data(fk_off + 20 * ccomps * dcomps);

            auto g_xxx_xyyyyyy = cbuffer.data(fk_off + 21 * ccomps * dcomps);

            auto g_xxx_xyyyyyz = cbuffer.data(fk_off + 22 * ccomps * dcomps);

            auto g_xxx_xyyyyzz = cbuffer.data(fk_off + 23 * ccomps * dcomps);

            auto g_xxx_xyyyzzz = cbuffer.data(fk_off + 24 * ccomps * dcomps);

            auto g_xxx_xyyzzzz = cbuffer.data(fk_off + 25 * ccomps * dcomps);

            auto g_xxx_xyzzzzz = cbuffer.data(fk_off + 26 * ccomps * dcomps);

            auto g_xxx_xzzzzzz = cbuffer.data(fk_off + 27 * ccomps * dcomps);

            auto g_xxx_yyyyyyy = cbuffer.data(fk_off + 28 * ccomps * dcomps);

            auto g_xxx_yyyyyyz = cbuffer.data(fk_off + 29 * ccomps * dcomps);

            auto g_xxx_yyyyyzz = cbuffer.data(fk_off + 30 * ccomps * dcomps);

            auto g_xxx_yyyyzzz = cbuffer.data(fk_off + 31 * ccomps * dcomps);

            auto g_xxx_yyyzzzz = cbuffer.data(fk_off + 32 * ccomps * dcomps);

            auto g_xxx_yyzzzzz = cbuffer.data(fk_off + 33 * ccomps * dcomps);

            auto g_xxx_yzzzzzz = cbuffer.data(fk_off + 34 * ccomps * dcomps);

            auto g_xxx_zzzzzzz = cbuffer.data(fk_off + 35 * ccomps * dcomps);

            auto g_xxy_xxxxxxx = cbuffer.data(fk_off + 36 * ccomps * dcomps);

            auto g_xxy_xxxxxxy = cbuffer.data(fk_off + 37 * ccomps * dcomps);

            auto g_xxy_xxxxxxz = cbuffer.data(fk_off + 38 * ccomps * dcomps);

            auto g_xxy_xxxxxyy = cbuffer.data(fk_off + 39 * ccomps * dcomps);

            auto g_xxy_xxxxxyz = cbuffer.data(fk_off + 40 * ccomps * dcomps);

            auto g_xxy_xxxxxzz = cbuffer.data(fk_off + 41 * ccomps * dcomps);

            auto g_xxy_xxxxyyy = cbuffer.data(fk_off + 42 * ccomps * dcomps);

            auto g_xxy_xxxxyyz = cbuffer.data(fk_off + 43 * ccomps * dcomps);

            auto g_xxy_xxxxyzz = cbuffer.data(fk_off + 44 * ccomps * dcomps);

            auto g_xxy_xxxxzzz = cbuffer.data(fk_off + 45 * ccomps * dcomps);

            auto g_xxy_xxxyyyy = cbuffer.data(fk_off + 46 * ccomps * dcomps);

            auto g_xxy_xxxyyyz = cbuffer.data(fk_off + 47 * ccomps * dcomps);

            auto g_xxy_xxxyyzz = cbuffer.data(fk_off + 48 * ccomps * dcomps);

            auto g_xxy_xxxyzzz = cbuffer.data(fk_off + 49 * ccomps * dcomps);

            auto g_xxy_xxxzzzz = cbuffer.data(fk_off + 50 * ccomps * dcomps);

            auto g_xxy_xxyyyyy = cbuffer.data(fk_off + 51 * ccomps * dcomps);

            auto g_xxy_xxyyyyz = cbuffer.data(fk_off + 52 * ccomps * dcomps);

            auto g_xxy_xxyyyzz = cbuffer.data(fk_off + 53 * ccomps * dcomps);

            auto g_xxy_xxyyzzz = cbuffer.data(fk_off + 54 * ccomps * dcomps);

            auto g_xxy_xxyzzzz = cbuffer.data(fk_off + 55 * ccomps * dcomps);

            auto g_xxy_xxzzzzz = cbuffer.data(fk_off + 56 * ccomps * dcomps);

            auto g_xxy_xyyyyyy = cbuffer.data(fk_off + 57 * ccomps * dcomps);

            auto g_xxy_xyyyyyz = cbuffer.data(fk_off + 58 * ccomps * dcomps);

            auto g_xxy_xyyyyzz = cbuffer.data(fk_off + 59 * ccomps * dcomps);

            auto g_xxy_xyyyzzz = cbuffer.data(fk_off + 60 * ccomps * dcomps);

            auto g_xxy_xyyzzzz = cbuffer.data(fk_off + 61 * ccomps * dcomps);

            auto g_xxy_xyzzzzz = cbuffer.data(fk_off + 62 * ccomps * dcomps);

            auto g_xxy_xzzzzzz = cbuffer.data(fk_off + 63 * ccomps * dcomps);

            auto g_xxy_yyyyyyy = cbuffer.data(fk_off + 64 * ccomps * dcomps);

            auto g_xxy_yyyyyyz = cbuffer.data(fk_off + 65 * ccomps * dcomps);

            auto g_xxy_yyyyyzz = cbuffer.data(fk_off + 66 * ccomps * dcomps);

            auto g_xxy_yyyyzzz = cbuffer.data(fk_off + 67 * ccomps * dcomps);

            auto g_xxy_yyyzzzz = cbuffer.data(fk_off + 68 * ccomps * dcomps);

            auto g_xxy_yyzzzzz = cbuffer.data(fk_off + 69 * ccomps * dcomps);

            auto g_xxy_yzzzzzz = cbuffer.data(fk_off + 70 * ccomps * dcomps);

            auto g_xxy_zzzzzzz = cbuffer.data(fk_off + 71 * ccomps * dcomps);

            auto g_xxz_xxxxxxx = cbuffer.data(fk_off + 72 * ccomps * dcomps);

            auto g_xxz_xxxxxxy = cbuffer.data(fk_off + 73 * ccomps * dcomps);

            auto g_xxz_xxxxxxz = cbuffer.data(fk_off + 74 * ccomps * dcomps);

            auto g_xxz_xxxxxyy = cbuffer.data(fk_off + 75 * ccomps * dcomps);

            auto g_xxz_xxxxxyz = cbuffer.data(fk_off + 76 * ccomps * dcomps);

            auto g_xxz_xxxxxzz = cbuffer.data(fk_off + 77 * ccomps * dcomps);

            auto g_xxz_xxxxyyy = cbuffer.data(fk_off + 78 * ccomps * dcomps);

            auto g_xxz_xxxxyyz = cbuffer.data(fk_off + 79 * ccomps * dcomps);

            auto g_xxz_xxxxyzz = cbuffer.data(fk_off + 80 * ccomps * dcomps);

            auto g_xxz_xxxxzzz = cbuffer.data(fk_off + 81 * ccomps * dcomps);

            auto g_xxz_xxxyyyy = cbuffer.data(fk_off + 82 * ccomps * dcomps);

            auto g_xxz_xxxyyyz = cbuffer.data(fk_off + 83 * ccomps * dcomps);

            auto g_xxz_xxxyyzz = cbuffer.data(fk_off + 84 * ccomps * dcomps);

            auto g_xxz_xxxyzzz = cbuffer.data(fk_off + 85 * ccomps * dcomps);

            auto g_xxz_xxxzzzz = cbuffer.data(fk_off + 86 * ccomps * dcomps);

            auto g_xxz_xxyyyyy = cbuffer.data(fk_off + 87 * ccomps * dcomps);

            auto g_xxz_xxyyyyz = cbuffer.data(fk_off + 88 * ccomps * dcomps);

            auto g_xxz_xxyyyzz = cbuffer.data(fk_off + 89 * ccomps * dcomps);

            auto g_xxz_xxyyzzz = cbuffer.data(fk_off + 90 * ccomps * dcomps);

            auto g_xxz_xxyzzzz = cbuffer.data(fk_off + 91 * ccomps * dcomps);

            auto g_xxz_xxzzzzz = cbuffer.data(fk_off + 92 * ccomps * dcomps);

            auto g_xxz_xyyyyyy = cbuffer.data(fk_off + 93 * ccomps * dcomps);

            auto g_xxz_xyyyyyz = cbuffer.data(fk_off + 94 * ccomps * dcomps);

            auto g_xxz_xyyyyzz = cbuffer.data(fk_off + 95 * ccomps * dcomps);

            auto g_xxz_xyyyzzz = cbuffer.data(fk_off + 96 * ccomps * dcomps);

            auto g_xxz_xyyzzzz = cbuffer.data(fk_off + 97 * ccomps * dcomps);

            auto g_xxz_xyzzzzz = cbuffer.data(fk_off + 98 * ccomps * dcomps);

            auto g_xxz_xzzzzzz = cbuffer.data(fk_off + 99 * ccomps * dcomps);

            auto g_xxz_yyyyyyy = cbuffer.data(fk_off + 100 * ccomps * dcomps);

            auto g_xxz_yyyyyyz = cbuffer.data(fk_off + 101 * ccomps * dcomps);

            auto g_xxz_yyyyyzz = cbuffer.data(fk_off + 102 * ccomps * dcomps);

            auto g_xxz_yyyyzzz = cbuffer.data(fk_off + 103 * ccomps * dcomps);

            auto g_xxz_yyyzzzz = cbuffer.data(fk_off + 104 * ccomps * dcomps);

            auto g_xxz_yyzzzzz = cbuffer.data(fk_off + 105 * ccomps * dcomps);

            auto g_xxz_yzzzzzz = cbuffer.data(fk_off + 106 * ccomps * dcomps);

            auto g_xxz_zzzzzzz = cbuffer.data(fk_off + 107 * ccomps * dcomps);

            auto g_xyy_xxxxxxx = cbuffer.data(fk_off + 108 * ccomps * dcomps);

            auto g_xyy_xxxxxxy = cbuffer.data(fk_off + 109 * ccomps * dcomps);

            auto g_xyy_xxxxxxz = cbuffer.data(fk_off + 110 * ccomps * dcomps);

            auto g_xyy_xxxxxyy = cbuffer.data(fk_off + 111 * ccomps * dcomps);

            auto g_xyy_xxxxxyz = cbuffer.data(fk_off + 112 * ccomps * dcomps);

            auto g_xyy_xxxxxzz = cbuffer.data(fk_off + 113 * ccomps * dcomps);

            auto g_xyy_xxxxyyy = cbuffer.data(fk_off + 114 * ccomps * dcomps);

            auto g_xyy_xxxxyyz = cbuffer.data(fk_off + 115 * ccomps * dcomps);

            auto g_xyy_xxxxyzz = cbuffer.data(fk_off + 116 * ccomps * dcomps);

            auto g_xyy_xxxxzzz = cbuffer.data(fk_off + 117 * ccomps * dcomps);

            auto g_xyy_xxxyyyy = cbuffer.data(fk_off + 118 * ccomps * dcomps);

            auto g_xyy_xxxyyyz = cbuffer.data(fk_off + 119 * ccomps * dcomps);

            auto g_xyy_xxxyyzz = cbuffer.data(fk_off + 120 * ccomps * dcomps);

            auto g_xyy_xxxyzzz = cbuffer.data(fk_off + 121 * ccomps * dcomps);

            auto g_xyy_xxxzzzz = cbuffer.data(fk_off + 122 * ccomps * dcomps);

            auto g_xyy_xxyyyyy = cbuffer.data(fk_off + 123 * ccomps * dcomps);

            auto g_xyy_xxyyyyz = cbuffer.data(fk_off + 124 * ccomps * dcomps);

            auto g_xyy_xxyyyzz = cbuffer.data(fk_off + 125 * ccomps * dcomps);

            auto g_xyy_xxyyzzz = cbuffer.data(fk_off + 126 * ccomps * dcomps);

            auto g_xyy_xxyzzzz = cbuffer.data(fk_off + 127 * ccomps * dcomps);

            auto g_xyy_xxzzzzz = cbuffer.data(fk_off + 128 * ccomps * dcomps);

            auto g_xyy_xyyyyyy = cbuffer.data(fk_off + 129 * ccomps * dcomps);

            auto g_xyy_xyyyyyz = cbuffer.data(fk_off + 130 * ccomps * dcomps);

            auto g_xyy_xyyyyzz = cbuffer.data(fk_off + 131 * ccomps * dcomps);

            auto g_xyy_xyyyzzz = cbuffer.data(fk_off + 132 * ccomps * dcomps);

            auto g_xyy_xyyzzzz = cbuffer.data(fk_off + 133 * ccomps * dcomps);

            auto g_xyy_xyzzzzz = cbuffer.data(fk_off + 134 * ccomps * dcomps);

            auto g_xyy_xzzzzzz = cbuffer.data(fk_off + 135 * ccomps * dcomps);

            auto g_xyy_yyyyyyy = cbuffer.data(fk_off + 136 * ccomps * dcomps);

            auto g_xyy_yyyyyyz = cbuffer.data(fk_off + 137 * ccomps * dcomps);

            auto g_xyy_yyyyyzz = cbuffer.data(fk_off + 138 * ccomps * dcomps);

            auto g_xyy_yyyyzzz = cbuffer.data(fk_off + 139 * ccomps * dcomps);

            auto g_xyy_yyyzzzz = cbuffer.data(fk_off + 140 * ccomps * dcomps);

            auto g_xyy_yyzzzzz = cbuffer.data(fk_off + 141 * ccomps * dcomps);

            auto g_xyy_yzzzzzz = cbuffer.data(fk_off + 142 * ccomps * dcomps);

            auto g_xyy_zzzzzzz = cbuffer.data(fk_off + 143 * ccomps * dcomps);

            auto g_xyz_xxxxxxx = cbuffer.data(fk_off + 144 * ccomps * dcomps);

            auto g_xyz_xxxxxxy = cbuffer.data(fk_off + 145 * ccomps * dcomps);

            auto g_xyz_xxxxxxz = cbuffer.data(fk_off + 146 * ccomps * dcomps);

            auto g_xyz_xxxxxyy = cbuffer.data(fk_off + 147 * ccomps * dcomps);

            auto g_xyz_xxxxxyz = cbuffer.data(fk_off + 148 * ccomps * dcomps);

            auto g_xyz_xxxxxzz = cbuffer.data(fk_off + 149 * ccomps * dcomps);

            auto g_xyz_xxxxyyy = cbuffer.data(fk_off + 150 * ccomps * dcomps);

            auto g_xyz_xxxxyyz = cbuffer.data(fk_off + 151 * ccomps * dcomps);

            auto g_xyz_xxxxyzz = cbuffer.data(fk_off + 152 * ccomps * dcomps);

            auto g_xyz_xxxxzzz = cbuffer.data(fk_off + 153 * ccomps * dcomps);

            auto g_xyz_xxxyyyy = cbuffer.data(fk_off + 154 * ccomps * dcomps);

            auto g_xyz_xxxyyyz = cbuffer.data(fk_off + 155 * ccomps * dcomps);

            auto g_xyz_xxxyyzz = cbuffer.data(fk_off + 156 * ccomps * dcomps);

            auto g_xyz_xxxyzzz = cbuffer.data(fk_off + 157 * ccomps * dcomps);

            auto g_xyz_xxxzzzz = cbuffer.data(fk_off + 158 * ccomps * dcomps);

            auto g_xyz_xxyyyyy = cbuffer.data(fk_off + 159 * ccomps * dcomps);

            auto g_xyz_xxyyyyz = cbuffer.data(fk_off + 160 * ccomps * dcomps);

            auto g_xyz_xxyyyzz = cbuffer.data(fk_off + 161 * ccomps * dcomps);

            auto g_xyz_xxyyzzz = cbuffer.data(fk_off + 162 * ccomps * dcomps);

            auto g_xyz_xxyzzzz = cbuffer.data(fk_off + 163 * ccomps * dcomps);

            auto g_xyz_xxzzzzz = cbuffer.data(fk_off + 164 * ccomps * dcomps);

            auto g_xyz_xyyyyyy = cbuffer.data(fk_off + 165 * ccomps * dcomps);

            auto g_xyz_xyyyyyz = cbuffer.data(fk_off + 166 * ccomps * dcomps);

            auto g_xyz_xyyyyzz = cbuffer.data(fk_off + 167 * ccomps * dcomps);

            auto g_xyz_xyyyzzz = cbuffer.data(fk_off + 168 * ccomps * dcomps);

            auto g_xyz_xyyzzzz = cbuffer.data(fk_off + 169 * ccomps * dcomps);

            auto g_xyz_xyzzzzz = cbuffer.data(fk_off + 170 * ccomps * dcomps);

            auto g_xyz_xzzzzzz = cbuffer.data(fk_off + 171 * ccomps * dcomps);

            auto g_xyz_yyyyyyy = cbuffer.data(fk_off + 172 * ccomps * dcomps);

            auto g_xyz_yyyyyyz = cbuffer.data(fk_off + 173 * ccomps * dcomps);

            auto g_xyz_yyyyyzz = cbuffer.data(fk_off + 174 * ccomps * dcomps);

            auto g_xyz_yyyyzzz = cbuffer.data(fk_off + 175 * ccomps * dcomps);

            auto g_xyz_yyyzzzz = cbuffer.data(fk_off + 176 * ccomps * dcomps);

            auto g_xyz_yyzzzzz = cbuffer.data(fk_off + 177 * ccomps * dcomps);

            auto g_xyz_yzzzzzz = cbuffer.data(fk_off + 178 * ccomps * dcomps);

            auto g_xyz_zzzzzzz = cbuffer.data(fk_off + 179 * ccomps * dcomps);

            auto g_xzz_xxxxxxx = cbuffer.data(fk_off + 180 * ccomps * dcomps);

            auto g_xzz_xxxxxxy = cbuffer.data(fk_off + 181 * ccomps * dcomps);

            auto g_xzz_xxxxxxz = cbuffer.data(fk_off + 182 * ccomps * dcomps);

            auto g_xzz_xxxxxyy = cbuffer.data(fk_off + 183 * ccomps * dcomps);

            auto g_xzz_xxxxxyz = cbuffer.data(fk_off + 184 * ccomps * dcomps);

            auto g_xzz_xxxxxzz = cbuffer.data(fk_off + 185 * ccomps * dcomps);

            auto g_xzz_xxxxyyy = cbuffer.data(fk_off + 186 * ccomps * dcomps);

            auto g_xzz_xxxxyyz = cbuffer.data(fk_off + 187 * ccomps * dcomps);

            auto g_xzz_xxxxyzz = cbuffer.data(fk_off + 188 * ccomps * dcomps);

            auto g_xzz_xxxxzzz = cbuffer.data(fk_off + 189 * ccomps * dcomps);

            auto g_xzz_xxxyyyy = cbuffer.data(fk_off + 190 * ccomps * dcomps);

            auto g_xzz_xxxyyyz = cbuffer.data(fk_off + 191 * ccomps * dcomps);

            auto g_xzz_xxxyyzz = cbuffer.data(fk_off + 192 * ccomps * dcomps);

            auto g_xzz_xxxyzzz = cbuffer.data(fk_off + 193 * ccomps * dcomps);

            auto g_xzz_xxxzzzz = cbuffer.data(fk_off + 194 * ccomps * dcomps);

            auto g_xzz_xxyyyyy = cbuffer.data(fk_off + 195 * ccomps * dcomps);

            auto g_xzz_xxyyyyz = cbuffer.data(fk_off + 196 * ccomps * dcomps);

            auto g_xzz_xxyyyzz = cbuffer.data(fk_off + 197 * ccomps * dcomps);

            auto g_xzz_xxyyzzz = cbuffer.data(fk_off + 198 * ccomps * dcomps);

            auto g_xzz_xxyzzzz = cbuffer.data(fk_off + 199 * ccomps * dcomps);

            auto g_xzz_xxzzzzz = cbuffer.data(fk_off + 200 * ccomps * dcomps);

            auto g_xzz_xyyyyyy = cbuffer.data(fk_off + 201 * ccomps * dcomps);

            auto g_xzz_xyyyyyz = cbuffer.data(fk_off + 202 * ccomps * dcomps);

            auto g_xzz_xyyyyzz = cbuffer.data(fk_off + 203 * ccomps * dcomps);

            auto g_xzz_xyyyzzz = cbuffer.data(fk_off + 204 * ccomps * dcomps);

            auto g_xzz_xyyzzzz = cbuffer.data(fk_off + 205 * ccomps * dcomps);

            auto g_xzz_xyzzzzz = cbuffer.data(fk_off + 206 * ccomps * dcomps);

            auto g_xzz_xzzzzzz = cbuffer.data(fk_off + 207 * ccomps * dcomps);

            auto g_xzz_yyyyyyy = cbuffer.data(fk_off + 208 * ccomps * dcomps);

            auto g_xzz_yyyyyyz = cbuffer.data(fk_off + 209 * ccomps * dcomps);

            auto g_xzz_yyyyyzz = cbuffer.data(fk_off + 210 * ccomps * dcomps);

            auto g_xzz_yyyyzzz = cbuffer.data(fk_off + 211 * ccomps * dcomps);

            auto g_xzz_yyyzzzz = cbuffer.data(fk_off + 212 * ccomps * dcomps);

            auto g_xzz_yyzzzzz = cbuffer.data(fk_off + 213 * ccomps * dcomps);

            auto g_xzz_yzzzzzz = cbuffer.data(fk_off + 214 * ccomps * dcomps);

            auto g_xzz_zzzzzzz = cbuffer.data(fk_off + 215 * ccomps * dcomps);

            auto g_yyy_xxxxxxx = cbuffer.data(fk_off + 216 * ccomps * dcomps);

            auto g_yyy_xxxxxxy = cbuffer.data(fk_off + 217 * ccomps * dcomps);

            auto g_yyy_xxxxxxz = cbuffer.data(fk_off + 218 * ccomps * dcomps);

            auto g_yyy_xxxxxyy = cbuffer.data(fk_off + 219 * ccomps * dcomps);

            auto g_yyy_xxxxxyz = cbuffer.data(fk_off + 220 * ccomps * dcomps);

            auto g_yyy_xxxxxzz = cbuffer.data(fk_off + 221 * ccomps * dcomps);

            auto g_yyy_xxxxyyy = cbuffer.data(fk_off + 222 * ccomps * dcomps);

            auto g_yyy_xxxxyyz = cbuffer.data(fk_off + 223 * ccomps * dcomps);

            auto g_yyy_xxxxyzz = cbuffer.data(fk_off + 224 * ccomps * dcomps);

            auto g_yyy_xxxxzzz = cbuffer.data(fk_off + 225 * ccomps * dcomps);

            auto g_yyy_xxxyyyy = cbuffer.data(fk_off + 226 * ccomps * dcomps);

            auto g_yyy_xxxyyyz = cbuffer.data(fk_off + 227 * ccomps * dcomps);

            auto g_yyy_xxxyyzz = cbuffer.data(fk_off + 228 * ccomps * dcomps);

            auto g_yyy_xxxyzzz = cbuffer.data(fk_off + 229 * ccomps * dcomps);

            auto g_yyy_xxxzzzz = cbuffer.data(fk_off + 230 * ccomps * dcomps);

            auto g_yyy_xxyyyyy = cbuffer.data(fk_off + 231 * ccomps * dcomps);

            auto g_yyy_xxyyyyz = cbuffer.data(fk_off + 232 * ccomps * dcomps);

            auto g_yyy_xxyyyzz = cbuffer.data(fk_off + 233 * ccomps * dcomps);

            auto g_yyy_xxyyzzz = cbuffer.data(fk_off + 234 * ccomps * dcomps);

            auto g_yyy_xxyzzzz = cbuffer.data(fk_off + 235 * ccomps * dcomps);

            auto g_yyy_xxzzzzz = cbuffer.data(fk_off + 236 * ccomps * dcomps);

            auto g_yyy_xyyyyyy = cbuffer.data(fk_off + 237 * ccomps * dcomps);

            auto g_yyy_xyyyyyz = cbuffer.data(fk_off + 238 * ccomps * dcomps);

            auto g_yyy_xyyyyzz = cbuffer.data(fk_off + 239 * ccomps * dcomps);

            auto g_yyy_xyyyzzz = cbuffer.data(fk_off + 240 * ccomps * dcomps);

            auto g_yyy_xyyzzzz = cbuffer.data(fk_off + 241 * ccomps * dcomps);

            auto g_yyy_xyzzzzz = cbuffer.data(fk_off + 242 * ccomps * dcomps);

            auto g_yyy_xzzzzzz = cbuffer.data(fk_off + 243 * ccomps * dcomps);

            auto g_yyy_yyyyyyy = cbuffer.data(fk_off + 244 * ccomps * dcomps);

            auto g_yyy_yyyyyyz = cbuffer.data(fk_off + 245 * ccomps * dcomps);

            auto g_yyy_yyyyyzz = cbuffer.data(fk_off + 246 * ccomps * dcomps);

            auto g_yyy_yyyyzzz = cbuffer.data(fk_off + 247 * ccomps * dcomps);

            auto g_yyy_yyyzzzz = cbuffer.data(fk_off + 248 * ccomps * dcomps);

            auto g_yyy_yyzzzzz = cbuffer.data(fk_off + 249 * ccomps * dcomps);

            auto g_yyy_yzzzzzz = cbuffer.data(fk_off + 250 * ccomps * dcomps);

            auto g_yyy_zzzzzzz = cbuffer.data(fk_off + 251 * ccomps * dcomps);

            auto g_yyz_xxxxxxx = cbuffer.data(fk_off + 252 * ccomps * dcomps);

            auto g_yyz_xxxxxxy = cbuffer.data(fk_off + 253 * ccomps * dcomps);

            auto g_yyz_xxxxxxz = cbuffer.data(fk_off + 254 * ccomps * dcomps);

            auto g_yyz_xxxxxyy = cbuffer.data(fk_off + 255 * ccomps * dcomps);

            auto g_yyz_xxxxxyz = cbuffer.data(fk_off + 256 * ccomps * dcomps);

            auto g_yyz_xxxxxzz = cbuffer.data(fk_off + 257 * ccomps * dcomps);

            auto g_yyz_xxxxyyy = cbuffer.data(fk_off + 258 * ccomps * dcomps);

            auto g_yyz_xxxxyyz = cbuffer.data(fk_off + 259 * ccomps * dcomps);

            auto g_yyz_xxxxyzz = cbuffer.data(fk_off + 260 * ccomps * dcomps);

            auto g_yyz_xxxxzzz = cbuffer.data(fk_off + 261 * ccomps * dcomps);

            auto g_yyz_xxxyyyy = cbuffer.data(fk_off + 262 * ccomps * dcomps);

            auto g_yyz_xxxyyyz = cbuffer.data(fk_off + 263 * ccomps * dcomps);

            auto g_yyz_xxxyyzz = cbuffer.data(fk_off + 264 * ccomps * dcomps);

            auto g_yyz_xxxyzzz = cbuffer.data(fk_off + 265 * ccomps * dcomps);

            auto g_yyz_xxxzzzz = cbuffer.data(fk_off + 266 * ccomps * dcomps);

            auto g_yyz_xxyyyyy = cbuffer.data(fk_off + 267 * ccomps * dcomps);

            auto g_yyz_xxyyyyz = cbuffer.data(fk_off + 268 * ccomps * dcomps);

            auto g_yyz_xxyyyzz = cbuffer.data(fk_off + 269 * ccomps * dcomps);

            auto g_yyz_xxyyzzz = cbuffer.data(fk_off + 270 * ccomps * dcomps);

            auto g_yyz_xxyzzzz = cbuffer.data(fk_off + 271 * ccomps * dcomps);

            auto g_yyz_xxzzzzz = cbuffer.data(fk_off + 272 * ccomps * dcomps);

            auto g_yyz_xyyyyyy = cbuffer.data(fk_off + 273 * ccomps * dcomps);

            auto g_yyz_xyyyyyz = cbuffer.data(fk_off + 274 * ccomps * dcomps);

            auto g_yyz_xyyyyzz = cbuffer.data(fk_off + 275 * ccomps * dcomps);

            auto g_yyz_xyyyzzz = cbuffer.data(fk_off + 276 * ccomps * dcomps);

            auto g_yyz_xyyzzzz = cbuffer.data(fk_off + 277 * ccomps * dcomps);

            auto g_yyz_xyzzzzz = cbuffer.data(fk_off + 278 * ccomps * dcomps);

            auto g_yyz_xzzzzzz = cbuffer.data(fk_off + 279 * ccomps * dcomps);

            auto g_yyz_yyyyyyy = cbuffer.data(fk_off + 280 * ccomps * dcomps);

            auto g_yyz_yyyyyyz = cbuffer.data(fk_off + 281 * ccomps * dcomps);

            auto g_yyz_yyyyyzz = cbuffer.data(fk_off + 282 * ccomps * dcomps);

            auto g_yyz_yyyyzzz = cbuffer.data(fk_off + 283 * ccomps * dcomps);

            auto g_yyz_yyyzzzz = cbuffer.data(fk_off + 284 * ccomps * dcomps);

            auto g_yyz_yyzzzzz = cbuffer.data(fk_off + 285 * ccomps * dcomps);

            auto g_yyz_yzzzzzz = cbuffer.data(fk_off + 286 * ccomps * dcomps);

            auto g_yyz_zzzzzzz = cbuffer.data(fk_off + 287 * ccomps * dcomps);

            auto g_yzz_xxxxxxx = cbuffer.data(fk_off + 288 * ccomps * dcomps);

            auto g_yzz_xxxxxxy = cbuffer.data(fk_off + 289 * ccomps * dcomps);

            auto g_yzz_xxxxxxz = cbuffer.data(fk_off + 290 * ccomps * dcomps);

            auto g_yzz_xxxxxyy = cbuffer.data(fk_off + 291 * ccomps * dcomps);

            auto g_yzz_xxxxxyz = cbuffer.data(fk_off + 292 * ccomps * dcomps);

            auto g_yzz_xxxxxzz = cbuffer.data(fk_off + 293 * ccomps * dcomps);

            auto g_yzz_xxxxyyy = cbuffer.data(fk_off + 294 * ccomps * dcomps);

            auto g_yzz_xxxxyyz = cbuffer.data(fk_off + 295 * ccomps * dcomps);

            auto g_yzz_xxxxyzz = cbuffer.data(fk_off + 296 * ccomps * dcomps);

            auto g_yzz_xxxxzzz = cbuffer.data(fk_off + 297 * ccomps * dcomps);

            auto g_yzz_xxxyyyy = cbuffer.data(fk_off + 298 * ccomps * dcomps);

            auto g_yzz_xxxyyyz = cbuffer.data(fk_off + 299 * ccomps * dcomps);

            auto g_yzz_xxxyyzz = cbuffer.data(fk_off + 300 * ccomps * dcomps);

            auto g_yzz_xxxyzzz = cbuffer.data(fk_off + 301 * ccomps * dcomps);

            auto g_yzz_xxxzzzz = cbuffer.data(fk_off + 302 * ccomps * dcomps);

            auto g_yzz_xxyyyyy = cbuffer.data(fk_off + 303 * ccomps * dcomps);

            auto g_yzz_xxyyyyz = cbuffer.data(fk_off + 304 * ccomps * dcomps);

            auto g_yzz_xxyyyzz = cbuffer.data(fk_off + 305 * ccomps * dcomps);

            auto g_yzz_xxyyzzz = cbuffer.data(fk_off + 306 * ccomps * dcomps);

            auto g_yzz_xxyzzzz = cbuffer.data(fk_off + 307 * ccomps * dcomps);

            auto g_yzz_xxzzzzz = cbuffer.data(fk_off + 308 * ccomps * dcomps);

            auto g_yzz_xyyyyyy = cbuffer.data(fk_off + 309 * ccomps * dcomps);

            auto g_yzz_xyyyyyz = cbuffer.data(fk_off + 310 * ccomps * dcomps);

            auto g_yzz_xyyyyzz = cbuffer.data(fk_off + 311 * ccomps * dcomps);

            auto g_yzz_xyyyzzz = cbuffer.data(fk_off + 312 * ccomps * dcomps);

            auto g_yzz_xyyzzzz = cbuffer.data(fk_off + 313 * ccomps * dcomps);

            auto g_yzz_xyzzzzz = cbuffer.data(fk_off + 314 * ccomps * dcomps);

            auto g_yzz_xzzzzzz = cbuffer.data(fk_off + 315 * ccomps * dcomps);

            auto g_yzz_yyyyyyy = cbuffer.data(fk_off + 316 * ccomps * dcomps);

            auto g_yzz_yyyyyyz = cbuffer.data(fk_off + 317 * ccomps * dcomps);

            auto g_yzz_yyyyyzz = cbuffer.data(fk_off + 318 * ccomps * dcomps);

            auto g_yzz_yyyyzzz = cbuffer.data(fk_off + 319 * ccomps * dcomps);

            auto g_yzz_yyyzzzz = cbuffer.data(fk_off + 320 * ccomps * dcomps);

            auto g_yzz_yyzzzzz = cbuffer.data(fk_off + 321 * ccomps * dcomps);

            auto g_yzz_yzzzzzz = cbuffer.data(fk_off + 322 * ccomps * dcomps);

            auto g_yzz_zzzzzzz = cbuffer.data(fk_off + 323 * ccomps * dcomps);

            auto g_zzz_xxxxxxx = cbuffer.data(fk_off + 324 * ccomps * dcomps);

            auto g_zzz_xxxxxxy = cbuffer.data(fk_off + 325 * ccomps * dcomps);

            auto g_zzz_xxxxxxz = cbuffer.data(fk_off + 326 * ccomps * dcomps);

            auto g_zzz_xxxxxyy = cbuffer.data(fk_off + 327 * ccomps * dcomps);

            auto g_zzz_xxxxxyz = cbuffer.data(fk_off + 328 * ccomps * dcomps);

            auto g_zzz_xxxxxzz = cbuffer.data(fk_off + 329 * ccomps * dcomps);

            auto g_zzz_xxxxyyy = cbuffer.data(fk_off + 330 * ccomps * dcomps);

            auto g_zzz_xxxxyyz = cbuffer.data(fk_off + 331 * ccomps * dcomps);

            auto g_zzz_xxxxyzz = cbuffer.data(fk_off + 332 * ccomps * dcomps);

            auto g_zzz_xxxxzzz = cbuffer.data(fk_off + 333 * ccomps * dcomps);

            auto g_zzz_xxxyyyy = cbuffer.data(fk_off + 334 * ccomps * dcomps);

            auto g_zzz_xxxyyyz = cbuffer.data(fk_off + 335 * ccomps * dcomps);

            auto g_zzz_xxxyyzz = cbuffer.data(fk_off + 336 * ccomps * dcomps);

            auto g_zzz_xxxyzzz = cbuffer.data(fk_off + 337 * ccomps * dcomps);

            auto g_zzz_xxxzzzz = cbuffer.data(fk_off + 338 * ccomps * dcomps);

            auto g_zzz_xxyyyyy = cbuffer.data(fk_off + 339 * ccomps * dcomps);

            auto g_zzz_xxyyyyz = cbuffer.data(fk_off + 340 * ccomps * dcomps);

            auto g_zzz_xxyyyzz = cbuffer.data(fk_off + 341 * ccomps * dcomps);

            auto g_zzz_xxyyzzz = cbuffer.data(fk_off + 342 * ccomps * dcomps);

            auto g_zzz_xxyzzzz = cbuffer.data(fk_off + 343 * ccomps * dcomps);

            auto g_zzz_xxzzzzz = cbuffer.data(fk_off + 344 * ccomps * dcomps);

            auto g_zzz_xyyyyyy = cbuffer.data(fk_off + 345 * ccomps * dcomps);

            auto g_zzz_xyyyyyz = cbuffer.data(fk_off + 346 * ccomps * dcomps);

            auto g_zzz_xyyyyzz = cbuffer.data(fk_off + 347 * ccomps * dcomps);

            auto g_zzz_xyyyzzz = cbuffer.data(fk_off + 348 * ccomps * dcomps);

            auto g_zzz_xyyzzzz = cbuffer.data(fk_off + 349 * ccomps * dcomps);

            auto g_zzz_xyzzzzz = cbuffer.data(fk_off + 350 * ccomps * dcomps);

            auto g_zzz_xzzzzzz = cbuffer.data(fk_off + 351 * ccomps * dcomps);

            auto g_zzz_yyyyyyy = cbuffer.data(fk_off + 352 * ccomps * dcomps);

            auto g_zzz_yyyyyyz = cbuffer.data(fk_off + 353 * ccomps * dcomps);

            auto g_zzz_yyyyyzz = cbuffer.data(fk_off + 354 * ccomps * dcomps);

            auto g_zzz_yyyyzzz = cbuffer.data(fk_off + 355 * ccomps * dcomps);

            auto g_zzz_yyyzzzz = cbuffer.data(fk_off + 356 * ccomps * dcomps);

            auto g_zzz_yyzzzzz = cbuffer.data(fk_off + 357 * ccomps * dcomps);

            auto g_zzz_yzzzzzz = cbuffer.data(fk_off + 358 * ccomps * dcomps);

            auto g_zzz_zzzzzzz = cbuffer.data(fk_off + 359 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FKSS

            const auto fk_geom_01_off = idx_geom_01_fkxx + i * dcomps + j;

            auto g_0_x_xxx_xxxxxxx = cbuffer.data(fk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxy = cbuffer.data(fk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxz = cbuffer.data(fk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxyy = cbuffer.data(fk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxyz = cbuffer.data(fk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxzz = cbuffer.data(fk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyyy = cbuffer.data(fk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyyz = cbuffer.data(fk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyzz = cbuffer.data(fk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxzzz = cbuffer.data(fk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyyy = cbuffer.data(fk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyyz = cbuffer.data(fk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyzz = cbuffer.data(fk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyzzz = cbuffer.data(fk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_xxxzzzz = cbuffer.data(fk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyyy = cbuffer.data(fk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyyz = cbuffer.data(fk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyzz = cbuffer.data(fk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyzzz = cbuffer.data(fk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxx_xxyzzzz = cbuffer.data(fk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxx_xxzzzzz = cbuffer.data(fk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyyy = cbuffer.data(fk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyyz = cbuffer.data(fk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyzz = cbuffer.data(fk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyzzz = cbuffer.data(fk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxx_xyyzzzz = cbuffer.data(fk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxx_xyzzzzz = cbuffer.data(fk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxx_xzzzzzz = cbuffer.data(fk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyyy = cbuffer.data(fk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyyz = cbuffer.data(fk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyzz = cbuffer.data(fk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyzzz = cbuffer.data(fk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxx_yyyzzzz = cbuffer.data(fk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxx_yyzzzzz = cbuffer.data(fk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxx_yzzzzzz = cbuffer.data(fk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxx_zzzzzzz = cbuffer.data(fk_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxx = cbuffer.data(fk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxy = cbuffer.data(fk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxz = cbuffer.data(fk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxyy = cbuffer.data(fk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxyz = cbuffer.data(fk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxzz = cbuffer.data(fk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyyy = cbuffer.data(fk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyyz = cbuffer.data(fk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyzz = cbuffer.data(fk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxzzz = cbuffer.data(fk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyyy = cbuffer.data(fk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyyz = cbuffer.data(fk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyzz = cbuffer.data(fk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyzzz = cbuffer.data(fk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxy_xxxzzzz = cbuffer.data(fk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyyy = cbuffer.data(fk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyyz = cbuffer.data(fk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyzz = cbuffer.data(fk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyzzz = cbuffer.data(fk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxy_xxyzzzz = cbuffer.data(fk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxy_xxzzzzz = cbuffer.data(fk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyyy = cbuffer.data(fk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyyz = cbuffer.data(fk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyzz = cbuffer.data(fk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyzzz = cbuffer.data(fk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxy_xyyzzzz = cbuffer.data(fk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxy_xyzzzzz = cbuffer.data(fk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxy_xzzzzzz = cbuffer.data(fk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyyy = cbuffer.data(fk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyyz = cbuffer.data(fk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyzz = cbuffer.data(fk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyzzz = cbuffer.data(fk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxy_yyyzzzz = cbuffer.data(fk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxy_yyzzzzz = cbuffer.data(fk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxy_yzzzzzz = cbuffer.data(fk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxy_zzzzzzz = cbuffer.data(fk_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxx = cbuffer.data(fk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxy = cbuffer.data(fk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxz = cbuffer.data(fk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxyy = cbuffer.data(fk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxyz = cbuffer.data(fk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxzz = cbuffer.data(fk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyyy = cbuffer.data(fk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyyz = cbuffer.data(fk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyzz = cbuffer.data(fk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxzzz = cbuffer.data(fk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyyy = cbuffer.data(fk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyyz = cbuffer.data(fk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyzz = cbuffer.data(fk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyzzz = cbuffer.data(fk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxz_xxxzzzz = cbuffer.data(fk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyyy = cbuffer.data(fk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyyz = cbuffer.data(fk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyzz = cbuffer.data(fk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyzzz = cbuffer.data(fk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxz_xxyzzzz = cbuffer.data(fk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxz_xxzzzzz = cbuffer.data(fk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyyy = cbuffer.data(fk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyyz = cbuffer.data(fk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyzz = cbuffer.data(fk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyzzz = cbuffer.data(fk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxz_xyyzzzz = cbuffer.data(fk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxz_xyzzzzz = cbuffer.data(fk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxz_xzzzzzz = cbuffer.data(fk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyyy = cbuffer.data(fk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyyz = cbuffer.data(fk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyzz = cbuffer.data(fk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyzzz = cbuffer.data(fk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxz_yyyzzzz = cbuffer.data(fk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxz_yyzzzzz = cbuffer.data(fk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxz_yzzzzzz = cbuffer.data(fk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxz_zzzzzzz = cbuffer.data(fk_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_yyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_yyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_yyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_yyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_yyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_yyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_yzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_yzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_yzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_yzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_yzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_yzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_yzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_zzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_zzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_zzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_zzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_zzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_zzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_zzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_zzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_zzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_zzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxx = cbuffer.data(fk_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxy = cbuffer.data(fk_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxz = cbuffer.data(fk_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxyy = cbuffer.data(fk_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxyz = cbuffer.data(fk_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxzz = cbuffer.data(fk_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyyy = cbuffer.data(fk_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyyz = cbuffer.data(fk_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyzz = cbuffer.data(fk_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxzzz = cbuffer.data(fk_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyyy = cbuffer.data(fk_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyyz = cbuffer.data(fk_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyzz = cbuffer.data(fk_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyzzz = cbuffer.data(fk_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxx_xxxzzzz = cbuffer.data(fk_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyyy = cbuffer.data(fk_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyyz = cbuffer.data(fk_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyzz = cbuffer.data(fk_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyzzz = cbuffer.data(fk_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxx_xxyzzzz = cbuffer.data(fk_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xxx_xxzzzzz = cbuffer.data(fk_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyyy = cbuffer.data(fk_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyyz = cbuffer.data(fk_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyzz = cbuffer.data(fk_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyzzz = cbuffer.data(fk_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxx_xyyzzzz = cbuffer.data(fk_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxx_xyzzzzz = cbuffer.data(fk_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxx_xzzzzzz = cbuffer.data(fk_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyyy = cbuffer.data(fk_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyyz = cbuffer.data(fk_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyzz = cbuffer.data(fk_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyzzz = cbuffer.data(fk_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxx_yyyzzzz = cbuffer.data(fk_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxx_yyzzzzz = cbuffer.data(fk_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxx_yzzzzzz = cbuffer.data(fk_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxx_zzzzzzz = cbuffer.data(fk_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxx = cbuffer.data(fk_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxy = cbuffer.data(fk_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxz = cbuffer.data(fk_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxyy = cbuffer.data(fk_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxyz = cbuffer.data(fk_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxzz = cbuffer.data(fk_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyyy = cbuffer.data(fk_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyyz = cbuffer.data(fk_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyzz = cbuffer.data(fk_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxzzz = cbuffer.data(fk_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyyy = cbuffer.data(fk_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyyz = cbuffer.data(fk_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyzz = cbuffer.data(fk_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyzzz = cbuffer.data(fk_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xxy_xxxzzzz = cbuffer.data(fk_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyyy = cbuffer.data(fk_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyyz = cbuffer.data(fk_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyzz = cbuffer.data(fk_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyzzz = cbuffer.data(fk_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxy_xxyzzzz = cbuffer.data(fk_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxy_xxzzzzz = cbuffer.data(fk_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyyy = cbuffer.data(fk_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyyz = cbuffer.data(fk_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyzz = cbuffer.data(fk_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyzzz = cbuffer.data(fk_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxy_xyyzzzz = cbuffer.data(fk_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxy_xyzzzzz = cbuffer.data(fk_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxy_xzzzzzz = cbuffer.data(fk_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyyy = cbuffer.data(fk_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyyz = cbuffer.data(fk_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyzz = cbuffer.data(fk_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyzzz = cbuffer.data(fk_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxy_yyyzzzz = cbuffer.data(fk_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxy_yyzzzzz = cbuffer.data(fk_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xxy_yzzzzzz = cbuffer.data(fk_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxy_zzzzzzz = cbuffer.data(fk_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxx = cbuffer.data(fk_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxy = cbuffer.data(fk_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxz = cbuffer.data(fk_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxyy = cbuffer.data(fk_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxyz = cbuffer.data(fk_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxzz = cbuffer.data(fk_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyyy = cbuffer.data(fk_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyyz = cbuffer.data(fk_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyzz = cbuffer.data(fk_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxzzz = cbuffer.data(fk_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyyy = cbuffer.data(fk_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyyz = cbuffer.data(fk_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyzz = cbuffer.data(fk_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyzzz = cbuffer.data(fk_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxz_xxxzzzz = cbuffer.data(fk_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyyy = cbuffer.data(fk_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyyz = cbuffer.data(fk_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyzz = cbuffer.data(fk_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyzzz = cbuffer.data(fk_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxz_xxyzzzz = cbuffer.data(fk_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxz_xxzzzzz = cbuffer.data(fk_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyyy = cbuffer.data(fk_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyyz = cbuffer.data(fk_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyzz = cbuffer.data(fk_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyzzz = cbuffer.data(fk_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxz_xyyzzzz = cbuffer.data(fk_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxz_xyzzzzz = cbuffer.data(fk_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxz_xzzzzzz = cbuffer.data(fk_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyyy = cbuffer.data(fk_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyyz = cbuffer.data(fk_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyzz = cbuffer.data(fk_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyzzz = cbuffer.data(fk_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxz_yyyzzzz = cbuffer.data(fk_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xxz_yyzzzzz = cbuffer.data(fk_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxz_yzzzzzz = cbuffer.data(fk_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxz_zzzzzzz = cbuffer.data(fk_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_yyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_yyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_yyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_yyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_yyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_yyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_yyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_yyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_yyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_yyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_yyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_yyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_yyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_yyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_yyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_yyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_yyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_yyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_yyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_yyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_yzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_yzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_yzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_yzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_yzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_yzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_yzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_yzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_yzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_yzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_zzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_zzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_zzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_zzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_zzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_zzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_zzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_zzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_zzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_zzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxx = cbuffer.data(fk_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxy = cbuffer.data(fk_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxz = cbuffer.data(fk_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxyy = cbuffer.data(fk_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxyz = cbuffer.data(fk_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxzz = cbuffer.data(fk_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyyy = cbuffer.data(fk_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyyz = cbuffer.data(fk_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyzz = cbuffer.data(fk_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxzzz = cbuffer.data(fk_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyyy = cbuffer.data(fk_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyyz = cbuffer.data(fk_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyzz = cbuffer.data(fk_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyzzz = cbuffer.data(fk_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xxx_xxxzzzz = cbuffer.data(fk_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyyy = cbuffer.data(fk_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyyz = cbuffer.data(fk_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyzz = cbuffer.data(fk_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyzzz = cbuffer.data(fk_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xxx_xxyzzzz = cbuffer.data(fk_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_xxx_xxzzzzz = cbuffer.data(fk_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyyy = cbuffer.data(fk_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyyz = cbuffer.data(fk_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyzz = cbuffer.data(fk_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyzzz = cbuffer.data(fk_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xxx_xyyzzzz = cbuffer.data(fk_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xxx_xyzzzzz = cbuffer.data(fk_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xxx_xzzzzzz = cbuffer.data(fk_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyyy = cbuffer.data(fk_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyyz = cbuffer.data(fk_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyzz = cbuffer.data(fk_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyzzz = cbuffer.data(fk_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xxx_yyyzzzz = cbuffer.data(fk_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xxx_yyzzzzz = cbuffer.data(fk_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xxx_yzzzzzz = cbuffer.data(fk_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xxx_zzzzzzz = cbuffer.data(fk_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxx = cbuffer.data(fk_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxy = cbuffer.data(fk_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxz = cbuffer.data(fk_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxyy = cbuffer.data(fk_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxyz = cbuffer.data(fk_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxzz = cbuffer.data(fk_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyyy = cbuffer.data(fk_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyyz = cbuffer.data(fk_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyzz = cbuffer.data(fk_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxzzz = cbuffer.data(fk_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyyy = cbuffer.data(fk_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyyz = cbuffer.data(fk_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyzz = cbuffer.data(fk_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyzzz = cbuffer.data(fk_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_xxy_xxxzzzz = cbuffer.data(fk_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyyy = cbuffer.data(fk_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyyz = cbuffer.data(fk_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyzz = cbuffer.data(fk_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyzzz = cbuffer.data(fk_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_xxy_xxyzzzz = cbuffer.data(fk_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_xxy_xxzzzzz = cbuffer.data(fk_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyyy = cbuffer.data(fk_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyyz = cbuffer.data(fk_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyzz = cbuffer.data(fk_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyzzz = cbuffer.data(fk_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_xxy_xyyzzzz = cbuffer.data(fk_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_xxy_xyzzzzz = cbuffer.data(fk_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_xxy_xzzzzzz = cbuffer.data(fk_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyyy = cbuffer.data(fk_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyyz = cbuffer.data(fk_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyzz = cbuffer.data(fk_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyzzz = cbuffer.data(fk_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_xxy_yyyzzzz = cbuffer.data(fk_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_xxy_yyzzzzz = cbuffer.data(fk_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_xxy_yzzzzzz = cbuffer.data(fk_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_xxy_zzzzzzz = cbuffer.data(fk_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxx = cbuffer.data(fk_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxy = cbuffer.data(fk_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxz = cbuffer.data(fk_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxyy = cbuffer.data(fk_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxyz = cbuffer.data(fk_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxzz = cbuffer.data(fk_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyyy = cbuffer.data(fk_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyyz = cbuffer.data(fk_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyzz = cbuffer.data(fk_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxzzz = cbuffer.data(fk_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyyy = cbuffer.data(fk_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyyz = cbuffer.data(fk_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyzz = cbuffer.data(fk_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyzzz = cbuffer.data(fk_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_xxz_xxxzzzz = cbuffer.data(fk_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyyy = cbuffer.data(fk_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyyz = cbuffer.data(fk_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyzz = cbuffer.data(fk_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyzzz = cbuffer.data(fk_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_xxz_xxyzzzz = cbuffer.data(fk_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_xxz_xxzzzzz = cbuffer.data(fk_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyyy = cbuffer.data(fk_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyyz = cbuffer.data(fk_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyzz = cbuffer.data(fk_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyzzz = cbuffer.data(fk_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_xxz_xyyzzzz = cbuffer.data(fk_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_xxz_xyzzzzz = cbuffer.data(fk_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_xxz_xzzzzzz = cbuffer.data(fk_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyyy = cbuffer.data(fk_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyyz = cbuffer.data(fk_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyzz = cbuffer.data(fk_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyzzz = cbuffer.data(fk_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_xxz_yyyzzzz = cbuffer.data(fk_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_xxz_yyzzzzz = cbuffer.data(fk_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_xxz_yzzzzzz = cbuffer.data(fk_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_xxz_zzzzzzz = cbuffer.data(fk_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_xyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_xyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_xyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_xyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_xyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_z_xyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_z_xyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_xyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_xyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_xyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_xyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_xyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_xyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_xyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_xyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_xyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_xyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_xyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_xyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_xyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_xzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxx = cbuffer.data(fk_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxy = cbuffer.data(fk_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxz = cbuffer.data(fk_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxyy = cbuffer.data(fk_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxyz = cbuffer.data(fk_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxzz = cbuffer.data(fk_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyyy = cbuffer.data(fk_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyyz = cbuffer.data(fk_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyzz = cbuffer.data(fk_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxzzz = cbuffer.data(fk_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyyy = cbuffer.data(fk_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyyz = cbuffer.data(fk_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyzz = cbuffer.data(fk_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyzzz = cbuffer.data(fk_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_z_yyy_xxxzzzz = cbuffer.data(fk_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyyy = cbuffer.data(fk_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyyz = cbuffer.data(fk_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyzz = cbuffer.data(fk_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyzzz = cbuffer.data(fk_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_yyy_xxyzzzz = cbuffer.data(fk_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_yyy_xxzzzzz = cbuffer.data(fk_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyyy = cbuffer.data(fk_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyyz = cbuffer.data(fk_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyzz = cbuffer.data(fk_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyzzz = cbuffer.data(fk_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_yyy_xyyzzzz = cbuffer.data(fk_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_yyy_xyzzzzz = cbuffer.data(fk_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_yyy_xzzzzzz = cbuffer.data(fk_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyyy = cbuffer.data(fk_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyyz = cbuffer.data(fk_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyzz = cbuffer.data(fk_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyzzz = cbuffer.data(fk_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_yyy_yyyzzzz = cbuffer.data(fk_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_yyy_yyzzzzz = cbuffer.data(fk_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_z_yyy_yzzzzzz = cbuffer.data(fk_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_yyy_zzzzzzz = cbuffer.data(fk_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxx = cbuffer.data(fk_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxy = cbuffer.data(fk_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxz = cbuffer.data(fk_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxyy = cbuffer.data(fk_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxyz = cbuffer.data(fk_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxzz = cbuffer.data(fk_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyyy = cbuffer.data(fk_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyyz = cbuffer.data(fk_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyzz = cbuffer.data(fk_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxzzz = cbuffer.data(fk_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyyy = cbuffer.data(fk_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyyz = cbuffer.data(fk_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyzz = cbuffer.data(fk_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyzzz = cbuffer.data(fk_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_yyz_xxxzzzz = cbuffer.data(fk_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyyy = cbuffer.data(fk_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyyz = cbuffer.data(fk_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyzz = cbuffer.data(fk_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyzzz = cbuffer.data(fk_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_yyz_xxyzzzz = cbuffer.data(fk_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_yyz_xxzzzzz = cbuffer.data(fk_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyyy = cbuffer.data(fk_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyyz = cbuffer.data(fk_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyzz = cbuffer.data(fk_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyzzz = cbuffer.data(fk_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_yyz_xyyzzzz = cbuffer.data(fk_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_yyz_xyzzzzz = cbuffer.data(fk_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_yyz_xzzzzzz = cbuffer.data(fk_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyyy = cbuffer.data(fk_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyyz = cbuffer.data(fk_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyzz = cbuffer.data(fk_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyzzz = cbuffer.data(fk_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_yyz_yyyzzzz = cbuffer.data(fk_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_z_yyz_yyzzzzz = cbuffer.data(fk_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_yyz_yzzzzzz = cbuffer.data(fk_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_yyz_zzzzzzz = cbuffer.data(fk_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_yzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_yzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_yzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_yzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_yzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_z_yzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_z_yzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_yzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_yzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_yzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxx = cbuffer.data(fk_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxy = cbuffer.data(fk_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxz = cbuffer.data(fk_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxyy = cbuffer.data(fk_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxyz = cbuffer.data(fk_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxzz = cbuffer.data(fk_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyyy = cbuffer.data(fk_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyyz = cbuffer.data(fk_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyzz = cbuffer.data(fk_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxzzz = cbuffer.data(fk_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyyy = cbuffer.data(fk_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyyz = cbuffer.data(fk_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyzz = cbuffer.data(fk_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyzzz = cbuffer.data(fk_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_zzz_xxxzzzz = cbuffer.data(fk_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyyy = cbuffer.data(fk_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyyz = cbuffer.data(fk_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyzz = cbuffer.data(fk_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyzzz = cbuffer.data(fk_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_zzz_xxyzzzz = cbuffer.data(fk_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_zzz_xxzzzzz = cbuffer.data(fk_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyyy = cbuffer.data(fk_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyyz = cbuffer.data(fk_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyzz = cbuffer.data(fk_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyzzz = cbuffer.data(fk_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_zzz_xyyzzzz = cbuffer.data(fk_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_z_zzz_xyzzzzz = cbuffer.data(fk_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_z_zzz_xzzzzzz = cbuffer.data(fk_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyyy = cbuffer.data(fk_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyyz = cbuffer.data(fk_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyzz = cbuffer.data(fk_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyzzz = cbuffer.data(fk_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_zzz_yyyzzzz = cbuffer.data(fk_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_zzz_yyzzzzz = cbuffer.data(fk_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_zzz_yzzzzzz = cbuffer.data(fk_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_zzz_zzzzzzz = cbuffer.data(fk_geom_01_off + 1079 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FLSS

            const auto fl_geom_01_off = idx_geom_01_flxx + i * dcomps + j;

            auto g_0_x_xxx_xxxxxxxx = cbuffer.data(fl_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxxz = cbuffer.data(fl_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxyy = cbuffer.data(fl_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxyz = cbuffer.data(fl_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxxzz = cbuffer.data(fl_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxyyy = cbuffer.data(fl_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxyyz = cbuffer.data(fl_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxyzz = cbuffer.data(fl_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxzzz = cbuffer.data(fl_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyyyy = cbuffer.data(fl_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyyyz = cbuffer.data(fl_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyyzz = cbuffer.data(fl_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyzzz = cbuffer.data(fl_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxzzzz = cbuffer.data(fl_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyyyy = cbuffer.data(fl_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyyyz = cbuffer.data(fl_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyyzz = cbuffer.data(fl_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyzzz = cbuffer.data(fl_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyzzzz = cbuffer.data(fl_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxx_xxxzzzzz = cbuffer.data(fl_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyyyy = cbuffer.data(fl_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyyyz = cbuffer.data(fl_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyyzz = cbuffer.data(fl_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyzzz = cbuffer.data(fl_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyzzzz = cbuffer.data(fl_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxx_xxyzzzzz = cbuffer.data(fl_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxx_xxzzzzzz = cbuffer.data(fl_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyyyy = cbuffer.data(fl_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyyyz = cbuffer.data(fl_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyyzz = cbuffer.data(fl_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyzzz = cbuffer.data(fl_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyzzzz = cbuffer.data(fl_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxx_xyyzzzzz = cbuffer.data(fl_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxx_xyzzzzzz = cbuffer.data(fl_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxx_xzzzzzzz = cbuffer.data(fl_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyyyy = cbuffer.data(fl_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyyyz = cbuffer.data(fl_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyyzz = cbuffer.data(fl_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyzzz = cbuffer.data(fl_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyzzzz = cbuffer.data(fl_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxx_yyyzzzzz = cbuffer.data(fl_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxx_yyzzzzzz = cbuffer.data(fl_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxx_yzzzzzzz = cbuffer.data(fl_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxx_zzzzzzzz = cbuffer.data(fl_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_yyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_yyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_yyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_yyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_yyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_yyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_yyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_yyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_yyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_yyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_yyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_yyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_yyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_yyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_yyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_yyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_yyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_yyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_yyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_yyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_yzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_yzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_yzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_yzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_yzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_yzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_yzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_yzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_yzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_yzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_zzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_zzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_zzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_zzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_zzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_zzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_zzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_zzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_zzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_zzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxxx = cbuffer.data(fl_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxxy = cbuffer.data(fl_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxxz = cbuffer.data(fl_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxyy = cbuffer.data(fl_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxyz = cbuffer.data(fl_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxxzz = cbuffer.data(fl_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxyyy = cbuffer.data(fl_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxyyz = cbuffer.data(fl_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxyzz = cbuffer.data(fl_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxzzz = cbuffer.data(fl_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyyyy = cbuffer.data(fl_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyyyz = cbuffer.data(fl_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyyzz = cbuffer.data(fl_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyzzz = cbuffer.data(fl_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxzzzz = cbuffer.data(fl_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyyyy = cbuffer.data(fl_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyyyz = cbuffer.data(fl_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyyzz = cbuffer.data(fl_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyzzz = cbuffer.data(fl_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyzzzz = cbuffer.data(fl_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xxx_xxxzzzzz = cbuffer.data(fl_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyyyy = cbuffer.data(fl_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyyyz = cbuffer.data(fl_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyyzz = cbuffer.data(fl_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyzzz = cbuffer.data(fl_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyzzzz = cbuffer.data(fl_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xxx_xxyzzzzz = cbuffer.data(fl_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xxx_xxzzzzzz = cbuffer.data(fl_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyyyy = cbuffer.data(fl_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyyyz = cbuffer.data(fl_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyyzz = cbuffer.data(fl_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyzzz = cbuffer.data(fl_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyzzzz = cbuffer.data(fl_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xxx_xyyzzzzz = cbuffer.data(fl_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xxx_xyzzzzzz = cbuffer.data(fl_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xxx_xzzzzzzz = cbuffer.data(fl_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyyyy = cbuffer.data(fl_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyyyz = cbuffer.data(fl_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyyzz = cbuffer.data(fl_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyzzz = cbuffer.data(fl_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyzzzz = cbuffer.data(fl_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xxx_yyyzzzzz = cbuffer.data(fl_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xxx_yyzzzzzz = cbuffer.data(fl_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xxx_yzzzzzzz = cbuffer.data(fl_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xxx_zzzzzzzz = cbuffer.data(fl_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xxy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xxy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xxy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xxy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xxy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xxy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xxy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xxy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xxy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xxy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_yyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_yyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_yyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_yyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_yyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_yyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_yyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_yyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_yyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_yyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_yyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_yyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_yyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_yyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_yyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_yyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_yyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_yyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_yyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_yyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_yzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_yzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_yzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_yzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_yzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_yzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_yzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_yzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_yzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_yzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_zzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_zzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_zzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_zzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_zzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_zzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_zzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_zzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_zzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_zzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxxx = cbuffer.data(fl_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxxy = cbuffer.data(fl_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxxz = cbuffer.data(fl_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxyy = cbuffer.data(fl_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxyz = cbuffer.data(fl_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxxzz = cbuffer.data(fl_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxyyy = cbuffer.data(fl_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxyyz = cbuffer.data(fl_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxyzz = cbuffer.data(fl_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxzzz = cbuffer.data(fl_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyyyy = cbuffer.data(fl_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyyyz = cbuffer.data(fl_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyyzz = cbuffer.data(fl_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyzzz = cbuffer.data(fl_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxzzzz = cbuffer.data(fl_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyyyy = cbuffer.data(fl_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyyyz = cbuffer.data(fl_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyyzz = cbuffer.data(fl_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyzzz = cbuffer.data(fl_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyzzzz = cbuffer.data(fl_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_xxx_xxxzzzzz = cbuffer.data(fl_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyyyy = cbuffer.data(fl_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyyyz = cbuffer.data(fl_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyyzz = cbuffer.data(fl_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyzzz = cbuffer.data(fl_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyzzzz = cbuffer.data(fl_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xxx_xxyzzzzz = cbuffer.data(fl_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xxx_xxzzzzzz = cbuffer.data(fl_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyyyy = cbuffer.data(fl_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyyyz = cbuffer.data(fl_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyyzz = cbuffer.data(fl_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyzzz = cbuffer.data(fl_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyzzzz = cbuffer.data(fl_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xxx_xyyzzzzz = cbuffer.data(fl_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xxx_xyzzzzzz = cbuffer.data(fl_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xxx_xzzzzzzz = cbuffer.data(fl_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyyyy = cbuffer.data(fl_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyyyz = cbuffer.data(fl_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyyzz = cbuffer.data(fl_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyzzz = cbuffer.data(fl_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyzzzz = cbuffer.data(fl_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_xxx_yyyzzzzz = cbuffer.data(fl_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_xxx_yyzzzzzz = cbuffer.data(fl_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_xxx_yzzzzzzz = cbuffer.data(fl_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_xxx_zzzzzzzz = cbuffer.data(fl_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_xxy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_xxy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_xxy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_xxy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_xxy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_z_xxy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_xxy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_z_xxy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_xxy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_xxy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_z_xxz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_xxz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_xxz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_xxz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_xxz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_xxz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_xxz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_xxz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_xxz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_xxz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_xyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_xyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_xyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_xyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_xyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_z_xyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_xyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_xyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_xyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_xyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_z_xyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_xyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_xyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_xyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_xyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_xzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_xzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_xzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_xzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_xzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_xzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_xzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_xzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_xzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_xzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_yyy_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_yyy_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_yyy_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_yyy_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_yyy_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_yyy_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_yyy_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_yyy_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_yyy_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_yyy_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_yyz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_yyz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_yyz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_yyz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_yyz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_yyz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_yyz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_yyz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_yyz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_yyz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1259 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_z_yzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_z_yzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_z_yzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_z_yzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_z_yzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_z_yzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_z_yzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_z_yzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_z_yzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_z_yzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxxx = cbuffer.data(fl_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxxy = cbuffer.data(fl_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxxz = cbuffer.data(fl_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxyy = cbuffer.data(fl_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxyz = cbuffer.data(fl_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxxzz = cbuffer.data(fl_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxyyy = cbuffer.data(fl_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxyyz = cbuffer.data(fl_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxyzz = cbuffer.data(fl_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxzzz = cbuffer.data(fl_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyyyy = cbuffer.data(fl_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyyyz = cbuffer.data(fl_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyyzz = cbuffer.data(fl_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyzzz = cbuffer.data(fl_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxzzzz = cbuffer.data(fl_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyyyy = cbuffer.data(fl_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyyyz = cbuffer.data(fl_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyyzz = cbuffer.data(fl_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyzzz = cbuffer.data(fl_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyzzzz = cbuffer.data(fl_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_z_zzz_xxxzzzzz = cbuffer.data(fl_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyyyy = cbuffer.data(fl_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyyyz = cbuffer.data(fl_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyyzz = cbuffer.data(fl_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyzzz = cbuffer.data(fl_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyzzzz = cbuffer.data(fl_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_z_zzz_xxyzzzzz = cbuffer.data(fl_geom_01_off + 1331 * ccomps * dcomps);

            auto g_0_z_zzz_xxzzzzzz = cbuffer.data(fl_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyyyy = cbuffer.data(fl_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyyyz = cbuffer.data(fl_geom_01_off + 1334 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyyzz = cbuffer.data(fl_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyzzz = cbuffer.data(fl_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyzzzz = cbuffer.data(fl_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_z_zzz_xyyzzzzz = cbuffer.data(fl_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_z_zzz_xyzzzzzz = cbuffer.data(fl_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_z_zzz_xzzzzzzz = cbuffer.data(fl_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyyyy = cbuffer.data(fl_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyyyz = cbuffer.data(fl_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyyzz = cbuffer.data(fl_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyzzz = cbuffer.data(fl_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyzzzz = cbuffer.data(fl_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_z_zzz_yyyzzzzz = cbuffer.data(fl_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_z_zzz_yyzzzzzz = cbuffer.data(fl_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_z_zzz_yzzzzzzz = cbuffer.data(fl_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_z_zzz_zzzzzzzz = cbuffer.data(fl_geom_01_off + 1349 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gkxx

            const auto gk_geom_01_off = idx_geom_01_gkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxxxx, g_0_x_xxx_xxxxxxxx, g_0_x_xxx_xxxxxxxy, g_0_x_xxx_xxxxxxxz, g_0_x_xxx_xxxxxxy, g_0_x_xxx_xxxxxxyy, g_0_x_xxx_xxxxxxyz, g_0_x_xxx_xxxxxxz, g_0_x_xxx_xxxxxxzz, g_0_x_xxx_xxxxxyy, g_0_x_xxx_xxxxxyyy, g_0_x_xxx_xxxxxyyz, g_0_x_xxx_xxxxxyz, g_0_x_xxx_xxxxxyzz, g_0_x_xxx_xxxxxzz, g_0_x_xxx_xxxxxzzz, g_0_x_xxx_xxxxyyy, g_0_x_xxx_xxxxyyyy, g_0_x_xxx_xxxxyyyz, g_0_x_xxx_xxxxyyz, g_0_x_xxx_xxxxyyzz, g_0_x_xxx_xxxxyzz, g_0_x_xxx_xxxxyzzz, g_0_x_xxx_xxxxzzz, g_0_x_xxx_xxxxzzzz, g_0_x_xxx_xxxyyyy, g_0_x_xxx_xxxyyyyy, g_0_x_xxx_xxxyyyyz, g_0_x_xxx_xxxyyyz, g_0_x_xxx_xxxyyyzz, g_0_x_xxx_xxxyyzz, g_0_x_xxx_xxxyyzzz, g_0_x_xxx_xxxyzzz, g_0_x_xxx_xxxyzzzz, g_0_x_xxx_xxxzzzz, g_0_x_xxx_xxxzzzzz, g_0_x_xxx_xxyyyyy, g_0_x_xxx_xxyyyyyy, g_0_x_xxx_xxyyyyyz, g_0_x_xxx_xxyyyyz, g_0_x_xxx_xxyyyyzz, g_0_x_xxx_xxyyyzz, g_0_x_xxx_xxyyyzzz, g_0_x_xxx_xxyyzzz, g_0_x_xxx_xxyyzzzz, g_0_x_xxx_xxyzzzz, g_0_x_xxx_xxyzzzzz, g_0_x_xxx_xxzzzzz, g_0_x_xxx_xxzzzzzz, g_0_x_xxx_xyyyyyy, g_0_x_xxx_xyyyyyyy, g_0_x_xxx_xyyyyyyz, g_0_x_xxx_xyyyyyz, g_0_x_xxx_xyyyyyzz, g_0_x_xxx_xyyyyzz, g_0_x_xxx_xyyyyzzz, g_0_x_xxx_xyyyzzz, g_0_x_xxx_xyyyzzzz, g_0_x_xxx_xyyzzzz, g_0_x_xxx_xyyzzzzz, g_0_x_xxx_xyzzzzz, g_0_x_xxx_xyzzzzzz, g_0_x_xxx_xzzzzzz, g_0_x_xxx_xzzzzzzz, g_0_x_xxx_yyyyyyy, g_0_x_xxx_yyyyyyz, g_0_x_xxx_yyyyyzz, g_0_x_xxx_yyyyzzz, g_0_x_xxx_yyyzzzz, g_0_x_xxx_yyzzzzz, g_0_x_xxx_yzzzzzz, g_0_x_xxx_zzzzzzz, g_0_x_xxxx_xxxxxxx, g_0_x_xxxx_xxxxxxy, g_0_x_xxxx_xxxxxxz, g_0_x_xxxx_xxxxxyy, g_0_x_xxxx_xxxxxyz, g_0_x_xxxx_xxxxxzz, g_0_x_xxxx_xxxxyyy, g_0_x_xxxx_xxxxyyz, g_0_x_xxxx_xxxxyzz, g_0_x_xxxx_xxxxzzz, g_0_x_xxxx_xxxyyyy, g_0_x_xxxx_xxxyyyz, g_0_x_xxxx_xxxyyzz, g_0_x_xxxx_xxxyzzz, g_0_x_xxxx_xxxzzzz, g_0_x_xxxx_xxyyyyy, g_0_x_xxxx_xxyyyyz, g_0_x_xxxx_xxyyyzz, g_0_x_xxxx_xxyyzzz, g_0_x_xxxx_xxyzzzz, g_0_x_xxxx_xxzzzzz, g_0_x_xxxx_xyyyyyy, g_0_x_xxxx_xyyyyyz, g_0_x_xxxx_xyyyyzz, g_0_x_xxxx_xyyyzzz, g_0_x_xxxx_xyyzzzz, g_0_x_xxxx_xyzzzzz, g_0_x_xxxx_xzzzzzz, g_0_x_xxxx_yyyyyyy, g_0_x_xxxx_yyyyyyz, g_0_x_xxxx_yyyyyzz, g_0_x_xxxx_yyyyzzz, g_0_x_xxxx_yyyzzzz, g_0_x_xxxx_yyzzzzz, g_0_x_xxxx_yzzzzzz, g_0_x_xxxx_zzzzzzz, g_xxx_xxxxxxx, g_xxx_xxxxxxy, g_xxx_xxxxxxz, g_xxx_xxxxxyy, g_xxx_xxxxxyz, g_xxx_xxxxxzz, g_xxx_xxxxyyy, g_xxx_xxxxyyz, g_xxx_xxxxyzz, g_xxx_xxxxzzz, g_xxx_xxxyyyy, g_xxx_xxxyyyz, g_xxx_xxxyyzz, g_xxx_xxxyzzz, g_xxx_xxxzzzz, g_xxx_xxyyyyy, g_xxx_xxyyyyz, g_xxx_xxyyyzz, g_xxx_xxyyzzz, g_xxx_xxyzzzz, g_xxx_xxzzzzz, g_xxx_xyyyyyy, g_xxx_xyyyyyz, g_xxx_xyyyyzz, g_xxx_xyyyzzz, g_xxx_xyyzzzz, g_xxx_xyzzzzz, g_xxx_xzzzzzz, g_xxx_yyyyyyy, g_xxx_yyyyyyz, g_xxx_yyyyyzz, g_xxx_yyyyzzz, g_xxx_yyyzzzz, g_xxx_yyzzzzz, g_xxx_yzzzzzz, g_xxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_xxxxxxx[k] = g_xxx_xxxxxxx[k] - g_0_x_xxx_xxxxxxx[k] * ab_x + g_0_x_xxx_xxxxxxxx[k];

                g_0_x_xxxx_xxxxxxy[k] = g_xxx_xxxxxxy[k] - g_0_x_xxx_xxxxxxy[k] * ab_x + g_0_x_xxx_xxxxxxxy[k];

                g_0_x_xxxx_xxxxxxz[k] = g_xxx_xxxxxxz[k] - g_0_x_xxx_xxxxxxz[k] * ab_x + g_0_x_xxx_xxxxxxxz[k];

                g_0_x_xxxx_xxxxxyy[k] = g_xxx_xxxxxyy[k] - g_0_x_xxx_xxxxxyy[k] * ab_x + g_0_x_xxx_xxxxxxyy[k];

                g_0_x_xxxx_xxxxxyz[k] = g_xxx_xxxxxyz[k] - g_0_x_xxx_xxxxxyz[k] * ab_x + g_0_x_xxx_xxxxxxyz[k];

                g_0_x_xxxx_xxxxxzz[k] = g_xxx_xxxxxzz[k] - g_0_x_xxx_xxxxxzz[k] * ab_x + g_0_x_xxx_xxxxxxzz[k];

                g_0_x_xxxx_xxxxyyy[k] = g_xxx_xxxxyyy[k] - g_0_x_xxx_xxxxyyy[k] * ab_x + g_0_x_xxx_xxxxxyyy[k];

                g_0_x_xxxx_xxxxyyz[k] = g_xxx_xxxxyyz[k] - g_0_x_xxx_xxxxyyz[k] * ab_x + g_0_x_xxx_xxxxxyyz[k];

                g_0_x_xxxx_xxxxyzz[k] = g_xxx_xxxxyzz[k] - g_0_x_xxx_xxxxyzz[k] * ab_x + g_0_x_xxx_xxxxxyzz[k];

                g_0_x_xxxx_xxxxzzz[k] = g_xxx_xxxxzzz[k] - g_0_x_xxx_xxxxzzz[k] * ab_x + g_0_x_xxx_xxxxxzzz[k];

                g_0_x_xxxx_xxxyyyy[k] = g_xxx_xxxyyyy[k] - g_0_x_xxx_xxxyyyy[k] * ab_x + g_0_x_xxx_xxxxyyyy[k];

                g_0_x_xxxx_xxxyyyz[k] = g_xxx_xxxyyyz[k] - g_0_x_xxx_xxxyyyz[k] * ab_x + g_0_x_xxx_xxxxyyyz[k];

                g_0_x_xxxx_xxxyyzz[k] = g_xxx_xxxyyzz[k] - g_0_x_xxx_xxxyyzz[k] * ab_x + g_0_x_xxx_xxxxyyzz[k];

                g_0_x_xxxx_xxxyzzz[k] = g_xxx_xxxyzzz[k] - g_0_x_xxx_xxxyzzz[k] * ab_x + g_0_x_xxx_xxxxyzzz[k];

                g_0_x_xxxx_xxxzzzz[k] = g_xxx_xxxzzzz[k] - g_0_x_xxx_xxxzzzz[k] * ab_x + g_0_x_xxx_xxxxzzzz[k];

                g_0_x_xxxx_xxyyyyy[k] = g_xxx_xxyyyyy[k] - g_0_x_xxx_xxyyyyy[k] * ab_x + g_0_x_xxx_xxxyyyyy[k];

                g_0_x_xxxx_xxyyyyz[k] = g_xxx_xxyyyyz[k] - g_0_x_xxx_xxyyyyz[k] * ab_x + g_0_x_xxx_xxxyyyyz[k];

                g_0_x_xxxx_xxyyyzz[k] = g_xxx_xxyyyzz[k] - g_0_x_xxx_xxyyyzz[k] * ab_x + g_0_x_xxx_xxxyyyzz[k];

                g_0_x_xxxx_xxyyzzz[k] = g_xxx_xxyyzzz[k] - g_0_x_xxx_xxyyzzz[k] * ab_x + g_0_x_xxx_xxxyyzzz[k];

                g_0_x_xxxx_xxyzzzz[k] = g_xxx_xxyzzzz[k] - g_0_x_xxx_xxyzzzz[k] * ab_x + g_0_x_xxx_xxxyzzzz[k];

                g_0_x_xxxx_xxzzzzz[k] = g_xxx_xxzzzzz[k] - g_0_x_xxx_xxzzzzz[k] * ab_x + g_0_x_xxx_xxxzzzzz[k];

                g_0_x_xxxx_xyyyyyy[k] = g_xxx_xyyyyyy[k] - g_0_x_xxx_xyyyyyy[k] * ab_x + g_0_x_xxx_xxyyyyyy[k];

                g_0_x_xxxx_xyyyyyz[k] = g_xxx_xyyyyyz[k] - g_0_x_xxx_xyyyyyz[k] * ab_x + g_0_x_xxx_xxyyyyyz[k];

                g_0_x_xxxx_xyyyyzz[k] = g_xxx_xyyyyzz[k] - g_0_x_xxx_xyyyyzz[k] * ab_x + g_0_x_xxx_xxyyyyzz[k];

                g_0_x_xxxx_xyyyzzz[k] = g_xxx_xyyyzzz[k] - g_0_x_xxx_xyyyzzz[k] * ab_x + g_0_x_xxx_xxyyyzzz[k];

                g_0_x_xxxx_xyyzzzz[k] = g_xxx_xyyzzzz[k] - g_0_x_xxx_xyyzzzz[k] * ab_x + g_0_x_xxx_xxyyzzzz[k];

                g_0_x_xxxx_xyzzzzz[k] = g_xxx_xyzzzzz[k] - g_0_x_xxx_xyzzzzz[k] * ab_x + g_0_x_xxx_xxyzzzzz[k];

                g_0_x_xxxx_xzzzzzz[k] = g_xxx_xzzzzzz[k] - g_0_x_xxx_xzzzzzz[k] * ab_x + g_0_x_xxx_xxzzzzzz[k];

                g_0_x_xxxx_yyyyyyy[k] = g_xxx_yyyyyyy[k] - g_0_x_xxx_yyyyyyy[k] * ab_x + g_0_x_xxx_xyyyyyyy[k];

                g_0_x_xxxx_yyyyyyz[k] = g_xxx_yyyyyyz[k] - g_0_x_xxx_yyyyyyz[k] * ab_x + g_0_x_xxx_xyyyyyyz[k];

                g_0_x_xxxx_yyyyyzz[k] = g_xxx_yyyyyzz[k] - g_0_x_xxx_yyyyyzz[k] * ab_x + g_0_x_xxx_xyyyyyzz[k];

                g_0_x_xxxx_yyyyzzz[k] = g_xxx_yyyyzzz[k] - g_0_x_xxx_yyyyzzz[k] * ab_x + g_0_x_xxx_xyyyyzzz[k];

                g_0_x_xxxx_yyyzzzz[k] = g_xxx_yyyzzzz[k] - g_0_x_xxx_yyyzzzz[k] * ab_x + g_0_x_xxx_xyyyzzzz[k];

                g_0_x_xxxx_yyzzzzz[k] = g_xxx_yyzzzzz[k] - g_0_x_xxx_yyzzzzz[k] * ab_x + g_0_x_xxx_xyyzzzzz[k];

                g_0_x_xxxx_yzzzzzz[k] = g_xxx_yzzzzzz[k] - g_0_x_xxx_yzzzzzz[k] * ab_x + g_0_x_xxx_xyzzzzzz[k];

                g_0_x_xxxx_zzzzzzz[k] = g_xxx_zzzzzzz[k] - g_0_x_xxx_zzzzzzz[k] * ab_x + g_0_x_xxx_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxxxx, g_0_x_xxx_xxxxxxxy, g_0_x_xxx_xxxxxxy, g_0_x_xxx_xxxxxxyy, g_0_x_xxx_xxxxxxyz, g_0_x_xxx_xxxxxxz, g_0_x_xxx_xxxxxyy, g_0_x_xxx_xxxxxyyy, g_0_x_xxx_xxxxxyyz, g_0_x_xxx_xxxxxyz, g_0_x_xxx_xxxxxyzz, g_0_x_xxx_xxxxxzz, g_0_x_xxx_xxxxyyy, g_0_x_xxx_xxxxyyyy, g_0_x_xxx_xxxxyyyz, g_0_x_xxx_xxxxyyz, g_0_x_xxx_xxxxyyzz, g_0_x_xxx_xxxxyzz, g_0_x_xxx_xxxxyzzz, g_0_x_xxx_xxxxzzz, g_0_x_xxx_xxxyyyy, g_0_x_xxx_xxxyyyyy, g_0_x_xxx_xxxyyyyz, g_0_x_xxx_xxxyyyz, g_0_x_xxx_xxxyyyzz, g_0_x_xxx_xxxyyzz, g_0_x_xxx_xxxyyzzz, g_0_x_xxx_xxxyzzz, g_0_x_xxx_xxxyzzzz, g_0_x_xxx_xxxzzzz, g_0_x_xxx_xxyyyyy, g_0_x_xxx_xxyyyyyy, g_0_x_xxx_xxyyyyyz, g_0_x_xxx_xxyyyyz, g_0_x_xxx_xxyyyyzz, g_0_x_xxx_xxyyyzz, g_0_x_xxx_xxyyyzzz, g_0_x_xxx_xxyyzzz, g_0_x_xxx_xxyyzzzz, g_0_x_xxx_xxyzzzz, g_0_x_xxx_xxyzzzzz, g_0_x_xxx_xxzzzzz, g_0_x_xxx_xyyyyyy, g_0_x_xxx_xyyyyyyy, g_0_x_xxx_xyyyyyyz, g_0_x_xxx_xyyyyyz, g_0_x_xxx_xyyyyyzz, g_0_x_xxx_xyyyyzz, g_0_x_xxx_xyyyyzzz, g_0_x_xxx_xyyyzzz, g_0_x_xxx_xyyyzzzz, g_0_x_xxx_xyyzzzz, g_0_x_xxx_xyyzzzzz, g_0_x_xxx_xyzzzzz, g_0_x_xxx_xyzzzzzz, g_0_x_xxx_xzzzzzz, g_0_x_xxx_yyyyyyy, g_0_x_xxx_yyyyyyyy, g_0_x_xxx_yyyyyyyz, g_0_x_xxx_yyyyyyz, g_0_x_xxx_yyyyyyzz, g_0_x_xxx_yyyyyzz, g_0_x_xxx_yyyyyzzz, g_0_x_xxx_yyyyzzz, g_0_x_xxx_yyyyzzzz, g_0_x_xxx_yyyzzzz, g_0_x_xxx_yyyzzzzz, g_0_x_xxx_yyzzzzz, g_0_x_xxx_yyzzzzzz, g_0_x_xxx_yzzzzzz, g_0_x_xxx_yzzzzzzz, g_0_x_xxx_zzzzzzz, g_0_x_xxxy_xxxxxxx, g_0_x_xxxy_xxxxxxy, g_0_x_xxxy_xxxxxxz, g_0_x_xxxy_xxxxxyy, g_0_x_xxxy_xxxxxyz, g_0_x_xxxy_xxxxxzz, g_0_x_xxxy_xxxxyyy, g_0_x_xxxy_xxxxyyz, g_0_x_xxxy_xxxxyzz, g_0_x_xxxy_xxxxzzz, g_0_x_xxxy_xxxyyyy, g_0_x_xxxy_xxxyyyz, g_0_x_xxxy_xxxyyzz, g_0_x_xxxy_xxxyzzz, g_0_x_xxxy_xxxzzzz, g_0_x_xxxy_xxyyyyy, g_0_x_xxxy_xxyyyyz, g_0_x_xxxy_xxyyyzz, g_0_x_xxxy_xxyyzzz, g_0_x_xxxy_xxyzzzz, g_0_x_xxxy_xxzzzzz, g_0_x_xxxy_xyyyyyy, g_0_x_xxxy_xyyyyyz, g_0_x_xxxy_xyyyyzz, g_0_x_xxxy_xyyyzzz, g_0_x_xxxy_xyyzzzz, g_0_x_xxxy_xyzzzzz, g_0_x_xxxy_xzzzzzz, g_0_x_xxxy_yyyyyyy, g_0_x_xxxy_yyyyyyz, g_0_x_xxxy_yyyyyzz, g_0_x_xxxy_yyyyzzz, g_0_x_xxxy_yyyzzzz, g_0_x_xxxy_yyzzzzz, g_0_x_xxxy_yzzzzzz, g_0_x_xxxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_xxxxxxx[k] = -g_0_x_xxx_xxxxxxx[k] * ab_y + g_0_x_xxx_xxxxxxxy[k];

                g_0_x_xxxy_xxxxxxy[k] = -g_0_x_xxx_xxxxxxy[k] * ab_y + g_0_x_xxx_xxxxxxyy[k];

                g_0_x_xxxy_xxxxxxz[k] = -g_0_x_xxx_xxxxxxz[k] * ab_y + g_0_x_xxx_xxxxxxyz[k];

                g_0_x_xxxy_xxxxxyy[k] = -g_0_x_xxx_xxxxxyy[k] * ab_y + g_0_x_xxx_xxxxxyyy[k];

                g_0_x_xxxy_xxxxxyz[k] = -g_0_x_xxx_xxxxxyz[k] * ab_y + g_0_x_xxx_xxxxxyyz[k];

                g_0_x_xxxy_xxxxxzz[k] = -g_0_x_xxx_xxxxxzz[k] * ab_y + g_0_x_xxx_xxxxxyzz[k];

                g_0_x_xxxy_xxxxyyy[k] = -g_0_x_xxx_xxxxyyy[k] * ab_y + g_0_x_xxx_xxxxyyyy[k];

                g_0_x_xxxy_xxxxyyz[k] = -g_0_x_xxx_xxxxyyz[k] * ab_y + g_0_x_xxx_xxxxyyyz[k];

                g_0_x_xxxy_xxxxyzz[k] = -g_0_x_xxx_xxxxyzz[k] * ab_y + g_0_x_xxx_xxxxyyzz[k];

                g_0_x_xxxy_xxxxzzz[k] = -g_0_x_xxx_xxxxzzz[k] * ab_y + g_0_x_xxx_xxxxyzzz[k];

                g_0_x_xxxy_xxxyyyy[k] = -g_0_x_xxx_xxxyyyy[k] * ab_y + g_0_x_xxx_xxxyyyyy[k];

                g_0_x_xxxy_xxxyyyz[k] = -g_0_x_xxx_xxxyyyz[k] * ab_y + g_0_x_xxx_xxxyyyyz[k];

                g_0_x_xxxy_xxxyyzz[k] = -g_0_x_xxx_xxxyyzz[k] * ab_y + g_0_x_xxx_xxxyyyzz[k];

                g_0_x_xxxy_xxxyzzz[k] = -g_0_x_xxx_xxxyzzz[k] * ab_y + g_0_x_xxx_xxxyyzzz[k];

                g_0_x_xxxy_xxxzzzz[k] = -g_0_x_xxx_xxxzzzz[k] * ab_y + g_0_x_xxx_xxxyzzzz[k];

                g_0_x_xxxy_xxyyyyy[k] = -g_0_x_xxx_xxyyyyy[k] * ab_y + g_0_x_xxx_xxyyyyyy[k];

                g_0_x_xxxy_xxyyyyz[k] = -g_0_x_xxx_xxyyyyz[k] * ab_y + g_0_x_xxx_xxyyyyyz[k];

                g_0_x_xxxy_xxyyyzz[k] = -g_0_x_xxx_xxyyyzz[k] * ab_y + g_0_x_xxx_xxyyyyzz[k];

                g_0_x_xxxy_xxyyzzz[k] = -g_0_x_xxx_xxyyzzz[k] * ab_y + g_0_x_xxx_xxyyyzzz[k];

                g_0_x_xxxy_xxyzzzz[k] = -g_0_x_xxx_xxyzzzz[k] * ab_y + g_0_x_xxx_xxyyzzzz[k];

                g_0_x_xxxy_xxzzzzz[k] = -g_0_x_xxx_xxzzzzz[k] * ab_y + g_0_x_xxx_xxyzzzzz[k];

                g_0_x_xxxy_xyyyyyy[k] = -g_0_x_xxx_xyyyyyy[k] * ab_y + g_0_x_xxx_xyyyyyyy[k];

                g_0_x_xxxy_xyyyyyz[k] = -g_0_x_xxx_xyyyyyz[k] * ab_y + g_0_x_xxx_xyyyyyyz[k];

                g_0_x_xxxy_xyyyyzz[k] = -g_0_x_xxx_xyyyyzz[k] * ab_y + g_0_x_xxx_xyyyyyzz[k];

                g_0_x_xxxy_xyyyzzz[k] = -g_0_x_xxx_xyyyzzz[k] * ab_y + g_0_x_xxx_xyyyyzzz[k];

                g_0_x_xxxy_xyyzzzz[k] = -g_0_x_xxx_xyyzzzz[k] * ab_y + g_0_x_xxx_xyyyzzzz[k];

                g_0_x_xxxy_xyzzzzz[k] = -g_0_x_xxx_xyzzzzz[k] * ab_y + g_0_x_xxx_xyyzzzzz[k];

                g_0_x_xxxy_xzzzzzz[k] = -g_0_x_xxx_xzzzzzz[k] * ab_y + g_0_x_xxx_xyzzzzzz[k];

                g_0_x_xxxy_yyyyyyy[k] = -g_0_x_xxx_yyyyyyy[k] * ab_y + g_0_x_xxx_yyyyyyyy[k];

                g_0_x_xxxy_yyyyyyz[k] = -g_0_x_xxx_yyyyyyz[k] * ab_y + g_0_x_xxx_yyyyyyyz[k];

                g_0_x_xxxy_yyyyyzz[k] = -g_0_x_xxx_yyyyyzz[k] * ab_y + g_0_x_xxx_yyyyyyzz[k];

                g_0_x_xxxy_yyyyzzz[k] = -g_0_x_xxx_yyyyzzz[k] * ab_y + g_0_x_xxx_yyyyyzzz[k];

                g_0_x_xxxy_yyyzzzz[k] = -g_0_x_xxx_yyyzzzz[k] * ab_y + g_0_x_xxx_yyyyzzzz[k];

                g_0_x_xxxy_yyzzzzz[k] = -g_0_x_xxx_yyzzzzz[k] * ab_y + g_0_x_xxx_yyyzzzzz[k];

                g_0_x_xxxy_yzzzzzz[k] = -g_0_x_xxx_yzzzzzz[k] * ab_y + g_0_x_xxx_yyzzzzzz[k];

                g_0_x_xxxy_zzzzzzz[k] = -g_0_x_xxx_zzzzzzz[k] * ab_y + g_0_x_xxx_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxxxx, g_0_x_xxx_xxxxxxxz, g_0_x_xxx_xxxxxxy, g_0_x_xxx_xxxxxxyz, g_0_x_xxx_xxxxxxz, g_0_x_xxx_xxxxxxzz, g_0_x_xxx_xxxxxyy, g_0_x_xxx_xxxxxyyz, g_0_x_xxx_xxxxxyz, g_0_x_xxx_xxxxxyzz, g_0_x_xxx_xxxxxzz, g_0_x_xxx_xxxxxzzz, g_0_x_xxx_xxxxyyy, g_0_x_xxx_xxxxyyyz, g_0_x_xxx_xxxxyyz, g_0_x_xxx_xxxxyyzz, g_0_x_xxx_xxxxyzz, g_0_x_xxx_xxxxyzzz, g_0_x_xxx_xxxxzzz, g_0_x_xxx_xxxxzzzz, g_0_x_xxx_xxxyyyy, g_0_x_xxx_xxxyyyyz, g_0_x_xxx_xxxyyyz, g_0_x_xxx_xxxyyyzz, g_0_x_xxx_xxxyyzz, g_0_x_xxx_xxxyyzzz, g_0_x_xxx_xxxyzzz, g_0_x_xxx_xxxyzzzz, g_0_x_xxx_xxxzzzz, g_0_x_xxx_xxxzzzzz, g_0_x_xxx_xxyyyyy, g_0_x_xxx_xxyyyyyz, g_0_x_xxx_xxyyyyz, g_0_x_xxx_xxyyyyzz, g_0_x_xxx_xxyyyzz, g_0_x_xxx_xxyyyzzz, g_0_x_xxx_xxyyzzz, g_0_x_xxx_xxyyzzzz, g_0_x_xxx_xxyzzzz, g_0_x_xxx_xxyzzzzz, g_0_x_xxx_xxzzzzz, g_0_x_xxx_xxzzzzzz, g_0_x_xxx_xyyyyyy, g_0_x_xxx_xyyyyyyz, g_0_x_xxx_xyyyyyz, g_0_x_xxx_xyyyyyzz, g_0_x_xxx_xyyyyzz, g_0_x_xxx_xyyyyzzz, g_0_x_xxx_xyyyzzz, g_0_x_xxx_xyyyzzzz, g_0_x_xxx_xyyzzzz, g_0_x_xxx_xyyzzzzz, g_0_x_xxx_xyzzzzz, g_0_x_xxx_xyzzzzzz, g_0_x_xxx_xzzzzzz, g_0_x_xxx_xzzzzzzz, g_0_x_xxx_yyyyyyy, g_0_x_xxx_yyyyyyyz, g_0_x_xxx_yyyyyyz, g_0_x_xxx_yyyyyyzz, g_0_x_xxx_yyyyyzz, g_0_x_xxx_yyyyyzzz, g_0_x_xxx_yyyyzzz, g_0_x_xxx_yyyyzzzz, g_0_x_xxx_yyyzzzz, g_0_x_xxx_yyyzzzzz, g_0_x_xxx_yyzzzzz, g_0_x_xxx_yyzzzzzz, g_0_x_xxx_yzzzzzz, g_0_x_xxx_yzzzzzzz, g_0_x_xxx_zzzzzzz, g_0_x_xxx_zzzzzzzz, g_0_x_xxxz_xxxxxxx, g_0_x_xxxz_xxxxxxy, g_0_x_xxxz_xxxxxxz, g_0_x_xxxz_xxxxxyy, g_0_x_xxxz_xxxxxyz, g_0_x_xxxz_xxxxxzz, g_0_x_xxxz_xxxxyyy, g_0_x_xxxz_xxxxyyz, g_0_x_xxxz_xxxxyzz, g_0_x_xxxz_xxxxzzz, g_0_x_xxxz_xxxyyyy, g_0_x_xxxz_xxxyyyz, g_0_x_xxxz_xxxyyzz, g_0_x_xxxz_xxxyzzz, g_0_x_xxxz_xxxzzzz, g_0_x_xxxz_xxyyyyy, g_0_x_xxxz_xxyyyyz, g_0_x_xxxz_xxyyyzz, g_0_x_xxxz_xxyyzzz, g_0_x_xxxz_xxyzzzz, g_0_x_xxxz_xxzzzzz, g_0_x_xxxz_xyyyyyy, g_0_x_xxxz_xyyyyyz, g_0_x_xxxz_xyyyyzz, g_0_x_xxxz_xyyyzzz, g_0_x_xxxz_xyyzzzz, g_0_x_xxxz_xyzzzzz, g_0_x_xxxz_xzzzzzz, g_0_x_xxxz_yyyyyyy, g_0_x_xxxz_yyyyyyz, g_0_x_xxxz_yyyyyzz, g_0_x_xxxz_yyyyzzz, g_0_x_xxxz_yyyzzzz, g_0_x_xxxz_yyzzzzz, g_0_x_xxxz_yzzzzzz, g_0_x_xxxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_xxxxxxx[k] = -g_0_x_xxx_xxxxxxx[k] * ab_z + g_0_x_xxx_xxxxxxxz[k];

                g_0_x_xxxz_xxxxxxy[k] = -g_0_x_xxx_xxxxxxy[k] * ab_z + g_0_x_xxx_xxxxxxyz[k];

                g_0_x_xxxz_xxxxxxz[k] = -g_0_x_xxx_xxxxxxz[k] * ab_z + g_0_x_xxx_xxxxxxzz[k];

                g_0_x_xxxz_xxxxxyy[k] = -g_0_x_xxx_xxxxxyy[k] * ab_z + g_0_x_xxx_xxxxxyyz[k];

                g_0_x_xxxz_xxxxxyz[k] = -g_0_x_xxx_xxxxxyz[k] * ab_z + g_0_x_xxx_xxxxxyzz[k];

                g_0_x_xxxz_xxxxxzz[k] = -g_0_x_xxx_xxxxxzz[k] * ab_z + g_0_x_xxx_xxxxxzzz[k];

                g_0_x_xxxz_xxxxyyy[k] = -g_0_x_xxx_xxxxyyy[k] * ab_z + g_0_x_xxx_xxxxyyyz[k];

                g_0_x_xxxz_xxxxyyz[k] = -g_0_x_xxx_xxxxyyz[k] * ab_z + g_0_x_xxx_xxxxyyzz[k];

                g_0_x_xxxz_xxxxyzz[k] = -g_0_x_xxx_xxxxyzz[k] * ab_z + g_0_x_xxx_xxxxyzzz[k];

                g_0_x_xxxz_xxxxzzz[k] = -g_0_x_xxx_xxxxzzz[k] * ab_z + g_0_x_xxx_xxxxzzzz[k];

                g_0_x_xxxz_xxxyyyy[k] = -g_0_x_xxx_xxxyyyy[k] * ab_z + g_0_x_xxx_xxxyyyyz[k];

                g_0_x_xxxz_xxxyyyz[k] = -g_0_x_xxx_xxxyyyz[k] * ab_z + g_0_x_xxx_xxxyyyzz[k];

                g_0_x_xxxz_xxxyyzz[k] = -g_0_x_xxx_xxxyyzz[k] * ab_z + g_0_x_xxx_xxxyyzzz[k];

                g_0_x_xxxz_xxxyzzz[k] = -g_0_x_xxx_xxxyzzz[k] * ab_z + g_0_x_xxx_xxxyzzzz[k];

                g_0_x_xxxz_xxxzzzz[k] = -g_0_x_xxx_xxxzzzz[k] * ab_z + g_0_x_xxx_xxxzzzzz[k];

                g_0_x_xxxz_xxyyyyy[k] = -g_0_x_xxx_xxyyyyy[k] * ab_z + g_0_x_xxx_xxyyyyyz[k];

                g_0_x_xxxz_xxyyyyz[k] = -g_0_x_xxx_xxyyyyz[k] * ab_z + g_0_x_xxx_xxyyyyzz[k];

                g_0_x_xxxz_xxyyyzz[k] = -g_0_x_xxx_xxyyyzz[k] * ab_z + g_0_x_xxx_xxyyyzzz[k];

                g_0_x_xxxz_xxyyzzz[k] = -g_0_x_xxx_xxyyzzz[k] * ab_z + g_0_x_xxx_xxyyzzzz[k];

                g_0_x_xxxz_xxyzzzz[k] = -g_0_x_xxx_xxyzzzz[k] * ab_z + g_0_x_xxx_xxyzzzzz[k];

                g_0_x_xxxz_xxzzzzz[k] = -g_0_x_xxx_xxzzzzz[k] * ab_z + g_0_x_xxx_xxzzzzzz[k];

                g_0_x_xxxz_xyyyyyy[k] = -g_0_x_xxx_xyyyyyy[k] * ab_z + g_0_x_xxx_xyyyyyyz[k];

                g_0_x_xxxz_xyyyyyz[k] = -g_0_x_xxx_xyyyyyz[k] * ab_z + g_0_x_xxx_xyyyyyzz[k];

                g_0_x_xxxz_xyyyyzz[k] = -g_0_x_xxx_xyyyyzz[k] * ab_z + g_0_x_xxx_xyyyyzzz[k];

                g_0_x_xxxz_xyyyzzz[k] = -g_0_x_xxx_xyyyzzz[k] * ab_z + g_0_x_xxx_xyyyzzzz[k];

                g_0_x_xxxz_xyyzzzz[k] = -g_0_x_xxx_xyyzzzz[k] * ab_z + g_0_x_xxx_xyyzzzzz[k];

                g_0_x_xxxz_xyzzzzz[k] = -g_0_x_xxx_xyzzzzz[k] * ab_z + g_0_x_xxx_xyzzzzzz[k];

                g_0_x_xxxz_xzzzzzz[k] = -g_0_x_xxx_xzzzzzz[k] * ab_z + g_0_x_xxx_xzzzzzzz[k];

                g_0_x_xxxz_yyyyyyy[k] = -g_0_x_xxx_yyyyyyy[k] * ab_z + g_0_x_xxx_yyyyyyyz[k];

                g_0_x_xxxz_yyyyyyz[k] = -g_0_x_xxx_yyyyyyz[k] * ab_z + g_0_x_xxx_yyyyyyzz[k];

                g_0_x_xxxz_yyyyyzz[k] = -g_0_x_xxx_yyyyyzz[k] * ab_z + g_0_x_xxx_yyyyyzzz[k];

                g_0_x_xxxz_yyyyzzz[k] = -g_0_x_xxx_yyyyzzz[k] * ab_z + g_0_x_xxx_yyyyzzzz[k];

                g_0_x_xxxz_yyyzzzz[k] = -g_0_x_xxx_yyyzzzz[k] * ab_z + g_0_x_xxx_yyyzzzzz[k];

                g_0_x_xxxz_yyzzzzz[k] = -g_0_x_xxx_yyzzzzz[k] * ab_z + g_0_x_xxx_yyzzzzzz[k];

                g_0_x_xxxz_yzzzzzz[k] = -g_0_x_xxx_yzzzzzz[k] * ab_z + g_0_x_xxx_yzzzzzzz[k];

                g_0_x_xxxz_zzzzzzz[k] = -g_0_x_xxx_zzzzzzz[k] * ab_z + g_0_x_xxx_zzzzzzzz[k];
            }

            /// Set up 108-144 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxy_xxxxxxx, g_0_x_xxy_xxxxxxxy, g_0_x_xxy_xxxxxxy, g_0_x_xxy_xxxxxxyy, g_0_x_xxy_xxxxxxyz, g_0_x_xxy_xxxxxxz, g_0_x_xxy_xxxxxyy, g_0_x_xxy_xxxxxyyy, g_0_x_xxy_xxxxxyyz, g_0_x_xxy_xxxxxyz, g_0_x_xxy_xxxxxyzz, g_0_x_xxy_xxxxxzz, g_0_x_xxy_xxxxyyy, g_0_x_xxy_xxxxyyyy, g_0_x_xxy_xxxxyyyz, g_0_x_xxy_xxxxyyz, g_0_x_xxy_xxxxyyzz, g_0_x_xxy_xxxxyzz, g_0_x_xxy_xxxxyzzz, g_0_x_xxy_xxxxzzz, g_0_x_xxy_xxxyyyy, g_0_x_xxy_xxxyyyyy, g_0_x_xxy_xxxyyyyz, g_0_x_xxy_xxxyyyz, g_0_x_xxy_xxxyyyzz, g_0_x_xxy_xxxyyzz, g_0_x_xxy_xxxyyzzz, g_0_x_xxy_xxxyzzz, g_0_x_xxy_xxxyzzzz, g_0_x_xxy_xxxzzzz, g_0_x_xxy_xxyyyyy, g_0_x_xxy_xxyyyyyy, g_0_x_xxy_xxyyyyyz, g_0_x_xxy_xxyyyyz, g_0_x_xxy_xxyyyyzz, g_0_x_xxy_xxyyyzz, g_0_x_xxy_xxyyyzzz, g_0_x_xxy_xxyyzzz, g_0_x_xxy_xxyyzzzz, g_0_x_xxy_xxyzzzz, g_0_x_xxy_xxyzzzzz, g_0_x_xxy_xxzzzzz, g_0_x_xxy_xyyyyyy, g_0_x_xxy_xyyyyyyy, g_0_x_xxy_xyyyyyyz, g_0_x_xxy_xyyyyyz, g_0_x_xxy_xyyyyyzz, g_0_x_xxy_xyyyyzz, g_0_x_xxy_xyyyyzzz, g_0_x_xxy_xyyyzzz, g_0_x_xxy_xyyyzzzz, g_0_x_xxy_xyyzzzz, g_0_x_xxy_xyyzzzzz, g_0_x_xxy_xyzzzzz, g_0_x_xxy_xyzzzzzz, g_0_x_xxy_xzzzzzz, g_0_x_xxy_yyyyyyy, g_0_x_xxy_yyyyyyyy, g_0_x_xxy_yyyyyyyz, g_0_x_xxy_yyyyyyz, g_0_x_xxy_yyyyyyzz, g_0_x_xxy_yyyyyzz, g_0_x_xxy_yyyyyzzz, g_0_x_xxy_yyyyzzz, g_0_x_xxy_yyyyzzzz, g_0_x_xxy_yyyzzzz, g_0_x_xxy_yyyzzzzz, g_0_x_xxy_yyzzzzz, g_0_x_xxy_yyzzzzzz, g_0_x_xxy_yzzzzzz, g_0_x_xxy_yzzzzzzz, g_0_x_xxy_zzzzzzz, g_0_x_xxyy_xxxxxxx, g_0_x_xxyy_xxxxxxy, g_0_x_xxyy_xxxxxxz, g_0_x_xxyy_xxxxxyy, g_0_x_xxyy_xxxxxyz, g_0_x_xxyy_xxxxxzz, g_0_x_xxyy_xxxxyyy, g_0_x_xxyy_xxxxyyz, g_0_x_xxyy_xxxxyzz, g_0_x_xxyy_xxxxzzz, g_0_x_xxyy_xxxyyyy, g_0_x_xxyy_xxxyyyz, g_0_x_xxyy_xxxyyzz, g_0_x_xxyy_xxxyzzz, g_0_x_xxyy_xxxzzzz, g_0_x_xxyy_xxyyyyy, g_0_x_xxyy_xxyyyyz, g_0_x_xxyy_xxyyyzz, g_0_x_xxyy_xxyyzzz, g_0_x_xxyy_xxyzzzz, g_0_x_xxyy_xxzzzzz, g_0_x_xxyy_xyyyyyy, g_0_x_xxyy_xyyyyyz, g_0_x_xxyy_xyyyyzz, g_0_x_xxyy_xyyyzzz, g_0_x_xxyy_xyyzzzz, g_0_x_xxyy_xyzzzzz, g_0_x_xxyy_xzzzzzz, g_0_x_xxyy_yyyyyyy, g_0_x_xxyy_yyyyyyz, g_0_x_xxyy_yyyyyzz, g_0_x_xxyy_yyyyzzz, g_0_x_xxyy_yyyzzzz, g_0_x_xxyy_yyzzzzz, g_0_x_xxyy_yzzzzzz, g_0_x_xxyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_xxxxxxx[k] = -g_0_x_xxy_xxxxxxx[k] * ab_y + g_0_x_xxy_xxxxxxxy[k];

                g_0_x_xxyy_xxxxxxy[k] = -g_0_x_xxy_xxxxxxy[k] * ab_y + g_0_x_xxy_xxxxxxyy[k];

                g_0_x_xxyy_xxxxxxz[k] = -g_0_x_xxy_xxxxxxz[k] * ab_y + g_0_x_xxy_xxxxxxyz[k];

                g_0_x_xxyy_xxxxxyy[k] = -g_0_x_xxy_xxxxxyy[k] * ab_y + g_0_x_xxy_xxxxxyyy[k];

                g_0_x_xxyy_xxxxxyz[k] = -g_0_x_xxy_xxxxxyz[k] * ab_y + g_0_x_xxy_xxxxxyyz[k];

                g_0_x_xxyy_xxxxxzz[k] = -g_0_x_xxy_xxxxxzz[k] * ab_y + g_0_x_xxy_xxxxxyzz[k];

                g_0_x_xxyy_xxxxyyy[k] = -g_0_x_xxy_xxxxyyy[k] * ab_y + g_0_x_xxy_xxxxyyyy[k];

                g_0_x_xxyy_xxxxyyz[k] = -g_0_x_xxy_xxxxyyz[k] * ab_y + g_0_x_xxy_xxxxyyyz[k];

                g_0_x_xxyy_xxxxyzz[k] = -g_0_x_xxy_xxxxyzz[k] * ab_y + g_0_x_xxy_xxxxyyzz[k];

                g_0_x_xxyy_xxxxzzz[k] = -g_0_x_xxy_xxxxzzz[k] * ab_y + g_0_x_xxy_xxxxyzzz[k];

                g_0_x_xxyy_xxxyyyy[k] = -g_0_x_xxy_xxxyyyy[k] * ab_y + g_0_x_xxy_xxxyyyyy[k];

                g_0_x_xxyy_xxxyyyz[k] = -g_0_x_xxy_xxxyyyz[k] * ab_y + g_0_x_xxy_xxxyyyyz[k];

                g_0_x_xxyy_xxxyyzz[k] = -g_0_x_xxy_xxxyyzz[k] * ab_y + g_0_x_xxy_xxxyyyzz[k];

                g_0_x_xxyy_xxxyzzz[k] = -g_0_x_xxy_xxxyzzz[k] * ab_y + g_0_x_xxy_xxxyyzzz[k];

                g_0_x_xxyy_xxxzzzz[k] = -g_0_x_xxy_xxxzzzz[k] * ab_y + g_0_x_xxy_xxxyzzzz[k];

                g_0_x_xxyy_xxyyyyy[k] = -g_0_x_xxy_xxyyyyy[k] * ab_y + g_0_x_xxy_xxyyyyyy[k];

                g_0_x_xxyy_xxyyyyz[k] = -g_0_x_xxy_xxyyyyz[k] * ab_y + g_0_x_xxy_xxyyyyyz[k];

                g_0_x_xxyy_xxyyyzz[k] = -g_0_x_xxy_xxyyyzz[k] * ab_y + g_0_x_xxy_xxyyyyzz[k];

                g_0_x_xxyy_xxyyzzz[k] = -g_0_x_xxy_xxyyzzz[k] * ab_y + g_0_x_xxy_xxyyyzzz[k];

                g_0_x_xxyy_xxyzzzz[k] = -g_0_x_xxy_xxyzzzz[k] * ab_y + g_0_x_xxy_xxyyzzzz[k];

                g_0_x_xxyy_xxzzzzz[k] = -g_0_x_xxy_xxzzzzz[k] * ab_y + g_0_x_xxy_xxyzzzzz[k];

                g_0_x_xxyy_xyyyyyy[k] = -g_0_x_xxy_xyyyyyy[k] * ab_y + g_0_x_xxy_xyyyyyyy[k];

                g_0_x_xxyy_xyyyyyz[k] = -g_0_x_xxy_xyyyyyz[k] * ab_y + g_0_x_xxy_xyyyyyyz[k];

                g_0_x_xxyy_xyyyyzz[k] = -g_0_x_xxy_xyyyyzz[k] * ab_y + g_0_x_xxy_xyyyyyzz[k];

                g_0_x_xxyy_xyyyzzz[k] = -g_0_x_xxy_xyyyzzz[k] * ab_y + g_0_x_xxy_xyyyyzzz[k];

                g_0_x_xxyy_xyyzzzz[k] = -g_0_x_xxy_xyyzzzz[k] * ab_y + g_0_x_xxy_xyyyzzzz[k];

                g_0_x_xxyy_xyzzzzz[k] = -g_0_x_xxy_xyzzzzz[k] * ab_y + g_0_x_xxy_xyyzzzzz[k];

                g_0_x_xxyy_xzzzzzz[k] = -g_0_x_xxy_xzzzzzz[k] * ab_y + g_0_x_xxy_xyzzzzzz[k];

                g_0_x_xxyy_yyyyyyy[k] = -g_0_x_xxy_yyyyyyy[k] * ab_y + g_0_x_xxy_yyyyyyyy[k];

                g_0_x_xxyy_yyyyyyz[k] = -g_0_x_xxy_yyyyyyz[k] * ab_y + g_0_x_xxy_yyyyyyyz[k];

                g_0_x_xxyy_yyyyyzz[k] = -g_0_x_xxy_yyyyyzz[k] * ab_y + g_0_x_xxy_yyyyyyzz[k];

                g_0_x_xxyy_yyyyzzz[k] = -g_0_x_xxy_yyyyzzz[k] * ab_y + g_0_x_xxy_yyyyyzzz[k];

                g_0_x_xxyy_yyyzzzz[k] = -g_0_x_xxy_yyyzzzz[k] * ab_y + g_0_x_xxy_yyyyzzzz[k];

                g_0_x_xxyy_yyzzzzz[k] = -g_0_x_xxy_yyzzzzz[k] * ab_y + g_0_x_xxy_yyyzzzzz[k];

                g_0_x_xxyy_yzzzzzz[k] = -g_0_x_xxy_yzzzzzz[k] * ab_y + g_0_x_xxy_yyzzzzzz[k];

                g_0_x_xxyy_zzzzzzz[k] = -g_0_x_xxy_zzzzzzz[k] * ab_y + g_0_x_xxy_yzzzzzzz[k];
            }

            /// Set up 144-180 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxyz_xxxxxxx, g_0_x_xxyz_xxxxxxy, g_0_x_xxyz_xxxxxxz, g_0_x_xxyz_xxxxxyy, g_0_x_xxyz_xxxxxyz, g_0_x_xxyz_xxxxxzz, g_0_x_xxyz_xxxxyyy, g_0_x_xxyz_xxxxyyz, g_0_x_xxyz_xxxxyzz, g_0_x_xxyz_xxxxzzz, g_0_x_xxyz_xxxyyyy, g_0_x_xxyz_xxxyyyz, g_0_x_xxyz_xxxyyzz, g_0_x_xxyz_xxxyzzz, g_0_x_xxyz_xxxzzzz, g_0_x_xxyz_xxyyyyy, g_0_x_xxyz_xxyyyyz, g_0_x_xxyz_xxyyyzz, g_0_x_xxyz_xxyyzzz, g_0_x_xxyz_xxyzzzz, g_0_x_xxyz_xxzzzzz, g_0_x_xxyz_xyyyyyy, g_0_x_xxyz_xyyyyyz, g_0_x_xxyz_xyyyyzz, g_0_x_xxyz_xyyyzzz, g_0_x_xxyz_xyyzzzz, g_0_x_xxyz_xyzzzzz, g_0_x_xxyz_xzzzzzz, g_0_x_xxyz_yyyyyyy, g_0_x_xxyz_yyyyyyz, g_0_x_xxyz_yyyyyzz, g_0_x_xxyz_yyyyzzz, g_0_x_xxyz_yyyzzzz, g_0_x_xxyz_yyzzzzz, g_0_x_xxyz_yzzzzzz, g_0_x_xxyz_zzzzzzz, g_0_x_xxz_xxxxxxx, g_0_x_xxz_xxxxxxxy, g_0_x_xxz_xxxxxxy, g_0_x_xxz_xxxxxxyy, g_0_x_xxz_xxxxxxyz, g_0_x_xxz_xxxxxxz, g_0_x_xxz_xxxxxyy, g_0_x_xxz_xxxxxyyy, g_0_x_xxz_xxxxxyyz, g_0_x_xxz_xxxxxyz, g_0_x_xxz_xxxxxyzz, g_0_x_xxz_xxxxxzz, g_0_x_xxz_xxxxyyy, g_0_x_xxz_xxxxyyyy, g_0_x_xxz_xxxxyyyz, g_0_x_xxz_xxxxyyz, g_0_x_xxz_xxxxyyzz, g_0_x_xxz_xxxxyzz, g_0_x_xxz_xxxxyzzz, g_0_x_xxz_xxxxzzz, g_0_x_xxz_xxxyyyy, g_0_x_xxz_xxxyyyyy, g_0_x_xxz_xxxyyyyz, g_0_x_xxz_xxxyyyz, g_0_x_xxz_xxxyyyzz, g_0_x_xxz_xxxyyzz, g_0_x_xxz_xxxyyzzz, g_0_x_xxz_xxxyzzz, g_0_x_xxz_xxxyzzzz, g_0_x_xxz_xxxzzzz, g_0_x_xxz_xxyyyyy, g_0_x_xxz_xxyyyyyy, g_0_x_xxz_xxyyyyyz, g_0_x_xxz_xxyyyyz, g_0_x_xxz_xxyyyyzz, g_0_x_xxz_xxyyyzz, g_0_x_xxz_xxyyyzzz, g_0_x_xxz_xxyyzzz, g_0_x_xxz_xxyyzzzz, g_0_x_xxz_xxyzzzz, g_0_x_xxz_xxyzzzzz, g_0_x_xxz_xxzzzzz, g_0_x_xxz_xyyyyyy, g_0_x_xxz_xyyyyyyy, g_0_x_xxz_xyyyyyyz, g_0_x_xxz_xyyyyyz, g_0_x_xxz_xyyyyyzz, g_0_x_xxz_xyyyyzz, g_0_x_xxz_xyyyyzzz, g_0_x_xxz_xyyyzzz, g_0_x_xxz_xyyyzzzz, g_0_x_xxz_xyyzzzz, g_0_x_xxz_xyyzzzzz, g_0_x_xxz_xyzzzzz, g_0_x_xxz_xyzzzzzz, g_0_x_xxz_xzzzzzz, g_0_x_xxz_yyyyyyy, g_0_x_xxz_yyyyyyyy, g_0_x_xxz_yyyyyyyz, g_0_x_xxz_yyyyyyz, g_0_x_xxz_yyyyyyzz, g_0_x_xxz_yyyyyzz, g_0_x_xxz_yyyyyzzz, g_0_x_xxz_yyyyzzz, g_0_x_xxz_yyyyzzzz, g_0_x_xxz_yyyzzzz, g_0_x_xxz_yyyzzzzz, g_0_x_xxz_yyzzzzz, g_0_x_xxz_yyzzzzzz, g_0_x_xxz_yzzzzzz, g_0_x_xxz_yzzzzzzz, g_0_x_xxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_xxxxxxx[k] = -g_0_x_xxz_xxxxxxx[k] * ab_y + g_0_x_xxz_xxxxxxxy[k];

                g_0_x_xxyz_xxxxxxy[k] = -g_0_x_xxz_xxxxxxy[k] * ab_y + g_0_x_xxz_xxxxxxyy[k];

                g_0_x_xxyz_xxxxxxz[k] = -g_0_x_xxz_xxxxxxz[k] * ab_y + g_0_x_xxz_xxxxxxyz[k];

                g_0_x_xxyz_xxxxxyy[k] = -g_0_x_xxz_xxxxxyy[k] * ab_y + g_0_x_xxz_xxxxxyyy[k];

                g_0_x_xxyz_xxxxxyz[k] = -g_0_x_xxz_xxxxxyz[k] * ab_y + g_0_x_xxz_xxxxxyyz[k];

                g_0_x_xxyz_xxxxxzz[k] = -g_0_x_xxz_xxxxxzz[k] * ab_y + g_0_x_xxz_xxxxxyzz[k];

                g_0_x_xxyz_xxxxyyy[k] = -g_0_x_xxz_xxxxyyy[k] * ab_y + g_0_x_xxz_xxxxyyyy[k];

                g_0_x_xxyz_xxxxyyz[k] = -g_0_x_xxz_xxxxyyz[k] * ab_y + g_0_x_xxz_xxxxyyyz[k];

                g_0_x_xxyz_xxxxyzz[k] = -g_0_x_xxz_xxxxyzz[k] * ab_y + g_0_x_xxz_xxxxyyzz[k];

                g_0_x_xxyz_xxxxzzz[k] = -g_0_x_xxz_xxxxzzz[k] * ab_y + g_0_x_xxz_xxxxyzzz[k];

                g_0_x_xxyz_xxxyyyy[k] = -g_0_x_xxz_xxxyyyy[k] * ab_y + g_0_x_xxz_xxxyyyyy[k];

                g_0_x_xxyz_xxxyyyz[k] = -g_0_x_xxz_xxxyyyz[k] * ab_y + g_0_x_xxz_xxxyyyyz[k];

                g_0_x_xxyz_xxxyyzz[k] = -g_0_x_xxz_xxxyyzz[k] * ab_y + g_0_x_xxz_xxxyyyzz[k];

                g_0_x_xxyz_xxxyzzz[k] = -g_0_x_xxz_xxxyzzz[k] * ab_y + g_0_x_xxz_xxxyyzzz[k];

                g_0_x_xxyz_xxxzzzz[k] = -g_0_x_xxz_xxxzzzz[k] * ab_y + g_0_x_xxz_xxxyzzzz[k];

                g_0_x_xxyz_xxyyyyy[k] = -g_0_x_xxz_xxyyyyy[k] * ab_y + g_0_x_xxz_xxyyyyyy[k];

                g_0_x_xxyz_xxyyyyz[k] = -g_0_x_xxz_xxyyyyz[k] * ab_y + g_0_x_xxz_xxyyyyyz[k];

                g_0_x_xxyz_xxyyyzz[k] = -g_0_x_xxz_xxyyyzz[k] * ab_y + g_0_x_xxz_xxyyyyzz[k];

                g_0_x_xxyz_xxyyzzz[k] = -g_0_x_xxz_xxyyzzz[k] * ab_y + g_0_x_xxz_xxyyyzzz[k];

                g_0_x_xxyz_xxyzzzz[k] = -g_0_x_xxz_xxyzzzz[k] * ab_y + g_0_x_xxz_xxyyzzzz[k];

                g_0_x_xxyz_xxzzzzz[k] = -g_0_x_xxz_xxzzzzz[k] * ab_y + g_0_x_xxz_xxyzzzzz[k];

                g_0_x_xxyz_xyyyyyy[k] = -g_0_x_xxz_xyyyyyy[k] * ab_y + g_0_x_xxz_xyyyyyyy[k];

                g_0_x_xxyz_xyyyyyz[k] = -g_0_x_xxz_xyyyyyz[k] * ab_y + g_0_x_xxz_xyyyyyyz[k];

                g_0_x_xxyz_xyyyyzz[k] = -g_0_x_xxz_xyyyyzz[k] * ab_y + g_0_x_xxz_xyyyyyzz[k];

                g_0_x_xxyz_xyyyzzz[k] = -g_0_x_xxz_xyyyzzz[k] * ab_y + g_0_x_xxz_xyyyyzzz[k];

                g_0_x_xxyz_xyyzzzz[k] = -g_0_x_xxz_xyyzzzz[k] * ab_y + g_0_x_xxz_xyyyzzzz[k];

                g_0_x_xxyz_xyzzzzz[k] = -g_0_x_xxz_xyzzzzz[k] * ab_y + g_0_x_xxz_xyyzzzzz[k];

                g_0_x_xxyz_xzzzzzz[k] = -g_0_x_xxz_xzzzzzz[k] * ab_y + g_0_x_xxz_xyzzzzzz[k];

                g_0_x_xxyz_yyyyyyy[k] = -g_0_x_xxz_yyyyyyy[k] * ab_y + g_0_x_xxz_yyyyyyyy[k];

                g_0_x_xxyz_yyyyyyz[k] = -g_0_x_xxz_yyyyyyz[k] * ab_y + g_0_x_xxz_yyyyyyyz[k];

                g_0_x_xxyz_yyyyyzz[k] = -g_0_x_xxz_yyyyyzz[k] * ab_y + g_0_x_xxz_yyyyyyzz[k];

                g_0_x_xxyz_yyyyzzz[k] = -g_0_x_xxz_yyyyzzz[k] * ab_y + g_0_x_xxz_yyyyyzzz[k];

                g_0_x_xxyz_yyyzzzz[k] = -g_0_x_xxz_yyyzzzz[k] * ab_y + g_0_x_xxz_yyyyzzzz[k];

                g_0_x_xxyz_yyzzzzz[k] = -g_0_x_xxz_yyzzzzz[k] * ab_y + g_0_x_xxz_yyyzzzzz[k];

                g_0_x_xxyz_yzzzzzz[k] = -g_0_x_xxz_yzzzzzz[k] * ab_y + g_0_x_xxz_yyzzzzzz[k];

                g_0_x_xxyz_zzzzzzz[k] = -g_0_x_xxz_zzzzzzz[k] * ab_y + g_0_x_xxz_yzzzzzzz[k];
            }

            /// Set up 180-216 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxz_xxxxxxx, g_0_x_xxz_xxxxxxxz, g_0_x_xxz_xxxxxxy, g_0_x_xxz_xxxxxxyz, g_0_x_xxz_xxxxxxz, g_0_x_xxz_xxxxxxzz, g_0_x_xxz_xxxxxyy, g_0_x_xxz_xxxxxyyz, g_0_x_xxz_xxxxxyz, g_0_x_xxz_xxxxxyzz, g_0_x_xxz_xxxxxzz, g_0_x_xxz_xxxxxzzz, g_0_x_xxz_xxxxyyy, g_0_x_xxz_xxxxyyyz, g_0_x_xxz_xxxxyyz, g_0_x_xxz_xxxxyyzz, g_0_x_xxz_xxxxyzz, g_0_x_xxz_xxxxyzzz, g_0_x_xxz_xxxxzzz, g_0_x_xxz_xxxxzzzz, g_0_x_xxz_xxxyyyy, g_0_x_xxz_xxxyyyyz, g_0_x_xxz_xxxyyyz, g_0_x_xxz_xxxyyyzz, g_0_x_xxz_xxxyyzz, g_0_x_xxz_xxxyyzzz, g_0_x_xxz_xxxyzzz, g_0_x_xxz_xxxyzzzz, g_0_x_xxz_xxxzzzz, g_0_x_xxz_xxxzzzzz, g_0_x_xxz_xxyyyyy, g_0_x_xxz_xxyyyyyz, g_0_x_xxz_xxyyyyz, g_0_x_xxz_xxyyyyzz, g_0_x_xxz_xxyyyzz, g_0_x_xxz_xxyyyzzz, g_0_x_xxz_xxyyzzz, g_0_x_xxz_xxyyzzzz, g_0_x_xxz_xxyzzzz, g_0_x_xxz_xxyzzzzz, g_0_x_xxz_xxzzzzz, g_0_x_xxz_xxzzzzzz, g_0_x_xxz_xyyyyyy, g_0_x_xxz_xyyyyyyz, g_0_x_xxz_xyyyyyz, g_0_x_xxz_xyyyyyzz, g_0_x_xxz_xyyyyzz, g_0_x_xxz_xyyyyzzz, g_0_x_xxz_xyyyzzz, g_0_x_xxz_xyyyzzzz, g_0_x_xxz_xyyzzzz, g_0_x_xxz_xyyzzzzz, g_0_x_xxz_xyzzzzz, g_0_x_xxz_xyzzzzzz, g_0_x_xxz_xzzzzzz, g_0_x_xxz_xzzzzzzz, g_0_x_xxz_yyyyyyy, g_0_x_xxz_yyyyyyyz, g_0_x_xxz_yyyyyyz, g_0_x_xxz_yyyyyyzz, g_0_x_xxz_yyyyyzz, g_0_x_xxz_yyyyyzzz, g_0_x_xxz_yyyyzzz, g_0_x_xxz_yyyyzzzz, g_0_x_xxz_yyyzzzz, g_0_x_xxz_yyyzzzzz, g_0_x_xxz_yyzzzzz, g_0_x_xxz_yyzzzzzz, g_0_x_xxz_yzzzzzz, g_0_x_xxz_yzzzzzzz, g_0_x_xxz_zzzzzzz, g_0_x_xxz_zzzzzzzz, g_0_x_xxzz_xxxxxxx, g_0_x_xxzz_xxxxxxy, g_0_x_xxzz_xxxxxxz, g_0_x_xxzz_xxxxxyy, g_0_x_xxzz_xxxxxyz, g_0_x_xxzz_xxxxxzz, g_0_x_xxzz_xxxxyyy, g_0_x_xxzz_xxxxyyz, g_0_x_xxzz_xxxxyzz, g_0_x_xxzz_xxxxzzz, g_0_x_xxzz_xxxyyyy, g_0_x_xxzz_xxxyyyz, g_0_x_xxzz_xxxyyzz, g_0_x_xxzz_xxxyzzz, g_0_x_xxzz_xxxzzzz, g_0_x_xxzz_xxyyyyy, g_0_x_xxzz_xxyyyyz, g_0_x_xxzz_xxyyyzz, g_0_x_xxzz_xxyyzzz, g_0_x_xxzz_xxyzzzz, g_0_x_xxzz_xxzzzzz, g_0_x_xxzz_xyyyyyy, g_0_x_xxzz_xyyyyyz, g_0_x_xxzz_xyyyyzz, g_0_x_xxzz_xyyyzzz, g_0_x_xxzz_xyyzzzz, g_0_x_xxzz_xyzzzzz, g_0_x_xxzz_xzzzzzz, g_0_x_xxzz_yyyyyyy, g_0_x_xxzz_yyyyyyz, g_0_x_xxzz_yyyyyzz, g_0_x_xxzz_yyyyzzz, g_0_x_xxzz_yyyzzzz, g_0_x_xxzz_yyzzzzz, g_0_x_xxzz_yzzzzzz, g_0_x_xxzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_xxxxxxx[k] = -g_0_x_xxz_xxxxxxx[k] * ab_z + g_0_x_xxz_xxxxxxxz[k];

                g_0_x_xxzz_xxxxxxy[k] = -g_0_x_xxz_xxxxxxy[k] * ab_z + g_0_x_xxz_xxxxxxyz[k];

                g_0_x_xxzz_xxxxxxz[k] = -g_0_x_xxz_xxxxxxz[k] * ab_z + g_0_x_xxz_xxxxxxzz[k];

                g_0_x_xxzz_xxxxxyy[k] = -g_0_x_xxz_xxxxxyy[k] * ab_z + g_0_x_xxz_xxxxxyyz[k];

                g_0_x_xxzz_xxxxxyz[k] = -g_0_x_xxz_xxxxxyz[k] * ab_z + g_0_x_xxz_xxxxxyzz[k];

                g_0_x_xxzz_xxxxxzz[k] = -g_0_x_xxz_xxxxxzz[k] * ab_z + g_0_x_xxz_xxxxxzzz[k];

                g_0_x_xxzz_xxxxyyy[k] = -g_0_x_xxz_xxxxyyy[k] * ab_z + g_0_x_xxz_xxxxyyyz[k];

                g_0_x_xxzz_xxxxyyz[k] = -g_0_x_xxz_xxxxyyz[k] * ab_z + g_0_x_xxz_xxxxyyzz[k];

                g_0_x_xxzz_xxxxyzz[k] = -g_0_x_xxz_xxxxyzz[k] * ab_z + g_0_x_xxz_xxxxyzzz[k];

                g_0_x_xxzz_xxxxzzz[k] = -g_0_x_xxz_xxxxzzz[k] * ab_z + g_0_x_xxz_xxxxzzzz[k];

                g_0_x_xxzz_xxxyyyy[k] = -g_0_x_xxz_xxxyyyy[k] * ab_z + g_0_x_xxz_xxxyyyyz[k];

                g_0_x_xxzz_xxxyyyz[k] = -g_0_x_xxz_xxxyyyz[k] * ab_z + g_0_x_xxz_xxxyyyzz[k];

                g_0_x_xxzz_xxxyyzz[k] = -g_0_x_xxz_xxxyyzz[k] * ab_z + g_0_x_xxz_xxxyyzzz[k];

                g_0_x_xxzz_xxxyzzz[k] = -g_0_x_xxz_xxxyzzz[k] * ab_z + g_0_x_xxz_xxxyzzzz[k];

                g_0_x_xxzz_xxxzzzz[k] = -g_0_x_xxz_xxxzzzz[k] * ab_z + g_0_x_xxz_xxxzzzzz[k];

                g_0_x_xxzz_xxyyyyy[k] = -g_0_x_xxz_xxyyyyy[k] * ab_z + g_0_x_xxz_xxyyyyyz[k];

                g_0_x_xxzz_xxyyyyz[k] = -g_0_x_xxz_xxyyyyz[k] * ab_z + g_0_x_xxz_xxyyyyzz[k];

                g_0_x_xxzz_xxyyyzz[k] = -g_0_x_xxz_xxyyyzz[k] * ab_z + g_0_x_xxz_xxyyyzzz[k];

                g_0_x_xxzz_xxyyzzz[k] = -g_0_x_xxz_xxyyzzz[k] * ab_z + g_0_x_xxz_xxyyzzzz[k];

                g_0_x_xxzz_xxyzzzz[k] = -g_0_x_xxz_xxyzzzz[k] * ab_z + g_0_x_xxz_xxyzzzzz[k];

                g_0_x_xxzz_xxzzzzz[k] = -g_0_x_xxz_xxzzzzz[k] * ab_z + g_0_x_xxz_xxzzzzzz[k];

                g_0_x_xxzz_xyyyyyy[k] = -g_0_x_xxz_xyyyyyy[k] * ab_z + g_0_x_xxz_xyyyyyyz[k];

                g_0_x_xxzz_xyyyyyz[k] = -g_0_x_xxz_xyyyyyz[k] * ab_z + g_0_x_xxz_xyyyyyzz[k];

                g_0_x_xxzz_xyyyyzz[k] = -g_0_x_xxz_xyyyyzz[k] * ab_z + g_0_x_xxz_xyyyyzzz[k];

                g_0_x_xxzz_xyyyzzz[k] = -g_0_x_xxz_xyyyzzz[k] * ab_z + g_0_x_xxz_xyyyzzzz[k];

                g_0_x_xxzz_xyyzzzz[k] = -g_0_x_xxz_xyyzzzz[k] * ab_z + g_0_x_xxz_xyyzzzzz[k];

                g_0_x_xxzz_xyzzzzz[k] = -g_0_x_xxz_xyzzzzz[k] * ab_z + g_0_x_xxz_xyzzzzzz[k];

                g_0_x_xxzz_xzzzzzz[k] = -g_0_x_xxz_xzzzzzz[k] * ab_z + g_0_x_xxz_xzzzzzzz[k];

                g_0_x_xxzz_yyyyyyy[k] = -g_0_x_xxz_yyyyyyy[k] * ab_z + g_0_x_xxz_yyyyyyyz[k];

                g_0_x_xxzz_yyyyyyz[k] = -g_0_x_xxz_yyyyyyz[k] * ab_z + g_0_x_xxz_yyyyyyzz[k];

                g_0_x_xxzz_yyyyyzz[k] = -g_0_x_xxz_yyyyyzz[k] * ab_z + g_0_x_xxz_yyyyyzzz[k];

                g_0_x_xxzz_yyyyzzz[k] = -g_0_x_xxz_yyyyzzz[k] * ab_z + g_0_x_xxz_yyyyzzzz[k];

                g_0_x_xxzz_yyyzzzz[k] = -g_0_x_xxz_yyyzzzz[k] * ab_z + g_0_x_xxz_yyyzzzzz[k];

                g_0_x_xxzz_yyzzzzz[k] = -g_0_x_xxz_yyzzzzz[k] * ab_z + g_0_x_xxz_yyzzzzzz[k];

                g_0_x_xxzz_yzzzzzz[k] = -g_0_x_xxz_yzzzzzz[k] * ab_z + g_0_x_xxz_yzzzzzzz[k];

                g_0_x_xxzz_zzzzzzz[k] = -g_0_x_xxz_zzzzzzz[k] * ab_z + g_0_x_xxz_zzzzzzzz[k];
            }

            /// Set up 216-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyy_xxxxxxx, g_0_x_xyy_xxxxxxxy, g_0_x_xyy_xxxxxxy, g_0_x_xyy_xxxxxxyy, g_0_x_xyy_xxxxxxyz, g_0_x_xyy_xxxxxxz, g_0_x_xyy_xxxxxyy, g_0_x_xyy_xxxxxyyy, g_0_x_xyy_xxxxxyyz, g_0_x_xyy_xxxxxyz, g_0_x_xyy_xxxxxyzz, g_0_x_xyy_xxxxxzz, g_0_x_xyy_xxxxyyy, g_0_x_xyy_xxxxyyyy, g_0_x_xyy_xxxxyyyz, g_0_x_xyy_xxxxyyz, g_0_x_xyy_xxxxyyzz, g_0_x_xyy_xxxxyzz, g_0_x_xyy_xxxxyzzz, g_0_x_xyy_xxxxzzz, g_0_x_xyy_xxxyyyy, g_0_x_xyy_xxxyyyyy, g_0_x_xyy_xxxyyyyz, g_0_x_xyy_xxxyyyz, g_0_x_xyy_xxxyyyzz, g_0_x_xyy_xxxyyzz, g_0_x_xyy_xxxyyzzz, g_0_x_xyy_xxxyzzz, g_0_x_xyy_xxxyzzzz, g_0_x_xyy_xxxzzzz, g_0_x_xyy_xxyyyyy, g_0_x_xyy_xxyyyyyy, g_0_x_xyy_xxyyyyyz, g_0_x_xyy_xxyyyyz, g_0_x_xyy_xxyyyyzz, g_0_x_xyy_xxyyyzz, g_0_x_xyy_xxyyyzzz, g_0_x_xyy_xxyyzzz, g_0_x_xyy_xxyyzzzz, g_0_x_xyy_xxyzzzz, g_0_x_xyy_xxyzzzzz, g_0_x_xyy_xxzzzzz, g_0_x_xyy_xyyyyyy, g_0_x_xyy_xyyyyyyy, g_0_x_xyy_xyyyyyyz, g_0_x_xyy_xyyyyyz, g_0_x_xyy_xyyyyyzz, g_0_x_xyy_xyyyyzz, g_0_x_xyy_xyyyyzzz, g_0_x_xyy_xyyyzzz, g_0_x_xyy_xyyyzzzz, g_0_x_xyy_xyyzzzz, g_0_x_xyy_xyyzzzzz, g_0_x_xyy_xyzzzzz, g_0_x_xyy_xyzzzzzz, g_0_x_xyy_xzzzzzz, g_0_x_xyy_yyyyyyy, g_0_x_xyy_yyyyyyyy, g_0_x_xyy_yyyyyyyz, g_0_x_xyy_yyyyyyz, g_0_x_xyy_yyyyyyzz, g_0_x_xyy_yyyyyzz, g_0_x_xyy_yyyyyzzz, g_0_x_xyy_yyyyzzz, g_0_x_xyy_yyyyzzzz, g_0_x_xyy_yyyzzzz, g_0_x_xyy_yyyzzzzz, g_0_x_xyy_yyzzzzz, g_0_x_xyy_yyzzzzzz, g_0_x_xyy_yzzzzzz, g_0_x_xyy_yzzzzzzz, g_0_x_xyy_zzzzzzz, g_0_x_xyyy_xxxxxxx, g_0_x_xyyy_xxxxxxy, g_0_x_xyyy_xxxxxxz, g_0_x_xyyy_xxxxxyy, g_0_x_xyyy_xxxxxyz, g_0_x_xyyy_xxxxxzz, g_0_x_xyyy_xxxxyyy, g_0_x_xyyy_xxxxyyz, g_0_x_xyyy_xxxxyzz, g_0_x_xyyy_xxxxzzz, g_0_x_xyyy_xxxyyyy, g_0_x_xyyy_xxxyyyz, g_0_x_xyyy_xxxyyzz, g_0_x_xyyy_xxxyzzz, g_0_x_xyyy_xxxzzzz, g_0_x_xyyy_xxyyyyy, g_0_x_xyyy_xxyyyyz, g_0_x_xyyy_xxyyyzz, g_0_x_xyyy_xxyyzzz, g_0_x_xyyy_xxyzzzz, g_0_x_xyyy_xxzzzzz, g_0_x_xyyy_xyyyyyy, g_0_x_xyyy_xyyyyyz, g_0_x_xyyy_xyyyyzz, g_0_x_xyyy_xyyyzzz, g_0_x_xyyy_xyyzzzz, g_0_x_xyyy_xyzzzzz, g_0_x_xyyy_xzzzzzz, g_0_x_xyyy_yyyyyyy, g_0_x_xyyy_yyyyyyz, g_0_x_xyyy_yyyyyzz, g_0_x_xyyy_yyyyzzz, g_0_x_xyyy_yyyzzzz, g_0_x_xyyy_yyzzzzz, g_0_x_xyyy_yzzzzzz, g_0_x_xyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_xxxxxxx[k] = -g_0_x_xyy_xxxxxxx[k] * ab_y + g_0_x_xyy_xxxxxxxy[k];

                g_0_x_xyyy_xxxxxxy[k] = -g_0_x_xyy_xxxxxxy[k] * ab_y + g_0_x_xyy_xxxxxxyy[k];

                g_0_x_xyyy_xxxxxxz[k] = -g_0_x_xyy_xxxxxxz[k] * ab_y + g_0_x_xyy_xxxxxxyz[k];

                g_0_x_xyyy_xxxxxyy[k] = -g_0_x_xyy_xxxxxyy[k] * ab_y + g_0_x_xyy_xxxxxyyy[k];

                g_0_x_xyyy_xxxxxyz[k] = -g_0_x_xyy_xxxxxyz[k] * ab_y + g_0_x_xyy_xxxxxyyz[k];

                g_0_x_xyyy_xxxxxzz[k] = -g_0_x_xyy_xxxxxzz[k] * ab_y + g_0_x_xyy_xxxxxyzz[k];

                g_0_x_xyyy_xxxxyyy[k] = -g_0_x_xyy_xxxxyyy[k] * ab_y + g_0_x_xyy_xxxxyyyy[k];

                g_0_x_xyyy_xxxxyyz[k] = -g_0_x_xyy_xxxxyyz[k] * ab_y + g_0_x_xyy_xxxxyyyz[k];

                g_0_x_xyyy_xxxxyzz[k] = -g_0_x_xyy_xxxxyzz[k] * ab_y + g_0_x_xyy_xxxxyyzz[k];

                g_0_x_xyyy_xxxxzzz[k] = -g_0_x_xyy_xxxxzzz[k] * ab_y + g_0_x_xyy_xxxxyzzz[k];

                g_0_x_xyyy_xxxyyyy[k] = -g_0_x_xyy_xxxyyyy[k] * ab_y + g_0_x_xyy_xxxyyyyy[k];

                g_0_x_xyyy_xxxyyyz[k] = -g_0_x_xyy_xxxyyyz[k] * ab_y + g_0_x_xyy_xxxyyyyz[k];

                g_0_x_xyyy_xxxyyzz[k] = -g_0_x_xyy_xxxyyzz[k] * ab_y + g_0_x_xyy_xxxyyyzz[k];

                g_0_x_xyyy_xxxyzzz[k] = -g_0_x_xyy_xxxyzzz[k] * ab_y + g_0_x_xyy_xxxyyzzz[k];

                g_0_x_xyyy_xxxzzzz[k] = -g_0_x_xyy_xxxzzzz[k] * ab_y + g_0_x_xyy_xxxyzzzz[k];

                g_0_x_xyyy_xxyyyyy[k] = -g_0_x_xyy_xxyyyyy[k] * ab_y + g_0_x_xyy_xxyyyyyy[k];

                g_0_x_xyyy_xxyyyyz[k] = -g_0_x_xyy_xxyyyyz[k] * ab_y + g_0_x_xyy_xxyyyyyz[k];

                g_0_x_xyyy_xxyyyzz[k] = -g_0_x_xyy_xxyyyzz[k] * ab_y + g_0_x_xyy_xxyyyyzz[k];

                g_0_x_xyyy_xxyyzzz[k] = -g_0_x_xyy_xxyyzzz[k] * ab_y + g_0_x_xyy_xxyyyzzz[k];

                g_0_x_xyyy_xxyzzzz[k] = -g_0_x_xyy_xxyzzzz[k] * ab_y + g_0_x_xyy_xxyyzzzz[k];

                g_0_x_xyyy_xxzzzzz[k] = -g_0_x_xyy_xxzzzzz[k] * ab_y + g_0_x_xyy_xxyzzzzz[k];

                g_0_x_xyyy_xyyyyyy[k] = -g_0_x_xyy_xyyyyyy[k] * ab_y + g_0_x_xyy_xyyyyyyy[k];

                g_0_x_xyyy_xyyyyyz[k] = -g_0_x_xyy_xyyyyyz[k] * ab_y + g_0_x_xyy_xyyyyyyz[k];

                g_0_x_xyyy_xyyyyzz[k] = -g_0_x_xyy_xyyyyzz[k] * ab_y + g_0_x_xyy_xyyyyyzz[k];

                g_0_x_xyyy_xyyyzzz[k] = -g_0_x_xyy_xyyyzzz[k] * ab_y + g_0_x_xyy_xyyyyzzz[k];

                g_0_x_xyyy_xyyzzzz[k] = -g_0_x_xyy_xyyzzzz[k] * ab_y + g_0_x_xyy_xyyyzzzz[k];

                g_0_x_xyyy_xyzzzzz[k] = -g_0_x_xyy_xyzzzzz[k] * ab_y + g_0_x_xyy_xyyzzzzz[k];

                g_0_x_xyyy_xzzzzzz[k] = -g_0_x_xyy_xzzzzzz[k] * ab_y + g_0_x_xyy_xyzzzzzz[k];

                g_0_x_xyyy_yyyyyyy[k] = -g_0_x_xyy_yyyyyyy[k] * ab_y + g_0_x_xyy_yyyyyyyy[k];

                g_0_x_xyyy_yyyyyyz[k] = -g_0_x_xyy_yyyyyyz[k] * ab_y + g_0_x_xyy_yyyyyyyz[k];

                g_0_x_xyyy_yyyyyzz[k] = -g_0_x_xyy_yyyyyzz[k] * ab_y + g_0_x_xyy_yyyyyyzz[k];

                g_0_x_xyyy_yyyyzzz[k] = -g_0_x_xyy_yyyyzzz[k] * ab_y + g_0_x_xyy_yyyyyzzz[k];

                g_0_x_xyyy_yyyzzzz[k] = -g_0_x_xyy_yyyzzzz[k] * ab_y + g_0_x_xyy_yyyyzzzz[k];

                g_0_x_xyyy_yyzzzzz[k] = -g_0_x_xyy_yyzzzzz[k] * ab_y + g_0_x_xyy_yyyzzzzz[k];

                g_0_x_xyyy_yzzzzzz[k] = -g_0_x_xyy_yzzzzzz[k] * ab_y + g_0_x_xyy_yyzzzzzz[k];

                g_0_x_xyyy_zzzzzzz[k] = -g_0_x_xyy_zzzzzzz[k] * ab_y + g_0_x_xyy_yzzzzzzz[k];
            }

            /// Set up 252-288 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyyz_xxxxxxx, g_0_x_xyyz_xxxxxxy, g_0_x_xyyz_xxxxxxz, g_0_x_xyyz_xxxxxyy, g_0_x_xyyz_xxxxxyz, g_0_x_xyyz_xxxxxzz, g_0_x_xyyz_xxxxyyy, g_0_x_xyyz_xxxxyyz, g_0_x_xyyz_xxxxyzz, g_0_x_xyyz_xxxxzzz, g_0_x_xyyz_xxxyyyy, g_0_x_xyyz_xxxyyyz, g_0_x_xyyz_xxxyyzz, g_0_x_xyyz_xxxyzzz, g_0_x_xyyz_xxxzzzz, g_0_x_xyyz_xxyyyyy, g_0_x_xyyz_xxyyyyz, g_0_x_xyyz_xxyyyzz, g_0_x_xyyz_xxyyzzz, g_0_x_xyyz_xxyzzzz, g_0_x_xyyz_xxzzzzz, g_0_x_xyyz_xyyyyyy, g_0_x_xyyz_xyyyyyz, g_0_x_xyyz_xyyyyzz, g_0_x_xyyz_xyyyzzz, g_0_x_xyyz_xyyzzzz, g_0_x_xyyz_xyzzzzz, g_0_x_xyyz_xzzzzzz, g_0_x_xyyz_yyyyyyy, g_0_x_xyyz_yyyyyyz, g_0_x_xyyz_yyyyyzz, g_0_x_xyyz_yyyyzzz, g_0_x_xyyz_yyyzzzz, g_0_x_xyyz_yyzzzzz, g_0_x_xyyz_yzzzzzz, g_0_x_xyyz_zzzzzzz, g_0_x_xyz_xxxxxxx, g_0_x_xyz_xxxxxxxy, g_0_x_xyz_xxxxxxy, g_0_x_xyz_xxxxxxyy, g_0_x_xyz_xxxxxxyz, g_0_x_xyz_xxxxxxz, g_0_x_xyz_xxxxxyy, g_0_x_xyz_xxxxxyyy, g_0_x_xyz_xxxxxyyz, g_0_x_xyz_xxxxxyz, g_0_x_xyz_xxxxxyzz, g_0_x_xyz_xxxxxzz, g_0_x_xyz_xxxxyyy, g_0_x_xyz_xxxxyyyy, g_0_x_xyz_xxxxyyyz, g_0_x_xyz_xxxxyyz, g_0_x_xyz_xxxxyyzz, g_0_x_xyz_xxxxyzz, g_0_x_xyz_xxxxyzzz, g_0_x_xyz_xxxxzzz, g_0_x_xyz_xxxyyyy, g_0_x_xyz_xxxyyyyy, g_0_x_xyz_xxxyyyyz, g_0_x_xyz_xxxyyyz, g_0_x_xyz_xxxyyyzz, g_0_x_xyz_xxxyyzz, g_0_x_xyz_xxxyyzzz, g_0_x_xyz_xxxyzzz, g_0_x_xyz_xxxyzzzz, g_0_x_xyz_xxxzzzz, g_0_x_xyz_xxyyyyy, g_0_x_xyz_xxyyyyyy, g_0_x_xyz_xxyyyyyz, g_0_x_xyz_xxyyyyz, g_0_x_xyz_xxyyyyzz, g_0_x_xyz_xxyyyzz, g_0_x_xyz_xxyyyzzz, g_0_x_xyz_xxyyzzz, g_0_x_xyz_xxyyzzzz, g_0_x_xyz_xxyzzzz, g_0_x_xyz_xxyzzzzz, g_0_x_xyz_xxzzzzz, g_0_x_xyz_xyyyyyy, g_0_x_xyz_xyyyyyyy, g_0_x_xyz_xyyyyyyz, g_0_x_xyz_xyyyyyz, g_0_x_xyz_xyyyyyzz, g_0_x_xyz_xyyyyzz, g_0_x_xyz_xyyyyzzz, g_0_x_xyz_xyyyzzz, g_0_x_xyz_xyyyzzzz, g_0_x_xyz_xyyzzzz, g_0_x_xyz_xyyzzzzz, g_0_x_xyz_xyzzzzz, g_0_x_xyz_xyzzzzzz, g_0_x_xyz_xzzzzzz, g_0_x_xyz_yyyyyyy, g_0_x_xyz_yyyyyyyy, g_0_x_xyz_yyyyyyyz, g_0_x_xyz_yyyyyyz, g_0_x_xyz_yyyyyyzz, g_0_x_xyz_yyyyyzz, g_0_x_xyz_yyyyyzzz, g_0_x_xyz_yyyyzzz, g_0_x_xyz_yyyyzzzz, g_0_x_xyz_yyyzzzz, g_0_x_xyz_yyyzzzzz, g_0_x_xyz_yyzzzzz, g_0_x_xyz_yyzzzzzz, g_0_x_xyz_yzzzzzz, g_0_x_xyz_yzzzzzzz, g_0_x_xyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_xxxxxxx[k] = -g_0_x_xyz_xxxxxxx[k] * ab_y + g_0_x_xyz_xxxxxxxy[k];

                g_0_x_xyyz_xxxxxxy[k] = -g_0_x_xyz_xxxxxxy[k] * ab_y + g_0_x_xyz_xxxxxxyy[k];

                g_0_x_xyyz_xxxxxxz[k] = -g_0_x_xyz_xxxxxxz[k] * ab_y + g_0_x_xyz_xxxxxxyz[k];

                g_0_x_xyyz_xxxxxyy[k] = -g_0_x_xyz_xxxxxyy[k] * ab_y + g_0_x_xyz_xxxxxyyy[k];

                g_0_x_xyyz_xxxxxyz[k] = -g_0_x_xyz_xxxxxyz[k] * ab_y + g_0_x_xyz_xxxxxyyz[k];

                g_0_x_xyyz_xxxxxzz[k] = -g_0_x_xyz_xxxxxzz[k] * ab_y + g_0_x_xyz_xxxxxyzz[k];

                g_0_x_xyyz_xxxxyyy[k] = -g_0_x_xyz_xxxxyyy[k] * ab_y + g_0_x_xyz_xxxxyyyy[k];

                g_0_x_xyyz_xxxxyyz[k] = -g_0_x_xyz_xxxxyyz[k] * ab_y + g_0_x_xyz_xxxxyyyz[k];

                g_0_x_xyyz_xxxxyzz[k] = -g_0_x_xyz_xxxxyzz[k] * ab_y + g_0_x_xyz_xxxxyyzz[k];

                g_0_x_xyyz_xxxxzzz[k] = -g_0_x_xyz_xxxxzzz[k] * ab_y + g_0_x_xyz_xxxxyzzz[k];

                g_0_x_xyyz_xxxyyyy[k] = -g_0_x_xyz_xxxyyyy[k] * ab_y + g_0_x_xyz_xxxyyyyy[k];

                g_0_x_xyyz_xxxyyyz[k] = -g_0_x_xyz_xxxyyyz[k] * ab_y + g_0_x_xyz_xxxyyyyz[k];

                g_0_x_xyyz_xxxyyzz[k] = -g_0_x_xyz_xxxyyzz[k] * ab_y + g_0_x_xyz_xxxyyyzz[k];

                g_0_x_xyyz_xxxyzzz[k] = -g_0_x_xyz_xxxyzzz[k] * ab_y + g_0_x_xyz_xxxyyzzz[k];

                g_0_x_xyyz_xxxzzzz[k] = -g_0_x_xyz_xxxzzzz[k] * ab_y + g_0_x_xyz_xxxyzzzz[k];

                g_0_x_xyyz_xxyyyyy[k] = -g_0_x_xyz_xxyyyyy[k] * ab_y + g_0_x_xyz_xxyyyyyy[k];

                g_0_x_xyyz_xxyyyyz[k] = -g_0_x_xyz_xxyyyyz[k] * ab_y + g_0_x_xyz_xxyyyyyz[k];

                g_0_x_xyyz_xxyyyzz[k] = -g_0_x_xyz_xxyyyzz[k] * ab_y + g_0_x_xyz_xxyyyyzz[k];

                g_0_x_xyyz_xxyyzzz[k] = -g_0_x_xyz_xxyyzzz[k] * ab_y + g_0_x_xyz_xxyyyzzz[k];

                g_0_x_xyyz_xxyzzzz[k] = -g_0_x_xyz_xxyzzzz[k] * ab_y + g_0_x_xyz_xxyyzzzz[k];

                g_0_x_xyyz_xxzzzzz[k] = -g_0_x_xyz_xxzzzzz[k] * ab_y + g_0_x_xyz_xxyzzzzz[k];

                g_0_x_xyyz_xyyyyyy[k] = -g_0_x_xyz_xyyyyyy[k] * ab_y + g_0_x_xyz_xyyyyyyy[k];

                g_0_x_xyyz_xyyyyyz[k] = -g_0_x_xyz_xyyyyyz[k] * ab_y + g_0_x_xyz_xyyyyyyz[k];

                g_0_x_xyyz_xyyyyzz[k] = -g_0_x_xyz_xyyyyzz[k] * ab_y + g_0_x_xyz_xyyyyyzz[k];

                g_0_x_xyyz_xyyyzzz[k] = -g_0_x_xyz_xyyyzzz[k] * ab_y + g_0_x_xyz_xyyyyzzz[k];

                g_0_x_xyyz_xyyzzzz[k] = -g_0_x_xyz_xyyzzzz[k] * ab_y + g_0_x_xyz_xyyyzzzz[k];

                g_0_x_xyyz_xyzzzzz[k] = -g_0_x_xyz_xyzzzzz[k] * ab_y + g_0_x_xyz_xyyzzzzz[k];

                g_0_x_xyyz_xzzzzzz[k] = -g_0_x_xyz_xzzzzzz[k] * ab_y + g_0_x_xyz_xyzzzzzz[k];

                g_0_x_xyyz_yyyyyyy[k] = -g_0_x_xyz_yyyyyyy[k] * ab_y + g_0_x_xyz_yyyyyyyy[k];

                g_0_x_xyyz_yyyyyyz[k] = -g_0_x_xyz_yyyyyyz[k] * ab_y + g_0_x_xyz_yyyyyyyz[k];

                g_0_x_xyyz_yyyyyzz[k] = -g_0_x_xyz_yyyyyzz[k] * ab_y + g_0_x_xyz_yyyyyyzz[k];

                g_0_x_xyyz_yyyyzzz[k] = -g_0_x_xyz_yyyyzzz[k] * ab_y + g_0_x_xyz_yyyyyzzz[k];

                g_0_x_xyyz_yyyzzzz[k] = -g_0_x_xyz_yyyzzzz[k] * ab_y + g_0_x_xyz_yyyyzzzz[k];

                g_0_x_xyyz_yyzzzzz[k] = -g_0_x_xyz_yyzzzzz[k] * ab_y + g_0_x_xyz_yyyzzzzz[k];

                g_0_x_xyyz_yzzzzzz[k] = -g_0_x_xyz_yzzzzzz[k] * ab_y + g_0_x_xyz_yyzzzzzz[k];

                g_0_x_xyyz_zzzzzzz[k] = -g_0_x_xyz_zzzzzzz[k] * ab_y + g_0_x_xyz_yzzzzzzz[k];
            }

            /// Set up 288-324 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyzz_xxxxxxx, g_0_x_xyzz_xxxxxxy, g_0_x_xyzz_xxxxxxz, g_0_x_xyzz_xxxxxyy, g_0_x_xyzz_xxxxxyz, g_0_x_xyzz_xxxxxzz, g_0_x_xyzz_xxxxyyy, g_0_x_xyzz_xxxxyyz, g_0_x_xyzz_xxxxyzz, g_0_x_xyzz_xxxxzzz, g_0_x_xyzz_xxxyyyy, g_0_x_xyzz_xxxyyyz, g_0_x_xyzz_xxxyyzz, g_0_x_xyzz_xxxyzzz, g_0_x_xyzz_xxxzzzz, g_0_x_xyzz_xxyyyyy, g_0_x_xyzz_xxyyyyz, g_0_x_xyzz_xxyyyzz, g_0_x_xyzz_xxyyzzz, g_0_x_xyzz_xxyzzzz, g_0_x_xyzz_xxzzzzz, g_0_x_xyzz_xyyyyyy, g_0_x_xyzz_xyyyyyz, g_0_x_xyzz_xyyyyzz, g_0_x_xyzz_xyyyzzz, g_0_x_xyzz_xyyzzzz, g_0_x_xyzz_xyzzzzz, g_0_x_xyzz_xzzzzzz, g_0_x_xyzz_yyyyyyy, g_0_x_xyzz_yyyyyyz, g_0_x_xyzz_yyyyyzz, g_0_x_xyzz_yyyyzzz, g_0_x_xyzz_yyyzzzz, g_0_x_xyzz_yyzzzzz, g_0_x_xyzz_yzzzzzz, g_0_x_xyzz_zzzzzzz, g_0_x_xzz_xxxxxxx, g_0_x_xzz_xxxxxxxy, g_0_x_xzz_xxxxxxy, g_0_x_xzz_xxxxxxyy, g_0_x_xzz_xxxxxxyz, g_0_x_xzz_xxxxxxz, g_0_x_xzz_xxxxxyy, g_0_x_xzz_xxxxxyyy, g_0_x_xzz_xxxxxyyz, g_0_x_xzz_xxxxxyz, g_0_x_xzz_xxxxxyzz, g_0_x_xzz_xxxxxzz, g_0_x_xzz_xxxxyyy, g_0_x_xzz_xxxxyyyy, g_0_x_xzz_xxxxyyyz, g_0_x_xzz_xxxxyyz, g_0_x_xzz_xxxxyyzz, g_0_x_xzz_xxxxyzz, g_0_x_xzz_xxxxyzzz, g_0_x_xzz_xxxxzzz, g_0_x_xzz_xxxyyyy, g_0_x_xzz_xxxyyyyy, g_0_x_xzz_xxxyyyyz, g_0_x_xzz_xxxyyyz, g_0_x_xzz_xxxyyyzz, g_0_x_xzz_xxxyyzz, g_0_x_xzz_xxxyyzzz, g_0_x_xzz_xxxyzzz, g_0_x_xzz_xxxyzzzz, g_0_x_xzz_xxxzzzz, g_0_x_xzz_xxyyyyy, g_0_x_xzz_xxyyyyyy, g_0_x_xzz_xxyyyyyz, g_0_x_xzz_xxyyyyz, g_0_x_xzz_xxyyyyzz, g_0_x_xzz_xxyyyzz, g_0_x_xzz_xxyyyzzz, g_0_x_xzz_xxyyzzz, g_0_x_xzz_xxyyzzzz, g_0_x_xzz_xxyzzzz, g_0_x_xzz_xxyzzzzz, g_0_x_xzz_xxzzzzz, g_0_x_xzz_xyyyyyy, g_0_x_xzz_xyyyyyyy, g_0_x_xzz_xyyyyyyz, g_0_x_xzz_xyyyyyz, g_0_x_xzz_xyyyyyzz, g_0_x_xzz_xyyyyzz, g_0_x_xzz_xyyyyzzz, g_0_x_xzz_xyyyzzz, g_0_x_xzz_xyyyzzzz, g_0_x_xzz_xyyzzzz, g_0_x_xzz_xyyzzzzz, g_0_x_xzz_xyzzzzz, g_0_x_xzz_xyzzzzzz, g_0_x_xzz_xzzzzzz, g_0_x_xzz_yyyyyyy, g_0_x_xzz_yyyyyyyy, g_0_x_xzz_yyyyyyyz, g_0_x_xzz_yyyyyyz, g_0_x_xzz_yyyyyyzz, g_0_x_xzz_yyyyyzz, g_0_x_xzz_yyyyyzzz, g_0_x_xzz_yyyyzzz, g_0_x_xzz_yyyyzzzz, g_0_x_xzz_yyyzzzz, g_0_x_xzz_yyyzzzzz, g_0_x_xzz_yyzzzzz, g_0_x_xzz_yyzzzzzz, g_0_x_xzz_yzzzzzz, g_0_x_xzz_yzzzzzzz, g_0_x_xzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_xxxxxxx[k] = -g_0_x_xzz_xxxxxxx[k] * ab_y + g_0_x_xzz_xxxxxxxy[k];

                g_0_x_xyzz_xxxxxxy[k] = -g_0_x_xzz_xxxxxxy[k] * ab_y + g_0_x_xzz_xxxxxxyy[k];

                g_0_x_xyzz_xxxxxxz[k] = -g_0_x_xzz_xxxxxxz[k] * ab_y + g_0_x_xzz_xxxxxxyz[k];

                g_0_x_xyzz_xxxxxyy[k] = -g_0_x_xzz_xxxxxyy[k] * ab_y + g_0_x_xzz_xxxxxyyy[k];

                g_0_x_xyzz_xxxxxyz[k] = -g_0_x_xzz_xxxxxyz[k] * ab_y + g_0_x_xzz_xxxxxyyz[k];

                g_0_x_xyzz_xxxxxzz[k] = -g_0_x_xzz_xxxxxzz[k] * ab_y + g_0_x_xzz_xxxxxyzz[k];

                g_0_x_xyzz_xxxxyyy[k] = -g_0_x_xzz_xxxxyyy[k] * ab_y + g_0_x_xzz_xxxxyyyy[k];

                g_0_x_xyzz_xxxxyyz[k] = -g_0_x_xzz_xxxxyyz[k] * ab_y + g_0_x_xzz_xxxxyyyz[k];

                g_0_x_xyzz_xxxxyzz[k] = -g_0_x_xzz_xxxxyzz[k] * ab_y + g_0_x_xzz_xxxxyyzz[k];

                g_0_x_xyzz_xxxxzzz[k] = -g_0_x_xzz_xxxxzzz[k] * ab_y + g_0_x_xzz_xxxxyzzz[k];

                g_0_x_xyzz_xxxyyyy[k] = -g_0_x_xzz_xxxyyyy[k] * ab_y + g_0_x_xzz_xxxyyyyy[k];

                g_0_x_xyzz_xxxyyyz[k] = -g_0_x_xzz_xxxyyyz[k] * ab_y + g_0_x_xzz_xxxyyyyz[k];

                g_0_x_xyzz_xxxyyzz[k] = -g_0_x_xzz_xxxyyzz[k] * ab_y + g_0_x_xzz_xxxyyyzz[k];

                g_0_x_xyzz_xxxyzzz[k] = -g_0_x_xzz_xxxyzzz[k] * ab_y + g_0_x_xzz_xxxyyzzz[k];

                g_0_x_xyzz_xxxzzzz[k] = -g_0_x_xzz_xxxzzzz[k] * ab_y + g_0_x_xzz_xxxyzzzz[k];

                g_0_x_xyzz_xxyyyyy[k] = -g_0_x_xzz_xxyyyyy[k] * ab_y + g_0_x_xzz_xxyyyyyy[k];

                g_0_x_xyzz_xxyyyyz[k] = -g_0_x_xzz_xxyyyyz[k] * ab_y + g_0_x_xzz_xxyyyyyz[k];

                g_0_x_xyzz_xxyyyzz[k] = -g_0_x_xzz_xxyyyzz[k] * ab_y + g_0_x_xzz_xxyyyyzz[k];

                g_0_x_xyzz_xxyyzzz[k] = -g_0_x_xzz_xxyyzzz[k] * ab_y + g_0_x_xzz_xxyyyzzz[k];

                g_0_x_xyzz_xxyzzzz[k] = -g_0_x_xzz_xxyzzzz[k] * ab_y + g_0_x_xzz_xxyyzzzz[k];

                g_0_x_xyzz_xxzzzzz[k] = -g_0_x_xzz_xxzzzzz[k] * ab_y + g_0_x_xzz_xxyzzzzz[k];

                g_0_x_xyzz_xyyyyyy[k] = -g_0_x_xzz_xyyyyyy[k] * ab_y + g_0_x_xzz_xyyyyyyy[k];

                g_0_x_xyzz_xyyyyyz[k] = -g_0_x_xzz_xyyyyyz[k] * ab_y + g_0_x_xzz_xyyyyyyz[k];

                g_0_x_xyzz_xyyyyzz[k] = -g_0_x_xzz_xyyyyzz[k] * ab_y + g_0_x_xzz_xyyyyyzz[k];

                g_0_x_xyzz_xyyyzzz[k] = -g_0_x_xzz_xyyyzzz[k] * ab_y + g_0_x_xzz_xyyyyzzz[k];

                g_0_x_xyzz_xyyzzzz[k] = -g_0_x_xzz_xyyzzzz[k] * ab_y + g_0_x_xzz_xyyyzzzz[k];

                g_0_x_xyzz_xyzzzzz[k] = -g_0_x_xzz_xyzzzzz[k] * ab_y + g_0_x_xzz_xyyzzzzz[k];

                g_0_x_xyzz_xzzzzzz[k] = -g_0_x_xzz_xzzzzzz[k] * ab_y + g_0_x_xzz_xyzzzzzz[k];

                g_0_x_xyzz_yyyyyyy[k] = -g_0_x_xzz_yyyyyyy[k] * ab_y + g_0_x_xzz_yyyyyyyy[k];

                g_0_x_xyzz_yyyyyyz[k] = -g_0_x_xzz_yyyyyyz[k] * ab_y + g_0_x_xzz_yyyyyyyz[k];

                g_0_x_xyzz_yyyyyzz[k] = -g_0_x_xzz_yyyyyzz[k] * ab_y + g_0_x_xzz_yyyyyyzz[k];

                g_0_x_xyzz_yyyyzzz[k] = -g_0_x_xzz_yyyyzzz[k] * ab_y + g_0_x_xzz_yyyyyzzz[k];

                g_0_x_xyzz_yyyzzzz[k] = -g_0_x_xzz_yyyzzzz[k] * ab_y + g_0_x_xzz_yyyyzzzz[k];

                g_0_x_xyzz_yyzzzzz[k] = -g_0_x_xzz_yyzzzzz[k] * ab_y + g_0_x_xzz_yyyzzzzz[k];

                g_0_x_xyzz_yzzzzzz[k] = -g_0_x_xzz_yzzzzzz[k] * ab_y + g_0_x_xzz_yyzzzzzz[k];

                g_0_x_xyzz_zzzzzzz[k] = -g_0_x_xzz_zzzzzzz[k] * ab_y + g_0_x_xzz_yzzzzzzz[k];
            }

            /// Set up 324-360 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xzz_xxxxxxx, g_0_x_xzz_xxxxxxxz, g_0_x_xzz_xxxxxxy, g_0_x_xzz_xxxxxxyz, g_0_x_xzz_xxxxxxz, g_0_x_xzz_xxxxxxzz, g_0_x_xzz_xxxxxyy, g_0_x_xzz_xxxxxyyz, g_0_x_xzz_xxxxxyz, g_0_x_xzz_xxxxxyzz, g_0_x_xzz_xxxxxzz, g_0_x_xzz_xxxxxzzz, g_0_x_xzz_xxxxyyy, g_0_x_xzz_xxxxyyyz, g_0_x_xzz_xxxxyyz, g_0_x_xzz_xxxxyyzz, g_0_x_xzz_xxxxyzz, g_0_x_xzz_xxxxyzzz, g_0_x_xzz_xxxxzzz, g_0_x_xzz_xxxxzzzz, g_0_x_xzz_xxxyyyy, g_0_x_xzz_xxxyyyyz, g_0_x_xzz_xxxyyyz, g_0_x_xzz_xxxyyyzz, g_0_x_xzz_xxxyyzz, g_0_x_xzz_xxxyyzzz, g_0_x_xzz_xxxyzzz, g_0_x_xzz_xxxyzzzz, g_0_x_xzz_xxxzzzz, g_0_x_xzz_xxxzzzzz, g_0_x_xzz_xxyyyyy, g_0_x_xzz_xxyyyyyz, g_0_x_xzz_xxyyyyz, g_0_x_xzz_xxyyyyzz, g_0_x_xzz_xxyyyzz, g_0_x_xzz_xxyyyzzz, g_0_x_xzz_xxyyzzz, g_0_x_xzz_xxyyzzzz, g_0_x_xzz_xxyzzzz, g_0_x_xzz_xxyzzzzz, g_0_x_xzz_xxzzzzz, g_0_x_xzz_xxzzzzzz, g_0_x_xzz_xyyyyyy, g_0_x_xzz_xyyyyyyz, g_0_x_xzz_xyyyyyz, g_0_x_xzz_xyyyyyzz, g_0_x_xzz_xyyyyzz, g_0_x_xzz_xyyyyzzz, g_0_x_xzz_xyyyzzz, g_0_x_xzz_xyyyzzzz, g_0_x_xzz_xyyzzzz, g_0_x_xzz_xyyzzzzz, g_0_x_xzz_xyzzzzz, g_0_x_xzz_xyzzzzzz, g_0_x_xzz_xzzzzzz, g_0_x_xzz_xzzzzzzz, g_0_x_xzz_yyyyyyy, g_0_x_xzz_yyyyyyyz, g_0_x_xzz_yyyyyyz, g_0_x_xzz_yyyyyyzz, g_0_x_xzz_yyyyyzz, g_0_x_xzz_yyyyyzzz, g_0_x_xzz_yyyyzzz, g_0_x_xzz_yyyyzzzz, g_0_x_xzz_yyyzzzz, g_0_x_xzz_yyyzzzzz, g_0_x_xzz_yyzzzzz, g_0_x_xzz_yyzzzzzz, g_0_x_xzz_yzzzzzz, g_0_x_xzz_yzzzzzzz, g_0_x_xzz_zzzzzzz, g_0_x_xzz_zzzzzzzz, g_0_x_xzzz_xxxxxxx, g_0_x_xzzz_xxxxxxy, g_0_x_xzzz_xxxxxxz, g_0_x_xzzz_xxxxxyy, g_0_x_xzzz_xxxxxyz, g_0_x_xzzz_xxxxxzz, g_0_x_xzzz_xxxxyyy, g_0_x_xzzz_xxxxyyz, g_0_x_xzzz_xxxxyzz, g_0_x_xzzz_xxxxzzz, g_0_x_xzzz_xxxyyyy, g_0_x_xzzz_xxxyyyz, g_0_x_xzzz_xxxyyzz, g_0_x_xzzz_xxxyzzz, g_0_x_xzzz_xxxzzzz, g_0_x_xzzz_xxyyyyy, g_0_x_xzzz_xxyyyyz, g_0_x_xzzz_xxyyyzz, g_0_x_xzzz_xxyyzzz, g_0_x_xzzz_xxyzzzz, g_0_x_xzzz_xxzzzzz, g_0_x_xzzz_xyyyyyy, g_0_x_xzzz_xyyyyyz, g_0_x_xzzz_xyyyyzz, g_0_x_xzzz_xyyyzzz, g_0_x_xzzz_xyyzzzz, g_0_x_xzzz_xyzzzzz, g_0_x_xzzz_xzzzzzz, g_0_x_xzzz_yyyyyyy, g_0_x_xzzz_yyyyyyz, g_0_x_xzzz_yyyyyzz, g_0_x_xzzz_yyyyzzz, g_0_x_xzzz_yyyzzzz, g_0_x_xzzz_yyzzzzz, g_0_x_xzzz_yzzzzzz, g_0_x_xzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_xxxxxxx[k] = -g_0_x_xzz_xxxxxxx[k] * ab_z + g_0_x_xzz_xxxxxxxz[k];

                g_0_x_xzzz_xxxxxxy[k] = -g_0_x_xzz_xxxxxxy[k] * ab_z + g_0_x_xzz_xxxxxxyz[k];

                g_0_x_xzzz_xxxxxxz[k] = -g_0_x_xzz_xxxxxxz[k] * ab_z + g_0_x_xzz_xxxxxxzz[k];

                g_0_x_xzzz_xxxxxyy[k] = -g_0_x_xzz_xxxxxyy[k] * ab_z + g_0_x_xzz_xxxxxyyz[k];

                g_0_x_xzzz_xxxxxyz[k] = -g_0_x_xzz_xxxxxyz[k] * ab_z + g_0_x_xzz_xxxxxyzz[k];

                g_0_x_xzzz_xxxxxzz[k] = -g_0_x_xzz_xxxxxzz[k] * ab_z + g_0_x_xzz_xxxxxzzz[k];

                g_0_x_xzzz_xxxxyyy[k] = -g_0_x_xzz_xxxxyyy[k] * ab_z + g_0_x_xzz_xxxxyyyz[k];

                g_0_x_xzzz_xxxxyyz[k] = -g_0_x_xzz_xxxxyyz[k] * ab_z + g_0_x_xzz_xxxxyyzz[k];

                g_0_x_xzzz_xxxxyzz[k] = -g_0_x_xzz_xxxxyzz[k] * ab_z + g_0_x_xzz_xxxxyzzz[k];

                g_0_x_xzzz_xxxxzzz[k] = -g_0_x_xzz_xxxxzzz[k] * ab_z + g_0_x_xzz_xxxxzzzz[k];

                g_0_x_xzzz_xxxyyyy[k] = -g_0_x_xzz_xxxyyyy[k] * ab_z + g_0_x_xzz_xxxyyyyz[k];

                g_0_x_xzzz_xxxyyyz[k] = -g_0_x_xzz_xxxyyyz[k] * ab_z + g_0_x_xzz_xxxyyyzz[k];

                g_0_x_xzzz_xxxyyzz[k] = -g_0_x_xzz_xxxyyzz[k] * ab_z + g_0_x_xzz_xxxyyzzz[k];

                g_0_x_xzzz_xxxyzzz[k] = -g_0_x_xzz_xxxyzzz[k] * ab_z + g_0_x_xzz_xxxyzzzz[k];

                g_0_x_xzzz_xxxzzzz[k] = -g_0_x_xzz_xxxzzzz[k] * ab_z + g_0_x_xzz_xxxzzzzz[k];

                g_0_x_xzzz_xxyyyyy[k] = -g_0_x_xzz_xxyyyyy[k] * ab_z + g_0_x_xzz_xxyyyyyz[k];

                g_0_x_xzzz_xxyyyyz[k] = -g_0_x_xzz_xxyyyyz[k] * ab_z + g_0_x_xzz_xxyyyyzz[k];

                g_0_x_xzzz_xxyyyzz[k] = -g_0_x_xzz_xxyyyzz[k] * ab_z + g_0_x_xzz_xxyyyzzz[k];

                g_0_x_xzzz_xxyyzzz[k] = -g_0_x_xzz_xxyyzzz[k] * ab_z + g_0_x_xzz_xxyyzzzz[k];

                g_0_x_xzzz_xxyzzzz[k] = -g_0_x_xzz_xxyzzzz[k] * ab_z + g_0_x_xzz_xxyzzzzz[k];

                g_0_x_xzzz_xxzzzzz[k] = -g_0_x_xzz_xxzzzzz[k] * ab_z + g_0_x_xzz_xxzzzzzz[k];

                g_0_x_xzzz_xyyyyyy[k] = -g_0_x_xzz_xyyyyyy[k] * ab_z + g_0_x_xzz_xyyyyyyz[k];

                g_0_x_xzzz_xyyyyyz[k] = -g_0_x_xzz_xyyyyyz[k] * ab_z + g_0_x_xzz_xyyyyyzz[k];

                g_0_x_xzzz_xyyyyzz[k] = -g_0_x_xzz_xyyyyzz[k] * ab_z + g_0_x_xzz_xyyyyzzz[k];

                g_0_x_xzzz_xyyyzzz[k] = -g_0_x_xzz_xyyyzzz[k] * ab_z + g_0_x_xzz_xyyyzzzz[k];

                g_0_x_xzzz_xyyzzzz[k] = -g_0_x_xzz_xyyzzzz[k] * ab_z + g_0_x_xzz_xyyzzzzz[k];

                g_0_x_xzzz_xyzzzzz[k] = -g_0_x_xzz_xyzzzzz[k] * ab_z + g_0_x_xzz_xyzzzzzz[k];

                g_0_x_xzzz_xzzzzzz[k] = -g_0_x_xzz_xzzzzzz[k] * ab_z + g_0_x_xzz_xzzzzzzz[k];

                g_0_x_xzzz_yyyyyyy[k] = -g_0_x_xzz_yyyyyyy[k] * ab_z + g_0_x_xzz_yyyyyyyz[k];

                g_0_x_xzzz_yyyyyyz[k] = -g_0_x_xzz_yyyyyyz[k] * ab_z + g_0_x_xzz_yyyyyyzz[k];

                g_0_x_xzzz_yyyyyzz[k] = -g_0_x_xzz_yyyyyzz[k] * ab_z + g_0_x_xzz_yyyyyzzz[k];

                g_0_x_xzzz_yyyyzzz[k] = -g_0_x_xzz_yyyyzzz[k] * ab_z + g_0_x_xzz_yyyyzzzz[k];

                g_0_x_xzzz_yyyzzzz[k] = -g_0_x_xzz_yyyzzzz[k] * ab_z + g_0_x_xzz_yyyzzzzz[k];

                g_0_x_xzzz_yyzzzzz[k] = -g_0_x_xzz_yyzzzzz[k] * ab_z + g_0_x_xzz_yyzzzzzz[k];

                g_0_x_xzzz_yzzzzzz[k] = -g_0_x_xzz_yzzzzzz[k] * ab_z + g_0_x_xzz_yzzzzzzz[k];

                g_0_x_xzzz_zzzzzzz[k] = -g_0_x_xzz_zzzzzzz[k] * ab_z + g_0_x_xzz_zzzzzzzz[k];
            }

            /// Set up 360-396 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyy_xxxxxxx, g_0_x_yyy_xxxxxxxy, g_0_x_yyy_xxxxxxy, g_0_x_yyy_xxxxxxyy, g_0_x_yyy_xxxxxxyz, g_0_x_yyy_xxxxxxz, g_0_x_yyy_xxxxxyy, g_0_x_yyy_xxxxxyyy, g_0_x_yyy_xxxxxyyz, g_0_x_yyy_xxxxxyz, g_0_x_yyy_xxxxxyzz, g_0_x_yyy_xxxxxzz, g_0_x_yyy_xxxxyyy, g_0_x_yyy_xxxxyyyy, g_0_x_yyy_xxxxyyyz, g_0_x_yyy_xxxxyyz, g_0_x_yyy_xxxxyyzz, g_0_x_yyy_xxxxyzz, g_0_x_yyy_xxxxyzzz, g_0_x_yyy_xxxxzzz, g_0_x_yyy_xxxyyyy, g_0_x_yyy_xxxyyyyy, g_0_x_yyy_xxxyyyyz, g_0_x_yyy_xxxyyyz, g_0_x_yyy_xxxyyyzz, g_0_x_yyy_xxxyyzz, g_0_x_yyy_xxxyyzzz, g_0_x_yyy_xxxyzzz, g_0_x_yyy_xxxyzzzz, g_0_x_yyy_xxxzzzz, g_0_x_yyy_xxyyyyy, g_0_x_yyy_xxyyyyyy, g_0_x_yyy_xxyyyyyz, g_0_x_yyy_xxyyyyz, g_0_x_yyy_xxyyyyzz, g_0_x_yyy_xxyyyzz, g_0_x_yyy_xxyyyzzz, g_0_x_yyy_xxyyzzz, g_0_x_yyy_xxyyzzzz, g_0_x_yyy_xxyzzzz, g_0_x_yyy_xxyzzzzz, g_0_x_yyy_xxzzzzz, g_0_x_yyy_xyyyyyy, g_0_x_yyy_xyyyyyyy, g_0_x_yyy_xyyyyyyz, g_0_x_yyy_xyyyyyz, g_0_x_yyy_xyyyyyzz, g_0_x_yyy_xyyyyzz, g_0_x_yyy_xyyyyzzz, g_0_x_yyy_xyyyzzz, g_0_x_yyy_xyyyzzzz, g_0_x_yyy_xyyzzzz, g_0_x_yyy_xyyzzzzz, g_0_x_yyy_xyzzzzz, g_0_x_yyy_xyzzzzzz, g_0_x_yyy_xzzzzzz, g_0_x_yyy_yyyyyyy, g_0_x_yyy_yyyyyyyy, g_0_x_yyy_yyyyyyyz, g_0_x_yyy_yyyyyyz, g_0_x_yyy_yyyyyyzz, g_0_x_yyy_yyyyyzz, g_0_x_yyy_yyyyyzzz, g_0_x_yyy_yyyyzzz, g_0_x_yyy_yyyyzzzz, g_0_x_yyy_yyyzzzz, g_0_x_yyy_yyyzzzzz, g_0_x_yyy_yyzzzzz, g_0_x_yyy_yyzzzzzz, g_0_x_yyy_yzzzzzz, g_0_x_yyy_yzzzzzzz, g_0_x_yyy_zzzzzzz, g_0_x_yyyy_xxxxxxx, g_0_x_yyyy_xxxxxxy, g_0_x_yyyy_xxxxxxz, g_0_x_yyyy_xxxxxyy, g_0_x_yyyy_xxxxxyz, g_0_x_yyyy_xxxxxzz, g_0_x_yyyy_xxxxyyy, g_0_x_yyyy_xxxxyyz, g_0_x_yyyy_xxxxyzz, g_0_x_yyyy_xxxxzzz, g_0_x_yyyy_xxxyyyy, g_0_x_yyyy_xxxyyyz, g_0_x_yyyy_xxxyyzz, g_0_x_yyyy_xxxyzzz, g_0_x_yyyy_xxxzzzz, g_0_x_yyyy_xxyyyyy, g_0_x_yyyy_xxyyyyz, g_0_x_yyyy_xxyyyzz, g_0_x_yyyy_xxyyzzz, g_0_x_yyyy_xxyzzzz, g_0_x_yyyy_xxzzzzz, g_0_x_yyyy_xyyyyyy, g_0_x_yyyy_xyyyyyz, g_0_x_yyyy_xyyyyzz, g_0_x_yyyy_xyyyzzz, g_0_x_yyyy_xyyzzzz, g_0_x_yyyy_xyzzzzz, g_0_x_yyyy_xzzzzzz, g_0_x_yyyy_yyyyyyy, g_0_x_yyyy_yyyyyyz, g_0_x_yyyy_yyyyyzz, g_0_x_yyyy_yyyyzzz, g_0_x_yyyy_yyyzzzz, g_0_x_yyyy_yyzzzzz, g_0_x_yyyy_yzzzzzz, g_0_x_yyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_xxxxxxx[k] = -g_0_x_yyy_xxxxxxx[k] * ab_y + g_0_x_yyy_xxxxxxxy[k];

                g_0_x_yyyy_xxxxxxy[k] = -g_0_x_yyy_xxxxxxy[k] * ab_y + g_0_x_yyy_xxxxxxyy[k];

                g_0_x_yyyy_xxxxxxz[k] = -g_0_x_yyy_xxxxxxz[k] * ab_y + g_0_x_yyy_xxxxxxyz[k];

                g_0_x_yyyy_xxxxxyy[k] = -g_0_x_yyy_xxxxxyy[k] * ab_y + g_0_x_yyy_xxxxxyyy[k];

                g_0_x_yyyy_xxxxxyz[k] = -g_0_x_yyy_xxxxxyz[k] * ab_y + g_0_x_yyy_xxxxxyyz[k];

                g_0_x_yyyy_xxxxxzz[k] = -g_0_x_yyy_xxxxxzz[k] * ab_y + g_0_x_yyy_xxxxxyzz[k];

                g_0_x_yyyy_xxxxyyy[k] = -g_0_x_yyy_xxxxyyy[k] * ab_y + g_0_x_yyy_xxxxyyyy[k];

                g_0_x_yyyy_xxxxyyz[k] = -g_0_x_yyy_xxxxyyz[k] * ab_y + g_0_x_yyy_xxxxyyyz[k];

                g_0_x_yyyy_xxxxyzz[k] = -g_0_x_yyy_xxxxyzz[k] * ab_y + g_0_x_yyy_xxxxyyzz[k];

                g_0_x_yyyy_xxxxzzz[k] = -g_0_x_yyy_xxxxzzz[k] * ab_y + g_0_x_yyy_xxxxyzzz[k];

                g_0_x_yyyy_xxxyyyy[k] = -g_0_x_yyy_xxxyyyy[k] * ab_y + g_0_x_yyy_xxxyyyyy[k];

                g_0_x_yyyy_xxxyyyz[k] = -g_0_x_yyy_xxxyyyz[k] * ab_y + g_0_x_yyy_xxxyyyyz[k];

                g_0_x_yyyy_xxxyyzz[k] = -g_0_x_yyy_xxxyyzz[k] * ab_y + g_0_x_yyy_xxxyyyzz[k];

                g_0_x_yyyy_xxxyzzz[k] = -g_0_x_yyy_xxxyzzz[k] * ab_y + g_0_x_yyy_xxxyyzzz[k];

                g_0_x_yyyy_xxxzzzz[k] = -g_0_x_yyy_xxxzzzz[k] * ab_y + g_0_x_yyy_xxxyzzzz[k];

                g_0_x_yyyy_xxyyyyy[k] = -g_0_x_yyy_xxyyyyy[k] * ab_y + g_0_x_yyy_xxyyyyyy[k];

                g_0_x_yyyy_xxyyyyz[k] = -g_0_x_yyy_xxyyyyz[k] * ab_y + g_0_x_yyy_xxyyyyyz[k];

                g_0_x_yyyy_xxyyyzz[k] = -g_0_x_yyy_xxyyyzz[k] * ab_y + g_0_x_yyy_xxyyyyzz[k];

                g_0_x_yyyy_xxyyzzz[k] = -g_0_x_yyy_xxyyzzz[k] * ab_y + g_0_x_yyy_xxyyyzzz[k];

                g_0_x_yyyy_xxyzzzz[k] = -g_0_x_yyy_xxyzzzz[k] * ab_y + g_0_x_yyy_xxyyzzzz[k];

                g_0_x_yyyy_xxzzzzz[k] = -g_0_x_yyy_xxzzzzz[k] * ab_y + g_0_x_yyy_xxyzzzzz[k];

                g_0_x_yyyy_xyyyyyy[k] = -g_0_x_yyy_xyyyyyy[k] * ab_y + g_0_x_yyy_xyyyyyyy[k];

                g_0_x_yyyy_xyyyyyz[k] = -g_0_x_yyy_xyyyyyz[k] * ab_y + g_0_x_yyy_xyyyyyyz[k];

                g_0_x_yyyy_xyyyyzz[k] = -g_0_x_yyy_xyyyyzz[k] * ab_y + g_0_x_yyy_xyyyyyzz[k];

                g_0_x_yyyy_xyyyzzz[k] = -g_0_x_yyy_xyyyzzz[k] * ab_y + g_0_x_yyy_xyyyyzzz[k];

                g_0_x_yyyy_xyyzzzz[k] = -g_0_x_yyy_xyyzzzz[k] * ab_y + g_0_x_yyy_xyyyzzzz[k];

                g_0_x_yyyy_xyzzzzz[k] = -g_0_x_yyy_xyzzzzz[k] * ab_y + g_0_x_yyy_xyyzzzzz[k];

                g_0_x_yyyy_xzzzzzz[k] = -g_0_x_yyy_xzzzzzz[k] * ab_y + g_0_x_yyy_xyzzzzzz[k];

                g_0_x_yyyy_yyyyyyy[k] = -g_0_x_yyy_yyyyyyy[k] * ab_y + g_0_x_yyy_yyyyyyyy[k];

                g_0_x_yyyy_yyyyyyz[k] = -g_0_x_yyy_yyyyyyz[k] * ab_y + g_0_x_yyy_yyyyyyyz[k];

                g_0_x_yyyy_yyyyyzz[k] = -g_0_x_yyy_yyyyyzz[k] * ab_y + g_0_x_yyy_yyyyyyzz[k];

                g_0_x_yyyy_yyyyzzz[k] = -g_0_x_yyy_yyyyzzz[k] * ab_y + g_0_x_yyy_yyyyyzzz[k];

                g_0_x_yyyy_yyyzzzz[k] = -g_0_x_yyy_yyyzzzz[k] * ab_y + g_0_x_yyy_yyyyzzzz[k];

                g_0_x_yyyy_yyzzzzz[k] = -g_0_x_yyy_yyzzzzz[k] * ab_y + g_0_x_yyy_yyyzzzzz[k];

                g_0_x_yyyy_yzzzzzz[k] = -g_0_x_yyy_yzzzzzz[k] * ab_y + g_0_x_yyy_yyzzzzzz[k];

                g_0_x_yyyy_zzzzzzz[k] = -g_0_x_yyy_zzzzzzz[k] * ab_y + g_0_x_yyy_yzzzzzzz[k];
            }

            /// Set up 396-432 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyyz_xxxxxxx, g_0_x_yyyz_xxxxxxy, g_0_x_yyyz_xxxxxxz, g_0_x_yyyz_xxxxxyy, g_0_x_yyyz_xxxxxyz, g_0_x_yyyz_xxxxxzz, g_0_x_yyyz_xxxxyyy, g_0_x_yyyz_xxxxyyz, g_0_x_yyyz_xxxxyzz, g_0_x_yyyz_xxxxzzz, g_0_x_yyyz_xxxyyyy, g_0_x_yyyz_xxxyyyz, g_0_x_yyyz_xxxyyzz, g_0_x_yyyz_xxxyzzz, g_0_x_yyyz_xxxzzzz, g_0_x_yyyz_xxyyyyy, g_0_x_yyyz_xxyyyyz, g_0_x_yyyz_xxyyyzz, g_0_x_yyyz_xxyyzzz, g_0_x_yyyz_xxyzzzz, g_0_x_yyyz_xxzzzzz, g_0_x_yyyz_xyyyyyy, g_0_x_yyyz_xyyyyyz, g_0_x_yyyz_xyyyyzz, g_0_x_yyyz_xyyyzzz, g_0_x_yyyz_xyyzzzz, g_0_x_yyyz_xyzzzzz, g_0_x_yyyz_xzzzzzz, g_0_x_yyyz_yyyyyyy, g_0_x_yyyz_yyyyyyz, g_0_x_yyyz_yyyyyzz, g_0_x_yyyz_yyyyzzz, g_0_x_yyyz_yyyzzzz, g_0_x_yyyz_yyzzzzz, g_0_x_yyyz_yzzzzzz, g_0_x_yyyz_zzzzzzz, g_0_x_yyz_xxxxxxx, g_0_x_yyz_xxxxxxxy, g_0_x_yyz_xxxxxxy, g_0_x_yyz_xxxxxxyy, g_0_x_yyz_xxxxxxyz, g_0_x_yyz_xxxxxxz, g_0_x_yyz_xxxxxyy, g_0_x_yyz_xxxxxyyy, g_0_x_yyz_xxxxxyyz, g_0_x_yyz_xxxxxyz, g_0_x_yyz_xxxxxyzz, g_0_x_yyz_xxxxxzz, g_0_x_yyz_xxxxyyy, g_0_x_yyz_xxxxyyyy, g_0_x_yyz_xxxxyyyz, g_0_x_yyz_xxxxyyz, g_0_x_yyz_xxxxyyzz, g_0_x_yyz_xxxxyzz, g_0_x_yyz_xxxxyzzz, g_0_x_yyz_xxxxzzz, g_0_x_yyz_xxxyyyy, g_0_x_yyz_xxxyyyyy, g_0_x_yyz_xxxyyyyz, g_0_x_yyz_xxxyyyz, g_0_x_yyz_xxxyyyzz, g_0_x_yyz_xxxyyzz, g_0_x_yyz_xxxyyzzz, g_0_x_yyz_xxxyzzz, g_0_x_yyz_xxxyzzzz, g_0_x_yyz_xxxzzzz, g_0_x_yyz_xxyyyyy, g_0_x_yyz_xxyyyyyy, g_0_x_yyz_xxyyyyyz, g_0_x_yyz_xxyyyyz, g_0_x_yyz_xxyyyyzz, g_0_x_yyz_xxyyyzz, g_0_x_yyz_xxyyyzzz, g_0_x_yyz_xxyyzzz, g_0_x_yyz_xxyyzzzz, g_0_x_yyz_xxyzzzz, g_0_x_yyz_xxyzzzzz, g_0_x_yyz_xxzzzzz, g_0_x_yyz_xyyyyyy, g_0_x_yyz_xyyyyyyy, g_0_x_yyz_xyyyyyyz, g_0_x_yyz_xyyyyyz, g_0_x_yyz_xyyyyyzz, g_0_x_yyz_xyyyyzz, g_0_x_yyz_xyyyyzzz, g_0_x_yyz_xyyyzzz, g_0_x_yyz_xyyyzzzz, g_0_x_yyz_xyyzzzz, g_0_x_yyz_xyyzzzzz, g_0_x_yyz_xyzzzzz, g_0_x_yyz_xyzzzzzz, g_0_x_yyz_xzzzzzz, g_0_x_yyz_yyyyyyy, g_0_x_yyz_yyyyyyyy, g_0_x_yyz_yyyyyyyz, g_0_x_yyz_yyyyyyz, g_0_x_yyz_yyyyyyzz, g_0_x_yyz_yyyyyzz, g_0_x_yyz_yyyyyzzz, g_0_x_yyz_yyyyzzz, g_0_x_yyz_yyyyzzzz, g_0_x_yyz_yyyzzzz, g_0_x_yyz_yyyzzzzz, g_0_x_yyz_yyzzzzz, g_0_x_yyz_yyzzzzzz, g_0_x_yyz_yzzzzzz, g_0_x_yyz_yzzzzzzz, g_0_x_yyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_xxxxxxx[k] = -g_0_x_yyz_xxxxxxx[k] * ab_y + g_0_x_yyz_xxxxxxxy[k];

                g_0_x_yyyz_xxxxxxy[k] = -g_0_x_yyz_xxxxxxy[k] * ab_y + g_0_x_yyz_xxxxxxyy[k];

                g_0_x_yyyz_xxxxxxz[k] = -g_0_x_yyz_xxxxxxz[k] * ab_y + g_0_x_yyz_xxxxxxyz[k];

                g_0_x_yyyz_xxxxxyy[k] = -g_0_x_yyz_xxxxxyy[k] * ab_y + g_0_x_yyz_xxxxxyyy[k];

                g_0_x_yyyz_xxxxxyz[k] = -g_0_x_yyz_xxxxxyz[k] * ab_y + g_0_x_yyz_xxxxxyyz[k];

                g_0_x_yyyz_xxxxxzz[k] = -g_0_x_yyz_xxxxxzz[k] * ab_y + g_0_x_yyz_xxxxxyzz[k];

                g_0_x_yyyz_xxxxyyy[k] = -g_0_x_yyz_xxxxyyy[k] * ab_y + g_0_x_yyz_xxxxyyyy[k];

                g_0_x_yyyz_xxxxyyz[k] = -g_0_x_yyz_xxxxyyz[k] * ab_y + g_0_x_yyz_xxxxyyyz[k];

                g_0_x_yyyz_xxxxyzz[k] = -g_0_x_yyz_xxxxyzz[k] * ab_y + g_0_x_yyz_xxxxyyzz[k];

                g_0_x_yyyz_xxxxzzz[k] = -g_0_x_yyz_xxxxzzz[k] * ab_y + g_0_x_yyz_xxxxyzzz[k];

                g_0_x_yyyz_xxxyyyy[k] = -g_0_x_yyz_xxxyyyy[k] * ab_y + g_0_x_yyz_xxxyyyyy[k];

                g_0_x_yyyz_xxxyyyz[k] = -g_0_x_yyz_xxxyyyz[k] * ab_y + g_0_x_yyz_xxxyyyyz[k];

                g_0_x_yyyz_xxxyyzz[k] = -g_0_x_yyz_xxxyyzz[k] * ab_y + g_0_x_yyz_xxxyyyzz[k];

                g_0_x_yyyz_xxxyzzz[k] = -g_0_x_yyz_xxxyzzz[k] * ab_y + g_0_x_yyz_xxxyyzzz[k];

                g_0_x_yyyz_xxxzzzz[k] = -g_0_x_yyz_xxxzzzz[k] * ab_y + g_0_x_yyz_xxxyzzzz[k];

                g_0_x_yyyz_xxyyyyy[k] = -g_0_x_yyz_xxyyyyy[k] * ab_y + g_0_x_yyz_xxyyyyyy[k];

                g_0_x_yyyz_xxyyyyz[k] = -g_0_x_yyz_xxyyyyz[k] * ab_y + g_0_x_yyz_xxyyyyyz[k];

                g_0_x_yyyz_xxyyyzz[k] = -g_0_x_yyz_xxyyyzz[k] * ab_y + g_0_x_yyz_xxyyyyzz[k];

                g_0_x_yyyz_xxyyzzz[k] = -g_0_x_yyz_xxyyzzz[k] * ab_y + g_0_x_yyz_xxyyyzzz[k];

                g_0_x_yyyz_xxyzzzz[k] = -g_0_x_yyz_xxyzzzz[k] * ab_y + g_0_x_yyz_xxyyzzzz[k];

                g_0_x_yyyz_xxzzzzz[k] = -g_0_x_yyz_xxzzzzz[k] * ab_y + g_0_x_yyz_xxyzzzzz[k];

                g_0_x_yyyz_xyyyyyy[k] = -g_0_x_yyz_xyyyyyy[k] * ab_y + g_0_x_yyz_xyyyyyyy[k];

                g_0_x_yyyz_xyyyyyz[k] = -g_0_x_yyz_xyyyyyz[k] * ab_y + g_0_x_yyz_xyyyyyyz[k];

                g_0_x_yyyz_xyyyyzz[k] = -g_0_x_yyz_xyyyyzz[k] * ab_y + g_0_x_yyz_xyyyyyzz[k];

                g_0_x_yyyz_xyyyzzz[k] = -g_0_x_yyz_xyyyzzz[k] * ab_y + g_0_x_yyz_xyyyyzzz[k];

                g_0_x_yyyz_xyyzzzz[k] = -g_0_x_yyz_xyyzzzz[k] * ab_y + g_0_x_yyz_xyyyzzzz[k];

                g_0_x_yyyz_xyzzzzz[k] = -g_0_x_yyz_xyzzzzz[k] * ab_y + g_0_x_yyz_xyyzzzzz[k];

                g_0_x_yyyz_xzzzzzz[k] = -g_0_x_yyz_xzzzzzz[k] * ab_y + g_0_x_yyz_xyzzzzzz[k];

                g_0_x_yyyz_yyyyyyy[k] = -g_0_x_yyz_yyyyyyy[k] * ab_y + g_0_x_yyz_yyyyyyyy[k];

                g_0_x_yyyz_yyyyyyz[k] = -g_0_x_yyz_yyyyyyz[k] * ab_y + g_0_x_yyz_yyyyyyyz[k];

                g_0_x_yyyz_yyyyyzz[k] = -g_0_x_yyz_yyyyyzz[k] * ab_y + g_0_x_yyz_yyyyyyzz[k];

                g_0_x_yyyz_yyyyzzz[k] = -g_0_x_yyz_yyyyzzz[k] * ab_y + g_0_x_yyz_yyyyyzzz[k];

                g_0_x_yyyz_yyyzzzz[k] = -g_0_x_yyz_yyyzzzz[k] * ab_y + g_0_x_yyz_yyyyzzzz[k];

                g_0_x_yyyz_yyzzzzz[k] = -g_0_x_yyz_yyzzzzz[k] * ab_y + g_0_x_yyz_yyyzzzzz[k];

                g_0_x_yyyz_yzzzzzz[k] = -g_0_x_yyz_yzzzzzz[k] * ab_y + g_0_x_yyz_yyzzzzzz[k];

                g_0_x_yyyz_zzzzzzz[k] = -g_0_x_yyz_zzzzzzz[k] * ab_y + g_0_x_yyz_yzzzzzzz[k];
            }

            /// Set up 432-468 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyzz_xxxxxxx, g_0_x_yyzz_xxxxxxy, g_0_x_yyzz_xxxxxxz, g_0_x_yyzz_xxxxxyy, g_0_x_yyzz_xxxxxyz, g_0_x_yyzz_xxxxxzz, g_0_x_yyzz_xxxxyyy, g_0_x_yyzz_xxxxyyz, g_0_x_yyzz_xxxxyzz, g_0_x_yyzz_xxxxzzz, g_0_x_yyzz_xxxyyyy, g_0_x_yyzz_xxxyyyz, g_0_x_yyzz_xxxyyzz, g_0_x_yyzz_xxxyzzz, g_0_x_yyzz_xxxzzzz, g_0_x_yyzz_xxyyyyy, g_0_x_yyzz_xxyyyyz, g_0_x_yyzz_xxyyyzz, g_0_x_yyzz_xxyyzzz, g_0_x_yyzz_xxyzzzz, g_0_x_yyzz_xxzzzzz, g_0_x_yyzz_xyyyyyy, g_0_x_yyzz_xyyyyyz, g_0_x_yyzz_xyyyyzz, g_0_x_yyzz_xyyyzzz, g_0_x_yyzz_xyyzzzz, g_0_x_yyzz_xyzzzzz, g_0_x_yyzz_xzzzzzz, g_0_x_yyzz_yyyyyyy, g_0_x_yyzz_yyyyyyz, g_0_x_yyzz_yyyyyzz, g_0_x_yyzz_yyyyzzz, g_0_x_yyzz_yyyzzzz, g_0_x_yyzz_yyzzzzz, g_0_x_yyzz_yzzzzzz, g_0_x_yyzz_zzzzzzz, g_0_x_yzz_xxxxxxx, g_0_x_yzz_xxxxxxxy, g_0_x_yzz_xxxxxxy, g_0_x_yzz_xxxxxxyy, g_0_x_yzz_xxxxxxyz, g_0_x_yzz_xxxxxxz, g_0_x_yzz_xxxxxyy, g_0_x_yzz_xxxxxyyy, g_0_x_yzz_xxxxxyyz, g_0_x_yzz_xxxxxyz, g_0_x_yzz_xxxxxyzz, g_0_x_yzz_xxxxxzz, g_0_x_yzz_xxxxyyy, g_0_x_yzz_xxxxyyyy, g_0_x_yzz_xxxxyyyz, g_0_x_yzz_xxxxyyz, g_0_x_yzz_xxxxyyzz, g_0_x_yzz_xxxxyzz, g_0_x_yzz_xxxxyzzz, g_0_x_yzz_xxxxzzz, g_0_x_yzz_xxxyyyy, g_0_x_yzz_xxxyyyyy, g_0_x_yzz_xxxyyyyz, g_0_x_yzz_xxxyyyz, g_0_x_yzz_xxxyyyzz, g_0_x_yzz_xxxyyzz, g_0_x_yzz_xxxyyzzz, g_0_x_yzz_xxxyzzz, g_0_x_yzz_xxxyzzzz, g_0_x_yzz_xxxzzzz, g_0_x_yzz_xxyyyyy, g_0_x_yzz_xxyyyyyy, g_0_x_yzz_xxyyyyyz, g_0_x_yzz_xxyyyyz, g_0_x_yzz_xxyyyyzz, g_0_x_yzz_xxyyyzz, g_0_x_yzz_xxyyyzzz, g_0_x_yzz_xxyyzzz, g_0_x_yzz_xxyyzzzz, g_0_x_yzz_xxyzzzz, g_0_x_yzz_xxyzzzzz, g_0_x_yzz_xxzzzzz, g_0_x_yzz_xyyyyyy, g_0_x_yzz_xyyyyyyy, g_0_x_yzz_xyyyyyyz, g_0_x_yzz_xyyyyyz, g_0_x_yzz_xyyyyyzz, g_0_x_yzz_xyyyyzz, g_0_x_yzz_xyyyyzzz, g_0_x_yzz_xyyyzzz, g_0_x_yzz_xyyyzzzz, g_0_x_yzz_xyyzzzz, g_0_x_yzz_xyyzzzzz, g_0_x_yzz_xyzzzzz, g_0_x_yzz_xyzzzzzz, g_0_x_yzz_xzzzzzz, g_0_x_yzz_yyyyyyy, g_0_x_yzz_yyyyyyyy, g_0_x_yzz_yyyyyyyz, g_0_x_yzz_yyyyyyz, g_0_x_yzz_yyyyyyzz, g_0_x_yzz_yyyyyzz, g_0_x_yzz_yyyyyzzz, g_0_x_yzz_yyyyzzz, g_0_x_yzz_yyyyzzzz, g_0_x_yzz_yyyzzzz, g_0_x_yzz_yyyzzzzz, g_0_x_yzz_yyzzzzz, g_0_x_yzz_yyzzzzzz, g_0_x_yzz_yzzzzzz, g_0_x_yzz_yzzzzzzz, g_0_x_yzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_xxxxxxx[k] = -g_0_x_yzz_xxxxxxx[k] * ab_y + g_0_x_yzz_xxxxxxxy[k];

                g_0_x_yyzz_xxxxxxy[k] = -g_0_x_yzz_xxxxxxy[k] * ab_y + g_0_x_yzz_xxxxxxyy[k];

                g_0_x_yyzz_xxxxxxz[k] = -g_0_x_yzz_xxxxxxz[k] * ab_y + g_0_x_yzz_xxxxxxyz[k];

                g_0_x_yyzz_xxxxxyy[k] = -g_0_x_yzz_xxxxxyy[k] * ab_y + g_0_x_yzz_xxxxxyyy[k];

                g_0_x_yyzz_xxxxxyz[k] = -g_0_x_yzz_xxxxxyz[k] * ab_y + g_0_x_yzz_xxxxxyyz[k];

                g_0_x_yyzz_xxxxxzz[k] = -g_0_x_yzz_xxxxxzz[k] * ab_y + g_0_x_yzz_xxxxxyzz[k];

                g_0_x_yyzz_xxxxyyy[k] = -g_0_x_yzz_xxxxyyy[k] * ab_y + g_0_x_yzz_xxxxyyyy[k];

                g_0_x_yyzz_xxxxyyz[k] = -g_0_x_yzz_xxxxyyz[k] * ab_y + g_0_x_yzz_xxxxyyyz[k];

                g_0_x_yyzz_xxxxyzz[k] = -g_0_x_yzz_xxxxyzz[k] * ab_y + g_0_x_yzz_xxxxyyzz[k];

                g_0_x_yyzz_xxxxzzz[k] = -g_0_x_yzz_xxxxzzz[k] * ab_y + g_0_x_yzz_xxxxyzzz[k];

                g_0_x_yyzz_xxxyyyy[k] = -g_0_x_yzz_xxxyyyy[k] * ab_y + g_0_x_yzz_xxxyyyyy[k];

                g_0_x_yyzz_xxxyyyz[k] = -g_0_x_yzz_xxxyyyz[k] * ab_y + g_0_x_yzz_xxxyyyyz[k];

                g_0_x_yyzz_xxxyyzz[k] = -g_0_x_yzz_xxxyyzz[k] * ab_y + g_0_x_yzz_xxxyyyzz[k];

                g_0_x_yyzz_xxxyzzz[k] = -g_0_x_yzz_xxxyzzz[k] * ab_y + g_0_x_yzz_xxxyyzzz[k];

                g_0_x_yyzz_xxxzzzz[k] = -g_0_x_yzz_xxxzzzz[k] * ab_y + g_0_x_yzz_xxxyzzzz[k];

                g_0_x_yyzz_xxyyyyy[k] = -g_0_x_yzz_xxyyyyy[k] * ab_y + g_0_x_yzz_xxyyyyyy[k];

                g_0_x_yyzz_xxyyyyz[k] = -g_0_x_yzz_xxyyyyz[k] * ab_y + g_0_x_yzz_xxyyyyyz[k];

                g_0_x_yyzz_xxyyyzz[k] = -g_0_x_yzz_xxyyyzz[k] * ab_y + g_0_x_yzz_xxyyyyzz[k];

                g_0_x_yyzz_xxyyzzz[k] = -g_0_x_yzz_xxyyzzz[k] * ab_y + g_0_x_yzz_xxyyyzzz[k];

                g_0_x_yyzz_xxyzzzz[k] = -g_0_x_yzz_xxyzzzz[k] * ab_y + g_0_x_yzz_xxyyzzzz[k];

                g_0_x_yyzz_xxzzzzz[k] = -g_0_x_yzz_xxzzzzz[k] * ab_y + g_0_x_yzz_xxyzzzzz[k];

                g_0_x_yyzz_xyyyyyy[k] = -g_0_x_yzz_xyyyyyy[k] * ab_y + g_0_x_yzz_xyyyyyyy[k];

                g_0_x_yyzz_xyyyyyz[k] = -g_0_x_yzz_xyyyyyz[k] * ab_y + g_0_x_yzz_xyyyyyyz[k];

                g_0_x_yyzz_xyyyyzz[k] = -g_0_x_yzz_xyyyyzz[k] * ab_y + g_0_x_yzz_xyyyyyzz[k];

                g_0_x_yyzz_xyyyzzz[k] = -g_0_x_yzz_xyyyzzz[k] * ab_y + g_0_x_yzz_xyyyyzzz[k];

                g_0_x_yyzz_xyyzzzz[k] = -g_0_x_yzz_xyyzzzz[k] * ab_y + g_0_x_yzz_xyyyzzzz[k];

                g_0_x_yyzz_xyzzzzz[k] = -g_0_x_yzz_xyzzzzz[k] * ab_y + g_0_x_yzz_xyyzzzzz[k];

                g_0_x_yyzz_xzzzzzz[k] = -g_0_x_yzz_xzzzzzz[k] * ab_y + g_0_x_yzz_xyzzzzzz[k];

                g_0_x_yyzz_yyyyyyy[k] = -g_0_x_yzz_yyyyyyy[k] * ab_y + g_0_x_yzz_yyyyyyyy[k];

                g_0_x_yyzz_yyyyyyz[k] = -g_0_x_yzz_yyyyyyz[k] * ab_y + g_0_x_yzz_yyyyyyyz[k];

                g_0_x_yyzz_yyyyyzz[k] = -g_0_x_yzz_yyyyyzz[k] * ab_y + g_0_x_yzz_yyyyyyzz[k];

                g_0_x_yyzz_yyyyzzz[k] = -g_0_x_yzz_yyyyzzz[k] * ab_y + g_0_x_yzz_yyyyyzzz[k];

                g_0_x_yyzz_yyyzzzz[k] = -g_0_x_yzz_yyyzzzz[k] * ab_y + g_0_x_yzz_yyyyzzzz[k];

                g_0_x_yyzz_yyzzzzz[k] = -g_0_x_yzz_yyzzzzz[k] * ab_y + g_0_x_yzz_yyyzzzzz[k];

                g_0_x_yyzz_yzzzzzz[k] = -g_0_x_yzz_yzzzzzz[k] * ab_y + g_0_x_yzz_yyzzzzzz[k];

                g_0_x_yyzz_zzzzzzz[k] = -g_0_x_yzz_zzzzzzz[k] * ab_y + g_0_x_yzz_yzzzzzzz[k];
            }

            /// Set up 468-504 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yzzz_xxxxxxx, g_0_x_yzzz_xxxxxxy, g_0_x_yzzz_xxxxxxz, g_0_x_yzzz_xxxxxyy, g_0_x_yzzz_xxxxxyz, g_0_x_yzzz_xxxxxzz, g_0_x_yzzz_xxxxyyy, g_0_x_yzzz_xxxxyyz, g_0_x_yzzz_xxxxyzz, g_0_x_yzzz_xxxxzzz, g_0_x_yzzz_xxxyyyy, g_0_x_yzzz_xxxyyyz, g_0_x_yzzz_xxxyyzz, g_0_x_yzzz_xxxyzzz, g_0_x_yzzz_xxxzzzz, g_0_x_yzzz_xxyyyyy, g_0_x_yzzz_xxyyyyz, g_0_x_yzzz_xxyyyzz, g_0_x_yzzz_xxyyzzz, g_0_x_yzzz_xxyzzzz, g_0_x_yzzz_xxzzzzz, g_0_x_yzzz_xyyyyyy, g_0_x_yzzz_xyyyyyz, g_0_x_yzzz_xyyyyzz, g_0_x_yzzz_xyyyzzz, g_0_x_yzzz_xyyzzzz, g_0_x_yzzz_xyzzzzz, g_0_x_yzzz_xzzzzzz, g_0_x_yzzz_yyyyyyy, g_0_x_yzzz_yyyyyyz, g_0_x_yzzz_yyyyyzz, g_0_x_yzzz_yyyyzzz, g_0_x_yzzz_yyyzzzz, g_0_x_yzzz_yyzzzzz, g_0_x_yzzz_yzzzzzz, g_0_x_yzzz_zzzzzzz, g_0_x_zzz_xxxxxxx, g_0_x_zzz_xxxxxxxy, g_0_x_zzz_xxxxxxy, g_0_x_zzz_xxxxxxyy, g_0_x_zzz_xxxxxxyz, g_0_x_zzz_xxxxxxz, g_0_x_zzz_xxxxxyy, g_0_x_zzz_xxxxxyyy, g_0_x_zzz_xxxxxyyz, g_0_x_zzz_xxxxxyz, g_0_x_zzz_xxxxxyzz, g_0_x_zzz_xxxxxzz, g_0_x_zzz_xxxxyyy, g_0_x_zzz_xxxxyyyy, g_0_x_zzz_xxxxyyyz, g_0_x_zzz_xxxxyyz, g_0_x_zzz_xxxxyyzz, g_0_x_zzz_xxxxyzz, g_0_x_zzz_xxxxyzzz, g_0_x_zzz_xxxxzzz, g_0_x_zzz_xxxyyyy, g_0_x_zzz_xxxyyyyy, g_0_x_zzz_xxxyyyyz, g_0_x_zzz_xxxyyyz, g_0_x_zzz_xxxyyyzz, g_0_x_zzz_xxxyyzz, g_0_x_zzz_xxxyyzzz, g_0_x_zzz_xxxyzzz, g_0_x_zzz_xxxyzzzz, g_0_x_zzz_xxxzzzz, g_0_x_zzz_xxyyyyy, g_0_x_zzz_xxyyyyyy, g_0_x_zzz_xxyyyyyz, g_0_x_zzz_xxyyyyz, g_0_x_zzz_xxyyyyzz, g_0_x_zzz_xxyyyzz, g_0_x_zzz_xxyyyzzz, g_0_x_zzz_xxyyzzz, g_0_x_zzz_xxyyzzzz, g_0_x_zzz_xxyzzzz, g_0_x_zzz_xxyzzzzz, g_0_x_zzz_xxzzzzz, g_0_x_zzz_xyyyyyy, g_0_x_zzz_xyyyyyyy, g_0_x_zzz_xyyyyyyz, g_0_x_zzz_xyyyyyz, g_0_x_zzz_xyyyyyzz, g_0_x_zzz_xyyyyzz, g_0_x_zzz_xyyyyzzz, g_0_x_zzz_xyyyzzz, g_0_x_zzz_xyyyzzzz, g_0_x_zzz_xyyzzzz, g_0_x_zzz_xyyzzzzz, g_0_x_zzz_xyzzzzz, g_0_x_zzz_xyzzzzzz, g_0_x_zzz_xzzzzzz, g_0_x_zzz_yyyyyyy, g_0_x_zzz_yyyyyyyy, g_0_x_zzz_yyyyyyyz, g_0_x_zzz_yyyyyyz, g_0_x_zzz_yyyyyyzz, g_0_x_zzz_yyyyyzz, g_0_x_zzz_yyyyyzzz, g_0_x_zzz_yyyyzzz, g_0_x_zzz_yyyyzzzz, g_0_x_zzz_yyyzzzz, g_0_x_zzz_yyyzzzzz, g_0_x_zzz_yyzzzzz, g_0_x_zzz_yyzzzzzz, g_0_x_zzz_yzzzzzz, g_0_x_zzz_yzzzzzzz, g_0_x_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_xxxxxxx[k] = -g_0_x_zzz_xxxxxxx[k] * ab_y + g_0_x_zzz_xxxxxxxy[k];

                g_0_x_yzzz_xxxxxxy[k] = -g_0_x_zzz_xxxxxxy[k] * ab_y + g_0_x_zzz_xxxxxxyy[k];

                g_0_x_yzzz_xxxxxxz[k] = -g_0_x_zzz_xxxxxxz[k] * ab_y + g_0_x_zzz_xxxxxxyz[k];

                g_0_x_yzzz_xxxxxyy[k] = -g_0_x_zzz_xxxxxyy[k] * ab_y + g_0_x_zzz_xxxxxyyy[k];

                g_0_x_yzzz_xxxxxyz[k] = -g_0_x_zzz_xxxxxyz[k] * ab_y + g_0_x_zzz_xxxxxyyz[k];

                g_0_x_yzzz_xxxxxzz[k] = -g_0_x_zzz_xxxxxzz[k] * ab_y + g_0_x_zzz_xxxxxyzz[k];

                g_0_x_yzzz_xxxxyyy[k] = -g_0_x_zzz_xxxxyyy[k] * ab_y + g_0_x_zzz_xxxxyyyy[k];

                g_0_x_yzzz_xxxxyyz[k] = -g_0_x_zzz_xxxxyyz[k] * ab_y + g_0_x_zzz_xxxxyyyz[k];

                g_0_x_yzzz_xxxxyzz[k] = -g_0_x_zzz_xxxxyzz[k] * ab_y + g_0_x_zzz_xxxxyyzz[k];

                g_0_x_yzzz_xxxxzzz[k] = -g_0_x_zzz_xxxxzzz[k] * ab_y + g_0_x_zzz_xxxxyzzz[k];

                g_0_x_yzzz_xxxyyyy[k] = -g_0_x_zzz_xxxyyyy[k] * ab_y + g_0_x_zzz_xxxyyyyy[k];

                g_0_x_yzzz_xxxyyyz[k] = -g_0_x_zzz_xxxyyyz[k] * ab_y + g_0_x_zzz_xxxyyyyz[k];

                g_0_x_yzzz_xxxyyzz[k] = -g_0_x_zzz_xxxyyzz[k] * ab_y + g_0_x_zzz_xxxyyyzz[k];

                g_0_x_yzzz_xxxyzzz[k] = -g_0_x_zzz_xxxyzzz[k] * ab_y + g_0_x_zzz_xxxyyzzz[k];

                g_0_x_yzzz_xxxzzzz[k] = -g_0_x_zzz_xxxzzzz[k] * ab_y + g_0_x_zzz_xxxyzzzz[k];

                g_0_x_yzzz_xxyyyyy[k] = -g_0_x_zzz_xxyyyyy[k] * ab_y + g_0_x_zzz_xxyyyyyy[k];

                g_0_x_yzzz_xxyyyyz[k] = -g_0_x_zzz_xxyyyyz[k] * ab_y + g_0_x_zzz_xxyyyyyz[k];

                g_0_x_yzzz_xxyyyzz[k] = -g_0_x_zzz_xxyyyzz[k] * ab_y + g_0_x_zzz_xxyyyyzz[k];

                g_0_x_yzzz_xxyyzzz[k] = -g_0_x_zzz_xxyyzzz[k] * ab_y + g_0_x_zzz_xxyyyzzz[k];

                g_0_x_yzzz_xxyzzzz[k] = -g_0_x_zzz_xxyzzzz[k] * ab_y + g_0_x_zzz_xxyyzzzz[k];

                g_0_x_yzzz_xxzzzzz[k] = -g_0_x_zzz_xxzzzzz[k] * ab_y + g_0_x_zzz_xxyzzzzz[k];

                g_0_x_yzzz_xyyyyyy[k] = -g_0_x_zzz_xyyyyyy[k] * ab_y + g_0_x_zzz_xyyyyyyy[k];

                g_0_x_yzzz_xyyyyyz[k] = -g_0_x_zzz_xyyyyyz[k] * ab_y + g_0_x_zzz_xyyyyyyz[k];

                g_0_x_yzzz_xyyyyzz[k] = -g_0_x_zzz_xyyyyzz[k] * ab_y + g_0_x_zzz_xyyyyyzz[k];

                g_0_x_yzzz_xyyyzzz[k] = -g_0_x_zzz_xyyyzzz[k] * ab_y + g_0_x_zzz_xyyyyzzz[k];

                g_0_x_yzzz_xyyzzzz[k] = -g_0_x_zzz_xyyzzzz[k] * ab_y + g_0_x_zzz_xyyyzzzz[k];

                g_0_x_yzzz_xyzzzzz[k] = -g_0_x_zzz_xyzzzzz[k] * ab_y + g_0_x_zzz_xyyzzzzz[k];

                g_0_x_yzzz_xzzzzzz[k] = -g_0_x_zzz_xzzzzzz[k] * ab_y + g_0_x_zzz_xyzzzzzz[k];

                g_0_x_yzzz_yyyyyyy[k] = -g_0_x_zzz_yyyyyyy[k] * ab_y + g_0_x_zzz_yyyyyyyy[k];

                g_0_x_yzzz_yyyyyyz[k] = -g_0_x_zzz_yyyyyyz[k] * ab_y + g_0_x_zzz_yyyyyyyz[k];

                g_0_x_yzzz_yyyyyzz[k] = -g_0_x_zzz_yyyyyzz[k] * ab_y + g_0_x_zzz_yyyyyyzz[k];

                g_0_x_yzzz_yyyyzzz[k] = -g_0_x_zzz_yyyyzzz[k] * ab_y + g_0_x_zzz_yyyyyzzz[k];

                g_0_x_yzzz_yyyzzzz[k] = -g_0_x_zzz_yyyzzzz[k] * ab_y + g_0_x_zzz_yyyyzzzz[k];

                g_0_x_yzzz_yyzzzzz[k] = -g_0_x_zzz_yyzzzzz[k] * ab_y + g_0_x_zzz_yyyzzzzz[k];

                g_0_x_yzzz_yzzzzzz[k] = -g_0_x_zzz_yzzzzzz[k] * ab_y + g_0_x_zzz_yyzzzzzz[k];

                g_0_x_yzzz_zzzzzzz[k] = -g_0_x_zzz_zzzzzzz[k] * ab_y + g_0_x_zzz_yzzzzzzz[k];
            }

            /// Set up 504-540 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_zzz_xxxxxxx, g_0_x_zzz_xxxxxxxz, g_0_x_zzz_xxxxxxy, g_0_x_zzz_xxxxxxyz, g_0_x_zzz_xxxxxxz, g_0_x_zzz_xxxxxxzz, g_0_x_zzz_xxxxxyy, g_0_x_zzz_xxxxxyyz, g_0_x_zzz_xxxxxyz, g_0_x_zzz_xxxxxyzz, g_0_x_zzz_xxxxxzz, g_0_x_zzz_xxxxxzzz, g_0_x_zzz_xxxxyyy, g_0_x_zzz_xxxxyyyz, g_0_x_zzz_xxxxyyz, g_0_x_zzz_xxxxyyzz, g_0_x_zzz_xxxxyzz, g_0_x_zzz_xxxxyzzz, g_0_x_zzz_xxxxzzz, g_0_x_zzz_xxxxzzzz, g_0_x_zzz_xxxyyyy, g_0_x_zzz_xxxyyyyz, g_0_x_zzz_xxxyyyz, g_0_x_zzz_xxxyyyzz, g_0_x_zzz_xxxyyzz, g_0_x_zzz_xxxyyzzz, g_0_x_zzz_xxxyzzz, g_0_x_zzz_xxxyzzzz, g_0_x_zzz_xxxzzzz, g_0_x_zzz_xxxzzzzz, g_0_x_zzz_xxyyyyy, g_0_x_zzz_xxyyyyyz, g_0_x_zzz_xxyyyyz, g_0_x_zzz_xxyyyyzz, g_0_x_zzz_xxyyyzz, g_0_x_zzz_xxyyyzzz, g_0_x_zzz_xxyyzzz, g_0_x_zzz_xxyyzzzz, g_0_x_zzz_xxyzzzz, g_0_x_zzz_xxyzzzzz, g_0_x_zzz_xxzzzzz, g_0_x_zzz_xxzzzzzz, g_0_x_zzz_xyyyyyy, g_0_x_zzz_xyyyyyyz, g_0_x_zzz_xyyyyyz, g_0_x_zzz_xyyyyyzz, g_0_x_zzz_xyyyyzz, g_0_x_zzz_xyyyyzzz, g_0_x_zzz_xyyyzzz, g_0_x_zzz_xyyyzzzz, g_0_x_zzz_xyyzzzz, g_0_x_zzz_xyyzzzzz, g_0_x_zzz_xyzzzzz, g_0_x_zzz_xyzzzzzz, g_0_x_zzz_xzzzzzz, g_0_x_zzz_xzzzzzzz, g_0_x_zzz_yyyyyyy, g_0_x_zzz_yyyyyyyz, g_0_x_zzz_yyyyyyz, g_0_x_zzz_yyyyyyzz, g_0_x_zzz_yyyyyzz, g_0_x_zzz_yyyyyzzz, g_0_x_zzz_yyyyzzz, g_0_x_zzz_yyyyzzzz, g_0_x_zzz_yyyzzzz, g_0_x_zzz_yyyzzzzz, g_0_x_zzz_yyzzzzz, g_0_x_zzz_yyzzzzzz, g_0_x_zzz_yzzzzzz, g_0_x_zzz_yzzzzzzz, g_0_x_zzz_zzzzzzz, g_0_x_zzz_zzzzzzzz, g_0_x_zzzz_xxxxxxx, g_0_x_zzzz_xxxxxxy, g_0_x_zzzz_xxxxxxz, g_0_x_zzzz_xxxxxyy, g_0_x_zzzz_xxxxxyz, g_0_x_zzzz_xxxxxzz, g_0_x_zzzz_xxxxyyy, g_0_x_zzzz_xxxxyyz, g_0_x_zzzz_xxxxyzz, g_0_x_zzzz_xxxxzzz, g_0_x_zzzz_xxxyyyy, g_0_x_zzzz_xxxyyyz, g_0_x_zzzz_xxxyyzz, g_0_x_zzzz_xxxyzzz, g_0_x_zzzz_xxxzzzz, g_0_x_zzzz_xxyyyyy, g_0_x_zzzz_xxyyyyz, g_0_x_zzzz_xxyyyzz, g_0_x_zzzz_xxyyzzz, g_0_x_zzzz_xxyzzzz, g_0_x_zzzz_xxzzzzz, g_0_x_zzzz_xyyyyyy, g_0_x_zzzz_xyyyyyz, g_0_x_zzzz_xyyyyzz, g_0_x_zzzz_xyyyzzz, g_0_x_zzzz_xyyzzzz, g_0_x_zzzz_xyzzzzz, g_0_x_zzzz_xzzzzzz, g_0_x_zzzz_yyyyyyy, g_0_x_zzzz_yyyyyyz, g_0_x_zzzz_yyyyyzz, g_0_x_zzzz_yyyyzzz, g_0_x_zzzz_yyyzzzz, g_0_x_zzzz_yyzzzzz, g_0_x_zzzz_yzzzzzz, g_0_x_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_xxxxxxx[k] = -g_0_x_zzz_xxxxxxx[k] * ab_z + g_0_x_zzz_xxxxxxxz[k];

                g_0_x_zzzz_xxxxxxy[k] = -g_0_x_zzz_xxxxxxy[k] * ab_z + g_0_x_zzz_xxxxxxyz[k];

                g_0_x_zzzz_xxxxxxz[k] = -g_0_x_zzz_xxxxxxz[k] * ab_z + g_0_x_zzz_xxxxxxzz[k];

                g_0_x_zzzz_xxxxxyy[k] = -g_0_x_zzz_xxxxxyy[k] * ab_z + g_0_x_zzz_xxxxxyyz[k];

                g_0_x_zzzz_xxxxxyz[k] = -g_0_x_zzz_xxxxxyz[k] * ab_z + g_0_x_zzz_xxxxxyzz[k];

                g_0_x_zzzz_xxxxxzz[k] = -g_0_x_zzz_xxxxxzz[k] * ab_z + g_0_x_zzz_xxxxxzzz[k];

                g_0_x_zzzz_xxxxyyy[k] = -g_0_x_zzz_xxxxyyy[k] * ab_z + g_0_x_zzz_xxxxyyyz[k];

                g_0_x_zzzz_xxxxyyz[k] = -g_0_x_zzz_xxxxyyz[k] * ab_z + g_0_x_zzz_xxxxyyzz[k];

                g_0_x_zzzz_xxxxyzz[k] = -g_0_x_zzz_xxxxyzz[k] * ab_z + g_0_x_zzz_xxxxyzzz[k];

                g_0_x_zzzz_xxxxzzz[k] = -g_0_x_zzz_xxxxzzz[k] * ab_z + g_0_x_zzz_xxxxzzzz[k];

                g_0_x_zzzz_xxxyyyy[k] = -g_0_x_zzz_xxxyyyy[k] * ab_z + g_0_x_zzz_xxxyyyyz[k];

                g_0_x_zzzz_xxxyyyz[k] = -g_0_x_zzz_xxxyyyz[k] * ab_z + g_0_x_zzz_xxxyyyzz[k];

                g_0_x_zzzz_xxxyyzz[k] = -g_0_x_zzz_xxxyyzz[k] * ab_z + g_0_x_zzz_xxxyyzzz[k];

                g_0_x_zzzz_xxxyzzz[k] = -g_0_x_zzz_xxxyzzz[k] * ab_z + g_0_x_zzz_xxxyzzzz[k];

                g_0_x_zzzz_xxxzzzz[k] = -g_0_x_zzz_xxxzzzz[k] * ab_z + g_0_x_zzz_xxxzzzzz[k];

                g_0_x_zzzz_xxyyyyy[k] = -g_0_x_zzz_xxyyyyy[k] * ab_z + g_0_x_zzz_xxyyyyyz[k];

                g_0_x_zzzz_xxyyyyz[k] = -g_0_x_zzz_xxyyyyz[k] * ab_z + g_0_x_zzz_xxyyyyzz[k];

                g_0_x_zzzz_xxyyyzz[k] = -g_0_x_zzz_xxyyyzz[k] * ab_z + g_0_x_zzz_xxyyyzzz[k];

                g_0_x_zzzz_xxyyzzz[k] = -g_0_x_zzz_xxyyzzz[k] * ab_z + g_0_x_zzz_xxyyzzzz[k];

                g_0_x_zzzz_xxyzzzz[k] = -g_0_x_zzz_xxyzzzz[k] * ab_z + g_0_x_zzz_xxyzzzzz[k];

                g_0_x_zzzz_xxzzzzz[k] = -g_0_x_zzz_xxzzzzz[k] * ab_z + g_0_x_zzz_xxzzzzzz[k];

                g_0_x_zzzz_xyyyyyy[k] = -g_0_x_zzz_xyyyyyy[k] * ab_z + g_0_x_zzz_xyyyyyyz[k];

                g_0_x_zzzz_xyyyyyz[k] = -g_0_x_zzz_xyyyyyz[k] * ab_z + g_0_x_zzz_xyyyyyzz[k];

                g_0_x_zzzz_xyyyyzz[k] = -g_0_x_zzz_xyyyyzz[k] * ab_z + g_0_x_zzz_xyyyyzzz[k];

                g_0_x_zzzz_xyyyzzz[k] = -g_0_x_zzz_xyyyzzz[k] * ab_z + g_0_x_zzz_xyyyzzzz[k];

                g_0_x_zzzz_xyyzzzz[k] = -g_0_x_zzz_xyyzzzz[k] * ab_z + g_0_x_zzz_xyyzzzzz[k];

                g_0_x_zzzz_xyzzzzz[k] = -g_0_x_zzz_xyzzzzz[k] * ab_z + g_0_x_zzz_xyzzzzzz[k];

                g_0_x_zzzz_xzzzzzz[k] = -g_0_x_zzz_xzzzzzz[k] * ab_z + g_0_x_zzz_xzzzzzzz[k];

                g_0_x_zzzz_yyyyyyy[k] = -g_0_x_zzz_yyyyyyy[k] * ab_z + g_0_x_zzz_yyyyyyyz[k];

                g_0_x_zzzz_yyyyyyz[k] = -g_0_x_zzz_yyyyyyz[k] * ab_z + g_0_x_zzz_yyyyyyzz[k];

                g_0_x_zzzz_yyyyyzz[k] = -g_0_x_zzz_yyyyyzz[k] * ab_z + g_0_x_zzz_yyyyyzzz[k];

                g_0_x_zzzz_yyyyzzz[k] = -g_0_x_zzz_yyyyzzz[k] * ab_z + g_0_x_zzz_yyyyzzzz[k];

                g_0_x_zzzz_yyyzzzz[k] = -g_0_x_zzz_yyyzzzz[k] * ab_z + g_0_x_zzz_yyyzzzzz[k];

                g_0_x_zzzz_yyzzzzz[k] = -g_0_x_zzz_yyzzzzz[k] * ab_z + g_0_x_zzz_yyzzzzzz[k];

                g_0_x_zzzz_yzzzzzz[k] = -g_0_x_zzz_yzzzzzz[k] * ab_z + g_0_x_zzz_yzzzzzzz[k];

                g_0_x_zzzz_zzzzzzz[k] = -g_0_x_zzz_zzzzzzz[k] * ab_z + g_0_x_zzz_zzzzzzzz[k];
            }

            /// Set up 540-576 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxx_xxxxxxx, g_0_y_xxx_xxxxxxxx, g_0_y_xxx_xxxxxxxy, g_0_y_xxx_xxxxxxxz, g_0_y_xxx_xxxxxxy, g_0_y_xxx_xxxxxxyy, g_0_y_xxx_xxxxxxyz, g_0_y_xxx_xxxxxxz, g_0_y_xxx_xxxxxxzz, g_0_y_xxx_xxxxxyy, g_0_y_xxx_xxxxxyyy, g_0_y_xxx_xxxxxyyz, g_0_y_xxx_xxxxxyz, g_0_y_xxx_xxxxxyzz, g_0_y_xxx_xxxxxzz, g_0_y_xxx_xxxxxzzz, g_0_y_xxx_xxxxyyy, g_0_y_xxx_xxxxyyyy, g_0_y_xxx_xxxxyyyz, g_0_y_xxx_xxxxyyz, g_0_y_xxx_xxxxyyzz, g_0_y_xxx_xxxxyzz, g_0_y_xxx_xxxxyzzz, g_0_y_xxx_xxxxzzz, g_0_y_xxx_xxxxzzzz, g_0_y_xxx_xxxyyyy, g_0_y_xxx_xxxyyyyy, g_0_y_xxx_xxxyyyyz, g_0_y_xxx_xxxyyyz, g_0_y_xxx_xxxyyyzz, g_0_y_xxx_xxxyyzz, g_0_y_xxx_xxxyyzzz, g_0_y_xxx_xxxyzzz, g_0_y_xxx_xxxyzzzz, g_0_y_xxx_xxxzzzz, g_0_y_xxx_xxxzzzzz, g_0_y_xxx_xxyyyyy, g_0_y_xxx_xxyyyyyy, g_0_y_xxx_xxyyyyyz, g_0_y_xxx_xxyyyyz, g_0_y_xxx_xxyyyyzz, g_0_y_xxx_xxyyyzz, g_0_y_xxx_xxyyyzzz, g_0_y_xxx_xxyyzzz, g_0_y_xxx_xxyyzzzz, g_0_y_xxx_xxyzzzz, g_0_y_xxx_xxyzzzzz, g_0_y_xxx_xxzzzzz, g_0_y_xxx_xxzzzzzz, g_0_y_xxx_xyyyyyy, g_0_y_xxx_xyyyyyyy, g_0_y_xxx_xyyyyyyz, g_0_y_xxx_xyyyyyz, g_0_y_xxx_xyyyyyzz, g_0_y_xxx_xyyyyzz, g_0_y_xxx_xyyyyzzz, g_0_y_xxx_xyyyzzz, g_0_y_xxx_xyyyzzzz, g_0_y_xxx_xyyzzzz, g_0_y_xxx_xyyzzzzz, g_0_y_xxx_xyzzzzz, g_0_y_xxx_xyzzzzzz, g_0_y_xxx_xzzzzzz, g_0_y_xxx_xzzzzzzz, g_0_y_xxx_yyyyyyy, g_0_y_xxx_yyyyyyz, g_0_y_xxx_yyyyyzz, g_0_y_xxx_yyyyzzz, g_0_y_xxx_yyyzzzz, g_0_y_xxx_yyzzzzz, g_0_y_xxx_yzzzzzz, g_0_y_xxx_zzzzzzz, g_0_y_xxxx_xxxxxxx, g_0_y_xxxx_xxxxxxy, g_0_y_xxxx_xxxxxxz, g_0_y_xxxx_xxxxxyy, g_0_y_xxxx_xxxxxyz, g_0_y_xxxx_xxxxxzz, g_0_y_xxxx_xxxxyyy, g_0_y_xxxx_xxxxyyz, g_0_y_xxxx_xxxxyzz, g_0_y_xxxx_xxxxzzz, g_0_y_xxxx_xxxyyyy, g_0_y_xxxx_xxxyyyz, g_0_y_xxxx_xxxyyzz, g_0_y_xxxx_xxxyzzz, g_0_y_xxxx_xxxzzzz, g_0_y_xxxx_xxyyyyy, g_0_y_xxxx_xxyyyyz, g_0_y_xxxx_xxyyyzz, g_0_y_xxxx_xxyyzzz, g_0_y_xxxx_xxyzzzz, g_0_y_xxxx_xxzzzzz, g_0_y_xxxx_xyyyyyy, g_0_y_xxxx_xyyyyyz, g_0_y_xxxx_xyyyyzz, g_0_y_xxxx_xyyyzzz, g_0_y_xxxx_xyyzzzz, g_0_y_xxxx_xyzzzzz, g_0_y_xxxx_xzzzzzz, g_0_y_xxxx_yyyyyyy, g_0_y_xxxx_yyyyyyz, g_0_y_xxxx_yyyyyzz, g_0_y_xxxx_yyyyzzz, g_0_y_xxxx_yyyzzzz, g_0_y_xxxx_yyzzzzz, g_0_y_xxxx_yzzzzzz, g_0_y_xxxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_xxxxxxx[k] = -g_0_y_xxx_xxxxxxx[k] * ab_x + g_0_y_xxx_xxxxxxxx[k];

                g_0_y_xxxx_xxxxxxy[k] = -g_0_y_xxx_xxxxxxy[k] * ab_x + g_0_y_xxx_xxxxxxxy[k];

                g_0_y_xxxx_xxxxxxz[k] = -g_0_y_xxx_xxxxxxz[k] * ab_x + g_0_y_xxx_xxxxxxxz[k];

                g_0_y_xxxx_xxxxxyy[k] = -g_0_y_xxx_xxxxxyy[k] * ab_x + g_0_y_xxx_xxxxxxyy[k];

                g_0_y_xxxx_xxxxxyz[k] = -g_0_y_xxx_xxxxxyz[k] * ab_x + g_0_y_xxx_xxxxxxyz[k];

                g_0_y_xxxx_xxxxxzz[k] = -g_0_y_xxx_xxxxxzz[k] * ab_x + g_0_y_xxx_xxxxxxzz[k];

                g_0_y_xxxx_xxxxyyy[k] = -g_0_y_xxx_xxxxyyy[k] * ab_x + g_0_y_xxx_xxxxxyyy[k];

                g_0_y_xxxx_xxxxyyz[k] = -g_0_y_xxx_xxxxyyz[k] * ab_x + g_0_y_xxx_xxxxxyyz[k];

                g_0_y_xxxx_xxxxyzz[k] = -g_0_y_xxx_xxxxyzz[k] * ab_x + g_0_y_xxx_xxxxxyzz[k];

                g_0_y_xxxx_xxxxzzz[k] = -g_0_y_xxx_xxxxzzz[k] * ab_x + g_0_y_xxx_xxxxxzzz[k];

                g_0_y_xxxx_xxxyyyy[k] = -g_0_y_xxx_xxxyyyy[k] * ab_x + g_0_y_xxx_xxxxyyyy[k];

                g_0_y_xxxx_xxxyyyz[k] = -g_0_y_xxx_xxxyyyz[k] * ab_x + g_0_y_xxx_xxxxyyyz[k];

                g_0_y_xxxx_xxxyyzz[k] = -g_0_y_xxx_xxxyyzz[k] * ab_x + g_0_y_xxx_xxxxyyzz[k];

                g_0_y_xxxx_xxxyzzz[k] = -g_0_y_xxx_xxxyzzz[k] * ab_x + g_0_y_xxx_xxxxyzzz[k];

                g_0_y_xxxx_xxxzzzz[k] = -g_0_y_xxx_xxxzzzz[k] * ab_x + g_0_y_xxx_xxxxzzzz[k];

                g_0_y_xxxx_xxyyyyy[k] = -g_0_y_xxx_xxyyyyy[k] * ab_x + g_0_y_xxx_xxxyyyyy[k];

                g_0_y_xxxx_xxyyyyz[k] = -g_0_y_xxx_xxyyyyz[k] * ab_x + g_0_y_xxx_xxxyyyyz[k];

                g_0_y_xxxx_xxyyyzz[k] = -g_0_y_xxx_xxyyyzz[k] * ab_x + g_0_y_xxx_xxxyyyzz[k];

                g_0_y_xxxx_xxyyzzz[k] = -g_0_y_xxx_xxyyzzz[k] * ab_x + g_0_y_xxx_xxxyyzzz[k];

                g_0_y_xxxx_xxyzzzz[k] = -g_0_y_xxx_xxyzzzz[k] * ab_x + g_0_y_xxx_xxxyzzzz[k];

                g_0_y_xxxx_xxzzzzz[k] = -g_0_y_xxx_xxzzzzz[k] * ab_x + g_0_y_xxx_xxxzzzzz[k];

                g_0_y_xxxx_xyyyyyy[k] = -g_0_y_xxx_xyyyyyy[k] * ab_x + g_0_y_xxx_xxyyyyyy[k];

                g_0_y_xxxx_xyyyyyz[k] = -g_0_y_xxx_xyyyyyz[k] * ab_x + g_0_y_xxx_xxyyyyyz[k];

                g_0_y_xxxx_xyyyyzz[k] = -g_0_y_xxx_xyyyyzz[k] * ab_x + g_0_y_xxx_xxyyyyzz[k];

                g_0_y_xxxx_xyyyzzz[k] = -g_0_y_xxx_xyyyzzz[k] * ab_x + g_0_y_xxx_xxyyyzzz[k];

                g_0_y_xxxx_xyyzzzz[k] = -g_0_y_xxx_xyyzzzz[k] * ab_x + g_0_y_xxx_xxyyzzzz[k];

                g_0_y_xxxx_xyzzzzz[k] = -g_0_y_xxx_xyzzzzz[k] * ab_x + g_0_y_xxx_xxyzzzzz[k];

                g_0_y_xxxx_xzzzzzz[k] = -g_0_y_xxx_xzzzzzz[k] * ab_x + g_0_y_xxx_xxzzzzzz[k];

                g_0_y_xxxx_yyyyyyy[k] = -g_0_y_xxx_yyyyyyy[k] * ab_x + g_0_y_xxx_xyyyyyyy[k];

                g_0_y_xxxx_yyyyyyz[k] = -g_0_y_xxx_yyyyyyz[k] * ab_x + g_0_y_xxx_xyyyyyyz[k];

                g_0_y_xxxx_yyyyyzz[k] = -g_0_y_xxx_yyyyyzz[k] * ab_x + g_0_y_xxx_xyyyyyzz[k];

                g_0_y_xxxx_yyyyzzz[k] = -g_0_y_xxx_yyyyzzz[k] * ab_x + g_0_y_xxx_xyyyyzzz[k];

                g_0_y_xxxx_yyyzzzz[k] = -g_0_y_xxx_yyyzzzz[k] * ab_x + g_0_y_xxx_xyyyzzzz[k];

                g_0_y_xxxx_yyzzzzz[k] = -g_0_y_xxx_yyzzzzz[k] * ab_x + g_0_y_xxx_xyyzzzzz[k];

                g_0_y_xxxx_yzzzzzz[k] = -g_0_y_xxx_yzzzzzz[k] * ab_x + g_0_y_xxx_xyzzzzzz[k];

                g_0_y_xxxx_zzzzzzz[k] = -g_0_y_xxx_zzzzzzz[k] * ab_x + g_0_y_xxx_xzzzzzzz[k];
            }

            /// Set up 576-612 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxxy_xxxxxxx, g_0_y_xxxy_xxxxxxy, g_0_y_xxxy_xxxxxxz, g_0_y_xxxy_xxxxxyy, g_0_y_xxxy_xxxxxyz, g_0_y_xxxy_xxxxxzz, g_0_y_xxxy_xxxxyyy, g_0_y_xxxy_xxxxyyz, g_0_y_xxxy_xxxxyzz, g_0_y_xxxy_xxxxzzz, g_0_y_xxxy_xxxyyyy, g_0_y_xxxy_xxxyyyz, g_0_y_xxxy_xxxyyzz, g_0_y_xxxy_xxxyzzz, g_0_y_xxxy_xxxzzzz, g_0_y_xxxy_xxyyyyy, g_0_y_xxxy_xxyyyyz, g_0_y_xxxy_xxyyyzz, g_0_y_xxxy_xxyyzzz, g_0_y_xxxy_xxyzzzz, g_0_y_xxxy_xxzzzzz, g_0_y_xxxy_xyyyyyy, g_0_y_xxxy_xyyyyyz, g_0_y_xxxy_xyyyyzz, g_0_y_xxxy_xyyyzzz, g_0_y_xxxy_xyyzzzz, g_0_y_xxxy_xyzzzzz, g_0_y_xxxy_xzzzzzz, g_0_y_xxxy_yyyyyyy, g_0_y_xxxy_yyyyyyz, g_0_y_xxxy_yyyyyzz, g_0_y_xxxy_yyyyzzz, g_0_y_xxxy_yyyzzzz, g_0_y_xxxy_yyzzzzz, g_0_y_xxxy_yzzzzzz, g_0_y_xxxy_zzzzzzz, g_0_y_xxy_xxxxxxx, g_0_y_xxy_xxxxxxxx, g_0_y_xxy_xxxxxxxy, g_0_y_xxy_xxxxxxxz, g_0_y_xxy_xxxxxxy, g_0_y_xxy_xxxxxxyy, g_0_y_xxy_xxxxxxyz, g_0_y_xxy_xxxxxxz, g_0_y_xxy_xxxxxxzz, g_0_y_xxy_xxxxxyy, g_0_y_xxy_xxxxxyyy, g_0_y_xxy_xxxxxyyz, g_0_y_xxy_xxxxxyz, g_0_y_xxy_xxxxxyzz, g_0_y_xxy_xxxxxzz, g_0_y_xxy_xxxxxzzz, g_0_y_xxy_xxxxyyy, g_0_y_xxy_xxxxyyyy, g_0_y_xxy_xxxxyyyz, g_0_y_xxy_xxxxyyz, g_0_y_xxy_xxxxyyzz, g_0_y_xxy_xxxxyzz, g_0_y_xxy_xxxxyzzz, g_0_y_xxy_xxxxzzz, g_0_y_xxy_xxxxzzzz, g_0_y_xxy_xxxyyyy, g_0_y_xxy_xxxyyyyy, g_0_y_xxy_xxxyyyyz, g_0_y_xxy_xxxyyyz, g_0_y_xxy_xxxyyyzz, g_0_y_xxy_xxxyyzz, g_0_y_xxy_xxxyyzzz, g_0_y_xxy_xxxyzzz, g_0_y_xxy_xxxyzzzz, g_0_y_xxy_xxxzzzz, g_0_y_xxy_xxxzzzzz, g_0_y_xxy_xxyyyyy, g_0_y_xxy_xxyyyyyy, g_0_y_xxy_xxyyyyyz, g_0_y_xxy_xxyyyyz, g_0_y_xxy_xxyyyyzz, g_0_y_xxy_xxyyyzz, g_0_y_xxy_xxyyyzzz, g_0_y_xxy_xxyyzzz, g_0_y_xxy_xxyyzzzz, g_0_y_xxy_xxyzzzz, g_0_y_xxy_xxyzzzzz, g_0_y_xxy_xxzzzzz, g_0_y_xxy_xxzzzzzz, g_0_y_xxy_xyyyyyy, g_0_y_xxy_xyyyyyyy, g_0_y_xxy_xyyyyyyz, g_0_y_xxy_xyyyyyz, g_0_y_xxy_xyyyyyzz, g_0_y_xxy_xyyyyzz, g_0_y_xxy_xyyyyzzz, g_0_y_xxy_xyyyzzz, g_0_y_xxy_xyyyzzzz, g_0_y_xxy_xyyzzzz, g_0_y_xxy_xyyzzzzz, g_0_y_xxy_xyzzzzz, g_0_y_xxy_xyzzzzzz, g_0_y_xxy_xzzzzzz, g_0_y_xxy_xzzzzzzz, g_0_y_xxy_yyyyyyy, g_0_y_xxy_yyyyyyz, g_0_y_xxy_yyyyyzz, g_0_y_xxy_yyyyzzz, g_0_y_xxy_yyyzzzz, g_0_y_xxy_yyzzzzz, g_0_y_xxy_yzzzzzz, g_0_y_xxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_xxxxxxx[k] = -g_0_y_xxy_xxxxxxx[k] * ab_x + g_0_y_xxy_xxxxxxxx[k];

                g_0_y_xxxy_xxxxxxy[k] = -g_0_y_xxy_xxxxxxy[k] * ab_x + g_0_y_xxy_xxxxxxxy[k];

                g_0_y_xxxy_xxxxxxz[k] = -g_0_y_xxy_xxxxxxz[k] * ab_x + g_0_y_xxy_xxxxxxxz[k];

                g_0_y_xxxy_xxxxxyy[k] = -g_0_y_xxy_xxxxxyy[k] * ab_x + g_0_y_xxy_xxxxxxyy[k];

                g_0_y_xxxy_xxxxxyz[k] = -g_0_y_xxy_xxxxxyz[k] * ab_x + g_0_y_xxy_xxxxxxyz[k];

                g_0_y_xxxy_xxxxxzz[k] = -g_0_y_xxy_xxxxxzz[k] * ab_x + g_0_y_xxy_xxxxxxzz[k];

                g_0_y_xxxy_xxxxyyy[k] = -g_0_y_xxy_xxxxyyy[k] * ab_x + g_0_y_xxy_xxxxxyyy[k];

                g_0_y_xxxy_xxxxyyz[k] = -g_0_y_xxy_xxxxyyz[k] * ab_x + g_0_y_xxy_xxxxxyyz[k];

                g_0_y_xxxy_xxxxyzz[k] = -g_0_y_xxy_xxxxyzz[k] * ab_x + g_0_y_xxy_xxxxxyzz[k];

                g_0_y_xxxy_xxxxzzz[k] = -g_0_y_xxy_xxxxzzz[k] * ab_x + g_0_y_xxy_xxxxxzzz[k];

                g_0_y_xxxy_xxxyyyy[k] = -g_0_y_xxy_xxxyyyy[k] * ab_x + g_0_y_xxy_xxxxyyyy[k];

                g_0_y_xxxy_xxxyyyz[k] = -g_0_y_xxy_xxxyyyz[k] * ab_x + g_0_y_xxy_xxxxyyyz[k];

                g_0_y_xxxy_xxxyyzz[k] = -g_0_y_xxy_xxxyyzz[k] * ab_x + g_0_y_xxy_xxxxyyzz[k];

                g_0_y_xxxy_xxxyzzz[k] = -g_0_y_xxy_xxxyzzz[k] * ab_x + g_0_y_xxy_xxxxyzzz[k];

                g_0_y_xxxy_xxxzzzz[k] = -g_0_y_xxy_xxxzzzz[k] * ab_x + g_0_y_xxy_xxxxzzzz[k];

                g_0_y_xxxy_xxyyyyy[k] = -g_0_y_xxy_xxyyyyy[k] * ab_x + g_0_y_xxy_xxxyyyyy[k];

                g_0_y_xxxy_xxyyyyz[k] = -g_0_y_xxy_xxyyyyz[k] * ab_x + g_0_y_xxy_xxxyyyyz[k];

                g_0_y_xxxy_xxyyyzz[k] = -g_0_y_xxy_xxyyyzz[k] * ab_x + g_0_y_xxy_xxxyyyzz[k];

                g_0_y_xxxy_xxyyzzz[k] = -g_0_y_xxy_xxyyzzz[k] * ab_x + g_0_y_xxy_xxxyyzzz[k];

                g_0_y_xxxy_xxyzzzz[k] = -g_0_y_xxy_xxyzzzz[k] * ab_x + g_0_y_xxy_xxxyzzzz[k];

                g_0_y_xxxy_xxzzzzz[k] = -g_0_y_xxy_xxzzzzz[k] * ab_x + g_0_y_xxy_xxxzzzzz[k];

                g_0_y_xxxy_xyyyyyy[k] = -g_0_y_xxy_xyyyyyy[k] * ab_x + g_0_y_xxy_xxyyyyyy[k];

                g_0_y_xxxy_xyyyyyz[k] = -g_0_y_xxy_xyyyyyz[k] * ab_x + g_0_y_xxy_xxyyyyyz[k];

                g_0_y_xxxy_xyyyyzz[k] = -g_0_y_xxy_xyyyyzz[k] * ab_x + g_0_y_xxy_xxyyyyzz[k];

                g_0_y_xxxy_xyyyzzz[k] = -g_0_y_xxy_xyyyzzz[k] * ab_x + g_0_y_xxy_xxyyyzzz[k];

                g_0_y_xxxy_xyyzzzz[k] = -g_0_y_xxy_xyyzzzz[k] * ab_x + g_0_y_xxy_xxyyzzzz[k];

                g_0_y_xxxy_xyzzzzz[k] = -g_0_y_xxy_xyzzzzz[k] * ab_x + g_0_y_xxy_xxyzzzzz[k];

                g_0_y_xxxy_xzzzzzz[k] = -g_0_y_xxy_xzzzzzz[k] * ab_x + g_0_y_xxy_xxzzzzzz[k];

                g_0_y_xxxy_yyyyyyy[k] = -g_0_y_xxy_yyyyyyy[k] * ab_x + g_0_y_xxy_xyyyyyyy[k];

                g_0_y_xxxy_yyyyyyz[k] = -g_0_y_xxy_yyyyyyz[k] * ab_x + g_0_y_xxy_xyyyyyyz[k];

                g_0_y_xxxy_yyyyyzz[k] = -g_0_y_xxy_yyyyyzz[k] * ab_x + g_0_y_xxy_xyyyyyzz[k];

                g_0_y_xxxy_yyyyzzz[k] = -g_0_y_xxy_yyyyzzz[k] * ab_x + g_0_y_xxy_xyyyyzzz[k];

                g_0_y_xxxy_yyyzzzz[k] = -g_0_y_xxy_yyyzzzz[k] * ab_x + g_0_y_xxy_xyyyzzzz[k];

                g_0_y_xxxy_yyzzzzz[k] = -g_0_y_xxy_yyzzzzz[k] * ab_x + g_0_y_xxy_xyyzzzzz[k];

                g_0_y_xxxy_yzzzzzz[k] = -g_0_y_xxy_yzzzzzz[k] * ab_x + g_0_y_xxy_xyzzzzzz[k];

                g_0_y_xxxy_zzzzzzz[k] = -g_0_y_xxy_zzzzzzz[k] * ab_x + g_0_y_xxy_xzzzzzzz[k];
            }

            /// Set up 612-648 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxxz_xxxxxxx, g_0_y_xxxz_xxxxxxy, g_0_y_xxxz_xxxxxxz, g_0_y_xxxz_xxxxxyy, g_0_y_xxxz_xxxxxyz, g_0_y_xxxz_xxxxxzz, g_0_y_xxxz_xxxxyyy, g_0_y_xxxz_xxxxyyz, g_0_y_xxxz_xxxxyzz, g_0_y_xxxz_xxxxzzz, g_0_y_xxxz_xxxyyyy, g_0_y_xxxz_xxxyyyz, g_0_y_xxxz_xxxyyzz, g_0_y_xxxz_xxxyzzz, g_0_y_xxxz_xxxzzzz, g_0_y_xxxz_xxyyyyy, g_0_y_xxxz_xxyyyyz, g_0_y_xxxz_xxyyyzz, g_0_y_xxxz_xxyyzzz, g_0_y_xxxz_xxyzzzz, g_0_y_xxxz_xxzzzzz, g_0_y_xxxz_xyyyyyy, g_0_y_xxxz_xyyyyyz, g_0_y_xxxz_xyyyyzz, g_0_y_xxxz_xyyyzzz, g_0_y_xxxz_xyyzzzz, g_0_y_xxxz_xyzzzzz, g_0_y_xxxz_xzzzzzz, g_0_y_xxxz_yyyyyyy, g_0_y_xxxz_yyyyyyz, g_0_y_xxxz_yyyyyzz, g_0_y_xxxz_yyyyzzz, g_0_y_xxxz_yyyzzzz, g_0_y_xxxz_yyzzzzz, g_0_y_xxxz_yzzzzzz, g_0_y_xxxz_zzzzzzz, g_0_y_xxz_xxxxxxx, g_0_y_xxz_xxxxxxxx, g_0_y_xxz_xxxxxxxy, g_0_y_xxz_xxxxxxxz, g_0_y_xxz_xxxxxxy, g_0_y_xxz_xxxxxxyy, g_0_y_xxz_xxxxxxyz, g_0_y_xxz_xxxxxxz, g_0_y_xxz_xxxxxxzz, g_0_y_xxz_xxxxxyy, g_0_y_xxz_xxxxxyyy, g_0_y_xxz_xxxxxyyz, g_0_y_xxz_xxxxxyz, g_0_y_xxz_xxxxxyzz, g_0_y_xxz_xxxxxzz, g_0_y_xxz_xxxxxzzz, g_0_y_xxz_xxxxyyy, g_0_y_xxz_xxxxyyyy, g_0_y_xxz_xxxxyyyz, g_0_y_xxz_xxxxyyz, g_0_y_xxz_xxxxyyzz, g_0_y_xxz_xxxxyzz, g_0_y_xxz_xxxxyzzz, g_0_y_xxz_xxxxzzz, g_0_y_xxz_xxxxzzzz, g_0_y_xxz_xxxyyyy, g_0_y_xxz_xxxyyyyy, g_0_y_xxz_xxxyyyyz, g_0_y_xxz_xxxyyyz, g_0_y_xxz_xxxyyyzz, g_0_y_xxz_xxxyyzz, g_0_y_xxz_xxxyyzzz, g_0_y_xxz_xxxyzzz, g_0_y_xxz_xxxyzzzz, g_0_y_xxz_xxxzzzz, g_0_y_xxz_xxxzzzzz, g_0_y_xxz_xxyyyyy, g_0_y_xxz_xxyyyyyy, g_0_y_xxz_xxyyyyyz, g_0_y_xxz_xxyyyyz, g_0_y_xxz_xxyyyyzz, g_0_y_xxz_xxyyyzz, g_0_y_xxz_xxyyyzzz, g_0_y_xxz_xxyyzzz, g_0_y_xxz_xxyyzzzz, g_0_y_xxz_xxyzzzz, g_0_y_xxz_xxyzzzzz, g_0_y_xxz_xxzzzzz, g_0_y_xxz_xxzzzzzz, g_0_y_xxz_xyyyyyy, g_0_y_xxz_xyyyyyyy, g_0_y_xxz_xyyyyyyz, g_0_y_xxz_xyyyyyz, g_0_y_xxz_xyyyyyzz, g_0_y_xxz_xyyyyzz, g_0_y_xxz_xyyyyzzz, g_0_y_xxz_xyyyzzz, g_0_y_xxz_xyyyzzzz, g_0_y_xxz_xyyzzzz, g_0_y_xxz_xyyzzzzz, g_0_y_xxz_xyzzzzz, g_0_y_xxz_xyzzzzzz, g_0_y_xxz_xzzzzzz, g_0_y_xxz_xzzzzzzz, g_0_y_xxz_yyyyyyy, g_0_y_xxz_yyyyyyz, g_0_y_xxz_yyyyyzz, g_0_y_xxz_yyyyzzz, g_0_y_xxz_yyyzzzz, g_0_y_xxz_yyzzzzz, g_0_y_xxz_yzzzzzz, g_0_y_xxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_xxxxxxx[k] = -g_0_y_xxz_xxxxxxx[k] * ab_x + g_0_y_xxz_xxxxxxxx[k];

                g_0_y_xxxz_xxxxxxy[k] = -g_0_y_xxz_xxxxxxy[k] * ab_x + g_0_y_xxz_xxxxxxxy[k];

                g_0_y_xxxz_xxxxxxz[k] = -g_0_y_xxz_xxxxxxz[k] * ab_x + g_0_y_xxz_xxxxxxxz[k];

                g_0_y_xxxz_xxxxxyy[k] = -g_0_y_xxz_xxxxxyy[k] * ab_x + g_0_y_xxz_xxxxxxyy[k];

                g_0_y_xxxz_xxxxxyz[k] = -g_0_y_xxz_xxxxxyz[k] * ab_x + g_0_y_xxz_xxxxxxyz[k];

                g_0_y_xxxz_xxxxxzz[k] = -g_0_y_xxz_xxxxxzz[k] * ab_x + g_0_y_xxz_xxxxxxzz[k];

                g_0_y_xxxz_xxxxyyy[k] = -g_0_y_xxz_xxxxyyy[k] * ab_x + g_0_y_xxz_xxxxxyyy[k];

                g_0_y_xxxz_xxxxyyz[k] = -g_0_y_xxz_xxxxyyz[k] * ab_x + g_0_y_xxz_xxxxxyyz[k];

                g_0_y_xxxz_xxxxyzz[k] = -g_0_y_xxz_xxxxyzz[k] * ab_x + g_0_y_xxz_xxxxxyzz[k];

                g_0_y_xxxz_xxxxzzz[k] = -g_0_y_xxz_xxxxzzz[k] * ab_x + g_0_y_xxz_xxxxxzzz[k];

                g_0_y_xxxz_xxxyyyy[k] = -g_0_y_xxz_xxxyyyy[k] * ab_x + g_0_y_xxz_xxxxyyyy[k];

                g_0_y_xxxz_xxxyyyz[k] = -g_0_y_xxz_xxxyyyz[k] * ab_x + g_0_y_xxz_xxxxyyyz[k];

                g_0_y_xxxz_xxxyyzz[k] = -g_0_y_xxz_xxxyyzz[k] * ab_x + g_0_y_xxz_xxxxyyzz[k];

                g_0_y_xxxz_xxxyzzz[k] = -g_0_y_xxz_xxxyzzz[k] * ab_x + g_0_y_xxz_xxxxyzzz[k];

                g_0_y_xxxz_xxxzzzz[k] = -g_0_y_xxz_xxxzzzz[k] * ab_x + g_0_y_xxz_xxxxzzzz[k];

                g_0_y_xxxz_xxyyyyy[k] = -g_0_y_xxz_xxyyyyy[k] * ab_x + g_0_y_xxz_xxxyyyyy[k];

                g_0_y_xxxz_xxyyyyz[k] = -g_0_y_xxz_xxyyyyz[k] * ab_x + g_0_y_xxz_xxxyyyyz[k];

                g_0_y_xxxz_xxyyyzz[k] = -g_0_y_xxz_xxyyyzz[k] * ab_x + g_0_y_xxz_xxxyyyzz[k];

                g_0_y_xxxz_xxyyzzz[k] = -g_0_y_xxz_xxyyzzz[k] * ab_x + g_0_y_xxz_xxxyyzzz[k];

                g_0_y_xxxz_xxyzzzz[k] = -g_0_y_xxz_xxyzzzz[k] * ab_x + g_0_y_xxz_xxxyzzzz[k];

                g_0_y_xxxz_xxzzzzz[k] = -g_0_y_xxz_xxzzzzz[k] * ab_x + g_0_y_xxz_xxxzzzzz[k];

                g_0_y_xxxz_xyyyyyy[k] = -g_0_y_xxz_xyyyyyy[k] * ab_x + g_0_y_xxz_xxyyyyyy[k];

                g_0_y_xxxz_xyyyyyz[k] = -g_0_y_xxz_xyyyyyz[k] * ab_x + g_0_y_xxz_xxyyyyyz[k];

                g_0_y_xxxz_xyyyyzz[k] = -g_0_y_xxz_xyyyyzz[k] * ab_x + g_0_y_xxz_xxyyyyzz[k];

                g_0_y_xxxz_xyyyzzz[k] = -g_0_y_xxz_xyyyzzz[k] * ab_x + g_0_y_xxz_xxyyyzzz[k];

                g_0_y_xxxz_xyyzzzz[k] = -g_0_y_xxz_xyyzzzz[k] * ab_x + g_0_y_xxz_xxyyzzzz[k];

                g_0_y_xxxz_xyzzzzz[k] = -g_0_y_xxz_xyzzzzz[k] * ab_x + g_0_y_xxz_xxyzzzzz[k];

                g_0_y_xxxz_xzzzzzz[k] = -g_0_y_xxz_xzzzzzz[k] * ab_x + g_0_y_xxz_xxzzzzzz[k];

                g_0_y_xxxz_yyyyyyy[k] = -g_0_y_xxz_yyyyyyy[k] * ab_x + g_0_y_xxz_xyyyyyyy[k];

                g_0_y_xxxz_yyyyyyz[k] = -g_0_y_xxz_yyyyyyz[k] * ab_x + g_0_y_xxz_xyyyyyyz[k];

                g_0_y_xxxz_yyyyyzz[k] = -g_0_y_xxz_yyyyyzz[k] * ab_x + g_0_y_xxz_xyyyyyzz[k];

                g_0_y_xxxz_yyyyzzz[k] = -g_0_y_xxz_yyyyzzz[k] * ab_x + g_0_y_xxz_xyyyyzzz[k];

                g_0_y_xxxz_yyyzzzz[k] = -g_0_y_xxz_yyyzzzz[k] * ab_x + g_0_y_xxz_xyyyzzzz[k];

                g_0_y_xxxz_yyzzzzz[k] = -g_0_y_xxz_yyzzzzz[k] * ab_x + g_0_y_xxz_xyyzzzzz[k];

                g_0_y_xxxz_yzzzzzz[k] = -g_0_y_xxz_yzzzzzz[k] * ab_x + g_0_y_xxz_xyzzzzzz[k];

                g_0_y_xxxz_zzzzzzz[k] = -g_0_y_xxz_zzzzzzz[k] * ab_x + g_0_y_xxz_xzzzzzzz[k];
            }

            /// Set up 648-684 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxyy_xxxxxxx, g_0_y_xxyy_xxxxxxy, g_0_y_xxyy_xxxxxxz, g_0_y_xxyy_xxxxxyy, g_0_y_xxyy_xxxxxyz, g_0_y_xxyy_xxxxxzz, g_0_y_xxyy_xxxxyyy, g_0_y_xxyy_xxxxyyz, g_0_y_xxyy_xxxxyzz, g_0_y_xxyy_xxxxzzz, g_0_y_xxyy_xxxyyyy, g_0_y_xxyy_xxxyyyz, g_0_y_xxyy_xxxyyzz, g_0_y_xxyy_xxxyzzz, g_0_y_xxyy_xxxzzzz, g_0_y_xxyy_xxyyyyy, g_0_y_xxyy_xxyyyyz, g_0_y_xxyy_xxyyyzz, g_0_y_xxyy_xxyyzzz, g_0_y_xxyy_xxyzzzz, g_0_y_xxyy_xxzzzzz, g_0_y_xxyy_xyyyyyy, g_0_y_xxyy_xyyyyyz, g_0_y_xxyy_xyyyyzz, g_0_y_xxyy_xyyyzzz, g_0_y_xxyy_xyyzzzz, g_0_y_xxyy_xyzzzzz, g_0_y_xxyy_xzzzzzz, g_0_y_xxyy_yyyyyyy, g_0_y_xxyy_yyyyyyz, g_0_y_xxyy_yyyyyzz, g_0_y_xxyy_yyyyzzz, g_0_y_xxyy_yyyzzzz, g_0_y_xxyy_yyzzzzz, g_0_y_xxyy_yzzzzzz, g_0_y_xxyy_zzzzzzz, g_0_y_xyy_xxxxxxx, g_0_y_xyy_xxxxxxxx, g_0_y_xyy_xxxxxxxy, g_0_y_xyy_xxxxxxxz, g_0_y_xyy_xxxxxxy, g_0_y_xyy_xxxxxxyy, g_0_y_xyy_xxxxxxyz, g_0_y_xyy_xxxxxxz, g_0_y_xyy_xxxxxxzz, g_0_y_xyy_xxxxxyy, g_0_y_xyy_xxxxxyyy, g_0_y_xyy_xxxxxyyz, g_0_y_xyy_xxxxxyz, g_0_y_xyy_xxxxxyzz, g_0_y_xyy_xxxxxzz, g_0_y_xyy_xxxxxzzz, g_0_y_xyy_xxxxyyy, g_0_y_xyy_xxxxyyyy, g_0_y_xyy_xxxxyyyz, g_0_y_xyy_xxxxyyz, g_0_y_xyy_xxxxyyzz, g_0_y_xyy_xxxxyzz, g_0_y_xyy_xxxxyzzz, g_0_y_xyy_xxxxzzz, g_0_y_xyy_xxxxzzzz, g_0_y_xyy_xxxyyyy, g_0_y_xyy_xxxyyyyy, g_0_y_xyy_xxxyyyyz, g_0_y_xyy_xxxyyyz, g_0_y_xyy_xxxyyyzz, g_0_y_xyy_xxxyyzz, g_0_y_xyy_xxxyyzzz, g_0_y_xyy_xxxyzzz, g_0_y_xyy_xxxyzzzz, g_0_y_xyy_xxxzzzz, g_0_y_xyy_xxxzzzzz, g_0_y_xyy_xxyyyyy, g_0_y_xyy_xxyyyyyy, g_0_y_xyy_xxyyyyyz, g_0_y_xyy_xxyyyyz, g_0_y_xyy_xxyyyyzz, g_0_y_xyy_xxyyyzz, g_0_y_xyy_xxyyyzzz, g_0_y_xyy_xxyyzzz, g_0_y_xyy_xxyyzzzz, g_0_y_xyy_xxyzzzz, g_0_y_xyy_xxyzzzzz, g_0_y_xyy_xxzzzzz, g_0_y_xyy_xxzzzzzz, g_0_y_xyy_xyyyyyy, g_0_y_xyy_xyyyyyyy, g_0_y_xyy_xyyyyyyz, g_0_y_xyy_xyyyyyz, g_0_y_xyy_xyyyyyzz, g_0_y_xyy_xyyyyzz, g_0_y_xyy_xyyyyzzz, g_0_y_xyy_xyyyzzz, g_0_y_xyy_xyyyzzzz, g_0_y_xyy_xyyzzzz, g_0_y_xyy_xyyzzzzz, g_0_y_xyy_xyzzzzz, g_0_y_xyy_xyzzzzzz, g_0_y_xyy_xzzzzzz, g_0_y_xyy_xzzzzzzz, g_0_y_xyy_yyyyyyy, g_0_y_xyy_yyyyyyz, g_0_y_xyy_yyyyyzz, g_0_y_xyy_yyyyzzz, g_0_y_xyy_yyyzzzz, g_0_y_xyy_yyzzzzz, g_0_y_xyy_yzzzzzz, g_0_y_xyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_xxxxxxx[k] = -g_0_y_xyy_xxxxxxx[k] * ab_x + g_0_y_xyy_xxxxxxxx[k];

                g_0_y_xxyy_xxxxxxy[k] = -g_0_y_xyy_xxxxxxy[k] * ab_x + g_0_y_xyy_xxxxxxxy[k];

                g_0_y_xxyy_xxxxxxz[k] = -g_0_y_xyy_xxxxxxz[k] * ab_x + g_0_y_xyy_xxxxxxxz[k];

                g_0_y_xxyy_xxxxxyy[k] = -g_0_y_xyy_xxxxxyy[k] * ab_x + g_0_y_xyy_xxxxxxyy[k];

                g_0_y_xxyy_xxxxxyz[k] = -g_0_y_xyy_xxxxxyz[k] * ab_x + g_0_y_xyy_xxxxxxyz[k];

                g_0_y_xxyy_xxxxxzz[k] = -g_0_y_xyy_xxxxxzz[k] * ab_x + g_0_y_xyy_xxxxxxzz[k];

                g_0_y_xxyy_xxxxyyy[k] = -g_0_y_xyy_xxxxyyy[k] * ab_x + g_0_y_xyy_xxxxxyyy[k];

                g_0_y_xxyy_xxxxyyz[k] = -g_0_y_xyy_xxxxyyz[k] * ab_x + g_0_y_xyy_xxxxxyyz[k];

                g_0_y_xxyy_xxxxyzz[k] = -g_0_y_xyy_xxxxyzz[k] * ab_x + g_0_y_xyy_xxxxxyzz[k];

                g_0_y_xxyy_xxxxzzz[k] = -g_0_y_xyy_xxxxzzz[k] * ab_x + g_0_y_xyy_xxxxxzzz[k];

                g_0_y_xxyy_xxxyyyy[k] = -g_0_y_xyy_xxxyyyy[k] * ab_x + g_0_y_xyy_xxxxyyyy[k];

                g_0_y_xxyy_xxxyyyz[k] = -g_0_y_xyy_xxxyyyz[k] * ab_x + g_0_y_xyy_xxxxyyyz[k];

                g_0_y_xxyy_xxxyyzz[k] = -g_0_y_xyy_xxxyyzz[k] * ab_x + g_0_y_xyy_xxxxyyzz[k];

                g_0_y_xxyy_xxxyzzz[k] = -g_0_y_xyy_xxxyzzz[k] * ab_x + g_0_y_xyy_xxxxyzzz[k];

                g_0_y_xxyy_xxxzzzz[k] = -g_0_y_xyy_xxxzzzz[k] * ab_x + g_0_y_xyy_xxxxzzzz[k];

                g_0_y_xxyy_xxyyyyy[k] = -g_0_y_xyy_xxyyyyy[k] * ab_x + g_0_y_xyy_xxxyyyyy[k];

                g_0_y_xxyy_xxyyyyz[k] = -g_0_y_xyy_xxyyyyz[k] * ab_x + g_0_y_xyy_xxxyyyyz[k];

                g_0_y_xxyy_xxyyyzz[k] = -g_0_y_xyy_xxyyyzz[k] * ab_x + g_0_y_xyy_xxxyyyzz[k];

                g_0_y_xxyy_xxyyzzz[k] = -g_0_y_xyy_xxyyzzz[k] * ab_x + g_0_y_xyy_xxxyyzzz[k];

                g_0_y_xxyy_xxyzzzz[k] = -g_0_y_xyy_xxyzzzz[k] * ab_x + g_0_y_xyy_xxxyzzzz[k];

                g_0_y_xxyy_xxzzzzz[k] = -g_0_y_xyy_xxzzzzz[k] * ab_x + g_0_y_xyy_xxxzzzzz[k];

                g_0_y_xxyy_xyyyyyy[k] = -g_0_y_xyy_xyyyyyy[k] * ab_x + g_0_y_xyy_xxyyyyyy[k];

                g_0_y_xxyy_xyyyyyz[k] = -g_0_y_xyy_xyyyyyz[k] * ab_x + g_0_y_xyy_xxyyyyyz[k];

                g_0_y_xxyy_xyyyyzz[k] = -g_0_y_xyy_xyyyyzz[k] * ab_x + g_0_y_xyy_xxyyyyzz[k];

                g_0_y_xxyy_xyyyzzz[k] = -g_0_y_xyy_xyyyzzz[k] * ab_x + g_0_y_xyy_xxyyyzzz[k];

                g_0_y_xxyy_xyyzzzz[k] = -g_0_y_xyy_xyyzzzz[k] * ab_x + g_0_y_xyy_xxyyzzzz[k];

                g_0_y_xxyy_xyzzzzz[k] = -g_0_y_xyy_xyzzzzz[k] * ab_x + g_0_y_xyy_xxyzzzzz[k];

                g_0_y_xxyy_xzzzzzz[k] = -g_0_y_xyy_xzzzzzz[k] * ab_x + g_0_y_xyy_xxzzzzzz[k];

                g_0_y_xxyy_yyyyyyy[k] = -g_0_y_xyy_yyyyyyy[k] * ab_x + g_0_y_xyy_xyyyyyyy[k];

                g_0_y_xxyy_yyyyyyz[k] = -g_0_y_xyy_yyyyyyz[k] * ab_x + g_0_y_xyy_xyyyyyyz[k];

                g_0_y_xxyy_yyyyyzz[k] = -g_0_y_xyy_yyyyyzz[k] * ab_x + g_0_y_xyy_xyyyyyzz[k];

                g_0_y_xxyy_yyyyzzz[k] = -g_0_y_xyy_yyyyzzz[k] * ab_x + g_0_y_xyy_xyyyyzzz[k];

                g_0_y_xxyy_yyyzzzz[k] = -g_0_y_xyy_yyyzzzz[k] * ab_x + g_0_y_xyy_xyyyzzzz[k];

                g_0_y_xxyy_yyzzzzz[k] = -g_0_y_xyy_yyzzzzz[k] * ab_x + g_0_y_xyy_xyyzzzzz[k];

                g_0_y_xxyy_yzzzzzz[k] = -g_0_y_xyy_yzzzzzz[k] * ab_x + g_0_y_xyy_xyzzzzzz[k];

                g_0_y_xxyy_zzzzzzz[k] = -g_0_y_xyy_zzzzzzz[k] * ab_x + g_0_y_xyy_xzzzzzzz[k];
            }

            /// Set up 684-720 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxyz_xxxxxxx, g_0_y_xxyz_xxxxxxy, g_0_y_xxyz_xxxxxxz, g_0_y_xxyz_xxxxxyy, g_0_y_xxyz_xxxxxyz, g_0_y_xxyz_xxxxxzz, g_0_y_xxyz_xxxxyyy, g_0_y_xxyz_xxxxyyz, g_0_y_xxyz_xxxxyzz, g_0_y_xxyz_xxxxzzz, g_0_y_xxyz_xxxyyyy, g_0_y_xxyz_xxxyyyz, g_0_y_xxyz_xxxyyzz, g_0_y_xxyz_xxxyzzz, g_0_y_xxyz_xxxzzzz, g_0_y_xxyz_xxyyyyy, g_0_y_xxyz_xxyyyyz, g_0_y_xxyz_xxyyyzz, g_0_y_xxyz_xxyyzzz, g_0_y_xxyz_xxyzzzz, g_0_y_xxyz_xxzzzzz, g_0_y_xxyz_xyyyyyy, g_0_y_xxyz_xyyyyyz, g_0_y_xxyz_xyyyyzz, g_0_y_xxyz_xyyyzzz, g_0_y_xxyz_xyyzzzz, g_0_y_xxyz_xyzzzzz, g_0_y_xxyz_xzzzzzz, g_0_y_xxyz_yyyyyyy, g_0_y_xxyz_yyyyyyz, g_0_y_xxyz_yyyyyzz, g_0_y_xxyz_yyyyzzz, g_0_y_xxyz_yyyzzzz, g_0_y_xxyz_yyzzzzz, g_0_y_xxyz_yzzzzzz, g_0_y_xxyz_zzzzzzz, g_0_y_xyz_xxxxxxx, g_0_y_xyz_xxxxxxxx, g_0_y_xyz_xxxxxxxy, g_0_y_xyz_xxxxxxxz, g_0_y_xyz_xxxxxxy, g_0_y_xyz_xxxxxxyy, g_0_y_xyz_xxxxxxyz, g_0_y_xyz_xxxxxxz, g_0_y_xyz_xxxxxxzz, g_0_y_xyz_xxxxxyy, g_0_y_xyz_xxxxxyyy, g_0_y_xyz_xxxxxyyz, g_0_y_xyz_xxxxxyz, g_0_y_xyz_xxxxxyzz, g_0_y_xyz_xxxxxzz, g_0_y_xyz_xxxxxzzz, g_0_y_xyz_xxxxyyy, g_0_y_xyz_xxxxyyyy, g_0_y_xyz_xxxxyyyz, g_0_y_xyz_xxxxyyz, g_0_y_xyz_xxxxyyzz, g_0_y_xyz_xxxxyzz, g_0_y_xyz_xxxxyzzz, g_0_y_xyz_xxxxzzz, g_0_y_xyz_xxxxzzzz, g_0_y_xyz_xxxyyyy, g_0_y_xyz_xxxyyyyy, g_0_y_xyz_xxxyyyyz, g_0_y_xyz_xxxyyyz, g_0_y_xyz_xxxyyyzz, g_0_y_xyz_xxxyyzz, g_0_y_xyz_xxxyyzzz, g_0_y_xyz_xxxyzzz, g_0_y_xyz_xxxyzzzz, g_0_y_xyz_xxxzzzz, g_0_y_xyz_xxxzzzzz, g_0_y_xyz_xxyyyyy, g_0_y_xyz_xxyyyyyy, g_0_y_xyz_xxyyyyyz, g_0_y_xyz_xxyyyyz, g_0_y_xyz_xxyyyyzz, g_0_y_xyz_xxyyyzz, g_0_y_xyz_xxyyyzzz, g_0_y_xyz_xxyyzzz, g_0_y_xyz_xxyyzzzz, g_0_y_xyz_xxyzzzz, g_0_y_xyz_xxyzzzzz, g_0_y_xyz_xxzzzzz, g_0_y_xyz_xxzzzzzz, g_0_y_xyz_xyyyyyy, g_0_y_xyz_xyyyyyyy, g_0_y_xyz_xyyyyyyz, g_0_y_xyz_xyyyyyz, g_0_y_xyz_xyyyyyzz, g_0_y_xyz_xyyyyzz, g_0_y_xyz_xyyyyzzz, g_0_y_xyz_xyyyzzz, g_0_y_xyz_xyyyzzzz, g_0_y_xyz_xyyzzzz, g_0_y_xyz_xyyzzzzz, g_0_y_xyz_xyzzzzz, g_0_y_xyz_xyzzzzzz, g_0_y_xyz_xzzzzzz, g_0_y_xyz_xzzzzzzz, g_0_y_xyz_yyyyyyy, g_0_y_xyz_yyyyyyz, g_0_y_xyz_yyyyyzz, g_0_y_xyz_yyyyzzz, g_0_y_xyz_yyyzzzz, g_0_y_xyz_yyzzzzz, g_0_y_xyz_yzzzzzz, g_0_y_xyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_xxxxxxx[k] = -g_0_y_xyz_xxxxxxx[k] * ab_x + g_0_y_xyz_xxxxxxxx[k];

                g_0_y_xxyz_xxxxxxy[k] = -g_0_y_xyz_xxxxxxy[k] * ab_x + g_0_y_xyz_xxxxxxxy[k];

                g_0_y_xxyz_xxxxxxz[k] = -g_0_y_xyz_xxxxxxz[k] * ab_x + g_0_y_xyz_xxxxxxxz[k];

                g_0_y_xxyz_xxxxxyy[k] = -g_0_y_xyz_xxxxxyy[k] * ab_x + g_0_y_xyz_xxxxxxyy[k];

                g_0_y_xxyz_xxxxxyz[k] = -g_0_y_xyz_xxxxxyz[k] * ab_x + g_0_y_xyz_xxxxxxyz[k];

                g_0_y_xxyz_xxxxxzz[k] = -g_0_y_xyz_xxxxxzz[k] * ab_x + g_0_y_xyz_xxxxxxzz[k];

                g_0_y_xxyz_xxxxyyy[k] = -g_0_y_xyz_xxxxyyy[k] * ab_x + g_0_y_xyz_xxxxxyyy[k];

                g_0_y_xxyz_xxxxyyz[k] = -g_0_y_xyz_xxxxyyz[k] * ab_x + g_0_y_xyz_xxxxxyyz[k];

                g_0_y_xxyz_xxxxyzz[k] = -g_0_y_xyz_xxxxyzz[k] * ab_x + g_0_y_xyz_xxxxxyzz[k];

                g_0_y_xxyz_xxxxzzz[k] = -g_0_y_xyz_xxxxzzz[k] * ab_x + g_0_y_xyz_xxxxxzzz[k];

                g_0_y_xxyz_xxxyyyy[k] = -g_0_y_xyz_xxxyyyy[k] * ab_x + g_0_y_xyz_xxxxyyyy[k];

                g_0_y_xxyz_xxxyyyz[k] = -g_0_y_xyz_xxxyyyz[k] * ab_x + g_0_y_xyz_xxxxyyyz[k];

                g_0_y_xxyz_xxxyyzz[k] = -g_0_y_xyz_xxxyyzz[k] * ab_x + g_0_y_xyz_xxxxyyzz[k];

                g_0_y_xxyz_xxxyzzz[k] = -g_0_y_xyz_xxxyzzz[k] * ab_x + g_0_y_xyz_xxxxyzzz[k];

                g_0_y_xxyz_xxxzzzz[k] = -g_0_y_xyz_xxxzzzz[k] * ab_x + g_0_y_xyz_xxxxzzzz[k];

                g_0_y_xxyz_xxyyyyy[k] = -g_0_y_xyz_xxyyyyy[k] * ab_x + g_0_y_xyz_xxxyyyyy[k];

                g_0_y_xxyz_xxyyyyz[k] = -g_0_y_xyz_xxyyyyz[k] * ab_x + g_0_y_xyz_xxxyyyyz[k];

                g_0_y_xxyz_xxyyyzz[k] = -g_0_y_xyz_xxyyyzz[k] * ab_x + g_0_y_xyz_xxxyyyzz[k];

                g_0_y_xxyz_xxyyzzz[k] = -g_0_y_xyz_xxyyzzz[k] * ab_x + g_0_y_xyz_xxxyyzzz[k];

                g_0_y_xxyz_xxyzzzz[k] = -g_0_y_xyz_xxyzzzz[k] * ab_x + g_0_y_xyz_xxxyzzzz[k];

                g_0_y_xxyz_xxzzzzz[k] = -g_0_y_xyz_xxzzzzz[k] * ab_x + g_0_y_xyz_xxxzzzzz[k];

                g_0_y_xxyz_xyyyyyy[k] = -g_0_y_xyz_xyyyyyy[k] * ab_x + g_0_y_xyz_xxyyyyyy[k];

                g_0_y_xxyz_xyyyyyz[k] = -g_0_y_xyz_xyyyyyz[k] * ab_x + g_0_y_xyz_xxyyyyyz[k];

                g_0_y_xxyz_xyyyyzz[k] = -g_0_y_xyz_xyyyyzz[k] * ab_x + g_0_y_xyz_xxyyyyzz[k];

                g_0_y_xxyz_xyyyzzz[k] = -g_0_y_xyz_xyyyzzz[k] * ab_x + g_0_y_xyz_xxyyyzzz[k];

                g_0_y_xxyz_xyyzzzz[k] = -g_0_y_xyz_xyyzzzz[k] * ab_x + g_0_y_xyz_xxyyzzzz[k];

                g_0_y_xxyz_xyzzzzz[k] = -g_0_y_xyz_xyzzzzz[k] * ab_x + g_0_y_xyz_xxyzzzzz[k];

                g_0_y_xxyz_xzzzzzz[k] = -g_0_y_xyz_xzzzzzz[k] * ab_x + g_0_y_xyz_xxzzzzzz[k];

                g_0_y_xxyz_yyyyyyy[k] = -g_0_y_xyz_yyyyyyy[k] * ab_x + g_0_y_xyz_xyyyyyyy[k];

                g_0_y_xxyz_yyyyyyz[k] = -g_0_y_xyz_yyyyyyz[k] * ab_x + g_0_y_xyz_xyyyyyyz[k];

                g_0_y_xxyz_yyyyyzz[k] = -g_0_y_xyz_yyyyyzz[k] * ab_x + g_0_y_xyz_xyyyyyzz[k];

                g_0_y_xxyz_yyyyzzz[k] = -g_0_y_xyz_yyyyzzz[k] * ab_x + g_0_y_xyz_xyyyyzzz[k];

                g_0_y_xxyz_yyyzzzz[k] = -g_0_y_xyz_yyyzzzz[k] * ab_x + g_0_y_xyz_xyyyzzzz[k];

                g_0_y_xxyz_yyzzzzz[k] = -g_0_y_xyz_yyzzzzz[k] * ab_x + g_0_y_xyz_xyyzzzzz[k];

                g_0_y_xxyz_yzzzzzz[k] = -g_0_y_xyz_yzzzzzz[k] * ab_x + g_0_y_xyz_xyzzzzzz[k];

                g_0_y_xxyz_zzzzzzz[k] = -g_0_y_xyz_zzzzzzz[k] * ab_x + g_0_y_xyz_xzzzzzzz[k];
            }

            /// Set up 720-756 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxzz_xxxxxxx, g_0_y_xxzz_xxxxxxy, g_0_y_xxzz_xxxxxxz, g_0_y_xxzz_xxxxxyy, g_0_y_xxzz_xxxxxyz, g_0_y_xxzz_xxxxxzz, g_0_y_xxzz_xxxxyyy, g_0_y_xxzz_xxxxyyz, g_0_y_xxzz_xxxxyzz, g_0_y_xxzz_xxxxzzz, g_0_y_xxzz_xxxyyyy, g_0_y_xxzz_xxxyyyz, g_0_y_xxzz_xxxyyzz, g_0_y_xxzz_xxxyzzz, g_0_y_xxzz_xxxzzzz, g_0_y_xxzz_xxyyyyy, g_0_y_xxzz_xxyyyyz, g_0_y_xxzz_xxyyyzz, g_0_y_xxzz_xxyyzzz, g_0_y_xxzz_xxyzzzz, g_0_y_xxzz_xxzzzzz, g_0_y_xxzz_xyyyyyy, g_0_y_xxzz_xyyyyyz, g_0_y_xxzz_xyyyyzz, g_0_y_xxzz_xyyyzzz, g_0_y_xxzz_xyyzzzz, g_0_y_xxzz_xyzzzzz, g_0_y_xxzz_xzzzzzz, g_0_y_xxzz_yyyyyyy, g_0_y_xxzz_yyyyyyz, g_0_y_xxzz_yyyyyzz, g_0_y_xxzz_yyyyzzz, g_0_y_xxzz_yyyzzzz, g_0_y_xxzz_yyzzzzz, g_0_y_xxzz_yzzzzzz, g_0_y_xxzz_zzzzzzz, g_0_y_xzz_xxxxxxx, g_0_y_xzz_xxxxxxxx, g_0_y_xzz_xxxxxxxy, g_0_y_xzz_xxxxxxxz, g_0_y_xzz_xxxxxxy, g_0_y_xzz_xxxxxxyy, g_0_y_xzz_xxxxxxyz, g_0_y_xzz_xxxxxxz, g_0_y_xzz_xxxxxxzz, g_0_y_xzz_xxxxxyy, g_0_y_xzz_xxxxxyyy, g_0_y_xzz_xxxxxyyz, g_0_y_xzz_xxxxxyz, g_0_y_xzz_xxxxxyzz, g_0_y_xzz_xxxxxzz, g_0_y_xzz_xxxxxzzz, g_0_y_xzz_xxxxyyy, g_0_y_xzz_xxxxyyyy, g_0_y_xzz_xxxxyyyz, g_0_y_xzz_xxxxyyz, g_0_y_xzz_xxxxyyzz, g_0_y_xzz_xxxxyzz, g_0_y_xzz_xxxxyzzz, g_0_y_xzz_xxxxzzz, g_0_y_xzz_xxxxzzzz, g_0_y_xzz_xxxyyyy, g_0_y_xzz_xxxyyyyy, g_0_y_xzz_xxxyyyyz, g_0_y_xzz_xxxyyyz, g_0_y_xzz_xxxyyyzz, g_0_y_xzz_xxxyyzz, g_0_y_xzz_xxxyyzzz, g_0_y_xzz_xxxyzzz, g_0_y_xzz_xxxyzzzz, g_0_y_xzz_xxxzzzz, g_0_y_xzz_xxxzzzzz, g_0_y_xzz_xxyyyyy, g_0_y_xzz_xxyyyyyy, g_0_y_xzz_xxyyyyyz, g_0_y_xzz_xxyyyyz, g_0_y_xzz_xxyyyyzz, g_0_y_xzz_xxyyyzz, g_0_y_xzz_xxyyyzzz, g_0_y_xzz_xxyyzzz, g_0_y_xzz_xxyyzzzz, g_0_y_xzz_xxyzzzz, g_0_y_xzz_xxyzzzzz, g_0_y_xzz_xxzzzzz, g_0_y_xzz_xxzzzzzz, g_0_y_xzz_xyyyyyy, g_0_y_xzz_xyyyyyyy, g_0_y_xzz_xyyyyyyz, g_0_y_xzz_xyyyyyz, g_0_y_xzz_xyyyyyzz, g_0_y_xzz_xyyyyzz, g_0_y_xzz_xyyyyzzz, g_0_y_xzz_xyyyzzz, g_0_y_xzz_xyyyzzzz, g_0_y_xzz_xyyzzzz, g_0_y_xzz_xyyzzzzz, g_0_y_xzz_xyzzzzz, g_0_y_xzz_xyzzzzzz, g_0_y_xzz_xzzzzzz, g_0_y_xzz_xzzzzzzz, g_0_y_xzz_yyyyyyy, g_0_y_xzz_yyyyyyz, g_0_y_xzz_yyyyyzz, g_0_y_xzz_yyyyzzz, g_0_y_xzz_yyyzzzz, g_0_y_xzz_yyzzzzz, g_0_y_xzz_yzzzzzz, g_0_y_xzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_xxxxxxx[k] = -g_0_y_xzz_xxxxxxx[k] * ab_x + g_0_y_xzz_xxxxxxxx[k];

                g_0_y_xxzz_xxxxxxy[k] = -g_0_y_xzz_xxxxxxy[k] * ab_x + g_0_y_xzz_xxxxxxxy[k];

                g_0_y_xxzz_xxxxxxz[k] = -g_0_y_xzz_xxxxxxz[k] * ab_x + g_0_y_xzz_xxxxxxxz[k];

                g_0_y_xxzz_xxxxxyy[k] = -g_0_y_xzz_xxxxxyy[k] * ab_x + g_0_y_xzz_xxxxxxyy[k];

                g_0_y_xxzz_xxxxxyz[k] = -g_0_y_xzz_xxxxxyz[k] * ab_x + g_0_y_xzz_xxxxxxyz[k];

                g_0_y_xxzz_xxxxxzz[k] = -g_0_y_xzz_xxxxxzz[k] * ab_x + g_0_y_xzz_xxxxxxzz[k];

                g_0_y_xxzz_xxxxyyy[k] = -g_0_y_xzz_xxxxyyy[k] * ab_x + g_0_y_xzz_xxxxxyyy[k];

                g_0_y_xxzz_xxxxyyz[k] = -g_0_y_xzz_xxxxyyz[k] * ab_x + g_0_y_xzz_xxxxxyyz[k];

                g_0_y_xxzz_xxxxyzz[k] = -g_0_y_xzz_xxxxyzz[k] * ab_x + g_0_y_xzz_xxxxxyzz[k];

                g_0_y_xxzz_xxxxzzz[k] = -g_0_y_xzz_xxxxzzz[k] * ab_x + g_0_y_xzz_xxxxxzzz[k];

                g_0_y_xxzz_xxxyyyy[k] = -g_0_y_xzz_xxxyyyy[k] * ab_x + g_0_y_xzz_xxxxyyyy[k];

                g_0_y_xxzz_xxxyyyz[k] = -g_0_y_xzz_xxxyyyz[k] * ab_x + g_0_y_xzz_xxxxyyyz[k];

                g_0_y_xxzz_xxxyyzz[k] = -g_0_y_xzz_xxxyyzz[k] * ab_x + g_0_y_xzz_xxxxyyzz[k];

                g_0_y_xxzz_xxxyzzz[k] = -g_0_y_xzz_xxxyzzz[k] * ab_x + g_0_y_xzz_xxxxyzzz[k];

                g_0_y_xxzz_xxxzzzz[k] = -g_0_y_xzz_xxxzzzz[k] * ab_x + g_0_y_xzz_xxxxzzzz[k];

                g_0_y_xxzz_xxyyyyy[k] = -g_0_y_xzz_xxyyyyy[k] * ab_x + g_0_y_xzz_xxxyyyyy[k];

                g_0_y_xxzz_xxyyyyz[k] = -g_0_y_xzz_xxyyyyz[k] * ab_x + g_0_y_xzz_xxxyyyyz[k];

                g_0_y_xxzz_xxyyyzz[k] = -g_0_y_xzz_xxyyyzz[k] * ab_x + g_0_y_xzz_xxxyyyzz[k];

                g_0_y_xxzz_xxyyzzz[k] = -g_0_y_xzz_xxyyzzz[k] * ab_x + g_0_y_xzz_xxxyyzzz[k];

                g_0_y_xxzz_xxyzzzz[k] = -g_0_y_xzz_xxyzzzz[k] * ab_x + g_0_y_xzz_xxxyzzzz[k];

                g_0_y_xxzz_xxzzzzz[k] = -g_0_y_xzz_xxzzzzz[k] * ab_x + g_0_y_xzz_xxxzzzzz[k];

                g_0_y_xxzz_xyyyyyy[k] = -g_0_y_xzz_xyyyyyy[k] * ab_x + g_0_y_xzz_xxyyyyyy[k];

                g_0_y_xxzz_xyyyyyz[k] = -g_0_y_xzz_xyyyyyz[k] * ab_x + g_0_y_xzz_xxyyyyyz[k];

                g_0_y_xxzz_xyyyyzz[k] = -g_0_y_xzz_xyyyyzz[k] * ab_x + g_0_y_xzz_xxyyyyzz[k];

                g_0_y_xxzz_xyyyzzz[k] = -g_0_y_xzz_xyyyzzz[k] * ab_x + g_0_y_xzz_xxyyyzzz[k];

                g_0_y_xxzz_xyyzzzz[k] = -g_0_y_xzz_xyyzzzz[k] * ab_x + g_0_y_xzz_xxyyzzzz[k];

                g_0_y_xxzz_xyzzzzz[k] = -g_0_y_xzz_xyzzzzz[k] * ab_x + g_0_y_xzz_xxyzzzzz[k];

                g_0_y_xxzz_xzzzzzz[k] = -g_0_y_xzz_xzzzzzz[k] * ab_x + g_0_y_xzz_xxzzzzzz[k];

                g_0_y_xxzz_yyyyyyy[k] = -g_0_y_xzz_yyyyyyy[k] * ab_x + g_0_y_xzz_xyyyyyyy[k];

                g_0_y_xxzz_yyyyyyz[k] = -g_0_y_xzz_yyyyyyz[k] * ab_x + g_0_y_xzz_xyyyyyyz[k];

                g_0_y_xxzz_yyyyyzz[k] = -g_0_y_xzz_yyyyyzz[k] * ab_x + g_0_y_xzz_xyyyyyzz[k];

                g_0_y_xxzz_yyyyzzz[k] = -g_0_y_xzz_yyyyzzz[k] * ab_x + g_0_y_xzz_xyyyyzzz[k];

                g_0_y_xxzz_yyyzzzz[k] = -g_0_y_xzz_yyyzzzz[k] * ab_x + g_0_y_xzz_xyyyzzzz[k];

                g_0_y_xxzz_yyzzzzz[k] = -g_0_y_xzz_yyzzzzz[k] * ab_x + g_0_y_xzz_xyyzzzzz[k];

                g_0_y_xxzz_yzzzzzz[k] = -g_0_y_xzz_yzzzzzz[k] * ab_x + g_0_y_xzz_xyzzzzzz[k];

                g_0_y_xxzz_zzzzzzz[k] = -g_0_y_xzz_zzzzzzz[k] * ab_x + g_0_y_xzz_xzzzzzzz[k];
            }

            /// Set up 756-792 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyyy_xxxxxxx, g_0_y_xyyy_xxxxxxy, g_0_y_xyyy_xxxxxxz, g_0_y_xyyy_xxxxxyy, g_0_y_xyyy_xxxxxyz, g_0_y_xyyy_xxxxxzz, g_0_y_xyyy_xxxxyyy, g_0_y_xyyy_xxxxyyz, g_0_y_xyyy_xxxxyzz, g_0_y_xyyy_xxxxzzz, g_0_y_xyyy_xxxyyyy, g_0_y_xyyy_xxxyyyz, g_0_y_xyyy_xxxyyzz, g_0_y_xyyy_xxxyzzz, g_0_y_xyyy_xxxzzzz, g_0_y_xyyy_xxyyyyy, g_0_y_xyyy_xxyyyyz, g_0_y_xyyy_xxyyyzz, g_0_y_xyyy_xxyyzzz, g_0_y_xyyy_xxyzzzz, g_0_y_xyyy_xxzzzzz, g_0_y_xyyy_xyyyyyy, g_0_y_xyyy_xyyyyyz, g_0_y_xyyy_xyyyyzz, g_0_y_xyyy_xyyyzzz, g_0_y_xyyy_xyyzzzz, g_0_y_xyyy_xyzzzzz, g_0_y_xyyy_xzzzzzz, g_0_y_xyyy_yyyyyyy, g_0_y_xyyy_yyyyyyz, g_0_y_xyyy_yyyyyzz, g_0_y_xyyy_yyyyzzz, g_0_y_xyyy_yyyzzzz, g_0_y_xyyy_yyzzzzz, g_0_y_xyyy_yzzzzzz, g_0_y_xyyy_zzzzzzz, g_0_y_yyy_xxxxxxx, g_0_y_yyy_xxxxxxxx, g_0_y_yyy_xxxxxxxy, g_0_y_yyy_xxxxxxxz, g_0_y_yyy_xxxxxxy, g_0_y_yyy_xxxxxxyy, g_0_y_yyy_xxxxxxyz, g_0_y_yyy_xxxxxxz, g_0_y_yyy_xxxxxxzz, g_0_y_yyy_xxxxxyy, g_0_y_yyy_xxxxxyyy, g_0_y_yyy_xxxxxyyz, g_0_y_yyy_xxxxxyz, g_0_y_yyy_xxxxxyzz, g_0_y_yyy_xxxxxzz, g_0_y_yyy_xxxxxzzz, g_0_y_yyy_xxxxyyy, g_0_y_yyy_xxxxyyyy, g_0_y_yyy_xxxxyyyz, g_0_y_yyy_xxxxyyz, g_0_y_yyy_xxxxyyzz, g_0_y_yyy_xxxxyzz, g_0_y_yyy_xxxxyzzz, g_0_y_yyy_xxxxzzz, g_0_y_yyy_xxxxzzzz, g_0_y_yyy_xxxyyyy, g_0_y_yyy_xxxyyyyy, g_0_y_yyy_xxxyyyyz, g_0_y_yyy_xxxyyyz, g_0_y_yyy_xxxyyyzz, g_0_y_yyy_xxxyyzz, g_0_y_yyy_xxxyyzzz, g_0_y_yyy_xxxyzzz, g_0_y_yyy_xxxyzzzz, g_0_y_yyy_xxxzzzz, g_0_y_yyy_xxxzzzzz, g_0_y_yyy_xxyyyyy, g_0_y_yyy_xxyyyyyy, g_0_y_yyy_xxyyyyyz, g_0_y_yyy_xxyyyyz, g_0_y_yyy_xxyyyyzz, g_0_y_yyy_xxyyyzz, g_0_y_yyy_xxyyyzzz, g_0_y_yyy_xxyyzzz, g_0_y_yyy_xxyyzzzz, g_0_y_yyy_xxyzzzz, g_0_y_yyy_xxyzzzzz, g_0_y_yyy_xxzzzzz, g_0_y_yyy_xxzzzzzz, g_0_y_yyy_xyyyyyy, g_0_y_yyy_xyyyyyyy, g_0_y_yyy_xyyyyyyz, g_0_y_yyy_xyyyyyz, g_0_y_yyy_xyyyyyzz, g_0_y_yyy_xyyyyzz, g_0_y_yyy_xyyyyzzz, g_0_y_yyy_xyyyzzz, g_0_y_yyy_xyyyzzzz, g_0_y_yyy_xyyzzzz, g_0_y_yyy_xyyzzzzz, g_0_y_yyy_xyzzzzz, g_0_y_yyy_xyzzzzzz, g_0_y_yyy_xzzzzzz, g_0_y_yyy_xzzzzzzz, g_0_y_yyy_yyyyyyy, g_0_y_yyy_yyyyyyz, g_0_y_yyy_yyyyyzz, g_0_y_yyy_yyyyzzz, g_0_y_yyy_yyyzzzz, g_0_y_yyy_yyzzzzz, g_0_y_yyy_yzzzzzz, g_0_y_yyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_xxxxxxx[k] = -g_0_y_yyy_xxxxxxx[k] * ab_x + g_0_y_yyy_xxxxxxxx[k];

                g_0_y_xyyy_xxxxxxy[k] = -g_0_y_yyy_xxxxxxy[k] * ab_x + g_0_y_yyy_xxxxxxxy[k];

                g_0_y_xyyy_xxxxxxz[k] = -g_0_y_yyy_xxxxxxz[k] * ab_x + g_0_y_yyy_xxxxxxxz[k];

                g_0_y_xyyy_xxxxxyy[k] = -g_0_y_yyy_xxxxxyy[k] * ab_x + g_0_y_yyy_xxxxxxyy[k];

                g_0_y_xyyy_xxxxxyz[k] = -g_0_y_yyy_xxxxxyz[k] * ab_x + g_0_y_yyy_xxxxxxyz[k];

                g_0_y_xyyy_xxxxxzz[k] = -g_0_y_yyy_xxxxxzz[k] * ab_x + g_0_y_yyy_xxxxxxzz[k];

                g_0_y_xyyy_xxxxyyy[k] = -g_0_y_yyy_xxxxyyy[k] * ab_x + g_0_y_yyy_xxxxxyyy[k];

                g_0_y_xyyy_xxxxyyz[k] = -g_0_y_yyy_xxxxyyz[k] * ab_x + g_0_y_yyy_xxxxxyyz[k];

                g_0_y_xyyy_xxxxyzz[k] = -g_0_y_yyy_xxxxyzz[k] * ab_x + g_0_y_yyy_xxxxxyzz[k];

                g_0_y_xyyy_xxxxzzz[k] = -g_0_y_yyy_xxxxzzz[k] * ab_x + g_0_y_yyy_xxxxxzzz[k];

                g_0_y_xyyy_xxxyyyy[k] = -g_0_y_yyy_xxxyyyy[k] * ab_x + g_0_y_yyy_xxxxyyyy[k];

                g_0_y_xyyy_xxxyyyz[k] = -g_0_y_yyy_xxxyyyz[k] * ab_x + g_0_y_yyy_xxxxyyyz[k];

                g_0_y_xyyy_xxxyyzz[k] = -g_0_y_yyy_xxxyyzz[k] * ab_x + g_0_y_yyy_xxxxyyzz[k];

                g_0_y_xyyy_xxxyzzz[k] = -g_0_y_yyy_xxxyzzz[k] * ab_x + g_0_y_yyy_xxxxyzzz[k];

                g_0_y_xyyy_xxxzzzz[k] = -g_0_y_yyy_xxxzzzz[k] * ab_x + g_0_y_yyy_xxxxzzzz[k];

                g_0_y_xyyy_xxyyyyy[k] = -g_0_y_yyy_xxyyyyy[k] * ab_x + g_0_y_yyy_xxxyyyyy[k];

                g_0_y_xyyy_xxyyyyz[k] = -g_0_y_yyy_xxyyyyz[k] * ab_x + g_0_y_yyy_xxxyyyyz[k];

                g_0_y_xyyy_xxyyyzz[k] = -g_0_y_yyy_xxyyyzz[k] * ab_x + g_0_y_yyy_xxxyyyzz[k];

                g_0_y_xyyy_xxyyzzz[k] = -g_0_y_yyy_xxyyzzz[k] * ab_x + g_0_y_yyy_xxxyyzzz[k];

                g_0_y_xyyy_xxyzzzz[k] = -g_0_y_yyy_xxyzzzz[k] * ab_x + g_0_y_yyy_xxxyzzzz[k];

                g_0_y_xyyy_xxzzzzz[k] = -g_0_y_yyy_xxzzzzz[k] * ab_x + g_0_y_yyy_xxxzzzzz[k];

                g_0_y_xyyy_xyyyyyy[k] = -g_0_y_yyy_xyyyyyy[k] * ab_x + g_0_y_yyy_xxyyyyyy[k];

                g_0_y_xyyy_xyyyyyz[k] = -g_0_y_yyy_xyyyyyz[k] * ab_x + g_0_y_yyy_xxyyyyyz[k];

                g_0_y_xyyy_xyyyyzz[k] = -g_0_y_yyy_xyyyyzz[k] * ab_x + g_0_y_yyy_xxyyyyzz[k];

                g_0_y_xyyy_xyyyzzz[k] = -g_0_y_yyy_xyyyzzz[k] * ab_x + g_0_y_yyy_xxyyyzzz[k];

                g_0_y_xyyy_xyyzzzz[k] = -g_0_y_yyy_xyyzzzz[k] * ab_x + g_0_y_yyy_xxyyzzzz[k];

                g_0_y_xyyy_xyzzzzz[k] = -g_0_y_yyy_xyzzzzz[k] * ab_x + g_0_y_yyy_xxyzzzzz[k];

                g_0_y_xyyy_xzzzzzz[k] = -g_0_y_yyy_xzzzzzz[k] * ab_x + g_0_y_yyy_xxzzzzzz[k];

                g_0_y_xyyy_yyyyyyy[k] = -g_0_y_yyy_yyyyyyy[k] * ab_x + g_0_y_yyy_xyyyyyyy[k];

                g_0_y_xyyy_yyyyyyz[k] = -g_0_y_yyy_yyyyyyz[k] * ab_x + g_0_y_yyy_xyyyyyyz[k];

                g_0_y_xyyy_yyyyyzz[k] = -g_0_y_yyy_yyyyyzz[k] * ab_x + g_0_y_yyy_xyyyyyzz[k];

                g_0_y_xyyy_yyyyzzz[k] = -g_0_y_yyy_yyyyzzz[k] * ab_x + g_0_y_yyy_xyyyyzzz[k];

                g_0_y_xyyy_yyyzzzz[k] = -g_0_y_yyy_yyyzzzz[k] * ab_x + g_0_y_yyy_xyyyzzzz[k];

                g_0_y_xyyy_yyzzzzz[k] = -g_0_y_yyy_yyzzzzz[k] * ab_x + g_0_y_yyy_xyyzzzzz[k];

                g_0_y_xyyy_yzzzzzz[k] = -g_0_y_yyy_yzzzzzz[k] * ab_x + g_0_y_yyy_xyzzzzzz[k];

                g_0_y_xyyy_zzzzzzz[k] = -g_0_y_yyy_zzzzzzz[k] * ab_x + g_0_y_yyy_xzzzzzzz[k];
            }

            /// Set up 792-828 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyyz_xxxxxxx, g_0_y_xyyz_xxxxxxy, g_0_y_xyyz_xxxxxxz, g_0_y_xyyz_xxxxxyy, g_0_y_xyyz_xxxxxyz, g_0_y_xyyz_xxxxxzz, g_0_y_xyyz_xxxxyyy, g_0_y_xyyz_xxxxyyz, g_0_y_xyyz_xxxxyzz, g_0_y_xyyz_xxxxzzz, g_0_y_xyyz_xxxyyyy, g_0_y_xyyz_xxxyyyz, g_0_y_xyyz_xxxyyzz, g_0_y_xyyz_xxxyzzz, g_0_y_xyyz_xxxzzzz, g_0_y_xyyz_xxyyyyy, g_0_y_xyyz_xxyyyyz, g_0_y_xyyz_xxyyyzz, g_0_y_xyyz_xxyyzzz, g_0_y_xyyz_xxyzzzz, g_0_y_xyyz_xxzzzzz, g_0_y_xyyz_xyyyyyy, g_0_y_xyyz_xyyyyyz, g_0_y_xyyz_xyyyyzz, g_0_y_xyyz_xyyyzzz, g_0_y_xyyz_xyyzzzz, g_0_y_xyyz_xyzzzzz, g_0_y_xyyz_xzzzzzz, g_0_y_xyyz_yyyyyyy, g_0_y_xyyz_yyyyyyz, g_0_y_xyyz_yyyyyzz, g_0_y_xyyz_yyyyzzz, g_0_y_xyyz_yyyzzzz, g_0_y_xyyz_yyzzzzz, g_0_y_xyyz_yzzzzzz, g_0_y_xyyz_zzzzzzz, g_0_y_yyz_xxxxxxx, g_0_y_yyz_xxxxxxxx, g_0_y_yyz_xxxxxxxy, g_0_y_yyz_xxxxxxxz, g_0_y_yyz_xxxxxxy, g_0_y_yyz_xxxxxxyy, g_0_y_yyz_xxxxxxyz, g_0_y_yyz_xxxxxxz, g_0_y_yyz_xxxxxxzz, g_0_y_yyz_xxxxxyy, g_0_y_yyz_xxxxxyyy, g_0_y_yyz_xxxxxyyz, g_0_y_yyz_xxxxxyz, g_0_y_yyz_xxxxxyzz, g_0_y_yyz_xxxxxzz, g_0_y_yyz_xxxxxzzz, g_0_y_yyz_xxxxyyy, g_0_y_yyz_xxxxyyyy, g_0_y_yyz_xxxxyyyz, g_0_y_yyz_xxxxyyz, g_0_y_yyz_xxxxyyzz, g_0_y_yyz_xxxxyzz, g_0_y_yyz_xxxxyzzz, g_0_y_yyz_xxxxzzz, g_0_y_yyz_xxxxzzzz, g_0_y_yyz_xxxyyyy, g_0_y_yyz_xxxyyyyy, g_0_y_yyz_xxxyyyyz, g_0_y_yyz_xxxyyyz, g_0_y_yyz_xxxyyyzz, g_0_y_yyz_xxxyyzz, g_0_y_yyz_xxxyyzzz, g_0_y_yyz_xxxyzzz, g_0_y_yyz_xxxyzzzz, g_0_y_yyz_xxxzzzz, g_0_y_yyz_xxxzzzzz, g_0_y_yyz_xxyyyyy, g_0_y_yyz_xxyyyyyy, g_0_y_yyz_xxyyyyyz, g_0_y_yyz_xxyyyyz, g_0_y_yyz_xxyyyyzz, g_0_y_yyz_xxyyyzz, g_0_y_yyz_xxyyyzzz, g_0_y_yyz_xxyyzzz, g_0_y_yyz_xxyyzzzz, g_0_y_yyz_xxyzzzz, g_0_y_yyz_xxyzzzzz, g_0_y_yyz_xxzzzzz, g_0_y_yyz_xxzzzzzz, g_0_y_yyz_xyyyyyy, g_0_y_yyz_xyyyyyyy, g_0_y_yyz_xyyyyyyz, g_0_y_yyz_xyyyyyz, g_0_y_yyz_xyyyyyzz, g_0_y_yyz_xyyyyzz, g_0_y_yyz_xyyyyzzz, g_0_y_yyz_xyyyzzz, g_0_y_yyz_xyyyzzzz, g_0_y_yyz_xyyzzzz, g_0_y_yyz_xyyzzzzz, g_0_y_yyz_xyzzzzz, g_0_y_yyz_xyzzzzzz, g_0_y_yyz_xzzzzzz, g_0_y_yyz_xzzzzzzz, g_0_y_yyz_yyyyyyy, g_0_y_yyz_yyyyyyz, g_0_y_yyz_yyyyyzz, g_0_y_yyz_yyyyzzz, g_0_y_yyz_yyyzzzz, g_0_y_yyz_yyzzzzz, g_0_y_yyz_yzzzzzz, g_0_y_yyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_xxxxxxx[k] = -g_0_y_yyz_xxxxxxx[k] * ab_x + g_0_y_yyz_xxxxxxxx[k];

                g_0_y_xyyz_xxxxxxy[k] = -g_0_y_yyz_xxxxxxy[k] * ab_x + g_0_y_yyz_xxxxxxxy[k];

                g_0_y_xyyz_xxxxxxz[k] = -g_0_y_yyz_xxxxxxz[k] * ab_x + g_0_y_yyz_xxxxxxxz[k];

                g_0_y_xyyz_xxxxxyy[k] = -g_0_y_yyz_xxxxxyy[k] * ab_x + g_0_y_yyz_xxxxxxyy[k];

                g_0_y_xyyz_xxxxxyz[k] = -g_0_y_yyz_xxxxxyz[k] * ab_x + g_0_y_yyz_xxxxxxyz[k];

                g_0_y_xyyz_xxxxxzz[k] = -g_0_y_yyz_xxxxxzz[k] * ab_x + g_0_y_yyz_xxxxxxzz[k];

                g_0_y_xyyz_xxxxyyy[k] = -g_0_y_yyz_xxxxyyy[k] * ab_x + g_0_y_yyz_xxxxxyyy[k];

                g_0_y_xyyz_xxxxyyz[k] = -g_0_y_yyz_xxxxyyz[k] * ab_x + g_0_y_yyz_xxxxxyyz[k];

                g_0_y_xyyz_xxxxyzz[k] = -g_0_y_yyz_xxxxyzz[k] * ab_x + g_0_y_yyz_xxxxxyzz[k];

                g_0_y_xyyz_xxxxzzz[k] = -g_0_y_yyz_xxxxzzz[k] * ab_x + g_0_y_yyz_xxxxxzzz[k];

                g_0_y_xyyz_xxxyyyy[k] = -g_0_y_yyz_xxxyyyy[k] * ab_x + g_0_y_yyz_xxxxyyyy[k];

                g_0_y_xyyz_xxxyyyz[k] = -g_0_y_yyz_xxxyyyz[k] * ab_x + g_0_y_yyz_xxxxyyyz[k];

                g_0_y_xyyz_xxxyyzz[k] = -g_0_y_yyz_xxxyyzz[k] * ab_x + g_0_y_yyz_xxxxyyzz[k];

                g_0_y_xyyz_xxxyzzz[k] = -g_0_y_yyz_xxxyzzz[k] * ab_x + g_0_y_yyz_xxxxyzzz[k];

                g_0_y_xyyz_xxxzzzz[k] = -g_0_y_yyz_xxxzzzz[k] * ab_x + g_0_y_yyz_xxxxzzzz[k];

                g_0_y_xyyz_xxyyyyy[k] = -g_0_y_yyz_xxyyyyy[k] * ab_x + g_0_y_yyz_xxxyyyyy[k];

                g_0_y_xyyz_xxyyyyz[k] = -g_0_y_yyz_xxyyyyz[k] * ab_x + g_0_y_yyz_xxxyyyyz[k];

                g_0_y_xyyz_xxyyyzz[k] = -g_0_y_yyz_xxyyyzz[k] * ab_x + g_0_y_yyz_xxxyyyzz[k];

                g_0_y_xyyz_xxyyzzz[k] = -g_0_y_yyz_xxyyzzz[k] * ab_x + g_0_y_yyz_xxxyyzzz[k];

                g_0_y_xyyz_xxyzzzz[k] = -g_0_y_yyz_xxyzzzz[k] * ab_x + g_0_y_yyz_xxxyzzzz[k];

                g_0_y_xyyz_xxzzzzz[k] = -g_0_y_yyz_xxzzzzz[k] * ab_x + g_0_y_yyz_xxxzzzzz[k];

                g_0_y_xyyz_xyyyyyy[k] = -g_0_y_yyz_xyyyyyy[k] * ab_x + g_0_y_yyz_xxyyyyyy[k];

                g_0_y_xyyz_xyyyyyz[k] = -g_0_y_yyz_xyyyyyz[k] * ab_x + g_0_y_yyz_xxyyyyyz[k];

                g_0_y_xyyz_xyyyyzz[k] = -g_0_y_yyz_xyyyyzz[k] * ab_x + g_0_y_yyz_xxyyyyzz[k];

                g_0_y_xyyz_xyyyzzz[k] = -g_0_y_yyz_xyyyzzz[k] * ab_x + g_0_y_yyz_xxyyyzzz[k];

                g_0_y_xyyz_xyyzzzz[k] = -g_0_y_yyz_xyyzzzz[k] * ab_x + g_0_y_yyz_xxyyzzzz[k];

                g_0_y_xyyz_xyzzzzz[k] = -g_0_y_yyz_xyzzzzz[k] * ab_x + g_0_y_yyz_xxyzzzzz[k];

                g_0_y_xyyz_xzzzzzz[k] = -g_0_y_yyz_xzzzzzz[k] * ab_x + g_0_y_yyz_xxzzzzzz[k];

                g_0_y_xyyz_yyyyyyy[k] = -g_0_y_yyz_yyyyyyy[k] * ab_x + g_0_y_yyz_xyyyyyyy[k];

                g_0_y_xyyz_yyyyyyz[k] = -g_0_y_yyz_yyyyyyz[k] * ab_x + g_0_y_yyz_xyyyyyyz[k];

                g_0_y_xyyz_yyyyyzz[k] = -g_0_y_yyz_yyyyyzz[k] * ab_x + g_0_y_yyz_xyyyyyzz[k];

                g_0_y_xyyz_yyyyzzz[k] = -g_0_y_yyz_yyyyzzz[k] * ab_x + g_0_y_yyz_xyyyyzzz[k];

                g_0_y_xyyz_yyyzzzz[k] = -g_0_y_yyz_yyyzzzz[k] * ab_x + g_0_y_yyz_xyyyzzzz[k];

                g_0_y_xyyz_yyzzzzz[k] = -g_0_y_yyz_yyzzzzz[k] * ab_x + g_0_y_yyz_xyyzzzzz[k];

                g_0_y_xyyz_yzzzzzz[k] = -g_0_y_yyz_yzzzzzz[k] * ab_x + g_0_y_yyz_xyzzzzzz[k];

                g_0_y_xyyz_zzzzzzz[k] = -g_0_y_yyz_zzzzzzz[k] * ab_x + g_0_y_yyz_xzzzzzzz[k];
            }

            /// Set up 828-864 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyzz_xxxxxxx, g_0_y_xyzz_xxxxxxy, g_0_y_xyzz_xxxxxxz, g_0_y_xyzz_xxxxxyy, g_0_y_xyzz_xxxxxyz, g_0_y_xyzz_xxxxxzz, g_0_y_xyzz_xxxxyyy, g_0_y_xyzz_xxxxyyz, g_0_y_xyzz_xxxxyzz, g_0_y_xyzz_xxxxzzz, g_0_y_xyzz_xxxyyyy, g_0_y_xyzz_xxxyyyz, g_0_y_xyzz_xxxyyzz, g_0_y_xyzz_xxxyzzz, g_0_y_xyzz_xxxzzzz, g_0_y_xyzz_xxyyyyy, g_0_y_xyzz_xxyyyyz, g_0_y_xyzz_xxyyyzz, g_0_y_xyzz_xxyyzzz, g_0_y_xyzz_xxyzzzz, g_0_y_xyzz_xxzzzzz, g_0_y_xyzz_xyyyyyy, g_0_y_xyzz_xyyyyyz, g_0_y_xyzz_xyyyyzz, g_0_y_xyzz_xyyyzzz, g_0_y_xyzz_xyyzzzz, g_0_y_xyzz_xyzzzzz, g_0_y_xyzz_xzzzzzz, g_0_y_xyzz_yyyyyyy, g_0_y_xyzz_yyyyyyz, g_0_y_xyzz_yyyyyzz, g_0_y_xyzz_yyyyzzz, g_0_y_xyzz_yyyzzzz, g_0_y_xyzz_yyzzzzz, g_0_y_xyzz_yzzzzzz, g_0_y_xyzz_zzzzzzz, g_0_y_yzz_xxxxxxx, g_0_y_yzz_xxxxxxxx, g_0_y_yzz_xxxxxxxy, g_0_y_yzz_xxxxxxxz, g_0_y_yzz_xxxxxxy, g_0_y_yzz_xxxxxxyy, g_0_y_yzz_xxxxxxyz, g_0_y_yzz_xxxxxxz, g_0_y_yzz_xxxxxxzz, g_0_y_yzz_xxxxxyy, g_0_y_yzz_xxxxxyyy, g_0_y_yzz_xxxxxyyz, g_0_y_yzz_xxxxxyz, g_0_y_yzz_xxxxxyzz, g_0_y_yzz_xxxxxzz, g_0_y_yzz_xxxxxzzz, g_0_y_yzz_xxxxyyy, g_0_y_yzz_xxxxyyyy, g_0_y_yzz_xxxxyyyz, g_0_y_yzz_xxxxyyz, g_0_y_yzz_xxxxyyzz, g_0_y_yzz_xxxxyzz, g_0_y_yzz_xxxxyzzz, g_0_y_yzz_xxxxzzz, g_0_y_yzz_xxxxzzzz, g_0_y_yzz_xxxyyyy, g_0_y_yzz_xxxyyyyy, g_0_y_yzz_xxxyyyyz, g_0_y_yzz_xxxyyyz, g_0_y_yzz_xxxyyyzz, g_0_y_yzz_xxxyyzz, g_0_y_yzz_xxxyyzzz, g_0_y_yzz_xxxyzzz, g_0_y_yzz_xxxyzzzz, g_0_y_yzz_xxxzzzz, g_0_y_yzz_xxxzzzzz, g_0_y_yzz_xxyyyyy, g_0_y_yzz_xxyyyyyy, g_0_y_yzz_xxyyyyyz, g_0_y_yzz_xxyyyyz, g_0_y_yzz_xxyyyyzz, g_0_y_yzz_xxyyyzz, g_0_y_yzz_xxyyyzzz, g_0_y_yzz_xxyyzzz, g_0_y_yzz_xxyyzzzz, g_0_y_yzz_xxyzzzz, g_0_y_yzz_xxyzzzzz, g_0_y_yzz_xxzzzzz, g_0_y_yzz_xxzzzzzz, g_0_y_yzz_xyyyyyy, g_0_y_yzz_xyyyyyyy, g_0_y_yzz_xyyyyyyz, g_0_y_yzz_xyyyyyz, g_0_y_yzz_xyyyyyzz, g_0_y_yzz_xyyyyzz, g_0_y_yzz_xyyyyzzz, g_0_y_yzz_xyyyzzz, g_0_y_yzz_xyyyzzzz, g_0_y_yzz_xyyzzzz, g_0_y_yzz_xyyzzzzz, g_0_y_yzz_xyzzzzz, g_0_y_yzz_xyzzzzzz, g_0_y_yzz_xzzzzzz, g_0_y_yzz_xzzzzzzz, g_0_y_yzz_yyyyyyy, g_0_y_yzz_yyyyyyz, g_0_y_yzz_yyyyyzz, g_0_y_yzz_yyyyzzz, g_0_y_yzz_yyyzzzz, g_0_y_yzz_yyzzzzz, g_0_y_yzz_yzzzzzz, g_0_y_yzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_xxxxxxx[k] = -g_0_y_yzz_xxxxxxx[k] * ab_x + g_0_y_yzz_xxxxxxxx[k];

                g_0_y_xyzz_xxxxxxy[k] = -g_0_y_yzz_xxxxxxy[k] * ab_x + g_0_y_yzz_xxxxxxxy[k];

                g_0_y_xyzz_xxxxxxz[k] = -g_0_y_yzz_xxxxxxz[k] * ab_x + g_0_y_yzz_xxxxxxxz[k];

                g_0_y_xyzz_xxxxxyy[k] = -g_0_y_yzz_xxxxxyy[k] * ab_x + g_0_y_yzz_xxxxxxyy[k];

                g_0_y_xyzz_xxxxxyz[k] = -g_0_y_yzz_xxxxxyz[k] * ab_x + g_0_y_yzz_xxxxxxyz[k];

                g_0_y_xyzz_xxxxxzz[k] = -g_0_y_yzz_xxxxxzz[k] * ab_x + g_0_y_yzz_xxxxxxzz[k];

                g_0_y_xyzz_xxxxyyy[k] = -g_0_y_yzz_xxxxyyy[k] * ab_x + g_0_y_yzz_xxxxxyyy[k];

                g_0_y_xyzz_xxxxyyz[k] = -g_0_y_yzz_xxxxyyz[k] * ab_x + g_0_y_yzz_xxxxxyyz[k];

                g_0_y_xyzz_xxxxyzz[k] = -g_0_y_yzz_xxxxyzz[k] * ab_x + g_0_y_yzz_xxxxxyzz[k];

                g_0_y_xyzz_xxxxzzz[k] = -g_0_y_yzz_xxxxzzz[k] * ab_x + g_0_y_yzz_xxxxxzzz[k];

                g_0_y_xyzz_xxxyyyy[k] = -g_0_y_yzz_xxxyyyy[k] * ab_x + g_0_y_yzz_xxxxyyyy[k];

                g_0_y_xyzz_xxxyyyz[k] = -g_0_y_yzz_xxxyyyz[k] * ab_x + g_0_y_yzz_xxxxyyyz[k];

                g_0_y_xyzz_xxxyyzz[k] = -g_0_y_yzz_xxxyyzz[k] * ab_x + g_0_y_yzz_xxxxyyzz[k];

                g_0_y_xyzz_xxxyzzz[k] = -g_0_y_yzz_xxxyzzz[k] * ab_x + g_0_y_yzz_xxxxyzzz[k];

                g_0_y_xyzz_xxxzzzz[k] = -g_0_y_yzz_xxxzzzz[k] * ab_x + g_0_y_yzz_xxxxzzzz[k];

                g_0_y_xyzz_xxyyyyy[k] = -g_0_y_yzz_xxyyyyy[k] * ab_x + g_0_y_yzz_xxxyyyyy[k];

                g_0_y_xyzz_xxyyyyz[k] = -g_0_y_yzz_xxyyyyz[k] * ab_x + g_0_y_yzz_xxxyyyyz[k];

                g_0_y_xyzz_xxyyyzz[k] = -g_0_y_yzz_xxyyyzz[k] * ab_x + g_0_y_yzz_xxxyyyzz[k];

                g_0_y_xyzz_xxyyzzz[k] = -g_0_y_yzz_xxyyzzz[k] * ab_x + g_0_y_yzz_xxxyyzzz[k];

                g_0_y_xyzz_xxyzzzz[k] = -g_0_y_yzz_xxyzzzz[k] * ab_x + g_0_y_yzz_xxxyzzzz[k];

                g_0_y_xyzz_xxzzzzz[k] = -g_0_y_yzz_xxzzzzz[k] * ab_x + g_0_y_yzz_xxxzzzzz[k];

                g_0_y_xyzz_xyyyyyy[k] = -g_0_y_yzz_xyyyyyy[k] * ab_x + g_0_y_yzz_xxyyyyyy[k];

                g_0_y_xyzz_xyyyyyz[k] = -g_0_y_yzz_xyyyyyz[k] * ab_x + g_0_y_yzz_xxyyyyyz[k];

                g_0_y_xyzz_xyyyyzz[k] = -g_0_y_yzz_xyyyyzz[k] * ab_x + g_0_y_yzz_xxyyyyzz[k];

                g_0_y_xyzz_xyyyzzz[k] = -g_0_y_yzz_xyyyzzz[k] * ab_x + g_0_y_yzz_xxyyyzzz[k];

                g_0_y_xyzz_xyyzzzz[k] = -g_0_y_yzz_xyyzzzz[k] * ab_x + g_0_y_yzz_xxyyzzzz[k];

                g_0_y_xyzz_xyzzzzz[k] = -g_0_y_yzz_xyzzzzz[k] * ab_x + g_0_y_yzz_xxyzzzzz[k];

                g_0_y_xyzz_xzzzzzz[k] = -g_0_y_yzz_xzzzzzz[k] * ab_x + g_0_y_yzz_xxzzzzzz[k];

                g_0_y_xyzz_yyyyyyy[k] = -g_0_y_yzz_yyyyyyy[k] * ab_x + g_0_y_yzz_xyyyyyyy[k];

                g_0_y_xyzz_yyyyyyz[k] = -g_0_y_yzz_yyyyyyz[k] * ab_x + g_0_y_yzz_xyyyyyyz[k];

                g_0_y_xyzz_yyyyyzz[k] = -g_0_y_yzz_yyyyyzz[k] * ab_x + g_0_y_yzz_xyyyyyzz[k];

                g_0_y_xyzz_yyyyzzz[k] = -g_0_y_yzz_yyyyzzz[k] * ab_x + g_0_y_yzz_xyyyyzzz[k];

                g_0_y_xyzz_yyyzzzz[k] = -g_0_y_yzz_yyyzzzz[k] * ab_x + g_0_y_yzz_xyyyzzzz[k];

                g_0_y_xyzz_yyzzzzz[k] = -g_0_y_yzz_yyzzzzz[k] * ab_x + g_0_y_yzz_xyyzzzzz[k];

                g_0_y_xyzz_yzzzzzz[k] = -g_0_y_yzz_yzzzzzz[k] * ab_x + g_0_y_yzz_xyzzzzzz[k];

                g_0_y_xyzz_zzzzzzz[k] = -g_0_y_yzz_zzzzzzz[k] * ab_x + g_0_y_yzz_xzzzzzzz[k];
            }

            /// Set up 864-900 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xzzz_xxxxxxx, g_0_y_xzzz_xxxxxxy, g_0_y_xzzz_xxxxxxz, g_0_y_xzzz_xxxxxyy, g_0_y_xzzz_xxxxxyz, g_0_y_xzzz_xxxxxzz, g_0_y_xzzz_xxxxyyy, g_0_y_xzzz_xxxxyyz, g_0_y_xzzz_xxxxyzz, g_0_y_xzzz_xxxxzzz, g_0_y_xzzz_xxxyyyy, g_0_y_xzzz_xxxyyyz, g_0_y_xzzz_xxxyyzz, g_0_y_xzzz_xxxyzzz, g_0_y_xzzz_xxxzzzz, g_0_y_xzzz_xxyyyyy, g_0_y_xzzz_xxyyyyz, g_0_y_xzzz_xxyyyzz, g_0_y_xzzz_xxyyzzz, g_0_y_xzzz_xxyzzzz, g_0_y_xzzz_xxzzzzz, g_0_y_xzzz_xyyyyyy, g_0_y_xzzz_xyyyyyz, g_0_y_xzzz_xyyyyzz, g_0_y_xzzz_xyyyzzz, g_0_y_xzzz_xyyzzzz, g_0_y_xzzz_xyzzzzz, g_0_y_xzzz_xzzzzzz, g_0_y_xzzz_yyyyyyy, g_0_y_xzzz_yyyyyyz, g_0_y_xzzz_yyyyyzz, g_0_y_xzzz_yyyyzzz, g_0_y_xzzz_yyyzzzz, g_0_y_xzzz_yyzzzzz, g_0_y_xzzz_yzzzzzz, g_0_y_xzzz_zzzzzzz, g_0_y_zzz_xxxxxxx, g_0_y_zzz_xxxxxxxx, g_0_y_zzz_xxxxxxxy, g_0_y_zzz_xxxxxxxz, g_0_y_zzz_xxxxxxy, g_0_y_zzz_xxxxxxyy, g_0_y_zzz_xxxxxxyz, g_0_y_zzz_xxxxxxz, g_0_y_zzz_xxxxxxzz, g_0_y_zzz_xxxxxyy, g_0_y_zzz_xxxxxyyy, g_0_y_zzz_xxxxxyyz, g_0_y_zzz_xxxxxyz, g_0_y_zzz_xxxxxyzz, g_0_y_zzz_xxxxxzz, g_0_y_zzz_xxxxxzzz, g_0_y_zzz_xxxxyyy, g_0_y_zzz_xxxxyyyy, g_0_y_zzz_xxxxyyyz, g_0_y_zzz_xxxxyyz, g_0_y_zzz_xxxxyyzz, g_0_y_zzz_xxxxyzz, g_0_y_zzz_xxxxyzzz, g_0_y_zzz_xxxxzzz, g_0_y_zzz_xxxxzzzz, g_0_y_zzz_xxxyyyy, g_0_y_zzz_xxxyyyyy, g_0_y_zzz_xxxyyyyz, g_0_y_zzz_xxxyyyz, g_0_y_zzz_xxxyyyzz, g_0_y_zzz_xxxyyzz, g_0_y_zzz_xxxyyzzz, g_0_y_zzz_xxxyzzz, g_0_y_zzz_xxxyzzzz, g_0_y_zzz_xxxzzzz, g_0_y_zzz_xxxzzzzz, g_0_y_zzz_xxyyyyy, g_0_y_zzz_xxyyyyyy, g_0_y_zzz_xxyyyyyz, g_0_y_zzz_xxyyyyz, g_0_y_zzz_xxyyyyzz, g_0_y_zzz_xxyyyzz, g_0_y_zzz_xxyyyzzz, g_0_y_zzz_xxyyzzz, g_0_y_zzz_xxyyzzzz, g_0_y_zzz_xxyzzzz, g_0_y_zzz_xxyzzzzz, g_0_y_zzz_xxzzzzz, g_0_y_zzz_xxzzzzzz, g_0_y_zzz_xyyyyyy, g_0_y_zzz_xyyyyyyy, g_0_y_zzz_xyyyyyyz, g_0_y_zzz_xyyyyyz, g_0_y_zzz_xyyyyyzz, g_0_y_zzz_xyyyyzz, g_0_y_zzz_xyyyyzzz, g_0_y_zzz_xyyyzzz, g_0_y_zzz_xyyyzzzz, g_0_y_zzz_xyyzzzz, g_0_y_zzz_xyyzzzzz, g_0_y_zzz_xyzzzzz, g_0_y_zzz_xyzzzzzz, g_0_y_zzz_xzzzzzz, g_0_y_zzz_xzzzzzzz, g_0_y_zzz_yyyyyyy, g_0_y_zzz_yyyyyyz, g_0_y_zzz_yyyyyzz, g_0_y_zzz_yyyyzzz, g_0_y_zzz_yyyzzzz, g_0_y_zzz_yyzzzzz, g_0_y_zzz_yzzzzzz, g_0_y_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_xxxxxxx[k] = -g_0_y_zzz_xxxxxxx[k] * ab_x + g_0_y_zzz_xxxxxxxx[k];

                g_0_y_xzzz_xxxxxxy[k] = -g_0_y_zzz_xxxxxxy[k] * ab_x + g_0_y_zzz_xxxxxxxy[k];

                g_0_y_xzzz_xxxxxxz[k] = -g_0_y_zzz_xxxxxxz[k] * ab_x + g_0_y_zzz_xxxxxxxz[k];

                g_0_y_xzzz_xxxxxyy[k] = -g_0_y_zzz_xxxxxyy[k] * ab_x + g_0_y_zzz_xxxxxxyy[k];

                g_0_y_xzzz_xxxxxyz[k] = -g_0_y_zzz_xxxxxyz[k] * ab_x + g_0_y_zzz_xxxxxxyz[k];

                g_0_y_xzzz_xxxxxzz[k] = -g_0_y_zzz_xxxxxzz[k] * ab_x + g_0_y_zzz_xxxxxxzz[k];

                g_0_y_xzzz_xxxxyyy[k] = -g_0_y_zzz_xxxxyyy[k] * ab_x + g_0_y_zzz_xxxxxyyy[k];

                g_0_y_xzzz_xxxxyyz[k] = -g_0_y_zzz_xxxxyyz[k] * ab_x + g_0_y_zzz_xxxxxyyz[k];

                g_0_y_xzzz_xxxxyzz[k] = -g_0_y_zzz_xxxxyzz[k] * ab_x + g_0_y_zzz_xxxxxyzz[k];

                g_0_y_xzzz_xxxxzzz[k] = -g_0_y_zzz_xxxxzzz[k] * ab_x + g_0_y_zzz_xxxxxzzz[k];

                g_0_y_xzzz_xxxyyyy[k] = -g_0_y_zzz_xxxyyyy[k] * ab_x + g_0_y_zzz_xxxxyyyy[k];

                g_0_y_xzzz_xxxyyyz[k] = -g_0_y_zzz_xxxyyyz[k] * ab_x + g_0_y_zzz_xxxxyyyz[k];

                g_0_y_xzzz_xxxyyzz[k] = -g_0_y_zzz_xxxyyzz[k] * ab_x + g_0_y_zzz_xxxxyyzz[k];

                g_0_y_xzzz_xxxyzzz[k] = -g_0_y_zzz_xxxyzzz[k] * ab_x + g_0_y_zzz_xxxxyzzz[k];

                g_0_y_xzzz_xxxzzzz[k] = -g_0_y_zzz_xxxzzzz[k] * ab_x + g_0_y_zzz_xxxxzzzz[k];

                g_0_y_xzzz_xxyyyyy[k] = -g_0_y_zzz_xxyyyyy[k] * ab_x + g_0_y_zzz_xxxyyyyy[k];

                g_0_y_xzzz_xxyyyyz[k] = -g_0_y_zzz_xxyyyyz[k] * ab_x + g_0_y_zzz_xxxyyyyz[k];

                g_0_y_xzzz_xxyyyzz[k] = -g_0_y_zzz_xxyyyzz[k] * ab_x + g_0_y_zzz_xxxyyyzz[k];

                g_0_y_xzzz_xxyyzzz[k] = -g_0_y_zzz_xxyyzzz[k] * ab_x + g_0_y_zzz_xxxyyzzz[k];

                g_0_y_xzzz_xxyzzzz[k] = -g_0_y_zzz_xxyzzzz[k] * ab_x + g_0_y_zzz_xxxyzzzz[k];

                g_0_y_xzzz_xxzzzzz[k] = -g_0_y_zzz_xxzzzzz[k] * ab_x + g_0_y_zzz_xxxzzzzz[k];

                g_0_y_xzzz_xyyyyyy[k] = -g_0_y_zzz_xyyyyyy[k] * ab_x + g_0_y_zzz_xxyyyyyy[k];

                g_0_y_xzzz_xyyyyyz[k] = -g_0_y_zzz_xyyyyyz[k] * ab_x + g_0_y_zzz_xxyyyyyz[k];

                g_0_y_xzzz_xyyyyzz[k] = -g_0_y_zzz_xyyyyzz[k] * ab_x + g_0_y_zzz_xxyyyyzz[k];

                g_0_y_xzzz_xyyyzzz[k] = -g_0_y_zzz_xyyyzzz[k] * ab_x + g_0_y_zzz_xxyyyzzz[k];

                g_0_y_xzzz_xyyzzzz[k] = -g_0_y_zzz_xyyzzzz[k] * ab_x + g_0_y_zzz_xxyyzzzz[k];

                g_0_y_xzzz_xyzzzzz[k] = -g_0_y_zzz_xyzzzzz[k] * ab_x + g_0_y_zzz_xxyzzzzz[k];

                g_0_y_xzzz_xzzzzzz[k] = -g_0_y_zzz_xzzzzzz[k] * ab_x + g_0_y_zzz_xxzzzzzz[k];

                g_0_y_xzzz_yyyyyyy[k] = -g_0_y_zzz_yyyyyyy[k] * ab_x + g_0_y_zzz_xyyyyyyy[k];

                g_0_y_xzzz_yyyyyyz[k] = -g_0_y_zzz_yyyyyyz[k] * ab_x + g_0_y_zzz_xyyyyyyz[k];

                g_0_y_xzzz_yyyyyzz[k] = -g_0_y_zzz_yyyyyzz[k] * ab_x + g_0_y_zzz_xyyyyyzz[k];

                g_0_y_xzzz_yyyyzzz[k] = -g_0_y_zzz_yyyyzzz[k] * ab_x + g_0_y_zzz_xyyyyzzz[k];

                g_0_y_xzzz_yyyzzzz[k] = -g_0_y_zzz_yyyzzzz[k] * ab_x + g_0_y_zzz_xyyyzzzz[k];

                g_0_y_xzzz_yyzzzzz[k] = -g_0_y_zzz_yyzzzzz[k] * ab_x + g_0_y_zzz_xyyzzzzz[k];

                g_0_y_xzzz_yzzzzzz[k] = -g_0_y_zzz_yzzzzzz[k] * ab_x + g_0_y_zzz_xyzzzzzz[k];

                g_0_y_xzzz_zzzzzzz[k] = -g_0_y_zzz_zzzzzzz[k] * ab_x + g_0_y_zzz_xzzzzzzz[k];
            }

            /// Set up 900-936 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyy_xxxxxxx, g_0_y_yyy_xxxxxxxy, g_0_y_yyy_xxxxxxy, g_0_y_yyy_xxxxxxyy, g_0_y_yyy_xxxxxxyz, g_0_y_yyy_xxxxxxz, g_0_y_yyy_xxxxxyy, g_0_y_yyy_xxxxxyyy, g_0_y_yyy_xxxxxyyz, g_0_y_yyy_xxxxxyz, g_0_y_yyy_xxxxxyzz, g_0_y_yyy_xxxxxzz, g_0_y_yyy_xxxxyyy, g_0_y_yyy_xxxxyyyy, g_0_y_yyy_xxxxyyyz, g_0_y_yyy_xxxxyyz, g_0_y_yyy_xxxxyyzz, g_0_y_yyy_xxxxyzz, g_0_y_yyy_xxxxyzzz, g_0_y_yyy_xxxxzzz, g_0_y_yyy_xxxyyyy, g_0_y_yyy_xxxyyyyy, g_0_y_yyy_xxxyyyyz, g_0_y_yyy_xxxyyyz, g_0_y_yyy_xxxyyyzz, g_0_y_yyy_xxxyyzz, g_0_y_yyy_xxxyyzzz, g_0_y_yyy_xxxyzzz, g_0_y_yyy_xxxyzzzz, g_0_y_yyy_xxxzzzz, g_0_y_yyy_xxyyyyy, g_0_y_yyy_xxyyyyyy, g_0_y_yyy_xxyyyyyz, g_0_y_yyy_xxyyyyz, g_0_y_yyy_xxyyyyzz, g_0_y_yyy_xxyyyzz, g_0_y_yyy_xxyyyzzz, g_0_y_yyy_xxyyzzz, g_0_y_yyy_xxyyzzzz, g_0_y_yyy_xxyzzzz, g_0_y_yyy_xxyzzzzz, g_0_y_yyy_xxzzzzz, g_0_y_yyy_xyyyyyy, g_0_y_yyy_xyyyyyyy, g_0_y_yyy_xyyyyyyz, g_0_y_yyy_xyyyyyz, g_0_y_yyy_xyyyyyzz, g_0_y_yyy_xyyyyzz, g_0_y_yyy_xyyyyzzz, g_0_y_yyy_xyyyzzz, g_0_y_yyy_xyyyzzzz, g_0_y_yyy_xyyzzzz, g_0_y_yyy_xyyzzzzz, g_0_y_yyy_xyzzzzz, g_0_y_yyy_xyzzzzzz, g_0_y_yyy_xzzzzzz, g_0_y_yyy_yyyyyyy, g_0_y_yyy_yyyyyyyy, g_0_y_yyy_yyyyyyyz, g_0_y_yyy_yyyyyyz, g_0_y_yyy_yyyyyyzz, g_0_y_yyy_yyyyyzz, g_0_y_yyy_yyyyyzzz, g_0_y_yyy_yyyyzzz, g_0_y_yyy_yyyyzzzz, g_0_y_yyy_yyyzzzz, g_0_y_yyy_yyyzzzzz, g_0_y_yyy_yyzzzzz, g_0_y_yyy_yyzzzzzz, g_0_y_yyy_yzzzzzz, g_0_y_yyy_yzzzzzzz, g_0_y_yyy_zzzzzzz, g_0_y_yyyy_xxxxxxx, g_0_y_yyyy_xxxxxxy, g_0_y_yyyy_xxxxxxz, g_0_y_yyyy_xxxxxyy, g_0_y_yyyy_xxxxxyz, g_0_y_yyyy_xxxxxzz, g_0_y_yyyy_xxxxyyy, g_0_y_yyyy_xxxxyyz, g_0_y_yyyy_xxxxyzz, g_0_y_yyyy_xxxxzzz, g_0_y_yyyy_xxxyyyy, g_0_y_yyyy_xxxyyyz, g_0_y_yyyy_xxxyyzz, g_0_y_yyyy_xxxyzzz, g_0_y_yyyy_xxxzzzz, g_0_y_yyyy_xxyyyyy, g_0_y_yyyy_xxyyyyz, g_0_y_yyyy_xxyyyzz, g_0_y_yyyy_xxyyzzz, g_0_y_yyyy_xxyzzzz, g_0_y_yyyy_xxzzzzz, g_0_y_yyyy_xyyyyyy, g_0_y_yyyy_xyyyyyz, g_0_y_yyyy_xyyyyzz, g_0_y_yyyy_xyyyzzz, g_0_y_yyyy_xyyzzzz, g_0_y_yyyy_xyzzzzz, g_0_y_yyyy_xzzzzzz, g_0_y_yyyy_yyyyyyy, g_0_y_yyyy_yyyyyyz, g_0_y_yyyy_yyyyyzz, g_0_y_yyyy_yyyyzzz, g_0_y_yyyy_yyyzzzz, g_0_y_yyyy_yyzzzzz, g_0_y_yyyy_yzzzzzz, g_0_y_yyyy_zzzzzzz, g_yyy_xxxxxxx, g_yyy_xxxxxxy, g_yyy_xxxxxxz, g_yyy_xxxxxyy, g_yyy_xxxxxyz, g_yyy_xxxxxzz, g_yyy_xxxxyyy, g_yyy_xxxxyyz, g_yyy_xxxxyzz, g_yyy_xxxxzzz, g_yyy_xxxyyyy, g_yyy_xxxyyyz, g_yyy_xxxyyzz, g_yyy_xxxyzzz, g_yyy_xxxzzzz, g_yyy_xxyyyyy, g_yyy_xxyyyyz, g_yyy_xxyyyzz, g_yyy_xxyyzzz, g_yyy_xxyzzzz, g_yyy_xxzzzzz, g_yyy_xyyyyyy, g_yyy_xyyyyyz, g_yyy_xyyyyzz, g_yyy_xyyyzzz, g_yyy_xyyzzzz, g_yyy_xyzzzzz, g_yyy_xzzzzzz, g_yyy_yyyyyyy, g_yyy_yyyyyyz, g_yyy_yyyyyzz, g_yyy_yyyyzzz, g_yyy_yyyzzzz, g_yyy_yyzzzzz, g_yyy_yzzzzzz, g_yyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_xxxxxxx[k] = g_yyy_xxxxxxx[k] - g_0_y_yyy_xxxxxxx[k] * ab_y + g_0_y_yyy_xxxxxxxy[k];

                g_0_y_yyyy_xxxxxxy[k] = g_yyy_xxxxxxy[k] - g_0_y_yyy_xxxxxxy[k] * ab_y + g_0_y_yyy_xxxxxxyy[k];

                g_0_y_yyyy_xxxxxxz[k] = g_yyy_xxxxxxz[k] - g_0_y_yyy_xxxxxxz[k] * ab_y + g_0_y_yyy_xxxxxxyz[k];

                g_0_y_yyyy_xxxxxyy[k] = g_yyy_xxxxxyy[k] - g_0_y_yyy_xxxxxyy[k] * ab_y + g_0_y_yyy_xxxxxyyy[k];

                g_0_y_yyyy_xxxxxyz[k] = g_yyy_xxxxxyz[k] - g_0_y_yyy_xxxxxyz[k] * ab_y + g_0_y_yyy_xxxxxyyz[k];

                g_0_y_yyyy_xxxxxzz[k] = g_yyy_xxxxxzz[k] - g_0_y_yyy_xxxxxzz[k] * ab_y + g_0_y_yyy_xxxxxyzz[k];

                g_0_y_yyyy_xxxxyyy[k] = g_yyy_xxxxyyy[k] - g_0_y_yyy_xxxxyyy[k] * ab_y + g_0_y_yyy_xxxxyyyy[k];

                g_0_y_yyyy_xxxxyyz[k] = g_yyy_xxxxyyz[k] - g_0_y_yyy_xxxxyyz[k] * ab_y + g_0_y_yyy_xxxxyyyz[k];

                g_0_y_yyyy_xxxxyzz[k] = g_yyy_xxxxyzz[k] - g_0_y_yyy_xxxxyzz[k] * ab_y + g_0_y_yyy_xxxxyyzz[k];

                g_0_y_yyyy_xxxxzzz[k] = g_yyy_xxxxzzz[k] - g_0_y_yyy_xxxxzzz[k] * ab_y + g_0_y_yyy_xxxxyzzz[k];

                g_0_y_yyyy_xxxyyyy[k] = g_yyy_xxxyyyy[k] - g_0_y_yyy_xxxyyyy[k] * ab_y + g_0_y_yyy_xxxyyyyy[k];

                g_0_y_yyyy_xxxyyyz[k] = g_yyy_xxxyyyz[k] - g_0_y_yyy_xxxyyyz[k] * ab_y + g_0_y_yyy_xxxyyyyz[k];

                g_0_y_yyyy_xxxyyzz[k] = g_yyy_xxxyyzz[k] - g_0_y_yyy_xxxyyzz[k] * ab_y + g_0_y_yyy_xxxyyyzz[k];

                g_0_y_yyyy_xxxyzzz[k] = g_yyy_xxxyzzz[k] - g_0_y_yyy_xxxyzzz[k] * ab_y + g_0_y_yyy_xxxyyzzz[k];

                g_0_y_yyyy_xxxzzzz[k] = g_yyy_xxxzzzz[k] - g_0_y_yyy_xxxzzzz[k] * ab_y + g_0_y_yyy_xxxyzzzz[k];

                g_0_y_yyyy_xxyyyyy[k] = g_yyy_xxyyyyy[k] - g_0_y_yyy_xxyyyyy[k] * ab_y + g_0_y_yyy_xxyyyyyy[k];

                g_0_y_yyyy_xxyyyyz[k] = g_yyy_xxyyyyz[k] - g_0_y_yyy_xxyyyyz[k] * ab_y + g_0_y_yyy_xxyyyyyz[k];

                g_0_y_yyyy_xxyyyzz[k] = g_yyy_xxyyyzz[k] - g_0_y_yyy_xxyyyzz[k] * ab_y + g_0_y_yyy_xxyyyyzz[k];

                g_0_y_yyyy_xxyyzzz[k] = g_yyy_xxyyzzz[k] - g_0_y_yyy_xxyyzzz[k] * ab_y + g_0_y_yyy_xxyyyzzz[k];

                g_0_y_yyyy_xxyzzzz[k] = g_yyy_xxyzzzz[k] - g_0_y_yyy_xxyzzzz[k] * ab_y + g_0_y_yyy_xxyyzzzz[k];

                g_0_y_yyyy_xxzzzzz[k] = g_yyy_xxzzzzz[k] - g_0_y_yyy_xxzzzzz[k] * ab_y + g_0_y_yyy_xxyzzzzz[k];

                g_0_y_yyyy_xyyyyyy[k] = g_yyy_xyyyyyy[k] - g_0_y_yyy_xyyyyyy[k] * ab_y + g_0_y_yyy_xyyyyyyy[k];

                g_0_y_yyyy_xyyyyyz[k] = g_yyy_xyyyyyz[k] - g_0_y_yyy_xyyyyyz[k] * ab_y + g_0_y_yyy_xyyyyyyz[k];

                g_0_y_yyyy_xyyyyzz[k] = g_yyy_xyyyyzz[k] - g_0_y_yyy_xyyyyzz[k] * ab_y + g_0_y_yyy_xyyyyyzz[k];

                g_0_y_yyyy_xyyyzzz[k] = g_yyy_xyyyzzz[k] - g_0_y_yyy_xyyyzzz[k] * ab_y + g_0_y_yyy_xyyyyzzz[k];

                g_0_y_yyyy_xyyzzzz[k] = g_yyy_xyyzzzz[k] - g_0_y_yyy_xyyzzzz[k] * ab_y + g_0_y_yyy_xyyyzzzz[k];

                g_0_y_yyyy_xyzzzzz[k] = g_yyy_xyzzzzz[k] - g_0_y_yyy_xyzzzzz[k] * ab_y + g_0_y_yyy_xyyzzzzz[k];

                g_0_y_yyyy_xzzzzzz[k] = g_yyy_xzzzzzz[k] - g_0_y_yyy_xzzzzzz[k] * ab_y + g_0_y_yyy_xyzzzzzz[k];

                g_0_y_yyyy_yyyyyyy[k] = g_yyy_yyyyyyy[k] - g_0_y_yyy_yyyyyyy[k] * ab_y + g_0_y_yyy_yyyyyyyy[k];

                g_0_y_yyyy_yyyyyyz[k] = g_yyy_yyyyyyz[k] - g_0_y_yyy_yyyyyyz[k] * ab_y + g_0_y_yyy_yyyyyyyz[k];

                g_0_y_yyyy_yyyyyzz[k] = g_yyy_yyyyyzz[k] - g_0_y_yyy_yyyyyzz[k] * ab_y + g_0_y_yyy_yyyyyyzz[k];

                g_0_y_yyyy_yyyyzzz[k] = g_yyy_yyyyzzz[k] - g_0_y_yyy_yyyyzzz[k] * ab_y + g_0_y_yyy_yyyyyzzz[k];

                g_0_y_yyyy_yyyzzzz[k] = g_yyy_yyyzzzz[k] - g_0_y_yyy_yyyzzzz[k] * ab_y + g_0_y_yyy_yyyyzzzz[k];

                g_0_y_yyyy_yyzzzzz[k] = g_yyy_yyzzzzz[k] - g_0_y_yyy_yyzzzzz[k] * ab_y + g_0_y_yyy_yyyzzzzz[k];

                g_0_y_yyyy_yzzzzzz[k] = g_yyy_yzzzzzz[k] - g_0_y_yyy_yzzzzzz[k] * ab_y + g_0_y_yyy_yyzzzzzz[k];

                g_0_y_yyyy_zzzzzzz[k] = g_yyy_zzzzzzz[k] - g_0_y_yyy_zzzzzzz[k] * ab_y + g_0_y_yyy_yzzzzzzz[k];
            }

            /// Set up 936-972 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyy_xxxxxxx, g_0_y_yyy_xxxxxxxz, g_0_y_yyy_xxxxxxy, g_0_y_yyy_xxxxxxyz, g_0_y_yyy_xxxxxxz, g_0_y_yyy_xxxxxxzz, g_0_y_yyy_xxxxxyy, g_0_y_yyy_xxxxxyyz, g_0_y_yyy_xxxxxyz, g_0_y_yyy_xxxxxyzz, g_0_y_yyy_xxxxxzz, g_0_y_yyy_xxxxxzzz, g_0_y_yyy_xxxxyyy, g_0_y_yyy_xxxxyyyz, g_0_y_yyy_xxxxyyz, g_0_y_yyy_xxxxyyzz, g_0_y_yyy_xxxxyzz, g_0_y_yyy_xxxxyzzz, g_0_y_yyy_xxxxzzz, g_0_y_yyy_xxxxzzzz, g_0_y_yyy_xxxyyyy, g_0_y_yyy_xxxyyyyz, g_0_y_yyy_xxxyyyz, g_0_y_yyy_xxxyyyzz, g_0_y_yyy_xxxyyzz, g_0_y_yyy_xxxyyzzz, g_0_y_yyy_xxxyzzz, g_0_y_yyy_xxxyzzzz, g_0_y_yyy_xxxzzzz, g_0_y_yyy_xxxzzzzz, g_0_y_yyy_xxyyyyy, g_0_y_yyy_xxyyyyyz, g_0_y_yyy_xxyyyyz, g_0_y_yyy_xxyyyyzz, g_0_y_yyy_xxyyyzz, g_0_y_yyy_xxyyyzzz, g_0_y_yyy_xxyyzzz, g_0_y_yyy_xxyyzzzz, g_0_y_yyy_xxyzzzz, g_0_y_yyy_xxyzzzzz, g_0_y_yyy_xxzzzzz, g_0_y_yyy_xxzzzzzz, g_0_y_yyy_xyyyyyy, g_0_y_yyy_xyyyyyyz, g_0_y_yyy_xyyyyyz, g_0_y_yyy_xyyyyyzz, g_0_y_yyy_xyyyyzz, g_0_y_yyy_xyyyyzzz, g_0_y_yyy_xyyyzzz, g_0_y_yyy_xyyyzzzz, g_0_y_yyy_xyyzzzz, g_0_y_yyy_xyyzzzzz, g_0_y_yyy_xyzzzzz, g_0_y_yyy_xyzzzzzz, g_0_y_yyy_xzzzzzz, g_0_y_yyy_xzzzzzzz, g_0_y_yyy_yyyyyyy, g_0_y_yyy_yyyyyyyz, g_0_y_yyy_yyyyyyz, g_0_y_yyy_yyyyyyzz, g_0_y_yyy_yyyyyzz, g_0_y_yyy_yyyyyzzz, g_0_y_yyy_yyyyzzz, g_0_y_yyy_yyyyzzzz, g_0_y_yyy_yyyzzzz, g_0_y_yyy_yyyzzzzz, g_0_y_yyy_yyzzzzz, g_0_y_yyy_yyzzzzzz, g_0_y_yyy_yzzzzzz, g_0_y_yyy_yzzzzzzz, g_0_y_yyy_zzzzzzz, g_0_y_yyy_zzzzzzzz, g_0_y_yyyz_xxxxxxx, g_0_y_yyyz_xxxxxxy, g_0_y_yyyz_xxxxxxz, g_0_y_yyyz_xxxxxyy, g_0_y_yyyz_xxxxxyz, g_0_y_yyyz_xxxxxzz, g_0_y_yyyz_xxxxyyy, g_0_y_yyyz_xxxxyyz, g_0_y_yyyz_xxxxyzz, g_0_y_yyyz_xxxxzzz, g_0_y_yyyz_xxxyyyy, g_0_y_yyyz_xxxyyyz, g_0_y_yyyz_xxxyyzz, g_0_y_yyyz_xxxyzzz, g_0_y_yyyz_xxxzzzz, g_0_y_yyyz_xxyyyyy, g_0_y_yyyz_xxyyyyz, g_0_y_yyyz_xxyyyzz, g_0_y_yyyz_xxyyzzz, g_0_y_yyyz_xxyzzzz, g_0_y_yyyz_xxzzzzz, g_0_y_yyyz_xyyyyyy, g_0_y_yyyz_xyyyyyz, g_0_y_yyyz_xyyyyzz, g_0_y_yyyz_xyyyzzz, g_0_y_yyyz_xyyzzzz, g_0_y_yyyz_xyzzzzz, g_0_y_yyyz_xzzzzzz, g_0_y_yyyz_yyyyyyy, g_0_y_yyyz_yyyyyyz, g_0_y_yyyz_yyyyyzz, g_0_y_yyyz_yyyyzzz, g_0_y_yyyz_yyyzzzz, g_0_y_yyyz_yyzzzzz, g_0_y_yyyz_yzzzzzz, g_0_y_yyyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_xxxxxxx[k] = -g_0_y_yyy_xxxxxxx[k] * ab_z + g_0_y_yyy_xxxxxxxz[k];

                g_0_y_yyyz_xxxxxxy[k] = -g_0_y_yyy_xxxxxxy[k] * ab_z + g_0_y_yyy_xxxxxxyz[k];

                g_0_y_yyyz_xxxxxxz[k] = -g_0_y_yyy_xxxxxxz[k] * ab_z + g_0_y_yyy_xxxxxxzz[k];

                g_0_y_yyyz_xxxxxyy[k] = -g_0_y_yyy_xxxxxyy[k] * ab_z + g_0_y_yyy_xxxxxyyz[k];

                g_0_y_yyyz_xxxxxyz[k] = -g_0_y_yyy_xxxxxyz[k] * ab_z + g_0_y_yyy_xxxxxyzz[k];

                g_0_y_yyyz_xxxxxzz[k] = -g_0_y_yyy_xxxxxzz[k] * ab_z + g_0_y_yyy_xxxxxzzz[k];

                g_0_y_yyyz_xxxxyyy[k] = -g_0_y_yyy_xxxxyyy[k] * ab_z + g_0_y_yyy_xxxxyyyz[k];

                g_0_y_yyyz_xxxxyyz[k] = -g_0_y_yyy_xxxxyyz[k] * ab_z + g_0_y_yyy_xxxxyyzz[k];

                g_0_y_yyyz_xxxxyzz[k] = -g_0_y_yyy_xxxxyzz[k] * ab_z + g_0_y_yyy_xxxxyzzz[k];

                g_0_y_yyyz_xxxxzzz[k] = -g_0_y_yyy_xxxxzzz[k] * ab_z + g_0_y_yyy_xxxxzzzz[k];

                g_0_y_yyyz_xxxyyyy[k] = -g_0_y_yyy_xxxyyyy[k] * ab_z + g_0_y_yyy_xxxyyyyz[k];

                g_0_y_yyyz_xxxyyyz[k] = -g_0_y_yyy_xxxyyyz[k] * ab_z + g_0_y_yyy_xxxyyyzz[k];

                g_0_y_yyyz_xxxyyzz[k] = -g_0_y_yyy_xxxyyzz[k] * ab_z + g_0_y_yyy_xxxyyzzz[k];

                g_0_y_yyyz_xxxyzzz[k] = -g_0_y_yyy_xxxyzzz[k] * ab_z + g_0_y_yyy_xxxyzzzz[k];

                g_0_y_yyyz_xxxzzzz[k] = -g_0_y_yyy_xxxzzzz[k] * ab_z + g_0_y_yyy_xxxzzzzz[k];

                g_0_y_yyyz_xxyyyyy[k] = -g_0_y_yyy_xxyyyyy[k] * ab_z + g_0_y_yyy_xxyyyyyz[k];

                g_0_y_yyyz_xxyyyyz[k] = -g_0_y_yyy_xxyyyyz[k] * ab_z + g_0_y_yyy_xxyyyyzz[k];

                g_0_y_yyyz_xxyyyzz[k] = -g_0_y_yyy_xxyyyzz[k] * ab_z + g_0_y_yyy_xxyyyzzz[k];

                g_0_y_yyyz_xxyyzzz[k] = -g_0_y_yyy_xxyyzzz[k] * ab_z + g_0_y_yyy_xxyyzzzz[k];

                g_0_y_yyyz_xxyzzzz[k] = -g_0_y_yyy_xxyzzzz[k] * ab_z + g_0_y_yyy_xxyzzzzz[k];

                g_0_y_yyyz_xxzzzzz[k] = -g_0_y_yyy_xxzzzzz[k] * ab_z + g_0_y_yyy_xxzzzzzz[k];

                g_0_y_yyyz_xyyyyyy[k] = -g_0_y_yyy_xyyyyyy[k] * ab_z + g_0_y_yyy_xyyyyyyz[k];

                g_0_y_yyyz_xyyyyyz[k] = -g_0_y_yyy_xyyyyyz[k] * ab_z + g_0_y_yyy_xyyyyyzz[k];

                g_0_y_yyyz_xyyyyzz[k] = -g_0_y_yyy_xyyyyzz[k] * ab_z + g_0_y_yyy_xyyyyzzz[k];

                g_0_y_yyyz_xyyyzzz[k] = -g_0_y_yyy_xyyyzzz[k] * ab_z + g_0_y_yyy_xyyyzzzz[k];

                g_0_y_yyyz_xyyzzzz[k] = -g_0_y_yyy_xyyzzzz[k] * ab_z + g_0_y_yyy_xyyzzzzz[k];

                g_0_y_yyyz_xyzzzzz[k] = -g_0_y_yyy_xyzzzzz[k] * ab_z + g_0_y_yyy_xyzzzzzz[k];

                g_0_y_yyyz_xzzzzzz[k] = -g_0_y_yyy_xzzzzzz[k] * ab_z + g_0_y_yyy_xzzzzzzz[k];

                g_0_y_yyyz_yyyyyyy[k] = -g_0_y_yyy_yyyyyyy[k] * ab_z + g_0_y_yyy_yyyyyyyz[k];

                g_0_y_yyyz_yyyyyyz[k] = -g_0_y_yyy_yyyyyyz[k] * ab_z + g_0_y_yyy_yyyyyyzz[k];

                g_0_y_yyyz_yyyyyzz[k] = -g_0_y_yyy_yyyyyzz[k] * ab_z + g_0_y_yyy_yyyyyzzz[k];

                g_0_y_yyyz_yyyyzzz[k] = -g_0_y_yyy_yyyyzzz[k] * ab_z + g_0_y_yyy_yyyyzzzz[k];

                g_0_y_yyyz_yyyzzzz[k] = -g_0_y_yyy_yyyzzzz[k] * ab_z + g_0_y_yyy_yyyzzzzz[k];

                g_0_y_yyyz_yyzzzzz[k] = -g_0_y_yyy_yyzzzzz[k] * ab_z + g_0_y_yyy_yyzzzzzz[k];

                g_0_y_yyyz_yzzzzzz[k] = -g_0_y_yyy_yzzzzzz[k] * ab_z + g_0_y_yyy_yzzzzzzz[k];

                g_0_y_yyyz_zzzzzzz[k] = -g_0_y_yyy_zzzzzzz[k] * ab_z + g_0_y_yyy_zzzzzzzz[k];
            }

            /// Set up 972-1008 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyz_xxxxxxx, g_0_y_yyz_xxxxxxxz, g_0_y_yyz_xxxxxxy, g_0_y_yyz_xxxxxxyz, g_0_y_yyz_xxxxxxz, g_0_y_yyz_xxxxxxzz, g_0_y_yyz_xxxxxyy, g_0_y_yyz_xxxxxyyz, g_0_y_yyz_xxxxxyz, g_0_y_yyz_xxxxxyzz, g_0_y_yyz_xxxxxzz, g_0_y_yyz_xxxxxzzz, g_0_y_yyz_xxxxyyy, g_0_y_yyz_xxxxyyyz, g_0_y_yyz_xxxxyyz, g_0_y_yyz_xxxxyyzz, g_0_y_yyz_xxxxyzz, g_0_y_yyz_xxxxyzzz, g_0_y_yyz_xxxxzzz, g_0_y_yyz_xxxxzzzz, g_0_y_yyz_xxxyyyy, g_0_y_yyz_xxxyyyyz, g_0_y_yyz_xxxyyyz, g_0_y_yyz_xxxyyyzz, g_0_y_yyz_xxxyyzz, g_0_y_yyz_xxxyyzzz, g_0_y_yyz_xxxyzzz, g_0_y_yyz_xxxyzzzz, g_0_y_yyz_xxxzzzz, g_0_y_yyz_xxxzzzzz, g_0_y_yyz_xxyyyyy, g_0_y_yyz_xxyyyyyz, g_0_y_yyz_xxyyyyz, g_0_y_yyz_xxyyyyzz, g_0_y_yyz_xxyyyzz, g_0_y_yyz_xxyyyzzz, g_0_y_yyz_xxyyzzz, g_0_y_yyz_xxyyzzzz, g_0_y_yyz_xxyzzzz, g_0_y_yyz_xxyzzzzz, g_0_y_yyz_xxzzzzz, g_0_y_yyz_xxzzzzzz, g_0_y_yyz_xyyyyyy, g_0_y_yyz_xyyyyyyz, g_0_y_yyz_xyyyyyz, g_0_y_yyz_xyyyyyzz, g_0_y_yyz_xyyyyzz, g_0_y_yyz_xyyyyzzz, g_0_y_yyz_xyyyzzz, g_0_y_yyz_xyyyzzzz, g_0_y_yyz_xyyzzzz, g_0_y_yyz_xyyzzzzz, g_0_y_yyz_xyzzzzz, g_0_y_yyz_xyzzzzzz, g_0_y_yyz_xzzzzzz, g_0_y_yyz_xzzzzzzz, g_0_y_yyz_yyyyyyy, g_0_y_yyz_yyyyyyyz, g_0_y_yyz_yyyyyyz, g_0_y_yyz_yyyyyyzz, g_0_y_yyz_yyyyyzz, g_0_y_yyz_yyyyyzzz, g_0_y_yyz_yyyyzzz, g_0_y_yyz_yyyyzzzz, g_0_y_yyz_yyyzzzz, g_0_y_yyz_yyyzzzzz, g_0_y_yyz_yyzzzzz, g_0_y_yyz_yyzzzzzz, g_0_y_yyz_yzzzzzz, g_0_y_yyz_yzzzzzzz, g_0_y_yyz_zzzzzzz, g_0_y_yyz_zzzzzzzz, g_0_y_yyzz_xxxxxxx, g_0_y_yyzz_xxxxxxy, g_0_y_yyzz_xxxxxxz, g_0_y_yyzz_xxxxxyy, g_0_y_yyzz_xxxxxyz, g_0_y_yyzz_xxxxxzz, g_0_y_yyzz_xxxxyyy, g_0_y_yyzz_xxxxyyz, g_0_y_yyzz_xxxxyzz, g_0_y_yyzz_xxxxzzz, g_0_y_yyzz_xxxyyyy, g_0_y_yyzz_xxxyyyz, g_0_y_yyzz_xxxyyzz, g_0_y_yyzz_xxxyzzz, g_0_y_yyzz_xxxzzzz, g_0_y_yyzz_xxyyyyy, g_0_y_yyzz_xxyyyyz, g_0_y_yyzz_xxyyyzz, g_0_y_yyzz_xxyyzzz, g_0_y_yyzz_xxyzzzz, g_0_y_yyzz_xxzzzzz, g_0_y_yyzz_xyyyyyy, g_0_y_yyzz_xyyyyyz, g_0_y_yyzz_xyyyyzz, g_0_y_yyzz_xyyyzzz, g_0_y_yyzz_xyyzzzz, g_0_y_yyzz_xyzzzzz, g_0_y_yyzz_xzzzzzz, g_0_y_yyzz_yyyyyyy, g_0_y_yyzz_yyyyyyz, g_0_y_yyzz_yyyyyzz, g_0_y_yyzz_yyyyzzz, g_0_y_yyzz_yyyzzzz, g_0_y_yyzz_yyzzzzz, g_0_y_yyzz_yzzzzzz, g_0_y_yyzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_xxxxxxx[k] = -g_0_y_yyz_xxxxxxx[k] * ab_z + g_0_y_yyz_xxxxxxxz[k];

                g_0_y_yyzz_xxxxxxy[k] = -g_0_y_yyz_xxxxxxy[k] * ab_z + g_0_y_yyz_xxxxxxyz[k];

                g_0_y_yyzz_xxxxxxz[k] = -g_0_y_yyz_xxxxxxz[k] * ab_z + g_0_y_yyz_xxxxxxzz[k];

                g_0_y_yyzz_xxxxxyy[k] = -g_0_y_yyz_xxxxxyy[k] * ab_z + g_0_y_yyz_xxxxxyyz[k];

                g_0_y_yyzz_xxxxxyz[k] = -g_0_y_yyz_xxxxxyz[k] * ab_z + g_0_y_yyz_xxxxxyzz[k];

                g_0_y_yyzz_xxxxxzz[k] = -g_0_y_yyz_xxxxxzz[k] * ab_z + g_0_y_yyz_xxxxxzzz[k];

                g_0_y_yyzz_xxxxyyy[k] = -g_0_y_yyz_xxxxyyy[k] * ab_z + g_0_y_yyz_xxxxyyyz[k];

                g_0_y_yyzz_xxxxyyz[k] = -g_0_y_yyz_xxxxyyz[k] * ab_z + g_0_y_yyz_xxxxyyzz[k];

                g_0_y_yyzz_xxxxyzz[k] = -g_0_y_yyz_xxxxyzz[k] * ab_z + g_0_y_yyz_xxxxyzzz[k];

                g_0_y_yyzz_xxxxzzz[k] = -g_0_y_yyz_xxxxzzz[k] * ab_z + g_0_y_yyz_xxxxzzzz[k];

                g_0_y_yyzz_xxxyyyy[k] = -g_0_y_yyz_xxxyyyy[k] * ab_z + g_0_y_yyz_xxxyyyyz[k];

                g_0_y_yyzz_xxxyyyz[k] = -g_0_y_yyz_xxxyyyz[k] * ab_z + g_0_y_yyz_xxxyyyzz[k];

                g_0_y_yyzz_xxxyyzz[k] = -g_0_y_yyz_xxxyyzz[k] * ab_z + g_0_y_yyz_xxxyyzzz[k];

                g_0_y_yyzz_xxxyzzz[k] = -g_0_y_yyz_xxxyzzz[k] * ab_z + g_0_y_yyz_xxxyzzzz[k];

                g_0_y_yyzz_xxxzzzz[k] = -g_0_y_yyz_xxxzzzz[k] * ab_z + g_0_y_yyz_xxxzzzzz[k];

                g_0_y_yyzz_xxyyyyy[k] = -g_0_y_yyz_xxyyyyy[k] * ab_z + g_0_y_yyz_xxyyyyyz[k];

                g_0_y_yyzz_xxyyyyz[k] = -g_0_y_yyz_xxyyyyz[k] * ab_z + g_0_y_yyz_xxyyyyzz[k];

                g_0_y_yyzz_xxyyyzz[k] = -g_0_y_yyz_xxyyyzz[k] * ab_z + g_0_y_yyz_xxyyyzzz[k];

                g_0_y_yyzz_xxyyzzz[k] = -g_0_y_yyz_xxyyzzz[k] * ab_z + g_0_y_yyz_xxyyzzzz[k];

                g_0_y_yyzz_xxyzzzz[k] = -g_0_y_yyz_xxyzzzz[k] * ab_z + g_0_y_yyz_xxyzzzzz[k];

                g_0_y_yyzz_xxzzzzz[k] = -g_0_y_yyz_xxzzzzz[k] * ab_z + g_0_y_yyz_xxzzzzzz[k];

                g_0_y_yyzz_xyyyyyy[k] = -g_0_y_yyz_xyyyyyy[k] * ab_z + g_0_y_yyz_xyyyyyyz[k];

                g_0_y_yyzz_xyyyyyz[k] = -g_0_y_yyz_xyyyyyz[k] * ab_z + g_0_y_yyz_xyyyyyzz[k];

                g_0_y_yyzz_xyyyyzz[k] = -g_0_y_yyz_xyyyyzz[k] * ab_z + g_0_y_yyz_xyyyyzzz[k];

                g_0_y_yyzz_xyyyzzz[k] = -g_0_y_yyz_xyyyzzz[k] * ab_z + g_0_y_yyz_xyyyzzzz[k];

                g_0_y_yyzz_xyyzzzz[k] = -g_0_y_yyz_xyyzzzz[k] * ab_z + g_0_y_yyz_xyyzzzzz[k];

                g_0_y_yyzz_xyzzzzz[k] = -g_0_y_yyz_xyzzzzz[k] * ab_z + g_0_y_yyz_xyzzzzzz[k];

                g_0_y_yyzz_xzzzzzz[k] = -g_0_y_yyz_xzzzzzz[k] * ab_z + g_0_y_yyz_xzzzzzzz[k];

                g_0_y_yyzz_yyyyyyy[k] = -g_0_y_yyz_yyyyyyy[k] * ab_z + g_0_y_yyz_yyyyyyyz[k];

                g_0_y_yyzz_yyyyyyz[k] = -g_0_y_yyz_yyyyyyz[k] * ab_z + g_0_y_yyz_yyyyyyzz[k];

                g_0_y_yyzz_yyyyyzz[k] = -g_0_y_yyz_yyyyyzz[k] * ab_z + g_0_y_yyz_yyyyyzzz[k];

                g_0_y_yyzz_yyyyzzz[k] = -g_0_y_yyz_yyyyzzz[k] * ab_z + g_0_y_yyz_yyyyzzzz[k];

                g_0_y_yyzz_yyyzzzz[k] = -g_0_y_yyz_yyyzzzz[k] * ab_z + g_0_y_yyz_yyyzzzzz[k];

                g_0_y_yyzz_yyzzzzz[k] = -g_0_y_yyz_yyzzzzz[k] * ab_z + g_0_y_yyz_yyzzzzzz[k];

                g_0_y_yyzz_yzzzzzz[k] = -g_0_y_yyz_yzzzzzz[k] * ab_z + g_0_y_yyz_yzzzzzzz[k];

                g_0_y_yyzz_zzzzzzz[k] = -g_0_y_yyz_zzzzzzz[k] * ab_z + g_0_y_yyz_zzzzzzzz[k];
            }

            /// Set up 1008-1044 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yzz_xxxxxxx, g_0_y_yzz_xxxxxxxz, g_0_y_yzz_xxxxxxy, g_0_y_yzz_xxxxxxyz, g_0_y_yzz_xxxxxxz, g_0_y_yzz_xxxxxxzz, g_0_y_yzz_xxxxxyy, g_0_y_yzz_xxxxxyyz, g_0_y_yzz_xxxxxyz, g_0_y_yzz_xxxxxyzz, g_0_y_yzz_xxxxxzz, g_0_y_yzz_xxxxxzzz, g_0_y_yzz_xxxxyyy, g_0_y_yzz_xxxxyyyz, g_0_y_yzz_xxxxyyz, g_0_y_yzz_xxxxyyzz, g_0_y_yzz_xxxxyzz, g_0_y_yzz_xxxxyzzz, g_0_y_yzz_xxxxzzz, g_0_y_yzz_xxxxzzzz, g_0_y_yzz_xxxyyyy, g_0_y_yzz_xxxyyyyz, g_0_y_yzz_xxxyyyz, g_0_y_yzz_xxxyyyzz, g_0_y_yzz_xxxyyzz, g_0_y_yzz_xxxyyzzz, g_0_y_yzz_xxxyzzz, g_0_y_yzz_xxxyzzzz, g_0_y_yzz_xxxzzzz, g_0_y_yzz_xxxzzzzz, g_0_y_yzz_xxyyyyy, g_0_y_yzz_xxyyyyyz, g_0_y_yzz_xxyyyyz, g_0_y_yzz_xxyyyyzz, g_0_y_yzz_xxyyyzz, g_0_y_yzz_xxyyyzzz, g_0_y_yzz_xxyyzzz, g_0_y_yzz_xxyyzzzz, g_0_y_yzz_xxyzzzz, g_0_y_yzz_xxyzzzzz, g_0_y_yzz_xxzzzzz, g_0_y_yzz_xxzzzzzz, g_0_y_yzz_xyyyyyy, g_0_y_yzz_xyyyyyyz, g_0_y_yzz_xyyyyyz, g_0_y_yzz_xyyyyyzz, g_0_y_yzz_xyyyyzz, g_0_y_yzz_xyyyyzzz, g_0_y_yzz_xyyyzzz, g_0_y_yzz_xyyyzzzz, g_0_y_yzz_xyyzzzz, g_0_y_yzz_xyyzzzzz, g_0_y_yzz_xyzzzzz, g_0_y_yzz_xyzzzzzz, g_0_y_yzz_xzzzzzz, g_0_y_yzz_xzzzzzzz, g_0_y_yzz_yyyyyyy, g_0_y_yzz_yyyyyyyz, g_0_y_yzz_yyyyyyz, g_0_y_yzz_yyyyyyzz, g_0_y_yzz_yyyyyzz, g_0_y_yzz_yyyyyzzz, g_0_y_yzz_yyyyzzz, g_0_y_yzz_yyyyzzzz, g_0_y_yzz_yyyzzzz, g_0_y_yzz_yyyzzzzz, g_0_y_yzz_yyzzzzz, g_0_y_yzz_yyzzzzzz, g_0_y_yzz_yzzzzzz, g_0_y_yzz_yzzzzzzz, g_0_y_yzz_zzzzzzz, g_0_y_yzz_zzzzzzzz, g_0_y_yzzz_xxxxxxx, g_0_y_yzzz_xxxxxxy, g_0_y_yzzz_xxxxxxz, g_0_y_yzzz_xxxxxyy, g_0_y_yzzz_xxxxxyz, g_0_y_yzzz_xxxxxzz, g_0_y_yzzz_xxxxyyy, g_0_y_yzzz_xxxxyyz, g_0_y_yzzz_xxxxyzz, g_0_y_yzzz_xxxxzzz, g_0_y_yzzz_xxxyyyy, g_0_y_yzzz_xxxyyyz, g_0_y_yzzz_xxxyyzz, g_0_y_yzzz_xxxyzzz, g_0_y_yzzz_xxxzzzz, g_0_y_yzzz_xxyyyyy, g_0_y_yzzz_xxyyyyz, g_0_y_yzzz_xxyyyzz, g_0_y_yzzz_xxyyzzz, g_0_y_yzzz_xxyzzzz, g_0_y_yzzz_xxzzzzz, g_0_y_yzzz_xyyyyyy, g_0_y_yzzz_xyyyyyz, g_0_y_yzzz_xyyyyzz, g_0_y_yzzz_xyyyzzz, g_0_y_yzzz_xyyzzzz, g_0_y_yzzz_xyzzzzz, g_0_y_yzzz_xzzzzzz, g_0_y_yzzz_yyyyyyy, g_0_y_yzzz_yyyyyyz, g_0_y_yzzz_yyyyyzz, g_0_y_yzzz_yyyyzzz, g_0_y_yzzz_yyyzzzz, g_0_y_yzzz_yyzzzzz, g_0_y_yzzz_yzzzzzz, g_0_y_yzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_xxxxxxx[k] = -g_0_y_yzz_xxxxxxx[k] * ab_z + g_0_y_yzz_xxxxxxxz[k];

                g_0_y_yzzz_xxxxxxy[k] = -g_0_y_yzz_xxxxxxy[k] * ab_z + g_0_y_yzz_xxxxxxyz[k];

                g_0_y_yzzz_xxxxxxz[k] = -g_0_y_yzz_xxxxxxz[k] * ab_z + g_0_y_yzz_xxxxxxzz[k];

                g_0_y_yzzz_xxxxxyy[k] = -g_0_y_yzz_xxxxxyy[k] * ab_z + g_0_y_yzz_xxxxxyyz[k];

                g_0_y_yzzz_xxxxxyz[k] = -g_0_y_yzz_xxxxxyz[k] * ab_z + g_0_y_yzz_xxxxxyzz[k];

                g_0_y_yzzz_xxxxxzz[k] = -g_0_y_yzz_xxxxxzz[k] * ab_z + g_0_y_yzz_xxxxxzzz[k];

                g_0_y_yzzz_xxxxyyy[k] = -g_0_y_yzz_xxxxyyy[k] * ab_z + g_0_y_yzz_xxxxyyyz[k];

                g_0_y_yzzz_xxxxyyz[k] = -g_0_y_yzz_xxxxyyz[k] * ab_z + g_0_y_yzz_xxxxyyzz[k];

                g_0_y_yzzz_xxxxyzz[k] = -g_0_y_yzz_xxxxyzz[k] * ab_z + g_0_y_yzz_xxxxyzzz[k];

                g_0_y_yzzz_xxxxzzz[k] = -g_0_y_yzz_xxxxzzz[k] * ab_z + g_0_y_yzz_xxxxzzzz[k];

                g_0_y_yzzz_xxxyyyy[k] = -g_0_y_yzz_xxxyyyy[k] * ab_z + g_0_y_yzz_xxxyyyyz[k];

                g_0_y_yzzz_xxxyyyz[k] = -g_0_y_yzz_xxxyyyz[k] * ab_z + g_0_y_yzz_xxxyyyzz[k];

                g_0_y_yzzz_xxxyyzz[k] = -g_0_y_yzz_xxxyyzz[k] * ab_z + g_0_y_yzz_xxxyyzzz[k];

                g_0_y_yzzz_xxxyzzz[k] = -g_0_y_yzz_xxxyzzz[k] * ab_z + g_0_y_yzz_xxxyzzzz[k];

                g_0_y_yzzz_xxxzzzz[k] = -g_0_y_yzz_xxxzzzz[k] * ab_z + g_0_y_yzz_xxxzzzzz[k];

                g_0_y_yzzz_xxyyyyy[k] = -g_0_y_yzz_xxyyyyy[k] * ab_z + g_0_y_yzz_xxyyyyyz[k];

                g_0_y_yzzz_xxyyyyz[k] = -g_0_y_yzz_xxyyyyz[k] * ab_z + g_0_y_yzz_xxyyyyzz[k];

                g_0_y_yzzz_xxyyyzz[k] = -g_0_y_yzz_xxyyyzz[k] * ab_z + g_0_y_yzz_xxyyyzzz[k];

                g_0_y_yzzz_xxyyzzz[k] = -g_0_y_yzz_xxyyzzz[k] * ab_z + g_0_y_yzz_xxyyzzzz[k];

                g_0_y_yzzz_xxyzzzz[k] = -g_0_y_yzz_xxyzzzz[k] * ab_z + g_0_y_yzz_xxyzzzzz[k];

                g_0_y_yzzz_xxzzzzz[k] = -g_0_y_yzz_xxzzzzz[k] * ab_z + g_0_y_yzz_xxzzzzzz[k];

                g_0_y_yzzz_xyyyyyy[k] = -g_0_y_yzz_xyyyyyy[k] * ab_z + g_0_y_yzz_xyyyyyyz[k];

                g_0_y_yzzz_xyyyyyz[k] = -g_0_y_yzz_xyyyyyz[k] * ab_z + g_0_y_yzz_xyyyyyzz[k];

                g_0_y_yzzz_xyyyyzz[k] = -g_0_y_yzz_xyyyyzz[k] * ab_z + g_0_y_yzz_xyyyyzzz[k];

                g_0_y_yzzz_xyyyzzz[k] = -g_0_y_yzz_xyyyzzz[k] * ab_z + g_0_y_yzz_xyyyzzzz[k];

                g_0_y_yzzz_xyyzzzz[k] = -g_0_y_yzz_xyyzzzz[k] * ab_z + g_0_y_yzz_xyyzzzzz[k];

                g_0_y_yzzz_xyzzzzz[k] = -g_0_y_yzz_xyzzzzz[k] * ab_z + g_0_y_yzz_xyzzzzzz[k];

                g_0_y_yzzz_xzzzzzz[k] = -g_0_y_yzz_xzzzzzz[k] * ab_z + g_0_y_yzz_xzzzzzzz[k];

                g_0_y_yzzz_yyyyyyy[k] = -g_0_y_yzz_yyyyyyy[k] * ab_z + g_0_y_yzz_yyyyyyyz[k];

                g_0_y_yzzz_yyyyyyz[k] = -g_0_y_yzz_yyyyyyz[k] * ab_z + g_0_y_yzz_yyyyyyzz[k];

                g_0_y_yzzz_yyyyyzz[k] = -g_0_y_yzz_yyyyyzz[k] * ab_z + g_0_y_yzz_yyyyyzzz[k];

                g_0_y_yzzz_yyyyzzz[k] = -g_0_y_yzz_yyyyzzz[k] * ab_z + g_0_y_yzz_yyyyzzzz[k];

                g_0_y_yzzz_yyyzzzz[k] = -g_0_y_yzz_yyyzzzz[k] * ab_z + g_0_y_yzz_yyyzzzzz[k];

                g_0_y_yzzz_yyzzzzz[k] = -g_0_y_yzz_yyzzzzz[k] * ab_z + g_0_y_yzz_yyzzzzzz[k];

                g_0_y_yzzz_yzzzzzz[k] = -g_0_y_yzz_yzzzzzz[k] * ab_z + g_0_y_yzz_yzzzzzzz[k];

                g_0_y_yzzz_zzzzzzz[k] = -g_0_y_yzz_zzzzzzz[k] * ab_z + g_0_y_yzz_zzzzzzzz[k];
            }

            /// Set up 1044-1080 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_zzz_xxxxxxx, g_0_y_zzz_xxxxxxxz, g_0_y_zzz_xxxxxxy, g_0_y_zzz_xxxxxxyz, g_0_y_zzz_xxxxxxz, g_0_y_zzz_xxxxxxzz, g_0_y_zzz_xxxxxyy, g_0_y_zzz_xxxxxyyz, g_0_y_zzz_xxxxxyz, g_0_y_zzz_xxxxxyzz, g_0_y_zzz_xxxxxzz, g_0_y_zzz_xxxxxzzz, g_0_y_zzz_xxxxyyy, g_0_y_zzz_xxxxyyyz, g_0_y_zzz_xxxxyyz, g_0_y_zzz_xxxxyyzz, g_0_y_zzz_xxxxyzz, g_0_y_zzz_xxxxyzzz, g_0_y_zzz_xxxxzzz, g_0_y_zzz_xxxxzzzz, g_0_y_zzz_xxxyyyy, g_0_y_zzz_xxxyyyyz, g_0_y_zzz_xxxyyyz, g_0_y_zzz_xxxyyyzz, g_0_y_zzz_xxxyyzz, g_0_y_zzz_xxxyyzzz, g_0_y_zzz_xxxyzzz, g_0_y_zzz_xxxyzzzz, g_0_y_zzz_xxxzzzz, g_0_y_zzz_xxxzzzzz, g_0_y_zzz_xxyyyyy, g_0_y_zzz_xxyyyyyz, g_0_y_zzz_xxyyyyz, g_0_y_zzz_xxyyyyzz, g_0_y_zzz_xxyyyzz, g_0_y_zzz_xxyyyzzz, g_0_y_zzz_xxyyzzz, g_0_y_zzz_xxyyzzzz, g_0_y_zzz_xxyzzzz, g_0_y_zzz_xxyzzzzz, g_0_y_zzz_xxzzzzz, g_0_y_zzz_xxzzzzzz, g_0_y_zzz_xyyyyyy, g_0_y_zzz_xyyyyyyz, g_0_y_zzz_xyyyyyz, g_0_y_zzz_xyyyyyzz, g_0_y_zzz_xyyyyzz, g_0_y_zzz_xyyyyzzz, g_0_y_zzz_xyyyzzz, g_0_y_zzz_xyyyzzzz, g_0_y_zzz_xyyzzzz, g_0_y_zzz_xyyzzzzz, g_0_y_zzz_xyzzzzz, g_0_y_zzz_xyzzzzzz, g_0_y_zzz_xzzzzzz, g_0_y_zzz_xzzzzzzz, g_0_y_zzz_yyyyyyy, g_0_y_zzz_yyyyyyyz, g_0_y_zzz_yyyyyyz, g_0_y_zzz_yyyyyyzz, g_0_y_zzz_yyyyyzz, g_0_y_zzz_yyyyyzzz, g_0_y_zzz_yyyyzzz, g_0_y_zzz_yyyyzzzz, g_0_y_zzz_yyyzzzz, g_0_y_zzz_yyyzzzzz, g_0_y_zzz_yyzzzzz, g_0_y_zzz_yyzzzzzz, g_0_y_zzz_yzzzzzz, g_0_y_zzz_yzzzzzzz, g_0_y_zzz_zzzzzzz, g_0_y_zzz_zzzzzzzz, g_0_y_zzzz_xxxxxxx, g_0_y_zzzz_xxxxxxy, g_0_y_zzzz_xxxxxxz, g_0_y_zzzz_xxxxxyy, g_0_y_zzzz_xxxxxyz, g_0_y_zzzz_xxxxxzz, g_0_y_zzzz_xxxxyyy, g_0_y_zzzz_xxxxyyz, g_0_y_zzzz_xxxxyzz, g_0_y_zzzz_xxxxzzz, g_0_y_zzzz_xxxyyyy, g_0_y_zzzz_xxxyyyz, g_0_y_zzzz_xxxyyzz, g_0_y_zzzz_xxxyzzz, g_0_y_zzzz_xxxzzzz, g_0_y_zzzz_xxyyyyy, g_0_y_zzzz_xxyyyyz, g_0_y_zzzz_xxyyyzz, g_0_y_zzzz_xxyyzzz, g_0_y_zzzz_xxyzzzz, g_0_y_zzzz_xxzzzzz, g_0_y_zzzz_xyyyyyy, g_0_y_zzzz_xyyyyyz, g_0_y_zzzz_xyyyyzz, g_0_y_zzzz_xyyyzzz, g_0_y_zzzz_xyyzzzz, g_0_y_zzzz_xyzzzzz, g_0_y_zzzz_xzzzzzz, g_0_y_zzzz_yyyyyyy, g_0_y_zzzz_yyyyyyz, g_0_y_zzzz_yyyyyzz, g_0_y_zzzz_yyyyzzz, g_0_y_zzzz_yyyzzzz, g_0_y_zzzz_yyzzzzz, g_0_y_zzzz_yzzzzzz, g_0_y_zzzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_xxxxxxx[k] = -g_0_y_zzz_xxxxxxx[k] * ab_z + g_0_y_zzz_xxxxxxxz[k];

                g_0_y_zzzz_xxxxxxy[k] = -g_0_y_zzz_xxxxxxy[k] * ab_z + g_0_y_zzz_xxxxxxyz[k];

                g_0_y_zzzz_xxxxxxz[k] = -g_0_y_zzz_xxxxxxz[k] * ab_z + g_0_y_zzz_xxxxxxzz[k];

                g_0_y_zzzz_xxxxxyy[k] = -g_0_y_zzz_xxxxxyy[k] * ab_z + g_0_y_zzz_xxxxxyyz[k];

                g_0_y_zzzz_xxxxxyz[k] = -g_0_y_zzz_xxxxxyz[k] * ab_z + g_0_y_zzz_xxxxxyzz[k];

                g_0_y_zzzz_xxxxxzz[k] = -g_0_y_zzz_xxxxxzz[k] * ab_z + g_0_y_zzz_xxxxxzzz[k];

                g_0_y_zzzz_xxxxyyy[k] = -g_0_y_zzz_xxxxyyy[k] * ab_z + g_0_y_zzz_xxxxyyyz[k];

                g_0_y_zzzz_xxxxyyz[k] = -g_0_y_zzz_xxxxyyz[k] * ab_z + g_0_y_zzz_xxxxyyzz[k];

                g_0_y_zzzz_xxxxyzz[k] = -g_0_y_zzz_xxxxyzz[k] * ab_z + g_0_y_zzz_xxxxyzzz[k];

                g_0_y_zzzz_xxxxzzz[k] = -g_0_y_zzz_xxxxzzz[k] * ab_z + g_0_y_zzz_xxxxzzzz[k];

                g_0_y_zzzz_xxxyyyy[k] = -g_0_y_zzz_xxxyyyy[k] * ab_z + g_0_y_zzz_xxxyyyyz[k];

                g_0_y_zzzz_xxxyyyz[k] = -g_0_y_zzz_xxxyyyz[k] * ab_z + g_0_y_zzz_xxxyyyzz[k];

                g_0_y_zzzz_xxxyyzz[k] = -g_0_y_zzz_xxxyyzz[k] * ab_z + g_0_y_zzz_xxxyyzzz[k];

                g_0_y_zzzz_xxxyzzz[k] = -g_0_y_zzz_xxxyzzz[k] * ab_z + g_0_y_zzz_xxxyzzzz[k];

                g_0_y_zzzz_xxxzzzz[k] = -g_0_y_zzz_xxxzzzz[k] * ab_z + g_0_y_zzz_xxxzzzzz[k];

                g_0_y_zzzz_xxyyyyy[k] = -g_0_y_zzz_xxyyyyy[k] * ab_z + g_0_y_zzz_xxyyyyyz[k];

                g_0_y_zzzz_xxyyyyz[k] = -g_0_y_zzz_xxyyyyz[k] * ab_z + g_0_y_zzz_xxyyyyzz[k];

                g_0_y_zzzz_xxyyyzz[k] = -g_0_y_zzz_xxyyyzz[k] * ab_z + g_0_y_zzz_xxyyyzzz[k];

                g_0_y_zzzz_xxyyzzz[k] = -g_0_y_zzz_xxyyzzz[k] * ab_z + g_0_y_zzz_xxyyzzzz[k];

                g_0_y_zzzz_xxyzzzz[k] = -g_0_y_zzz_xxyzzzz[k] * ab_z + g_0_y_zzz_xxyzzzzz[k];

                g_0_y_zzzz_xxzzzzz[k] = -g_0_y_zzz_xxzzzzz[k] * ab_z + g_0_y_zzz_xxzzzzzz[k];

                g_0_y_zzzz_xyyyyyy[k] = -g_0_y_zzz_xyyyyyy[k] * ab_z + g_0_y_zzz_xyyyyyyz[k];

                g_0_y_zzzz_xyyyyyz[k] = -g_0_y_zzz_xyyyyyz[k] * ab_z + g_0_y_zzz_xyyyyyzz[k];

                g_0_y_zzzz_xyyyyzz[k] = -g_0_y_zzz_xyyyyzz[k] * ab_z + g_0_y_zzz_xyyyyzzz[k];

                g_0_y_zzzz_xyyyzzz[k] = -g_0_y_zzz_xyyyzzz[k] * ab_z + g_0_y_zzz_xyyyzzzz[k];

                g_0_y_zzzz_xyyzzzz[k] = -g_0_y_zzz_xyyzzzz[k] * ab_z + g_0_y_zzz_xyyzzzzz[k];

                g_0_y_zzzz_xyzzzzz[k] = -g_0_y_zzz_xyzzzzz[k] * ab_z + g_0_y_zzz_xyzzzzzz[k];

                g_0_y_zzzz_xzzzzzz[k] = -g_0_y_zzz_xzzzzzz[k] * ab_z + g_0_y_zzz_xzzzzzzz[k];

                g_0_y_zzzz_yyyyyyy[k] = -g_0_y_zzz_yyyyyyy[k] * ab_z + g_0_y_zzz_yyyyyyyz[k];

                g_0_y_zzzz_yyyyyyz[k] = -g_0_y_zzz_yyyyyyz[k] * ab_z + g_0_y_zzz_yyyyyyzz[k];

                g_0_y_zzzz_yyyyyzz[k] = -g_0_y_zzz_yyyyyzz[k] * ab_z + g_0_y_zzz_yyyyyzzz[k];

                g_0_y_zzzz_yyyyzzz[k] = -g_0_y_zzz_yyyyzzz[k] * ab_z + g_0_y_zzz_yyyyzzzz[k];

                g_0_y_zzzz_yyyzzzz[k] = -g_0_y_zzz_yyyzzzz[k] * ab_z + g_0_y_zzz_yyyzzzzz[k];

                g_0_y_zzzz_yyzzzzz[k] = -g_0_y_zzz_yyzzzzz[k] * ab_z + g_0_y_zzz_yyzzzzzz[k];

                g_0_y_zzzz_yzzzzzz[k] = -g_0_y_zzz_yzzzzzz[k] * ab_z + g_0_y_zzz_yzzzzzzz[k];

                g_0_y_zzzz_zzzzzzz[k] = -g_0_y_zzz_zzzzzzz[k] * ab_z + g_0_y_zzz_zzzzzzzz[k];
            }

            /// Set up 1080-1116 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxx_xxxxxxx, g_0_z_xxx_xxxxxxxx, g_0_z_xxx_xxxxxxxy, g_0_z_xxx_xxxxxxxz, g_0_z_xxx_xxxxxxy, g_0_z_xxx_xxxxxxyy, g_0_z_xxx_xxxxxxyz, g_0_z_xxx_xxxxxxz, g_0_z_xxx_xxxxxxzz, g_0_z_xxx_xxxxxyy, g_0_z_xxx_xxxxxyyy, g_0_z_xxx_xxxxxyyz, g_0_z_xxx_xxxxxyz, g_0_z_xxx_xxxxxyzz, g_0_z_xxx_xxxxxzz, g_0_z_xxx_xxxxxzzz, g_0_z_xxx_xxxxyyy, g_0_z_xxx_xxxxyyyy, g_0_z_xxx_xxxxyyyz, g_0_z_xxx_xxxxyyz, g_0_z_xxx_xxxxyyzz, g_0_z_xxx_xxxxyzz, g_0_z_xxx_xxxxyzzz, g_0_z_xxx_xxxxzzz, g_0_z_xxx_xxxxzzzz, g_0_z_xxx_xxxyyyy, g_0_z_xxx_xxxyyyyy, g_0_z_xxx_xxxyyyyz, g_0_z_xxx_xxxyyyz, g_0_z_xxx_xxxyyyzz, g_0_z_xxx_xxxyyzz, g_0_z_xxx_xxxyyzzz, g_0_z_xxx_xxxyzzz, g_0_z_xxx_xxxyzzzz, g_0_z_xxx_xxxzzzz, g_0_z_xxx_xxxzzzzz, g_0_z_xxx_xxyyyyy, g_0_z_xxx_xxyyyyyy, g_0_z_xxx_xxyyyyyz, g_0_z_xxx_xxyyyyz, g_0_z_xxx_xxyyyyzz, g_0_z_xxx_xxyyyzz, g_0_z_xxx_xxyyyzzz, g_0_z_xxx_xxyyzzz, g_0_z_xxx_xxyyzzzz, g_0_z_xxx_xxyzzzz, g_0_z_xxx_xxyzzzzz, g_0_z_xxx_xxzzzzz, g_0_z_xxx_xxzzzzzz, g_0_z_xxx_xyyyyyy, g_0_z_xxx_xyyyyyyy, g_0_z_xxx_xyyyyyyz, g_0_z_xxx_xyyyyyz, g_0_z_xxx_xyyyyyzz, g_0_z_xxx_xyyyyzz, g_0_z_xxx_xyyyyzzz, g_0_z_xxx_xyyyzzz, g_0_z_xxx_xyyyzzzz, g_0_z_xxx_xyyzzzz, g_0_z_xxx_xyyzzzzz, g_0_z_xxx_xyzzzzz, g_0_z_xxx_xyzzzzzz, g_0_z_xxx_xzzzzzz, g_0_z_xxx_xzzzzzzz, g_0_z_xxx_yyyyyyy, g_0_z_xxx_yyyyyyz, g_0_z_xxx_yyyyyzz, g_0_z_xxx_yyyyzzz, g_0_z_xxx_yyyzzzz, g_0_z_xxx_yyzzzzz, g_0_z_xxx_yzzzzzz, g_0_z_xxx_zzzzzzz, g_0_z_xxxx_xxxxxxx, g_0_z_xxxx_xxxxxxy, g_0_z_xxxx_xxxxxxz, g_0_z_xxxx_xxxxxyy, g_0_z_xxxx_xxxxxyz, g_0_z_xxxx_xxxxxzz, g_0_z_xxxx_xxxxyyy, g_0_z_xxxx_xxxxyyz, g_0_z_xxxx_xxxxyzz, g_0_z_xxxx_xxxxzzz, g_0_z_xxxx_xxxyyyy, g_0_z_xxxx_xxxyyyz, g_0_z_xxxx_xxxyyzz, g_0_z_xxxx_xxxyzzz, g_0_z_xxxx_xxxzzzz, g_0_z_xxxx_xxyyyyy, g_0_z_xxxx_xxyyyyz, g_0_z_xxxx_xxyyyzz, g_0_z_xxxx_xxyyzzz, g_0_z_xxxx_xxyzzzz, g_0_z_xxxx_xxzzzzz, g_0_z_xxxx_xyyyyyy, g_0_z_xxxx_xyyyyyz, g_0_z_xxxx_xyyyyzz, g_0_z_xxxx_xyyyzzz, g_0_z_xxxx_xyyzzzz, g_0_z_xxxx_xyzzzzz, g_0_z_xxxx_xzzzzzz, g_0_z_xxxx_yyyyyyy, g_0_z_xxxx_yyyyyyz, g_0_z_xxxx_yyyyyzz, g_0_z_xxxx_yyyyzzz, g_0_z_xxxx_yyyzzzz, g_0_z_xxxx_yyzzzzz, g_0_z_xxxx_yzzzzzz, g_0_z_xxxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_xxxxxxx[k] = -g_0_z_xxx_xxxxxxx[k] * ab_x + g_0_z_xxx_xxxxxxxx[k];

                g_0_z_xxxx_xxxxxxy[k] = -g_0_z_xxx_xxxxxxy[k] * ab_x + g_0_z_xxx_xxxxxxxy[k];

                g_0_z_xxxx_xxxxxxz[k] = -g_0_z_xxx_xxxxxxz[k] * ab_x + g_0_z_xxx_xxxxxxxz[k];

                g_0_z_xxxx_xxxxxyy[k] = -g_0_z_xxx_xxxxxyy[k] * ab_x + g_0_z_xxx_xxxxxxyy[k];

                g_0_z_xxxx_xxxxxyz[k] = -g_0_z_xxx_xxxxxyz[k] * ab_x + g_0_z_xxx_xxxxxxyz[k];

                g_0_z_xxxx_xxxxxzz[k] = -g_0_z_xxx_xxxxxzz[k] * ab_x + g_0_z_xxx_xxxxxxzz[k];

                g_0_z_xxxx_xxxxyyy[k] = -g_0_z_xxx_xxxxyyy[k] * ab_x + g_0_z_xxx_xxxxxyyy[k];

                g_0_z_xxxx_xxxxyyz[k] = -g_0_z_xxx_xxxxyyz[k] * ab_x + g_0_z_xxx_xxxxxyyz[k];

                g_0_z_xxxx_xxxxyzz[k] = -g_0_z_xxx_xxxxyzz[k] * ab_x + g_0_z_xxx_xxxxxyzz[k];

                g_0_z_xxxx_xxxxzzz[k] = -g_0_z_xxx_xxxxzzz[k] * ab_x + g_0_z_xxx_xxxxxzzz[k];

                g_0_z_xxxx_xxxyyyy[k] = -g_0_z_xxx_xxxyyyy[k] * ab_x + g_0_z_xxx_xxxxyyyy[k];

                g_0_z_xxxx_xxxyyyz[k] = -g_0_z_xxx_xxxyyyz[k] * ab_x + g_0_z_xxx_xxxxyyyz[k];

                g_0_z_xxxx_xxxyyzz[k] = -g_0_z_xxx_xxxyyzz[k] * ab_x + g_0_z_xxx_xxxxyyzz[k];

                g_0_z_xxxx_xxxyzzz[k] = -g_0_z_xxx_xxxyzzz[k] * ab_x + g_0_z_xxx_xxxxyzzz[k];

                g_0_z_xxxx_xxxzzzz[k] = -g_0_z_xxx_xxxzzzz[k] * ab_x + g_0_z_xxx_xxxxzzzz[k];

                g_0_z_xxxx_xxyyyyy[k] = -g_0_z_xxx_xxyyyyy[k] * ab_x + g_0_z_xxx_xxxyyyyy[k];

                g_0_z_xxxx_xxyyyyz[k] = -g_0_z_xxx_xxyyyyz[k] * ab_x + g_0_z_xxx_xxxyyyyz[k];

                g_0_z_xxxx_xxyyyzz[k] = -g_0_z_xxx_xxyyyzz[k] * ab_x + g_0_z_xxx_xxxyyyzz[k];

                g_0_z_xxxx_xxyyzzz[k] = -g_0_z_xxx_xxyyzzz[k] * ab_x + g_0_z_xxx_xxxyyzzz[k];

                g_0_z_xxxx_xxyzzzz[k] = -g_0_z_xxx_xxyzzzz[k] * ab_x + g_0_z_xxx_xxxyzzzz[k];

                g_0_z_xxxx_xxzzzzz[k] = -g_0_z_xxx_xxzzzzz[k] * ab_x + g_0_z_xxx_xxxzzzzz[k];

                g_0_z_xxxx_xyyyyyy[k] = -g_0_z_xxx_xyyyyyy[k] * ab_x + g_0_z_xxx_xxyyyyyy[k];

                g_0_z_xxxx_xyyyyyz[k] = -g_0_z_xxx_xyyyyyz[k] * ab_x + g_0_z_xxx_xxyyyyyz[k];

                g_0_z_xxxx_xyyyyzz[k] = -g_0_z_xxx_xyyyyzz[k] * ab_x + g_0_z_xxx_xxyyyyzz[k];

                g_0_z_xxxx_xyyyzzz[k] = -g_0_z_xxx_xyyyzzz[k] * ab_x + g_0_z_xxx_xxyyyzzz[k];

                g_0_z_xxxx_xyyzzzz[k] = -g_0_z_xxx_xyyzzzz[k] * ab_x + g_0_z_xxx_xxyyzzzz[k];

                g_0_z_xxxx_xyzzzzz[k] = -g_0_z_xxx_xyzzzzz[k] * ab_x + g_0_z_xxx_xxyzzzzz[k];

                g_0_z_xxxx_xzzzzzz[k] = -g_0_z_xxx_xzzzzzz[k] * ab_x + g_0_z_xxx_xxzzzzzz[k];

                g_0_z_xxxx_yyyyyyy[k] = -g_0_z_xxx_yyyyyyy[k] * ab_x + g_0_z_xxx_xyyyyyyy[k];

                g_0_z_xxxx_yyyyyyz[k] = -g_0_z_xxx_yyyyyyz[k] * ab_x + g_0_z_xxx_xyyyyyyz[k];

                g_0_z_xxxx_yyyyyzz[k] = -g_0_z_xxx_yyyyyzz[k] * ab_x + g_0_z_xxx_xyyyyyzz[k];

                g_0_z_xxxx_yyyyzzz[k] = -g_0_z_xxx_yyyyzzz[k] * ab_x + g_0_z_xxx_xyyyyzzz[k];

                g_0_z_xxxx_yyyzzzz[k] = -g_0_z_xxx_yyyzzzz[k] * ab_x + g_0_z_xxx_xyyyzzzz[k];

                g_0_z_xxxx_yyzzzzz[k] = -g_0_z_xxx_yyzzzzz[k] * ab_x + g_0_z_xxx_xyyzzzzz[k];

                g_0_z_xxxx_yzzzzzz[k] = -g_0_z_xxx_yzzzzzz[k] * ab_x + g_0_z_xxx_xyzzzzzz[k];

                g_0_z_xxxx_zzzzzzz[k] = -g_0_z_xxx_zzzzzzz[k] * ab_x + g_0_z_xxx_xzzzzzzz[k];
            }

            /// Set up 1116-1152 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxxy_xxxxxxx, g_0_z_xxxy_xxxxxxy, g_0_z_xxxy_xxxxxxz, g_0_z_xxxy_xxxxxyy, g_0_z_xxxy_xxxxxyz, g_0_z_xxxy_xxxxxzz, g_0_z_xxxy_xxxxyyy, g_0_z_xxxy_xxxxyyz, g_0_z_xxxy_xxxxyzz, g_0_z_xxxy_xxxxzzz, g_0_z_xxxy_xxxyyyy, g_0_z_xxxy_xxxyyyz, g_0_z_xxxy_xxxyyzz, g_0_z_xxxy_xxxyzzz, g_0_z_xxxy_xxxzzzz, g_0_z_xxxy_xxyyyyy, g_0_z_xxxy_xxyyyyz, g_0_z_xxxy_xxyyyzz, g_0_z_xxxy_xxyyzzz, g_0_z_xxxy_xxyzzzz, g_0_z_xxxy_xxzzzzz, g_0_z_xxxy_xyyyyyy, g_0_z_xxxy_xyyyyyz, g_0_z_xxxy_xyyyyzz, g_0_z_xxxy_xyyyzzz, g_0_z_xxxy_xyyzzzz, g_0_z_xxxy_xyzzzzz, g_0_z_xxxy_xzzzzzz, g_0_z_xxxy_yyyyyyy, g_0_z_xxxy_yyyyyyz, g_0_z_xxxy_yyyyyzz, g_0_z_xxxy_yyyyzzz, g_0_z_xxxy_yyyzzzz, g_0_z_xxxy_yyzzzzz, g_0_z_xxxy_yzzzzzz, g_0_z_xxxy_zzzzzzz, g_0_z_xxy_xxxxxxx, g_0_z_xxy_xxxxxxxx, g_0_z_xxy_xxxxxxxy, g_0_z_xxy_xxxxxxxz, g_0_z_xxy_xxxxxxy, g_0_z_xxy_xxxxxxyy, g_0_z_xxy_xxxxxxyz, g_0_z_xxy_xxxxxxz, g_0_z_xxy_xxxxxxzz, g_0_z_xxy_xxxxxyy, g_0_z_xxy_xxxxxyyy, g_0_z_xxy_xxxxxyyz, g_0_z_xxy_xxxxxyz, g_0_z_xxy_xxxxxyzz, g_0_z_xxy_xxxxxzz, g_0_z_xxy_xxxxxzzz, g_0_z_xxy_xxxxyyy, g_0_z_xxy_xxxxyyyy, g_0_z_xxy_xxxxyyyz, g_0_z_xxy_xxxxyyz, g_0_z_xxy_xxxxyyzz, g_0_z_xxy_xxxxyzz, g_0_z_xxy_xxxxyzzz, g_0_z_xxy_xxxxzzz, g_0_z_xxy_xxxxzzzz, g_0_z_xxy_xxxyyyy, g_0_z_xxy_xxxyyyyy, g_0_z_xxy_xxxyyyyz, g_0_z_xxy_xxxyyyz, g_0_z_xxy_xxxyyyzz, g_0_z_xxy_xxxyyzz, g_0_z_xxy_xxxyyzzz, g_0_z_xxy_xxxyzzz, g_0_z_xxy_xxxyzzzz, g_0_z_xxy_xxxzzzz, g_0_z_xxy_xxxzzzzz, g_0_z_xxy_xxyyyyy, g_0_z_xxy_xxyyyyyy, g_0_z_xxy_xxyyyyyz, g_0_z_xxy_xxyyyyz, g_0_z_xxy_xxyyyyzz, g_0_z_xxy_xxyyyzz, g_0_z_xxy_xxyyyzzz, g_0_z_xxy_xxyyzzz, g_0_z_xxy_xxyyzzzz, g_0_z_xxy_xxyzzzz, g_0_z_xxy_xxyzzzzz, g_0_z_xxy_xxzzzzz, g_0_z_xxy_xxzzzzzz, g_0_z_xxy_xyyyyyy, g_0_z_xxy_xyyyyyyy, g_0_z_xxy_xyyyyyyz, g_0_z_xxy_xyyyyyz, g_0_z_xxy_xyyyyyzz, g_0_z_xxy_xyyyyzz, g_0_z_xxy_xyyyyzzz, g_0_z_xxy_xyyyzzz, g_0_z_xxy_xyyyzzzz, g_0_z_xxy_xyyzzzz, g_0_z_xxy_xyyzzzzz, g_0_z_xxy_xyzzzzz, g_0_z_xxy_xyzzzzzz, g_0_z_xxy_xzzzzzz, g_0_z_xxy_xzzzzzzz, g_0_z_xxy_yyyyyyy, g_0_z_xxy_yyyyyyz, g_0_z_xxy_yyyyyzz, g_0_z_xxy_yyyyzzz, g_0_z_xxy_yyyzzzz, g_0_z_xxy_yyzzzzz, g_0_z_xxy_yzzzzzz, g_0_z_xxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_xxxxxxx[k] = -g_0_z_xxy_xxxxxxx[k] * ab_x + g_0_z_xxy_xxxxxxxx[k];

                g_0_z_xxxy_xxxxxxy[k] = -g_0_z_xxy_xxxxxxy[k] * ab_x + g_0_z_xxy_xxxxxxxy[k];

                g_0_z_xxxy_xxxxxxz[k] = -g_0_z_xxy_xxxxxxz[k] * ab_x + g_0_z_xxy_xxxxxxxz[k];

                g_0_z_xxxy_xxxxxyy[k] = -g_0_z_xxy_xxxxxyy[k] * ab_x + g_0_z_xxy_xxxxxxyy[k];

                g_0_z_xxxy_xxxxxyz[k] = -g_0_z_xxy_xxxxxyz[k] * ab_x + g_0_z_xxy_xxxxxxyz[k];

                g_0_z_xxxy_xxxxxzz[k] = -g_0_z_xxy_xxxxxzz[k] * ab_x + g_0_z_xxy_xxxxxxzz[k];

                g_0_z_xxxy_xxxxyyy[k] = -g_0_z_xxy_xxxxyyy[k] * ab_x + g_0_z_xxy_xxxxxyyy[k];

                g_0_z_xxxy_xxxxyyz[k] = -g_0_z_xxy_xxxxyyz[k] * ab_x + g_0_z_xxy_xxxxxyyz[k];

                g_0_z_xxxy_xxxxyzz[k] = -g_0_z_xxy_xxxxyzz[k] * ab_x + g_0_z_xxy_xxxxxyzz[k];

                g_0_z_xxxy_xxxxzzz[k] = -g_0_z_xxy_xxxxzzz[k] * ab_x + g_0_z_xxy_xxxxxzzz[k];

                g_0_z_xxxy_xxxyyyy[k] = -g_0_z_xxy_xxxyyyy[k] * ab_x + g_0_z_xxy_xxxxyyyy[k];

                g_0_z_xxxy_xxxyyyz[k] = -g_0_z_xxy_xxxyyyz[k] * ab_x + g_0_z_xxy_xxxxyyyz[k];

                g_0_z_xxxy_xxxyyzz[k] = -g_0_z_xxy_xxxyyzz[k] * ab_x + g_0_z_xxy_xxxxyyzz[k];

                g_0_z_xxxy_xxxyzzz[k] = -g_0_z_xxy_xxxyzzz[k] * ab_x + g_0_z_xxy_xxxxyzzz[k];

                g_0_z_xxxy_xxxzzzz[k] = -g_0_z_xxy_xxxzzzz[k] * ab_x + g_0_z_xxy_xxxxzzzz[k];

                g_0_z_xxxy_xxyyyyy[k] = -g_0_z_xxy_xxyyyyy[k] * ab_x + g_0_z_xxy_xxxyyyyy[k];

                g_0_z_xxxy_xxyyyyz[k] = -g_0_z_xxy_xxyyyyz[k] * ab_x + g_0_z_xxy_xxxyyyyz[k];

                g_0_z_xxxy_xxyyyzz[k] = -g_0_z_xxy_xxyyyzz[k] * ab_x + g_0_z_xxy_xxxyyyzz[k];

                g_0_z_xxxy_xxyyzzz[k] = -g_0_z_xxy_xxyyzzz[k] * ab_x + g_0_z_xxy_xxxyyzzz[k];

                g_0_z_xxxy_xxyzzzz[k] = -g_0_z_xxy_xxyzzzz[k] * ab_x + g_0_z_xxy_xxxyzzzz[k];

                g_0_z_xxxy_xxzzzzz[k] = -g_0_z_xxy_xxzzzzz[k] * ab_x + g_0_z_xxy_xxxzzzzz[k];

                g_0_z_xxxy_xyyyyyy[k] = -g_0_z_xxy_xyyyyyy[k] * ab_x + g_0_z_xxy_xxyyyyyy[k];

                g_0_z_xxxy_xyyyyyz[k] = -g_0_z_xxy_xyyyyyz[k] * ab_x + g_0_z_xxy_xxyyyyyz[k];

                g_0_z_xxxy_xyyyyzz[k] = -g_0_z_xxy_xyyyyzz[k] * ab_x + g_0_z_xxy_xxyyyyzz[k];

                g_0_z_xxxy_xyyyzzz[k] = -g_0_z_xxy_xyyyzzz[k] * ab_x + g_0_z_xxy_xxyyyzzz[k];

                g_0_z_xxxy_xyyzzzz[k] = -g_0_z_xxy_xyyzzzz[k] * ab_x + g_0_z_xxy_xxyyzzzz[k];

                g_0_z_xxxy_xyzzzzz[k] = -g_0_z_xxy_xyzzzzz[k] * ab_x + g_0_z_xxy_xxyzzzzz[k];

                g_0_z_xxxy_xzzzzzz[k] = -g_0_z_xxy_xzzzzzz[k] * ab_x + g_0_z_xxy_xxzzzzzz[k];

                g_0_z_xxxy_yyyyyyy[k] = -g_0_z_xxy_yyyyyyy[k] * ab_x + g_0_z_xxy_xyyyyyyy[k];

                g_0_z_xxxy_yyyyyyz[k] = -g_0_z_xxy_yyyyyyz[k] * ab_x + g_0_z_xxy_xyyyyyyz[k];

                g_0_z_xxxy_yyyyyzz[k] = -g_0_z_xxy_yyyyyzz[k] * ab_x + g_0_z_xxy_xyyyyyzz[k];

                g_0_z_xxxy_yyyyzzz[k] = -g_0_z_xxy_yyyyzzz[k] * ab_x + g_0_z_xxy_xyyyyzzz[k];

                g_0_z_xxxy_yyyzzzz[k] = -g_0_z_xxy_yyyzzzz[k] * ab_x + g_0_z_xxy_xyyyzzzz[k];

                g_0_z_xxxy_yyzzzzz[k] = -g_0_z_xxy_yyzzzzz[k] * ab_x + g_0_z_xxy_xyyzzzzz[k];

                g_0_z_xxxy_yzzzzzz[k] = -g_0_z_xxy_yzzzzzz[k] * ab_x + g_0_z_xxy_xyzzzzzz[k];

                g_0_z_xxxy_zzzzzzz[k] = -g_0_z_xxy_zzzzzzz[k] * ab_x + g_0_z_xxy_xzzzzzzz[k];
            }

            /// Set up 1152-1188 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxxz_xxxxxxx, g_0_z_xxxz_xxxxxxy, g_0_z_xxxz_xxxxxxz, g_0_z_xxxz_xxxxxyy, g_0_z_xxxz_xxxxxyz, g_0_z_xxxz_xxxxxzz, g_0_z_xxxz_xxxxyyy, g_0_z_xxxz_xxxxyyz, g_0_z_xxxz_xxxxyzz, g_0_z_xxxz_xxxxzzz, g_0_z_xxxz_xxxyyyy, g_0_z_xxxz_xxxyyyz, g_0_z_xxxz_xxxyyzz, g_0_z_xxxz_xxxyzzz, g_0_z_xxxz_xxxzzzz, g_0_z_xxxz_xxyyyyy, g_0_z_xxxz_xxyyyyz, g_0_z_xxxz_xxyyyzz, g_0_z_xxxz_xxyyzzz, g_0_z_xxxz_xxyzzzz, g_0_z_xxxz_xxzzzzz, g_0_z_xxxz_xyyyyyy, g_0_z_xxxz_xyyyyyz, g_0_z_xxxz_xyyyyzz, g_0_z_xxxz_xyyyzzz, g_0_z_xxxz_xyyzzzz, g_0_z_xxxz_xyzzzzz, g_0_z_xxxz_xzzzzzz, g_0_z_xxxz_yyyyyyy, g_0_z_xxxz_yyyyyyz, g_0_z_xxxz_yyyyyzz, g_0_z_xxxz_yyyyzzz, g_0_z_xxxz_yyyzzzz, g_0_z_xxxz_yyzzzzz, g_0_z_xxxz_yzzzzzz, g_0_z_xxxz_zzzzzzz, g_0_z_xxz_xxxxxxx, g_0_z_xxz_xxxxxxxx, g_0_z_xxz_xxxxxxxy, g_0_z_xxz_xxxxxxxz, g_0_z_xxz_xxxxxxy, g_0_z_xxz_xxxxxxyy, g_0_z_xxz_xxxxxxyz, g_0_z_xxz_xxxxxxz, g_0_z_xxz_xxxxxxzz, g_0_z_xxz_xxxxxyy, g_0_z_xxz_xxxxxyyy, g_0_z_xxz_xxxxxyyz, g_0_z_xxz_xxxxxyz, g_0_z_xxz_xxxxxyzz, g_0_z_xxz_xxxxxzz, g_0_z_xxz_xxxxxzzz, g_0_z_xxz_xxxxyyy, g_0_z_xxz_xxxxyyyy, g_0_z_xxz_xxxxyyyz, g_0_z_xxz_xxxxyyz, g_0_z_xxz_xxxxyyzz, g_0_z_xxz_xxxxyzz, g_0_z_xxz_xxxxyzzz, g_0_z_xxz_xxxxzzz, g_0_z_xxz_xxxxzzzz, g_0_z_xxz_xxxyyyy, g_0_z_xxz_xxxyyyyy, g_0_z_xxz_xxxyyyyz, g_0_z_xxz_xxxyyyz, g_0_z_xxz_xxxyyyzz, g_0_z_xxz_xxxyyzz, g_0_z_xxz_xxxyyzzz, g_0_z_xxz_xxxyzzz, g_0_z_xxz_xxxyzzzz, g_0_z_xxz_xxxzzzz, g_0_z_xxz_xxxzzzzz, g_0_z_xxz_xxyyyyy, g_0_z_xxz_xxyyyyyy, g_0_z_xxz_xxyyyyyz, g_0_z_xxz_xxyyyyz, g_0_z_xxz_xxyyyyzz, g_0_z_xxz_xxyyyzz, g_0_z_xxz_xxyyyzzz, g_0_z_xxz_xxyyzzz, g_0_z_xxz_xxyyzzzz, g_0_z_xxz_xxyzzzz, g_0_z_xxz_xxyzzzzz, g_0_z_xxz_xxzzzzz, g_0_z_xxz_xxzzzzzz, g_0_z_xxz_xyyyyyy, g_0_z_xxz_xyyyyyyy, g_0_z_xxz_xyyyyyyz, g_0_z_xxz_xyyyyyz, g_0_z_xxz_xyyyyyzz, g_0_z_xxz_xyyyyzz, g_0_z_xxz_xyyyyzzz, g_0_z_xxz_xyyyzzz, g_0_z_xxz_xyyyzzzz, g_0_z_xxz_xyyzzzz, g_0_z_xxz_xyyzzzzz, g_0_z_xxz_xyzzzzz, g_0_z_xxz_xyzzzzzz, g_0_z_xxz_xzzzzzz, g_0_z_xxz_xzzzzzzz, g_0_z_xxz_yyyyyyy, g_0_z_xxz_yyyyyyz, g_0_z_xxz_yyyyyzz, g_0_z_xxz_yyyyzzz, g_0_z_xxz_yyyzzzz, g_0_z_xxz_yyzzzzz, g_0_z_xxz_yzzzzzz, g_0_z_xxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_xxxxxxx[k] = -g_0_z_xxz_xxxxxxx[k] * ab_x + g_0_z_xxz_xxxxxxxx[k];

                g_0_z_xxxz_xxxxxxy[k] = -g_0_z_xxz_xxxxxxy[k] * ab_x + g_0_z_xxz_xxxxxxxy[k];

                g_0_z_xxxz_xxxxxxz[k] = -g_0_z_xxz_xxxxxxz[k] * ab_x + g_0_z_xxz_xxxxxxxz[k];

                g_0_z_xxxz_xxxxxyy[k] = -g_0_z_xxz_xxxxxyy[k] * ab_x + g_0_z_xxz_xxxxxxyy[k];

                g_0_z_xxxz_xxxxxyz[k] = -g_0_z_xxz_xxxxxyz[k] * ab_x + g_0_z_xxz_xxxxxxyz[k];

                g_0_z_xxxz_xxxxxzz[k] = -g_0_z_xxz_xxxxxzz[k] * ab_x + g_0_z_xxz_xxxxxxzz[k];

                g_0_z_xxxz_xxxxyyy[k] = -g_0_z_xxz_xxxxyyy[k] * ab_x + g_0_z_xxz_xxxxxyyy[k];

                g_0_z_xxxz_xxxxyyz[k] = -g_0_z_xxz_xxxxyyz[k] * ab_x + g_0_z_xxz_xxxxxyyz[k];

                g_0_z_xxxz_xxxxyzz[k] = -g_0_z_xxz_xxxxyzz[k] * ab_x + g_0_z_xxz_xxxxxyzz[k];

                g_0_z_xxxz_xxxxzzz[k] = -g_0_z_xxz_xxxxzzz[k] * ab_x + g_0_z_xxz_xxxxxzzz[k];

                g_0_z_xxxz_xxxyyyy[k] = -g_0_z_xxz_xxxyyyy[k] * ab_x + g_0_z_xxz_xxxxyyyy[k];

                g_0_z_xxxz_xxxyyyz[k] = -g_0_z_xxz_xxxyyyz[k] * ab_x + g_0_z_xxz_xxxxyyyz[k];

                g_0_z_xxxz_xxxyyzz[k] = -g_0_z_xxz_xxxyyzz[k] * ab_x + g_0_z_xxz_xxxxyyzz[k];

                g_0_z_xxxz_xxxyzzz[k] = -g_0_z_xxz_xxxyzzz[k] * ab_x + g_0_z_xxz_xxxxyzzz[k];

                g_0_z_xxxz_xxxzzzz[k] = -g_0_z_xxz_xxxzzzz[k] * ab_x + g_0_z_xxz_xxxxzzzz[k];

                g_0_z_xxxz_xxyyyyy[k] = -g_0_z_xxz_xxyyyyy[k] * ab_x + g_0_z_xxz_xxxyyyyy[k];

                g_0_z_xxxz_xxyyyyz[k] = -g_0_z_xxz_xxyyyyz[k] * ab_x + g_0_z_xxz_xxxyyyyz[k];

                g_0_z_xxxz_xxyyyzz[k] = -g_0_z_xxz_xxyyyzz[k] * ab_x + g_0_z_xxz_xxxyyyzz[k];

                g_0_z_xxxz_xxyyzzz[k] = -g_0_z_xxz_xxyyzzz[k] * ab_x + g_0_z_xxz_xxxyyzzz[k];

                g_0_z_xxxz_xxyzzzz[k] = -g_0_z_xxz_xxyzzzz[k] * ab_x + g_0_z_xxz_xxxyzzzz[k];

                g_0_z_xxxz_xxzzzzz[k] = -g_0_z_xxz_xxzzzzz[k] * ab_x + g_0_z_xxz_xxxzzzzz[k];

                g_0_z_xxxz_xyyyyyy[k] = -g_0_z_xxz_xyyyyyy[k] * ab_x + g_0_z_xxz_xxyyyyyy[k];

                g_0_z_xxxz_xyyyyyz[k] = -g_0_z_xxz_xyyyyyz[k] * ab_x + g_0_z_xxz_xxyyyyyz[k];

                g_0_z_xxxz_xyyyyzz[k] = -g_0_z_xxz_xyyyyzz[k] * ab_x + g_0_z_xxz_xxyyyyzz[k];

                g_0_z_xxxz_xyyyzzz[k] = -g_0_z_xxz_xyyyzzz[k] * ab_x + g_0_z_xxz_xxyyyzzz[k];

                g_0_z_xxxz_xyyzzzz[k] = -g_0_z_xxz_xyyzzzz[k] * ab_x + g_0_z_xxz_xxyyzzzz[k];

                g_0_z_xxxz_xyzzzzz[k] = -g_0_z_xxz_xyzzzzz[k] * ab_x + g_0_z_xxz_xxyzzzzz[k];

                g_0_z_xxxz_xzzzzzz[k] = -g_0_z_xxz_xzzzzzz[k] * ab_x + g_0_z_xxz_xxzzzzzz[k];

                g_0_z_xxxz_yyyyyyy[k] = -g_0_z_xxz_yyyyyyy[k] * ab_x + g_0_z_xxz_xyyyyyyy[k];

                g_0_z_xxxz_yyyyyyz[k] = -g_0_z_xxz_yyyyyyz[k] * ab_x + g_0_z_xxz_xyyyyyyz[k];

                g_0_z_xxxz_yyyyyzz[k] = -g_0_z_xxz_yyyyyzz[k] * ab_x + g_0_z_xxz_xyyyyyzz[k];

                g_0_z_xxxz_yyyyzzz[k] = -g_0_z_xxz_yyyyzzz[k] * ab_x + g_0_z_xxz_xyyyyzzz[k];

                g_0_z_xxxz_yyyzzzz[k] = -g_0_z_xxz_yyyzzzz[k] * ab_x + g_0_z_xxz_xyyyzzzz[k];

                g_0_z_xxxz_yyzzzzz[k] = -g_0_z_xxz_yyzzzzz[k] * ab_x + g_0_z_xxz_xyyzzzzz[k];

                g_0_z_xxxz_yzzzzzz[k] = -g_0_z_xxz_yzzzzzz[k] * ab_x + g_0_z_xxz_xyzzzzzz[k];

                g_0_z_xxxz_zzzzzzz[k] = -g_0_z_xxz_zzzzzzz[k] * ab_x + g_0_z_xxz_xzzzzzzz[k];
            }

            /// Set up 1188-1224 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxyy_xxxxxxx, g_0_z_xxyy_xxxxxxy, g_0_z_xxyy_xxxxxxz, g_0_z_xxyy_xxxxxyy, g_0_z_xxyy_xxxxxyz, g_0_z_xxyy_xxxxxzz, g_0_z_xxyy_xxxxyyy, g_0_z_xxyy_xxxxyyz, g_0_z_xxyy_xxxxyzz, g_0_z_xxyy_xxxxzzz, g_0_z_xxyy_xxxyyyy, g_0_z_xxyy_xxxyyyz, g_0_z_xxyy_xxxyyzz, g_0_z_xxyy_xxxyzzz, g_0_z_xxyy_xxxzzzz, g_0_z_xxyy_xxyyyyy, g_0_z_xxyy_xxyyyyz, g_0_z_xxyy_xxyyyzz, g_0_z_xxyy_xxyyzzz, g_0_z_xxyy_xxyzzzz, g_0_z_xxyy_xxzzzzz, g_0_z_xxyy_xyyyyyy, g_0_z_xxyy_xyyyyyz, g_0_z_xxyy_xyyyyzz, g_0_z_xxyy_xyyyzzz, g_0_z_xxyy_xyyzzzz, g_0_z_xxyy_xyzzzzz, g_0_z_xxyy_xzzzzzz, g_0_z_xxyy_yyyyyyy, g_0_z_xxyy_yyyyyyz, g_0_z_xxyy_yyyyyzz, g_0_z_xxyy_yyyyzzz, g_0_z_xxyy_yyyzzzz, g_0_z_xxyy_yyzzzzz, g_0_z_xxyy_yzzzzzz, g_0_z_xxyy_zzzzzzz, g_0_z_xyy_xxxxxxx, g_0_z_xyy_xxxxxxxx, g_0_z_xyy_xxxxxxxy, g_0_z_xyy_xxxxxxxz, g_0_z_xyy_xxxxxxy, g_0_z_xyy_xxxxxxyy, g_0_z_xyy_xxxxxxyz, g_0_z_xyy_xxxxxxz, g_0_z_xyy_xxxxxxzz, g_0_z_xyy_xxxxxyy, g_0_z_xyy_xxxxxyyy, g_0_z_xyy_xxxxxyyz, g_0_z_xyy_xxxxxyz, g_0_z_xyy_xxxxxyzz, g_0_z_xyy_xxxxxzz, g_0_z_xyy_xxxxxzzz, g_0_z_xyy_xxxxyyy, g_0_z_xyy_xxxxyyyy, g_0_z_xyy_xxxxyyyz, g_0_z_xyy_xxxxyyz, g_0_z_xyy_xxxxyyzz, g_0_z_xyy_xxxxyzz, g_0_z_xyy_xxxxyzzz, g_0_z_xyy_xxxxzzz, g_0_z_xyy_xxxxzzzz, g_0_z_xyy_xxxyyyy, g_0_z_xyy_xxxyyyyy, g_0_z_xyy_xxxyyyyz, g_0_z_xyy_xxxyyyz, g_0_z_xyy_xxxyyyzz, g_0_z_xyy_xxxyyzz, g_0_z_xyy_xxxyyzzz, g_0_z_xyy_xxxyzzz, g_0_z_xyy_xxxyzzzz, g_0_z_xyy_xxxzzzz, g_0_z_xyy_xxxzzzzz, g_0_z_xyy_xxyyyyy, g_0_z_xyy_xxyyyyyy, g_0_z_xyy_xxyyyyyz, g_0_z_xyy_xxyyyyz, g_0_z_xyy_xxyyyyzz, g_0_z_xyy_xxyyyzz, g_0_z_xyy_xxyyyzzz, g_0_z_xyy_xxyyzzz, g_0_z_xyy_xxyyzzzz, g_0_z_xyy_xxyzzzz, g_0_z_xyy_xxyzzzzz, g_0_z_xyy_xxzzzzz, g_0_z_xyy_xxzzzzzz, g_0_z_xyy_xyyyyyy, g_0_z_xyy_xyyyyyyy, g_0_z_xyy_xyyyyyyz, g_0_z_xyy_xyyyyyz, g_0_z_xyy_xyyyyyzz, g_0_z_xyy_xyyyyzz, g_0_z_xyy_xyyyyzzz, g_0_z_xyy_xyyyzzz, g_0_z_xyy_xyyyzzzz, g_0_z_xyy_xyyzzzz, g_0_z_xyy_xyyzzzzz, g_0_z_xyy_xyzzzzz, g_0_z_xyy_xyzzzzzz, g_0_z_xyy_xzzzzzz, g_0_z_xyy_xzzzzzzz, g_0_z_xyy_yyyyyyy, g_0_z_xyy_yyyyyyz, g_0_z_xyy_yyyyyzz, g_0_z_xyy_yyyyzzz, g_0_z_xyy_yyyzzzz, g_0_z_xyy_yyzzzzz, g_0_z_xyy_yzzzzzz, g_0_z_xyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_xxxxxxx[k] = -g_0_z_xyy_xxxxxxx[k] * ab_x + g_0_z_xyy_xxxxxxxx[k];

                g_0_z_xxyy_xxxxxxy[k] = -g_0_z_xyy_xxxxxxy[k] * ab_x + g_0_z_xyy_xxxxxxxy[k];

                g_0_z_xxyy_xxxxxxz[k] = -g_0_z_xyy_xxxxxxz[k] * ab_x + g_0_z_xyy_xxxxxxxz[k];

                g_0_z_xxyy_xxxxxyy[k] = -g_0_z_xyy_xxxxxyy[k] * ab_x + g_0_z_xyy_xxxxxxyy[k];

                g_0_z_xxyy_xxxxxyz[k] = -g_0_z_xyy_xxxxxyz[k] * ab_x + g_0_z_xyy_xxxxxxyz[k];

                g_0_z_xxyy_xxxxxzz[k] = -g_0_z_xyy_xxxxxzz[k] * ab_x + g_0_z_xyy_xxxxxxzz[k];

                g_0_z_xxyy_xxxxyyy[k] = -g_0_z_xyy_xxxxyyy[k] * ab_x + g_0_z_xyy_xxxxxyyy[k];

                g_0_z_xxyy_xxxxyyz[k] = -g_0_z_xyy_xxxxyyz[k] * ab_x + g_0_z_xyy_xxxxxyyz[k];

                g_0_z_xxyy_xxxxyzz[k] = -g_0_z_xyy_xxxxyzz[k] * ab_x + g_0_z_xyy_xxxxxyzz[k];

                g_0_z_xxyy_xxxxzzz[k] = -g_0_z_xyy_xxxxzzz[k] * ab_x + g_0_z_xyy_xxxxxzzz[k];

                g_0_z_xxyy_xxxyyyy[k] = -g_0_z_xyy_xxxyyyy[k] * ab_x + g_0_z_xyy_xxxxyyyy[k];

                g_0_z_xxyy_xxxyyyz[k] = -g_0_z_xyy_xxxyyyz[k] * ab_x + g_0_z_xyy_xxxxyyyz[k];

                g_0_z_xxyy_xxxyyzz[k] = -g_0_z_xyy_xxxyyzz[k] * ab_x + g_0_z_xyy_xxxxyyzz[k];

                g_0_z_xxyy_xxxyzzz[k] = -g_0_z_xyy_xxxyzzz[k] * ab_x + g_0_z_xyy_xxxxyzzz[k];

                g_0_z_xxyy_xxxzzzz[k] = -g_0_z_xyy_xxxzzzz[k] * ab_x + g_0_z_xyy_xxxxzzzz[k];

                g_0_z_xxyy_xxyyyyy[k] = -g_0_z_xyy_xxyyyyy[k] * ab_x + g_0_z_xyy_xxxyyyyy[k];

                g_0_z_xxyy_xxyyyyz[k] = -g_0_z_xyy_xxyyyyz[k] * ab_x + g_0_z_xyy_xxxyyyyz[k];

                g_0_z_xxyy_xxyyyzz[k] = -g_0_z_xyy_xxyyyzz[k] * ab_x + g_0_z_xyy_xxxyyyzz[k];

                g_0_z_xxyy_xxyyzzz[k] = -g_0_z_xyy_xxyyzzz[k] * ab_x + g_0_z_xyy_xxxyyzzz[k];

                g_0_z_xxyy_xxyzzzz[k] = -g_0_z_xyy_xxyzzzz[k] * ab_x + g_0_z_xyy_xxxyzzzz[k];

                g_0_z_xxyy_xxzzzzz[k] = -g_0_z_xyy_xxzzzzz[k] * ab_x + g_0_z_xyy_xxxzzzzz[k];

                g_0_z_xxyy_xyyyyyy[k] = -g_0_z_xyy_xyyyyyy[k] * ab_x + g_0_z_xyy_xxyyyyyy[k];

                g_0_z_xxyy_xyyyyyz[k] = -g_0_z_xyy_xyyyyyz[k] * ab_x + g_0_z_xyy_xxyyyyyz[k];

                g_0_z_xxyy_xyyyyzz[k] = -g_0_z_xyy_xyyyyzz[k] * ab_x + g_0_z_xyy_xxyyyyzz[k];

                g_0_z_xxyy_xyyyzzz[k] = -g_0_z_xyy_xyyyzzz[k] * ab_x + g_0_z_xyy_xxyyyzzz[k];

                g_0_z_xxyy_xyyzzzz[k] = -g_0_z_xyy_xyyzzzz[k] * ab_x + g_0_z_xyy_xxyyzzzz[k];

                g_0_z_xxyy_xyzzzzz[k] = -g_0_z_xyy_xyzzzzz[k] * ab_x + g_0_z_xyy_xxyzzzzz[k];

                g_0_z_xxyy_xzzzzzz[k] = -g_0_z_xyy_xzzzzzz[k] * ab_x + g_0_z_xyy_xxzzzzzz[k];

                g_0_z_xxyy_yyyyyyy[k] = -g_0_z_xyy_yyyyyyy[k] * ab_x + g_0_z_xyy_xyyyyyyy[k];

                g_0_z_xxyy_yyyyyyz[k] = -g_0_z_xyy_yyyyyyz[k] * ab_x + g_0_z_xyy_xyyyyyyz[k];

                g_0_z_xxyy_yyyyyzz[k] = -g_0_z_xyy_yyyyyzz[k] * ab_x + g_0_z_xyy_xyyyyyzz[k];

                g_0_z_xxyy_yyyyzzz[k] = -g_0_z_xyy_yyyyzzz[k] * ab_x + g_0_z_xyy_xyyyyzzz[k];

                g_0_z_xxyy_yyyzzzz[k] = -g_0_z_xyy_yyyzzzz[k] * ab_x + g_0_z_xyy_xyyyzzzz[k];

                g_0_z_xxyy_yyzzzzz[k] = -g_0_z_xyy_yyzzzzz[k] * ab_x + g_0_z_xyy_xyyzzzzz[k];

                g_0_z_xxyy_yzzzzzz[k] = -g_0_z_xyy_yzzzzzz[k] * ab_x + g_0_z_xyy_xyzzzzzz[k];

                g_0_z_xxyy_zzzzzzz[k] = -g_0_z_xyy_zzzzzzz[k] * ab_x + g_0_z_xyy_xzzzzzzz[k];
            }

            /// Set up 1224-1260 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxyz_xxxxxxx, g_0_z_xxyz_xxxxxxy, g_0_z_xxyz_xxxxxxz, g_0_z_xxyz_xxxxxyy, g_0_z_xxyz_xxxxxyz, g_0_z_xxyz_xxxxxzz, g_0_z_xxyz_xxxxyyy, g_0_z_xxyz_xxxxyyz, g_0_z_xxyz_xxxxyzz, g_0_z_xxyz_xxxxzzz, g_0_z_xxyz_xxxyyyy, g_0_z_xxyz_xxxyyyz, g_0_z_xxyz_xxxyyzz, g_0_z_xxyz_xxxyzzz, g_0_z_xxyz_xxxzzzz, g_0_z_xxyz_xxyyyyy, g_0_z_xxyz_xxyyyyz, g_0_z_xxyz_xxyyyzz, g_0_z_xxyz_xxyyzzz, g_0_z_xxyz_xxyzzzz, g_0_z_xxyz_xxzzzzz, g_0_z_xxyz_xyyyyyy, g_0_z_xxyz_xyyyyyz, g_0_z_xxyz_xyyyyzz, g_0_z_xxyz_xyyyzzz, g_0_z_xxyz_xyyzzzz, g_0_z_xxyz_xyzzzzz, g_0_z_xxyz_xzzzzzz, g_0_z_xxyz_yyyyyyy, g_0_z_xxyz_yyyyyyz, g_0_z_xxyz_yyyyyzz, g_0_z_xxyz_yyyyzzz, g_0_z_xxyz_yyyzzzz, g_0_z_xxyz_yyzzzzz, g_0_z_xxyz_yzzzzzz, g_0_z_xxyz_zzzzzzz, g_0_z_xyz_xxxxxxx, g_0_z_xyz_xxxxxxxx, g_0_z_xyz_xxxxxxxy, g_0_z_xyz_xxxxxxxz, g_0_z_xyz_xxxxxxy, g_0_z_xyz_xxxxxxyy, g_0_z_xyz_xxxxxxyz, g_0_z_xyz_xxxxxxz, g_0_z_xyz_xxxxxxzz, g_0_z_xyz_xxxxxyy, g_0_z_xyz_xxxxxyyy, g_0_z_xyz_xxxxxyyz, g_0_z_xyz_xxxxxyz, g_0_z_xyz_xxxxxyzz, g_0_z_xyz_xxxxxzz, g_0_z_xyz_xxxxxzzz, g_0_z_xyz_xxxxyyy, g_0_z_xyz_xxxxyyyy, g_0_z_xyz_xxxxyyyz, g_0_z_xyz_xxxxyyz, g_0_z_xyz_xxxxyyzz, g_0_z_xyz_xxxxyzz, g_0_z_xyz_xxxxyzzz, g_0_z_xyz_xxxxzzz, g_0_z_xyz_xxxxzzzz, g_0_z_xyz_xxxyyyy, g_0_z_xyz_xxxyyyyy, g_0_z_xyz_xxxyyyyz, g_0_z_xyz_xxxyyyz, g_0_z_xyz_xxxyyyzz, g_0_z_xyz_xxxyyzz, g_0_z_xyz_xxxyyzzz, g_0_z_xyz_xxxyzzz, g_0_z_xyz_xxxyzzzz, g_0_z_xyz_xxxzzzz, g_0_z_xyz_xxxzzzzz, g_0_z_xyz_xxyyyyy, g_0_z_xyz_xxyyyyyy, g_0_z_xyz_xxyyyyyz, g_0_z_xyz_xxyyyyz, g_0_z_xyz_xxyyyyzz, g_0_z_xyz_xxyyyzz, g_0_z_xyz_xxyyyzzz, g_0_z_xyz_xxyyzzz, g_0_z_xyz_xxyyzzzz, g_0_z_xyz_xxyzzzz, g_0_z_xyz_xxyzzzzz, g_0_z_xyz_xxzzzzz, g_0_z_xyz_xxzzzzzz, g_0_z_xyz_xyyyyyy, g_0_z_xyz_xyyyyyyy, g_0_z_xyz_xyyyyyyz, g_0_z_xyz_xyyyyyz, g_0_z_xyz_xyyyyyzz, g_0_z_xyz_xyyyyzz, g_0_z_xyz_xyyyyzzz, g_0_z_xyz_xyyyzzz, g_0_z_xyz_xyyyzzzz, g_0_z_xyz_xyyzzzz, g_0_z_xyz_xyyzzzzz, g_0_z_xyz_xyzzzzz, g_0_z_xyz_xyzzzzzz, g_0_z_xyz_xzzzzzz, g_0_z_xyz_xzzzzzzz, g_0_z_xyz_yyyyyyy, g_0_z_xyz_yyyyyyz, g_0_z_xyz_yyyyyzz, g_0_z_xyz_yyyyzzz, g_0_z_xyz_yyyzzzz, g_0_z_xyz_yyzzzzz, g_0_z_xyz_yzzzzzz, g_0_z_xyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_xxxxxxx[k] = -g_0_z_xyz_xxxxxxx[k] * ab_x + g_0_z_xyz_xxxxxxxx[k];

                g_0_z_xxyz_xxxxxxy[k] = -g_0_z_xyz_xxxxxxy[k] * ab_x + g_0_z_xyz_xxxxxxxy[k];

                g_0_z_xxyz_xxxxxxz[k] = -g_0_z_xyz_xxxxxxz[k] * ab_x + g_0_z_xyz_xxxxxxxz[k];

                g_0_z_xxyz_xxxxxyy[k] = -g_0_z_xyz_xxxxxyy[k] * ab_x + g_0_z_xyz_xxxxxxyy[k];

                g_0_z_xxyz_xxxxxyz[k] = -g_0_z_xyz_xxxxxyz[k] * ab_x + g_0_z_xyz_xxxxxxyz[k];

                g_0_z_xxyz_xxxxxzz[k] = -g_0_z_xyz_xxxxxzz[k] * ab_x + g_0_z_xyz_xxxxxxzz[k];

                g_0_z_xxyz_xxxxyyy[k] = -g_0_z_xyz_xxxxyyy[k] * ab_x + g_0_z_xyz_xxxxxyyy[k];

                g_0_z_xxyz_xxxxyyz[k] = -g_0_z_xyz_xxxxyyz[k] * ab_x + g_0_z_xyz_xxxxxyyz[k];

                g_0_z_xxyz_xxxxyzz[k] = -g_0_z_xyz_xxxxyzz[k] * ab_x + g_0_z_xyz_xxxxxyzz[k];

                g_0_z_xxyz_xxxxzzz[k] = -g_0_z_xyz_xxxxzzz[k] * ab_x + g_0_z_xyz_xxxxxzzz[k];

                g_0_z_xxyz_xxxyyyy[k] = -g_0_z_xyz_xxxyyyy[k] * ab_x + g_0_z_xyz_xxxxyyyy[k];

                g_0_z_xxyz_xxxyyyz[k] = -g_0_z_xyz_xxxyyyz[k] * ab_x + g_0_z_xyz_xxxxyyyz[k];

                g_0_z_xxyz_xxxyyzz[k] = -g_0_z_xyz_xxxyyzz[k] * ab_x + g_0_z_xyz_xxxxyyzz[k];

                g_0_z_xxyz_xxxyzzz[k] = -g_0_z_xyz_xxxyzzz[k] * ab_x + g_0_z_xyz_xxxxyzzz[k];

                g_0_z_xxyz_xxxzzzz[k] = -g_0_z_xyz_xxxzzzz[k] * ab_x + g_0_z_xyz_xxxxzzzz[k];

                g_0_z_xxyz_xxyyyyy[k] = -g_0_z_xyz_xxyyyyy[k] * ab_x + g_0_z_xyz_xxxyyyyy[k];

                g_0_z_xxyz_xxyyyyz[k] = -g_0_z_xyz_xxyyyyz[k] * ab_x + g_0_z_xyz_xxxyyyyz[k];

                g_0_z_xxyz_xxyyyzz[k] = -g_0_z_xyz_xxyyyzz[k] * ab_x + g_0_z_xyz_xxxyyyzz[k];

                g_0_z_xxyz_xxyyzzz[k] = -g_0_z_xyz_xxyyzzz[k] * ab_x + g_0_z_xyz_xxxyyzzz[k];

                g_0_z_xxyz_xxyzzzz[k] = -g_0_z_xyz_xxyzzzz[k] * ab_x + g_0_z_xyz_xxxyzzzz[k];

                g_0_z_xxyz_xxzzzzz[k] = -g_0_z_xyz_xxzzzzz[k] * ab_x + g_0_z_xyz_xxxzzzzz[k];

                g_0_z_xxyz_xyyyyyy[k] = -g_0_z_xyz_xyyyyyy[k] * ab_x + g_0_z_xyz_xxyyyyyy[k];

                g_0_z_xxyz_xyyyyyz[k] = -g_0_z_xyz_xyyyyyz[k] * ab_x + g_0_z_xyz_xxyyyyyz[k];

                g_0_z_xxyz_xyyyyzz[k] = -g_0_z_xyz_xyyyyzz[k] * ab_x + g_0_z_xyz_xxyyyyzz[k];

                g_0_z_xxyz_xyyyzzz[k] = -g_0_z_xyz_xyyyzzz[k] * ab_x + g_0_z_xyz_xxyyyzzz[k];

                g_0_z_xxyz_xyyzzzz[k] = -g_0_z_xyz_xyyzzzz[k] * ab_x + g_0_z_xyz_xxyyzzzz[k];

                g_0_z_xxyz_xyzzzzz[k] = -g_0_z_xyz_xyzzzzz[k] * ab_x + g_0_z_xyz_xxyzzzzz[k];

                g_0_z_xxyz_xzzzzzz[k] = -g_0_z_xyz_xzzzzzz[k] * ab_x + g_0_z_xyz_xxzzzzzz[k];

                g_0_z_xxyz_yyyyyyy[k] = -g_0_z_xyz_yyyyyyy[k] * ab_x + g_0_z_xyz_xyyyyyyy[k];

                g_0_z_xxyz_yyyyyyz[k] = -g_0_z_xyz_yyyyyyz[k] * ab_x + g_0_z_xyz_xyyyyyyz[k];

                g_0_z_xxyz_yyyyyzz[k] = -g_0_z_xyz_yyyyyzz[k] * ab_x + g_0_z_xyz_xyyyyyzz[k];

                g_0_z_xxyz_yyyyzzz[k] = -g_0_z_xyz_yyyyzzz[k] * ab_x + g_0_z_xyz_xyyyyzzz[k];

                g_0_z_xxyz_yyyzzzz[k] = -g_0_z_xyz_yyyzzzz[k] * ab_x + g_0_z_xyz_xyyyzzzz[k];

                g_0_z_xxyz_yyzzzzz[k] = -g_0_z_xyz_yyzzzzz[k] * ab_x + g_0_z_xyz_xyyzzzzz[k];

                g_0_z_xxyz_yzzzzzz[k] = -g_0_z_xyz_yzzzzzz[k] * ab_x + g_0_z_xyz_xyzzzzzz[k];

                g_0_z_xxyz_zzzzzzz[k] = -g_0_z_xyz_zzzzzzz[k] * ab_x + g_0_z_xyz_xzzzzzzz[k];
            }

            /// Set up 1260-1296 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxzz_xxxxxxx, g_0_z_xxzz_xxxxxxy, g_0_z_xxzz_xxxxxxz, g_0_z_xxzz_xxxxxyy, g_0_z_xxzz_xxxxxyz, g_0_z_xxzz_xxxxxzz, g_0_z_xxzz_xxxxyyy, g_0_z_xxzz_xxxxyyz, g_0_z_xxzz_xxxxyzz, g_0_z_xxzz_xxxxzzz, g_0_z_xxzz_xxxyyyy, g_0_z_xxzz_xxxyyyz, g_0_z_xxzz_xxxyyzz, g_0_z_xxzz_xxxyzzz, g_0_z_xxzz_xxxzzzz, g_0_z_xxzz_xxyyyyy, g_0_z_xxzz_xxyyyyz, g_0_z_xxzz_xxyyyzz, g_0_z_xxzz_xxyyzzz, g_0_z_xxzz_xxyzzzz, g_0_z_xxzz_xxzzzzz, g_0_z_xxzz_xyyyyyy, g_0_z_xxzz_xyyyyyz, g_0_z_xxzz_xyyyyzz, g_0_z_xxzz_xyyyzzz, g_0_z_xxzz_xyyzzzz, g_0_z_xxzz_xyzzzzz, g_0_z_xxzz_xzzzzzz, g_0_z_xxzz_yyyyyyy, g_0_z_xxzz_yyyyyyz, g_0_z_xxzz_yyyyyzz, g_0_z_xxzz_yyyyzzz, g_0_z_xxzz_yyyzzzz, g_0_z_xxzz_yyzzzzz, g_0_z_xxzz_yzzzzzz, g_0_z_xxzz_zzzzzzz, g_0_z_xzz_xxxxxxx, g_0_z_xzz_xxxxxxxx, g_0_z_xzz_xxxxxxxy, g_0_z_xzz_xxxxxxxz, g_0_z_xzz_xxxxxxy, g_0_z_xzz_xxxxxxyy, g_0_z_xzz_xxxxxxyz, g_0_z_xzz_xxxxxxz, g_0_z_xzz_xxxxxxzz, g_0_z_xzz_xxxxxyy, g_0_z_xzz_xxxxxyyy, g_0_z_xzz_xxxxxyyz, g_0_z_xzz_xxxxxyz, g_0_z_xzz_xxxxxyzz, g_0_z_xzz_xxxxxzz, g_0_z_xzz_xxxxxzzz, g_0_z_xzz_xxxxyyy, g_0_z_xzz_xxxxyyyy, g_0_z_xzz_xxxxyyyz, g_0_z_xzz_xxxxyyz, g_0_z_xzz_xxxxyyzz, g_0_z_xzz_xxxxyzz, g_0_z_xzz_xxxxyzzz, g_0_z_xzz_xxxxzzz, g_0_z_xzz_xxxxzzzz, g_0_z_xzz_xxxyyyy, g_0_z_xzz_xxxyyyyy, g_0_z_xzz_xxxyyyyz, g_0_z_xzz_xxxyyyz, g_0_z_xzz_xxxyyyzz, g_0_z_xzz_xxxyyzz, g_0_z_xzz_xxxyyzzz, g_0_z_xzz_xxxyzzz, g_0_z_xzz_xxxyzzzz, g_0_z_xzz_xxxzzzz, g_0_z_xzz_xxxzzzzz, g_0_z_xzz_xxyyyyy, g_0_z_xzz_xxyyyyyy, g_0_z_xzz_xxyyyyyz, g_0_z_xzz_xxyyyyz, g_0_z_xzz_xxyyyyzz, g_0_z_xzz_xxyyyzz, g_0_z_xzz_xxyyyzzz, g_0_z_xzz_xxyyzzz, g_0_z_xzz_xxyyzzzz, g_0_z_xzz_xxyzzzz, g_0_z_xzz_xxyzzzzz, g_0_z_xzz_xxzzzzz, g_0_z_xzz_xxzzzzzz, g_0_z_xzz_xyyyyyy, g_0_z_xzz_xyyyyyyy, g_0_z_xzz_xyyyyyyz, g_0_z_xzz_xyyyyyz, g_0_z_xzz_xyyyyyzz, g_0_z_xzz_xyyyyzz, g_0_z_xzz_xyyyyzzz, g_0_z_xzz_xyyyzzz, g_0_z_xzz_xyyyzzzz, g_0_z_xzz_xyyzzzz, g_0_z_xzz_xyyzzzzz, g_0_z_xzz_xyzzzzz, g_0_z_xzz_xyzzzzzz, g_0_z_xzz_xzzzzzz, g_0_z_xzz_xzzzzzzz, g_0_z_xzz_yyyyyyy, g_0_z_xzz_yyyyyyz, g_0_z_xzz_yyyyyzz, g_0_z_xzz_yyyyzzz, g_0_z_xzz_yyyzzzz, g_0_z_xzz_yyzzzzz, g_0_z_xzz_yzzzzzz, g_0_z_xzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_xxxxxxx[k] = -g_0_z_xzz_xxxxxxx[k] * ab_x + g_0_z_xzz_xxxxxxxx[k];

                g_0_z_xxzz_xxxxxxy[k] = -g_0_z_xzz_xxxxxxy[k] * ab_x + g_0_z_xzz_xxxxxxxy[k];

                g_0_z_xxzz_xxxxxxz[k] = -g_0_z_xzz_xxxxxxz[k] * ab_x + g_0_z_xzz_xxxxxxxz[k];

                g_0_z_xxzz_xxxxxyy[k] = -g_0_z_xzz_xxxxxyy[k] * ab_x + g_0_z_xzz_xxxxxxyy[k];

                g_0_z_xxzz_xxxxxyz[k] = -g_0_z_xzz_xxxxxyz[k] * ab_x + g_0_z_xzz_xxxxxxyz[k];

                g_0_z_xxzz_xxxxxzz[k] = -g_0_z_xzz_xxxxxzz[k] * ab_x + g_0_z_xzz_xxxxxxzz[k];

                g_0_z_xxzz_xxxxyyy[k] = -g_0_z_xzz_xxxxyyy[k] * ab_x + g_0_z_xzz_xxxxxyyy[k];

                g_0_z_xxzz_xxxxyyz[k] = -g_0_z_xzz_xxxxyyz[k] * ab_x + g_0_z_xzz_xxxxxyyz[k];

                g_0_z_xxzz_xxxxyzz[k] = -g_0_z_xzz_xxxxyzz[k] * ab_x + g_0_z_xzz_xxxxxyzz[k];

                g_0_z_xxzz_xxxxzzz[k] = -g_0_z_xzz_xxxxzzz[k] * ab_x + g_0_z_xzz_xxxxxzzz[k];

                g_0_z_xxzz_xxxyyyy[k] = -g_0_z_xzz_xxxyyyy[k] * ab_x + g_0_z_xzz_xxxxyyyy[k];

                g_0_z_xxzz_xxxyyyz[k] = -g_0_z_xzz_xxxyyyz[k] * ab_x + g_0_z_xzz_xxxxyyyz[k];

                g_0_z_xxzz_xxxyyzz[k] = -g_0_z_xzz_xxxyyzz[k] * ab_x + g_0_z_xzz_xxxxyyzz[k];

                g_0_z_xxzz_xxxyzzz[k] = -g_0_z_xzz_xxxyzzz[k] * ab_x + g_0_z_xzz_xxxxyzzz[k];

                g_0_z_xxzz_xxxzzzz[k] = -g_0_z_xzz_xxxzzzz[k] * ab_x + g_0_z_xzz_xxxxzzzz[k];

                g_0_z_xxzz_xxyyyyy[k] = -g_0_z_xzz_xxyyyyy[k] * ab_x + g_0_z_xzz_xxxyyyyy[k];

                g_0_z_xxzz_xxyyyyz[k] = -g_0_z_xzz_xxyyyyz[k] * ab_x + g_0_z_xzz_xxxyyyyz[k];

                g_0_z_xxzz_xxyyyzz[k] = -g_0_z_xzz_xxyyyzz[k] * ab_x + g_0_z_xzz_xxxyyyzz[k];

                g_0_z_xxzz_xxyyzzz[k] = -g_0_z_xzz_xxyyzzz[k] * ab_x + g_0_z_xzz_xxxyyzzz[k];

                g_0_z_xxzz_xxyzzzz[k] = -g_0_z_xzz_xxyzzzz[k] * ab_x + g_0_z_xzz_xxxyzzzz[k];

                g_0_z_xxzz_xxzzzzz[k] = -g_0_z_xzz_xxzzzzz[k] * ab_x + g_0_z_xzz_xxxzzzzz[k];

                g_0_z_xxzz_xyyyyyy[k] = -g_0_z_xzz_xyyyyyy[k] * ab_x + g_0_z_xzz_xxyyyyyy[k];

                g_0_z_xxzz_xyyyyyz[k] = -g_0_z_xzz_xyyyyyz[k] * ab_x + g_0_z_xzz_xxyyyyyz[k];

                g_0_z_xxzz_xyyyyzz[k] = -g_0_z_xzz_xyyyyzz[k] * ab_x + g_0_z_xzz_xxyyyyzz[k];

                g_0_z_xxzz_xyyyzzz[k] = -g_0_z_xzz_xyyyzzz[k] * ab_x + g_0_z_xzz_xxyyyzzz[k];

                g_0_z_xxzz_xyyzzzz[k] = -g_0_z_xzz_xyyzzzz[k] * ab_x + g_0_z_xzz_xxyyzzzz[k];

                g_0_z_xxzz_xyzzzzz[k] = -g_0_z_xzz_xyzzzzz[k] * ab_x + g_0_z_xzz_xxyzzzzz[k];

                g_0_z_xxzz_xzzzzzz[k] = -g_0_z_xzz_xzzzzzz[k] * ab_x + g_0_z_xzz_xxzzzzzz[k];

                g_0_z_xxzz_yyyyyyy[k] = -g_0_z_xzz_yyyyyyy[k] * ab_x + g_0_z_xzz_xyyyyyyy[k];

                g_0_z_xxzz_yyyyyyz[k] = -g_0_z_xzz_yyyyyyz[k] * ab_x + g_0_z_xzz_xyyyyyyz[k];

                g_0_z_xxzz_yyyyyzz[k] = -g_0_z_xzz_yyyyyzz[k] * ab_x + g_0_z_xzz_xyyyyyzz[k];

                g_0_z_xxzz_yyyyzzz[k] = -g_0_z_xzz_yyyyzzz[k] * ab_x + g_0_z_xzz_xyyyyzzz[k];

                g_0_z_xxzz_yyyzzzz[k] = -g_0_z_xzz_yyyzzzz[k] * ab_x + g_0_z_xzz_xyyyzzzz[k];

                g_0_z_xxzz_yyzzzzz[k] = -g_0_z_xzz_yyzzzzz[k] * ab_x + g_0_z_xzz_xyyzzzzz[k];

                g_0_z_xxzz_yzzzzzz[k] = -g_0_z_xzz_yzzzzzz[k] * ab_x + g_0_z_xzz_xyzzzzzz[k];

                g_0_z_xxzz_zzzzzzz[k] = -g_0_z_xzz_zzzzzzz[k] * ab_x + g_0_z_xzz_xzzzzzzz[k];
            }

            /// Set up 1296-1332 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyyy_xxxxxxx, g_0_z_xyyy_xxxxxxy, g_0_z_xyyy_xxxxxxz, g_0_z_xyyy_xxxxxyy, g_0_z_xyyy_xxxxxyz, g_0_z_xyyy_xxxxxzz, g_0_z_xyyy_xxxxyyy, g_0_z_xyyy_xxxxyyz, g_0_z_xyyy_xxxxyzz, g_0_z_xyyy_xxxxzzz, g_0_z_xyyy_xxxyyyy, g_0_z_xyyy_xxxyyyz, g_0_z_xyyy_xxxyyzz, g_0_z_xyyy_xxxyzzz, g_0_z_xyyy_xxxzzzz, g_0_z_xyyy_xxyyyyy, g_0_z_xyyy_xxyyyyz, g_0_z_xyyy_xxyyyzz, g_0_z_xyyy_xxyyzzz, g_0_z_xyyy_xxyzzzz, g_0_z_xyyy_xxzzzzz, g_0_z_xyyy_xyyyyyy, g_0_z_xyyy_xyyyyyz, g_0_z_xyyy_xyyyyzz, g_0_z_xyyy_xyyyzzz, g_0_z_xyyy_xyyzzzz, g_0_z_xyyy_xyzzzzz, g_0_z_xyyy_xzzzzzz, g_0_z_xyyy_yyyyyyy, g_0_z_xyyy_yyyyyyz, g_0_z_xyyy_yyyyyzz, g_0_z_xyyy_yyyyzzz, g_0_z_xyyy_yyyzzzz, g_0_z_xyyy_yyzzzzz, g_0_z_xyyy_yzzzzzz, g_0_z_xyyy_zzzzzzz, g_0_z_yyy_xxxxxxx, g_0_z_yyy_xxxxxxxx, g_0_z_yyy_xxxxxxxy, g_0_z_yyy_xxxxxxxz, g_0_z_yyy_xxxxxxy, g_0_z_yyy_xxxxxxyy, g_0_z_yyy_xxxxxxyz, g_0_z_yyy_xxxxxxz, g_0_z_yyy_xxxxxxzz, g_0_z_yyy_xxxxxyy, g_0_z_yyy_xxxxxyyy, g_0_z_yyy_xxxxxyyz, g_0_z_yyy_xxxxxyz, g_0_z_yyy_xxxxxyzz, g_0_z_yyy_xxxxxzz, g_0_z_yyy_xxxxxzzz, g_0_z_yyy_xxxxyyy, g_0_z_yyy_xxxxyyyy, g_0_z_yyy_xxxxyyyz, g_0_z_yyy_xxxxyyz, g_0_z_yyy_xxxxyyzz, g_0_z_yyy_xxxxyzz, g_0_z_yyy_xxxxyzzz, g_0_z_yyy_xxxxzzz, g_0_z_yyy_xxxxzzzz, g_0_z_yyy_xxxyyyy, g_0_z_yyy_xxxyyyyy, g_0_z_yyy_xxxyyyyz, g_0_z_yyy_xxxyyyz, g_0_z_yyy_xxxyyyzz, g_0_z_yyy_xxxyyzz, g_0_z_yyy_xxxyyzzz, g_0_z_yyy_xxxyzzz, g_0_z_yyy_xxxyzzzz, g_0_z_yyy_xxxzzzz, g_0_z_yyy_xxxzzzzz, g_0_z_yyy_xxyyyyy, g_0_z_yyy_xxyyyyyy, g_0_z_yyy_xxyyyyyz, g_0_z_yyy_xxyyyyz, g_0_z_yyy_xxyyyyzz, g_0_z_yyy_xxyyyzz, g_0_z_yyy_xxyyyzzz, g_0_z_yyy_xxyyzzz, g_0_z_yyy_xxyyzzzz, g_0_z_yyy_xxyzzzz, g_0_z_yyy_xxyzzzzz, g_0_z_yyy_xxzzzzz, g_0_z_yyy_xxzzzzzz, g_0_z_yyy_xyyyyyy, g_0_z_yyy_xyyyyyyy, g_0_z_yyy_xyyyyyyz, g_0_z_yyy_xyyyyyz, g_0_z_yyy_xyyyyyzz, g_0_z_yyy_xyyyyzz, g_0_z_yyy_xyyyyzzz, g_0_z_yyy_xyyyzzz, g_0_z_yyy_xyyyzzzz, g_0_z_yyy_xyyzzzz, g_0_z_yyy_xyyzzzzz, g_0_z_yyy_xyzzzzz, g_0_z_yyy_xyzzzzzz, g_0_z_yyy_xzzzzzz, g_0_z_yyy_xzzzzzzz, g_0_z_yyy_yyyyyyy, g_0_z_yyy_yyyyyyz, g_0_z_yyy_yyyyyzz, g_0_z_yyy_yyyyzzz, g_0_z_yyy_yyyzzzz, g_0_z_yyy_yyzzzzz, g_0_z_yyy_yzzzzzz, g_0_z_yyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_xxxxxxx[k] = -g_0_z_yyy_xxxxxxx[k] * ab_x + g_0_z_yyy_xxxxxxxx[k];

                g_0_z_xyyy_xxxxxxy[k] = -g_0_z_yyy_xxxxxxy[k] * ab_x + g_0_z_yyy_xxxxxxxy[k];

                g_0_z_xyyy_xxxxxxz[k] = -g_0_z_yyy_xxxxxxz[k] * ab_x + g_0_z_yyy_xxxxxxxz[k];

                g_0_z_xyyy_xxxxxyy[k] = -g_0_z_yyy_xxxxxyy[k] * ab_x + g_0_z_yyy_xxxxxxyy[k];

                g_0_z_xyyy_xxxxxyz[k] = -g_0_z_yyy_xxxxxyz[k] * ab_x + g_0_z_yyy_xxxxxxyz[k];

                g_0_z_xyyy_xxxxxzz[k] = -g_0_z_yyy_xxxxxzz[k] * ab_x + g_0_z_yyy_xxxxxxzz[k];

                g_0_z_xyyy_xxxxyyy[k] = -g_0_z_yyy_xxxxyyy[k] * ab_x + g_0_z_yyy_xxxxxyyy[k];

                g_0_z_xyyy_xxxxyyz[k] = -g_0_z_yyy_xxxxyyz[k] * ab_x + g_0_z_yyy_xxxxxyyz[k];

                g_0_z_xyyy_xxxxyzz[k] = -g_0_z_yyy_xxxxyzz[k] * ab_x + g_0_z_yyy_xxxxxyzz[k];

                g_0_z_xyyy_xxxxzzz[k] = -g_0_z_yyy_xxxxzzz[k] * ab_x + g_0_z_yyy_xxxxxzzz[k];

                g_0_z_xyyy_xxxyyyy[k] = -g_0_z_yyy_xxxyyyy[k] * ab_x + g_0_z_yyy_xxxxyyyy[k];

                g_0_z_xyyy_xxxyyyz[k] = -g_0_z_yyy_xxxyyyz[k] * ab_x + g_0_z_yyy_xxxxyyyz[k];

                g_0_z_xyyy_xxxyyzz[k] = -g_0_z_yyy_xxxyyzz[k] * ab_x + g_0_z_yyy_xxxxyyzz[k];

                g_0_z_xyyy_xxxyzzz[k] = -g_0_z_yyy_xxxyzzz[k] * ab_x + g_0_z_yyy_xxxxyzzz[k];

                g_0_z_xyyy_xxxzzzz[k] = -g_0_z_yyy_xxxzzzz[k] * ab_x + g_0_z_yyy_xxxxzzzz[k];

                g_0_z_xyyy_xxyyyyy[k] = -g_0_z_yyy_xxyyyyy[k] * ab_x + g_0_z_yyy_xxxyyyyy[k];

                g_0_z_xyyy_xxyyyyz[k] = -g_0_z_yyy_xxyyyyz[k] * ab_x + g_0_z_yyy_xxxyyyyz[k];

                g_0_z_xyyy_xxyyyzz[k] = -g_0_z_yyy_xxyyyzz[k] * ab_x + g_0_z_yyy_xxxyyyzz[k];

                g_0_z_xyyy_xxyyzzz[k] = -g_0_z_yyy_xxyyzzz[k] * ab_x + g_0_z_yyy_xxxyyzzz[k];

                g_0_z_xyyy_xxyzzzz[k] = -g_0_z_yyy_xxyzzzz[k] * ab_x + g_0_z_yyy_xxxyzzzz[k];

                g_0_z_xyyy_xxzzzzz[k] = -g_0_z_yyy_xxzzzzz[k] * ab_x + g_0_z_yyy_xxxzzzzz[k];

                g_0_z_xyyy_xyyyyyy[k] = -g_0_z_yyy_xyyyyyy[k] * ab_x + g_0_z_yyy_xxyyyyyy[k];

                g_0_z_xyyy_xyyyyyz[k] = -g_0_z_yyy_xyyyyyz[k] * ab_x + g_0_z_yyy_xxyyyyyz[k];

                g_0_z_xyyy_xyyyyzz[k] = -g_0_z_yyy_xyyyyzz[k] * ab_x + g_0_z_yyy_xxyyyyzz[k];

                g_0_z_xyyy_xyyyzzz[k] = -g_0_z_yyy_xyyyzzz[k] * ab_x + g_0_z_yyy_xxyyyzzz[k];

                g_0_z_xyyy_xyyzzzz[k] = -g_0_z_yyy_xyyzzzz[k] * ab_x + g_0_z_yyy_xxyyzzzz[k];

                g_0_z_xyyy_xyzzzzz[k] = -g_0_z_yyy_xyzzzzz[k] * ab_x + g_0_z_yyy_xxyzzzzz[k];

                g_0_z_xyyy_xzzzzzz[k] = -g_0_z_yyy_xzzzzzz[k] * ab_x + g_0_z_yyy_xxzzzzzz[k];

                g_0_z_xyyy_yyyyyyy[k] = -g_0_z_yyy_yyyyyyy[k] * ab_x + g_0_z_yyy_xyyyyyyy[k];

                g_0_z_xyyy_yyyyyyz[k] = -g_0_z_yyy_yyyyyyz[k] * ab_x + g_0_z_yyy_xyyyyyyz[k];

                g_0_z_xyyy_yyyyyzz[k] = -g_0_z_yyy_yyyyyzz[k] * ab_x + g_0_z_yyy_xyyyyyzz[k];

                g_0_z_xyyy_yyyyzzz[k] = -g_0_z_yyy_yyyyzzz[k] * ab_x + g_0_z_yyy_xyyyyzzz[k];

                g_0_z_xyyy_yyyzzzz[k] = -g_0_z_yyy_yyyzzzz[k] * ab_x + g_0_z_yyy_xyyyzzzz[k];

                g_0_z_xyyy_yyzzzzz[k] = -g_0_z_yyy_yyzzzzz[k] * ab_x + g_0_z_yyy_xyyzzzzz[k];

                g_0_z_xyyy_yzzzzzz[k] = -g_0_z_yyy_yzzzzzz[k] * ab_x + g_0_z_yyy_xyzzzzzz[k];

                g_0_z_xyyy_zzzzzzz[k] = -g_0_z_yyy_zzzzzzz[k] * ab_x + g_0_z_yyy_xzzzzzzz[k];
            }

            /// Set up 1332-1368 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyyz_xxxxxxx, g_0_z_xyyz_xxxxxxy, g_0_z_xyyz_xxxxxxz, g_0_z_xyyz_xxxxxyy, g_0_z_xyyz_xxxxxyz, g_0_z_xyyz_xxxxxzz, g_0_z_xyyz_xxxxyyy, g_0_z_xyyz_xxxxyyz, g_0_z_xyyz_xxxxyzz, g_0_z_xyyz_xxxxzzz, g_0_z_xyyz_xxxyyyy, g_0_z_xyyz_xxxyyyz, g_0_z_xyyz_xxxyyzz, g_0_z_xyyz_xxxyzzz, g_0_z_xyyz_xxxzzzz, g_0_z_xyyz_xxyyyyy, g_0_z_xyyz_xxyyyyz, g_0_z_xyyz_xxyyyzz, g_0_z_xyyz_xxyyzzz, g_0_z_xyyz_xxyzzzz, g_0_z_xyyz_xxzzzzz, g_0_z_xyyz_xyyyyyy, g_0_z_xyyz_xyyyyyz, g_0_z_xyyz_xyyyyzz, g_0_z_xyyz_xyyyzzz, g_0_z_xyyz_xyyzzzz, g_0_z_xyyz_xyzzzzz, g_0_z_xyyz_xzzzzzz, g_0_z_xyyz_yyyyyyy, g_0_z_xyyz_yyyyyyz, g_0_z_xyyz_yyyyyzz, g_0_z_xyyz_yyyyzzz, g_0_z_xyyz_yyyzzzz, g_0_z_xyyz_yyzzzzz, g_0_z_xyyz_yzzzzzz, g_0_z_xyyz_zzzzzzz, g_0_z_yyz_xxxxxxx, g_0_z_yyz_xxxxxxxx, g_0_z_yyz_xxxxxxxy, g_0_z_yyz_xxxxxxxz, g_0_z_yyz_xxxxxxy, g_0_z_yyz_xxxxxxyy, g_0_z_yyz_xxxxxxyz, g_0_z_yyz_xxxxxxz, g_0_z_yyz_xxxxxxzz, g_0_z_yyz_xxxxxyy, g_0_z_yyz_xxxxxyyy, g_0_z_yyz_xxxxxyyz, g_0_z_yyz_xxxxxyz, g_0_z_yyz_xxxxxyzz, g_0_z_yyz_xxxxxzz, g_0_z_yyz_xxxxxzzz, g_0_z_yyz_xxxxyyy, g_0_z_yyz_xxxxyyyy, g_0_z_yyz_xxxxyyyz, g_0_z_yyz_xxxxyyz, g_0_z_yyz_xxxxyyzz, g_0_z_yyz_xxxxyzz, g_0_z_yyz_xxxxyzzz, g_0_z_yyz_xxxxzzz, g_0_z_yyz_xxxxzzzz, g_0_z_yyz_xxxyyyy, g_0_z_yyz_xxxyyyyy, g_0_z_yyz_xxxyyyyz, g_0_z_yyz_xxxyyyz, g_0_z_yyz_xxxyyyzz, g_0_z_yyz_xxxyyzz, g_0_z_yyz_xxxyyzzz, g_0_z_yyz_xxxyzzz, g_0_z_yyz_xxxyzzzz, g_0_z_yyz_xxxzzzz, g_0_z_yyz_xxxzzzzz, g_0_z_yyz_xxyyyyy, g_0_z_yyz_xxyyyyyy, g_0_z_yyz_xxyyyyyz, g_0_z_yyz_xxyyyyz, g_0_z_yyz_xxyyyyzz, g_0_z_yyz_xxyyyzz, g_0_z_yyz_xxyyyzzz, g_0_z_yyz_xxyyzzz, g_0_z_yyz_xxyyzzzz, g_0_z_yyz_xxyzzzz, g_0_z_yyz_xxyzzzzz, g_0_z_yyz_xxzzzzz, g_0_z_yyz_xxzzzzzz, g_0_z_yyz_xyyyyyy, g_0_z_yyz_xyyyyyyy, g_0_z_yyz_xyyyyyyz, g_0_z_yyz_xyyyyyz, g_0_z_yyz_xyyyyyzz, g_0_z_yyz_xyyyyzz, g_0_z_yyz_xyyyyzzz, g_0_z_yyz_xyyyzzz, g_0_z_yyz_xyyyzzzz, g_0_z_yyz_xyyzzzz, g_0_z_yyz_xyyzzzzz, g_0_z_yyz_xyzzzzz, g_0_z_yyz_xyzzzzzz, g_0_z_yyz_xzzzzzz, g_0_z_yyz_xzzzzzzz, g_0_z_yyz_yyyyyyy, g_0_z_yyz_yyyyyyz, g_0_z_yyz_yyyyyzz, g_0_z_yyz_yyyyzzz, g_0_z_yyz_yyyzzzz, g_0_z_yyz_yyzzzzz, g_0_z_yyz_yzzzzzz, g_0_z_yyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_xxxxxxx[k] = -g_0_z_yyz_xxxxxxx[k] * ab_x + g_0_z_yyz_xxxxxxxx[k];

                g_0_z_xyyz_xxxxxxy[k] = -g_0_z_yyz_xxxxxxy[k] * ab_x + g_0_z_yyz_xxxxxxxy[k];

                g_0_z_xyyz_xxxxxxz[k] = -g_0_z_yyz_xxxxxxz[k] * ab_x + g_0_z_yyz_xxxxxxxz[k];

                g_0_z_xyyz_xxxxxyy[k] = -g_0_z_yyz_xxxxxyy[k] * ab_x + g_0_z_yyz_xxxxxxyy[k];

                g_0_z_xyyz_xxxxxyz[k] = -g_0_z_yyz_xxxxxyz[k] * ab_x + g_0_z_yyz_xxxxxxyz[k];

                g_0_z_xyyz_xxxxxzz[k] = -g_0_z_yyz_xxxxxzz[k] * ab_x + g_0_z_yyz_xxxxxxzz[k];

                g_0_z_xyyz_xxxxyyy[k] = -g_0_z_yyz_xxxxyyy[k] * ab_x + g_0_z_yyz_xxxxxyyy[k];

                g_0_z_xyyz_xxxxyyz[k] = -g_0_z_yyz_xxxxyyz[k] * ab_x + g_0_z_yyz_xxxxxyyz[k];

                g_0_z_xyyz_xxxxyzz[k] = -g_0_z_yyz_xxxxyzz[k] * ab_x + g_0_z_yyz_xxxxxyzz[k];

                g_0_z_xyyz_xxxxzzz[k] = -g_0_z_yyz_xxxxzzz[k] * ab_x + g_0_z_yyz_xxxxxzzz[k];

                g_0_z_xyyz_xxxyyyy[k] = -g_0_z_yyz_xxxyyyy[k] * ab_x + g_0_z_yyz_xxxxyyyy[k];

                g_0_z_xyyz_xxxyyyz[k] = -g_0_z_yyz_xxxyyyz[k] * ab_x + g_0_z_yyz_xxxxyyyz[k];

                g_0_z_xyyz_xxxyyzz[k] = -g_0_z_yyz_xxxyyzz[k] * ab_x + g_0_z_yyz_xxxxyyzz[k];

                g_0_z_xyyz_xxxyzzz[k] = -g_0_z_yyz_xxxyzzz[k] * ab_x + g_0_z_yyz_xxxxyzzz[k];

                g_0_z_xyyz_xxxzzzz[k] = -g_0_z_yyz_xxxzzzz[k] * ab_x + g_0_z_yyz_xxxxzzzz[k];

                g_0_z_xyyz_xxyyyyy[k] = -g_0_z_yyz_xxyyyyy[k] * ab_x + g_0_z_yyz_xxxyyyyy[k];

                g_0_z_xyyz_xxyyyyz[k] = -g_0_z_yyz_xxyyyyz[k] * ab_x + g_0_z_yyz_xxxyyyyz[k];

                g_0_z_xyyz_xxyyyzz[k] = -g_0_z_yyz_xxyyyzz[k] * ab_x + g_0_z_yyz_xxxyyyzz[k];

                g_0_z_xyyz_xxyyzzz[k] = -g_0_z_yyz_xxyyzzz[k] * ab_x + g_0_z_yyz_xxxyyzzz[k];

                g_0_z_xyyz_xxyzzzz[k] = -g_0_z_yyz_xxyzzzz[k] * ab_x + g_0_z_yyz_xxxyzzzz[k];

                g_0_z_xyyz_xxzzzzz[k] = -g_0_z_yyz_xxzzzzz[k] * ab_x + g_0_z_yyz_xxxzzzzz[k];

                g_0_z_xyyz_xyyyyyy[k] = -g_0_z_yyz_xyyyyyy[k] * ab_x + g_0_z_yyz_xxyyyyyy[k];

                g_0_z_xyyz_xyyyyyz[k] = -g_0_z_yyz_xyyyyyz[k] * ab_x + g_0_z_yyz_xxyyyyyz[k];

                g_0_z_xyyz_xyyyyzz[k] = -g_0_z_yyz_xyyyyzz[k] * ab_x + g_0_z_yyz_xxyyyyzz[k];

                g_0_z_xyyz_xyyyzzz[k] = -g_0_z_yyz_xyyyzzz[k] * ab_x + g_0_z_yyz_xxyyyzzz[k];

                g_0_z_xyyz_xyyzzzz[k] = -g_0_z_yyz_xyyzzzz[k] * ab_x + g_0_z_yyz_xxyyzzzz[k];

                g_0_z_xyyz_xyzzzzz[k] = -g_0_z_yyz_xyzzzzz[k] * ab_x + g_0_z_yyz_xxyzzzzz[k];

                g_0_z_xyyz_xzzzzzz[k] = -g_0_z_yyz_xzzzzzz[k] * ab_x + g_0_z_yyz_xxzzzzzz[k];

                g_0_z_xyyz_yyyyyyy[k] = -g_0_z_yyz_yyyyyyy[k] * ab_x + g_0_z_yyz_xyyyyyyy[k];

                g_0_z_xyyz_yyyyyyz[k] = -g_0_z_yyz_yyyyyyz[k] * ab_x + g_0_z_yyz_xyyyyyyz[k];

                g_0_z_xyyz_yyyyyzz[k] = -g_0_z_yyz_yyyyyzz[k] * ab_x + g_0_z_yyz_xyyyyyzz[k];

                g_0_z_xyyz_yyyyzzz[k] = -g_0_z_yyz_yyyyzzz[k] * ab_x + g_0_z_yyz_xyyyyzzz[k];

                g_0_z_xyyz_yyyzzzz[k] = -g_0_z_yyz_yyyzzzz[k] * ab_x + g_0_z_yyz_xyyyzzzz[k];

                g_0_z_xyyz_yyzzzzz[k] = -g_0_z_yyz_yyzzzzz[k] * ab_x + g_0_z_yyz_xyyzzzzz[k];

                g_0_z_xyyz_yzzzzzz[k] = -g_0_z_yyz_yzzzzzz[k] * ab_x + g_0_z_yyz_xyzzzzzz[k];

                g_0_z_xyyz_zzzzzzz[k] = -g_0_z_yyz_zzzzzzz[k] * ab_x + g_0_z_yyz_xzzzzzzz[k];
            }

            /// Set up 1368-1404 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyzz_xxxxxxx, g_0_z_xyzz_xxxxxxy, g_0_z_xyzz_xxxxxxz, g_0_z_xyzz_xxxxxyy, g_0_z_xyzz_xxxxxyz, g_0_z_xyzz_xxxxxzz, g_0_z_xyzz_xxxxyyy, g_0_z_xyzz_xxxxyyz, g_0_z_xyzz_xxxxyzz, g_0_z_xyzz_xxxxzzz, g_0_z_xyzz_xxxyyyy, g_0_z_xyzz_xxxyyyz, g_0_z_xyzz_xxxyyzz, g_0_z_xyzz_xxxyzzz, g_0_z_xyzz_xxxzzzz, g_0_z_xyzz_xxyyyyy, g_0_z_xyzz_xxyyyyz, g_0_z_xyzz_xxyyyzz, g_0_z_xyzz_xxyyzzz, g_0_z_xyzz_xxyzzzz, g_0_z_xyzz_xxzzzzz, g_0_z_xyzz_xyyyyyy, g_0_z_xyzz_xyyyyyz, g_0_z_xyzz_xyyyyzz, g_0_z_xyzz_xyyyzzz, g_0_z_xyzz_xyyzzzz, g_0_z_xyzz_xyzzzzz, g_0_z_xyzz_xzzzzzz, g_0_z_xyzz_yyyyyyy, g_0_z_xyzz_yyyyyyz, g_0_z_xyzz_yyyyyzz, g_0_z_xyzz_yyyyzzz, g_0_z_xyzz_yyyzzzz, g_0_z_xyzz_yyzzzzz, g_0_z_xyzz_yzzzzzz, g_0_z_xyzz_zzzzzzz, g_0_z_yzz_xxxxxxx, g_0_z_yzz_xxxxxxxx, g_0_z_yzz_xxxxxxxy, g_0_z_yzz_xxxxxxxz, g_0_z_yzz_xxxxxxy, g_0_z_yzz_xxxxxxyy, g_0_z_yzz_xxxxxxyz, g_0_z_yzz_xxxxxxz, g_0_z_yzz_xxxxxxzz, g_0_z_yzz_xxxxxyy, g_0_z_yzz_xxxxxyyy, g_0_z_yzz_xxxxxyyz, g_0_z_yzz_xxxxxyz, g_0_z_yzz_xxxxxyzz, g_0_z_yzz_xxxxxzz, g_0_z_yzz_xxxxxzzz, g_0_z_yzz_xxxxyyy, g_0_z_yzz_xxxxyyyy, g_0_z_yzz_xxxxyyyz, g_0_z_yzz_xxxxyyz, g_0_z_yzz_xxxxyyzz, g_0_z_yzz_xxxxyzz, g_0_z_yzz_xxxxyzzz, g_0_z_yzz_xxxxzzz, g_0_z_yzz_xxxxzzzz, g_0_z_yzz_xxxyyyy, g_0_z_yzz_xxxyyyyy, g_0_z_yzz_xxxyyyyz, g_0_z_yzz_xxxyyyz, g_0_z_yzz_xxxyyyzz, g_0_z_yzz_xxxyyzz, g_0_z_yzz_xxxyyzzz, g_0_z_yzz_xxxyzzz, g_0_z_yzz_xxxyzzzz, g_0_z_yzz_xxxzzzz, g_0_z_yzz_xxxzzzzz, g_0_z_yzz_xxyyyyy, g_0_z_yzz_xxyyyyyy, g_0_z_yzz_xxyyyyyz, g_0_z_yzz_xxyyyyz, g_0_z_yzz_xxyyyyzz, g_0_z_yzz_xxyyyzz, g_0_z_yzz_xxyyyzzz, g_0_z_yzz_xxyyzzz, g_0_z_yzz_xxyyzzzz, g_0_z_yzz_xxyzzzz, g_0_z_yzz_xxyzzzzz, g_0_z_yzz_xxzzzzz, g_0_z_yzz_xxzzzzzz, g_0_z_yzz_xyyyyyy, g_0_z_yzz_xyyyyyyy, g_0_z_yzz_xyyyyyyz, g_0_z_yzz_xyyyyyz, g_0_z_yzz_xyyyyyzz, g_0_z_yzz_xyyyyzz, g_0_z_yzz_xyyyyzzz, g_0_z_yzz_xyyyzzz, g_0_z_yzz_xyyyzzzz, g_0_z_yzz_xyyzzzz, g_0_z_yzz_xyyzzzzz, g_0_z_yzz_xyzzzzz, g_0_z_yzz_xyzzzzzz, g_0_z_yzz_xzzzzzz, g_0_z_yzz_xzzzzzzz, g_0_z_yzz_yyyyyyy, g_0_z_yzz_yyyyyyz, g_0_z_yzz_yyyyyzz, g_0_z_yzz_yyyyzzz, g_0_z_yzz_yyyzzzz, g_0_z_yzz_yyzzzzz, g_0_z_yzz_yzzzzzz, g_0_z_yzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_xxxxxxx[k] = -g_0_z_yzz_xxxxxxx[k] * ab_x + g_0_z_yzz_xxxxxxxx[k];

                g_0_z_xyzz_xxxxxxy[k] = -g_0_z_yzz_xxxxxxy[k] * ab_x + g_0_z_yzz_xxxxxxxy[k];

                g_0_z_xyzz_xxxxxxz[k] = -g_0_z_yzz_xxxxxxz[k] * ab_x + g_0_z_yzz_xxxxxxxz[k];

                g_0_z_xyzz_xxxxxyy[k] = -g_0_z_yzz_xxxxxyy[k] * ab_x + g_0_z_yzz_xxxxxxyy[k];

                g_0_z_xyzz_xxxxxyz[k] = -g_0_z_yzz_xxxxxyz[k] * ab_x + g_0_z_yzz_xxxxxxyz[k];

                g_0_z_xyzz_xxxxxzz[k] = -g_0_z_yzz_xxxxxzz[k] * ab_x + g_0_z_yzz_xxxxxxzz[k];

                g_0_z_xyzz_xxxxyyy[k] = -g_0_z_yzz_xxxxyyy[k] * ab_x + g_0_z_yzz_xxxxxyyy[k];

                g_0_z_xyzz_xxxxyyz[k] = -g_0_z_yzz_xxxxyyz[k] * ab_x + g_0_z_yzz_xxxxxyyz[k];

                g_0_z_xyzz_xxxxyzz[k] = -g_0_z_yzz_xxxxyzz[k] * ab_x + g_0_z_yzz_xxxxxyzz[k];

                g_0_z_xyzz_xxxxzzz[k] = -g_0_z_yzz_xxxxzzz[k] * ab_x + g_0_z_yzz_xxxxxzzz[k];

                g_0_z_xyzz_xxxyyyy[k] = -g_0_z_yzz_xxxyyyy[k] * ab_x + g_0_z_yzz_xxxxyyyy[k];

                g_0_z_xyzz_xxxyyyz[k] = -g_0_z_yzz_xxxyyyz[k] * ab_x + g_0_z_yzz_xxxxyyyz[k];

                g_0_z_xyzz_xxxyyzz[k] = -g_0_z_yzz_xxxyyzz[k] * ab_x + g_0_z_yzz_xxxxyyzz[k];

                g_0_z_xyzz_xxxyzzz[k] = -g_0_z_yzz_xxxyzzz[k] * ab_x + g_0_z_yzz_xxxxyzzz[k];

                g_0_z_xyzz_xxxzzzz[k] = -g_0_z_yzz_xxxzzzz[k] * ab_x + g_0_z_yzz_xxxxzzzz[k];

                g_0_z_xyzz_xxyyyyy[k] = -g_0_z_yzz_xxyyyyy[k] * ab_x + g_0_z_yzz_xxxyyyyy[k];

                g_0_z_xyzz_xxyyyyz[k] = -g_0_z_yzz_xxyyyyz[k] * ab_x + g_0_z_yzz_xxxyyyyz[k];

                g_0_z_xyzz_xxyyyzz[k] = -g_0_z_yzz_xxyyyzz[k] * ab_x + g_0_z_yzz_xxxyyyzz[k];

                g_0_z_xyzz_xxyyzzz[k] = -g_0_z_yzz_xxyyzzz[k] * ab_x + g_0_z_yzz_xxxyyzzz[k];

                g_0_z_xyzz_xxyzzzz[k] = -g_0_z_yzz_xxyzzzz[k] * ab_x + g_0_z_yzz_xxxyzzzz[k];

                g_0_z_xyzz_xxzzzzz[k] = -g_0_z_yzz_xxzzzzz[k] * ab_x + g_0_z_yzz_xxxzzzzz[k];

                g_0_z_xyzz_xyyyyyy[k] = -g_0_z_yzz_xyyyyyy[k] * ab_x + g_0_z_yzz_xxyyyyyy[k];

                g_0_z_xyzz_xyyyyyz[k] = -g_0_z_yzz_xyyyyyz[k] * ab_x + g_0_z_yzz_xxyyyyyz[k];

                g_0_z_xyzz_xyyyyzz[k] = -g_0_z_yzz_xyyyyzz[k] * ab_x + g_0_z_yzz_xxyyyyzz[k];

                g_0_z_xyzz_xyyyzzz[k] = -g_0_z_yzz_xyyyzzz[k] * ab_x + g_0_z_yzz_xxyyyzzz[k];

                g_0_z_xyzz_xyyzzzz[k] = -g_0_z_yzz_xyyzzzz[k] * ab_x + g_0_z_yzz_xxyyzzzz[k];

                g_0_z_xyzz_xyzzzzz[k] = -g_0_z_yzz_xyzzzzz[k] * ab_x + g_0_z_yzz_xxyzzzzz[k];

                g_0_z_xyzz_xzzzzzz[k] = -g_0_z_yzz_xzzzzzz[k] * ab_x + g_0_z_yzz_xxzzzzzz[k];

                g_0_z_xyzz_yyyyyyy[k] = -g_0_z_yzz_yyyyyyy[k] * ab_x + g_0_z_yzz_xyyyyyyy[k];

                g_0_z_xyzz_yyyyyyz[k] = -g_0_z_yzz_yyyyyyz[k] * ab_x + g_0_z_yzz_xyyyyyyz[k];

                g_0_z_xyzz_yyyyyzz[k] = -g_0_z_yzz_yyyyyzz[k] * ab_x + g_0_z_yzz_xyyyyyzz[k];

                g_0_z_xyzz_yyyyzzz[k] = -g_0_z_yzz_yyyyzzz[k] * ab_x + g_0_z_yzz_xyyyyzzz[k];

                g_0_z_xyzz_yyyzzzz[k] = -g_0_z_yzz_yyyzzzz[k] * ab_x + g_0_z_yzz_xyyyzzzz[k];

                g_0_z_xyzz_yyzzzzz[k] = -g_0_z_yzz_yyzzzzz[k] * ab_x + g_0_z_yzz_xyyzzzzz[k];

                g_0_z_xyzz_yzzzzzz[k] = -g_0_z_yzz_yzzzzzz[k] * ab_x + g_0_z_yzz_xyzzzzzz[k];

                g_0_z_xyzz_zzzzzzz[k] = -g_0_z_yzz_zzzzzzz[k] * ab_x + g_0_z_yzz_xzzzzzzz[k];
            }

            /// Set up 1404-1440 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xzzz_xxxxxxx, g_0_z_xzzz_xxxxxxy, g_0_z_xzzz_xxxxxxz, g_0_z_xzzz_xxxxxyy, g_0_z_xzzz_xxxxxyz, g_0_z_xzzz_xxxxxzz, g_0_z_xzzz_xxxxyyy, g_0_z_xzzz_xxxxyyz, g_0_z_xzzz_xxxxyzz, g_0_z_xzzz_xxxxzzz, g_0_z_xzzz_xxxyyyy, g_0_z_xzzz_xxxyyyz, g_0_z_xzzz_xxxyyzz, g_0_z_xzzz_xxxyzzz, g_0_z_xzzz_xxxzzzz, g_0_z_xzzz_xxyyyyy, g_0_z_xzzz_xxyyyyz, g_0_z_xzzz_xxyyyzz, g_0_z_xzzz_xxyyzzz, g_0_z_xzzz_xxyzzzz, g_0_z_xzzz_xxzzzzz, g_0_z_xzzz_xyyyyyy, g_0_z_xzzz_xyyyyyz, g_0_z_xzzz_xyyyyzz, g_0_z_xzzz_xyyyzzz, g_0_z_xzzz_xyyzzzz, g_0_z_xzzz_xyzzzzz, g_0_z_xzzz_xzzzzzz, g_0_z_xzzz_yyyyyyy, g_0_z_xzzz_yyyyyyz, g_0_z_xzzz_yyyyyzz, g_0_z_xzzz_yyyyzzz, g_0_z_xzzz_yyyzzzz, g_0_z_xzzz_yyzzzzz, g_0_z_xzzz_yzzzzzz, g_0_z_xzzz_zzzzzzz, g_0_z_zzz_xxxxxxx, g_0_z_zzz_xxxxxxxx, g_0_z_zzz_xxxxxxxy, g_0_z_zzz_xxxxxxxz, g_0_z_zzz_xxxxxxy, g_0_z_zzz_xxxxxxyy, g_0_z_zzz_xxxxxxyz, g_0_z_zzz_xxxxxxz, g_0_z_zzz_xxxxxxzz, g_0_z_zzz_xxxxxyy, g_0_z_zzz_xxxxxyyy, g_0_z_zzz_xxxxxyyz, g_0_z_zzz_xxxxxyz, g_0_z_zzz_xxxxxyzz, g_0_z_zzz_xxxxxzz, g_0_z_zzz_xxxxxzzz, g_0_z_zzz_xxxxyyy, g_0_z_zzz_xxxxyyyy, g_0_z_zzz_xxxxyyyz, g_0_z_zzz_xxxxyyz, g_0_z_zzz_xxxxyyzz, g_0_z_zzz_xxxxyzz, g_0_z_zzz_xxxxyzzz, g_0_z_zzz_xxxxzzz, g_0_z_zzz_xxxxzzzz, g_0_z_zzz_xxxyyyy, g_0_z_zzz_xxxyyyyy, g_0_z_zzz_xxxyyyyz, g_0_z_zzz_xxxyyyz, g_0_z_zzz_xxxyyyzz, g_0_z_zzz_xxxyyzz, g_0_z_zzz_xxxyyzzz, g_0_z_zzz_xxxyzzz, g_0_z_zzz_xxxyzzzz, g_0_z_zzz_xxxzzzz, g_0_z_zzz_xxxzzzzz, g_0_z_zzz_xxyyyyy, g_0_z_zzz_xxyyyyyy, g_0_z_zzz_xxyyyyyz, g_0_z_zzz_xxyyyyz, g_0_z_zzz_xxyyyyzz, g_0_z_zzz_xxyyyzz, g_0_z_zzz_xxyyyzzz, g_0_z_zzz_xxyyzzz, g_0_z_zzz_xxyyzzzz, g_0_z_zzz_xxyzzzz, g_0_z_zzz_xxyzzzzz, g_0_z_zzz_xxzzzzz, g_0_z_zzz_xxzzzzzz, g_0_z_zzz_xyyyyyy, g_0_z_zzz_xyyyyyyy, g_0_z_zzz_xyyyyyyz, g_0_z_zzz_xyyyyyz, g_0_z_zzz_xyyyyyzz, g_0_z_zzz_xyyyyzz, g_0_z_zzz_xyyyyzzz, g_0_z_zzz_xyyyzzz, g_0_z_zzz_xyyyzzzz, g_0_z_zzz_xyyzzzz, g_0_z_zzz_xyyzzzzz, g_0_z_zzz_xyzzzzz, g_0_z_zzz_xyzzzzzz, g_0_z_zzz_xzzzzzz, g_0_z_zzz_xzzzzzzz, g_0_z_zzz_yyyyyyy, g_0_z_zzz_yyyyyyz, g_0_z_zzz_yyyyyzz, g_0_z_zzz_yyyyzzz, g_0_z_zzz_yyyzzzz, g_0_z_zzz_yyzzzzz, g_0_z_zzz_yzzzzzz, g_0_z_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_xxxxxxx[k] = -g_0_z_zzz_xxxxxxx[k] * ab_x + g_0_z_zzz_xxxxxxxx[k];

                g_0_z_xzzz_xxxxxxy[k] = -g_0_z_zzz_xxxxxxy[k] * ab_x + g_0_z_zzz_xxxxxxxy[k];

                g_0_z_xzzz_xxxxxxz[k] = -g_0_z_zzz_xxxxxxz[k] * ab_x + g_0_z_zzz_xxxxxxxz[k];

                g_0_z_xzzz_xxxxxyy[k] = -g_0_z_zzz_xxxxxyy[k] * ab_x + g_0_z_zzz_xxxxxxyy[k];

                g_0_z_xzzz_xxxxxyz[k] = -g_0_z_zzz_xxxxxyz[k] * ab_x + g_0_z_zzz_xxxxxxyz[k];

                g_0_z_xzzz_xxxxxzz[k] = -g_0_z_zzz_xxxxxzz[k] * ab_x + g_0_z_zzz_xxxxxxzz[k];

                g_0_z_xzzz_xxxxyyy[k] = -g_0_z_zzz_xxxxyyy[k] * ab_x + g_0_z_zzz_xxxxxyyy[k];

                g_0_z_xzzz_xxxxyyz[k] = -g_0_z_zzz_xxxxyyz[k] * ab_x + g_0_z_zzz_xxxxxyyz[k];

                g_0_z_xzzz_xxxxyzz[k] = -g_0_z_zzz_xxxxyzz[k] * ab_x + g_0_z_zzz_xxxxxyzz[k];

                g_0_z_xzzz_xxxxzzz[k] = -g_0_z_zzz_xxxxzzz[k] * ab_x + g_0_z_zzz_xxxxxzzz[k];

                g_0_z_xzzz_xxxyyyy[k] = -g_0_z_zzz_xxxyyyy[k] * ab_x + g_0_z_zzz_xxxxyyyy[k];

                g_0_z_xzzz_xxxyyyz[k] = -g_0_z_zzz_xxxyyyz[k] * ab_x + g_0_z_zzz_xxxxyyyz[k];

                g_0_z_xzzz_xxxyyzz[k] = -g_0_z_zzz_xxxyyzz[k] * ab_x + g_0_z_zzz_xxxxyyzz[k];

                g_0_z_xzzz_xxxyzzz[k] = -g_0_z_zzz_xxxyzzz[k] * ab_x + g_0_z_zzz_xxxxyzzz[k];

                g_0_z_xzzz_xxxzzzz[k] = -g_0_z_zzz_xxxzzzz[k] * ab_x + g_0_z_zzz_xxxxzzzz[k];

                g_0_z_xzzz_xxyyyyy[k] = -g_0_z_zzz_xxyyyyy[k] * ab_x + g_0_z_zzz_xxxyyyyy[k];

                g_0_z_xzzz_xxyyyyz[k] = -g_0_z_zzz_xxyyyyz[k] * ab_x + g_0_z_zzz_xxxyyyyz[k];

                g_0_z_xzzz_xxyyyzz[k] = -g_0_z_zzz_xxyyyzz[k] * ab_x + g_0_z_zzz_xxxyyyzz[k];

                g_0_z_xzzz_xxyyzzz[k] = -g_0_z_zzz_xxyyzzz[k] * ab_x + g_0_z_zzz_xxxyyzzz[k];

                g_0_z_xzzz_xxyzzzz[k] = -g_0_z_zzz_xxyzzzz[k] * ab_x + g_0_z_zzz_xxxyzzzz[k];

                g_0_z_xzzz_xxzzzzz[k] = -g_0_z_zzz_xxzzzzz[k] * ab_x + g_0_z_zzz_xxxzzzzz[k];

                g_0_z_xzzz_xyyyyyy[k] = -g_0_z_zzz_xyyyyyy[k] * ab_x + g_0_z_zzz_xxyyyyyy[k];

                g_0_z_xzzz_xyyyyyz[k] = -g_0_z_zzz_xyyyyyz[k] * ab_x + g_0_z_zzz_xxyyyyyz[k];

                g_0_z_xzzz_xyyyyzz[k] = -g_0_z_zzz_xyyyyzz[k] * ab_x + g_0_z_zzz_xxyyyyzz[k];

                g_0_z_xzzz_xyyyzzz[k] = -g_0_z_zzz_xyyyzzz[k] * ab_x + g_0_z_zzz_xxyyyzzz[k];

                g_0_z_xzzz_xyyzzzz[k] = -g_0_z_zzz_xyyzzzz[k] * ab_x + g_0_z_zzz_xxyyzzzz[k];

                g_0_z_xzzz_xyzzzzz[k] = -g_0_z_zzz_xyzzzzz[k] * ab_x + g_0_z_zzz_xxyzzzzz[k];

                g_0_z_xzzz_xzzzzzz[k] = -g_0_z_zzz_xzzzzzz[k] * ab_x + g_0_z_zzz_xxzzzzzz[k];

                g_0_z_xzzz_yyyyyyy[k] = -g_0_z_zzz_yyyyyyy[k] * ab_x + g_0_z_zzz_xyyyyyyy[k];

                g_0_z_xzzz_yyyyyyz[k] = -g_0_z_zzz_yyyyyyz[k] * ab_x + g_0_z_zzz_xyyyyyyz[k];

                g_0_z_xzzz_yyyyyzz[k] = -g_0_z_zzz_yyyyyzz[k] * ab_x + g_0_z_zzz_xyyyyyzz[k];

                g_0_z_xzzz_yyyyzzz[k] = -g_0_z_zzz_yyyyzzz[k] * ab_x + g_0_z_zzz_xyyyyzzz[k];

                g_0_z_xzzz_yyyzzzz[k] = -g_0_z_zzz_yyyzzzz[k] * ab_x + g_0_z_zzz_xyyyzzzz[k];

                g_0_z_xzzz_yyzzzzz[k] = -g_0_z_zzz_yyzzzzz[k] * ab_x + g_0_z_zzz_xyyzzzzz[k];

                g_0_z_xzzz_yzzzzzz[k] = -g_0_z_zzz_yzzzzzz[k] * ab_x + g_0_z_zzz_xyzzzzzz[k];

                g_0_z_xzzz_zzzzzzz[k] = -g_0_z_zzz_zzzzzzz[k] * ab_x + g_0_z_zzz_xzzzzzzz[k];
            }

            /// Set up 1440-1476 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyy_xxxxxxx, g_0_z_yyy_xxxxxxxy, g_0_z_yyy_xxxxxxy, g_0_z_yyy_xxxxxxyy, g_0_z_yyy_xxxxxxyz, g_0_z_yyy_xxxxxxz, g_0_z_yyy_xxxxxyy, g_0_z_yyy_xxxxxyyy, g_0_z_yyy_xxxxxyyz, g_0_z_yyy_xxxxxyz, g_0_z_yyy_xxxxxyzz, g_0_z_yyy_xxxxxzz, g_0_z_yyy_xxxxyyy, g_0_z_yyy_xxxxyyyy, g_0_z_yyy_xxxxyyyz, g_0_z_yyy_xxxxyyz, g_0_z_yyy_xxxxyyzz, g_0_z_yyy_xxxxyzz, g_0_z_yyy_xxxxyzzz, g_0_z_yyy_xxxxzzz, g_0_z_yyy_xxxyyyy, g_0_z_yyy_xxxyyyyy, g_0_z_yyy_xxxyyyyz, g_0_z_yyy_xxxyyyz, g_0_z_yyy_xxxyyyzz, g_0_z_yyy_xxxyyzz, g_0_z_yyy_xxxyyzzz, g_0_z_yyy_xxxyzzz, g_0_z_yyy_xxxyzzzz, g_0_z_yyy_xxxzzzz, g_0_z_yyy_xxyyyyy, g_0_z_yyy_xxyyyyyy, g_0_z_yyy_xxyyyyyz, g_0_z_yyy_xxyyyyz, g_0_z_yyy_xxyyyyzz, g_0_z_yyy_xxyyyzz, g_0_z_yyy_xxyyyzzz, g_0_z_yyy_xxyyzzz, g_0_z_yyy_xxyyzzzz, g_0_z_yyy_xxyzzzz, g_0_z_yyy_xxyzzzzz, g_0_z_yyy_xxzzzzz, g_0_z_yyy_xyyyyyy, g_0_z_yyy_xyyyyyyy, g_0_z_yyy_xyyyyyyz, g_0_z_yyy_xyyyyyz, g_0_z_yyy_xyyyyyzz, g_0_z_yyy_xyyyyzz, g_0_z_yyy_xyyyyzzz, g_0_z_yyy_xyyyzzz, g_0_z_yyy_xyyyzzzz, g_0_z_yyy_xyyzzzz, g_0_z_yyy_xyyzzzzz, g_0_z_yyy_xyzzzzz, g_0_z_yyy_xyzzzzzz, g_0_z_yyy_xzzzzzz, g_0_z_yyy_yyyyyyy, g_0_z_yyy_yyyyyyyy, g_0_z_yyy_yyyyyyyz, g_0_z_yyy_yyyyyyz, g_0_z_yyy_yyyyyyzz, g_0_z_yyy_yyyyyzz, g_0_z_yyy_yyyyyzzz, g_0_z_yyy_yyyyzzz, g_0_z_yyy_yyyyzzzz, g_0_z_yyy_yyyzzzz, g_0_z_yyy_yyyzzzzz, g_0_z_yyy_yyzzzzz, g_0_z_yyy_yyzzzzzz, g_0_z_yyy_yzzzzzz, g_0_z_yyy_yzzzzzzz, g_0_z_yyy_zzzzzzz, g_0_z_yyyy_xxxxxxx, g_0_z_yyyy_xxxxxxy, g_0_z_yyyy_xxxxxxz, g_0_z_yyyy_xxxxxyy, g_0_z_yyyy_xxxxxyz, g_0_z_yyyy_xxxxxzz, g_0_z_yyyy_xxxxyyy, g_0_z_yyyy_xxxxyyz, g_0_z_yyyy_xxxxyzz, g_0_z_yyyy_xxxxzzz, g_0_z_yyyy_xxxyyyy, g_0_z_yyyy_xxxyyyz, g_0_z_yyyy_xxxyyzz, g_0_z_yyyy_xxxyzzz, g_0_z_yyyy_xxxzzzz, g_0_z_yyyy_xxyyyyy, g_0_z_yyyy_xxyyyyz, g_0_z_yyyy_xxyyyzz, g_0_z_yyyy_xxyyzzz, g_0_z_yyyy_xxyzzzz, g_0_z_yyyy_xxzzzzz, g_0_z_yyyy_xyyyyyy, g_0_z_yyyy_xyyyyyz, g_0_z_yyyy_xyyyyzz, g_0_z_yyyy_xyyyzzz, g_0_z_yyyy_xyyzzzz, g_0_z_yyyy_xyzzzzz, g_0_z_yyyy_xzzzzzz, g_0_z_yyyy_yyyyyyy, g_0_z_yyyy_yyyyyyz, g_0_z_yyyy_yyyyyzz, g_0_z_yyyy_yyyyzzz, g_0_z_yyyy_yyyzzzz, g_0_z_yyyy_yyzzzzz, g_0_z_yyyy_yzzzzzz, g_0_z_yyyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_xxxxxxx[k] = -g_0_z_yyy_xxxxxxx[k] * ab_y + g_0_z_yyy_xxxxxxxy[k];

                g_0_z_yyyy_xxxxxxy[k] = -g_0_z_yyy_xxxxxxy[k] * ab_y + g_0_z_yyy_xxxxxxyy[k];

                g_0_z_yyyy_xxxxxxz[k] = -g_0_z_yyy_xxxxxxz[k] * ab_y + g_0_z_yyy_xxxxxxyz[k];

                g_0_z_yyyy_xxxxxyy[k] = -g_0_z_yyy_xxxxxyy[k] * ab_y + g_0_z_yyy_xxxxxyyy[k];

                g_0_z_yyyy_xxxxxyz[k] = -g_0_z_yyy_xxxxxyz[k] * ab_y + g_0_z_yyy_xxxxxyyz[k];

                g_0_z_yyyy_xxxxxzz[k] = -g_0_z_yyy_xxxxxzz[k] * ab_y + g_0_z_yyy_xxxxxyzz[k];

                g_0_z_yyyy_xxxxyyy[k] = -g_0_z_yyy_xxxxyyy[k] * ab_y + g_0_z_yyy_xxxxyyyy[k];

                g_0_z_yyyy_xxxxyyz[k] = -g_0_z_yyy_xxxxyyz[k] * ab_y + g_0_z_yyy_xxxxyyyz[k];

                g_0_z_yyyy_xxxxyzz[k] = -g_0_z_yyy_xxxxyzz[k] * ab_y + g_0_z_yyy_xxxxyyzz[k];

                g_0_z_yyyy_xxxxzzz[k] = -g_0_z_yyy_xxxxzzz[k] * ab_y + g_0_z_yyy_xxxxyzzz[k];

                g_0_z_yyyy_xxxyyyy[k] = -g_0_z_yyy_xxxyyyy[k] * ab_y + g_0_z_yyy_xxxyyyyy[k];

                g_0_z_yyyy_xxxyyyz[k] = -g_0_z_yyy_xxxyyyz[k] * ab_y + g_0_z_yyy_xxxyyyyz[k];

                g_0_z_yyyy_xxxyyzz[k] = -g_0_z_yyy_xxxyyzz[k] * ab_y + g_0_z_yyy_xxxyyyzz[k];

                g_0_z_yyyy_xxxyzzz[k] = -g_0_z_yyy_xxxyzzz[k] * ab_y + g_0_z_yyy_xxxyyzzz[k];

                g_0_z_yyyy_xxxzzzz[k] = -g_0_z_yyy_xxxzzzz[k] * ab_y + g_0_z_yyy_xxxyzzzz[k];

                g_0_z_yyyy_xxyyyyy[k] = -g_0_z_yyy_xxyyyyy[k] * ab_y + g_0_z_yyy_xxyyyyyy[k];

                g_0_z_yyyy_xxyyyyz[k] = -g_0_z_yyy_xxyyyyz[k] * ab_y + g_0_z_yyy_xxyyyyyz[k];

                g_0_z_yyyy_xxyyyzz[k] = -g_0_z_yyy_xxyyyzz[k] * ab_y + g_0_z_yyy_xxyyyyzz[k];

                g_0_z_yyyy_xxyyzzz[k] = -g_0_z_yyy_xxyyzzz[k] * ab_y + g_0_z_yyy_xxyyyzzz[k];

                g_0_z_yyyy_xxyzzzz[k] = -g_0_z_yyy_xxyzzzz[k] * ab_y + g_0_z_yyy_xxyyzzzz[k];

                g_0_z_yyyy_xxzzzzz[k] = -g_0_z_yyy_xxzzzzz[k] * ab_y + g_0_z_yyy_xxyzzzzz[k];

                g_0_z_yyyy_xyyyyyy[k] = -g_0_z_yyy_xyyyyyy[k] * ab_y + g_0_z_yyy_xyyyyyyy[k];

                g_0_z_yyyy_xyyyyyz[k] = -g_0_z_yyy_xyyyyyz[k] * ab_y + g_0_z_yyy_xyyyyyyz[k];

                g_0_z_yyyy_xyyyyzz[k] = -g_0_z_yyy_xyyyyzz[k] * ab_y + g_0_z_yyy_xyyyyyzz[k];

                g_0_z_yyyy_xyyyzzz[k] = -g_0_z_yyy_xyyyzzz[k] * ab_y + g_0_z_yyy_xyyyyzzz[k];

                g_0_z_yyyy_xyyzzzz[k] = -g_0_z_yyy_xyyzzzz[k] * ab_y + g_0_z_yyy_xyyyzzzz[k];

                g_0_z_yyyy_xyzzzzz[k] = -g_0_z_yyy_xyzzzzz[k] * ab_y + g_0_z_yyy_xyyzzzzz[k];

                g_0_z_yyyy_xzzzzzz[k] = -g_0_z_yyy_xzzzzzz[k] * ab_y + g_0_z_yyy_xyzzzzzz[k];

                g_0_z_yyyy_yyyyyyy[k] = -g_0_z_yyy_yyyyyyy[k] * ab_y + g_0_z_yyy_yyyyyyyy[k];

                g_0_z_yyyy_yyyyyyz[k] = -g_0_z_yyy_yyyyyyz[k] * ab_y + g_0_z_yyy_yyyyyyyz[k];

                g_0_z_yyyy_yyyyyzz[k] = -g_0_z_yyy_yyyyyzz[k] * ab_y + g_0_z_yyy_yyyyyyzz[k];

                g_0_z_yyyy_yyyyzzz[k] = -g_0_z_yyy_yyyyzzz[k] * ab_y + g_0_z_yyy_yyyyyzzz[k];

                g_0_z_yyyy_yyyzzzz[k] = -g_0_z_yyy_yyyzzzz[k] * ab_y + g_0_z_yyy_yyyyzzzz[k];

                g_0_z_yyyy_yyzzzzz[k] = -g_0_z_yyy_yyzzzzz[k] * ab_y + g_0_z_yyy_yyyzzzzz[k];

                g_0_z_yyyy_yzzzzzz[k] = -g_0_z_yyy_yzzzzzz[k] * ab_y + g_0_z_yyy_yyzzzzzz[k];

                g_0_z_yyyy_zzzzzzz[k] = -g_0_z_yyy_zzzzzzz[k] * ab_y + g_0_z_yyy_yzzzzzzz[k];
            }

            /// Set up 1476-1512 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyyz_xxxxxxx, g_0_z_yyyz_xxxxxxy, g_0_z_yyyz_xxxxxxz, g_0_z_yyyz_xxxxxyy, g_0_z_yyyz_xxxxxyz, g_0_z_yyyz_xxxxxzz, g_0_z_yyyz_xxxxyyy, g_0_z_yyyz_xxxxyyz, g_0_z_yyyz_xxxxyzz, g_0_z_yyyz_xxxxzzz, g_0_z_yyyz_xxxyyyy, g_0_z_yyyz_xxxyyyz, g_0_z_yyyz_xxxyyzz, g_0_z_yyyz_xxxyzzz, g_0_z_yyyz_xxxzzzz, g_0_z_yyyz_xxyyyyy, g_0_z_yyyz_xxyyyyz, g_0_z_yyyz_xxyyyzz, g_0_z_yyyz_xxyyzzz, g_0_z_yyyz_xxyzzzz, g_0_z_yyyz_xxzzzzz, g_0_z_yyyz_xyyyyyy, g_0_z_yyyz_xyyyyyz, g_0_z_yyyz_xyyyyzz, g_0_z_yyyz_xyyyzzz, g_0_z_yyyz_xyyzzzz, g_0_z_yyyz_xyzzzzz, g_0_z_yyyz_xzzzzzz, g_0_z_yyyz_yyyyyyy, g_0_z_yyyz_yyyyyyz, g_0_z_yyyz_yyyyyzz, g_0_z_yyyz_yyyyzzz, g_0_z_yyyz_yyyzzzz, g_0_z_yyyz_yyzzzzz, g_0_z_yyyz_yzzzzzz, g_0_z_yyyz_zzzzzzz, g_0_z_yyz_xxxxxxx, g_0_z_yyz_xxxxxxxy, g_0_z_yyz_xxxxxxy, g_0_z_yyz_xxxxxxyy, g_0_z_yyz_xxxxxxyz, g_0_z_yyz_xxxxxxz, g_0_z_yyz_xxxxxyy, g_0_z_yyz_xxxxxyyy, g_0_z_yyz_xxxxxyyz, g_0_z_yyz_xxxxxyz, g_0_z_yyz_xxxxxyzz, g_0_z_yyz_xxxxxzz, g_0_z_yyz_xxxxyyy, g_0_z_yyz_xxxxyyyy, g_0_z_yyz_xxxxyyyz, g_0_z_yyz_xxxxyyz, g_0_z_yyz_xxxxyyzz, g_0_z_yyz_xxxxyzz, g_0_z_yyz_xxxxyzzz, g_0_z_yyz_xxxxzzz, g_0_z_yyz_xxxyyyy, g_0_z_yyz_xxxyyyyy, g_0_z_yyz_xxxyyyyz, g_0_z_yyz_xxxyyyz, g_0_z_yyz_xxxyyyzz, g_0_z_yyz_xxxyyzz, g_0_z_yyz_xxxyyzzz, g_0_z_yyz_xxxyzzz, g_0_z_yyz_xxxyzzzz, g_0_z_yyz_xxxzzzz, g_0_z_yyz_xxyyyyy, g_0_z_yyz_xxyyyyyy, g_0_z_yyz_xxyyyyyz, g_0_z_yyz_xxyyyyz, g_0_z_yyz_xxyyyyzz, g_0_z_yyz_xxyyyzz, g_0_z_yyz_xxyyyzzz, g_0_z_yyz_xxyyzzz, g_0_z_yyz_xxyyzzzz, g_0_z_yyz_xxyzzzz, g_0_z_yyz_xxyzzzzz, g_0_z_yyz_xxzzzzz, g_0_z_yyz_xyyyyyy, g_0_z_yyz_xyyyyyyy, g_0_z_yyz_xyyyyyyz, g_0_z_yyz_xyyyyyz, g_0_z_yyz_xyyyyyzz, g_0_z_yyz_xyyyyzz, g_0_z_yyz_xyyyyzzz, g_0_z_yyz_xyyyzzz, g_0_z_yyz_xyyyzzzz, g_0_z_yyz_xyyzzzz, g_0_z_yyz_xyyzzzzz, g_0_z_yyz_xyzzzzz, g_0_z_yyz_xyzzzzzz, g_0_z_yyz_xzzzzzz, g_0_z_yyz_yyyyyyy, g_0_z_yyz_yyyyyyyy, g_0_z_yyz_yyyyyyyz, g_0_z_yyz_yyyyyyz, g_0_z_yyz_yyyyyyzz, g_0_z_yyz_yyyyyzz, g_0_z_yyz_yyyyyzzz, g_0_z_yyz_yyyyzzz, g_0_z_yyz_yyyyzzzz, g_0_z_yyz_yyyzzzz, g_0_z_yyz_yyyzzzzz, g_0_z_yyz_yyzzzzz, g_0_z_yyz_yyzzzzzz, g_0_z_yyz_yzzzzzz, g_0_z_yyz_yzzzzzzz, g_0_z_yyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_xxxxxxx[k] = -g_0_z_yyz_xxxxxxx[k] * ab_y + g_0_z_yyz_xxxxxxxy[k];

                g_0_z_yyyz_xxxxxxy[k] = -g_0_z_yyz_xxxxxxy[k] * ab_y + g_0_z_yyz_xxxxxxyy[k];

                g_0_z_yyyz_xxxxxxz[k] = -g_0_z_yyz_xxxxxxz[k] * ab_y + g_0_z_yyz_xxxxxxyz[k];

                g_0_z_yyyz_xxxxxyy[k] = -g_0_z_yyz_xxxxxyy[k] * ab_y + g_0_z_yyz_xxxxxyyy[k];

                g_0_z_yyyz_xxxxxyz[k] = -g_0_z_yyz_xxxxxyz[k] * ab_y + g_0_z_yyz_xxxxxyyz[k];

                g_0_z_yyyz_xxxxxzz[k] = -g_0_z_yyz_xxxxxzz[k] * ab_y + g_0_z_yyz_xxxxxyzz[k];

                g_0_z_yyyz_xxxxyyy[k] = -g_0_z_yyz_xxxxyyy[k] * ab_y + g_0_z_yyz_xxxxyyyy[k];

                g_0_z_yyyz_xxxxyyz[k] = -g_0_z_yyz_xxxxyyz[k] * ab_y + g_0_z_yyz_xxxxyyyz[k];

                g_0_z_yyyz_xxxxyzz[k] = -g_0_z_yyz_xxxxyzz[k] * ab_y + g_0_z_yyz_xxxxyyzz[k];

                g_0_z_yyyz_xxxxzzz[k] = -g_0_z_yyz_xxxxzzz[k] * ab_y + g_0_z_yyz_xxxxyzzz[k];

                g_0_z_yyyz_xxxyyyy[k] = -g_0_z_yyz_xxxyyyy[k] * ab_y + g_0_z_yyz_xxxyyyyy[k];

                g_0_z_yyyz_xxxyyyz[k] = -g_0_z_yyz_xxxyyyz[k] * ab_y + g_0_z_yyz_xxxyyyyz[k];

                g_0_z_yyyz_xxxyyzz[k] = -g_0_z_yyz_xxxyyzz[k] * ab_y + g_0_z_yyz_xxxyyyzz[k];

                g_0_z_yyyz_xxxyzzz[k] = -g_0_z_yyz_xxxyzzz[k] * ab_y + g_0_z_yyz_xxxyyzzz[k];

                g_0_z_yyyz_xxxzzzz[k] = -g_0_z_yyz_xxxzzzz[k] * ab_y + g_0_z_yyz_xxxyzzzz[k];

                g_0_z_yyyz_xxyyyyy[k] = -g_0_z_yyz_xxyyyyy[k] * ab_y + g_0_z_yyz_xxyyyyyy[k];

                g_0_z_yyyz_xxyyyyz[k] = -g_0_z_yyz_xxyyyyz[k] * ab_y + g_0_z_yyz_xxyyyyyz[k];

                g_0_z_yyyz_xxyyyzz[k] = -g_0_z_yyz_xxyyyzz[k] * ab_y + g_0_z_yyz_xxyyyyzz[k];

                g_0_z_yyyz_xxyyzzz[k] = -g_0_z_yyz_xxyyzzz[k] * ab_y + g_0_z_yyz_xxyyyzzz[k];

                g_0_z_yyyz_xxyzzzz[k] = -g_0_z_yyz_xxyzzzz[k] * ab_y + g_0_z_yyz_xxyyzzzz[k];

                g_0_z_yyyz_xxzzzzz[k] = -g_0_z_yyz_xxzzzzz[k] * ab_y + g_0_z_yyz_xxyzzzzz[k];

                g_0_z_yyyz_xyyyyyy[k] = -g_0_z_yyz_xyyyyyy[k] * ab_y + g_0_z_yyz_xyyyyyyy[k];

                g_0_z_yyyz_xyyyyyz[k] = -g_0_z_yyz_xyyyyyz[k] * ab_y + g_0_z_yyz_xyyyyyyz[k];

                g_0_z_yyyz_xyyyyzz[k] = -g_0_z_yyz_xyyyyzz[k] * ab_y + g_0_z_yyz_xyyyyyzz[k];

                g_0_z_yyyz_xyyyzzz[k] = -g_0_z_yyz_xyyyzzz[k] * ab_y + g_0_z_yyz_xyyyyzzz[k];

                g_0_z_yyyz_xyyzzzz[k] = -g_0_z_yyz_xyyzzzz[k] * ab_y + g_0_z_yyz_xyyyzzzz[k];

                g_0_z_yyyz_xyzzzzz[k] = -g_0_z_yyz_xyzzzzz[k] * ab_y + g_0_z_yyz_xyyzzzzz[k];

                g_0_z_yyyz_xzzzzzz[k] = -g_0_z_yyz_xzzzzzz[k] * ab_y + g_0_z_yyz_xyzzzzzz[k];

                g_0_z_yyyz_yyyyyyy[k] = -g_0_z_yyz_yyyyyyy[k] * ab_y + g_0_z_yyz_yyyyyyyy[k];

                g_0_z_yyyz_yyyyyyz[k] = -g_0_z_yyz_yyyyyyz[k] * ab_y + g_0_z_yyz_yyyyyyyz[k];

                g_0_z_yyyz_yyyyyzz[k] = -g_0_z_yyz_yyyyyzz[k] * ab_y + g_0_z_yyz_yyyyyyzz[k];

                g_0_z_yyyz_yyyyzzz[k] = -g_0_z_yyz_yyyyzzz[k] * ab_y + g_0_z_yyz_yyyyyzzz[k];

                g_0_z_yyyz_yyyzzzz[k] = -g_0_z_yyz_yyyzzzz[k] * ab_y + g_0_z_yyz_yyyyzzzz[k];

                g_0_z_yyyz_yyzzzzz[k] = -g_0_z_yyz_yyzzzzz[k] * ab_y + g_0_z_yyz_yyyzzzzz[k];

                g_0_z_yyyz_yzzzzzz[k] = -g_0_z_yyz_yzzzzzz[k] * ab_y + g_0_z_yyz_yyzzzzzz[k];

                g_0_z_yyyz_zzzzzzz[k] = -g_0_z_yyz_zzzzzzz[k] * ab_y + g_0_z_yyz_yzzzzzzz[k];
            }

            /// Set up 1512-1548 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyzz_xxxxxxx, g_0_z_yyzz_xxxxxxy, g_0_z_yyzz_xxxxxxz, g_0_z_yyzz_xxxxxyy, g_0_z_yyzz_xxxxxyz, g_0_z_yyzz_xxxxxzz, g_0_z_yyzz_xxxxyyy, g_0_z_yyzz_xxxxyyz, g_0_z_yyzz_xxxxyzz, g_0_z_yyzz_xxxxzzz, g_0_z_yyzz_xxxyyyy, g_0_z_yyzz_xxxyyyz, g_0_z_yyzz_xxxyyzz, g_0_z_yyzz_xxxyzzz, g_0_z_yyzz_xxxzzzz, g_0_z_yyzz_xxyyyyy, g_0_z_yyzz_xxyyyyz, g_0_z_yyzz_xxyyyzz, g_0_z_yyzz_xxyyzzz, g_0_z_yyzz_xxyzzzz, g_0_z_yyzz_xxzzzzz, g_0_z_yyzz_xyyyyyy, g_0_z_yyzz_xyyyyyz, g_0_z_yyzz_xyyyyzz, g_0_z_yyzz_xyyyzzz, g_0_z_yyzz_xyyzzzz, g_0_z_yyzz_xyzzzzz, g_0_z_yyzz_xzzzzzz, g_0_z_yyzz_yyyyyyy, g_0_z_yyzz_yyyyyyz, g_0_z_yyzz_yyyyyzz, g_0_z_yyzz_yyyyzzz, g_0_z_yyzz_yyyzzzz, g_0_z_yyzz_yyzzzzz, g_0_z_yyzz_yzzzzzz, g_0_z_yyzz_zzzzzzz, g_0_z_yzz_xxxxxxx, g_0_z_yzz_xxxxxxxy, g_0_z_yzz_xxxxxxy, g_0_z_yzz_xxxxxxyy, g_0_z_yzz_xxxxxxyz, g_0_z_yzz_xxxxxxz, g_0_z_yzz_xxxxxyy, g_0_z_yzz_xxxxxyyy, g_0_z_yzz_xxxxxyyz, g_0_z_yzz_xxxxxyz, g_0_z_yzz_xxxxxyzz, g_0_z_yzz_xxxxxzz, g_0_z_yzz_xxxxyyy, g_0_z_yzz_xxxxyyyy, g_0_z_yzz_xxxxyyyz, g_0_z_yzz_xxxxyyz, g_0_z_yzz_xxxxyyzz, g_0_z_yzz_xxxxyzz, g_0_z_yzz_xxxxyzzz, g_0_z_yzz_xxxxzzz, g_0_z_yzz_xxxyyyy, g_0_z_yzz_xxxyyyyy, g_0_z_yzz_xxxyyyyz, g_0_z_yzz_xxxyyyz, g_0_z_yzz_xxxyyyzz, g_0_z_yzz_xxxyyzz, g_0_z_yzz_xxxyyzzz, g_0_z_yzz_xxxyzzz, g_0_z_yzz_xxxyzzzz, g_0_z_yzz_xxxzzzz, g_0_z_yzz_xxyyyyy, g_0_z_yzz_xxyyyyyy, g_0_z_yzz_xxyyyyyz, g_0_z_yzz_xxyyyyz, g_0_z_yzz_xxyyyyzz, g_0_z_yzz_xxyyyzz, g_0_z_yzz_xxyyyzzz, g_0_z_yzz_xxyyzzz, g_0_z_yzz_xxyyzzzz, g_0_z_yzz_xxyzzzz, g_0_z_yzz_xxyzzzzz, g_0_z_yzz_xxzzzzz, g_0_z_yzz_xyyyyyy, g_0_z_yzz_xyyyyyyy, g_0_z_yzz_xyyyyyyz, g_0_z_yzz_xyyyyyz, g_0_z_yzz_xyyyyyzz, g_0_z_yzz_xyyyyzz, g_0_z_yzz_xyyyyzzz, g_0_z_yzz_xyyyzzz, g_0_z_yzz_xyyyzzzz, g_0_z_yzz_xyyzzzz, g_0_z_yzz_xyyzzzzz, g_0_z_yzz_xyzzzzz, g_0_z_yzz_xyzzzzzz, g_0_z_yzz_xzzzzzz, g_0_z_yzz_yyyyyyy, g_0_z_yzz_yyyyyyyy, g_0_z_yzz_yyyyyyyz, g_0_z_yzz_yyyyyyz, g_0_z_yzz_yyyyyyzz, g_0_z_yzz_yyyyyzz, g_0_z_yzz_yyyyyzzz, g_0_z_yzz_yyyyzzz, g_0_z_yzz_yyyyzzzz, g_0_z_yzz_yyyzzzz, g_0_z_yzz_yyyzzzzz, g_0_z_yzz_yyzzzzz, g_0_z_yzz_yyzzzzzz, g_0_z_yzz_yzzzzzz, g_0_z_yzz_yzzzzzzz, g_0_z_yzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_xxxxxxx[k] = -g_0_z_yzz_xxxxxxx[k] * ab_y + g_0_z_yzz_xxxxxxxy[k];

                g_0_z_yyzz_xxxxxxy[k] = -g_0_z_yzz_xxxxxxy[k] * ab_y + g_0_z_yzz_xxxxxxyy[k];

                g_0_z_yyzz_xxxxxxz[k] = -g_0_z_yzz_xxxxxxz[k] * ab_y + g_0_z_yzz_xxxxxxyz[k];

                g_0_z_yyzz_xxxxxyy[k] = -g_0_z_yzz_xxxxxyy[k] * ab_y + g_0_z_yzz_xxxxxyyy[k];

                g_0_z_yyzz_xxxxxyz[k] = -g_0_z_yzz_xxxxxyz[k] * ab_y + g_0_z_yzz_xxxxxyyz[k];

                g_0_z_yyzz_xxxxxzz[k] = -g_0_z_yzz_xxxxxzz[k] * ab_y + g_0_z_yzz_xxxxxyzz[k];

                g_0_z_yyzz_xxxxyyy[k] = -g_0_z_yzz_xxxxyyy[k] * ab_y + g_0_z_yzz_xxxxyyyy[k];

                g_0_z_yyzz_xxxxyyz[k] = -g_0_z_yzz_xxxxyyz[k] * ab_y + g_0_z_yzz_xxxxyyyz[k];

                g_0_z_yyzz_xxxxyzz[k] = -g_0_z_yzz_xxxxyzz[k] * ab_y + g_0_z_yzz_xxxxyyzz[k];

                g_0_z_yyzz_xxxxzzz[k] = -g_0_z_yzz_xxxxzzz[k] * ab_y + g_0_z_yzz_xxxxyzzz[k];

                g_0_z_yyzz_xxxyyyy[k] = -g_0_z_yzz_xxxyyyy[k] * ab_y + g_0_z_yzz_xxxyyyyy[k];

                g_0_z_yyzz_xxxyyyz[k] = -g_0_z_yzz_xxxyyyz[k] * ab_y + g_0_z_yzz_xxxyyyyz[k];

                g_0_z_yyzz_xxxyyzz[k] = -g_0_z_yzz_xxxyyzz[k] * ab_y + g_0_z_yzz_xxxyyyzz[k];

                g_0_z_yyzz_xxxyzzz[k] = -g_0_z_yzz_xxxyzzz[k] * ab_y + g_0_z_yzz_xxxyyzzz[k];

                g_0_z_yyzz_xxxzzzz[k] = -g_0_z_yzz_xxxzzzz[k] * ab_y + g_0_z_yzz_xxxyzzzz[k];

                g_0_z_yyzz_xxyyyyy[k] = -g_0_z_yzz_xxyyyyy[k] * ab_y + g_0_z_yzz_xxyyyyyy[k];

                g_0_z_yyzz_xxyyyyz[k] = -g_0_z_yzz_xxyyyyz[k] * ab_y + g_0_z_yzz_xxyyyyyz[k];

                g_0_z_yyzz_xxyyyzz[k] = -g_0_z_yzz_xxyyyzz[k] * ab_y + g_0_z_yzz_xxyyyyzz[k];

                g_0_z_yyzz_xxyyzzz[k] = -g_0_z_yzz_xxyyzzz[k] * ab_y + g_0_z_yzz_xxyyyzzz[k];

                g_0_z_yyzz_xxyzzzz[k] = -g_0_z_yzz_xxyzzzz[k] * ab_y + g_0_z_yzz_xxyyzzzz[k];

                g_0_z_yyzz_xxzzzzz[k] = -g_0_z_yzz_xxzzzzz[k] * ab_y + g_0_z_yzz_xxyzzzzz[k];

                g_0_z_yyzz_xyyyyyy[k] = -g_0_z_yzz_xyyyyyy[k] * ab_y + g_0_z_yzz_xyyyyyyy[k];

                g_0_z_yyzz_xyyyyyz[k] = -g_0_z_yzz_xyyyyyz[k] * ab_y + g_0_z_yzz_xyyyyyyz[k];

                g_0_z_yyzz_xyyyyzz[k] = -g_0_z_yzz_xyyyyzz[k] * ab_y + g_0_z_yzz_xyyyyyzz[k];

                g_0_z_yyzz_xyyyzzz[k] = -g_0_z_yzz_xyyyzzz[k] * ab_y + g_0_z_yzz_xyyyyzzz[k];

                g_0_z_yyzz_xyyzzzz[k] = -g_0_z_yzz_xyyzzzz[k] * ab_y + g_0_z_yzz_xyyyzzzz[k];

                g_0_z_yyzz_xyzzzzz[k] = -g_0_z_yzz_xyzzzzz[k] * ab_y + g_0_z_yzz_xyyzzzzz[k];

                g_0_z_yyzz_xzzzzzz[k] = -g_0_z_yzz_xzzzzzz[k] * ab_y + g_0_z_yzz_xyzzzzzz[k];

                g_0_z_yyzz_yyyyyyy[k] = -g_0_z_yzz_yyyyyyy[k] * ab_y + g_0_z_yzz_yyyyyyyy[k];

                g_0_z_yyzz_yyyyyyz[k] = -g_0_z_yzz_yyyyyyz[k] * ab_y + g_0_z_yzz_yyyyyyyz[k];

                g_0_z_yyzz_yyyyyzz[k] = -g_0_z_yzz_yyyyyzz[k] * ab_y + g_0_z_yzz_yyyyyyzz[k];

                g_0_z_yyzz_yyyyzzz[k] = -g_0_z_yzz_yyyyzzz[k] * ab_y + g_0_z_yzz_yyyyyzzz[k];

                g_0_z_yyzz_yyyzzzz[k] = -g_0_z_yzz_yyyzzzz[k] * ab_y + g_0_z_yzz_yyyyzzzz[k];

                g_0_z_yyzz_yyzzzzz[k] = -g_0_z_yzz_yyzzzzz[k] * ab_y + g_0_z_yzz_yyyzzzzz[k];

                g_0_z_yyzz_yzzzzzz[k] = -g_0_z_yzz_yzzzzzz[k] * ab_y + g_0_z_yzz_yyzzzzzz[k];

                g_0_z_yyzz_zzzzzzz[k] = -g_0_z_yzz_zzzzzzz[k] * ab_y + g_0_z_yzz_yzzzzzzz[k];
            }

            /// Set up 1548-1584 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yzzz_xxxxxxx, g_0_z_yzzz_xxxxxxy, g_0_z_yzzz_xxxxxxz, g_0_z_yzzz_xxxxxyy, g_0_z_yzzz_xxxxxyz, g_0_z_yzzz_xxxxxzz, g_0_z_yzzz_xxxxyyy, g_0_z_yzzz_xxxxyyz, g_0_z_yzzz_xxxxyzz, g_0_z_yzzz_xxxxzzz, g_0_z_yzzz_xxxyyyy, g_0_z_yzzz_xxxyyyz, g_0_z_yzzz_xxxyyzz, g_0_z_yzzz_xxxyzzz, g_0_z_yzzz_xxxzzzz, g_0_z_yzzz_xxyyyyy, g_0_z_yzzz_xxyyyyz, g_0_z_yzzz_xxyyyzz, g_0_z_yzzz_xxyyzzz, g_0_z_yzzz_xxyzzzz, g_0_z_yzzz_xxzzzzz, g_0_z_yzzz_xyyyyyy, g_0_z_yzzz_xyyyyyz, g_0_z_yzzz_xyyyyzz, g_0_z_yzzz_xyyyzzz, g_0_z_yzzz_xyyzzzz, g_0_z_yzzz_xyzzzzz, g_0_z_yzzz_xzzzzzz, g_0_z_yzzz_yyyyyyy, g_0_z_yzzz_yyyyyyz, g_0_z_yzzz_yyyyyzz, g_0_z_yzzz_yyyyzzz, g_0_z_yzzz_yyyzzzz, g_0_z_yzzz_yyzzzzz, g_0_z_yzzz_yzzzzzz, g_0_z_yzzz_zzzzzzz, g_0_z_zzz_xxxxxxx, g_0_z_zzz_xxxxxxxy, g_0_z_zzz_xxxxxxy, g_0_z_zzz_xxxxxxyy, g_0_z_zzz_xxxxxxyz, g_0_z_zzz_xxxxxxz, g_0_z_zzz_xxxxxyy, g_0_z_zzz_xxxxxyyy, g_0_z_zzz_xxxxxyyz, g_0_z_zzz_xxxxxyz, g_0_z_zzz_xxxxxyzz, g_0_z_zzz_xxxxxzz, g_0_z_zzz_xxxxyyy, g_0_z_zzz_xxxxyyyy, g_0_z_zzz_xxxxyyyz, g_0_z_zzz_xxxxyyz, g_0_z_zzz_xxxxyyzz, g_0_z_zzz_xxxxyzz, g_0_z_zzz_xxxxyzzz, g_0_z_zzz_xxxxzzz, g_0_z_zzz_xxxyyyy, g_0_z_zzz_xxxyyyyy, g_0_z_zzz_xxxyyyyz, g_0_z_zzz_xxxyyyz, g_0_z_zzz_xxxyyyzz, g_0_z_zzz_xxxyyzz, g_0_z_zzz_xxxyyzzz, g_0_z_zzz_xxxyzzz, g_0_z_zzz_xxxyzzzz, g_0_z_zzz_xxxzzzz, g_0_z_zzz_xxyyyyy, g_0_z_zzz_xxyyyyyy, g_0_z_zzz_xxyyyyyz, g_0_z_zzz_xxyyyyz, g_0_z_zzz_xxyyyyzz, g_0_z_zzz_xxyyyzz, g_0_z_zzz_xxyyyzzz, g_0_z_zzz_xxyyzzz, g_0_z_zzz_xxyyzzzz, g_0_z_zzz_xxyzzzz, g_0_z_zzz_xxyzzzzz, g_0_z_zzz_xxzzzzz, g_0_z_zzz_xyyyyyy, g_0_z_zzz_xyyyyyyy, g_0_z_zzz_xyyyyyyz, g_0_z_zzz_xyyyyyz, g_0_z_zzz_xyyyyyzz, g_0_z_zzz_xyyyyzz, g_0_z_zzz_xyyyyzzz, g_0_z_zzz_xyyyzzz, g_0_z_zzz_xyyyzzzz, g_0_z_zzz_xyyzzzz, g_0_z_zzz_xyyzzzzz, g_0_z_zzz_xyzzzzz, g_0_z_zzz_xyzzzzzz, g_0_z_zzz_xzzzzzz, g_0_z_zzz_yyyyyyy, g_0_z_zzz_yyyyyyyy, g_0_z_zzz_yyyyyyyz, g_0_z_zzz_yyyyyyz, g_0_z_zzz_yyyyyyzz, g_0_z_zzz_yyyyyzz, g_0_z_zzz_yyyyyzzz, g_0_z_zzz_yyyyzzz, g_0_z_zzz_yyyyzzzz, g_0_z_zzz_yyyzzzz, g_0_z_zzz_yyyzzzzz, g_0_z_zzz_yyzzzzz, g_0_z_zzz_yyzzzzzz, g_0_z_zzz_yzzzzzz, g_0_z_zzz_yzzzzzzz, g_0_z_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_xxxxxxx[k] = -g_0_z_zzz_xxxxxxx[k] * ab_y + g_0_z_zzz_xxxxxxxy[k];

                g_0_z_yzzz_xxxxxxy[k] = -g_0_z_zzz_xxxxxxy[k] * ab_y + g_0_z_zzz_xxxxxxyy[k];

                g_0_z_yzzz_xxxxxxz[k] = -g_0_z_zzz_xxxxxxz[k] * ab_y + g_0_z_zzz_xxxxxxyz[k];

                g_0_z_yzzz_xxxxxyy[k] = -g_0_z_zzz_xxxxxyy[k] * ab_y + g_0_z_zzz_xxxxxyyy[k];

                g_0_z_yzzz_xxxxxyz[k] = -g_0_z_zzz_xxxxxyz[k] * ab_y + g_0_z_zzz_xxxxxyyz[k];

                g_0_z_yzzz_xxxxxzz[k] = -g_0_z_zzz_xxxxxzz[k] * ab_y + g_0_z_zzz_xxxxxyzz[k];

                g_0_z_yzzz_xxxxyyy[k] = -g_0_z_zzz_xxxxyyy[k] * ab_y + g_0_z_zzz_xxxxyyyy[k];

                g_0_z_yzzz_xxxxyyz[k] = -g_0_z_zzz_xxxxyyz[k] * ab_y + g_0_z_zzz_xxxxyyyz[k];

                g_0_z_yzzz_xxxxyzz[k] = -g_0_z_zzz_xxxxyzz[k] * ab_y + g_0_z_zzz_xxxxyyzz[k];

                g_0_z_yzzz_xxxxzzz[k] = -g_0_z_zzz_xxxxzzz[k] * ab_y + g_0_z_zzz_xxxxyzzz[k];

                g_0_z_yzzz_xxxyyyy[k] = -g_0_z_zzz_xxxyyyy[k] * ab_y + g_0_z_zzz_xxxyyyyy[k];

                g_0_z_yzzz_xxxyyyz[k] = -g_0_z_zzz_xxxyyyz[k] * ab_y + g_0_z_zzz_xxxyyyyz[k];

                g_0_z_yzzz_xxxyyzz[k] = -g_0_z_zzz_xxxyyzz[k] * ab_y + g_0_z_zzz_xxxyyyzz[k];

                g_0_z_yzzz_xxxyzzz[k] = -g_0_z_zzz_xxxyzzz[k] * ab_y + g_0_z_zzz_xxxyyzzz[k];

                g_0_z_yzzz_xxxzzzz[k] = -g_0_z_zzz_xxxzzzz[k] * ab_y + g_0_z_zzz_xxxyzzzz[k];

                g_0_z_yzzz_xxyyyyy[k] = -g_0_z_zzz_xxyyyyy[k] * ab_y + g_0_z_zzz_xxyyyyyy[k];

                g_0_z_yzzz_xxyyyyz[k] = -g_0_z_zzz_xxyyyyz[k] * ab_y + g_0_z_zzz_xxyyyyyz[k];

                g_0_z_yzzz_xxyyyzz[k] = -g_0_z_zzz_xxyyyzz[k] * ab_y + g_0_z_zzz_xxyyyyzz[k];

                g_0_z_yzzz_xxyyzzz[k] = -g_0_z_zzz_xxyyzzz[k] * ab_y + g_0_z_zzz_xxyyyzzz[k];

                g_0_z_yzzz_xxyzzzz[k] = -g_0_z_zzz_xxyzzzz[k] * ab_y + g_0_z_zzz_xxyyzzzz[k];

                g_0_z_yzzz_xxzzzzz[k] = -g_0_z_zzz_xxzzzzz[k] * ab_y + g_0_z_zzz_xxyzzzzz[k];

                g_0_z_yzzz_xyyyyyy[k] = -g_0_z_zzz_xyyyyyy[k] * ab_y + g_0_z_zzz_xyyyyyyy[k];

                g_0_z_yzzz_xyyyyyz[k] = -g_0_z_zzz_xyyyyyz[k] * ab_y + g_0_z_zzz_xyyyyyyz[k];

                g_0_z_yzzz_xyyyyzz[k] = -g_0_z_zzz_xyyyyzz[k] * ab_y + g_0_z_zzz_xyyyyyzz[k];

                g_0_z_yzzz_xyyyzzz[k] = -g_0_z_zzz_xyyyzzz[k] * ab_y + g_0_z_zzz_xyyyyzzz[k];

                g_0_z_yzzz_xyyzzzz[k] = -g_0_z_zzz_xyyzzzz[k] * ab_y + g_0_z_zzz_xyyyzzzz[k];

                g_0_z_yzzz_xyzzzzz[k] = -g_0_z_zzz_xyzzzzz[k] * ab_y + g_0_z_zzz_xyyzzzzz[k];

                g_0_z_yzzz_xzzzzzz[k] = -g_0_z_zzz_xzzzzzz[k] * ab_y + g_0_z_zzz_xyzzzzzz[k];

                g_0_z_yzzz_yyyyyyy[k] = -g_0_z_zzz_yyyyyyy[k] * ab_y + g_0_z_zzz_yyyyyyyy[k];

                g_0_z_yzzz_yyyyyyz[k] = -g_0_z_zzz_yyyyyyz[k] * ab_y + g_0_z_zzz_yyyyyyyz[k];

                g_0_z_yzzz_yyyyyzz[k] = -g_0_z_zzz_yyyyyzz[k] * ab_y + g_0_z_zzz_yyyyyyzz[k];

                g_0_z_yzzz_yyyyzzz[k] = -g_0_z_zzz_yyyyzzz[k] * ab_y + g_0_z_zzz_yyyyyzzz[k];

                g_0_z_yzzz_yyyzzzz[k] = -g_0_z_zzz_yyyzzzz[k] * ab_y + g_0_z_zzz_yyyyzzzz[k];

                g_0_z_yzzz_yyzzzzz[k] = -g_0_z_zzz_yyzzzzz[k] * ab_y + g_0_z_zzz_yyyzzzzz[k];

                g_0_z_yzzz_yzzzzzz[k] = -g_0_z_zzz_yzzzzzz[k] * ab_y + g_0_z_zzz_yyzzzzzz[k];

                g_0_z_yzzz_zzzzzzz[k] = -g_0_z_zzz_zzzzzzz[k] * ab_y + g_0_z_zzz_yzzzzzzz[k];
            }

            /// Set up 1584-1620 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_zzz_xxxxxxx, g_0_z_zzz_xxxxxxxz, g_0_z_zzz_xxxxxxy, g_0_z_zzz_xxxxxxyz, g_0_z_zzz_xxxxxxz, g_0_z_zzz_xxxxxxzz, g_0_z_zzz_xxxxxyy, g_0_z_zzz_xxxxxyyz, g_0_z_zzz_xxxxxyz, g_0_z_zzz_xxxxxyzz, g_0_z_zzz_xxxxxzz, g_0_z_zzz_xxxxxzzz, g_0_z_zzz_xxxxyyy, g_0_z_zzz_xxxxyyyz, g_0_z_zzz_xxxxyyz, g_0_z_zzz_xxxxyyzz, g_0_z_zzz_xxxxyzz, g_0_z_zzz_xxxxyzzz, g_0_z_zzz_xxxxzzz, g_0_z_zzz_xxxxzzzz, g_0_z_zzz_xxxyyyy, g_0_z_zzz_xxxyyyyz, g_0_z_zzz_xxxyyyz, g_0_z_zzz_xxxyyyzz, g_0_z_zzz_xxxyyzz, g_0_z_zzz_xxxyyzzz, g_0_z_zzz_xxxyzzz, g_0_z_zzz_xxxyzzzz, g_0_z_zzz_xxxzzzz, g_0_z_zzz_xxxzzzzz, g_0_z_zzz_xxyyyyy, g_0_z_zzz_xxyyyyyz, g_0_z_zzz_xxyyyyz, g_0_z_zzz_xxyyyyzz, g_0_z_zzz_xxyyyzz, g_0_z_zzz_xxyyyzzz, g_0_z_zzz_xxyyzzz, g_0_z_zzz_xxyyzzzz, g_0_z_zzz_xxyzzzz, g_0_z_zzz_xxyzzzzz, g_0_z_zzz_xxzzzzz, g_0_z_zzz_xxzzzzzz, g_0_z_zzz_xyyyyyy, g_0_z_zzz_xyyyyyyz, g_0_z_zzz_xyyyyyz, g_0_z_zzz_xyyyyyzz, g_0_z_zzz_xyyyyzz, g_0_z_zzz_xyyyyzzz, g_0_z_zzz_xyyyzzz, g_0_z_zzz_xyyyzzzz, g_0_z_zzz_xyyzzzz, g_0_z_zzz_xyyzzzzz, g_0_z_zzz_xyzzzzz, g_0_z_zzz_xyzzzzzz, g_0_z_zzz_xzzzzzz, g_0_z_zzz_xzzzzzzz, g_0_z_zzz_yyyyyyy, g_0_z_zzz_yyyyyyyz, g_0_z_zzz_yyyyyyz, g_0_z_zzz_yyyyyyzz, g_0_z_zzz_yyyyyzz, g_0_z_zzz_yyyyyzzz, g_0_z_zzz_yyyyzzz, g_0_z_zzz_yyyyzzzz, g_0_z_zzz_yyyzzzz, g_0_z_zzz_yyyzzzzz, g_0_z_zzz_yyzzzzz, g_0_z_zzz_yyzzzzzz, g_0_z_zzz_yzzzzzz, g_0_z_zzz_yzzzzzzz, g_0_z_zzz_zzzzzzz, g_0_z_zzz_zzzzzzzz, g_0_z_zzzz_xxxxxxx, g_0_z_zzzz_xxxxxxy, g_0_z_zzzz_xxxxxxz, g_0_z_zzzz_xxxxxyy, g_0_z_zzzz_xxxxxyz, g_0_z_zzzz_xxxxxzz, g_0_z_zzzz_xxxxyyy, g_0_z_zzzz_xxxxyyz, g_0_z_zzzz_xxxxyzz, g_0_z_zzzz_xxxxzzz, g_0_z_zzzz_xxxyyyy, g_0_z_zzzz_xxxyyyz, g_0_z_zzzz_xxxyyzz, g_0_z_zzzz_xxxyzzz, g_0_z_zzzz_xxxzzzz, g_0_z_zzzz_xxyyyyy, g_0_z_zzzz_xxyyyyz, g_0_z_zzzz_xxyyyzz, g_0_z_zzzz_xxyyzzz, g_0_z_zzzz_xxyzzzz, g_0_z_zzzz_xxzzzzz, g_0_z_zzzz_xyyyyyy, g_0_z_zzzz_xyyyyyz, g_0_z_zzzz_xyyyyzz, g_0_z_zzzz_xyyyzzz, g_0_z_zzzz_xyyzzzz, g_0_z_zzzz_xyzzzzz, g_0_z_zzzz_xzzzzzz, g_0_z_zzzz_yyyyyyy, g_0_z_zzzz_yyyyyyz, g_0_z_zzzz_yyyyyzz, g_0_z_zzzz_yyyyzzz, g_0_z_zzzz_yyyzzzz, g_0_z_zzzz_yyzzzzz, g_0_z_zzzz_yzzzzzz, g_0_z_zzzz_zzzzzzz, g_zzz_xxxxxxx, g_zzz_xxxxxxy, g_zzz_xxxxxxz, g_zzz_xxxxxyy, g_zzz_xxxxxyz, g_zzz_xxxxxzz, g_zzz_xxxxyyy, g_zzz_xxxxyyz, g_zzz_xxxxyzz, g_zzz_xxxxzzz, g_zzz_xxxyyyy, g_zzz_xxxyyyz, g_zzz_xxxyyzz, g_zzz_xxxyzzz, g_zzz_xxxzzzz, g_zzz_xxyyyyy, g_zzz_xxyyyyz, g_zzz_xxyyyzz, g_zzz_xxyyzzz, g_zzz_xxyzzzz, g_zzz_xxzzzzz, g_zzz_xyyyyyy, g_zzz_xyyyyyz, g_zzz_xyyyyzz, g_zzz_xyyyzzz, g_zzz_xyyzzzz, g_zzz_xyzzzzz, g_zzz_xzzzzzz, g_zzz_yyyyyyy, g_zzz_yyyyyyz, g_zzz_yyyyyzz, g_zzz_yyyyzzz, g_zzz_yyyzzzz, g_zzz_yyzzzzz, g_zzz_yzzzzzz, g_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_xxxxxxx[k] = g_zzz_xxxxxxx[k] - g_0_z_zzz_xxxxxxx[k] * ab_z + g_0_z_zzz_xxxxxxxz[k];

                g_0_z_zzzz_xxxxxxy[k] = g_zzz_xxxxxxy[k] - g_0_z_zzz_xxxxxxy[k] * ab_z + g_0_z_zzz_xxxxxxyz[k];

                g_0_z_zzzz_xxxxxxz[k] = g_zzz_xxxxxxz[k] - g_0_z_zzz_xxxxxxz[k] * ab_z + g_0_z_zzz_xxxxxxzz[k];

                g_0_z_zzzz_xxxxxyy[k] = g_zzz_xxxxxyy[k] - g_0_z_zzz_xxxxxyy[k] * ab_z + g_0_z_zzz_xxxxxyyz[k];

                g_0_z_zzzz_xxxxxyz[k] = g_zzz_xxxxxyz[k] - g_0_z_zzz_xxxxxyz[k] * ab_z + g_0_z_zzz_xxxxxyzz[k];

                g_0_z_zzzz_xxxxxzz[k] = g_zzz_xxxxxzz[k] - g_0_z_zzz_xxxxxzz[k] * ab_z + g_0_z_zzz_xxxxxzzz[k];

                g_0_z_zzzz_xxxxyyy[k] = g_zzz_xxxxyyy[k] - g_0_z_zzz_xxxxyyy[k] * ab_z + g_0_z_zzz_xxxxyyyz[k];

                g_0_z_zzzz_xxxxyyz[k] = g_zzz_xxxxyyz[k] - g_0_z_zzz_xxxxyyz[k] * ab_z + g_0_z_zzz_xxxxyyzz[k];

                g_0_z_zzzz_xxxxyzz[k] = g_zzz_xxxxyzz[k] - g_0_z_zzz_xxxxyzz[k] * ab_z + g_0_z_zzz_xxxxyzzz[k];

                g_0_z_zzzz_xxxxzzz[k] = g_zzz_xxxxzzz[k] - g_0_z_zzz_xxxxzzz[k] * ab_z + g_0_z_zzz_xxxxzzzz[k];

                g_0_z_zzzz_xxxyyyy[k] = g_zzz_xxxyyyy[k] - g_0_z_zzz_xxxyyyy[k] * ab_z + g_0_z_zzz_xxxyyyyz[k];

                g_0_z_zzzz_xxxyyyz[k] = g_zzz_xxxyyyz[k] - g_0_z_zzz_xxxyyyz[k] * ab_z + g_0_z_zzz_xxxyyyzz[k];

                g_0_z_zzzz_xxxyyzz[k] = g_zzz_xxxyyzz[k] - g_0_z_zzz_xxxyyzz[k] * ab_z + g_0_z_zzz_xxxyyzzz[k];

                g_0_z_zzzz_xxxyzzz[k] = g_zzz_xxxyzzz[k] - g_0_z_zzz_xxxyzzz[k] * ab_z + g_0_z_zzz_xxxyzzzz[k];

                g_0_z_zzzz_xxxzzzz[k] = g_zzz_xxxzzzz[k] - g_0_z_zzz_xxxzzzz[k] * ab_z + g_0_z_zzz_xxxzzzzz[k];

                g_0_z_zzzz_xxyyyyy[k] = g_zzz_xxyyyyy[k] - g_0_z_zzz_xxyyyyy[k] * ab_z + g_0_z_zzz_xxyyyyyz[k];

                g_0_z_zzzz_xxyyyyz[k] = g_zzz_xxyyyyz[k] - g_0_z_zzz_xxyyyyz[k] * ab_z + g_0_z_zzz_xxyyyyzz[k];

                g_0_z_zzzz_xxyyyzz[k] = g_zzz_xxyyyzz[k] - g_0_z_zzz_xxyyyzz[k] * ab_z + g_0_z_zzz_xxyyyzzz[k];

                g_0_z_zzzz_xxyyzzz[k] = g_zzz_xxyyzzz[k] - g_0_z_zzz_xxyyzzz[k] * ab_z + g_0_z_zzz_xxyyzzzz[k];

                g_0_z_zzzz_xxyzzzz[k] = g_zzz_xxyzzzz[k] - g_0_z_zzz_xxyzzzz[k] * ab_z + g_0_z_zzz_xxyzzzzz[k];

                g_0_z_zzzz_xxzzzzz[k] = g_zzz_xxzzzzz[k] - g_0_z_zzz_xxzzzzz[k] * ab_z + g_0_z_zzz_xxzzzzzz[k];

                g_0_z_zzzz_xyyyyyy[k] = g_zzz_xyyyyyy[k] - g_0_z_zzz_xyyyyyy[k] * ab_z + g_0_z_zzz_xyyyyyyz[k];

                g_0_z_zzzz_xyyyyyz[k] = g_zzz_xyyyyyz[k] - g_0_z_zzz_xyyyyyz[k] * ab_z + g_0_z_zzz_xyyyyyzz[k];

                g_0_z_zzzz_xyyyyzz[k] = g_zzz_xyyyyzz[k] - g_0_z_zzz_xyyyyzz[k] * ab_z + g_0_z_zzz_xyyyyzzz[k];

                g_0_z_zzzz_xyyyzzz[k] = g_zzz_xyyyzzz[k] - g_0_z_zzz_xyyyzzz[k] * ab_z + g_0_z_zzz_xyyyzzzz[k];

                g_0_z_zzzz_xyyzzzz[k] = g_zzz_xyyzzzz[k] - g_0_z_zzz_xyyzzzz[k] * ab_z + g_0_z_zzz_xyyzzzzz[k];

                g_0_z_zzzz_xyzzzzz[k] = g_zzz_xyzzzzz[k] - g_0_z_zzz_xyzzzzz[k] * ab_z + g_0_z_zzz_xyzzzzzz[k];

                g_0_z_zzzz_xzzzzzz[k] = g_zzz_xzzzzzz[k] - g_0_z_zzz_xzzzzzz[k] * ab_z + g_0_z_zzz_xzzzzzzz[k];

                g_0_z_zzzz_yyyyyyy[k] = g_zzz_yyyyyyy[k] - g_0_z_zzz_yyyyyyy[k] * ab_z + g_0_z_zzz_yyyyyyyz[k];

                g_0_z_zzzz_yyyyyyz[k] = g_zzz_yyyyyyz[k] - g_0_z_zzz_yyyyyyz[k] * ab_z + g_0_z_zzz_yyyyyyzz[k];

                g_0_z_zzzz_yyyyyzz[k] = g_zzz_yyyyyzz[k] - g_0_z_zzz_yyyyyzz[k] * ab_z + g_0_z_zzz_yyyyyzzz[k];

                g_0_z_zzzz_yyyyzzz[k] = g_zzz_yyyyzzz[k] - g_0_z_zzz_yyyyzzz[k] * ab_z + g_0_z_zzz_yyyyzzzz[k];

                g_0_z_zzzz_yyyzzzz[k] = g_zzz_yyyzzzz[k] - g_0_z_zzz_yyyzzzz[k] * ab_z + g_0_z_zzz_yyyzzzzz[k];

                g_0_z_zzzz_yyzzzzz[k] = g_zzz_yyzzzzz[k] - g_0_z_zzz_yyzzzzz[k] * ab_z + g_0_z_zzz_yyzzzzzz[k];

                g_0_z_zzzz_yzzzzzz[k] = g_zzz_yzzzzzz[k] - g_0_z_zzz_yzzzzzz[k] * ab_z + g_0_z_zzz_yzzzzzzz[k];

                g_0_z_zzzz_zzzzzzz[k] = g_zzz_zzzzzzz[k] - g_0_z_zzz_zzzzzzz[k] * ab_z + g_0_z_zzz_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

