#include "ElectronRepulsionGeom0100ContrRecFKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_fkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_fkxx,
                                            const size_t idx_dkxx,
                                            const size_t idx_geom_01_dkxx,
                                            const size_t idx_geom_01_dlxx,
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
            /// Set up components of auxilary buffer : DKSS

            const auto dk_off = idx_dkxx + i * dcomps + j;

            auto g_xx_xxxxxxx = cbuffer.data(dk_off + 0 * ccomps * dcomps);

            auto g_xx_xxxxxxy = cbuffer.data(dk_off + 1 * ccomps * dcomps);

            auto g_xx_xxxxxxz = cbuffer.data(dk_off + 2 * ccomps * dcomps);

            auto g_xx_xxxxxyy = cbuffer.data(dk_off + 3 * ccomps * dcomps);

            auto g_xx_xxxxxyz = cbuffer.data(dk_off + 4 * ccomps * dcomps);

            auto g_xx_xxxxxzz = cbuffer.data(dk_off + 5 * ccomps * dcomps);

            auto g_xx_xxxxyyy = cbuffer.data(dk_off + 6 * ccomps * dcomps);

            auto g_xx_xxxxyyz = cbuffer.data(dk_off + 7 * ccomps * dcomps);

            auto g_xx_xxxxyzz = cbuffer.data(dk_off + 8 * ccomps * dcomps);

            auto g_xx_xxxxzzz = cbuffer.data(dk_off + 9 * ccomps * dcomps);

            auto g_xx_xxxyyyy = cbuffer.data(dk_off + 10 * ccomps * dcomps);

            auto g_xx_xxxyyyz = cbuffer.data(dk_off + 11 * ccomps * dcomps);

            auto g_xx_xxxyyzz = cbuffer.data(dk_off + 12 * ccomps * dcomps);

            auto g_xx_xxxyzzz = cbuffer.data(dk_off + 13 * ccomps * dcomps);

            auto g_xx_xxxzzzz = cbuffer.data(dk_off + 14 * ccomps * dcomps);

            auto g_xx_xxyyyyy = cbuffer.data(dk_off + 15 * ccomps * dcomps);

            auto g_xx_xxyyyyz = cbuffer.data(dk_off + 16 * ccomps * dcomps);

            auto g_xx_xxyyyzz = cbuffer.data(dk_off + 17 * ccomps * dcomps);

            auto g_xx_xxyyzzz = cbuffer.data(dk_off + 18 * ccomps * dcomps);

            auto g_xx_xxyzzzz = cbuffer.data(dk_off + 19 * ccomps * dcomps);

            auto g_xx_xxzzzzz = cbuffer.data(dk_off + 20 * ccomps * dcomps);

            auto g_xx_xyyyyyy = cbuffer.data(dk_off + 21 * ccomps * dcomps);

            auto g_xx_xyyyyyz = cbuffer.data(dk_off + 22 * ccomps * dcomps);

            auto g_xx_xyyyyzz = cbuffer.data(dk_off + 23 * ccomps * dcomps);

            auto g_xx_xyyyzzz = cbuffer.data(dk_off + 24 * ccomps * dcomps);

            auto g_xx_xyyzzzz = cbuffer.data(dk_off + 25 * ccomps * dcomps);

            auto g_xx_xyzzzzz = cbuffer.data(dk_off + 26 * ccomps * dcomps);

            auto g_xx_xzzzzzz = cbuffer.data(dk_off + 27 * ccomps * dcomps);

            auto g_xx_yyyyyyy = cbuffer.data(dk_off + 28 * ccomps * dcomps);

            auto g_xx_yyyyyyz = cbuffer.data(dk_off + 29 * ccomps * dcomps);

            auto g_xx_yyyyyzz = cbuffer.data(dk_off + 30 * ccomps * dcomps);

            auto g_xx_yyyyzzz = cbuffer.data(dk_off + 31 * ccomps * dcomps);

            auto g_xx_yyyzzzz = cbuffer.data(dk_off + 32 * ccomps * dcomps);

            auto g_xx_yyzzzzz = cbuffer.data(dk_off + 33 * ccomps * dcomps);

            auto g_xx_yzzzzzz = cbuffer.data(dk_off + 34 * ccomps * dcomps);

            auto g_xx_zzzzzzz = cbuffer.data(dk_off + 35 * ccomps * dcomps);

            auto g_xy_xxxxxxx = cbuffer.data(dk_off + 36 * ccomps * dcomps);

            auto g_xy_xxxxxxy = cbuffer.data(dk_off + 37 * ccomps * dcomps);

            auto g_xy_xxxxxxz = cbuffer.data(dk_off + 38 * ccomps * dcomps);

            auto g_xy_xxxxxyy = cbuffer.data(dk_off + 39 * ccomps * dcomps);

            auto g_xy_xxxxxyz = cbuffer.data(dk_off + 40 * ccomps * dcomps);

            auto g_xy_xxxxxzz = cbuffer.data(dk_off + 41 * ccomps * dcomps);

            auto g_xy_xxxxyyy = cbuffer.data(dk_off + 42 * ccomps * dcomps);

            auto g_xy_xxxxyyz = cbuffer.data(dk_off + 43 * ccomps * dcomps);

            auto g_xy_xxxxyzz = cbuffer.data(dk_off + 44 * ccomps * dcomps);

            auto g_xy_xxxxzzz = cbuffer.data(dk_off + 45 * ccomps * dcomps);

            auto g_xy_xxxyyyy = cbuffer.data(dk_off + 46 * ccomps * dcomps);

            auto g_xy_xxxyyyz = cbuffer.data(dk_off + 47 * ccomps * dcomps);

            auto g_xy_xxxyyzz = cbuffer.data(dk_off + 48 * ccomps * dcomps);

            auto g_xy_xxxyzzz = cbuffer.data(dk_off + 49 * ccomps * dcomps);

            auto g_xy_xxxzzzz = cbuffer.data(dk_off + 50 * ccomps * dcomps);

            auto g_xy_xxyyyyy = cbuffer.data(dk_off + 51 * ccomps * dcomps);

            auto g_xy_xxyyyyz = cbuffer.data(dk_off + 52 * ccomps * dcomps);

            auto g_xy_xxyyyzz = cbuffer.data(dk_off + 53 * ccomps * dcomps);

            auto g_xy_xxyyzzz = cbuffer.data(dk_off + 54 * ccomps * dcomps);

            auto g_xy_xxyzzzz = cbuffer.data(dk_off + 55 * ccomps * dcomps);

            auto g_xy_xxzzzzz = cbuffer.data(dk_off + 56 * ccomps * dcomps);

            auto g_xy_xyyyyyy = cbuffer.data(dk_off + 57 * ccomps * dcomps);

            auto g_xy_xyyyyyz = cbuffer.data(dk_off + 58 * ccomps * dcomps);

            auto g_xy_xyyyyzz = cbuffer.data(dk_off + 59 * ccomps * dcomps);

            auto g_xy_xyyyzzz = cbuffer.data(dk_off + 60 * ccomps * dcomps);

            auto g_xy_xyyzzzz = cbuffer.data(dk_off + 61 * ccomps * dcomps);

            auto g_xy_xyzzzzz = cbuffer.data(dk_off + 62 * ccomps * dcomps);

            auto g_xy_xzzzzzz = cbuffer.data(dk_off + 63 * ccomps * dcomps);

            auto g_xy_yyyyyyy = cbuffer.data(dk_off + 64 * ccomps * dcomps);

            auto g_xy_yyyyyyz = cbuffer.data(dk_off + 65 * ccomps * dcomps);

            auto g_xy_yyyyyzz = cbuffer.data(dk_off + 66 * ccomps * dcomps);

            auto g_xy_yyyyzzz = cbuffer.data(dk_off + 67 * ccomps * dcomps);

            auto g_xy_yyyzzzz = cbuffer.data(dk_off + 68 * ccomps * dcomps);

            auto g_xy_yyzzzzz = cbuffer.data(dk_off + 69 * ccomps * dcomps);

            auto g_xy_yzzzzzz = cbuffer.data(dk_off + 70 * ccomps * dcomps);

            auto g_xy_zzzzzzz = cbuffer.data(dk_off + 71 * ccomps * dcomps);

            auto g_xz_xxxxxxx = cbuffer.data(dk_off + 72 * ccomps * dcomps);

            auto g_xz_xxxxxxy = cbuffer.data(dk_off + 73 * ccomps * dcomps);

            auto g_xz_xxxxxxz = cbuffer.data(dk_off + 74 * ccomps * dcomps);

            auto g_xz_xxxxxyy = cbuffer.data(dk_off + 75 * ccomps * dcomps);

            auto g_xz_xxxxxyz = cbuffer.data(dk_off + 76 * ccomps * dcomps);

            auto g_xz_xxxxxzz = cbuffer.data(dk_off + 77 * ccomps * dcomps);

            auto g_xz_xxxxyyy = cbuffer.data(dk_off + 78 * ccomps * dcomps);

            auto g_xz_xxxxyyz = cbuffer.data(dk_off + 79 * ccomps * dcomps);

            auto g_xz_xxxxyzz = cbuffer.data(dk_off + 80 * ccomps * dcomps);

            auto g_xz_xxxxzzz = cbuffer.data(dk_off + 81 * ccomps * dcomps);

            auto g_xz_xxxyyyy = cbuffer.data(dk_off + 82 * ccomps * dcomps);

            auto g_xz_xxxyyyz = cbuffer.data(dk_off + 83 * ccomps * dcomps);

            auto g_xz_xxxyyzz = cbuffer.data(dk_off + 84 * ccomps * dcomps);

            auto g_xz_xxxyzzz = cbuffer.data(dk_off + 85 * ccomps * dcomps);

            auto g_xz_xxxzzzz = cbuffer.data(dk_off + 86 * ccomps * dcomps);

            auto g_xz_xxyyyyy = cbuffer.data(dk_off + 87 * ccomps * dcomps);

            auto g_xz_xxyyyyz = cbuffer.data(dk_off + 88 * ccomps * dcomps);

            auto g_xz_xxyyyzz = cbuffer.data(dk_off + 89 * ccomps * dcomps);

            auto g_xz_xxyyzzz = cbuffer.data(dk_off + 90 * ccomps * dcomps);

            auto g_xz_xxyzzzz = cbuffer.data(dk_off + 91 * ccomps * dcomps);

            auto g_xz_xxzzzzz = cbuffer.data(dk_off + 92 * ccomps * dcomps);

            auto g_xz_xyyyyyy = cbuffer.data(dk_off + 93 * ccomps * dcomps);

            auto g_xz_xyyyyyz = cbuffer.data(dk_off + 94 * ccomps * dcomps);

            auto g_xz_xyyyyzz = cbuffer.data(dk_off + 95 * ccomps * dcomps);

            auto g_xz_xyyyzzz = cbuffer.data(dk_off + 96 * ccomps * dcomps);

            auto g_xz_xyyzzzz = cbuffer.data(dk_off + 97 * ccomps * dcomps);

            auto g_xz_xyzzzzz = cbuffer.data(dk_off + 98 * ccomps * dcomps);

            auto g_xz_xzzzzzz = cbuffer.data(dk_off + 99 * ccomps * dcomps);

            auto g_xz_yyyyyyy = cbuffer.data(dk_off + 100 * ccomps * dcomps);

            auto g_xz_yyyyyyz = cbuffer.data(dk_off + 101 * ccomps * dcomps);

            auto g_xz_yyyyyzz = cbuffer.data(dk_off + 102 * ccomps * dcomps);

            auto g_xz_yyyyzzz = cbuffer.data(dk_off + 103 * ccomps * dcomps);

            auto g_xz_yyyzzzz = cbuffer.data(dk_off + 104 * ccomps * dcomps);

            auto g_xz_yyzzzzz = cbuffer.data(dk_off + 105 * ccomps * dcomps);

            auto g_xz_yzzzzzz = cbuffer.data(dk_off + 106 * ccomps * dcomps);

            auto g_xz_zzzzzzz = cbuffer.data(dk_off + 107 * ccomps * dcomps);

            auto g_yy_xxxxxxx = cbuffer.data(dk_off + 108 * ccomps * dcomps);

            auto g_yy_xxxxxxy = cbuffer.data(dk_off + 109 * ccomps * dcomps);

            auto g_yy_xxxxxxz = cbuffer.data(dk_off + 110 * ccomps * dcomps);

            auto g_yy_xxxxxyy = cbuffer.data(dk_off + 111 * ccomps * dcomps);

            auto g_yy_xxxxxyz = cbuffer.data(dk_off + 112 * ccomps * dcomps);

            auto g_yy_xxxxxzz = cbuffer.data(dk_off + 113 * ccomps * dcomps);

            auto g_yy_xxxxyyy = cbuffer.data(dk_off + 114 * ccomps * dcomps);

            auto g_yy_xxxxyyz = cbuffer.data(dk_off + 115 * ccomps * dcomps);

            auto g_yy_xxxxyzz = cbuffer.data(dk_off + 116 * ccomps * dcomps);

            auto g_yy_xxxxzzz = cbuffer.data(dk_off + 117 * ccomps * dcomps);

            auto g_yy_xxxyyyy = cbuffer.data(dk_off + 118 * ccomps * dcomps);

            auto g_yy_xxxyyyz = cbuffer.data(dk_off + 119 * ccomps * dcomps);

            auto g_yy_xxxyyzz = cbuffer.data(dk_off + 120 * ccomps * dcomps);

            auto g_yy_xxxyzzz = cbuffer.data(dk_off + 121 * ccomps * dcomps);

            auto g_yy_xxxzzzz = cbuffer.data(dk_off + 122 * ccomps * dcomps);

            auto g_yy_xxyyyyy = cbuffer.data(dk_off + 123 * ccomps * dcomps);

            auto g_yy_xxyyyyz = cbuffer.data(dk_off + 124 * ccomps * dcomps);

            auto g_yy_xxyyyzz = cbuffer.data(dk_off + 125 * ccomps * dcomps);

            auto g_yy_xxyyzzz = cbuffer.data(dk_off + 126 * ccomps * dcomps);

            auto g_yy_xxyzzzz = cbuffer.data(dk_off + 127 * ccomps * dcomps);

            auto g_yy_xxzzzzz = cbuffer.data(dk_off + 128 * ccomps * dcomps);

            auto g_yy_xyyyyyy = cbuffer.data(dk_off + 129 * ccomps * dcomps);

            auto g_yy_xyyyyyz = cbuffer.data(dk_off + 130 * ccomps * dcomps);

            auto g_yy_xyyyyzz = cbuffer.data(dk_off + 131 * ccomps * dcomps);

            auto g_yy_xyyyzzz = cbuffer.data(dk_off + 132 * ccomps * dcomps);

            auto g_yy_xyyzzzz = cbuffer.data(dk_off + 133 * ccomps * dcomps);

            auto g_yy_xyzzzzz = cbuffer.data(dk_off + 134 * ccomps * dcomps);

            auto g_yy_xzzzzzz = cbuffer.data(dk_off + 135 * ccomps * dcomps);

            auto g_yy_yyyyyyy = cbuffer.data(dk_off + 136 * ccomps * dcomps);

            auto g_yy_yyyyyyz = cbuffer.data(dk_off + 137 * ccomps * dcomps);

            auto g_yy_yyyyyzz = cbuffer.data(dk_off + 138 * ccomps * dcomps);

            auto g_yy_yyyyzzz = cbuffer.data(dk_off + 139 * ccomps * dcomps);

            auto g_yy_yyyzzzz = cbuffer.data(dk_off + 140 * ccomps * dcomps);

            auto g_yy_yyzzzzz = cbuffer.data(dk_off + 141 * ccomps * dcomps);

            auto g_yy_yzzzzzz = cbuffer.data(dk_off + 142 * ccomps * dcomps);

            auto g_yy_zzzzzzz = cbuffer.data(dk_off + 143 * ccomps * dcomps);

            auto g_yz_xxxxxxx = cbuffer.data(dk_off + 144 * ccomps * dcomps);

            auto g_yz_xxxxxxy = cbuffer.data(dk_off + 145 * ccomps * dcomps);

            auto g_yz_xxxxxxz = cbuffer.data(dk_off + 146 * ccomps * dcomps);

            auto g_yz_xxxxxyy = cbuffer.data(dk_off + 147 * ccomps * dcomps);

            auto g_yz_xxxxxyz = cbuffer.data(dk_off + 148 * ccomps * dcomps);

            auto g_yz_xxxxxzz = cbuffer.data(dk_off + 149 * ccomps * dcomps);

            auto g_yz_xxxxyyy = cbuffer.data(dk_off + 150 * ccomps * dcomps);

            auto g_yz_xxxxyyz = cbuffer.data(dk_off + 151 * ccomps * dcomps);

            auto g_yz_xxxxyzz = cbuffer.data(dk_off + 152 * ccomps * dcomps);

            auto g_yz_xxxxzzz = cbuffer.data(dk_off + 153 * ccomps * dcomps);

            auto g_yz_xxxyyyy = cbuffer.data(dk_off + 154 * ccomps * dcomps);

            auto g_yz_xxxyyyz = cbuffer.data(dk_off + 155 * ccomps * dcomps);

            auto g_yz_xxxyyzz = cbuffer.data(dk_off + 156 * ccomps * dcomps);

            auto g_yz_xxxyzzz = cbuffer.data(dk_off + 157 * ccomps * dcomps);

            auto g_yz_xxxzzzz = cbuffer.data(dk_off + 158 * ccomps * dcomps);

            auto g_yz_xxyyyyy = cbuffer.data(dk_off + 159 * ccomps * dcomps);

            auto g_yz_xxyyyyz = cbuffer.data(dk_off + 160 * ccomps * dcomps);

            auto g_yz_xxyyyzz = cbuffer.data(dk_off + 161 * ccomps * dcomps);

            auto g_yz_xxyyzzz = cbuffer.data(dk_off + 162 * ccomps * dcomps);

            auto g_yz_xxyzzzz = cbuffer.data(dk_off + 163 * ccomps * dcomps);

            auto g_yz_xxzzzzz = cbuffer.data(dk_off + 164 * ccomps * dcomps);

            auto g_yz_xyyyyyy = cbuffer.data(dk_off + 165 * ccomps * dcomps);

            auto g_yz_xyyyyyz = cbuffer.data(dk_off + 166 * ccomps * dcomps);

            auto g_yz_xyyyyzz = cbuffer.data(dk_off + 167 * ccomps * dcomps);

            auto g_yz_xyyyzzz = cbuffer.data(dk_off + 168 * ccomps * dcomps);

            auto g_yz_xyyzzzz = cbuffer.data(dk_off + 169 * ccomps * dcomps);

            auto g_yz_xyzzzzz = cbuffer.data(dk_off + 170 * ccomps * dcomps);

            auto g_yz_xzzzzzz = cbuffer.data(dk_off + 171 * ccomps * dcomps);

            auto g_yz_yyyyyyy = cbuffer.data(dk_off + 172 * ccomps * dcomps);

            auto g_yz_yyyyyyz = cbuffer.data(dk_off + 173 * ccomps * dcomps);

            auto g_yz_yyyyyzz = cbuffer.data(dk_off + 174 * ccomps * dcomps);

            auto g_yz_yyyyzzz = cbuffer.data(dk_off + 175 * ccomps * dcomps);

            auto g_yz_yyyzzzz = cbuffer.data(dk_off + 176 * ccomps * dcomps);

            auto g_yz_yyzzzzz = cbuffer.data(dk_off + 177 * ccomps * dcomps);

            auto g_yz_yzzzzzz = cbuffer.data(dk_off + 178 * ccomps * dcomps);

            auto g_yz_zzzzzzz = cbuffer.data(dk_off + 179 * ccomps * dcomps);

            auto g_zz_xxxxxxx = cbuffer.data(dk_off + 180 * ccomps * dcomps);

            auto g_zz_xxxxxxy = cbuffer.data(dk_off + 181 * ccomps * dcomps);

            auto g_zz_xxxxxxz = cbuffer.data(dk_off + 182 * ccomps * dcomps);

            auto g_zz_xxxxxyy = cbuffer.data(dk_off + 183 * ccomps * dcomps);

            auto g_zz_xxxxxyz = cbuffer.data(dk_off + 184 * ccomps * dcomps);

            auto g_zz_xxxxxzz = cbuffer.data(dk_off + 185 * ccomps * dcomps);

            auto g_zz_xxxxyyy = cbuffer.data(dk_off + 186 * ccomps * dcomps);

            auto g_zz_xxxxyyz = cbuffer.data(dk_off + 187 * ccomps * dcomps);

            auto g_zz_xxxxyzz = cbuffer.data(dk_off + 188 * ccomps * dcomps);

            auto g_zz_xxxxzzz = cbuffer.data(dk_off + 189 * ccomps * dcomps);

            auto g_zz_xxxyyyy = cbuffer.data(dk_off + 190 * ccomps * dcomps);

            auto g_zz_xxxyyyz = cbuffer.data(dk_off + 191 * ccomps * dcomps);

            auto g_zz_xxxyyzz = cbuffer.data(dk_off + 192 * ccomps * dcomps);

            auto g_zz_xxxyzzz = cbuffer.data(dk_off + 193 * ccomps * dcomps);

            auto g_zz_xxxzzzz = cbuffer.data(dk_off + 194 * ccomps * dcomps);

            auto g_zz_xxyyyyy = cbuffer.data(dk_off + 195 * ccomps * dcomps);

            auto g_zz_xxyyyyz = cbuffer.data(dk_off + 196 * ccomps * dcomps);

            auto g_zz_xxyyyzz = cbuffer.data(dk_off + 197 * ccomps * dcomps);

            auto g_zz_xxyyzzz = cbuffer.data(dk_off + 198 * ccomps * dcomps);

            auto g_zz_xxyzzzz = cbuffer.data(dk_off + 199 * ccomps * dcomps);

            auto g_zz_xxzzzzz = cbuffer.data(dk_off + 200 * ccomps * dcomps);

            auto g_zz_xyyyyyy = cbuffer.data(dk_off + 201 * ccomps * dcomps);

            auto g_zz_xyyyyyz = cbuffer.data(dk_off + 202 * ccomps * dcomps);

            auto g_zz_xyyyyzz = cbuffer.data(dk_off + 203 * ccomps * dcomps);

            auto g_zz_xyyyzzz = cbuffer.data(dk_off + 204 * ccomps * dcomps);

            auto g_zz_xyyzzzz = cbuffer.data(dk_off + 205 * ccomps * dcomps);

            auto g_zz_xyzzzzz = cbuffer.data(dk_off + 206 * ccomps * dcomps);

            auto g_zz_xzzzzzz = cbuffer.data(dk_off + 207 * ccomps * dcomps);

            auto g_zz_yyyyyyy = cbuffer.data(dk_off + 208 * ccomps * dcomps);

            auto g_zz_yyyyyyz = cbuffer.data(dk_off + 209 * ccomps * dcomps);

            auto g_zz_yyyyyzz = cbuffer.data(dk_off + 210 * ccomps * dcomps);

            auto g_zz_yyyyzzz = cbuffer.data(dk_off + 211 * ccomps * dcomps);

            auto g_zz_yyyzzzz = cbuffer.data(dk_off + 212 * ccomps * dcomps);

            auto g_zz_yyzzzzz = cbuffer.data(dk_off + 213 * ccomps * dcomps);

            auto g_zz_yzzzzzz = cbuffer.data(dk_off + 214 * ccomps * dcomps);

            auto g_zz_zzzzzzz = cbuffer.data(dk_off + 215 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DKSS

            const auto dk_geom_01_off = idx_geom_01_dkxx + i * dcomps + j;

            auto g_0_x_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 647 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DLSS

            const auto dl_geom_01_off = idx_geom_01_dlxx + i * dcomps + j;

            auto g_0_x_xx_xxxxxxxx = cbuffer.data(dl_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxxy = cbuffer.data(dl_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxxz = cbuffer.data(dl_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxyy = cbuffer.data(dl_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxyz = cbuffer.data(dl_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxzz = cbuffer.data(dl_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyyy = cbuffer.data(dl_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyyz = cbuffer.data(dl_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyzz = cbuffer.data(dl_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxzzz = cbuffer.data(dl_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyyy = cbuffer.data(dl_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyyz = cbuffer.data(dl_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyzz = cbuffer.data(dl_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyzzz = cbuffer.data(dl_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_xxxxzzzz = cbuffer.data(dl_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyyy = cbuffer.data(dl_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyyz = cbuffer.data(dl_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyzz = cbuffer.data(dl_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyzzz = cbuffer.data(dl_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xx_xxxyzzzz = cbuffer.data(dl_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xx_xxxzzzzz = cbuffer.data(dl_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyyy = cbuffer.data(dl_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyyz = cbuffer.data(dl_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyzz = cbuffer.data(dl_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyzzz = cbuffer.data(dl_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xx_xxyyzzzz = cbuffer.data(dl_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xx_xxyzzzzz = cbuffer.data(dl_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xx_xxzzzzzz = cbuffer.data(dl_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyyy = cbuffer.data(dl_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyyz = cbuffer.data(dl_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyzz = cbuffer.data(dl_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyzzz = cbuffer.data(dl_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xx_xyyyzzzz = cbuffer.data(dl_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xx_xyyzzzzz = cbuffer.data(dl_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xx_xyzzzzzz = cbuffer.data(dl_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xx_xzzzzzzz = cbuffer.data(dl_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyyy = cbuffer.data(dl_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyyz = cbuffer.data(dl_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyzz = cbuffer.data(dl_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyzzz = cbuffer.data(dl_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xx_yyyyzzzz = cbuffer.data(dl_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xx_yyyzzzzz = cbuffer.data(dl_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xx_yyzzzzzz = cbuffer.data(dl_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xx_yzzzzzzz = cbuffer.data(dl_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xx_zzzzzzzz = cbuffer.data(dl_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_yz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_yz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_yz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_yz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_yz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_yz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_yz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_yz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_yz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_yz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_zz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_zz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_zz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_zz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_zz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_zz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_zz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_zz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_zz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_zz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_zz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_zz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_zz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_zz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_zz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxxx = cbuffer.data(dl_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxxy = cbuffer.data(dl_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxxz = cbuffer.data(dl_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxyy = cbuffer.data(dl_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxyz = cbuffer.data(dl_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxzz = cbuffer.data(dl_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyyy = cbuffer.data(dl_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyyz = cbuffer.data(dl_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyzz = cbuffer.data(dl_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxzzz = cbuffer.data(dl_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyyy = cbuffer.data(dl_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyyz = cbuffer.data(dl_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyzz = cbuffer.data(dl_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyzzz = cbuffer.data(dl_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xx_xxxxzzzz = cbuffer.data(dl_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyyy = cbuffer.data(dl_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyyz = cbuffer.data(dl_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyzz = cbuffer.data(dl_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyzzz = cbuffer.data(dl_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xx_xxxyzzzz = cbuffer.data(dl_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xx_xxxzzzzz = cbuffer.data(dl_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyyy = cbuffer.data(dl_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyyz = cbuffer.data(dl_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyzz = cbuffer.data(dl_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyzzz = cbuffer.data(dl_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xx_xxyyzzzz = cbuffer.data(dl_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xx_xxyzzzzz = cbuffer.data(dl_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xx_xxzzzzzz = cbuffer.data(dl_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyyy = cbuffer.data(dl_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyyz = cbuffer.data(dl_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyzz = cbuffer.data(dl_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyzzz = cbuffer.data(dl_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xx_xyyyzzzz = cbuffer.data(dl_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xx_xyyzzzzz = cbuffer.data(dl_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xx_xyzzzzzz = cbuffer.data(dl_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xx_xzzzzzzz = cbuffer.data(dl_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyyy = cbuffer.data(dl_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyyz = cbuffer.data(dl_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyzz = cbuffer.data(dl_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyzzz = cbuffer.data(dl_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xx_yyyyzzzz = cbuffer.data(dl_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xx_yyyzzzzz = cbuffer.data(dl_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xx_yyzzzzzz = cbuffer.data(dl_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xx_yzzzzzzz = cbuffer.data(dl_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xx_zzzzzzzz = cbuffer.data(dl_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_yy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_yy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_yy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_yy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_yy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_yy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_yy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_yy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_yy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_yy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_yy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_yy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_yy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_yy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_yy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_yz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_yz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_yz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_yz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_yz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_yz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_yz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_yz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_yz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_yz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_yz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_yz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_yz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_yz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_yz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_zz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_zz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_zz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_zz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_zz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_zz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_zz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_zz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_zz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_zz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_zz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_zz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_zz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_zz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_zz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxxx = cbuffer.data(dl_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxxy = cbuffer.data(dl_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxxz = cbuffer.data(dl_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxyy = cbuffer.data(dl_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxyz = cbuffer.data(dl_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxzz = cbuffer.data(dl_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyyy = cbuffer.data(dl_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyyz = cbuffer.data(dl_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyzz = cbuffer.data(dl_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxzzz = cbuffer.data(dl_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyyy = cbuffer.data(dl_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyyz = cbuffer.data(dl_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyzz = cbuffer.data(dl_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyzzz = cbuffer.data(dl_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_xx_xxxxzzzz = cbuffer.data(dl_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyyy = cbuffer.data(dl_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyyz = cbuffer.data(dl_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyzz = cbuffer.data(dl_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyzzz = cbuffer.data(dl_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_xx_xxxyzzzz = cbuffer.data(dl_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xx_xxxzzzzz = cbuffer.data(dl_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyyy = cbuffer.data(dl_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyyz = cbuffer.data(dl_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyzz = cbuffer.data(dl_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyzzz = cbuffer.data(dl_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xx_xxyyzzzz = cbuffer.data(dl_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xx_xxyzzzzz = cbuffer.data(dl_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xx_xxzzzzzz = cbuffer.data(dl_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyyy = cbuffer.data(dl_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyyz = cbuffer.data(dl_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyzz = cbuffer.data(dl_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyzzz = cbuffer.data(dl_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xx_xyyyzzzz = cbuffer.data(dl_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xx_xyyzzzzz = cbuffer.data(dl_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xx_xyzzzzzz = cbuffer.data(dl_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xx_xzzzzzzz = cbuffer.data(dl_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyyy = cbuffer.data(dl_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyyz = cbuffer.data(dl_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyzz = cbuffer.data(dl_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyzzz = cbuffer.data(dl_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xx_yyyyzzzz = cbuffer.data(dl_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xx_yyyzzzzz = cbuffer.data(dl_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_xx_yyzzzzzz = cbuffer.data(dl_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xx_yzzzzzzz = cbuffer.data(dl_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xx_zzzzzzzz = cbuffer.data(dl_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_xy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_xy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_xy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_xy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_xy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_xy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_xy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_xy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_xy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_xy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_xy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_xy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_xy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_xy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxxx = cbuffer.data(dl_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxxy = cbuffer.data(dl_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxxz = cbuffer.data(dl_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxyy = cbuffer.data(dl_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxyz = cbuffer.data(dl_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxzz = cbuffer.data(dl_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyyy = cbuffer.data(dl_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyyz = cbuffer.data(dl_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyzz = cbuffer.data(dl_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxzzz = cbuffer.data(dl_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyyy = cbuffer.data(dl_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyyz = cbuffer.data(dl_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyzz = cbuffer.data(dl_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyzzz = cbuffer.data(dl_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_yy_xxxxzzzz = cbuffer.data(dl_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyyy = cbuffer.data(dl_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyyz = cbuffer.data(dl_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyzz = cbuffer.data(dl_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyzzz = cbuffer.data(dl_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_yy_xxxyzzzz = cbuffer.data(dl_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_yy_xxxzzzzz = cbuffer.data(dl_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyyy = cbuffer.data(dl_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyyz = cbuffer.data(dl_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyzz = cbuffer.data(dl_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyzzz = cbuffer.data(dl_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_yy_xxyyzzzz = cbuffer.data(dl_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_yy_xxyzzzzz = cbuffer.data(dl_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_yy_xxzzzzzz = cbuffer.data(dl_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyyy = cbuffer.data(dl_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyyz = cbuffer.data(dl_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyzz = cbuffer.data(dl_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyzzz = cbuffer.data(dl_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_yy_xyyyzzzz = cbuffer.data(dl_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_yy_xyyzzzzz = cbuffer.data(dl_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_yy_xyzzzzzz = cbuffer.data(dl_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_yy_xzzzzzzz = cbuffer.data(dl_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyyy = cbuffer.data(dl_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyyz = cbuffer.data(dl_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyzz = cbuffer.data(dl_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyzzz = cbuffer.data(dl_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_yy_yyyyzzzz = cbuffer.data(dl_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_yy_yyyzzzzz = cbuffer.data(dl_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_yy_yyzzzzzz = cbuffer.data(dl_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_yy_yzzzzzzz = cbuffer.data(dl_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_yy_zzzzzzzz = cbuffer.data(dl_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_yz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_yz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_yz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_yz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_yz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_yz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_yz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_yz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_yz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_yz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_yz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_yz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_yz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_yz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_yz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxxx = cbuffer.data(dl_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxxy = cbuffer.data(dl_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxxz = cbuffer.data(dl_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxyy = cbuffer.data(dl_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxyz = cbuffer.data(dl_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxzz = cbuffer.data(dl_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyyy = cbuffer.data(dl_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyyz = cbuffer.data(dl_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyzz = cbuffer.data(dl_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxzzz = cbuffer.data(dl_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyyy = cbuffer.data(dl_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyyz = cbuffer.data(dl_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyzz = cbuffer.data(dl_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyzzz = cbuffer.data(dl_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_zz_xxxxzzzz = cbuffer.data(dl_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyyy = cbuffer.data(dl_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyyz = cbuffer.data(dl_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyzz = cbuffer.data(dl_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyzzz = cbuffer.data(dl_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_zz_xxxyzzzz = cbuffer.data(dl_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_zz_xxxzzzzz = cbuffer.data(dl_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyyy = cbuffer.data(dl_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyyz = cbuffer.data(dl_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyzz = cbuffer.data(dl_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyzzz = cbuffer.data(dl_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_zz_xxyyzzzz = cbuffer.data(dl_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_zz_xxyzzzzz = cbuffer.data(dl_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_zz_xxzzzzzz = cbuffer.data(dl_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyyy = cbuffer.data(dl_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyyz = cbuffer.data(dl_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyzz = cbuffer.data(dl_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyzzz = cbuffer.data(dl_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_zz_xyyyzzzz = cbuffer.data(dl_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_zz_xyyzzzzz = cbuffer.data(dl_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_zz_xyzzzzzz = cbuffer.data(dl_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_zz_xzzzzzzz = cbuffer.data(dl_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyyy = cbuffer.data(dl_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyyz = cbuffer.data(dl_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyzz = cbuffer.data(dl_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyzzz = cbuffer.data(dl_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_zz_yyyyzzzz = cbuffer.data(dl_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_zz_yyyzzzzz = cbuffer.data(dl_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_zz_yyzzzzzz = cbuffer.data(dl_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_zz_yzzzzzzz = cbuffer.data(dl_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_zz_zzzzzzzz = cbuffer.data(dl_geom_01_off + 809 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fkxx

            const auto fk_geom_01_off = idx_geom_01_fkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxxxxxx, g_0_x_xx_xxxxxxxx, g_0_x_xx_xxxxxxxy, g_0_x_xx_xxxxxxxz, g_0_x_xx_xxxxxxy, g_0_x_xx_xxxxxxyy, g_0_x_xx_xxxxxxyz, g_0_x_xx_xxxxxxz, g_0_x_xx_xxxxxxzz, g_0_x_xx_xxxxxyy, g_0_x_xx_xxxxxyyy, g_0_x_xx_xxxxxyyz, g_0_x_xx_xxxxxyz, g_0_x_xx_xxxxxyzz, g_0_x_xx_xxxxxzz, g_0_x_xx_xxxxxzzz, g_0_x_xx_xxxxyyy, g_0_x_xx_xxxxyyyy, g_0_x_xx_xxxxyyyz, g_0_x_xx_xxxxyyz, g_0_x_xx_xxxxyyzz, g_0_x_xx_xxxxyzz, g_0_x_xx_xxxxyzzz, g_0_x_xx_xxxxzzz, g_0_x_xx_xxxxzzzz, g_0_x_xx_xxxyyyy, g_0_x_xx_xxxyyyyy, g_0_x_xx_xxxyyyyz, g_0_x_xx_xxxyyyz, g_0_x_xx_xxxyyyzz, g_0_x_xx_xxxyyzz, g_0_x_xx_xxxyyzzz, g_0_x_xx_xxxyzzz, g_0_x_xx_xxxyzzzz, g_0_x_xx_xxxzzzz, g_0_x_xx_xxxzzzzz, g_0_x_xx_xxyyyyy, g_0_x_xx_xxyyyyyy, g_0_x_xx_xxyyyyyz, g_0_x_xx_xxyyyyz, g_0_x_xx_xxyyyyzz, g_0_x_xx_xxyyyzz, g_0_x_xx_xxyyyzzz, g_0_x_xx_xxyyzzz, g_0_x_xx_xxyyzzzz, g_0_x_xx_xxyzzzz, g_0_x_xx_xxyzzzzz, g_0_x_xx_xxzzzzz, g_0_x_xx_xxzzzzzz, g_0_x_xx_xyyyyyy, g_0_x_xx_xyyyyyyy, g_0_x_xx_xyyyyyyz, g_0_x_xx_xyyyyyz, g_0_x_xx_xyyyyyzz, g_0_x_xx_xyyyyzz, g_0_x_xx_xyyyyzzz, g_0_x_xx_xyyyzzz, g_0_x_xx_xyyyzzzz, g_0_x_xx_xyyzzzz, g_0_x_xx_xyyzzzzz, g_0_x_xx_xyzzzzz, g_0_x_xx_xyzzzzzz, g_0_x_xx_xzzzzzz, g_0_x_xx_xzzzzzzz, g_0_x_xx_yyyyyyy, g_0_x_xx_yyyyyyz, g_0_x_xx_yyyyyzz, g_0_x_xx_yyyyzzz, g_0_x_xx_yyyzzzz, g_0_x_xx_yyzzzzz, g_0_x_xx_yzzzzzz, g_0_x_xx_zzzzzzz, g_0_x_xxx_xxxxxxx, g_0_x_xxx_xxxxxxy, g_0_x_xxx_xxxxxxz, g_0_x_xxx_xxxxxyy, g_0_x_xxx_xxxxxyz, g_0_x_xxx_xxxxxzz, g_0_x_xxx_xxxxyyy, g_0_x_xxx_xxxxyyz, g_0_x_xxx_xxxxyzz, g_0_x_xxx_xxxxzzz, g_0_x_xxx_xxxyyyy, g_0_x_xxx_xxxyyyz, g_0_x_xxx_xxxyyzz, g_0_x_xxx_xxxyzzz, g_0_x_xxx_xxxzzzz, g_0_x_xxx_xxyyyyy, g_0_x_xxx_xxyyyyz, g_0_x_xxx_xxyyyzz, g_0_x_xxx_xxyyzzz, g_0_x_xxx_xxyzzzz, g_0_x_xxx_xxzzzzz, g_0_x_xxx_xyyyyyy, g_0_x_xxx_xyyyyyz, g_0_x_xxx_xyyyyzz, g_0_x_xxx_xyyyzzz, g_0_x_xxx_xyyzzzz, g_0_x_xxx_xyzzzzz, g_0_x_xxx_xzzzzzz, g_0_x_xxx_yyyyyyy, g_0_x_xxx_yyyyyyz, g_0_x_xxx_yyyyyzz, g_0_x_xxx_yyyyzzz, g_0_x_xxx_yyyzzzz, g_0_x_xxx_yyzzzzz, g_0_x_xxx_yzzzzzz, g_0_x_xxx_zzzzzzz, g_xx_xxxxxxx, g_xx_xxxxxxy, g_xx_xxxxxxz, g_xx_xxxxxyy, g_xx_xxxxxyz, g_xx_xxxxxzz, g_xx_xxxxyyy, g_xx_xxxxyyz, g_xx_xxxxyzz, g_xx_xxxxzzz, g_xx_xxxyyyy, g_xx_xxxyyyz, g_xx_xxxyyzz, g_xx_xxxyzzz, g_xx_xxxzzzz, g_xx_xxyyyyy, g_xx_xxyyyyz, g_xx_xxyyyzz, g_xx_xxyyzzz, g_xx_xxyzzzz, g_xx_xxzzzzz, g_xx_xyyyyyy, g_xx_xyyyyyz, g_xx_xyyyyzz, g_xx_xyyyzzz, g_xx_xyyzzzz, g_xx_xyzzzzz, g_xx_xzzzzzz, g_xx_yyyyyyy, g_xx_yyyyyyz, g_xx_yyyyyzz, g_xx_yyyyzzz, g_xx_yyyzzzz, g_xx_yyzzzzz, g_xx_yzzzzzz, g_xx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxx_xxxxxxx[k] = g_xx_xxxxxxx[k] - g_0_x_xx_xxxxxxx[k] * ab_x + g_0_x_xx_xxxxxxxx[k];

                g_0_x_xxx_xxxxxxy[k] = g_xx_xxxxxxy[k] - g_0_x_xx_xxxxxxy[k] * ab_x + g_0_x_xx_xxxxxxxy[k];

                g_0_x_xxx_xxxxxxz[k] = g_xx_xxxxxxz[k] - g_0_x_xx_xxxxxxz[k] * ab_x + g_0_x_xx_xxxxxxxz[k];

                g_0_x_xxx_xxxxxyy[k] = g_xx_xxxxxyy[k] - g_0_x_xx_xxxxxyy[k] * ab_x + g_0_x_xx_xxxxxxyy[k];

                g_0_x_xxx_xxxxxyz[k] = g_xx_xxxxxyz[k] - g_0_x_xx_xxxxxyz[k] * ab_x + g_0_x_xx_xxxxxxyz[k];

                g_0_x_xxx_xxxxxzz[k] = g_xx_xxxxxzz[k] - g_0_x_xx_xxxxxzz[k] * ab_x + g_0_x_xx_xxxxxxzz[k];

                g_0_x_xxx_xxxxyyy[k] = g_xx_xxxxyyy[k] - g_0_x_xx_xxxxyyy[k] * ab_x + g_0_x_xx_xxxxxyyy[k];

                g_0_x_xxx_xxxxyyz[k] = g_xx_xxxxyyz[k] - g_0_x_xx_xxxxyyz[k] * ab_x + g_0_x_xx_xxxxxyyz[k];

                g_0_x_xxx_xxxxyzz[k] = g_xx_xxxxyzz[k] - g_0_x_xx_xxxxyzz[k] * ab_x + g_0_x_xx_xxxxxyzz[k];

                g_0_x_xxx_xxxxzzz[k] = g_xx_xxxxzzz[k] - g_0_x_xx_xxxxzzz[k] * ab_x + g_0_x_xx_xxxxxzzz[k];

                g_0_x_xxx_xxxyyyy[k] = g_xx_xxxyyyy[k] - g_0_x_xx_xxxyyyy[k] * ab_x + g_0_x_xx_xxxxyyyy[k];

                g_0_x_xxx_xxxyyyz[k] = g_xx_xxxyyyz[k] - g_0_x_xx_xxxyyyz[k] * ab_x + g_0_x_xx_xxxxyyyz[k];

                g_0_x_xxx_xxxyyzz[k] = g_xx_xxxyyzz[k] - g_0_x_xx_xxxyyzz[k] * ab_x + g_0_x_xx_xxxxyyzz[k];

                g_0_x_xxx_xxxyzzz[k] = g_xx_xxxyzzz[k] - g_0_x_xx_xxxyzzz[k] * ab_x + g_0_x_xx_xxxxyzzz[k];

                g_0_x_xxx_xxxzzzz[k] = g_xx_xxxzzzz[k] - g_0_x_xx_xxxzzzz[k] * ab_x + g_0_x_xx_xxxxzzzz[k];

                g_0_x_xxx_xxyyyyy[k] = g_xx_xxyyyyy[k] - g_0_x_xx_xxyyyyy[k] * ab_x + g_0_x_xx_xxxyyyyy[k];

                g_0_x_xxx_xxyyyyz[k] = g_xx_xxyyyyz[k] - g_0_x_xx_xxyyyyz[k] * ab_x + g_0_x_xx_xxxyyyyz[k];

                g_0_x_xxx_xxyyyzz[k] = g_xx_xxyyyzz[k] - g_0_x_xx_xxyyyzz[k] * ab_x + g_0_x_xx_xxxyyyzz[k];

                g_0_x_xxx_xxyyzzz[k] = g_xx_xxyyzzz[k] - g_0_x_xx_xxyyzzz[k] * ab_x + g_0_x_xx_xxxyyzzz[k];

                g_0_x_xxx_xxyzzzz[k] = g_xx_xxyzzzz[k] - g_0_x_xx_xxyzzzz[k] * ab_x + g_0_x_xx_xxxyzzzz[k];

                g_0_x_xxx_xxzzzzz[k] = g_xx_xxzzzzz[k] - g_0_x_xx_xxzzzzz[k] * ab_x + g_0_x_xx_xxxzzzzz[k];

                g_0_x_xxx_xyyyyyy[k] = g_xx_xyyyyyy[k] - g_0_x_xx_xyyyyyy[k] * ab_x + g_0_x_xx_xxyyyyyy[k];

                g_0_x_xxx_xyyyyyz[k] = g_xx_xyyyyyz[k] - g_0_x_xx_xyyyyyz[k] * ab_x + g_0_x_xx_xxyyyyyz[k];

                g_0_x_xxx_xyyyyzz[k] = g_xx_xyyyyzz[k] - g_0_x_xx_xyyyyzz[k] * ab_x + g_0_x_xx_xxyyyyzz[k];

                g_0_x_xxx_xyyyzzz[k] = g_xx_xyyyzzz[k] - g_0_x_xx_xyyyzzz[k] * ab_x + g_0_x_xx_xxyyyzzz[k];

                g_0_x_xxx_xyyzzzz[k] = g_xx_xyyzzzz[k] - g_0_x_xx_xyyzzzz[k] * ab_x + g_0_x_xx_xxyyzzzz[k];

                g_0_x_xxx_xyzzzzz[k] = g_xx_xyzzzzz[k] - g_0_x_xx_xyzzzzz[k] * ab_x + g_0_x_xx_xxyzzzzz[k];

                g_0_x_xxx_xzzzzzz[k] = g_xx_xzzzzzz[k] - g_0_x_xx_xzzzzzz[k] * ab_x + g_0_x_xx_xxzzzzzz[k];

                g_0_x_xxx_yyyyyyy[k] = g_xx_yyyyyyy[k] - g_0_x_xx_yyyyyyy[k] * ab_x + g_0_x_xx_xyyyyyyy[k];

                g_0_x_xxx_yyyyyyz[k] = g_xx_yyyyyyz[k] - g_0_x_xx_yyyyyyz[k] * ab_x + g_0_x_xx_xyyyyyyz[k];

                g_0_x_xxx_yyyyyzz[k] = g_xx_yyyyyzz[k] - g_0_x_xx_yyyyyzz[k] * ab_x + g_0_x_xx_xyyyyyzz[k];

                g_0_x_xxx_yyyyzzz[k] = g_xx_yyyyzzz[k] - g_0_x_xx_yyyyzzz[k] * ab_x + g_0_x_xx_xyyyyzzz[k];

                g_0_x_xxx_yyyzzzz[k] = g_xx_yyyzzzz[k] - g_0_x_xx_yyyzzzz[k] * ab_x + g_0_x_xx_xyyyzzzz[k];

                g_0_x_xxx_yyzzzzz[k] = g_xx_yyzzzzz[k] - g_0_x_xx_yyzzzzz[k] * ab_x + g_0_x_xx_xyyzzzzz[k];

                g_0_x_xxx_yzzzzzz[k] = g_xx_yzzzzzz[k] - g_0_x_xx_yzzzzzz[k] * ab_x + g_0_x_xx_xyzzzzzz[k];

                g_0_x_xxx_zzzzzzz[k] = g_xx_zzzzzzz[k] - g_0_x_xx_zzzzzzz[k] * ab_x + g_0_x_xx_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxxxxxx, g_0_x_xx_xxxxxxxy, g_0_x_xx_xxxxxxy, g_0_x_xx_xxxxxxyy, g_0_x_xx_xxxxxxyz, g_0_x_xx_xxxxxxz, g_0_x_xx_xxxxxyy, g_0_x_xx_xxxxxyyy, g_0_x_xx_xxxxxyyz, g_0_x_xx_xxxxxyz, g_0_x_xx_xxxxxyzz, g_0_x_xx_xxxxxzz, g_0_x_xx_xxxxyyy, g_0_x_xx_xxxxyyyy, g_0_x_xx_xxxxyyyz, g_0_x_xx_xxxxyyz, g_0_x_xx_xxxxyyzz, g_0_x_xx_xxxxyzz, g_0_x_xx_xxxxyzzz, g_0_x_xx_xxxxzzz, g_0_x_xx_xxxyyyy, g_0_x_xx_xxxyyyyy, g_0_x_xx_xxxyyyyz, g_0_x_xx_xxxyyyz, g_0_x_xx_xxxyyyzz, g_0_x_xx_xxxyyzz, g_0_x_xx_xxxyyzzz, g_0_x_xx_xxxyzzz, g_0_x_xx_xxxyzzzz, g_0_x_xx_xxxzzzz, g_0_x_xx_xxyyyyy, g_0_x_xx_xxyyyyyy, g_0_x_xx_xxyyyyyz, g_0_x_xx_xxyyyyz, g_0_x_xx_xxyyyyzz, g_0_x_xx_xxyyyzz, g_0_x_xx_xxyyyzzz, g_0_x_xx_xxyyzzz, g_0_x_xx_xxyyzzzz, g_0_x_xx_xxyzzzz, g_0_x_xx_xxyzzzzz, g_0_x_xx_xxzzzzz, g_0_x_xx_xyyyyyy, g_0_x_xx_xyyyyyyy, g_0_x_xx_xyyyyyyz, g_0_x_xx_xyyyyyz, g_0_x_xx_xyyyyyzz, g_0_x_xx_xyyyyzz, g_0_x_xx_xyyyyzzz, g_0_x_xx_xyyyzzz, g_0_x_xx_xyyyzzzz, g_0_x_xx_xyyzzzz, g_0_x_xx_xyyzzzzz, g_0_x_xx_xyzzzzz, g_0_x_xx_xyzzzzzz, g_0_x_xx_xzzzzzz, g_0_x_xx_yyyyyyy, g_0_x_xx_yyyyyyyy, g_0_x_xx_yyyyyyyz, g_0_x_xx_yyyyyyz, g_0_x_xx_yyyyyyzz, g_0_x_xx_yyyyyzz, g_0_x_xx_yyyyyzzz, g_0_x_xx_yyyyzzz, g_0_x_xx_yyyyzzzz, g_0_x_xx_yyyzzzz, g_0_x_xx_yyyzzzzz, g_0_x_xx_yyzzzzz, g_0_x_xx_yyzzzzzz, g_0_x_xx_yzzzzzz, g_0_x_xx_yzzzzzzz, g_0_x_xx_zzzzzzz, g_0_x_xxy_xxxxxxx, g_0_x_xxy_xxxxxxy, g_0_x_xxy_xxxxxxz, g_0_x_xxy_xxxxxyy, g_0_x_xxy_xxxxxyz, g_0_x_xxy_xxxxxzz, g_0_x_xxy_xxxxyyy, g_0_x_xxy_xxxxyyz, g_0_x_xxy_xxxxyzz, g_0_x_xxy_xxxxzzz, g_0_x_xxy_xxxyyyy, g_0_x_xxy_xxxyyyz, g_0_x_xxy_xxxyyzz, g_0_x_xxy_xxxyzzz, g_0_x_xxy_xxxzzzz, g_0_x_xxy_xxyyyyy, g_0_x_xxy_xxyyyyz, g_0_x_xxy_xxyyyzz, g_0_x_xxy_xxyyzzz, g_0_x_xxy_xxyzzzz, g_0_x_xxy_xxzzzzz, g_0_x_xxy_xyyyyyy, g_0_x_xxy_xyyyyyz, g_0_x_xxy_xyyyyzz, g_0_x_xxy_xyyyzzz, g_0_x_xxy_xyyzzzz, g_0_x_xxy_xyzzzzz, g_0_x_xxy_xzzzzzz, g_0_x_xxy_yyyyyyy, g_0_x_xxy_yyyyyyz, g_0_x_xxy_yyyyyzz, g_0_x_xxy_yyyyzzz, g_0_x_xxy_yyyzzzz, g_0_x_xxy_yyzzzzz, g_0_x_xxy_yzzzzzz, g_0_x_xxy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxy_xxxxxxx[k] = -g_0_x_xx_xxxxxxx[k] * ab_y + g_0_x_xx_xxxxxxxy[k];

                g_0_x_xxy_xxxxxxy[k] = -g_0_x_xx_xxxxxxy[k] * ab_y + g_0_x_xx_xxxxxxyy[k];

                g_0_x_xxy_xxxxxxz[k] = -g_0_x_xx_xxxxxxz[k] * ab_y + g_0_x_xx_xxxxxxyz[k];

                g_0_x_xxy_xxxxxyy[k] = -g_0_x_xx_xxxxxyy[k] * ab_y + g_0_x_xx_xxxxxyyy[k];

                g_0_x_xxy_xxxxxyz[k] = -g_0_x_xx_xxxxxyz[k] * ab_y + g_0_x_xx_xxxxxyyz[k];

                g_0_x_xxy_xxxxxzz[k] = -g_0_x_xx_xxxxxzz[k] * ab_y + g_0_x_xx_xxxxxyzz[k];

                g_0_x_xxy_xxxxyyy[k] = -g_0_x_xx_xxxxyyy[k] * ab_y + g_0_x_xx_xxxxyyyy[k];

                g_0_x_xxy_xxxxyyz[k] = -g_0_x_xx_xxxxyyz[k] * ab_y + g_0_x_xx_xxxxyyyz[k];

                g_0_x_xxy_xxxxyzz[k] = -g_0_x_xx_xxxxyzz[k] * ab_y + g_0_x_xx_xxxxyyzz[k];

                g_0_x_xxy_xxxxzzz[k] = -g_0_x_xx_xxxxzzz[k] * ab_y + g_0_x_xx_xxxxyzzz[k];

                g_0_x_xxy_xxxyyyy[k] = -g_0_x_xx_xxxyyyy[k] * ab_y + g_0_x_xx_xxxyyyyy[k];

                g_0_x_xxy_xxxyyyz[k] = -g_0_x_xx_xxxyyyz[k] * ab_y + g_0_x_xx_xxxyyyyz[k];

                g_0_x_xxy_xxxyyzz[k] = -g_0_x_xx_xxxyyzz[k] * ab_y + g_0_x_xx_xxxyyyzz[k];

                g_0_x_xxy_xxxyzzz[k] = -g_0_x_xx_xxxyzzz[k] * ab_y + g_0_x_xx_xxxyyzzz[k];

                g_0_x_xxy_xxxzzzz[k] = -g_0_x_xx_xxxzzzz[k] * ab_y + g_0_x_xx_xxxyzzzz[k];

                g_0_x_xxy_xxyyyyy[k] = -g_0_x_xx_xxyyyyy[k] * ab_y + g_0_x_xx_xxyyyyyy[k];

                g_0_x_xxy_xxyyyyz[k] = -g_0_x_xx_xxyyyyz[k] * ab_y + g_0_x_xx_xxyyyyyz[k];

                g_0_x_xxy_xxyyyzz[k] = -g_0_x_xx_xxyyyzz[k] * ab_y + g_0_x_xx_xxyyyyzz[k];

                g_0_x_xxy_xxyyzzz[k] = -g_0_x_xx_xxyyzzz[k] * ab_y + g_0_x_xx_xxyyyzzz[k];

                g_0_x_xxy_xxyzzzz[k] = -g_0_x_xx_xxyzzzz[k] * ab_y + g_0_x_xx_xxyyzzzz[k];

                g_0_x_xxy_xxzzzzz[k] = -g_0_x_xx_xxzzzzz[k] * ab_y + g_0_x_xx_xxyzzzzz[k];

                g_0_x_xxy_xyyyyyy[k] = -g_0_x_xx_xyyyyyy[k] * ab_y + g_0_x_xx_xyyyyyyy[k];

                g_0_x_xxy_xyyyyyz[k] = -g_0_x_xx_xyyyyyz[k] * ab_y + g_0_x_xx_xyyyyyyz[k];

                g_0_x_xxy_xyyyyzz[k] = -g_0_x_xx_xyyyyzz[k] * ab_y + g_0_x_xx_xyyyyyzz[k];

                g_0_x_xxy_xyyyzzz[k] = -g_0_x_xx_xyyyzzz[k] * ab_y + g_0_x_xx_xyyyyzzz[k];

                g_0_x_xxy_xyyzzzz[k] = -g_0_x_xx_xyyzzzz[k] * ab_y + g_0_x_xx_xyyyzzzz[k];

                g_0_x_xxy_xyzzzzz[k] = -g_0_x_xx_xyzzzzz[k] * ab_y + g_0_x_xx_xyyzzzzz[k];

                g_0_x_xxy_xzzzzzz[k] = -g_0_x_xx_xzzzzzz[k] * ab_y + g_0_x_xx_xyzzzzzz[k];

                g_0_x_xxy_yyyyyyy[k] = -g_0_x_xx_yyyyyyy[k] * ab_y + g_0_x_xx_yyyyyyyy[k];

                g_0_x_xxy_yyyyyyz[k] = -g_0_x_xx_yyyyyyz[k] * ab_y + g_0_x_xx_yyyyyyyz[k];

                g_0_x_xxy_yyyyyzz[k] = -g_0_x_xx_yyyyyzz[k] * ab_y + g_0_x_xx_yyyyyyzz[k];

                g_0_x_xxy_yyyyzzz[k] = -g_0_x_xx_yyyyzzz[k] * ab_y + g_0_x_xx_yyyyyzzz[k];

                g_0_x_xxy_yyyzzzz[k] = -g_0_x_xx_yyyzzzz[k] * ab_y + g_0_x_xx_yyyyzzzz[k];

                g_0_x_xxy_yyzzzzz[k] = -g_0_x_xx_yyzzzzz[k] * ab_y + g_0_x_xx_yyyzzzzz[k];

                g_0_x_xxy_yzzzzzz[k] = -g_0_x_xx_yzzzzzz[k] * ab_y + g_0_x_xx_yyzzzzzz[k];

                g_0_x_xxy_zzzzzzz[k] = -g_0_x_xx_zzzzzzz[k] * ab_y + g_0_x_xx_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxxxxxx, g_0_x_xx_xxxxxxxz, g_0_x_xx_xxxxxxy, g_0_x_xx_xxxxxxyz, g_0_x_xx_xxxxxxz, g_0_x_xx_xxxxxxzz, g_0_x_xx_xxxxxyy, g_0_x_xx_xxxxxyyz, g_0_x_xx_xxxxxyz, g_0_x_xx_xxxxxyzz, g_0_x_xx_xxxxxzz, g_0_x_xx_xxxxxzzz, g_0_x_xx_xxxxyyy, g_0_x_xx_xxxxyyyz, g_0_x_xx_xxxxyyz, g_0_x_xx_xxxxyyzz, g_0_x_xx_xxxxyzz, g_0_x_xx_xxxxyzzz, g_0_x_xx_xxxxzzz, g_0_x_xx_xxxxzzzz, g_0_x_xx_xxxyyyy, g_0_x_xx_xxxyyyyz, g_0_x_xx_xxxyyyz, g_0_x_xx_xxxyyyzz, g_0_x_xx_xxxyyzz, g_0_x_xx_xxxyyzzz, g_0_x_xx_xxxyzzz, g_0_x_xx_xxxyzzzz, g_0_x_xx_xxxzzzz, g_0_x_xx_xxxzzzzz, g_0_x_xx_xxyyyyy, g_0_x_xx_xxyyyyyz, g_0_x_xx_xxyyyyz, g_0_x_xx_xxyyyyzz, g_0_x_xx_xxyyyzz, g_0_x_xx_xxyyyzzz, g_0_x_xx_xxyyzzz, g_0_x_xx_xxyyzzzz, g_0_x_xx_xxyzzzz, g_0_x_xx_xxyzzzzz, g_0_x_xx_xxzzzzz, g_0_x_xx_xxzzzzzz, g_0_x_xx_xyyyyyy, g_0_x_xx_xyyyyyyz, g_0_x_xx_xyyyyyz, g_0_x_xx_xyyyyyzz, g_0_x_xx_xyyyyzz, g_0_x_xx_xyyyyzzz, g_0_x_xx_xyyyzzz, g_0_x_xx_xyyyzzzz, g_0_x_xx_xyyzzzz, g_0_x_xx_xyyzzzzz, g_0_x_xx_xyzzzzz, g_0_x_xx_xyzzzzzz, g_0_x_xx_xzzzzzz, g_0_x_xx_xzzzzzzz, g_0_x_xx_yyyyyyy, g_0_x_xx_yyyyyyyz, g_0_x_xx_yyyyyyz, g_0_x_xx_yyyyyyzz, g_0_x_xx_yyyyyzz, g_0_x_xx_yyyyyzzz, g_0_x_xx_yyyyzzz, g_0_x_xx_yyyyzzzz, g_0_x_xx_yyyzzzz, g_0_x_xx_yyyzzzzz, g_0_x_xx_yyzzzzz, g_0_x_xx_yyzzzzzz, g_0_x_xx_yzzzzzz, g_0_x_xx_yzzzzzzz, g_0_x_xx_zzzzzzz, g_0_x_xx_zzzzzzzz, g_0_x_xxz_xxxxxxx, g_0_x_xxz_xxxxxxy, g_0_x_xxz_xxxxxxz, g_0_x_xxz_xxxxxyy, g_0_x_xxz_xxxxxyz, g_0_x_xxz_xxxxxzz, g_0_x_xxz_xxxxyyy, g_0_x_xxz_xxxxyyz, g_0_x_xxz_xxxxyzz, g_0_x_xxz_xxxxzzz, g_0_x_xxz_xxxyyyy, g_0_x_xxz_xxxyyyz, g_0_x_xxz_xxxyyzz, g_0_x_xxz_xxxyzzz, g_0_x_xxz_xxxzzzz, g_0_x_xxz_xxyyyyy, g_0_x_xxz_xxyyyyz, g_0_x_xxz_xxyyyzz, g_0_x_xxz_xxyyzzz, g_0_x_xxz_xxyzzzz, g_0_x_xxz_xxzzzzz, g_0_x_xxz_xyyyyyy, g_0_x_xxz_xyyyyyz, g_0_x_xxz_xyyyyzz, g_0_x_xxz_xyyyzzz, g_0_x_xxz_xyyzzzz, g_0_x_xxz_xyzzzzz, g_0_x_xxz_xzzzzzz, g_0_x_xxz_yyyyyyy, g_0_x_xxz_yyyyyyz, g_0_x_xxz_yyyyyzz, g_0_x_xxz_yyyyzzz, g_0_x_xxz_yyyzzzz, g_0_x_xxz_yyzzzzz, g_0_x_xxz_yzzzzzz, g_0_x_xxz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxz_xxxxxxx[k] = -g_0_x_xx_xxxxxxx[k] * ab_z + g_0_x_xx_xxxxxxxz[k];

                g_0_x_xxz_xxxxxxy[k] = -g_0_x_xx_xxxxxxy[k] * ab_z + g_0_x_xx_xxxxxxyz[k];

                g_0_x_xxz_xxxxxxz[k] = -g_0_x_xx_xxxxxxz[k] * ab_z + g_0_x_xx_xxxxxxzz[k];

                g_0_x_xxz_xxxxxyy[k] = -g_0_x_xx_xxxxxyy[k] * ab_z + g_0_x_xx_xxxxxyyz[k];

                g_0_x_xxz_xxxxxyz[k] = -g_0_x_xx_xxxxxyz[k] * ab_z + g_0_x_xx_xxxxxyzz[k];

                g_0_x_xxz_xxxxxzz[k] = -g_0_x_xx_xxxxxzz[k] * ab_z + g_0_x_xx_xxxxxzzz[k];

                g_0_x_xxz_xxxxyyy[k] = -g_0_x_xx_xxxxyyy[k] * ab_z + g_0_x_xx_xxxxyyyz[k];

                g_0_x_xxz_xxxxyyz[k] = -g_0_x_xx_xxxxyyz[k] * ab_z + g_0_x_xx_xxxxyyzz[k];

                g_0_x_xxz_xxxxyzz[k] = -g_0_x_xx_xxxxyzz[k] * ab_z + g_0_x_xx_xxxxyzzz[k];

                g_0_x_xxz_xxxxzzz[k] = -g_0_x_xx_xxxxzzz[k] * ab_z + g_0_x_xx_xxxxzzzz[k];

                g_0_x_xxz_xxxyyyy[k] = -g_0_x_xx_xxxyyyy[k] * ab_z + g_0_x_xx_xxxyyyyz[k];

                g_0_x_xxz_xxxyyyz[k] = -g_0_x_xx_xxxyyyz[k] * ab_z + g_0_x_xx_xxxyyyzz[k];

                g_0_x_xxz_xxxyyzz[k] = -g_0_x_xx_xxxyyzz[k] * ab_z + g_0_x_xx_xxxyyzzz[k];

                g_0_x_xxz_xxxyzzz[k] = -g_0_x_xx_xxxyzzz[k] * ab_z + g_0_x_xx_xxxyzzzz[k];

                g_0_x_xxz_xxxzzzz[k] = -g_0_x_xx_xxxzzzz[k] * ab_z + g_0_x_xx_xxxzzzzz[k];

                g_0_x_xxz_xxyyyyy[k] = -g_0_x_xx_xxyyyyy[k] * ab_z + g_0_x_xx_xxyyyyyz[k];

                g_0_x_xxz_xxyyyyz[k] = -g_0_x_xx_xxyyyyz[k] * ab_z + g_0_x_xx_xxyyyyzz[k];

                g_0_x_xxz_xxyyyzz[k] = -g_0_x_xx_xxyyyzz[k] * ab_z + g_0_x_xx_xxyyyzzz[k];

                g_0_x_xxz_xxyyzzz[k] = -g_0_x_xx_xxyyzzz[k] * ab_z + g_0_x_xx_xxyyzzzz[k];

                g_0_x_xxz_xxyzzzz[k] = -g_0_x_xx_xxyzzzz[k] * ab_z + g_0_x_xx_xxyzzzzz[k];

                g_0_x_xxz_xxzzzzz[k] = -g_0_x_xx_xxzzzzz[k] * ab_z + g_0_x_xx_xxzzzzzz[k];

                g_0_x_xxz_xyyyyyy[k] = -g_0_x_xx_xyyyyyy[k] * ab_z + g_0_x_xx_xyyyyyyz[k];

                g_0_x_xxz_xyyyyyz[k] = -g_0_x_xx_xyyyyyz[k] * ab_z + g_0_x_xx_xyyyyyzz[k];

                g_0_x_xxz_xyyyyzz[k] = -g_0_x_xx_xyyyyzz[k] * ab_z + g_0_x_xx_xyyyyzzz[k];

                g_0_x_xxz_xyyyzzz[k] = -g_0_x_xx_xyyyzzz[k] * ab_z + g_0_x_xx_xyyyzzzz[k];

                g_0_x_xxz_xyyzzzz[k] = -g_0_x_xx_xyyzzzz[k] * ab_z + g_0_x_xx_xyyzzzzz[k];

                g_0_x_xxz_xyzzzzz[k] = -g_0_x_xx_xyzzzzz[k] * ab_z + g_0_x_xx_xyzzzzzz[k];

                g_0_x_xxz_xzzzzzz[k] = -g_0_x_xx_xzzzzzz[k] * ab_z + g_0_x_xx_xzzzzzzz[k];

                g_0_x_xxz_yyyyyyy[k] = -g_0_x_xx_yyyyyyy[k] * ab_z + g_0_x_xx_yyyyyyyz[k];

                g_0_x_xxz_yyyyyyz[k] = -g_0_x_xx_yyyyyyz[k] * ab_z + g_0_x_xx_yyyyyyzz[k];

                g_0_x_xxz_yyyyyzz[k] = -g_0_x_xx_yyyyyzz[k] * ab_z + g_0_x_xx_yyyyyzzz[k];

                g_0_x_xxz_yyyyzzz[k] = -g_0_x_xx_yyyyzzz[k] * ab_z + g_0_x_xx_yyyyzzzz[k];

                g_0_x_xxz_yyyzzzz[k] = -g_0_x_xx_yyyzzzz[k] * ab_z + g_0_x_xx_yyyzzzzz[k];

                g_0_x_xxz_yyzzzzz[k] = -g_0_x_xx_yyzzzzz[k] * ab_z + g_0_x_xx_yyzzzzzz[k];

                g_0_x_xxz_yzzzzzz[k] = -g_0_x_xx_yzzzzzz[k] * ab_z + g_0_x_xx_yzzzzzzz[k];

                g_0_x_xxz_zzzzzzz[k] = -g_0_x_xx_zzzzzzz[k] * ab_z + g_0_x_xx_zzzzzzzz[k];
            }

            /// Set up 108-144 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xy_xxxxxxx, g_0_x_xy_xxxxxxxy, g_0_x_xy_xxxxxxy, g_0_x_xy_xxxxxxyy, g_0_x_xy_xxxxxxyz, g_0_x_xy_xxxxxxz, g_0_x_xy_xxxxxyy, g_0_x_xy_xxxxxyyy, g_0_x_xy_xxxxxyyz, g_0_x_xy_xxxxxyz, g_0_x_xy_xxxxxyzz, g_0_x_xy_xxxxxzz, g_0_x_xy_xxxxyyy, g_0_x_xy_xxxxyyyy, g_0_x_xy_xxxxyyyz, g_0_x_xy_xxxxyyz, g_0_x_xy_xxxxyyzz, g_0_x_xy_xxxxyzz, g_0_x_xy_xxxxyzzz, g_0_x_xy_xxxxzzz, g_0_x_xy_xxxyyyy, g_0_x_xy_xxxyyyyy, g_0_x_xy_xxxyyyyz, g_0_x_xy_xxxyyyz, g_0_x_xy_xxxyyyzz, g_0_x_xy_xxxyyzz, g_0_x_xy_xxxyyzzz, g_0_x_xy_xxxyzzz, g_0_x_xy_xxxyzzzz, g_0_x_xy_xxxzzzz, g_0_x_xy_xxyyyyy, g_0_x_xy_xxyyyyyy, g_0_x_xy_xxyyyyyz, g_0_x_xy_xxyyyyz, g_0_x_xy_xxyyyyzz, g_0_x_xy_xxyyyzz, g_0_x_xy_xxyyyzzz, g_0_x_xy_xxyyzzz, g_0_x_xy_xxyyzzzz, g_0_x_xy_xxyzzzz, g_0_x_xy_xxyzzzzz, g_0_x_xy_xxzzzzz, g_0_x_xy_xyyyyyy, g_0_x_xy_xyyyyyyy, g_0_x_xy_xyyyyyyz, g_0_x_xy_xyyyyyz, g_0_x_xy_xyyyyyzz, g_0_x_xy_xyyyyzz, g_0_x_xy_xyyyyzzz, g_0_x_xy_xyyyzzz, g_0_x_xy_xyyyzzzz, g_0_x_xy_xyyzzzz, g_0_x_xy_xyyzzzzz, g_0_x_xy_xyzzzzz, g_0_x_xy_xyzzzzzz, g_0_x_xy_xzzzzzz, g_0_x_xy_yyyyyyy, g_0_x_xy_yyyyyyyy, g_0_x_xy_yyyyyyyz, g_0_x_xy_yyyyyyz, g_0_x_xy_yyyyyyzz, g_0_x_xy_yyyyyzz, g_0_x_xy_yyyyyzzz, g_0_x_xy_yyyyzzz, g_0_x_xy_yyyyzzzz, g_0_x_xy_yyyzzzz, g_0_x_xy_yyyzzzzz, g_0_x_xy_yyzzzzz, g_0_x_xy_yyzzzzzz, g_0_x_xy_yzzzzzz, g_0_x_xy_yzzzzzzz, g_0_x_xy_zzzzzzz, g_0_x_xyy_xxxxxxx, g_0_x_xyy_xxxxxxy, g_0_x_xyy_xxxxxxz, g_0_x_xyy_xxxxxyy, g_0_x_xyy_xxxxxyz, g_0_x_xyy_xxxxxzz, g_0_x_xyy_xxxxyyy, g_0_x_xyy_xxxxyyz, g_0_x_xyy_xxxxyzz, g_0_x_xyy_xxxxzzz, g_0_x_xyy_xxxyyyy, g_0_x_xyy_xxxyyyz, g_0_x_xyy_xxxyyzz, g_0_x_xyy_xxxyzzz, g_0_x_xyy_xxxzzzz, g_0_x_xyy_xxyyyyy, g_0_x_xyy_xxyyyyz, g_0_x_xyy_xxyyyzz, g_0_x_xyy_xxyyzzz, g_0_x_xyy_xxyzzzz, g_0_x_xyy_xxzzzzz, g_0_x_xyy_xyyyyyy, g_0_x_xyy_xyyyyyz, g_0_x_xyy_xyyyyzz, g_0_x_xyy_xyyyzzz, g_0_x_xyy_xyyzzzz, g_0_x_xyy_xyzzzzz, g_0_x_xyy_xzzzzzz, g_0_x_xyy_yyyyyyy, g_0_x_xyy_yyyyyyz, g_0_x_xyy_yyyyyzz, g_0_x_xyy_yyyyzzz, g_0_x_xyy_yyyzzzz, g_0_x_xyy_yyzzzzz, g_0_x_xyy_yzzzzzz, g_0_x_xyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyy_xxxxxxx[k] = -g_0_x_xy_xxxxxxx[k] * ab_y + g_0_x_xy_xxxxxxxy[k];

                g_0_x_xyy_xxxxxxy[k] = -g_0_x_xy_xxxxxxy[k] * ab_y + g_0_x_xy_xxxxxxyy[k];

                g_0_x_xyy_xxxxxxz[k] = -g_0_x_xy_xxxxxxz[k] * ab_y + g_0_x_xy_xxxxxxyz[k];

                g_0_x_xyy_xxxxxyy[k] = -g_0_x_xy_xxxxxyy[k] * ab_y + g_0_x_xy_xxxxxyyy[k];

                g_0_x_xyy_xxxxxyz[k] = -g_0_x_xy_xxxxxyz[k] * ab_y + g_0_x_xy_xxxxxyyz[k];

                g_0_x_xyy_xxxxxzz[k] = -g_0_x_xy_xxxxxzz[k] * ab_y + g_0_x_xy_xxxxxyzz[k];

                g_0_x_xyy_xxxxyyy[k] = -g_0_x_xy_xxxxyyy[k] * ab_y + g_0_x_xy_xxxxyyyy[k];

                g_0_x_xyy_xxxxyyz[k] = -g_0_x_xy_xxxxyyz[k] * ab_y + g_0_x_xy_xxxxyyyz[k];

                g_0_x_xyy_xxxxyzz[k] = -g_0_x_xy_xxxxyzz[k] * ab_y + g_0_x_xy_xxxxyyzz[k];

                g_0_x_xyy_xxxxzzz[k] = -g_0_x_xy_xxxxzzz[k] * ab_y + g_0_x_xy_xxxxyzzz[k];

                g_0_x_xyy_xxxyyyy[k] = -g_0_x_xy_xxxyyyy[k] * ab_y + g_0_x_xy_xxxyyyyy[k];

                g_0_x_xyy_xxxyyyz[k] = -g_0_x_xy_xxxyyyz[k] * ab_y + g_0_x_xy_xxxyyyyz[k];

                g_0_x_xyy_xxxyyzz[k] = -g_0_x_xy_xxxyyzz[k] * ab_y + g_0_x_xy_xxxyyyzz[k];

                g_0_x_xyy_xxxyzzz[k] = -g_0_x_xy_xxxyzzz[k] * ab_y + g_0_x_xy_xxxyyzzz[k];

                g_0_x_xyy_xxxzzzz[k] = -g_0_x_xy_xxxzzzz[k] * ab_y + g_0_x_xy_xxxyzzzz[k];

                g_0_x_xyy_xxyyyyy[k] = -g_0_x_xy_xxyyyyy[k] * ab_y + g_0_x_xy_xxyyyyyy[k];

                g_0_x_xyy_xxyyyyz[k] = -g_0_x_xy_xxyyyyz[k] * ab_y + g_0_x_xy_xxyyyyyz[k];

                g_0_x_xyy_xxyyyzz[k] = -g_0_x_xy_xxyyyzz[k] * ab_y + g_0_x_xy_xxyyyyzz[k];

                g_0_x_xyy_xxyyzzz[k] = -g_0_x_xy_xxyyzzz[k] * ab_y + g_0_x_xy_xxyyyzzz[k];

                g_0_x_xyy_xxyzzzz[k] = -g_0_x_xy_xxyzzzz[k] * ab_y + g_0_x_xy_xxyyzzzz[k];

                g_0_x_xyy_xxzzzzz[k] = -g_0_x_xy_xxzzzzz[k] * ab_y + g_0_x_xy_xxyzzzzz[k];

                g_0_x_xyy_xyyyyyy[k] = -g_0_x_xy_xyyyyyy[k] * ab_y + g_0_x_xy_xyyyyyyy[k];

                g_0_x_xyy_xyyyyyz[k] = -g_0_x_xy_xyyyyyz[k] * ab_y + g_0_x_xy_xyyyyyyz[k];

                g_0_x_xyy_xyyyyzz[k] = -g_0_x_xy_xyyyyzz[k] * ab_y + g_0_x_xy_xyyyyyzz[k];

                g_0_x_xyy_xyyyzzz[k] = -g_0_x_xy_xyyyzzz[k] * ab_y + g_0_x_xy_xyyyyzzz[k];

                g_0_x_xyy_xyyzzzz[k] = -g_0_x_xy_xyyzzzz[k] * ab_y + g_0_x_xy_xyyyzzzz[k];

                g_0_x_xyy_xyzzzzz[k] = -g_0_x_xy_xyzzzzz[k] * ab_y + g_0_x_xy_xyyzzzzz[k];

                g_0_x_xyy_xzzzzzz[k] = -g_0_x_xy_xzzzzzz[k] * ab_y + g_0_x_xy_xyzzzzzz[k];

                g_0_x_xyy_yyyyyyy[k] = -g_0_x_xy_yyyyyyy[k] * ab_y + g_0_x_xy_yyyyyyyy[k];

                g_0_x_xyy_yyyyyyz[k] = -g_0_x_xy_yyyyyyz[k] * ab_y + g_0_x_xy_yyyyyyyz[k];

                g_0_x_xyy_yyyyyzz[k] = -g_0_x_xy_yyyyyzz[k] * ab_y + g_0_x_xy_yyyyyyzz[k];

                g_0_x_xyy_yyyyzzz[k] = -g_0_x_xy_yyyyzzz[k] * ab_y + g_0_x_xy_yyyyyzzz[k];

                g_0_x_xyy_yyyzzzz[k] = -g_0_x_xy_yyyzzzz[k] * ab_y + g_0_x_xy_yyyyzzzz[k];

                g_0_x_xyy_yyzzzzz[k] = -g_0_x_xy_yyzzzzz[k] * ab_y + g_0_x_xy_yyyzzzzz[k];

                g_0_x_xyy_yzzzzzz[k] = -g_0_x_xy_yzzzzzz[k] * ab_y + g_0_x_xy_yyzzzzzz[k];

                g_0_x_xyy_zzzzzzz[k] = -g_0_x_xy_zzzzzzz[k] * ab_y + g_0_x_xy_yzzzzzzz[k];
            }

            /// Set up 144-180 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyz_xxxxxxx, g_0_x_xyz_xxxxxxy, g_0_x_xyz_xxxxxxz, g_0_x_xyz_xxxxxyy, g_0_x_xyz_xxxxxyz, g_0_x_xyz_xxxxxzz, g_0_x_xyz_xxxxyyy, g_0_x_xyz_xxxxyyz, g_0_x_xyz_xxxxyzz, g_0_x_xyz_xxxxzzz, g_0_x_xyz_xxxyyyy, g_0_x_xyz_xxxyyyz, g_0_x_xyz_xxxyyzz, g_0_x_xyz_xxxyzzz, g_0_x_xyz_xxxzzzz, g_0_x_xyz_xxyyyyy, g_0_x_xyz_xxyyyyz, g_0_x_xyz_xxyyyzz, g_0_x_xyz_xxyyzzz, g_0_x_xyz_xxyzzzz, g_0_x_xyz_xxzzzzz, g_0_x_xyz_xyyyyyy, g_0_x_xyz_xyyyyyz, g_0_x_xyz_xyyyyzz, g_0_x_xyz_xyyyzzz, g_0_x_xyz_xyyzzzz, g_0_x_xyz_xyzzzzz, g_0_x_xyz_xzzzzzz, g_0_x_xyz_yyyyyyy, g_0_x_xyz_yyyyyyz, g_0_x_xyz_yyyyyzz, g_0_x_xyz_yyyyzzz, g_0_x_xyz_yyyzzzz, g_0_x_xyz_yyzzzzz, g_0_x_xyz_yzzzzzz, g_0_x_xyz_zzzzzzz, g_0_x_xz_xxxxxxx, g_0_x_xz_xxxxxxxy, g_0_x_xz_xxxxxxy, g_0_x_xz_xxxxxxyy, g_0_x_xz_xxxxxxyz, g_0_x_xz_xxxxxxz, g_0_x_xz_xxxxxyy, g_0_x_xz_xxxxxyyy, g_0_x_xz_xxxxxyyz, g_0_x_xz_xxxxxyz, g_0_x_xz_xxxxxyzz, g_0_x_xz_xxxxxzz, g_0_x_xz_xxxxyyy, g_0_x_xz_xxxxyyyy, g_0_x_xz_xxxxyyyz, g_0_x_xz_xxxxyyz, g_0_x_xz_xxxxyyzz, g_0_x_xz_xxxxyzz, g_0_x_xz_xxxxyzzz, g_0_x_xz_xxxxzzz, g_0_x_xz_xxxyyyy, g_0_x_xz_xxxyyyyy, g_0_x_xz_xxxyyyyz, g_0_x_xz_xxxyyyz, g_0_x_xz_xxxyyyzz, g_0_x_xz_xxxyyzz, g_0_x_xz_xxxyyzzz, g_0_x_xz_xxxyzzz, g_0_x_xz_xxxyzzzz, g_0_x_xz_xxxzzzz, g_0_x_xz_xxyyyyy, g_0_x_xz_xxyyyyyy, g_0_x_xz_xxyyyyyz, g_0_x_xz_xxyyyyz, g_0_x_xz_xxyyyyzz, g_0_x_xz_xxyyyzz, g_0_x_xz_xxyyyzzz, g_0_x_xz_xxyyzzz, g_0_x_xz_xxyyzzzz, g_0_x_xz_xxyzzzz, g_0_x_xz_xxyzzzzz, g_0_x_xz_xxzzzzz, g_0_x_xz_xyyyyyy, g_0_x_xz_xyyyyyyy, g_0_x_xz_xyyyyyyz, g_0_x_xz_xyyyyyz, g_0_x_xz_xyyyyyzz, g_0_x_xz_xyyyyzz, g_0_x_xz_xyyyyzzz, g_0_x_xz_xyyyzzz, g_0_x_xz_xyyyzzzz, g_0_x_xz_xyyzzzz, g_0_x_xz_xyyzzzzz, g_0_x_xz_xyzzzzz, g_0_x_xz_xyzzzzzz, g_0_x_xz_xzzzzzz, g_0_x_xz_yyyyyyy, g_0_x_xz_yyyyyyyy, g_0_x_xz_yyyyyyyz, g_0_x_xz_yyyyyyz, g_0_x_xz_yyyyyyzz, g_0_x_xz_yyyyyzz, g_0_x_xz_yyyyyzzz, g_0_x_xz_yyyyzzz, g_0_x_xz_yyyyzzzz, g_0_x_xz_yyyzzzz, g_0_x_xz_yyyzzzzz, g_0_x_xz_yyzzzzz, g_0_x_xz_yyzzzzzz, g_0_x_xz_yzzzzzz, g_0_x_xz_yzzzzzzz, g_0_x_xz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyz_xxxxxxx[k] = -g_0_x_xz_xxxxxxx[k] * ab_y + g_0_x_xz_xxxxxxxy[k];

                g_0_x_xyz_xxxxxxy[k] = -g_0_x_xz_xxxxxxy[k] * ab_y + g_0_x_xz_xxxxxxyy[k];

                g_0_x_xyz_xxxxxxz[k] = -g_0_x_xz_xxxxxxz[k] * ab_y + g_0_x_xz_xxxxxxyz[k];

                g_0_x_xyz_xxxxxyy[k] = -g_0_x_xz_xxxxxyy[k] * ab_y + g_0_x_xz_xxxxxyyy[k];

                g_0_x_xyz_xxxxxyz[k] = -g_0_x_xz_xxxxxyz[k] * ab_y + g_0_x_xz_xxxxxyyz[k];

                g_0_x_xyz_xxxxxzz[k] = -g_0_x_xz_xxxxxzz[k] * ab_y + g_0_x_xz_xxxxxyzz[k];

                g_0_x_xyz_xxxxyyy[k] = -g_0_x_xz_xxxxyyy[k] * ab_y + g_0_x_xz_xxxxyyyy[k];

                g_0_x_xyz_xxxxyyz[k] = -g_0_x_xz_xxxxyyz[k] * ab_y + g_0_x_xz_xxxxyyyz[k];

                g_0_x_xyz_xxxxyzz[k] = -g_0_x_xz_xxxxyzz[k] * ab_y + g_0_x_xz_xxxxyyzz[k];

                g_0_x_xyz_xxxxzzz[k] = -g_0_x_xz_xxxxzzz[k] * ab_y + g_0_x_xz_xxxxyzzz[k];

                g_0_x_xyz_xxxyyyy[k] = -g_0_x_xz_xxxyyyy[k] * ab_y + g_0_x_xz_xxxyyyyy[k];

                g_0_x_xyz_xxxyyyz[k] = -g_0_x_xz_xxxyyyz[k] * ab_y + g_0_x_xz_xxxyyyyz[k];

                g_0_x_xyz_xxxyyzz[k] = -g_0_x_xz_xxxyyzz[k] * ab_y + g_0_x_xz_xxxyyyzz[k];

                g_0_x_xyz_xxxyzzz[k] = -g_0_x_xz_xxxyzzz[k] * ab_y + g_0_x_xz_xxxyyzzz[k];

                g_0_x_xyz_xxxzzzz[k] = -g_0_x_xz_xxxzzzz[k] * ab_y + g_0_x_xz_xxxyzzzz[k];

                g_0_x_xyz_xxyyyyy[k] = -g_0_x_xz_xxyyyyy[k] * ab_y + g_0_x_xz_xxyyyyyy[k];

                g_0_x_xyz_xxyyyyz[k] = -g_0_x_xz_xxyyyyz[k] * ab_y + g_0_x_xz_xxyyyyyz[k];

                g_0_x_xyz_xxyyyzz[k] = -g_0_x_xz_xxyyyzz[k] * ab_y + g_0_x_xz_xxyyyyzz[k];

                g_0_x_xyz_xxyyzzz[k] = -g_0_x_xz_xxyyzzz[k] * ab_y + g_0_x_xz_xxyyyzzz[k];

                g_0_x_xyz_xxyzzzz[k] = -g_0_x_xz_xxyzzzz[k] * ab_y + g_0_x_xz_xxyyzzzz[k];

                g_0_x_xyz_xxzzzzz[k] = -g_0_x_xz_xxzzzzz[k] * ab_y + g_0_x_xz_xxyzzzzz[k];

                g_0_x_xyz_xyyyyyy[k] = -g_0_x_xz_xyyyyyy[k] * ab_y + g_0_x_xz_xyyyyyyy[k];

                g_0_x_xyz_xyyyyyz[k] = -g_0_x_xz_xyyyyyz[k] * ab_y + g_0_x_xz_xyyyyyyz[k];

                g_0_x_xyz_xyyyyzz[k] = -g_0_x_xz_xyyyyzz[k] * ab_y + g_0_x_xz_xyyyyyzz[k];

                g_0_x_xyz_xyyyzzz[k] = -g_0_x_xz_xyyyzzz[k] * ab_y + g_0_x_xz_xyyyyzzz[k];

                g_0_x_xyz_xyyzzzz[k] = -g_0_x_xz_xyyzzzz[k] * ab_y + g_0_x_xz_xyyyzzzz[k];

                g_0_x_xyz_xyzzzzz[k] = -g_0_x_xz_xyzzzzz[k] * ab_y + g_0_x_xz_xyyzzzzz[k];

                g_0_x_xyz_xzzzzzz[k] = -g_0_x_xz_xzzzzzz[k] * ab_y + g_0_x_xz_xyzzzzzz[k];

                g_0_x_xyz_yyyyyyy[k] = -g_0_x_xz_yyyyyyy[k] * ab_y + g_0_x_xz_yyyyyyyy[k];

                g_0_x_xyz_yyyyyyz[k] = -g_0_x_xz_yyyyyyz[k] * ab_y + g_0_x_xz_yyyyyyyz[k];

                g_0_x_xyz_yyyyyzz[k] = -g_0_x_xz_yyyyyzz[k] * ab_y + g_0_x_xz_yyyyyyzz[k];

                g_0_x_xyz_yyyyzzz[k] = -g_0_x_xz_yyyyzzz[k] * ab_y + g_0_x_xz_yyyyyzzz[k];

                g_0_x_xyz_yyyzzzz[k] = -g_0_x_xz_yyyzzzz[k] * ab_y + g_0_x_xz_yyyyzzzz[k];

                g_0_x_xyz_yyzzzzz[k] = -g_0_x_xz_yyzzzzz[k] * ab_y + g_0_x_xz_yyyzzzzz[k];

                g_0_x_xyz_yzzzzzz[k] = -g_0_x_xz_yzzzzzz[k] * ab_y + g_0_x_xz_yyzzzzzz[k];

                g_0_x_xyz_zzzzzzz[k] = -g_0_x_xz_zzzzzzz[k] * ab_y + g_0_x_xz_yzzzzzzz[k];
            }

            /// Set up 180-216 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xz_xxxxxxx, g_0_x_xz_xxxxxxxz, g_0_x_xz_xxxxxxy, g_0_x_xz_xxxxxxyz, g_0_x_xz_xxxxxxz, g_0_x_xz_xxxxxxzz, g_0_x_xz_xxxxxyy, g_0_x_xz_xxxxxyyz, g_0_x_xz_xxxxxyz, g_0_x_xz_xxxxxyzz, g_0_x_xz_xxxxxzz, g_0_x_xz_xxxxxzzz, g_0_x_xz_xxxxyyy, g_0_x_xz_xxxxyyyz, g_0_x_xz_xxxxyyz, g_0_x_xz_xxxxyyzz, g_0_x_xz_xxxxyzz, g_0_x_xz_xxxxyzzz, g_0_x_xz_xxxxzzz, g_0_x_xz_xxxxzzzz, g_0_x_xz_xxxyyyy, g_0_x_xz_xxxyyyyz, g_0_x_xz_xxxyyyz, g_0_x_xz_xxxyyyzz, g_0_x_xz_xxxyyzz, g_0_x_xz_xxxyyzzz, g_0_x_xz_xxxyzzz, g_0_x_xz_xxxyzzzz, g_0_x_xz_xxxzzzz, g_0_x_xz_xxxzzzzz, g_0_x_xz_xxyyyyy, g_0_x_xz_xxyyyyyz, g_0_x_xz_xxyyyyz, g_0_x_xz_xxyyyyzz, g_0_x_xz_xxyyyzz, g_0_x_xz_xxyyyzzz, g_0_x_xz_xxyyzzz, g_0_x_xz_xxyyzzzz, g_0_x_xz_xxyzzzz, g_0_x_xz_xxyzzzzz, g_0_x_xz_xxzzzzz, g_0_x_xz_xxzzzzzz, g_0_x_xz_xyyyyyy, g_0_x_xz_xyyyyyyz, g_0_x_xz_xyyyyyz, g_0_x_xz_xyyyyyzz, g_0_x_xz_xyyyyzz, g_0_x_xz_xyyyyzzz, g_0_x_xz_xyyyzzz, g_0_x_xz_xyyyzzzz, g_0_x_xz_xyyzzzz, g_0_x_xz_xyyzzzzz, g_0_x_xz_xyzzzzz, g_0_x_xz_xyzzzzzz, g_0_x_xz_xzzzzzz, g_0_x_xz_xzzzzzzz, g_0_x_xz_yyyyyyy, g_0_x_xz_yyyyyyyz, g_0_x_xz_yyyyyyz, g_0_x_xz_yyyyyyzz, g_0_x_xz_yyyyyzz, g_0_x_xz_yyyyyzzz, g_0_x_xz_yyyyzzz, g_0_x_xz_yyyyzzzz, g_0_x_xz_yyyzzzz, g_0_x_xz_yyyzzzzz, g_0_x_xz_yyzzzzz, g_0_x_xz_yyzzzzzz, g_0_x_xz_yzzzzzz, g_0_x_xz_yzzzzzzz, g_0_x_xz_zzzzzzz, g_0_x_xz_zzzzzzzz, g_0_x_xzz_xxxxxxx, g_0_x_xzz_xxxxxxy, g_0_x_xzz_xxxxxxz, g_0_x_xzz_xxxxxyy, g_0_x_xzz_xxxxxyz, g_0_x_xzz_xxxxxzz, g_0_x_xzz_xxxxyyy, g_0_x_xzz_xxxxyyz, g_0_x_xzz_xxxxyzz, g_0_x_xzz_xxxxzzz, g_0_x_xzz_xxxyyyy, g_0_x_xzz_xxxyyyz, g_0_x_xzz_xxxyyzz, g_0_x_xzz_xxxyzzz, g_0_x_xzz_xxxzzzz, g_0_x_xzz_xxyyyyy, g_0_x_xzz_xxyyyyz, g_0_x_xzz_xxyyyzz, g_0_x_xzz_xxyyzzz, g_0_x_xzz_xxyzzzz, g_0_x_xzz_xxzzzzz, g_0_x_xzz_xyyyyyy, g_0_x_xzz_xyyyyyz, g_0_x_xzz_xyyyyzz, g_0_x_xzz_xyyyzzz, g_0_x_xzz_xyyzzzz, g_0_x_xzz_xyzzzzz, g_0_x_xzz_xzzzzzz, g_0_x_xzz_yyyyyyy, g_0_x_xzz_yyyyyyz, g_0_x_xzz_yyyyyzz, g_0_x_xzz_yyyyzzz, g_0_x_xzz_yyyzzzz, g_0_x_xzz_yyzzzzz, g_0_x_xzz_yzzzzzz, g_0_x_xzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzz_xxxxxxx[k] = -g_0_x_xz_xxxxxxx[k] * ab_z + g_0_x_xz_xxxxxxxz[k];

                g_0_x_xzz_xxxxxxy[k] = -g_0_x_xz_xxxxxxy[k] * ab_z + g_0_x_xz_xxxxxxyz[k];

                g_0_x_xzz_xxxxxxz[k] = -g_0_x_xz_xxxxxxz[k] * ab_z + g_0_x_xz_xxxxxxzz[k];

                g_0_x_xzz_xxxxxyy[k] = -g_0_x_xz_xxxxxyy[k] * ab_z + g_0_x_xz_xxxxxyyz[k];

                g_0_x_xzz_xxxxxyz[k] = -g_0_x_xz_xxxxxyz[k] * ab_z + g_0_x_xz_xxxxxyzz[k];

                g_0_x_xzz_xxxxxzz[k] = -g_0_x_xz_xxxxxzz[k] * ab_z + g_0_x_xz_xxxxxzzz[k];

                g_0_x_xzz_xxxxyyy[k] = -g_0_x_xz_xxxxyyy[k] * ab_z + g_0_x_xz_xxxxyyyz[k];

                g_0_x_xzz_xxxxyyz[k] = -g_0_x_xz_xxxxyyz[k] * ab_z + g_0_x_xz_xxxxyyzz[k];

                g_0_x_xzz_xxxxyzz[k] = -g_0_x_xz_xxxxyzz[k] * ab_z + g_0_x_xz_xxxxyzzz[k];

                g_0_x_xzz_xxxxzzz[k] = -g_0_x_xz_xxxxzzz[k] * ab_z + g_0_x_xz_xxxxzzzz[k];

                g_0_x_xzz_xxxyyyy[k] = -g_0_x_xz_xxxyyyy[k] * ab_z + g_0_x_xz_xxxyyyyz[k];

                g_0_x_xzz_xxxyyyz[k] = -g_0_x_xz_xxxyyyz[k] * ab_z + g_0_x_xz_xxxyyyzz[k];

                g_0_x_xzz_xxxyyzz[k] = -g_0_x_xz_xxxyyzz[k] * ab_z + g_0_x_xz_xxxyyzzz[k];

                g_0_x_xzz_xxxyzzz[k] = -g_0_x_xz_xxxyzzz[k] * ab_z + g_0_x_xz_xxxyzzzz[k];

                g_0_x_xzz_xxxzzzz[k] = -g_0_x_xz_xxxzzzz[k] * ab_z + g_0_x_xz_xxxzzzzz[k];

                g_0_x_xzz_xxyyyyy[k] = -g_0_x_xz_xxyyyyy[k] * ab_z + g_0_x_xz_xxyyyyyz[k];

                g_0_x_xzz_xxyyyyz[k] = -g_0_x_xz_xxyyyyz[k] * ab_z + g_0_x_xz_xxyyyyzz[k];

                g_0_x_xzz_xxyyyzz[k] = -g_0_x_xz_xxyyyzz[k] * ab_z + g_0_x_xz_xxyyyzzz[k];

                g_0_x_xzz_xxyyzzz[k] = -g_0_x_xz_xxyyzzz[k] * ab_z + g_0_x_xz_xxyyzzzz[k];

                g_0_x_xzz_xxyzzzz[k] = -g_0_x_xz_xxyzzzz[k] * ab_z + g_0_x_xz_xxyzzzzz[k];

                g_0_x_xzz_xxzzzzz[k] = -g_0_x_xz_xxzzzzz[k] * ab_z + g_0_x_xz_xxzzzzzz[k];

                g_0_x_xzz_xyyyyyy[k] = -g_0_x_xz_xyyyyyy[k] * ab_z + g_0_x_xz_xyyyyyyz[k];

                g_0_x_xzz_xyyyyyz[k] = -g_0_x_xz_xyyyyyz[k] * ab_z + g_0_x_xz_xyyyyyzz[k];

                g_0_x_xzz_xyyyyzz[k] = -g_0_x_xz_xyyyyzz[k] * ab_z + g_0_x_xz_xyyyyzzz[k];

                g_0_x_xzz_xyyyzzz[k] = -g_0_x_xz_xyyyzzz[k] * ab_z + g_0_x_xz_xyyyzzzz[k];

                g_0_x_xzz_xyyzzzz[k] = -g_0_x_xz_xyyzzzz[k] * ab_z + g_0_x_xz_xyyzzzzz[k];

                g_0_x_xzz_xyzzzzz[k] = -g_0_x_xz_xyzzzzz[k] * ab_z + g_0_x_xz_xyzzzzzz[k];

                g_0_x_xzz_xzzzzzz[k] = -g_0_x_xz_xzzzzzz[k] * ab_z + g_0_x_xz_xzzzzzzz[k];

                g_0_x_xzz_yyyyyyy[k] = -g_0_x_xz_yyyyyyy[k] * ab_z + g_0_x_xz_yyyyyyyz[k];

                g_0_x_xzz_yyyyyyz[k] = -g_0_x_xz_yyyyyyz[k] * ab_z + g_0_x_xz_yyyyyyzz[k];

                g_0_x_xzz_yyyyyzz[k] = -g_0_x_xz_yyyyyzz[k] * ab_z + g_0_x_xz_yyyyyzzz[k];

                g_0_x_xzz_yyyyzzz[k] = -g_0_x_xz_yyyyzzz[k] * ab_z + g_0_x_xz_yyyyzzzz[k];

                g_0_x_xzz_yyyzzzz[k] = -g_0_x_xz_yyyzzzz[k] * ab_z + g_0_x_xz_yyyzzzzz[k];

                g_0_x_xzz_yyzzzzz[k] = -g_0_x_xz_yyzzzzz[k] * ab_z + g_0_x_xz_yyzzzzzz[k];

                g_0_x_xzz_yzzzzzz[k] = -g_0_x_xz_yzzzzzz[k] * ab_z + g_0_x_xz_yzzzzzzz[k];

                g_0_x_xzz_zzzzzzz[k] = -g_0_x_xz_zzzzzzz[k] * ab_z + g_0_x_xz_zzzzzzzz[k];
            }

            /// Set up 216-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yy_xxxxxxx, g_0_x_yy_xxxxxxxy, g_0_x_yy_xxxxxxy, g_0_x_yy_xxxxxxyy, g_0_x_yy_xxxxxxyz, g_0_x_yy_xxxxxxz, g_0_x_yy_xxxxxyy, g_0_x_yy_xxxxxyyy, g_0_x_yy_xxxxxyyz, g_0_x_yy_xxxxxyz, g_0_x_yy_xxxxxyzz, g_0_x_yy_xxxxxzz, g_0_x_yy_xxxxyyy, g_0_x_yy_xxxxyyyy, g_0_x_yy_xxxxyyyz, g_0_x_yy_xxxxyyz, g_0_x_yy_xxxxyyzz, g_0_x_yy_xxxxyzz, g_0_x_yy_xxxxyzzz, g_0_x_yy_xxxxzzz, g_0_x_yy_xxxyyyy, g_0_x_yy_xxxyyyyy, g_0_x_yy_xxxyyyyz, g_0_x_yy_xxxyyyz, g_0_x_yy_xxxyyyzz, g_0_x_yy_xxxyyzz, g_0_x_yy_xxxyyzzz, g_0_x_yy_xxxyzzz, g_0_x_yy_xxxyzzzz, g_0_x_yy_xxxzzzz, g_0_x_yy_xxyyyyy, g_0_x_yy_xxyyyyyy, g_0_x_yy_xxyyyyyz, g_0_x_yy_xxyyyyz, g_0_x_yy_xxyyyyzz, g_0_x_yy_xxyyyzz, g_0_x_yy_xxyyyzzz, g_0_x_yy_xxyyzzz, g_0_x_yy_xxyyzzzz, g_0_x_yy_xxyzzzz, g_0_x_yy_xxyzzzzz, g_0_x_yy_xxzzzzz, g_0_x_yy_xyyyyyy, g_0_x_yy_xyyyyyyy, g_0_x_yy_xyyyyyyz, g_0_x_yy_xyyyyyz, g_0_x_yy_xyyyyyzz, g_0_x_yy_xyyyyzz, g_0_x_yy_xyyyyzzz, g_0_x_yy_xyyyzzz, g_0_x_yy_xyyyzzzz, g_0_x_yy_xyyzzzz, g_0_x_yy_xyyzzzzz, g_0_x_yy_xyzzzzz, g_0_x_yy_xyzzzzzz, g_0_x_yy_xzzzzzz, g_0_x_yy_yyyyyyy, g_0_x_yy_yyyyyyyy, g_0_x_yy_yyyyyyyz, g_0_x_yy_yyyyyyz, g_0_x_yy_yyyyyyzz, g_0_x_yy_yyyyyzz, g_0_x_yy_yyyyyzzz, g_0_x_yy_yyyyzzz, g_0_x_yy_yyyyzzzz, g_0_x_yy_yyyzzzz, g_0_x_yy_yyyzzzzz, g_0_x_yy_yyzzzzz, g_0_x_yy_yyzzzzzz, g_0_x_yy_yzzzzzz, g_0_x_yy_yzzzzzzz, g_0_x_yy_zzzzzzz, g_0_x_yyy_xxxxxxx, g_0_x_yyy_xxxxxxy, g_0_x_yyy_xxxxxxz, g_0_x_yyy_xxxxxyy, g_0_x_yyy_xxxxxyz, g_0_x_yyy_xxxxxzz, g_0_x_yyy_xxxxyyy, g_0_x_yyy_xxxxyyz, g_0_x_yyy_xxxxyzz, g_0_x_yyy_xxxxzzz, g_0_x_yyy_xxxyyyy, g_0_x_yyy_xxxyyyz, g_0_x_yyy_xxxyyzz, g_0_x_yyy_xxxyzzz, g_0_x_yyy_xxxzzzz, g_0_x_yyy_xxyyyyy, g_0_x_yyy_xxyyyyz, g_0_x_yyy_xxyyyzz, g_0_x_yyy_xxyyzzz, g_0_x_yyy_xxyzzzz, g_0_x_yyy_xxzzzzz, g_0_x_yyy_xyyyyyy, g_0_x_yyy_xyyyyyz, g_0_x_yyy_xyyyyzz, g_0_x_yyy_xyyyzzz, g_0_x_yyy_xyyzzzz, g_0_x_yyy_xyzzzzz, g_0_x_yyy_xzzzzzz, g_0_x_yyy_yyyyyyy, g_0_x_yyy_yyyyyyz, g_0_x_yyy_yyyyyzz, g_0_x_yyy_yyyyzzz, g_0_x_yyy_yyyzzzz, g_0_x_yyy_yyzzzzz, g_0_x_yyy_yzzzzzz, g_0_x_yyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyy_xxxxxxx[k] = -g_0_x_yy_xxxxxxx[k] * ab_y + g_0_x_yy_xxxxxxxy[k];

                g_0_x_yyy_xxxxxxy[k] = -g_0_x_yy_xxxxxxy[k] * ab_y + g_0_x_yy_xxxxxxyy[k];

                g_0_x_yyy_xxxxxxz[k] = -g_0_x_yy_xxxxxxz[k] * ab_y + g_0_x_yy_xxxxxxyz[k];

                g_0_x_yyy_xxxxxyy[k] = -g_0_x_yy_xxxxxyy[k] * ab_y + g_0_x_yy_xxxxxyyy[k];

                g_0_x_yyy_xxxxxyz[k] = -g_0_x_yy_xxxxxyz[k] * ab_y + g_0_x_yy_xxxxxyyz[k];

                g_0_x_yyy_xxxxxzz[k] = -g_0_x_yy_xxxxxzz[k] * ab_y + g_0_x_yy_xxxxxyzz[k];

                g_0_x_yyy_xxxxyyy[k] = -g_0_x_yy_xxxxyyy[k] * ab_y + g_0_x_yy_xxxxyyyy[k];

                g_0_x_yyy_xxxxyyz[k] = -g_0_x_yy_xxxxyyz[k] * ab_y + g_0_x_yy_xxxxyyyz[k];

                g_0_x_yyy_xxxxyzz[k] = -g_0_x_yy_xxxxyzz[k] * ab_y + g_0_x_yy_xxxxyyzz[k];

                g_0_x_yyy_xxxxzzz[k] = -g_0_x_yy_xxxxzzz[k] * ab_y + g_0_x_yy_xxxxyzzz[k];

                g_0_x_yyy_xxxyyyy[k] = -g_0_x_yy_xxxyyyy[k] * ab_y + g_0_x_yy_xxxyyyyy[k];

                g_0_x_yyy_xxxyyyz[k] = -g_0_x_yy_xxxyyyz[k] * ab_y + g_0_x_yy_xxxyyyyz[k];

                g_0_x_yyy_xxxyyzz[k] = -g_0_x_yy_xxxyyzz[k] * ab_y + g_0_x_yy_xxxyyyzz[k];

                g_0_x_yyy_xxxyzzz[k] = -g_0_x_yy_xxxyzzz[k] * ab_y + g_0_x_yy_xxxyyzzz[k];

                g_0_x_yyy_xxxzzzz[k] = -g_0_x_yy_xxxzzzz[k] * ab_y + g_0_x_yy_xxxyzzzz[k];

                g_0_x_yyy_xxyyyyy[k] = -g_0_x_yy_xxyyyyy[k] * ab_y + g_0_x_yy_xxyyyyyy[k];

                g_0_x_yyy_xxyyyyz[k] = -g_0_x_yy_xxyyyyz[k] * ab_y + g_0_x_yy_xxyyyyyz[k];

                g_0_x_yyy_xxyyyzz[k] = -g_0_x_yy_xxyyyzz[k] * ab_y + g_0_x_yy_xxyyyyzz[k];

                g_0_x_yyy_xxyyzzz[k] = -g_0_x_yy_xxyyzzz[k] * ab_y + g_0_x_yy_xxyyyzzz[k];

                g_0_x_yyy_xxyzzzz[k] = -g_0_x_yy_xxyzzzz[k] * ab_y + g_0_x_yy_xxyyzzzz[k];

                g_0_x_yyy_xxzzzzz[k] = -g_0_x_yy_xxzzzzz[k] * ab_y + g_0_x_yy_xxyzzzzz[k];

                g_0_x_yyy_xyyyyyy[k] = -g_0_x_yy_xyyyyyy[k] * ab_y + g_0_x_yy_xyyyyyyy[k];

                g_0_x_yyy_xyyyyyz[k] = -g_0_x_yy_xyyyyyz[k] * ab_y + g_0_x_yy_xyyyyyyz[k];

                g_0_x_yyy_xyyyyzz[k] = -g_0_x_yy_xyyyyzz[k] * ab_y + g_0_x_yy_xyyyyyzz[k];

                g_0_x_yyy_xyyyzzz[k] = -g_0_x_yy_xyyyzzz[k] * ab_y + g_0_x_yy_xyyyyzzz[k];

                g_0_x_yyy_xyyzzzz[k] = -g_0_x_yy_xyyzzzz[k] * ab_y + g_0_x_yy_xyyyzzzz[k];

                g_0_x_yyy_xyzzzzz[k] = -g_0_x_yy_xyzzzzz[k] * ab_y + g_0_x_yy_xyyzzzzz[k];

                g_0_x_yyy_xzzzzzz[k] = -g_0_x_yy_xzzzzzz[k] * ab_y + g_0_x_yy_xyzzzzzz[k];

                g_0_x_yyy_yyyyyyy[k] = -g_0_x_yy_yyyyyyy[k] * ab_y + g_0_x_yy_yyyyyyyy[k];

                g_0_x_yyy_yyyyyyz[k] = -g_0_x_yy_yyyyyyz[k] * ab_y + g_0_x_yy_yyyyyyyz[k];

                g_0_x_yyy_yyyyyzz[k] = -g_0_x_yy_yyyyyzz[k] * ab_y + g_0_x_yy_yyyyyyzz[k];

                g_0_x_yyy_yyyyzzz[k] = -g_0_x_yy_yyyyzzz[k] * ab_y + g_0_x_yy_yyyyyzzz[k];

                g_0_x_yyy_yyyzzzz[k] = -g_0_x_yy_yyyzzzz[k] * ab_y + g_0_x_yy_yyyyzzzz[k];

                g_0_x_yyy_yyzzzzz[k] = -g_0_x_yy_yyzzzzz[k] * ab_y + g_0_x_yy_yyyzzzzz[k];

                g_0_x_yyy_yzzzzzz[k] = -g_0_x_yy_yzzzzzz[k] * ab_y + g_0_x_yy_yyzzzzzz[k];

                g_0_x_yyy_zzzzzzz[k] = -g_0_x_yy_zzzzzzz[k] * ab_y + g_0_x_yy_yzzzzzzz[k];
            }

            /// Set up 252-288 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyz_xxxxxxx, g_0_x_yyz_xxxxxxy, g_0_x_yyz_xxxxxxz, g_0_x_yyz_xxxxxyy, g_0_x_yyz_xxxxxyz, g_0_x_yyz_xxxxxzz, g_0_x_yyz_xxxxyyy, g_0_x_yyz_xxxxyyz, g_0_x_yyz_xxxxyzz, g_0_x_yyz_xxxxzzz, g_0_x_yyz_xxxyyyy, g_0_x_yyz_xxxyyyz, g_0_x_yyz_xxxyyzz, g_0_x_yyz_xxxyzzz, g_0_x_yyz_xxxzzzz, g_0_x_yyz_xxyyyyy, g_0_x_yyz_xxyyyyz, g_0_x_yyz_xxyyyzz, g_0_x_yyz_xxyyzzz, g_0_x_yyz_xxyzzzz, g_0_x_yyz_xxzzzzz, g_0_x_yyz_xyyyyyy, g_0_x_yyz_xyyyyyz, g_0_x_yyz_xyyyyzz, g_0_x_yyz_xyyyzzz, g_0_x_yyz_xyyzzzz, g_0_x_yyz_xyzzzzz, g_0_x_yyz_xzzzzzz, g_0_x_yyz_yyyyyyy, g_0_x_yyz_yyyyyyz, g_0_x_yyz_yyyyyzz, g_0_x_yyz_yyyyzzz, g_0_x_yyz_yyyzzzz, g_0_x_yyz_yyzzzzz, g_0_x_yyz_yzzzzzz, g_0_x_yyz_zzzzzzz, g_0_x_yz_xxxxxxx, g_0_x_yz_xxxxxxxy, g_0_x_yz_xxxxxxy, g_0_x_yz_xxxxxxyy, g_0_x_yz_xxxxxxyz, g_0_x_yz_xxxxxxz, g_0_x_yz_xxxxxyy, g_0_x_yz_xxxxxyyy, g_0_x_yz_xxxxxyyz, g_0_x_yz_xxxxxyz, g_0_x_yz_xxxxxyzz, g_0_x_yz_xxxxxzz, g_0_x_yz_xxxxyyy, g_0_x_yz_xxxxyyyy, g_0_x_yz_xxxxyyyz, g_0_x_yz_xxxxyyz, g_0_x_yz_xxxxyyzz, g_0_x_yz_xxxxyzz, g_0_x_yz_xxxxyzzz, g_0_x_yz_xxxxzzz, g_0_x_yz_xxxyyyy, g_0_x_yz_xxxyyyyy, g_0_x_yz_xxxyyyyz, g_0_x_yz_xxxyyyz, g_0_x_yz_xxxyyyzz, g_0_x_yz_xxxyyzz, g_0_x_yz_xxxyyzzz, g_0_x_yz_xxxyzzz, g_0_x_yz_xxxyzzzz, g_0_x_yz_xxxzzzz, g_0_x_yz_xxyyyyy, g_0_x_yz_xxyyyyyy, g_0_x_yz_xxyyyyyz, g_0_x_yz_xxyyyyz, g_0_x_yz_xxyyyyzz, g_0_x_yz_xxyyyzz, g_0_x_yz_xxyyyzzz, g_0_x_yz_xxyyzzz, g_0_x_yz_xxyyzzzz, g_0_x_yz_xxyzzzz, g_0_x_yz_xxyzzzzz, g_0_x_yz_xxzzzzz, g_0_x_yz_xyyyyyy, g_0_x_yz_xyyyyyyy, g_0_x_yz_xyyyyyyz, g_0_x_yz_xyyyyyz, g_0_x_yz_xyyyyyzz, g_0_x_yz_xyyyyzz, g_0_x_yz_xyyyyzzz, g_0_x_yz_xyyyzzz, g_0_x_yz_xyyyzzzz, g_0_x_yz_xyyzzzz, g_0_x_yz_xyyzzzzz, g_0_x_yz_xyzzzzz, g_0_x_yz_xyzzzzzz, g_0_x_yz_xzzzzzz, g_0_x_yz_yyyyyyy, g_0_x_yz_yyyyyyyy, g_0_x_yz_yyyyyyyz, g_0_x_yz_yyyyyyz, g_0_x_yz_yyyyyyzz, g_0_x_yz_yyyyyzz, g_0_x_yz_yyyyyzzz, g_0_x_yz_yyyyzzz, g_0_x_yz_yyyyzzzz, g_0_x_yz_yyyzzzz, g_0_x_yz_yyyzzzzz, g_0_x_yz_yyzzzzz, g_0_x_yz_yyzzzzzz, g_0_x_yz_yzzzzzz, g_0_x_yz_yzzzzzzz, g_0_x_yz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyz_xxxxxxx[k] = -g_0_x_yz_xxxxxxx[k] * ab_y + g_0_x_yz_xxxxxxxy[k];

                g_0_x_yyz_xxxxxxy[k] = -g_0_x_yz_xxxxxxy[k] * ab_y + g_0_x_yz_xxxxxxyy[k];

                g_0_x_yyz_xxxxxxz[k] = -g_0_x_yz_xxxxxxz[k] * ab_y + g_0_x_yz_xxxxxxyz[k];

                g_0_x_yyz_xxxxxyy[k] = -g_0_x_yz_xxxxxyy[k] * ab_y + g_0_x_yz_xxxxxyyy[k];

                g_0_x_yyz_xxxxxyz[k] = -g_0_x_yz_xxxxxyz[k] * ab_y + g_0_x_yz_xxxxxyyz[k];

                g_0_x_yyz_xxxxxzz[k] = -g_0_x_yz_xxxxxzz[k] * ab_y + g_0_x_yz_xxxxxyzz[k];

                g_0_x_yyz_xxxxyyy[k] = -g_0_x_yz_xxxxyyy[k] * ab_y + g_0_x_yz_xxxxyyyy[k];

                g_0_x_yyz_xxxxyyz[k] = -g_0_x_yz_xxxxyyz[k] * ab_y + g_0_x_yz_xxxxyyyz[k];

                g_0_x_yyz_xxxxyzz[k] = -g_0_x_yz_xxxxyzz[k] * ab_y + g_0_x_yz_xxxxyyzz[k];

                g_0_x_yyz_xxxxzzz[k] = -g_0_x_yz_xxxxzzz[k] * ab_y + g_0_x_yz_xxxxyzzz[k];

                g_0_x_yyz_xxxyyyy[k] = -g_0_x_yz_xxxyyyy[k] * ab_y + g_0_x_yz_xxxyyyyy[k];

                g_0_x_yyz_xxxyyyz[k] = -g_0_x_yz_xxxyyyz[k] * ab_y + g_0_x_yz_xxxyyyyz[k];

                g_0_x_yyz_xxxyyzz[k] = -g_0_x_yz_xxxyyzz[k] * ab_y + g_0_x_yz_xxxyyyzz[k];

                g_0_x_yyz_xxxyzzz[k] = -g_0_x_yz_xxxyzzz[k] * ab_y + g_0_x_yz_xxxyyzzz[k];

                g_0_x_yyz_xxxzzzz[k] = -g_0_x_yz_xxxzzzz[k] * ab_y + g_0_x_yz_xxxyzzzz[k];

                g_0_x_yyz_xxyyyyy[k] = -g_0_x_yz_xxyyyyy[k] * ab_y + g_0_x_yz_xxyyyyyy[k];

                g_0_x_yyz_xxyyyyz[k] = -g_0_x_yz_xxyyyyz[k] * ab_y + g_0_x_yz_xxyyyyyz[k];

                g_0_x_yyz_xxyyyzz[k] = -g_0_x_yz_xxyyyzz[k] * ab_y + g_0_x_yz_xxyyyyzz[k];

                g_0_x_yyz_xxyyzzz[k] = -g_0_x_yz_xxyyzzz[k] * ab_y + g_0_x_yz_xxyyyzzz[k];

                g_0_x_yyz_xxyzzzz[k] = -g_0_x_yz_xxyzzzz[k] * ab_y + g_0_x_yz_xxyyzzzz[k];

                g_0_x_yyz_xxzzzzz[k] = -g_0_x_yz_xxzzzzz[k] * ab_y + g_0_x_yz_xxyzzzzz[k];

                g_0_x_yyz_xyyyyyy[k] = -g_0_x_yz_xyyyyyy[k] * ab_y + g_0_x_yz_xyyyyyyy[k];

                g_0_x_yyz_xyyyyyz[k] = -g_0_x_yz_xyyyyyz[k] * ab_y + g_0_x_yz_xyyyyyyz[k];

                g_0_x_yyz_xyyyyzz[k] = -g_0_x_yz_xyyyyzz[k] * ab_y + g_0_x_yz_xyyyyyzz[k];

                g_0_x_yyz_xyyyzzz[k] = -g_0_x_yz_xyyyzzz[k] * ab_y + g_0_x_yz_xyyyyzzz[k];

                g_0_x_yyz_xyyzzzz[k] = -g_0_x_yz_xyyzzzz[k] * ab_y + g_0_x_yz_xyyyzzzz[k];

                g_0_x_yyz_xyzzzzz[k] = -g_0_x_yz_xyzzzzz[k] * ab_y + g_0_x_yz_xyyzzzzz[k];

                g_0_x_yyz_xzzzzzz[k] = -g_0_x_yz_xzzzzzz[k] * ab_y + g_0_x_yz_xyzzzzzz[k];

                g_0_x_yyz_yyyyyyy[k] = -g_0_x_yz_yyyyyyy[k] * ab_y + g_0_x_yz_yyyyyyyy[k];

                g_0_x_yyz_yyyyyyz[k] = -g_0_x_yz_yyyyyyz[k] * ab_y + g_0_x_yz_yyyyyyyz[k];

                g_0_x_yyz_yyyyyzz[k] = -g_0_x_yz_yyyyyzz[k] * ab_y + g_0_x_yz_yyyyyyzz[k];

                g_0_x_yyz_yyyyzzz[k] = -g_0_x_yz_yyyyzzz[k] * ab_y + g_0_x_yz_yyyyyzzz[k];

                g_0_x_yyz_yyyzzzz[k] = -g_0_x_yz_yyyzzzz[k] * ab_y + g_0_x_yz_yyyyzzzz[k];

                g_0_x_yyz_yyzzzzz[k] = -g_0_x_yz_yyzzzzz[k] * ab_y + g_0_x_yz_yyyzzzzz[k];

                g_0_x_yyz_yzzzzzz[k] = -g_0_x_yz_yzzzzzz[k] * ab_y + g_0_x_yz_yyzzzzzz[k];

                g_0_x_yyz_zzzzzzz[k] = -g_0_x_yz_zzzzzzz[k] * ab_y + g_0_x_yz_yzzzzzzz[k];
            }

            /// Set up 288-324 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yzz_xxxxxxx, g_0_x_yzz_xxxxxxy, g_0_x_yzz_xxxxxxz, g_0_x_yzz_xxxxxyy, g_0_x_yzz_xxxxxyz, g_0_x_yzz_xxxxxzz, g_0_x_yzz_xxxxyyy, g_0_x_yzz_xxxxyyz, g_0_x_yzz_xxxxyzz, g_0_x_yzz_xxxxzzz, g_0_x_yzz_xxxyyyy, g_0_x_yzz_xxxyyyz, g_0_x_yzz_xxxyyzz, g_0_x_yzz_xxxyzzz, g_0_x_yzz_xxxzzzz, g_0_x_yzz_xxyyyyy, g_0_x_yzz_xxyyyyz, g_0_x_yzz_xxyyyzz, g_0_x_yzz_xxyyzzz, g_0_x_yzz_xxyzzzz, g_0_x_yzz_xxzzzzz, g_0_x_yzz_xyyyyyy, g_0_x_yzz_xyyyyyz, g_0_x_yzz_xyyyyzz, g_0_x_yzz_xyyyzzz, g_0_x_yzz_xyyzzzz, g_0_x_yzz_xyzzzzz, g_0_x_yzz_xzzzzzz, g_0_x_yzz_yyyyyyy, g_0_x_yzz_yyyyyyz, g_0_x_yzz_yyyyyzz, g_0_x_yzz_yyyyzzz, g_0_x_yzz_yyyzzzz, g_0_x_yzz_yyzzzzz, g_0_x_yzz_yzzzzzz, g_0_x_yzz_zzzzzzz, g_0_x_zz_xxxxxxx, g_0_x_zz_xxxxxxxy, g_0_x_zz_xxxxxxy, g_0_x_zz_xxxxxxyy, g_0_x_zz_xxxxxxyz, g_0_x_zz_xxxxxxz, g_0_x_zz_xxxxxyy, g_0_x_zz_xxxxxyyy, g_0_x_zz_xxxxxyyz, g_0_x_zz_xxxxxyz, g_0_x_zz_xxxxxyzz, g_0_x_zz_xxxxxzz, g_0_x_zz_xxxxyyy, g_0_x_zz_xxxxyyyy, g_0_x_zz_xxxxyyyz, g_0_x_zz_xxxxyyz, g_0_x_zz_xxxxyyzz, g_0_x_zz_xxxxyzz, g_0_x_zz_xxxxyzzz, g_0_x_zz_xxxxzzz, g_0_x_zz_xxxyyyy, g_0_x_zz_xxxyyyyy, g_0_x_zz_xxxyyyyz, g_0_x_zz_xxxyyyz, g_0_x_zz_xxxyyyzz, g_0_x_zz_xxxyyzz, g_0_x_zz_xxxyyzzz, g_0_x_zz_xxxyzzz, g_0_x_zz_xxxyzzzz, g_0_x_zz_xxxzzzz, g_0_x_zz_xxyyyyy, g_0_x_zz_xxyyyyyy, g_0_x_zz_xxyyyyyz, g_0_x_zz_xxyyyyz, g_0_x_zz_xxyyyyzz, g_0_x_zz_xxyyyzz, g_0_x_zz_xxyyyzzz, g_0_x_zz_xxyyzzz, g_0_x_zz_xxyyzzzz, g_0_x_zz_xxyzzzz, g_0_x_zz_xxyzzzzz, g_0_x_zz_xxzzzzz, g_0_x_zz_xyyyyyy, g_0_x_zz_xyyyyyyy, g_0_x_zz_xyyyyyyz, g_0_x_zz_xyyyyyz, g_0_x_zz_xyyyyyzz, g_0_x_zz_xyyyyzz, g_0_x_zz_xyyyyzzz, g_0_x_zz_xyyyzzz, g_0_x_zz_xyyyzzzz, g_0_x_zz_xyyzzzz, g_0_x_zz_xyyzzzzz, g_0_x_zz_xyzzzzz, g_0_x_zz_xyzzzzzz, g_0_x_zz_xzzzzzz, g_0_x_zz_yyyyyyy, g_0_x_zz_yyyyyyyy, g_0_x_zz_yyyyyyyz, g_0_x_zz_yyyyyyz, g_0_x_zz_yyyyyyzz, g_0_x_zz_yyyyyzz, g_0_x_zz_yyyyyzzz, g_0_x_zz_yyyyzzz, g_0_x_zz_yyyyzzzz, g_0_x_zz_yyyzzzz, g_0_x_zz_yyyzzzzz, g_0_x_zz_yyzzzzz, g_0_x_zz_yyzzzzzz, g_0_x_zz_yzzzzzz, g_0_x_zz_yzzzzzzz, g_0_x_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzz_xxxxxxx[k] = -g_0_x_zz_xxxxxxx[k] * ab_y + g_0_x_zz_xxxxxxxy[k];

                g_0_x_yzz_xxxxxxy[k] = -g_0_x_zz_xxxxxxy[k] * ab_y + g_0_x_zz_xxxxxxyy[k];

                g_0_x_yzz_xxxxxxz[k] = -g_0_x_zz_xxxxxxz[k] * ab_y + g_0_x_zz_xxxxxxyz[k];

                g_0_x_yzz_xxxxxyy[k] = -g_0_x_zz_xxxxxyy[k] * ab_y + g_0_x_zz_xxxxxyyy[k];

                g_0_x_yzz_xxxxxyz[k] = -g_0_x_zz_xxxxxyz[k] * ab_y + g_0_x_zz_xxxxxyyz[k];

                g_0_x_yzz_xxxxxzz[k] = -g_0_x_zz_xxxxxzz[k] * ab_y + g_0_x_zz_xxxxxyzz[k];

                g_0_x_yzz_xxxxyyy[k] = -g_0_x_zz_xxxxyyy[k] * ab_y + g_0_x_zz_xxxxyyyy[k];

                g_0_x_yzz_xxxxyyz[k] = -g_0_x_zz_xxxxyyz[k] * ab_y + g_0_x_zz_xxxxyyyz[k];

                g_0_x_yzz_xxxxyzz[k] = -g_0_x_zz_xxxxyzz[k] * ab_y + g_0_x_zz_xxxxyyzz[k];

                g_0_x_yzz_xxxxzzz[k] = -g_0_x_zz_xxxxzzz[k] * ab_y + g_0_x_zz_xxxxyzzz[k];

                g_0_x_yzz_xxxyyyy[k] = -g_0_x_zz_xxxyyyy[k] * ab_y + g_0_x_zz_xxxyyyyy[k];

                g_0_x_yzz_xxxyyyz[k] = -g_0_x_zz_xxxyyyz[k] * ab_y + g_0_x_zz_xxxyyyyz[k];

                g_0_x_yzz_xxxyyzz[k] = -g_0_x_zz_xxxyyzz[k] * ab_y + g_0_x_zz_xxxyyyzz[k];

                g_0_x_yzz_xxxyzzz[k] = -g_0_x_zz_xxxyzzz[k] * ab_y + g_0_x_zz_xxxyyzzz[k];

                g_0_x_yzz_xxxzzzz[k] = -g_0_x_zz_xxxzzzz[k] * ab_y + g_0_x_zz_xxxyzzzz[k];

                g_0_x_yzz_xxyyyyy[k] = -g_0_x_zz_xxyyyyy[k] * ab_y + g_0_x_zz_xxyyyyyy[k];

                g_0_x_yzz_xxyyyyz[k] = -g_0_x_zz_xxyyyyz[k] * ab_y + g_0_x_zz_xxyyyyyz[k];

                g_0_x_yzz_xxyyyzz[k] = -g_0_x_zz_xxyyyzz[k] * ab_y + g_0_x_zz_xxyyyyzz[k];

                g_0_x_yzz_xxyyzzz[k] = -g_0_x_zz_xxyyzzz[k] * ab_y + g_0_x_zz_xxyyyzzz[k];

                g_0_x_yzz_xxyzzzz[k] = -g_0_x_zz_xxyzzzz[k] * ab_y + g_0_x_zz_xxyyzzzz[k];

                g_0_x_yzz_xxzzzzz[k] = -g_0_x_zz_xxzzzzz[k] * ab_y + g_0_x_zz_xxyzzzzz[k];

                g_0_x_yzz_xyyyyyy[k] = -g_0_x_zz_xyyyyyy[k] * ab_y + g_0_x_zz_xyyyyyyy[k];

                g_0_x_yzz_xyyyyyz[k] = -g_0_x_zz_xyyyyyz[k] * ab_y + g_0_x_zz_xyyyyyyz[k];

                g_0_x_yzz_xyyyyzz[k] = -g_0_x_zz_xyyyyzz[k] * ab_y + g_0_x_zz_xyyyyyzz[k];

                g_0_x_yzz_xyyyzzz[k] = -g_0_x_zz_xyyyzzz[k] * ab_y + g_0_x_zz_xyyyyzzz[k];

                g_0_x_yzz_xyyzzzz[k] = -g_0_x_zz_xyyzzzz[k] * ab_y + g_0_x_zz_xyyyzzzz[k];

                g_0_x_yzz_xyzzzzz[k] = -g_0_x_zz_xyzzzzz[k] * ab_y + g_0_x_zz_xyyzzzzz[k];

                g_0_x_yzz_xzzzzzz[k] = -g_0_x_zz_xzzzzzz[k] * ab_y + g_0_x_zz_xyzzzzzz[k];

                g_0_x_yzz_yyyyyyy[k] = -g_0_x_zz_yyyyyyy[k] * ab_y + g_0_x_zz_yyyyyyyy[k];

                g_0_x_yzz_yyyyyyz[k] = -g_0_x_zz_yyyyyyz[k] * ab_y + g_0_x_zz_yyyyyyyz[k];

                g_0_x_yzz_yyyyyzz[k] = -g_0_x_zz_yyyyyzz[k] * ab_y + g_0_x_zz_yyyyyyzz[k];

                g_0_x_yzz_yyyyzzz[k] = -g_0_x_zz_yyyyzzz[k] * ab_y + g_0_x_zz_yyyyyzzz[k];

                g_0_x_yzz_yyyzzzz[k] = -g_0_x_zz_yyyzzzz[k] * ab_y + g_0_x_zz_yyyyzzzz[k];

                g_0_x_yzz_yyzzzzz[k] = -g_0_x_zz_yyzzzzz[k] * ab_y + g_0_x_zz_yyyzzzzz[k];

                g_0_x_yzz_yzzzzzz[k] = -g_0_x_zz_yzzzzzz[k] * ab_y + g_0_x_zz_yyzzzzzz[k];

                g_0_x_yzz_zzzzzzz[k] = -g_0_x_zz_zzzzzzz[k] * ab_y + g_0_x_zz_yzzzzzzz[k];
            }

            /// Set up 324-360 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_zz_xxxxxxx, g_0_x_zz_xxxxxxxz, g_0_x_zz_xxxxxxy, g_0_x_zz_xxxxxxyz, g_0_x_zz_xxxxxxz, g_0_x_zz_xxxxxxzz, g_0_x_zz_xxxxxyy, g_0_x_zz_xxxxxyyz, g_0_x_zz_xxxxxyz, g_0_x_zz_xxxxxyzz, g_0_x_zz_xxxxxzz, g_0_x_zz_xxxxxzzz, g_0_x_zz_xxxxyyy, g_0_x_zz_xxxxyyyz, g_0_x_zz_xxxxyyz, g_0_x_zz_xxxxyyzz, g_0_x_zz_xxxxyzz, g_0_x_zz_xxxxyzzz, g_0_x_zz_xxxxzzz, g_0_x_zz_xxxxzzzz, g_0_x_zz_xxxyyyy, g_0_x_zz_xxxyyyyz, g_0_x_zz_xxxyyyz, g_0_x_zz_xxxyyyzz, g_0_x_zz_xxxyyzz, g_0_x_zz_xxxyyzzz, g_0_x_zz_xxxyzzz, g_0_x_zz_xxxyzzzz, g_0_x_zz_xxxzzzz, g_0_x_zz_xxxzzzzz, g_0_x_zz_xxyyyyy, g_0_x_zz_xxyyyyyz, g_0_x_zz_xxyyyyz, g_0_x_zz_xxyyyyzz, g_0_x_zz_xxyyyzz, g_0_x_zz_xxyyyzzz, g_0_x_zz_xxyyzzz, g_0_x_zz_xxyyzzzz, g_0_x_zz_xxyzzzz, g_0_x_zz_xxyzzzzz, g_0_x_zz_xxzzzzz, g_0_x_zz_xxzzzzzz, g_0_x_zz_xyyyyyy, g_0_x_zz_xyyyyyyz, g_0_x_zz_xyyyyyz, g_0_x_zz_xyyyyyzz, g_0_x_zz_xyyyyzz, g_0_x_zz_xyyyyzzz, g_0_x_zz_xyyyzzz, g_0_x_zz_xyyyzzzz, g_0_x_zz_xyyzzzz, g_0_x_zz_xyyzzzzz, g_0_x_zz_xyzzzzz, g_0_x_zz_xyzzzzzz, g_0_x_zz_xzzzzzz, g_0_x_zz_xzzzzzzz, g_0_x_zz_yyyyyyy, g_0_x_zz_yyyyyyyz, g_0_x_zz_yyyyyyz, g_0_x_zz_yyyyyyzz, g_0_x_zz_yyyyyzz, g_0_x_zz_yyyyyzzz, g_0_x_zz_yyyyzzz, g_0_x_zz_yyyyzzzz, g_0_x_zz_yyyzzzz, g_0_x_zz_yyyzzzzz, g_0_x_zz_yyzzzzz, g_0_x_zz_yyzzzzzz, g_0_x_zz_yzzzzzz, g_0_x_zz_yzzzzzzz, g_0_x_zz_zzzzzzz, g_0_x_zz_zzzzzzzz, g_0_x_zzz_xxxxxxx, g_0_x_zzz_xxxxxxy, g_0_x_zzz_xxxxxxz, g_0_x_zzz_xxxxxyy, g_0_x_zzz_xxxxxyz, g_0_x_zzz_xxxxxzz, g_0_x_zzz_xxxxyyy, g_0_x_zzz_xxxxyyz, g_0_x_zzz_xxxxyzz, g_0_x_zzz_xxxxzzz, g_0_x_zzz_xxxyyyy, g_0_x_zzz_xxxyyyz, g_0_x_zzz_xxxyyzz, g_0_x_zzz_xxxyzzz, g_0_x_zzz_xxxzzzz, g_0_x_zzz_xxyyyyy, g_0_x_zzz_xxyyyyz, g_0_x_zzz_xxyyyzz, g_0_x_zzz_xxyyzzz, g_0_x_zzz_xxyzzzz, g_0_x_zzz_xxzzzzz, g_0_x_zzz_xyyyyyy, g_0_x_zzz_xyyyyyz, g_0_x_zzz_xyyyyzz, g_0_x_zzz_xyyyzzz, g_0_x_zzz_xyyzzzz, g_0_x_zzz_xyzzzzz, g_0_x_zzz_xzzzzzz, g_0_x_zzz_yyyyyyy, g_0_x_zzz_yyyyyyz, g_0_x_zzz_yyyyyzz, g_0_x_zzz_yyyyzzz, g_0_x_zzz_yyyzzzz, g_0_x_zzz_yyzzzzz, g_0_x_zzz_yzzzzzz, g_0_x_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzz_xxxxxxx[k] = -g_0_x_zz_xxxxxxx[k] * ab_z + g_0_x_zz_xxxxxxxz[k];

                g_0_x_zzz_xxxxxxy[k] = -g_0_x_zz_xxxxxxy[k] * ab_z + g_0_x_zz_xxxxxxyz[k];

                g_0_x_zzz_xxxxxxz[k] = -g_0_x_zz_xxxxxxz[k] * ab_z + g_0_x_zz_xxxxxxzz[k];

                g_0_x_zzz_xxxxxyy[k] = -g_0_x_zz_xxxxxyy[k] * ab_z + g_0_x_zz_xxxxxyyz[k];

                g_0_x_zzz_xxxxxyz[k] = -g_0_x_zz_xxxxxyz[k] * ab_z + g_0_x_zz_xxxxxyzz[k];

                g_0_x_zzz_xxxxxzz[k] = -g_0_x_zz_xxxxxzz[k] * ab_z + g_0_x_zz_xxxxxzzz[k];

                g_0_x_zzz_xxxxyyy[k] = -g_0_x_zz_xxxxyyy[k] * ab_z + g_0_x_zz_xxxxyyyz[k];

                g_0_x_zzz_xxxxyyz[k] = -g_0_x_zz_xxxxyyz[k] * ab_z + g_0_x_zz_xxxxyyzz[k];

                g_0_x_zzz_xxxxyzz[k] = -g_0_x_zz_xxxxyzz[k] * ab_z + g_0_x_zz_xxxxyzzz[k];

                g_0_x_zzz_xxxxzzz[k] = -g_0_x_zz_xxxxzzz[k] * ab_z + g_0_x_zz_xxxxzzzz[k];

                g_0_x_zzz_xxxyyyy[k] = -g_0_x_zz_xxxyyyy[k] * ab_z + g_0_x_zz_xxxyyyyz[k];

                g_0_x_zzz_xxxyyyz[k] = -g_0_x_zz_xxxyyyz[k] * ab_z + g_0_x_zz_xxxyyyzz[k];

                g_0_x_zzz_xxxyyzz[k] = -g_0_x_zz_xxxyyzz[k] * ab_z + g_0_x_zz_xxxyyzzz[k];

                g_0_x_zzz_xxxyzzz[k] = -g_0_x_zz_xxxyzzz[k] * ab_z + g_0_x_zz_xxxyzzzz[k];

                g_0_x_zzz_xxxzzzz[k] = -g_0_x_zz_xxxzzzz[k] * ab_z + g_0_x_zz_xxxzzzzz[k];

                g_0_x_zzz_xxyyyyy[k] = -g_0_x_zz_xxyyyyy[k] * ab_z + g_0_x_zz_xxyyyyyz[k];

                g_0_x_zzz_xxyyyyz[k] = -g_0_x_zz_xxyyyyz[k] * ab_z + g_0_x_zz_xxyyyyzz[k];

                g_0_x_zzz_xxyyyzz[k] = -g_0_x_zz_xxyyyzz[k] * ab_z + g_0_x_zz_xxyyyzzz[k];

                g_0_x_zzz_xxyyzzz[k] = -g_0_x_zz_xxyyzzz[k] * ab_z + g_0_x_zz_xxyyzzzz[k];

                g_0_x_zzz_xxyzzzz[k] = -g_0_x_zz_xxyzzzz[k] * ab_z + g_0_x_zz_xxyzzzzz[k];

                g_0_x_zzz_xxzzzzz[k] = -g_0_x_zz_xxzzzzz[k] * ab_z + g_0_x_zz_xxzzzzzz[k];

                g_0_x_zzz_xyyyyyy[k] = -g_0_x_zz_xyyyyyy[k] * ab_z + g_0_x_zz_xyyyyyyz[k];

                g_0_x_zzz_xyyyyyz[k] = -g_0_x_zz_xyyyyyz[k] * ab_z + g_0_x_zz_xyyyyyzz[k];

                g_0_x_zzz_xyyyyzz[k] = -g_0_x_zz_xyyyyzz[k] * ab_z + g_0_x_zz_xyyyyzzz[k];

                g_0_x_zzz_xyyyzzz[k] = -g_0_x_zz_xyyyzzz[k] * ab_z + g_0_x_zz_xyyyzzzz[k];

                g_0_x_zzz_xyyzzzz[k] = -g_0_x_zz_xyyzzzz[k] * ab_z + g_0_x_zz_xyyzzzzz[k];

                g_0_x_zzz_xyzzzzz[k] = -g_0_x_zz_xyzzzzz[k] * ab_z + g_0_x_zz_xyzzzzzz[k];

                g_0_x_zzz_xzzzzzz[k] = -g_0_x_zz_xzzzzzz[k] * ab_z + g_0_x_zz_xzzzzzzz[k];

                g_0_x_zzz_yyyyyyy[k] = -g_0_x_zz_yyyyyyy[k] * ab_z + g_0_x_zz_yyyyyyyz[k];

                g_0_x_zzz_yyyyyyz[k] = -g_0_x_zz_yyyyyyz[k] * ab_z + g_0_x_zz_yyyyyyzz[k];

                g_0_x_zzz_yyyyyzz[k] = -g_0_x_zz_yyyyyzz[k] * ab_z + g_0_x_zz_yyyyyzzz[k];

                g_0_x_zzz_yyyyzzz[k] = -g_0_x_zz_yyyyzzz[k] * ab_z + g_0_x_zz_yyyyzzzz[k];

                g_0_x_zzz_yyyzzzz[k] = -g_0_x_zz_yyyzzzz[k] * ab_z + g_0_x_zz_yyyzzzzz[k];

                g_0_x_zzz_yyzzzzz[k] = -g_0_x_zz_yyzzzzz[k] * ab_z + g_0_x_zz_yyzzzzzz[k];

                g_0_x_zzz_yzzzzzz[k] = -g_0_x_zz_yzzzzzz[k] * ab_z + g_0_x_zz_yzzzzzzz[k];

                g_0_x_zzz_zzzzzzz[k] = -g_0_x_zz_zzzzzzz[k] * ab_z + g_0_x_zz_zzzzzzzz[k];
            }

            /// Set up 360-396 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xx_xxxxxxx, g_0_y_xx_xxxxxxxx, g_0_y_xx_xxxxxxxy, g_0_y_xx_xxxxxxxz, g_0_y_xx_xxxxxxy, g_0_y_xx_xxxxxxyy, g_0_y_xx_xxxxxxyz, g_0_y_xx_xxxxxxz, g_0_y_xx_xxxxxxzz, g_0_y_xx_xxxxxyy, g_0_y_xx_xxxxxyyy, g_0_y_xx_xxxxxyyz, g_0_y_xx_xxxxxyz, g_0_y_xx_xxxxxyzz, g_0_y_xx_xxxxxzz, g_0_y_xx_xxxxxzzz, g_0_y_xx_xxxxyyy, g_0_y_xx_xxxxyyyy, g_0_y_xx_xxxxyyyz, g_0_y_xx_xxxxyyz, g_0_y_xx_xxxxyyzz, g_0_y_xx_xxxxyzz, g_0_y_xx_xxxxyzzz, g_0_y_xx_xxxxzzz, g_0_y_xx_xxxxzzzz, g_0_y_xx_xxxyyyy, g_0_y_xx_xxxyyyyy, g_0_y_xx_xxxyyyyz, g_0_y_xx_xxxyyyz, g_0_y_xx_xxxyyyzz, g_0_y_xx_xxxyyzz, g_0_y_xx_xxxyyzzz, g_0_y_xx_xxxyzzz, g_0_y_xx_xxxyzzzz, g_0_y_xx_xxxzzzz, g_0_y_xx_xxxzzzzz, g_0_y_xx_xxyyyyy, g_0_y_xx_xxyyyyyy, g_0_y_xx_xxyyyyyz, g_0_y_xx_xxyyyyz, g_0_y_xx_xxyyyyzz, g_0_y_xx_xxyyyzz, g_0_y_xx_xxyyyzzz, g_0_y_xx_xxyyzzz, g_0_y_xx_xxyyzzzz, g_0_y_xx_xxyzzzz, g_0_y_xx_xxyzzzzz, g_0_y_xx_xxzzzzz, g_0_y_xx_xxzzzzzz, g_0_y_xx_xyyyyyy, g_0_y_xx_xyyyyyyy, g_0_y_xx_xyyyyyyz, g_0_y_xx_xyyyyyz, g_0_y_xx_xyyyyyzz, g_0_y_xx_xyyyyzz, g_0_y_xx_xyyyyzzz, g_0_y_xx_xyyyzzz, g_0_y_xx_xyyyzzzz, g_0_y_xx_xyyzzzz, g_0_y_xx_xyyzzzzz, g_0_y_xx_xyzzzzz, g_0_y_xx_xyzzzzzz, g_0_y_xx_xzzzzzz, g_0_y_xx_xzzzzzzz, g_0_y_xx_yyyyyyy, g_0_y_xx_yyyyyyz, g_0_y_xx_yyyyyzz, g_0_y_xx_yyyyzzz, g_0_y_xx_yyyzzzz, g_0_y_xx_yyzzzzz, g_0_y_xx_yzzzzzz, g_0_y_xx_zzzzzzz, g_0_y_xxx_xxxxxxx, g_0_y_xxx_xxxxxxy, g_0_y_xxx_xxxxxxz, g_0_y_xxx_xxxxxyy, g_0_y_xxx_xxxxxyz, g_0_y_xxx_xxxxxzz, g_0_y_xxx_xxxxyyy, g_0_y_xxx_xxxxyyz, g_0_y_xxx_xxxxyzz, g_0_y_xxx_xxxxzzz, g_0_y_xxx_xxxyyyy, g_0_y_xxx_xxxyyyz, g_0_y_xxx_xxxyyzz, g_0_y_xxx_xxxyzzz, g_0_y_xxx_xxxzzzz, g_0_y_xxx_xxyyyyy, g_0_y_xxx_xxyyyyz, g_0_y_xxx_xxyyyzz, g_0_y_xxx_xxyyzzz, g_0_y_xxx_xxyzzzz, g_0_y_xxx_xxzzzzz, g_0_y_xxx_xyyyyyy, g_0_y_xxx_xyyyyyz, g_0_y_xxx_xyyyyzz, g_0_y_xxx_xyyyzzz, g_0_y_xxx_xyyzzzz, g_0_y_xxx_xyzzzzz, g_0_y_xxx_xzzzzzz, g_0_y_xxx_yyyyyyy, g_0_y_xxx_yyyyyyz, g_0_y_xxx_yyyyyzz, g_0_y_xxx_yyyyzzz, g_0_y_xxx_yyyzzzz, g_0_y_xxx_yyzzzzz, g_0_y_xxx_yzzzzzz, g_0_y_xxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxx_xxxxxxx[k] = -g_0_y_xx_xxxxxxx[k] * ab_x + g_0_y_xx_xxxxxxxx[k];

                g_0_y_xxx_xxxxxxy[k] = -g_0_y_xx_xxxxxxy[k] * ab_x + g_0_y_xx_xxxxxxxy[k];

                g_0_y_xxx_xxxxxxz[k] = -g_0_y_xx_xxxxxxz[k] * ab_x + g_0_y_xx_xxxxxxxz[k];

                g_0_y_xxx_xxxxxyy[k] = -g_0_y_xx_xxxxxyy[k] * ab_x + g_0_y_xx_xxxxxxyy[k];

                g_0_y_xxx_xxxxxyz[k] = -g_0_y_xx_xxxxxyz[k] * ab_x + g_0_y_xx_xxxxxxyz[k];

                g_0_y_xxx_xxxxxzz[k] = -g_0_y_xx_xxxxxzz[k] * ab_x + g_0_y_xx_xxxxxxzz[k];

                g_0_y_xxx_xxxxyyy[k] = -g_0_y_xx_xxxxyyy[k] * ab_x + g_0_y_xx_xxxxxyyy[k];

                g_0_y_xxx_xxxxyyz[k] = -g_0_y_xx_xxxxyyz[k] * ab_x + g_0_y_xx_xxxxxyyz[k];

                g_0_y_xxx_xxxxyzz[k] = -g_0_y_xx_xxxxyzz[k] * ab_x + g_0_y_xx_xxxxxyzz[k];

                g_0_y_xxx_xxxxzzz[k] = -g_0_y_xx_xxxxzzz[k] * ab_x + g_0_y_xx_xxxxxzzz[k];

                g_0_y_xxx_xxxyyyy[k] = -g_0_y_xx_xxxyyyy[k] * ab_x + g_0_y_xx_xxxxyyyy[k];

                g_0_y_xxx_xxxyyyz[k] = -g_0_y_xx_xxxyyyz[k] * ab_x + g_0_y_xx_xxxxyyyz[k];

                g_0_y_xxx_xxxyyzz[k] = -g_0_y_xx_xxxyyzz[k] * ab_x + g_0_y_xx_xxxxyyzz[k];

                g_0_y_xxx_xxxyzzz[k] = -g_0_y_xx_xxxyzzz[k] * ab_x + g_0_y_xx_xxxxyzzz[k];

                g_0_y_xxx_xxxzzzz[k] = -g_0_y_xx_xxxzzzz[k] * ab_x + g_0_y_xx_xxxxzzzz[k];

                g_0_y_xxx_xxyyyyy[k] = -g_0_y_xx_xxyyyyy[k] * ab_x + g_0_y_xx_xxxyyyyy[k];

                g_0_y_xxx_xxyyyyz[k] = -g_0_y_xx_xxyyyyz[k] * ab_x + g_0_y_xx_xxxyyyyz[k];

                g_0_y_xxx_xxyyyzz[k] = -g_0_y_xx_xxyyyzz[k] * ab_x + g_0_y_xx_xxxyyyzz[k];

                g_0_y_xxx_xxyyzzz[k] = -g_0_y_xx_xxyyzzz[k] * ab_x + g_0_y_xx_xxxyyzzz[k];

                g_0_y_xxx_xxyzzzz[k] = -g_0_y_xx_xxyzzzz[k] * ab_x + g_0_y_xx_xxxyzzzz[k];

                g_0_y_xxx_xxzzzzz[k] = -g_0_y_xx_xxzzzzz[k] * ab_x + g_0_y_xx_xxxzzzzz[k];

                g_0_y_xxx_xyyyyyy[k] = -g_0_y_xx_xyyyyyy[k] * ab_x + g_0_y_xx_xxyyyyyy[k];

                g_0_y_xxx_xyyyyyz[k] = -g_0_y_xx_xyyyyyz[k] * ab_x + g_0_y_xx_xxyyyyyz[k];

                g_0_y_xxx_xyyyyzz[k] = -g_0_y_xx_xyyyyzz[k] * ab_x + g_0_y_xx_xxyyyyzz[k];

                g_0_y_xxx_xyyyzzz[k] = -g_0_y_xx_xyyyzzz[k] * ab_x + g_0_y_xx_xxyyyzzz[k];

                g_0_y_xxx_xyyzzzz[k] = -g_0_y_xx_xyyzzzz[k] * ab_x + g_0_y_xx_xxyyzzzz[k];

                g_0_y_xxx_xyzzzzz[k] = -g_0_y_xx_xyzzzzz[k] * ab_x + g_0_y_xx_xxyzzzzz[k];

                g_0_y_xxx_xzzzzzz[k] = -g_0_y_xx_xzzzzzz[k] * ab_x + g_0_y_xx_xxzzzzzz[k];

                g_0_y_xxx_yyyyyyy[k] = -g_0_y_xx_yyyyyyy[k] * ab_x + g_0_y_xx_xyyyyyyy[k];

                g_0_y_xxx_yyyyyyz[k] = -g_0_y_xx_yyyyyyz[k] * ab_x + g_0_y_xx_xyyyyyyz[k];

                g_0_y_xxx_yyyyyzz[k] = -g_0_y_xx_yyyyyzz[k] * ab_x + g_0_y_xx_xyyyyyzz[k];

                g_0_y_xxx_yyyyzzz[k] = -g_0_y_xx_yyyyzzz[k] * ab_x + g_0_y_xx_xyyyyzzz[k];

                g_0_y_xxx_yyyzzzz[k] = -g_0_y_xx_yyyzzzz[k] * ab_x + g_0_y_xx_xyyyzzzz[k];

                g_0_y_xxx_yyzzzzz[k] = -g_0_y_xx_yyzzzzz[k] * ab_x + g_0_y_xx_xyyzzzzz[k];

                g_0_y_xxx_yzzzzzz[k] = -g_0_y_xx_yzzzzzz[k] * ab_x + g_0_y_xx_xyzzzzzz[k];

                g_0_y_xxx_zzzzzzz[k] = -g_0_y_xx_zzzzzzz[k] * ab_x + g_0_y_xx_xzzzzzzz[k];
            }

            /// Set up 396-432 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxy_xxxxxxx, g_0_y_xxy_xxxxxxy, g_0_y_xxy_xxxxxxz, g_0_y_xxy_xxxxxyy, g_0_y_xxy_xxxxxyz, g_0_y_xxy_xxxxxzz, g_0_y_xxy_xxxxyyy, g_0_y_xxy_xxxxyyz, g_0_y_xxy_xxxxyzz, g_0_y_xxy_xxxxzzz, g_0_y_xxy_xxxyyyy, g_0_y_xxy_xxxyyyz, g_0_y_xxy_xxxyyzz, g_0_y_xxy_xxxyzzz, g_0_y_xxy_xxxzzzz, g_0_y_xxy_xxyyyyy, g_0_y_xxy_xxyyyyz, g_0_y_xxy_xxyyyzz, g_0_y_xxy_xxyyzzz, g_0_y_xxy_xxyzzzz, g_0_y_xxy_xxzzzzz, g_0_y_xxy_xyyyyyy, g_0_y_xxy_xyyyyyz, g_0_y_xxy_xyyyyzz, g_0_y_xxy_xyyyzzz, g_0_y_xxy_xyyzzzz, g_0_y_xxy_xyzzzzz, g_0_y_xxy_xzzzzzz, g_0_y_xxy_yyyyyyy, g_0_y_xxy_yyyyyyz, g_0_y_xxy_yyyyyzz, g_0_y_xxy_yyyyzzz, g_0_y_xxy_yyyzzzz, g_0_y_xxy_yyzzzzz, g_0_y_xxy_yzzzzzz, g_0_y_xxy_zzzzzzz, g_0_y_xy_xxxxxxx, g_0_y_xy_xxxxxxxx, g_0_y_xy_xxxxxxxy, g_0_y_xy_xxxxxxxz, g_0_y_xy_xxxxxxy, g_0_y_xy_xxxxxxyy, g_0_y_xy_xxxxxxyz, g_0_y_xy_xxxxxxz, g_0_y_xy_xxxxxxzz, g_0_y_xy_xxxxxyy, g_0_y_xy_xxxxxyyy, g_0_y_xy_xxxxxyyz, g_0_y_xy_xxxxxyz, g_0_y_xy_xxxxxyzz, g_0_y_xy_xxxxxzz, g_0_y_xy_xxxxxzzz, g_0_y_xy_xxxxyyy, g_0_y_xy_xxxxyyyy, g_0_y_xy_xxxxyyyz, g_0_y_xy_xxxxyyz, g_0_y_xy_xxxxyyzz, g_0_y_xy_xxxxyzz, g_0_y_xy_xxxxyzzz, g_0_y_xy_xxxxzzz, g_0_y_xy_xxxxzzzz, g_0_y_xy_xxxyyyy, g_0_y_xy_xxxyyyyy, g_0_y_xy_xxxyyyyz, g_0_y_xy_xxxyyyz, g_0_y_xy_xxxyyyzz, g_0_y_xy_xxxyyzz, g_0_y_xy_xxxyyzzz, g_0_y_xy_xxxyzzz, g_0_y_xy_xxxyzzzz, g_0_y_xy_xxxzzzz, g_0_y_xy_xxxzzzzz, g_0_y_xy_xxyyyyy, g_0_y_xy_xxyyyyyy, g_0_y_xy_xxyyyyyz, g_0_y_xy_xxyyyyz, g_0_y_xy_xxyyyyzz, g_0_y_xy_xxyyyzz, g_0_y_xy_xxyyyzzz, g_0_y_xy_xxyyzzz, g_0_y_xy_xxyyzzzz, g_0_y_xy_xxyzzzz, g_0_y_xy_xxyzzzzz, g_0_y_xy_xxzzzzz, g_0_y_xy_xxzzzzzz, g_0_y_xy_xyyyyyy, g_0_y_xy_xyyyyyyy, g_0_y_xy_xyyyyyyz, g_0_y_xy_xyyyyyz, g_0_y_xy_xyyyyyzz, g_0_y_xy_xyyyyzz, g_0_y_xy_xyyyyzzz, g_0_y_xy_xyyyzzz, g_0_y_xy_xyyyzzzz, g_0_y_xy_xyyzzzz, g_0_y_xy_xyyzzzzz, g_0_y_xy_xyzzzzz, g_0_y_xy_xyzzzzzz, g_0_y_xy_xzzzzzz, g_0_y_xy_xzzzzzzz, g_0_y_xy_yyyyyyy, g_0_y_xy_yyyyyyz, g_0_y_xy_yyyyyzz, g_0_y_xy_yyyyzzz, g_0_y_xy_yyyzzzz, g_0_y_xy_yyzzzzz, g_0_y_xy_yzzzzzz, g_0_y_xy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxy_xxxxxxx[k] = -g_0_y_xy_xxxxxxx[k] * ab_x + g_0_y_xy_xxxxxxxx[k];

                g_0_y_xxy_xxxxxxy[k] = -g_0_y_xy_xxxxxxy[k] * ab_x + g_0_y_xy_xxxxxxxy[k];

                g_0_y_xxy_xxxxxxz[k] = -g_0_y_xy_xxxxxxz[k] * ab_x + g_0_y_xy_xxxxxxxz[k];

                g_0_y_xxy_xxxxxyy[k] = -g_0_y_xy_xxxxxyy[k] * ab_x + g_0_y_xy_xxxxxxyy[k];

                g_0_y_xxy_xxxxxyz[k] = -g_0_y_xy_xxxxxyz[k] * ab_x + g_0_y_xy_xxxxxxyz[k];

                g_0_y_xxy_xxxxxzz[k] = -g_0_y_xy_xxxxxzz[k] * ab_x + g_0_y_xy_xxxxxxzz[k];

                g_0_y_xxy_xxxxyyy[k] = -g_0_y_xy_xxxxyyy[k] * ab_x + g_0_y_xy_xxxxxyyy[k];

                g_0_y_xxy_xxxxyyz[k] = -g_0_y_xy_xxxxyyz[k] * ab_x + g_0_y_xy_xxxxxyyz[k];

                g_0_y_xxy_xxxxyzz[k] = -g_0_y_xy_xxxxyzz[k] * ab_x + g_0_y_xy_xxxxxyzz[k];

                g_0_y_xxy_xxxxzzz[k] = -g_0_y_xy_xxxxzzz[k] * ab_x + g_0_y_xy_xxxxxzzz[k];

                g_0_y_xxy_xxxyyyy[k] = -g_0_y_xy_xxxyyyy[k] * ab_x + g_0_y_xy_xxxxyyyy[k];

                g_0_y_xxy_xxxyyyz[k] = -g_0_y_xy_xxxyyyz[k] * ab_x + g_0_y_xy_xxxxyyyz[k];

                g_0_y_xxy_xxxyyzz[k] = -g_0_y_xy_xxxyyzz[k] * ab_x + g_0_y_xy_xxxxyyzz[k];

                g_0_y_xxy_xxxyzzz[k] = -g_0_y_xy_xxxyzzz[k] * ab_x + g_0_y_xy_xxxxyzzz[k];

                g_0_y_xxy_xxxzzzz[k] = -g_0_y_xy_xxxzzzz[k] * ab_x + g_0_y_xy_xxxxzzzz[k];

                g_0_y_xxy_xxyyyyy[k] = -g_0_y_xy_xxyyyyy[k] * ab_x + g_0_y_xy_xxxyyyyy[k];

                g_0_y_xxy_xxyyyyz[k] = -g_0_y_xy_xxyyyyz[k] * ab_x + g_0_y_xy_xxxyyyyz[k];

                g_0_y_xxy_xxyyyzz[k] = -g_0_y_xy_xxyyyzz[k] * ab_x + g_0_y_xy_xxxyyyzz[k];

                g_0_y_xxy_xxyyzzz[k] = -g_0_y_xy_xxyyzzz[k] * ab_x + g_0_y_xy_xxxyyzzz[k];

                g_0_y_xxy_xxyzzzz[k] = -g_0_y_xy_xxyzzzz[k] * ab_x + g_0_y_xy_xxxyzzzz[k];

                g_0_y_xxy_xxzzzzz[k] = -g_0_y_xy_xxzzzzz[k] * ab_x + g_0_y_xy_xxxzzzzz[k];

                g_0_y_xxy_xyyyyyy[k] = -g_0_y_xy_xyyyyyy[k] * ab_x + g_0_y_xy_xxyyyyyy[k];

                g_0_y_xxy_xyyyyyz[k] = -g_0_y_xy_xyyyyyz[k] * ab_x + g_0_y_xy_xxyyyyyz[k];

                g_0_y_xxy_xyyyyzz[k] = -g_0_y_xy_xyyyyzz[k] * ab_x + g_0_y_xy_xxyyyyzz[k];

                g_0_y_xxy_xyyyzzz[k] = -g_0_y_xy_xyyyzzz[k] * ab_x + g_0_y_xy_xxyyyzzz[k];

                g_0_y_xxy_xyyzzzz[k] = -g_0_y_xy_xyyzzzz[k] * ab_x + g_0_y_xy_xxyyzzzz[k];

                g_0_y_xxy_xyzzzzz[k] = -g_0_y_xy_xyzzzzz[k] * ab_x + g_0_y_xy_xxyzzzzz[k];

                g_0_y_xxy_xzzzzzz[k] = -g_0_y_xy_xzzzzzz[k] * ab_x + g_0_y_xy_xxzzzzzz[k];

                g_0_y_xxy_yyyyyyy[k] = -g_0_y_xy_yyyyyyy[k] * ab_x + g_0_y_xy_xyyyyyyy[k];

                g_0_y_xxy_yyyyyyz[k] = -g_0_y_xy_yyyyyyz[k] * ab_x + g_0_y_xy_xyyyyyyz[k];

                g_0_y_xxy_yyyyyzz[k] = -g_0_y_xy_yyyyyzz[k] * ab_x + g_0_y_xy_xyyyyyzz[k];

                g_0_y_xxy_yyyyzzz[k] = -g_0_y_xy_yyyyzzz[k] * ab_x + g_0_y_xy_xyyyyzzz[k];

                g_0_y_xxy_yyyzzzz[k] = -g_0_y_xy_yyyzzzz[k] * ab_x + g_0_y_xy_xyyyzzzz[k];

                g_0_y_xxy_yyzzzzz[k] = -g_0_y_xy_yyzzzzz[k] * ab_x + g_0_y_xy_xyyzzzzz[k];

                g_0_y_xxy_yzzzzzz[k] = -g_0_y_xy_yzzzzzz[k] * ab_x + g_0_y_xy_xyzzzzzz[k];

                g_0_y_xxy_zzzzzzz[k] = -g_0_y_xy_zzzzzzz[k] * ab_x + g_0_y_xy_xzzzzzzz[k];
            }

            /// Set up 432-468 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxz_xxxxxxx, g_0_y_xxz_xxxxxxy, g_0_y_xxz_xxxxxxz, g_0_y_xxz_xxxxxyy, g_0_y_xxz_xxxxxyz, g_0_y_xxz_xxxxxzz, g_0_y_xxz_xxxxyyy, g_0_y_xxz_xxxxyyz, g_0_y_xxz_xxxxyzz, g_0_y_xxz_xxxxzzz, g_0_y_xxz_xxxyyyy, g_0_y_xxz_xxxyyyz, g_0_y_xxz_xxxyyzz, g_0_y_xxz_xxxyzzz, g_0_y_xxz_xxxzzzz, g_0_y_xxz_xxyyyyy, g_0_y_xxz_xxyyyyz, g_0_y_xxz_xxyyyzz, g_0_y_xxz_xxyyzzz, g_0_y_xxz_xxyzzzz, g_0_y_xxz_xxzzzzz, g_0_y_xxz_xyyyyyy, g_0_y_xxz_xyyyyyz, g_0_y_xxz_xyyyyzz, g_0_y_xxz_xyyyzzz, g_0_y_xxz_xyyzzzz, g_0_y_xxz_xyzzzzz, g_0_y_xxz_xzzzzzz, g_0_y_xxz_yyyyyyy, g_0_y_xxz_yyyyyyz, g_0_y_xxz_yyyyyzz, g_0_y_xxz_yyyyzzz, g_0_y_xxz_yyyzzzz, g_0_y_xxz_yyzzzzz, g_0_y_xxz_yzzzzzz, g_0_y_xxz_zzzzzzz, g_0_y_xz_xxxxxxx, g_0_y_xz_xxxxxxxx, g_0_y_xz_xxxxxxxy, g_0_y_xz_xxxxxxxz, g_0_y_xz_xxxxxxy, g_0_y_xz_xxxxxxyy, g_0_y_xz_xxxxxxyz, g_0_y_xz_xxxxxxz, g_0_y_xz_xxxxxxzz, g_0_y_xz_xxxxxyy, g_0_y_xz_xxxxxyyy, g_0_y_xz_xxxxxyyz, g_0_y_xz_xxxxxyz, g_0_y_xz_xxxxxyzz, g_0_y_xz_xxxxxzz, g_0_y_xz_xxxxxzzz, g_0_y_xz_xxxxyyy, g_0_y_xz_xxxxyyyy, g_0_y_xz_xxxxyyyz, g_0_y_xz_xxxxyyz, g_0_y_xz_xxxxyyzz, g_0_y_xz_xxxxyzz, g_0_y_xz_xxxxyzzz, g_0_y_xz_xxxxzzz, g_0_y_xz_xxxxzzzz, g_0_y_xz_xxxyyyy, g_0_y_xz_xxxyyyyy, g_0_y_xz_xxxyyyyz, g_0_y_xz_xxxyyyz, g_0_y_xz_xxxyyyzz, g_0_y_xz_xxxyyzz, g_0_y_xz_xxxyyzzz, g_0_y_xz_xxxyzzz, g_0_y_xz_xxxyzzzz, g_0_y_xz_xxxzzzz, g_0_y_xz_xxxzzzzz, g_0_y_xz_xxyyyyy, g_0_y_xz_xxyyyyyy, g_0_y_xz_xxyyyyyz, g_0_y_xz_xxyyyyz, g_0_y_xz_xxyyyyzz, g_0_y_xz_xxyyyzz, g_0_y_xz_xxyyyzzz, g_0_y_xz_xxyyzzz, g_0_y_xz_xxyyzzzz, g_0_y_xz_xxyzzzz, g_0_y_xz_xxyzzzzz, g_0_y_xz_xxzzzzz, g_0_y_xz_xxzzzzzz, g_0_y_xz_xyyyyyy, g_0_y_xz_xyyyyyyy, g_0_y_xz_xyyyyyyz, g_0_y_xz_xyyyyyz, g_0_y_xz_xyyyyyzz, g_0_y_xz_xyyyyzz, g_0_y_xz_xyyyyzzz, g_0_y_xz_xyyyzzz, g_0_y_xz_xyyyzzzz, g_0_y_xz_xyyzzzz, g_0_y_xz_xyyzzzzz, g_0_y_xz_xyzzzzz, g_0_y_xz_xyzzzzzz, g_0_y_xz_xzzzzzz, g_0_y_xz_xzzzzzzz, g_0_y_xz_yyyyyyy, g_0_y_xz_yyyyyyz, g_0_y_xz_yyyyyzz, g_0_y_xz_yyyyzzz, g_0_y_xz_yyyzzzz, g_0_y_xz_yyzzzzz, g_0_y_xz_yzzzzzz, g_0_y_xz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxz_xxxxxxx[k] = -g_0_y_xz_xxxxxxx[k] * ab_x + g_0_y_xz_xxxxxxxx[k];

                g_0_y_xxz_xxxxxxy[k] = -g_0_y_xz_xxxxxxy[k] * ab_x + g_0_y_xz_xxxxxxxy[k];

                g_0_y_xxz_xxxxxxz[k] = -g_0_y_xz_xxxxxxz[k] * ab_x + g_0_y_xz_xxxxxxxz[k];

                g_0_y_xxz_xxxxxyy[k] = -g_0_y_xz_xxxxxyy[k] * ab_x + g_0_y_xz_xxxxxxyy[k];

                g_0_y_xxz_xxxxxyz[k] = -g_0_y_xz_xxxxxyz[k] * ab_x + g_0_y_xz_xxxxxxyz[k];

                g_0_y_xxz_xxxxxzz[k] = -g_0_y_xz_xxxxxzz[k] * ab_x + g_0_y_xz_xxxxxxzz[k];

                g_0_y_xxz_xxxxyyy[k] = -g_0_y_xz_xxxxyyy[k] * ab_x + g_0_y_xz_xxxxxyyy[k];

                g_0_y_xxz_xxxxyyz[k] = -g_0_y_xz_xxxxyyz[k] * ab_x + g_0_y_xz_xxxxxyyz[k];

                g_0_y_xxz_xxxxyzz[k] = -g_0_y_xz_xxxxyzz[k] * ab_x + g_0_y_xz_xxxxxyzz[k];

                g_0_y_xxz_xxxxzzz[k] = -g_0_y_xz_xxxxzzz[k] * ab_x + g_0_y_xz_xxxxxzzz[k];

                g_0_y_xxz_xxxyyyy[k] = -g_0_y_xz_xxxyyyy[k] * ab_x + g_0_y_xz_xxxxyyyy[k];

                g_0_y_xxz_xxxyyyz[k] = -g_0_y_xz_xxxyyyz[k] * ab_x + g_0_y_xz_xxxxyyyz[k];

                g_0_y_xxz_xxxyyzz[k] = -g_0_y_xz_xxxyyzz[k] * ab_x + g_0_y_xz_xxxxyyzz[k];

                g_0_y_xxz_xxxyzzz[k] = -g_0_y_xz_xxxyzzz[k] * ab_x + g_0_y_xz_xxxxyzzz[k];

                g_0_y_xxz_xxxzzzz[k] = -g_0_y_xz_xxxzzzz[k] * ab_x + g_0_y_xz_xxxxzzzz[k];

                g_0_y_xxz_xxyyyyy[k] = -g_0_y_xz_xxyyyyy[k] * ab_x + g_0_y_xz_xxxyyyyy[k];

                g_0_y_xxz_xxyyyyz[k] = -g_0_y_xz_xxyyyyz[k] * ab_x + g_0_y_xz_xxxyyyyz[k];

                g_0_y_xxz_xxyyyzz[k] = -g_0_y_xz_xxyyyzz[k] * ab_x + g_0_y_xz_xxxyyyzz[k];

                g_0_y_xxz_xxyyzzz[k] = -g_0_y_xz_xxyyzzz[k] * ab_x + g_0_y_xz_xxxyyzzz[k];

                g_0_y_xxz_xxyzzzz[k] = -g_0_y_xz_xxyzzzz[k] * ab_x + g_0_y_xz_xxxyzzzz[k];

                g_0_y_xxz_xxzzzzz[k] = -g_0_y_xz_xxzzzzz[k] * ab_x + g_0_y_xz_xxxzzzzz[k];

                g_0_y_xxz_xyyyyyy[k] = -g_0_y_xz_xyyyyyy[k] * ab_x + g_0_y_xz_xxyyyyyy[k];

                g_0_y_xxz_xyyyyyz[k] = -g_0_y_xz_xyyyyyz[k] * ab_x + g_0_y_xz_xxyyyyyz[k];

                g_0_y_xxz_xyyyyzz[k] = -g_0_y_xz_xyyyyzz[k] * ab_x + g_0_y_xz_xxyyyyzz[k];

                g_0_y_xxz_xyyyzzz[k] = -g_0_y_xz_xyyyzzz[k] * ab_x + g_0_y_xz_xxyyyzzz[k];

                g_0_y_xxz_xyyzzzz[k] = -g_0_y_xz_xyyzzzz[k] * ab_x + g_0_y_xz_xxyyzzzz[k];

                g_0_y_xxz_xyzzzzz[k] = -g_0_y_xz_xyzzzzz[k] * ab_x + g_0_y_xz_xxyzzzzz[k];

                g_0_y_xxz_xzzzzzz[k] = -g_0_y_xz_xzzzzzz[k] * ab_x + g_0_y_xz_xxzzzzzz[k];

                g_0_y_xxz_yyyyyyy[k] = -g_0_y_xz_yyyyyyy[k] * ab_x + g_0_y_xz_xyyyyyyy[k];

                g_0_y_xxz_yyyyyyz[k] = -g_0_y_xz_yyyyyyz[k] * ab_x + g_0_y_xz_xyyyyyyz[k];

                g_0_y_xxz_yyyyyzz[k] = -g_0_y_xz_yyyyyzz[k] * ab_x + g_0_y_xz_xyyyyyzz[k];

                g_0_y_xxz_yyyyzzz[k] = -g_0_y_xz_yyyyzzz[k] * ab_x + g_0_y_xz_xyyyyzzz[k];

                g_0_y_xxz_yyyzzzz[k] = -g_0_y_xz_yyyzzzz[k] * ab_x + g_0_y_xz_xyyyzzzz[k];

                g_0_y_xxz_yyzzzzz[k] = -g_0_y_xz_yyzzzzz[k] * ab_x + g_0_y_xz_xyyzzzzz[k];

                g_0_y_xxz_yzzzzzz[k] = -g_0_y_xz_yzzzzzz[k] * ab_x + g_0_y_xz_xyzzzzzz[k];

                g_0_y_xxz_zzzzzzz[k] = -g_0_y_xz_zzzzzzz[k] * ab_x + g_0_y_xz_xzzzzzzz[k];
            }

            /// Set up 468-504 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyy_xxxxxxx, g_0_y_xyy_xxxxxxy, g_0_y_xyy_xxxxxxz, g_0_y_xyy_xxxxxyy, g_0_y_xyy_xxxxxyz, g_0_y_xyy_xxxxxzz, g_0_y_xyy_xxxxyyy, g_0_y_xyy_xxxxyyz, g_0_y_xyy_xxxxyzz, g_0_y_xyy_xxxxzzz, g_0_y_xyy_xxxyyyy, g_0_y_xyy_xxxyyyz, g_0_y_xyy_xxxyyzz, g_0_y_xyy_xxxyzzz, g_0_y_xyy_xxxzzzz, g_0_y_xyy_xxyyyyy, g_0_y_xyy_xxyyyyz, g_0_y_xyy_xxyyyzz, g_0_y_xyy_xxyyzzz, g_0_y_xyy_xxyzzzz, g_0_y_xyy_xxzzzzz, g_0_y_xyy_xyyyyyy, g_0_y_xyy_xyyyyyz, g_0_y_xyy_xyyyyzz, g_0_y_xyy_xyyyzzz, g_0_y_xyy_xyyzzzz, g_0_y_xyy_xyzzzzz, g_0_y_xyy_xzzzzzz, g_0_y_xyy_yyyyyyy, g_0_y_xyy_yyyyyyz, g_0_y_xyy_yyyyyzz, g_0_y_xyy_yyyyzzz, g_0_y_xyy_yyyzzzz, g_0_y_xyy_yyzzzzz, g_0_y_xyy_yzzzzzz, g_0_y_xyy_zzzzzzz, g_0_y_yy_xxxxxxx, g_0_y_yy_xxxxxxxx, g_0_y_yy_xxxxxxxy, g_0_y_yy_xxxxxxxz, g_0_y_yy_xxxxxxy, g_0_y_yy_xxxxxxyy, g_0_y_yy_xxxxxxyz, g_0_y_yy_xxxxxxz, g_0_y_yy_xxxxxxzz, g_0_y_yy_xxxxxyy, g_0_y_yy_xxxxxyyy, g_0_y_yy_xxxxxyyz, g_0_y_yy_xxxxxyz, g_0_y_yy_xxxxxyzz, g_0_y_yy_xxxxxzz, g_0_y_yy_xxxxxzzz, g_0_y_yy_xxxxyyy, g_0_y_yy_xxxxyyyy, g_0_y_yy_xxxxyyyz, g_0_y_yy_xxxxyyz, g_0_y_yy_xxxxyyzz, g_0_y_yy_xxxxyzz, g_0_y_yy_xxxxyzzz, g_0_y_yy_xxxxzzz, g_0_y_yy_xxxxzzzz, g_0_y_yy_xxxyyyy, g_0_y_yy_xxxyyyyy, g_0_y_yy_xxxyyyyz, g_0_y_yy_xxxyyyz, g_0_y_yy_xxxyyyzz, g_0_y_yy_xxxyyzz, g_0_y_yy_xxxyyzzz, g_0_y_yy_xxxyzzz, g_0_y_yy_xxxyzzzz, g_0_y_yy_xxxzzzz, g_0_y_yy_xxxzzzzz, g_0_y_yy_xxyyyyy, g_0_y_yy_xxyyyyyy, g_0_y_yy_xxyyyyyz, g_0_y_yy_xxyyyyz, g_0_y_yy_xxyyyyzz, g_0_y_yy_xxyyyzz, g_0_y_yy_xxyyyzzz, g_0_y_yy_xxyyzzz, g_0_y_yy_xxyyzzzz, g_0_y_yy_xxyzzzz, g_0_y_yy_xxyzzzzz, g_0_y_yy_xxzzzzz, g_0_y_yy_xxzzzzzz, g_0_y_yy_xyyyyyy, g_0_y_yy_xyyyyyyy, g_0_y_yy_xyyyyyyz, g_0_y_yy_xyyyyyz, g_0_y_yy_xyyyyyzz, g_0_y_yy_xyyyyzz, g_0_y_yy_xyyyyzzz, g_0_y_yy_xyyyzzz, g_0_y_yy_xyyyzzzz, g_0_y_yy_xyyzzzz, g_0_y_yy_xyyzzzzz, g_0_y_yy_xyzzzzz, g_0_y_yy_xyzzzzzz, g_0_y_yy_xzzzzzz, g_0_y_yy_xzzzzzzz, g_0_y_yy_yyyyyyy, g_0_y_yy_yyyyyyz, g_0_y_yy_yyyyyzz, g_0_y_yy_yyyyzzz, g_0_y_yy_yyyzzzz, g_0_y_yy_yyzzzzz, g_0_y_yy_yzzzzzz, g_0_y_yy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyy_xxxxxxx[k] = -g_0_y_yy_xxxxxxx[k] * ab_x + g_0_y_yy_xxxxxxxx[k];

                g_0_y_xyy_xxxxxxy[k] = -g_0_y_yy_xxxxxxy[k] * ab_x + g_0_y_yy_xxxxxxxy[k];

                g_0_y_xyy_xxxxxxz[k] = -g_0_y_yy_xxxxxxz[k] * ab_x + g_0_y_yy_xxxxxxxz[k];

                g_0_y_xyy_xxxxxyy[k] = -g_0_y_yy_xxxxxyy[k] * ab_x + g_0_y_yy_xxxxxxyy[k];

                g_0_y_xyy_xxxxxyz[k] = -g_0_y_yy_xxxxxyz[k] * ab_x + g_0_y_yy_xxxxxxyz[k];

                g_0_y_xyy_xxxxxzz[k] = -g_0_y_yy_xxxxxzz[k] * ab_x + g_0_y_yy_xxxxxxzz[k];

                g_0_y_xyy_xxxxyyy[k] = -g_0_y_yy_xxxxyyy[k] * ab_x + g_0_y_yy_xxxxxyyy[k];

                g_0_y_xyy_xxxxyyz[k] = -g_0_y_yy_xxxxyyz[k] * ab_x + g_0_y_yy_xxxxxyyz[k];

                g_0_y_xyy_xxxxyzz[k] = -g_0_y_yy_xxxxyzz[k] * ab_x + g_0_y_yy_xxxxxyzz[k];

                g_0_y_xyy_xxxxzzz[k] = -g_0_y_yy_xxxxzzz[k] * ab_x + g_0_y_yy_xxxxxzzz[k];

                g_0_y_xyy_xxxyyyy[k] = -g_0_y_yy_xxxyyyy[k] * ab_x + g_0_y_yy_xxxxyyyy[k];

                g_0_y_xyy_xxxyyyz[k] = -g_0_y_yy_xxxyyyz[k] * ab_x + g_0_y_yy_xxxxyyyz[k];

                g_0_y_xyy_xxxyyzz[k] = -g_0_y_yy_xxxyyzz[k] * ab_x + g_0_y_yy_xxxxyyzz[k];

                g_0_y_xyy_xxxyzzz[k] = -g_0_y_yy_xxxyzzz[k] * ab_x + g_0_y_yy_xxxxyzzz[k];

                g_0_y_xyy_xxxzzzz[k] = -g_0_y_yy_xxxzzzz[k] * ab_x + g_0_y_yy_xxxxzzzz[k];

                g_0_y_xyy_xxyyyyy[k] = -g_0_y_yy_xxyyyyy[k] * ab_x + g_0_y_yy_xxxyyyyy[k];

                g_0_y_xyy_xxyyyyz[k] = -g_0_y_yy_xxyyyyz[k] * ab_x + g_0_y_yy_xxxyyyyz[k];

                g_0_y_xyy_xxyyyzz[k] = -g_0_y_yy_xxyyyzz[k] * ab_x + g_0_y_yy_xxxyyyzz[k];

                g_0_y_xyy_xxyyzzz[k] = -g_0_y_yy_xxyyzzz[k] * ab_x + g_0_y_yy_xxxyyzzz[k];

                g_0_y_xyy_xxyzzzz[k] = -g_0_y_yy_xxyzzzz[k] * ab_x + g_0_y_yy_xxxyzzzz[k];

                g_0_y_xyy_xxzzzzz[k] = -g_0_y_yy_xxzzzzz[k] * ab_x + g_0_y_yy_xxxzzzzz[k];

                g_0_y_xyy_xyyyyyy[k] = -g_0_y_yy_xyyyyyy[k] * ab_x + g_0_y_yy_xxyyyyyy[k];

                g_0_y_xyy_xyyyyyz[k] = -g_0_y_yy_xyyyyyz[k] * ab_x + g_0_y_yy_xxyyyyyz[k];

                g_0_y_xyy_xyyyyzz[k] = -g_0_y_yy_xyyyyzz[k] * ab_x + g_0_y_yy_xxyyyyzz[k];

                g_0_y_xyy_xyyyzzz[k] = -g_0_y_yy_xyyyzzz[k] * ab_x + g_0_y_yy_xxyyyzzz[k];

                g_0_y_xyy_xyyzzzz[k] = -g_0_y_yy_xyyzzzz[k] * ab_x + g_0_y_yy_xxyyzzzz[k];

                g_0_y_xyy_xyzzzzz[k] = -g_0_y_yy_xyzzzzz[k] * ab_x + g_0_y_yy_xxyzzzzz[k];

                g_0_y_xyy_xzzzzzz[k] = -g_0_y_yy_xzzzzzz[k] * ab_x + g_0_y_yy_xxzzzzzz[k];

                g_0_y_xyy_yyyyyyy[k] = -g_0_y_yy_yyyyyyy[k] * ab_x + g_0_y_yy_xyyyyyyy[k];

                g_0_y_xyy_yyyyyyz[k] = -g_0_y_yy_yyyyyyz[k] * ab_x + g_0_y_yy_xyyyyyyz[k];

                g_0_y_xyy_yyyyyzz[k] = -g_0_y_yy_yyyyyzz[k] * ab_x + g_0_y_yy_xyyyyyzz[k];

                g_0_y_xyy_yyyyzzz[k] = -g_0_y_yy_yyyyzzz[k] * ab_x + g_0_y_yy_xyyyyzzz[k];

                g_0_y_xyy_yyyzzzz[k] = -g_0_y_yy_yyyzzzz[k] * ab_x + g_0_y_yy_xyyyzzzz[k];

                g_0_y_xyy_yyzzzzz[k] = -g_0_y_yy_yyzzzzz[k] * ab_x + g_0_y_yy_xyyzzzzz[k];

                g_0_y_xyy_yzzzzzz[k] = -g_0_y_yy_yzzzzzz[k] * ab_x + g_0_y_yy_xyzzzzzz[k];

                g_0_y_xyy_zzzzzzz[k] = -g_0_y_yy_zzzzzzz[k] * ab_x + g_0_y_yy_xzzzzzzz[k];
            }

            /// Set up 504-540 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyz_xxxxxxx, g_0_y_xyz_xxxxxxy, g_0_y_xyz_xxxxxxz, g_0_y_xyz_xxxxxyy, g_0_y_xyz_xxxxxyz, g_0_y_xyz_xxxxxzz, g_0_y_xyz_xxxxyyy, g_0_y_xyz_xxxxyyz, g_0_y_xyz_xxxxyzz, g_0_y_xyz_xxxxzzz, g_0_y_xyz_xxxyyyy, g_0_y_xyz_xxxyyyz, g_0_y_xyz_xxxyyzz, g_0_y_xyz_xxxyzzz, g_0_y_xyz_xxxzzzz, g_0_y_xyz_xxyyyyy, g_0_y_xyz_xxyyyyz, g_0_y_xyz_xxyyyzz, g_0_y_xyz_xxyyzzz, g_0_y_xyz_xxyzzzz, g_0_y_xyz_xxzzzzz, g_0_y_xyz_xyyyyyy, g_0_y_xyz_xyyyyyz, g_0_y_xyz_xyyyyzz, g_0_y_xyz_xyyyzzz, g_0_y_xyz_xyyzzzz, g_0_y_xyz_xyzzzzz, g_0_y_xyz_xzzzzzz, g_0_y_xyz_yyyyyyy, g_0_y_xyz_yyyyyyz, g_0_y_xyz_yyyyyzz, g_0_y_xyz_yyyyzzz, g_0_y_xyz_yyyzzzz, g_0_y_xyz_yyzzzzz, g_0_y_xyz_yzzzzzz, g_0_y_xyz_zzzzzzz, g_0_y_yz_xxxxxxx, g_0_y_yz_xxxxxxxx, g_0_y_yz_xxxxxxxy, g_0_y_yz_xxxxxxxz, g_0_y_yz_xxxxxxy, g_0_y_yz_xxxxxxyy, g_0_y_yz_xxxxxxyz, g_0_y_yz_xxxxxxz, g_0_y_yz_xxxxxxzz, g_0_y_yz_xxxxxyy, g_0_y_yz_xxxxxyyy, g_0_y_yz_xxxxxyyz, g_0_y_yz_xxxxxyz, g_0_y_yz_xxxxxyzz, g_0_y_yz_xxxxxzz, g_0_y_yz_xxxxxzzz, g_0_y_yz_xxxxyyy, g_0_y_yz_xxxxyyyy, g_0_y_yz_xxxxyyyz, g_0_y_yz_xxxxyyz, g_0_y_yz_xxxxyyzz, g_0_y_yz_xxxxyzz, g_0_y_yz_xxxxyzzz, g_0_y_yz_xxxxzzz, g_0_y_yz_xxxxzzzz, g_0_y_yz_xxxyyyy, g_0_y_yz_xxxyyyyy, g_0_y_yz_xxxyyyyz, g_0_y_yz_xxxyyyz, g_0_y_yz_xxxyyyzz, g_0_y_yz_xxxyyzz, g_0_y_yz_xxxyyzzz, g_0_y_yz_xxxyzzz, g_0_y_yz_xxxyzzzz, g_0_y_yz_xxxzzzz, g_0_y_yz_xxxzzzzz, g_0_y_yz_xxyyyyy, g_0_y_yz_xxyyyyyy, g_0_y_yz_xxyyyyyz, g_0_y_yz_xxyyyyz, g_0_y_yz_xxyyyyzz, g_0_y_yz_xxyyyzz, g_0_y_yz_xxyyyzzz, g_0_y_yz_xxyyzzz, g_0_y_yz_xxyyzzzz, g_0_y_yz_xxyzzzz, g_0_y_yz_xxyzzzzz, g_0_y_yz_xxzzzzz, g_0_y_yz_xxzzzzzz, g_0_y_yz_xyyyyyy, g_0_y_yz_xyyyyyyy, g_0_y_yz_xyyyyyyz, g_0_y_yz_xyyyyyz, g_0_y_yz_xyyyyyzz, g_0_y_yz_xyyyyzz, g_0_y_yz_xyyyyzzz, g_0_y_yz_xyyyzzz, g_0_y_yz_xyyyzzzz, g_0_y_yz_xyyzzzz, g_0_y_yz_xyyzzzzz, g_0_y_yz_xyzzzzz, g_0_y_yz_xyzzzzzz, g_0_y_yz_xzzzzzz, g_0_y_yz_xzzzzzzz, g_0_y_yz_yyyyyyy, g_0_y_yz_yyyyyyz, g_0_y_yz_yyyyyzz, g_0_y_yz_yyyyzzz, g_0_y_yz_yyyzzzz, g_0_y_yz_yyzzzzz, g_0_y_yz_yzzzzzz, g_0_y_yz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyz_xxxxxxx[k] = -g_0_y_yz_xxxxxxx[k] * ab_x + g_0_y_yz_xxxxxxxx[k];

                g_0_y_xyz_xxxxxxy[k] = -g_0_y_yz_xxxxxxy[k] * ab_x + g_0_y_yz_xxxxxxxy[k];

                g_0_y_xyz_xxxxxxz[k] = -g_0_y_yz_xxxxxxz[k] * ab_x + g_0_y_yz_xxxxxxxz[k];

                g_0_y_xyz_xxxxxyy[k] = -g_0_y_yz_xxxxxyy[k] * ab_x + g_0_y_yz_xxxxxxyy[k];

                g_0_y_xyz_xxxxxyz[k] = -g_0_y_yz_xxxxxyz[k] * ab_x + g_0_y_yz_xxxxxxyz[k];

                g_0_y_xyz_xxxxxzz[k] = -g_0_y_yz_xxxxxzz[k] * ab_x + g_0_y_yz_xxxxxxzz[k];

                g_0_y_xyz_xxxxyyy[k] = -g_0_y_yz_xxxxyyy[k] * ab_x + g_0_y_yz_xxxxxyyy[k];

                g_0_y_xyz_xxxxyyz[k] = -g_0_y_yz_xxxxyyz[k] * ab_x + g_0_y_yz_xxxxxyyz[k];

                g_0_y_xyz_xxxxyzz[k] = -g_0_y_yz_xxxxyzz[k] * ab_x + g_0_y_yz_xxxxxyzz[k];

                g_0_y_xyz_xxxxzzz[k] = -g_0_y_yz_xxxxzzz[k] * ab_x + g_0_y_yz_xxxxxzzz[k];

                g_0_y_xyz_xxxyyyy[k] = -g_0_y_yz_xxxyyyy[k] * ab_x + g_0_y_yz_xxxxyyyy[k];

                g_0_y_xyz_xxxyyyz[k] = -g_0_y_yz_xxxyyyz[k] * ab_x + g_0_y_yz_xxxxyyyz[k];

                g_0_y_xyz_xxxyyzz[k] = -g_0_y_yz_xxxyyzz[k] * ab_x + g_0_y_yz_xxxxyyzz[k];

                g_0_y_xyz_xxxyzzz[k] = -g_0_y_yz_xxxyzzz[k] * ab_x + g_0_y_yz_xxxxyzzz[k];

                g_0_y_xyz_xxxzzzz[k] = -g_0_y_yz_xxxzzzz[k] * ab_x + g_0_y_yz_xxxxzzzz[k];

                g_0_y_xyz_xxyyyyy[k] = -g_0_y_yz_xxyyyyy[k] * ab_x + g_0_y_yz_xxxyyyyy[k];

                g_0_y_xyz_xxyyyyz[k] = -g_0_y_yz_xxyyyyz[k] * ab_x + g_0_y_yz_xxxyyyyz[k];

                g_0_y_xyz_xxyyyzz[k] = -g_0_y_yz_xxyyyzz[k] * ab_x + g_0_y_yz_xxxyyyzz[k];

                g_0_y_xyz_xxyyzzz[k] = -g_0_y_yz_xxyyzzz[k] * ab_x + g_0_y_yz_xxxyyzzz[k];

                g_0_y_xyz_xxyzzzz[k] = -g_0_y_yz_xxyzzzz[k] * ab_x + g_0_y_yz_xxxyzzzz[k];

                g_0_y_xyz_xxzzzzz[k] = -g_0_y_yz_xxzzzzz[k] * ab_x + g_0_y_yz_xxxzzzzz[k];

                g_0_y_xyz_xyyyyyy[k] = -g_0_y_yz_xyyyyyy[k] * ab_x + g_0_y_yz_xxyyyyyy[k];

                g_0_y_xyz_xyyyyyz[k] = -g_0_y_yz_xyyyyyz[k] * ab_x + g_0_y_yz_xxyyyyyz[k];

                g_0_y_xyz_xyyyyzz[k] = -g_0_y_yz_xyyyyzz[k] * ab_x + g_0_y_yz_xxyyyyzz[k];

                g_0_y_xyz_xyyyzzz[k] = -g_0_y_yz_xyyyzzz[k] * ab_x + g_0_y_yz_xxyyyzzz[k];

                g_0_y_xyz_xyyzzzz[k] = -g_0_y_yz_xyyzzzz[k] * ab_x + g_0_y_yz_xxyyzzzz[k];

                g_0_y_xyz_xyzzzzz[k] = -g_0_y_yz_xyzzzzz[k] * ab_x + g_0_y_yz_xxyzzzzz[k];

                g_0_y_xyz_xzzzzzz[k] = -g_0_y_yz_xzzzzzz[k] * ab_x + g_0_y_yz_xxzzzzzz[k];

                g_0_y_xyz_yyyyyyy[k] = -g_0_y_yz_yyyyyyy[k] * ab_x + g_0_y_yz_xyyyyyyy[k];

                g_0_y_xyz_yyyyyyz[k] = -g_0_y_yz_yyyyyyz[k] * ab_x + g_0_y_yz_xyyyyyyz[k];

                g_0_y_xyz_yyyyyzz[k] = -g_0_y_yz_yyyyyzz[k] * ab_x + g_0_y_yz_xyyyyyzz[k];

                g_0_y_xyz_yyyyzzz[k] = -g_0_y_yz_yyyyzzz[k] * ab_x + g_0_y_yz_xyyyyzzz[k];

                g_0_y_xyz_yyyzzzz[k] = -g_0_y_yz_yyyzzzz[k] * ab_x + g_0_y_yz_xyyyzzzz[k];

                g_0_y_xyz_yyzzzzz[k] = -g_0_y_yz_yyzzzzz[k] * ab_x + g_0_y_yz_xyyzzzzz[k];

                g_0_y_xyz_yzzzzzz[k] = -g_0_y_yz_yzzzzzz[k] * ab_x + g_0_y_yz_xyzzzzzz[k];

                g_0_y_xyz_zzzzzzz[k] = -g_0_y_yz_zzzzzzz[k] * ab_x + g_0_y_yz_xzzzzzzz[k];
            }

            /// Set up 540-576 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xzz_xxxxxxx, g_0_y_xzz_xxxxxxy, g_0_y_xzz_xxxxxxz, g_0_y_xzz_xxxxxyy, g_0_y_xzz_xxxxxyz, g_0_y_xzz_xxxxxzz, g_0_y_xzz_xxxxyyy, g_0_y_xzz_xxxxyyz, g_0_y_xzz_xxxxyzz, g_0_y_xzz_xxxxzzz, g_0_y_xzz_xxxyyyy, g_0_y_xzz_xxxyyyz, g_0_y_xzz_xxxyyzz, g_0_y_xzz_xxxyzzz, g_0_y_xzz_xxxzzzz, g_0_y_xzz_xxyyyyy, g_0_y_xzz_xxyyyyz, g_0_y_xzz_xxyyyzz, g_0_y_xzz_xxyyzzz, g_0_y_xzz_xxyzzzz, g_0_y_xzz_xxzzzzz, g_0_y_xzz_xyyyyyy, g_0_y_xzz_xyyyyyz, g_0_y_xzz_xyyyyzz, g_0_y_xzz_xyyyzzz, g_0_y_xzz_xyyzzzz, g_0_y_xzz_xyzzzzz, g_0_y_xzz_xzzzzzz, g_0_y_xzz_yyyyyyy, g_0_y_xzz_yyyyyyz, g_0_y_xzz_yyyyyzz, g_0_y_xzz_yyyyzzz, g_0_y_xzz_yyyzzzz, g_0_y_xzz_yyzzzzz, g_0_y_xzz_yzzzzzz, g_0_y_xzz_zzzzzzz, g_0_y_zz_xxxxxxx, g_0_y_zz_xxxxxxxx, g_0_y_zz_xxxxxxxy, g_0_y_zz_xxxxxxxz, g_0_y_zz_xxxxxxy, g_0_y_zz_xxxxxxyy, g_0_y_zz_xxxxxxyz, g_0_y_zz_xxxxxxz, g_0_y_zz_xxxxxxzz, g_0_y_zz_xxxxxyy, g_0_y_zz_xxxxxyyy, g_0_y_zz_xxxxxyyz, g_0_y_zz_xxxxxyz, g_0_y_zz_xxxxxyzz, g_0_y_zz_xxxxxzz, g_0_y_zz_xxxxxzzz, g_0_y_zz_xxxxyyy, g_0_y_zz_xxxxyyyy, g_0_y_zz_xxxxyyyz, g_0_y_zz_xxxxyyz, g_0_y_zz_xxxxyyzz, g_0_y_zz_xxxxyzz, g_0_y_zz_xxxxyzzz, g_0_y_zz_xxxxzzz, g_0_y_zz_xxxxzzzz, g_0_y_zz_xxxyyyy, g_0_y_zz_xxxyyyyy, g_0_y_zz_xxxyyyyz, g_0_y_zz_xxxyyyz, g_0_y_zz_xxxyyyzz, g_0_y_zz_xxxyyzz, g_0_y_zz_xxxyyzzz, g_0_y_zz_xxxyzzz, g_0_y_zz_xxxyzzzz, g_0_y_zz_xxxzzzz, g_0_y_zz_xxxzzzzz, g_0_y_zz_xxyyyyy, g_0_y_zz_xxyyyyyy, g_0_y_zz_xxyyyyyz, g_0_y_zz_xxyyyyz, g_0_y_zz_xxyyyyzz, g_0_y_zz_xxyyyzz, g_0_y_zz_xxyyyzzz, g_0_y_zz_xxyyzzz, g_0_y_zz_xxyyzzzz, g_0_y_zz_xxyzzzz, g_0_y_zz_xxyzzzzz, g_0_y_zz_xxzzzzz, g_0_y_zz_xxzzzzzz, g_0_y_zz_xyyyyyy, g_0_y_zz_xyyyyyyy, g_0_y_zz_xyyyyyyz, g_0_y_zz_xyyyyyz, g_0_y_zz_xyyyyyzz, g_0_y_zz_xyyyyzz, g_0_y_zz_xyyyyzzz, g_0_y_zz_xyyyzzz, g_0_y_zz_xyyyzzzz, g_0_y_zz_xyyzzzz, g_0_y_zz_xyyzzzzz, g_0_y_zz_xyzzzzz, g_0_y_zz_xyzzzzzz, g_0_y_zz_xzzzzzz, g_0_y_zz_xzzzzzzz, g_0_y_zz_yyyyyyy, g_0_y_zz_yyyyyyz, g_0_y_zz_yyyyyzz, g_0_y_zz_yyyyzzz, g_0_y_zz_yyyzzzz, g_0_y_zz_yyzzzzz, g_0_y_zz_yzzzzzz, g_0_y_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzz_xxxxxxx[k] = -g_0_y_zz_xxxxxxx[k] * ab_x + g_0_y_zz_xxxxxxxx[k];

                g_0_y_xzz_xxxxxxy[k] = -g_0_y_zz_xxxxxxy[k] * ab_x + g_0_y_zz_xxxxxxxy[k];

                g_0_y_xzz_xxxxxxz[k] = -g_0_y_zz_xxxxxxz[k] * ab_x + g_0_y_zz_xxxxxxxz[k];

                g_0_y_xzz_xxxxxyy[k] = -g_0_y_zz_xxxxxyy[k] * ab_x + g_0_y_zz_xxxxxxyy[k];

                g_0_y_xzz_xxxxxyz[k] = -g_0_y_zz_xxxxxyz[k] * ab_x + g_0_y_zz_xxxxxxyz[k];

                g_0_y_xzz_xxxxxzz[k] = -g_0_y_zz_xxxxxzz[k] * ab_x + g_0_y_zz_xxxxxxzz[k];

                g_0_y_xzz_xxxxyyy[k] = -g_0_y_zz_xxxxyyy[k] * ab_x + g_0_y_zz_xxxxxyyy[k];

                g_0_y_xzz_xxxxyyz[k] = -g_0_y_zz_xxxxyyz[k] * ab_x + g_0_y_zz_xxxxxyyz[k];

                g_0_y_xzz_xxxxyzz[k] = -g_0_y_zz_xxxxyzz[k] * ab_x + g_0_y_zz_xxxxxyzz[k];

                g_0_y_xzz_xxxxzzz[k] = -g_0_y_zz_xxxxzzz[k] * ab_x + g_0_y_zz_xxxxxzzz[k];

                g_0_y_xzz_xxxyyyy[k] = -g_0_y_zz_xxxyyyy[k] * ab_x + g_0_y_zz_xxxxyyyy[k];

                g_0_y_xzz_xxxyyyz[k] = -g_0_y_zz_xxxyyyz[k] * ab_x + g_0_y_zz_xxxxyyyz[k];

                g_0_y_xzz_xxxyyzz[k] = -g_0_y_zz_xxxyyzz[k] * ab_x + g_0_y_zz_xxxxyyzz[k];

                g_0_y_xzz_xxxyzzz[k] = -g_0_y_zz_xxxyzzz[k] * ab_x + g_0_y_zz_xxxxyzzz[k];

                g_0_y_xzz_xxxzzzz[k] = -g_0_y_zz_xxxzzzz[k] * ab_x + g_0_y_zz_xxxxzzzz[k];

                g_0_y_xzz_xxyyyyy[k] = -g_0_y_zz_xxyyyyy[k] * ab_x + g_0_y_zz_xxxyyyyy[k];

                g_0_y_xzz_xxyyyyz[k] = -g_0_y_zz_xxyyyyz[k] * ab_x + g_0_y_zz_xxxyyyyz[k];

                g_0_y_xzz_xxyyyzz[k] = -g_0_y_zz_xxyyyzz[k] * ab_x + g_0_y_zz_xxxyyyzz[k];

                g_0_y_xzz_xxyyzzz[k] = -g_0_y_zz_xxyyzzz[k] * ab_x + g_0_y_zz_xxxyyzzz[k];

                g_0_y_xzz_xxyzzzz[k] = -g_0_y_zz_xxyzzzz[k] * ab_x + g_0_y_zz_xxxyzzzz[k];

                g_0_y_xzz_xxzzzzz[k] = -g_0_y_zz_xxzzzzz[k] * ab_x + g_0_y_zz_xxxzzzzz[k];

                g_0_y_xzz_xyyyyyy[k] = -g_0_y_zz_xyyyyyy[k] * ab_x + g_0_y_zz_xxyyyyyy[k];

                g_0_y_xzz_xyyyyyz[k] = -g_0_y_zz_xyyyyyz[k] * ab_x + g_0_y_zz_xxyyyyyz[k];

                g_0_y_xzz_xyyyyzz[k] = -g_0_y_zz_xyyyyzz[k] * ab_x + g_0_y_zz_xxyyyyzz[k];

                g_0_y_xzz_xyyyzzz[k] = -g_0_y_zz_xyyyzzz[k] * ab_x + g_0_y_zz_xxyyyzzz[k];

                g_0_y_xzz_xyyzzzz[k] = -g_0_y_zz_xyyzzzz[k] * ab_x + g_0_y_zz_xxyyzzzz[k];

                g_0_y_xzz_xyzzzzz[k] = -g_0_y_zz_xyzzzzz[k] * ab_x + g_0_y_zz_xxyzzzzz[k];

                g_0_y_xzz_xzzzzzz[k] = -g_0_y_zz_xzzzzzz[k] * ab_x + g_0_y_zz_xxzzzzzz[k];

                g_0_y_xzz_yyyyyyy[k] = -g_0_y_zz_yyyyyyy[k] * ab_x + g_0_y_zz_xyyyyyyy[k];

                g_0_y_xzz_yyyyyyz[k] = -g_0_y_zz_yyyyyyz[k] * ab_x + g_0_y_zz_xyyyyyyz[k];

                g_0_y_xzz_yyyyyzz[k] = -g_0_y_zz_yyyyyzz[k] * ab_x + g_0_y_zz_xyyyyyzz[k];

                g_0_y_xzz_yyyyzzz[k] = -g_0_y_zz_yyyyzzz[k] * ab_x + g_0_y_zz_xyyyyzzz[k];

                g_0_y_xzz_yyyzzzz[k] = -g_0_y_zz_yyyzzzz[k] * ab_x + g_0_y_zz_xyyyzzzz[k];

                g_0_y_xzz_yyzzzzz[k] = -g_0_y_zz_yyzzzzz[k] * ab_x + g_0_y_zz_xyyzzzzz[k];

                g_0_y_xzz_yzzzzzz[k] = -g_0_y_zz_yzzzzzz[k] * ab_x + g_0_y_zz_xyzzzzzz[k];

                g_0_y_xzz_zzzzzzz[k] = -g_0_y_zz_zzzzzzz[k] * ab_x + g_0_y_zz_xzzzzzzz[k];
            }

            /// Set up 576-612 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yy_xxxxxxx, g_0_y_yy_xxxxxxxy, g_0_y_yy_xxxxxxy, g_0_y_yy_xxxxxxyy, g_0_y_yy_xxxxxxyz, g_0_y_yy_xxxxxxz, g_0_y_yy_xxxxxyy, g_0_y_yy_xxxxxyyy, g_0_y_yy_xxxxxyyz, g_0_y_yy_xxxxxyz, g_0_y_yy_xxxxxyzz, g_0_y_yy_xxxxxzz, g_0_y_yy_xxxxyyy, g_0_y_yy_xxxxyyyy, g_0_y_yy_xxxxyyyz, g_0_y_yy_xxxxyyz, g_0_y_yy_xxxxyyzz, g_0_y_yy_xxxxyzz, g_0_y_yy_xxxxyzzz, g_0_y_yy_xxxxzzz, g_0_y_yy_xxxyyyy, g_0_y_yy_xxxyyyyy, g_0_y_yy_xxxyyyyz, g_0_y_yy_xxxyyyz, g_0_y_yy_xxxyyyzz, g_0_y_yy_xxxyyzz, g_0_y_yy_xxxyyzzz, g_0_y_yy_xxxyzzz, g_0_y_yy_xxxyzzzz, g_0_y_yy_xxxzzzz, g_0_y_yy_xxyyyyy, g_0_y_yy_xxyyyyyy, g_0_y_yy_xxyyyyyz, g_0_y_yy_xxyyyyz, g_0_y_yy_xxyyyyzz, g_0_y_yy_xxyyyzz, g_0_y_yy_xxyyyzzz, g_0_y_yy_xxyyzzz, g_0_y_yy_xxyyzzzz, g_0_y_yy_xxyzzzz, g_0_y_yy_xxyzzzzz, g_0_y_yy_xxzzzzz, g_0_y_yy_xyyyyyy, g_0_y_yy_xyyyyyyy, g_0_y_yy_xyyyyyyz, g_0_y_yy_xyyyyyz, g_0_y_yy_xyyyyyzz, g_0_y_yy_xyyyyzz, g_0_y_yy_xyyyyzzz, g_0_y_yy_xyyyzzz, g_0_y_yy_xyyyzzzz, g_0_y_yy_xyyzzzz, g_0_y_yy_xyyzzzzz, g_0_y_yy_xyzzzzz, g_0_y_yy_xyzzzzzz, g_0_y_yy_xzzzzzz, g_0_y_yy_yyyyyyy, g_0_y_yy_yyyyyyyy, g_0_y_yy_yyyyyyyz, g_0_y_yy_yyyyyyz, g_0_y_yy_yyyyyyzz, g_0_y_yy_yyyyyzz, g_0_y_yy_yyyyyzzz, g_0_y_yy_yyyyzzz, g_0_y_yy_yyyyzzzz, g_0_y_yy_yyyzzzz, g_0_y_yy_yyyzzzzz, g_0_y_yy_yyzzzzz, g_0_y_yy_yyzzzzzz, g_0_y_yy_yzzzzzz, g_0_y_yy_yzzzzzzz, g_0_y_yy_zzzzzzz, g_0_y_yyy_xxxxxxx, g_0_y_yyy_xxxxxxy, g_0_y_yyy_xxxxxxz, g_0_y_yyy_xxxxxyy, g_0_y_yyy_xxxxxyz, g_0_y_yyy_xxxxxzz, g_0_y_yyy_xxxxyyy, g_0_y_yyy_xxxxyyz, g_0_y_yyy_xxxxyzz, g_0_y_yyy_xxxxzzz, g_0_y_yyy_xxxyyyy, g_0_y_yyy_xxxyyyz, g_0_y_yyy_xxxyyzz, g_0_y_yyy_xxxyzzz, g_0_y_yyy_xxxzzzz, g_0_y_yyy_xxyyyyy, g_0_y_yyy_xxyyyyz, g_0_y_yyy_xxyyyzz, g_0_y_yyy_xxyyzzz, g_0_y_yyy_xxyzzzz, g_0_y_yyy_xxzzzzz, g_0_y_yyy_xyyyyyy, g_0_y_yyy_xyyyyyz, g_0_y_yyy_xyyyyzz, g_0_y_yyy_xyyyzzz, g_0_y_yyy_xyyzzzz, g_0_y_yyy_xyzzzzz, g_0_y_yyy_xzzzzzz, g_0_y_yyy_yyyyyyy, g_0_y_yyy_yyyyyyz, g_0_y_yyy_yyyyyzz, g_0_y_yyy_yyyyzzz, g_0_y_yyy_yyyzzzz, g_0_y_yyy_yyzzzzz, g_0_y_yyy_yzzzzzz, g_0_y_yyy_zzzzzzz, g_yy_xxxxxxx, g_yy_xxxxxxy, g_yy_xxxxxxz, g_yy_xxxxxyy, g_yy_xxxxxyz, g_yy_xxxxxzz, g_yy_xxxxyyy, g_yy_xxxxyyz, g_yy_xxxxyzz, g_yy_xxxxzzz, g_yy_xxxyyyy, g_yy_xxxyyyz, g_yy_xxxyyzz, g_yy_xxxyzzz, g_yy_xxxzzzz, g_yy_xxyyyyy, g_yy_xxyyyyz, g_yy_xxyyyzz, g_yy_xxyyzzz, g_yy_xxyzzzz, g_yy_xxzzzzz, g_yy_xyyyyyy, g_yy_xyyyyyz, g_yy_xyyyyzz, g_yy_xyyyzzz, g_yy_xyyzzzz, g_yy_xyzzzzz, g_yy_xzzzzzz, g_yy_yyyyyyy, g_yy_yyyyyyz, g_yy_yyyyyzz, g_yy_yyyyzzz, g_yy_yyyzzzz, g_yy_yyzzzzz, g_yy_yzzzzzz, g_yy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyy_xxxxxxx[k] = g_yy_xxxxxxx[k] - g_0_y_yy_xxxxxxx[k] * ab_y + g_0_y_yy_xxxxxxxy[k];

                g_0_y_yyy_xxxxxxy[k] = g_yy_xxxxxxy[k] - g_0_y_yy_xxxxxxy[k] * ab_y + g_0_y_yy_xxxxxxyy[k];

                g_0_y_yyy_xxxxxxz[k] = g_yy_xxxxxxz[k] - g_0_y_yy_xxxxxxz[k] * ab_y + g_0_y_yy_xxxxxxyz[k];

                g_0_y_yyy_xxxxxyy[k] = g_yy_xxxxxyy[k] - g_0_y_yy_xxxxxyy[k] * ab_y + g_0_y_yy_xxxxxyyy[k];

                g_0_y_yyy_xxxxxyz[k] = g_yy_xxxxxyz[k] - g_0_y_yy_xxxxxyz[k] * ab_y + g_0_y_yy_xxxxxyyz[k];

                g_0_y_yyy_xxxxxzz[k] = g_yy_xxxxxzz[k] - g_0_y_yy_xxxxxzz[k] * ab_y + g_0_y_yy_xxxxxyzz[k];

                g_0_y_yyy_xxxxyyy[k] = g_yy_xxxxyyy[k] - g_0_y_yy_xxxxyyy[k] * ab_y + g_0_y_yy_xxxxyyyy[k];

                g_0_y_yyy_xxxxyyz[k] = g_yy_xxxxyyz[k] - g_0_y_yy_xxxxyyz[k] * ab_y + g_0_y_yy_xxxxyyyz[k];

                g_0_y_yyy_xxxxyzz[k] = g_yy_xxxxyzz[k] - g_0_y_yy_xxxxyzz[k] * ab_y + g_0_y_yy_xxxxyyzz[k];

                g_0_y_yyy_xxxxzzz[k] = g_yy_xxxxzzz[k] - g_0_y_yy_xxxxzzz[k] * ab_y + g_0_y_yy_xxxxyzzz[k];

                g_0_y_yyy_xxxyyyy[k] = g_yy_xxxyyyy[k] - g_0_y_yy_xxxyyyy[k] * ab_y + g_0_y_yy_xxxyyyyy[k];

                g_0_y_yyy_xxxyyyz[k] = g_yy_xxxyyyz[k] - g_0_y_yy_xxxyyyz[k] * ab_y + g_0_y_yy_xxxyyyyz[k];

                g_0_y_yyy_xxxyyzz[k] = g_yy_xxxyyzz[k] - g_0_y_yy_xxxyyzz[k] * ab_y + g_0_y_yy_xxxyyyzz[k];

                g_0_y_yyy_xxxyzzz[k] = g_yy_xxxyzzz[k] - g_0_y_yy_xxxyzzz[k] * ab_y + g_0_y_yy_xxxyyzzz[k];

                g_0_y_yyy_xxxzzzz[k] = g_yy_xxxzzzz[k] - g_0_y_yy_xxxzzzz[k] * ab_y + g_0_y_yy_xxxyzzzz[k];

                g_0_y_yyy_xxyyyyy[k] = g_yy_xxyyyyy[k] - g_0_y_yy_xxyyyyy[k] * ab_y + g_0_y_yy_xxyyyyyy[k];

                g_0_y_yyy_xxyyyyz[k] = g_yy_xxyyyyz[k] - g_0_y_yy_xxyyyyz[k] * ab_y + g_0_y_yy_xxyyyyyz[k];

                g_0_y_yyy_xxyyyzz[k] = g_yy_xxyyyzz[k] - g_0_y_yy_xxyyyzz[k] * ab_y + g_0_y_yy_xxyyyyzz[k];

                g_0_y_yyy_xxyyzzz[k] = g_yy_xxyyzzz[k] - g_0_y_yy_xxyyzzz[k] * ab_y + g_0_y_yy_xxyyyzzz[k];

                g_0_y_yyy_xxyzzzz[k] = g_yy_xxyzzzz[k] - g_0_y_yy_xxyzzzz[k] * ab_y + g_0_y_yy_xxyyzzzz[k];

                g_0_y_yyy_xxzzzzz[k] = g_yy_xxzzzzz[k] - g_0_y_yy_xxzzzzz[k] * ab_y + g_0_y_yy_xxyzzzzz[k];

                g_0_y_yyy_xyyyyyy[k] = g_yy_xyyyyyy[k] - g_0_y_yy_xyyyyyy[k] * ab_y + g_0_y_yy_xyyyyyyy[k];

                g_0_y_yyy_xyyyyyz[k] = g_yy_xyyyyyz[k] - g_0_y_yy_xyyyyyz[k] * ab_y + g_0_y_yy_xyyyyyyz[k];

                g_0_y_yyy_xyyyyzz[k] = g_yy_xyyyyzz[k] - g_0_y_yy_xyyyyzz[k] * ab_y + g_0_y_yy_xyyyyyzz[k];

                g_0_y_yyy_xyyyzzz[k] = g_yy_xyyyzzz[k] - g_0_y_yy_xyyyzzz[k] * ab_y + g_0_y_yy_xyyyyzzz[k];

                g_0_y_yyy_xyyzzzz[k] = g_yy_xyyzzzz[k] - g_0_y_yy_xyyzzzz[k] * ab_y + g_0_y_yy_xyyyzzzz[k];

                g_0_y_yyy_xyzzzzz[k] = g_yy_xyzzzzz[k] - g_0_y_yy_xyzzzzz[k] * ab_y + g_0_y_yy_xyyzzzzz[k];

                g_0_y_yyy_xzzzzzz[k] = g_yy_xzzzzzz[k] - g_0_y_yy_xzzzzzz[k] * ab_y + g_0_y_yy_xyzzzzzz[k];

                g_0_y_yyy_yyyyyyy[k] = g_yy_yyyyyyy[k] - g_0_y_yy_yyyyyyy[k] * ab_y + g_0_y_yy_yyyyyyyy[k];

                g_0_y_yyy_yyyyyyz[k] = g_yy_yyyyyyz[k] - g_0_y_yy_yyyyyyz[k] * ab_y + g_0_y_yy_yyyyyyyz[k];

                g_0_y_yyy_yyyyyzz[k] = g_yy_yyyyyzz[k] - g_0_y_yy_yyyyyzz[k] * ab_y + g_0_y_yy_yyyyyyzz[k];

                g_0_y_yyy_yyyyzzz[k] = g_yy_yyyyzzz[k] - g_0_y_yy_yyyyzzz[k] * ab_y + g_0_y_yy_yyyyyzzz[k];

                g_0_y_yyy_yyyzzzz[k] = g_yy_yyyzzzz[k] - g_0_y_yy_yyyzzzz[k] * ab_y + g_0_y_yy_yyyyzzzz[k];

                g_0_y_yyy_yyzzzzz[k] = g_yy_yyzzzzz[k] - g_0_y_yy_yyzzzzz[k] * ab_y + g_0_y_yy_yyyzzzzz[k];

                g_0_y_yyy_yzzzzzz[k] = g_yy_yzzzzzz[k] - g_0_y_yy_yzzzzzz[k] * ab_y + g_0_y_yy_yyzzzzzz[k];

                g_0_y_yyy_zzzzzzz[k] = g_yy_zzzzzzz[k] - g_0_y_yy_zzzzzzz[k] * ab_y + g_0_y_yy_yzzzzzzz[k];
            }

            /// Set up 612-648 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yy_xxxxxxx, g_0_y_yy_xxxxxxxz, g_0_y_yy_xxxxxxy, g_0_y_yy_xxxxxxyz, g_0_y_yy_xxxxxxz, g_0_y_yy_xxxxxxzz, g_0_y_yy_xxxxxyy, g_0_y_yy_xxxxxyyz, g_0_y_yy_xxxxxyz, g_0_y_yy_xxxxxyzz, g_0_y_yy_xxxxxzz, g_0_y_yy_xxxxxzzz, g_0_y_yy_xxxxyyy, g_0_y_yy_xxxxyyyz, g_0_y_yy_xxxxyyz, g_0_y_yy_xxxxyyzz, g_0_y_yy_xxxxyzz, g_0_y_yy_xxxxyzzz, g_0_y_yy_xxxxzzz, g_0_y_yy_xxxxzzzz, g_0_y_yy_xxxyyyy, g_0_y_yy_xxxyyyyz, g_0_y_yy_xxxyyyz, g_0_y_yy_xxxyyyzz, g_0_y_yy_xxxyyzz, g_0_y_yy_xxxyyzzz, g_0_y_yy_xxxyzzz, g_0_y_yy_xxxyzzzz, g_0_y_yy_xxxzzzz, g_0_y_yy_xxxzzzzz, g_0_y_yy_xxyyyyy, g_0_y_yy_xxyyyyyz, g_0_y_yy_xxyyyyz, g_0_y_yy_xxyyyyzz, g_0_y_yy_xxyyyzz, g_0_y_yy_xxyyyzzz, g_0_y_yy_xxyyzzz, g_0_y_yy_xxyyzzzz, g_0_y_yy_xxyzzzz, g_0_y_yy_xxyzzzzz, g_0_y_yy_xxzzzzz, g_0_y_yy_xxzzzzzz, g_0_y_yy_xyyyyyy, g_0_y_yy_xyyyyyyz, g_0_y_yy_xyyyyyz, g_0_y_yy_xyyyyyzz, g_0_y_yy_xyyyyzz, g_0_y_yy_xyyyyzzz, g_0_y_yy_xyyyzzz, g_0_y_yy_xyyyzzzz, g_0_y_yy_xyyzzzz, g_0_y_yy_xyyzzzzz, g_0_y_yy_xyzzzzz, g_0_y_yy_xyzzzzzz, g_0_y_yy_xzzzzzz, g_0_y_yy_xzzzzzzz, g_0_y_yy_yyyyyyy, g_0_y_yy_yyyyyyyz, g_0_y_yy_yyyyyyz, g_0_y_yy_yyyyyyzz, g_0_y_yy_yyyyyzz, g_0_y_yy_yyyyyzzz, g_0_y_yy_yyyyzzz, g_0_y_yy_yyyyzzzz, g_0_y_yy_yyyzzzz, g_0_y_yy_yyyzzzzz, g_0_y_yy_yyzzzzz, g_0_y_yy_yyzzzzzz, g_0_y_yy_yzzzzzz, g_0_y_yy_yzzzzzzz, g_0_y_yy_zzzzzzz, g_0_y_yy_zzzzzzzz, g_0_y_yyz_xxxxxxx, g_0_y_yyz_xxxxxxy, g_0_y_yyz_xxxxxxz, g_0_y_yyz_xxxxxyy, g_0_y_yyz_xxxxxyz, g_0_y_yyz_xxxxxzz, g_0_y_yyz_xxxxyyy, g_0_y_yyz_xxxxyyz, g_0_y_yyz_xxxxyzz, g_0_y_yyz_xxxxzzz, g_0_y_yyz_xxxyyyy, g_0_y_yyz_xxxyyyz, g_0_y_yyz_xxxyyzz, g_0_y_yyz_xxxyzzz, g_0_y_yyz_xxxzzzz, g_0_y_yyz_xxyyyyy, g_0_y_yyz_xxyyyyz, g_0_y_yyz_xxyyyzz, g_0_y_yyz_xxyyzzz, g_0_y_yyz_xxyzzzz, g_0_y_yyz_xxzzzzz, g_0_y_yyz_xyyyyyy, g_0_y_yyz_xyyyyyz, g_0_y_yyz_xyyyyzz, g_0_y_yyz_xyyyzzz, g_0_y_yyz_xyyzzzz, g_0_y_yyz_xyzzzzz, g_0_y_yyz_xzzzzzz, g_0_y_yyz_yyyyyyy, g_0_y_yyz_yyyyyyz, g_0_y_yyz_yyyyyzz, g_0_y_yyz_yyyyzzz, g_0_y_yyz_yyyzzzz, g_0_y_yyz_yyzzzzz, g_0_y_yyz_yzzzzzz, g_0_y_yyz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyz_xxxxxxx[k] = -g_0_y_yy_xxxxxxx[k] * ab_z + g_0_y_yy_xxxxxxxz[k];

                g_0_y_yyz_xxxxxxy[k] = -g_0_y_yy_xxxxxxy[k] * ab_z + g_0_y_yy_xxxxxxyz[k];

                g_0_y_yyz_xxxxxxz[k] = -g_0_y_yy_xxxxxxz[k] * ab_z + g_0_y_yy_xxxxxxzz[k];

                g_0_y_yyz_xxxxxyy[k] = -g_0_y_yy_xxxxxyy[k] * ab_z + g_0_y_yy_xxxxxyyz[k];

                g_0_y_yyz_xxxxxyz[k] = -g_0_y_yy_xxxxxyz[k] * ab_z + g_0_y_yy_xxxxxyzz[k];

                g_0_y_yyz_xxxxxzz[k] = -g_0_y_yy_xxxxxzz[k] * ab_z + g_0_y_yy_xxxxxzzz[k];

                g_0_y_yyz_xxxxyyy[k] = -g_0_y_yy_xxxxyyy[k] * ab_z + g_0_y_yy_xxxxyyyz[k];

                g_0_y_yyz_xxxxyyz[k] = -g_0_y_yy_xxxxyyz[k] * ab_z + g_0_y_yy_xxxxyyzz[k];

                g_0_y_yyz_xxxxyzz[k] = -g_0_y_yy_xxxxyzz[k] * ab_z + g_0_y_yy_xxxxyzzz[k];

                g_0_y_yyz_xxxxzzz[k] = -g_0_y_yy_xxxxzzz[k] * ab_z + g_0_y_yy_xxxxzzzz[k];

                g_0_y_yyz_xxxyyyy[k] = -g_0_y_yy_xxxyyyy[k] * ab_z + g_0_y_yy_xxxyyyyz[k];

                g_0_y_yyz_xxxyyyz[k] = -g_0_y_yy_xxxyyyz[k] * ab_z + g_0_y_yy_xxxyyyzz[k];

                g_0_y_yyz_xxxyyzz[k] = -g_0_y_yy_xxxyyzz[k] * ab_z + g_0_y_yy_xxxyyzzz[k];

                g_0_y_yyz_xxxyzzz[k] = -g_0_y_yy_xxxyzzz[k] * ab_z + g_0_y_yy_xxxyzzzz[k];

                g_0_y_yyz_xxxzzzz[k] = -g_0_y_yy_xxxzzzz[k] * ab_z + g_0_y_yy_xxxzzzzz[k];

                g_0_y_yyz_xxyyyyy[k] = -g_0_y_yy_xxyyyyy[k] * ab_z + g_0_y_yy_xxyyyyyz[k];

                g_0_y_yyz_xxyyyyz[k] = -g_0_y_yy_xxyyyyz[k] * ab_z + g_0_y_yy_xxyyyyzz[k];

                g_0_y_yyz_xxyyyzz[k] = -g_0_y_yy_xxyyyzz[k] * ab_z + g_0_y_yy_xxyyyzzz[k];

                g_0_y_yyz_xxyyzzz[k] = -g_0_y_yy_xxyyzzz[k] * ab_z + g_0_y_yy_xxyyzzzz[k];

                g_0_y_yyz_xxyzzzz[k] = -g_0_y_yy_xxyzzzz[k] * ab_z + g_0_y_yy_xxyzzzzz[k];

                g_0_y_yyz_xxzzzzz[k] = -g_0_y_yy_xxzzzzz[k] * ab_z + g_0_y_yy_xxzzzzzz[k];

                g_0_y_yyz_xyyyyyy[k] = -g_0_y_yy_xyyyyyy[k] * ab_z + g_0_y_yy_xyyyyyyz[k];

                g_0_y_yyz_xyyyyyz[k] = -g_0_y_yy_xyyyyyz[k] * ab_z + g_0_y_yy_xyyyyyzz[k];

                g_0_y_yyz_xyyyyzz[k] = -g_0_y_yy_xyyyyzz[k] * ab_z + g_0_y_yy_xyyyyzzz[k];

                g_0_y_yyz_xyyyzzz[k] = -g_0_y_yy_xyyyzzz[k] * ab_z + g_0_y_yy_xyyyzzzz[k];

                g_0_y_yyz_xyyzzzz[k] = -g_0_y_yy_xyyzzzz[k] * ab_z + g_0_y_yy_xyyzzzzz[k];

                g_0_y_yyz_xyzzzzz[k] = -g_0_y_yy_xyzzzzz[k] * ab_z + g_0_y_yy_xyzzzzzz[k];

                g_0_y_yyz_xzzzzzz[k] = -g_0_y_yy_xzzzzzz[k] * ab_z + g_0_y_yy_xzzzzzzz[k];

                g_0_y_yyz_yyyyyyy[k] = -g_0_y_yy_yyyyyyy[k] * ab_z + g_0_y_yy_yyyyyyyz[k];

                g_0_y_yyz_yyyyyyz[k] = -g_0_y_yy_yyyyyyz[k] * ab_z + g_0_y_yy_yyyyyyzz[k];

                g_0_y_yyz_yyyyyzz[k] = -g_0_y_yy_yyyyyzz[k] * ab_z + g_0_y_yy_yyyyyzzz[k];

                g_0_y_yyz_yyyyzzz[k] = -g_0_y_yy_yyyyzzz[k] * ab_z + g_0_y_yy_yyyyzzzz[k];

                g_0_y_yyz_yyyzzzz[k] = -g_0_y_yy_yyyzzzz[k] * ab_z + g_0_y_yy_yyyzzzzz[k];

                g_0_y_yyz_yyzzzzz[k] = -g_0_y_yy_yyzzzzz[k] * ab_z + g_0_y_yy_yyzzzzzz[k];

                g_0_y_yyz_yzzzzzz[k] = -g_0_y_yy_yzzzzzz[k] * ab_z + g_0_y_yy_yzzzzzzz[k];

                g_0_y_yyz_zzzzzzz[k] = -g_0_y_yy_zzzzzzz[k] * ab_z + g_0_y_yy_zzzzzzzz[k];
            }

            /// Set up 648-684 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yz_xxxxxxx, g_0_y_yz_xxxxxxxz, g_0_y_yz_xxxxxxy, g_0_y_yz_xxxxxxyz, g_0_y_yz_xxxxxxz, g_0_y_yz_xxxxxxzz, g_0_y_yz_xxxxxyy, g_0_y_yz_xxxxxyyz, g_0_y_yz_xxxxxyz, g_0_y_yz_xxxxxyzz, g_0_y_yz_xxxxxzz, g_0_y_yz_xxxxxzzz, g_0_y_yz_xxxxyyy, g_0_y_yz_xxxxyyyz, g_0_y_yz_xxxxyyz, g_0_y_yz_xxxxyyzz, g_0_y_yz_xxxxyzz, g_0_y_yz_xxxxyzzz, g_0_y_yz_xxxxzzz, g_0_y_yz_xxxxzzzz, g_0_y_yz_xxxyyyy, g_0_y_yz_xxxyyyyz, g_0_y_yz_xxxyyyz, g_0_y_yz_xxxyyyzz, g_0_y_yz_xxxyyzz, g_0_y_yz_xxxyyzzz, g_0_y_yz_xxxyzzz, g_0_y_yz_xxxyzzzz, g_0_y_yz_xxxzzzz, g_0_y_yz_xxxzzzzz, g_0_y_yz_xxyyyyy, g_0_y_yz_xxyyyyyz, g_0_y_yz_xxyyyyz, g_0_y_yz_xxyyyyzz, g_0_y_yz_xxyyyzz, g_0_y_yz_xxyyyzzz, g_0_y_yz_xxyyzzz, g_0_y_yz_xxyyzzzz, g_0_y_yz_xxyzzzz, g_0_y_yz_xxyzzzzz, g_0_y_yz_xxzzzzz, g_0_y_yz_xxzzzzzz, g_0_y_yz_xyyyyyy, g_0_y_yz_xyyyyyyz, g_0_y_yz_xyyyyyz, g_0_y_yz_xyyyyyzz, g_0_y_yz_xyyyyzz, g_0_y_yz_xyyyyzzz, g_0_y_yz_xyyyzzz, g_0_y_yz_xyyyzzzz, g_0_y_yz_xyyzzzz, g_0_y_yz_xyyzzzzz, g_0_y_yz_xyzzzzz, g_0_y_yz_xyzzzzzz, g_0_y_yz_xzzzzzz, g_0_y_yz_xzzzzzzz, g_0_y_yz_yyyyyyy, g_0_y_yz_yyyyyyyz, g_0_y_yz_yyyyyyz, g_0_y_yz_yyyyyyzz, g_0_y_yz_yyyyyzz, g_0_y_yz_yyyyyzzz, g_0_y_yz_yyyyzzz, g_0_y_yz_yyyyzzzz, g_0_y_yz_yyyzzzz, g_0_y_yz_yyyzzzzz, g_0_y_yz_yyzzzzz, g_0_y_yz_yyzzzzzz, g_0_y_yz_yzzzzzz, g_0_y_yz_yzzzzzzz, g_0_y_yz_zzzzzzz, g_0_y_yz_zzzzzzzz, g_0_y_yzz_xxxxxxx, g_0_y_yzz_xxxxxxy, g_0_y_yzz_xxxxxxz, g_0_y_yzz_xxxxxyy, g_0_y_yzz_xxxxxyz, g_0_y_yzz_xxxxxzz, g_0_y_yzz_xxxxyyy, g_0_y_yzz_xxxxyyz, g_0_y_yzz_xxxxyzz, g_0_y_yzz_xxxxzzz, g_0_y_yzz_xxxyyyy, g_0_y_yzz_xxxyyyz, g_0_y_yzz_xxxyyzz, g_0_y_yzz_xxxyzzz, g_0_y_yzz_xxxzzzz, g_0_y_yzz_xxyyyyy, g_0_y_yzz_xxyyyyz, g_0_y_yzz_xxyyyzz, g_0_y_yzz_xxyyzzz, g_0_y_yzz_xxyzzzz, g_0_y_yzz_xxzzzzz, g_0_y_yzz_xyyyyyy, g_0_y_yzz_xyyyyyz, g_0_y_yzz_xyyyyzz, g_0_y_yzz_xyyyzzz, g_0_y_yzz_xyyzzzz, g_0_y_yzz_xyzzzzz, g_0_y_yzz_xzzzzzz, g_0_y_yzz_yyyyyyy, g_0_y_yzz_yyyyyyz, g_0_y_yzz_yyyyyzz, g_0_y_yzz_yyyyzzz, g_0_y_yzz_yyyzzzz, g_0_y_yzz_yyzzzzz, g_0_y_yzz_yzzzzzz, g_0_y_yzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzz_xxxxxxx[k] = -g_0_y_yz_xxxxxxx[k] * ab_z + g_0_y_yz_xxxxxxxz[k];

                g_0_y_yzz_xxxxxxy[k] = -g_0_y_yz_xxxxxxy[k] * ab_z + g_0_y_yz_xxxxxxyz[k];

                g_0_y_yzz_xxxxxxz[k] = -g_0_y_yz_xxxxxxz[k] * ab_z + g_0_y_yz_xxxxxxzz[k];

                g_0_y_yzz_xxxxxyy[k] = -g_0_y_yz_xxxxxyy[k] * ab_z + g_0_y_yz_xxxxxyyz[k];

                g_0_y_yzz_xxxxxyz[k] = -g_0_y_yz_xxxxxyz[k] * ab_z + g_0_y_yz_xxxxxyzz[k];

                g_0_y_yzz_xxxxxzz[k] = -g_0_y_yz_xxxxxzz[k] * ab_z + g_0_y_yz_xxxxxzzz[k];

                g_0_y_yzz_xxxxyyy[k] = -g_0_y_yz_xxxxyyy[k] * ab_z + g_0_y_yz_xxxxyyyz[k];

                g_0_y_yzz_xxxxyyz[k] = -g_0_y_yz_xxxxyyz[k] * ab_z + g_0_y_yz_xxxxyyzz[k];

                g_0_y_yzz_xxxxyzz[k] = -g_0_y_yz_xxxxyzz[k] * ab_z + g_0_y_yz_xxxxyzzz[k];

                g_0_y_yzz_xxxxzzz[k] = -g_0_y_yz_xxxxzzz[k] * ab_z + g_0_y_yz_xxxxzzzz[k];

                g_0_y_yzz_xxxyyyy[k] = -g_0_y_yz_xxxyyyy[k] * ab_z + g_0_y_yz_xxxyyyyz[k];

                g_0_y_yzz_xxxyyyz[k] = -g_0_y_yz_xxxyyyz[k] * ab_z + g_0_y_yz_xxxyyyzz[k];

                g_0_y_yzz_xxxyyzz[k] = -g_0_y_yz_xxxyyzz[k] * ab_z + g_0_y_yz_xxxyyzzz[k];

                g_0_y_yzz_xxxyzzz[k] = -g_0_y_yz_xxxyzzz[k] * ab_z + g_0_y_yz_xxxyzzzz[k];

                g_0_y_yzz_xxxzzzz[k] = -g_0_y_yz_xxxzzzz[k] * ab_z + g_0_y_yz_xxxzzzzz[k];

                g_0_y_yzz_xxyyyyy[k] = -g_0_y_yz_xxyyyyy[k] * ab_z + g_0_y_yz_xxyyyyyz[k];

                g_0_y_yzz_xxyyyyz[k] = -g_0_y_yz_xxyyyyz[k] * ab_z + g_0_y_yz_xxyyyyzz[k];

                g_0_y_yzz_xxyyyzz[k] = -g_0_y_yz_xxyyyzz[k] * ab_z + g_0_y_yz_xxyyyzzz[k];

                g_0_y_yzz_xxyyzzz[k] = -g_0_y_yz_xxyyzzz[k] * ab_z + g_0_y_yz_xxyyzzzz[k];

                g_0_y_yzz_xxyzzzz[k] = -g_0_y_yz_xxyzzzz[k] * ab_z + g_0_y_yz_xxyzzzzz[k];

                g_0_y_yzz_xxzzzzz[k] = -g_0_y_yz_xxzzzzz[k] * ab_z + g_0_y_yz_xxzzzzzz[k];

                g_0_y_yzz_xyyyyyy[k] = -g_0_y_yz_xyyyyyy[k] * ab_z + g_0_y_yz_xyyyyyyz[k];

                g_0_y_yzz_xyyyyyz[k] = -g_0_y_yz_xyyyyyz[k] * ab_z + g_0_y_yz_xyyyyyzz[k];

                g_0_y_yzz_xyyyyzz[k] = -g_0_y_yz_xyyyyzz[k] * ab_z + g_0_y_yz_xyyyyzzz[k];

                g_0_y_yzz_xyyyzzz[k] = -g_0_y_yz_xyyyzzz[k] * ab_z + g_0_y_yz_xyyyzzzz[k];

                g_0_y_yzz_xyyzzzz[k] = -g_0_y_yz_xyyzzzz[k] * ab_z + g_0_y_yz_xyyzzzzz[k];

                g_0_y_yzz_xyzzzzz[k] = -g_0_y_yz_xyzzzzz[k] * ab_z + g_0_y_yz_xyzzzzzz[k];

                g_0_y_yzz_xzzzzzz[k] = -g_0_y_yz_xzzzzzz[k] * ab_z + g_0_y_yz_xzzzzzzz[k];

                g_0_y_yzz_yyyyyyy[k] = -g_0_y_yz_yyyyyyy[k] * ab_z + g_0_y_yz_yyyyyyyz[k];

                g_0_y_yzz_yyyyyyz[k] = -g_0_y_yz_yyyyyyz[k] * ab_z + g_0_y_yz_yyyyyyzz[k];

                g_0_y_yzz_yyyyyzz[k] = -g_0_y_yz_yyyyyzz[k] * ab_z + g_0_y_yz_yyyyyzzz[k];

                g_0_y_yzz_yyyyzzz[k] = -g_0_y_yz_yyyyzzz[k] * ab_z + g_0_y_yz_yyyyzzzz[k];

                g_0_y_yzz_yyyzzzz[k] = -g_0_y_yz_yyyzzzz[k] * ab_z + g_0_y_yz_yyyzzzzz[k];

                g_0_y_yzz_yyzzzzz[k] = -g_0_y_yz_yyzzzzz[k] * ab_z + g_0_y_yz_yyzzzzzz[k];

                g_0_y_yzz_yzzzzzz[k] = -g_0_y_yz_yzzzzzz[k] * ab_z + g_0_y_yz_yzzzzzzz[k];

                g_0_y_yzz_zzzzzzz[k] = -g_0_y_yz_zzzzzzz[k] * ab_z + g_0_y_yz_zzzzzzzz[k];
            }

            /// Set up 684-720 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_zz_xxxxxxx, g_0_y_zz_xxxxxxxz, g_0_y_zz_xxxxxxy, g_0_y_zz_xxxxxxyz, g_0_y_zz_xxxxxxz, g_0_y_zz_xxxxxxzz, g_0_y_zz_xxxxxyy, g_0_y_zz_xxxxxyyz, g_0_y_zz_xxxxxyz, g_0_y_zz_xxxxxyzz, g_0_y_zz_xxxxxzz, g_0_y_zz_xxxxxzzz, g_0_y_zz_xxxxyyy, g_0_y_zz_xxxxyyyz, g_0_y_zz_xxxxyyz, g_0_y_zz_xxxxyyzz, g_0_y_zz_xxxxyzz, g_0_y_zz_xxxxyzzz, g_0_y_zz_xxxxzzz, g_0_y_zz_xxxxzzzz, g_0_y_zz_xxxyyyy, g_0_y_zz_xxxyyyyz, g_0_y_zz_xxxyyyz, g_0_y_zz_xxxyyyzz, g_0_y_zz_xxxyyzz, g_0_y_zz_xxxyyzzz, g_0_y_zz_xxxyzzz, g_0_y_zz_xxxyzzzz, g_0_y_zz_xxxzzzz, g_0_y_zz_xxxzzzzz, g_0_y_zz_xxyyyyy, g_0_y_zz_xxyyyyyz, g_0_y_zz_xxyyyyz, g_0_y_zz_xxyyyyzz, g_0_y_zz_xxyyyzz, g_0_y_zz_xxyyyzzz, g_0_y_zz_xxyyzzz, g_0_y_zz_xxyyzzzz, g_0_y_zz_xxyzzzz, g_0_y_zz_xxyzzzzz, g_0_y_zz_xxzzzzz, g_0_y_zz_xxzzzzzz, g_0_y_zz_xyyyyyy, g_0_y_zz_xyyyyyyz, g_0_y_zz_xyyyyyz, g_0_y_zz_xyyyyyzz, g_0_y_zz_xyyyyzz, g_0_y_zz_xyyyyzzz, g_0_y_zz_xyyyzzz, g_0_y_zz_xyyyzzzz, g_0_y_zz_xyyzzzz, g_0_y_zz_xyyzzzzz, g_0_y_zz_xyzzzzz, g_0_y_zz_xyzzzzzz, g_0_y_zz_xzzzzzz, g_0_y_zz_xzzzzzzz, g_0_y_zz_yyyyyyy, g_0_y_zz_yyyyyyyz, g_0_y_zz_yyyyyyz, g_0_y_zz_yyyyyyzz, g_0_y_zz_yyyyyzz, g_0_y_zz_yyyyyzzz, g_0_y_zz_yyyyzzz, g_0_y_zz_yyyyzzzz, g_0_y_zz_yyyzzzz, g_0_y_zz_yyyzzzzz, g_0_y_zz_yyzzzzz, g_0_y_zz_yyzzzzzz, g_0_y_zz_yzzzzzz, g_0_y_zz_yzzzzzzz, g_0_y_zz_zzzzzzz, g_0_y_zz_zzzzzzzz, g_0_y_zzz_xxxxxxx, g_0_y_zzz_xxxxxxy, g_0_y_zzz_xxxxxxz, g_0_y_zzz_xxxxxyy, g_0_y_zzz_xxxxxyz, g_0_y_zzz_xxxxxzz, g_0_y_zzz_xxxxyyy, g_0_y_zzz_xxxxyyz, g_0_y_zzz_xxxxyzz, g_0_y_zzz_xxxxzzz, g_0_y_zzz_xxxyyyy, g_0_y_zzz_xxxyyyz, g_0_y_zzz_xxxyyzz, g_0_y_zzz_xxxyzzz, g_0_y_zzz_xxxzzzz, g_0_y_zzz_xxyyyyy, g_0_y_zzz_xxyyyyz, g_0_y_zzz_xxyyyzz, g_0_y_zzz_xxyyzzz, g_0_y_zzz_xxyzzzz, g_0_y_zzz_xxzzzzz, g_0_y_zzz_xyyyyyy, g_0_y_zzz_xyyyyyz, g_0_y_zzz_xyyyyzz, g_0_y_zzz_xyyyzzz, g_0_y_zzz_xyyzzzz, g_0_y_zzz_xyzzzzz, g_0_y_zzz_xzzzzzz, g_0_y_zzz_yyyyyyy, g_0_y_zzz_yyyyyyz, g_0_y_zzz_yyyyyzz, g_0_y_zzz_yyyyzzz, g_0_y_zzz_yyyzzzz, g_0_y_zzz_yyzzzzz, g_0_y_zzz_yzzzzzz, g_0_y_zzz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzz_xxxxxxx[k] = -g_0_y_zz_xxxxxxx[k] * ab_z + g_0_y_zz_xxxxxxxz[k];

                g_0_y_zzz_xxxxxxy[k] = -g_0_y_zz_xxxxxxy[k] * ab_z + g_0_y_zz_xxxxxxyz[k];

                g_0_y_zzz_xxxxxxz[k] = -g_0_y_zz_xxxxxxz[k] * ab_z + g_0_y_zz_xxxxxxzz[k];

                g_0_y_zzz_xxxxxyy[k] = -g_0_y_zz_xxxxxyy[k] * ab_z + g_0_y_zz_xxxxxyyz[k];

                g_0_y_zzz_xxxxxyz[k] = -g_0_y_zz_xxxxxyz[k] * ab_z + g_0_y_zz_xxxxxyzz[k];

                g_0_y_zzz_xxxxxzz[k] = -g_0_y_zz_xxxxxzz[k] * ab_z + g_0_y_zz_xxxxxzzz[k];

                g_0_y_zzz_xxxxyyy[k] = -g_0_y_zz_xxxxyyy[k] * ab_z + g_0_y_zz_xxxxyyyz[k];

                g_0_y_zzz_xxxxyyz[k] = -g_0_y_zz_xxxxyyz[k] * ab_z + g_0_y_zz_xxxxyyzz[k];

                g_0_y_zzz_xxxxyzz[k] = -g_0_y_zz_xxxxyzz[k] * ab_z + g_0_y_zz_xxxxyzzz[k];

                g_0_y_zzz_xxxxzzz[k] = -g_0_y_zz_xxxxzzz[k] * ab_z + g_0_y_zz_xxxxzzzz[k];

                g_0_y_zzz_xxxyyyy[k] = -g_0_y_zz_xxxyyyy[k] * ab_z + g_0_y_zz_xxxyyyyz[k];

                g_0_y_zzz_xxxyyyz[k] = -g_0_y_zz_xxxyyyz[k] * ab_z + g_0_y_zz_xxxyyyzz[k];

                g_0_y_zzz_xxxyyzz[k] = -g_0_y_zz_xxxyyzz[k] * ab_z + g_0_y_zz_xxxyyzzz[k];

                g_0_y_zzz_xxxyzzz[k] = -g_0_y_zz_xxxyzzz[k] * ab_z + g_0_y_zz_xxxyzzzz[k];

                g_0_y_zzz_xxxzzzz[k] = -g_0_y_zz_xxxzzzz[k] * ab_z + g_0_y_zz_xxxzzzzz[k];

                g_0_y_zzz_xxyyyyy[k] = -g_0_y_zz_xxyyyyy[k] * ab_z + g_0_y_zz_xxyyyyyz[k];

                g_0_y_zzz_xxyyyyz[k] = -g_0_y_zz_xxyyyyz[k] * ab_z + g_0_y_zz_xxyyyyzz[k];

                g_0_y_zzz_xxyyyzz[k] = -g_0_y_zz_xxyyyzz[k] * ab_z + g_0_y_zz_xxyyyzzz[k];

                g_0_y_zzz_xxyyzzz[k] = -g_0_y_zz_xxyyzzz[k] * ab_z + g_0_y_zz_xxyyzzzz[k];

                g_0_y_zzz_xxyzzzz[k] = -g_0_y_zz_xxyzzzz[k] * ab_z + g_0_y_zz_xxyzzzzz[k];

                g_0_y_zzz_xxzzzzz[k] = -g_0_y_zz_xxzzzzz[k] * ab_z + g_0_y_zz_xxzzzzzz[k];

                g_0_y_zzz_xyyyyyy[k] = -g_0_y_zz_xyyyyyy[k] * ab_z + g_0_y_zz_xyyyyyyz[k];

                g_0_y_zzz_xyyyyyz[k] = -g_0_y_zz_xyyyyyz[k] * ab_z + g_0_y_zz_xyyyyyzz[k];

                g_0_y_zzz_xyyyyzz[k] = -g_0_y_zz_xyyyyzz[k] * ab_z + g_0_y_zz_xyyyyzzz[k];

                g_0_y_zzz_xyyyzzz[k] = -g_0_y_zz_xyyyzzz[k] * ab_z + g_0_y_zz_xyyyzzzz[k];

                g_0_y_zzz_xyyzzzz[k] = -g_0_y_zz_xyyzzzz[k] * ab_z + g_0_y_zz_xyyzzzzz[k];

                g_0_y_zzz_xyzzzzz[k] = -g_0_y_zz_xyzzzzz[k] * ab_z + g_0_y_zz_xyzzzzzz[k];

                g_0_y_zzz_xzzzzzz[k] = -g_0_y_zz_xzzzzzz[k] * ab_z + g_0_y_zz_xzzzzzzz[k];

                g_0_y_zzz_yyyyyyy[k] = -g_0_y_zz_yyyyyyy[k] * ab_z + g_0_y_zz_yyyyyyyz[k];

                g_0_y_zzz_yyyyyyz[k] = -g_0_y_zz_yyyyyyz[k] * ab_z + g_0_y_zz_yyyyyyzz[k];

                g_0_y_zzz_yyyyyzz[k] = -g_0_y_zz_yyyyyzz[k] * ab_z + g_0_y_zz_yyyyyzzz[k];

                g_0_y_zzz_yyyyzzz[k] = -g_0_y_zz_yyyyzzz[k] * ab_z + g_0_y_zz_yyyyzzzz[k];

                g_0_y_zzz_yyyzzzz[k] = -g_0_y_zz_yyyzzzz[k] * ab_z + g_0_y_zz_yyyzzzzz[k];

                g_0_y_zzz_yyzzzzz[k] = -g_0_y_zz_yyzzzzz[k] * ab_z + g_0_y_zz_yyzzzzzz[k];

                g_0_y_zzz_yzzzzzz[k] = -g_0_y_zz_yzzzzzz[k] * ab_z + g_0_y_zz_yzzzzzzz[k];

                g_0_y_zzz_zzzzzzz[k] = -g_0_y_zz_zzzzzzz[k] * ab_z + g_0_y_zz_zzzzzzzz[k];
            }

            /// Set up 720-756 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xx_xxxxxxx, g_0_z_xx_xxxxxxxx, g_0_z_xx_xxxxxxxy, g_0_z_xx_xxxxxxxz, g_0_z_xx_xxxxxxy, g_0_z_xx_xxxxxxyy, g_0_z_xx_xxxxxxyz, g_0_z_xx_xxxxxxz, g_0_z_xx_xxxxxxzz, g_0_z_xx_xxxxxyy, g_0_z_xx_xxxxxyyy, g_0_z_xx_xxxxxyyz, g_0_z_xx_xxxxxyz, g_0_z_xx_xxxxxyzz, g_0_z_xx_xxxxxzz, g_0_z_xx_xxxxxzzz, g_0_z_xx_xxxxyyy, g_0_z_xx_xxxxyyyy, g_0_z_xx_xxxxyyyz, g_0_z_xx_xxxxyyz, g_0_z_xx_xxxxyyzz, g_0_z_xx_xxxxyzz, g_0_z_xx_xxxxyzzz, g_0_z_xx_xxxxzzz, g_0_z_xx_xxxxzzzz, g_0_z_xx_xxxyyyy, g_0_z_xx_xxxyyyyy, g_0_z_xx_xxxyyyyz, g_0_z_xx_xxxyyyz, g_0_z_xx_xxxyyyzz, g_0_z_xx_xxxyyzz, g_0_z_xx_xxxyyzzz, g_0_z_xx_xxxyzzz, g_0_z_xx_xxxyzzzz, g_0_z_xx_xxxzzzz, g_0_z_xx_xxxzzzzz, g_0_z_xx_xxyyyyy, g_0_z_xx_xxyyyyyy, g_0_z_xx_xxyyyyyz, g_0_z_xx_xxyyyyz, g_0_z_xx_xxyyyyzz, g_0_z_xx_xxyyyzz, g_0_z_xx_xxyyyzzz, g_0_z_xx_xxyyzzz, g_0_z_xx_xxyyzzzz, g_0_z_xx_xxyzzzz, g_0_z_xx_xxyzzzzz, g_0_z_xx_xxzzzzz, g_0_z_xx_xxzzzzzz, g_0_z_xx_xyyyyyy, g_0_z_xx_xyyyyyyy, g_0_z_xx_xyyyyyyz, g_0_z_xx_xyyyyyz, g_0_z_xx_xyyyyyzz, g_0_z_xx_xyyyyzz, g_0_z_xx_xyyyyzzz, g_0_z_xx_xyyyzzz, g_0_z_xx_xyyyzzzz, g_0_z_xx_xyyzzzz, g_0_z_xx_xyyzzzzz, g_0_z_xx_xyzzzzz, g_0_z_xx_xyzzzzzz, g_0_z_xx_xzzzzzz, g_0_z_xx_xzzzzzzz, g_0_z_xx_yyyyyyy, g_0_z_xx_yyyyyyz, g_0_z_xx_yyyyyzz, g_0_z_xx_yyyyzzz, g_0_z_xx_yyyzzzz, g_0_z_xx_yyzzzzz, g_0_z_xx_yzzzzzz, g_0_z_xx_zzzzzzz, g_0_z_xxx_xxxxxxx, g_0_z_xxx_xxxxxxy, g_0_z_xxx_xxxxxxz, g_0_z_xxx_xxxxxyy, g_0_z_xxx_xxxxxyz, g_0_z_xxx_xxxxxzz, g_0_z_xxx_xxxxyyy, g_0_z_xxx_xxxxyyz, g_0_z_xxx_xxxxyzz, g_0_z_xxx_xxxxzzz, g_0_z_xxx_xxxyyyy, g_0_z_xxx_xxxyyyz, g_0_z_xxx_xxxyyzz, g_0_z_xxx_xxxyzzz, g_0_z_xxx_xxxzzzz, g_0_z_xxx_xxyyyyy, g_0_z_xxx_xxyyyyz, g_0_z_xxx_xxyyyzz, g_0_z_xxx_xxyyzzz, g_0_z_xxx_xxyzzzz, g_0_z_xxx_xxzzzzz, g_0_z_xxx_xyyyyyy, g_0_z_xxx_xyyyyyz, g_0_z_xxx_xyyyyzz, g_0_z_xxx_xyyyzzz, g_0_z_xxx_xyyzzzz, g_0_z_xxx_xyzzzzz, g_0_z_xxx_xzzzzzz, g_0_z_xxx_yyyyyyy, g_0_z_xxx_yyyyyyz, g_0_z_xxx_yyyyyzz, g_0_z_xxx_yyyyzzz, g_0_z_xxx_yyyzzzz, g_0_z_xxx_yyzzzzz, g_0_z_xxx_yzzzzzz, g_0_z_xxx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxx_xxxxxxx[k] = -g_0_z_xx_xxxxxxx[k] * ab_x + g_0_z_xx_xxxxxxxx[k];

                g_0_z_xxx_xxxxxxy[k] = -g_0_z_xx_xxxxxxy[k] * ab_x + g_0_z_xx_xxxxxxxy[k];

                g_0_z_xxx_xxxxxxz[k] = -g_0_z_xx_xxxxxxz[k] * ab_x + g_0_z_xx_xxxxxxxz[k];

                g_0_z_xxx_xxxxxyy[k] = -g_0_z_xx_xxxxxyy[k] * ab_x + g_0_z_xx_xxxxxxyy[k];

                g_0_z_xxx_xxxxxyz[k] = -g_0_z_xx_xxxxxyz[k] * ab_x + g_0_z_xx_xxxxxxyz[k];

                g_0_z_xxx_xxxxxzz[k] = -g_0_z_xx_xxxxxzz[k] * ab_x + g_0_z_xx_xxxxxxzz[k];

                g_0_z_xxx_xxxxyyy[k] = -g_0_z_xx_xxxxyyy[k] * ab_x + g_0_z_xx_xxxxxyyy[k];

                g_0_z_xxx_xxxxyyz[k] = -g_0_z_xx_xxxxyyz[k] * ab_x + g_0_z_xx_xxxxxyyz[k];

                g_0_z_xxx_xxxxyzz[k] = -g_0_z_xx_xxxxyzz[k] * ab_x + g_0_z_xx_xxxxxyzz[k];

                g_0_z_xxx_xxxxzzz[k] = -g_0_z_xx_xxxxzzz[k] * ab_x + g_0_z_xx_xxxxxzzz[k];

                g_0_z_xxx_xxxyyyy[k] = -g_0_z_xx_xxxyyyy[k] * ab_x + g_0_z_xx_xxxxyyyy[k];

                g_0_z_xxx_xxxyyyz[k] = -g_0_z_xx_xxxyyyz[k] * ab_x + g_0_z_xx_xxxxyyyz[k];

                g_0_z_xxx_xxxyyzz[k] = -g_0_z_xx_xxxyyzz[k] * ab_x + g_0_z_xx_xxxxyyzz[k];

                g_0_z_xxx_xxxyzzz[k] = -g_0_z_xx_xxxyzzz[k] * ab_x + g_0_z_xx_xxxxyzzz[k];

                g_0_z_xxx_xxxzzzz[k] = -g_0_z_xx_xxxzzzz[k] * ab_x + g_0_z_xx_xxxxzzzz[k];

                g_0_z_xxx_xxyyyyy[k] = -g_0_z_xx_xxyyyyy[k] * ab_x + g_0_z_xx_xxxyyyyy[k];

                g_0_z_xxx_xxyyyyz[k] = -g_0_z_xx_xxyyyyz[k] * ab_x + g_0_z_xx_xxxyyyyz[k];

                g_0_z_xxx_xxyyyzz[k] = -g_0_z_xx_xxyyyzz[k] * ab_x + g_0_z_xx_xxxyyyzz[k];

                g_0_z_xxx_xxyyzzz[k] = -g_0_z_xx_xxyyzzz[k] * ab_x + g_0_z_xx_xxxyyzzz[k];

                g_0_z_xxx_xxyzzzz[k] = -g_0_z_xx_xxyzzzz[k] * ab_x + g_0_z_xx_xxxyzzzz[k];

                g_0_z_xxx_xxzzzzz[k] = -g_0_z_xx_xxzzzzz[k] * ab_x + g_0_z_xx_xxxzzzzz[k];

                g_0_z_xxx_xyyyyyy[k] = -g_0_z_xx_xyyyyyy[k] * ab_x + g_0_z_xx_xxyyyyyy[k];

                g_0_z_xxx_xyyyyyz[k] = -g_0_z_xx_xyyyyyz[k] * ab_x + g_0_z_xx_xxyyyyyz[k];

                g_0_z_xxx_xyyyyzz[k] = -g_0_z_xx_xyyyyzz[k] * ab_x + g_0_z_xx_xxyyyyzz[k];

                g_0_z_xxx_xyyyzzz[k] = -g_0_z_xx_xyyyzzz[k] * ab_x + g_0_z_xx_xxyyyzzz[k];

                g_0_z_xxx_xyyzzzz[k] = -g_0_z_xx_xyyzzzz[k] * ab_x + g_0_z_xx_xxyyzzzz[k];

                g_0_z_xxx_xyzzzzz[k] = -g_0_z_xx_xyzzzzz[k] * ab_x + g_0_z_xx_xxyzzzzz[k];

                g_0_z_xxx_xzzzzzz[k] = -g_0_z_xx_xzzzzzz[k] * ab_x + g_0_z_xx_xxzzzzzz[k];

                g_0_z_xxx_yyyyyyy[k] = -g_0_z_xx_yyyyyyy[k] * ab_x + g_0_z_xx_xyyyyyyy[k];

                g_0_z_xxx_yyyyyyz[k] = -g_0_z_xx_yyyyyyz[k] * ab_x + g_0_z_xx_xyyyyyyz[k];

                g_0_z_xxx_yyyyyzz[k] = -g_0_z_xx_yyyyyzz[k] * ab_x + g_0_z_xx_xyyyyyzz[k];

                g_0_z_xxx_yyyyzzz[k] = -g_0_z_xx_yyyyzzz[k] * ab_x + g_0_z_xx_xyyyyzzz[k];

                g_0_z_xxx_yyyzzzz[k] = -g_0_z_xx_yyyzzzz[k] * ab_x + g_0_z_xx_xyyyzzzz[k];

                g_0_z_xxx_yyzzzzz[k] = -g_0_z_xx_yyzzzzz[k] * ab_x + g_0_z_xx_xyyzzzzz[k];

                g_0_z_xxx_yzzzzzz[k] = -g_0_z_xx_yzzzzzz[k] * ab_x + g_0_z_xx_xyzzzzzz[k];

                g_0_z_xxx_zzzzzzz[k] = -g_0_z_xx_zzzzzzz[k] * ab_x + g_0_z_xx_xzzzzzzz[k];
            }

            /// Set up 756-792 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxy_xxxxxxx, g_0_z_xxy_xxxxxxy, g_0_z_xxy_xxxxxxz, g_0_z_xxy_xxxxxyy, g_0_z_xxy_xxxxxyz, g_0_z_xxy_xxxxxzz, g_0_z_xxy_xxxxyyy, g_0_z_xxy_xxxxyyz, g_0_z_xxy_xxxxyzz, g_0_z_xxy_xxxxzzz, g_0_z_xxy_xxxyyyy, g_0_z_xxy_xxxyyyz, g_0_z_xxy_xxxyyzz, g_0_z_xxy_xxxyzzz, g_0_z_xxy_xxxzzzz, g_0_z_xxy_xxyyyyy, g_0_z_xxy_xxyyyyz, g_0_z_xxy_xxyyyzz, g_0_z_xxy_xxyyzzz, g_0_z_xxy_xxyzzzz, g_0_z_xxy_xxzzzzz, g_0_z_xxy_xyyyyyy, g_0_z_xxy_xyyyyyz, g_0_z_xxy_xyyyyzz, g_0_z_xxy_xyyyzzz, g_0_z_xxy_xyyzzzz, g_0_z_xxy_xyzzzzz, g_0_z_xxy_xzzzzzz, g_0_z_xxy_yyyyyyy, g_0_z_xxy_yyyyyyz, g_0_z_xxy_yyyyyzz, g_0_z_xxy_yyyyzzz, g_0_z_xxy_yyyzzzz, g_0_z_xxy_yyzzzzz, g_0_z_xxy_yzzzzzz, g_0_z_xxy_zzzzzzz, g_0_z_xy_xxxxxxx, g_0_z_xy_xxxxxxxx, g_0_z_xy_xxxxxxxy, g_0_z_xy_xxxxxxxz, g_0_z_xy_xxxxxxy, g_0_z_xy_xxxxxxyy, g_0_z_xy_xxxxxxyz, g_0_z_xy_xxxxxxz, g_0_z_xy_xxxxxxzz, g_0_z_xy_xxxxxyy, g_0_z_xy_xxxxxyyy, g_0_z_xy_xxxxxyyz, g_0_z_xy_xxxxxyz, g_0_z_xy_xxxxxyzz, g_0_z_xy_xxxxxzz, g_0_z_xy_xxxxxzzz, g_0_z_xy_xxxxyyy, g_0_z_xy_xxxxyyyy, g_0_z_xy_xxxxyyyz, g_0_z_xy_xxxxyyz, g_0_z_xy_xxxxyyzz, g_0_z_xy_xxxxyzz, g_0_z_xy_xxxxyzzz, g_0_z_xy_xxxxzzz, g_0_z_xy_xxxxzzzz, g_0_z_xy_xxxyyyy, g_0_z_xy_xxxyyyyy, g_0_z_xy_xxxyyyyz, g_0_z_xy_xxxyyyz, g_0_z_xy_xxxyyyzz, g_0_z_xy_xxxyyzz, g_0_z_xy_xxxyyzzz, g_0_z_xy_xxxyzzz, g_0_z_xy_xxxyzzzz, g_0_z_xy_xxxzzzz, g_0_z_xy_xxxzzzzz, g_0_z_xy_xxyyyyy, g_0_z_xy_xxyyyyyy, g_0_z_xy_xxyyyyyz, g_0_z_xy_xxyyyyz, g_0_z_xy_xxyyyyzz, g_0_z_xy_xxyyyzz, g_0_z_xy_xxyyyzzz, g_0_z_xy_xxyyzzz, g_0_z_xy_xxyyzzzz, g_0_z_xy_xxyzzzz, g_0_z_xy_xxyzzzzz, g_0_z_xy_xxzzzzz, g_0_z_xy_xxzzzzzz, g_0_z_xy_xyyyyyy, g_0_z_xy_xyyyyyyy, g_0_z_xy_xyyyyyyz, g_0_z_xy_xyyyyyz, g_0_z_xy_xyyyyyzz, g_0_z_xy_xyyyyzz, g_0_z_xy_xyyyyzzz, g_0_z_xy_xyyyzzz, g_0_z_xy_xyyyzzzz, g_0_z_xy_xyyzzzz, g_0_z_xy_xyyzzzzz, g_0_z_xy_xyzzzzz, g_0_z_xy_xyzzzzzz, g_0_z_xy_xzzzzzz, g_0_z_xy_xzzzzzzz, g_0_z_xy_yyyyyyy, g_0_z_xy_yyyyyyz, g_0_z_xy_yyyyyzz, g_0_z_xy_yyyyzzz, g_0_z_xy_yyyzzzz, g_0_z_xy_yyzzzzz, g_0_z_xy_yzzzzzz, g_0_z_xy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxy_xxxxxxx[k] = -g_0_z_xy_xxxxxxx[k] * ab_x + g_0_z_xy_xxxxxxxx[k];

                g_0_z_xxy_xxxxxxy[k] = -g_0_z_xy_xxxxxxy[k] * ab_x + g_0_z_xy_xxxxxxxy[k];

                g_0_z_xxy_xxxxxxz[k] = -g_0_z_xy_xxxxxxz[k] * ab_x + g_0_z_xy_xxxxxxxz[k];

                g_0_z_xxy_xxxxxyy[k] = -g_0_z_xy_xxxxxyy[k] * ab_x + g_0_z_xy_xxxxxxyy[k];

                g_0_z_xxy_xxxxxyz[k] = -g_0_z_xy_xxxxxyz[k] * ab_x + g_0_z_xy_xxxxxxyz[k];

                g_0_z_xxy_xxxxxzz[k] = -g_0_z_xy_xxxxxzz[k] * ab_x + g_0_z_xy_xxxxxxzz[k];

                g_0_z_xxy_xxxxyyy[k] = -g_0_z_xy_xxxxyyy[k] * ab_x + g_0_z_xy_xxxxxyyy[k];

                g_0_z_xxy_xxxxyyz[k] = -g_0_z_xy_xxxxyyz[k] * ab_x + g_0_z_xy_xxxxxyyz[k];

                g_0_z_xxy_xxxxyzz[k] = -g_0_z_xy_xxxxyzz[k] * ab_x + g_0_z_xy_xxxxxyzz[k];

                g_0_z_xxy_xxxxzzz[k] = -g_0_z_xy_xxxxzzz[k] * ab_x + g_0_z_xy_xxxxxzzz[k];

                g_0_z_xxy_xxxyyyy[k] = -g_0_z_xy_xxxyyyy[k] * ab_x + g_0_z_xy_xxxxyyyy[k];

                g_0_z_xxy_xxxyyyz[k] = -g_0_z_xy_xxxyyyz[k] * ab_x + g_0_z_xy_xxxxyyyz[k];

                g_0_z_xxy_xxxyyzz[k] = -g_0_z_xy_xxxyyzz[k] * ab_x + g_0_z_xy_xxxxyyzz[k];

                g_0_z_xxy_xxxyzzz[k] = -g_0_z_xy_xxxyzzz[k] * ab_x + g_0_z_xy_xxxxyzzz[k];

                g_0_z_xxy_xxxzzzz[k] = -g_0_z_xy_xxxzzzz[k] * ab_x + g_0_z_xy_xxxxzzzz[k];

                g_0_z_xxy_xxyyyyy[k] = -g_0_z_xy_xxyyyyy[k] * ab_x + g_0_z_xy_xxxyyyyy[k];

                g_0_z_xxy_xxyyyyz[k] = -g_0_z_xy_xxyyyyz[k] * ab_x + g_0_z_xy_xxxyyyyz[k];

                g_0_z_xxy_xxyyyzz[k] = -g_0_z_xy_xxyyyzz[k] * ab_x + g_0_z_xy_xxxyyyzz[k];

                g_0_z_xxy_xxyyzzz[k] = -g_0_z_xy_xxyyzzz[k] * ab_x + g_0_z_xy_xxxyyzzz[k];

                g_0_z_xxy_xxyzzzz[k] = -g_0_z_xy_xxyzzzz[k] * ab_x + g_0_z_xy_xxxyzzzz[k];

                g_0_z_xxy_xxzzzzz[k] = -g_0_z_xy_xxzzzzz[k] * ab_x + g_0_z_xy_xxxzzzzz[k];

                g_0_z_xxy_xyyyyyy[k] = -g_0_z_xy_xyyyyyy[k] * ab_x + g_0_z_xy_xxyyyyyy[k];

                g_0_z_xxy_xyyyyyz[k] = -g_0_z_xy_xyyyyyz[k] * ab_x + g_0_z_xy_xxyyyyyz[k];

                g_0_z_xxy_xyyyyzz[k] = -g_0_z_xy_xyyyyzz[k] * ab_x + g_0_z_xy_xxyyyyzz[k];

                g_0_z_xxy_xyyyzzz[k] = -g_0_z_xy_xyyyzzz[k] * ab_x + g_0_z_xy_xxyyyzzz[k];

                g_0_z_xxy_xyyzzzz[k] = -g_0_z_xy_xyyzzzz[k] * ab_x + g_0_z_xy_xxyyzzzz[k];

                g_0_z_xxy_xyzzzzz[k] = -g_0_z_xy_xyzzzzz[k] * ab_x + g_0_z_xy_xxyzzzzz[k];

                g_0_z_xxy_xzzzzzz[k] = -g_0_z_xy_xzzzzzz[k] * ab_x + g_0_z_xy_xxzzzzzz[k];

                g_0_z_xxy_yyyyyyy[k] = -g_0_z_xy_yyyyyyy[k] * ab_x + g_0_z_xy_xyyyyyyy[k];

                g_0_z_xxy_yyyyyyz[k] = -g_0_z_xy_yyyyyyz[k] * ab_x + g_0_z_xy_xyyyyyyz[k];

                g_0_z_xxy_yyyyyzz[k] = -g_0_z_xy_yyyyyzz[k] * ab_x + g_0_z_xy_xyyyyyzz[k];

                g_0_z_xxy_yyyyzzz[k] = -g_0_z_xy_yyyyzzz[k] * ab_x + g_0_z_xy_xyyyyzzz[k];

                g_0_z_xxy_yyyzzzz[k] = -g_0_z_xy_yyyzzzz[k] * ab_x + g_0_z_xy_xyyyzzzz[k];

                g_0_z_xxy_yyzzzzz[k] = -g_0_z_xy_yyzzzzz[k] * ab_x + g_0_z_xy_xyyzzzzz[k];

                g_0_z_xxy_yzzzzzz[k] = -g_0_z_xy_yzzzzzz[k] * ab_x + g_0_z_xy_xyzzzzzz[k];

                g_0_z_xxy_zzzzzzz[k] = -g_0_z_xy_zzzzzzz[k] * ab_x + g_0_z_xy_xzzzzzzz[k];
            }

            /// Set up 792-828 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxz_xxxxxxx, g_0_z_xxz_xxxxxxy, g_0_z_xxz_xxxxxxz, g_0_z_xxz_xxxxxyy, g_0_z_xxz_xxxxxyz, g_0_z_xxz_xxxxxzz, g_0_z_xxz_xxxxyyy, g_0_z_xxz_xxxxyyz, g_0_z_xxz_xxxxyzz, g_0_z_xxz_xxxxzzz, g_0_z_xxz_xxxyyyy, g_0_z_xxz_xxxyyyz, g_0_z_xxz_xxxyyzz, g_0_z_xxz_xxxyzzz, g_0_z_xxz_xxxzzzz, g_0_z_xxz_xxyyyyy, g_0_z_xxz_xxyyyyz, g_0_z_xxz_xxyyyzz, g_0_z_xxz_xxyyzzz, g_0_z_xxz_xxyzzzz, g_0_z_xxz_xxzzzzz, g_0_z_xxz_xyyyyyy, g_0_z_xxz_xyyyyyz, g_0_z_xxz_xyyyyzz, g_0_z_xxz_xyyyzzz, g_0_z_xxz_xyyzzzz, g_0_z_xxz_xyzzzzz, g_0_z_xxz_xzzzzzz, g_0_z_xxz_yyyyyyy, g_0_z_xxz_yyyyyyz, g_0_z_xxz_yyyyyzz, g_0_z_xxz_yyyyzzz, g_0_z_xxz_yyyzzzz, g_0_z_xxz_yyzzzzz, g_0_z_xxz_yzzzzzz, g_0_z_xxz_zzzzzzz, g_0_z_xz_xxxxxxx, g_0_z_xz_xxxxxxxx, g_0_z_xz_xxxxxxxy, g_0_z_xz_xxxxxxxz, g_0_z_xz_xxxxxxy, g_0_z_xz_xxxxxxyy, g_0_z_xz_xxxxxxyz, g_0_z_xz_xxxxxxz, g_0_z_xz_xxxxxxzz, g_0_z_xz_xxxxxyy, g_0_z_xz_xxxxxyyy, g_0_z_xz_xxxxxyyz, g_0_z_xz_xxxxxyz, g_0_z_xz_xxxxxyzz, g_0_z_xz_xxxxxzz, g_0_z_xz_xxxxxzzz, g_0_z_xz_xxxxyyy, g_0_z_xz_xxxxyyyy, g_0_z_xz_xxxxyyyz, g_0_z_xz_xxxxyyz, g_0_z_xz_xxxxyyzz, g_0_z_xz_xxxxyzz, g_0_z_xz_xxxxyzzz, g_0_z_xz_xxxxzzz, g_0_z_xz_xxxxzzzz, g_0_z_xz_xxxyyyy, g_0_z_xz_xxxyyyyy, g_0_z_xz_xxxyyyyz, g_0_z_xz_xxxyyyz, g_0_z_xz_xxxyyyzz, g_0_z_xz_xxxyyzz, g_0_z_xz_xxxyyzzz, g_0_z_xz_xxxyzzz, g_0_z_xz_xxxyzzzz, g_0_z_xz_xxxzzzz, g_0_z_xz_xxxzzzzz, g_0_z_xz_xxyyyyy, g_0_z_xz_xxyyyyyy, g_0_z_xz_xxyyyyyz, g_0_z_xz_xxyyyyz, g_0_z_xz_xxyyyyzz, g_0_z_xz_xxyyyzz, g_0_z_xz_xxyyyzzz, g_0_z_xz_xxyyzzz, g_0_z_xz_xxyyzzzz, g_0_z_xz_xxyzzzz, g_0_z_xz_xxyzzzzz, g_0_z_xz_xxzzzzz, g_0_z_xz_xxzzzzzz, g_0_z_xz_xyyyyyy, g_0_z_xz_xyyyyyyy, g_0_z_xz_xyyyyyyz, g_0_z_xz_xyyyyyz, g_0_z_xz_xyyyyyzz, g_0_z_xz_xyyyyzz, g_0_z_xz_xyyyyzzz, g_0_z_xz_xyyyzzz, g_0_z_xz_xyyyzzzz, g_0_z_xz_xyyzzzz, g_0_z_xz_xyyzzzzz, g_0_z_xz_xyzzzzz, g_0_z_xz_xyzzzzzz, g_0_z_xz_xzzzzzz, g_0_z_xz_xzzzzzzz, g_0_z_xz_yyyyyyy, g_0_z_xz_yyyyyyz, g_0_z_xz_yyyyyzz, g_0_z_xz_yyyyzzz, g_0_z_xz_yyyzzzz, g_0_z_xz_yyzzzzz, g_0_z_xz_yzzzzzz, g_0_z_xz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxz_xxxxxxx[k] = -g_0_z_xz_xxxxxxx[k] * ab_x + g_0_z_xz_xxxxxxxx[k];

                g_0_z_xxz_xxxxxxy[k] = -g_0_z_xz_xxxxxxy[k] * ab_x + g_0_z_xz_xxxxxxxy[k];

                g_0_z_xxz_xxxxxxz[k] = -g_0_z_xz_xxxxxxz[k] * ab_x + g_0_z_xz_xxxxxxxz[k];

                g_0_z_xxz_xxxxxyy[k] = -g_0_z_xz_xxxxxyy[k] * ab_x + g_0_z_xz_xxxxxxyy[k];

                g_0_z_xxz_xxxxxyz[k] = -g_0_z_xz_xxxxxyz[k] * ab_x + g_0_z_xz_xxxxxxyz[k];

                g_0_z_xxz_xxxxxzz[k] = -g_0_z_xz_xxxxxzz[k] * ab_x + g_0_z_xz_xxxxxxzz[k];

                g_0_z_xxz_xxxxyyy[k] = -g_0_z_xz_xxxxyyy[k] * ab_x + g_0_z_xz_xxxxxyyy[k];

                g_0_z_xxz_xxxxyyz[k] = -g_0_z_xz_xxxxyyz[k] * ab_x + g_0_z_xz_xxxxxyyz[k];

                g_0_z_xxz_xxxxyzz[k] = -g_0_z_xz_xxxxyzz[k] * ab_x + g_0_z_xz_xxxxxyzz[k];

                g_0_z_xxz_xxxxzzz[k] = -g_0_z_xz_xxxxzzz[k] * ab_x + g_0_z_xz_xxxxxzzz[k];

                g_0_z_xxz_xxxyyyy[k] = -g_0_z_xz_xxxyyyy[k] * ab_x + g_0_z_xz_xxxxyyyy[k];

                g_0_z_xxz_xxxyyyz[k] = -g_0_z_xz_xxxyyyz[k] * ab_x + g_0_z_xz_xxxxyyyz[k];

                g_0_z_xxz_xxxyyzz[k] = -g_0_z_xz_xxxyyzz[k] * ab_x + g_0_z_xz_xxxxyyzz[k];

                g_0_z_xxz_xxxyzzz[k] = -g_0_z_xz_xxxyzzz[k] * ab_x + g_0_z_xz_xxxxyzzz[k];

                g_0_z_xxz_xxxzzzz[k] = -g_0_z_xz_xxxzzzz[k] * ab_x + g_0_z_xz_xxxxzzzz[k];

                g_0_z_xxz_xxyyyyy[k] = -g_0_z_xz_xxyyyyy[k] * ab_x + g_0_z_xz_xxxyyyyy[k];

                g_0_z_xxz_xxyyyyz[k] = -g_0_z_xz_xxyyyyz[k] * ab_x + g_0_z_xz_xxxyyyyz[k];

                g_0_z_xxz_xxyyyzz[k] = -g_0_z_xz_xxyyyzz[k] * ab_x + g_0_z_xz_xxxyyyzz[k];

                g_0_z_xxz_xxyyzzz[k] = -g_0_z_xz_xxyyzzz[k] * ab_x + g_0_z_xz_xxxyyzzz[k];

                g_0_z_xxz_xxyzzzz[k] = -g_0_z_xz_xxyzzzz[k] * ab_x + g_0_z_xz_xxxyzzzz[k];

                g_0_z_xxz_xxzzzzz[k] = -g_0_z_xz_xxzzzzz[k] * ab_x + g_0_z_xz_xxxzzzzz[k];

                g_0_z_xxz_xyyyyyy[k] = -g_0_z_xz_xyyyyyy[k] * ab_x + g_0_z_xz_xxyyyyyy[k];

                g_0_z_xxz_xyyyyyz[k] = -g_0_z_xz_xyyyyyz[k] * ab_x + g_0_z_xz_xxyyyyyz[k];

                g_0_z_xxz_xyyyyzz[k] = -g_0_z_xz_xyyyyzz[k] * ab_x + g_0_z_xz_xxyyyyzz[k];

                g_0_z_xxz_xyyyzzz[k] = -g_0_z_xz_xyyyzzz[k] * ab_x + g_0_z_xz_xxyyyzzz[k];

                g_0_z_xxz_xyyzzzz[k] = -g_0_z_xz_xyyzzzz[k] * ab_x + g_0_z_xz_xxyyzzzz[k];

                g_0_z_xxz_xyzzzzz[k] = -g_0_z_xz_xyzzzzz[k] * ab_x + g_0_z_xz_xxyzzzzz[k];

                g_0_z_xxz_xzzzzzz[k] = -g_0_z_xz_xzzzzzz[k] * ab_x + g_0_z_xz_xxzzzzzz[k];

                g_0_z_xxz_yyyyyyy[k] = -g_0_z_xz_yyyyyyy[k] * ab_x + g_0_z_xz_xyyyyyyy[k];

                g_0_z_xxz_yyyyyyz[k] = -g_0_z_xz_yyyyyyz[k] * ab_x + g_0_z_xz_xyyyyyyz[k];

                g_0_z_xxz_yyyyyzz[k] = -g_0_z_xz_yyyyyzz[k] * ab_x + g_0_z_xz_xyyyyyzz[k];

                g_0_z_xxz_yyyyzzz[k] = -g_0_z_xz_yyyyzzz[k] * ab_x + g_0_z_xz_xyyyyzzz[k];

                g_0_z_xxz_yyyzzzz[k] = -g_0_z_xz_yyyzzzz[k] * ab_x + g_0_z_xz_xyyyzzzz[k];

                g_0_z_xxz_yyzzzzz[k] = -g_0_z_xz_yyzzzzz[k] * ab_x + g_0_z_xz_xyyzzzzz[k];

                g_0_z_xxz_yzzzzzz[k] = -g_0_z_xz_yzzzzzz[k] * ab_x + g_0_z_xz_xyzzzzzz[k];

                g_0_z_xxz_zzzzzzz[k] = -g_0_z_xz_zzzzzzz[k] * ab_x + g_0_z_xz_xzzzzzzz[k];
            }

            /// Set up 828-864 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyy_xxxxxxx, g_0_z_xyy_xxxxxxy, g_0_z_xyy_xxxxxxz, g_0_z_xyy_xxxxxyy, g_0_z_xyy_xxxxxyz, g_0_z_xyy_xxxxxzz, g_0_z_xyy_xxxxyyy, g_0_z_xyy_xxxxyyz, g_0_z_xyy_xxxxyzz, g_0_z_xyy_xxxxzzz, g_0_z_xyy_xxxyyyy, g_0_z_xyy_xxxyyyz, g_0_z_xyy_xxxyyzz, g_0_z_xyy_xxxyzzz, g_0_z_xyy_xxxzzzz, g_0_z_xyy_xxyyyyy, g_0_z_xyy_xxyyyyz, g_0_z_xyy_xxyyyzz, g_0_z_xyy_xxyyzzz, g_0_z_xyy_xxyzzzz, g_0_z_xyy_xxzzzzz, g_0_z_xyy_xyyyyyy, g_0_z_xyy_xyyyyyz, g_0_z_xyy_xyyyyzz, g_0_z_xyy_xyyyzzz, g_0_z_xyy_xyyzzzz, g_0_z_xyy_xyzzzzz, g_0_z_xyy_xzzzzzz, g_0_z_xyy_yyyyyyy, g_0_z_xyy_yyyyyyz, g_0_z_xyy_yyyyyzz, g_0_z_xyy_yyyyzzz, g_0_z_xyy_yyyzzzz, g_0_z_xyy_yyzzzzz, g_0_z_xyy_yzzzzzz, g_0_z_xyy_zzzzzzz, g_0_z_yy_xxxxxxx, g_0_z_yy_xxxxxxxx, g_0_z_yy_xxxxxxxy, g_0_z_yy_xxxxxxxz, g_0_z_yy_xxxxxxy, g_0_z_yy_xxxxxxyy, g_0_z_yy_xxxxxxyz, g_0_z_yy_xxxxxxz, g_0_z_yy_xxxxxxzz, g_0_z_yy_xxxxxyy, g_0_z_yy_xxxxxyyy, g_0_z_yy_xxxxxyyz, g_0_z_yy_xxxxxyz, g_0_z_yy_xxxxxyzz, g_0_z_yy_xxxxxzz, g_0_z_yy_xxxxxzzz, g_0_z_yy_xxxxyyy, g_0_z_yy_xxxxyyyy, g_0_z_yy_xxxxyyyz, g_0_z_yy_xxxxyyz, g_0_z_yy_xxxxyyzz, g_0_z_yy_xxxxyzz, g_0_z_yy_xxxxyzzz, g_0_z_yy_xxxxzzz, g_0_z_yy_xxxxzzzz, g_0_z_yy_xxxyyyy, g_0_z_yy_xxxyyyyy, g_0_z_yy_xxxyyyyz, g_0_z_yy_xxxyyyz, g_0_z_yy_xxxyyyzz, g_0_z_yy_xxxyyzz, g_0_z_yy_xxxyyzzz, g_0_z_yy_xxxyzzz, g_0_z_yy_xxxyzzzz, g_0_z_yy_xxxzzzz, g_0_z_yy_xxxzzzzz, g_0_z_yy_xxyyyyy, g_0_z_yy_xxyyyyyy, g_0_z_yy_xxyyyyyz, g_0_z_yy_xxyyyyz, g_0_z_yy_xxyyyyzz, g_0_z_yy_xxyyyzz, g_0_z_yy_xxyyyzzz, g_0_z_yy_xxyyzzz, g_0_z_yy_xxyyzzzz, g_0_z_yy_xxyzzzz, g_0_z_yy_xxyzzzzz, g_0_z_yy_xxzzzzz, g_0_z_yy_xxzzzzzz, g_0_z_yy_xyyyyyy, g_0_z_yy_xyyyyyyy, g_0_z_yy_xyyyyyyz, g_0_z_yy_xyyyyyz, g_0_z_yy_xyyyyyzz, g_0_z_yy_xyyyyzz, g_0_z_yy_xyyyyzzz, g_0_z_yy_xyyyzzz, g_0_z_yy_xyyyzzzz, g_0_z_yy_xyyzzzz, g_0_z_yy_xyyzzzzz, g_0_z_yy_xyzzzzz, g_0_z_yy_xyzzzzzz, g_0_z_yy_xzzzzzz, g_0_z_yy_xzzzzzzz, g_0_z_yy_yyyyyyy, g_0_z_yy_yyyyyyz, g_0_z_yy_yyyyyzz, g_0_z_yy_yyyyzzz, g_0_z_yy_yyyzzzz, g_0_z_yy_yyzzzzz, g_0_z_yy_yzzzzzz, g_0_z_yy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyy_xxxxxxx[k] = -g_0_z_yy_xxxxxxx[k] * ab_x + g_0_z_yy_xxxxxxxx[k];

                g_0_z_xyy_xxxxxxy[k] = -g_0_z_yy_xxxxxxy[k] * ab_x + g_0_z_yy_xxxxxxxy[k];

                g_0_z_xyy_xxxxxxz[k] = -g_0_z_yy_xxxxxxz[k] * ab_x + g_0_z_yy_xxxxxxxz[k];

                g_0_z_xyy_xxxxxyy[k] = -g_0_z_yy_xxxxxyy[k] * ab_x + g_0_z_yy_xxxxxxyy[k];

                g_0_z_xyy_xxxxxyz[k] = -g_0_z_yy_xxxxxyz[k] * ab_x + g_0_z_yy_xxxxxxyz[k];

                g_0_z_xyy_xxxxxzz[k] = -g_0_z_yy_xxxxxzz[k] * ab_x + g_0_z_yy_xxxxxxzz[k];

                g_0_z_xyy_xxxxyyy[k] = -g_0_z_yy_xxxxyyy[k] * ab_x + g_0_z_yy_xxxxxyyy[k];

                g_0_z_xyy_xxxxyyz[k] = -g_0_z_yy_xxxxyyz[k] * ab_x + g_0_z_yy_xxxxxyyz[k];

                g_0_z_xyy_xxxxyzz[k] = -g_0_z_yy_xxxxyzz[k] * ab_x + g_0_z_yy_xxxxxyzz[k];

                g_0_z_xyy_xxxxzzz[k] = -g_0_z_yy_xxxxzzz[k] * ab_x + g_0_z_yy_xxxxxzzz[k];

                g_0_z_xyy_xxxyyyy[k] = -g_0_z_yy_xxxyyyy[k] * ab_x + g_0_z_yy_xxxxyyyy[k];

                g_0_z_xyy_xxxyyyz[k] = -g_0_z_yy_xxxyyyz[k] * ab_x + g_0_z_yy_xxxxyyyz[k];

                g_0_z_xyy_xxxyyzz[k] = -g_0_z_yy_xxxyyzz[k] * ab_x + g_0_z_yy_xxxxyyzz[k];

                g_0_z_xyy_xxxyzzz[k] = -g_0_z_yy_xxxyzzz[k] * ab_x + g_0_z_yy_xxxxyzzz[k];

                g_0_z_xyy_xxxzzzz[k] = -g_0_z_yy_xxxzzzz[k] * ab_x + g_0_z_yy_xxxxzzzz[k];

                g_0_z_xyy_xxyyyyy[k] = -g_0_z_yy_xxyyyyy[k] * ab_x + g_0_z_yy_xxxyyyyy[k];

                g_0_z_xyy_xxyyyyz[k] = -g_0_z_yy_xxyyyyz[k] * ab_x + g_0_z_yy_xxxyyyyz[k];

                g_0_z_xyy_xxyyyzz[k] = -g_0_z_yy_xxyyyzz[k] * ab_x + g_0_z_yy_xxxyyyzz[k];

                g_0_z_xyy_xxyyzzz[k] = -g_0_z_yy_xxyyzzz[k] * ab_x + g_0_z_yy_xxxyyzzz[k];

                g_0_z_xyy_xxyzzzz[k] = -g_0_z_yy_xxyzzzz[k] * ab_x + g_0_z_yy_xxxyzzzz[k];

                g_0_z_xyy_xxzzzzz[k] = -g_0_z_yy_xxzzzzz[k] * ab_x + g_0_z_yy_xxxzzzzz[k];

                g_0_z_xyy_xyyyyyy[k] = -g_0_z_yy_xyyyyyy[k] * ab_x + g_0_z_yy_xxyyyyyy[k];

                g_0_z_xyy_xyyyyyz[k] = -g_0_z_yy_xyyyyyz[k] * ab_x + g_0_z_yy_xxyyyyyz[k];

                g_0_z_xyy_xyyyyzz[k] = -g_0_z_yy_xyyyyzz[k] * ab_x + g_0_z_yy_xxyyyyzz[k];

                g_0_z_xyy_xyyyzzz[k] = -g_0_z_yy_xyyyzzz[k] * ab_x + g_0_z_yy_xxyyyzzz[k];

                g_0_z_xyy_xyyzzzz[k] = -g_0_z_yy_xyyzzzz[k] * ab_x + g_0_z_yy_xxyyzzzz[k];

                g_0_z_xyy_xyzzzzz[k] = -g_0_z_yy_xyzzzzz[k] * ab_x + g_0_z_yy_xxyzzzzz[k];

                g_0_z_xyy_xzzzzzz[k] = -g_0_z_yy_xzzzzzz[k] * ab_x + g_0_z_yy_xxzzzzzz[k];

                g_0_z_xyy_yyyyyyy[k] = -g_0_z_yy_yyyyyyy[k] * ab_x + g_0_z_yy_xyyyyyyy[k];

                g_0_z_xyy_yyyyyyz[k] = -g_0_z_yy_yyyyyyz[k] * ab_x + g_0_z_yy_xyyyyyyz[k];

                g_0_z_xyy_yyyyyzz[k] = -g_0_z_yy_yyyyyzz[k] * ab_x + g_0_z_yy_xyyyyyzz[k];

                g_0_z_xyy_yyyyzzz[k] = -g_0_z_yy_yyyyzzz[k] * ab_x + g_0_z_yy_xyyyyzzz[k];

                g_0_z_xyy_yyyzzzz[k] = -g_0_z_yy_yyyzzzz[k] * ab_x + g_0_z_yy_xyyyzzzz[k];

                g_0_z_xyy_yyzzzzz[k] = -g_0_z_yy_yyzzzzz[k] * ab_x + g_0_z_yy_xyyzzzzz[k];

                g_0_z_xyy_yzzzzzz[k] = -g_0_z_yy_yzzzzzz[k] * ab_x + g_0_z_yy_xyzzzzzz[k];

                g_0_z_xyy_zzzzzzz[k] = -g_0_z_yy_zzzzzzz[k] * ab_x + g_0_z_yy_xzzzzzzz[k];
            }

            /// Set up 864-900 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyz_xxxxxxx, g_0_z_xyz_xxxxxxy, g_0_z_xyz_xxxxxxz, g_0_z_xyz_xxxxxyy, g_0_z_xyz_xxxxxyz, g_0_z_xyz_xxxxxzz, g_0_z_xyz_xxxxyyy, g_0_z_xyz_xxxxyyz, g_0_z_xyz_xxxxyzz, g_0_z_xyz_xxxxzzz, g_0_z_xyz_xxxyyyy, g_0_z_xyz_xxxyyyz, g_0_z_xyz_xxxyyzz, g_0_z_xyz_xxxyzzz, g_0_z_xyz_xxxzzzz, g_0_z_xyz_xxyyyyy, g_0_z_xyz_xxyyyyz, g_0_z_xyz_xxyyyzz, g_0_z_xyz_xxyyzzz, g_0_z_xyz_xxyzzzz, g_0_z_xyz_xxzzzzz, g_0_z_xyz_xyyyyyy, g_0_z_xyz_xyyyyyz, g_0_z_xyz_xyyyyzz, g_0_z_xyz_xyyyzzz, g_0_z_xyz_xyyzzzz, g_0_z_xyz_xyzzzzz, g_0_z_xyz_xzzzzzz, g_0_z_xyz_yyyyyyy, g_0_z_xyz_yyyyyyz, g_0_z_xyz_yyyyyzz, g_0_z_xyz_yyyyzzz, g_0_z_xyz_yyyzzzz, g_0_z_xyz_yyzzzzz, g_0_z_xyz_yzzzzzz, g_0_z_xyz_zzzzzzz, g_0_z_yz_xxxxxxx, g_0_z_yz_xxxxxxxx, g_0_z_yz_xxxxxxxy, g_0_z_yz_xxxxxxxz, g_0_z_yz_xxxxxxy, g_0_z_yz_xxxxxxyy, g_0_z_yz_xxxxxxyz, g_0_z_yz_xxxxxxz, g_0_z_yz_xxxxxxzz, g_0_z_yz_xxxxxyy, g_0_z_yz_xxxxxyyy, g_0_z_yz_xxxxxyyz, g_0_z_yz_xxxxxyz, g_0_z_yz_xxxxxyzz, g_0_z_yz_xxxxxzz, g_0_z_yz_xxxxxzzz, g_0_z_yz_xxxxyyy, g_0_z_yz_xxxxyyyy, g_0_z_yz_xxxxyyyz, g_0_z_yz_xxxxyyz, g_0_z_yz_xxxxyyzz, g_0_z_yz_xxxxyzz, g_0_z_yz_xxxxyzzz, g_0_z_yz_xxxxzzz, g_0_z_yz_xxxxzzzz, g_0_z_yz_xxxyyyy, g_0_z_yz_xxxyyyyy, g_0_z_yz_xxxyyyyz, g_0_z_yz_xxxyyyz, g_0_z_yz_xxxyyyzz, g_0_z_yz_xxxyyzz, g_0_z_yz_xxxyyzzz, g_0_z_yz_xxxyzzz, g_0_z_yz_xxxyzzzz, g_0_z_yz_xxxzzzz, g_0_z_yz_xxxzzzzz, g_0_z_yz_xxyyyyy, g_0_z_yz_xxyyyyyy, g_0_z_yz_xxyyyyyz, g_0_z_yz_xxyyyyz, g_0_z_yz_xxyyyyzz, g_0_z_yz_xxyyyzz, g_0_z_yz_xxyyyzzz, g_0_z_yz_xxyyzzz, g_0_z_yz_xxyyzzzz, g_0_z_yz_xxyzzzz, g_0_z_yz_xxyzzzzz, g_0_z_yz_xxzzzzz, g_0_z_yz_xxzzzzzz, g_0_z_yz_xyyyyyy, g_0_z_yz_xyyyyyyy, g_0_z_yz_xyyyyyyz, g_0_z_yz_xyyyyyz, g_0_z_yz_xyyyyyzz, g_0_z_yz_xyyyyzz, g_0_z_yz_xyyyyzzz, g_0_z_yz_xyyyzzz, g_0_z_yz_xyyyzzzz, g_0_z_yz_xyyzzzz, g_0_z_yz_xyyzzzzz, g_0_z_yz_xyzzzzz, g_0_z_yz_xyzzzzzz, g_0_z_yz_xzzzzzz, g_0_z_yz_xzzzzzzz, g_0_z_yz_yyyyyyy, g_0_z_yz_yyyyyyz, g_0_z_yz_yyyyyzz, g_0_z_yz_yyyyzzz, g_0_z_yz_yyyzzzz, g_0_z_yz_yyzzzzz, g_0_z_yz_yzzzzzz, g_0_z_yz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyz_xxxxxxx[k] = -g_0_z_yz_xxxxxxx[k] * ab_x + g_0_z_yz_xxxxxxxx[k];

                g_0_z_xyz_xxxxxxy[k] = -g_0_z_yz_xxxxxxy[k] * ab_x + g_0_z_yz_xxxxxxxy[k];

                g_0_z_xyz_xxxxxxz[k] = -g_0_z_yz_xxxxxxz[k] * ab_x + g_0_z_yz_xxxxxxxz[k];

                g_0_z_xyz_xxxxxyy[k] = -g_0_z_yz_xxxxxyy[k] * ab_x + g_0_z_yz_xxxxxxyy[k];

                g_0_z_xyz_xxxxxyz[k] = -g_0_z_yz_xxxxxyz[k] * ab_x + g_0_z_yz_xxxxxxyz[k];

                g_0_z_xyz_xxxxxzz[k] = -g_0_z_yz_xxxxxzz[k] * ab_x + g_0_z_yz_xxxxxxzz[k];

                g_0_z_xyz_xxxxyyy[k] = -g_0_z_yz_xxxxyyy[k] * ab_x + g_0_z_yz_xxxxxyyy[k];

                g_0_z_xyz_xxxxyyz[k] = -g_0_z_yz_xxxxyyz[k] * ab_x + g_0_z_yz_xxxxxyyz[k];

                g_0_z_xyz_xxxxyzz[k] = -g_0_z_yz_xxxxyzz[k] * ab_x + g_0_z_yz_xxxxxyzz[k];

                g_0_z_xyz_xxxxzzz[k] = -g_0_z_yz_xxxxzzz[k] * ab_x + g_0_z_yz_xxxxxzzz[k];

                g_0_z_xyz_xxxyyyy[k] = -g_0_z_yz_xxxyyyy[k] * ab_x + g_0_z_yz_xxxxyyyy[k];

                g_0_z_xyz_xxxyyyz[k] = -g_0_z_yz_xxxyyyz[k] * ab_x + g_0_z_yz_xxxxyyyz[k];

                g_0_z_xyz_xxxyyzz[k] = -g_0_z_yz_xxxyyzz[k] * ab_x + g_0_z_yz_xxxxyyzz[k];

                g_0_z_xyz_xxxyzzz[k] = -g_0_z_yz_xxxyzzz[k] * ab_x + g_0_z_yz_xxxxyzzz[k];

                g_0_z_xyz_xxxzzzz[k] = -g_0_z_yz_xxxzzzz[k] * ab_x + g_0_z_yz_xxxxzzzz[k];

                g_0_z_xyz_xxyyyyy[k] = -g_0_z_yz_xxyyyyy[k] * ab_x + g_0_z_yz_xxxyyyyy[k];

                g_0_z_xyz_xxyyyyz[k] = -g_0_z_yz_xxyyyyz[k] * ab_x + g_0_z_yz_xxxyyyyz[k];

                g_0_z_xyz_xxyyyzz[k] = -g_0_z_yz_xxyyyzz[k] * ab_x + g_0_z_yz_xxxyyyzz[k];

                g_0_z_xyz_xxyyzzz[k] = -g_0_z_yz_xxyyzzz[k] * ab_x + g_0_z_yz_xxxyyzzz[k];

                g_0_z_xyz_xxyzzzz[k] = -g_0_z_yz_xxyzzzz[k] * ab_x + g_0_z_yz_xxxyzzzz[k];

                g_0_z_xyz_xxzzzzz[k] = -g_0_z_yz_xxzzzzz[k] * ab_x + g_0_z_yz_xxxzzzzz[k];

                g_0_z_xyz_xyyyyyy[k] = -g_0_z_yz_xyyyyyy[k] * ab_x + g_0_z_yz_xxyyyyyy[k];

                g_0_z_xyz_xyyyyyz[k] = -g_0_z_yz_xyyyyyz[k] * ab_x + g_0_z_yz_xxyyyyyz[k];

                g_0_z_xyz_xyyyyzz[k] = -g_0_z_yz_xyyyyzz[k] * ab_x + g_0_z_yz_xxyyyyzz[k];

                g_0_z_xyz_xyyyzzz[k] = -g_0_z_yz_xyyyzzz[k] * ab_x + g_0_z_yz_xxyyyzzz[k];

                g_0_z_xyz_xyyzzzz[k] = -g_0_z_yz_xyyzzzz[k] * ab_x + g_0_z_yz_xxyyzzzz[k];

                g_0_z_xyz_xyzzzzz[k] = -g_0_z_yz_xyzzzzz[k] * ab_x + g_0_z_yz_xxyzzzzz[k];

                g_0_z_xyz_xzzzzzz[k] = -g_0_z_yz_xzzzzzz[k] * ab_x + g_0_z_yz_xxzzzzzz[k];

                g_0_z_xyz_yyyyyyy[k] = -g_0_z_yz_yyyyyyy[k] * ab_x + g_0_z_yz_xyyyyyyy[k];

                g_0_z_xyz_yyyyyyz[k] = -g_0_z_yz_yyyyyyz[k] * ab_x + g_0_z_yz_xyyyyyyz[k];

                g_0_z_xyz_yyyyyzz[k] = -g_0_z_yz_yyyyyzz[k] * ab_x + g_0_z_yz_xyyyyyzz[k];

                g_0_z_xyz_yyyyzzz[k] = -g_0_z_yz_yyyyzzz[k] * ab_x + g_0_z_yz_xyyyyzzz[k];

                g_0_z_xyz_yyyzzzz[k] = -g_0_z_yz_yyyzzzz[k] * ab_x + g_0_z_yz_xyyyzzzz[k];

                g_0_z_xyz_yyzzzzz[k] = -g_0_z_yz_yyzzzzz[k] * ab_x + g_0_z_yz_xyyzzzzz[k];

                g_0_z_xyz_yzzzzzz[k] = -g_0_z_yz_yzzzzzz[k] * ab_x + g_0_z_yz_xyzzzzzz[k];

                g_0_z_xyz_zzzzzzz[k] = -g_0_z_yz_zzzzzzz[k] * ab_x + g_0_z_yz_xzzzzzzz[k];
            }

            /// Set up 900-936 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xzz_xxxxxxx, g_0_z_xzz_xxxxxxy, g_0_z_xzz_xxxxxxz, g_0_z_xzz_xxxxxyy, g_0_z_xzz_xxxxxyz, g_0_z_xzz_xxxxxzz, g_0_z_xzz_xxxxyyy, g_0_z_xzz_xxxxyyz, g_0_z_xzz_xxxxyzz, g_0_z_xzz_xxxxzzz, g_0_z_xzz_xxxyyyy, g_0_z_xzz_xxxyyyz, g_0_z_xzz_xxxyyzz, g_0_z_xzz_xxxyzzz, g_0_z_xzz_xxxzzzz, g_0_z_xzz_xxyyyyy, g_0_z_xzz_xxyyyyz, g_0_z_xzz_xxyyyzz, g_0_z_xzz_xxyyzzz, g_0_z_xzz_xxyzzzz, g_0_z_xzz_xxzzzzz, g_0_z_xzz_xyyyyyy, g_0_z_xzz_xyyyyyz, g_0_z_xzz_xyyyyzz, g_0_z_xzz_xyyyzzz, g_0_z_xzz_xyyzzzz, g_0_z_xzz_xyzzzzz, g_0_z_xzz_xzzzzzz, g_0_z_xzz_yyyyyyy, g_0_z_xzz_yyyyyyz, g_0_z_xzz_yyyyyzz, g_0_z_xzz_yyyyzzz, g_0_z_xzz_yyyzzzz, g_0_z_xzz_yyzzzzz, g_0_z_xzz_yzzzzzz, g_0_z_xzz_zzzzzzz, g_0_z_zz_xxxxxxx, g_0_z_zz_xxxxxxxx, g_0_z_zz_xxxxxxxy, g_0_z_zz_xxxxxxxz, g_0_z_zz_xxxxxxy, g_0_z_zz_xxxxxxyy, g_0_z_zz_xxxxxxyz, g_0_z_zz_xxxxxxz, g_0_z_zz_xxxxxxzz, g_0_z_zz_xxxxxyy, g_0_z_zz_xxxxxyyy, g_0_z_zz_xxxxxyyz, g_0_z_zz_xxxxxyz, g_0_z_zz_xxxxxyzz, g_0_z_zz_xxxxxzz, g_0_z_zz_xxxxxzzz, g_0_z_zz_xxxxyyy, g_0_z_zz_xxxxyyyy, g_0_z_zz_xxxxyyyz, g_0_z_zz_xxxxyyz, g_0_z_zz_xxxxyyzz, g_0_z_zz_xxxxyzz, g_0_z_zz_xxxxyzzz, g_0_z_zz_xxxxzzz, g_0_z_zz_xxxxzzzz, g_0_z_zz_xxxyyyy, g_0_z_zz_xxxyyyyy, g_0_z_zz_xxxyyyyz, g_0_z_zz_xxxyyyz, g_0_z_zz_xxxyyyzz, g_0_z_zz_xxxyyzz, g_0_z_zz_xxxyyzzz, g_0_z_zz_xxxyzzz, g_0_z_zz_xxxyzzzz, g_0_z_zz_xxxzzzz, g_0_z_zz_xxxzzzzz, g_0_z_zz_xxyyyyy, g_0_z_zz_xxyyyyyy, g_0_z_zz_xxyyyyyz, g_0_z_zz_xxyyyyz, g_0_z_zz_xxyyyyzz, g_0_z_zz_xxyyyzz, g_0_z_zz_xxyyyzzz, g_0_z_zz_xxyyzzz, g_0_z_zz_xxyyzzzz, g_0_z_zz_xxyzzzz, g_0_z_zz_xxyzzzzz, g_0_z_zz_xxzzzzz, g_0_z_zz_xxzzzzzz, g_0_z_zz_xyyyyyy, g_0_z_zz_xyyyyyyy, g_0_z_zz_xyyyyyyz, g_0_z_zz_xyyyyyz, g_0_z_zz_xyyyyyzz, g_0_z_zz_xyyyyzz, g_0_z_zz_xyyyyzzz, g_0_z_zz_xyyyzzz, g_0_z_zz_xyyyzzzz, g_0_z_zz_xyyzzzz, g_0_z_zz_xyyzzzzz, g_0_z_zz_xyzzzzz, g_0_z_zz_xyzzzzzz, g_0_z_zz_xzzzzzz, g_0_z_zz_xzzzzzzz, g_0_z_zz_yyyyyyy, g_0_z_zz_yyyyyyz, g_0_z_zz_yyyyyzz, g_0_z_zz_yyyyzzz, g_0_z_zz_yyyzzzz, g_0_z_zz_yyzzzzz, g_0_z_zz_yzzzzzz, g_0_z_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzz_xxxxxxx[k] = -g_0_z_zz_xxxxxxx[k] * ab_x + g_0_z_zz_xxxxxxxx[k];

                g_0_z_xzz_xxxxxxy[k] = -g_0_z_zz_xxxxxxy[k] * ab_x + g_0_z_zz_xxxxxxxy[k];

                g_0_z_xzz_xxxxxxz[k] = -g_0_z_zz_xxxxxxz[k] * ab_x + g_0_z_zz_xxxxxxxz[k];

                g_0_z_xzz_xxxxxyy[k] = -g_0_z_zz_xxxxxyy[k] * ab_x + g_0_z_zz_xxxxxxyy[k];

                g_0_z_xzz_xxxxxyz[k] = -g_0_z_zz_xxxxxyz[k] * ab_x + g_0_z_zz_xxxxxxyz[k];

                g_0_z_xzz_xxxxxzz[k] = -g_0_z_zz_xxxxxzz[k] * ab_x + g_0_z_zz_xxxxxxzz[k];

                g_0_z_xzz_xxxxyyy[k] = -g_0_z_zz_xxxxyyy[k] * ab_x + g_0_z_zz_xxxxxyyy[k];

                g_0_z_xzz_xxxxyyz[k] = -g_0_z_zz_xxxxyyz[k] * ab_x + g_0_z_zz_xxxxxyyz[k];

                g_0_z_xzz_xxxxyzz[k] = -g_0_z_zz_xxxxyzz[k] * ab_x + g_0_z_zz_xxxxxyzz[k];

                g_0_z_xzz_xxxxzzz[k] = -g_0_z_zz_xxxxzzz[k] * ab_x + g_0_z_zz_xxxxxzzz[k];

                g_0_z_xzz_xxxyyyy[k] = -g_0_z_zz_xxxyyyy[k] * ab_x + g_0_z_zz_xxxxyyyy[k];

                g_0_z_xzz_xxxyyyz[k] = -g_0_z_zz_xxxyyyz[k] * ab_x + g_0_z_zz_xxxxyyyz[k];

                g_0_z_xzz_xxxyyzz[k] = -g_0_z_zz_xxxyyzz[k] * ab_x + g_0_z_zz_xxxxyyzz[k];

                g_0_z_xzz_xxxyzzz[k] = -g_0_z_zz_xxxyzzz[k] * ab_x + g_0_z_zz_xxxxyzzz[k];

                g_0_z_xzz_xxxzzzz[k] = -g_0_z_zz_xxxzzzz[k] * ab_x + g_0_z_zz_xxxxzzzz[k];

                g_0_z_xzz_xxyyyyy[k] = -g_0_z_zz_xxyyyyy[k] * ab_x + g_0_z_zz_xxxyyyyy[k];

                g_0_z_xzz_xxyyyyz[k] = -g_0_z_zz_xxyyyyz[k] * ab_x + g_0_z_zz_xxxyyyyz[k];

                g_0_z_xzz_xxyyyzz[k] = -g_0_z_zz_xxyyyzz[k] * ab_x + g_0_z_zz_xxxyyyzz[k];

                g_0_z_xzz_xxyyzzz[k] = -g_0_z_zz_xxyyzzz[k] * ab_x + g_0_z_zz_xxxyyzzz[k];

                g_0_z_xzz_xxyzzzz[k] = -g_0_z_zz_xxyzzzz[k] * ab_x + g_0_z_zz_xxxyzzzz[k];

                g_0_z_xzz_xxzzzzz[k] = -g_0_z_zz_xxzzzzz[k] * ab_x + g_0_z_zz_xxxzzzzz[k];

                g_0_z_xzz_xyyyyyy[k] = -g_0_z_zz_xyyyyyy[k] * ab_x + g_0_z_zz_xxyyyyyy[k];

                g_0_z_xzz_xyyyyyz[k] = -g_0_z_zz_xyyyyyz[k] * ab_x + g_0_z_zz_xxyyyyyz[k];

                g_0_z_xzz_xyyyyzz[k] = -g_0_z_zz_xyyyyzz[k] * ab_x + g_0_z_zz_xxyyyyzz[k];

                g_0_z_xzz_xyyyzzz[k] = -g_0_z_zz_xyyyzzz[k] * ab_x + g_0_z_zz_xxyyyzzz[k];

                g_0_z_xzz_xyyzzzz[k] = -g_0_z_zz_xyyzzzz[k] * ab_x + g_0_z_zz_xxyyzzzz[k];

                g_0_z_xzz_xyzzzzz[k] = -g_0_z_zz_xyzzzzz[k] * ab_x + g_0_z_zz_xxyzzzzz[k];

                g_0_z_xzz_xzzzzzz[k] = -g_0_z_zz_xzzzzzz[k] * ab_x + g_0_z_zz_xxzzzzzz[k];

                g_0_z_xzz_yyyyyyy[k] = -g_0_z_zz_yyyyyyy[k] * ab_x + g_0_z_zz_xyyyyyyy[k];

                g_0_z_xzz_yyyyyyz[k] = -g_0_z_zz_yyyyyyz[k] * ab_x + g_0_z_zz_xyyyyyyz[k];

                g_0_z_xzz_yyyyyzz[k] = -g_0_z_zz_yyyyyzz[k] * ab_x + g_0_z_zz_xyyyyyzz[k];

                g_0_z_xzz_yyyyzzz[k] = -g_0_z_zz_yyyyzzz[k] * ab_x + g_0_z_zz_xyyyyzzz[k];

                g_0_z_xzz_yyyzzzz[k] = -g_0_z_zz_yyyzzzz[k] * ab_x + g_0_z_zz_xyyyzzzz[k];

                g_0_z_xzz_yyzzzzz[k] = -g_0_z_zz_yyzzzzz[k] * ab_x + g_0_z_zz_xyyzzzzz[k];

                g_0_z_xzz_yzzzzzz[k] = -g_0_z_zz_yzzzzzz[k] * ab_x + g_0_z_zz_xyzzzzzz[k];

                g_0_z_xzz_zzzzzzz[k] = -g_0_z_zz_zzzzzzz[k] * ab_x + g_0_z_zz_xzzzzzzz[k];
            }

            /// Set up 936-972 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yy_xxxxxxx, g_0_z_yy_xxxxxxxy, g_0_z_yy_xxxxxxy, g_0_z_yy_xxxxxxyy, g_0_z_yy_xxxxxxyz, g_0_z_yy_xxxxxxz, g_0_z_yy_xxxxxyy, g_0_z_yy_xxxxxyyy, g_0_z_yy_xxxxxyyz, g_0_z_yy_xxxxxyz, g_0_z_yy_xxxxxyzz, g_0_z_yy_xxxxxzz, g_0_z_yy_xxxxyyy, g_0_z_yy_xxxxyyyy, g_0_z_yy_xxxxyyyz, g_0_z_yy_xxxxyyz, g_0_z_yy_xxxxyyzz, g_0_z_yy_xxxxyzz, g_0_z_yy_xxxxyzzz, g_0_z_yy_xxxxzzz, g_0_z_yy_xxxyyyy, g_0_z_yy_xxxyyyyy, g_0_z_yy_xxxyyyyz, g_0_z_yy_xxxyyyz, g_0_z_yy_xxxyyyzz, g_0_z_yy_xxxyyzz, g_0_z_yy_xxxyyzzz, g_0_z_yy_xxxyzzz, g_0_z_yy_xxxyzzzz, g_0_z_yy_xxxzzzz, g_0_z_yy_xxyyyyy, g_0_z_yy_xxyyyyyy, g_0_z_yy_xxyyyyyz, g_0_z_yy_xxyyyyz, g_0_z_yy_xxyyyyzz, g_0_z_yy_xxyyyzz, g_0_z_yy_xxyyyzzz, g_0_z_yy_xxyyzzz, g_0_z_yy_xxyyzzzz, g_0_z_yy_xxyzzzz, g_0_z_yy_xxyzzzzz, g_0_z_yy_xxzzzzz, g_0_z_yy_xyyyyyy, g_0_z_yy_xyyyyyyy, g_0_z_yy_xyyyyyyz, g_0_z_yy_xyyyyyz, g_0_z_yy_xyyyyyzz, g_0_z_yy_xyyyyzz, g_0_z_yy_xyyyyzzz, g_0_z_yy_xyyyzzz, g_0_z_yy_xyyyzzzz, g_0_z_yy_xyyzzzz, g_0_z_yy_xyyzzzzz, g_0_z_yy_xyzzzzz, g_0_z_yy_xyzzzzzz, g_0_z_yy_xzzzzzz, g_0_z_yy_yyyyyyy, g_0_z_yy_yyyyyyyy, g_0_z_yy_yyyyyyyz, g_0_z_yy_yyyyyyz, g_0_z_yy_yyyyyyzz, g_0_z_yy_yyyyyzz, g_0_z_yy_yyyyyzzz, g_0_z_yy_yyyyzzz, g_0_z_yy_yyyyzzzz, g_0_z_yy_yyyzzzz, g_0_z_yy_yyyzzzzz, g_0_z_yy_yyzzzzz, g_0_z_yy_yyzzzzzz, g_0_z_yy_yzzzzzz, g_0_z_yy_yzzzzzzz, g_0_z_yy_zzzzzzz, g_0_z_yyy_xxxxxxx, g_0_z_yyy_xxxxxxy, g_0_z_yyy_xxxxxxz, g_0_z_yyy_xxxxxyy, g_0_z_yyy_xxxxxyz, g_0_z_yyy_xxxxxzz, g_0_z_yyy_xxxxyyy, g_0_z_yyy_xxxxyyz, g_0_z_yyy_xxxxyzz, g_0_z_yyy_xxxxzzz, g_0_z_yyy_xxxyyyy, g_0_z_yyy_xxxyyyz, g_0_z_yyy_xxxyyzz, g_0_z_yyy_xxxyzzz, g_0_z_yyy_xxxzzzz, g_0_z_yyy_xxyyyyy, g_0_z_yyy_xxyyyyz, g_0_z_yyy_xxyyyzz, g_0_z_yyy_xxyyzzz, g_0_z_yyy_xxyzzzz, g_0_z_yyy_xxzzzzz, g_0_z_yyy_xyyyyyy, g_0_z_yyy_xyyyyyz, g_0_z_yyy_xyyyyzz, g_0_z_yyy_xyyyzzz, g_0_z_yyy_xyyzzzz, g_0_z_yyy_xyzzzzz, g_0_z_yyy_xzzzzzz, g_0_z_yyy_yyyyyyy, g_0_z_yyy_yyyyyyz, g_0_z_yyy_yyyyyzz, g_0_z_yyy_yyyyzzz, g_0_z_yyy_yyyzzzz, g_0_z_yyy_yyzzzzz, g_0_z_yyy_yzzzzzz, g_0_z_yyy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyy_xxxxxxx[k] = -g_0_z_yy_xxxxxxx[k] * ab_y + g_0_z_yy_xxxxxxxy[k];

                g_0_z_yyy_xxxxxxy[k] = -g_0_z_yy_xxxxxxy[k] * ab_y + g_0_z_yy_xxxxxxyy[k];

                g_0_z_yyy_xxxxxxz[k] = -g_0_z_yy_xxxxxxz[k] * ab_y + g_0_z_yy_xxxxxxyz[k];

                g_0_z_yyy_xxxxxyy[k] = -g_0_z_yy_xxxxxyy[k] * ab_y + g_0_z_yy_xxxxxyyy[k];

                g_0_z_yyy_xxxxxyz[k] = -g_0_z_yy_xxxxxyz[k] * ab_y + g_0_z_yy_xxxxxyyz[k];

                g_0_z_yyy_xxxxxzz[k] = -g_0_z_yy_xxxxxzz[k] * ab_y + g_0_z_yy_xxxxxyzz[k];

                g_0_z_yyy_xxxxyyy[k] = -g_0_z_yy_xxxxyyy[k] * ab_y + g_0_z_yy_xxxxyyyy[k];

                g_0_z_yyy_xxxxyyz[k] = -g_0_z_yy_xxxxyyz[k] * ab_y + g_0_z_yy_xxxxyyyz[k];

                g_0_z_yyy_xxxxyzz[k] = -g_0_z_yy_xxxxyzz[k] * ab_y + g_0_z_yy_xxxxyyzz[k];

                g_0_z_yyy_xxxxzzz[k] = -g_0_z_yy_xxxxzzz[k] * ab_y + g_0_z_yy_xxxxyzzz[k];

                g_0_z_yyy_xxxyyyy[k] = -g_0_z_yy_xxxyyyy[k] * ab_y + g_0_z_yy_xxxyyyyy[k];

                g_0_z_yyy_xxxyyyz[k] = -g_0_z_yy_xxxyyyz[k] * ab_y + g_0_z_yy_xxxyyyyz[k];

                g_0_z_yyy_xxxyyzz[k] = -g_0_z_yy_xxxyyzz[k] * ab_y + g_0_z_yy_xxxyyyzz[k];

                g_0_z_yyy_xxxyzzz[k] = -g_0_z_yy_xxxyzzz[k] * ab_y + g_0_z_yy_xxxyyzzz[k];

                g_0_z_yyy_xxxzzzz[k] = -g_0_z_yy_xxxzzzz[k] * ab_y + g_0_z_yy_xxxyzzzz[k];

                g_0_z_yyy_xxyyyyy[k] = -g_0_z_yy_xxyyyyy[k] * ab_y + g_0_z_yy_xxyyyyyy[k];

                g_0_z_yyy_xxyyyyz[k] = -g_0_z_yy_xxyyyyz[k] * ab_y + g_0_z_yy_xxyyyyyz[k];

                g_0_z_yyy_xxyyyzz[k] = -g_0_z_yy_xxyyyzz[k] * ab_y + g_0_z_yy_xxyyyyzz[k];

                g_0_z_yyy_xxyyzzz[k] = -g_0_z_yy_xxyyzzz[k] * ab_y + g_0_z_yy_xxyyyzzz[k];

                g_0_z_yyy_xxyzzzz[k] = -g_0_z_yy_xxyzzzz[k] * ab_y + g_0_z_yy_xxyyzzzz[k];

                g_0_z_yyy_xxzzzzz[k] = -g_0_z_yy_xxzzzzz[k] * ab_y + g_0_z_yy_xxyzzzzz[k];

                g_0_z_yyy_xyyyyyy[k] = -g_0_z_yy_xyyyyyy[k] * ab_y + g_0_z_yy_xyyyyyyy[k];

                g_0_z_yyy_xyyyyyz[k] = -g_0_z_yy_xyyyyyz[k] * ab_y + g_0_z_yy_xyyyyyyz[k];

                g_0_z_yyy_xyyyyzz[k] = -g_0_z_yy_xyyyyzz[k] * ab_y + g_0_z_yy_xyyyyyzz[k];

                g_0_z_yyy_xyyyzzz[k] = -g_0_z_yy_xyyyzzz[k] * ab_y + g_0_z_yy_xyyyyzzz[k];

                g_0_z_yyy_xyyzzzz[k] = -g_0_z_yy_xyyzzzz[k] * ab_y + g_0_z_yy_xyyyzzzz[k];

                g_0_z_yyy_xyzzzzz[k] = -g_0_z_yy_xyzzzzz[k] * ab_y + g_0_z_yy_xyyzzzzz[k];

                g_0_z_yyy_xzzzzzz[k] = -g_0_z_yy_xzzzzzz[k] * ab_y + g_0_z_yy_xyzzzzzz[k];

                g_0_z_yyy_yyyyyyy[k] = -g_0_z_yy_yyyyyyy[k] * ab_y + g_0_z_yy_yyyyyyyy[k];

                g_0_z_yyy_yyyyyyz[k] = -g_0_z_yy_yyyyyyz[k] * ab_y + g_0_z_yy_yyyyyyyz[k];

                g_0_z_yyy_yyyyyzz[k] = -g_0_z_yy_yyyyyzz[k] * ab_y + g_0_z_yy_yyyyyyzz[k];

                g_0_z_yyy_yyyyzzz[k] = -g_0_z_yy_yyyyzzz[k] * ab_y + g_0_z_yy_yyyyyzzz[k];

                g_0_z_yyy_yyyzzzz[k] = -g_0_z_yy_yyyzzzz[k] * ab_y + g_0_z_yy_yyyyzzzz[k];

                g_0_z_yyy_yyzzzzz[k] = -g_0_z_yy_yyzzzzz[k] * ab_y + g_0_z_yy_yyyzzzzz[k];

                g_0_z_yyy_yzzzzzz[k] = -g_0_z_yy_yzzzzzz[k] * ab_y + g_0_z_yy_yyzzzzzz[k];

                g_0_z_yyy_zzzzzzz[k] = -g_0_z_yy_zzzzzzz[k] * ab_y + g_0_z_yy_yzzzzzzz[k];
            }

            /// Set up 972-1008 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyz_xxxxxxx, g_0_z_yyz_xxxxxxy, g_0_z_yyz_xxxxxxz, g_0_z_yyz_xxxxxyy, g_0_z_yyz_xxxxxyz, g_0_z_yyz_xxxxxzz, g_0_z_yyz_xxxxyyy, g_0_z_yyz_xxxxyyz, g_0_z_yyz_xxxxyzz, g_0_z_yyz_xxxxzzz, g_0_z_yyz_xxxyyyy, g_0_z_yyz_xxxyyyz, g_0_z_yyz_xxxyyzz, g_0_z_yyz_xxxyzzz, g_0_z_yyz_xxxzzzz, g_0_z_yyz_xxyyyyy, g_0_z_yyz_xxyyyyz, g_0_z_yyz_xxyyyzz, g_0_z_yyz_xxyyzzz, g_0_z_yyz_xxyzzzz, g_0_z_yyz_xxzzzzz, g_0_z_yyz_xyyyyyy, g_0_z_yyz_xyyyyyz, g_0_z_yyz_xyyyyzz, g_0_z_yyz_xyyyzzz, g_0_z_yyz_xyyzzzz, g_0_z_yyz_xyzzzzz, g_0_z_yyz_xzzzzzz, g_0_z_yyz_yyyyyyy, g_0_z_yyz_yyyyyyz, g_0_z_yyz_yyyyyzz, g_0_z_yyz_yyyyzzz, g_0_z_yyz_yyyzzzz, g_0_z_yyz_yyzzzzz, g_0_z_yyz_yzzzzzz, g_0_z_yyz_zzzzzzz, g_0_z_yz_xxxxxxx, g_0_z_yz_xxxxxxxy, g_0_z_yz_xxxxxxy, g_0_z_yz_xxxxxxyy, g_0_z_yz_xxxxxxyz, g_0_z_yz_xxxxxxz, g_0_z_yz_xxxxxyy, g_0_z_yz_xxxxxyyy, g_0_z_yz_xxxxxyyz, g_0_z_yz_xxxxxyz, g_0_z_yz_xxxxxyzz, g_0_z_yz_xxxxxzz, g_0_z_yz_xxxxyyy, g_0_z_yz_xxxxyyyy, g_0_z_yz_xxxxyyyz, g_0_z_yz_xxxxyyz, g_0_z_yz_xxxxyyzz, g_0_z_yz_xxxxyzz, g_0_z_yz_xxxxyzzz, g_0_z_yz_xxxxzzz, g_0_z_yz_xxxyyyy, g_0_z_yz_xxxyyyyy, g_0_z_yz_xxxyyyyz, g_0_z_yz_xxxyyyz, g_0_z_yz_xxxyyyzz, g_0_z_yz_xxxyyzz, g_0_z_yz_xxxyyzzz, g_0_z_yz_xxxyzzz, g_0_z_yz_xxxyzzzz, g_0_z_yz_xxxzzzz, g_0_z_yz_xxyyyyy, g_0_z_yz_xxyyyyyy, g_0_z_yz_xxyyyyyz, g_0_z_yz_xxyyyyz, g_0_z_yz_xxyyyyzz, g_0_z_yz_xxyyyzz, g_0_z_yz_xxyyyzzz, g_0_z_yz_xxyyzzz, g_0_z_yz_xxyyzzzz, g_0_z_yz_xxyzzzz, g_0_z_yz_xxyzzzzz, g_0_z_yz_xxzzzzz, g_0_z_yz_xyyyyyy, g_0_z_yz_xyyyyyyy, g_0_z_yz_xyyyyyyz, g_0_z_yz_xyyyyyz, g_0_z_yz_xyyyyyzz, g_0_z_yz_xyyyyzz, g_0_z_yz_xyyyyzzz, g_0_z_yz_xyyyzzz, g_0_z_yz_xyyyzzzz, g_0_z_yz_xyyzzzz, g_0_z_yz_xyyzzzzz, g_0_z_yz_xyzzzzz, g_0_z_yz_xyzzzzzz, g_0_z_yz_xzzzzzz, g_0_z_yz_yyyyyyy, g_0_z_yz_yyyyyyyy, g_0_z_yz_yyyyyyyz, g_0_z_yz_yyyyyyz, g_0_z_yz_yyyyyyzz, g_0_z_yz_yyyyyzz, g_0_z_yz_yyyyyzzz, g_0_z_yz_yyyyzzz, g_0_z_yz_yyyyzzzz, g_0_z_yz_yyyzzzz, g_0_z_yz_yyyzzzzz, g_0_z_yz_yyzzzzz, g_0_z_yz_yyzzzzzz, g_0_z_yz_yzzzzzz, g_0_z_yz_yzzzzzzz, g_0_z_yz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyz_xxxxxxx[k] = -g_0_z_yz_xxxxxxx[k] * ab_y + g_0_z_yz_xxxxxxxy[k];

                g_0_z_yyz_xxxxxxy[k] = -g_0_z_yz_xxxxxxy[k] * ab_y + g_0_z_yz_xxxxxxyy[k];

                g_0_z_yyz_xxxxxxz[k] = -g_0_z_yz_xxxxxxz[k] * ab_y + g_0_z_yz_xxxxxxyz[k];

                g_0_z_yyz_xxxxxyy[k] = -g_0_z_yz_xxxxxyy[k] * ab_y + g_0_z_yz_xxxxxyyy[k];

                g_0_z_yyz_xxxxxyz[k] = -g_0_z_yz_xxxxxyz[k] * ab_y + g_0_z_yz_xxxxxyyz[k];

                g_0_z_yyz_xxxxxzz[k] = -g_0_z_yz_xxxxxzz[k] * ab_y + g_0_z_yz_xxxxxyzz[k];

                g_0_z_yyz_xxxxyyy[k] = -g_0_z_yz_xxxxyyy[k] * ab_y + g_0_z_yz_xxxxyyyy[k];

                g_0_z_yyz_xxxxyyz[k] = -g_0_z_yz_xxxxyyz[k] * ab_y + g_0_z_yz_xxxxyyyz[k];

                g_0_z_yyz_xxxxyzz[k] = -g_0_z_yz_xxxxyzz[k] * ab_y + g_0_z_yz_xxxxyyzz[k];

                g_0_z_yyz_xxxxzzz[k] = -g_0_z_yz_xxxxzzz[k] * ab_y + g_0_z_yz_xxxxyzzz[k];

                g_0_z_yyz_xxxyyyy[k] = -g_0_z_yz_xxxyyyy[k] * ab_y + g_0_z_yz_xxxyyyyy[k];

                g_0_z_yyz_xxxyyyz[k] = -g_0_z_yz_xxxyyyz[k] * ab_y + g_0_z_yz_xxxyyyyz[k];

                g_0_z_yyz_xxxyyzz[k] = -g_0_z_yz_xxxyyzz[k] * ab_y + g_0_z_yz_xxxyyyzz[k];

                g_0_z_yyz_xxxyzzz[k] = -g_0_z_yz_xxxyzzz[k] * ab_y + g_0_z_yz_xxxyyzzz[k];

                g_0_z_yyz_xxxzzzz[k] = -g_0_z_yz_xxxzzzz[k] * ab_y + g_0_z_yz_xxxyzzzz[k];

                g_0_z_yyz_xxyyyyy[k] = -g_0_z_yz_xxyyyyy[k] * ab_y + g_0_z_yz_xxyyyyyy[k];

                g_0_z_yyz_xxyyyyz[k] = -g_0_z_yz_xxyyyyz[k] * ab_y + g_0_z_yz_xxyyyyyz[k];

                g_0_z_yyz_xxyyyzz[k] = -g_0_z_yz_xxyyyzz[k] * ab_y + g_0_z_yz_xxyyyyzz[k];

                g_0_z_yyz_xxyyzzz[k] = -g_0_z_yz_xxyyzzz[k] * ab_y + g_0_z_yz_xxyyyzzz[k];

                g_0_z_yyz_xxyzzzz[k] = -g_0_z_yz_xxyzzzz[k] * ab_y + g_0_z_yz_xxyyzzzz[k];

                g_0_z_yyz_xxzzzzz[k] = -g_0_z_yz_xxzzzzz[k] * ab_y + g_0_z_yz_xxyzzzzz[k];

                g_0_z_yyz_xyyyyyy[k] = -g_0_z_yz_xyyyyyy[k] * ab_y + g_0_z_yz_xyyyyyyy[k];

                g_0_z_yyz_xyyyyyz[k] = -g_0_z_yz_xyyyyyz[k] * ab_y + g_0_z_yz_xyyyyyyz[k];

                g_0_z_yyz_xyyyyzz[k] = -g_0_z_yz_xyyyyzz[k] * ab_y + g_0_z_yz_xyyyyyzz[k];

                g_0_z_yyz_xyyyzzz[k] = -g_0_z_yz_xyyyzzz[k] * ab_y + g_0_z_yz_xyyyyzzz[k];

                g_0_z_yyz_xyyzzzz[k] = -g_0_z_yz_xyyzzzz[k] * ab_y + g_0_z_yz_xyyyzzzz[k];

                g_0_z_yyz_xyzzzzz[k] = -g_0_z_yz_xyzzzzz[k] * ab_y + g_0_z_yz_xyyzzzzz[k];

                g_0_z_yyz_xzzzzzz[k] = -g_0_z_yz_xzzzzzz[k] * ab_y + g_0_z_yz_xyzzzzzz[k];

                g_0_z_yyz_yyyyyyy[k] = -g_0_z_yz_yyyyyyy[k] * ab_y + g_0_z_yz_yyyyyyyy[k];

                g_0_z_yyz_yyyyyyz[k] = -g_0_z_yz_yyyyyyz[k] * ab_y + g_0_z_yz_yyyyyyyz[k];

                g_0_z_yyz_yyyyyzz[k] = -g_0_z_yz_yyyyyzz[k] * ab_y + g_0_z_yz_yyyyyyzz[k];

                g_0_z_yyz_yyyyzzz[k] = -g_0_z_yz_yyyyzzz[k] * ab_y + g_0_z_yz_yyyyyzzz[k];

                g_0_z_yyz_yyyzzzz[k] = -g_0_z_yz_yyyzzzz[k] * ab_y + g_0_z_yz_yyyyzzzz[k];

                g_0_z_yyz_yyzzzzz[k] = -g_0_z_yz_yyzzzzz[k] * ab_y + g_0_z_yz_yyyzzzzz[k];

                g_0_z_yyz_yzzzzzz[k] = -g_0_z_yz_yzzzzzz[k] * ab_y + g_0_z_yz_yyzzzzzz[k];

                g_0_z_yyz_zzzzzzz[k] = -g_0_z_yz_zzzzzzz[k] * ab_y + g_0_z_yz_yzzzzzzz[k];
            }

            /// Set up 1008-1044 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yzz_xxxxxxx, g_0_z_yzz_xxxxxxy, g_0_z_yzz_xxxxxxz, g_0_z_yzz_xxxxxyy, g_0_z_yzz_xxxxxyz, g_0_z_yzz_xxxxxzz, g_0_z_yzz_xxxxyyy, g_0_z_yzz_xxxxyyz, g_0_z_yzz_xxxxyzz, g_0_z_yzz_xxxxzzz, g_0_z_yzz_xxxyyyy, g_0_z_yzz_xxxyyyz, g_0_z_yzz_xxxyyzz, g_0_z_yzz_xxxyzzz, g_0_z_yzz_xxxzzzz, g_0_z_yzz_xxyyyyy, g_0_z_yzz_xxyyyyz, g_0_z_yzz_xxyyyzz, g_0_z_yzz_xxyyzzz, g_0_z_yzz_xxyzzzz, g_0_z_yzz_xxzzzzz, g_0_z_yzz_xyyyyyy, g_0_z_yzz_xyyyyyz, g_0_z_yzz_xyyyyzz, g_0_z_yzz_xyyyzzz, g_0_z_yzz_xyyzzzz, g_0_z_yzz_xyzzzzz, g_0_z_yzz_xzzzzzz, g_0_z_yzz_yyyyyyy, g_0_z_yzz_yyyyyyz, g_0_z_yzz_yyyyyzz, g_0_z_yzz_yyyyzzz, g_0_z_yzz_yyyzzzz, g_0_z_yzz_yyzzzzz, g_0_z_yzz_yzzzzzz, g_0_z_yzz_zzzzzzz, g_0_z_zz_xxxxxxx, g_0_z_zz_xxxxxxxy, g_0_z_zz_xxxxxxy, g_0_z_zz_xxxxxxyy, g_0_z_zz_xxxxxxyz, g_0_z_zz_xxxxxxz, g_0_z_zz_xxxxxyy, g_0_z_zz_xxxxxyyy, g_0_z_zz_xxxxxyyz, g_0_z_zz_xxxxxyz, g_0_z_zz_xxxxxyzz, g_0_z_zz_xxxxxzz, g_0_z_zz_xxxxyyy, g_0_z_zz_xxxxyyyy, g_0_z_zz_xxxxyyyz, g_0_z_zz_xxxxyyz, g_0_z_zz_xxxxyyzz, g_0_z_zz_xxxxyzz, g_0_z_zz_xxxxyzzz, g_0_z_zz_xxxxzzz, g_0_z_zz_xxxyyyy, g_0_z_zz_xxxyyyyy, g_0_z_zz_xxxyyyyz, g_0_z_zz_xxxyyyz, g_0_z_zz_xxxyyyzz, g_0_z_zz_xxxyyzz, g_0_z_zz_xxxyyzzz, g_0_z_zz_xxxyzzz, g_0_z_zz_xxxyzzzz, g_0_z_zz_xxxzzzz, g_0_z_zz_xxyyyyy, g_0_z_zz_xxyyyyyy, g_0_z_zz_xxyyyyyz, g_0_z_zz_xxyyyyz, g_0_z_zz_xxyyyyzz, g_0_z_zz_xxyyyzz, g_0_z_zz_xxyyyzzz, g_0_z_zz_xxyyzzz, g_0_z_zz_xxyyzzzz, g_0_z_zz_xxyzzzz, g_0_z_zz_xxyzzzzz, g_0_z_zz_xxzzzzz, g_0_z_zz_xyyyyyy, g_0_z_zz_xyyyyyyy, g_0_z_zz_xyyyyyyz, g_0_z_zz_xyyyyyz, g_0_z_zz_xyyyyyzz, g_0_z_zz_xyyyyzz, g_0_z_zz_xyyyyzzz, g_0_z_zz_xyyyzzz, g_0_z_zz_xyyyzzzz, g_0_z_zz_xyyzzzz, g_0_z_zz_xyyzzzzz, g_0_z_zz_xyzzzzz, g_0_z_zz_xyzzzzzz, g_0_z_zz_xzzzzzz, g_0_z_zz_yyyyyyy, g_0_z_zz_yyyyyyyy, g_0_z_zz_yyyyyyyz, g_0_z_zz_yyyyyyz, g_0_z_zz_yyyyyyzz, g_0_z_zz_yyyyyzz, g_0_z_zz_yyyyyzzz, g_0_z_zz_yyyyzzz, g_0_z_zz_yyyyzzzz, g_0_z_zz_yyyzzzz, g_0_z_zz_yyyzzzzz, g_0_z_zz_yyzzzzz, g_0_z_zz_yyzzzzzz, g_0_z_zz_yzzzzzz, g_0_z_zz_yzzzzzzz, g_0_z_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzz_xxxxxxx[k] = -g_0_z_zz_xxxxxxx[k] * ab_y + g_0_z_zz_xxxxxxxy[k];

                g_0_z_yzz_xxxxxxy[k] = -g_0_z_zz_xxxxxxy[k] * ab_y + g_0_z_zz_xxxxxxyy[k];

                g_0_z_yzz_xxxxxxz[k] = -g_0_z_zz_xxxxxxz[k] * ab_y + g_0_z_zz_xxxxxxyz[k];

                g_0_z_yzz_xxxxxyy[k] = -g_0_z_zz_xxxxxyy[k] * ab_y + g_0_z_zz_xxxxxyyy[k];

                g_0_z_yzz_xxxxxyz[k] = -g_0_z_zz_xxxxxyz[k] * ab_y + g_0_z_zz_xxxxxyyz[k];

                g_0_z_yzz_xxxxxzz[k] = -g_0_z_zz_xxxxxzz[k] * ab_y + g_0_z_zz_xxxxxyzz[k];

                g_0_z_yzz_xxxxyyy[k] = -g_0_z_zz_xxxxyyy[k] * ab_y + g_0_z_zz_xxxxyyyy[k];

                g_0_z_yzz_xxxxyyz[k] = -g_0_z_zz_xxxxyyz[k] * ab_y + g_0_z_zz_xxxxyyyz[k];

                g_0_z_yzz_xxxxyzz[k] = -g_0_z_zz_xxxxyzz[k] * ab_y + g_0_z_zz_xxxxyyzz[k];

                g_0_z_yzz_xxxxzzz[k] = -g_0_z_zz_xxxxzzz[k] * ab_y + g_0_z_zz_xxxxyzzz[k];

                g_0_z_yzz_xxxyyyy[k] = -g_0_z_zz_xxxyyyy[k] * ab_y + g_0_z_zz_xxxyyyyy[k];

                g_0_z_yzz_xxxyyyz[k] = -g_0_z_zz_xxxyyyz[k] * ab_y + g_0_z_zz_xxxyyyyz[k];

                g_0_z_yzz_xxxyyzz[k] = -g_0_z_zz_xxxyyzz[k] * ab_y + g_0_z_zz_xxxyyyzz[k];

                g_0_z_yzz_xxxyzzz[k] = -g_0_z_zz_xxxyzzz[k] * ab_y + g_0_z_zz_xxxyyzzz[k];

                g_0_z_yzz_xxxzzzz[k] = -g_0_z_zz_xxxzzzz[k] * ab_y + g_0_z_zz_xxxyzzzz[k];

                g_0_z_yzz_xxyyyyy[k] = -g_0_z_zz_xxyyyyy[k] * ab_y + g_0_z_zz_xxyyyyyy[k];

                g_0_z_yzz_xxyyyyz[k] = -g_0_z_zz_xxyyyyz[k] * ab_y + g_0_z_zz_xxyyyyyz[k];

                g_0_z_yzz_xxyyyzz[k] = -g_0_z_zz_xxyyyzz[k] * ab_y + g_0_z_zz_xxyyyyzz[k];

                g_0_z_yzz_xxyyzzz[k] = -g_0_z_zz_xxyyzzz[k] * ab_y + g_0_z_zz_xxyyyzzz[k];

                g_0_z_yzz_xxyzzzz[k] = -g_0_z_zz_xxyzzzz[k] * ab_y + g_0_z_zz_xxyyzzzz[k];

                g_0_z_yzz_xxzzzzz[k] = -g_0_z_zz_xxzzzzz[k] * ab_y + g_0_z_zz_xxyzzzzz[k];

                g_0_z_yzz_xyyyyyy[k] = -g_0_z_zz_xyyyyyy[k] * ab_y + g_0_z_zz_xyyyyyyy[k];

                g_0_z_yzz_xyyyyyz[k] = -g_0_z_zz_xyyyyyz[k] * ab_y + g_0_z_zz_xyyyyyyz[k];

                g_0_z_yzz_xyyyyzz[k] = -g_0_z_zz_xyyyyzz[k] * ab_y + g_0_z_zz_xyyyyyzz[k];

                g_0_z_yzz_xyyyzzz[k] = -g_0_z_zz_xyyyzzz[k] * ab_y + g_0_z_zz_xyyyyzzz[k];

                g_0_z_yzz_xyyzzzz[k] = -g_0_z_zz_xyyzzzz[k] * ab_y + g_0_z_zz_xyyyzzzz[k];

                g_0_z_yzz_xyzzzzz[k] = -g_0_z_zz_xyzzzzz[k] * ab_y + g_0_z_zz_xyyzzzzz[k];

                g_0_z_yzz_xzzzzzz[k] = -g_0_z_zz_xzzzzzz[k] * ab_y + g_0_z_zz_xyzzzzzz[k];

                g_0_z_yzz_yyyyyyy[k] = -g_0_z_zz_yyyyyyy[k] * ab_y + g_0_z_zz_yyyyyyyy[k];

                g_0_z_yzz_yyyyyyz[k] = -g_0_z_zz_yyyyyyz[k] * ab_y + g_0_z_zz_yyyyyyyz[k];

                g_0_z_yzz_yyyyyzz[k] = -g_0_z_zz_yyyyyzz[k] * ab_y + g_0_z_zz_yyyyyyzz[k];

                g_0_z_yzz_yyyyzzz[k] = -g_0_z_zz_yyyyzzz[k] * ab_y + g_0_z_zz_yyyyyzzz[k];

                g_0_z_yzz_yyyzzzz[k] = -g_0_z_zz_yyyzzzz[k] * ab_y + g_0_z_zz_yyyyzzzz[k];

                g_0_z_yzz_yyzzzzz[k] = -g_0_z_zz_yyzzzzz[k] * ab_y + g_0_z_zz_yyyzzzzz[k];

                g_0_z_yzz_yzzzzzz[k] = -g_0_z_zz_yzzzzzz[k] * ab_y + g_0_z_zz_yyzzzzzz[k];

                g_0_z_yzz_zzzzzzz[k] = -g_0_z_zz_zzzzzzz[k] * ab_y + g_0_z_zz_yzzzzzzz[k];
            }

            /// Set up 1044-1080 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_zz_xxxxxxx, g_0_z_zz_xxxxxxxz, g_0_z_zz_xxxxxxy, g_0_z_zz_xxxxxxyz, g_0_z_zz_xxxxxxz, g_0_z_zz_xxxxxxzz, g_0_z_zz_xxxxxyy, g_0_z_zz_xxxxxyyz, g_0_z_zz_xxxxxyz, g_0_z_zz_xxxxxyzz, g_0_z_zz_xxxxxzz, g_0_z_zz_xxxxxzzz, g_0_z_zz_xxxxyyy, g_0_z_zz_xxxxyyyz, g_0_z_zz_xxxxyyz, g_0_z_zz_xxxxyyzz, g_0_z_zz_xxxxyzz, g_0_z_zz_xxxxyzzz, g_0_z_zz_xxxxzzz, g_0_z_zz_xxxxzzzz, g_0_z_zz_xxxyyyy, g_0_z_zz_xxxyyyyz, g_0_z_zz_xxxyyyz, g_0_z_zz_xxxyyyzz, g_0_z_zz_xxxyyzz, g_0_z_zz_xxxyyzzz, g_0_z_zz_xxxyzzz, g_0_z_zz_xxxyzzzz, g_0_z_zz_xxxzzzz, g_0_z_zz_xxxzzzzz, g_0_z_zz_xxyyyyy, g_0_z_zz_xxyyyyyz, g_0_z_zz_xxyyyyz, g_0_z_zz_xxyyyyzz, g_0_z_zz_xxyyyzz, g_0_z_zz_xxyyyzzz, g_0_z_zz_xxyyzzz, g_0_z_zz_xxyyzzzz, g_0_z_zz_xxyzzzz, g_0_z_zz_xxyzzzzz, g_0_z_zz_xxzzzzz, g_0_z_zz_xxzzzzzz, g_0_z_zz_xyyyyyy, g_0_z_zz_xyyyyyyz, g_0_z_zz_xyyyyyz, g_0_z_zz_xyyyyyzz, g_0_z_zz_xyyyyzz, g_0_z_zz_xyyyyzzz, g_0_z_zz_xyyyzzz, g_0_z_zz_xyyyzzzz, g_0_z_zz_xyyzzzz, g_0_z_zz_xyyzzzzz, g_0_z_zz_xyzzzzz, g_0_z_zz_xyzzzzzz, g_0_z_zz_xzzzzzz, g_0_z_zz_xzzzzzzz, g_0_z_zz_yyyyyyy, g_0_z_zz_yyyyyyyz, g_0_z_zz_yyyyyyz, g_0_z_zz_yyyyyyzz, g_0_z_zz_yyyyyzz, g_0_z_zz_yyyyyzzz, g_0_z_zz_yyyyzzz, g_0_z_zz_yyyyzzzz, g_0_z_zz_yyyzzzz, g_0_z_zz_yyyzzzzz, g_0_z_zz_yyzzzzz, g_0_z_zz_yyzzzzzz, g_0_z_zz_yzzzzzz, g_0_z_zz_yzzzzzzz, g_0_z_zz_zzzzzzz, g_0_z_zz_zzzzzzzz, g_0_z_zzz_xxxxxxx, g_0_z_zzz_xxxxxxy, g_0_z_zzz_xxxxxxz, g_0_z_zzz_xxxxxyy, g_0_z_zzz_xxxxxyz, g_0_z_zzz_xxxxxzz, g_0_z_zzz_xxxxyyy, g_0_z_zzz_xxxxyyz, g_0_z_zzz_xxxxyzz, g_0_z_zzz_xxxxzzz, g_0_z_zzz_xxxyyyy, g_0_z_zzz_xxxyyyz, g_0_z_zzz_xxxyyzz, g_0_z_zzz_xxxyzzz, g_0_z_zzz_xxxzzzz, g_0_z_zzz_xxyyyyy, g_0_z_zzz_xxyyyyz, g_0_z_zzz_xxyyyzz, g_0_z_zzz_xxyyzzz, g_0_z_zzz_xxyzzzz, g_0_z_zzz_xxzzzzz, g_0_z_zzz_xyyyyyy, g_0_z_zzz_xyyyyyz, g_0_z_zzz_xyyyyzz, g_0_z_zzz_xyyyzzz, g_0_z_zzz_xyyzzzz, g_0_z_zzz_xyzzzzz, g_0_z_zzz_xzzzzzz, g_0_z_zzz_yyyyyyy, g_0_z_zzz_yyyyyyz, g_0_z_zzz_yyyyyzz, g_0_z_zzz_yyyyzzz, g_0_z_zzz_yyyzzzz, g_0_z_zzz_yyzzzzz, g_0_z_zzz_yzzzzzz, g_0_z_zzz_zzzzzzz, g_zz_xxxxxxx, g_zz_xxxxxxy, g_zz_xxxxxxz, g_zz_xxxxxyy, g_zz_xxxxxyz, g_zz_xxxxxzz, g_zz_xxxxyyy, g_zz_xxxxyyz, g_zz_xxxxyzz, g_zz_xxxxzzz, g_zz_xxxyyyy, g_zz_xxxyyyz, g_zz_xxxyyzz, g_zz_xxxyzzz, g_zz_xxxzzzz, g_zz_xxyyyyy, g_zz_xxyyyyz, g_zz_xxyyyzz, g_zz_xxyyzzz, g_zz_xxyzzzz, g_zz_xxzzzzz, g_zz_xyyyyyy, g_zz_xyyyyyz, g_zz_xyyyyzz, g_zz_xyyyzzz, g_zz_xyyzzzz, g_zz_xyzzzzz, g_zz_xzzzzzz, g_zz_yyyyyyy, g_zz_yyyyyyz, g_zz_yyyyyzz, g_zz_yyyyzzz, g_zz_yyyzzzz, g_zz_yyzzzzz, g_zz_yzzzzzz, g_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzz_xxxxxxx[k] = g_zz_xxxxxxx[k] - g_0_z_zz_xxxxxxx[k] * ab_z + g_0_z_zz_xxxxxxxz[k];

                g_0_z_zzz_xxxxxxy[k] = g_zz_xxxxxxy[k] - g_0_z_zz_xxxxxxy[k] * ab_z + g_0_z_zz_xxxxxxyz[k];

                g_0_z_zzz_xxxxxxz[k] = g_zz_xxxxxxz[k] - g_0_z_zz_xxxxxxz[k] * ab_z + g_0_z_zz_xxxxxxzz[k];

                g_0_z_zzz_xxxxxyy[k] = g_zz_xxxxxyy[k] - g_0_z_zz_xxxxxyy[k] * ab_z + g_0_z_zz_xxxxxyyz[k];

                g_0_z_zzz_xxxxxyz[k] = g_zz_xxxxxyz[k] - g_0_z_zz_xxxxxyz[k] * ab_z + g_0_z_zz_xxxxxyzz[k];

                g_0_z_zzz_xxxxxzz[k] = g_zz_xxxxxzz[k] - g_0_z_zz_xxxxxzz[k] * ab_z + g_0_z_zz_xxxxxzzz[k];

                g_0_z_zzz_xxxxyyy[k] = g_zz_xxxxyyy[k] - g_0_z_zz_xxxxyyy[k] * ab_z + g_0_z_zz_xxxxyyyz[k];

                g_0_z_zzz_xxxxyyz[k] = g_zz_xxxxyyz[k] - g_0_z_zz_xxxxyyz[k] * ab_z + g_0_z_zz_xxxxyyzz[k];

                g_0_z_zzz_xxxxyzz[k] = g_zz_xxxxyzz[k] - g_0_z_zz_xxxxyzz[k] * ab_z + g_0_z_zz_xxxxyzzz[k];

                g_0_z_zzz_xxxxzzz[k] = g_zz_xxxxzzz[k] - g_0_z_zz_xxxxzzz[k] * ab_z + g_0_z_zz_xxxxzzzz[k];

                g_0_z_zzz_xxxyyyy[k] = g_zz_xxxyyyy[k] - g_0_z_zz_xxxyyyy[k] * ab_z + g_0_z_zz_xxxyyyyz[k];

                g_0_z_zzz_xxxyyyz[k] = g_zz_xxxyyyz[k] - g_0_z_zz_xxxyyyz[k] * ab_z + g_0_z_zz_xxxyyyzz[k];

                g_0_z_zzz_xxxyyzz[k] = g_zz_xxxyyzz[k] - g_0_z_zz_xxxyyzz[k] * ab_z + g_0_z_zz_xxxyyzzz[k];

                g_0_z_zzz_xxxyzzz[k] = g_zz_xxxyzzz[k] - g_0_z_zz_xxxyzzz[k] * ab_z + g_0_z_zz_xxxyzzzz[k];

                g_0_z_zzz_xxxzzzz[k] = g_zz_xxxzzzz[k] - g_0_z_zz_xxxzzzz[k] * ab_z + g_0_z_zz_xxxzzzzz[k];

                g_0_z_zzz_xxyyyyy[k] = g_zz_xxyyyyy[k] - g_0_z_zz_xxyyyyy[k] * ab_z + g_0_z_zz_xxyyyyyz[k];

                g_0_z_zzz_xxyyyyz[k] = g_zz_xxyyyyz[k] - g_0_z_zz_xxyyyyz[k] * ab_z + g_0_z_zz_xxyyyyzz[k];

                g_0_z_zzz_xxyyyzz[k] = g_zz_xxyyyzz[k] - g_0_z_zz_xxyyyzz[k] * ab_z + g_0_z_zz_xxyyyzzz[k];

                g_0_z_zzz_xxyyzzz[k] = g_zz_xxyyzzz[k] - g_0_z_zz_xxyyzzz[k] * ab_z + g_0_z_zz_xxyyzzzz[k];

                g_0_z_zzz_xxyzzzz[k] = g_zz_xxyzzzz[k] - g_0_z_zz_xxyzzzz[k] * ab_z + g_0_z_zz_xxyzzzzz[k];

                g_0_z_zzz_xxzzzzz[k] = g_zz_xxzzzzz[k] - g_0_z_zz_xxzzzzz[k] * ab_z + g_0_z_zz_xxzzzzzz[k];

                g_0_z_zzz_xyyyyyy[k] = g_zz_xyyyyyy[k] - g_0_z_zz_xyyyyyy[k] * ab_z + g_0_z_zz_xyyyyyyz[k];

                g_0_z_zzz_xyyyyyz[k] = g_zz_xyyyyyz[k] - g_0_z_zz_xyyyyyz[k] * ab_z + g_0_z_zz_xyyyyyzz[k];

                g_0_z_zzz_xyyyyzz[k] = g_zz_xyyyyzz[k] - g_0_z_zz_xyyyyzz[k] * ab_z + g_0_z_zz_xyyyyzzz[k];

                g_0_z_zzz_xyyyzzz[k] = g_zz_xyyyzzz[k] - g_0_z_zz_xyyyzzz[k] * ab_z + g_0_z_zz_xyyyzzzz[k];

                g_0_z_zzz_xyyzzzz[k] = g_zz_xyyzzzz[k] - g_0_z_zz_xyyzzzz[k] * ab_z + g_0_z_zz_xyyzzzzz[k];

                g_0_z_zzz_xyzzzzz[k] = g_zz_xyzzzzz[k] - g_0_z_zz_xyzzzzz[k] * ab_z + g_0_z_zz_xyzzzzzz[k];

                g_0_z_zzz_xzzzzzz[k] = g_zz_xzzzzzz[k] - g_0_z_zz_xzzzzzz[k] * ab_z + g_0_z_zz_xzzzzzzz[k];

                g_0_z_zzz_yyyyyyy[k] = g_zz_yyyyyyy[k] - g_0_z_zz_yyyyyyy[k] * ab_z + g_0_z_zz_yyyyyyyz[k];

                g_0_z_zzz_yyyyyyz[k] = g_zz_yyyyyyz[k] - g_0_z_zz_yyyyyyz[k] * ab_z + g_0_z_zz_yyyyyyzz[k];

                g_0_z_zzz_yyyyyzz[k] = g_zz_yyyyyzz[k] - g_0_z_zz_yyyyyzz[k] * ab_z + g_0_z_zz_yyyyyzzz[k];

                g_0_z_zzz_yyyyzzz[k] = g_zz_yyyyzzz[k] - g_0_z_zz_yyyyzzz[k] * ab_z + g_0_z_zz_yyyyzzzz[k];

                g_0_z_zzz_yyyzzzz[k] = g_zz_yyyzzzz[k] - g_0_z_zz_yyyzzzz[k] * ab_z + g_0_z_zz_yyyzzzzz[k];

                g_0_z_zzz_yyzzzzz[k] = g_zz_yyzzzzz[k] - g_0_z_zz_yyzzzzz[k] * ab_z + g_0_z_zz_yyzzzzzz[k];

                g_0_z_zzz_yzzzzzz[k] = g_zz_yzzzzzz[k] - g_0_z_zz_yzzzzzz[k] * ab_z + g_0_z_zz_yzzzzzzz[k];

                g_0_z_zzz_zzzzzzz[k] = g_zz_zzzzzzz[k] - g_0_z_zz_zzzzzzz[k] * ab_z + g_0_z_zz_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

