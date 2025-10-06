#include "ElectronRepulsionGeom0100ContrRecGHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_ghxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ghxx,
                                            const size_t idx_fhxx,
                                            const size_t idx_geom_01_fhxx,
                                            const size_t idx_geom_01_fixx,
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
            /// Set up components of auxilary buffer : FHSS

            const auto fh_off = idx_fhxx + i * dcomps + j;

            auto g_xxx_xxxxx = cbuffer.data(fh_off + 0 * ccomps * dcomps);

            auto g_xxx_xxxxy = cbuffer.data(fh_off + 1 * ccomps * dcomps);

            auto g_xxx_xxxxz = cbuffer.data(fh_off + 2 * ccomps * dcomps);

            auto g_xxx_xxxyy = cbuffer.data(fh_off + 3 * ccomps * dcomps);

            auto g_xxx_xxxyz = cbuffer.data(fh_off + 4 * ccomps * dcomps);

            auto g_xxx_xxxzz = cbuffer.data(fh_off + 5 * ccomps * dcomps);

            auto g_xxx_xxyyy = cbuffer.data(fh_off + 6 * ccomps * dcomps);

            auto g_xxx_xxyyz = cbuffer.data(fh_off + 7 * ccomps * dcomps);

            auto g_xxx_xxyzz = cbuffer.data(fh_off + 8 * ccomps * dcomps);

            auto g_xxx_xxzzz = cbuffer.data(fh_off + 9 * ccomps * dcomps);

            auto g_xxx_xyyyy = cbuffer.data(fh_off + 10 * ccomps * dcomps);

            auto g_xxx_xyyyz = cbuffer.data(fh_off + 11 * ccomps * dcomps);

            auto g_xxx_xyyzz = cbuffer.data(fh_off + 12 * ccomps * dcomps);

            auto g_xxx_xyzzz = cbuffer.data(fh_off + 13 * ccomps * dcomps);

            auto g_xxx_xzzzz = cbuffer.data(fh_off + 14 * ccomps * dcomps);

            auto g_xxx_yyyyy = cbuffer.data(fh_off + 15 * ccomps * dcomps);

            auto g_xxx_yyyyz = cbuffer.data(fh_off + 16 * ccomps * dcomps);

            auto g_xxx_yyyzz = cbuffer.data(fh_off + 17 * ccomps * dcomps);

            auto g_xxx_yyzzz = cbuffer.data(fh_off + 18 * ccomps * dcomps);

            auto g_xxx_yzzzz = cbuffer.data(fh_off + 19 * ccomps * dcomps);

            auto g_xxx_zzzzz = cbuffer.data(fh_off + 20 * ccomps * dcomps);

            auto g_xxy_xxxxx = cbuffer.data(fh_off + 21 * ccomps * dcomps);

            auto g_xxy_xxxxy = cbuffer.data(fh_off + 22 * ccomps * dcomps);

            auto g_xxy_xxxxz = cbuffer.data(fh_off + 23 * ccomps * dcomps);

            auto g_xxy_xxxyy = cbuffer.data(fh_off + 24 * ccomps * dcomps);

            auto g_xxy_xxxyz = cbuffer.data(fh_off + 25 * ccomps * dcomps);

            auto g_xxy_xxxzz = cbuffer.data(fh_off + 26 * ccomps * dcomps);

            auto g_xxy_xxyyy = cbuffer.data(fh_off + 27 * ccomps * dcomps);

            auto g_xxy_xxyyz = cbuffer.data(fh_off + 28 * ccomps * dcomps);

            auto g_xxy_xxyzz = cbuffer.data(fh_off + 29 * ccomps * dcomps);

            auto g_xxy_xxzzz = cbuffer.data(fh_off + 30 * ccomps * dcomps);

            auto g_xxy_xyyyy = cbuffer.data(fh_off + 31 * ccomps * dcomps);

            auto g_xxy_xyyyz = cbuffer.data(fh_off + 32 * ccomps * dcomps);

            auto g_xxy_xyyzz = cbuffer.data(fh_off + 33 * ccomps * dcomps);

            auto g_xxy_xyzzz = cbuffer.data(fh_off + 34 * ccomps * dcomps);

            auto g_xxy_xzzzz = cbuffer.data(fh_off + 35 * ccomps * dcomps);

            auto g_xxy_yyyyy = cbuffer.data(fh_off + 36 * ccomps * dcomps);

            auto g_xxy_yyyyz = cbuffer.data(fh_off + 37 * ccomps * dcomps);

            auto g_xxy_yyyzz = cbuffer.data(fh_off + 38 * ccomps * dcomps);

            auto g_xxy_yyzzz = cbuffer.data(fh_off + 39 * ccomps * dcomps);

            auto g_xxy_yzzzz = cbuffer.data(fh_off + 40 * ccomps * dcomps);

            auto g_xxy_zzzzz = cbuffer.data(fh_off + 41 * ccomps * dcomps);

            auto g_xxz_xxxxx = cbuffer.data(fh_off + 42 * ccomps * dcomps);

            auto g_xxz_xxxxy = cbuffer.data(fh_off + 43 * ccomps * dcomps);

            auto g_xxz_xxxxz = cbuffer.data(fh_off + 44 * ccomps * dcomps);

            auto g_xxz_xxxyy = cbuffer.data(fh_off + 45 * ccomps * dcomps);

            auto g_xxz_xxxyz = cbuffer.data(fh_off + 46 * ccomps * dcomps);

            auto g_xxz_xxxzz = cbuffer.data(fh_off + 47 * ccomps * dcomps);

            auto g_xxz_xxyyy = cbuffer.data(fh_off + 48 * ccomps * dcomps);

            auto g_xxz_xxyyz = cbuffer.data(fh_off + 49 * ccomps * dcomps);

            auto g_xxz_xxyzz = cbuffer.data(fh_off + 50 * ccomps * dcomps);

            auto g_xxz_xxzzz = cbuffer.data(fh_off + 51 * ccomps * dcomps);

            auto g_xxz_xyyyy = cbuffer.data(fh_off + 52 * ccomps * dcomps);

            auto g_xxz_xyyyz = cbuffer.data(fh_off + 53 * ccomps * dcomps);

            auto g_xxz_xyyzz = cbuffer.data(fh_off + 54 * ccomps * dcomps);

            auto g_xxz_xyzzz = cbuffer.data(fh_off + 55 * ccomps * dcomps);

            auto g_xxz_xzzzz = cbuffer.data(fh_off + 56 * ccomps * dcomps);

            auto g_xxz_yyyyy = cbuffer.data(fh_off + 57 * ccomps * dcomps);

            auto g_xxz_yyyyz = cbuffer.data(fh_off + 58 * ccomps * dcomps);

            auto g_xxz_yyyzz = cbuffer.data(fh_off + 59 * ccomps * dcomps);

            auto g_xxz_yyzzz = cbuffer.data(fh_off + 60 * ccomps * dcomps);

            auto g_xxz_yzzzz = cbuffer.data(fh_off + 61 * ccomps * dcomps);

            auto g_xxz_zzzzz = cbuffer.data(fh_off + 62 * ccomps * dcomps);

            auto g_xyy_xxxxx = cbuffer.data(fh_off + 63 * ccomps * dcomps);

            auto g_xyy_xxxxy = cbuffer.data(fh_off + 64 * ccomps * dcomps);

            auto g_xyy_xxxxz = cbuffer.data(fh_off + 65 * ccomps * dcomps);

            auto g_xyy_xxxyy = cbuffer.data(fh_off + 66 * ccomps * dcomps);

            auto g_xyy_xxxyz = cbuffer.data(fh_off + 67 * ccomps * dcomps);

            auto g_xyy_xxxzz = cbuffer.data(fh_off + 68 * ccomps * dcomps);

            auto g_xyy_xxyyy = cbuffer.data(fh_off + 69 * ccomps * dcomps);

            auto g_xyy_xxyyz = cbuffer.data(fh_off + 70 * ccomps * dcomps);

            auto g_xyy_xxyzz = cbuffer.data(fh_off + 71 * ccomps * dcomps);

            auto g_xyy_xxzzz = cbuffer.data(fh_off + 72 * ccomps * dcomps);

            auto g_xyy_xyyyy = cbuffer.data(fh_off + 73 * ccomps * dcomps);

            auto g_xyy_xyyyz = cbuffer.data(fh_off + 74 * ccomps * dcomps);

            auto g_xyy_xyyzz = cbuffer.data(fh_off + 75 * ccomps * dcomps);

            auto g_xyy_xyzzz = cbuffer.data(fh_off + 76 * ccomps * dcomps);

            auto g_xyy_xzzzz = cbuffer.data(fh_off + 77 * ccomps * dcomps);

            auto g_xyy_yyyyy = cbuffer.data(fh_off + 78 * ccomps * dcomps);

            auto g_xyy_yyyyz = cbuffer.data(fh_off + 79 * ccomps * dcomps);

            auto g_xyy_yyyzz = cbuffer.data(fh_off + 80 * ccomps * dcomps);

            auto g_xyy_yyzzz = cbuffer.data(fh_off + 81 * ccomps * dcomps);

            auto g_xyy_yzzzz = cbuffer.data(fh_off + 82 * ccomps * dcomps);

            auto g_xyy_zzzzz = cbuffer.data(fh_off + 83 * ccomps * dcomps);

            auto g_xyz_xxxxx = cbuffer.data(fh_off + 84 * ccomps * dcomps);

            auto g_xyz_xxxxy = cbuffer.data(fh_off + 85 * ccomps * dcomps);

            auto g_xyz_xxxxz = cbuffer.data(fh_off + 86 * ccomps * dcomps);

            auto g_xyz_xxxyy = cbuffer.data(fh_off + 87 * ccomps * dcomps);

            auto g_xyz_xxxyz = cbuffer.data(fh_off + 88 * ccomps * dcomps);

            auto g_xyz_xxxzz = cbuffer.data(fh_off + 89 * ccomps * dcomps);

            auto g_xyz_xxyyy = cbuffer.data(fh_off + 90 * ccomps * dcomps);

            auto g_xyz_xxyyz = cbuffer.data(fh_off + 91 * ccomps * dcomps);

            auto g_xyz_xxyzz = cbuffer.data(fh_off + 92 * ccomps * dcomps);

            auto g_xyz_xxzzz = cbuffer.data(fh_off + 93 * ccomps * dcomps);

            auto g_xyz_xyyyy = cbuffer.data(fh_off + 94 * ccomps * dcomps);

            auto g_xyz_xyyyz = cbuffer.data(fh_off + 95 * ccomps * dcomps);

            auto g_xyz_xyyzz = cbuffer.data(fh_off + 96 * ccomps * dcomps);

            auto g_xyz_xyzzz = cbuffer.data(fh_off + 97 * ccomps * dcomps);

            auto g_xyz_xzzzz = cbuffer.data(fh_off + 98 * ccomps * dcomps);

            auto g_xyz_yyyyy = cbuffer.data(fh_off + 99 * ccomps * dcomps);

            auto g_xyz_yyyyz = cbuffer.data(fh_off + 100 * ccomps * dcomps);

            auto g_xyz_yyyzz = cbuffer.data(fh_off + 101 * ccomps * dcomps);

            auto g_xyz_yyzzz = cbuffer.data(fh_off + 102 * ccomps * dcomps);

            auto g_xyz_yzzzz = cbuffer.data(fh_off + 103 * ccomps * dcomps);

            auto g_xyz_zzzzz = cbuffer.data(fh_off + 104 * ccomps * dcomps);

            auto g_xzz_xxxxx = cbuffer.data(fh_off + 105 * ccomps * dcomps);

            auto g_xzz_xxxxy = cbuffer.data(fh_off + 106 * ccomps * dcomps);

            auto g_xzz_xxxxz = cbuffer.data(fh_off + 107 * ccomps * dcomps);

            auto g_xzz_xxxyy = cbuffer.data(fh_off + 108 * ccomps * dcomps);

            auto g_xzz_xxxyz = cbuffer.data(fh_off + 109 * ccomps * dcomps);

            auto g_xzz_xxxzz = cbuffer.data(fh_off + 110 * ccomps * dcomps);

            auto g_xzz_xxyyy = cbuffer.data(fh_off + 111 * ccomps * dcomps);

            auto g_xzz_xxyyz = cbuffer.data(fh_off + 112 * ccomps * dcomps);

            auto g_xzz_xxyzz = cbuffer.data(fh_off + 113 * ccomps * dcomps);

            auto g_xzz_xxzzz = cbuffer.data(fh_off + 114 * ccomps * dcomps);

            auto g_xzz_xyyyy = cbuffer.data(fh_off + 115 * ccomps * dcomps);

            auto g_xzz_xyyyz = cbuffer.data(fh_off + 116 * ccomps * dcomps);

            auto g_xzz_xyyzz = cbuffer.data(fh_off + 117 * ccomps * dcomps);

            auto g_xzz_xyzzz = cbuffer.data(fh_off + 118 * ccomps * dcomps);

            auto g_xzz_xzzzz = cbuffer.data(fh_off + 119 * ccomps * dcomps);

            auto g_xzz_yyyyy = cbuffer.data(fh_off + 120 * ccomps * dcomps);

            auto g_xzz_yyyyz = cbuffer.data(fh_off + 121 * ccomps * dcomps);

            auto g_xzz_yyyzz = cbuffer.data(fh_off + 122 * ccomps * dcomps);

            auto g_xzz_yyzzz = cbuffer.data(fh_off + 123 * ccomps * dcomps);

            auto g_xzz_yzzzz = cbuffer.data(fh_off + 124 * ccomps * dcomps);

            auto g_xzz_zzzzz = cbuffer.data(fh_off + 125 * ccomps * dcomps);

            auto g_yyy_xxxxx = cbuffer.data(fh_off + 126 * ccomps * dcomps);

            auto g_yyy_xxxxy = cbuffer.data(fh_off + 127 * ccomps * dcomps);

            auto g_yyy_xxxxz = cbuffer.data(fh_off + 128 * ccomps * dcomps);

            auto g_yyy_xxxyy = cbuffer.data(fh_off + 129 * ccomps * dcomps);

            auto g_yyy_xxxyz = cbuffer.data(fh_off + 130 * ccomps * dcomps);

            auto g_yyy_xxxzz = cbuffer.data(fh_off + 131 * ccomps * dcomps);

            auto g_yyy_xxyyy = cbuffer.data(fh_off + 132 * ccomps * dcomps);

            auto g_yyy_xxyyz = cbuffer.data(fh_off + 133 * ccomps * dcomps);

            auto g_yyy_xxyzz = cbuffer.data(fh_off + 134 * ccomps * dcomps);

            auto g_yyy_xxzzz = cbuffer.data(fh_off + 135 * ccomps * dcomps);

            auto g_yyy_xyyyy = cbuffer.data(fh_off + 136 * ccomps * dcomps);

            auto g_yyy_xyyyz = cbuffer.data(fh_off + 137 * ccomps * dcomps);

            auto g_yyy_xyyzz = cbuffer.data(fh_off + 138 * ccomps * dcomps);

            auto g_yyy_xyzzz = cbuffer.data(fh_off + 139 * ccomps * dcomps);

            auto g_yyy_xzzzz = cbuffer.data(fh_off + 140 * ccomps * dcomps);

            auto g_yyy_yyyyy = cbuffer.data(fh_off + 141 * ccomps * dcomps);

            auto g_yyy_yyyyz = cbuffer.data(fh_off + 142 * ccomps * dcomps);

            auto g_yyy_yyyzz = cbuffer.data(fh_off + 143 * ccomps * dcomps);

            auto g_yyy_yyzzz = cbuffer.data(fh_off + 144 * ccomps * dcomps);

            auto g_yyy_yzzzz = cbuffer.data(fh_off + 145 * ccomps * dcomps);

            auto g_yyy_zzzzz = cbuffer.data(fh_off + 146 * ccomps * dcomps);

            auto g_yyz_xxxxx = cbuffer.data(fh_off + 147 * ccomps * dcomps);

            auto g_yyz_xxxxy = cbuffer.data(fh_off + 148 * ccomps * dcomps);

            auto g_yyz_xxxxz = cbuffer.data(fh_off + 149 * ccomps * dcomps);

            auto g_yyz_xxxyy = cbuffer.data(fh_off + 150 * ccomps * dcomps);

            auto g_yyz_xxxyz = cbuffer.data(fh_off + 151 * ccomps * dcomps);

            auto g_yyz_xxxzz = cbuffer.data(fh_off + 152 * ccomps * dcomps);

            auto g_yyz_xxyyy = cbuffer.data(fh_off + 153 * ccomps * dcomps);

            auto g_yyz_xxyyz = cbuffer.data(fh_off + 154 * ccomps * dcomps);

            auto g_yyz_xxyzz = cbuffer.data(fh_off + 155 * ccomps * dcomps);

            auto g_yyz_xxzzz = cbuffer.data(fh_off + 156 * ccomps * dcomps);

            auto g_yyz_xyyyy = cbuffer.data(fh_off + 157 * ccomps * dcomps);

            auto g_yyz_xyyyz = cbuffer.data(fh_off + 158 * ccomps * dcomps);

            auto g_yyz_xyyzz = cbuffer.data(fh_off + 159 * ccomps * dcomps);

            auto g_yyz_xyzzz = cbuffer.data(fh_off + 160 * ccomps * dcomps);

            auto g_yyz_xzzzz = cbuffer.data(fh_off + 161 * ccomps * dcomps);

            auto g_yyz_yyyyy = cbuffer.data(fh_off + 162 * ccomps * dcomps);

            auto g_yyz_yyyyz = cbuffer.data(fh_off + 163 * ccomps * dcomps);

            auto g_yyz_yyyzz = cbuffer.data(fh_off + 164 * ccomps * dcomps);

            auto g_yyz_yyzzz = cbuffer.data(fh_off + 165 * ccomps * dcomps);

            auto g_yyz_yzzzz = cbuffer.data(fh_off + 166 * ccomps * dcomps);

            auto g_yyz_zzzzz = cbuffer.data(fh_off + 167 * ccomps * dcomps);

            auto g_yzz_xxxxx = cbuffer.data(fh_off + 168 * ccomps * dcomps);

            auto g_yzz_xxxxy = cbuffer.data(fh_off + 169 * ccomps * dcomps);

            auto g_yzz_xxxxz = cbuffer.data(fh_off + 170 * ccomps * dcomps);

            auto g_yzz_xxxyy = cbuffer.data(fh_off + 171 * ccomps * dcomps);

            auto g_yzz_xxxyz = cbuffer.data(fh_off + 172 * ccomps * dcomps);

            auto g_yzz_xxxzz = cbuffer.data(fh_off + 173 * ccomps * dcomps);

            auto g_yzz_xxyyy = cbuffer.data(fh_off + 174 * ccomps * dcomps);

            auto g_yzz_xxyyz = cbuffer.data(fh_off + 175 * ccomps * dcomps);

            auto g_yzz_xxyzz = cbuffer.data(fh_off + 176 * ccomps * dcomps);

            auto g_yzz_xxzzz = cbuffer.data(fh_off + 177 * ccomps * dcomps);

            auto g_yzz_xyyyy = cbuffer.data(fh_off + 178 * ccomps * dcomps);

            auto g_yzz_xyyyz = cbuffer.data(fh_off + 179 * ccomps * dcomps);

            auto g_yzz_xyyzz = cbuffer.data(fh_off + 180 * ccomps * dcomps);

            auto g_yzz_xyzzz = cbuffer.data(fh_off + 181 * ccomps * dcomps);

            auto g_yzz_xzzzz = cbuffer.data(fh_off + 182 * ccomps * dcomps);

            auto g_yzz_yyyyy = cbuffer.data(fh_off + 183 * ccomps * dcomps);

            auto g_yzz_yyyyz = cbuffer.data(fh_off + 184 * ccomps * dcomps);

            auto g_yzz_yyyzz = cbuffer.data(fh_off + 185 * ccomps * dcomps);

            auto g_yzz_yyzzz = cbuffer.data(fh_off + 186 * ccomps * dcomps);

            auto g_yzz_yzzzz = cbuffer.data(fh_off + 187 * ccomps * dcomps);

            auto g_yzz_zzzzz = cbuffer.data(fh_off + 188 * ccomps * dcomps);

            auto g_zzz_xxxxx = cbuffer.data(fh_off + 189 * ccomps * dcomps);

            auto g_zzz_xxxxy = cbuffer.data(fh_off + 190 * ccomps * dcomps);

            auto g_zzz_xxxxz = cbuffer.data(fh_off + 191 * ccomps * dcomps);

            auto g_zzz_xxxyy = cbuffer.data(fh_off + 192 * ccomps * dcomps);

            auto g_zzz_xxxyz = cbuffer.data(fh_off + 193 * ccomps * dcomps);

            auto g_zzz_xxxzz = cbuffer.data(fh_off + 194 * ccomps * dcomps);

            auto g_zzz_xxyyy = cbuffer.data(fh_off + 195 * ccomps * dcomps);

            auto g_zzz_xxyyz = cbuffer.data(fh_off + 196 * ccomps * dcomps);

            auto g_zzz_xxyzz = cbuffer.data(fh_off + 197 * ccomps * dcomps);

            auto g_zzz_xxzzz = cbuffer.data(fh_off + 198 * ccomps * dcomps);

            auto g_zzz_xyyyy = cbuffer.data(fh_off + 199 * ccomps * dcomps);

            auto g_zzz_xyyyz = cbuffer.data(fh_off + 200 * ccomps * dcomps);

            auto g_zzz_xyyzz = cbuffer.data(fh_off + 201 * ccomps * dcomps);

            auto g_zzz_xyzzz = cbuffer.data(fh_off + 202 * ccomps * dcomps);

            auto g_zzz_xzzzz = cbuffer.data(fh_off + 203 * ccomps * dcomps);

            auto g_zzz_yyyyy = cbuffer.data(fh_off + 204 * ccomps * dcomps);

            auto g_zzz_yyyyz = cbuffer.data(fh_off + 205 * ccomps * dcomps);

            auto g_zzz_yyyzz = cbuffer.data(fh_off + 206 * ccomps * dcomps);

            auto g_zzz_yyzzz = cbuffer.data(fh_off + 207 * ccomps * dcomps);

            auto g_zzz_yzzzz = cbuffer.data(fh_off + 208 * ccomps * dcomps);

            auto g_zzz_zzzzz = cbuffer.data(fh_off + 209 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FHSS

            const auto fh_geom_01_off = idx_geom_01_fhxx + i * dcomps + j;

            auto g_0_x_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxx = cbuffer.data(fh_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxy = cbuffer.data(fh_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxz = cbuffer.data(fh_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyy = cbuffer.data(fh_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyz = cbuffer.data(fh_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxx_xxxzz = cbuffer.data(fh_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyy = cbuffer.data(fh_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyz = cbuffer.data(fh_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xxx_xxyzz = cbuffer.data(fh_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xxx_xxzzz = cbuffer.data(fh_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyy = cbuffer.data(fh_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyz = cbuffer.data(fh_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xxx_xyyzz = cbuffer.data(fh_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xxx_xyzzz = cbuffer.data(fh_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xxx_xzzzz = cbuffer.data(fh_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyy = cbuffer.data(fh_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyz = cbuffer.data(fh_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xxx_yyyzz = cbuffer.data(fh_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xxx_yyzzz = cbuffer.data(fh_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xxx_yzzzz = cbuffer.data(fh_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xxx_zzzzz = cbuffer.data(fh_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxx = cbuffer.data(fh_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxy = cbuffer.data(fh_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxz = cbuffer.data(fh_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyy = cbuffer.data(fh_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyz = cbuffer.data(fh_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xxy_xxxzz = cbuffer.data(fh_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyy = cbuffer.data(fh_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyz = cbuffer.data(fh_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xxy_xxyzz = cbuffer.data(fh_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xxy_xxzzz = cbuffer.data(fh_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyy = cbuffer.data(fh_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyz = cbuffer.data(fh_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xxy_xyyzz = cbuffer.data(fh_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xxy_xyzzz = cbuffer.data(fh_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xxy_xzzzz = cbuffer.data(fh_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyy = cbuffer.data(fh_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyz = cbuffer.data(fh_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xxy_yyyzz = cbuffer.data(fh_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xxy_yyzzz = cbuffer.data(fh_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xxy_yzzzz = cbuffer.data(fh_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xxy_zzzzz = cbuffer.data(fh_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxx = cbuffer.data(fh_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxy = cbuffer.data(fh_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxz = cbuffer.data(fh_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyy = cbuffer.data(fh_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyz = cbuffer.data(fh_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xxz_xxxzz = cbuffer.data(fh_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyy = cbuffer.data(fh_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyz = cbuffer.data(fh_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xxz_xxyzz = cbuffer.data(fh_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xxz_xxzzz = cbuffer.data(fh_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyy = cbuffer.data(fh_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyz = cbuffer.data(fh_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xxz_xyyzz = cbuffer.data(fh_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xxz_xyzzz = cbuffer.data(fh_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xxz_xzzzz = cbuffer.data(fh_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyy = cbuffer.data(fh_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyz = cbuffer.data(fh_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xxz_yyyzz = cbuffer.data(fh_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_xxz_yyzzz = cbuffer.data(fh_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xxz_yzzzz = cbuffer.data(fh_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xxz_zzzzz = cbuffer.data(fh_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxx = cbuffer.data(fh_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxy = cbuffer.data(fh_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxz = cbuffer.data(fh_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyy = cbuffer.data(fh_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyz = cbuffer.data(fh_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xyy_xxxzz = cbuffer.data(fh_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyy = cbuffer.data(fh_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyz = cbuffer.data(fh_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xyy_xxyzz = cbuffer.data(fh_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xyy_xxzzz = cbuffer.data(fh_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyy = cbuffer.data(fh_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyz = cbuffer.data(fh_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xyy_xyyzz = cbuffer.data(fh_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xyy_xyzzz = cbuffer.data(fh_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xyy_xzzzz = cbuffer.data(fh_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyy = cbuffer.data(fh_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyz = cbuffer.data(fh_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xyy_yyyzz = cbuffer.data(fh_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xyy_yyzzz = cbuffer.data(fh_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xyy_yzzzz = cbuffer.data(fh_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xyy_zzzzz = cbuffer.data(fh_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxx = cbuffer.data(fh_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxy = cbuffer.data(fh_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxz = cbuffer.data(fh_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyy = cbuffer.data(fh_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyz = cbuffer.data(fh_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xyz_xxxzz = cbuffer.data(fh_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyy = cbuffer.data(fh_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyz = cbuffer.data(fh_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xyz_xxyzz = cbuffer.data(fh_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xyz_xxzzz = cbuffer.data(fh_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyy = cbuffer.data(fh_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyz = cbuffer.data(fh_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xyz_xyyzz = cbuffer.data(fh_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xyz_xyzzz = cbuffer.data(fh_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xyz_xzzzz = cbuffer.data(fh_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyy = cbuffer.data(fh_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyz = cbuffer.data(fh_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xyz_yyyzz = cbuffer.data(fh_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xyz_yyzzz = cbuffer.data(fh_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xyz_yzzzz = cbuffer.data(fh_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xyz_zzzzz = cbuffer.data(fh_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxx = cbuffer.data(fh_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxy = cbuffer.data(fh_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxz = cbuffer.data(fh_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyy = cbuffer.data(fh_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyz = cbuffer.data(fh_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xzz_xxxzz = cbuffer.data(fh_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyy = cbuffer.data(fh_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyz = cbuffer.data(fh_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xzz_xxyzz = cbuffer.data(fh_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xzz_xxzzz = cbuffer.data(fh_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyy = cbuffer.data(fh_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyz = cbuffer.data(fh_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xzz_xyyzz = cbuffer.data(fh_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xzz_xyzzz = cbuffer.data(fh_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xzz_xzzzz = cbuffer.data(fh_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyy = cbuffer.data(fh_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyz = cbuffer.data(fh_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xzz_yyyzz = cbuffer.data(fh_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xzz_yyzzz = cbuffer.data(fh_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xzz_yzzzz = cbuffer.data(fh_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xzz_zzzzz = cbuffer.data(fh_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxx = cbuffer.data(fh_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxy = cbuffer.data(fh_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxz = cbuffer.data(fh_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyy = cbuffer.data(fh_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyz = cbuffer.data(fh_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_yyy_xxxzz = cbuffer.data(fh_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyy = cbuffer.data(fh_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyz = cbuffer.data(fh_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_yyy_xxyzz = cbuffer.data(fh_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_yyy_xxzzz = cbuffer.data(fh_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyy = cbuffer.data(fh_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyz = cbuffer.data(fh_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_yyy_xyyzz = cbuffer.data(fh_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_yyy_xyzzz = cbuffer.data(fh_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_yyy_xzzzz = cbuffer.data(fh_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyy = cbuffer.data(fh_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyz = cbuffer.data(fh_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_yyy_yyyzz = cbuffer.data(fh_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_yyy_yyzzz = cbuffer.data(fh_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_yyy_yzzzz = cbuffer.data(fh_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_yyy_zzzzz = cbuffer.data(fh_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxx = cbuffer.data(fh_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxy = cbuffer.data(fh_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxz = cbuffer.data(fh_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyy = cbuffer.data(fh_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyz = cbuffer.data(fh_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_yyz_xxxzz = cbuffer.data(fh_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyy = cbuffer.data(fh_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyz = cbuffer.data(fh_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_yyz_xxyzz = cbuffer.data(fh_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_yyz_xxzzz = cbuffer.data(fh_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyy = cbuffer.data(fh_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyz = cbuffer.data(fh_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_yyz_xyyzz = cbuffer.data(fh_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_yyz_xyzzz = cbuffer.data(fh_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_yyz_xzzzz = cbuffer.data(fh_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyy = cbuffer.data(fh_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyz = cbuffer.data(fh_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_yyz_yyyzz = cbuffer.data(fh_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_yyz_yyzzz = cbuffer.data(fh_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_yyz_yzzzz = cbuffer.data(fh_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_yyz_zzzzz = cbuffer.data(fh_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxx = cbuffer.data(fh_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxy = cbuffer.data(fh_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxz = cbuffer.data(fh_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyy = cbuffer.data(fh_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyz = cbuffer.data(fh_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_yzz_xxxzz = cbuffer.data(fh_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyy = cbuffer.data(fh_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyz = cbuffer.data(fh_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_yzz_xxyzz = cbuffer.data(fh_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_yzz_xxzzz = cbuffer.data(fh_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyy = cbuffer.data(fh_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyz = cbuffer.data(fh_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_yzz_xyyzz = cbuffer.data(fh_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yzz_xyzzz = cbuffer.data(fh_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yzz_xzzzz = cbuffer.data(fh_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyy = cbuffer.data(fh_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyz = cbuffer.data(fh_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yzz_yyyzz = cbuffer.data(fh_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yzz_yyzzz = cbuffer.data(fh_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yzz_yzzzz = cbuffer.data(fh_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yzz_zzzzz = cbuffer.data(fh_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxx = cbuffer.data(fh_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxy = cbuffer.data(fh_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxz = cbuffer.data(fh_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyy = cbuffer.data(fh_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyz = cbuffer.data(fh_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_zzz_xxxzz = cbuffer.data(fh_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyy = cbuffer.data(fh_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyz = cbuffer.data(fh_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_zzz_xxyzz = cbuffer.data(fh_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_zzz_xxzzz = cbuffer.data(fh_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyy = cbuffer.data(fh_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyz = cbuffer.data(fh_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_zzz_xyyzz = cbuffer.data(fh_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_zzz_xyzzz = cbuffer.data(fh_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_zzz_xzzzz = cbuffer.data(fh_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyy = cbuffer.data(fh_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyz = cbuffer.data(fh_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_zzz_yyyzz = cbuffer.data(fh_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_zzz_yyzzz = cbuffer.data(fh_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_zzz_yzzzz = cbuffer.data(fh_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_zzz_zzzzz = cbuffer.data(fh_geom_01_off + 629 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FISS

            const auto fi_geom_01_off = idx_geom_01_fixx + i * dcomps + j;

            auto g_0_x_xxx_xxxxxx = cbuffer.data(fi_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxy = cbuffer.data(fi_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxxz = cbuffer.data(fi_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyy = cbuffer.data(fi_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxyz = cbuffer.data(fi_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xxxxzz = cbuffer.data(fi_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyy = cbuffer.data(fi_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyyz = cbuffer.data(fi_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_xxxyzz = cbuffer.data(fi_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_xxxzzz = cbuffer.data(fi_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyy = cbuffer.data(fi_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyyz = cbuffer.data(fi_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxx_xxyyzz = cbuffer.data(fi_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxx_xxyzzz = cbuffer.data(fi_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxx_xxzzzz = cbuffer.data(fi_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyy = cbuffer.data(fi_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyyz = cbuffer.data(fi_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxx_xyyyzz = cbuffer.data(fi_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxx_xyyzzz = cbuffer.data(fi_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxx_xyzzzz = cbuffer.data(fi_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxx_xzzzzz = cbuffer.data(fi_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyy = cbuffer.data(fi_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyyz = cbuffer.data(fi_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxx_yyyyzz = cbuffer.data(fi_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxx_yyyzzz = cbuffer.data(fi_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxx_yyzzzz = cbuffer.data(fi_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxx_yzzzzz = cbuffer.data(fi_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxx_zzzzzz = cbuffer.data(fi_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxx = cbuffer.data(fi_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxy = cbuffer.data(fi_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxxz = cbuffer.data(fi_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyy = cbuffer.data(fi_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxyz = cbuffer.data(fi_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxy_xxxxzz = cbuffer.data(fi_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyy = cbuffer.data(fi_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyyz = cbuffer.data(fi_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxy_xxxyzz = cbuffer.data(fi_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxy_xxxzzz = cbuffer.data(fi_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyy = cbuffer.data(fi_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyyz = cbuffer.data(fi_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxy_xxyyzz = cbuffer.data(fi_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxy_xxyzzz = cbuffer.data(fi_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxy_xxzzzz = cbuffer.data(fi_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyy = cbuffer.data(fi_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyyz = cbuffer.data(fi_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxy_xyyyzz = cbuffer.data(fi_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxy_xyyzzz = cbuffer.data(fi_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxy_xyzzzz = cbuffer.data(fi_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxy_xzzzzz = cbuffer.data(fi_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyy = cbuffer.data(fi_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyyz = cbuffer.data(fi_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxy_yyyyzz = cbuffer.data(fi_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxy_yyyzzz = cbuffer.data(fi_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxy_yyzzzz = cbuffer.data(fi_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxy_yzzzzz = cbuffer.data(fi_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxy_zzzzzz = cbuffer.data(fi_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxx = cbuffer.data(fi_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxy = cbuffer.data(fi_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxxz = cbuffer.data(fi_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyy = cbuffer.data(fi_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxyz = cbuffer.data(fi_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxz_xxxxzz = cbuffer.data(fi_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyy = cbuffer.data(fi_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyyz = cbuffer.data(fi_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxz_xxxyzz = cbuffer.data(fi_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxz_xxxzzz = cbuffer.data(fi_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyy = cbuffer.data(fi_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyyz = cbuffer.data(fi_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxz_xxyyzz = cbuffer.data(fi_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxz_xxyzzz = cbuffer.data(fi_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxz_xxzzzz = cbuffer.data(fi_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyy = cbuffer.data(fi_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyyz = cbuffer.data(fi_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxz_xyyyzz = cbuffer.data(fi_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxz_xyyzzz = cbuffer.data(fi_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxz_xyzzzz = cbuffer.data(fi_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxz_xzzzzz = cbuffer.data(fi_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyy = cbuffer.data(fi_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyyz = cbuffer.data(fi_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxz_yyyyzz = cbuffer.data(fi_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxz_yyyzzz = cbuffer.data(fi_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxz_yyzzzz = cbuffer.data(fi_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxz_yzzzzz = cbuffer.data(fi_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxz_zzzzzz = cbuffer.data(fi_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxx = cbuffer.data(fi_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxy = cbuffer.data(fi_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxxz = cbuffer.data(fi_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyy = cbuffer.data(fi_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxyz = cbuffer.data(fi_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xyy_xxxxzz = cbuffer.data(fi_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyy = cbuffer.data(fi_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyyz = cbuffer.data(fi_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyy_xxxyzz = cbuffer.data(fi_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyy_xxxzzz = cbuffer.data(fi_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyy = cbuffer.data(fi_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyyz = cbuffer.data(fi_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyy_xxyyzz = cbuffer.data(fi_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyy_xxyzzz = cbuffer.data(fi_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyy_xxzzzz = cbuffer.data(fi_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyy = cbuffer.data(fi_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyyz = cbuffer.data(fi_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyy_xyyyzz = cbuffer.data(fi_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyy_xyyzzz = cbuffer.data(fi_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyy_xyzzzz = cbuffer.data(fi_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyy_xzzzzz = cbuffer.data(fi_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyy = cbuffer.data(fi_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyyz = cbuffer.data(fi_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyy_yyyyzz = cbuffer.data(fi_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyy_yyyzzz = cbuffer.data(fi_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyy_yyzzzz = cbuffer.data(fi_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyy_yzzzzz = cbuffer.data(fi_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyy_zzzzzz = cbuffer.data(fi_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxx = cbuffer.data(fi_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxy = cbuffer.data(fi_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxxz = cbuffer.data(fi_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyy = cbuffer.data(fi_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxyz = cbuffer.data(fi_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyz_xxxxzz = cbuffer.data(fi_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyy = cbuffer.data(fi_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyyz = cbuffer.data(fi_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xyz_xxxyzz = cbuffer.data(fi_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xyz_xxxzzz = cbuffer.data(fi_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyy = cbuffer.data(fi_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyyz = cbuffer.data(fi_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xyz_xxyyzz = cbuffer.data(fi_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xyz_xxyzzz = cbuffer.data(fi_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xyz_xxzzzz = cbuffer.data(fi_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyy = cbuffer.data(fi_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyyz = cbuffer.data(fi_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyz_xyyyzz = cbuffer.data(fi_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyz_xyyzzz = cbuffer.data(fi_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyz_xyzzzz = cbuffer.data(fi_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xyz_xzzzzz = cbuffer.data(fi_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyy = cbuffer.data(fi_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyyz = cbuffer.data(fi_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyz_yyyyzz = cbuffer.data(fi_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyz_yyyzzz = cbuffer.data(fi_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyz_yyzzzz = cbuffer.data(fi_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xyz_yzzzzz = cbuffer.data(fi_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyz_zzzzzz = cbuffer.data(fi_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxx = cbuffer.data(fi_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxy = cbuffer.data(fi_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxxz = cbuffer.data(fi_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyy = cbuffer.data(fi_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxyz = cbuffer.data(fi_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xzz_xxxxzz = cbuffer.data(fi_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyy = cbuffer.data(fi_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyyz = cbuffer.data(fi_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xzz_xxxyzz = cbuffer.data(fi_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xzz_xxxzzz = cbuffer.data(fi_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyy = cbuffer.data(fi_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyyz = cbuffer.data(fi_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xzz_xxyyzz = cbuffer.data(fi_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xzz_xxyzzz = cbuffer.data(fi_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xzz_xxzzzz = cbuffer.data(fi_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyy = cbuffer.data(fi_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyyz = cbuffer.data(fi_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xzz_xyyyzz = cbuffer.data(fi_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xzz_xyyzzz = cbuffer.data(fi_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xzz_xyzzzz = cbuffer.data(fi_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xzz_xzzzzz = cbuffer.data(fi_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyy = cbuffer.data(fi_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyyz = cbuffer.data(fi_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xzz_yyyyzz = cbuffer.data(fi_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xzz_yyyzzz = cbuffer.data(fi_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xzz_yyzzzz = cbuffer.data(fi_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xzz_yzzzzz = cbuffer.data(fi_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xzz_zzzzzz = cbuffer.data(fi_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxx = cbuffer.data(fi_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxy = cbuffer.data(fi_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxxz = cbuffer.data(fi_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyy = cbuffer.data(fi_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxyz = cbuffer.data(fi_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yyy_xxxxzz = cbuffer.data(fi_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyy = cbuffer.data(fi_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyyz = cbuffer.data(fi_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yyy_xxxyzz = cbuffer.data(fi_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yyy_xxxzzz = cbuffer.data(fi_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyy = cbuffer.data(fi_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyyz = cbuffer.data(fi_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_yyy_xxyyzz = cbuffer.data(fi_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yyy_xxyzzz = cbuffer.data(fi_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yyy_xxzzzz = cbuffer.data(fi_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyy = cbuffer.data(fi_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyyz = cbuffer.data(fi_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yyy_xyyyzz = cbuffer.data(fi_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yyy_xyyzzz = cbuffer.data(fi_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yyy_xyzzzz = cbuffer.data(fi_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yyy_xzzzzz = cbuffer.data(fi_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyy = cbuffer.data(fi_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyyz = cbuffer.data(fi_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_yyy_yyyyzz = cbuffer.data(fi_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_yyy_yyyzzz = cbuffer.data(fi_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_yyy_yyzzzz = cbuffer.data(fi_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_yyy_yzzzzz = cbuffer.data(fi_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_yyy_zzzzzz = cbuffer.data(fi_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxx = cbuffer.data(fi_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxy = cbuffer.data(fi_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxxz = cbuffer.data(fi_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyy = cbuffer.data(fi_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxyz = cbuffer.data(fi_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_yyz_xxxxzz = cbuffer.data(fi_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyy = cbuffer.data(fi_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyyz = cbuffer.data(fi_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_yyz_xxxyzz = cbuffer.data(fi_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_yyz_xxxzzz = cbuffer.data(fi_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyy = cbuffer.data(fi_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyyz = cbuffer.data(fi_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_yyz_xxyyzz = cbuffer.data(fi_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_yyz_xxyzzz = cbuffer.data(fi_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_yyz_xxzzzz = cbuffer.data(fi_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyy = cbuffer.data(fi_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyyz = cbuffer.data(fi_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_yyz_xyyyzz = cbuffer.data(fi_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_yyz_xyyzzz = cbuffer.data(fi_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_yyz_xyzzzz = cbuffer.data(fi_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yyz_xzzzzz = cbuffer.data(fi_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyy = cbuffer.data(fi_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyyz = cbuffer.data(fi_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yyz_yyyyzz = cbuffer.data(fi_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yyz_yyyzzz = cbuffer.data(fi_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yyz_yyzzzz = cbuffer.data(fi_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yyz_yzzzzz = cbuffer.data(fi_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yyz_zzzzzz = cbuffer.data(fi_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxx = cbuffer.data(fi_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxy = cbuffer.data(fi_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxxz = cbuffer.data(fi_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyy = cbuffer.data(fi_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxyz = cbuffer.data(fi_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yzz_xxxxzz = cbuffer.data(fi_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyy = cbuffer.data(fi_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyyz = cbuffer.data(fi_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yzz_xxxyzz = cbuffer.data(fi_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yzz_xxxzzz = cbuffer.data(fi_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyy = cbuffer.data(fi_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyyz = cbuffer.data(fi_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yzz_xxyyzz = cbuffer.data(fi_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yzz_xxyzzz = cbuffer.data(fi_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yzz_xxzzzz = cbuffer.data(fi_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyy = cbuffer.data(fi_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyyz = cbuffer.data(fi_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yzz_xyyyzz = cbuffer.data(fi_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yzz_xyyzzz = cbuffer.data(fi_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yzz_xyzzzz = cbuffer.data(fi_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yzz_xzzzzz = cbuffer.data(fi_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyy = cbuffer.data(fi_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyyz = cbuffer.data(fi_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yzz_yyyyzz = cbuffer.data(fi_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yzz_yyyzzz = cbuffer.data(fi_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yzz_yyzzzz = cbuffer.data(fi_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yzz_yzzzzz = cbuffer.data(fi_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yzz_zzzzzz = cbuffer.data(fi_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxx = cbuffer.data(fi_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxy = cbuffer.data(fi_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxxz = cbuffer.data(fi_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyy = cbuffer.data(fi_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxyz = cbuffer.data(fi_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_zzz_xxxxzz = cbuffer.data(fi_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyy = cbuffer.data(fi_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyyz = cbuffer.data(fi_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_zzz_xxxyzz = cbuffer.data(fi_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_zzz_xxxzzz = cbuffer.data(fi_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyy = cbuffer.data(fi_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyyz = cbuffer.data(fi_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_zzz_xxyyzz = cbuffer.data(fi_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_zzz_xxyzzz = cbuffer.data(fi_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_zzz_xxzzzz = cbuffer.data(fi_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyy = cbuffer.data(fi_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyyz = cbuffer.data(fi_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_zzz_xyyyzz = cbuffer.data(fi_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_zzz_xyyzzz = cbuffer.data(fi_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_zzz_xyzzzz = cbuffer.data(fi_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_zzz_xzzzzz = cbuffer.data(fi_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyy = cbuffer.data(fi_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyyz = cbuffer.data(fi_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_zzz_yyyyzz = cbuffer.data(fi_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_zzz_yyyzzz = cbuffer.data(fi_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_zzz_yyzzzz = cbuffer.data(fi_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_zzz_yzzzzz = cbuffer.data(fi_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_zzz_zzzzzz = cbuffer.data(fi_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxx = cbuffer.data(fi_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxy = cbuffer.data(fi_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxxz = cbuffer.data(fi_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyy = cbuffer.data(fi_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxyz = cbuffer.data(fi_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xxx_xxxxzz = cbuffer.data(fi_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyy = cbuffer.data(fi_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyyz = cbuffer.data(fi_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xxx_xxxyzz = cbuffer.data(fi_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxx_xxxzzz = cbuffer.data(fi_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyy = cbuffer.data(fi_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyyz = cbuffer.data(fi_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxx_xxyyzz = cbuffer.data(fi_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxx_xxyzzz = cbuffer.data(fi_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xxx_xxzzzz = cbuffer.data(fi_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyy = cbuffer.data(fi_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyyz = cbuffer.data(fi_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxx_xyyyzz = cbuffer.data(fi_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxx_xyyzzz = cbuffer.data(fi_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxx_xyzzzz = cbuffer.data(fi_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xxx_xzzzzz = cbuffer.data(fi_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyy = cbuffer.data(fi_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyyz = cbuffer.data(fi_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxx_yyyyzz = cbuffer.data(fi_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxx_yyyzzz = cbuffer.data(fi_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxx_yyzzzz = cbuffer.data(fi_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xxx_yzzzzz = cbuffer.data(fi_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxx_zzzzzz = cbuffer.data(fi_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxx = cbuffer.data(fi_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxy = cbuffer.data(fi_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxxz = cbuffer.data(fi_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyy = cbuffer.data(fi_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxyz = cbuffer.data(fi_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xxy_xxxxzz = cbuffer.data(fi_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyy = cbuffer.data(fi_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyyz = cbuffer.data(fi_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxy_xxxyzz = cbuffer.data(fi_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxy_xxxzzz = cbuffer.data(fi_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyy = cbuffer.data(fi_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyyz = cbuffer.data(fi_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxy_xxyyzz = cbuffer.data(fi_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxy_xxyzzz = cbuffer.data(fi_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxy_xxzzzz = cbuffer.data(fi_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyy = cbuffer.data(fi_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyyz = cbuffer.data(fi_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxy_xyyyzz = cbuffer.data(fi_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxy_xyyzzz = cbuffer.data(fi_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxy_xyzzzz = cbuffer.data(fi_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxy_xzzzzz = cbuffer.data(fi_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyy = cbuffer.data(fi_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyyz = cbuffer.data(fi_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxy_yyyyzz = cbuffer.data(fi_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxy_yyyzzz = cbuffer.data(fi_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxy_yyzzzz = cbuffer.data(fi_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxy_yzzzzz = cbuffer.data(fi_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxy_zzzzzz = cbuffer.data(fi_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxx = cbuffer.data(fi_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxy = cbuffer.data(fi_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxxz = cbuffer.data(fi_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyy = cbuffer.data(fi_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxyz = cbuffer.data(fi_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxz_xxxxzz = cbuffer.data(fi_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyy = cbuffer.data(fi_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyyz = cbuffer.data(fi_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xxz_xxxyzz = cbuffer.data(fi_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xxz_xxxzzz = cbuffer.data(fi_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyy = cbuffer.data(fi_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyyz = cbuffer.data(fi_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xxz_xxyyzz = cbuffer.data(fi_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xxz_xxyzzz = cbuffer.data(fi_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xxz_xxzzzz = cbuffer.data(fi_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyy = cbuffer.data(fi_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyyz = cbuffer.data(fi_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xxz_xyyyzz = cbuffer.data(fi_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xxz_xyyzzz = cbuffer.data(fi_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xxz_xyzzzz = cbuffer.data(fi_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xxz_xzzzzz = cbuffer.data(fi_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyy = cbuffer.data(fi_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyyz = cbuffer.data(fi_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xxz_yyyyzz = cbuffer.data(fi_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxz_yyyzzz = cbuffer.data(fi_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxz_yyzzzz = cbuffer.data(fi_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxz_yzzzzz = cbuffer.data(fi_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxz_zzzzzz = cbuffer.data(fi_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxx = cbuffer.data(fi_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxy = cbuffer.data(fi_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxxz = cbuffer.data(fi_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyy = cbuffer.data(fi_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxyz = cbuffer.data(fi_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xyy_xxxxzz = cbuffer.data(fi_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyy = cbuffer.data(fi_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyyz = cbuffer.data(fi_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xyy_xxxyzz = cbuffer.data(fi_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xyy_xxxzzz = cbuffer.data(fi_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyy = cbuffer.data(fi_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyyz = cbuffer.data(fi_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xyy_xxyyzz = cbuffer.data(fi_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xyy_xxyzzz = cbuffer.data(fi_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xyy_xxzzzz = cbuffer.data(fi_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyy = cbuffer.data(fi_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyyz = cbuffer.data(fi_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xyy_xyyyzz = cbuffer.data(fi_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xyy_xyyzzz = cbuffer.data(fi_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xyy_xyzzzz = cbuffer.data(fi_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xyy_xzzzzz = cbuffer.data(fi_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyy = cbuffer.data(fi_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyyz = cbuffer.data(fi_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xyy_yyyyzz = cbuffer.data(fi_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xyy_yyyzzz = cbuffer.data(fi_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xyy_yyzzzz = cbuffer.data(fi_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xyy_yzzzzz = cbuffer.data(fi_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xyy_zzzzzz = cbuffer.data(fi_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxx = cbuffer.data(fi_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxy = cbuffer.data(fi_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxxz = cbuffer.data(fi_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyy = cbuffer.data(fi_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxyz = cbuffer.data(fi_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xyz_xxxxzz = cbuffer.data(fi_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyy = cbuffer.data(fi_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyyz = cbuffer.data(fi_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xyz_xxxyzz = cbuffer.data(fi_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xyz_xxxzzz = cbuffer.data(fi_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyy = cbuffer.data(fi_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyyz = cbuffer.data(fi_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xyz_xxyyzz = cbuffer.data(fi_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xyz_xxyzzz = cbuffer.data(fi_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xyz_xxzzzz = cbuffer.data(fi_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyy = cbuffer.data(fi_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyyz = cbuffer.data(fi_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xyz_xyyyzz = cbuffer.data(fi_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xyz_xyyzzz = cbuffer.data(fi_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xyz_xyzzzz = cbuffer.data(fi_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xyz_xzzzzz = cbuffer.data(fi_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyy = cbuffer.data(fi_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyyz = cbuffer.data(fi_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xyz_yyyyzz = cbuffer.data(fi_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xyz_yyyzzz = cbuffer.data(fi_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xyz_yyzzzz = cbuffer.data(fi_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xyz_yzzzzz = cbuffer.data(fi_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xyz_zzzzzz = cbuffer.data(fi_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxx = cbuffer.data(fi_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxy = cbuffer.data(fi_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxxz = cbuffer.data(fi_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyy = cbuffer.data(fi_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxyz = cbuffer.data(fi_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xzz_xxxxzz = cbuffer.data(fi_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyy = cbuffer.data(fi_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyyz = cbuffer.data(fi_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xzz_xxxyzz = cbuffer.data(fi_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xzz_xxxzzz = cbuffer.data(fi_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyy = cbuffer.data(fi_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyyz = cbuffer.data(fi_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xzz_xxyyzz = cbuffer.data(fi_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xzz_xxyzzz = cbuffer.data(fi_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xzz_xxzzzz = cbuffer.data(fi_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyy = cbuffer.data(fi_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyyz = cbuffer.data(fi_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xzz_xyyyzz = cbuffer.data(fi_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xzz_xyyzzz = cbuffer.data(fi_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xzz_xyzzzz = cbuffer.data(fi_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xzz_xzzzzz = cbuffer.data(fi_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyy = cbuffer.data(fi_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyyz = cbuffer.data(fi_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xzz_yyyyzz = cbuffer.data(fi_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xzz_yyyzzz = cbuffer.data(fi_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xzz_yyzzzz = cbuffer.data(fi_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xzz_yzzzzz = cbuffer.data(fi_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xzz_zzzzzz = cbuffer.data(fi_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxx = cbuffer.data(fi_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxy = cbuffer.data(fi_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxxz = cbuffer.data(fi_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyy = cbuffer.data(fi_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxyz = cbuffer.data(fi_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_yyy_xxxxzz = cbuffer.data(fi_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyy = cbuffer.data(fi_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyyz = cbuffer.data(fi_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_yyy_xxxyzz = cbuffer.data(fi_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_yyy_xxxzzz = cbuffer.data(fi_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyy = cbuffer.data(fi_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyyz = cbuffer.data(fi_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_yyy_xxyyzz = cbuffer.data(fi_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_yyy_xxyzzz = cbuffer.data(fi_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_yyy_xxzzzz = cbuffer.data(fi_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyy = cbuffer.data(fi_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyyz = cbuffer.data(fi_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_yyy_xyyyzz = cbuffer.data(fi_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_yyy_xyyzzz = cbuffer.data(fi_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_yyy_xyzzzz = cbuffer.data(fi_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_yyy_xzzzzz = cbuffer.data(fi_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyy = cbuffer.data(fi_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyyz = cbuffer.data(fi_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_yyy_yyyyzz = cbuffer.data(fi_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_yyy_yyyzzz = cbuffer.data(fi_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_yyy_yyzzzz = cbuffer.data(fi_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_yyy_yzzzzz = cbuffer.data(fi_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_yyy_zzzzzz = cbuffer.data(fi_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxx = cbuffer.data(fi_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxy = cbuffer.data(fi_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxxz = cbuffer.data(fi_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyy = cbuffer.data(fi_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxyz = cbuffer.data(fi_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_yyz_xxxxzz = cbuffer.data(fi_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyy = cbuffer.data(fi_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyyz = cbuffer.data(fi_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_yyz_xxxyzz = cbuffer.data(fi_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_yyz_xxxzzz = cbuffer.data(fi_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyy = cbuffer.data(fi_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyyz = cbuffer.data(fi_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_yyz_xxyyzz = cbuffer.data(fi_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_yyz_xxyzzz = cbuffer.data(fi_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_yyz_xxzzzz = cbuffer.data(fi_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyy = cbuffer.data(fi_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyyz = cbuffer.data(fi_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_yyz_xyyyzz = cbuffer.data(fi_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_yyz_xyyzzz = cbuffer.data(fi_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_yyz_xyzzzz = cbuffer.data(fi_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_yyz_xzzzzz = cbuffer.data(fi_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyy = cbuffer.data(fi_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyyz = cbuffer.data(fi_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_yyz_yyyyzz = cbuffer.data(fi_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_yyz_yyyzzz = cbuffer.data(fi_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_yyz_yyzzzz = cbuffer.data(fi_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_yyz_yzzzzz = cbuffer.data(fi_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_yyz_zzzzzz = cbuffer.data(fi_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxx = cbuffer.data(fi_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxy = cbuffer.data(fi_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxxz = cbuffer.data(fi_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyy = cbuffer.data(fi_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxyz = cbuffer.data(fi_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_yzz_xxxxzz = cbuffer.data(fi_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyy = cbuffer.data(fi_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyyz = cbuffer.data(fi_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_yzz_xxxyzz = cbuffer.data(fi_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_yzz_xxxzzz = cbuffer.data(fi_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyy = cbuffer.data(fi_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyyz = cbuffer.data(fi_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_yzz_xxyyzz = cbuffer.data(fi_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_yzz_xxyzzz = cbuffer.data(fi_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_yzz_xxzzzz = cbuffer.data(fi_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyy = cbuffer.data(fi_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyyz = cbuffer.data(fi_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_yzz_xyyyzz = cbuffer.data(fi_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_yzz_xyyzzz = cbuffer.data(fi_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_yzz_xyzzzz = cbuffer.data(fi_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_yzz_xzzzzz = cbuffer.data(fi_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyy = cbuffer.data(fi_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyyz = cbuffer.data(fi_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_yzz_yyyyzz = cbuffer.data(fi_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_yzz_yyyzzz = cbuffer.data(fi_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_yzz_yyzzzz = cbuffer.data(fi_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_yzz_yzzzzz = cbuffer.data(fi_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_yzz_zzzzzz = cbuffer.data(fi_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxx = cbuffer.data(fi_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxy = cbuffer.data(fi_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxxz = cbuffer.data(fi_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyy = cbuffer.data(fi_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxyz = cbuffer.data(fi_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_zzz_xxxxzz = cbuffer.data(fi_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyy = cbuffer.data(fi_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyyz = cbuffer.data(fi_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_zzz_xxxyzz = cbuffer.data(fi_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_zzz_xxxzzz = cbuffer.data(fi_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyy = cbuffer.data(fi_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyyz = cbuffer.data(fi_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_zzz_xxyyzz = cbuffer.data(fi_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_zzz_xxyzzz = cbuffer.data(fi_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_zzz_xxzzzz = cbuffer.data(fi_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyy = cbuffer.data(fi_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyyz = cbuffer.data(fi_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_zzz_xyyyzz = cbuffer.data(fi_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_zzz_xyyzzz = cbuffer.data(fi_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_zzz_xyzzzz = cbuffer.data(fi_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_zzz_xzzzzz = cbuffer.data(fi_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyy = cbuffer.data(fi_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyyz = cbuffer.data(fi_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_zzz_yyyyzz = cbuffer.data(fi_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_zzz_yyyzzz = cbuffer.data(fi_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_zzz_yyzzzz = cbuffer.data(fi_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_zzz_yzzzzz = cbuffer.data(fi_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_zzz_zzzzzz = cbuffer.data(fi_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxx = cbuffer.data(fi_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxy = cbuffer.data(fi_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxxz = cbuffer.data(fi_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyy = cbuffer.data(fi_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxyz = cbuffer.data(fi_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xxx_xxxxzz = cbuffer.data(fi_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyy = cbuffer.data(fi_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyyz = cbuffer.data(fi_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xxx_xxxyzz = cbuffer.data(fi_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xxx_xxxzzz = cbuffer.data(fi_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyy = cbuffer.data(fi_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyyz = cbuffer.data(fi_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xxx_xxyyzz = cbuffer.data(fi_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xxx_xxyzzz = cbuffer.data(fi_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xxx_xxzzzz = cbuffer.data(fi_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyy = cbuffer.data(fi_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyyz = cbuffer.data(fi_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xxx_xyyyzz = cbuffer.data(fi_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xxx_xyyzzz = cbuffer.data(fi_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xxx_xyzzzz = cbuffer.data(fi_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xxx_xzzzzz = cbuffer.data(fi_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyy = cbuffer.data(fi_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyyz = cbuffer.data(fi_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xxx_yyyyzz = cbuffer.data(fi_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xxx_yyyzzz = cbuffer.data(fi_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_xxx_yyzzzz = cbuffer.data(fi_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xxx_yzzzzz = cbuffer.data(fi_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xxx_zzzzzz = cbuffer.data(fi_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxx = cbuffer.data(fi_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxy = cbuffer.data(fi_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxxz = cbuffer.data(fi_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyy = cbuffer.data(fi_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxyz = cbuffer.data(fi_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xxy_xxxxzz = cbuffer.data(fi_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyy = cbuffer.data(fi_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyyz = cbuffer.data(fi_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xxy_xxxyzz = cbuffer.data(fi_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xxy_xxxzzz = cbuffer.data(fi_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyy = cbuffer.data(fi_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyyz = cbuffer.data(fi_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_xxy_xxyyzz = cbuffer.data(fi_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_xxy_xxyzzz = cbuffer.data(fi_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_xxy_xxzzzz = cbuffer.data(fi_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyy = cbuffer.data(fi_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyyz = cbuffer.data(fi_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_xxy_xyyyzz = cbuffer.data(fi_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_xxy_xyyzzz = cbuffer.data(fi_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_xxy_xyzzzz = cbuffer.data(fi_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_xxy_xzzzzz = cbuffer.data(fi_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyy = cbuffer.data(fi_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyyz = cbuffer.data(fi_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_xxy_yyyyzz = cbuffer.data(fi_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_xxy_yyyzzz = cbuffer.data(fi_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_xxy_yyzzzz = cbuffer.data(fi_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_xxy_yzzzzz = cbuffer.data(fi_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_xxy_zzzzzz = cbuffer.data(fi_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxx = cbuffer.data(fi_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxy = cbuffer.data(fi_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxxz = cbuffer.data(fi_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyy = cbuffer.data(fi_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxyz = cbuffer.data(fi_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_xxz_xxxxzz = cbuffer.data(fi_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyy = cbuffer.data(fi_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyyz = cbuffer.data(fi_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_xxz_xxxyzz = cbuffer.data(fi_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_xxz_xxxzzz = cbuffer.data(fi_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyy = cbuffer.data(fi_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyyz = cbuffer.data(fi_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_xxz_xxyyzz = cbuffer.data(fi_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_xxz_xxyzzz = cbuffer.data(fi_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xxz_xxzzzz = cbuffer.data(fi_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyy = cbuffer.data(fi_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyyz = cbuffer.data(fi_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xxz_xyyyzz = cbuffer.data(fi_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xxz_xyyzzz = cbuffer.data(fi_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xxz_xyzzzz = cbuffer.data(fi_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xxz_xzzzzz = cbuffer.data(fi_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyy = cbuffer.data(fi_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyyz = cbuffer.data(fi_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xxz_yyyyzz = cbuffer.data(fi_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xxz_yyyzzz = cbuffer.data(fi_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xxz_yyzzzz = cbuffer.data(fi_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xxz_yzzzzz = cbuffer.data(fi_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xxz_zzzzzz = cbuffer.data(fi_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxx = cbuffer.data(fi_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxy = cbuffer.data(fi_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxxz = cbuffer.data(fi_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyy = cbuffer.data(fi_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxyz = cbuffer.data(fi_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xyy_xxxxzz = cbuffer.data(fi_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyy = cbuffer.data(fi_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyyz = cbuffer.data(fi_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xyy_xxxyzz = cbuffer.data(fi_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xyy_xxxzzz = cbuffer.data(fi_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyy = cbuffer.data(fi_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyyz = cbuffer.data(fi_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xyy_xxyyzz = cbuffer.data(fi_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xyy_xxyzzz = cbuffer.data(fi_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xyy_xxzzzz = cbuffer.data(fi_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyy = cbuffer.data(fi_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyyz = cbuffer.data(fi_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xyy_xyyyzz = cbuffer.data(fi_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xyy_xyyzzz = cbuffer.data(fi_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xyy_xyzzzz = cbuffer.data(fi_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xyy_xzzzzz = cbuffer.data(fi_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyy = cbuffer.data(fi_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyyz = cbuffer.data(fi_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xyy_yyyyzz = cbuffer.data(fi_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xyy_yyyzzz = cbuffer.data(fi_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xyy_yyzzzz = cbuffer.data(fi_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xyy_yzzzzz = cbuffer.data(fi_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xyy_zzzzzz = cbuffer.data(fi_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxx = cbuffer.data(fi_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxy = cbuffer.data(fi_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxxz = cbuffer.data(fi_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyy = cbuffer.data(fi_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxyz = cbuffer.data(fi_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_xyz_xxxxzz = cbuffer.data(fi_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyy = cbuffer.data(fi_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyyz = cbuffer.data(fi_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_xyz_xxxyzz = cbuffer.data(fi_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_xyz_xxxzzz = cbuffer.data(fi_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyy = cbuffer.data(fi_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyyz = cbuffer.data(fi_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_xyz_xxyyzz = cbuffer.data(fi_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_xyz_xxyzzz = cbuffer.data(fi_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_xyz_xxzzzz = cbuffer.data(fi_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyy = cbuffer.data(fi_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyyz = cbuffer.data(fi_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_xyz_xyyyzz = cbuffer.data(fi_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_xyz_xyyzzz = cbuffer.data(fi_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_xyz_xyzzzz = cbuffer.data(fi_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_xyz_xzzzzz = cbuffer.data(fi_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyy = cbuffer.data(fi_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyyz = cbuffer.data(fi_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_xyz_yyyyzz = cbuffer.data(fi_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_xyz_yyyzzz = cbuffer.data(fi_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_xyz_yyzzzz = cbuffer.data(fi_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_xyz_yzzzzz = cbuffer.data(fi_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_xyz_zzzzzz = cbuffer.data(fi_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxx = cbuffer.data(fi_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxy = cbuffer.data(fi_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxxz = cbuffer.data(fi_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyy = cbuffer.data(fi_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxyz = cbuffer.data(fi_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_xzz_xxxxzz = cbuffer.data(fi_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyy = cbuffer.data(fi_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyyz = cbuffer.data(fi_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_xzz_xxxyzz = cbuffer.data(fi_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_xzz_xxxzzz = cbuffer.data(fi_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyy = cbuffer.data(fi_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyyz = cbuffer.data(fi_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_xzz_xxyyzz = cbuffer.data(fi_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_xzz_xxyzzz = cbuffer.data(fi_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_xzz_xxzzzz = cbuffer.data(fi_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyy = cbuffer.data(fi_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyyz = cbuffer.data(fi_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_xzz_xyyyzz = cbuffer.data(fi_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_xzz_xyyzzz = cbuffer.data(fi_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_xzz_xyzzzz = cbuffer.data(fi_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xzz_xzzzzz = cbuffer.data(fi_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyy = cbuffer.data(fi_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyyz = cbuffer.data(fi_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xzz_yyyyzz = cbuffer.data(fi_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xzz_yyyzzz = cbuffer.data(fi_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xzz_yyzzzz = cbuffer.data(fi_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xzz_yzzzzz = cbuffer.data(fi_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xzz_zzzzzz = cbuffer.data(fi_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxx = cbuffer.data(fi_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxy = cbuffer.data(fi_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxxz = cbuffer.data(fi_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyy = cbuffer.data(fi_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxyz = cbuffer.data(fi_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_yyy_xxxxzz = cbuffer.data(fi_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyy = cbuffer.data(fi_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyyz = cbuffer.data(fi_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_yyy_xxxyzz = cbuffer.data(fi_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_yyy_xxxzzz = cbuffer.data(fi_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyy = cbuffer.data(fi_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyyz = cbuffer.data(fi_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_yyy_xxyyzz = cbuffer.data(fi_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_yyy_xxyzzz = cbuffer.data(fi_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_yyy_xxzzzz = cbuffer.data(fi_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyy = cbuffer.data(fi_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyyz = cbuffer.data(fi_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_yyy_xyyyzz = cbuffer.data(fi_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_yyy_xyyzzz = cbuffer.data(fi_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_yyy_xyzzzz = cbuffer.data(fi_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_yyy_xzzzzz = cbuffer.data(fi_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyy = cbuffer.data(fi_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyyz = cbuffer.data(fi_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_yyy_yyyyzz = cbuffer.data(fi_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_yyy_yyyzzz = cbuffer.data(fi_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_yyy_yyzzzz = cbuffer.data(fi_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_yyy_yzzzzz = cbuffer.data(fi_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_yyy_zzzzzz = cbuffer.data(fi_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxx = cbuffer.data(fi_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxy = cbuffer.data(fi_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxxz = cbuffer.data(fi_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyy = cbuffer.data(fi_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxyz = cbuffer.data(fi_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_yyz_xxxxzz = cbuffer.data(fi_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyy = cbuffer.data(fi_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyyz = cbuffer.data(fi_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_yyz_xxxyzz = cbuffer.data(fi_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_yyz_xxxzzz = cbuffer.data(fi_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyy = cbuffer.data(fi_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyyz = cbuffer.data(fi_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_yyz_xxyyzz = cbuffer.data(fi_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_yyz_xxyzzz = cbuffer.data(fi_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_yyz_xxzzzz = cbuffer.data(fi_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyy = cbuffer.data(fi_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyyz = cbuffer.data(fi_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_yyz_xyyyzz = cbuffer.data(fi_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_yyz_xyyzzz = cbuffer.data(fi_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_yyz_xyzzzz = cbuffer.data(fi_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_yyz_xzzzzz = cbuffer.data(fi_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyy = cbuffer.data(fi_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyyz = cbuffer.data(fi_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_yyz_yyyyzz = cbuffer.data(fi_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_yyz_yyyzzz = cbuffer.data(fi_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_yyz_yyzzzz = cbuffer.data(fi_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_yyz_yzzzzz = cbuffer.data(fi_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_yyz_zzzzzz = cbuffer.data(fi_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxx = cbuffer.data(fi_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxy = cbuffer.data(fi_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxxz = cbuffer.data(fi_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyy = cbuffer.data(fi_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxyz = cbuffer.data(fi_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_yzz_xxxxzz = cbuffer.data(fi_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyy = cbuffer.data(fi_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyyz = cbuffer.data(fi_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_yzz_xxxyzz = cbuffer.data(fi_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_yzz_xxxzzz = cbuffer.data(fi_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyy = cbuffer.data(fi_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyyz = cbuffer.data(fi_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_yzz_xxyyzz = cbuffer.data(fi_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_yzz_xxyzzz = cbuffer.data(fi_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_yzz_xxzzzz = cbuffer.data(fi_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyy = cbuffer.data(fi_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyyz = cbuffer.data(fi_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_yzz_xyyyzz = cbuffer.data(fi_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_yzz_xyyzzz = cbuffer.data(fi_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_yzz_xyzzzz = cbuffer.data(fi_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_yzz_xzzzzz = cbuffer.data(fi_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyy = cbuffer.data(fi_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyyz = cbuffer.data(fi_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_yzz_yyyyzz = cbuffer.data(fi_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_yzz_yyyzzz = cbuffer.data(fi_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_yzz_yyzzzz = cbuffer.data(fi_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_yzz_yzzzzz = cbuffer.data(fi_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_yzz_zzzzzz = cbuffer.data(fi_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxx = cbuffer.data(fi_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxy = cbuffer.data(fi_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxxz = cbuffer.data(fi_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyy = cbuffer.data(fi_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxyz = cbuffer.data(fi_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_zzz_xxxxzz = cbuffer.data(fi_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyy = cbuffer.data(fi_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyyz = cbuffer.data(fi_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_zzz_xxxyzz = cbuffer.data(fi_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_zzz_xxxzzz = cbuffer.data(fi_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyy = cbuffer.data(fi_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyyz = cbuffer.data(fi_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_zzz_xxyyzz = cbuffer.data(fi_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_zzz_xxyzzz = cbuffer.data(fi_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_zzz_xxzzzz = cbuffer.data(fi_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyy = cbuffer.data(fi_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyyz = cbuffer.data(fi_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_zzz_xyyyzz = cbuffer.data(fi_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_zzz_xyyzzz = cbuffer.data(fi_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_zzz_xyzzzz = cbuffer.data(fi_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_zzz_xzzzzz = cbuffer.data(fi_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyy = cbuffer.data(fi_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyyz = cbuffer.data(fi_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_zzz_yyyyzz = cbuffer.data(fi_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_zzz_yyyzzz = cbuffer.data(fi_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_zzz_yyzzzz = cbuffer.data(fi_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_zzz_yzzzzz = cbuffer.data(fi_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_zzz_zzzzzz = cbuffer.data(fi_geom_01_off + 839 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ghxx

            const auto gh_geom_01_off = idx_geom_01_ghxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxx, g_0_x_xxx_xxxxxx, g_0_x_xxx_xxxxxy, g_0_x_xxx_xxxxxz, g_0_x_xxx_xxxxy, g_0_x_xxx_xxxxyy, g_0_x_xxx_xxxxyz, g_0_x_xxx_xxxxz, g_0_x_xxx_xxxxzz, g_0_x_xxx_xxxyy, g_0_x_xxx_xxxyyy, g_0_x_xxx_xxxyyz, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxyzz, g_0_x_xxx_xxxzz, g_0_x_xxx_xxxzzz, g_0_x_xxx_xxyyy, g_0_x_xxx_xxyyyy, g_0_x_xxx_xxyyyz, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyyzz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxyzzz, g_0_x_xxx_xxzzz, g_0_x_xxx_xxzzzz, g_0_x_xxx_xyyyy, g_0_x_xxx_xyyyyy, g_0_x_xxx_xyyyyz, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyyzz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyyzzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xyzzzz, g_0_x_xxx_xzzzz, g_0_x_xxx_xzzzzz, g_0_x_xxx_yyyyy, g_0_x_xxx_yyyyz, g_0_x_xxx_yyyzz, g_0_x_xxx_yyzzz, g_0_x_xxx_yzzzz, g_0_x_xxx_zzzzz, g_0_x_xxxx_xxxxx, g_0_x_xxxx_xxxxy, g_0_x_xxxx_xxxxz, g_0_x_xxxx_xxxyy, g_0_x_xxxx_xxxyz, g_0_x_xxxx_xxxzz, g_0_x_xxxx_xxyyy, g_0_x_xxxx_xxyyz, g_0_x_xxxx_xxyzz, g_0_x_xxxx_xxzzz, g_0_x_xxxx_xyyyy, g_0_x_xxxx_xyyyz, g_0_x_xxxx_xyyzz, g_0_x_xxxx_xyzzz, g_0_x_xxxx_xzzzz, g_0_x_xxxx_yyyyy, g_0_x_xxxx_yyyyz, g_0_x_xxxx_yyyzz, g_0_x_xxxx_yyzzz, g_0_x_xxxx_yzzzz, g_0_x_xxxx_zzzzz, g_xxx_xxxxx, g_xxx_xxxxy, g_xxx_xxxxz, g_xxx_xxxyy, g_xxx_xxxyz, g_xxx_xxxzz, g_xxx_xxyyy, g_xxx_xxyyz, g_xxx_xxyzz, g_xxx_xxzzz, g_xxx_xyyyy, g_xxx_xyyyz, g_xxx_xyyzz, g_xxx_xyzzz, g_xxx_xzzzz, g_xxx_yyyyy, g_xxx_yyyyz, g_xxx_yyyzz, g_xxx_yyzzz, g_xxx_yzzzz, g_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_xxxxx[k] = g_xxx_xxxxx[k] - g_0_x_xxx_xxxxx[k] * ab_x + g_0_x_xxx_xxxxxx[k];

                g_0_x_xxxx_xxxxy[k] = g_xxx_xxxxy[k] - g_0_x_xxx_xxxxy[k] * ab_x + g_0_x_xxx_xxxxxy[k];

                g_0_x_xxxx_xxxxz[k] = g_xxx_xxxxz[k] - g_0_x_xxx_xxxxz[k] * ab_x + g_0_x_xxx_xxxxxz[k];

                g_0_x_xxxx_xxxyy[k] = g_xxx_xxxyy[k] - g_0_x_xxx_xxxyy[k] * ab_x + g_0_x_xxx_xxxxyy[k];

                g_0_x_xxxx_xxxyz[k] = g_xxx_xxxyz[k] - g_0_x_xxx_xxxyz[k] * ab_x + g_0_x_xxx_xxxxyz[k];

                g_0_x_xxxx_xxxzz[k] = g_xxx_xxxzz[k] - g_0_x_xxx_xxxzz[k] * ab_x + g_0_x_xxx_xxxxzz[k];

                g_0_x_xxxx_xxyyy[k] = g_xxx_xxyyy[k] - g_0_x_xxx_xxyyy[k] * ab_x + g_0_x_xxx_xxxyyy[k];

                g_0_x_xxxx_xxyyz[k] = g_xxx_xxyyz[k] - g_0_x_xxx_xxyyz[k] * ab_x + g_0_x_xxx_xxxyyz[k];

                g_0_x_xxxx_xxyzz[k] = g_xxx_xxyzz[k] - g_0_x_xxx_xxyzz[k] * ab_x + g_0_x_xxx_xxxyzz[k];

                g_0_x_xxxx_xxzzz[k] = g_xxx_xxzzz[k] - g_0_x_xxx_xxzzz[k] * ab_x + g_0_x_xxx_xxxzzz[k];

                g_0_x_xxxx_xyyyy[k] = g_xxx_xyyyy[k] - g_0_x_xxx_xyyyy[k] * ab_x + g_0_x_xxx_xxyyyy[k];

                g_0_x_xxxx_xyyyz[k] = g_xxx_xyyyz[k] - g_0_x_xxx_xyyyz[k] * ab_x + g_0_x_xxx_xxyyyz[k];

                g_0_x_xxxx_xyyzz[k] = g_xxx_xyyzz[k] - g_0_x_xxx_xyyzz[k] * ab_x + g_0_x_xxx_xxyyzz[k];

                g_0_x_xxxx_xyzzz[k] = g_xxx_xyzzz[k] - g_0_x_xxx_xyzzz[k] * ab_x + g_0_x_xxx_xxyzzz[k];

                g_0_x_xxxx_xzzzz[k] = g_xxx_xzzzz[k] - g_0_x_xxx_xzzzz[k] * ab_x + g_0_x_xxx_xxzzzz[k];

                g_0_x_xxxx_yyyyy[k] = g_xxx_yyyyy[k] - g_0_x_xxx_yyyyy[k] * ab_x + g_0_x_xxx_xyyyyy[k];

                g_0_x_xxxx_yyyyz[k] = g_xxx_yyyyz[k] - g_0_x_xxx_yyyyz[k] * ab_x + g_0_x_xxx_xyyyyz[k];

                g_0_x_xxxx_yyyzz[k] = g_xxx_yyyzz[k] - g_0_x_xxx_yyyzz[k] * ab_x + g_0_x_xxx_xyyyzz[k];

                g_0_x_xxxx_yyzzz[k] = g_xxx_yyzzz[k] - g_0_x_xxx_yyzzz[k] * ab_x + g_0_x_xxx_xyyzzz[k];

                g_0_x_xxxx_yzzzz[k] = g_xxx_yzzzz[k] - g_0_x_xxx_yzzzz[k] * ab_x + g_0_x_xxx_xyzzzz[k];

                g_0_x_xxxx_zzzzz[k] = g_xxx_zzzzz[k] - g_0_x_xxx_zzzzz[k] * ab_x + g_0_x_xxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxx, g_0_x_xxx_xxxxxy, g_0_x_xxx_xxxxy, g_0_x_xxx_xxxxyy, g_0_x_xxx_xxxxyz, g_0_x_xxx_xxxxz, g_0_x_xxx_xxxyy, g_0_x_xxx_xxxyyy, g_0_x_xxx_xxxyyz, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxyzz, g_0_x_xxx_xxxzz, g_0_x_xxx_xxyyy, g_0_x_xxx_xxyyyy, g_0_x_xxx_xxyyyz, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyyzz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxyzzz, g_0_x_xxx_xxzzz, g_0_x_xxx_xyyyy, g_0_x_xxx_xyyyyy, g_0_x_xxx_xyyyyz, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyyzz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyyzzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xyzzzz, g_0_x_xxx_xzzzz, g_0_x_xxx_yyyyy, g_0_x_xxx_yyyyyy, g_0_x_xxx_yyyyyz, g_0_x_xxx_yyyyz, g_0_x_xxx_yyyyzz, g_0_x_xxx_yyyzz, g_0_x_xxx_yyyzzz, g_0_x_xxx_yyzzz, g_0_x_xxx_yyzzzz, g_0_x_xxx_yzzzz, g_0_x_xxx_yzzzzz, g_0_x_xxx_zzzzz, g_0_x_xxxy_xxxxx, g_0_x_xxxy_xxxxy, g_0_x_xxxy_xxxxz, g_0_x_xxxy_xxxyy, g_0_x_xxxy_xxxyz, g_0_x_xxxy_xxxzz, g_0_x_xxxy_xxyyy, g_0_x_xxxy_xxyyz, g_0_x_xxxy_xxyzz, g_0_x_xxxy_xxzzz, g_0_x_xxxy_xyyyy, g_0_x_xxxy_xyyyz, g_0_x_xxxy_xyyzz, g_0_x_xxxy_xyzzz, g_0_x_xxxy_xzzzz, g_0_x_xxxy_yyyyy, g_0_x_xxxy_yyyyz, g_0_x_xxxy_yyyzz, g_0_x_xxxy_yyzzz, g_0_x_xxxy_yzzzz, g_0_x_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_xxxxx[k] = -g_0_x_xxx_xxxxx[k] * ab_y + g_0_x_xxx_xxxxxy[k];

                g_0_x_xxxy_xxxxy[k] = -g_0_x_xxx_xxxxy[k] * ab_y + g_0_x_xxx_xxxxyy[k];

                g_0_x_xxxy_xxxxz[k] = -g_0_x_xxx_xxxxz[k] * ab_y + g_0_x_xxx_xxxxyz[k];

                g_0_x_xxxy_xxxyy[k] = -g_0_x_xxx_xxxyy[k] * ab_y + g_0_x_xxx_xxxyyy[k];

                g_0_x_xxxy_xxxyz[k] = -g_0_x_xxx_xxxyz[k] * ab_y + g_0_x_xxx_xxxyyz[k];

                g_0_x_xxxy_xxxzz[k] = -g_0_x_xxx_xxxzz[k] * ab_y + g_0_x_xxx_xxxyzz[k];

                g_0_x_xxxy_xxyyy[k] = -g_0_x_xxx_xxyyy[k] * ab_y + g_0_x_xxx_xxyyyy[k];

                g_0_x_xxxy_xxyyz[k] = -g_0_x_xxx_xxyyz[k] * ab_y + g_0_x_xxx_xxyyyz[k];

                g_0_x_xxxy_xxyzz[k] = -g_0_x_xxx_xxyzz[k] * ab_y + g_0_x_xxx_xxyyzz[k];

                g_0_x_xxxy_xxzzz[k] = -g_0_x_xxx_xxzzz[k] * ab_y + g_0_x_xxx_xxyzzz[k];

                g_0_x_xxxy_xyyyy[k] = -g_0_x_xxx_xyyyy[k] * ab_y + g_0_x_xxx_xyyyyy[k];

                g_0_x_xxxy_xyyyz[k] = -g_0_x_xxx_xyyyz[k] * ab_y + g_0_x_xxx_xyyyyz[k];

                g_0_x_xxxy_xyyzz[k] = -g_0_x_xxx_xyyzz[k] * ab_y + g_0_x_xxx_xyyyzz[k];

                g_0_x_xxxy_xyzzz[k] = -g_0_x_xxx_xyzzz[k] * ab_y + g_0_x_xxx_xyyzzz[k];

                g_0_x_xxxy_xzzzz[k] = -g_0_x_xxx_xzzzz[k] * ab_y + g_0_x_xxx_xyzzzz[k];

                g_0_x_xxxy_yyyyy[k] = -g_0_x_xxx_yyyyy[k] * ab_y + g_0_x_xxx_yyyyyy[k];

                g_0_x_xxxy_yyyyz[k] = -g_0_x_xxx_yyyyz[k] * ab_y + g_0_x_xxx_yyyyyz[k];

                g_0_x_xxxy_yyyzz[k] = -g_0_x_xxx_yyyzz[k] * ab_y + g_0_x_xxx_yyyyzz[k];

                g_0_x_xxxy_yyzzz[k] = -g_0_x_xxx_yyzzz[k] * ab_y + g_0_x_xxx_yyyzzz[k];

                g_0_x_xxxy_yzzzz[k] = -g_0_x_xxx_yzzzz[k] * ab_y + g_0_x_xxx_yyzzzz[k];

                g_0_x_xxxy_zzzzz[k] = -g_0_x_xxx_zzzzz[k] * ab_y + g_0_x_xxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxx_xxxxx, g_0_x_xxx_xxxxxz, g_0_x_xxx_xxxxy, g_0_x_xxx_xxxxyz, g_0_x_xxx_xxxxz, g_0_x_xxx_xxxxzz, g_0_x_xxx_xxxyy, g_0_x_xxx_xxxyyz, g_0_x_xxx_xxxyz, g_0_x_xxx_xxxyzz, g_0_x_xxx_xxxzz, g_0_x_xxx_xxxzzz, g_0_x_xxx_xxyyy, g_0_x_xxx_xxyyyz, g_0_x_xxx_xxyyz, g_0_x_xxx_xxyyzz, g_0_x_xxx_xxyzz, g_0_x_xxx_xxyzzz, g_0_x_xxx_xxzzz, g_0_x_xxx_xxzzzz, g_0_x_xxx_xyyyy, g_0_x_xxx_xyyyyz, g_0_x_xxx_xyyyz, g_0_x_xxx_xyyyzz, g_0_x_xxx_xyyzz, g_0_x_xxx_xyyzzz, g_0_x_xxx_xyzzz, g_0_x_xxx_xyzzzz, g_0_x_xxx_xzzzz, g_0_x_xxx_xzzzzz, g_0_x_xxx_yyyyy, g_0_x_xxx_yyyyyz, g_0_x_xxx_yyyyz, g_0_x_xxx_yyyyzz, g_0_x_xxx_yyyzz, g_0_x_xxx_yyyzzz, g_0_x_xxx_yyzzz, g_0_x_xxx_yyzzzz, g_0_x_xxx_yzzzz, g_0_x_xxx_yzzzzz, g_0_x_xxx_zzzzz, g_0_x_xxx_zzzzzz, g_0_x_xxxz_xxxxx, g_0_x_xxxz_xxxxy, g_0_x_xxxz_xxxxz, g_0_x_xxxz_xxxyy, g_0_x_xxxz_xxxyz, g_0_x_xxxz_xxxzz, g_0_x_xxxz_xxyyy, g_0_x_xxxz_xxyyz, g_0_x_xxxz_xxyzz, g_0_x_xxxz_xxzzz, g_0_x_xxxz_xyyyy, g_0_x_xxxz_xyyyz, g_0_x_xxxz_xyyzz, g_0_x_xxxz_xyzzz, g_0_x_xxxz_xzzzz, g_0_x_xxxz_yyyyy, g_0_x_xxxz_yyyyz, g_0_x_xxxz_yyyzz, g_0_x_xxxz_yyzzz, g_0_x_xxxz_yzzzz, g_0_x_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_xxxxx[k] = -g_0_x_xxx_xxxxx[k] * ab_z + g_0_x_xxx_xxxxxz[k];

                g_0_x_xxxz_xxxxy[k] = -g_0_x_xxx_xxxxy[k] * ab_z + g_0_x_xxx_xxxxyz[k];

                g_0_x_xxxz_xxxxz[k] = -g_0_x_xxx_xxxxz[k] * ab_z + g_0_x_xxx_xxxxzz[k];

                g_0_x_xxxz_xxxyy[k] = -g_0_x_xxx_xxxyy[k] * ab_z + g_0_x_xxx_xxxyyz[k];

                g_0_x_xxxz_xxxyz[k] = -g_0_x_xxx_xxxyz[k] * ab_z + g_0_x_xxx_xxxyzz[k];

                g_0_x_xxxz_xxxzz[k] = -g_0_x_xxx_xxxzz[k] * ab_z + g_0_x_xxx_xxxzzz[k];

                g_0_x_xxxz_xxyyy[k] = -g_0_x_xxx_xxyyy[k] * ab_z + g_0_x_xxx_xxyyyz[k];

                g_0_x_xxxz_xxyyz[k] = -g_0_x_xxx_xxyyz[k] * ab_z + g_0_x_xxx_xxyyzz[k];

                g_0_x_xxxz_xxyzz[k] = -g_0_x_xxx_xxyzz[k] * ab_z + g_0_x_xxx_xxyzzz[k];

                g_0_x_xxxz_xxzzz[k] = -g_0_x_xxx_xxzzz[k] * ab_z + g_0_x_xxx_xxzzzz[k];

                g_0_x_xxxz_xyyyy[k] = -g_0_x_xxx_xyyyy[k] * ab_z + g_0_x_xxx_xyyyyz[k];

                g_0_x_xxxz_xyyyz[k] = -g_0_x_xxx_xyyyz[k] * ab_z + g_0_x_xxx_xyyyzz[k];

                g_0_x_xxxz_xyyzz[k] = -g_0_x_xxx_xyyzz[k] * ab_z + g_0_x_xxx_xyyzzz[k];

                g_0_x_xxxz_xyzzz[k] = -g_0_x_xxx_xyzzz[k] * ab_z + g_0_x_xxx_xyzzzz[k];

                g_0_x_xxxz_xzzzz[k] = -g_0_x_xxx_xzzzz[k] * ab_z + g_0_x_xxx_xzzzzz[k];

                g_0_x_xxxz_yyyyy[k] = -g_0_x_xxx_yyyyy[k] * ab_z + g_0_x_xxx_yyyyyz[k];

                g_0_x_xxxz_yyyyz[k] = -g_0_x_xxx_yyyyz[k] * ab_z + g_0_x_xxx_yyyyzz[k];

                g_0_x_xxxz_yyyzz[k] = -g_0_x_xxx_yyyzz[k] * ab_z + g_0_x_xxx_yyyzzz[k];

                g_0_x_xxxz_yyzzz[k] = -g_0_x_xxx_yyzzz[k] * ab_z + g_0_x_xxx_yyzzzz[k];

                g_0_x_xxxz_yzzzz[k] = -g_0_x_xxx_yzzzz[k] * ab_z + g_0_x_xxx_yzzzzz[k];

                g_0_x_xxxz_zzzzz[k] = -g_0_x_xxx_zzzzz[k] * ab_z + g_0_x_xxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxy_xxxxx, g_0_x_xxy_xxxxxy, g_0_x_xxy_xxxxy, g_0_x_xxy_xxxxyy, g_0_x_xxy_xxxxyz, g_0_x_xxy_xxxxz, g_0_x_xxy_xxxyy, g_0_x_xxy_xxxyyy, g_0_x_xxy_xxxyyz, g_0_x_xxy_xxxyz, g_0_x_xxy_xxxyzz, g_0_x_xxy_xxxzz, g_0_x_xxy_xxyyy, g_0_x_xxy_xxyyyy, g_0_x_xxy_xxyyyz, g_0_x_xxy_xxyyz, g_0_x_xxy_xxyyzz, g_0_x_xxy_xxyzz, g_0_x_xxy_xxyzzz, g_0_x_xxy_xxzzz, g_0_x_xxy_xyyyy, g_0_x_xxy_xyyyyy, g_0_x_xxy_xyyyyz, g_0_x_xxy_xyyyz, g_0_x_xxy_xyyyzz, g_0_x_xxy_xyyzz, g_0_x_xxy_xyyzzz, g_0_x_xxy_xyzzz, g_0_x_xxy_xyzzzz, g_0_x_xxy_xzzzz, g_0_x_xxy_yyyyy, g_0_x_xxy_yyyyyy, g_0_x_xxy_yyyyyz, g_0_x_xxy_yyyyz, g_0_x_xxy_yyyyzz, g_0_x_xxy_yyyzz, g_0_x_xxy_yyyzzz, g_0_x_xxy_yyzzz, g_0_x_xxy_yyzzzz, g_0_x_xxy_yzzzz, g_0_x_xxy_yzzzzz, g_0_x_xxy_zzzzz, g_0_x_xxyy_xxxxx, g_0_x_xxyy_xxxxy, g_0_x_xxyy_xxxxz, g_0_x_xxyy_xxxyy, g_0_x_xxyy_xxxyz, g_0_x_xxyy_xxxzz, g_0_x_xxyy_xxyyy, g_0_x_xxyy_xxyyz, g_0_x_xxyy_xxyzz, g_0_x_xxyy_xxzzz, g_0_x_xxyy_xyyyy, g_0_x_xxyy_xyyyz, g_0_x_xxyy_xyyzz, g_0_x_xxyy_xyzzz, g_0_x_xxyy_xzzzz, g_0_x_xxyy_yyyyy, g_0_x_xxyy_yyyyz, g_0_x_xxyy_yyyzz, g_0_x_xxyy_yyzzz, g_0_x_xxyy_yzzzz, g_0_x_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_xxxxx[k] = -g_0_x_xxy_xxxxx[k] * ab_y + g_0_x_xxy_xxxxxy[k];

                g_0_x_xxyy_xxxxy[k] = -g_0_x_xxy_xxxxy[k] * ab_y + g_0_x_xxy_xxxxyy[k];

                g_0_x_xxyy_xxxxz[k] = -g_0_x_xxy_xxxxz[k] * ab_y + g_0_x_xxy_xxxxyz[k];

                g_0_x_xxyy_xxxyy[k] = -g_0_x_xxy_xxxyy[k] * ab_y + g_0_x_xxy_xxxyyy[k];

                g_0_x_xxyy_xxxyz[k] = -g_0_x_xxy_xxxyz[k] * ab_y + g_0_x_xxy_xxxyyz[k];

                g_0_x_xxyy_xxxzz[k] = -g_0_x_xxy_xxxzz[k] * ab_y + g_0_x_xxy_xxxyzz[k];

                g_0_x_xxyy_xxyyy[k] = -g_0_x_xxy_xxyyy[k] * ab_y + g_0_x_xxy_xxyyyy[k];

                g_0_x_xxyy_xxyyz[k] = -g_0_x_xxy_xxyyz[k] * ab_y + g_0_x_xxy_xxyyyz[k];

                g_0_x_xxyy_xxyzz[k] = -g_0_x_xxy_xxyzz[k] * ab_y + g_0_x_xxy_xxyyzz[k];

                g_0_x_xxyy_xxzzz[k] = -g_0_x_xxy_xxzzz[k] * ab_y + g_0_x_xxy_xxyzzz[k];

                g_0_x_xxyy_xyyyy[k] = -g_0_x_xxy_xyyyy[k] * ab_y + g_0_x_xxy_xyyyyy[k];

                g_0_x_xxyy_xyyyz[k] = -g_0_x_xxy_xyyyz[k] * ab_y + g_0_x_xxy_xyyyyz[k];

                g_0_x_xxyy_xyyzz[k] = -g_0_x_xxy_xyyzz[k] * ab_y + g_0_x_xxy_xyyyzz[k];

                g_0_x_xxyy_xyzzz[k] = -g_0_x_xxy_xyzzz[k] * ab_y + g_0_x_xxy_xyyzzz[k];

                g_0_x_xxyy_xzzzz[k] = -g_0_x_xxy_xzzzz[k] * ab_y + g_0_x_xxy_xyzzzz[k];

                g_0_x_xxyy_yyyyy[k] = -g_0_x_xxy_yyyyy[k] * ab_y + g_0_x_xxy_yyyyyy[k];

                g_0_x_xxyy_yyyyz[k] = -g_0_x_xxy_yyyyz[k] * ab_y + g_0_x_xxy_yyyyyz[k];

                g_0_x_xxyy_yyyzz[k] = -g_0_x_xxy_yyyzz[k] * ab_y + g_0_x_xxy_yyyyzz[k];

                g_0_x_xxyy_yyzzz[k] = -g_0_x_xxy_yyzzz[k] * ab_y + g_0_x_xxy_yyyzzz[k];

                g_0_x_xxyy_yzzzz[k] = -g_0_x_xxy_yzzzz[k] * ab_y + g_0_x_xxy_yyzzzz[k];

                g_0_x_xxyy_zzzzz[k] = -g_0_x_xxy_zzzzz[k] * ab_y + g_0_x_xxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxyz_xxxxx, g_0_x_xxyz_xxxxy, g_0_x_xxyz_xxxxz, g_0_x_xxyz_xxxyy, g_0_x_xxyz_xxxyz, g_0_x_xxyz_xxxzz, g_0_x_xxyz_xxyyy, g_0_x_xxyz_xxyyz, g_0_x_xxyz_xxyzz, g_0_x_xxyz_xxzzz, g_0_x_xxyz_xyyyy, g_0_x_xxyz_xyyyz, g_0_x_xxyz_xyyzz, g_0_x_xxyz_xyzzz, g_0_x_xxyz_xzzzz, g_0_x_xxyz_yyyyy, g_0_x_xxyz_yyyyz, g_0_x_xxyz_yyyzz, g_0_x_xxyz_yyzzz, g_0_x_xxyz_yzzzz, g_0_x_xxyz_zzzzz, g_0_x_xxz_xxxxx, g_0_x_xxz_xxxxxy, g_0_x_xxz_xxxxy, g_0_x_xxz_xxxxyy, g_0_x_xxz_xxxxyz, g_0_x_xxz_xxxxz, g_0_x_xxz_xxxyy, g_0_x_xxz_xxxyyy, g_0_x_xxz_xxxyyz, g_0_x_xxz_xxxyz, g_0_x_xxz_xxxyzz, g_0_x_xxz_xxxzz, g_0_x_xxz_xxyyy, g_0_x_xxz_xxyyyy, g_0_x_xxz_xxyyyz, g_0_x_xxz_xxyyz, g_0_x_xxz_xxyyzz, g_0_x_xxz_xxyzz, g_0_x_xxz_xxyzzz, g_0_x_xxz_xxzzz, g_0_x_xxz_xyyyy, g_0_x_xxz_xyyyyy, g_0_x_xxz_xyyyyz, g_0_x_xxz_xyyyz, g_0_x_xxz_xyyyzz, g_0_x_xxz_xyyzz, g_0_x_xxz_xyyzzz, g_0_x_xxz_xyzzz, g_0_x_xxz_xyzzzz, g_0_x_xxz_xzzzz, g_0_x_xxz_yyyyy, g_0_x_xxz_yyyyyy, g_0_x_xxz_yyyyyz, g_0_x_xxz_yyyyz, g_0_x_xxz_yyyyzz, g_0_x_xxz_yyyzz, g_0_x_xxz_yyyzzz, g_0_x_xxz_yyzzz, g_0_x_xxz_yyzzzz, g_0_x_xxz_yzzzz, g_0_x_xxz_yzzzzz, g_0_x_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_xxxxx[k] = -g_0_x_xxz_xxxxx[k] * ab_y + g_0_x_xxz_xxxxxy[k];

                g_0_x_xxyz_xxxxy[k] = -g_0_x_xxz_xxxxy[k] * ab_y + g_0_x_xxz_xxxxyy[k];

                g_0_x_xxyz_xxxxz[k] = -g_0_x_xxz_xxxxz[k] * ab_y + g_0_x_xxz_xxxxyz[k];

                g_0_x_xxyz_xxxyy[k] = -g_0_x_xxz_xxxyy[k] * ab_y + g_0_x_xxz_xxxyyy[k];

                g_0_x_xxyz_xxxyz[k] = -g_0_x_xxz_xxxyz[k] * ab_y + g_0_x_xxz_xxxyyz[k];

                g_0_x_xxyz_xxxzz[k] = -g_0_x_xxz_xxxzz[k] * ab_y + g_0_x_xxz_xxxyzz[k];

                g_0_x_xxyz_xxyyy[k] = -g_0_x_xxz_xxyyy[k] * ab_y + g_0_x_xxz_xxyyyy[k];

                g_0_x_xxyz_xxyyz[k] = -g_0_x_xxz_xxyyz[k] * ab_y + g_0_x_xxz_xxyyyz[k];

                g_0_x_xxyz_xxyzz[k] = -g_0_x_xxz_xxyzz[k] * ab_y + g_0_x_xxz_xxyyzz[k];

                g_0_x_xxyz_xxzzz[k] = -g_0_x_xxz_xxzzz[k] * ab_y + g_0_x_xxz_xxyzzz[k];

                g_0_x_xxyz_xyyyy[k] = -g_0_x_xxz_xyyyy[k] * ab_y + g_0_x_xxz_xyyyyy[k];

                g_0_x_xxyz_xyyyz[k] = -g_0_x_xxz_xyyyz[k] * ab_y + g_0_x_xxz_xyyyyz[k];

                g_0_x_xxyz_xyyzz[k] = -g_0_x_xxz_xyyzz[k] * ab_y + g_0_x_xxz_xyyyzz[k];

                g_0_x_xxyz_xyzzz[k] = -g_0_x_xxz_xyzzz[k] * ab_y + g_0_x_xxz_xyyzzz[k];

                g_0_x_xxyz_xzzzz[k] = -g_0_x_xxz_xzzzz[k] * ab_y + g_0_x_xxz_xyzzzz[k];

                g_0_x_xxyz_yyyyy[k] = -g_0_x_xxz_yyyyy[k] * ab_y + g_0_x_xxz_yyyyyy[k];

                g_0_x_xxyz_yyyyz[k] = -g_0_x_xxz_yyyyz[k] * ab_y + g_0_x_xxz_yyyyyz[k];

                g_0_x_xxyz_yyyzz[k] = -g_0_x_xxz_yyyzz[k] * ab_y + g_0_x_xxz_yyyyzz[k];

                g_0_x_xxyz_yyzzz[k] = -g_0_x_xxz_yyzzz[k] * ab_y + g_0_x_xxz_yyyzzz[k];

                g_0_x_xxyz_yzzzz[k] = -g_0_x_xxz_yzzzz[k] * ab_y + g_0_x_xxz_yyzzzz[k];

                g_0_x_xxyz_zzzzz[k] = -g_0_x_xxz_zzzzz[k] * ab_y + g_0_x_xxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xxz_xxxxx, g_0_x_xxz_xxxxxz, g_0_x_xxz_xxxxy, g_0_x_xxz_xxxxyz, g_0_x_xxz_xxxxz, g_0_x_xxz_xxxxzz, g_0_x_xxz_xxxyy, g_0_x_xxz_xxxyyz, g_0_x_xxz_xxxyz, g_0_x_xxz_xxxyzz, g_0_x_xxz_xxxzz, g_0_x_xxz_xxxzzz, g_0_x_xxz_xxyyy, g_0_x_xxz_xxyyyz, g_0_x_xxz_xxyyz, g_0_x_xxz_xxyyzz, g_0_x_xxz_xxyzz, g_0_x_xxz_xxyzzz, g_0_x_xxz_xxzzz, g_0_x_xxz_xxzzzz, g_0_x_xxz_xyyyy, g_0_x_xxz_xyyyyz, g_0_x_xxz_xyyyz, g_0_x_xxz_xyyyzz, g_0_x_xxz_xyyzz, g_0_x_xxz_xyyzzz, g_0_x_xxz_xyzzz, g_0_x_xxz_xyzzzz, g_0_x_xxz_xzzzz, g_0_x_xxz_xzzzzz, g_0_x_xxz_yyyyy, g_0_x_xxz_yyyyyz, g_0_x_xxz_yyyyz, g_0_x_xxz_yyyyzz, g_0_x_xxz_yyyzz, g_0_x_xxz_yyyzzz, g_0_x_xxz_yyzzz, g_0_x_xxz_yyzzzz, g_0_x_xxz_yzzzz, g_0_x_xxz_yzzzzz, g_0_x_xxz_zzzzz, g_0_x_xxz_zzzzzz, g_0_x_xxzz_xxxxx, g_0_x_xxzz_xxxxy, g_0_x_xxzz_xxxxz, g_0_x_xxzz_xxxyy, g_0_x_xxzz_xxxyz, g_0_x_xxzz_xxxzz, g_0_x_xxzz_xxyyy, g_0_x_xxzz_xxyyz, g_0_x_xxzz_xxyzz, g_0_x_xxzz_xxzzz, g_0_x_xxzz_xyyyy, g_0_x_xxzz_xyyyz, g_0_x_xxzz_xyyzz, g_0_x_xxzz_xyzzz, g_0_x_xxzz_xzzzz, g_0_x_xxzz_yyyyy, g_0_x_xxzz_yyyyz, g_0_x_xxzz_yyyzz, g_0_x_xxzz_yyzzz, g_0_x_xxzz_yzzzz, g_0_x_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_xxxxx[k] = -g_0_x_xxz_xxxxx[k] * ab_z + g_0_x_xxz_xxxxxz[k];

                g_0_x_xxzz_xxxxy[k] = -g_0_x_xxz_xxxxy[k] * ab_z + g_0_x_xxz_xxxxyz[k];

                g_0_x_xxzz_xxxxz[k] = -g_0_x_xxz_xxxxz[k] * ab_z + g_0_x_xxz_xxxxzz[k];

                g_0_x_xxzz_xxxyy[k] = -g_0_x_xxz_xxxyy[k] * ab_z + g_0_x_xxz_xxxyyz[k];

                g_0_x_xxzz_xxxyz[k] = -g_0_x_xxz_xxxyz[k] * ab_z + g_0_x_xxz_xxxyzz[k];

                g_0_x_xxzz_xxxzz[k] = -g_0_x_xxz_xxxzz[k] * ab_z + g_0_x_xxz_xxxzzz[k];

                g_0_x_xxzz_xxyyy[k] = -g_0_x_xxz_xxyyy[k] * ab_z + g_0_x_xxz_xxyyyz[k];

                g_0_x_xxzz_xxyyz[k] = -g_0_x_xxz_xxyyz[k] * ab_z + g_0_x_xxz_xxyyzz[k];

                g_0_x_xxzz_xxyzz[k] = -g_0_x_xxz_xxyzz[k] * ab_z + g_0_x_xxz_xxyzzz[k];

                g_0_x_xxzz_xxzzz[k] = -g_0_x_xxz_xxzzz[k] * ab_z + g_0_x_xxz_xxzzzz[k];

                g_0_x_xxzz_xyyyy[k] = -g_0_x_xxz_xyyyy[k] * ab_z + g_0_x_xxz_xyyyyz[k];

                g_0_x_xxzz_xyyyz[k] = -g_0_x_xxz_xyyyz[k] * ab_z + g_0_x_xxz_xyyyzz[k];

                g_0_x_xxzz_xyyzz[k] = -g_0_x_xxz_xyyzz[k] * ab_z + g_0_x_xxz_xyyzzz[k];

                g_0_x_xxzz_xyzzz[k] = -g_0_x_xxz_xyzzz[k] * ab_z + g_0_x_xxz_xyzzzz[k];

                g_0_x_xxzz_xzzzz[k] = -g_0_x_xxz_xzzzz[k] * ab_z + g_0_x_xxz_xzzzzz[k];

                g_0_x_xxzz_yyyyy[k] = -g_0_x_xxz_yyyyy[k] * ab_z + g_0_x_xxz_yyyyyz[k];

                g_0_x_xxzz_yyyyz[k] = -g_0_x_xxz_yyyyz[k] * ab_z + g_0_x_xxz_yyyyzz[k];

                g_0_x_xxzz_yyyzz[k] = -g_0_x_xxz_yyyzz[k] * ab_z + g_0_x_xxz_yyyzzz[k];

                g_0_x_xxzz_yyzzz[k] = -g_0_x_xxz_yyzzz[k] * ab_z + g_0_x_xxz_yyzzzz[k];

                g_0_x_xxzz_yzzzz[k] = -g_0_x_xxz_yzzzz[k] * ab_z + g_0_x_xxz_yzzzzz[k];

                g_0_x_xxzz_zzzzz[k] = -g_0_x_xxz_zzzzz[k] * ab_z + g_0_x_xxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyy_xxxxx, g_0_x_xyy_xxxxxy, g_0_x_xyy_xxxxy, g_0_x_xyy_xxxxyy, g_0_x_xyy_xxxxyz, g_0_x_xyy_xxxxz, g_0_x_xyy_xxxyy, g_0_x_xyy_xxxyyy, g_0_x_xyy_xxxyyz, g_0_x_xyy_xxxyz, g_0_x_xyy_xxxyzz, g_0_x_xyy_xxxzz, g_0_x_xyy_xxyyy, g_0_x_xyy_xxyyyy, g_0_x_xyy_xxyyyz, g_0_x_xyy_xxyyz, g_0_x_xyy_xxyyzz, g_0_x_xyy_xxyzz, g_0_x_xyy_xxyzzz, g_0_x_xyy_xxzzz, g_0_x_xyy_xyyyy, g_0_x_xyy_xyyyyy, g_0_x_xyy_xyyyyz, g_0_x_xyy_xyyyz, g_0_x_xyy_xyyyzz, g_0_x_xyy_xyyzz, g_0_x_xyy_xyyzzz, g_0_x_xyy_xyzzz, g_0_x_xyy_xyzzzz, g_0_x_xyy_xzzzz, g_0_x_xyy_yyyyy, g_0_x_xyy_yyyyyy, g_0_x_xyy_yyyyyz, g_0_x_xyy_yyyyz, g_0_x_xyy_yyyyzz, g_0_x_xyy_yyyzz, g_0_x_xyy_yyyzzz, g_0_x_xyy_yyzzz, g_0_x_xyy_yyzzzz, g_0_x_xyy_yzzzz, g_0_x_xyy_yzzzzz, g_0_x_xyy_zzzzz, g_0_x_xyyy_xxxxx, g_0_x_xyyy_xxxxy, g_0_x_xyyy_xxxxz, g_0_x_xyyy_xxxyy, g_0_x_xyyy_xxxyz, g_0_x_xyyy_xxxzz, g_0_x_xyyy_xxyyy, g_0_x_xyyy_xxyyz, g_0_x_xyyy_xxyzz, g_0_x_xyyy_xxzzz, g_0_x_xyyy_xyyyy, g_0_x_xyyy_xyyyz, g_0_x_xyyy_xyyzz, g_0_x_xyyy_xyzzz, g_0_x_xyyy_xzzzz, g_0_x_xyyy_yyyyy, g_0_x_xyyy_yyyyz, g_0_x_xyyy_yyyzz, g_0_x_xyyy_yyzzz, g_0_x_xyyy_yzzzz, g_0_x_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_xxxxx[k] = -g_0_x_xyy_xxxxx[k] * ab_y + g_0_x_xyy_xxxxxy[k];

                g_0_x_xyyy_xxxxy[k] = -g_0_x_xyy_xxxxy[k] * ab_y + g_0_x_xyy_xxxxyy[k];

                g_0_x_xyyy_xxxxz[k] = -g_0_x_xyy_xxxxz[k] * ab_y + g_0_x_xyy_xxxxyz[k];

                g_0_x_xyyy_xxxyy[k] = -g_0_x_xyy_xxxyy[k] * ab_y + g_0_x_xyy_xxxyyy[k];

                g_0_x_xyyy_xxxyz[k] = -g_0_x_xyy_xxxyz[k] * ab_y + g_0_x_xyy_xxxyyz[k];

                g_0_x_xyyy_xxxzz[k] = -g_0_x_xyy_xxxzz[k] * ab_y + g_0_x_xyy_xxxyzz[k];

                g_0_x_xyyy_xxyyy[k] = -g_0_x_xyy_xxyyy[k] * ab_y + g_0_x_xyy_xxyyyy[k];

                g_0_x_xyyy_xxyyz[k] = -g_0_x_xyy_xxyyz[k] * ab_y + g_0_x_xyy_xxyyyz[k];

                g_0_x_xyyy_xxyzz[k] = -g_0_x_xyy_xxyzz[k] * ab_y + g_0_x_xyy_xxyyzz[k];

                g_0_x_xyyy_xxzzz[k] = -g_0_x_xyy_xxzzz[k] * ab_y + g_0_x_xyy_xxyzzz[k];

                g_0_x_xyyy_xyyyy[k] = -g_0_x_xyy_xyyyy[k] * ab_y + g_0_x_xyy_xyyyyy[k];

                g_0_x_xyyy_xyyyz[k] = -g_0_x_xyy_xyyyz[k] * ab_y + g_0_x_xyy_xyyyyz[k];

                g_0_x_xyyy_xyyzz[k] = -g_0_x_xyy_xyyzz[k] * ab_y + g_0_x_xyy_xyyyzz[k];

                g_0_x_xyyy_xyzzz[k] = -g_0_x_xyy_xyzzz[k] * ab_y + g_0_x_xyy_xyyzzz[k];

                g_0_x_xyyy_xzzzz[k] = -g_0_x_xyy_xzzzz[k] * ab_y + g_0_x_xyy_xyzzzz[k];

                g_0_x_xyyy_yyyyy[k] = -g_0_x_xyy_yyyyy[k] * ab_y + g_0_x_xyy_yyyyyy[k];

                g_0_x_xyyy_yyyyz[k] = -g_0_x_xyy_yyyyz[k] * ab_y + g_0_x_xyy_yyyyyz[k];

                g_0_x_xyyy_yyyzz[k] = -g_0_x_xyy_yyyzz[k] * ab_y + g_0_x_xyy_yyyyzz[k];

                g_0_x_xyyy_yyzzz[k] = -g_0_x_xyy_yyzzz[k] * ab_y + g_0_x_xyy_yyyzzz[k];

                g_0_x_xyyy_yzzzz[k] = -g_0_x_xyy_yzzzz[k] * ab_y + g_0_x_xyy_yyzzzz[k];

                g_0_x_xyyy_zzzzz[k] = -g_0_x_xyy_zzzzz[k] * ab_y + g_0_x_xyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyyz_xxxxx, g_0_x_xyyz_xxxxy, g_0_x_xyyz_xxxxz, g_0_x_xyyz_xxxyy, g_0_x_xyyz_xxxyz, g_0_x_xyyz_xxxzz, g_0_x_xyyz_xxyyy, g_0_x_xyyz_xxyyz, g_0_x_xyyz_xxyzz, g_0_x_xyyz_xxzzz, g_0_x_xyyz_xyyyy, g_0_x_xyyz_xyyyz, g_0_x_xyyz_xyyzz, g_0_x_xyyz_xyzzz, g_0_x_xyyz_xzzzz, g_0_x_xyyz_yyyyy, g_0_x_xyyz_yyyyz, g_0_x_xyyz_yyyzz, g_0_x_xyyz_yyzzz, g_0_x_xyyz_yzzzz, g_0_x_xyyz_zzzzz, g_0_x_xyz_xxxxx, g_0_x_xyz_xxxxxy, g_0_x_xyz_xxxxy, g_0_x_xyz_xxxxyy, g_0_x_xyz_xxxxyz, g_0_x_xyz_xxxxz, g_0_x_xyz_xxxyy, g_0_x_xyz_xxxyyy, g_0_x_xyz_xxxyyz, g_0_x_xyz_xxxyz, g_0_x_xyz_xxxyzz, g_0_x_xyz_xxxzz, g_0_x_xyz_xxyyy, g_0_x_xyz_xxyyyy, g_0_x_xyz_xxyyyz, g_0_x_xyz_xxyyz, g_0_x_xyz_xxyyzz, g_0_x_xyz_xxyzz, g_0_x_xyz_xxyzzz, g_0_x_xyz_xxzzz, g_0_x_xyz_xyyyy, g_0_x_xyz_xyyyyy, g_0_x_xyz_xyyyyz, g_0_x_xyz_xyyyz, g_0_x_xyz_xyyyzz, g_0_x_xyz_xyyzz, g_0_x_xyz_xyyzzz, g_0_x_xyz_xyzzz, g_0_x_xyz_xyzzzz, g_0_x_xyz_xzzzz, g_0_x_xyz_yyyyy, g_0_x_xyz_yyyyyy, g_0_x_xyz_yyyyyz, g_0_x_xyz_yyyyz, g_0_x_xyz_yyyyzz, g_0_x_xyz_yyyzz, g_0_x_xyz_yyyzzz, g_0_x_xyz_yyzzz, g_0_x_xyz_yyzzzz, g_0_x_xyz_yzzzz, g_0_x_xyz_yzzzzz, g_0_x_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_xxxxx[k] = -g_0_x_xyz_xxxxx[k] * ab_y + g_0_x_xyz_xxxxxy[k];

                g_0_x_xyyz_xxxxy[k] = -g_0_x_xyz_xxxxy[k] * ab_y + g_0_x_xyz_xxxxyy[k];

                g_0_x_xyyz_xxxxz[k] = -g_0_x_xyz_xxxxz[k] * ab_y + g_0_x_xyz_xxxxyz[k];

                g_0_x_xyyz_xxxyy[k] = -g_0_x_xyz_xxxyy[k] * ab_y + g_0_x_xyz_xxxyyy[k];

                g_0_x_xyyz_xxxyz[k] = -g_0_x_xyz_xxxyz[k] * ab_y + g_0_x_xyz_xxxyyz[k];

                g_0_x_xyyz_xxxzz[k] = -g_0_x_xyz_xxxzz[k] * ab_y + g_0_x_xyz_xxxyzz[k];

                g_0_x_xyyz_xxyyy[k] = -g_0_x_xyz_xxyyy[k] * ab_y + g_0_x_xyz_xxyyyy[k];

                g_0_x_xyyz_xxyyz[k] = -g_0_x_xyz_xxyyz[k] * ab_y + g_0_x_xyz_xxyyyz[k];

                g_0_x_xyyz_xxyzz[k] = -g_0_x_xyz_xxyzz[k] * ab_y + g_0_x_xyz_xxyyzz[k];

                g_0_x_xyyz_xxzzz[k] = -g_0_x_xyz_xxzzz[k] * ab_y + g_0_x_xyz_xxyzzz[k];

                g_0_x_xyyz_xyyyy[k] = -g_0_x_xyz_xyyyy[k] * ab_y + g_0_x_xyz_xyyyyy[k];

                g_0_x_xyyz_xyyyz[k] = -g_0_x_xyz_xyyyz[k] * ab_y + g_0_x_xyz_xyyyyz[k];

                g_0_x_xyyz_xyyzz[k] = -g_0_x_xyz_xyyzz[k] * ab_y + g_0_x_xyz_xyyyzz[k];

                g_0_x_xyyz_xyzzz[k] = -g_0_x_xyz_xyzzz[k] * ab_y + g_0_x_xyz_xyyzzz[k];

                g_0_x_xyyz_xzzzz[k] = -g_0_x_xyz_xzzzz[k] * ab_y + g_0_x_xyz_xyzzzz[k];

                g_0_x_xyyz_yyyyy[k] = -g_0_x_xyz_yyyyy[k] * ab_y + g_0_x_xyz_yyyyyy[k];

                g_0_x_xyyz_yyyyz[k] = -g_0_x_xyz_yyyyz[k] * ab_y + g_0_x_xyz_yyyyyz[k];

                g_0_x_xyyz_yyyzz[k] = -g_0_x_xyz_yyyzz[k] * ab_y + g_0_x_xyz_yyyyzz[k];

                g_0_x_xyyz_yyzzz[k] = -g_0_x_xyz_yyzzz[k] * ab_y + g_0_x_xyz_yyyzzz[k];

                g_0_x_xyyz_yzzzz[k] = -g_0_x_xyz_yzzzz[k] * ab_y + g_0_x_xyz_yyzzzz[k];

                g_0_x_xyyz_zzzzz[k] = -g_0_x_xyz_zzzzz[k] * ab_y + g_0_x_xyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyzz_xxxxx, g_0_x_xyzz_xxxxy, g_0_x_xyzz_xxxxz, g_0_x_xyzz_xxxyy, g_0_x_xyzz_xxxyz, g_0_x_xyzz_xxxzz, g_0_x_xyzz_xxyyy, g_0_x_xyzz_xxyyz, g_0_x_xyzz_xxyzz, g_0_x_xyzz_xxzzz, g_0_x_xyzz_xyyyy, g_0_x_xyzz_xyyyz, g_0_x_xyzz_xyyzz, g_0_x_xyzz_xyzzz, g_0_x_xyzz_xzzzz, g_0_x_xyzz_yyyyy, g_0_x_xyzz_yyyyz, g_0_x_xyzz_yyyzz, g_0_x_xyzz_yyzzz, g_0_x_xyzz_yzzzz, g_0_x_xyzz_zzzzz, g_0_x_xzz_xxxxx, g_0_x_xzz_xxxxxy, g_0_x_xzz_xxxxy, g_0_x_xzz_xxxxyy, g_0_x_xzz_xxxxyz, g_0_x_xzz_xxxxz, g_0_x_xzz_xxxyy, g_0_x_xzz_xxxyyy, g_0_x_xzz_xxxyyz, g_0_x_xzz_xxxyz, g_0_x_xzz_xxxyzz, g_0_x_xzz_xxxzz, g_0_x_xzz_xxyyy, g_0_x_xzz_xxyyyy, g_0_x_xzz_xxyyyz, g_0_x_xzz_xxyyz, g_0_x_xzz_xxyyzz, g_0_x_xzz_xxyzz, g_0_x_xzz_xxyzzz, g_0_x_xzz_xxzzz, g_0_x_xzz_xyyyy, g_0_x_xzz_xyyyyy, g_0_x_xzz_xyyyyz, g_0_x_xzz_xyyyz, g_0_x_xzz_xyyyzz, g_0_x_xzz_xyyzz, g_0_x_xzz_xyyzzz, g_0_x_xzz_xyzzz, g_0_x_xzz_xyzzzz, g_0_x_xzz_xzzzz, g_0_x_xzz_yyyyy, g_0_x_xzz_yyyyyy, g_0_x_xzz_yyyyyz, g_0_x_xzz_yyyyz, g_0_x_xzz_yyyyzz, g_0_x_xzz_yyyzz, g_0_x_xzz_yyyzzz, g_0_x_xzz_yyzzz, g_0_x_xzz_yyzzzz, g_0_x_xzz_yzzzz, g_0_x_xzz_yzzzzz, g_0_x_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_xxxxx[k] = -g_0_x_xzz_xxxxx[k] * ab_y + g_0_x_xzz_xxxxxy[k];

                g_0_x_xyzz_xxxxy[k] = -g_0_x_xzz_xxxxy[k] * ab_y + g_0_x_xzz_xxxxyy[k];

                g_0_x_xyzz_xxxxz[k] = -g_0_x_xzz_xxxxz[k] * ab_y + g_0_x_xzz_xxxxyz[k];

                g_0_x_xyzz_xxxyy[k] = -g_0_x_xzz_xxxyy[k] * ab_y + g_0_x_xzz_xxxyyy[k];

                g_0_x_xyzz_xxxyz[k] = -g_0_x_xzz_xxxyz[k] * ab_y + g_0_x_xzz_xxxyyz[k];

                g_0_x_xyzz_xxxzz[k] = -g_0_x_xzz_xxxzz[k] * ab_y + g_0_x_xzz_xxxyzz[k];

                g_0_x_xyzz_xxyyy[k] = -g_0_x_xzz_xxyyy[k] * ab_y + g_0_x_xzz_xxyyyy[k];

                g_0_x_xyzz_xxyyz[k] = -g_0_x_xzz_xxyyz[k] * ab_y + g_0_x_xzz_xxyyyz[k];

                g_0_x_xyzz_xxyzz[k] = -g_0_x_xzz_xxyzz[k] * ab_y + g_0_x_xzz_xxyyzz[k];

                g_0_x_xyzz_xxzzz[k] = -g_0_x_xzz_xxzzz[k] * ab_y + g_0_x_xzz_xxyzzz[k];

                g_0_x_xyzz_xyyyy[k] = -g_0_x_xzz_xyyyy[k] * ab_y + g_0_x_xzz_xyyyyy[k];

                g_0_x_xyzz_xyyyz[k] = -g_0_x_xzz_xyyyz[k] * ab_y + g_0_x_xzz_xyyyyz[k];

                g_0_x_xyzz_xyyzz[k] = -g_0_x_xzz_xyyzz[k] * ab_y + g_0_x_xzz_xyyyzz[k];

                g_0_x_xyzz_xyzzz[k] = -g_0_x_xzz_xyzzz[k] * ab_y + g_0_x_xzz_xyyzzz[k];

                g_0_x_xyzz_xzzzz[k] = -g_0_x_xzz_xzzzz[k] * ab_y + g_0_x_xzz_xyzzzz[k];

                g_0_x_xyzz_yyyyy[k] = -g_0_x_xzz_yyyyy[k] * ab_y + g_0_x_xzz_yyyyyy[k];

                g_0_x_xyzz_yyyyz[k] = -g_0_x_xzz_yyyyz[k] * ab_y + g_0_x_xzz_yyyyyz[k];

                g_0_x_xyzz_yyyzz[k] = -g_0_x_xzz_yyyzz[k] * ab_y + g_0_x_xzz_yyyyzz[k];

                g_0_x_xyzz_yyzzz[k] = -g_0_x_xzz_yyzzz[k] * ab_y + g_0_x_xzz_yyyzzz[k];

                g_0_x_xyzz_yzzzz[k] = -g_0_x_xzz_yzzzz[k] * ab_y + g_0_x_xzz_yyzzzz[k];

                g_0_x_xyzz_zzzzz[k] = -g_0_x_xzz_zzzzz[k] * ab_y + g_0_x_xzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xzz_xxxxx, g_0_x_xzz_xxxxxz, g_0_x_xzz_xxxxy, g_0_x_xzz_xxxxyz, g_0_x_xzz_xxxxz, g_0_x_xzz_xxxxzz, g_0_x_xzz_xxxyy, g_0_x_xzz_xxxyyz, g_0_x_xzz_xxxyz, g_0_x_xzz_xxxyzz, g_0_x_xzz_xxxzz, g_0_x_xzz_xxxzzz, g_0_x_xzz_xxyyy, g_0_x_xzz_xxyyyz, g_0_x_xzz_xxyyz, g_0_x_xzz_xxyyzz, g_0_x_xzz_xxyzz, g_0_x_xzz_xxyzzz, g_0_x_xzz_xxzzz, g_0_x_xzz_xxzzzz, g_0_x_xzz_xyyyy, g_0_x_xzz_xyyyyz, g_0_x_xzz_xyyyz, g_0_x_xzz_xyyyzz, g_0_x_xzz_xyyzz, g_0_x_xzz_xyyzzz, g_0_x_xzz_xyzzz, g_0_x_xzz_xyzzzz, g_0_x_xzz_xzzzz, g_0_x_xzz_xzzzzz, g_0_x_xzz_yyyyy, g_0_x_xzz_yyyyyz, g_0_x_xzz_yyyyz, g_0_x_xzz_yyyyzz, g_0_x_xzz_yyyzz, g_0_x_xzz_yyyzzz, g_0_x_xzz_yyzzz, g_0_x_xzz_yyzzzz, g_0_x_xzz_yzzzz, g_0_x_xzz_yzzzzz, g_0_x_xzz_zzzzz, g_0_x_xzz_zzzzzz, g_0_x_xzzz_xxxxx, g_0_x_xzzz_xxxxy, g_0_x_xzzz_xxxxz, g_0_x_xzzz_xxxyy, g_0_x_xzzz_xxxyz, g_0_x_xzzz_xxxzz, g_0_x_xzzz_xxyyy, g_0_x_xzzz_xxyyz, g_0_x_xzzz_xxyzz, g_0_x_xzzz_xxzzz, g_0_x_xzzz_xyyyy, g_0_x_xzzz_xyyyz, g_0_x_xzzz_xyyzz, g_0_x_xzzz_xyzzz, g_0_x_xzzz_xzzzz, g_0_x_xzzz_yyyyy, g_0_x_xzzz_yyyyz, g_0_x_xzzz_yyyzz, g_0_x_xzzz_yyzzz, g_0_x_xzzz_yzzzz, g_0_x_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_xxxxx[k] = -g_0_x_xzz_xxxxx[k] * ab_z + g_0_x_xzz_xxxxxz[k];

                g_0_x_xzzz_xxxxy[k] = -g_0_x_xzz_xxxxy[k] * ab_z + g_0_x_xzz_xxxxyz[k];

                g_0_x_xzzz_xxxxz[k] = -g_0_x_xzz_xxxxz[k] * ab_z + g_0_x_xzz_xxxxzz[k];

                g_0_x_xzzz_xxxyy[k] = -g_0_x_xzz_xxxyy[k] * ab_z + g_0_x_xzz_xxxyyz[k];

                g_0_x_xzzz_xxxyz[k] = -g_0_x_xzz_xxxyz[k] * ab_z + g_0_x_xzz_xxxyzz[k];

                g_0_x_xzzz_xxxzz[k] = -g_0_x_xzz_xxxzz[k] * ab_z + g_0_x_xzz_xxxzzz[k];

                g_0_x_xzzz_xxyyy[k] = -g_0_x_xzz_xxyyy[k] * ab_z + g_0_x_xzz_xxyyyz[k];

                g_0_x_xzzz_xxyyz[k] = -g_0_x_xzz_xxyyz[k] * ab_z + g_0_x_xzz_xxyyzz[k];

                g_0_x_xzzz_xxyzz[k] = -g_0_x_xzz_xxyzz[k] * ab_z + g_0_x_xzz_xxyzzz[k];

                g_0_x_xzzz_xxzzz[k] = -g_0_x_xzz_xxzzz[k] * ab_z + g_0_x_xzz_xxzzzz[k];

                g_0_x_xzzz_xyyyy[k] = -g_0_x_xzz_xyyyy[k] * ab_z + g_0_x_xzz_xyyyyz[k];

                g_0_x_xzzz_xyyyz[k] = -g_0_x_xzz_xyyyz[k] * ab_z + g_0_x_xzz_xyyyzz[k];

                g_0_x_xzzz_xyyzz[k] = -g_0_x_xzz_xyyzz[k] * ab_z + g_0_x_xzz_xyyzzz[k];

                g_0_x_xzzz_xyzzz[k] = -g_0_x_xzz_xyzzz[k] * ab_z + g_0_x_xzz_xyzzzz[k];

                g_0_x_xzzz_xzzzz[k] = -g_0_x_xzz_xzzzz[k] * ab_z + g_0_x_xzz_xzzzzz[k];

                g_0_x_xzzz_yyyyy[k] = -g_0_x_xzz_yyyyy[k] * ab_z + g_0_x_xzz_yyyyyz[k];

                g_0_x_xzzz_yyyyz[k] = -g_0_x_xzz_yyyyz[k] * ab_z + g_0_x_xzz_yyyyzz[k];

                g_0_x_xzzz_yyyzz[k] = -g_0_x_xzz_yyyzz[k] * ab_z + g_0_x_xzz_yyyzzz[k];

                g_0_x_xzzz_yyzzz[k] = -g_0_x_xzz_yyzzz[k] * ab_z + g_0_x_xzz_yyzzzz[k];

                g_0_x_xzzz_yzzzz[k] = -g_0_x_xzz_yzzzz[k] * ab_z + g_0_x_xzz_yzzzzz[k];

                g_0_x_xzzz_zzzzz[k] = -g_0_x_xzz_zzzzz[k] * ab_z + g_0_x_xzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyy_xxxxx, g_0_x_yyy_xxxxxy, g_0_x_yyy_xxxxy, g_0_x_yyy_xxxxyy, g_0_x_yyy_xxxxyz, g_0_x_yyy_xxxxz, g_0_x_yyy_xxxyy, g_0_x_yyy_xxxyyy, g_0_x_yyy_xxxyyz, g_0_x_yyy_xxxyz, g_0_x_yyy_xxxyzz, g_0_x_yyy_xxxzz, g_0_x_yyy_xxyyy, g_0_x_yyy_xxyyyy, g_0_x_yyy_xxyyyz, g_0_x_yyy_xxyyz, g_0_x_yyy_xxyyzz, g_0_x_yyy_xxyzz, g_0_x_yyy_xxyzzz, g_0_x_yyy_xxzzz, g_0_x_yyy_xyyyy, g_0_x_yyy_xyyyyy, g_0_x_yyy_xyyyyz, g_0_x_yyy_xyyyz, g_0_x_yyy_xyyyzz, g_0_x_yyy_xyyzz, g_0_x_yyy_xyyzzz, g_0_x_yyy_xyzzz, g_0_x_yyy_xyzzzz, g_0_x_yyy_xzzzz, g_0_x_yyy_yyyyy, g_0_x_yyy_yyyyyy, g_0_x_yyy_yyyyyz, g_0_x_yyy_yyyyz, g_0_x_yyy_yyyyzz, g_0_x_yyy_yyyzz, g_0_x_yyy_yyyzzz, g_0_x_yyy_yyzzz, g_0_x_yyy_yyzzzz, g_0_x_yyy_yzzzz, g_0_x_yyy_yzzzzz, g_0_x_yyy_zzzzz, g_0_x_yyyy_xxxxx, g_0_x_yyyy_xxxxy, g_0_x_yyyy_xxxxz, g_0_x_yyyy_xxxyy, g_0_x_yyyy_xxxyz, g_0_x_yyyy_xxxzz, g_0_x_yyyy_xxyyy, g_0_x_yyyy_xxyyz, g_0_x_yyyy_xxyzz, g_0_x_yyyy_xxzzz, g_0_x_yyyy_xyyyy, g_0_x_yyyy_xyyyz, g_0_x_yyyy_xyyzz, g_0_x_yyyy_xyzzz, g_0_x_yyyy_xzzzz, g_0_x_yyyy_yyyyy, g_0_x_yyyy_yyyyz, g_0_x_yyyy_yyyzz, g_0_x_yyyy_yyzzz, g_0_x_yyyy_yzzzz, g_0_x_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_xxxxx[k] = -g_0_x_yyy_xxxxx[k] * ab_y + g_0_x_yyy_xxxxxy[k];

                g_0_x_yyyy_xxxxy[k] = -g_0_x_yyy_xxxxy[k] * ab_y + g_0_x_yyy_xxxxyy[k];

                g_0_x_yyyy_xxxxz[k] = -g_0_x_yyy_xxxxz[k] * ab_y + g_0_x_yyy_xxxxyz[k];

                g_0_x_yyyy_xxxyy[k] = -g_0_x_yyy_xxxyy[k] * ab_y + g_0_x_yyy_xxxyyy[k];

                g_0_x_yyyy_xxxyz[k] = -g_0_x_yyy_xxxyz[k] * ab_y + g_0_x_yyy_xxxyyz[k];

                g_0_x_yyyy_xxxzz[k] = -g_0_x_yyy_xxxzz[k] * ab_y + g_0_x_yyy_xxxyzz[k];

                g_0_x_yyyy_xxyyy[k] = -g_0_x_yyy_xxyyy[k] * ab_y + g_0_x_yyy_xxyyyy[k];

                g_0_x_yyyy_xxyyz[k] = -g_0_x_yyy_xxyyz[k] * ab_y + g_0_x_yyy_xxyyyz[k];

                g_0_x_yyyy_xxyzz[k] = -g_0_x_yyy_xxyzz[k] * ab_y + g_0_x_yyy_xxyyzz[k];

                g_0_x_yyyy_xxzzz[k] = -g_0_x_yyy_xxzzz[k] * ab_y + g_0_x_yyy_xxyzzz[k];

                g_0_x_yyyy_xyyyy[k] = -g_0_x_yyy_xyyyy[k] * ab_y + g_0_x_yyy_xyyyyy[k];

                g_0_x_yyyy_xyyyz[k] = -g_0_x_yyy_xyyyz[k] * ab_y + g_0_x_yyy_xyyyyz[k];

                g_0_x_yyyy_xyyzz[k] = -g_0_x_yyy_xyyzz[k] * ab_y + g_0_x_yyy_xyyyzz[k];

                g_0_x_yyyy_xyzzz[k] = -g_0_x_yyy_xyzzz[k] * ab_y + g_0_x_yyy_xyyzzz[k];

                g_0_x_yyyy_xzzzz[k] = -g_0_x_yyy_xzzzz[k] * ab_y + g_0_x_yyy_xyzzzz[k];

                g_0_x_yyyy_yyyyy[k] = -g_0_x_yyy_yyyyy[k] * ab_y + g_0_x_yyy_yyyyyy[k];

                g_0_x_yyyy_yyyyz[k] = -g_0_x_yyy_yyyyz[k] * ab_y + g_0_x_yyy_yyyyyz[k];

                g_0_x_yyyy_yyyzz[k] = -g_0_x_yyy_yyyzz[k] * ab_y + g_0_x_yyy_yyyyzz[k];

                g_0_x_yyyy_yyzzz[k] = -g_0_x_yyy_yyzzz[k] * ab_y + g_0_x_yyy_yyyzzz[k];

                g_0_x_yyyy_yzzzz[k] = -g_0_x_yyy_yzzzz[k] * ab_y + g_0_x_yyy_yyzzzz[k];

                g_0_x_yyyy_zzzzz[k] = -g_0_x_yyy_zzzzz[k] * ab_y + g_0_x_yyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyyz_xxxxx, g_0_x_yyyz_xxxxy, g_0_x_yyyz_xxxxz, g_0_x_yyyz_xxxyy, g_0_x_yyyz_xxxyz, g_0_x_yyyz_xxxzz, g_0_x_yyyz_xxyyy, g_0_x_yyyz_xxyyz, g_0_x_yyyz_xxyzz, g_0_x_yyyz_xxzzz, g_0_x_yyyz_xyyyy, g_0_x_yyyz_xyyyz, g_0_x_yyyz_xyyzz, g_0_x_yyyz_xyzzz, g_0_x_yyyz_xzzzz, g_0_x_yyyz_yyyyy, g_0_x_yyyz_yyyyz, g_0_x_yyyz_yyyzz, g_0_x_yyyz_yyzzz, g_0_x_yyyz_yzzzz, g_0_x_yyyz_zzzzz, g_0_x_yyz_xxxxx, g_0_x_yyz_xxxxxy, g_0_x_yyz_xxxxy, g_0_x_yyz_xxxxyy, g_0_x_yyz_xxxxyz, g_0_x_yyz_xxxxz, g_0_x_yyz_xxxyy, g_0_x_yyz_xxxyyy, g_0_x_yyz_xxxyyz, g_0_x_yyz_xxxyz, g_0_x_yyz_xxxyzz, g_0_x_yyz_xxxzz, g_0_x_yyz_xxyyy, g_0_x_yyz_xxyyyy, g_0_x_yyz_xxyyyz, g_0_x_yyz_xxyyz, g_0_x_yyz_xxyyzz, g_0_x_yyz_xxyzz, g_0_x_yyz_xxyzzz, g_0_x_yyz_xxzzz, g_0_x_yyz_xyyyy, g_0_x_yyz_xyyyyy, g_0_x_yyz_xyyyyz, g_0_x_yyz_xyyyz, g_0_x_yyz_xyyyzz, g_0_x_yyz_xyyzz, g_0_x_yyz_xyyzzz, g_0_x_yyz_xyzzz, g_0_x_yyz_xyzzzz, g_0_x_yyz_xzzzz, g_0_x_yyz_yyyyy, g_0_x_yyz_yyyyyy, g_0_x_yyz_yyyyyz, g_0_x_yyz_yyyyz, g_0_x_yyz_yyyyzz, g_0_x_yyz_yyyzz, g_0_x_yyz_yyyzzz, g_0_x_yyz_yyzzz, g_0_x_yyz_yyzzzz, g_0_x_yyz_yzzzz, g_0_x_yyz_yzzzzz, g_0_x_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_xxxxx[k] = -g_0_x_yyz_xxxxx[k] * ab_y + g_0_x_yyz_xxxxxy[k];

                g_0_x_yyyz_xxxxy[k] = -g_0_x_yyz_xxxxy[k] * ab_y + g_0_x_yyz_xxxxyy[k];

                g_0_x_yyyz_xxxxz[k] = -g_0_x_yyz_xxxxz[k] * ab_y + g_0_x_yyz_xxxxyz[k];

                g_0_x_yyyz_xxxyy[k] = -g_0_x_yyz_xxxyy[k] * ab_y + g_0_x_yyz_xxxyyy[k];

                g_0_x_yyyz_xxxyz[k] = -g_0_x_yyz_xxxyz[k] * ab_y + g_0_x_yyz_xxxyyz[k];

                g_0_x_yyyz_xxxzz[k] = -g_0_x_yyz_xxxzz[k] * ab_y + g_0_x_yyz_xxxyzz[k];

                g_0_x_yyyz_xxyyy[k] = -g_0_x_yyz_xxyyy[k] * ab_y + g_0_x_yyz_xxyyyy[k];

                g_0_x_yyyz_xxyyz[k] = -g_0_x_yyz_xxyyz[k] * ab_y + g_0_x_yyz_xxyyyz[k];

                g_0_x_yyyz_xxyzz[k] = -g_0_x_yyz_xxyzz[k] * ab_y + g_0_x_yyz_xxyyzz[k];

                g_0_x_yyyz_xxzzz[k] = -g_0_x_yyz_xxzzz[k] * ab_y + g_0_x_yyz_xxyzzz[k];

                g_0_x_yyyz_xyyyy[k] = -g_0_x_yyz_xyyyy[k] * ab_y + g_0_x_yyz_xyyyyy[k];

                g_0_x_yyyz_xyyyz[k] = -g_0_x_yyz_xyyyz[k] * ab_y + g_0_x_yyz_xyyyyz[k];

                g_0_x_yyyz_xyyzz[k] = -g_0_x_yyz_xyyzz[k] * ab_y + g_0_x_yyz_xyyyzz[k];

                g_0_x_yyyz_xyzzz[k] = -g_0_x_yyz_xyzzz[k] * ab_y + g_0_x_yyz_xyyzzz[k];

                g_0_x_yyyz_xzzzz[k] = -g_0_x_yyz_xzzzz[k] * ab_y + g_0_x_yyz_xyzzzz[k];

                g_0_x_yyyz_yyyyy[k] = -g_0_x_yyz_yyyyy[k] * ab_y + g_0_x_yyz_yyyyyy[k];

                g_0_x_yyyz_yyyyz[k] = -g_0_x_yyz_yyyyz[k] * ab_y + g_0_x_yyz_yyyyyz[k];

                g_0_x_yyyz_yyyzz[k] = -g_0_x_yyz_yyyzz[k] * ab_y + g_0_x_yyz_yyyyzz[k];

                g_0_x_yyyz_yyzzz[k] = -g_0_x_yyz_yyzzz[k] * ab_y + g_0_x_yyz_yyyzzz[k];

                g_0_x_yyyz_yzzzz[k] = -g_0_x_yyz_yzzzz[k] * ab_y + g_0_x_yyz_yyzzzz[k];

                g_0_x_yyyz_zzzzz[k] = -g_0_x_yyz_zzzzz[k] * ab_y + g_0_x_yyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyzz_xxxxx, g_0_x_yyzz_xxxxy, g_0_x_yyzz_xxxxz, g_0_x_yyzz_xxxyy, g_0_x_yyzz_xxxyz, g_0_x_yyzz_xxxzz, g_0_x_yyzz_xxyyy, g_0_x_yyzz_xxyyz, g_0_x_yyzz_xxyzz, g_0_x_yyzz_xxzzz, g_0_x_yyzz_xyyyy, g_0_x_yyzz_xyyyz, g_0_x_yyzz_xyyzz, g_0_x_yyzz_xyzzz, g_0_x_yyzz_xzzzz, g_0_x_yyzz_yyyyy, g_0_x_yyzz_yyyyz, g_0_x_yyzz_yyyzz, g_0_x_yyzz_yyzzz, g_0_x_yyzz_yzzzz, g_0_x_yyzz_zzzzz, g_0_x_yzz_xxxxx, g_0_x_yzz_xxxxxy, g_0_x_yzz_xxxxy, g_0_x_yzz_xxxxyy, g_0_x_yzz_xxxxyz, g_0_x_yzz_xxxxz, g_0_x_yzz_xxxyy, g_0_x_yzz_xxxyyy, g_0_x_yzz_xxxyyz, g_0_x_yzz_xxxyz, g_0_x_yzz_xxxyzz, g_0_x_yzz_xxxzz, g_0_x_yzz_xxyyy, g_0_x_yzz_xxyyyy, g_0_x_yzz_xxyyyz, g_0_x_yzz_xxyyz, g_0_x_yzz_xxyyzz, g_0_x_yzz_xxyzz, g_0_x_yzz_xxyzzz, g_0_x_yzz_xxzzz, g_0_x_yzz_xyyyy, g_0_x_yzz_xyyyyy, g_0_x_yzz_xyyyyz, g_0_x_yzz_xyyyz, g_0_x_yzz_xyyyzz, g_0_x_yzz_xyyzz, g_0_x_yzz_xyyzzz, g_0_x_yzz_xyzzz, g_0_x_yzz_xyzzzz, g_0_x_yzz_xzzzz, g_0_x_yzz_yyyyy, g_0_x_yzz_yyyyyy, g_0_x_yzz_yyyyyz, g_0_x_yzz_yyyyz, g_0_x_yzz_yyyyzz, g_0_x_yzz_yyyzz, g_0_x_yzz_yyyzzz, g_0_x_yzz_yyzzz, g_0_x_yzz_yyzzzz, g_0_x_yzz_yzzzz, g_0_x_yzz_yzzzzz, g_0_x_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_xxxxx[k] = -g_0_x_yzz_xxxxx[k] * ab_y + g_0_x_yzz_xxxxxy[k];

                g_0_x_yyzz_xxxxy[k] = -g_0_x_yzz_xxxxy[k] * ab_y + g_0_x_yzz_xxxxyy[k];

                g_0_x_yyzz_xxxxz[k] = -g_0_x_yzz_xxxxz[k] * ab_y + g_0_x_yzz_xxxxyz[k];

                g_0_x_yyzz_xxxyy[k] = -g_0_x_yzz_xxxyy[k] * ab_y + g_0_x_yzz_xxxyyy[k];

                g_0_x_yyzz_xxxyz[k] = -g_0_x_yzz_xxxyz[k] * ab_y + g_0_x_yzz_xxxyyz[k];

                g_0_x_yyzz_xxxzz[k] = -g_0_x_yzz_xxxzz[k] * ab_y + g_0_x_yzz_xxxyzz[k];

                g_0_x_yyzz_xxyyy[k] = -g_0_x_yzz_xxyyy[k] * ab_y + g_0_x_yzz_xxyyyy[k];

                g_0_x_yyzz_xxyyz[k] = -g_0_x_yzz_xxyyz[k] * ab_y + g_0_x_yzz_xxyyyz[k];

                g_0_x_yyzz_xxyzz[k] = -g_0_x_yzz_xxyzz[k] * ab_y + g_0_x_yzz_xxyyzz[k];

                g_0_x_yyzz_xxzzz[k] = -g_0_x_yzz_xxzzz[k] * ab_y + g_0_x_yzz_xxyzzz[k];

                g_0_x_yyzz_xyyyy[k] = -g_0_x_yzz_xyyyy[k] * ab_y + g_0_x_yzz_xyyyyy[k];

                g_0_x_yyzz_xyyyz[k] = -g_0_x_yzz_xyyyz[k] * ab_y + g_0_x_yzz_xyyyyz[k];

                g_0_x_yyzz_xyyzz[k] = -g_0_x_yzz_xyyzz[k] * ab_y + g_0_x_yzz_xyyyzz[k];

                g_0_x_yyzz_xyzzz[k] = -g_0_x_yzz_xyzzz[k] * ab_y + g_0_x_yzz_xyyzzz[k];

                g_0_x_yyzz_xzzzz[k] = -g_0_x_yzz_xzzzz[k] * ab_y + g_0_x_yzz_xyzzzz[k];

                g_0_x_yyzz_yyyyy[k] = -g_0_x_yzz_yyyyy[k] * ab_y + g_0_x_yzz_yyyyyy[k];

                g_0_x_yyzz_yyyyz[k] = -g_0_x_yzz_yyyyz[k] * ab_y + g_0_x_yzz_yyyyyz[k];

                g_0_x_yyzz_yyyzz[k] = -g_0_x_yzz_yyyzz[k] * ab_y + g_0_x_yzz_yyyyzz[k];

                g_0_x_yyzz_yyzzz[k] = -g_0_x_yzz_yyzzz[k] * ab_y + g_0_x_yzz_yyyzzz[k];

                g_0_x_yyzz_yzzzz[k] = -g_0_x_yzz_yzzzz[k] * ab_y + g_0_x_yzz_yyzzzz[k];

                g_0_x_yyzz_zzzzz[k] = -g_0_x_yzz_zzzzz[k] * ab_y + g_0_x_yzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yzzz_xxxxx, g_0_x_yzzz_xxxxy, g_0_x_yzzz_xxxxz, g_0_x_yzzz_xxxyy, g_0_x_yzzz_xxxyz, g_0_x_yzzz_xxxzz, g_0_x_yzzz_xxyyy, g_0_x_yzzz_xxyyz, g_0_x_yzzz_xxyzz, g_0_x_yzzz_xxzzz, g_0_x_yzzz_xyyyy, g_0_x_yzzz_xyyyz, g_0_x_yzzz_xyyzz, g_0_x_yzzz_xyzzz, g_0_x_yzzz_xzzzz, g_0_x_yzzz_yyyyy, g_0_x_yzzz_yyyyz, g_0_x_yzzz_yyyzz, g_0_x_yzzz_yyzzz, g_0_x_yzzz_yzzzz, g_0_x_yzzz_zzzzz, g_0_x_zzz_xxxxx, g_0_x_zzz_xxxxxy, g_0_x_zzz_xxxxy, g_0_x_zzz_xxxxyy, g_0_x_zzz_xxxxyz, g_0_x_zzz_xxxxz, g_0_x_zzz_xxxyy, g_0_x_zzz_xxxyyy, g_0_x_zzz_xxxyyz, g_0_x_zzz_xxxyz, g_0_x_zzz_xxxyzz, g_0_x_zzz_xxxzz, g_0_x_zzz_xxyyy, g_0_x_zzz_xxyyyy, g_0_x_zzz_xxyyyz, g_0_x_zzz_xxyyz, g_0_x_zzz_xxyyzz, g_0_x_zzz_xxyzz, g_0_x_zzz_xxyzzz, g_0_x_zzz_xxzzz, g_0_x_zzz_xyyyy, g_0_x_zzz_xyyyyy, g_0_x_zzz_xyyyyz, g_0_x_zzz_xyyyz, g_0_x_zzz_xyyyzz, g_0_x_zzz_xyyzz, g_0_x_zzz_xyyzzz, g_0_x_zzz_xyzzz, g_0_x_zzz_xyzzzz, g_0_x_zzz_xzzzz, g_0_x_zzz_yyyyy, g_0_x_zzz_yyyyyy, g_0_x_zzz_yyyyyz, g_0_x_zzz_yyyyz, g_0_x_zzz_yyyyzz, g_0_x_zzz_yyyzz, g_0_x_zzz_yyyzzz, g_0_x_zzz_yyzzz, g_0_x_zzz_yyzzzz, g_0_x_zzz_yzzzz, g_0_x_zzz_yzzzzz, g_0_x_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_xxxxx[k] = -g_0_x_zzz_xxxxx[k] * ab_y + g_0_x_zzz_xxxxxy[k];

                g_0_x_yzzz_xxxxy[k] = -g_0_x_zzz_xxxxy[k] * ab_y + g_0_x_zzz_xxxxyy[k];

                g_0_x_yzzz_xxxxz[k] = -g_0_x_zzz_xxxxz[k] * ab_y + g_0_x_zzz_xxxxyz[k];

                g_0_x_yzzz_xxxyy[k] = -g_0_x_zzz_xxxyy[k] * ab_y + g_0_x_zzz_xxxyyy[k];

                g_0_x_yzzz_xxxyz[k] = -g_0_x_zzz_xxxyz[k] * ab_y + g_0_x_zzz_xxxyyz[k];

                g_0_x_yzzz_xxxzz[k] = -g_0_x_zzz_xxxzz[k] * ab_y + g_0_x_zzz_xxxyzz[k];

                g_0_x_yzzz_xxyyy[k] = -g_0_x_zzz_xxyyy[k] * ab_y + g_0_x_zzz_xxyyyy[k];

                g_0_x_yzzz_xxyyz[k] = -g_0_x_zzz_xxyyz[k] * ab_y + g_0_x_zzz_xxyyyz[k];

                g_0_x_yzzz_xxyzz[k] = -g_0_x_zzz_xxyzz[k] * ab_y + g_0_x_zzz_xxyyzz[k];

                g_0_x_yzzz_xxzzz[k] = -g_0_x_zzz_xxzzz[k] * ab_y + g_0_x_zzz_xxyzzz[k];

                g_0_x_yzzz_xyyyy[k] = -g_0_x_zzz_xyyyy[k] * ab_y + g_0_x_zzz_xyyyyy[k];

                g_0_x_yzzz_xyyyz[k] = -g_0_x_zzz_xyyyz[k] * ab_y + g_0_x_zzz_xyyyyz[k];

                g_0_x_yzzz_xyyzz[k] = -g_0_x_zzz_xyyzz[k] * ab_y + g_0_x_zzz_xyyyzz[k];

                g_0_x_yzzz_xyzzz[k] = -g_0_x_zzz_xyzzz[k] * ab_y + g_0_x_zzz_xyyzzz[k];

                g_0_x_yzzz_xzzzz[k] = -g_0_x_zzz_xzzzz[k] * ab_y + g_0_x_zzz_xyzzzz[k];

                g_0_x_yzzz_yyyyy[k] = -g_0_x_zzz_yyyyy[k] * ab_y + g_0_x_zzz_yyyyyy[k];

                g_0_x_yzzz_yyyyz[k] = -g_0_x_zzz_yyyyz[k] * ab_y + g_0_x_zzz_yyyyyz[k];

                g_0_x_yzzz_yyyzz[k] = -g_0_x_zzz_yyyzz[k] * ab_y + g_0_x_zzz_yyyyzz[k];

                g_0_x_yzzz_yyzzz[k] = -g_0_x_zzz_yyzzz[k] * ab_y + g_0_x_zzz_yyyzzz[k];

                g_0_x_yzzz_yzzzz[k] = -g_0_x_zzz_yzzzz[k] * ab_y + g_0_x_zzz_yyzzzz[k];

                g_0_x_yzzz_zzzzz[k] = -g_0_x_zzz_zzzzz[k] * ab_y + g_0_x_zzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_zzz_xxxxx, g_0_x_zzz_xxxxxz, g_0_x_zzz_xxxxy, g_0_x_zzz_xxxxyz, g_0_x_zzz_xxxxz, g_0_x_zzz_xxxxzz, g_0_x_zzz_xxxyy, g_0_x_zzz_xxxyyz, g_0_x_zzz_xxxyz, g_0_x_zzz_xxxyzz, g_0_x_zzz_xxxzz, g_0_x_zzz_xxxzzz, g_0_x_zzz_xxyyy, g_0_x_zzz_xxyyyz, g_0_x_zzz_xxyyz, g_0_x_zzz_xxyyzz, g_0_x_zzz_xxyzz, g_0_x_zzz_xxyzzz, g_0_x_zzz_xxzzz, g_0_x_zzz_xxzzzz, g_0_x_zzz_xyyyy, g_0_x_zzz_xyyyyz, g_0_x_zzz_xyyyz, g_0_x_zzz_xyyyzz, g_0_x_zzz_xyyzz, g_0_x_zzz_xyyzzz, g_0_x_zzz_xyzzz, g_0_x_zzz_xyzzzz, g_0_x_zzz_xzzzz, g_0_x_zzz_xzzzzz, g_0_x_zzz_yyyyy, g_0_x_zzz_yyyyyz, g_0_x_zzz_yyyyz, g_0_x_zzz_yyyyzz, g_0_x_zzz_yyyzz, g_0_x_zzz_yyyzzz, g_0_x_zzz_yyzzz, g_0_x_zzz_yyzzzz, g_0_x_zzz_yzzzz, g_0_x_zzz_yzzzzz, g_0_x_zzz_zzzzz, g_0_x_zzz_zzzzzz, g_0_x_zzzz_xxxxx, g_0_x_zzzz_xxxxy, g_0_x_zzzz_xxxxz, g_0_x_zzzz_xxxyy, g_0_x_zzzz_xxxyz, g_0_x_zzzz_xxxzz, g_0_x_zzzz_xxyyy, g_0_x_zzzz_xxyyz, g_0_x_zzzz_xxyzz, g_0_x_zzzz_xxzzz, g_0_x_zzzz_xyyyy, g_0_x_zzzz_xyyyz, g_0_x_zzzz_xyyzz, g_0_x_zzzz_xyzzz, g_0_x_zzzz_xzzzz, g_0_x_zzzz_yyyyy, g_0_x_zzzz_yyyyz, g_0_x_zzzz_yyyzz, g_0_x_zzzz_yyzzz, g_0_x_zzzz_yzzzz, g_0_x_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_xxxxx[k] = -g_0_x_zzz_xxxxx[k] * ab_z + g_0_x_zzz_xxxxxz[k];

                g_0_x_zzzz_xxxxy[k] = -g_0_x_zzz_xxxxy[k] * ab_z + g_0_x_zzz_xxxxyz[k];

                g_0_x_zzzz_xxxxz[k] = -g_0_x_zzz_xxxxz[k] * ab_z + g_0_x_zzz_xxxxzz[k];

                g_0_x_zzzz_xxxyy[k] = -g_0_x_zzz_xxxyy[k] * ab_z + g_0_x_zzz_xxxyyz[k];

                g_0_x_zzzz_xxxyz[k] = -g_0_x_zzz_xxxyz[k] * ab_z + g_0_x_zzz_xxxyzz[k];

                g_0_x_zzzz_xxxzz[k] = -g_0_x_zzz_xxxzz[k] * ab_z + g_0_x_zzz_xxxzzz[k];

                g_0_x_zzzz_xxyyy[k] = -g_0_x_zzz_xxyyy[k] * ab_z + g_0_x_zzz_xxyyyz[k];

                g_0_x_zzzz_xxyyz[k] = -g_0_x_zzz_xxyyz[k] * ab_z + g_0_x_zzz_xxyyzz[k];

                g_0_x_zzzz_xxyzz[k] = -g_0_x_zzz_xxyzz[k] * ab_z + g_0_x_zzz_xxyzzz[k];

                g_0_x_zzzz_xxzzz[k] = -g_0_x_zzz_xxzzz[k] * ab_z + g_0_x_zzz_xxzzzz[k];

                g_0_x_zzzz_xyyyy[k] = -g_0_x_zzz_xyyyy[k] * ab_z + g_0_x_zzz_xyyyyz[k];

                g_0_x_zzzz_xyyyz[k] = -g_0_x_zzz_xyyyz[k] * ab_z + g_0_x_zzz_xyyyzz[k];

                g_0_x_zzzz_xyyzz[k] = -g_0_x_zzz_xyyzz[k] * ab_z + g_0_x_zzz_xyyzzz[k];

                g_0_x_zzzz_xyzzz[k] = -g_0_x_zzz_xyzzz[k] * ab_z + g_0_x_zzz_xyzzzz[k];

                g_0_x_zzzz_xzzzz[k] = -g_0_x_zzz_xzzzz[k] * ab_z + g_0_x_zzz_xzzzzz[k];

                g_0_x_zzzz_yyyyy[k] = -g_0_x_zzz_yyyyy[k] * ab_z + g_0_x_zzz_yyyyyz[k];

                g_0_x_zzzz_yyyyz[k] = -g_0_x_zzz_yyyyz[k] * ab_z + g_0_x_zzz_yyyyzz[k];

                g_0_x_zzzz_yyyzz[k] = -g_0_x_zzz_yyyzz[k] * ab_z + g_0_x_zzz_yyyzzz[k];

                g_0_x_zzzz_yyzzz[k] = -g_0_x_zzz_yyzzz[k] * ab_z + g_0_x_zzz_yyzzzz[k];

                g_0_x_zzzz_yzzzz[k] = -g_0_x_zzz_yzzzz[k] * ab_z + g_0_x_zzz_yzzzzz[k];

                g_0_x_zzzz_zzzzz[k] = -g_0_x_zzz_zzzzz[k] * ab_z + g_0_x_zzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxx_xxxxx, g_0_y_xxx_xxxxxx, g_0_y_xxx_xxxxxy, g_0_y_xxx_xxxxxz, g_0_y_xxx_xxxxy, g_0_y_xxx_xxxxyy, g_0_y_xxx_xxxxyz, g_0_y_xxx_xxxxz, g_0_y_xxx_xxxxzz, g_0_y_xxx_xxxyy, g_0_y_xxx_xxxyyy, g_0_y_xxx_xxxyyz, g_0_y_xxx_xxxyz, g_0_y_xxx_xxxyzz, g_0_y_xxx_xxxzz, g_0_y_xxx_xxxzzz, g_0_y_xxx_xxyyy, g_0_y_xxx_xxyyyy, g_0_y_xxx_xxyyyz, g_0_y_xxx_xxyyz, g_0_y_xxx_xxyyzz, g_0_y_xxx_xxyzz, g_0_y_xxx_xxyzzz, g_0_y_xxx_xxzzz, g_0_y_xxx_xxzzzz, g_0_y_xxx_xyyyy, g_0_y_xxx_xyyyyy, g_0_y_xxx_xyyyyz, g_0_y_xxx_xyyyz, g_0_y_xxx_xyyyzz, g_0_y_xxx_xyyzz, g_0_y_xxx_xyyzzz, g_0_y_xxx_xyzzz, g_0_y_xxx_xyzzzz, g_0_y_xxx_xzzzz, g_0_y_xxx_xzzzzz, g_0_y_xxx_yyyyy, g_0_y_xxx_yyyyz, g_0_y_xxx_yyyzz, g_0_y_xxx_yyzzz, g_0_y_xxx_yzzzz, g_0_y_xxx_zzzzz, g_0_y_xxxx_xxxxx, g_0_y_xxxx_xxxxy, g_0_y_xxxx_xxxxz, g_0_y_xxxx_xxxyy, g_0_y_xxxx_xxxyz, g_0_y_xxxx_xxxzz, g_0_y_xxxx_xxyyy, g_0_y_xxxx_xxyyz, g_0_y_xxxx_xxyzz, g_0_y_xxxx_xxzzz, g_0_y_xxxx_xyyyy, g_0_y_xxxx_xyyyz, g_0_y_xxxx_xyyzz, g_0_y_xxxx_xyzzz, g_0_y_xxxx_xzzzz, g_0_y_xxxx_yyyyy, g_0_y_xxxx_yyyyz, g_0_y_xxxx_yyyzz, g_0_y_xxxx_yyzzz, g_0_y_xxxx_yzzzz, g_0_y_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_xxxxx[k] = -g_0_y_xxx_xxxxx[k] * ab_x + g_0_y_xxx_xxxxxx[k];

                g_0_y_xxxx_xxxxy[k] = -g_0_y_xxx_xxxxy[k] * ab_x + g_0_y_xxx_xxxxxy[k];

                g_0_y_xxxx_xxxxz[k] = -g_0_y_xxx_xxxxz[k] * ab_x + g_0_y_xxx_xxxxxz[k];

                g_0_y_xxxx_xxxyy[k] = -g_0_y_xxx_xxxyy[k] * ab_x + g_0_y_xxx_xxxxyy[k];

                g_0_y_xxxx_xxxyz[k] = -g_0_y_xxx_xxxyz[k] * ab_x + g_0_y_xxx_xxxxyz[k];

                g_0_y_xxxx_xxxzz[k] = -g_0_y_xxx_xxxzz[k] * ab_x + g_0_y_xxx_xxxxzz[k];

                g_0_y_xxxx_xxyyy[k] = -g_0_y_xxx_xxyyy[k] * ab_x + g_0_y_xxx_xxxyyy[k];

                g_0_y_xxxx_xxyyz[k] = -g_0_y_xxx_xxyyz[k] * ab_x + g_0_y_xxx_xxxyyz[k];

                g_0_y_xxxx_xxyzz[k] = -g_0_y_xxx_xxyzz[k] * ab_x + g_0_y_xxx_xxxyzz[k];

                g_0_y_xxxx_xxzzz[k] = -g_0_y_xxx_xxzzz[k] * ab_x + g_0_y_xxx_xxxzzz[k];

                g_0_y_xxxx_xyyyy[k] = -g_0_y_xxx_xyyyy[k] * ab_x + g_0_y_xxx_xxyyyy[k];

                g_0_y_xxxx_xyyyz[k] = -g_0_y_xxx_xyyyz[k] * ab_x + g_0_y_xxx_xxyyyz[k];

                g_0_y_xxxx_xyyzz[k] = -g_0_y_xxx_xyyzz[k] * ab_x + g_0_y_xxx_xxyyzz[k];

                g_0_y_xxxx_xyzzz[k] = -g_0_y_xxx_xyzzz[k] * ab_x + g_0_y_xxx_xxyzzz[k];

                g_0_y_xxxx_xzzzz[k] = -g_0_y_xxx_xzzzz[k] * ab_x + g_0_y_xxx_xxzzzz[k];

                g_0_y_xxxx_yyyyy[k] = -g_0_y_xxx_yyyyy[k] * ab_x + g_0_y_xxx_xyyyyy[k];

                g_0_y_xxxx_yyyyz[k] = -g_0_y_xxx_yyyyz[k] * ab_x + g_0_y_xxx_xyyyyz[k];

                g_0_y_xxxx_yyyzz[k] = -g_0_y_xxx_yyyzz[k] * ab_x + g_0_y_xxx_xyyyzz[k];

                g_0_y_xxxx_yyzzz[k] = -g_0_y_xxx_yyzzz[k] * ab_x + g_0_y_xxx_xyyzzz[k];

                g_0_y_xxxx_yzzzz[k] = -g_0_y_xxx_yzzzz[k] * ab_x + g_0_y_xxx_xyzzzz[k];

                g_0_y_xxxx_zzzzz[k] = -g_0_y_xxx_zzzzz[k] * ab_x + g_0_y_xxx_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxxy_xxxxx, g_0_y_xxxy_xxxxy, g_0_y_xxxy_xxxxz, g_0_y_xxxy_xxxyy, g_0_y_xxxy_xxxyz, g_0_y_xxxy_xxxzz, g_0_y_xxxy_xxyyy, g_0_y_xxxy_xxyyz, g_0_y_xxxy_xxyzz, g_0_y_xxxy_xxzzz, g_0_y_xxxy_xyyyy, g_0_y_xxxy_xyyyz, g_0_y_xxxy_xyyzz, g_0_y_xxxy_xyzzz, g_0_y_xxxy_xzzzz, g_0_y_xxxy_yyyyy, g_0_y_xxxy_yyyyz, g_0_y_xxxy_yyyzz, g_0_y_xxxy_yyzzz, g_0_y_xxxy_yzzzz, g_0_y_xxxy_zzzzz, g_0_y_xxy_xxxxx, g_0_y_xxy_xxxxxx, g_0_y_xxy_xxxxxy, g_0_y_xxy_xxxxxz, g_0_y_xxy_xxxxy, g_0_y_xxy_xxxxyy, g_0_y_xxy_xxxxyz, g_0_y_xxy_xxxxz, g_0_y_xxy_xxxxzz, g_0_y_xxy_xxxyy, g_0_y_xxy_xxxyyy, g_0_y_xxy_xxxyyz, g_0_y_xxy_xxxyz, g_0_y_xxy_xxxyzz, g_0_y_xxy_xxxzz, g_0_y_xxy_xxxzzz, g_0_y_xxy_xxyyy, g_0_y_xxy_xxyyyy, g_0_y_xxy_xxyyyz, g_0_y_xxy_xxyyz, g_0_y_xxy_xxyyzz, g_0_y_xxy_xxyzz, g_0_y_xxy_xxyzzz, g_0_y_xxy_xxzzz, g_0_y_xxy_xxzzzz, g_0_y_xxy_xyyyy, g_0_y_xxy_xyyyyy, g_0_y_xxy_xyyyyz, g_0_y_xxy_xyyyz, g_0_y_xxy_xyyyzz, g_0_y_xxy_xyyzz, g_0_y_xxy_xyyzzz, g_0_y_xxy_xyzzz, g_0_y_xxy_xyzzzz, g_0_y_xxy_xzzzz, g_0_y_xxy_xzzzzz, g_0_y_xxy_yyyyy, g_0_y_xxy_yyyyz, g_0_y_xxy_yyyzz, g_0_y_xxy_yyzzz, g_0_y_xxy_yzzzz, g_0_y_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_xxxxx[k] = -g_0_y_xxy_xxxxx[k] * ab_x + g_0_y_xxy_xxxxxx[k];

                g_0_y_xxxy_xxxxy[k] = -g_0_y_xxy_xxxxy[k] * ab_x + g_0_y_xxy_xxxxxy[k];

                g_0_y_xxxy_xxxxz[k] = -g_0_y_xxy_xxxxz[k] * ab_x + g_0_y_xxy_xxxxxz[k];

                g_0_y_xxxy_xxxyy[k] = -g_0_y_xxy_xxxyy[k] * ab_x + g_0_y_xxy_xxxxyy[k];

                g_0_y_xxxy_xxxyz[k] = -g_0_y_xxy_xxxyz[k] * ab_x + g_0_y_xxy_xxxxyz[k];

                g_0_y_xxxy_xxxzz[k] = -g_0_y_xxy_xxxzz[k] * ab_x + g_0_y_xxy_xxxxzz[k];

                g_0_y_xxxy_xxyyy[k] = -g_0_y_xxy_xxyyy[k] * ab_x + g_0_y_xxy_xxxyyy[k];

                g_0_y_xxxy_xxyyz[k] = -g_0_y_xxy_xxyyz[k] * ab_x + g_0_y_xxy_xxxyyz[k];

                g_0_y_xxxy_xxyzz[k] = -g_0_y_xxy_xxyzz[k] * ab_x + g_0_y_xxy_xxxyzz[k];

                g_0_y_xxxy_xxzzz[k] = -g_0_y_xxy_xxzzz[k] * ab_x + g_0_y_xxy_xxxzzz[k];

                g_0_y_xxxy_xyyyy[k] = -g_0_y_xxy_xyyyy[k] * ab_x + g_0_y_xxy_xxyyyy[k];

                g_0_y_xxxy_xyyyz[k] = -g_0_y_xxy_xyyyz[k] * ab_x + g_0_y_xxy_xxyyyz[k];

                g_0_y_xxxy_xyyzz[k] = -g_0_y_xxy_xyyzz[k] * ab_x + g_0_y_xxy_xxyyzz[k];

                g_0_y_xxxy_xyzzz[k] = -g_0_y_xxy_xyzzz[k] * ab_x + g_0_y_xxy_xxyzzz[k];

                g_0_y_xxxy_xzzzz[k] = -g_0_y_xxy_xzzzz[k] * ab_x + g_0_y_xxy_xxzzzz[k];

                g_0_y_xxxy_yyyyy[k] = -g_0_y_xxy_yyyyy[k] * ab_x + g_0_y_xxy_xyyyyy[k];

                g_0_y_xxxy_yyyyz[k] = -g_0_y_xxy_yyyyz[k] * ab_x + g_0_y_xxy_xyyyyz[k];

                g_0_y_xxxy_yyyzz[k] = -g_0_y_xxy_yyyzz[k] * ab_x + g_0_y_xxy_xyyyzz[k];

                g_0_y_xxxy_yyzzz[k] = -g_0_y_xxy_yyzzz[k] * ab_x + g_0_y_xxy_xyyzzz[k];

                g_0_y_xxxy_yzzzz[k] = -g_0_y_xxy_yzzzz[k] * ab_x + g_0_y_xxy_xyzzzz[k];

                g_0_y_xxxy_zzzzz[k] = -g_0_y_xxy_zzzzz[k] * ab_x + g_0_y_xxy_xzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxxz_xxxxx, g_0_y_xxxz_xxxxy, g_0_y_xxxz_xxxxz, g_0_y_xxxz_xxxyy, g_0_y_xxxz_xxxyz, g_0_y_xxxz_xxxzz, g_0_y_xxxz_xxyyy, g_0_y_xxxz_xxyyz, g_0_y_xxxz_xxyzz, g_0_y_xxxz_xxzzz, g_0_y_xxxz_xyyyy, g_0_y_xxxz_xyyyz, g_0_y_xxxz_xyyzz, g_0_y_xxxz_xyzzz, g_0_y_xxxz_xzzzz, g_0_y_xxxz_yyyyy, g_0_y_xxxz_yyyyz, g_0_y_xxxz_yyyzz, g_0_y_xxxz_yyzzz, g_0_y_xxxz_yzzzz, g_0_y_xxxz_zzzzz, g_0_y_xxz_xxxxx, g_0_y_xxz_xxxxxx, g_0_y_xxz_xxxxxy, g_0_y_xxz_xxxxxz, g_0_y_xxz_xxxxy, g_0_y_xxz_xxxxyy, g_0_y_xxz_xxxxyz, g_0_y_xxz_xxxxz, g_0_y_xxz_xxxxzz, g_0_y_xxz_xxxyy, g_0_y_xxz_xxxyyy, g_0_y_xxz_xxxyyz, g_0_y_xxz_xxxyz, g_0_y_xxz_xxxyzz, g_0_y_xxz_xxxzz, g_0_y_xxz_xxxzzz, g_0_y_xxz_xxyyy, g_0_y_xxz_xxyyyy, g_0_y_xxz_xxyyyz, g_0_y_xxz_xxyyz, g_0_y_xxz_xxyyzz, g_0_y_xxz_xxyzz, g_0_y_xxz_xxyzzz, g_0_y_xxz_xxzzz, g_0_y_xxz_xxzzzz, g_0_y_xxz_xyyyy, g_0_y_xxz_xyyyyy, g_0_y_xxz_xyyyyz, g_0_y_xxz_xyyyz, g_0_y_xxz_xyyyzz, g_0_y_xxz_xyyzz, g_0_y_xxz_xyyzzz, g_0_y_xxz_xyzzz, g_0_y_xxz_xyzzzz, g_0_y_xxz_xzzzz, g_0_y_xxz_xzzzzz, g_0_y_xxz_yyyyy, g_0_y_xxz_yyyyz, g_0_y_xxz_yyyzz, g_0_y_xxz_yyzzz, g_0_y_xxz_yzzzz, g_0_y_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_xxxxx[k] = -g_0_y_xxz_xxxxx[k] * ab_x + g_0_y_xxz_xxxxxx[k];

                g_0_y_xxxz_xxxxy[k] = -g_0_y_xxz_xxxxy[k] * ab_x + g_0_y_xxz_xxxxxy[k];

                g_0_y_xxxz_xxxxz[k] = -g_0_y_xxz_xxxxz[k] * ab_x + g_0_y_xxz_xxxxxz[k];

                g_0_y_xxxz_xxxyy[k] = -g_0_y_xxz_xxxyy[k] * ab_x + g_0_y_xxz_xxxxyy[k];

                g_0_y_xxxz_xxxyz[k] = -g_0_y_xxz_xxxyz[k] * ab_x + g_0_y_xxz_xxxxyz[k];

                g_0_y_xxxz_xxxzz[k] = -g_0_y_xxz_xxxzz[k] * ab_x + g_0_y_xxz_xxxxzz[k];

                g_0_y_xxxz_xxyyy[k] = -g_0_y_xxz_xxyyy[k] * ab_x + g_0_y_xxz_xxxyyy[k];

                g_0_y_xxxz_xxyyz[k] = -g_0_y_xxz_xxyyz[k] * ab_x + g_0_y_xxz_xxxyyz[k];

                g_0_y_xxxz_xxyzz[k] = -g_0_y_xxz_xxyzz[k] * ab_x + g_0_y_xxz_xxxyzz[k];

                g_0_y_xxxz_xxzzz[k] = -g_0_y_xxz_xxzzz[k] * ab_x + g_0_y_xxz_xxxzzz[k];

                g_0_y_xxxz_xyyyy[k] = -g_0_y_xxz_xyyyy[k] * ab_x + g_0_y_xxz_xxyyyy[k];

                g_0_y_xxxz_xyyyz[k] = -g_0_y_xxz_xyyyz[k] * ab_x + g_0_y_xxz_xxyyyz[k];

                g_0_y_xxxz_xyyzz[k] = -g_0_y_xxz_xyyzz[k] * ab_x + g_0_y_xxz_xxyyzz[k];

                g_0_y_xxxz_xyzzz[k] = -g_0_y_xxz_xyzzz[k] * ab_x + g_0_y_xxz_xxyzzz[k];

                g_0_y_xxxz_xzzzz[k] = -g_0_y_xxz_xzzzz[k] * ab_x + g_0_y_xxz_xxzzzz[k];

                g_0_y_xxxz_yyyyy[k] = -g_0_y_xxz_yyyyy[k] * ab_x + g_0_y_xxz_xyyyyy[k];

                g_0_y_xxxz_yyyyz[k] = -g_0_y_xxz_yyyyz[k] * ab_x + g_0_y_xxz_xyyyyz[k];

                g_0_y_xxxz_yyyzz[k] = -g_0_y_xxz_yyyzz[k] * ab_x + g_0_y_xxz_xyyyzz[k];

                g_0_y_xxxz_yyzzz[k] = -g_0_y_xxz_yyzzz[k] * ab_x + g_0_y_xxz_xyyzzz[k];

                g_0_y_xxxz_yzzzz[k] = -g_0_y_xxz_yzzzz[k] * ab_x + g_0_y_xxz_xyzzzz[k];

                g_0_y_xxxz_zzzzz[k] = -g_0_y_xxz_zzzzz[k] * ab_x + g_0_y_xxz_xzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxyy_xxxxx, g_0_y_xxyy_xxxxy, g_0_y_xxyy_xxxxz, g_0_y_xxyy_xxxyy, g_0_y_xxyy_xxxyz, g_0_y_xxyy_xxxzz, g_0_y_xxyy_xxyyy, g_0_y_xxyy_xxyyz, g_0_y_xxyy_xxyzz, g_0_y_xxyy_xxzzz, g_0_y_xxyy_xyyyy, g_0_y_xxyy_xyyyz, g_0_y_xxyy_xyyzz, g_0_y_xxyy_xyzzz, g_0_y_xxyy_xzzzz, g_0_y_xxyy_yyyyy, g_0_y_xxyy_yyyyz, g_0_y_xxyy_yyyzz, g_0_y_xxyy_yyzzz, g_0_y_xxyy_yzzzz, g_0_y_xxyy_zzzzz, g_0_y_xyy_xxxxx, g_0_y_xyy_xxxxxx, g_0_y_xyy_xxxxxy, g_0_y_xyy_xxxxxz, g_0_y_xyy_xxxxy, g_0_y_xyy_xxxxyy, g_0_y_xyy_xxxxyz, g_0_y_xyy_xxxxz, g_0_y_xyy_xxxxzz, g_0_y_xyy_xxxyy, g_0_y_xyy_xxxyyy, g_0_y_xyy_xxxyyz, g_0_y_xyy_xxxyz, g_0_y_xyy_xxxyzz, g_0_y_xyy_xxxzz, g_0_y_xyy_xxxzzz, g_0_y_xyy_xxyyy, g_0_y_xyy_xxyyyy, g_0_y_xyy_xxyyyz, g_0_y_xyy_xxyyz, g_0_y_xyy_xxyyzz, g_0_y_xyy_xxyzz, g_0_y_xyy_xxyzzz, g_0_y_xyy_xxzzz, g_0_y_xyy_xxzzzz, g_0_y_xyy_xyyyy, g_0_y_xyy_xyyyyy, g_0_y_xyy_xyyyyz, g_0_y_xyy_xyyyz, g_0_y_xyy_xyyyzz, g_0_y_xyy_xyyzz, g_0_y_xyy_xyyzzz, g_0_y_xyy_xyzzz, g_0_y_xyy_xyzzzz, g_0_y_xyy_xzzzz, g_0_y_xyy_xzzzzz, g_0_y_xyy_yyyyy, g_0_y_xyy_yyyyz, g_0_y_xyy_yyyzz, g_0_y_xyy_yyzzz, g_0_y_xyy_yzzzz, g_0_y_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_xxxxx[k] = -g_0_y_xyy_xxxxx[k] * ab_x + g_0_y_xyy_xxxxxx[k];

                g_0_y_xxyy_xxxxy[k] = -g_0_y_xyy_xxxxy[k] * ab_x + g_0_y_xyy_xxxxxy[k];

                g_0_y_xxyy_xxxxz[k] = -g_0_y_xyy_xxxxz[k] * ab_x + g_0_y_xyy_xxxxxz[k];

                g_0_y_xxyy_xxxyy[k] = -g_0_y_xyy_xxxyy[k] * ab_x + g_0_y_xyy_xxxxyy[k];

                g_0_y_xxyy_xxxyz[k] = -g_0_y_xyy_xxxyz[k] * ab_x + g_0_y_xyy_xxxxyz[k];

                g_0_y_xxyy_xxxzz[k] = -g_0_y_xyy_xxxzz[k] * ab_x + g_0_y_xyy_xxxxzz[k];

                g_0_y_xxyy_xxyyy[k] = -g_0_y_xyy_xxyyy[k] * ab_x + g_0_y_xyy_xxxyyy[k];

                g_0_y_xxyy_xxyyz[k] = -g_0_y_xyy_xxyyz[k] * ab_x + g_0_y_xyy_xxxyyz[k];

                g_0_y_xxyy_xxyzz[k] = -g_0_y_xyy_xxyzz[k] * ab_x + g_0_y_xyy_xxxyzz[k];

                g_0_y_xxyy_xxzzz[k] = -g_0_y_xyy_xxzzz[k] * ab_x + g_0_y_xyy_xxxzzz[k];

                g_0_y_xxyy_xyyyy[k] = -g_0_y_xyy_xyyyy[k] * ab_x + g_0_y_xyy_xxyyyy[k];

                g_0_y_xxyy_xyyyz[k] = -g_0_y_xyy_xyyyz[k] * ab_x + g_0_y_xyy_xxyyyz[k];

                g_0_y_xxyy_xyyzz[k] = -g_0_y_xyy_xyyzz[k] * ab_x + g_0_y_xyy_xxyyzz[k];

                g_0_y_xxyy_xyzzz[k] = -g_0_y_xyy_xyzzz[k] * ab_x + g_0_y_xyy_xxyzzz[k];

                g_0_y_xxyy_xzzzz[k] = -g_0_y_xyy_xzzzz[k] * ab_x + g_0_y_xyy_xxzzzz[k];

                g_0_y_xxyy_yyyyy[k] = -g_0_y_xyy_yyyyy[k] * ab_x + g_0_y_xyy_xyyyyy[k];

                g_0_y_xxyy_yyyyz[k] = -g_0_y_xyy_yyyyz[k] * ab_x + g_0_y_xyy_xyyyyz[k];

                g_0_y_xxyy_yyyzz[k] = -g_0_y_xyy_yyyzz[k] * ab_x + g_0_y_xyy_xyyyzz[k];

                g_0_y_xxyy_yyzzz[k] = -g_0_y_xyy_yyzzz[k] * ab_x + g_0_y_xyy_xyyzzz[k];

                g_0_y_xxyy_yzzzz[k] = -g_0_y_xyy_yzzzz[k] * ab_x + g_0_y_xyy_xyzzzz[k];

                g_0_y_xxyy_zzzzz[k] = -g_0_y_xyy_zzzzz[k] * ab_x + g_0_y_xyy_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxyz_xxxxx, g_0_y_xxyz_xxxxy, g_0_y_xxyz_xxxxz, g_0_y_xxyz_xxxyy, g_0_y_xxyz_xxxyz, g_0_y_xxyz_xxxzz, g_0_y_xxyz_xxyyy, g_0_y_xxyz_xxyyz, g_0_y_xxyz_xxyzz, g_0_y_xxyz_xxzzz, g_0_y_xxyz_xyyyy, g_0_y_xxyz_xyyyz, g_0_y_xxyz_xyyzz, g_0_y_xxyz_xyzzz, g_0_y_xxyz_xzzzz, g_0_y_xxyz_yyyyy, g_0_y_xxyz_yyyyz, g_0_y_xxyz_yyyzz, g_0_y_xxyz_yyzzz, g_0_y_xxyz_yzzzz, g_0_y_xxyz_zzzzz, g_0_y_xyz_xxxxx, g_0_y_xyz_xxxxxx, g_0_y_xyz_xxxxxy, g_0_y_xyz_xxxxxz, g_0_y_xyz_xxxxy, g_0_y_xyz_xxxxyy, g_0_y_xyz_xxxxyz, g_0_y_xyz_xxxxz, g_0_y_xyz_xxxxzz, g_0_y_xyz_xxxyy, g_0_y_xyz_xxxyyy, g_0_y_xyz_xxxyyz, g_0_y_xyz_xxxyz, g_0_y_xyz_xxxyzz, g_0_y_xyz_xxxzz, g_0_y_xyz_xxxzzz, g_0_y_xyz_xxyyy, g_0_y_xyz_xxyyyy, g_0_y_xyz_xxyyyz, g_0_y_xyz_xxyyz, g_0_y_xyz_xxyyzz, g_0_y_xyz_xxyzz, g_0_y_xyz_xxyzzz, g_0_y_xyz_xxzzz, g_0_y_xyz_xxzzzz, g_0_y_xyz_xyyyy, g_0_y_xyz_xyyyyy, g_0_y_xyz_xyyyyz, g_0_y_xyz_xyyyz, g_0_y_xyz_xyyyzz, g_0_y_xyz_xyyzz, g_0_y_xyz_xyyzzz, g_0_y_xyz_xyzzz, g_0_y_xyz_xyzzzz, g_0_y_xyz_xzzzz, g_0_y_xyz_xzzzzz, g_0_y_xyz_yyyyy, g_0_y_xyz_yyyyz, g_0_y_xyz_yyyzz, g_0_y_xyz_yyzzz, g_0_y_xyz_yzzzz, g_0_y_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_xxxxx[k] = -g_0_y_xyz_xxxxx[k] * ab_x + g_0_y_xyz_xxxxxx[k];

                g_0_y_xxyz_xxxxy[k] = -g_0_y_xyz_xxxxy[k] * ab_x + g_0_y_xyz_xxxxxy[k];

                g_0_y_xxyz_xxxxz[k] = -g_0_y_xyz_xxxxz[k] * ab_x + g_0_y_xyz_xxxxxz[k];

                g_0_y_xxyz_xxxyy[k] = -g_0_y_xyz_xxxyy[k] * ab_x + g_0_y_xyz_xxxxyy[k];

                g_0_y_xxyz_xxxyz[k] = -g_0_y_xyz_xxxyz[k] * ab_x + g_0_y_xyz_xxxxyz[k];

                g_0_y_xxyz_xxxzz[k] = -g_0_y_xyz_xxxzz[k] * ab_x + g_0_y_xyz_xxxxzz[k];

                g_0_y_xxyz_xxyyy[k] = -g_0_y_xyz_xxyyy[k] * ab_x + g_0_y_xyz_xxxyyy[k];

                g_0_y_xxyz_xxyyz[k] = -g_0_y_xyz_xxyyz[k] * ab_x + g_0_y_xyz_xxxyyz[k];

                g_0_y_xxyz_xxyzz[k] = -g_0_y_xyz_xxyzz[k] * ab_x + g_0_y_xyz_xxxyzz[k];

                g_0_y_xxyz_xxzzz[k] = -g_0_y_xyz_xxzzz[k] * ab_x + g_0_y_xyz_xxxzzz[k];

                g_0_y_xxyz_xyyyy[k] = -g_0_y_xyz_xyyyy[k] * ab_x + g_0_y_xyz_xxyyyy[k];

                g_0_y_xxyz_xyyyz[k] = -g_0_y_xyz_xyyyz[k] * ab_x + g_0_y_xyz_xxyyyz[k];

                g_0_y_xxyz_xyyzz[k] = -g_0_y_xyz_xyyzz[k] * ab_x + g_0_y_xyz_xxyyzz[k];

                g_0_y_xxyz_xyzzz[k] = -g_0_y_xyz_xyzzz[k] * ab_x + g_0_y_xyz_xxyzzz[k];

                g_0_y_xxyz_xzzzz[k] = -g_0_y_xyz_xzzzz[k] * ab_x + g_0_y_xyz_xxzzzz[k];

                g_0_y_xxyz_yyyyy[k] = -g_0_y_xyz_yyyyy[k] * ab_x + g_0_y_xyz_xyyyyy[k];

                g_0_y_xxyz_yyyyz[k] = -g_0_y_xyz_yyyyz[k] * ab_x + g_0_y_xyz_xyyyyz[k];

                g_0_y_xxyz_yyyzz[k] = -g_0_y_xyz_yyyzz[k] * ab_x + g_0_y_xyz_xyyyzz[k];

                g_0_y_xxyz_yyzzz[k] = -g_0_y_xyz_yyzzz[k] * ab_x + g_0_y_xyz_xyyzzz[k];

                g_0_y_xxyz_yzzzz[k] = -g_0_y_xyz_yzzzz[k] * ab_x + g_0_y_xyz_xyzzzz[k];

                g_0_y_xxyz_zzzzz[k] = -g_0_y_xyz_zzzzz[k] * ab_x + g_0_y_xyz_xzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxzz_xxxxx, g_0_y_xxzz_xxxxy, g_0_y_xxzz_xxxxz, g_0_y_xxzz_xxxyy, g_0_y_xxzz_xxxyz, g_0_y_xxzz_xxxzz, g_0_y_xxzz_xxyyy, g_0_y_xxzz_xxyyz, g_0_y_xxzz_xxyzz, g_0_y_xxzz_xxzzz, g_0_y_xxzz_xyyyy, g_0_y_xxzz_xyyyz, g_0_y_xxzz_xyyzz, g_0_y_xxzz_xyzzz, g_0_y_xxzz_xzzzz, g_0_y_xxzz_yyyyy, g_0_y_xxzz_yyyyz, g_0_y_xxzz_yyyzz, g_0_y_xxzz_yyzzz, g_0_y_xxzz_yzzzz, g_0_y_xxzz_zzzzz, g_0_y_xzz_xxxxx, g_0_y_xzz_xxxxxx, g_0_y_xzz_xxxxxy, g_0_y_xzz_xxxxxz, g_0_y_xzz_xxxxy, g_0_y_xzz_xxxxyy, g_0_y_xzz_xxxxyz, g_0_y_xzz_xxxxz, g_0_y_xzz_xxxxzz, g_0_y_xzz_xxxyy, g_0_y_xzz_xxxyyy, g_0_y_xzz_xxxyyz, g_0_y_xzz_xxxyz, g_0_y_xzz_xxxyzz, g_0_y_xzz_xxxzz, g_0_y_xzz_xxxzzz, g_0_y_xzz_xxyyy, g_0_y_xzz_xxyyyy, g_0_y_xzz_xxyyyz, g_0_y_xzz_xxyyz, g_0_y_xzz_xxyyzz, g_0_y_xzz_xxyzz, g_0_y_xzz_xxyzzz, g_0_y_xzz_xxzzz, g_0_y_xzz_xxzzzz, g_0_y_xzz_xyyyy, g_0_y_xzz_xyyyyy, g_0_y_xzz_xyyyyz, g_0_y_xzz_xyyyz, g_0_y_xzz_xyyyzz, g_0_y_xzz_xyyzz, g_0_y_xzz_xyyzzz, g_0_y_xzz_xyzzz, g_0_y_xzz_xyzzzz, g_0_y_xzz_xzzzz, g_0_y_xzz_xzzzzz, g_0_y_xzz_yyyyy, g_0_y_xzz_yyyyz, g_0_y_xzz_yyyzz, g_0_y_xzz_yyzzz, g_0_y_xzz_yzzzz, g_0_y_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_xxxxx[k] = -g_0_y_xzz_xxxxx[k] * ab_x + g_0_y_xzz_xxxxxx[k];

                g_0_y_xxzz_xxxxy[k] = -g_0_y_xzz_xxxxy[k] * ab_x + g_0_y_xzz_xxxxxy[k];

                g_0_y_xxzz_xxxxz[k] = -g_0_y_xzz_xxxxz[k] * ab_x + g_0_y_xzz_xxxxxz[k];

                g_0_y_xxzz_xxxyy[k] = -g_0_y_xzz_xxxyy[k] * ab_x + g_0_y_xzz_xxxxyy[k];

                g_0_y_xxzz_xxxyz[k] = -g_0_y_xzz_xxxyz[k] * ab_x + g_0_y_xzz_xxxxyz[k];

                g_0_y_xxzz_xxxzz[k] = -g_0_y_xzz_xxxzz[k] * ab_x + g_0_y_xzz_xxxxzz[k];

                g_0_y_xxzz_xxyyy[k] = -g_0_y_xzz_xxyyy[k] * ab_x + g_0_y_xzz_xxxyyy[k];

                g_0_y_xxzz_xxyyz[k] = -g_0_y_xzz_xxyyz[k] * ab_x + g_0_y_xzz_xxxyyz[k];

                g_0_y_xxzz_xxyzz[k] = -g_0_y_xzz_xxyzz[k] * ab_x + g_0_y_xzz_xxxyzz[k];

                g_0_y_xxzz_xxzzz[k] = -g_0_y_xzz_xxzzz[k] * ab_x + g_0_y_xzz_xxxzzz[k];

                g_0_y_xxzz_xyyyy[k] = -g_0_y_xzz_xyyyy[k] * ab_x + g_0_y_xzz_xxyyyy[k];

                g_0_y_xxzz_xyyyz[k] = -g_0_y_xzz_xyyyz[k] * ab_x + g_0_y_xzz_xxyyyz[k];

                g_0_y_xxzz_xyyzz[k] = -g_0_y_xzz_xyyzz[k] * ab_x + g_0_y_xzz_xxyyzz[k];

                g_0_y_xxzz_xyzzz[k] = -g_0_y_xzz_xyzzz[k] * ab_x + g_0_y_xzz_xxyzzz[k];

                g_0_y_xxzz_xzzzz[k] = -g_0_y_xzz_xzzzz[k] * ab_x + g_0_y_xzz_xxzzzz[k];

                g_0_y_xxzz_yyyyy[k] = -g_0_y_xzz_yyyyy[k] * ab_x + g_0_y_xzz_xyyyyy[k];

                g_0_y_xxzz_yyyyz[k] = -g_0_y_xzz_yyyyz[k] * ab_x + g_0_y_xzz_xyyyyz[k];

                g_0_y_xxzz_yyyzz[k] = -g_0_y_xzz_yyyzz[k] * ab_x + g_0_y_xzz_xyyyzz[k];

                g_0_y_xxzz_yyzzz[k] = -g_0_y_xzz_yyzzz[k] * ab_x + g_0_y_xzz_xyyzzz[k];

                g_0_y_xxzz_yzzzz[k] = -g_0_y_xzz_yzzzz[k] * ab_x + g_0_y_xzz_xyzzzz[k];

                g_0_y_xxzz_zzzzz[k] = -g_0_y_xzz_zzzzz[k] * ab_x + g_0_y_xzz_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyyy_xxxxx, g_0_y_xyyy_xxxxy, g_0_y_xyyy_xxxxz, g_0_y_xyyy_xxxyy, g_0_y_xyyy_xxxyz, g_0_y_xyyy_xxxzz, g_0_y_xyyy_xxyyy, g_0_y_xyyy_xxyyz, g_0_y_xyyy_xxyzz, g_0_y_xyyy_xxzzz, g_0_y_xyyy_xyyyy, g_0_y_xyyy_xyyyz, g_0_y_xyyy_xyyzz, g_0_y_xyyy_xyzzz, g_0_y_xyyy_xzzzz, g_0_y_xyyy_yyyyy, g_0_y_xyyy_yyyyz, g_0_y_xyyy_yyyzz, g_0_y_xyyy_yyzzz, g_0_y_xyyy_yzzzz, g_0_y_xyyy_zzzzz, g_0_y_yyy_xxxxx, g_0_y_yyy_xxxxxx, g_0_y_yyy_xxxxxy, g_0_y_yyy_xxxxxz, g_0_y_yyy_xxxxy, g_0_y_yyy_xxxxyy, g_0_y_yyy_xxxxyz, g_0_y_yyy_xxxxz, g_0_y_yyy_xxxxzz, g_0_y_yyy_xxxyy, g_0_y_yyy_xxxyyy, g_0_y_yyy_xxxyyz, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxyzz, g_0_y_yyy_xxxzz, g_0_y_yyy_xxxzzz, g_0_y_yyy_xxyyy, g_0_y_yyy_xxyyyy, g_0_y_yyy_xxyyyz, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyyzz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxyzzz, g_0_y_yyy_xxzzz, g_0_y_yyy_xxzzzz, g_0_y_yyy_xyyyy, g_0_y_yyy_xyyyyy, g_0_y_yyy_xyyyyz, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyyzz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyyzzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xyzzzz, g_0_y_yyy_xzzzz, g_0_y_yyy_xzzzzz, g_0_y_yyy_yyyyy, g_0_y_yyy_yyyyz, g_0_y_yyy_yyyzz, g_0_y_yyy_yyzzz, g_0_y_yyy_yzzzz, g_0_y_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_xxxxx[k] = -g_0_y_yyy_xxxxx[k] * ab_x + g_0_y_yyy_xxxxxx[k];

                g_0_y_xyyy_xxxxy[k] = -g_0_y_yyy_xxxxy[k] * ab_x + g_0_y_yyy_xxxxxy[k];

                g_0_y_xyyy_xxxxz[k] = -g_0_y_yyy_xxxxz[k] * ab_x + g_0_y_yyy_xxxxxz[k];

                g_0_y_xyyy_xxxyy[k] = -g_0_y_yyy_xxxyy[k] * ab_x + g_0_y_yyy_xxxxyy[k];

                g_0_y_xyyy_xxxyz[k] = -g_0_y_yyy_xxxyz[k] * ab_x + g_0_y_yyy_xxxxyz[k];

                g_0_y_xyyy_xxxzz[k] = -g_0_y_yyy_xxxzz[k] * ab_x + g_0_y_yyy_xxxxzz[k];

                g_0_y_xyyy_xxyyy[k] = -g_0_y_yyy_xxyyy[k] * ab_x + g_0_y_yyy_xxxyyy[k];

                g_0_y_xyyy_xxyyz[k] = -g_0_y_yyy_xxyyz[k] * ab_x + g_0_y_yyy_xxxyyz[k];

                g_0_y_xyyy_xxyzz[k] = -g_0_y_yyy_xxyzz[k] * ab_x + g_0_y_yyy_xxxyzz[k];

                g_0_y_xyyy_xxzzz[k] = -g_0_y_yyy_xxzzz[k] * ab_x + g_0_y_yyy_xxxzzz[k];

                g_0_y_xyyy_xyyyy[k] = -g_0_y_yyy_xyyyy[k] * ab_x + g_0_y_yyy_xxyyyy[k];

                g_0_y_xyyy_xyyyz[k] = -g_0_y_yyy_xyyyz[k] * ab_x + g_0_y_yyy_xxyyyz[k];

                g_0_y_xyyy_xyyzz[k] = -g_0_y_yyy_xyyzz[k] * ab_x + g_0_y_yyy_xxyyzz[k];

                g_0_y_xyyy_xyzzz[k] = -g_0_y_yyy_xyzzz[k] * ab_x + g_0_y_yyy_xxyzzz[k];

                g_0_y_xyyy_xzzzz[k] = -g_0_y_yyy_xzzzz[k] * ab_x + g_0_y_yyy_xxzzzz[k];

                g_0_y_xyyy_yyyyy[k] = -g_0_y_yyy_yyyyy[k] * ab_x + g_0_y_yyy_xyyyyy[k];

                g_0_y_xyyy_yyyyz[k] = -g_0_y_yyy_yyyyz[k] * ab_x + g_0_y_yyy_xyyyyz[k];

                g_0_y_xyyy_yyyzz[k] = -g_0_y_yyy_yyyzz[k] * ab_x + g_0_y_yyy_xyyyzz[k];

                g_0_y_xyyy_yyzzz[k] = -g_0_y_yyy_yyzzz[k] * ab_x + g_0_y_yyy_xyyzzz[k];

                g_0_y_xyyy_yzzzz[k] = -g_0_y_yyy_yzzzz[k] * ab_x + g_0_y_yyy_xyzzzz[k];

                g_0_y_xyyy_zzzzz[k] = -g_0_y_yyy_zzzzz[k] * ab_x + g_0_y_yyy_xzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyyz_xxxxx, g_0_y_xyyz_xxxxy, g_0_y_xyyz_xxxxz, g_0_y_xyyz_xxxyy, g_0_y_xyyz_xxxyz, g_0_y_xyyz_xxxzz, g_0_y_xyyz_xxyyy, g_0_y_xyyz_xxyyz, g_0_y_xyyz_xxyzz, g_0_y_xyyz_xxzzz, g_0_y_xyyz_xyyyy, g_0_y_xyyz_xyyyz, g_0_y_xyyz_xyyzz, g_0_y_xyyz_xyzzz, g_0_y_xyyz_xzzzz, g_0_y_xyyz_yyyyy, g_0_y_xyyz_yyyyz, g_0_y_xyyz_yyyzz, g_0_y_xyyz_yyzzz, g_0_y_xyyz_yzzzz, g_0_y_xyyz_zzzzz, g_0_y_yyz_xxxxx, g_0_y_yyz_xxxxxx, g_0_y_yyz_xxxxxy, g_0_y_yyz_xxxxxz, g_0_y_yyz_xxxxy, g_0_y_yyz_xxxxyy, g_0_y_yyz_xxxxyz, g_0_y_yyz_xxxxz, g_0_y_yyz_xxxxzz, g_0_y_yyz_xxxyy, g_0_y_yyz_xxxyyy, g_0_y_yyz_xxxyyz, g_0_y_yyz_xxxyz, g_0_y_yyz_xxxyzz, g_0_y_yyz_xxxzz, g_0_y_yyz_xxxzzz, g_0_y_yyz_xxyyy, g_0_y_yyz_xxyyyy, g_0_y_yyz_xxyyyz, g_0_y_yyz_xxyyz, g_0_y_yyz_xxyyzz, g_0_y_yyz_xxyzz, g_0_y_yyz_xxyzzz, g_0_y_yyz_xxzzz, g_0_y_yyz_xxzzzz, g_0_y_yyz_xyyyy, g_0_y_yyz_xyyyyy, g_0_y_yyz_xyyyyz, g_0_y_yyz_xyyyz, g_0_y_yyz_xyyyzz, g_0_y_yyz_xyyzz, g_0_y_yyz_xyyzzz, g_0_y_yyz_xyzzz, g_0_y_yyz_xyzzzz, g_0_y_yyz_xzzzz, g_0_y_yyz_xzzzzz, g_0_y_yyz_yyyyy, g_0_y_yyz_yyyyz, g_0_y_yyz_yyyzz, g_0_y_yyz_yyzzz, g_0_y_yyz_yzzzz, g_0_y_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_xxxxx[k] = -g_0_y_yyz_xxxxx[k] * ab_x + g_0_y_yyz_xxxxxx[k];

                g_0_y_xyyz_xxxxy[k] = -g_0_y_yyz_xxxxy[k] * ab_x + g_0_y_yyz_xxxxxy[k];

                g_0_y_xyyz_xxxxz[k] = -g_0_y_yyz_xxxxz[k] * ab_x + g_0_y_yyz_xxxxxz[k];

                g_0_y_xyyz_xxxyy[k] = -g_0_y_yyz_xxxyy[k] * ab_x + g_0_y_yyz_xxxxyy[k];

                g_0_y_xyyz_xxxyz[k] = -g_0_y_yyz_xxxyz[k] * ab_x + g_0_y_yyz_xxxxyz[k];

                g_0_y_xyyz_xxxzz[k] = -g_0_y_yyz_xxxzz[k] * ab_x + g_0_y_yyz_xxxxzz[k];

                g_0_y_xyyz_xxyyy[k] = -g_0_y_yyz_xxyyy[k] * ab_x + g_0_y_yyz_xxxyyy[k];

                g_0_y_xyyz_xxyyz[k] = -g_0_y_yyz_xxyyz[k] * ab_x + g_0_y_yyz_xxxyyz[k];

                g_0_y_xyyz_xxyzz[k] = -g_0_y_yyz_xxyzz[k] * ab_x + g_0_y_yyz_xxxyzz[k];

                g_0_y_xyyz_xxzzz[k] = -g_0_y_yyz_xxzzz[k] * ab_x + g_0_y_yyz_xxxzzz[k];

                g_0_y_xyyz_xyyyy[k] = -g_0_y_yyz_xyyyy[k] * ab_x + g_0_y_yyz_xxyyyy[k];

                g_0_y_xyyz_xyyyz[k] = -g_0_y_yyz_xyyyz[k] * ab_x + g_0_y_yyz_xxyyyz[k];

                g_0_y_xyyz_xyyzz[k] = -g_0_y_yyz_xyyzz[k] * ab_x + g_0_y_yyz_xxyyzz[k];

                g_0_y_xyyz_xyzzz[k] = -g_0_y_yyz_xyzzz[k] * ab_x + g_0_y_yyz_xxyzzz[k];

                g_0_y_xyyz_xzzzz[k] = -g_0_y_yyz_xzzzz[k] * ab_x + g_0_y_yyz_xxzzzz[k];

                g_0_y_xyyz_yyyyy[k] = -g_0_y_yyz_yyyyy[k] * ab_x + g_0_y_yyz_xyyyyy[k];

                g_0_y_xyyz_yyyyz[k] = -g_0_y_yyz_yyyyz[k] * ab_x + g_0_y_yyz_xyyyyz[k];

                g_0_y_xyyz_yyyzz[k] = -g_0_y_yyz_yyyzz[k] * ab_x + g_0_y_yyz_xyyyzz[k];

                g_0_y_xyyz_yyzzz[k] = -g_0_y_yyz_yyzzz[k] * ab_x + g_0_y_yyz_xyyzzz[k];

                g_0_y_xyyz_yzzzz[k] = -g_0_y_yyz_yzzzz[k] * ab_x + g_0_y_yyz_xyzzzz[k];

                g_0_y_xyyz_zzzzz[k] = -g_0_y_yyz_zzzzz[k] * ab_x + g_0_y_yyz_xzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyzz_xxxxx, g_0_y_xyzz_xxxxy, g_0_y_xyzz_xxxxz, g_0_y_xyzz_xxxyy, g_0_y_xyzz_xxxyz, g_0_y_xyzz_xxxzz, g_0_y_xyzz_xxyyy, g_0_y_xyzz_xxyyz, g_0_y_xyzz_xxyzz, g_0_y_xyzz_xxzzz, g_0_y_xyzz_xyyyy, g_0_y_xyzz_xyyyz, g_0_y_xyzz_xyyzz, g_0_y_xyzz_xyzzz, g_0_y_xyzz_xzzzz, g_0_y_xyzz_yyyyy, g_0_y_xyzz_yyyyz, g_0_y_xyzz_yyyzz, g_0_y_xyzz_yyzzz, g_0_y_xyzz_yzzzz, g_0_y_xyzz_zzzzz, g_0_y_yzz_xxxxx, g_0_y_yzz_xxxxxx, g_0_y_yzz_xxxxxy, g_0_y_yzz_xxxxxz, g_0_y_yzz_xxxxy, g_0_y_yzz_xxxxyy, g_0_y_yzz_xxxxyz, g_0_y_yzz_xxxxz, g_0_y_yzz_xxxxzz, g_0_y_yzz_xxxyy, g_0_y_yzz_xxxyyy, g_0_y_yzz_xxxyyz, g_0_y_yzz_xxxyz, g_0_y_yzz_xxxyzz, g_0_y_yzz_xxxzz, g_0_y_yzz_xxxzzz, g_0_y_yzz_xxyyy, g_0_y_yzz_xxyyyy, g_0_y_yzz_xxyyyz, g_0_y_yzz_xxyyz, g_0_y_yzz_xxyyzz, g_0_y_yzz_xxyzz, g_0_y_yzz_xxyzzz, g_0_y_yzz_xxzzz, g_0_y_yzz_xxzzzz, g_0_y_yzz_xyyyy, g_0_y_yzz_xyyyyy, g_0_y_yzz_xyyyyz, g_0_y_yzz_xyyyz, g_0_y_yzz_xyyyzz, g_0_y_yzz_xyyzz, g_0_y_yzz_xyyzzz, g_0_y_yzz_xyzzz, g_0_y_yzz_xyzzzz, g_0_y_yzz_xzzzz, g_0_y_yzz_xzzzzz, g_0_y_yzz_yyyyy, g_0_y_yzz_yyyyz, g_0_y_yzz_yyyzz, g_0_y_yzz_yyzzz, g_0_y_yzz_yzzzz, g_0_y_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_xxxxx[k] = -g_0_y_yzz_xxxxx[k] * ab_x + g_0_y_yzz_xxxxxx[k];

                g_0_y_xyzz_xxxxy[k] = -g_0_y_yzz_xxxxy[k] * ab_x + g_0_y_yzz_xxxxxy[k];

                g_0_y_xyzz_xxxxz[k] = -g_0_y_yzz_xxxxz[k] * ab_x + g_0_y_yzz_xxxxxz[k];

                g_0_y_xyzz_xxxyy[k] = -g_0_y_yzz_xxxyy[k] * ab_x + g_0_y_yzz_xxxxyy[k];

                g_0_y_xyzz_xxxyz[k] = -g_0_y_yzz_xxxyz[k] * ab_x + g_0_y_yzz_xxxxyz[k];

                g_0_y_xyzz_xxxzz[k] = -g_0_y_yzz_xxxzz[k] * ab_x + g_0_y_yzz_xxxxzz[k];

                g_0_y_xyzz_xxyyy[k] = -g_0_y_yzz_xxyyy[k] * ab_x + g_0_y_yzz_xxxyyy[k];

                g_0_y_xyzz_xxyyz[k] = -g_0_y_yzz_xxyyz[k] * ab_x + g_0_y_yzz_xxxyyz[k];

                g_0_y_xyzz_xxyzz[k] = -g_0_y_yzz_xxyzz[k] * ab_x + g_0_y_yzz_xxxyzz[k];

                g_0_y_xyzz_xxzzz[k] = -g_0_y_yzz_xxzzz[k] * ab_x + g_0_y_yzz_xxxzzz[k];

                g_0_y_xyzz_xyyyy[k] = -g_0_y_yzz_xyyyy[k] * ab_x + g_0_y_yzz_xxyyyy[k];

                g_0_y_xyzz_xyyyz[k] = -g_0_y_yzz_xyyyz[k] * ab_x + g_0_y_yzz_xxyyyz[k];

                g_0_y_xyzz_xyyzz[k] = -g_0_y_yzz_xyyzz[k] * ab_x + g_0_y_yzz_xxyyzz[k];

                g_0_y_xyzz_xyzzz[k] = -g_0_y_yzz_xyzzz[k] * ab_x + g_0_y_yzz_xxyzzz[k];

                g_0_y_xyzz_xzzzz[k] = -g_0_y_yzz_xzzzz[k] * ab_x + g_0_y_yzz_xxzzzz[k];

                g_0_y_xyzz_yyyyy[k] = -g_0_y_yzz_yyyyy[k] * ab_x + g_0_y_yzz_xyyyyy[k];

                g_0_y_xyzz_yyyyz[k] = -g_0_y_yzz_yyyyz[k] * ab_x + g_0_y_yzz_xyyyyz[k];

                g_0_y_xyzz_yyyzz[k] = -g_0_y_yzz_yyyzz[k] * ab_x + g_0_y_yzz_xyyyzz[k];

                g_0_y_xyzz_yyzzz[k] = -g_0_y_yzz_yyzzz[k] * ab_x + g_0_y_yzz_xyyzzz[k];

                g_0_y_xyzz_yzzzz[k] = -g_0_y_yzz_yzzzz[k] * ab_x + g_0_y_yzz_xyzzzz[k];

                g_0_y_xyzz_zzzzz[k] = -g_0_y_yzz_zzzzz[k] * ab_x + g_0_y_yzz_xzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xzzz_xxxxx, g_0_y_xzzz_xxxxy, g_0_y_xzzz_xxxxz, g_0_y_xzzz_xxxyy, g_0_y_xzzz_xxxyz, g_0_y_xzzz_xxxzz, g_0_y_xzzz_xxyyy, g_0_y_xzzz_xxyyz, g_0_y_xzzz_xxyzz, g_0_y_xzzz_xxzzz, g_0_y_xzzz_xyyyy, g_0_y_xzzz_xyyyz, g_0_y_xzzz_xyyzz, g_0_y_xzzz_xyzzz, g_0_y_xzzz_xzzzz, g_0_y_xzzz_yyyyy, g_0_y_xzzz_yyyyz, g_0_y_xzzz_yyyzz, g_0_y_xzzz_yyzzz, g_0_y_xzzz_yzzzz, g_0_y_xzzz_zzzzz, g_0_y_zzz_xxxxx, g_0_y_zzz_xxxxxx, g_0_y_zzz_xxxxxy, g_0_y_zzz_xxxxxz, g_0_y_zzz_xxxxy, g_0_y_zzz_xxxxyy, g_0_y_zzz_xxxxyz, g_0_y_zzz_xxxxz, g_0_y_zzz_xxxxzz, g_0_y_zzz_xxxyy, g_0_y_zzz_xxxyyy, g_0_y_zzz_xxxyyz, g_0_y_zzz_xxxyz, g_0_y_zzz_xxxyzz, g_0_y_zzz_xxxzz, g_0_y_zzz_xxxzzz, g_0_y_zzz_xxyyy, g_0_y_zzz_xxyyyy, g_0_y_zzz_xxyyyz, g_0_y_zzz_xxyyz, g_0_y_zzz_xxyyzz, g_0_y_zzz_xxyzz, g_0_y_zzz_xxyzzz, g_0_y_zzz_xxzzz, g_0_y_zzz_xxzzzz, g_0_y_zzz_xyyyy, g_0_y_zzz_xyyyyy, g_0_y_zzz_xyyyyz, g_0_y_zzz_xyyyz, g_0_y_zzz_xyyyzz, g_0_y_zzz_xyyzz, g_0_y_zzz_xyyzzz, g_0_y_zzz_xyzzz, g_0_y_zzz_xyzzzz, g_0_y_zzz_xzzzz, g_0_y_zzz_xzzzzz, g_0_y_zzz_yyyyy, g_0_y_zzz_yyyyz, g_0_y_zzz_yyyzz, g_0_y_zzz_yyzzz, g_0_y_zzz_yzzzz, g_0_y_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_xxxxx[k] = -g_0_y_zzz_xxxxx[k] * ab_x + g_0_y_zzz_xxxxxx[k];

                g_0_y_xzzz_xxxxy[k] = -g_0_y_zzz_xxxxy[k] * ab_x + g_0_y_zzz_xxxxxy[k];

                g_0_y_xzzz_xxxxz[k] = -g_0_y_zzz_xxxxz[k] * ab_x + g_0_y_zzz_xxxxxz[k];

                g_0_y_xzzz_xxxyy[k] = -g_0_y_zzz_xxxyy[k] * ab_x + g_0_y_zzz_xxxxyy[k];

                g_0_y_xzzz_xxxyz[k] = -g_0_y_zzz_xxxyz[k] * ab_x + g_0_y_zzz_xxxxyz[k];

                g_0_y_xzzz_xxxzz[k] = -g_0_y_zzz_xxxzz[k] * ab_x + g_0_y_zzz_xxxxzz[k];

                g_0_y_xzzz_xxyyy[k] = -g_0_y_zzz_xxyyy[k] * ab_x + g_0_y_zzz_xxxyyy[k];

                g_0_y_xzzz_xxyyz[k] = -g_0_y_zzz_xxyyz[k] * ab_x + g_0_y_zzz_xxxyyz[k];

                g_0_y_xzzz_xxyzz[k] = -g_0_y_zzz_xxyzz[k] * ab_x + g_0_y_zzz_xxxyzz[k];

                g_0_y_xzzz_xxzzz[k] = -g_0_y_zzz_xxzzz[k] * ab_x + g_0_y_zzz_xxxzzz[k];

                g_0_y_xzzz_xyyyy[k] = -g_0_y_zzz_xyyyy[k] * ab_x + g_0_y_zzz_xxyyyy[k];

                g_0_y_xzzz_xyyyz[k] = -g_0_y_zzz_xyyyz[k] * ab_x + g_0_y_zzz_xxyyyz[k];

                g_0_y_xzzz_xyyzz[k] = -g_0_y_zzz_xyyzz[k] * ab_x + g_0_y_zzz_xxyyzz[k];

                g_0_y_xzzz_xyzzz[k] = -g_0_y_zzz_xyzzz[k] * ab_x + g_0_y_zzz_xxyzzz[k];

                g_0_y_xzzz_xzzzz[k] = -g_0_y_zzz_xzzzz[k] * ab_x + g_0_y_zzz_xxzzzz[k];

                g_0_y_xzzz_yyyyy[k] = -g_0_y_zzz_yyyyy[k] * ab_x + g_0_y_zzz_xyyyyy[k];

                g_0_y_xzzz_yyyyz[k] = -g_0_y_zzz_yyyyz[k] * ab_x + g_0_y_zzz_xyyyyz[k];

                g_0_y_xzzz_yyyzz[k] = -g_0_y_zzz_yyyzz[k] * ab_x + g_0_y_zzz_xyyyzz[k];

                g_0_y_xzzz_yyzzz[k] = -g_0_y_zzz_yyzzz[k] * ab_x + g_0_y_zzz_xyyzzz[k];

                g_0_y_xzzz_yzzzz[k] = -g_0_y_zzz_yzzzz[k] * ab_x + g_0_y_zzz_xyzzzz[k];

                g_0_y_xzzz_zzzzz[k] = -g_0_y_zzz_zzzzz[k] * ab_x + g_0_y_zzz_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyy_xxxxx, g_0_y_yyy_xxxxxy, g_0_y_yyy_xxxxy, g_0_y_yyy_xxxxyy, g_0_y_yyy_xxxxyz, g_0_y_yyy_xxxxz, g_0_y_yyy_xxxyy, g_0_y_yyy_xxxyyy, g_0_y_yyy_xxxyyz, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxyzz, g_0_y_yyy_xxxzz, g_0_y_yyy_xxyyy, g_0_y_yyy_xxyyyy, g_0_y_yyy_xxyyyz, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyyzz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxyzzz, g_0_y_yyy_xxzzz, g_0_y_yyy_xyyyy, g_0_y_yyy_xyyyyy, g_0_y_yyy_xyyyyz, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyyzz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyyzzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xyzzzz, g_0_y_yyy_xzzzz, g_0_y_yyy_yyyyy, g_0_y_yyy_yyyyyy, g_0_y_yyy_yyyyyz, g_0_y_yyy_yyyyz, g_0_y_yyy_yyyyzz, g_0_y_yyy_yyyzz, g_0_y_yyy_yyyzzz, g_0_y_yyy_yyzzz, g_0_y_yyy_yyzzzz, g_0_y_yyy_yzzzz, g_0_y_yyy_yzzzzz, g_0_y_yyy_zzzzz, g_0_y_yyyy_xxxxx, g_0_y_yyyy_xxxxy, g_0_y_yyyy_xxxxz, g_0_y_yyyy_xxxyy, g_0_y_yyyy_xxxyz, g_0_y_yyyy_xxxzz, g_0_y_yyyy_xxyyy, g_0_y_yyyy_xxyyz, g_0_y_yyyy_xxyzz, g_0_y_yyyy_xxzzz, g_0_y_yyyy_xyyyy, g_0_y_yyyy_xyyyz, g_0_y_yyyy_xyyzz, g_0_y_yyyy_xyzzz, g_0_y_yyyy_xzzzz, g_0_y_yyyy_yyyyy, g_0_y_yyyy_yyyyz, g_0_y_yyyy_yyyzz, g_0_y_yyyy_yyzzz, g_0_y_yyyy_yzzzz, g_0_y_yyyy_zzzzz, g_yyy_xxxxx, g_yyy_xxxxy, g_yyy_xxxxz, g_yyy_xxxyy, g_yyy_xxxyz, g_yyy_xxxzz, g_yyy_xxyyy, g_yyy_xxyyz, g_yyy_xxyzz, g_yyy_xxzzz, g_yyy_xyyyy, g_yyy_xyyyz, g_yyy_xyyzz, g_yyy_xyzzz, g_yyy_xzzzz, g_yyy_yyyyy, g_yyy_yyyyz, g_yyy_yyyzz, g_yyy_yyzzz, g_yyy_yzzzz, g_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_xxxxx[k] = g_yyy_xxxxx[k] - g_0_y_yyy_xxxxx[k] * ab_y + g_0_y_yyy_xxxxxy[k];

                g_0_y_yyyy_xxxxy[k] = g_yyy_xxxxy[k] - g_0_y_yyy_xxxxy[k] * ab_y + g_0_y_yyy_xxxxyy[k];

                g_0_y_yyyy_xxxxz[k] = g_yyy_xxxxz[k] - g_0_y_yyy_xxxxz[k] * ab_y + g_0_y_yyy_xxxxyz[k];

                g_0_y_yyyy_xxxyy[k] = g_yyy_xxxyy[k] - g_0_y_yyy_xxxyy[k] * ab_y + g_0_y_yyy_xxxyyy[k];

                g_0_y_yyyy_xxxyz[k] = g_yyy_xxxyz[k] - g_0_y_yyy_xxxyz[k] * ab_y + g_0_y_yyy_xxxyyz[k];

                g_0_y_yyyy_xxxzz[k] = g_yyy_xxxzz[k] - g_0_y_yyy_xxxzz[k] * ab_y + g_0_y_yyy_xxxyzz[k];

                g_0_y_yyyy_xxyyy[k] = g_yyy_xxyyy[k] - g_0_y_yyy_xxyyy[k] * ab_y + g_0_y_yyy_xxyyyy[k];

                g_0_y_yyyy_xxyyz[k] = g_yyy_xxyyz[k] - g_0_y_yyy_xxyyz[k] * ab_y + g_0_y_yyy_xxyyyz[k];

                g_0_y_yyyy_xxyzz[k] = g_yyy_xxyzz[k] - g_0_y_yyy_xxyzz[k] * ab_y + g_0_y_yyy_xxyyzz[k];

                g_0_y_yyyy_xxzzz[k] = g_yyy_xxzzz[k] - g_0_y_yyy_xxzzz[k] * ab_y + g_0_y_yyy_xxyzzz[k];

                g_0_y_yyyy_xyyyy[k] = g_yyy_xyyyy[k] - g_0_y_yyy_xyyyy[k] * ab_y + g_0_y_yyy_xyyyyy[k];

                g_0_y_yyyy_xyyyz[k] = g_yyy_xyyyz[k] - g_0_y_yyy_xyyyz[k] * ab_y + g_0_y_yyy_xyyyyz[k];

                g_0_y_yyyy_xyyzz[k] = g_yyy_xyyzz[k] - g_0_y_yyy_xyyzz[k] * ab_y + g_0_y_yyy_xyyyzz[k];

                g_0_y_yyyy_xyzzz[k] = g_yyy_xyzzz[k] - g_0_y_yyy_xyzzz[k] * ab_y + g_0_y_yyy_xyyzzz[k];

                g_0_y_yyyy_xzzzz[k] = g_yyy_xzzzz[k] - g_0_y_yyy_xzzzz[k] * ab_y + g_0_y_yyy_xyzzzz[k];

                g_0_y_yyyy_yyyyy[k] = g_yyy_yyyyy[k] - g_0_y_yyy_yyyyy[k] * ab_y + g_0_y_yyy_yyyyyy[k];

                g_0_y_yyyy_yyyyz[k] = g_yyy_yyyyz[k] - g_0_y_yyy_yyyyz[k] * ab_y + g_0_y_yyy_yyyyyz[k];

                g_0_y_yyyy_yyyzz[k] = g_yyy_yyyzz[k] - g_0_y_yyy_yyyzz[k] * ab_y + g_0_y_yyy_yyyyzz[k];

                g_0_y_yyyy_yyzzz[k] = g_yyy_yyzzz[k] - g_0_y_yyy_yyzzz[k] * ab_y + g_0_y_yyy_yyyzzz[k];

                g_0_y_yyyy_yzzzz[k] = g_yyy_yzzzz[k] - g_0_y_yyy_yzzzz[k] * ab_y + g_0_y_yyy_yyzzzz[k];

                g_0_y_yyyy_zzzzz[k] = g_yyy_zzzzz[k] - g_0_y_yyy_zzzzz[k] * ab_y + g_0_y_yyy_yzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyy_xxxxx, g_0_y_yyy_xxxxxz, g_0_y_yyy_xxxxy, g_0_y_yyy_xxxxyz, g_0_y_yyy_xxxxz, g_0_y_yyy_xxxxzz, g_0_y_yyy_xxxyy, g_0_y_yyy_xxxyyz, g_0_y_yyy_xxxyz, g_0_y_yyy_xxxyzz, g_0_y_yyy_xxxzz, g_0_y_yyy_xxxzzz, g_0_y_yyy_xxyyy, g_0_y_yyy_xxyyyz, g_0_y_yyy_xxyyz, g_0_y_yyy_xxyyzz, g_0_y_yyy_xxyzz, g_0_y_yyy_xxyzzz, g_0_y_yyy_xxzzz, g_0_y_yyy_xxzzzz, g_0_y_yyy_xyyyy, g_0_y_yyy_xyyyyz, g_0_y_yyy_xyyyz, g_0_y_yyy_xyyyzz, g_0_y_yyy_xyyzz, g_0_y_yyy_xyyzzz, g_0_y_yyy_xyzzz, g_0_y_yyy_xyzzzz, g_0_y_yyy_xzzzz, g_0_y_yyy_xzzzzz, g_0_y_yyy_yyyyy, g_0_y_yyy_yyyyyz, g_0_y_yyy_yyyyz, g_0_y_yyy_yyyyzz, g_0_y_yyy_yyyzz, g_0_y_yyy_yyyzzz, g_0_y_yyy_yyzzz, g_0_y_yyy_yyzzzz, g_0_y_yyy_yzzzz, g_0_y_yyy_yzzzzz, g_0_y_yyy_zzzzz, g_0_y_yyy_zzzzzz, g_0_y_yyyz_xxxxx, g_0_y_yyyz_xxxxy, g_0_y_yyyz_xxxxz, g_0_y_yyyz_xxxyy, g_0_y_yyyz_xxxyz, g_0_y_yyyz_xxxzz, g_0_y_yyyz_xxyyy, g_0_y_yyyz_xxyyz, g_0_y_yyyz_xxyzz, g_0_y_yyyz_xxzzz, g_0_y_yyyz_xyyyy, g_0_y_yyyz_xyyyz, g_0_y_yyyz_xyyzz, g_0_y_yyyz_xyzzz, g_0_y_yyyz_xzzzz, g_0_y_yyyz_yyyyy, g_0_y_yyyz_yyyyz, g_0_y_yyyz_yyyzz, g_0_y_yyyz_yyzzz, g_0_y_yyyz_yzzzz, g_0_y_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_xxxxx[k] = -g_0_y_yyy_xxxxx[k] * ab_z + g_0_y_yyy_xxxxxz[k];

                g_0_y_yyyz_xxxxy[k] = -g_0_y_yyy_xxxxy[k] * ab_z + g_0_y_yyy_xxxxyz[k];

                g_0_y_yyyz_xxxxz[k] = -g_0_y_yyy_xxxxz[k] * ab_z + g_0_y_yyy_xxxxzz[k];

                g_0_y_yyyz_xxxyy[k] = -g_0_y_yyy_xxxyy[k] * ab_z + g_0_y_yyy_xxxyyz[k];

                g_0_y_yyyz_xxxyz[k] = -g_0_y_yyy_xxxyz[k] * ab_z + g_0_y_yyy_xxxyzz[k];

                g_0_y_yyyz_xxxzz[k] = -g_0_y_yyy_xxxzz[k] * ab_z + g_0_y_yyy_xxxzzz[k];

                g_0_y_yyyz_xxyyy[k] = -g_0_y_yyy_xxyyy[k] * ab_z + g_0_y_yyy_xxyyyz[k];

                g_0_y_yyyz_xxyyz[k] = -g_0_y_yyy_xxyyz[k] * ab_z + g_0_y_yyy_xxyyzz[k];

                g_0_y_yyyz_xxyzz[k] = -g_0_y_yyy_xxyzz[k] * ab_z + g_0_y_yyy_xxyzzz[k];

                g_0_y_yyyz_xxzzz[k] = -g_0_y_yyy_xxzzz[k] * ab_z + g_0_y_yyy_xxzzzz[k];

                g_0_y_yyyz_xyyyy[k] = -g_0_y_yyy_xyyyy[k] * ab_z + g_0_y_yyy_xyyyyz[k];

                g_0_y_yyyz_xyyyz[k] = -g_0_y_yyy_xyyyz[k] * ab_z + g_0_y_yyy_xyyyzz[k];

                g_0_y_yyyz_xyyzz[k] = -g_0_y_yyy_xyyzz[k] * ab_z + g_0_y_yyy_xyyzzz[k];

                g_0_y_yyyz_xyzzz[k] = -g_0_y_yyy_xyzzz[k] * ab_z + g_0_y_yyy_xyzzzz[k];

                g_0_y_yyyz_xzzzz[k] = -g_0_y_yyy_xzzzz[k] * ab_z + g_0_y_yyy_xzzzzz[k];

                g_0_y_yyyz_yyyyy[k] = -g_0_y_yyy_yyyyy[k] * ab_z + g_0_y_yyy_yyyyyz[k];

                g_0_y_yyyz_yyyyz[k] = -g_0_y_yyy_yyyyz[k] * ab_z + g_0_y_yyy_yyyyzz[k];

                g_0_y_yyyz_yyyzz[k] = -g_0_y_yyy_yyyzz[k] * ab_z + g_0_y_yyy_yyyzzz[k];

                g_0_y_yyyz_yyzzz[k] = -g_0_y_yyy_yyzzz[k] * ab_z + g_0_y_yyy_yyzzzz[k];

                g_0_y_yyyz_yzzzz[k] = -g_0_y_yyy_yzzzz[k] * ab_z + g_0_y_yyy_yzzzzz[k];

                g_0_y_yyyz_zzzzz[k] = -g_0_y_yyy_zzzzz[k] * ab_z + g_0_y_yyy_zzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yyz_xxxxx, g_0_y_yyz_xxxxxz, g_0_y_yyz_xxxxy, g_0_y_yyz_xxxxyz, g_0_y_yyz_xxxxz, g_0_y_yyz_xxxxzz, g_0_y_yyz_xxxyy, g_0_y_yyz_xxxyyz, g_0_y_yyz_xxxyz, g_0_y_yyz_xxxyzz, g_0_y_yyz_xxxzz, g_0_y_yyz_xxxzzz, g_0_y_yyz_xxyyy, g_0_y_yyz_xxyyyz, g_0_y_yyz_xxyyz, g_0_y_yyz_xxyyzz, g_0_y_yyz_xxyzz, g_0_y_yyz_xxyzzz, g_0_y_yyz_xxzzz, g_0_y_yyz_xxzzzz, g_0_y_yyz_xyyyy, g_0_y_yyz_xyyyyz, g_0_y_yyz_xyyyz, g_0_y_yyz_xyyyzz, g_0_y_yyz_xyyzz, g_0_y_yyz_xyyzzz, g_0_y_yyz_xyzzz, g_0_y_yyz_xyzzzz, g_0_y_yyz_xzzzz, g_0_y_yyz_xzzzzz, g_0_y_yyz_yyyyy, g_0_y_yyz_yyyyyz, g_0_y_yyz_yyyyz, g_0_y_yyz_yyyyzz, g_0_y_yyz_yyyzz, g_0_y_yyz_yyyzzz, g_0_y_yyz_yyzzz, g_0_y_yyz_yyzzzz, g_0_y_yyz_yzzzz, g_0_y_yyz_yzzzzz, g_0_y_yyz_zzzzz, g_0_y_yyz_zzzzzz, g_0_y_yyzz_xxxxx, g_0_y_yyzz_xxxxy, g_0_y_yyzz_xxxxz, g_0_y_yyzz_xxxyy, g_0_y_yyzz_xxxyz, g_0_y_yyzz_xxxzz, g_0_y_yyzz_xxyyy, g_0_y_yyzz_xxyyz, g_0_y_yyzz_xxyzz, g_0_y_yyzz_xxzzz, g_0_y_yyzz_xyyyy, g_0_y_yyzz_xyyyz, g_0_y_yyzz_xyyzz, g_0_y_yyzz_xyzzz, g_0_y_yyzz_xzzzz, g_0_y_yyzz_yyyyy, g_0_y_yyzz_yyyyz, g_0_y_yyzz_yyyzz, g_0_y_yyzz_yyzzz, g_0_y_yyzz_yzzzz, g_0_y_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_xxxxx[k] = -g_0_y_yyz_xxxxx[k] * ab_z + g_0_y_yyz_xxxxxz[k];

                g_0_y_yyzz_xxxxy[k] = -g_0_y_yyz_xxxxy[k] * ab_z + g_0_y_yyz_xxxxyz[k];

                g_0_y_yyzz_xxxxz[k] = -g_0_y_yyz_xxxxz[k] * ab_z + g_0_y_yyz_xxxxzz[k];

                g_0_y_yyzz_xxxyy[k] = -g_0_y_yyz_xxxyy[k] * ab_z + g_0_y_yyz_xxxyyz[k];

                g_0_y_yyzz_xxxyz[k] = -g_0_y_yyz_xxxyz[k] * ab_z + g_0_y_yyz_xxxyzz[k];

                g_0_y_yyzz_xxxzz[k] = -g_0_y_yyz_xxxzz[k] * ab_z + g_0_y_yyz_xxxzzz[k];

                g_0_y_yyzz_xxyyy[k] = -g_0_y_yyz_xxyyy[k] * ab_z + g_0_y_yyz_xxyyyz[k];

                g_0_y_yyzz_xxyyz[k] = -g_0_y_yyz_xxyyz[k] * ab_z + g_0_y_yyz_xxyyzz[k];

                g_0_y_yyzz_xxyzz[k] = -g_0_y_yyz_xxyzz[k] * ab_z + g_0_y_yyz_xxyzzz[k];

                g_0_y_yyzz_xxzzz[k] = -g_0_y_yyz_xxzzz[k] * ab_z + g_0_y_yyz_xxzzzz[k];

                g_0_y_yyzz_xyyyy[k] = -g_0_y_yyz_xyyyy[k] * ab_z + g_0_y_yyz_xyyyyz[k];

                g_0_y_yyzz_xyyyz[k] = -g_0_y_yyz_xyyyz[k] * ab_z + g_0_y_yyz_xyyyzz[k];

                g_0_y_yyzz_xyyzz[k] = -g_0_y_yyz_xyyzz[k] * ab_z + g_0_y_yyz_xyyzzz[k];

                g_0_y_yyzz_xyzzz[k] = -g_0_y_yyz_xyzzz[k] * ab_z + g_0_y_yyz_xyzzzz[k];

                g_0_y_yyzz_xzzzz[k] = -g_0_y_yyz_xzzzz[k] * ab_z + g_0_y_yyz_xzzzzz[k];

                g_0_y_yyzz_yyyyy[k] = -g_0_y_yyz_yyyyy[k] * ab_z + g_0_y_yyz_yyyyyz[k];

                g_0_y_yyzz_yyyyz[k] = -g_0_y_yyz_yyyyz[k] * ab_z + g_0_y_yyz_yyyyzz[k];

                g_0_y_yyzz_yyyzz[k] = -g_0_y_yyz_yyyzz[k] * ab_z + g_0_y_yyz_yyyzzz[k];

                g_0_y_yyzz_yyzzz[k] = -g_0_y_yyz_yyzzz[k] * ab_z + g_0_y_yyz_yyzzzz[k];

                g_0_y_yyzz_yzzzz[k] = -g_0_y_yyz_yzzzz[k] * ab_z + g_0_y_yyz_yzzzzz[k];

                g_0_y_yyzz_zzzzz[k] = -g_0_y_yyz_zzzzz[k] * ab_z + g_0_y_yyz_zzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yzz_xxxxx, g_0_y_yzz_xxxxxz, g_0_y_yzz_xxxxy, g_0_y_yzz_xxxxyz, g_0_y_yzz_xxxxz, g_0_y_yzz_xxxxzz, g_0_y_yzz_xxxyy, g_0_y_yzz_xxxyyz, g_0_y_yzz_xxxyz, g_0_y_yzz_xxxyzz, g_0_y_yzz_xxxzz, g_0_y_yzz_xxxzzz, g_0_y_yzz_xxyyy, g_0_y_yzz_xxyyyz, g_0_y_yzz_xxyyz, g_0_y_yzz_xxyyzz, g_0_y_yzz_xxyzz, g_0_y_yzz_xxyzzz, g_0_y_yzz_xxzzz, g_0_y_yzz_xxzzzz, g_0_y_yzz_xyyyy, g_0_y_yzz_xyyyyz, g_0_y_yzz_xyyyz, g_0_y_yzz_xyyyzz, g_0_y_yzz_xyyzz, g_0_y_yzz_xyyzzz, g_0_y_yzz_xyzzz, g_0_y_yzz_xyzzzz, g_0_y_yzz_xzzzz, g_0_y_yzz_xzzzzz, g_0_y_yzz_yyyyy, g_0_y_yzz_yyyyyz, g_0_y_yzz_yyyyz, g_0_y_yzz_yyyyzz, g_0_y_yzz_yyyzz, g_0_y_yzz_yyyzzz, g_0_y_yzz_yyzzz, g_0_y_yzz_yyzzzz, g_0_y_yzz_yzzzz, g_0_y_yzz_yzzzzz, g_0_y_yzz_zzzzz, g_0_y_yzz_zzzzzz, g_0_y_yzzz_xxxxx, g_0_y_yzzz_xxxxy, g_0_y_yzzz_xxxxz, g_0_y_yzzz_xxxyy, g_0_y_yzzz_xxxyz, g_0_y_yzzz_xxxzz, g_0_y_yzzz_xxyyy, g_0_y_yzzz_xxyyz, g_0_y_yzzz_xxyzz, g_0_y_yzzz_xxzzz, g_0_y_yzzz_xyyyy, g_0_y_yzzz_xyyyz, g_0_y_yzzz_xyyzz, g_0_y_yzzz_xyzzz, g_0_y_yzzz_xzzzz, g_0_y_yzzz_yyyyy, g_0_y_yzzz_yyyyz, g_0_y_yzzz_yyyzz, g_0_y_yzzz_yyzzz, g_0_y_yzzz_yzzzz, g_0_y_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_xxxxx[k] = -g_0_y_yzz_xxxxx[k] * ab_z + g_0_y_yzz_xxxxxz[k];

                g_0_y_yzzz_xxxxy[k] = -g_0_y_yzz_xxxxy[k] * ab_z + g_0_y_yzz_xxxxyz[k];

                g_0_y_yzzz_xxxxz[k] = -g_0_y_yzz_xxxxz[k] * ab_z + g_0_y_yzz_xxxxzz[k];

                g_0_y_yzzz_xxxyy[k] = -g_0_y_yzz_xxxyy[k] * ab_z + g_0_y_yzz_xxxyyz[k];

                g_0_y_yzzz_xxxyz[k] = -g_0_y_yzz_xxxyz[k] * ab_z + g_0_y_yzz_xxxyzz[k];

                g_0_y_yzzz_xxxzz[k] = -g_0_y_yzz_xxxzz[k] * ab_z + g_0_y_yzz_xxxzzz[k];

                g_0_y_yzzz_xxyyy[k] = -g_0_y_yzz_xxyyy[k] * ab_z + g_0_y_yzz_xxyyyz[k];

                g_0_y_yzzz_xxyyz[k] = -g_0_y_yzz_xxyyz[k] * ab_z + g_0_y_yzz_xxyyzz[k];

                g_0_y_yzzz_xxyzz[k] = -g_0_y_yzz_xxyzz[k] * ab_z + g_0_y_yzz_xxyzzz[k];

                g_0_y_yzzz_xxzzz[k] = -g_0_y_yzz_xxzzz[k] * ab_z + g_0_y_yzz_xxzzzz[k];

                g_0_y_yzzz_xyyyy[k] = -g_0_y_yzz_xyyyy[k] * ab_z + g_0_y_yzz_xyyyyz[k];

                g_0_y_yzzz_xyyyz[k] = -g_0_y_yzz_xyyyz[k] * ab_z + g_0_y_yzz_xyyyzz[k];

                g_0_y_yzzz_xyyzz[k] = -g_0_y_yzz_xyyzz[k] * ab_z + g_0_y_yzz_xyyzzz[k];

                g_0_y_yzzz_xyzzz[k] = -g_0_y_yzz_xyzzz[k] * ab_z + g_0_y_yzz_xyzzzz[k];

                g_0_y_yzzz_xzzzz[k] = -g_0_y_yzz_xzzzz[k] * ab_z + g_0_y_yzz_xzzzzz[k];

                g_0_y_yzzz_yyyyy[k] = -g_0_y_yzz_yyyyy[k] * ab_z + g_0_y_yzz_yyyyyz[k];

                g_0_y_yzzz_yyyyz[k] = -g_0_y_yzz_yyyyz[k] * ab_z + g_0_y_yzz_yyyyzz[k];

                g_0_y_yzzz_yyyzz[k] = -g_0_y_yzz_yyyzz[k] * ab_z + g_0_y_yzz_yyyzzz[k];

                g_0_y_yzzz_yyzzz[k] = -g_0_y_yzz_yyzzz[k] * ab_z + g_0_y_yzz_yyzzzz[k];

                g_0_y_yzzz_yzzzz[k] = -g_0_y_yzz_yzzzz[k] * ab_z + g_0_y_yzz_yzzzzz[k];

                g_0_y_yzzz_zzzzz[k] = -g_0_y_yzz_zzzzz[k] * ab_z + g_0_y_yzz_zzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_zzz_xxxxx, g_0_y_zzz_xxxxxz, g_0_y_zzz_xxxxy, g_0_y_zzz_xxxxyz, g_0_y_zzz_xxxxz, g_0_y_zzz_xxxxzz, g_0_y_zzz_xxxyy, g_0_y_zzz_xxxyyz, g_0_y_zzz_xxxyz, g_0_y_zzz_xxxyzz, g_0_y_zzz_xxxzz, g_0_y_zzz_xxxzzz, g_0_y_zzz_xxyyy, g_0_y_zzz_xxyyyz, g_0_y_zzz_xxyyz, g_0_y_zzz_xxyyzz, g_0_y_zzz_xxyzz, g_0_y_zzz_xxyzzz, g_0_y_zzz_xxzzz, g_0_y_zzz_xxzzzz, g_0_y_zzz_xyyyy, g_0_y_zzz_xyyyyz, g_0_y_zzz_xyyyz, g_0_y_zzz_xyyyzz, g_0_y_zzz_xyyzz, g_0_y_zzz_xyyzzz, g_0_y_zzz_xyzzz, g_0_y_zzz_xyzzzz, g_0_y_zzz_xzzzz, g_0_y_zzz_xzzzzz, g_0_y_zzz_yyyyy, g_0_y_zzz_yyyyyz, g_0_y_zzz_yyyyz, g_0_y_zzz_yyyyzz, g_0_y_zzz_yyyzz, g_0_y_zzz_yyyzzz, g_0_y_zzz_yyzzz, g_0_y_zzz_yyzzzz, g_0_y_zzz_yzzzz, g_0_y_zzz_yzzzzz, g_0_y_zzz_zzzzz, g_0_y_zzz_zzzzzz, g_0_y_zzzz_xxxxx, g_0_y_zzzz_xxxxy, g_0_y_zzzz_xxxxz, g_0_y_zzzz_xxxyy, g_0_y_zzzz_xxxyz, g_0_y_zzzz_xxxzz, g_0_y_zzzz_xxyyy, g_0_y_zzzz_xxyyz, g_0_y_zzzz_xxyzz, g_0_y_zzzz_xxzzz, g_0_y_zzzz_xyyyy, g_0_y_zzzz_xyyyz, g_0_y_zzzz_xyyzz, g_0_y_zzzz_xyzzz, g_0_y_zzzz_xzzzz, g_0_y_zzzz_yyyyy, g_0_y_zzzz_yyyyz, g_0_y_zzzz_yyyzz, g_0_y_zzzz_yyzzz, g_0_y_zzzz_yzzzz, g_0_y_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_xxxxx[k] = -g_0_y_zzz_xxxxx[k] * ab_z + g_0_y_zzz_xxxxxz[k];

                g_0_y_zzzz_xxxxy[k] = -g_0_y_zzz_xxxxy[k] * ab_z + g_0_y_zzz_xxxxyz[k];

                g_0_y_zzzz_xxxxz[k] = -g_0_y_zzz_xxxxz[k] * ab_z + g_0_y_zzz_xxxxzz[k];

                g_0_y_zzzz_xxxyy[k] = -g_0_y_zzz_xxxyy[k] * ab_z + g_0_y_zzz_xxxyyz[k];

                g_0_y_zzzz_xxxyz[k] = -g_0_y_zzz_xxxyz[k] * ab_z + g_0_y_zzz_xxxyzz[k];

                g_0_y_zzzz_xxxzz[k] = -g_0_y_zzz_xxxzz[k] * ab_z + g_0_y_zzz_xxxzzz[k];

                g_0_y_zzzz_xxyyy[k] = -g_0_y_zzz_xxyyy[k] * ab_z + g_0_y_zzz_xxyyyz[k];

                g_0_y_zzzz_xxyyz[k] = -g_0_y_zzz_xxyyz[k] * ab_z + g_0_y_zzz_xxyyzz[k];

                g_0_y_zzzz_xxyzz[k] = -g_0_y_zzz_xxyzz[k] * ab_z + g_0_y_zzz_xxyzzz[k];

                g_0_y_zzzz_xxzzz[k] = -g_0_y_zzz_xxzzz[k] * ab_z + g_0_y_zzz_xxzzzz[k];

                g_0_y_zzzz_xyyyy[k] = -g_0_y_zzz_xyyyy[k] * ab_z + g_0_y_zzz_xyyyyz[k];

                g_0_y_zzzz_xyyyz[k] = -g_0_y_zzz_xyyyz[k] * ab_z + g_0_y_zzz_xyyyzz[k];

                g_0_y_zzzz_xyyzz[k] = -g_0_y_zzz_xyyzz[k] * ab_z + g_0_y_zzz_xyyzzz[k];

                g_0_y_zzzz_xyzzz[k] = -g_0_y_zzz_xyzzz[k] * ab_z + g_0_y_zzz_xyzzzz[k];

                g_0_y_zzzz_xzzzz[k] = -g_0_y_zzz_xzzzz[k] * ab_z + g_0_y_zzz_xzzzzz[k];

                g_0_y_zzzz_yyyyy[k] = -g_0_y_zzz_yyyyy[k] * ab_z + g_0_y_zzz_yyyyyz[k];

                g_0_y_zzzz_yyyyz[k] = -g_0_y_zzz_yyyyz[k] * ab_z + g_0_y_zzz_yyyyzz[k];

                g_0_y_zzzz_yyyzz[k] = -g_0_y_zzz_yyyzz[k] * ab_z + g_0_y_zzz_yyyzzz[k];

                g_0_y_zzzz_yyzzz[k] = -g_0_y_zzz_yyzzz[k] * ab_z + g_0_y_zzz_yyzzzz[k];

                g_0_y_zzzz_yzzzz[k] = -g_0_y_zzz_yzzzz[k] * ab_z + g_0_y_zzz_yzzzzz[k];

                g_0_y_zzzz_zzzzz[k] = -g_0_y_zzz_zzzzz[k] * ab_z + g_0_y_zzz_zzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxx_xxxxx, g_0_z_xxx_xxxxxx, g_0_z_xxx_xxxxxy, g_0_z_xxx_xxxxxz, g_0_z_xxx_xxxxy, g_0_z_xxx_xxxxyy, g_0_z_xxx_xxxxyz, g_0_z_xxx_xxxxz, g_0_z_xxx_xxxxzz, g_0_z_xxx_xxxyy, g_0_z_xxx_xxxyyy, g_0_z_xxx_xxxyyz, g_0_z_xxx_xxxyz, g_0_z_xxx_xxxyzz, g_0_z_xxx_xxxzz, g_0_z_xxx_xxxzzz, g_0_z_xxx_xxyyy, g_0_z_xxx_xxyyyy, g_0_z_xxx_xxyyyz, g_0_z_xxx_xxyyz, g_0_z_xxx_xxyyzz, g_0_z_xxx_xxyzz, g_0_z_xxx_xxyzzz, g_0_z_xxx_xxzzz, g_0_z_xxx_xxzzzz, g_0_z_xxx_xyyyy, g_0_z_xxx_xyyyyy, g_0_z_xxx_xyyyyz, g_0_z_xxx_xyyyz, g_0_z_xxx_xyyyzz, g_0_z_xxx_xyyzz, g_0_z_xxx_xyyzzz, g_0_z_xxx_xyzzz, g_0_z_xxx_xyzzzz, g_0_z_xxx_xzzzz, g_0_z_xxx_xzzzzz, g_0_z_xxx_yyyyy, g_0_z_xxx_yyyyz, g_0_z_xxx_yyyzz, g_0_z_xxx_yyzzz, g_0_z_xxx_yzzzz, g_0_z_xxx_zzzzz, g_0_z_xxxx_xxxxx, g_0_z_xxxx_xxxxy, g_0_z_xxxx_xxxxz, g_0_z_xxxx_xxxyy, g_0_z_xxxx_xxxyz, g_0_z_xxxx_xxxzz, g_0_z_xxxx_xxyyy, g_0_z_xxxx_xxyyz, g_0_z_xxxx_xxyzz, g_0_z_xxxx_xxzzz, g_0_z_xxxx_xyyyy, g_0_z_xxxx_xyyyz, g_0_z_xxxx_xyyzz, g_0_z_xxxx_xyzzz, g_0_z_xxxx_xzzzz, g_0_z_xxxx_yyyyy, g_0_z_xxxx_yyyyz, g_0_z_xxxx_yyyzz, g_0_z_xxxx_yyzzz, g_0_z_xxxx_yzzzz, g_0_z_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_xxxxx[k] = -g_0_z_xxx_xxxxx[k] * ab_x + g_0_z_xxx_xxxxxx[k];

                g_0_z_xxxx_xxxxy[k] = -g_0_z_xxx_xxxxy[k] * ab_x + g_0_z_xxx_xxxxxy[k];

                g_0_z_xxxx_xxxxz[k] = -g_0_z_xxx_xxxxz[k] * ab_x + g_0_z_xxx_xxxxxz[k];

                g_0_z_xxxx_xxxyy[k] = -g_0_z_xxx_xxxyy[k] * ab_x + g_0_z_xxx_xxxxyy[k];

                g_0_z_xxxx_xxxyz[k] = -g_0_z_xxx_xxxyz[k] * ab_x + g_0_z_xxx_xxxxyz[k];

                g_0_z_xxxx_xxxzz[k] = -g_0_z_xxx_xxxzz[k] * ab_x + g_0_z_xxx_xxxxzz[k];

                g_0_z_xxxx_xxyyy[k] = -g_0_z_xxx_xxyyy[k] * ab_x + g_0_z_xxx_xxxyyy[k];

                g_0_z_xxxx_xxyyz[k] = -g_0_z_xxx_xxyyz[k] * ab_x + g_0_z_xxx_xxxyyz[k];

                g_0_z_xxxx_xxyzz[k] = -g_0_z_xxx_xxyzz[k] * ab_x + g_0_z_xxx_xxxyzz[k];

                g_0_z_xxxx_xxzzz[k] = -g_0_z_xxx_xxzzz[k] * ab_x + g_0_z_xxx_xxxzzz[k];

                g_0_z_xxxx_xyyyy[k] = -g_0_z_xxx_xyyyy[k] * ab_x + g_0_z_xxx_xxyyyy[k];

                g_0_z_xxxx_xyyyz[k] = -g_0_z_xxx_xyyyz[k] * ab_x + g_0_z_xxx_xxyyyz[k];

                g_0_z_xxxx_xyyzz[k] = -g_0_z_xxx_xyyzz[k] * ab_x + g_0_z_xxx_xxyyzz[k];

                g_0_z_xxxx_xyzzz[k] = -g_0_z_xxx_xyzzz[k] * ab_x + g_0_z_xxx_xxyzzz[k];

                g_0_z_xxxx_xzzzz[k] = -g_0_z_xxx_xzzzz[k] * ab_x + g_0_z_xxx_xxzzzz[k];

                g_0_z_xxxx_yyyyy[k] = -g_0_z_xxx_yyyyy[k] * ab_x + g_0_z_xxx_xyyyyy[k];

                g_0_z_xxxx_yyyyz[k] = -g_0_z_xxx_yyyyz[k] * ab_x + g_0_z_xxx_xyyyyz[k];

                g_0_z_xxxx_yyyzz[k] = -g_0_z_xxx_yyyzz[k] * ab_x + g_0_z_xxx_xyyyzz[k];

                g_0_z_xxxx_yyzzz[k] = -g_0_z_xxx_yyzzz[k] * ab_x + g_0_z_xxx_xyyzzz[k];

                g_0_z_xxxx_yzzzz[k] = -g_0_z_xxx_yzzzz[k] * ab_x + g_0_z_xxx_xyzzzz[k];

                g_0_z_xxxx_zzzzz[k] = -g_0_z_xxx_zzzzz[k] * ab_x + g_0_z_xxx_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxxy_xxxxx, g_0_z_xxxy_xxxxy, g_0_z_xxxy_xxxxz, g_0_z_xxxy_xxxyy, g_0_z_xxxy_xxxyz, g_0_z_xxxy_xxxzz, g_0_z_xxxy_xxyyy, g_0_z_xxxy_xxyyz, g_0_z_xxxy_xxyzz, g_0_z_xxxy_xxzzz, g_0_z_xxxy_xyyyy, g_0_z_xxxy_xyyyz, g_0_z_xxxy_xyyzz, g_0_z_xxxy_xyzzz, g_0_z_xxxy_xzzzz, g_0_z_xxxy_yyyyy, g_0_z_xxxy_yyyyz, g_0_z_xxxy_yyyzz, g_0_z_xxxy_yyzzz, g_0_z_xxxy_yzzzz, g_0_z_xxxy_zzzzz, g_0_z_xxy_xxxxx, g_0_z_xxy_xxxxxx, g_0_z_xxy_xxxxxy, g_0_z_xxy_xxxxxz, g_0_z_xxy_xxxxy, g_0_z_xxy_xxxxyy, g_0_z_xxy_xxxxyz, g_0_z_xxy_xxxxz, g_0_z_xxy_xxxxzz, g_0_z_xxy_xxxyy, g_0_z_xxy_xxxyyy, g_0_z_xxy_xxxyyz, g_0_z_xxy_xxxyz, g_0_z_xxy_xxxyzz, g_0_z_xxy_xxxzz, g_0_z_xxy_xxxzzz, g_0_z_xxy_xxyyy, g_0_z_xxy_xxyyyy, g_0_z_xxy_xxyyyz, g_0_z_xxy_xxyyz, g_0_z_xxy_xxyyzz, g_0_z_xxy_xxyzz, g_0_z_xxy_xxyzzz, g_0_z_xxy_xxzzz, g_0_z_xxy_xxzzzz, g_0_z_xxy_xyyyy, g_0_z_xxy_xyyyyy, g_0_z_xxy_xyyyyz, g_0_z_xxy_xyyyz, g_0_z_xxy_xyyyzz, g_0_z_xxy_xyyzz, g_0_z_xxy_xyyzzz, g_0_z_xxy_xyzzz, g_0_z_xxy_xyzzzz, g_0_z_xxy_xzzzz, g_0_z_xxy_xzzzzz, g_0_z_xxy_yyyyy, g_0_z_xxy_yyyyz, g_0_z_xxy_yyyzz, g_0_z_xxy_yyzzz, g_0_z_xxy_yzzzz, g_0_z_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_xxxxx[k] = -g_0_z_xxy_xxxxx[k] * ab_x + g_0_z_xxy_xxxxxx[k];

                g_0_z_xxxy_xxxxy[k] = -g_0_z_xxy_xxxxy[k] * ab_x + g_0_z_xxy_xxxxxy[k];

                g_0_z_xxxy_xxxxz[k] = -g_0_z_xxy_xxxxz[k] * ab_x + g_0_z_xxy_xxxxxz[k];

                g_0_z_xxxy_xxxyy[k] = -g_0_z_xxy_xxxyy[k] * ab_x + g_0_z_xxy_xxxxyy[k];

                g_0_z_xxxy_xxxyz[k] = -g_0_z_xxy_xxxyz[k] * ab_x + g_0_z_xxy_xxxxyz[k];

                g_0_z_xxxy_xxxzz[k] = -g_0_z_xxy_xxxzz[k] * ab_x + g_0_z_xxy_xxxxzz[k];

                g_0_z_xxxy_xxyyy[k] = -g_0_z_xxy_xxyyy[k] * ab_x + g_0_z_xxy_xxxyyy[k];

                g_0_z_xxxy_xxyyz[k] = -g_0_z_xxy_xxyyz[k] * ab_x + g_0_z_xxy_xxxyyz[k];

                g_0_z_xxxy_xxyzz[k] = -g_0_z_xxy_xxyzz[k] * ab_x + g_0_z_xxy_xxxyzz[k];

                g_0_z_xxxy_xxzzz[k] = -g_0_z_xxy_xxzzz[k] * ab_x + g_0_z_xxy_xxxzzz[k];

                g_0_z_xxxy_xyyyy[k] = -g_0_z_xxy_xyyyy[k] * ab_x + g_0_z_xxy_xxyyyy[k];

                g_0_z_xxxy_xyyyz[k] = -g_0_z_xxy_xyyyz[k] * ab_x + g_0_z_xxy_xxyyyz[k];

                g_0_z_xxxy_xyyzz[k] = -g_0_z_xxy_xyyzz[k] * ab_x + g_0_z_xxy_xxyyzz[k];

                g_0_z_xxxy_xyzzz[k] = -g_0_z_xxy_xyzzz[k] * ab_x + g_0_z_xxy_xxyzzz[k];

                g_0_z_xxxy_xzzzz[k] = -g_0_z_xxy_xzzzz[k] * ab_x + g_0_z_xxy_xxzzzz[k];

                g_0_z_xxxy_yyyyy[k] = -g_0_z_xxy_yyyyy[k] * ab_x + g_0_z_xxy_xyyyyy[k];

                g_0_z_xxxy_yyyyz[k] = -g_0_z_xxy_yyyyz[k] * ab_x + g_0_z_xxy_xyyyyz[k];

                g_0_z_xxxy_yyyzz[k] = -g_0_z_xxy_yyyzz[k] * ab_x + g_0_z_xxy_xyyyzz[k];

                g_0_z_xxxy_yyzzz[k] = -g_0_z_xxy_yyzzz[k] * ab_x + g_0_z_xxy_xyyzzz[k];

                g_0_z_xxxy_yzzzz[k] = -g_0_z_xxy_yzzzz[k] * ab_x + g_0_z_xxy_xyzzzz[k];

                g_0_z_xxxy_zzzzz[k] = -g_0_z_xxy_zzzzz[k] * ab_x + g_0_z_xxy_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxxz_xxxxx, g_0_z_xxxz_xxxxy, g_0_z_xxxz_xxxxz, g_0_z_xxxz_xxxyy, g_0_z_xxxz_xxxyz, g_0_z_xxxz_xxxzz, g_0_z_xxxz_xxyyy, g_0_z_xxxz_xxyyz, g_0_z_xxxz_xxyzz, g_0_z_xxxz_xxzzz, g_0_z_xxxz_xyyyy, g_0_z_xxxz_xyyyz, g_0_z_xxxz_xyyzz, g_0_z_xxxz_xyzzz, g_0_z_xxxz_xzzzz, g_0_z_xxxz_yyyyy, g_0_z_xxxz_yyyyz, g_0_z_xxxz_yyyzz, g_0_z_xxxz_yyzzz, g_0_z_xxxz_yzzzz, g_0_z_xxxz_zzzzz, g_0_z_xxz_xxxxx, g_0_z_xxz_xxxxxx, g_0_z_xxz_xxxxxy, g_0_z_xxz_xxxxxz, g_0_z_xxz_xxxxy, g_0_z_xxz_xxxxyy, g_0_z_xxz_xxxxyz, g_0_z_xxz_xxxxz, g_0_z_xxz_xxxxzz, g_0_z_xxz_xxxyy, g_0_z_xxz_xxxyyy, g_0_z_xxz_xxxyyz, g_0_z_xxz_xxxyz, g_0_z_xxz_xxxyzz, g_0_z_xxz_xxxzz, g_0_z_xxz_xxxzzz, g_0_z_xxz_xxyyy, g_0_z_xxz_xxyyyy, g_0_z_xxz_xxyyyz, g_0_z_xxz_xxyyz, g_0_z_xxz_xxyyzz, g_0_z_xxz_xxyzz, g_0_z_xxz_xxyzzz, g_0_z_xxz_xxzzz, g_0_z_xxz_xxzzzz, g_0_z_xxz_xyyyy, g_0_z_xxz_xyyyyy, g_0_z_xxz_xyyyyz, g_0_z_xxz_xyyyz, g_0_z_xxz_xyyyzz, g_0_z_xxz_xyyzz, g_0_z_xxz_xyyzzz, g_0_z_xxz_xyzzz, g_0_z_xxz_xyzzzz, g_0_z_xxz_xzzzz, g_0_z_xxz_xzzzzz, g_0_z_xxz_yyyyy, g_0_z_xxz_yyyyz, g_0_z_xxz_yyyzz, g_0_z_xxz_yyzzz, g_0_z_xxz_yzzzz, g_0_z_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_xxxxx[k] = -g_0_z_xxz_xxxxx[k] * ab_x + g_0_z_xxz_xxxxxx[k];

                g_0_z_xxxz_xxxxy[k] = -g_0_z_xxz_xxxxy[k] * ab_x + g_0_z_xxz_xxxxxy[k];

                g_0_z_xxxz_xxxxz[k] = -g_0_z_xxz_xxxxz[k] * ab_x + g_0_z_xxz_xxxxxz[k];

                g_0_z_xxxz_xxxyy[k] = -g_0_z_xxz_xxxyy[k] * ab_x + g_0_z_xxz_xxxxyy[k];

                g_0_z_xxxz_xxxyz[k] = -g_0_z_xxz_xxxyz[k] * ab_x + g_0_z_xxz_xxxxyz[k];

                g_0_z_xxxz_xxxzz[k] = -g_0_z_xxz_xxxzz[k] * ab_x + g_0_z_xxz_xxxxzz[k];

                g_0_z_xxxz_xxyyy[k] = -g_0_z_xxz_xxyyy[k] * ab_x + g_0_z_xxz_xxxyyy[k];

                g_0_z_xxxz_xxyyz[k] = -g_0_z_xxz_xxyyz[k] * ab_x + g_0_z_xxz_xxxyyz[k];

                g_0_z_xxxz_xxyzz[k] = -g_0_z_xxz_xxyzz[k] * ab_x + g_0_z_xxz_xxxyzz[k];

                g_0_z_xxxz_xxzzz[k] = -g_0_z_xxz_xxzzz[k] * ab_x + g_0_z_xxz_xxxzzz[k];

                g_0_z_xxxz_xyyyy[k] = -g_0_z_xxz_xyyyy[k] * ab_x + g_0_z_xxz_xxyyyy[k];

                g_0_z_xxxz_xyyyz[k] = -g_0_z_xxz_xyyyz[k] * ab_x + g_0_z_xxz_xxyyyz[k];

                g_0_z_xxxz_xyyzz[k] = -g_0_z_xxz_xyyzz[k] * ab_x + g_0_z_xxz_xxyyzz[k];

                g_0_z_xxxz_xyzzz[k] = -g_0_z_xxz_xyzzz[k] * ab_x + g_0_z_xxz_xxyzzz[k];

                g_0_z_xxxz_xzzzz[k] = -g_0_z_xxz_xzzzz[k] * ab_x + g_0_z_xxz_xxzzzz[k];

                g_0_z_xxxz_yyyyy[k] = -g_0_z_xxz_yyyyy[k] * ab_x + g_0_z_xxz_xyyyyy[k];

                g_0_z_xxxz_yyyyz[k] = -g_0_z_xxz_yyyyz[k] * ab_x + g_0_z_xxz_xyyyyz[k];

                g_0_z_xxxz_yyyzz[k] = -g_0_z_xxz_yyyzz[k] * ab_x + g_0_z_xxz_xyyyzz[k];

                g_0_z_xxxz_yyzzz[k] = -g_0_z_xxz_yyzzz[k] * ab_x + g_0_z_xxz_xyyzzz[k];

                g_0_z_xxxz_yzzzz[k] = -g_0_z_xxz_yzzzz[k] * ab_x + g_0_z_xxz_xyzzzz[k];

                g_0_z_xxxz_zzzzz[k] = -g_0_z_xxz_zzzzz[k] * ab_x + g_0_z_xxz_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxyy_xxxxx, g_0_z_xxyy_xxxxy, g_0_z_xxyy_xxxxz, g_0_z_xxyy_xxxyy, g_0_z_xxyy_xxxyz, g_0_z_xxyy_xxxzz, g_0_z_xxyy_xxyyy, g_0_z_xxyy_xxyyz, g_0_z_xxyy_xxyzz, g_0_z_xxyy_xxzzz, g_0_z_xxyy_xyyyy, g_0_z_xxyy_xyyyz, g_0_z_xxyy_xyyzz, g_0_z_xxyy_xyzzz, g_0_z_xxyy_xzzzz, g_0_z_xxyy_yyyyy, g_0_z_xxyy_yyyyz, g_0_z_xxyy_yyyzz, g_0_z_xxyy_yyzzz, g_0_z_xxyy_yzzzz, g_0_z_xxyy_zzzzz, g_0_z_xyy_xxxxx, g_0_z_xyy_xxxxxx, g_0_z_xyy_xxxxxy, g_0_z_xyy_xxxxxz, g_0_z_xyy_xxxxy, g_0_z_xyy_xxxxyy, g_0_z_xyy_xxxxyz, g_0_z_xyy_xxxxz, g_0_z_xyy_xxxxzz, g_0_z_xyy_xxxyy, g_0_z_xyy_xxxyyy, g_0_z_xyy_xxxyyz, g_0_z_xyy_xxxyz, g_0_z_xyy_xxxyzz, g_0_z_xyy_xxxzz, g_0_z_xyy_xxxzzz, g_0_z_xyy_xxyyy, g_0_z_xyy_xxyyyy, g_0_z_xyy_xxyyyz, g_0_z_xyy_xxyyz, g_0_z_xyy_xxyyzz, g_0_z_xyy_xxyzz, g_0_z_xyy_xxyzzz, g_0_z_xyy_xxzzz, g_0_z_xyy_xxzzzz, g_0_z_xyy_xyyyy, g_0_z_xyy_xyyyyy, g_0_z_xyy_xyyyyz, g_0_z_xyy_xyyyz, g_0_z_xyy_xyyyzz, g_0_z_xyy_xyyzz, g_0_z_xyy_xyyzzz, g_0_z_xyy_xyzzz, g_0_z_xyy_xyzzzz, g_0_z_xyy_xzzzz, g_0_z_xyy_xzzzzz, g_0_z_xyy_yyyyy, g_0_z_xyy_yyyyz, g_0_z_xyy_yyyzz, g_0_z_xyy_yyzzz, g_0_z_xyy_yzzzz, g_0_z_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_xxxxx[k] = -g_0_z_xyy_xxxxx[k] * ab_x + g_0_z_xyy_xxxxxx[k];

                g_0_z_xxyy_xxxxy[k] = -g_0_z_xyy_xxxxy[k] * ab_x + g_0_z_xyy_xxxxxy[k];

                g_0_z_xxyy_xxxxz[k] = -g_0_z_xyy_xxxxz[k] * ab_x + g_0_z_xyy_xxxxxz[k];

                g_0_z_xxyy_xxxyy[k] = -g_0_z_xyy_xxxyy[k] * ab_x + g_0_z_xyy_xxxxyy[k];

                g_0_z_xxyy_xxxyz[k] = -g_0_z_xyy_xxxyz[k] * ab_x + g_0_z_xyy_xxxxyz[k];

                g_0_z_xxyy_xxxzz[k] = -g_0_z_xyy_xxxzz[k] * ab_x + g_0_z_xyy_xxxxzz[k];

                g_0_z_xxyy_xxyyy[k] = -g_0_z_xyy_xxyyy[k] * ab_x + g_0_z_xyy_xxxyyy[k];

                g_0_z_xxyy_xxyyz[k] = -g_0_z_xyy_xxyyz[k] * ab_x + g_0_z_xyy_xxxyyz[k];

                g_0_z_xxyy_xxyzz[k] = -g_0_z_xyy_xxyzz[k] * ab_x + g_0_z_xyy_xxxyzz[k];

                g_0_z_xxyy_xxzzz[k] = -g_0_z_xyy_xxzzz[k] * ab_x + g_0_z_xyy_xxxzzz[k];

                g_0_z_xxyy_xyyyy[k] = -g_0_z_xyy_xyyyy[k] * ab_x + g_0_z_xyy_xxyyyy[k];

                g_0_z_xxyy_xyyyz[k] = -g_0_z_xyy_xyyyz[k] * ab_x + g_0_z_xyy_xxyyyz[k];

                g_0_z_xxyy_xyyzz[k] = -g_0_z_xyy_xyyzz[k] * ab_x + g_0_z_xyy_xxyyzz[k];

                g_0_z_xxyy_xyzzz[k] = -g_0_z_xyy_xyzzz[k] * ab_x + g_0_z_xyy_xxyzzz[k];

                g_0_z_xxyy_xzzzz[k] = -g_0_z_xyy_xzzzz[k] * ab_x + g_0_z_xyy_xxzzzz[k];

                g_0_z_xxyy_yyyyy[k] = -g_0_z_xyy_yyyyy[k] * ab_x + g_0_z_xyy_xyyyyy[k];

                g_0_z_xxyy_yyyyz[k] = -g_0_z_xyy_yyyyz[k] * ab_x + g_0_z_xyy_xyyyyz[k];

                g_0_z_xxyy_yyyzz[k] = -g_0_z_xyy_yyyzz[k] * ab_x + g_0_z_xyy_xyyyzz[k];

                g_0_z_xxyy_yyzzz[k] = -g_0_z_xyy_yyzzz[k] * ab_x + g_0_z_xyy_xyyzzz[k];

                g_0_z_xxyy_yzzzz[k] = -g_0_z_xyy_yzzzz[k] * ab_x + g_0_z_xyy_xyzzzz[k];

                g_0_z_xxyy_zzzzz[k] = -g_0_z_xyy_zzzzz[k] * ab_x + g_0_z_xyy_xzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxyz_xxxxx, g_0_z_xxyz_xxxxy, g_0_z_xxyz_xxxxz, g_0_z_xxyz_xxxyy, g_0_z_xxyz_xxxyz, g_0_z_xxyz_xxxzz, g_0_z_xxyz_xxyyy, g_0_z_xxyz_xxyyz, g_0_z_xxyz_xxyzz, g_0_z_xxyz_xxzzz, g_0_z_xxyz_xyyyy, g_0_z_xxyz_xyyyz, g_0_z_xxyz_xyyzz, g_0_z_xxyz_xyzzz, g_0_z_xxyz_xzzzz, g_0_z_xxyz_yyyyy, g_0_z_xxyz_yyyyz, g_0_z_xxyz_yyyzz, g_0_z_xxyz_yyzzz, g_0_z_xxyz_yzzzz, g_0_z_xxyz_zzzzz, g_0_z_xyz_xxxxx, g_0_z_xyz_xxxxxx, g_0_z_xyz_xxxxxy, g_0_z_xyz_xxxxxz, g_0_z_xyz_xxxxy, g_0_z_xyz_xxxxyy, g_0_z_xyz_xxxxyz, g_0_z_xyz_xxxxz, g_0_z_xyz_xxxxzz, g_0_z_xyz_xxxyy, g_0_z_xyz_xxxyyy, g_0_z_xyz_xxxyyz, g_0_z_xyz_xxxyz, g_0_z_xyz_xxxyzz, g_0_z_xyz_xxxzz, g_0_z_xyz_xxxzzz, g_0_z_xyz_xxyyy, g_0_z_xyz_xxyyyy, g_0_z_xyz_xxyyyz, g_0_z_xyz_xxyyz, g_0_z_xyz_xxyyzz, g_0_z_xyz_xxyzz, g_0_z_xyz_xxyzzz, g_0_z_xyz_xxzzz, g_0_z_xyz_xxzzzz, g_0_z_xyz_xyyyy, g_0_z_xyz_xyyyyy, g_0_z_xyz_xyyyyz, g_0_z_xyz_xyyyz, g_0_z_xyz_xyyyzz, g_0_z_xyz_xyyzz, g_0_z_xyz_xyyzzz, g_0_z_xyz_xyzzz, g_0_z_xyz_xyzzzz, g_0_z_xyz_xzzzz, g_0_z_xyz_xzzzzz, g_0_z_xyz_yyyyy, g_0_z_xyz_yyyyz, g_0_z_xyz_yyyzz, g_0_z_xyz_yyzzz, g_0_z_xyz_yzzzz, g_0_z_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_xxxxx[k] = -g_0_z_xyz_xxxxx[k] * ab_x + g_0_z_xyz_xxxxxx[k];

                g_0_z_xxyz_xxxxy[k] = -g_0_z_xyz_xxxxy[k] * ab_x + g_0_z_xyz_xxxxxy[k];

                g_0_z_xxyz_xxxxz[k] = -g_0_z_xyz_xxxxz[k] * ab_x + g_0_z_xyz_xxxxxz[k];

                g_0_z_xxyz_xxxyy[k] = -g_0_z_xyz_xxxyy[k] * ab_x + g_0_z_xyz_xxxxyy[k];

                g_0_z_xxyz_xxxyz[k] = -g_0_z_xyz_xxxyz[k] * ab_x + g_0_z_xyz_xxxxyz[k];

                g_0_z_xxyz_xxxzz[k] = -g_0_z_xyz_xxxzz[k] * ab_x + g_0_z_xyz_xxxxzz[k];

                g_0_z_xxyz_xxyyy[k] = -g_0_z_xyz_xxyyy[k] * ab_x + g_0_z_xyz_xxxyyy[k];

                g_0_z_xxyz_xxyyz[k] = -g_0_z_xyz_xxyyz[k] * ab_x + g_0_z_xyz_xxxyyz[k];

                g_0_z_xxyz_xxyzz[k] = -g_0_z_xyz_xxyzz[k] * ab_x + g_0_z_xyz_xxxyzz[k];

                g_0_z_xxyz_xxzzz[k] = -g_0_z_xyz_xxzzz[k] * ab_x + g_0_z_xyz_xxxzzz[k];

                g_0_z_xxyz_xyyyy[k] = -g_0_z_xyz_xyyyy[k] * ab_x + g_0_z_xyz_xxyyyy[k];

                g_0_z_xxyz_xyyyz[k] = -g_0_z_xyz_xyyyz[k] * ab_x + g_0_z_xyz_xxyyyz[k];

                g_0_z_xxyz_xyyzz[k] = -g_0_z_xyz_xyyzz[k] * ab_x + g_0_z_xyz_xxyyzz[k];

                g_0_z_xxyz_xyzzz[k] = -g_0_z_xyz_xyzzz[k] * ab_x + g_0_z_xyz_xxyzzz[k];

                g_0_z_xxyz_xzzzz[k] = -g_0_z_xyz_xzzzz[k] * ab_x + g_0_z_xyz_xxzzzz[k];

                g_0_z_xxyz_yyyyy[k] = -g_0_z_xyz_yyyyy[k] * ab_x + g_0_z_xyz_xyyyyy[k];

                g_0_z_xxyz_yyyyz[k] = -g_0_z_xyz_yyyyz[k] * ab_x + g_0_z_xyz_xyyyyz[k];

                g_0_z_xxyz_yyyzz[k] = -g_0_z_xyz_yyyzz[k] * ab_x + g_0_z_xyz_xyyyzz[k];

                g_0_z_xxyz_yyzzz[k] = -g_0_z_xyz_yyzzz[k] * ab_x + g_0_z_xyz_xyyzzz[k];

                g_0_z_xxyz_yzzzz[k] = -g_0_z_xyz_yzzzz[k] * ab_x + g_0_z_xyz_xyzzzz[k];

                g_0_z_xxyz_zzzzz[k] = -g_0_z_xyz_zzzzz[k] * ab_x + g_0_z_xyz_xzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxzz_xxxxx, g_0_z_xxzz_xxxxy, g_0_z_xxzz_xxxxz, g_0_z_xxzz_xxxyy, g_0_z_xxzz_xxxyz, g_0_z_xxzz_xxxzz, g_0_z_xxzz_xxyyy, g_0_z_xxzz_xxyyz, g_0_z_xxzz_xxyzz, g_0_z_xxzz_xxzzz, g_0_z_xxzz_xyyyy, g_0_z_xxzz_xyyyz, g_0_z_xxzz_xyyzz, g_0_z_xxzz_xyzzz, g_0_z_xxzz_xzzzz, g_0_z_xxzz_yyyyy, g_0_z_xxzz_yyyyz, g_0_z_xxzz_yyyzz, g_0_z_xxzz_yyzzz, g_0_z_xxzz_yzzzz, g_0_z_xxzz_zzzzz, g_0_z_xzz_xxxxx, g_0_z_xzz_xxxxxx, g_0_z_xzz_xxxxxy, g_0_z_xzz_xxxxxz, g_0_z_xzz_xxxxy, g_0_z_xzz_xxxxyy, g_0_z_xzz_xxxxyz, g_0_z_xzz_xxxxz, g_0_z_xzz_xxxxzz, g_0_z_xzz_xxxyy, g_0_z_xzz_xxxyyy, g_0_z_xzz_xxxyyz, g_0_z_xzz_xxxyz, g_0_z_xzz_xxxyzz, g_0_z_xzz_xxxzz, g_0_z_xzz_xxxzzz, g_0_z_xzz_xxyyy, g_0_z_xzz_xxyyyy, g_0_z_xzz_xxyyyz, g_0_z_xzz_xxyyz, g_0_z_xzz_xxyyzz, g_0_z_xzz_xxyzz, g_0_z_xzz_xxyzzz, g_0_z_xzz_xxzzz, g_0_z_xzz_xxzzzz, g_0_z_xzz_xyyyy, g_0_z_xzz_xyyyyy, g_0_z_xzz_xyyyyz, g_0_z_xzz_xyyyz, g_0_z_xzz_xyyyzz, g_0_z_xzz_xyyzz, g_0_z_xzz_xyyzzz, g_0_z_xzz_xyzzz, g_0_z_xzz_xyzzzz, g_0_z_xzz_xzzzz, g_0_z_xzz_xzzzzz, g_0_z_xzz_yyyyy, g_0_z_xzz_yyyyz, g_0_z_xzz_yyyzz, g_0_z_xzz_yyzzz, g_0_z_xzz_yzzzz, g_0_z_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_xxxxx[k] = -g_0_z_xzz_xxxxx[k] * ab_x + g_0_z_xzz_xxxxxx[k];

                g_0_z_xxzz_xxxxy[k] = -g_0_z_xzz_xxxxy[k] * ab_x + g_0_z_xzz_xxxxxy[k];

                g_0_z_xxzz_xxxxz[k] = -g_0_z_xzz_xxxxz[k] * ab_x + g_0_z_xzz_xxxxxz[k];

                g_0_z_xxzz_xxxyy[k] = -g_0_z_xzz_xxxyy[k] * ab_x + g_0_z_xzz_xxxxyy[k];

                g_0_z_xxzz_xxxyz[k] = -g_0_z_xzz_xxxyz[k] * ab_x + g_0_z_xzz_xxxxyz[k];

                g_0_z_xxzz_xxxzz[k] = -g_0_z_xzz_xxxzz[k] * ab_x + g_0_z_xzz_xxxxzz[k];

                g_0_z_xxzz_xxyyy[k] = -g_0_z_xzz_xxyyy[k] * ab_x + g_0_z_xzz_xxxyyy[k];

                g_0_z_xxzz_xxyyz[k] = -g_0_z_xzz_xxyyz[k] * ab_x + g_0_z_xzz_xxxyyz[k];

                g_0_z_xxzz_xxyzz[k] = -g_0_z_xzz_xxyzz[k] * ab_x + g_0_z_xzz_xxxyzz[k];

                g_0_z_xxzz_xxzzz[k] = -g_0_z_xzz_xxzzz[k] * ab_x + g_0_z_xzz_xxxzzz[k];

                g_0_z_xxzz_xyyyy[k] = -g_0_z_xzz_xyyyy[k] * ab_x + g_0_z_xzz_xxyyyy[k];

                g_0_z_xxzz_xyyyz[k] = -g_0_z_xzz_xyyyz[k] * ab_x + g_0_z_xzz_xxyyyz[k];

                g_0_z_xxzz_xyyzz[k] = -g_0_z_xzz_xyyzz[k] * ab_x + g_0_z_xzz_xxyyzz[k];

                g_0_z_xxzz_xyzzz[k] = -g_0_z_xzz_xyzzz[k] * ab_x + g_0_z_xzz_xxyzzz[k];

                g_0_z_xxzz_xzzzz[k] = -g_0_z_xzz_xzzzz[k] * ab_x + g_0_z_xzz_xxzzzz[k];

                g_0_z_xxzz_yyyyy[k] = -g_0_z_xzz_yyyyy[k] * ab_x + g_0_z_xzz_xyyyyy[k];

                g_0_z_xxzz_yyyyz[k] = -g_0_z_xzz_yyyyz[k] * ab_x + g_0_z_xzz_xyyyyz[k];

                g_0_z_xxzz_yyyzz[k] = -g_0_z_xzz_yyyzz[k] * ab_x + g_0_z_xzz_xyyyzz[k];

                g_0_z_xxzz_yyzzz[k] = -g_0_z_xzz_yyzzz[k] * ab_x + g_0_z_xzz_xyyzzz[k];

                g_0_z_xxzz_yzzzz[k] = -g_0_z_xzz_yzzzz[k] * ab_x + g_0_z_xzz_xyzzzz[k];

                g_0_z_xxzz_zzzzz[k] = -g_0_z_xzz_zzzzz[k] * ab_x + g_0_z_xzz_xzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyyy_xxxxx, g_0_z_xyyy_xxxxy, g_0_z_xyyy_xxxxz, g_0_z_xyyy_xxxyy, g_0_z_xyyy_xxxyz, g_0_z_xyyy_xxxzz, g_0_z_xyyy_xxyyy, g_0_z_xyyy_xxyyz, g_0_z_xyyy_xxyzz, g_0_z_xyyy_xxzzz, g_0_z_xyyy_xyyyy, g_0_z_xyyy_xyyyz, g_0_z_xyyy_xyyzz, g_0_z_xyyy_xyzzz, g_0_z_xyyy_xzzzz, g_0_z_xyyy_yyyyy, g_0_z_xyyy_yyyyz, g_0_z_xyyy_yyyzz, g_0_z_xyyy_yyzzz, g_0_z_xyyy_yzzzz, g_0_z_xyyy_zzzzz, g_0_z_yyy_xxxxx, g_0_z_yyy_xxxxxx, g_0_z_yyy_xxxxxy, g_0_z_yyy_xxxxxz, g_0_z_yyy_xxxxy, g_0_z_yyy_xxxxyy, g_0_z_yyy_xxxxyz, g_0_z_yyy_xxxxz, g_0_z_yyy_xxxxzz, g_0_z_yyy_xxxyy, g_0_z_yyy_xxxyyy, g_0_z_yyy_xxxyyz, g_0_z_yyy_xxxyz, g_0_z_yyy_xxxyzz, g_0_z_yyy_xxxzz, g_0_z_yyy_xxxzzz, g_0_z_yyy_xxyyy, g_0_z_yyy_xxyyyy, g_0_z_yyy_xxyyyz, g_0_z_yyy_xxyyz, g_0_z_yyy_xxyyzz, g_0_z_yyy_xxyzz, g_0_z_yyy_xxyzzz, g_0_z_yyy_xxzzz, g_0_z_yyy_xxzzzz, g_0_z_yyy_xyyyy, g_0_z_yyy_xyyyyy, g_0_z_yyy_xyyyyz, g_0_z_yyy_xyyyz, g_0_z_yyy_xyyyzz, g_0_z_yyy_xyyzz, g_0_z_yyy_xyyzzz, g_0_z_yyy_xyzzz, g_0_z_yyy_xyzzzz, g_0_z_yyy_xzzzz, g_0_z_yyy_xzzzzz, g_0_z_yyy_yyyyy, g_0_z_yyy_yyyyz, g_0_z_yyy_yyyzz, g_0_z_yyy_yyzzz, g_0_z_yyy_yzzzz, g_0_z_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_xxxxx[k] = -g_0_z_yyy_xxxxx[k] * ab_x + g_0_z_yyy_xxxxxx[k];

                g_0_z_xyyy_xxxxy[k] = -g_0_z_yyy_xxxxy[k] * ab_x + g_0_z_yyy_xxxxxy[k];

                g_0_z_xyyy_xxxxz[k] = -g_0_z_yyy_xxxxz[k] * ab_x + g_0_z_yyy_xxxxxz[k];

                g_0_z_xyyy_xxxyy[k] = -g_0_z_yyy_xxxyy[k] * ab_x + g_0_z_yyy_xxxxyy[k];

                g_0_z_xyyy_xxxyz[k] = -g_0_z_yyy_xxxyz[k] * ab_x + g_0_z_yyy_xxxxyz[k];

                g_0_z_xyyy_xxxzz[k] = -g_0_z_yyy_xxxzz[k] * ab_x + g_0_z_yyy_xxxxzz[k];

                g_0_z_xyyy_xxyyy[k] = -g_0_z_yyy_xxyyy[k] * ab_x + g_0_z_yyy_xxxyyy[k];

                g_0_z_xyyy_xxyyz[k] = -g_0_z_yyy_xxyyz[k] * ab_x + g_0_z_yyy_xxxyyz[k];

                g_0_z_xyyy_xxyzz[k] = -g_0_z_yyy_xxyzz[k] * ab_x + g_0_z_yyy_xxxyzz[k];

                g_0_z_xyyy_xxzzz[k] = -g_0_z_yyy_xxzzz[k] * ab_x + g_0_z_yyy_xxxzzz[k];

                g_0_z_xyyy_xyyyy[k] = -g_0_z_yyy_xyyyy[k] * ab_x + g_0_z_yyy_xxyyyy[k];

                g_0_z_xyyy_xyyyz[k] = -g_0_z_yyy_xyyyz[k] * ab_x + g_0_z_yyy_xxyyyz[k];

                g_0_z_xyyy_xyyzz[k] = -g_0_z_yyy_xyyzz[k] * ab_x + g_0_z_yyy_xxyyzz[k];

                g_0_z_xyyy_xyzzz[k] = -g_0_z_yyy_xyzzz[k] * ab_x + g_0_z_yyy_xxyzzz[k];

                g_0_z_xyyy_xzzzz[k] = -g_0_z_yyy_xzzzz[k] * ab_x + g_0_z_yyy_xxzzzz[k];

                g_0_z_xyyy_yyyyy[k] = -g_0_z_yyy_yyyyy[k] * ab_x + g_0_z_yyy_xyyyyy[k];

                g_0_z_xyyy_yyyyz[k] = -g_0_z_yyy_yyyyz[k] * ab_x + g_0_z_yyy_xyyyyz[k];

                g_0_z_xyyy_yyyzz[k] = -g_0_z_yyy_yyyzz[k] * ab_x + g_0_z_yyy_xyyyzz[k];

                g_0_z_xyyy_yyzzz[k] = -g_0_z_yyy_yyzzz[k] * ab_x + g_0_z_yyy_xyyzzz[k];

                g_0_z_xyyy_yzzzz[k] = -g_0_z_yyy_yzzzz[k] * ab_x + g_0_z_yyy_xyzzzz[k];

                g_0_z_xyyy_zzzzz[k] = -g_0_z_yyy_zzzzz[k] * ab_x + g_0_z_yyy_xzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyyz_xxxxx, g_0_z_xyyz_xxxxy, g_0_z_xyyz_xxxxz, g_0_z_xyyz_xxxyy, g_0_z_xyyz_xxxyz, g_0_z_xyyz_xxxzz, g_0_z_xyyz_xxyyy, g_0_z_xyyz_xxyyz, g_0_z_xyyz_xxyzz, g_0_z_xyyz_xxzzz, g_0_z_xyyz_xyyyy, g_0_z_xyyz_xyyyz, g_0_z_xyyz_xyyzz, g_0_z_xyyz_xyzzz, g_0_z_xyyz_xzzzz, g_0_z_xyyz_yyyyy, g_0_z_xyyz_yyyyz, g_0_z_xyyz_yyyzz, g_0_z_xyyz_yyzzz, g_0_z_xyyz_yzzzz, g_0_z_xyyz_zzzzz, g_0_z_yyz_xxxxx, g_0_z_yyz_xxxxxx, g_0_z_yyz_xxxxxy, g_0_z_yyz_xxxxxz, g_0_z_yyz_xxxxy, g_0_z_yyz_xxxxyy, g_0_z_yyz_xxxxyz, g_0_z_yyz_xxxxz, g_0_z_yyz_xxxxzz, g_0_z_yyz_xxxyy, g_0_z_yyz_xxxyyy, g_0_z_yyz_xxxyyz, g_0_z_yyz_xxxyz, g_0_z_yyz_xxxyzz, g_0_z_yyz_xxxzz, g_0_z_yyz_xxxzzz, g_0_z_yyz_xxyyy, g_0_z_yyz_xxyyyy, g_0_z_yyz_xxyyyz, g_0_z_yyz_xxyyz, g_0_z_yyz_xxyyzz, g_0_z_yyz_xxyzz, g_0_z_yyz_xxyzzz, g_0_z_yyz_xxzzz, g_0_z_yyz_xxzzzz, g_0_z_yyz_xyyyy, g_0_z_yyz_xyyyyy, g_0_z_yyz_xyyyyz, g_0_z_yyz_xyyyz, g_0_z_yyz_xyyyzz, g_0_z_yyz_xyyzz, g_0_z_yyz_xyyzzz, g_0_z_yyz_xyzzz, g_0_z_yyz_xyzzzz, g_0_z_yyz_xzzzz, g_0_z_yyz_xzzzzz, g_0_z_yyz_yyyyy, g_0_z_yyz_yyyyz, g_0_z_yyz_yyyzz, g_0_z_yyz_yyzzz, g_0_z_yyz_yzzzz, g_0_z_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_xxxxx[k] = -g_0_z_yyz_xxxxx[k] * ab_x + g_0_z_yyz_xxxxxx[k];

                g_0_z_xyyz_xxxxy[k] = -g_0_z_yyz_xxxxy[k] * ab_x + g_0_z_yyz_xxxxxy[k];

                g_0_z_xyyz_xxxxz[k] = -g_0_z_yyz_xxxxz[k] * ab_x + g_0_z_yyz_xxxxxz[k];

                g_0_z_xyyz_xxxyy[k] = -g_0_z_yyz_xxxyy[k] * ab_x + g_0_z_yyz_xxxxyy[k];

                g_0_z_xyyz_xxxyz[k] = -g_0_z_yyz_xxxyz[k] * ab_x + g_0_z_yyz_xxxxyz[k];

                g_0_z_xyyz_xxxzz[k] = -g_0_z_yyz_xxxzz[k] * ab_x + g_0_z_yyz_xxxxzz[k];

                g_0_z_xyyz_xxyyy[k] = -g_0_z_yyz_xxyyy[k] * ab_x + g_0_z_yyz_xxxyyy[k];

                g_0_z_xyyz_xxyyz[k] = -g_0_z_yyz_xxyyz[k] * ab_x + g_0_z_yyz_xxxyyz[k];

                g_0_z_xyyz_xxyzz[k] = -g_0_z_yyz_xxyzz[k] * ab_x + g_0_z_yyz_xxxyzz[k];

                g_0_z_xyyz_xxzzz[k] = -g_0_z_yyz_xxzzz[k] * ab_x + g_0_z_yyz_xxxzzz[k];

                g_0_z_xyyz_xyyyy[k] = -g_0_z_yyz_xyyyy[k] * ab_x + g_0_z_yyz_xxyyyy[k];

                g_0_z_xyyz_xyyyz[k] = -g_0_z_yyz_xyyyz[k] * ab_x + g_0_z_yyz_xxyyyz[k];

                g_0_z_xyyz_xyyzz[k] = -g_0_z_yyz_xyyzz[k] * ab_x + g_0_z_yyz_xxyyzz[k];

                g_0_z_xyyz_xyzzz[k] = -g_0_z_yyz_xyzzz[k] * ab_x + g_0_z_yyz_xxyzzz[k];

                g_0_z_xyyz_xzzzz[k] = -g_0_z_yyz_xzzzz[k] * ab_x + g_0_z_yyz_xxzzzz[k];

                g_0_z_xyyz_yyyyy[k] = -g_0_z_yyz_yyyyy[k] * ab_x + g_0_z_yyz_xyyyyy[k];

                g_0_z_xyyz_yyyyz[k] = -g_0_z_yyz_yyyyz[k] * ab_x + g_0_z_yyz_xyyyyz[k];

                g_0_z_xyyz_yyyzz[k] = -g_0_z_yyz_yyyzz[k] * ab_x + g_0_z_yyz_xyyyzz[k];

                g_0_z_xyyz_yyzzz[k] = -g_0_z_yyz_yyzzz[k] * ab_x + g_0_z_yyz_xyyzzz[k];

                g_0_z_xyyz_yzzzz[k] = -g_0_z_yyz_yzzzz[k] * ab_x + g_0_z_yyz_xyzzzz[k];

                g_0_z_xyyz_zzzzz[k] = -g_0_z_yyz_zzzzz[k] * ab_x + g_0_z_yyz_xzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyzz_xxxxx, g_0_z_xyzz_xxxxy, g_0_z_xyzz_xxxxz, g_0_z_xyzz_xxxyy, g_0_z_xyzz_xxxyz, g_0_z_xyzz_xxxzz, g_0_z_xyzz_xxyyy, g_0_z_xyzz_xxyyz, g_0_z_xyzz_xxyzz, g_0_z_xyzz_xxzzz, g_0_z_xyzz_xyyyy, g_0_z_xyzz_xyyyz, g_0_z_xyzz_xyyzz, g_0_z_xyzz_xyzzz, g_0_z_xyzz_xzzzz, g_0_z_xyzz_yyyyy, g_0_z_xyzz_yyyyz, g_0_z_xyzz_yyyzz, g_0_z_xyzz_yyzzz, g_0_z_xyzz_yzzzz, g_0_z_xyzz_zzzzz, g_0_z_yzz_xxxxx, g_0_z_yzz_xxxxxx, g_0_z_yzz_xxxxxy, g_0_z_yzz_xxxxxz, g_0_z_yzz_xxxxy, g_0_z_yzz_xxxxyy, g_0_z_yzz_xxxxyz, g_0_z_yzz_xxxxz, g_0_z_yzz_xxxxzz, g_0_z_yzz_xxxyy, g_0_z_yzz_xxxyyy, g_0_z_yzz_xxxyyz, g_0_z_yzz_xxxyz, g_0_z_yzz_xxxyzz, g_0_z_yzz_xxxzz, g_0_z_yzz_xxxzzz, g_0_z_yzz_xxyyy, g_0_z_yzz_xxyyyy, g_0_z_yzz_xxyyyz, g_0_z_yzz_xxyyz, g_0_z_yzz_xxyyzz, g_0_z_yzz_xxyzz, g_0_z_yzz_xxyzzz, g_0_z_yzz_xxzzz, g_0_z_yzz_xxzzzz, g_0_z_yzz_xyyyy, g_0_z_yzz_xyyyyy, g_0_z_yzz_xyyyyz, g_0_z_yzz_xyyyz, g_0_z_yzz_xyyyzz, g_0_z_yzz_xyyzz, g_0_z_yzz_xyyzzz, g_0_z_yzz_xyzzz, g_0_z_yzz_xyzzzz, g_0_z_yzz_xzzzz, g_0_z_yzz_xzzzzz, g_0_z_yzz_yyyyy, g_0_z_yzz_yyyyz, g_0_z_yzz_yyyzz, g_0_z_yzz_yyzzz, g_0_z_yzz_yzzzz, g_0_z_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_xxxxx[k] = -g_0_z_yzz_xxxxx[k] * ab_x + g_0_z_yzz_xxxxxx[k];

                g_0_z_xyzz_xxxxy[k] = -g_0_z_yzz_xxxxy[k] * ab_x + g_0_z_yzz_xxxxxy[k];

                g_0_z_xyzz_xxxxz[k] = -g_0_z_yzz_xxxxz[k] * ab_x + g_0_z_yzz_xxxxxz[k];

                g_0_z_xyzz_xxxyy[k] = -g_0_z_yzz_xxxyy[k] * ab_x + g_0_z_yzz_xxxxyy[k];

                g_0_z_xyzz_xxxyz[k] = -g_0_z_yzz_xxxyz[k] * ab_x + g_0_z_yzz_xxxxyz[k];

                g_0_z_xyzz_xxxzz[k] = -g_0_z_yzz_xxxzz[k] * ab_x + g_0_z_yzz_xxxxzz[k];

                g_0_z_xyzz_xxyyy[k] = -g_0_z_yzz_xxyyy[k] * ab_x + g_0_z_yzz_xxxyyy[k];

                g_0_z_xyzz_xxyyz[k] = -g_0_z_yzz_xxyyz[k] * ab_x + g_0_z_yzz_xxxyyz[k];

                g_0_z_xyzz_xxyzz[k] = -g_0_z_yzz_xxyzz[k] * ab_x + g_0_z_yzz_xxxyzz[k];

                g_0_z_xyzz_xxzzz[k] = -g_0_z_yzz_xxzzz[k] * ab_x + g_0_z_yzz_xxxzzz[k];

                g_0_z_xyzz_xyyyy[k] = -g_0_z_yzz_xyyyy[k] * ab_x + g_0_z_yzz_xxyyyy[k];

                g_0_z_xyzz_xyyyz[k] = -g_0_z_yzz_xyyyz[k] * ab_x + g_0_z_yzz_xxyyyz[k];

                g_0_z_xyzz_xyyzz[k] = -g_0_z_yzz_xyyzz[k] * ab_x + g_0_z_yzz_xxyyzz[k];

                g_0_z_xyzz_xyzzz[k] = -g_0_z_yzz_xyzzz[k] * ab_x + g_0_z_yzz_xxyzzz[k];

                g_0_z_xyzz_xzzzz[k] = -g_0_z_yzz_xzzzz[k] * ab_x + g_0_z_yzz_xxzzzz[k];

                g_0_z_xyzz_yyyyy[k] = -g_0_z_yzz_yyyyy[k] * ab_x + g_0_z_yzz_xyyyyy[k];

                g_0_z_xyzz_yyyyz[k] = -g_0_z_yzz_yyyyz[k] * ab_x + g_0_z_yzz_xyyyyz[k];

                g_0_z_xyzz_yyyzz[k] = -g_0_z_yzz_yyyzz[k] * ab_x + g_0_z_yzz_xyyyzz[k];

                g_0_z_xyzz_yyzzz[k] = -g_0_z_yzz_yyzzz[k] * ab_x + g_0_z_yzz_xyyzzz[k];

                g_0_z_xyzz_yzzzz[k] = -g_0_z_yzz_yzzzz[k] * ab_x + g_0_z_yzz_xyzzzz[k];

                g_0_z_xyzz_zzzzz[k] = -g_0_z_yzz_zzzzz[k] * ab_x + g_0_z_yzz_xzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xzzz_xxxxx, g_0_z_xzzz_xxxxy, g_0_z_xzzz_xxxxz, g_0_z_xzzz_xxxyy, g_0_z_xzzz_xxxyz, g_0_z_xzzz_xxxzz, g_0_z_xzzz_xxyyy, g_0_z_xzzz_xxyyz, g_0_z_xzzz_xxyzz, g_0_z_xzzz_xxzzz, g_0_z_xzzz_xyyyy, g_0_z_xzzz_xyyyz, g_0_z_xzzz_xyyzz, g_0_z_xzzz_xyzzz, g_0_z_xzzz_xzzzz, g_0_z_xzzz_yyyyy, g_0_z_xzzz_yyyyz, g_0_z_xzzz_yyyzz, g_0_z_xzzz_yyzzz, g_0_z_xzzz_yzzzz, g_0_z_xzzz_zzzzz, g_0_z_zzz_xxxxx, g_0_z_zzz_xxxxxx, g_0_z_zzz_xxxxxy, g_0_z_zzz_xxxxxz, g_0_z_zzz_xxxxy, g_0_z_zzz_xxxxyy, g_0_z_zzz_xxxxyz, g_0_z_zzz_xxxxz, g_0_z_zzz_xxxxzz, g_0_z_zzz_xxxyy, g_0_z_zzz_xxxyyy, g_0_z_zzz_xxxyyz, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxyzz, g_0_z_zzz_xxxzz, g_0_z_zzz_xxxzzz, g_0_z_zzz_xxyyy, g_0_z_zzz_xxyyyy, g_0_z_zzz_xxyyyz, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyyzz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxyzzz, g_0_z_zzz_xxzzz, g_0_z_zzz_xxzzzz, g_0_z_zzz_xyyyy, g_0_z_zzz_xyyyyy, g_0_z_zzz_xyyyyz, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyyzz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyyzzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xyzzzz, g_0_z_zzz_xzzzz, g_0_z_zzz_xzzzzz, g_0_z_zzz_yyyyy, g_0_z_zzz_yyyyz, g_0_z_zzz_yyyzz, g_0_z_zzz_yyzzz, g_0_z_zzz_yzzzz, g_0_z_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_xxxxx[k] = -g_0_z_zzz_xxxxx[k] * ab_x + g_0_z_zzz_xxxxxx[k];

                g_0_z_xzzz_xxxxy[k] = -g_0_z_zzz_xxxxy[k] * ab_x + g_0_z_zzz_xxxxxy[k];

                g_0_z_xzzz_xxxxz[k] = -g_0_z_zzz_xxxxz[k] * ab_x + g_0_z_zzz_xxxxxz[k];

                g_0_z_xzzz_xxxyy[k] = -g_0_z_zzz_xxxyy[k] * ab_x + g_0_z_zzz_xxxxyy[k];

                g_0_z_xzzz_xxxyz[k] = -g_0_z_zzz_xxxyz[k] * ab_x + g_0_z_zzz_xxxxyz[k];

                g_0_z_xzzz_xxxzz[k] = -g_0_z_zzz_xxxzz[k] * ab_x + g_0_z_zzz_xxxxzz[k];

                g_0_z_xzzz_xxyyy[k] = -g_0_z_zzz_xxyyy[k] * ab_x + g_0_z_zzz_xxxyyy[k];

                g_0_z_xzzz_xxyyz[k] = -g_0_z_zzz_xxyyz[k] * ab_x + g_0_z_zzz_xxxyyz[k];

                g_0_z_xzzz_xxyzz[k] = -g_0_z_zzz_xxyzz[k] * ab_x + g_0_z_zzz_xxxyzz[k];

                g_0_z_xzzz_xxzzz[k] = -g_0_z_zzz_xxzzz[k] * ab_x + g_0_z_zzz_xxxzzz[k];

                g_0_z_xzzz_xyyyy[k] = -g_0_z_zzz_xyyyy[k] * ab_x + g_0_z_zzz_xxyyyy[k];

                g_0_z_xzzz_xyyyz[k] = -g_0_z_zzz_xyyyz[k] * ab_x + g_0_z_zzz_xxyyyz[k];

                g_0_z_xzzz_xyyzz[k] = -g_0_z_zzz_xyyzz[k] * ab_x + g_0_z_zzz_xxyyzz[k];

                g_0_z_xzzz_xyzzz[k] = -g_0_z_zzz_xyzzz[k] * ab_x + g_0_z_zzz_xxyzzz[k];

                g_0_z_xzzz_xzzzz[k] = -g_0_z_zzz_xzzzz[k] * ab_x + g_0_z_zzz_xxzzzz[k];

                g_0_z_xzzz_yyyyy[k] = -g_0_z_zzz_yyyyy[k] * ab_x + g_0_z_zzz_xyyyyy[k];

                g_0_z_xzzz_yyyyz[k] = -g_0_z_zzz_yyyyz[k] * ab_x + g_0_z_zzz_xyyyyz[k];

                g_0_z_xzzz_yyyzz[k] = -g_0_z_zzz_yyyzz[k] * ab_x + g_0_z_zzz_xyyyzz[k];

                g_0_z_xzzz_yyzzz[k] = -g_0_z_zzz_yyzzz[k] * ab_x + g_0_z_zzz_xyyzzz[k];

                g_0_z_xzzz_yzzzz[k] = -g_0_z_zzz_yzzzz[k] * ab_x + g_0_z_zzz_xyzzzz[k];

                g_0_z_xzzz_zzzzz[k] = -g_0_z_zzz_zzzzz[k] * ab_x + g_0_z_zzz_xzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyy_xxxxx, g_0_z_yyy_xxxxxy, g_0_z_yyy_xxxxy, g_0_z_yyy_xxxxyy, g_0_z_yyy_xxxxyz, g_0_z_yyy_xxxxz, g_0_z_yyy_xxxyy, g_0_z_yyy_xxxyyy, g_0_z_yyy_xxxyyz, g_0_z_yyy_xxxyz, g_0_z_yyy_xxxyzz, g_0_z_yyy_xxxzz, g_0_z_yyy_xxyyy, g_0_z_yyy_xxyyyy, g_0_z_yyy_xxyyyz, g_0_z_yyy_xxyyz, g_0_z_yyy_xxyyzz, g_0_z_yyy_xxyzz, g_0_z_yyy_xxyzzz, g_0_z_yyy_xxzzz, g_0_z_yyy_xyyyy, g_0_z_yyy_xyyyyy, g_0_z_yyy_xyyyyz, g_0_z_yyy_xyyyz, g_0_z_yyy_xyyyzz, g_0_z_yyy_xyyzz, g_0_z_yyy_xyyzzz, g_0_z_yyy_xyzzz, g_0_z_yyy_xyzzzz, g_0_z_yyy_xzzzz, g_0_z_yyy_yyyyy, g_0_z_yyy_yyyyyy, g_0_z_yyy_yyyyyz, g_0_z_yyy_yyyyz, g_0_z_yyy_yyyyzz, g_0_z_yyy_yyyzz, g_0_z_yyy_yyyzzz, g_0_z_yyy_yyzzz, g_0_z_yyy_yyzzzz, g_0_z_yyy_yzzzz, g_0_z_yyy_yzzzzz, g_0_z_yyy_zzzzz, g_0_z_yyyy_xxxxx, g_0_z_yyyy_xxxxy, g_0_z_yyyy_xxxxz, g_0_z_yyyy_xxxyy, g_0_z_yyyy_xxxyz, g_0_z_yyyy_xxxzz, g_0_z_yyyy_xxyyy, g_0_z_yyyy_xxyyz, g_0_z_yyyy_xxyzz, g_0_z_yyyy_xxzzz, g_0_z_yyyy_xyyyy, g_0_z_yyyy_xyyyz, g_0_z_yyyy_xyyzz, g_0_z_yyyy_xyzzz, g_0_z_yyyy_xzzzz, g_0_z_yyyy_yyyyy, g_0_z_yyyy_yyyyz, g_0_z_yyyy_yyyzz, g_0_z_yyyy_yyzzz, g_0_z_yyyy_yzzzz, g_0_z_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_xxxxx[k] = -g_0_z_yyy_xxxxx[k] * ab_y + g_0_z_yyy_xxxxxy[k];

                g_0_z_yyyy_xxxxy[k] = -g_0_z_yyy_xxxxy[k] * ab_y + g_0_z_yyy_xxxxyy[k];

                g_0_z_yyyy_xxxxz[k] = -g_0_z_yyy_xxxxz[k] * ab_y + g_0_z_yyy_xxxxyz[k];

                g_0_z_yyyy_xxxyy[k] = -g_0_z_yyy_xxxyy[k] * ab_y + g_0_z_yyy_xxxyyy[k];

                g_0_z_yyyy_xxxyz[k] = -g_0_z_yyy_xxxyz[k] * ab_y + g_0_z_yyy_xxxyyz[k];

                g_0_z_yyyy_xxxzz[k] = -g_0_z_yyy_xxxzz[k] * ab_y + g_0_z_yyy_xxxyzz[k];

                g_0_z_yyyy_xxyyy[k] = -g_0_z_yyy_xxyyy[k] * ab_y + g_0_z_yyy_xxyyyy[k];

                g_0_z_yyyy_xxyyz[k] = -g_0_z_yyy_xxyyz[k] * ab_y + g_0_z_yyy_xxyyyz[k];

                g_0_z_yyyy_xxyzz[k] = -g_0_z_yyy_xxyzz[k] * ab_y + g_0_z_yyy_xxyyzz[k];

                g_0_z_yyyy_xxzzz[k] = -g_0_z_yyy_xxzzz[k] * ab_y + g_0_z_yyy_xxyzzz[k];

                g_0_z_yyyy_xyyyy[k] = -g_0_z_yyy_xyyyy[k] * ab_y + g_0_z_yyy_xyyyyy[k];

                g_0_z_yyyy_xyyyz[k] = -g_0_z_yyy_xyyyz[k] * ab_y + g_0_z_yyy_xyyyyz[k];

                g_0_z_yyyy_xyyzz[k] = -g_0_z_yyy_xyyzz[k] * ab_y + g_0_z_yyy_xyyyzz[k];

                g_0_z_yyyy_xyzzz[k] = -g_0_z_yyy_xyzzz[k] * ab_y + g_0_z_yyy_xyyzzz[k];

                g_0_z_yyyy_xzzzz[k] = -g_0_z_yyy_xzzzz[k] * ab_y + g_0_z_yyy_xyzzzz[k];

                g_0_z_yyyy_yyyyy[k] = -g_0_z_yyy_yyyyy[k] * ab_y + g_0_z_yyy_yyyyyy[k];

                g_0_z_yyyy_yyyyz[k] = -g_0_z_yyy_yyyyz[k] * ab_y + g_0_z_yyy_yyyyyz[k];

                g_0_z_yyyy_yyyzz[k] = -g_0_z_yyy_yyyzz[k] * ab_y + g_0_z_yyy_yyyyzz[k];

                g_0_z_yyyy_yyzzz[k] = -g_0_z_yyy_yyzzz[k] * ab_y + g_0_z_yyy_yyyzzz[k];

                g_0_z_yyyy_yzzzz[k] = -g_0_z_yyy_yzzzz[k] * ab_y + g_0_z_yyy_yyzzzz[k];

                g_0_z_yyyy_zzzzz[k] = -g_0_z_yyy_zzzzz[k] * ab_y + g_0_z_yyy_yzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyyz_xxxxx, g_0_z_yyyz_xxxxy, g_0_z_yyyz_xxxxz, g_0_z_yyyz_xxxyy, g_0_z_yyyz_xxxyz, g_0_z_yyyz_xxxzz, g_0_z_yyyz_xxyyy, g_0_z_yyyz_xxyyz, g_0_z_yyyz_xxyzz, g_0_z_yyyz_xxzzz, g_0_z_yyyz_xyyyy, g_0_z_yyyz_xyyyz, g_0_z_yyyz_xyyzz, g_0_z_yyyz_xyzzz, g_0_z_yyyz_xzzzz, g_0_z_yyyz_yyyyy, g_0_z_yyyz_yyyyz, g_0_z_yyyz_yyyzz, g_0_z_yyyz_yyzzz, g_0_z_yyyz_yzzzz, g_0_z_yyyz_zzzzz, g_0_z_yyz_xxxxx, g_0_z_yyz_xxxxxy, g_0_z_yyz_xxxxy, g_0_z_yyz_xxxxyy, g_0_z_yyz_xxxxyz, g_0_z_yyz_xxxxz, g_0_z_yyz_xxxyy, g_0_z_yyz_xxxyyy, g_0_z_yyz_xxxyyz, g_0_z_yyz_xxxyz, g_0_z_yyz_xxxyzz, g_0_z_yyz_xxxzz, g_0_z_yyz_xxyyy, g_0_z_yyz_xxyyyy, g_0_z_yyz_xxyyyz, g_0_z_yyz_xxyyz, g_0_z_yyz_xxyyzz, g_0_z_yyz_xxyzz, g_0_z_yyz_xxyzzz, g_0_z_yyz_xxzzz, g_0_z_yyz_xyyyy, g_0_z_yyz_xyyyyy, g_0_z_yyz_xyyyyz, g_0_z_yyz_xyyyz, g_0_z_yyz_xyyyzz, g_0_z_yyz_xyyzz, g_0_z_yyz_xyyzzz, g_0_z_yyz_xyzzz, g_0_z_yyz_xyzzzz, g_0_z_yyz_xzzzz, g_0_z_yyz_yyyyy, g_0_z_yyz_yyyyyy, g_0_z_yyz_yyyyyz, g_0_z_yyz_yyyyz, g_0_z_yyz_yyyyzz, g_0_z_yyz_yyyzz, g_0_z_yyz_yyyzzz, g_0_z_yyz_yyzzz, g_0_z_yyz_yyzzzz, g_0_z_yyz_yzzzz, g_0_z_yyz_yzzzzz, g_0_z_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_xxxxx[k] = -g_0_z_yyz_xxxxx[k] * ab_y + g_0_z_yyz_xxxxxy[k];

                g_0_z_yyyz_xxxxy[k] = -g_0_z_yyz_xxxxy[k] * ab_y + g_0_z_yyz_xxxxyy[k];

                g_0_z_yyyz_xxxxz[k] = -g_0_z_yyz_xxxxz[k] * ab_y + g_0_z_yyz_xxxxyz[k];

                g_0_z_yyyz_xxxyy[k] = -g_0_z_yyz_xxxyy[k] * ab_y + g_0_z_yyz_xxxyyy[k];

                g_0_z_yyyz_xxxyz[k] = -g_0_z_yyz_xxxyz[k] * ab_y + g_0_z_yyz_xxxyyz[k];

                g_0_z_yyyz_xxxzz[k] = -g_0_z_yyz_xxxzz[k] * ab_y + g_0_z_yyz_xxxyzz[k];

                g_0_z_yyyz_xxyyy[k] = -g_0_z_yyz_xxyyy[k] * ab_y + g_0_z_yyz_xxyyyy[k];

                g_0_z_yyyz_xxyyz[k] = -g_0_z_yyz_xxyyz[k] * ab_y + g_0_z_yyz_xxyyyz[k];

                g_0_z_yyyz_xxyzz[k] = -g_0_z_yyz_xxyzz[k] * ab_y + g_0_z_yyz_xxyyzz[k];

                g_0_z_yyyz_xxzzz[k] = -g_0_z_yyz_xxzzz[k] * ab_y + g_0_z_yyz_xxyzzz[k];

                g_0_z_yyyz_xyyyy[k] = -g_0_z_yyz_xyyyy[k] * ab_y + g_0_z_yyz_xyyyyy[k];

                g_0_z_yyyz_xyyyz[k] = -g_0_z_yyz_xyyyz[k] * ab_y + g_0_z_yyz_xyyyyz[k];

                g_0_z_yyyz_xyyzz[k] = -g_0_z_yyz_xyyzz[k] * ab_y + g_0_z_yyz_xyyyzz[k];

                g_0_z_yyyz_xyzzz[k] = -g_0_z_yyz_xyzzz[k] * ab_y + g_0_z_yyz_xyyzzz[k];

                g_0_z_yyyz_xzzzz[k] = -g_0_z_yyz_xzzzz[k] * ab_y + g_0_z_yyz_xyzzzz[k];

                g_0_z_yyyz_yyyyy[k] = -g_0_z_yyz_yyyyy[k] * ab_y + g_0_z_yyz_yyyyyy[k];

                g_0_z_yyyz_yyyyz[k] = -g_0_z_yyz_yyyyz[k] * ab_y + g_0_z_yyz_yyyyyz[k];

                g_0_z_yyyz_yyyzz[k] = -g_0_z_yyz_yyyzz[k] * ab_y + g_0_z_yyz_yyyyzz[k];

                g_0_z_yyyz_yyzzz[k] = -g_0_z_yyz_yyzzz[k] * ab_y + g_0_z_yyz_yyyzzz[k];

                g_0_z_yyyz_yzzzz[k] = -g_0_z_yyz_yzzzz[k] * ab_y + g_0_z_yyz_yyzzzz[k];

                g_0_z_yyyz_zzzzz[k] = -g_0_z_yyz_zzzzz[k] * ab_y + g_0_z_yyz_yzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyzz_xxxxx, g_0_z_yyzz_xxxxy, g_0_z_yyzz_xxxxz, g_0_z_yyzz_xxxyy, g_0_z_yyzz_xxxyz, g_0_z_yyzz_xxxzz, g_0_z_yyzz_xxyyy, g_0_z_yyzz_xxyyz, g_0_z_yyzz_xxyzz, g_0_z_yyzz_xxzzz, g_0_z_yyzz_xyyyy, g_0_z_yyzz_xyyyz, g_0_z_yyzz_xyyzz, g_0_z_yyzz_xyzzz, g_0_z_yyzz_xzzzz, g_0_z_yyzz_yyyyy, g_0_z_yyzz_yyyyz, g_0_z_yyzz_yyyzz, g_0_z_yyzz_yyzzz, g_0_z_yyzz_yzzzz, g_0_z_yyzz_zzzzz, g_0_z_yzz_xxxxx, g_0_z_yzz_xxxxxy, g_0_z_yzz_xxxxy, g_0_z_yzz_xxxxyy, g_0_z_yzz_xxxxyz, g_0_z_yzz_xxxxz, g_0_z_yzz_xxxyy, g_0_z_yzz_xxxyyy, g_0_z_yzz_xxxyyz, g_0_z_yzz_xxxyz, g_0_z_yzz_xxxyzz, g_0_z_yzz_xxxzz, g_0_z_yzz_xxyyy, g_0_z_yzz_xxyyyy, g_0_z_yzz_xxyyyz, g_0_z_yzz_xxyyz, g_0_z_yzz_xxyyzz, g_0_z_yzz_xxyzz, g_0_z_yzz_xxyzzz, g_0_z_yzz_xxzzz, g_0_z_yzz_xyyyy, g_0_z_yzz_xyyyyy, g_0_z_yzz_xyyyyz, g_0_z_yzz_xyyyz, g_0_z_yzz_xyyyzz, g_0_z_yzz_xyyzz, g_0_z_yzz_xyyzzz, g_0_z_yzz_xyzzz, g_0_z_yzz_xyzzzz, g_0_z_yzz_xzzzz, g_0_z_yzz_yyyyy, g_0_z_yzz_yyyyyy, g_0_z_yzz_yyyyyz, g_0_z_yzz_yyyyz, g_0_z_yzz_yyyyzz, g_0_z_yzz_yyyzz, g_0_z_yzz_yyyzzz, g_0_z_yzz_yyzzz, g_0_z_yzz_yyzzzz, g_0_z_yzz_yzzzz, g_0_z_yzz_yzzzzz, g_0_z_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_xxxxx[k] = -g_0_z_yzz_xxxxx[k] * ab_y + g_0_z_yzz_xxxxxy[k];

                g_0_z_yyzz_xxxxy[k] = -g_0_z_yzz_xxxxy[k] * ab_y + g_0_z_yzz_xxxxyy[k];

                g_0_z_yyzz_xxxxz[k] = -g_0_z_yzz_xxxxz[k] * ab_y + g_0_z_yzz_xxxxyz[k];

                g_0_z_yyzz_xxxyy[k] = -g_0_z_yzz_xxxyy[k] * ab_y + g_0_z_yzz_xxxyyy[k];

                g_0_z_yyzz_xxxyz[k] = -g_0_z_yzz_xxxyz[k] * ab_y + g_0_z_yzz_xxxyyz[k];

                g_0_z_yyzz_xxxzz[k] = -g_0_z_yzz_xxxzz[k] * ab_y + g_0_z_yzz_xxxyzz[k];

                g_0_z_yyzz_xxyyy[k] = -g_0_z_yzz_xxyyy[k] * ab_y + g_0_z_yzz_xxyyyy[k];

                g_0_z_yyzz_xxyyz[k] = -g_0_z_yzz_xxyyz[k] * ab_y + g_0_z_yzz_xxyyyz[k];

                g_0_z_yyzz_xxyzz[k] = -g_0_z_yzz_xxyzz[k] * ab_y + g_0_z_yzz_xxyyzz[k];

                g_0_z_yyzz_xxzzz[k] = -g_0_z_yzz_xxzzz[k] * ab_y + g_0_z_yzz_xxyzzz[k];

                g_0_z_yyzz_xyyyy[k] = -g_0_z_yzz_xyyyy[k] * ab_y + g_0_z_yzz_xyyyyy[k];

                g_0_z_yyzz_xyyyz[k] = -g_0_z_yzz_xyyyz[k] * ab_y + g_0_z_yzz_xyyyyz[k];

                g_0_z_yyzz_xyyzz[k] = -g_0_z_yzz_xyyzz[k] * ab_y + g_0_z_yzz_xyyyzz[k];

                g_0_z_yyzz_xyzzz[k] = -g_0_z_yzz_xyzzz[k] * ab_y + g_0_z_yzz_xyyzzz[k];

                g_0_z_yyzz_xzzzz[k] = -g_0_z_yzz_xzzzz[k] * ab_y + g_0_z_yzz_xyzzzz[k];

                g_0_z_yyzz_yyyyy[k] = -g_0_z_yzz_yyyyy[k] * ab_y + g_0_z_yzz_yyyyyy[k];

                g_0_z_yyzz_yyyyz[k] = -g_0_z_yzz_yyyyz[k] * ab_y + g_0_z_yzz_yyyyyz[k];

                g_0_z_yyzz_yyyzz[k] = -g_0_z_yzz_yyyzz[k] * ab_y + g_0_z_yzz_yyyyzz[k];

                g_0_z_yyzz_yyzzz[k] = -g_0_z_yzz_yyzzz[k] * ab_y + g_0_z_yzz_yyyzzz[k];

                g_0_z_yyzz_yzzzz[k] = -g_0_z_yzz_yzzzz[k] * ab_y + g_0_z_yzz_yyzzzz[k];

                g_0_z_yyzz_zzzzz[k] = -g_0_z_yzz_zzzzz[k] * ab_y + g_0_z_yzz_yzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yzzz_xxxxx, g_0_z_yzzz_xxxxy, g_0_z_yzzz_xxxxz, g_0_z_yzzz_xxxyy, g_0_z_yzzz_xxxyz, g_0_z_yzzz_xxxzz, g_0_z_yzzz_xxyyy, g_0_z_yzzz_xxyyz, g_0_z_yzzz_xxyzz, g_0_z_yzzz_xxzzz, g_0_z_yzzz_xyyyy, g_0_z_yzzz_xyyyz, g_0_z_yzzz_xyyzz, g_0_z_yzzz_xyzzz, g_0_z_yzzz_xzzzz, g_0_z_yzzz_yyyyy, g_0_z_yzzz_yyyyz, g_0_z_yzzz_yyyzz, g_0_z_yzzz_yyzzz, g_0_z_yzzz_yzzzz, g_0_z_yzzz_zzzzz, g_0_z_zzz_xxxxx, g_0_z_zzz_xxxxxy, g_0_z_zzz_xxxxy, g_0_z_zzz_xxxxyy, g_0_z_zzz_xxxxyz, g_0_z_zzz_xxxxz, g_0_z_zzz_xxxyy, g_0_z_zzz_xxxyyy, g_0_z_zzz_xxxyyz, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxyzz, g_0_z_zzz_xxxzz, g_0_z_zzz_xxyyy, g_0_z_zzz_xxyyyy, g_0_z_zzz_xxyyyz, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyyzz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxyzzz, g_0_z_zzz_xxzzz, g_0_z_zzz_xyyyy, g_0_z_zzz_xyyyyy, g_0_z_zzz_xyyyyz, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyyzz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyyzzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xyzzzz, g_0_z_zzz_xzzzz, g_0_z_zzz_yyyyy, g_0_z_zzz_yyyyyy, g_0_z_zzz_yyyyyz, g_0_z_zzz_yyyyz, g_0_z_zzz_yyyyzz, g_0_z_zzz_yyyzz, g_0_z_zzz_yyyzzz, g_0_z_zzz_yyzzz, g_0_z_zzz_yyzzzz, g_0_z_zzz_yzzzz, g_0_z_zzz_yzzzzz, g_0_z_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_xxxxx[k] = -g_0_z_zzz_xxxxx[k] * ab_y + g_0_z_zzz_xxxxxy[k];

                g_0_z_yzzz_xxxxy[k] = -g_0_z_zzz_xxxxy[k] * ab_y + g_0_z_zzz_xxxxyy[k];

                g_0_z_yzzz_xxxxz[k] = -g_0_z_zzz_xxxxz[k] * ab_y + g_0_z_zzz_xxxxyz[k];

                g_0_z_yzzz_xxxyy[k] = -g_0_z_zzz_xxxyy[k] * ab_y + g_0_z_zzz_xxxyyy[k];

                g_0_z_yzzz_xxxyz[k] = -g_0_z_zzz_xxxyz[k] * ab_y + g_0_z_zzz_xxxyyz[k];

                g_0_z_yzzz_xxxzz[k] = -g_0_z_zzz_xxxzz[k] * ab_y + g_0_z_zzz_xxxyzz[k];

                g_0_z_yzzz_xxyyy[k] = -g_0_z_zzz_xxyyy[k] * ab_y + g_0_z_zzz_xxyyyy[k];

                g_0_z_yzzz_xxyyz[k] = -g_0_z_zzz_xxyyz[k] * ab_y + g_0_z_zzz_xxyyyz[k];

                g_0_z_yzzz_xxyzz[k] = -g_0_z_zzz_xxyzz[k] * ab_y + g_0_z_zzz_xxyyzz[k];

                g_0_z_yzzz_xxzzz[k] = -g_0_z_zzz_xxzzz[k] * ab_y + g_0_z_zzz_xxyzzz[k];

                g_0_z_yzzz_xyyyy[k] = -g_0_z_zzz_xyyyy[k] * ab_y + g_0_z_zzz_xyyyyy[k];

                g_0_z_yzzz_xyyyz[k] = -g_0_z_zzz_xyyyz[k] * ab_y + g_0_z_zzz_xyyyyz[k];

                g_0_z_yzzz_xyyzz[k] = -g_0_z_zzz_xyyzz[k] * ab_y + g_0_z_zzz_xyyyzz[k];

                g_0_z_yzzz_xyzzz[k] = -g_0_z_zzz_xyzzz[k] * ab_y + g_0_z_zzz_xyyzzz[k];

                g_0_z_yzzz_xzzzz[k] = -g_0_z_zzz_xzzzz[k] * ab_y + g_0_z_zzz_xyzzzz[k];

                g_0_z_yzzz_yyyyy[k] = -g_0_z_zzz_yyyyy[k] * ab_y + g_0_z_zzz_yyyyyy[k];

                g_0_z_yzzz_yyyyz[k] = -g_0_z_zzz_yyyyz[k] * ab_y + g_0_z_zzz_yyyyyz[k];

                g_0_z_yzzz_yyyzz[k] = -g_0_z_zzz_yyyzz[k] * ab_y + g_0_z_zzz_yyyyzz[k];

                g_0_z_yzzz_yyzzz[k] = -g_0_z_zzz_yyzzz[k] * ab_y + g_0_z_zzz_yyyzzz[k];

                g_0_z_yzzz_yzzzz[k] = -g_0_z_zzz_yzzzz[k] * ab_y + g_0_z_zzz_yyzzzz[k];

                g_0_z_yzzz_zzzzz[k] = -g_0_z_zzz_zzzzz[k] * ab_y + g_0_z_zzz_yzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_zzz_xxxxx, g_0_z_zzz_xxxxxz, g_0_z_zzz_xxxxy, g_0_z_zzz_xxxxyz, g_0_z_zzz_xxxxz, g_0_z_zzz_xxxxzz, g_0_z_zzz_xxxyy, g_0_z_zzz_xxxyyz, g_0_z_zzz_xxxyz, g_0_z_zzz_xxxyzz, g_0_z_zzz_xxxzz, g_0_z_zzz_xxxzzz, g_0_z_zzz_xxyyy, g_0_z_zzz_xxyyyz, g_0_z_zzz_xxyyz, g_0_z_zzz_xxyyzz, g_0_z_zzz_xxyzz, g_0_z_zzz_xxyzzz, g_0_z_zzz_xxzzz, g_0_z_zzz_xxzzzz, g_0_z_zzz_xyyyy, g_0_z_zzz_xyyyyz, g_0_z_zzz_xyyyz, g_0_z_zzz_xyyyzz, g_0_z_zzz_xyyzz, g_0_z_zzz_xyyzzz, g_0_z_zzz_xyzzz, g_0_z_zzz_xyzzzz, g_0_z_zzz_xzzzz, g_0_z_zzz_xzzzzz, g_0_z_zzz_yyyyy, g_0_z_zzz_yyyyyz, g_0_z_zzz_yyyyz, g_0_z_zzz_yyyyzz, g_0_z_zzz_yyyzz, g_0_z_zzz_yyyzzz, g_0_z_zzz_yyzzz, g_0_z_zzz_yyzzzz, g_0_z_zzz_yzzzz, g_0_z_zzz_yzzzzz, g_0_z_zzz_zzzzz, g_0_z_zzz_zzzzzz, g_0_z_zzzz_xxxxx, g_0_z_zzzz_xxxxy, g_0_z_zzzz_xxxxz, g_0_z_zzzz_xxxyy, g_0_z_zzzz_xxxyz, g_0_z_zzzz_xxxzz, g_0_z_zzzz_xxyyy, g_0_z_zzzz_xxyyz, g_0_z_zzzz_xxyzz, g_0_z_zzzz_xxzzz, g_0_z_zzzz_xyyyy, g_0_z_zzzz_xyyyz, g_0_z_zzzz_xyyzz, g_0_z_zzzz_xyzzz, g_0_z_zzzz_xzzzz, g_0_z_zzzz_yyyyy, g_0_z_zzzz_yyyyz, g_0_z_zzzz_yyyzz, g_0_z_zzzz_yyzzz, g_0_z_zzzz_yzzzz, g_0_z_zzzz_zzzzz, g_zzz_xxxxx, g_zzz_xxxxy, g_zzz_xxxxz, g_zzz_xxxyy, g_zzz_xxxyz, g_zzz_xxxzz, g_zzz_xxyyy, g_zzz_xxyyz, g_zzz_xxyzz, g_zzz_xxzzz, g_zzz_xyyyy, g_zzz_xyyyz, g_zzz_xyyzz, g_zzz_xyzzz, g_zzz_xzzzz, g_zzz_yyyyy, g_zzz_yyyyz, g_zzz_yyyzz, g_zzz_yyzzz, g_zzz_yzzzz, g_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_xxxxx[k] = g_zzz_xxxxx[k] - g_0_z_zzz_xxxxx[k] * ab_z + g_0_z_zzz_xxxxxz[k];

                g_0_z_zzzz_xxxxy[k] = g_zzz_xxxxy[k] - g_0_z_zzz_xxxxy[k] * ab_z + g_0_z_zzz_xxxxyz[k];

                g_0_z_zzzz_xxxxz[k] = g_zzz_xxxxz[k] - g_0_z_zzz_xxxxz[k] * ab_z + g_0_z_zzz_xxxxzz[k];

                g_0_z_zzzz_xxxyy[k] = g_zzz_xxxyy[k] - g_0_z_zzz_xxxyy[k] * ab_z + g_0_z_zzz_xxxyyz[k];

                g_0_z_zzzz_xxxyz[k] = g_zzz_xxxyz[k] - g_0_z_zzz_xxxyz[k] * ab_z + g_0_z_zzz_xxxyzz[k];

                g_0_z_zzzz_xxxzz[k] = g_zzz_xxxzz[k] - g_0_z_zzz_xxxzz[k] * ab_z + g_0_z_zzz_xxxzzz[k];

                g_0_z_zzzz_xxyyy[k] = g_zzz_xxyyy[k] - g_0_z_zzz_xxyyy[k] * ab_z + g_0_z_zzz_xxyyyz[k];

                g_0_z_zzzz_xxyyz[k] = g_zzz_xxyyz[k] - g_0_z_zzz_xxyyz[k] * ab_z + g_0_z_zzz_xxyyzz[k];

                g_0_z_zzzz_xxyzz[k] = g_zzz_xxyzz[k] - g_0_z_zzz_xxyzz[k] * ab_z + g_0_z_zzz_xxyzzz[k];

                g_0_z_zzzz_xxzzz[k] = g_zzz_xxzzz[k] - g_0_z_zzz_xxzzz[k] * ab_z + g_0_z_zzz_xxzzzz[k];

                g_0_z_zzzz_xyyyy[k] = g_zzz_xyyyy[k] - g_0_z_zzz_xyyyy[k] * ab_z + g_0_z_zzz_xyyyyz[k];

                g_0_z_zzzz_xyyyz[k] = g_zzz_xyyyz[k] - g_0_z_zzz_xyyyz[k] * ab_z + g_0_z_zzz_xyyyzz[k];

                g_0_z_zzzz_xyyzz[k] = g_zzz_xyyzz[k] - g_0_z_zzz_xyyzz[k] * ab_z + g_0_z_zzz_xyyzzz[k];

                g_0_z_zzzz_xyzzz[k] = g_zzz_xyzzz[k] - g_0_z_zzz_xyzzz[k] * ab_z + g_0_z_zzz_xyzzzz[k];

                g_0_z_zzzz_xzzzz[k] = g_zzz_xzzzz[k] - g_0_z_zzz_xzzzz[k] * ab_z + g_0_z_zzz_xzzzzz[k];

                g_0_z_zzzz_yyyyy[k] = g_zzz_yyyyy[k] - g_0_z_zzz_yyyyy[k] * ab_z + g_0_z_zzz_yyyyyz[k];

                g_0_z_zzzz_yyyyz[k] = g_zzz_yyyyz[k] - g_0_z_zzz_yyyyz[k] * ab_z + g_0_z_zzz_yyyyzz[k];

                g_0_z_zzzz_yyyzz[k] = g_zzz_yyyzz[k] - g_0_z_zzz_yyyzz[k] * ab_z + g_0_z_zzz_yyyzzz[k];

                g_0_z_zzzz_yyzzz[k] = g_zzz_yyzzz[k] - g_0_z_zzz_yyzzz[k] * ab_z + g_0_z_zzz_yyzzzz[k];

                g_0_z_zzzz_yzzzz[k] = g_zzz_yzzzz[k] - g_0_z_zzz_yzzzz[k] * ab_z + g_0_z_zzz_yzzzzz[k];

                g_0_z_zzzz_zzzzz[k] = g_zzz_zzzzz[k] - g_0_z_zzz_zzzzz[k] * ab_z + g_0_z_zzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

