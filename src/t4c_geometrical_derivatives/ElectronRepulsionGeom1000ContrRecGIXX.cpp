#include "ElectronRepulsionGeom1000ContrRecGIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_gixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gixx,
                                            const size_t idx_fixx,
                                            const size_t idx_geom_10_fixx,
                                            const size_t idx_geom_10_fkxx,
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
            /// Set up components of auxilary buffer : FISS

            const auto fi_off = idx_fixx + i * dcomps + j;

            auto g_xxx_xxxxxx = cbuffer.data(fi_off + 0 * ccomps * dcomps);

            auto g_xxx_xxxxxy = cbuffer.data(fi_off + 1 * ccomps * dcomps);

            auto g_xxx_xxxxxz = cbuffer.data(fi_off + 2 * ccomps * dcomps);

            auto g_xxx_xxxxyy = cbuffer.data(fi_off + 3 * ccomps * dcomps);

            auto g_xxx_xxxxyz = cbuffer.data(fi_off + 4 * ccomps * dcomps);

            auto g_xxx_xxxxzz = cbuffer.data(fi_off + 5 * ccomps * dcomps);

            auto g_xxx_xxxyyy = cbuffer.data(fi_off + 6 * ccomps * dcomps);

            auto g_xxx_xxxyyz = cbuffer.data(fi_off + 7 * ccomps * dcomps);

            auto g_xxx_xxxyzz = cbuffer.data(fi_off + 8 * ccomps * dcomps);

            auto g_xxx_xxxzzz = cbuffer.data(fi_off + 9 * ccomps * dcomps);

            auto g_xxx_xxyyyy = cbuffer.data(fi_off + 10 * ccomps * dcomps);

            auto g_xxx_xxyyyz = cbuffer.data(fi_off + 11 * ccomps * dcomps);

            auto g_xxx_xxyyzz = cbuffer.data(fi_off + 12 * ccomps * dcomps);

            auto g_xxx_xxyzzz = cbuffer.data(fi_off + 13 * ccomps * dcomps);

            auto g_xxx_xxzzzz = cbuffer.data(fi_off + 14 * ccomps * dcomps);

            auto g_xxx_xyyyyy = cbuffer.data(fi_off + 15 * ccomps * dcomps);

            auto g_xxx_xyyyyz = cbuffer.data(fi_off + 16 * ccomps * dcomps);

            auto g_xxx_xyyyzz = cbuffer.data(fi_off + 17 * ccomps * dcomps);

            auto g_xxx_xyyzzz = cbuffer.data(fi_off + 18 * ccomps * dcomps);

            auto g_xxx_xyzzzz = cbuffer.data(fi_off + 19 * ccomps * dcomps);

            auto g_xxx_xzzzzz = cbuffer.data(fi_off + 20 * ccomps * dcomps);

            auto g_xxx_yyyyyy = cbuffer.data(fi_off + 21 * ccomps * dcomps);

            auto g_xxx_yyyyyz = cbuffer.data(fi_off + 22 * ccomps * dcomps);

            auto g_xxx_yyyyzz = cbuffer.data(fi_off + 23 * ccomps * dcomps);

            auto g_xxx_yyyzzz = cbuffer.data(fi_off + 24 * ccomps * dcomps);

            auto g_xxx_yyzzzz = cbuffer.data(fi_off + 25 * ccomps * dcomps);

            auto g_xxx_yzzzzz = cbuffer.data(fi_off + 26 * ccomps * dcomps);

            auto g_xxx_zzzzzz = cbuffer.data(fi_off + 27 * ccomps * dcomps);

            auto g_xxy_xxxxxx = cbuffer.data(fi_off + 28 * ccomps * dcomps);

            auto g_xxy_xxxxxy = cbuffer.data(fi_off + 29 * ccomps * dcomps);

            auto g_xxy_xxxxxz = cbuffer.data(fi_off + 30 * ccomps * dcomps);

            auto g_xxy_xxxxyy = cbuffer.data(fi_off + 31 * ccomps * dcomps);

            auto g_xxy_xxxxyz = cbuffer.data(fi_off + 32 * ccomps * dcomps);

            auto g_xxy_xxxxzz = cbuffer.data(fi_off + 33 * ccomps * dcomps);

            auto g_xxy_xxxyyy = cbuffer.data(fi_off + 34 * ccomps * dcomps);

            auto g_xxy_xxxyyz = cbuffer.data(fi_off + 35 * ccomps * dcomps);

            auto g_xxy_xxxyzz = cbuffer.data(fi_off + 36 * ccomps * dcomps);

            auto g_xxy_xxxzzz = cbuffer.data(fi_off + 37 * ccomps * dcomps);

            auto g_xxy_xxyyyy = cbuffer.data(fi_off + 38 * ccomps * dcomps);

            auto g_xxy_xxyyyz = cbuffer.data(fi_off + 39 * ccomps * dcomps);

            auto g_xxy_xxyyzz = cbuffer.data(fi_off + 40 * ccomps * dcomps);

            auto g_xxy_xxyzzz = cbuffer.data(fi_off + 41 * ccomps * dcomps);

            auto g_xxy_xxzzzz = cbuffer.data(fi_off + 42 * ccomps * dcomps);

            auto g_xxy_xyyyyy = cbuffer.data(fi_off + 43 * ccomps * dcomps);

            auto g_xxy_xyyyyz = cbuffer.data(fi_off + 44 * ccomps * dcomps);

            auto g_xxy_xyyyzz = cbuffer.data(fi_off + 45 * ccomps * dcomps);

            auto g_xxy_xyyzzz = cbuffer.data(fi_off + 46 * ccomps * dcomps);

            auto g_xxy_xyzzzz = cbuffer.data(fi_off + 47 * ccomps * dcomps);

            auto g_xxy_xzzzzz = cbuffer.data(fi_off + 48 * ccomps * dcomps);

            auto g_xxy_yyyyyy = cbuffer.data(fi_off + 49 * ccomps * dcomps);

            auto g_xxy_yyyyyz = cbuffer.data(fi_off + 50 * ccomps * dcomps);

            auto g_xxy_yyyyzz = cbuffer.data(fi_off + 51 * ccomps * dcomps);

            auto g_xxy_yyyzzz = cbuffer.data(fi_off + 52 * ccomps * dcomps);

            auto g_xxy_yyzzzz = cbuffer.data(fi_off + 53 * ccomps * dcomps);

            auto g_xxy_yzzzzz = cbuffer.data(fi_off + 54 * ccomps * dcomps);

            auto g_xxy_zzzzzz = cbuffer.data(fi_off + 55 * ccomps * dcomps);

            auto g_xxz_xxxxxx = cbuffer.data(fi_off + 56 * ccomps * dcomps);

            auto g_xxz_xxxxxy = cbuffer.data(fi_off + 57 * ccomps * dcomps);

            auto g_xxz_xxxxxz = cbuffer.data(fi_off + 58 * ccomps * dcomps);

            auto g_xxz_xxxxyy = cbuffer.data(fi_off + 59 * ccomps * dcomps);

            auto g_xxz_xxxxyz = cbuffer.data(fi_off + 60 * ccomps * dcomps);

            auto g_xxz_xxxxzz = cbuffer.data(fi_off + 61 * ccomps * dcomps);

            auto g_xxz_xxxyyy = cbuffer.data(fi_off + 62 * ccomps * dcomps);

            auto g_xxz_xxxyyz = cbuffer.data(fi_off + 63 * ccomps * dcomps);

            auto g_xxz_xxxyzz = cbuffer.data(fi_off + 64 * ccomps * dcomps);

            auto g_xxz_xxxzzz = cbuffer.data(fi_off + 65 * ccomps * dcomps);

            auto g_xxz_xxyyyy = cbuffer.data(fi_off + 66 * ccomps * dcomps);

            auto g_xxz_xxyyyz = cbuffer.data(fi_off + 67 * ccomps * dcomps);

            auto g_xxz_xxyyzz = cbuffer.data(fi_off + 68 * ccomps * dcomps);

            auto g_xxz_xxyzzz = cbuffer.data(fi_off + 69 * ccomps * dcomps);

            auto g_xxz_xxzzzz = cbuffer.data(fi_off + 70 * ccomps * dcomps);

            auto g_xxz_xyyyyy = cbuffer.data(fi_off + 71 * ccomps * dcomps);

            auto g_xxz_xyyyyz = cbuffer.data(fi_off + 72 * ccomps * dcomps);

            auto g_xxz_xyyyzz = cbuffer.data(fi_off + 73 * ccomps * dcomps);

            auto g_xxz_xyyzzz = cbuffer.data(fi_off + 74 * ccomps * dcomps);

            auto g_xxz_xyzzzz = cbuffer.data(fi_off + 75 * ccomps * dcomps);

            auto g_xxz_xzzzzz = cbuffer.data(fi_off + 76 * ccomps * dcomps);

            auto g_xxz_yyyyyy = cbuffer.data(fi_off + 77 * ccomps * dcomps);

            auto g_xxz_yyyyyz = cbuffer.data(fi_off + 78 * ccomps * dcomps);

            auto g_xxz_yyyyzz = cbuffer.data(fi_off + 79 * ccomps * dcomps);

            auto g_xxz_yyyzzz = cbuffer.data(fi_off + 80 * ccomps * dcomps);

            auto g_xxz_yyzzzz = cbuffer.data(fi_off + 81 * ccomps * dcomps);

            auto g_xxz_yzzzzz = cbuffer.data(fi_off + 82 * ccomps * dcomps);

            auto g_xxz_zzzzzz = cbuffer.data(fi_off + 83 * ccomps * dcomps);

            auto g_xyy_xxxxxx = cbuffer.data(fi_off + 84 * ccomps * dcomps);

            auto g_xyy_xxxxxy = cbuffer.data(fi_off + 85 * ccomps * dcomps);

            auto g_xyy_xxxxxz = cbuffer.data(fi_off + 86 * ccomps * dcomps);

            auto g_xyy_xxxxyy = cbuffer.data(fi_off + 87 * ccomps * dcomps);

            auto g_xyy_xxxxyz = cbuffer.data(fi_off + 88 * ccomps * dcomps);

            auto g_xyy_xxxxzz = cbuffer.data(fi_off + 89 * ccomps * dcomps);

            auto g_xyy_xxxyyy = cbuffer.data(fi_off + 90 * ccomps * dcomps);

            auto g_xyy_xxxyyz = cbuffer.data(fi_off + 91 * ccomps * dcomps);

            auto g_xyy_xxxyzz = cbuffer.data(fi_off + 92 * ccomps * dcomps);

            auto g_xyy_xxxzzz = cbuffer.data(fi_off + 93 * ccomps * dcomps);

            auto g_xyy_xxyyyy = cbuffer.data(fi_off + 94 * ccomps * dcomps);

            auto g_xyy_xxyyyz = cbuffer.data(fi_off + 95 * ccomps * dcomps);

            auto g_xyy_xxyyzz = cbuffer.data(fi_off + 96 * ccomps * dcomps);

            auto g_xyy_xxyzzz = cbuffer.data(fi_off + 97 * ccomps * dcomps);

            auto g_xyy_xxzzzz = cbuffer.data(fi_off + 98 * ccomps * dcomps);

            auto g_xyy_xyyyyy = cbuffer.data(fi_off + 99 * ccomps * dcomps);

            auto g_xyy_xyyyyz = cbuffer.data(fi_off + 100 * ccomps * dcomps);

            auto g_xyy_xyyyzz = cbuffer.data(fi_off + 101 * ccomps * dcomps);

            auto g_xyy_xyyzzz = cbuffer.data(fi_off + 102 * ccomps * dcomps);

            auto g_xyy_xyzzzz = cbuffer.data(fi_off + 103 * ccomps * dcomps);

            auto g_xyy_xzzzzz = cbuffer.data(fi_off + 104 * ccomps * dcomps);

            auto g_xyy_yyyyyy = cbuffer.data(fi_off + 105 * ccomps * dcomps);

            auto g_xyy_yyyyyz = cbuffer.data(fi_off + 106 * ccomps * dcomps);

            auto g_xyy_yyyyzz = cbuffer.data(fi_off + 107 * ccomps * dcomps);

            auto g_xyy_yyyzzz = cbuffer.data(fi_off + 108 * ccomps * dcomps);

            auto g_xyy_yyzzzz = cbuffer.data(fi_off + 109 * ccomps * dcomps);

            auto g_xyy_yzzzzz = cbuffer.data(fi_off + 110 * ccomps * dcomps);

            auto g_xyy_zzzzzz = cbuffer.data(fi_off + 111 * ccomps * dcomps);

            auto g_xyz_xxxxxx = cbuffer.data(fi_off + 112 * ccomps * dcomps);

            auto g_xyz_xxxxxy = cbuffer.data(fi_off + 113 * ccomps * dcomps);

            auto g_xyz_xxxxxz = cbuffer.data(fi_off + 114 * ccomps * dcomps);

            auto g_xyz_xxxxyy = cbuffer.data(fi_off + 115 * ccomps * dcomps);

            auto g_xyz_xxxxyz = cbuffer.data(fi_off + 116 * ccomps * dcomps);

            auto g_xyz_xxxxzz = cbuffer.data(fi_off + 117 * ccomps * dcomps);

            auto g_xyz_xxxyyy = cbuffer.data(fi_off + 118 * ccomps * dcomps);

            auto g_xyz_xxxyyz = cbuffer.data(fi_off + 119 * ccomps * dcomps);

            auto g_xyz_xxxyzz = cbuffer.data(fi_off + 120 * ccomps * dcomps);

            auto g_xyz_xxxzzz = cbuffer.data(fi_off + 121 * ccomps * dcomps);

            auto g_xyz_xxyyyy = cbuffer.data(fi_off + 122 * ccomps * dcomps);

            auto g_xyz_xxyyyz = cbuffer.data(fi_off + 123 * ccomps * dcomps);

            auto g_xyz_xxyyzz = cbuffer.data(fi_off + 124 * ccomps * dcomps);

            auto g_xyz_xxyzzz = cbuffer.data(fi_off + 125 * ccomps * dcomps);

            auto g_xyz_xxzzzz = cbuffer.data(fi_off + 126 * ccomps * dcomps);

            auto g_xyz_xyyyyy = cbuffer.data(fi_off + 127 * ccomps * dcomps);

            auto g_xyz_xyyyyz = cbuffer.data(fi_off + 128 * ccomps * dcomps);

            auto g_xyz_xyyyzz = cbuffer.data(fi_off + 129 * ccomps * dcomps);

            auto g_xyz_xyyzzz = cbuffer.data(fi_off + 130 * ccomps * dcomps);

            auto g_xyz_xyzzzz = cbuffer.data(fi_off + 131 * ccomps * dcomps);

            auto g_xyz_xzzzzz = cbuffer.data(fi_off + 132 * ccomps * dcomps);

            auto g_xyz_yyyyyy = cbuffer.data(fi_off + 133 * ccomps * dcomps);

            auto g_xyz_yyyyyz = cbuffer.data(fi_off + 134 * ccomps * dcomps);

            auto g_xyz_yyyyzz = cbuffer.data(fi_off + 135 * ccomps * dcomps);

            auto g_xyz_yyyzzz = cbuffer.data(fi_off + 136 * ccomps * dcomps);

            auto g_xyz_yyzzzz = cbuffer.data(fi_off + 137 * ccomps * dcomps);

            auto g_xyz_yzzzzz = cbuffer.data(fi_off + 138 * ccomps * dcomps);

            auto g_xyz_zzzzzz = cbuffer.data(fi_off + 139 * ccomps * dcomps);

            auto g_xzz_xxxxxx = cbuffer.data(fi_off + 140 * ccomps * dcomps);

            auto g_xzz_xxxxxy = cbuffer.data(fi_off + 141 * ccomps * dcomps);

            auto g_xzz_xxxxxz = cbuffer.data(fi_off + 142 * ccomps * dcomps);

            auto g_xzz_xxxxyy = cbuffer.data(fi_off + 143 * ccomps * dcomps);

            auto g_xzz_xxxxyz = cbuffer.data(fi_off + 144 * ccomps * dcomps);

            auto g_xzz_xxxxzz = cbuffer.data(fi_off + 145 * ccomps * dcomps);

            auto g_xzz_xxxyyy = cbuffer.data(fi_off + 146 * ccomps * dcomps);

            auto g_xzz_xxxyyz = cbuffer.data(fi_off + 147 * ccomps * dcomps);

            auto g_xzz_xxxyzz = cbuffer.data(fi_off + 148 * ccomps * dcomps);

            auto g_xzz_xxxzzz = cbuffer.data(fi_off + 149 * ccomps * dcomps);

            auto g_xzz_xxyyyy = cbuffer.data(fi_off + 150 * ccomps * dcomps);

            auto g_xzz_xxyyyz = cbuffer.data(fi_off + 151 * ccomps * dcomps);

            auto g_xzz_xxyyzz = cbuffer.data(fi_off + 152 * ccomps * dcomps);

            auto g_xzz_xxyzzz = cbuffer.data(fi_off + 153 * ccomps * dcomps);

            auto g_xzz_xxzzzz = cbuffer.data(fi_off + 154 * ccomps * dcomps);

            auto g_xzz_xyyyyy = cbuffer.data(fi_off + 155 * ccomps * dcomps);

            auto g_xzz_xyyyyz = cbuffer.data(fi_off + 156 * ccomps * dcomps);

            auto g_xzz_xyyyzz = cbuffer.data(fi_off + 157 * ccomps * dcomps);

            auto g_xzz_xyyzzz = cbuffer.data(fi_off + 158 * ccomps * dcomps);

            auto g_xzz_xyzzzz = cbuffer.data(fi_off + 159 * ccomps * dcomps);

            auto g_xzz_xzzzzz = cbuffer.data(fi_off + 160 * ccomps * dcomps);

            auto g_xzz_yyyyyy = cbuffer.data(fi_off + 161 * ccomps * dcomps);

            auto g_xzz_yyyyyz = cbuffer.data(fi_off + 162 * ccomps * dcomps);

            auto g_xzz_yyyyzz = cbuffer.data(fi_off + 163 * ccomps * dcomps);

            auto g_xzz_yyyzzz = cbuffer.data(fi_off + 164 * ccomps * dcomps);

            auto g_xzz_yyzzzz = cbuffer.data(fi_off + 165 * ccomps * dcomps);

            auto g_xzz_yzzzzz = cbuffer.data(fi_off + 166 * ccomps * dcomps);

            auto g_xzz_zzzzzz = cbuffer.data(fi_off + 167 * ccomps * dcomps);

            auto g_yyy_xxxxxx = cbuffer.data(fi_off + 168 * ccomps * dcomps);

            auto g_yyy_xxxxxy = cbuffer.data(fi_off + 169 * ccomps * dcomps);

            auto g_yyy_xxxxxz = cbuffer.data(fi_off + 170 * ccomps * dcomps);

            auto g_yyy_xxxxyy = cbuffer.data(fi_off + 171 * ccomps * dcomps);

            auto g_yyy_xxxxyz = cbuffer.data(fi_off + 172 * ccomps * dcomps);

            auto g_yyy_xxxxzz = cbuffer.data(fi_off + 173 * ccomps * dcomps);

            auto g_yyy_xxxyyy = cbuffer.data(fi_off + 174 * ccomps * dcomps);

            auto g_yyy_xxxyyz = cbuffer.data(fi_off + 175 * ccomps * dcomps);

            auto g_yyy_xxxyzz = cbuffer.data(fi_off + 176 * ccomps * dcomps);

            auto g_yyy_xxxzzz = cbuffer.data(fi_off + 177 * ccomps * dcomps);

            auto g_yyy_xxyyyy = cbuffer.data(fi_off + 178 * ccomps * dcomps);

            auto g_yyy_xxyyyz = cbuffer.data(fi_off + 179 * ccomps * dcomps);

            auto g_yyy_xxyyzz = cbuffer.data(fi_off + 180 * ccomps * dcomps);

            auto g_yyy_xxyzzz = cbuffer.data(fi_off + 181 * ccomps * dcomps);

            auto g_yyy_xxzzzz = cbuffer.data(fi_off + 182 * ccomps * dcomps);

            auto g_yyy_xyyyyy = cbuffer.data(fi_off + 183 * ccomps * dcomps);

            auto g_yyy_xyyyyz = cbuffer.data(fi_off + 184 * ccomps * dcomps);

            auto g_yyy_xyyyzz = cbuffer.data(fi_off + 185 * ccomps * dcomps);

            auto g_yyy_xyyzzz = cbuffer.data(fi_off + 186 * ccomps * dcomps);

            auto g_yyy_xyzzzz = cbuffer.data(fi_off + 187 * ccomps * dcomps);

            auto g_yyy_xzzzzz = cbuffer.data(fi_off + 188 * ccomps * dcomps);

            auto g_yyy_yyyyyy = cbuffer.data(fi_off + 189 * ccomps * dcomps);

            auto g_yyy_yyyyyz = cbuffer.data(fi_off + 190 * ccomps * dcomps);

            auto g_yyy_yyyyzz = cbuffer.data(fi_off + 191 * ccomps * dcomps);

            auto g_yyy_yyyzzz = cbuffer.data(fi_off + 192 * ccomps * dcomps);

            auto g_yyy_yyzzzz = cbuffer.data(fi_off + 193 * ccomps * dcomps);

            auto g_yyy_yzzzzz = cbuffer.data(fi_off + 194 * ccomps * dcomps);

            auto g_yyy_zzzzzz = cbuffer.data(fi_off + 195 * ccomps * dcomps);

            auto g_yyz_xxxxxx = cbuffer.data(fi_off + 196 * ccomps * dcomps);

            auto g_yyz_xxxxxy = cbuffer.data(fi_off + 197 * ccomps * dcomps);

            auto g_yyz_xxxxxz = cbuffer.data(fi_off + 198 * ccomps * dcomps);

            auto g_yyz_xxxxyy = cbuffer.data(fi_off + 199 * ccomps * dcomps);

            auto g_yyz_xxxxyz = cbuffer.data(fi_off + 200 * ccomps * dcomps);

            auto g_yyz_xxxxzz = cbuffer.data(fi_off + 201 * ccomps * dcomps);

            auto g_yyz_xxxyyy = cbuffer.data(fi_off + 202 * ccomps * dcomps);

            auto g_yyz_xxxyyz = cbuffer.data(fi_off + 203 * ccomps * dcomps);

            auto g_yyz_xxxyzz = cbuffer.data(fi_off + 204 * ccomps * dcomps);

            auto g_yyz_xxxzzz = cbuffer.data(fi_off + 205 * ccomps * dcomps);

            auto g_yyz_xxyyyy = cbuffer.data(fi_off + 206 * ccomps * dcomps);

            auto g_yyz_xxyyyz = cbuffer.data(fi_off + 207 * ccomps * dcomps);

            auto g_yyz_xxyyzz = cbuffer.data(fi_off + 208 * ccomps * dcomps);

            auto g_yyz_xxyzzz = cbuffer.data(fi_off + 209 * ccomps * dcomps);

            auto g_yyz_xxzzzz = cbuffer.data(fi_off + 210 * ccomps * dcomps);

            auto g_yyz_xyyyyy = cbuffer.data(fi_off + 211 * ccomps * dcomps);

            auto g_yyz_xyyyyz = cbuffer.data(fi_off + 212 * ccomps * dcomps);

            auto g_yyz_xyyyzz = cbuffer.data(fi_off + 213 * ccomps * dcomps);

            auto g_yyz_xyyzzz = cbuffer.data(fi_off + 214 * ccomps * dcomps);

            auto g_yyz_xyzzzz = cbuffer.data(fi_off + 215 * ccomps * dcomps);

            auto g_yyz_xzzzzz = cbuffer.data(fi_off + 216 * ccomps * dcomps);

            auto g_yyz_yyyyyy = cbuffer.data(fi_off + 217 * ccomps * dcomps);

            auto g_yyz_yyyyyz = cbuffer.data(fi_off + 218 * ccomps * dcomps);

            auto g_yyz_yyyyzz = cbuffer.data(fi_off + 219 * ccomps * dcomps);

            auto g_yyz_yyyzzz = cbuffer.data(fi_off + 220 * ccomps * dcomps);

            auto g_yyz_yyzzzz = cbuffer.data(fi_off + 221 * ccomps * dcomps);

            auto g_yyz_yzzzzz = cbuffer.data(fi_off + 222 * ccomps * dcomps);

            auto g_yyz_zzzzzz = cbuffer.data(fi_off + 223 * ccomps * dcomps);

            auto g_yzz_xxxxxx = cbuffer.data(fi_off + 224 * ccomps * dcomps);

            auto g_yzz_xxxxxy = cbuffer.data(fi_off + 225 * ccomps * dcomps);

            auto g_yzz_xxxxxz = cbuffer.data(fi_off + 226 * ccomps * dcomps);

            auto g_yzz_xxxxyy = cbuffer.data(fi_off + 227 * ccomps * dcomps);

            auto g_yzz_xxxxyz = cbuffer.data(fi_off + 228 * ccomps * dcomps);

            auto g_yzz_xxxxzz = cbuffer.data(fi_off + 229 * ccomps * dcomps);

            auto g_yzz_xxxyyy = cbuffer.data(fi_off + 230 * ccomps * dcomps);

            auto g_yzz_xxxyyz = cbuffer.data(fi_off + 231 * ccomps * dcomps);

            auto g_yzz_xxxyzz = cbuffer.data(fi_off + 232 * ccomps * dcomps);

            auto g_yzz_xxxzzz = cbuffer.data(fi_off + 233 * ccomps * dcomps);

            auto g_yzz_xxyyyy = cbuffer.data(fi_off + 234 * ccomps * dcomps);

            auto g_yzz_xxyyyz = cbuffer.data(fi_off + 235 * ccomps * dcomps);

            auto g_yzz_xxyyzz = cbuffer.data(fi_off + 236 * ccomps * dcomps);

            auto g_yzz_xxyzzz = cbuffer.data(fi_off + 237 * ccomps * dcomps);

            auto g_yzz_xxzzzz = cbuffer.data(fi_off + 238 * ccomps * dcomps);

            auto g_yzz_xyyyyy = cbuffer.data(fi_off + 239 * ccomps * dcomps);

            auto g_yzz_xyyyyz = cbuffer.data(fi_off + 240 * ccomps * dcomps);

            auto g_yzz_xyyyzz = cbuffer.data(fi_off + 241 * ccomps * dcomps);

            auto g_yzz_xyyzzz = cbuffer.data(fi_off + 242 * ccomps * dcomps);

            auto g_yzz_xyzzzz = cbuffer.data(fi_off + 243 * ccomps * dcomps);

            auto g_yzz_xzzzzz = cbuffer.data(fi_off + 244 * ccomps * dcomps);

            auto g_yzz_yyyyyy = cbuffer.data(fi_off + 245 * ccomps * dcomps);

            auto g_yzz_yyyyyz = cbuffer.data(fi_off + 246 * ccomps * dcomps);

            auto g_yzz_yyyyzz = cbuffer.data(fi_off + 247 * ccomps * dcomps);

            auto g_yzz_yyyzzz = cbuffer.data(fi_off + 248 * ccomps * dcomps);

            auto g_yzz_yyzzzz = cbuffer.data(fi_off + 249 * ccomps * dcomps);

            auto g_yzz_yzzzzz = cbuffer.data(fi_off + 250 * ccomps * dcomps);

            auto g_yzz_zzzzzz = cbuffer.data(fi_off + 251 * ccomps * dcomps);

            auto g_zzz_xxxxxx = cbuffer.data(fi_off + 252 * ccomps * dcomps);

            auto g_zzz_xxxxxy = cbuffer.data(fi_off + 253 * ccomps * dcomps);

            auto g_zzz_xxxxxz = cbuffer.data(fi_off + 254 * ccomps * dcomps);

            auto g_zzz_xxxxyy = cbuffer.data(fi_off + 255 * ccomps * dcomps);

            auto g_zzz_xxxxyz = cbuffer.data(fi_off + 256 * ccomps * dcomps);

            auto g_zzz_xxxxzz = cbuffer.data(fi_off + 257 * ccomps * dcomps);

            auto g_zzz_xxxyyy = cbuffer.data(fi_off + 258 * ccomps * dcomps);

            auto g_zzz_xxxyyz = cbuffer.data(fi_off + 259 * ccomps * dcomps);

            auto g_zzz_xxxyzz = cbuffer.data(fi_off + 260 * ccomps * dcomps);

            auto g_zzz_xxxzzz = cbuffer.data(fi_off + 261 * ccomps * dcomps);

            auto g_zzz_xxyyyy = cbuffer.data(fi_off + 262 * ccomps * dcomps);

            auto g_zzz_xxyyyz = cbuffer.data(fi_off + 263 * ccomps * dcomps);

            auto g_zzz_xxyyzz = cbuffer.data(fi_off + 264 * ccomps * dcomps);

            auto g_zzz_xxyzzz = cbuffer.data(fi_off + 265 * ccomps * dcomps);

            auto g_zzz_xxzzzz = cbuffer.data(fi_off + 266 * ccomps * dcomps);

            auto g_zzz_xyyyyy = cbuffer.data(fi_off + 267 * ccomps * dcomps);

            auto g_zzz_xyyyyz = cbuffer.data(fi_off + 268 * ccomps * dcomps);

            auto g_zzz_xyyyzz = cbuffer.data(fi_off + 269 * ccomps * dcomps);

            auto g_zzz_xyyzzz = cbuffer.data(fi_off + 270 * ccomps * dcomps);

            auto g_zzz_xyzzzz = cbuffer.data(fi_off + 271 * ccomps * dcomps);

            auto g_zzz_xzzzzz = cbuffer.data(fi_off + 272 * ccomps * dcomps);

            auto g_zzz_yyyyyy = cbuffer.data(fi_off + 273 * ccomps * dcomps);

            auto g_zzz_yyyyyz = cbuffer.data(fi_off + 274 * ccomps * dcomps);

            auto g_zzz_yyyyzz = cbuffer.data(fi_off + 275 * ccomps * dcomps);

            auto g_zzz_yyyzzz = cbuffer.data(fi_off + 276 * ccomps * dcomps);

            auto g_zzz_yyzzzz = cbuffer.data(fi_off + 277 * ccomps * dcomps);

            auto g_zzz_yzzzzz = cbuffer.data(fi_off + 278 * ccomps * dcomps);

            auto g_zzz_zzzzzz = cbuffer.data(fi_off + 279 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FISS

            const auto fi_geom_10_off = idx_geom_10_fixx + i * dcomps + j;

            auto g_x_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 659 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 674 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 675 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 676 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 677 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 678 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 679 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 680 * ccomps * dcomps);

            auto g_z_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 681 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 682 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 683 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 684 * ccomps * dcomps);

            auto g_z_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 685 * ccomps * dcomps);

            auto g_z_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 686 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 687 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 688 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 689 * ccomps * dcomps);

            auto g_z_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 690 * ccomps * dcomps);

            auto g_z_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 691 * ccomps * dcomps);

            auto g_z_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 692 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 693 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 694 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 695 * ccomps * dcomps);

            auto g_z_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 696 * ccomps * dcomps);

            auto g_z_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 697 * ccomps * dcomps);

            auto g_z_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 698 * ccomps * dcomps);

            auto g_z_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 699 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 700 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 701 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 702 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 703 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 704 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 705 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 706 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 707 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 708 * ccomps * dcomps);

            auto g_z_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 709 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 710 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 711 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 712 * ccomps * dcomps);

            auto g_z_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 713 * ccomps * dcomps);

            auto g_z_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 714 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 715 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 716 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 717 * ccomps * dcomps);

            auto g_z_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 718 * ccomps * dcomps);

            auto g_z_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 719 * ccomps * dcomps);

            auto g_z_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 727 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 729 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 739 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 749 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 755 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 759 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 769 * ccomps * dcomps);

            auto g_z_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 779 * ccomps * dcomps);

            auto g_z_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 783 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 789 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 799 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 809 * ccomps * dcomps);

            auto g_z_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 811 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 819 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 829 * ccomps * dcomps);

            auto g_z_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 839 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FKSS

            const auto fk_geom_10_off = idx_geom_10_fkxx + i * dcomps + j;

            auto g_x_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxxx = cbuffer.data(fk_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxxz = cbuffer.data(fk_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxzz = cbuffer.data(fk_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxzzz = cbuffer.data(fk_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxxx = cbuffer.data(fk_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxxz = cbuffer.data(fk_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxzz = cbuffer.data(fk_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxzzz = cbuffer.data(fk_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 719 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 727 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 729 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 739 * ccomps * dcomps);

            auto g_z_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 749 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 755 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxxx = cbuffer.data(fk_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxxz = cbuffer.data(fk_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 759 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxzz = cbuffer.data(fk_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxzzz = cbuffer.data(fk_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 769 * ccomps * dcomps);

            auto g_z_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 779 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 783 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 789 * ccomps * dcomps);

            auto g_z_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 799 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 809 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 811 * ccomps * dcomps);

            auto g_z_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 819 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 829 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 839 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 854 * ccomps * dcomps);

            auto g_z_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 867 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 869 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 944 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 951 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 959 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 974 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 979 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 989 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 1004 * ccomps * dcomps);

            auto g_z_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 1007 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 1019 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 1034 * ccomps * dcomps);

            auto g_z_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 1035 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 1049 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 1063 * ccomps * dcomps);

            auto g_z_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 1064 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 1079 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gixx

            const auto gi_geom_10_off = idx_geom_10_gixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxx, g_x_0_xxx_xxxxxxy, g_x_0_xxx_xxxxxxz, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxxzz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxxzzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxxzzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xxzzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_xzzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzzz, g_xxx_xxxxxx, g_xxx_xxxxxy, g_xxx_xxxxxz, g_xxx_xxxxyy, g_xxx_xxxxyz, g_xxx_xxxxzz, g_xxx_xxxyyy, g_xxx_xxxyyz, g_xxx_xxxyzz, g_xxx_xxxzzz, g_xxx_xxyyyy, g_xxx_xxyyyz, g_xxx_xxyyzz, g_xxx_xxyzzz, g_xxx_xxzzzz, g_xxx_xyyyyy, g_xxx_xyyyyz, g_xxx_xyyyzz, g_xxx_xyyzzz, g_xxx_xyzzzz, g_xxx_xzzzzz, g_xxx_yyyyyy, g_xxx_yyyyyz, g_xxx_yyyyzz, g_xxx_yyyzzz, g_xxx_yyzzzz, g_xxx_yzzzzz, g_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_xxxxxx[k] = -g_xxx_xxxxxx[k] - g_x_0_xxx_xxxxxx[k] * ab_x + g_x_0_xxx_xxxxxxx[k];

                g_x_0_xxxx_xxxxxy[k] = -g_xxx_xxxxxy[k] - g_x_0_xxx_xxxxxy[k] * ab_x + g_x_0_xxx_xxxxxxy[k];

                g_x_0_xxxx_xxxxxz[k] = -g_xxx_xxxxxz[k] - g_x_0_xxx_xxxxxz[k] * ab_x + g_x_0_xxx_xxxxxxz[k];

                g_x_0_xxxx_xxxxyy[k] = -g_xxx_xxxxyy[k] - g_x_0_xxx_xxxxyy[k] * ab_x + g_x_0_xxx_xxxxxyy[k];

                g_x_0_xxxx_xxxxyz[k] = -g_xxx_xxxxyz[k] - g_x_0_xxx_xxxxyz[k] * ab_x + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxx_xxxxzz[k] = -g_xxx_xxxxzz[k] - g_x_0_xxx_xxxxzz[k] * ab_x + g_x_0_xxx_xxxxxzz[k];

                g_x_0_xxxx_xxxyyy[k] = -g_xxx_xxxyyy[k] - g_x_0_xxx_xxxyyy[k] * ab_x + g_x_0_xxx_xxxxyyy[k];

                g_x_0_xxxx_xxxyyz[k] = -g_xxx_xxxyyz[k] - g_x_0_xxx_xxxyyz[k] * ab_x + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxx_xxxyzz[k] = -g_xxx_xxxyzz[k] - g_x_0_xxx_xxxyzz[k] * ab_x + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxx_xxxzzz[k] = -g_xxx_xxxzzz[k] - g_x_0_xxx_xxxzzz[k] * ab_x + g_x_0_xxx_xxxxzzz[k];

                g_x_0_xxxx_xxyyyy[k] = -g_xxx_xxyyyy[k] - g_x_0_xxx_xxyyyy[k] * ab_x + g_x_0_xxx_xxxyyyy[k];

                g_x_0_xxxx_xxyyyz[k] = -g_xxx_xxyyyz[k] - g_x_0_xxx_xxyyyz[k] * ab_x + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxx_xxyyzz[k] = -g_xxx_xxyyzz[k] - g_x_0_xxx_xxyyzz[k] * ab_x + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxx_xxyzzz[k] = -g_xxx_xxyzzz[k] - g_x_0_xxx_xxyzzz[k] * ab_x + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxx_xxzzzz[k] = -g_xxx_xxzzzz[k] - g_x_0_xxx_xxzzzz[k] * ab_x + g_x_0_xxx_xxxzzzz[k];

                g_x_0_xxxx_xyyyyy[k] = -g_xxx_xyyyyy[k] - g_x_0_xxx_xyyyyy[k] * ab_x + g_x_0_xxx_xxyyyyy[k];

                g_x_0_xxxx_xyyyyz[k] = -g_xxx_xyyyyz[k] - g_x_0_xxx_xyyyyz[k] * ab_x + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxx_xyyyzz[k] = -g_xxx_xyyyzz[k] - g_x_0_xxx_xyyyzz[k] * ab_x + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxx_xyyzzz[k] = -g_xxx_xyyzzz[k] - g_x_0_xxx_xyyzzz[k] * ab_x + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxx_xyzzzz[k] = -g_xxx_xyzzzz[k] - g_x_0_xxx_xyzzzz[k] * ab_x + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxx_xzzzzz[k] = -g_xxx_xzzzzz[k] - g_x_0_xxx_xzzzzz[k] * ab_x + g_x_0_xxx_xxzzzzz[k];

                g_x_0_xxxx_yyyyyy[k] = -g_xxx_yyyyyy[k] - g_x_0_xxx_yyyyyy[k] * ab_x + g_x_0_xxx_xyyyyyy[k];

                g_x_0_xxxx_yyyyyz[k] = -g_xxx_yyyyyz[k] - g_x_0_xxx_yyyyyz[k] * ab_x + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxx_yyyyzz[k] = -g_xxx_yyyyzz[k] - g_x_0_xxx_yyyyzz[k] * ab_x + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxx_yyyzzz[k] = -g_xxx_yyyzzz[k] - g_x_0_xxx_yyyzzz[k] * ab_x + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxx_yyzzzz[k] = -g_xxx_yyzzzz[k] - g_x_0_xxx_yyzzzz[k] * ab_x + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxx_yzzzzz[k] = -g_xxx_yzzzzz[k] - g_x_0_xxx_yzzzzz[k] * ab_x + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxx_zzzzzz[k] = -g_xxx_zzzzzz[k] - g_x_0_xxx_zzzzzz[k] * ab_x + g_x_0_xxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxy, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyyy, g_x_0_xxx_yyyyyyz, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyyzz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyyzzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyyzzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yyzzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_yzzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxxy_xxxxxx, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxxz, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxzz, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxzzz, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxzzzz, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xzzzzz, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_xxxxxx[k] = -g_x_0_xxx_xxxxxx[k] * ab_y + g_x_0_xxx_xxxxxxy[k];

                g_x_0_xxxy_xxxxxy[k] = -g_x_0_xxx_xxxxxy[k] * ab_y + g_x_0_xxx_xxxxxyy[k];

                g_x_0_xxxy_xxxxxz[k] = -g_x_0_xxx_xxxxxz[k] * ab_y + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxy_xxxxyy[k] = -g_x_0_xxx_xxxxyy[k] * ab_y + g_x_0_xxx_xxxxyyy[k];

                g_x_0_xxxy_xxxxyz[k] = -g_x_0_xxx_xxxxyz[k] * ab_y + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxy_xxxxzz[k] = -g_x_0_xxx_xxxxzz[k] * ab_y + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxy_xxxyyy[k] = -g_x_0_xxx_xxxyyy[k] * ab_y + g_x_0_xxx_xxxyyyy[k];

                g_x_0_xxxy_xxxyyz[k] = -g_x_0_xxx_xxxyyz[k] * ab_y + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxy_xxxyzz[k] = -g_x_0_xxx_xxxyzz[k] * ab_y + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxy_xxxzzz[k] = -g_x_0_xxx_xxxzzz[k] * ab_y + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxy_xxyyyy[k] = -g_x_0_xxx_xxyyyy[k] * ab_y + g_x_0_xxx_xxyyyyy[k];

                g_x_0_xxxy_xxyyyz[k] = -g_x_0_xxx_xxyyyz[k] * ab_y + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxy_xxyyzz[k] = -g_x_0_xxx_xxyyzz[k] * ab_y + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxy_xxyzzz[k] = -g_x_0_xxx_xxyzzz[k] * ab_y + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxy_xxzzzz[k] = -g_x_0_xxx_xxzzzz[k] * ab_y + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxy_xyyyyy[k] = -g_x_0_xxx_xyyyyy[k] * ab_y + g_x_0_xxx_xyyyyyy[k];

                g_x_0_xxxy_xyyyyz[k] = -g_x_0_xxx_xyyyyz[k] * ab_y + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxy_xyyyzz[k] = -g_x_0_xxx_xyyyzz[k] * ab_y + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxy_xyyzzz[k] = -g_x_0_xxx_xyyzzz[k] * ab_y + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxy_xyzzzz[k] = -g_x_0_xxx_xyzzzz[k] * ab_y + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxy_xzzzzz[k] = -g_x_0_xxx_xzzzzz[k] * ab_y + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxy_yyyyyy[k] = -g_x_0_xxx_yyyyyy[k] * ab_y + g_x_0_xxx_yyyyyyy[k];

                g_x_0_xxxy_yyyyyz[k] = -g_x_0_xxx_yyyyyz[k] * ab_y + g_x_0_xxx_yyyyyyz[k];

                g_x_0_xxxy_yyyyzz[k] = -g_x_0_xxx_yyyyzz[k] * ab_y + g_x_0_xxx_yyyyyzz[k];

                g_x_0_xxxy_yyyzzz[k] = -g_x_0_xxx_yyyzzz[k] * ab_y + g_x_0_xxx_yyyyzzz[k];

                g_x_0_xxxy_yyzzzz[k] = -g_x_0_xxx_yyzzzz[k] * ab_y + g_x_0_xxx_yyyzzzz[k];

                g_x_0_xxxy_yzzzzz[k] = -g_x_0_xxx_yzzzzz[k] * ab_y + g_x_0_xxx_yyzzzzz[k];

                g_x_0_xxxy_zzzzzz[k] = -g_x_0_xxx_zzzzzz[k] * ab_y + g_x_0_xxx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxz, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxxzz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxxzzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxxzzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xxzzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_xzzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyyz, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyyzz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyyzzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyyzzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yyzzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_yzzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxx_zzzzzzz, g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_xxxxxx[k] = -g_x_0_xxx_xxxxxx[k] * ab_z + g_x_0_xxx_xxxxxxz[k];

                g_x_0_xxxz_xxxxxy[k] = -g_x_0_xxx_xxxxxy[k] * ab_z + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxz_xxxxxz[k] = -g_x_0_xxx_xxxxxz[k] * ab_z + g_x_0_xxx_xxxxxzz[k];

                g_x_0_xxxz_xxxxyy[k] = -g_x_0_xxx_xxxxyy[k] * ab_z + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxz_xxxxyz[k] = -g_x_0_xxx_xxxxyz[k] * ab_z + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxz_xxxxzz[k] = -g_x_0_xxx_xxxxzz[k] * ab_z + g_x_0_xxx_xxxxzzz[k];

                g_x_0_xxxz_xxxyyy[k] = -g_x_0_xxx_xxxyyy[k] * ab_z + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxz_xxxyyz[k] = -g_x_0_xxx_xxxyyz[k] * ab_z + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxz_xxxyzz[k] = -g_x_0_xxx_xxxyzz[k] * ab_z + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxz_xxxzzz[k] = -g_x_0_xxx_xxxzzz[k] * ab_z + g_x_0_xxx_xxxzzzz[k];

                g_x_0_xxxz_xxyyyy[k] = -g_x_0_xxx_xxyyyy[k] * ab_z + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxz_xxyyyz[k] = -g_x_0_xxx_xxyyyz[k] * ab_z + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxz_xxyyzz[k] = -g_x_0_xxx_xxyyzz[k] * ab_z + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxz_xxyzzz[k] = -g_x_0_xxx_xxyzzz[k] * ab_z + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxz_xxzzzz[k] = -g_x_0_xxx_xxzzzz[k] * ab_z + g_x_0_xxx_xxzzzzz[k];

                g_x_0_xxxz_xyyyyy[k] = -g_x_0_xxx_xyyyyy[k] * ab_z + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxz_xyyyyz[k] = -g_x_0_xxx_xyyyyz[k] * ab_z + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxz_xyyyzz[k] = -g_x_0_xxx_xyyyzz[k] * ab_z + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxz_xyyzzz[k] = -g_x_0_xxx_xyyzzz[k] * ab_z + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxz_xyzzzz[k] = -g_x_0_xxx_xyzzzz[k] * ab_z + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxz_xzzzzz[k] = -g_x_0_xxx_xzzzzz[k] * ab_z + g_x_0_xxx_xzzzzzz[k];

                g_x_0_xxxz_yyyyyy[k] = -g_x_0_xxx_yyyyyy[k] * ab_z + g_x_0_xxx_yyyyyyz[k];

                g_x_0_xxxz_yyyyyz[k] = -g_x_0_xxx_yyyyyz[k] * ab_z + g_x_0_xxx_yyyyyzz[k];

                g_x_0_xxxz_yyyyzz[k] = -g_x_0_xxx_yyyyzz[k] * ab_z + g_x_0_xxx_yyyyzzz[k];

                g_x_0_xxxz_yyyzzz[k] = -g_x_0_xxx_yyyzzz[k] * ab_z + g_x_0_xxx_yyyzzzz[k];

                g_x_0_xxxz_yyzzzz[k] = -g_x_0_xxx_yyzzzz[k] * ab_z + g_x_0_xxx_yyzzzzz[k];

                g_x_0_xxxz_yzzzzz[k] = -g_x_0_xxx_yzzzzz[k] * ab_z + g_x_0_xxx_yzzzzzz[k];

                g_x_0_xxxz_zzzzzz[k] = -g_x_0_xxx_zzzzzz[k] * ab_z + g_x_0_xxx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxy_xxxxxx, g_x_0_xxy_xxxxxxy, g_x_0_xxy_xxxxxy, g_x_0_xxy_xxxxxyy, g_x_0_xxy_xxxxxyz, g_x_0_xxy_xxxxxz, g_x_0_xxy_xxxxyy, g_x_0_xxy_xxxxyyy, g_x_0_xxy_xxxxyyz, g_x_0_xxy_xxxxyz, g_x_0_xxy_xxxxyzz, g_x_0_xxy_xxxxzz, g_x_0_xxy_xxxyyy, g_x_0_xxy_xxxyyyy, g_x_0_xxy_xxxyyyz, g_x_0_xxy_xxxyyz, g_x_0_xxy_xxxyyzz, g_x_0_xxy_xxxyzz, g_x_0_xxy_xxxyzzz, g_x_0_xxy_xxxzzz, g_x_0_xxy_xxyyyy, g_x_0_xxy_xxyyyyy, g_x_0_xxy_xxyyyyz, g_x_0_xxy_xxyyyz, g_x_0_xxy_xxyyyzz, g_x_0_xxy_xxyyzz, g_x_0_xxy_xxyyzzz, g_x_0_xxy_xxyzzz, g_x_0_xxy_xxyzzzz, g_x_0_xxy_xxzzzz, g_x_0_xxy_xyyyyy, g_x_0_xxy_xyyyyyy, g_x_0_xxy_xyyyyyz, g_x_0_xxy_xyyyyz, g_x_0_xxy_xyyyyzz, g_x_0_xxy_xyyyzz, g_x_0_xxy_xyyyzzz, g_x_0_xxy_xyyzzz, g_x_0_xxy_xyyzzzz, g_x_0_xxy_xyzzzz, g_x_0_xxy_xyzzzzz, g_x_0_xxy_xzzzzz, g_x_0_xxy_yyyyyy, g_x_0_xxy_yyyyyyy, g_x_0_xxy_yyyyyyz, g_x_0_xxy_yyyyyz, g_x_0_xxy_yyyyyzz, g_x_0_xxy_yyyyzz, g_x_0_xxy_yyyyzzz, g_x_0_xxy_yyyzzz, g_x_0_xxy_yyyzzzz, g_x_0_xxy_yyzzzz, g_x_0_xxy_yyzzzzz, g_x_0_xxy_yzzzzz, g_x_0_xxy_yzzzzzz, g_x_0_xxy_zzzzzz, g_x_0_xxyy_xxxxxx, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxxz, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxzz, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxzzz, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxzzzz, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xzzzzz, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_xxxxxx[k] = -g_x_0_xxy_xxxxxx[k] * ab_y + g_x_0_xxy_xxxxxxy[k];

                g_x_0_xxyy_xxxxxy[k] = -g_x_0_xxy_xxxxxy[k] * ab_y + g_x_0_xxy_xxxxxyy[k];

                g_x_0_xxyy_xxxxxz[k] = -g_x_0_xxy_xxxxxz[k] * ab_y + g_x_0_xxy_xxxxxyz[k];

                g_x_0_xxyy_xxxxyy[k] = -g_x_0_xxy_xxxxyy[k] * ab_y + g_x_0_xxy_xxxxyyy[k];

                g_x_0_xxyy_xxxxyz[k] = -g_x_0_xxy_xxxxyz[k] * ab_y + g_x_0_xxy_xxxxyyz[k];

                g_x_0_xxyy_xxxxzz[k] = -g_x_0_xxy_xxxxzz[k] * ab_y + g_x_0_xxy_xxxxyzz[k];

                g_x_0_xxyy_xxxyyy[k] = -g_x_0_xxy_xxxyyy[k] * ab_y + g_x_0_xxy_xxxyyyy[k];

                g_x_0_xxyy_xxxyyz[k] = -g_x_0_xxy_xxxyyz[k] * ab_y + g_x_0_xxy_xxxyyyz[k];

                g_x_0_xxyy_xxxyzz[k] = -g_x_0_xxy_xxxyzz[k] * ab_y + g_x_0_xxy_xxxyyzz[k];

                g_x_0_xxyy_xxxzzz[k] = -g_x_0_xxy_xxxzzz[k] * ab_y + g_x_0_xxy_xxxyzzz[k];

                g_x_0_xxyy_xxyyyy[k] = -g_x_0_xxy_xxyyyy[k] * ab_y + g_x_0_xxy_xxyyyyy[k];

                g_x_0_xxyy_xxyyyz[k] = -g_x_0_xxy_xxyyyz[k] * ab_y + g_x_0_xxy_xxyyyyz[k];

                g_x_0_xxyy_xxyyzz[k] = -g_x_0_xxy_xxyyzz[k] * ab_y + g_x_0_xxy_xxyyyzz[k];

                g_x_0_xxyy_xxyzzz[k] = -g_x_0_xxy_xxyzzz[k] * ab_y + g_x_0_xxy_xxyyzzz[k];

                g_x_0_xxyy_xxzzzz[k] = -g_x_0_xxy_xxzzzz[k] * ab_y + g_x_0_xxy_xxyzzzz[k];

                g_x_0_xxyy_xyyyyy[k] = -g_x_0_xxy_xyyyyy[k] * ab_y + g_x_0_xxy_xyyyyyy[k];

                g_x_0_xxyy_xyyyyz[k] = -g_x_0_xxy_xyyyyz[k] * ab_y + g_x_0_xxy_xyyyyyz[k];

                g_x_0_xxyy_xyyyzz[k] = -g_x_0_xxy_xyyyzz[k] * ab_y + g_x_0_xxy_xyyyyzz[k];

                g_x_0_xxyy_xyyzzz[k] = -g_x_0_xxy_xyyzzz[k] * ab_y + g_x_0_xxy_xyyyzzz[k];

                g_x_0_xxyy_xyzzzz[k] = -g_x_0_xxy_xyzzzz[k] * ab_y + g_x_0_xxy_xyyzzzz[k];

                g_x_0_xxyy_xzzzzz[k] = -g_x_0_xxy_xzzzzz[k] * ab_y + g_x_0_xxy_xyzzzzz[k];

                g_x_0_xxyy_yyyyyy[k] = -g_x_0_xxy_yyyyyy[k] * ab_y + g_x_0_xxy_yyyyyyy[k];

                g_x_0_xxyy_yyyyyz[k] = -g_x_0_xxy_yyyyyz[k] * ab_y + g_x_0_xxy_yyyyyyz[k];

                g_x_0_xxyy_yyyyzz[k] = -g_x_0_xxy_yyyyzz[k] * ab_y + g_x_0_xxy_yyyyyzz[k];

                g_x_0_xxyy_yyyzzz[k] = -g_x_0_xxy_yyyzzz[k] * ab_y + g_x_0_xxy_yyyyzzz[k];

                g_x_0_xxyy_yyzzzz[k] = -g_x_0_xxy_yyzzzz[k] * ab_y + g_x_0_xxy_yyyzzzz[k];

                g_x_0_xxyy_yzzzzz[k] = -g_x_0_xxy_yzzzzz[k] * ab_y + g_x_0_xxy_yyzzzzz[k];

                g_x_0_xxyy_zzzzzz[k] = -g_x_0_xxy_zzzzzz[k] * ab_y + g_x_0_xxy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyz_xxxxxx, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxxz, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxzz, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxzzz, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxzzzz, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xzzzzz, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_zzzzzz, g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxxy, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxyy, g_x_0_xxz_xxxxxyz, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyyy, g_x_0_xxz_xxxxyyz, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxyzz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyyy, g_x_0_xxz_xxxyyyz, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyyzz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxyzzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyyy, g_x_0_xxz_xxyyyyz, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyyzz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyyzzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxyzzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyyy, g_x_0_xxz_xyyyyyz, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyyzz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyyzzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyyzzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xyzzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyyy, g_x_0_xxz_yyyyyyz, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyyzz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyyzzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyyzzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yyzzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_yzzzzzz, g_x_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_xxxxxx[k] = -g_x_0_xxz_xxxxxx[k] * ab_y + g_x_0_xxz_xxxxxxy[k];

                g_x_0_xxyz_xxxxxy[k] = -g_x_0_xxz_xxxxxy[k] * ab_y + g_x_0_xxz_xxxxxyy[k];

                g_x_0_xxyz_xxxxxz[k] = -g_x_0_xxz_xxxxxz[k] * ab_y + g_x_0_xxz_xxxxxyz[k];

                g_x_0_xxyz_xxxxyy[k] = -g_x_0_xxz_xxxxyy[k] * ab_y + g_x_0_xxz_xxxxyyy[k];

                g_x_0_xxyz_xxxxyz[k] = -g_x_0_xxz_xxxxyz[k] * ab_y + g_x_0_xxz_xxxxyyz[k];

                g_x_0_xxyz_xxxxzz[k] = -g_x_0_xxz_xxxxzz[k] * ab_y + g_x_0_xxz_xxxxyzz[k];

                g_x_0_xxyz_xxxyyy[k] = -g_x_0_xxz_xxxyyy[k] * ab_y + g_x_0_xxz_xxxyyyy[k];

                g_x_0_xxyz_xxxyyz[k] = -g_x_0_xxz_xxxyyz[k] * ab_y + g_x_0_xxz_xxxyyyz[k];

                g_x_0_xxyz_xxxyzz[k] = -g_x_0_xxz_xxxyzz[k] * ab_y + g_x_0_xxz_xxxyyzz[k];

                g_x_0_xxyz_xxxzzz[k] = -g_x_0_xxz_xxxzzz[k] * ab_y + g_x_0_xxz_xxxyzzz[k];

                g_x_0_xxyz_xxyyyy[k] = -g_x_0_xxz_xxyyyy[k] * ab_y + g_x_0_xxz_xxyyyyy[k];

                g_x_0_xxyz_xxyyyz[k] = -g_x_0_xxz_xxyyyz[k] * ab_y + g_x_0_xxz_xxyyyyz[k];

                g_x_0_xxyz_xxyyzz[k] = -g_x_0_xxz_xxyyzz[k] * ab_y + g_x_0_xxz_xxyyyzz[k];

                g_x_0_xxyz_xxyzzz[k] = -g_x_0_xxz_xxyzzz[k] * ab_y + g_x_0_xxz_xxyyzzz[k];

                g_x_0_xxyz_xxzzzz[k] = -g_x_0_xxz_xxzzzz[k] * ab_y + g_x_0_xxz_xxyzzzz[k];

                g_x_0_xxyz_xyyyyy[k] = -g_x_0_xxz_xyyyyy[k] * ab_y + g_x_0_xxz_xyyyyyy[k];

                g_x_0_xxyz_xyyyyz[k] = -g_x_0_xxz_xyyyyz[k] * ab_y + g_x_0_xxz_xyyyyyz[k];

                g_x_0_xxyz_xyyyzz[k] = -g_x_0_xxz_xyyyzz[k] * ab_y + g_x_0_xxz_xyyyyzz[k];

                g_x_0_xxyz_xyyzzz[k] = -g_x_0_xxz_xyyzzz[k] * ab_y + g_x_0_xxz_xyyyzzz[k];

                g_x_0_xxyz_xyzzzz[k] = -g_x_0_xxz_xyzzzz[k] * ab_y + g_x_0_xxz_xyyzzzz[k];

                g_x_0_xxyz_xzzzzz[k] = -g_x_0_xxz_xzzzzz[k] * ab_y + g_x_0_xxz_xyzzzzz[k];

                g_x_0_xxyz_yyyyyy[k] = -g_x_0_xxz_yyyyyy[k] * ab_y + g_x_0_xxz_yyyyyyy[k];

                g_x_0_xxyz_yyyyyz[k] = -g_x_0_xxz_yyyyyz[k] * ab_y + g_x_0_xxz_yyyyyyz[k];

                g_x_0_xxyz_yyyyzz[k] = -g_x_0_xxz_yyyyzz[k] * ab_y + g_x_0_xxz_yyyyyzz[k];

                g_x_0_xxyz_yyyzzz[k] = -g_x_0_xxz_yyyzzz[k] * ab_y + g_x_0_xxz_yyyyzzz[k];

                g_x_0_xxyz_yyzzzz[k] = -g_x_0_xxz_yyzzzz[k] * ab_y + g_x_0_xxz_yyyzzzz[k];

                g_x_0_xxyz_yzzzzz[k] = -g_x_0_xxz_yzzzzz[k] * ab_y + g_x_0_xxz_yyzzzzz[k];

                g_x_0_xxyz_zzzzzz[k] = -g_x_0_xxz_zzzzzz[k] * ab_y + g_x_0_xxz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxxz, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxyz, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxxzz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyyz, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxyzz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxxzzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyyz, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyyzz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxyzzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxxzzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyyz, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyyzz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyyzzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxyzzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xxzzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyyz, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyyzz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyyzzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyyzzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xyzzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_xzzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyyz, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyyzz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyyzzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyyzzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yyzzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_yzzzzzz, g_x_0_xxz_zzzzzz, g_x_0_xxz_zzzzzzz, g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_xxxxxx[k] = -g_x_0_xxz_xxxxxx[k] * ab_z + g_x_0_xxz_xxxxxxz[k];

                g_x_0_xxzz_xxxxxy[k] = -g_x_0_xxz_xxxxxy[k] * ab_z + g_x_0_xxz_xxxxxyz[k];

                g_x_0_xxzz_xxxxxz[k] = -g_x_0_xxz_xxxxxz[k] * ab_z + g_x_0_xxz_xxxxxzz[k];

                g_x_0_xxzz_xxxxyy[k] = -g_x_0_xxz_xxxxyy[k] * ab_z + g_x_0_xxz_xxxxyyz[k];

                g_x_0_xxzz_xxxxyz[k] = -g_x_0_xxz_xxxxyz[k] * ab_z + g_x_0_xxz_xxxxyzz[k];

                g_x_0_xxzz_xxxxzz[k] = -g_x_0_xxz_xxxxzz[k] * ab_z + g_x_0_xxz_xxxxzzz[k];

                g_x_0_xxzz_xxxyyy[k] = -g_x_0_xxz_xxxyyy[k] * ab_z + g_x_0_xxz_xxxyyyz[k];

                g_x_0_xxzz_xxxyyz[k] = -g_x_0_xxz_xxxyyz[k] * ab_z + g_x_0_xxz_xxxyyzz[k];

                g_x_0_xxzz_xxxyzz[k] = -g_x_0_xxz_xxxyzz[k] * ab_z + g_x_0_xxz_xxxyzzz[k];

                g_x_0_xxzz_xxxzzz[k] = -g_x_0_xxz_xxxzzz[k] * ab_z + g_x_0_xxz_xxxzzzz[k];

                g_x_0_xxzz_xxyyyy[k] = -g_x_0_xxz_xxyyyy[k] * ab_z + g_x_0_xxz_xxyyyyz[k];

                g_x_0_xxzz_xxyyyz[k] = -g_x_0_xxz_xxyyyz[k] * ab_z + g_x_0_xxz_xxyyyzz[k];

                g_x_0_xxzz_xxyyzz[k] = -g_x_0_xxz_xxyyzz[k] * ab_z + g_x_0_xxz_xxyyzzz[k];

                g_x_0_xxzz_xxyzzz[k] = -g_x_0_xxz_xxyzzz[k] * ab_z + g_x_0_xxz_xxyzzzz[k];

                g_x_0_xxzz_xxzzzz[k] = -g_x_0_xxz_xxzzzz[k] * ab_z + g_x_0_xxz_xxzzzzz[k];

                g_x_0_xxzz_xyyyyy[k] = -g_x_0_xxz_xyyyyy[k] * ab_z + g_x_0_xxz_xyyyyyz[k];

                g_x_0_xxzz_xyyyyz[k] = -g_x_0_xxz_xyyyyz[k] * ab_z + g_x_0_xxz_xyyyyzz[k];

                g_x_0_xxzz_xyyyzz[k] = -g_x_0_xxz_xyyyzz[k] * ab_z + g_x_0_xxz_xyyyzzz[k];

                g_x_0_xxzz_xyyzzz[k] = -g_x_0_xxz_xyyzzz[k] * ab_z + g_x_0_xxz_xyyzzzz[k];

                g_x_0_xxzz_xyzzzz[k] = -g_x_0_xxz_xyzzzz[k] * ab_z + g_x_0_xxz_xyzzzzz[k];

                g_x_0_xxzz_xzzzzz[k] = -g_x_0_xxz_xzzzzz[k] * ab_z + g_x_0_xxz_xzzzzzz[k];

                g_x_0_xxzz_yyyyyy[k] = -g_x_0_xxz_yyyyyy[k] * ab_z + g_x_0_xxz_yyyyyyz[k];

                g_x_0_xxzz_yyyyyz[k] = -g_x_0_xxz_yyyyyz[k] * ab_z + g_x_0_xxz_yyyyyzz[k];

                g_x_0_xxzz_yyyyzz[k] = -g_x_0_xxz_yyyyzz[k] * ab_z + g_x_0_xxz_yyyyzzz[k];

                g_x_0_xxzz_yyyzzz[k] = -g_x_0_xxz_yyyzzz[k] * ab_z + g_x_0_xxz_yyyzzzz[k];

                g_x_0_xxzz_yyzzzz[k] = -g_x_0_xxz_yyzzzz[k] * ab_z + g_x_0_xxz_yyzzzzz[k];

                g_x_0_xxzz_yzzzzz[k] = -g_x_0_xxz_yzzzzz[k] * ab_z + g_x_0_xxz_yzzzzzz[k];

                g_x_0_xxzz_zzzzzz[k] = -g_x_0_xxz_zzzzzz[k] * ab_z + g_x_0_xxz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyy_xxxxxx, g_x_0_xyy_xxxxxxy, g_x_0_xyy_xxxxxy, g_x_0_xyy_xxxxxyy, g_x_0_xyy_xxxxxyz, g_x_0_xyy_xxxxxz, g_x_0_xyy_xxxxyy, g_x_0_xyy_xxxxyyy, g_x_0_xyy_xxxxyyz, g_x_0_xyy_xxxxyz, g_x_0_xyy_xxxxyzz, g_x_0_xyy_xxxxzz, g_x_0_xyy_xxxyyy, g_x_0_xyy_xxxyyyy, g_x_0_xyy_xxxyyyz, g_x_0_xyy_xxxyyz, g_x_0_xyy_xxxyyzz, g_x_0_xyy_xxxyzz, g_x_0_xyy_xxxyzzz, g_x_0_xyy_xxxzzz, g_x_0_xyy_xxyyyy, g_x_0_xyy_xxyyyyy, g_x_0_xyy_xxyyyyz, g_x_0_xyy_xxyyyz, g_x_0_xyy_xxyyyzz, g_x_0_xyy_xxyyzz, g_x_0_xyy_xxyyzzz, g_x_0_xyy_xxyzzz, g_x_0_xyy_xxyzzzz, g_x_0_xyy_xxzzzz, g_x_0_xyy_xyyyyy, g_x_0_xyy_xyyyyyy, g_x_0_xyy_xyyyyyz, g_x_0_xyy_xyyyyz, g_x_0_xyy_xyyyyzz, g_x_0_xyy_xyyyzz, g_x_0_xyy_xyyyzzz, g_x_0_xyy_xyyzzz, g_x_0_xyy_xyyzzzz, g_x_0_xyy_xyzzzz, g_x_0_xyy_xyzzzzz, g_x_0_xyy_xzzzzz, g_x_0_xyy_yyyyyy, g_x_0_xyy_yyyyyyy, g_x_0_xyy_yyyyyyz, g_x_0_xyy_yyyyyz, g_x_0_xyy_yyyyyzz, g_x_0_xyy_yyyyzz, g_x_0_xyy_yyyyzzz, g_x_0_xyy_yyyzzz, g_x_0_xyy_yyyzzzz, g_x_0_xyy_yyzzzz, g_x_0_xyy_yyzzzzz, g_x_0_xyy_yzzzzz, g_x_0_xyy_yzzzzzz, g_x_0_xyy_zzzzzz, g_x_0_xyyy_xxxxxx, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxxz, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxzz, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxzzz, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxzzzz, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xzzzzz, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_xxxxxx[k] = -g_x_0_xyy_xxxxxx[k] * ab_y + g_x_0_xyy_xxxxxxy[k];

                g_x_0_xyyy_xxxxxy[k] = -g_x_0_xyy_xxxxxy[k] * ab_y + g_x_0_xyy_xxxxxyy[k];

                g_x_0_xyyy_xxxxxz[k] = -g_x_0_xyy_xxxxxz[k] * ab_y + g_x_0_xyy_xxxxxyz[k];

                g_x_0_xyyy_xxxxyy[k] = -g_x_0_xyy_xxxxyy[k] * ab_y + g_x_0_xyy_xxxxyyy[k];

                g_x_0_xyyy_xxxxyz[k] = -g_x_0_xyy_xxxxyz[k] * ab_y + g_x_0_xyy_xxxxyyz[k];

                g_x_0_xyyy_xxxxzz[k] = -g_x_0_xyy_xxxxzz[k] * ab_y + g_x_0_xyy_xxxxyzz[k];

                g_x_0_xyyy_xxxyyy[k] = -g_x_0_xyy_xxxyyy[k] * ab_y + g_x_0_xyy_xxxyyyy[k];

                g_x_0_xyyy_xxxyyz[k] = -g_x_0_xyy_xxxyyz[k] * ab_y + g_x_0_xyy_xxxyyyz[k];

                g_x_0_xyyy_xxxyzz[k] = -g_x_0_xyy_xxxyzz[k] * ab_y + g_x_0_xyy_xxxyyzz[k];

                g_x_0_xyyy_xxxzzz[k] = -g_x_0_xyy_xxxzzz[k] * ab_y + g_x_0_xyy_xxxyzzz[k];

                g_x_0_xyyy_xxyyyy[k] = -g_x_0_xyy_xxyyyy[k] * ab_y + g_x_0_xyy_xxyyyyy[k];

                g_x_0_xyyy_xxyyyz[k] = -g_x_0_xyy_xxyyyz[k] * ab_y + g_x_0_xyy_xxyyyyz[k];

                g_x_0_xyyy_xxyyzz[k] = -g_x_0_xyy_xxyyzz[k] * ab_y + g_x_0_xyy_xxyyyzz[k];

                g_x_0_xyyy_xxyzzz[k] = -g_x_0_xyy_xxyzzz[k] * ab_y + g_x_0_xyy_xxyyzzz[k];

                g_x_0_xyyy_xxzzzz[k] = -g_x_0_xyy_xxzzzz[k] * ab_y + g_x_0_xyy_xxyzzzz[k];

                g_x_0_xyyy_xyyyyy[k] = -g_x_0_xyy_xyyyyy[k] * ab_y + g_x_0_xyy_xyyyyyy[k];

                g_x_0_xyyy_xyyyyz[k] = -g_x_0_xyy_xyyyyz[k] * ab_y + g_x_0_xyy_xyyyyyz[k];

                g_x_0_xyyy_xyyyzz[k] = -g_x_0_xyy_xyyyzz[k] * ab_y + g_x_0_xyy_xyyyyzz[k];

                g_x_0_xyyy_xyyzzz[k] = -g_x_0_xyy_xyyzzz[k] * ab_y + g_x_0_xyy_xyyyzzz[k];

                g_x_0_xyyy_xyzzzz[k] = -g_x_0_xyy_xyzzzz[k] * ab_y + g_x_0_xyy_xyyzzzz[k];

                g_x_0_xyyy_xzzzzz[k] = -g_x_0_xyy_xzzzzz[k] * ab_y + g_x_0_xyy_xyzzzzz[k];

                g_x_0_xyyy_yyyyyy[k] = -g_x_0_xyy_yyyyyy[k] * ab_y + g_x_0_xyy_yyyyyyy[k];

                g_x_0_xyyy_yyyyyz[k] = -g_x_0_xyy_yyyyyz[k] * ab_y + g_x_0_xyy_yyyyyyz[k];

                g_x_0_xyyy_yyyyzz[k] = -g_x_0_xyy_yyyyzz[k] * ab_y + g_x_0_xyy_yyyyyzz[k];

                g_x_0_xyyy_yyyzzz[k] = -g_x_0_xyy_yyyzzz[k] * ab_y + g_x_0_xyy_yyyyzzz[k];

                g_x_0_xyyy_yyzzzz[k] = -g_x_0_xyy_yyzzzz[k] * ab_y + g_x_0_xyy_yyyzzzz[k];

                g_x_0_xyyy_yzzzzz[k] = -g_x_0_xyy_yzzzzz[k] * ab_y + g_x_0_xyy_yyzzzzz[k];

                g_x_0_xyyy_zzzzzz[k] = -g_x_0_xyy_zzzzzz[k] * ab_y + g_x_0_xyy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyz_xxxxxx, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxxz, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxzz, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxzzz, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxzzzz, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xzzzzz, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_zzzzzz, g_x_0_xyz_xxxxxx, g_x_0_xyz_xxxxxxy, g_x_0_xyz_xxxxxy, g_x_0_xyz_xxxxxyy, g_x_0_xyz_xxxxxyz, g_x_0_xyz_xxxxxz, g_x_0_xyz_xxxxyy, g_x_0_xyz_xxxxyyy, g_x_0_xyz_xxxxyyz, g_x_0_xyz_xxxxyz, g_x_0_xyz_xxxxyzz, g_x_0_xyz_xxxxzz, g_x_0_xyz_xxxyyy, g_x_0_xyz_xxxyyyy, g_x_0_xyz_xxxyyyz, g_x_0_xyz_xxxyyz, g_x_0_xyz_xxxyyzz, g_x_0_xyz_xxxyzz, g_x_0_xyz_xxxyzzz, g_x_0_xyz_xxxzzz, g_x_0_xyz_xxyyyy, g_x_0_xyz_xxyyyyy, g_x_0_xyz_xxyyyyz, g_x_0_xyz_xxyyyz, g_x_0_xyz_xxyyyzz, g_x_0_xyz_xxyyzz, g_x_0_xyz_xxyyzzz, g_x_0_xyz_xxyzzz, g_x_0_xyz_xxyzzzz, g_x_0_xyz_xxzzzz, g_x_0_xyz_xyyyyy, g_x_0_xyz_xyyyyyy, g_x_0_xyz_xyyyyyz, g_x_0_xyz_xyyyyz, g_x_0_xyz_xyyyyzz, g_x_0_xyz_xyyyzz, g_x_0_xyz_xyyyzzz, g_x_0_xyz_xyyzzz, g_x_0_xyz_xyyzzzz, g_x_0_xyz_xyzzzz, g_x_0_xyz_xyzzzzz, g_x_0_xyz_xzzzzz, g_x_0_xyz_yyyyyy, g_x_0_xyz_yyyyyyy, g_x_0_xyz_yyyyyyz, g_x_0_xyz_yyyyyz, g_x_0_xyz_yyyyyzz, g_x_0_xyz_yyyyzz, g_x_0_xyz_yyyyzzz, g_x_0_xyz_yyyzzz, g_x_0_xyz_yyyzzzz, g_x_0_xyz_yyzzzz, g_x_0_xyz_yyzzzzz, g_x_0_xyz_yzzzzz, g_x_0_xyz_yzzzzzz, g_x_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_xxxxxx[k] = -g_x_0_xyz_xxxxxx[k] * ab_y + g_x_0_xyz_xxxxxxy[k];

                g_x_0_xyyz_xxxxxy[k] = -g_x_0_xyz_xxxxxy[k] * ab_y + g_x_0_xyz_xxxxxyy[k];

                g_x_0_xyyz_xxxxxz[k] = -g_x_0_xyz_xxxxxz[k] * ab_y + g_x_0_xyz_xxxxxyz[k];

                g_x_0_xyyz_xxxxyy[k] = -g_x_0_xyz_xxxxyy[k] * ab_y + g_x_0_xyz_xxxxyyy[k];

                g_x_0_xyyz_xxxxyz[k] = -g_x_0_xyz_xxxxyz[k] * ab_y + g_x_0_xyz_xxxxyyz[k];

                g_x_0_xyyz_xxxxzz[k] = -g_x_0_xyz_xxxxzz[k] * ab_y + g_x_0_xyz_xxxxyzz[k];

                g_x_0_xyyz_xxxyyy[k] = -g_x_0_xyz_xxxyyy[k] * ab_y + g_x_0_xyz_xxxyyyy[k];

                g_x_0_xyyz_xxxyyz[k] = -g_x_0_xyz_xxxyyz[k] * ab_y + g_x_0_xyz_xxxyyyz[k];

                g_x_0_xyyz_xxxyzz[k] = -g_x_0_xyz_xxxyzz[k] * ab_y + g_x_0_xyz_xxxyyzz[k];

                g_x_0_xyyz_xxxzzz[k] = -g_x_0_xyz_xxxzzz[k] * ab_y + g_x_0_xyz_xxxyzzz[k];

                g_x_0_xyyz_xxyyyy[k] = -g_x_0_xyz_xxyyyy[k] * ab_y + g_x_0_xyz_xxyyyyy[k];

                g_x_0_xyyz_xxyyyz[k] = -g_x_0_xyz_xxyyyz[k] * ab_y + g_x_0_xyz_xxyyyyz[k];

                g_x_0_xyyz_xxyyzz[k] = -g_x_0_xyz_xxyyzz[k] * ab_y + g_x_0_xyz_xxyyyzz[k];

                g_x_0_xyyz_xxyzzz[k] = -g_x_0_xyz_xxyzzz[k] * ab_y + g_x_0_xyz_xxyyzzz[k];

                g_x_0_xyyz_xxzzzz[k] = -g_x_0_xyz_xxzzzz[k] * ab_y + g_x_0_xyz_xxyzzzz[k];

                g_x_0_xyyz_xyyyyy[k] = -g_x_0_xyz_xyyyyy[k] * ab_y + g_x_0_xyz_xyyyyyy[k];

                g_x_0_xyyz_xyyyyz[k] = -g_x_0_xyz_xyyyyz[k] * ab_y + g_x_0_xyz_xyyyyyz[k];

                g_x_0_xyyz_xyyyzz[k] = -g_x_0_xyz_xyyyzz[k] * ab_y + g_x_0_xyz_xyyyyzz[k];

                g_x_0_xyyz_xyyzzz[k] = -g_x_0_xyz_xyyzzz[k] * ab_y + g_x_0_xyz_xyyyzzz[k];

                g_x_0_xyyz_xyzzzz[k] = -g_x_0_xyz_xyzzzz[k] * ab_y + g_x_0_xyz_xyyzzzz[k];

                g_x_0_xyyz_xzzzzz[k] = -g_x_0_xyz_xzzzzz[k] * ab_y + g_x_0_xyz_xyzzzzz[k];

                g_x_0_xyyz_yyyyyy[k] = -g_x_0_xyz_yyyyyy[k] * ab_y + g_x_0_xyz_yyyyyyy[k];

                g_x_0_xyyz_yyyyyz[k] = -g_x_0_xyz_yyyyyz[k] * ab_y + g_x_0_xyz_yyyyyyz[k];

                g_x_0_xyyz_yyyyzz[k] = -g_x_0_xyz_yyyyzz[k] * ab_y + g_x_0_xyz_yyyyyzz[k];

                g_x_0_xyyz_yyyzzz[k] = -g_x_0_xyz_yyyzzz[k] * ab_y + g_x_0_xyz_yyyyzzz[k];

                g_x_0_xyyz_yyzzzz[k] = -g_x_0_xyz_yyzzzz[k] * ab_y + g_x_0_xyz_yyyzzzz[k];

                g_x_0_xyyz_yzzzzz[k] = -g_x_0_xyz_yzzzzz[k] * ab_y + g_x_0_xyz_yyzzzzz[k];

                g_x_0_xyyz_zzzzzz[k] = -g_x_0_xyz_zzzzzz[k] * ab_y + g_x_0_xyz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzz_xxxxxx, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxxz, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxzz, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxzzz, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxzzzz, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xzzzzz, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_zzzzzz, g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxxy, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxyy, g_x_0_xzz_xxxxxyz, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyyy, g_x_0_xzz_xxxxyyz, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxyzz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyyy, g_x_0_xzz_xxxyyyz, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyyzz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxyzzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyyy, g_x_0_xzz_xxyyyyz, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyyzz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyyzzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxyzzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyyy, g_x_0_xzz_xyyyyyz, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyyzz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyyzzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyyzzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xyzzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyyy, g_x_0_xzz_yyyyyyz, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyyzz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyyzzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyyzzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yyzzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_yzzzzzz, g_x_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_xxxxxx[k] = -g_x_0_xzz_xxxxxx[k] * ab_y + g_x_0_xzz_xxxxxxy[k];

                g_x_0_xyzz_xxxxxy[k] = -g_x_0_xzz_xxxxxy[k] * ab_y + g_x_0_xzz_xxxxxyy[k];

                g_x_0_xyzz_xxxxxz[k] = -g_x_0_xzz_xxxxxz[k] * ab_y + g_x_0_xzz_xxxxxyz[k];

                g_x_0_xyzz_xxxxyy[k] = -g_x_0_xzz_xxxxyy[k] * ab_y + g_x_0_xzz_xxxxyyy[k];

                g_x_0_xyzz_xxxxyz[k] = -g_x_0_xzz_xxxxyz[k] * ab_y + g_x_0_xzz_xxxxyyz[k];

                g_x_0_xyzz_xxxxzz[k] = -g_x_0_xzz_xxxxzz[k] * ab_y + g_x_0_xzz_xxxxyzz[k];

                g_x_0_xyzz_xxxyyy[k] = -g_x_0_xzz_xxxyyy[k] * ab_y + g_x_0_xzz_xxxyyyy[k];

                g_x_0_xyzz_xxxyyz[k] = -g_x_0_xzz_xxxyyz[k] * ab_y + g_x_0_xzz_xxxyyyz[k];

                g_x_0_xyzz_xxxyzz[k] = -g_x_0_xzz_xxxyzz[k] * ab_y + g_x_0_xzz_xxxyyzz[k];

                g_x_0_xyzz_xxxzzz[k] = -g_x_0_xzz_xxxzzz[k] * ab_y + g_x_0_xzz_xxxyzzz[k];

                g_x_0_xyzz_xxyyyy[k] = -g_x_0_xzz_xxyyyy[k] * ab_y + g_x_0_xzz_xxyyyyy[k];

                g_x_0_xyzz_xxyyyz[k] = -g_x_0_xzz_xxyyyz[k] * ab_y + g_x_0_xzz_xxyyyyz[k];

                g_x_0_xyzz_xxyyzz[k] = -g_x_0_xzz_xxyyzz[k] * ab_y + g_x_0_xzz_xxyyyzz[k];

                g_x_0_xyzz_xxyzzz[k] = -g_x_0_xzz_xxyzzz[k] * ab_y + g_x_0_xzz_xxyyzzz[k];

                g_x_0_xyzz_xxzzzz[k] = -g_x_0_xzz_xxzzzz[k] * ab_y + g_x_0_xzz_xxyzzzz[k];

                g_x_0_xyzz_xyyyyy[k] = -g_x_0_xzz_xyyyyy[k] * ab_y + g_x_0_xzz_xyyyyyy[k];

                g_x_0_xyzz_xyyyyz[k] = -g_x_0_xzz_xyyyyz[k] * ab_y + g_x_0_xzz_xyyyyyz[k];

                g_x_0_xyzz_xyyyzz[k] = -g_x_0_xzz_xyyyzz[k] * ab_y + g_x_0_xzz_xyyyyzz[k];

                g_x_0_xyzz_xyyzzz[k] = -g_x_0_xzz_xyyzzz[k] * ab_y + g_x_0_xzz_xyyyzzz[k];

                g_x_0_xyzz_xyzzzz[k] = -g_x_0_xzz_xyzzzz[k] * ab_y + g_x_0_xzz_xyyzzzz[k];

                g_x_0_xyzz_xzzzzz[k] = -g_x_0_xzz_xzzzzz[k] * ab_y + g_x_0_xzz_xyzzzzz[k];

                g_x_0_xyzz_yyyyyy[k] = -g_x_0_xzz_yyyyyy[k] * ab_y + g_x_0_xzz_yyyyyyy[k];

                g_x_0_xyzz_yyyyyz[k] = -g_x_0_xzz_yyyyyz[k] * ab_y + g_x_0_xzz_yyyyyyz[k];

                g_x_0_xyzz_yyyyzz[k] = -g_x_0_xzz_yyyyzz[k] * ab_y + g_x_0_xzz_yyyyyzz[k];

                g_x_0_xyzz_yyyzzz[k] = -g_x_0_xzz_yyyzzz[k] * ab_y + g_x_0_xzz_yyyyzzz[k];

                g_x_0_xyzz_yyzzzz[k] = -g_x_0_xzz_yyzzzz[k] * ab_y + g_x_0_xzz_yyyzzzz[k];

                g_x_0_xyzz_yzzzzz[k] = -g_x_0_xzz_yzzzzz[k] * ab_y + g_x_0_xzz_yyzzzzz[k];

                g_x_0_xyzz_zzzzzz[k] = -g_x_0_xzz_zzzzzz[k] * ab_y + g_x_0_xzz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxxz, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxyz, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxxzz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyyz, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxyzz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxxzzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyyz, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyyzz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxyzzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxxzzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyyz, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyyzz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyyzzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxyzzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xxzzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyyz, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyyzz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyyzzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyyzzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xyzzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_xzzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyyz, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyyzz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyyzzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyyzzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yyzzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_yzzzzzz, g_x_0_xzz_zzzzzz, g_x_0_xzz_zzzzzzz, g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_xxxxxx[k] = -g_x_0_xzz_xxxxxx[k] * ab_z + g_x_0_xzz_xxxxxxz[k];

                g_x_0_xzzz_xxxxxy[k] = -g_x_0_xzz_xxxxxy[k] * ab_z + g_x_0_xzz_xxxxxyz[k];

                g_x_0_xzzz_xxxxxz[k] = -g_x_0_xzz_xxxxxz[k] * ab_z + g_x_0_xzz_xxxxxzz[k];

                g_x_0_xzzz_xxxxyy[k] = -g_x_0_xzz_xxxxyy[k] * ab_z + g_x_0_xzz_xxxxyyz[k];

                g_x_0_xzzz_xxxxyz[k] = -g_x_0_xzz_xxxxyz[k] * ab_z + g_x_0_xzz_xxxxyzz[k];

                g_x_0_xzzz_xxxxzz[k] = -g_x_0_xzz_xxxxzz[k] * ab_z + g_x_0_xzz_xxxxzzz[k];

                g_x_0_xzzz_xxxyyy[k] = -g_x_0_xzz_xxxyyy[k] * ab_z + g_x_0_xzz_xxxyyyz[k];

                g_x_0_xzzz_xxxyyz[k] = -g_x_0_xzz_xxxyyz[k] * ab_z + g_x_0_xzz_xxxyyzz[k];

                g_x_0_xzzz_xxxyzz[k] = -g_x_0_xzz_xxxyzz[k] * ab_z + g_x_0_xzz_xxxyzzz[k];

                g_x_0_xzzz_xxxzzz[k] = -g_x_0_xzz_xxxzzz[k] * ab_z + g_x_0_xzz_xxxzzzz[k];

                g_x_0_xzzz_xxyyyy[k] = -g_x_0_xzz_xxyyyy[k] * ab_z + g_x_0_xzz_xxyyyyz[k];

                g_x_0_xzzz_xxyyyz[k] = -g_x_0_xzz_xxyyyz[k] * ab_z + g_x_0_xzz_xxyyyzz[k];

                g_x_0_xzzz_xxyyzz[k] = -g_x_0_xzz_xxyyzz[k] * ab_z + g_x_0_xzz_xxyyzzz[k];

                g_x_0_xzzz_xxyzzz[k] = -g_x_0_xzz_xxyzzz[k] * ab_z + g_x_0_xzz_xxyzzzz[k];

                g_x_0_xzzz_xxzzzz[k] = -g_x_0_xzz_xxzzzz[k] * ab_z + g_x_0_xzz_xxzzzzz[k];

                g_x_0_xzzz_xyyyyy[k] = -g_x_0_xzz_xyyyyy[k] * ab_z + g_x_0_xzz_xyyyyyz[k];

                g_x_0_xzzz_xyyyyz[k] = -g_x_0_xzz_xyyyyz[k] * ab_z + g_x_0_xzz_xyyyyzz[k];

                g_x_0_xzzz_xyyyzz[k] = -g_x_0_xzz_xyyyzz[k] * ab_z + g_x_0_xzz_xyyyzzz[k];

                g_x_0_xzzz_xyyzzz[k] = -g_x_0_xzz_xyyzzz[k] * ab_z + g_x_0_xzz_xyyzzzz[k];

                g_x_0_xzzz_xyzzzz[k] = -g_x_0_xzz_xyzzzz[k] * ab_z + g_x_0_xzz_xyzzzzz[k];

                g_x_0_xzzz_xzzzzz[k] = -g_x_0_xzz_xzzzzz[k] * ab_z + g_x_0_xzz_xzzzzzz[k];

                g_x_0_xzzz_yyyyyy[k] = -g_x_0_xzz_yyyyyy[k] * ab_z + g_x_0_xzz_yyyyyyz[k];

                g_x_0_xzzz_yyyyyz[k] = -g_x_0_xzz_yyyyyz[k] * ab_z + g_x_0_xzz_yyyyyzz[k];

                g_x_0_xzzz_yyyyzz[k] = -g_x_0_xzz_yyyyzz[k] * ab_z + g_x_0_xzz_yyyyzzz[k];

                g_x_0_xzzz_yyyzzz[k] = -g_x_0_xzz_yyyzzz[k] * ab_z + g_x_0_xzz_yyyzzzz[k];

                g_x_0_xzzz_yyzzzz[k] = -g_x_0_xzz_yyzzzz[k] * ab_z + g_x_0_xzz_yyzzzzz[k];

                g_x_0_xzzz_yzzzzz[k] = -g_x_0_xzz_yzzzzz[k] * ab_z + g_x_0_xzz_yzzzzzz[k];

                g_x_0_xzzz_zzzzzz[k] = -g_x_0_xzz_zzzzzz[k] * ab_z + g_x_0_xzz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_xxxxxx, g_x_0_yyy_xxxxxxy, g_x_0_yyy_xxxxxy, g_x_0_yyy_xxxxxyy, g_x_0_yyy_xxxxxyz, g_x_0_yyy_xxxxxz, g_x_0_yyy_xxxxyy, g_x_0_yyy_xxxxyyy, g_x_0_yyy_xxxxyyz, g_x_0_yyy_xxxxyz, g_x_0_yyy_xxxxyzz, g_x_0_yyy_xxxxzz, g_x_0_yyy_xxxyyy, g_x_0_yyy_xxxyyyy, g_x_0_yyy_xxxyyyz, g_x_0_yyy_xxxyyz, g_x_0_yyy_xxxyyzz, g_x_0_yyy_xxxyzz, g_x_0_yyy_xxxyzzz, g_x_0_yyy_xxxzzz, g_x_0_yyy_xxyyyy, g_x_0_yyy_xxyyyyy, g_x_0_yyy_xxyyyyz, g_x_0_yyy_xxyyyz, g_x_0_yyy_xxyyyzz, g_x_0_yyy_xxyyzz, g_x_0_yyy_xxyyzzz, g_x_0_yyy_xxyzzz, g_x_0_yyy_xxyzzzz, g_x_0_yyy_xxzzzz, g_x_0_yyy_xyyyyy, g_x_0_yyy_xyyyyyy, g_x_0_yyy_xyyyyyz, g_x_0_yyy_xyyyyz, g_x_0_yyy_xyyyyzz, g_x_0_yyy_xyyyzz, g_x_0_yyy_xyyyzzz, g_x_0_yyy_xyyzzz, g_x_0_yyy_xyyzzzz, g_x_0_yyy_xyzzzz, g_x_0_yyy_xyzzzzz, g_x_0_yyy_xzzzzz, g_x_0_yyy_yyyyyy, g_x_0_yyy_yyyyyyy, g_x_0_yyy_yyyyyyz, g_x_0_yyy_yyyyyz, g_x_0_yyy_yyyyyzz, g_x_0_yyy_yyyyzz, g_x_0_yyy_yyyyzzz, g_x_0_yyy_yyyzzz, g_x_0_yyy_yyyzzzz, g_x_0_yyy_yyzzzz, g_x_0_yyy_yyzzzzz, g_x_0_yyy_yzzzzz, g_x_0_yyy_yzzzzzz, g_x_0_yyy_zzzzzz, g_x_0_yyyy_xxxxxx, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxxz, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxzz, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxzzz, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxzzzz, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xzzzzz, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_xxxxxx[k] = -g_x_0_yyy_xxxxxx[k] * ab_y + g_x_0_yyy_xxxxxxy[k];

                g_x_0_yyyy_xxxxxy[k] = -g_x_0_yyy_xxxxxy[k] * ab_y + g_x_0_yyy_xxxxxyy[k];

                g_x_0_yyyy_xxxxxz[k] = -g_x_0_yyy_xxxxxz[k] * ab_y + g_x_0_yyy_xxxxxyz[k];

                g_x_0_yyyy_xxxxyy[k] = -g_x_0_yyy_xxxxyy[k] * ab_y + g_x_0_yyy_xxxxyyy[k];

                g_x_0_yyyy_xxxxyz[k] = -g_x_0_yyy_xxxxyz[k] * ab_y + g_x_0_yyy_xxxxyyz[k];

                g_x_0_yyyy_xxxxzz[k] = -g_x_0_yyy_xxxxzz[k] * ab_y + g_x_0_yyy_xxxxyzz[k];

                g_x_0_yyyy_xxxyyy[k] = -g_x_0_yyy_xxxyyy[k] * ab_y + g_x_0_yyy_xxxyyyy[k];

                g_x_0_yyyy_xxxyyz[k] = -g_x_0_yyy_xxxyyz[k] * ab_y + g_x_0_yyy_xxxyyyz[k];

                g_x_0_yyyy_xxxyzz[k] = -g_x_0_yyy_xxxyzz[k] * ab_y + g_x_0_yyy_xxxyyzz[k];

                g_x_0_yyyy_xxxzzz[k] = -g_x_0_yyy_xxxzzz[k] * ab_y + g_x_0_yyy_xxxyzzz[k];

                g_x_0_yyyy_xxyyyy[k] = -g_x_0_yyy_xxyyyy[k] * ab_y + g_x_0_yyy_xxyyyyy[k];

                g_x_0_yyyy_xxyyyz[k] = -g_x_0_yyy_xxyyyz[k] * ab_y + g_x_0_yyy_xxyyyyz[k];

                g_x_0_yyyy_xxyyzz[k] = -g_x_0_yyy_xxyyzz[k] * ab_y + g_x_0_yyy_xxyyyzz[k];

                g_x_0_yyyy_xxyzzz[k] = -g_x_0_yyy_xxyzzz[k] * ab_y + g_x_0_yyy_xxyyzzz[k];

                g_x_0_yyyy_xxzzzz[k] = -g_x_0_yyy_xxzzzz[k] * ab_y + g_x_0_yyy_xxyzzzz[k];

                g_x_0_yyyy_xyyyyy[k] = -g_x_0_yyy_xyyyyy[k] * ab_y + g_x_0_yyy_xyyyyyy[k];

                g_x_0_yyyy_xyyyyz[k] = -g_x_0_yyy_xyyyyz[k] * ab_y + g_x_0_yyy_xyyyyyz[k];

                g_x_0_yyyy_xyyyzz[k] = -g_x_0_yyy_xyyyzz[k] * ab_y + g_x_0_yyy_xyyyyzz[k];

                g_x_0_yyyy_xyyzzz[k] = -g_x_0_yyy_xyyzzz[k] * ab_y + g_x_0_yyy_xyyyzzz[k];

                g_x_0_yyyy_xyzzzz[k] = -g_x_0_yyy_xyzzzz[k] * ab_y + g_x_0_yyy_xyyzzzz[k];

                g_x_0_yyyy_xzzzzz[k] = -g_x_0_yyy_xzzzzz[k] * ab_y + g_x_0_yyy_xyzzzzz[k];

                g_x_0_yyyy_yyyyyy[k] = -g_x_0_yyy_yyyyyy[k] * ab_y + g_x_0_yyy_yyyyyyy[k];

                g_x_0_yyyy_yyyyyz[k] = -g_x_0_yyy_yyyyyz[k] * ab_y + g_x_0_yyy_yyyyyyz[k];

                g_x_0_yyyy_yyyyzz[k] = -g_x_0_yyy_yyyyzz[k] * ab_y + g_x_0_yyy_yyyyyzz[k];

                g_x_0_yyyy_yyyzzz[k] = -g_x_0_yyy_yyyzzz[k] * ab_y + g_x_0_yyy_yyyyzzz[k];

                g_x_0_yyyy_yyzzzz[k] = -g_x_0_yyy_yyzzzz[k] * ab_y + g_x_0_yyy_yyyzzzz[k];

                g_x_0_yyyy_yzzzzz[k] = -g_x_0_yyy_yzzzzz[k] * ab_y + g_x_0_yyy_yyzzzzz[k];

                g_x_0_yyyy_zzzzzz[k] = -g_x_0_yyy_zzzzzz[k] * ab_y + g_x_0_yyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyz_xxxxxx, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxxz, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxzz, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxzzz, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxzzzz, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xzzzzz, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_zzzzzz, g_x_0_yyz_xxxxxx, g_x_0_yyz_xxxxxxy, g_x_0_yyz_xxxxxy, g_x_0_yyz_xxxxxyy, g_x_0_yyz_xxxxxyz, g_x_0_yyz_xxxxxz, g_x_0_yyz_xxxxyy, g_x_0_yyz_xxxxyyy, g_x_0_yyz_xxxxyyz, g_x_0_yyz_xxxxyz, g_x_0_yyz_xxxxyzz, g_x_0_yyz_xxxxzz, g_x_0_yyz_xxxyyy, g_x_0_yyz_xxxyyyy, g_x_0_yyz_xxxyyyz, g_x_0_yyz_xxxyyz, g_x_0_yyz_xxxyyzz, g_x_0_yyz_xxxyzz, g_x_0_yyz_xxxyzzz, g_x_0_yyz_xxxzzz, g_x_0_yyz_xxyyyy, g_x_0_yyz_xxyyyyy, g_x_0_yyz_xxyyyyz, g_x_0_yyz_xxyyyz, g_x_0_yyz_xxyyyzz, g_x_0_yyz_xxyyzz, g_x_0_yyz_xxyyzzz, g_x_0_yyz_xxyzzz, g_x_0_yyz_xxyzzzz, g_x_0_yyz_xxzzzz, g_x_0_yyz_xyyyyy, g_x_0_yyz_xyyyyyy, g_x_0_yyz_xyyyyyz, g_x_0_yyz_xyyyyz, g_x_0_yyz_xyyyyzz, g_x_0_yyz_xyyyzz, g_x_0_yyz_xyyyzzz, g_x_0_yyz_xyyzzz, g_x_0_yyz_xyyzzzz, g_x_0_yyz_xyzzzz, g_x_0_yyz_xyzzzzz, g_x_0_yyz_xzzzzz, g_x_0_yyz_yyyyyy, g_x_0_yyz_yyyyyyy, g_x_0_yyz_yyyyyyz, g_x_0_yyz_yyyyyz, g_x_0_yyz_yyyyyzz, g_x_0_yyz_yyyyzz, g_x_0_yyz_yyyyzzz, g_x_0_yyz_yyyzzz, g_x_0_yyz_yyyzzzz, g_x_0_yyz_yyzzzz, g_x_0_yyz_yyzzzzz, g_x_0_yyz_yzzzzz, g_x_0_yyz_yzzzzzz, g_x_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_xxxxxx[k] = -g_x_0_yyz_xxxxxx[k] * ab_y + g_x_0_yyz_xxxxxxy[k];

                g_x_0_yyyz_xxxxxy[k] = -g_x_0_yyz_xxxxxy[k] * ab_y + g_x_0_yyz_xxxxxyy[k];

                g_x_0_yyyz_xxxxxz[k] = -g_x_0_yyz_xxxxxz[k] * ab_y + g_x_0_yyz_xxxxxyz[k];

                g_x_0_yyyz_xxxxyy[k] = -g_x_0_yyz_xxxxyy[k] * ab_y + g_x_0_yyz_xxxxyyy[k];

                g_x_0_yyyz_xxxxyz[k] = -g_x_0_yyz_xxxxyz[k] * ab_y + g_x_0_yyz_xxxxyyz[k];

                g_x_0_yyyz_xxxxzz[k] = -g_x_0_yyz_xxxxzz[k] * ab_y + g_x_0_yyz_xxxxyzz[k];

                g_x_0_yyyz_xxxyyy[k] = -g_x_0_yyz_xxxyyy[k] * ab_y + g_x_0_yyz_xxxyyyy[k];

                g_x_0_yyyz_xxxyyz[k] = -g_x_0_yyz_xxxyyz[k] * ab_y + g_x_0_yyz_xxxyyyz[k];

                g_x_0_yyyz_xxxyzz[k] = -g_x_0_yyz_xxxyzz[k] * ab_y + g_x_0_yyz_xxxyyzz[k];

                g_x_0_yyyz_xxxzzz[k] = -g_x_0_yyz_xxxzzz[k] * ab_y + g_x_0_yyz_xxxyzzz[k];

                g_x_0_yyyz_xxyyyy[k] = -g_x_0_yyz_xxyyyy[k] * ab_y + g_x_0_yyz_xxyyyyy[k];

                g_x_0_yyyz_xxyyyz[k] = -g_x_0_yyz_xxyyyz[k] * ab_y + g_x_0_yyz_xxyyyyz[k];

                g_x_0_yyyz_xxyyzz[k] = -g_x_0_yyz_xxyyzz[k] * ab_y + g_x_0_yyz_xxyyyzz[k];

                g_x_0_yyyz_xxyzzz[k] = -g_x_0_yyz_xxyzzz[k] * ab_y + g_x_0_yyz_xxyyzzz[k];

                g_x_0_yyyz_xxzzzz[k] = -g_x_0_yyz_xxzzzz[k] * ab_y + g_x_0_yyz_xxyzzzz[k];

                g_x_0_yyyz_xyyyyy[k] = -g_x_0_yyz_xyyyyy[k] * ab_y + g_x_0_yyz_xyyyyyy[k];

                g_x_0_yyyz_xyyyyz[k] = -g_x_0_yyz_xyyyyz[k] * ab_y + g_x_0_yyz_xyyyyyz[k];

                g_x_0_yyyz_xyyyzz[k] = -g_x_0_yyz_xyyyzz[k] * ab_y + g_x_0_yyz_xyyyyzz[k];

                g_x_0_yyyz_xyyzzz[k] = -g_x_0_yyz_xyyzzz[k] * ab_y + g_x_0_yyz_xyyyzzz[k];

                g_x_0_yyyz_xyzzzz[k] = -g_x_0_yyz_xyzzzz[k] * ab_y + g_x_0_yyz_xyyzzzz[k];

                g_x_0_yyyz_xzzzzz[k] = -g_x_0_yyz_xzzzzz[k] * ab_y + g_x_0_yyz_xyzzzzz[k];

                g_x_0_yyyz_yyyyyy[k] = -g_x_0_yyz_yyyyyy[k] * ab_y + g_x_0_yyz_yyyyyyy[k];

                g_x_0_yyyz_yyyyyz[k] = -g_x_0_yyz_yyyyyz[k] * ab_y + g_x_0_yyz_yyyyyyz[k];

                g_x_0_yyyz_yyyyzz[k] = -g_x_0_yyz_yyyyzz[k] * ab_y + g_x_0_yyz_yyyyyzz[k];

                g_x_0_yyyz_yyyzzz[k] = -g_x_0_yyz_yyyzzz[k] * ab_y + g_x_0_yyz_yyyyzzz[k];

                g_x_0_yyyz_yyzzzz[k] = -g_x_0_yyz_yyzzzz[k] * ab_y + g_x_0_yyz_yyyzzzz[k];

                g_x_0_yyyz_yzzzzz[k] = -g_x_0_yyz_yzzzzz[k] * ab_y + g_x_0_yyz_yyzzzzz[k];

                g_x_0_yyyz_zzzzzz[k] = -g_x_0_yyz_zzzzzz[k] * ab_y + g_x_0_yyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzz_xxxxxx, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxxz, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxzz, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxzzz, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxzzzz, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xzzzzz, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_zzzzzz, g_x_0_yzz_xxxxxx, g_x_0_yzz_xxxxxxy, g_x_0_yzz_xxxxxy, g_x_0_yzz_xxxxxyy, g_x_0_yzz_xxxxxyz, g_x_0_yzz_xxxxxz, g_x_0_yzz_xxxxyy, g_x_0_yzz_xxxxyyy, g_x_0_yzz_xxxxyyz, g_x_0_yzz_xxxxyz, g_x_0_yzz_xxxxyzz, g_x_0_yzz_xxxxzz, g_x_0_yzz_xxxyyy, g_x_0_yzz_xxxyyyy, g_x_0_yzz_xxxyyyz, g_x_0_yzz_xxxyyz, g_x_0_yzz_xxxyyzz, g_x_0_yzz_xxxyzz, g_x_0_yzz_xxxyzzz, g_x_0_yzz_xxxzzz, g_x_0_yzz_xxyyyy, g_x_0_yzz_xxyyyyy, g_x_0_yzz_xxyyyyz, g_x_0_yzz_xxyyyz, g_x_0_yzz_xxyyyzz, g_x_0_yzz_xxyyzz, g_x_0_yzz_xxyyzzz, g_x_0_yzz_xxyzzz, g_x_0_yzz_xxyzzzz, g_x_0_yzz_xxzzzz, g_x_0_yzz_xyyyyy, g_x_0_yzz_xyyyyyy, g_x_0_yzz_xyyyyyz, g_x_0_yzz_xyyyyz, g_x_0_yzz_xyyyyzz, g_x_0_yzz_xyyyzz, g_x_0_yzz_xyyyzzz, g_x_0_yzz_xyyzzz, g_x_0_yzz_xyyzzzz, g_x_0_yzz_xyzzzz, g_x_0_yzz_xyzzzzz, g_x_0_yzz_xzzzzz, g_x_0_yzz_yyyyyy, g_x_0_yzz_yyyyyyy, g_x_0_yzz_yyyyyyz, g_x_0_yzz_yyyyyz, g_x_0_yzz_yyyyyzz, g_x_0_yzz_yyyyzz, g_x_0_yzz_yyyyzzz, g_x_0_yzz_yyyzzz, g_x_0_yzz_yyyzzzz, g_x_0_yzz_yyzzzz, g_x_0_yzz_yyzzzzz, g_x_0_yzz_yzzzzz, g_x_0_yzz_yzzzzzz, g_x_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_xxxxxx[k] = -g_x_0_yzz_xxxxxx[k] * ab_y + g_x_0_yzz_xxxxxxy[k];

                g_x_0_yyzz_xxxxxy[k] = -g_x_0_yzz_xxxxxy[k] * ab_y + g_x_0_yzz_xxxxxyy[k];

                g_x_0_yyzz_xxxxxz[k] = -g_x_0_yzz_xxxxxz[k] * ab_y + g_x_0_yzz_xxxxxyz[k];

                g_x_0_yyzz_xxxxyy[k] = -g_x_0_yzz_xxxxyy[k] * ab_y + g_x_0_yzz_xxxxyyy[k];

                g_x_0_yyzz_xxxxyz[k] = -g_x_0_yzz_xxxxyz[k] * ab_y + g_x_0_yzz_xxxxyyz[k];

                g_x_0_yyzz_xxxxzz[k] = -g_x_0_yzz_xxxxzz[k] * ab_y + g_x_0_yzz_xxxxyzz[k];

                g_x_0_yyzz_xxxyyy[k] = -g_x_0_yzz_xxxyyy[k] * ab_y + g_x_0_yzz_xxxyyyy[k];

                g_x_0_yyzz_xxxyyz[k] = -g_x_0_yzz_xxxyyz[k] * ab_y + g_x_0_yzz_xxxyyyz[k];

                g_x_0_yyzz_xxxyzz[k] = -g_x_0_yzz_xxxyzz[k] * ab_y + g_x_0_yzz_xxxyyzz[k];

                g_x_0_yyzz_xxxzzz[k] = -g_x_0_yzz_xxxzzz[k] * ab_y + g_x_0_yzz_xxxyzzz[k];

                g_x_0_yyzz_xxyyyy[k] = -g_x_0_yzz_xxyyyy[k] * ab_y + g_x_0_yzz_xxyyyyy[k];

                g_x_0_yyzz_xxyyyz[k] = -g_x_0_yzz_xxyyyz[k] * ab_y + g_x_0_yzz_xxyyyyz[k];

                g_x_0_yyzz_xxyyzz[k] = -g_x_0_yzz_xxyyzz[k] * ab_y + g_x_0_yzz_xxyyyzz[k];

                g_x_0_yyzz_xxyzzz[k] = -g_x_0_yzz_xxyzzz[k] * ab_y + g_x_0_yzz_xxyyzzz[k];

                g_x_0_yyzz_xxzzzz[k] = -g_x_0_yzz_xxzzzz[k] * ab_y + g_x_0_yzz_xxyzzzz[k];

                g_x_0_yyzz_xyyyyy[k] = -g_x_0_yzz_xyyyyy[k] * ab_y + g_x_0_yzz_xyyyyyy[k];

                g_x_0_yyzz_xyyyyz[k] = -g_x_0_yzz_xyyyyz[k] * ab_y + g_x_0_yzz_xyyyyyz[k];

                g_x_0_yyzz_xyyyzz[k] = -g_x_0_yzz_xyyyzz[k] * ab_y + g_x_0_yzz_xyyyyzz[k];

                g_x_0_yyzz_xyyzzz[k] = -g_x_0_yzz_xyyzzz[k] * ab_y + g_x_0_yzz_xyyyzzz[k];

                g_x_0_yyzz_xyzzzz[k] = -g_x_0_yzz_xyzzzz[k] * ab_y + g_x_0_yzz_xyyzzzz[k];

                g_x_0_yyzz_xzzzzz[k] = -g_x_0_yzz_xzzzzz[k] * ab_y + g_x_0_yzz_xyzzzzz[k];

                g_x_0_yyzz_yyyyyy[k] = -g_x_0_yzz_yyyyyy[k] * ab_y + g_x_0_yzz_yyyyyyy[k];

                g_x_0_yyzz_yyyyyz[k] = -g_x_0_yzz_yyyyyz[k] * ab_y + g_x_0_yzz_yyyyyyz[k];

                g_x_0_yyzz_yyyyzz[k] = -g_x_0_yzz_yyyyzz[k] * ab_y + g_x_0_yzz_yyyyyzz[k];

                g_x_0_yyzz_yyyzzz[k] = -g_x_0_yzz_yyyzzz[k] * ab_y + g_x_0_yzz_yyyyzzz[k];

                g_x_0_yyzz_yyzzzz[k] = -g_x_0_yzz_yyzzzz[k] * ab_y + g_x_0_yzz_yyyzzzz[k];

                g_x_0_yyzz_yzzzzz[k] = -g_x_0_yzz_yzzzzz[k] * ab_y + g_x_0_yzz_yyzzzzz[k];

                g_x_0_yyzz_zzzzzz[k] = -g_x_0_yzz_zzzzzz[k] * ab_y + g_x_0_yzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzz_xxxxxx, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxxz, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxzz, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxzzz, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxzzzz, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xzzzzz, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_zzzzzz, g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxxy, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxyy, g_x_0_zzz_xxxxxyz, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyyy, g_x_0_zzz_xxxxyyz, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxyzz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyyy, g_x_0_zzz_xxxyyyz, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyyzz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxyzzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyyy, g_x_0_zzz_xxyyyyz, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyyzz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyyzzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxyzzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyyy, g_x_0_zzz_xyyyyyz, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyyzz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyyzzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyyzzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xyzzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyyy, g_x_0_zzz_yyyyyyz, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyyzz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyyzzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyyzzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yyzzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_yzzzzzz, g_x_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_xxxxxx[k] = -g_x_0_zzz_xxxxxx[k] * ab_y + g_x_0_zzz_xxxxxxy[k];

                g_x_0_yzzz_xxxxxy[k] = -g_x_0_zzz_xxxxxy[k] * ab_y + g_x_0_zzz_xxxxxyy[k];

                g_x_0_yzzz_xxxxxz[k] = -g_x_0_zzz_xxxxxz[k] * ab_y + g_x_0_zzz_xxxxxyz[k];

                g_x_0_yzzz_xxxxyy[k] = -g_x_0_zzz_xxxxyy[k] * ab_y + g_x_0_zzz_xxxxyyy[k];

                g_x_0_yzzz_xxxxyz[k] = -g_x_0_zzz_xxxxyz[k] * ab_y + g_x_0_zzz_xxxxyyz[k];

                g_x_0_yzzz_xxxxzz[k] = -g_x_0_zzz_xxxxzz[k] * ab_y + g_x_0_zzz_xxxxyzz[k];

                g_x_0_yzzz_xxxyyy[k] = -g_x_0_zzz_xxxyyy[k] * ab_y + g_x_0_zzz_xxxyyyy[k];

                g_x_0_yzzz_xxxyyz[k] = -g_x_0_zzz_xxxyyz[k] * ab_y + g_x_0_zzz_xxxyyyz[k];

                g_x_0_yzzz_xxxyzz[k] = -g_x_0_zzz_xxxyzz[k] * ab_y + g_x_0_zzz_xxxyyzz[k];

                g_x_0_yzzz_xxxzzz[k] = -g_x_0_zzz_xxxzzz[k] * ab_y + g_x_0_zzz_xxxyzzz[k];

                g_x_0_yzzz_xxyyyy[k] = -g_x_0_zzz_xxyyyy[k] * ab_y + g_x_0_zzz_xxyyyyy[k];

                g_x_0_yzzz_xxyyyz[k] = -g_x_0_zzz_xxyyyz[k] * ab_y + g_x_0_zzz_xxyyyyz[k];

                g_x_0_yzzz_xxyyzz[k] = -g_x_0_zzz_xxyyzz[k] * ab_y + g_x_0_zzz_xxyyyzz[k];

                g_x_0_yzzz_xxyzzz[k] = -g_x_0_zzz_xxyzzz[k] * ab_y + g_x_0_zzz_xxyyzzz[k];

                g_x_0_yzzz_xxzzzz[k] = -g_x_0_zzz_xxzzzz[k] * ab_y + g_x_0_zzz_xxyzzzz[k];

                g_x_0_yzzz_xyyyyy[k] = -g_x_0_zzz_xyyyyy[k] * ab_y + g_x_0_zzz_xyyyyyy[k];

                g_x_0_yzzz_xyyyyz[k] = -g_x_0_zzz_xyyyyz[k] * ab_y + g_x_0_zzz_xyyyyyz[k];

                g_x_0_yzzz_xyyyzz[k] = -g_x_0_zzz_xyyyzz[k] * ab_y + g_x_0_zzz_xyyyyzz[k];

                g_x_0_yzzz_xyyzzz[k] = -g_x_0_zzz_xyyzzz[k] * ab_y + g_x_0_zzz_xyyyzzz[k];

                g_x_0_yzzz_xyzzzz[k] = -g_x_0_zzz_xyzzzz[k] * ab_y + g_x_0_zzz_xyyzzzz[k];

                g_x_0_yzzz_xzzzzz[k] = -g_x_0_zzz_xzzzzz[k] * ab_y + g_x_0_zzz_xyzzzzz[k];

                g_x_0_yzzz_yyyyyy[k] = -g_x_0_zzz_yyyyyy[k] * ab_y + g_x_0_zzz_yyyyyyy[k];

                g_x_0_yzzz_yyyyyz[k] = -g_x_0_zzz_yyyyyz[k] * ab_y + g_x_0_zzz_yyyyyyz[k];

                g_x_0_yzzz_yyyyzz[k] = -g_x_0_zzz_yyyyzz[k] * ab_y + g_x_0_zzz_yyyyyzz[k];

                g_x_0_yzzz_yyyzzz[k] = -g_x_0_zzz_yyyzzz[k] * ab_y + g_x_0_zzz_yyyyzzz[k];

                g_x_0_yzzz_yyzzzz[k] = -g_x_0_zzz_yyzzzz[k] * ab_y + g_x_0_zzz_yyyzzzz[k];

                g_x_0_yzzz_yzzzzz[k] = -g_x_0_zzz_yzzzzz[k] * ab_y + g_x_0_zzz_yyzzzzz[k];

                g_x_0_yzzz_zzzzzz[k] = -g_x_0_zzz_zzzzzz[k] * ab_y + g_x_0_zzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxxz, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxyz, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxxzz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyyz, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxyzz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxxzzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyyz, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyyzz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxyzzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxxzzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyyz, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyyzz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyyzzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxyzzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xxzzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyyz, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyyzz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyyzzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyyzzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xyzzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_xzzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyyz, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyyzz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyyzzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyyzzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yyzzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_yzzzzzz, g_x_0_zzz_zzzzzz, g_x_0_zzz_zzzzzzz, g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_xxxxxx[k] = -g_x_0_zzz_xxxxxx[k] * ab_z + g_x_0_zzz_xxxxxxz[k];

                g_x_0_zzzz_xxxxxy[k] = -g_x_0_zzz_xxxxxy[k] * ab_z + g_x_0_zzz_xxxxxyz[k];

                g_x_0_zzzz_xxxxxz[k] = -g_x_0_zzz_xxxxxz[k] * ab_z + g_x_0_zzz_xxxxxzz[k];

                g_x_0_zzzz_xxxxyy[k] = -g_x_0_zzz_xxxxyy[k] * ab_z + g_x_0_zzz_xxxxyyz[k];

                g_x_0_zzzz_xxxxyz[k] = -g_x_0_zzz_xxxxyz[k] * ab_z + g_x_0_zzz_xxxxyzz[k];

                g_x_0_zzzz_xxxxzz[k] = -g_x_0_zzz_xxxxzz[k] * ab_z + g_x_0_zzz_xxxxzzz[k];

                g_x_0_zzzz_xxxyyy[k] = -g_x_0_zzz_xxxyyy[k] * ab_z + g_x_0_zzz_xxxyyyz[k];

                g_x_0_zzzz_xxxyyz[k] = -g_x_0_zzz_xxxyyz[k] * ab_z + g_x_0_zzz_xxxyyzz[k];

                g_x_0_zzzz_xxxyzz[k] = -g_x_0_zzz_xxxyzz[k] * ab_z + g_x_0_zzz_xxxyzzz[k];

                g_x_0_zzzz_xxxzzz[k] = -g_x_0_zzz_xxxzzz[k] * ab_z + g_x_0_zzz_xxxzzzz[k];

                g_x_0_zzzz_xxyyyy[k] = -g_x_0_zzz_xxyyyy[k] * ab_z + g_x_0_zzz_xxyyyyz[k];

                g_x_0_zzzz_xxyyyz[k] = -g_x_0_zzz_xxyyyz[k] * ab_z + g_x_0_zzz_xxyyyzz[k];

                g_x_0_zzzz_xxyyzz[k] = -g_x_0_zzz_xxyyzz[k] * ab_z + g_x_0_zzz_xxyyzzz[k];

                g_x_0_zzzz_xxyzzz[k] = -g_x_0_zzz_xxyzzz[k] * ab_z + g_x_0_zzz_xxyzzzz[k];

                g_x_0_zzzz_xxzzzz[k] = -g_x_0_zzz_xxzzzz[k] * ab_z + g_x_0_zzz_xxzzzzz[k];

                g_x_0_zzzz_xyyyyy[k] = -g_x_0_zzz_xyyyyy[k] * ab_z + g_x_0_zzz_xyyyyyz[k];

                g_x_0_zzzz_xyyyyz[k] = -g_x_0_zzz_xyyyyz[k] * ab_z + g_x_0_zzz_xyyyyzz[k];

                g_x_0_zzzz_xyyyzz[k] = -g_x_0_zzz_xyyyzz[k] * ab_z + g_x_0_zzz_xyyyzzz[k];

                g_x_0_zzzz_xyyzzz[k] = -g_x_0_zzz_xyyzzz[k] * ab_z + g_x_0_zzz_xyyzzzz[k];

                g_x_0_zzzz_xyzzzz[k] = -g_x_0_zzz_xyzzzz[k] * ab_z + g_x_0_zzz_xyzzzzz[k];

                g_x_0_zzzz_xzzzzz[k] = -g_x_0_zzz_xzzzzz[k] * ab_z + g_x_0_zzz_xzzzzzz[k];

                g_x_0_zzzz_yyyyyy[k] = -g_x_0_zzz_yyyyyy[k] * ab_z + g_x_0_zzz_yyyyyyz[k];

                g_x_0_zzzz_yyyyyz[k] = -g_x_0_zzz_yyyyyz[k] * ab_z + g_x_0_zzz_yyyyyzz[k];

                g_x_0_zzzz_yyyyzz[k] = -g_x_0_zzz_yyyyzz[k] * ab_z + g_x_0_zzz_yyyyzzz[k];

                g_x_0_zzzz_yyyzzz[k] = -g_x_0_zzz_yyyzzz[k] * ab_z + g_x_0_zzz_yyyzzzz[k];

                g_x_0_zzzz_yyzzzz[k] = -g_x_0_zzz_yyzzzz[k] * ab_z + g_x_0_zzz_yyzzzzz[k];

                g_x_0_zzzz_yzzzzz[k] = -g_x_0_zzz_yzzzzz[k] * ab_z + g_x_0_zzz_yzzzzzz[k];

                g_x_0_zzzz_zzzzzz[k] = -g_x_0_zzz_zzzzzz[k] * ab_z + g_x_0_zzz_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_xxxxxx, g_y_0_xxx_xxxxxxx, g_y_0_xxx_xxxxxxy, g_y_0_xxx_xxxxxxz, g_y_0_xxx_xxxxxy, g_y_0_xxx_xxxxxyy, g_y_0_xxx_xxxxxyz, g_y_0_xxx_xxxxxz, g_y_0_xxx_xxxxxzz, g_y_0_xxx_xxxxyy, g_y_0_xxx_xxxxyyy, g_y_0_xxx_xxxxyyz, g_y_0_xxx_xxxxyz, g_y_0_xxx_xxxxyzz, g_y_0_xxx_xxxxzz, g_y_0_xxx_xxxxzzz, g_y_0_xxx_xxxyyy, g_y_0_xxx_xxxyyyy, g_y_0_xxx_xxxyyyz, g_y_0_xxx_xxxyyz, g_y_0_xxx_xxxyyzz, g_y_0_xxx_xxxyzz, g_y_0_xxx_xxxyzzz, g_y_0_xxx_xxxzzz, g_y_0_xxx_xxxzzzz, g_y_0_xxx_xxyyyy, g_y_0_xxx_xxyyyyy, g_y_0_xxx_xxyyyyz, g_y_0_xxx_xxyyyz, g_y_0_xxx_xxyyyzz, g_y_0_xxx_xxyyzz, g_y_0_xxx_xxyyzzz, g_y_0_xxx_xxyzzz, g_y_0_xxx_xxyzzzz, g_y_0_xxx_xxzzzz, g_y_0_xxx_xxzzzzz, g_y_0_xxx_xyyyyy, g_y_0_xxx_xyyyyyy, g_y_0_xxx_xyyyyyz, g_y_0_xxx_xyyyyz, g_y_0_xxx_xyyyyzz, g_y_0_xxx_xyyyzz, g_y_0_xxx_xyyyzzz, g_y_0_xxx_xyyzzz, g_y_0_xxx_xyyzzzz, g_y_0_xxx_xyzzzz, g_y_0_xxx_xyzzzzz, g_y_0_xxx_xzzzzz, g_y_0_xxx_xzzzzzz, g_y_0_xxx_yyyyyy, g_y_0_xxx_yyyyyz, g_y_0_xxx_yyyyzz, g_y_0_xxx_yyyzzz, g_y_0_xxx_yyzzzz, g_y_0_xxx_yzzzzz, g_y_0_xxx_zzzzzz, g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_yyyyyy, g_y_0_xxxx_yyyyyz, g_y_0_xxxx_yyyyzz, g_y_0_xxxx_yyyzzz, g_y_0_xxxx_yyzzzz, g_y_0_xxxx_yzzzzz, g_y_0_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_xxxxxx[k] = -g_y_0_xxx_xxxxxx[k] * ab_x + g_y_0_xxx_xxxxxxx[k];

                g_y_0_xxxx_xxxxxy[k] = -g_y_0_xxx_xxxxxy[k] * ab_x + g_y_0_xxx_xxxxxxy[k];

                g_y_0_xxxx_xxxxxz[k] = -g_y_0_xxx_xxxxxz[k] * ab_x + g_y_0_xxx_xxxxxxz[k];

                g_y_0_xxxx_xxxxyy[k] = -g_y_0_xxx_xxxxyy[k] * ab_x + g_y_0_xxx_xxxxxyy[k];

                g_y_0_xxxx_xxxxyz[k] = -g_y_0_xxx_xxxxyz[k] * ab_x + g_y_0_xxx_xxxxxyz[k];

                g_y_0_xxxx_xxxxzz[k] = -g_y_0_xxx_xxxxzz[k] * ab_x + g_y_0_xxx_xxxxxzz[k];

                g_y_0_xxxx_xxxyyy[k] = -g_y_0_xxx_xxxyyy[k] * ab_x + g_y_0_xxx_xxxxyyy[k];

                g_y_0_xxxx_xxxyyz[k] = -g_y_0_xxx_xxxyyz[k] * ab_x + g_y_0_xxx_xxxxyyz[k];

                g_y_0_xxxx_xxxyzz[k] = -g_y_0_xxx_xxxyzz[k] * ab_x + g_y_0_xxx_xxxxyzz[k];

                g_y_0_xxxx_xxxzzz[k] = -g_y_0_xxx_xxxzzz[k] * ab_x + g_y_0_xxx_xxxxzzz[k];

                g_y_0_xxxx_xxyyyy[k] = -g_y_0_xxx_xxyyyy[k] * ab_x + g_y_0_xxx_xxxyyyy[k];

                g_y_0_xxxx_xxyyyz[k] = -g_y_0_xxx_xxyyyz[k] * ab_x + g_y_0_xxx_xxxyyyz[k];

                g_y_0_xxxx_xxyyzz[k] = -g_y_0_xxx_xxyyzz[k] * ab_x + g_y_0_xxx_xxxyyzz[k];

                g_y_0_xxxx_xxyzzz[k] = -g_y_0_xxx_xxyzzz[k] * ab_x + g_y_0_xxx_xxxyzzz[k];

                g_y_0_xxxx_xxzzzz[k] = -g_y_0_xxx_xxzzzz[k] * ab_x + g_y_0_xxx_xxxzzzz[k];

                g_y_0_xxxx_xyyyyy[k] = -g_y_0_xxx_xyyyyy[k] * ab_x + g_y_0_xxx_xxyyyyy[k];

                g_y_0_xxxx_xyyyyz[k] = -g_y_0_xxx_xyyyyz[k] * ab_x + g_y_0_xxx_xxyyyyz[k];

                g_y_0_xxxx_xyyyzz[k] = -g_y_0_xxx_xyyyzz[k] * ab_x + g_y_0_xxx_xxyyyzz[k];

                g_y_0_xxxx_xyyzzz[k] = -g_y_0_xxx_xyyzzz[k] * ab_x + g_y_0_xxx_xxyyzzz[k];

                g_y_0_xxxx_xyzzzz[k] = -g_y_0_xxx_xyzzzz[k] * ab_x + g_y_0_xxx_xxyzzzz[k];

                g_y_0_xxxx_xzzzzz[k] = -g_y_0_xxx_xzzzzz[k] * ab_x + g_y_0_xxx_xxzzzzz[k];

                g_y_0_xxxx_yyyyyy[k] = -g_y_0_xxx_yyyyyy[k] * ab_x + g_y_0_xxx_xyyyyyy[k];

                g_y_0_xxxx_yyyyyz[k] = -g_y_0_xxx_yyyyyz[k] * ab_x + g_y_0_xxx_xyyyyyz[k];

                g_y_0_xxxx_yyyyzz[k] = -g_y_0_xxx_yyyyzz[k] * ab_x + g_y_0_xxx_xyyyyzz[k];

                g_y_0_xxxx_yyyzzz[k] = -g_y_0_xxx_yyyzzz[k] * ab_x + g_y_0_xxx_xyyyzzz[k];

                g_y_0_xxxx_yyzzzz[k] = -g_y_0_xxx_yyzzzz[k] * ab_x + g_y_0_xxx_xyyzzzz[k];

                g_y_0_xxxx_yzzzzz[k] = -g_y_0_xxx_yzzzzz[k] * ab_x + g_y_0_xxx_xyzzzzz[k];

                g_y_0_xxxx_zzzzzz[k] = -g_y_0_xxx_zzzzzz[k] * ab_x + g_y_0_xxx_xzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_yyyyyy, g_y_0_xxxy_yyyyyz, g_y_0_xxxy_yyyyzz, g_y_0_xxxy_yyyzzz, g_y_0_xxxy_yyzzzz, g_y_0_xxxy_yzzzzz, g_y_0_xxxy_zzzzzz, g_y_0_xxy_xxxxxx, g_y_0_xxy_xxxxxxx, g_y_0_xxy_xxxxxxy, g_y_0_xxy_xxxxxxz, g_y_0_xxy_xxxxxy, g_y_0_xxy_xxxxxyy, g_y_0_xxy_xxxxxyz, g_y_0_xxy_xxxxxz, g_y_0_xxy_xxxxxzz, g_y_0_xxy_xxxxyy, g_y_0_xxy_xxxxyyy, g_y_0_xxy_xxxxyyz, g_y_0_xxy_xxxxyz, g_y_0_xxy_xxxxyzz, g_y_0_xxy_xxxxzz, g_y_0_xxy_xxxxzzz, g_y_0_xxy_xxxyyy, g_y_0_xxy_xxxyyyy, g_y_0_xxy_xxxyyyz, g_y_0_xxy_xxxyyz, g_y_0_xxy_xxxyyzz, g_y_0_xxy_xxxyzz, g_y_0_xxy_xxxyzzz, g_y_0_xxy_xxxzzz, g_y_0_xxy_xxxzzzz, g_y_0_xxy_xxyyyy, g_y_0_xxy_xxyyyyy, g_y_0_xxy_xxyyyyz, g_y_0_xxy_xxyyyz, g_y_0_xxy_xxyyyzz, g_y_0_xxy_xxyyzz, g_y_0_xxy_xxyyzzz, g_y_0_xxy_xxyzzz, g_y_0_xxy_xxyzzzz, g_y_0_xxy_xxzzzz, g_y_0_xxy_xxzzzzz, g_y_0_xxy_xyyyyy, g_y_0_xxy_xyyyyyy, g_y_0_xxy_xyyyyyz, g_y_0_xxy_xyyyyz, g_y_0_xxy_xyyyyzz, g_y_0_xxy_xyyyzz, g_y_0_xxy_xyyyzzz, g_y_0_xxy_xyyzzz, g_y_0_xxy_xyyzzzz, g_y_0_xxy_xyzzzz, g_y_0_xxy_xyzzzzz, g_y_0_xxy_xzzzzz, g_y_0_xxy_xzzzzzz, g_y_0_xxy_yyyyyy, g_y_0_xxy_yyyyyz, g_y_0_xxy_yyyyzz, g_y_0_xxy_yyyzzz, g_y_0_xxy_yyzzzz, g_y_0_xxy_yzzzzz, g_y_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_xxxxxx[k] = -g_y_0_xxy_xxxxxx[k] * ab_x + g_y_0_xxy_xxxxxxx[k];

                g_y_0_xxxy_xxxxxy[k] = -g_y_0_xxy_xxxxxy[k] * ab_x + g_y_0_xxy_xxxxxxy[k];

                g_y_0_xxxy_xxxxxz[k] = -g_y_0_xxy_xxxxxz[k] * ab_x + g_y_0_xxy_xxxxxxz[k];

                g_y_0_xxxy_xxxxyy[k] = -g_y_0_xxy_xxxxyy[k] * ab_x + g_y_0_xxy_xxxxxyy[k];

                g_y_0_xxxy_xxxxyz[k] = -g_y_0_xxy_xxxxyz[k] * ab_x + g_y_0_xxy_xxxxxyz[k];

                g_y_0_xxxy_xxxxzz[k] = -g_y_0_xxy_xxxxzz[k] * ab_x + g_y_0_xxy_xxxxxzz[k];

                g_y_0_xxxy_xxxyyy[k] = -g_y_0_xxy_xxxyyy[k] * ab_x + g_y_0_xxy_xxxxyyy[k];

                g_y_0_xxxy_xxxyyz[k] = -g_y_0_xxy_xxxyyz[k] * ab_x + g_y_0_xxy_xxxxyyz[k];

                g_y_0_xxxy_xxxyzz[k] = -g_y_0_xxy_xxxyzz[k] * ab_x + g_y_0_xxy_xxxxyzz[k];

                g_y_0_xxxy_xxxzzz[k] = -g_y_0_xxy_xxxzzz[k] * ab_x + g_y_0_xxy_xxxxzzz[k];

                g_y_0_xxxy_xxyyyy[k] = -g_y_0_xxy_xxyyyy[k] * ab_x + g_y_0_xxy_xxxyyyy[k];

                g_y_0_xxxy_xxyyyz[k] = -g_y_0_xxy_xxyyyz[k] * ab_x + g_y_0_xxy_xxxyyyz[k];

                g_y_0_xxxy_xxyyzz[k] = -g_y_0_xxy_xxyyzz[k] * ab_x + g_y_0_xxy_xxxyyzz[k];

                g_y_0_xxxy_xxyzzz[k] = -g_y_0_xxy_xxyzzz[k] * ab_x + g_y_0_xxy_xxxyzzz[k];

                g_y_0_xxxy_xxzzzz[k] = -g_y_0_xxy_xxzzzz[k] * ab_x + g_y_0_xxy_xxxzzzz[k];

                g_y_0_xxxy_xyyyyy[k] = -g_y_0_xxy_xyyyyy[k] * ab_x + g_y_0_xxy_xxyyyyy[k];

                g_y_0_xxxy_xyyyyz[k] = -g_y_0_xxy_xyyyyz[k] * ab_x + g_y_0_xxy_xxyyyyz[k];

                g_y_0_xxxy_xyyyzz[k] = -g_y_0_xxy_xyyyzz[k] * ab_x + g_y_0_xxy_xxyyyzz[k];

                g_y_0_xxxy_xyyzzz[k] = -g_y_0_xxy_xyyzzz[k] * ab_x + g_y_0_xxy_xxyyzzz[k];

                g_y_0_xxxy_xyzzzz[k] = -g_y_0_xxy_xyzzzz[k] * ab_x + g_y_0_xxy_xxyzzzz[k];

                g_y_0_xxxy_xzzzzz[k] = -g_y_0_xxy_xzzzzz[k] * ab_x + g_y_0_xxy_xxzzzzz[k];

                g_y_0_xxxy_yyyyyy[k] = -g_y_0_xxy_yyyyyy[k] * ab_x + g_y_0_xxy_xyyyyyy[k];

                g_y_0_xxxy_yyyyyz[k] = -g_y_0_xxy_yyyyyz[k] * ab_x + g_y_0_xxy_xyyyyyz[k];

                g_y_0_xxxy_yyyyzz[k] = -g_y_0_xxy_yyyyzz[k] * ab_x + g_y_0_xxy_xyyyyzz[k];

                g_y_0_xxxy_yyyzzz[k] = -g_y_0_xxy_yyyzzz[k] * ab_x + g_y_0_xxy_xyyyzzz[k];

                g_y_0_xxxy_yyzzzz[k] = -g_y_0_xxy_yyzzzz[k] * ab_x + g_y_0_xxy_xyyzzzz[k];

                g_y_0_xxxy_yzzzzz[k] = -g_y_0_xxy_yzzzzz[k] * ab_x + g_y_0_xxy_xyzzzzz[k];

                g_y_0_xxxy_zzzzzz[k] = -g_y_0_xxy_zzzzzz[k] * ab_x + g_y_0_xxy_xzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_yyyyyy, g_y_0_xxxz_yyyyyz, g_y_0_xxxz_yyyyzz, g_y_0_xxxz_yyyzzz, g_y_0_xxxz_yyzzzz, g_y_0_xxxz_yzzzzz, g_y_0_xxxz_zzzzzz, g_y_0_xxz_xxxxxx, g_y_0_xxz_xxxxxxx, g_y_0_xxz_xxxxxxy, g_y_0_xxz_xxxxxxz, g_y_0_xxz_xxxxxy, g_y_0_xxz_xxxxxyy, g_y_0_xxz_xxxxxyz, g_y_0_xxz_xxxxxz, g_y_0_xxz_xxxxxzz, g_y_0_xxz_xxxxyy, g_y_0_xxz_xxxxyyy, g_y_0_xxz_xxxxyyz, g_y_0_xxz_xxxxyz, g_y_0_xxz_xxxxyzz, g_y_0_xxz_xxxxzz, g_y_0_xxz_xxxxzzz, g_y_0_xxz_xxxyyy, g_y_0_xxz_xxxyyyy, g_y_0_xxz_xxxyyyz, g_y_0_xxz_xxxyyz, g_y_0_xxz_xxxyyzz, g_y_0_xxz_xxxyzz, g_y_0_xxz_xxxyzzz, g_y_0_xxz_xxxzzz, g_y_0_xxz_xxxzzzz, g_y_0_xxz_xxyyyy, g_y_0_xxz_xxyyyyy, g_y_0_xxz_xxyyyyz, g_y_0_xxz_xxyyyz, g_y_0_xxz_xxyyyzz, g_y_0_xxz_xxyyzz, g_y_0_xxz_xxyyzzz, g_y_0_xxz_xxyzzz, g_y_0_xxz_xxyzzzz, g_y_0_xxz_xxzzzz, g_y_0_xxz_xxzzzzz, g_y_0_xxz_xyyyyy, g_y_0_xxz_xyyyyyy, g_y_0_xxz_xyyyyyz, g_y_0_xxz_xyyyyz, g_y_0_xxz_xyyyyzz, g_y_0_xxz_xyyyzz, g_y_0_xxz_xyyyzzz, g_y_0_xxz_xyyzzz, g_y_0_xxz_xyyzzzz, g_y_0_xxz_xyzzzz, g_y_0_xxz_xyzzzzz, g_y_0_xxz_xzzzzz, g_y_0_xxz_xzzzzzz, g_y_0_xxz_yyyyyy, g_y_0_xxz_yyyyyz, g_y_0_xxz_yyyyzz, g_y_0_xxz_yyyzzz, g_y_0_xxz_yyzzzz, g_y_0_xxz_yzzzzz, g_y_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_xxxxxx[k] = -g_y_0_xxz_xxxxxx[k] * ab_x + g_y_0_xxz_xxxxxxx[k];

                g_y_0_xxxz_xxxxxy[k] = -g_y_0_xxz_xxxxxy[k] * ab_x + g_y_0_xxz_xxxxxxy[k];

                g_y_0_xxxz_xxxxxz[k] = -g_y_0_xxz_xxxxxz[k] * ab_x + g_y_0_xxz_xxxxxxz[k];

                g_y_0_xxxz_xxxxyy[k] = -g_y_0_xxz_xxxxyy[k] * ab_x + g_y_0_xxz_xxxxxyy[k];

                g_y_0_xxxz_xxxxyz[k] = -g_y_0_xxz_xxxxyz[k] * ab_x + g_y_0_xxz_xxxxxyz[k];

                g_y_0_xxxz_xxxxzz[k] = -g_y_0_xxz_xxxxzz[k] * ab_x + g_y_0_xxz_xxxxxzz[k];

                g_y_0_xxxz_xxxyyy[k] = -g_y_0_xxz_xxxyyy[k] * ab_x + g_y_0_xxz_xxxxyyy[k];

                g_y_0_xxxz_xxxyyz[k] = -g_y_0_xxz_xxxyyz[k] * ab_x + g_y_0_xxz_xxxxyyz[k];

                g_y_0_xxxz_xxxyzz[k] = -g_y_0_xxz_xxxyzz[k] * ab_x + g_y_0_xxz_xxxxyzz[k];

                g_y_0_xxxz_xxxzzz[k] = -g_y_0_xxz_xxxzzz[k] * ab_x + g_y_0_xxz_xxxxzzz[k];

                g_y_0_xxxz_xxyyyy[k] = -g_y_0_xxz_xxyyyy[k] * ab_x + g_y_0_xxz_xxxyyyy[k];

                g_y_0_xxxz_xxyyyz[k] = -g_y_0_xxz_xxyyyz[k] * ab_x + g_y_0_xxz_xxxyyyz[k];

                g_y_0_xxxz_xxyyzz[k] = -g_y_0_xxz_xxyyzz[k] * ab_x + g_y_0_xxz_xxxyyzz[k];

                g_y_0_xxxz_xxyzzz[k] = -g_y_0_xxz_xxyzzz[k] * ab_x + g_y_0_xxz_xxxyzzz[k];

                g_y_0_xxxz_xxzzzz[k] = -g_y_0_xxz_xxzzzz[k] * ab_x + g_y_0_xxz_xxxzzzz[k];

                g_y_0_xxxz_xyyyyy[k] = -g_y_0_xxz_xyyyyy[k] * ab_x + g_y_0_xxz_xxyyyyy[k];

                g_y_0_xxxz_xyyyyz[k] = -g_y_0_xxz_xyyyyz[k] * ab_x + g_y_0_xxz_xxyyyyz[k];

                g_y_0_xxxz_xyyyzz[k] = -g_y_0_xxz_xyyyzz[k] * ab_x + g_y_0_xxz_xxyyyzz[k];

                g_y_0_xxxz_xyyzzz[k] = -g_y_0_xxz_xyyzzz[k] * ab_x + g_y_0_xxz_xxyyzzz[k];

                g_y_0_xxxz_xyzzzz[k] = -g_y_0_xxz_xyzzzz[k] * ab_x + g_y_0_xxz_xxyzzzz[k];

                g_y_0_xxxz_xzzzzz[k] = -g_y_0_xxz_xzzzzz[k] * ab_x + g_y_0_xxz_xxzzzzz[k];

                g_y_0_xxxz_yyyyyy[k] = -g_y_0_xxz_yyyyyy[k] * ab_x + g_y_0_xxz_xyyyyyy[k];

                g_y_0_xxxz_yyyyyz[k] = -g_y_0_xxz_yyyyyz[k] * ab_x + g_y_0_xxz_xyyyyyz[k];

                g_y_0_xxxz_yyyyzz[k] = -g_y_0_xxz_yyyyzz[k] * ab_x + g_y_0_xxz_xyyyyzz[k];

                g_y_0_xxxz_yyyzzz[k] = -g_y_0_xxz_yyyzzz[k] * ab_x + g_y_0_xxz_xyyyzzz[k];

                g_y_0_xxxz_yyzzzz[k] = -g_y_0_xxz_yyzzzz[k] * ab_x + g_y_0_xxz_xyyzzzz[k];

                g_y_0_xxxz_yzzzzz[k] = -g_y_0_xxz_yzzzzz[k] * ab_x + g_y_0_xxz_xyzzzzz[k];

                g_y_0_xxxz_zzzzzz[k] = -g_y_0_xxz_zzzzzz[k] * ab_x + g_y_0_xxz_xzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_yyyyyy, g_y_0_xxyy_yyyyyz, g_y_0_xxyy_yyyyzz, g_y_0_xxyy_yyyzzz, g_y_0_xxyy_yyzzzz, g_y_0_xxyy_yzzzzz, g_y_0_xxyy_zzzzzz, g_y_0_xyy_xxxxxx, g_y_0_xyy_xxxxxxx, g_y_0_xyy_xxxxxxy, g_y_0_xyy_xxxxxxz, g_y_0_xyy_xxxxxy, g_y_0_xyy_xxxxxyy, g_y_0_xyy_xxxxxyz, g_y_0_xyy_xxxxxz, g_y_0_xyy_xxxxxzz, g_y_0_xyy_xxxxyy, g_y_0_xyy_xxxxyyy, g_y_0_xyy_xxxxyyz, g_y_0_xyy_xxxxyz, g_y_0_xyy_xxxxyzz, g_y_0_xyy_xxxxzz, g_y_0_xyy_xxxxzzz, g_y_0_xyy_xxxyyy, g_y_0_xyy_xxxyyyy, g_y_0_xyy_xxxyyyz, g_y_0_xyy_xxxyyz, g_y_0_xyy_xxxyyzz, g_y_0_xyy_xxxyzz, g_y_0_xyy_xxxyzzz, g_y_0_xyy_xxxzzz, g_y_0_xyy_xxxzzzz, g_y_0_xyy_xxyyyy, g_y_0_xyy_xxyyyyy, g_y_0_xyy_xxyyyyz, g_y_0_xyy_xxyyyz, g_y_0_xyy_xxyyyzz, g_y_0_xyy_xxyyzz, g_y_0_xyy_xxyyzzz, g_y_0_xyy_xxyzzz, g_y_0_xyy_xxyzzzz, g_y_0_xyy_xxzzzz, g_y_0_xyy_xxzzzzz, g_y_0_xyy_xyyyyy, g_y_0_xyy_xyyyyyy, g_y_0_xyy_xyyyyyz, g_y_0_xyy_xyyyyz, g_y_0_xyy_xyyyyzz, g_y_0_xyy_xyyyzz, g_y_0_xyy_xyyyzzz, g_y_0_xyy_xyyzzz, g_y_0_xyy_xyyzzzz, g_y_0_xyy_xyzzzz, g_y_0_xyy_xyzzzzz, g_y_0_xyy_xzzzzz, g_y_0_xyy_xzzzzzz, g_y_0_xyy_yyyyyy, g_y_0_xyy_yyyyyz, g_y_0_xyy_yyyyzz, g_y_0_xyy_yyyzzz, g_y_0_xyy_yyzzzz, g_y_0_xyy_yzzzzz, g_y_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_xxxxxx[k] = -g_y_0_xyy_xxxxxx[k] * ab_x + g_y_0_xyy_xxxxxxx[k];

                g_y_0_xxyy_xxxxxy[k] = -g_y_0_xyy_xxxxxy[k] * ab_x + g_y_0_xyy_xxxxxxy[k];

                g_y_0_xxyy_xxxxxz[k] = -g_y_0_xyy_xxxxxz[k] * ab_x + g_y_0_xyy_xxxxxxz[k];

                g_y_0_xxyy_xxxxyy[k] = -g_y_0_xyy_xxxxyy[k] * ab_x + g_y_0_xyy_xxxxxyy[k];

                g_y_0_xxyy_xxxxyz[k] = -g_y_0_xyy_xxxxyz[k] * ab_x + g_y_0_xyy_xxxxxyz[k];

                g_y_0_xxyy_xxxxzz[k] = -g_y_0_xyy_xxxxzz[k] * ab_x + g_y_0_xyy_xxxxxzz[k];

                g_y_0_xxyy_xxxyyy[k] = -g_y_0_xyy_xxxyyy[k] * ab_x + g_y_0_xyy_xxxxyyy[k];

                g_y_0_xxyy_xxxyyz[k] = -g_y_0_xyy_xxxyyz[k] * ab_x + g_y_0_xyy_xxxxyyz[k];

                g_y_0_xxyy_xxxyzz[k] = -g_y_0_xyy_xxxyzz[k] * ab_x + g_y_0_xyy_xxxxyzz[k];

                g_y_0_xxyy_xxxzzz[k] = -g_y_0_xyy_xxxzzz[k] * ab_x + g_y_0_xyy_xxxxzzz[k];

                g_y_0_xxyy_xxyyyy[k] = -g_y_0_xyy_xxyyyy[k] * ab_x + g_y_0_xyy_xxxyyyy[k];

                g_y_0_xxyy_xxyyyz[k] = -g_y_0_xyy_xxyyyz[k] * ab_x + g_y_0_xyy_xxxyyyz[k];

                g_y_0_xxyy_xxyyzz[k] = -g_y_0_xyy_xxyyzz[k] * ab_x + g_y_0_xyy_xxxyyzz[k];

                g_y_0_xxyy_xxyzzz[k] = -g_y_0_xyy_xxyzzz[k] * ab_x + g_y_0_xyy_xxxyzzz[k];

                g_y_0_xxyy_xxzzzz[k] = -g_y_0_xyy_xxzzzz[k] * ab_x + g_y_0_xyy_xxxzzzz[k];

                g_y_0_xxyy_xyyyyy[k] = -g_y_0_xyy_xyyyyy[k] * ab_x + g_y_0_xyy_xxyyyyy[k];

                g_y_0_xxyy_xyyyyz[k] = -g_y_0_xyy_xyyyyz[k] * ab_x + g_y_0_xyy_xxyyyyz[k];

                g_y_0_xxyy_xyyyzz[k] = -g_y_0_xyy_xyyyzz[k] * ab_x + g_y_0_xyy_xxyyyzz[k];

                g_y_0_xxyy_xyyzzz[k] = -g_y_0_xyy_xyyzzz[k] * ab_x + g_y_0_xyy_xxyyzzz[k];

                g_y_0_xxyy_xyzzzz[k] = -g_y_0_xyy_xyzzzz[k] * ab_x + g_y_0_xyy_xxyzzzz[k];

                g_y_0_xxyy_xzzzzz[k] = -g_y_0_xyy_xzzzzz[k] * ab_x + g_y_0_xyy_xxzzzzz[k];

                g_y_0_xxyy_yyyyyy[k] = -g_y_0_xyy_yyyyyy[k] * ab_x + g_y_0_xyy_xyyyyyy[k];

                g_y_0_xxyy_yyyyyz[k] = -g_y_0_xyy_yyyyyz[k] * ab_x + g_y_0_xyy_xyyyyyz[k];

                g_y_0_xxyy_yyyyzz[k] = -g_y_0_xyy_yyyyzz[k] * ab_x + g_y_0_xyy_xyyyyzz[k];

                g_y_0_xxyy_yyyzzz[k] = -g_y_0_xyy_yyyzzz[k] * ab_x + g_y_0_xyy_xyyyzzz[k];

                g_y_0_xxyy_yyzzzz[k] = -g_y_0_xyy_yyzzzz[k] * ab_x + g_y_0_xyy_xyyzzzz[k];

                g_y_0_xxyy_yzzzzz[k] = -g_y_0_xyy_yzzzzz[k] * ab_x + g_y_0_xyy_xyzzzzz[k];

                g_y_0_xxyy_zzzzzz[k] = -g_y_0_xyy_zzzzzz[k] * ab_x + g_y_0_xyy_xzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_yyyyyy, g_y_0_xxyz_yyyyyz, g_y_0_xxyz_yyyyzz, g_y_0_xxyz_yyyzzz, g_y_0_xxyz_yyzzzz, g_y_0_xxyz_yzzzzz, g_y_0_xxyz_zzzzzz, g_y_0_xyz_xxxxxx, g_y_0_xyz_xxxxxxx, g_y_0_xyz_xxxxxxy, g_y_0_xyz_xxxxxxz, g_y_0_xyz_xxxxxy, g_y_0_xyz_xxxxxyy, g_y_0_xyz_xxxxxyz, g_y_0_xyz_xxxxxz, g_y_0_xyz_xxxxxzz, g_y_0_xyz_xxxxyy, g_y_0_xyz_xxxxyyy, g_y_0_xyz_xxxxyyz, g_y_0_xyz_xxxxyz, g_y_0_xyz_xxxxyzz, g_y_0_xyz_xxxxzz, g_y_0_xyz_xxxxzzz, g_y_0_xyz_xxxyyy, g_y_0_xyz_xxxyyyy, g_y_0_xyz_xxxyyyz, g_y_0_xyz_xxxyyz, g_y_0_xyz_xxxyyzz, g_y_0_xyz_xxxyzz, g_y_0_xyz_xxxyzzz, g_y_0_xyz_xxxzzz, g_y_0_xyz_xxxzzzz, g_y_0_xyz_xxyyyy, g_y_0_xyz_xxyyyyy, g_y_0_xyz_xxyyyyz, g_y_0_xyz_xxyyyz, g_y_0_xyz_xxyyyzz, g_y_0_xyz_xxyyzz, g_y_0_xyz_xxyyzzz, g_y_0_xyz_xxyzzz, g_y_0_xyz_xxyzzzz, g_y_0_xyz_xxzzzz, g_y_0_xyz_xxzzzzz, g_y_0_xyz_xyyyyy, g_y_0_xyz_xyyyyyy, g_y_0_xyz_xyyyyyz, g_y_0_xyz_xyyyyz, g_y_0_xyz_xyyyyzz, g_y_0_xyz_xyyyzz, g_y_0_xyz_xyyyzzz, g_y_0_xyz_xyyzzz, g_y_0_xyz_xyyzzzz, g_y_0_xyz_xyzzzz, g_y_0_xyz_xyzzzzz, g_y_0_xyz_xzzzzz, g_y_0_xyz_xzzzzzz, g_y_0_xyz_yyyyyy, g_y_0_xyz_yyyyyz, g_y_0_xyz_yyyyzz, g_y_0_xyz_yyyzzz, g_y_0_xyz_yyzzzz, g_y_0_xyz_yzzzzz, g_y_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_xxxxxx[k] = -g_y_0_xyz_xxxxxx[k] * ab_x + g_y_0_xyz_xxxxxxx[k];

                g_y_0_xxyz_xxxxxy[k] = -g_y_0_xyz_xxxxxy[k] * ab_x + g_y_0_xyz_xxxxxxy[k];

                g_y_0_xxyz_xxxxxz[k] = -g_y_0_xyz_xxxxxz[k] * ab_x + g_y_0_xyz_xxxxxxz[k];

                g_y_0_xxyz_xxxxyy[k] = -g_y_0_xyz_xxxxyy[k] * ab_x + g_y_0_xyz_xxxxxyy[k];

                g_y_0_xxyz_xxxxyz[k] = -g_y_0_xyz_xxxxyz[k] * ab_x + g_y_0_xyz_xxxxxyz[k];

                g_y_0_xxyz_xxxxzz[k] = -g_y_0_xyz_xxxxzz[k] * ab_x + g_y_0_xyz_xxxxxzz[k];

                g_y_0_xxyz_xxxyyy[k] = -g_y_0_xyz_xxxyyy[k] * ab_x + g_y_0_xyz_xxxxyyy[k];

                g_y_0_xxyz_xxxyyz[k] = -g_y_0_xyz_xxxyyz[k] * ab_x + g_y_0_xyz_xxxxyyz[k];

                g_y_0_xxyz_xxxyzz[k] = -g_y_0_xyz_xxxyzz[k] * ab_x + g_y_0_xyz_xxxxyzz[k];

                g_y_0_xxyz_xxxzzz[k] = -g_y_0_xyz_xxxzzz[k] * ab_x + g_y_0_xyz_xxxxzzz[k];

                g_y_0_xxyz_xxyyyy[k] = -g_y_0_xyz_xxyyyy[k] * ab_x + g_y_0_xyz_xxxyyyy[k];

                g_y_0_xxyz_xxyyyz[k] = -g_y_0_xyz_xxyyyz[k] * ab_x + g_y_0_xyz_xxxyyyz[k];

                g_y_0_xxyz_xxyyzz[k] = -g_y_0_xyz_xxyyzz[k] * ab_x + g_y_0_xyz_xxxyyzz[k];

                g_y_0_xxyz_xxyzzz[k] = -g_y_0_xyz_xxyzzz[k] * ab_x + g_y_0_xyz_xxxyzzz[k];

                g_y_0_xxyz_xxzzzz[k] = -g_y_0_xyz_xxzzzz[k] * ab_x + g_y_0_xyz_xxxzzzz[k];

                g_y_0_xxyz_xyyyyy[k] = -g_y_0_xyz_xyyyyy[k] * ab_x + g_y_0_xyz_xxyyyyy[k];

                g_y_0_xxyz_xyyyyz[k] = -g_y_0_xyz_xyyyyz[k] * ab_x + g_y_0_xyz_xxyyyyz[k];

                g_y_0_xxyz_xyyyzz[k] = -g_y_0_xyz_xyyyzz[k] * ab_x + g_y_0_xyz_xxyyyzz[k];

                g_y_0_xxyz_xyyzzz[k] = -g_y_0_xyz_xyyzzz[k] * ab_x + g_y_0_xyz_xxyyzzz[k];

                g_y_0_xxyz_xyzzzz[k] = -g_y_0_xyz_xyzzzz[k] * ab_x + g_y_0_xyz_xxyzzzz[k];

                g_y_0_xxyz_xzzzzz[k] = -g_y_0_xyz_xzzzzz[k] * ab_x + g_y_0_xyz_xxzzzzz[k];

                g_y_0_xxyz_yyyyyy[k] = -g_y_0_xyz_yyyyyy[k] * ab_x + g_y_0_xyz_xyyyyyy[k];

                g_y_0_xxyz_yyyyyz[k] = -g_y_0_xyz_yyyyyz[k] * ab_x + g_y_0_xyz_xyyyyyz[k];

                g_y_0_xxyz_yyyyzz[k] = -g_y_0_xyz_yyyyzz[k] * ab_x + g_y_0_xyz_xyyyyzz[k];

                g_y_0_xxyz_yyyzzz[k] = -g_y_0_xyz_yyyzzz[k] * ab_x + g_y_0_xyz_xyyyzzz[k];

                g_y_0_xxyz_yyzzzz[k] = -g_y_0_xyz_yyzzzz[k] * ab_x + g_y_0_xyz_xyyzzzz[k];

                g_y_0_xxyz_yzzzzz[k] = -g_y_0_xyz_yzzzzz[k] * ab_x + g_y_0_xyz_xyzzzzz[k];

                g_y_0_xxyz_zzzzzz[k] = -g_y_0_xyz_zzzzzz[k] * ab_x + g_y_0_xyz_xzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_yyyyyy, g_y_0_xxzz_yyyyyz, g_y_0_xxzz_yyyyzz, g_y_0_xxzz_yyyzzz, g_y_0_xxzz_yyzzzz, g_y_0_xxzz_yzzzzz, g_y_0_xxzz_zzzzzz, g_y_0_xzz_xxxxxx, g_y_0_xzz_xxxxxxx, g_y_0_xzz_xxxxxxy, g_y_0_xzz_xxxxxxz, g_y_0_xzz_xxxxxy, g_y_0_xzz_xxxxxyy, g_y_0_xzz_xxxxxyz, g_y_0_xzz_xxxxxz, g_y_0_xzz_xxxxxzz, g_y_0_xzz_xxxxyy, g_y_0_xzz_xxxxyyy, g_y_0_xzz_xxxxyyz, g_y_0_xzz_xxxxyz, g_y_0_xzz_xxxxyzz, g_y_0_xzz_xxxxzz, g_y_0_xzz_xxxxzzz, g_y_0_xzz_xxxyyy, g_y_0_xzz_xxxyyyy, g_y_0_xzz_xxxyyyz, g_y_0_xzz_xxxyyz, g_y_0_xzz_xxxyyzz, g_y_0_xzz_xxxyzz, g_y_0_xzz_xxxyzzz, g_y_0_xzz_xxxzzz, g_y_0_xzz_xxxzzzz, g_y_0_xzz_xxyyyy, g_y_0_xzz_xxyyyyy, g_y_0_xzz_xxyyyyz, g_y_0_xzz_xxyyyz, g_y_0_xzz_xxyyyzz, g_y_0_xzz_xxyyzz, g_y_0_xzz_xxyyzzz, g_y_0_xzz_xxyzzz, g_y_0_xzz_xxyzzzz, g_y_0_xzz_xxzzzz, g_y_0_xzz_xxzzzzz, g_y_0_xzz_xyyyyy, g_y_0_xzz_xyyyyyy, g_y_0_xzz_xyyyyyz, g_y_0_xzz_xyyyyz, g_y_0_xzz_xyyyyzz, g_y_0_xzz_xyyyzz, g_y_0_xzz_xyyyzzz, g_y_0_xzz_xyyzzz, g_y_0_xzz_xyyzzzz, g_y_0_xzz_xyzzzz, g_y_0_xzz_xyzzzzz, g_y_0_xzz_xzzzzz, g_y_0_xzz_xzzzzzz, g_y_0_xzz_yyyyyy, g_y_0_xzz_yyyyyz, g_y_0_xzz_yyyyzz, g_y_0_xzz_yyyzzz, g_y_0_xzz_yyzzzz, g_y_0_xzz_yzzzzz, g_y_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_xxxxxx[k] = -g_y_0_xzz_xxxxxx[k] * ab_x + g_y_0_xzz_xxxxxxx[k];

                g_y_0_xxzz_xxxxxy[k] = -g_y_0_xzz_xxxxxy[k] * ab_x + g_y_0_xzz_xxxxxxy[k];

                g_y_0_xxzz_xxxxxz[k] = -g_y_0_xzz_xxxxxz[k] * ab_x + g_y_0_xzz_xxxxxxz[k];

                g_y_0_xxzz_xxxxyy[k] = -g_y_0_xzz_xxxxyy[k] * ab_x + g_y_0_xzz_xxxxxyy[k];

                g_y_0_xxzz_xxxxyz[k] = -g_y_0_xzz_xxxxyz[k] * ab_x + g_y_0_xzz_xxxxxyz[k];

                g_y_0_xxzz_xxxxzz[k] = -g_y_0_xzz_xxxxzz[k] * ab_x + g_y_0_xzz_xxxxxzz[k];

                g_y_0_xxzz_xxxyyy[k] = -g_y_0_xzz_xxxyyy[k] * ab_x + g_y_0_xzz_xxxxyyy[k];

                g_y_0_xxzz_xxxyyz[k] = -g_y_0_xzz_xxxyyz[k] * ab_x + g_y_0_xzz_xxxxyyz[k];

                g_y_0_xxzz_xxxyzz[k] = -g_y_0_xzz_xxxyzz[k] * ab_x + g_y_0_xzz_xxxxyzz[k];

                g_y_0_xxzz_xxxzzz[k] = -g_y_0_xzz_xxxzzz[k] * ab_x + g_y_0_xzz_xxxxzzz[k];

                g_y_0_xxzz_xxyyyy[k] = -g_y_0_xzz_xxyyyy[k] * ab_x + g_y_0_xzz_xxxyyyy[k];

                g_y_0_xxzz_xxyyyz[k] = -g_y_0_xzz_xxyyyz[k] * ab_x + g_y_0_xzz_xxxyyyz[k];

                g_y_0_xxzz_xxyyzz[k] = -g_y_0_xzz_xxyyzz[k] * ab_x + g_y_0_xzz_xxxyyzz[k];

                g_y_0_xxzz_xxyzzz[k] = -g_y_0_xzz_xxyzzz[k] * ab_x + g_y_0_xzz_xxxyzzz[k];

                g_y_0_xxzz_xxzzzz[k] = -g_y_0_xzz_xxzzzz[k] * ab_x + g_y_0_xzz_xxxzzzz[k];

                g_y_0_xxzz_xyyyyy[k] = -g_y_0_xzz_xyyyyy[k] * ab_x + g_y_0_xzz_xxyyyyy[k];

                g_y_0_xxzz_xyyyyz[k] = -g_y_0_xzz_xyyyyz[k] * ab_x + g_y_0_xzz_xxyyyyz[k];

                g_y_0_xxzz_xyyyzz[k] = -g_y_0_xzz_xyyyzz[k] * ab_x + g_y_0_xzz_xxyyyzz[k];

                g_y_0_xxzz_xyyzzz[k] = -g_y_0_xzz_xyyzzz[k] * ab_x + g_y_0_xzz_xxyyzzz[k];

                g_y_0_xxzz_xyzzzz[k] = -g_y_0_xzz_xyzzzz[k] * ab_x + g_y_0_xzz_xxyzzzz[k];

                g_y_0_xxzz_xzzzzz[k] = -g_y_0_xzz_xzzzzz[k] * ab_x + g_y_0_xzz_xxzzzzz[k];

                g_y_0_xxzz_yyyyyy[k] = -g_y_0_xzz_yyyyyy[k] * ab_x + g_y_0_xzz_xyyyyyy[k];

                g_y_0_xxzz_yyyyyz[k] = -g_y_0_xzz_yyyyyz[k] * ab_x + g_y_0_xzz_xyyyyyz[k];

                g_y_0_xxzz_yyyyzz[k] = -g_y_0_xzz_yyyyzz[k] * ab_x + g_y_0_xzz_xyyyyzz[k];

                g_y_0_xxzz_yyyzzz[k] = -g_y_0_xzz_yyyzzz[k] * ab_x + g_y_0_xzz_xyyyzzz[k];

                g_y_0_xxzz_yyzzzz[k] = -g_y_0_xzz_yyzzzz[k] * ab_x + g_y_0_xzz_xyyzzzz[k];

                g_y_0_xxzz_yzzzzz[k] = -g_y_0_xzz_yzzzzz[k] * ab_x + g_y_0_xzz_xyzzzzz[k];

                g_y_0_xxzz_zzzzzz[k] = -g_y_0_xzz_zzzzzz[k] * ab_x + g_y_0_xzz_xzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_yyyyyy, g_y_0_xyyy_yyyyyz, g_y_0_xyyy_yyyyzz, g_y_0_xyyy_yyyzzz, g_y_0_xyyy_yyzzzz, g_y_0_xyyy_yzzzzz, g_y_0_xyyy_zzzzzz, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxx, g_y_0_yyy_xxxxxxy, g_y_0_yyy_xxxxxxz, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxxzz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxxzzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxxzzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xxzzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_xzzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_xxxxxx[k] = -g_y_0_yyy_xxxxxx[k] * ab_x + g_y_0_yyy_xxxxxxx[k];

                g_y_0_xyyy_xxxxxy[k] = -g_y_0_yyy_xxxxxy[k] * ab_x + g_y_0_yyy_xxxxxxy[k];

                g_y_0_xyyy_xxxxxz[k] = -g_y_0_yyy_xxxxxz[k] * ab_x + g_y_0_yyy_xxxxxxz[k];

                g_y_0_xyyy_xxxxyy[k] = -g_y_0_yyy_xxxxyy[k] * ab_x + g_y_0_yyy_xxxxxyy[k];

                g_y_0_xyyy_xxxxyz[k] = -g_y_0_yyy_xxxxyz[k] * ab_x + g_y_0_yyy_xxxxxyz[k];

                g_y_0_xyyy_xxxxzz[k] = -g_y_0_yyy_xxxxzz[k] * ab_x + g_y_0_yyy_xxxxxzz[k];

                g_y_0_xyyy_xxxyyy[k] = -g_y_0_yyy_xxxyyy[k] * ab_x + g_y_0_yyy_xxxxyyy[k];

                g_y_0_xyyy_xxxyyz[k] = -g_y_0_yyy_xxxyyz[k] * ab_x + g_y_0_yyy_xxxxyyz[k];

                g_y_0_xyyy_xxxyzz[k] = -g_y_0_yyy_xxxyzz[k] * ab_x + g_y_0_yyy_xxxxyzz[k];

                g_y_0_xyyy_xxxzzz[k] = -g_y_0_yyy_xxxzzz[k] * ab_x + g_y_0_yyy_xxxxzzz[k];

                g_y_0_xyyy_xxyyyy[k] = -g_y_0_yyy_xxyyyy[k] * ab_x + g_y_0_yyy_xxxyyyy[k];

                g_y_0_xyyy_xxyyyz[k] = -g_y_0_yyy_xxyyyz[k] * ab_x + g_y_0_yyy_xxxyyyz[k];

                g_y_0_xyyy_xxyyzz[k] = -g_y_0_yyy_xxyyzz[k] * ab_x + g_y_0_yyy_xxxyyzz[k];

                g_y_0_xyyy_xxyzzz[k] = -g_y_0_yyy_xxyzzz[k] * ab_x + g_y_0_yyy_xxxyzzz[k];

                g_y_0_xyyy_xxzzzz[k] = -g_y_0_yyy_xxzzzz[k] * ab_x + g_y_0_yyy_xxxzzzz[k];

                g_y_0_xyyy_xyyyyy[k] = -g_y_0_yyy_xyyyyy[k] * ab_x + g_y_0_yyy_xxyyyyy[k];

                g_y_0_xyyy_xyyyyz[k] = -g_y_0_yyy_xyyyyz[k] * ab_x + g_y_0_yyy_xxyyyyz[k];

                g_y_0_xyyy_xyyyzz[k] = -g_y_0_yyy_xyyyzz[k] * ab_x + g_y_0_yyy_xxyyyzz[k];

                g_y_0_xyyy_xyyzzz[k] = -g_y_0_yyy_xyyzzz[k] * ab_x + g_y_0_yyy_xxyyzzz[k];

                g_y_0_xyyy_xyzzzz[k] = -g_y_0_yyy_xyzzzz[k] * ab_x + g_y_0_yyy_xxyzzzz[k];

                g_y_0_xyyy_xzzzzz[k] = -g_y_0_yyy_xzzzzz[k] * ab_x + g_y_0_yyy_xxzzzzz[k];

                g_y_0_xyyy_yyyyyy[k] = -g_y_0_yyy_yyyyyy[k] * ab_x + g_y_0_yyy_xyyyyyy[k];

                g_y_0_xyyy_yyyyyz[k] = -g_y_0_yyy_yyyyyz[k] * ab_x + g_y_0_yyy_xyyyyyz[k];

                g_y_0_xyyy_yyyyzz[k] = -g_y_0_yyy_yyyyzz[k] * ab_x + g_y_0_yyy_xyyyyzz[k];

                g_y_0_xyyy_yyyzzz[k] = -g_y_0_yyy_yyyzzz[k] * ab_x + g_y_0_yyy_xyyyzzz[k];

                g_y_0_xyyy_yyzzzz[k] = -g_y_0_yyy_yyzzzz[k] * ab_x + g_y_0_yyy_xyyzzzz[k];

                g_y_0_xyyy_yzzzzz[k] = -g_y_0_yyy_yzzzzz[k] * ab_x + g_y_0_yyy_xyzzzzz[k];

                g_y_0_xyyy_zzzzzz[k] = -g_y_0_yyy_zzzzzz[k] * ab_x + g_y_0_yyy_xzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_yyyyyy, g_y_0_xyyz_yyyyyz, g_y_0_xyyz_yyyyzz, g_y_0_xyyz_yyyzzz, g_y_0_xyyz_yyzzzz, g_y_0_xyyz_yzzzzz, g_y_0_xyyz_zzzzzz, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxxx, g_y_0_yyz_xxxxxxy, g_y_0_yyz_xxxxxxz, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxyy, g_y_0_yyz_xxxxxyz, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxxzz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyyy, g_y_0_yyz_xxxxyyz, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxyzz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxxzzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyyy, g_y_0_yyz_xxxyyyz, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyyzz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxyzzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxxzzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyyy, g_y_0_yyz_xxyyyyz, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyyzz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyyzzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxyzzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xxzzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyyy, g_y_0_yyz_xyyyyyz, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyyzz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyyzzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyyzzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xyzzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_xzzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_xxxxxx[k] = -g_y_0_yyz_xxxxxx[k] * ab_x + g_y_0_yyz_xxxxxxx[k];

                g_y_0_xyyz_xxxxxy[k] = -g_y_0_yyz_xxxxxy[k] * ab_x + g_y_0_yyz_xxxxxxy[k];

                g_y_0_xyyz_xxxxxz[k] = -g_y_0_yyz_xxxxxz[k] * ab_x + g_y_0_yyz_xxxxxxz[k];

                g_y_0_xyyz_xxxxyy[k] = -g_y_0_yyz_xxxxyy[k] * ab_x + g_y_0_yyz_xxxxxyy[k];

                g_y_0_xyyz_xxxxyz[k] = -g_y_0_yyz_xxxxyz[k] * ab_x + g_y_0_yyz_xxxxxyz[k];

                g_y_0_xyyz_xxxxzz[k] = -g_y_0_yyz_xxxxzz[k] * ab_x + g_y_0_yyz_xxxxxzz[k];

                g_y_0_xyyz_xxxyyy[k] = -g_y_0_yyz_xxxyyy[k] * ab_x + g_y_0_yyz_xxxxyyy[k];

                g_y_0_xyyz_xxxyyz[k] = -g_y_0_yyz_xxxyyz[k] * ab_x + g_y_0_yyz_xxxxyyz[k];

                g_y_0_xyyz_xxxyzz[k] = -g_y_0_yyz_xxxyzz[k] * ab_x + g_y_0_yyz_xxxxyzz[k];

                g_y_0_xyyz_xxxzzz[k] = -g_y_0_yyz_xxxzzz[k] * ab_x + g_y_0_yyz_xxxxzzz[k];

                g_y_0_xyyz_xxyyyy[k] = -g_y_0_yyz_xxyyyy[k] * ab_x + g_y_0_yyz_xxxyyyy[k];

                g_y_0_xyyz_xxyyyz[k] = -g_y_0_yyz_xxyyyz[k] * ab_x + g_y_0_yyz_xxxyyyz[k];

                g_y_0_xyyz_xxyyzz[k] = -g_y_0_yyz_xxyyzz[k] * ab_x + g_y_0_yyz_xxxyyzz[k];

                g_y_0_xyyz_xxyzzz[k] = -g_y_0_yyz_xxyzzz[k] * ab_x + g_y_0_yyz_xxxyzzz[k];

                g_y_0_xyyz_xxzzzz[k] = -g_y_0_yyz_xxzzzz[k] * ab_x + g_y_0_yyz_xxxzzzz[k];

                g_y_0_xyyz_xyyyyy[k] = -g_y_0_yyz_xyyyyy[k] * ab_x + g_y_0_yyz_xxyyyyy[k];

                g_y_0_xyyz_xyyyyz[k] = -g_y_0_yyz_xyyyyz[k] * ab_x + g_y_0_yyz_xxyyyyz[k];

                g_y_0_xyyz_xyyyzz[k] = -g_y_0_yyz_xyyyzz[k] * ab_x + g_y_0_yyz_xxyyyzz[k];

                g_y_0_xyyz_xyyzzz[k] = -g_y_0_yyz_xyyzzz[k] * ab_x + g_y_0_yyz_xxyyzzz[k];

                g_y_0_xyyz_xyzzzz[k] = -g_y_0_yyz_xyzzzz[k] * ab_x + g_y_0_yyz_xxyzzzz[k];

                g_y_0_xyyz_xzzzzz[k] = -g_y_0_yyz_xzzzzz[k] * ab_x + g_y_0_yyz_xxzzzzz[k];

                g_y_0_xyyz_yyyyyy[k] = -g_y_0_yyz_yyyyyy[k] * ab_x + g_y_0_yyz_xyyyyyy[k];

                g_y_0_xyyz_yyyyyz[k] = -g_y_0_yyz_yyyyyz[k] * ab_x + g_y_0_yyz_xyyyyyz[k];

                g_y_0_xyyz_yyyyzz[k] = -g_y_0_yyz_yyyyzz[k] * ab_x + g_y_0_yyz_xyyyyzz[k];

                g_y_0_xyyz_yyyzzz[k] = -g_y_0_yyz_yyyzzz[k] * ab_x + g_y_0_yyz_xyyyzzz[k];

                g_y_0_xyyz_yyzzzz[k] = -g_y_0_yyz_yyzzzz[k] * ab_x + g_y_0_yyz_xyyzzzz[k];

                g_y_0_xyyz_yzzzzz[k] = -g_y_0_yyz_yzzzzz[k] * ab_x + g_y_0_yyz_xyzzzzz[k];

                g_y_0_xyyz_zzzzzz[k] = -g_y_0_yyz_zzzzzz[k] * ab_x + g_y_0_yyz_xzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_yyyyyy, g_y_0_xyzz_yyyyyz, g_y_0_xyzz_yyyyzz, g_y_0_xyzz_yyyzzz, g_y_0_xyzz_yyzzzz, g_y_0_xyzz_yzzzzz, g_y_0_xyzz_zzzzzz, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxxx, g_y_0_yzz_xxxxxxy, g_y_0_yzz_xxxxxxz, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxyy, g_y_0_yzz_xxxxxyz, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxxzz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyyy, g_y_0_yzz_xxxxyyz, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxyzz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxxzzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyyy, g_y_0_yzz_xxxyyyz, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyyzz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxyzzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxxzzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyyy, g_y_0_yzz_xxyyyyz, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyyzz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyyzzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxyzzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xxzzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyyy, g_y_0_yzz_xyyyyyz, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyyzz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyyzzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyyzzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xyzzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_xzzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_xxxxxx[k] = -g_y_0_yzz_xxxxxx[k] * ab_x + g_y_0_yzz_xxxxxxx[k];

                g_y_0_xyzz_xxxxxy[k] = -g_y_0_yzz_xxxxxy[k] * ab_x + g_y_0_yzz_xxxxxxy[k];

                g_y_0_xyzz_xxxxxz[k] = -g_y_0_yzz_xxxxxz[k] * ab_x + g_y_0_yzz_xxxxxxz[k];

                g_y_0_xyzz_xxxxyy[k] = -g_y_0_yzz_xxxxyy[k] * ab_x + g_y_0_yzz_xxxxxyy[k];

                g_y_0_xyzz_xxxxyz[k] = -g_y_0_yzz_xxxxyz[k] * ab_x + g_y_0_yzz_xxxxxyz[k];

                g_y_0_xyzz_xxxxzz[k] = -g_y_0_yzz_xxxxzz[k] * ab_x + g_y_0_yzz_xxxxxzz[k];

                g_y_0_xyzz_xxxyyy[k] = -g_y_0_yzz_xxxyyy[k] * ab_x + g_y_0_yzz_xxxxyyy[k];

                g_y_0_xyzz_xxxyyz[k] = -g_y_0_yzz_xxxyyz[k] * ab_x + g_y_0_yzz_xxxxyyz[k];

                g_y_0_xyzz_xxxyzz[k] = -g_y_0_yzz_xxxyzz[k] * ab_x + g_y_0_yzz_xxxxyzz[k];

                g_y_0_xyzz_xxxzzz[k] = -g_y_0_yzz_xxxzzz[k] * ab_x + g_y_0_yzz_xxxxzzz[k];

                g_y_0_xyzz_xxyyyy[k] = -g_y_0_yzz_xxyyyy[k] * ab_x + g_y_0_yzz_xxxyyyy[k];

                g_y_0_xyzz_xxyyyz[k] = -g_y_0_yzz_xxyyyz[k] * ab_x + g_y_0_yzz_xxxyyyz[k];

                g_y_0_xyzz_xxyyzz[k] = -g_y_0_yzz_xxyyzz[k] * ab_x + g_y_0_yzz_xxxyyzz[k];

                g_y_0_xyzz_xxyzzz[k] = -g_y_0_yzz_xxyzzz[k] * ab_x + g_y_0_yzz_xxxyzzz[k];

                g_y_0_xyzz_xxzzzz[k] = -g_y_0_yzz_xxzzzz[k] * ab_x + g_y_0_yzz_xxxzzzz[k];

                g_y_0_xyzz_xyyyyy[k] = -g_y_0_yzz_xyyyyy[k] * ab_x + g_y_0_yzz_xxyyyyy[k];

                g_y_0_xyzz_xyyyyz[k] = -g_y_0_yzz_xyyyyz[k] * ab_x + g_y_0_yzz_xxyyyyz[k];

                g_y_0_xyzz_xyyyzz[k] = -g_y_0_yzz_xyyyzz[k] * ab_x + g_y_0_yzz_xxyyyzz[k];

                g_y_0_xyzz_xyyzzz[k] = -g_y_0_yzz_xyyzzz[k] * ab_x + g_y_0_yzz_xxyyzzz[k];

                g_y_0_xyzz_xyzzzz[k] = -g_y_0_yzz_xyzzzz[k] * ab_x + g_y_0_yzz_xxyzzzz[k];

                g_y_0_xyzz_xzzzzz[k] = -g_y_0_yzz_xzzzzz[k] * ab_x + g_y_0_yzz_xxzzzzz[k];

                g_y_0_xyzz_yyyyyy[k] = -g_y_0_yzz_yyyyyy[k] * ab_x + g_y_0_yzz_xyyyyyy[k];

                g_y_0_xyzz_yyyyyz[k] = -g_y_0_yzz_yyyyyz[k] * ab_x + g_y_0_yzz_xyyyyyz[k];

                g_y_0_xyzz_yyyyzz[k] = -g_y_0_yzz_yyyyzz[k] * ab_x + g_y_0_yzz_xyyyyzz[k];

                g_y_0_xyzz_yyyzzz[k] = -g_y_0_yzz_yyyzzz[k] * ab_x + g_y_0_yzz_xyyyzzz[k];

                g_y_0_xyzz_yyzzzz[k] = -g_y_0_yzz_yyzzzz[k] * ab_x + g_y_0_yzz_xyyzzzz[k];

                g_y_0_xyzz_yzzzzz[k] = -g_y_0_yzz_yzzzzz[k] * ab_x + g_y_0_yzz_xyzzzzz[k];

                g_y_0_xyzz_zzzzzz[k] = -g_y_0_yzz_zzzzzz[k] * ab_x + g_y_0_yzz_xzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_yyyyyy, g_y_0_xzzz_yyyyyz, g_y_0_xzzz_yyyyzz, g_y_0_xzzz_yyyzzz, g_y_0_xzzz_yyzzzz, g_y_0_xzzz_yzzzzz, g_y_0_xzzz_zzzzzz, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxxx, g_y_0_zzz_xxxxxxy, g_y_0_zzz_xxxxxxz, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxyy, g_y_0_zzz_xxxxxyz, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxxzz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyyy, g_y_0_zzz_xxxxyyz, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxyzz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxxzzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyyy, g_y_0_zzz_xxxyyyz, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyyzz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxyzzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxxzzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyyy, g_y_0_zzz_xxyyyyz, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyyzz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyyzzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxyzzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xxzzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyyy, g_y_0_zzz_xyyyyyz, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyyzz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyyzzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyyzzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xyzzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_xzzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_xxxxxx[k] = -g_y_0_zzz_xxxxxx[k] * ab_x + g_y_0_zzz_xxxxxxx[k];

                g_y_0_xzzz_xxxxxy[k] = -g_y_0_zzz_xxxxxy[k] * ab_x + g_y_0_zzz_xxxxxxy[k];

                g_y_0_xzzz_xxxxxz[k] = -g_y_0_zzz_xxxxxz[k] * ab_x + g_y_0_zzz_xxxxxxz[k];

                g_y_0_xzzz_xxxxyy[k] = -g_y_0_zzz_xxxxyy[k] * ab_x + g_y_0_zzz_xxxxxyy[k];

                g_y_0_xzzz_xxxxyz[k] = -g_y_0_zzz_xxxxyz[k] * ab_x + g_y_0_zzz_xxxxxyz[k];

                g_y_0_xzzz_xxxxzz[k] = -g_y_0_zzz_xxxxzz[k] * ab_x + g_y_0_zzz_xxxxxzz[k];

                g_y_0_xzzz_xxxyyy[k] = -g_y_0_zzz_xxxyyy[k] * ab_x + g_y_0_zzz_xxxxyyy[k];

                g_y_0_xzzz_xxxyyz[k] = -g_y_0_zzz_xxxyyz[k] * ab_x + g_y_0_zzz_xxxxyyz[k];

                g_y_0_xzzz_xxxyzz[k] = -g_y_0_zzz_xxxyzz[k] * ab_x + g_y_0_zzz_xxxxyzz[k];

                g_y_0_xzzz_xxxzzz[k] = -g_y_0_zzz_xxxzzz[k] * ab_x + g_y_0_zzz_xxxxzzz[k];

                g_y_0_xzzz_xxyyyy[k] = -g_y_0_zzz_xxyyyy[k] * ab_x + g_y_0_zzz_xxxyyyy[k];

                g_y_0_xzzz_xxyyyz[k] = -g_y_0_zzz_xxyyyz[k] * ab_x + g_y_0_zzz_xxxyyyz[k];

                g_y_0_xzzz_xxyyzz[k] = -g_y_0_zzz_xxyyzz[k] * ab_x + g_y_0_zzz_xxxyyzz[k];

                g_y_0_xzzz_xxyzzz[k] = -g_y_0_zzz_xxyzzz[k] * ab_x + g_y_0_zzz_xxxyzzz[k];

                g_y_0_xzzz_xxzzzz[k] = -g_y_0_zzz_xxzzzz[k] * ab_x + g_y_0_zzz_xxxzzzz[k];

                g_y_0_xzzz_xyyyyy[k] = -g_y_0_zzz_xyyyyy[k] * ab_x + g_y_0_zzz_xxyyyyy[k];

                g_y_0_xzzz_xyyyyz[k] = -g_y_0_zzz_xyyyyz[k] * ab_x + g_y_0_zzz_xxyyyyz[k];

                g_y_0_xzzz_xyyyzz[k] = -g_y_0_zzz_xyyyzz[k] * ab_x + g_y_0_zzz_xxyyyzz[k];

                g_y_0_xzzz_xyyzzz[k] = -g_y_0_zzz_xyyzzz[k] * ab_x + g_y_0_zzz_xxyyzzz[k];

                g_y_0_xzzz_xyzzzz[k] = -g_y_0_zzz_xyzzzz[k] * ab_x + g_y_0_zzz_xxyzzzz[k];

                g_y_0_xzzz_xzzzzz[k] = -g_y_0_zzz_xzzzzz[k] * ab_x + g_y_0_zzz_xxzzzzz[k];

                g_y_0_xzzz_yyyyyy[k] = -g_y_0_zzz_yyyyyy[k] * ab_x + g_y_0_zzz_xyyyyyy[k];

                g_y_0_xzzz_yyyyyz[k] = -g_y_0_zzz_yyyyyz[k] * ab_x + g_y_0_zzz_xyyyyyz[k];

                g_y_0_xzzz_yyyyzz[k] = -g_y_0_zzz_yyyyzz[k] * ab_x + g_y_0_zzz_xyyyyzz[k];

                g_y_0_xzzz_yyyzzz[k] = -g_y_0_zzz_yyyzzz[k] * ab_x + g_y_0_zzz_xyyyzzz[k];

                g_y_0_xzzz_yyzzzz[k] = -g_y_0_zzz_yyzzzz[k] * ab_x + g_y_0_zzz_xyyzzzz[k];

                g_y_0_xzzz_yzzzzz[k] = -g_y_0_zzz_yzzzzz[k] * ab_x + g_y_0_zzz_xyzzzzz[k];

                g_y_0_xzzz_zzzzzz[k] = -g_y_0_zzz_zzzzzz[k] * ab_x + g_y_0_zzz_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxy, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyyy, g_y_0_yyy_yyyyyyz, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyyzz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyyzzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyyzzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yyzzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_yzzzzzz, g_y_0_yyy_zzzzzz, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzzz, g_yyy_xxxxxx, g_yyy_xxxxxy, g_yyy_xxxxxz, g_yyy_xxxxyy, g_yyy_xxxxyz, g_yyy_xxxxzz, g_yyy_xxxyyy, g_yyy_xxxyyz, g_yyy_xxxyzz, g_yyy_xxxzzz, g_yyy_xxyyyy, g_yyy_xxyyyz, g_yyy_xxyyzz, g_yyy_xxyzzz, g_yyy_xxzzzz, g_yyy_xyyyyy, g_yyy_xyyyyz, g_yyy_xyyyzz, g_yyy_xyyzzz, g_yyy_xyzzzz, g_yyy_xzzzzz, g_yyy_yyyyyy, g_yyy_yyyyyz, g_yyy_yyyyzz, g_yyy_yyyzzz, g_yyy_yyzzzz, g_yyy_yzzzzz, g_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_xxxxxx[k] = -g_yyy_xxxxxx[k] - g_y_0_yyy_xxxxxx[k] * ab_y + g_y_0_yyy_xxxxxxy[k];

                g_y_0_yyyy_xxxxxy[k] = -g_yyy_xxxxxy[k] - g_y_0_yyy_xxxxxy[k] * ab_y + g_y_0_yyy_xxxxxyy[k];

                g_y_0_yyyy_xxxxxz[k] = -g_yyy_xxxxxz[k] - g_y_0_yyy_xxxxxz[k] * ab_y + g_y_0_yyy_xxxxxyz[k];

                g_y_0_yyyy_xxxxyy[k] = -g_yyy_xxxxyy[k] - g_y_0_yyy_xxxxyy[k] * ab_y + g_y_0_yyy_xxxxyyy[k];

                g_y_0_yyyy_xxxxyz[k] = -g_yyy_xxxxyz[k] - g_y_0_yyy_xxxxyz[k] * ab_y + g_y_0_yyy_xxxxyyz[k];

                g_y_0_yyyy_xxxxzz[k] = -g_yyy_xxxxzz[k] - g_y_0_yyy_xxxxzz[k] * ab_y + g_y_0_yyy_xxxxyzz[k];

                g_y_0_yyyy_xxxyyy[k] = -g_yyy_xxxyyy[k] - g_y_0_yyy_xxxyyy[k] * ab_y + g_y_0_yyy_xxxyyyy[k];

                g_y_0_yyyy_xxxyyz[k] = -g_yyy_xxxyyz[k] - g_y_0_yyy_xxxyyz[k] * ab_y + g_y_0_yyy_xxxyyyz[k];

                g_y_0_yyyy_xxxyzz[k] = -g_yyy_xxxyzz[k] - g_y_0_yyy_xxxyzz[k] * ab_y + g_y_0_yyy_xxxyyzz[k];

                g_y_0_yyyy_xxxzzz[k] = -g_yyy_xxxzzz[k] - g_y_0_yyy_xxxzzz[k] * ab_y + g_y_0_yyy_xxxyzzz[k];

                g_y_0_yyyy_xxyyyy[k] = -g_yyy_xxyyyy[k] - g_y_0_yyy_xxyyyy[k] * ab_y + g_y_0_yyy_xxyyyyy[k];

                g_y_0_yyyy_xxyyyz[k] = -g_yyy_xxyyyz[k] - g_y_0_yyy_xxyyyz[k] * ab_y + g_y_0_yyy_xxyyyyz[k];

                g_y_0_yyyy_xxyyzz[k] = -g_yyy_xxyyzz[k] - g_y_0_yyy_xxyyzz[k] * ab_y + g_y_0_yyy_xxyyyzz[k];

                g_y_0_yyyy_xxyzzz[k] = -g_yyy_xxyzzz[k] - g_y_0_yyy_xxyzzz[k] * ab_y + g_y_0_yyy_xxyyzzz[k];

                g_y_0_yyyy_xxzzzz[k] = -g_yyy_xxzzzz[k] - g_y_0_yyy_xxzzzz[k] * ab_y + g_y_0_yyy_xxyzzzz[k];

                g_y_0_yyyy_xyyyyy[k] = -g_yyy_xyyyyy[k] - g_y_0_yyy_xyyyyy[k] * ab_y + g_y_0_yyy_xyyyyyy[k];

                g_y_0_yyyy_xyyyyz[k] = -g_yyy_xyyyyz[k] - g_y_0_yyy_xyyyyz[k] * ab_y + g_y_0_yyy_xyyyyyz[k];

                g_y_0_yyyy_xyyyzz[k] = -g_yyy_xyyyzz[k] - g_y_0_yyy_xyyyzz[k] * ab_y + g_y_0_yyy_xyyyyzz[k];

                g_y_0_yyyy_xyyzzz[k] = -g_yyy_xyyzzz[k] - g_y_0_yyy_xyyzzz[k] * ab_y + g_y_0_yyy_xyyyzzz[k];

                g_y_0_yyyy_xyzzzz[k] = -g_yyy_xyzzzz[k] - g_y_0_yyy_xyzzzz[k] * ab_y + g_y_0_yyy_xyyzzzz[k];

                g_y_0_yyyy_xzzzzz[k] = -g_yyy_xzzzzz[k] - g_y_0_yyy_xzzzzz[k] * ab_y + g_y_0_yyy_xyzzzzz[k];

                g_y_0_yyyy_yyyyyy[k] = -g_yyy_yyyyyy[k] - g_y_0_yyy_yyyyyy[k] * ab_y + g_y_0_yyy_yyyyyyy[k];

                g_y_0_yyyy_yyyyyz[k] = -g_yyy_yyyyyz[k] - g_y_0_yyy_yyyyyz[k] * ab_y + g_y_0_yyy_yyyyyyz[k];

                g_y_0_yyyy_yyyyzz[k] = -g_yyy_yyyyzz[k] - g_y_0_yyy_yyyyzz[k] * ab_y + g_y_0_yyy_yyyyyzz[k];

                g_y_0_yyyy_yyyzzz[k] = -g_yyy_yyyzzz[k] - g_y_0_yyy_yyyzzz[k] * ab_y + g_y_0_yyy_yyyyzzz[k];

                g_y_0_yyyy_yyzzzz[k] = -g_yyy_yyzzzz[k] - g_y_0_yyy_yyzzzz[k] * ab_y + g_y_0_yyy_yyyzzzz[k];

                g_y_0_yyyy_yzzzzz[k] = -g_yyy_yzzzzz[k] - g_y_0_yyy_yzzzzz[k] * ab_y + g_y_0_yyy_yyzzzzz[k];

                g_y_0_yyyy_zzzzzz[k] = -g_yyy_zzzzzz[k] - g_y_0_yyy_zzzzzz[k] * ab_y + g_y_0_yyy_yzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxz, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxxzz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxxzzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxxzzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xxzzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_xzzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyyz, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyyzz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyyzzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyyzzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yyzzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_yzzzzzz, g_y_0_yyy_zzzzzz, g_y_0_yyy_zzzzzzz, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_xxxxxx[k] = -g_y_0_yyy_xxxxxx[k] * ab_z + g_y_0_yyy_xxxxxxz[k];

                g_y_0_yyyz_xxxxxy[k] = -g_y_0_yyy_xxxxxy[k] * ab_z + g_y_0_yyy_xxxxxyz[k];

                g_y_0_yyyz_xxxxxz[k] = -g_y_0_yyy_xxxxxz[k] * ab_z + g_y_0_yyy_xxxxxzz[k];

                g_y_0_yyyz_xxxxyy[k] = -g_y_0_yyy_xxxxyy[k] * ab_z + g_y_0_yyy_xxxxyyz[k];

                g_y_0_yyyz_xxxxyz[k] = -g_y_0_yyy_xxxxyz[k] * ab_z + g_y_0_yyy_xxxxyzz[k];

                g_y_0_yyyz_xxxxzz[k] = -g_y_0_yyy_xxxxzz[k] * ab_z + g_y_0_yyy_xxxxzzz[k];

                g_y_0_yyyz_xxxyyy[k] = -g_y_0_yyy_xxxyyy[k] * ab_z + g_y_0_yyy_xxxyyyz[k];

                g_y_0_yyyz_xxxyyz[k] = -g_y_0_yyy_xxxyyz[k] * ab_z + g_y_0_yyy_xxxyyzz[k];

                g_y_0_yyyz_xxxyzz[k] = -g_y_0_yyy_xxxyzz[k] * ab_z + g_y_0_yyy_xxxyzzz[k];

                g_y_0_yyyz_xxxzzz[k] = -g_y_0_yyy_xxxzzz[k] * ab_z + g_y_0_yyy_xxxzzzz[k];

                g_y_0_yyyz_xxyyyy[k] = -g_y_0_yyy_xxyyyy[k] * ab_z + g_y_0_yyy_xxyyyyz[k];

                g_y_0_yyyz_xxyyyz[k] = -g_y_0_yyy_xxyyyz[k] * ab_z + g_y_0_yyy_xxyyyzz[k];

                g_y_0_yyyz_xxyyzz[k] = -g_y_0_yyy_xxyyzz[k] * ab_z + g_y_0_yyy_xxyyzzz[k];

                g_y_0_yyyz_xxyzzz[k] = -g_y_0_yyy_xxyzzz[k] * ab_z + g_y_0_yyy_xxyzzzz[k];

                g_y_0_yyyz_xxzzzz[k] = -g_y_0_yyy_xxzzzz[k] * ab_z + g_y_0_yyy_xxzzzzz[k];

                g_y_0_yyyz_xyyyyy[k] = -g_y_0_yyy_xyyyyy[k] * ab_z + g_y_0_yyy_xyyyyyz[k];

                g_y_0_yyyz_xyyyyz[k] = -g_y_0_yyy_xyyyyz[k] * ab_z + g_y_0_yyy_xyyyyzz[k];

                g_y_0_yyyz_xyyyzz[k] = -g_y_0_yyy_xyyyzz[k] * ab_z + g_y_0_yyy_xyyyzzz[k];

                g_y_0_yyyz_xyyzzz[k] = -g_y_0_yyy_xyyzzz[k] * ab_z + g_y_0_yyy_xyyzzzz[k];

                g_y_0_yyyz_xyzzzz[k] = -g_y_0_yyy_xyzzzz[k] * ab_z + g_y_0_yyy_xyzzzzz[k];

                g_y_0_yyyz_xzzzzz[k] = -g_y_0_yyy_xzzzzz[k] * ab_z + g_y_0_yyy_xzzzzzz[k];

                g_y_0_yyyz_yyyyyy[k] = -g_y_0_yyy_yyyyyy[k] * ab_z + g_y_0_yyy_yyyyyyz[k];

                g_y_0_yyyz_yyyyyz[k] = -g_y_0_yyy_yyyyyz[k] * ab_z + g_y_0_yyy_yyyyyzz[k];

                g_y_0_yyyz_yyyyzz[k] = -g_y_0_yyy_yyyyzz[k] * ab_z + g_y_0_yyy_yyyyzzz[k];

                g_y_0_yyyz_yyyzzz[k] = -g_y_0_yyy_yyyzzz[k] * ab_z + g_y_0_yyy_yyyzzzz[k];

                g_y_0_yyyz_yyzzzz[k] = -g_y_0_yyy_yyzzzz[k] * ab_z + g_y_0_yyy_yyzzzzz[k];

                g_y_0_yyyz_yzzzzz[k] = -g_y_0_yyy_yzzzzz[k] * ab_z + g_y_0_yyy_yzzzzzz[k];

                g_y_0_yyyz_zzzzzz[k] = -g_y_0_yyy_zzzzzz[k] * ab_z + g_y_0_yyy_zzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxxz, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxyz, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxxzz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyyz, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxyzz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxxzzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyyz, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyyzz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxyzzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxxzzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyyz, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyyzz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyyzzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxyzzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xxzzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyyz, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyyzz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyyzzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyyzzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xyzzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_xzzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyyz, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyyzz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyyzzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyyzzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yyzzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_yzzzzzz, g_y_0_yyz_zzzzzz, g_y_0_yyz_zzzzzzz, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_xxxxxx[k] = -g_y_0_yyz_xxxxxx[k] * ab_z + g_y_0_yyz_xxxxxxz[k];

                g_y_0_yyzz_xxxxxy[k] = -g_y_0_yyz_xxxxxy[k] * ab_z + g_y_0_yyz_xxxxxyz[k];

                g_y_0_yyzz_xxxxxz[k] = -g_y_0_yyz_xxxxxz[k] * ab_z + g_y_0_yyz_xxxxxzz[k];

                g_y_0_yyzz_xxxxyy[k] = -g_y_0_yyz_xxxxyy[k] * ab_z + g_y_0_yyz_xxxxyyz[k];

                g_y_0_yyzz_xxxxyz[k] = -g_y_0_yyz_xxxxyz[k] * ab_z + g_y_0_yyz_xxxxyzz[k];

                g_y_0_yyzz_xxxxzz[k] = -g_y_0_yyz_xxxxzz[k] * ab_z + g_y_0_yyz_xxxxzzz[k];

                g_y_0_yyzz_xxxyyy[k] = -g_y_0_yyz_xxxyyy[k] * ab_z + g_y_0_yyz_xxxyyyz[k];

                g_y_0_yyzz_xxxyyz[k] = -g_y_0_yyz_xxxyyz[k] * ab_z + g_y_0_yyz_xxxyyzz[k];

                g_y_0_yyzz_xxxyzz[k] = -g_y_0_yyz_xxxyzz[k] * ab_z + g_y_0_yyz_xxxyzzz[k];

                g_y_0_yyzz_xxxzzz[k] = -g_y_0_yyz_xxxzzz[k] * ab_z + g_y_0_yyz_xxxzzzz[k];

                g_y_0_yyzz_xxyyyy[k] = -g_y_0_yyz_xxyyyy[k] * ab_z + g_y_0_yyz_xxyyyyz[k];

                g_y_0_yyzz_xxyyyz[k] = -g_y_0_yyz_xxyyyz[k] * ab_z + g_y_0_yyz_xxyyyzz[k];

                g_y_0_yyzz_xxyyzz[k] = -g_y_0_yyz_xxyyzz[k] * ab_z + g_y_0_yyz_xxyyzzz[k];

                g_y_0_yyzz_xxyzzz[k] = -g_y_0_yyz_xxyzzz[k] * ab_z + g_y_0_yyz_xxyzzzz[k];

                g_y_0_yyzz_xxzzzz[k] = -g_y_0_yyz_xxzzzz[k] * ab_z + g_y_0_yyz_xxzzzzz[k];

                g_y_0_yyzz_xyyyyy[k] = -g_y_0_yyz_xyyyyy[k] * ab_z + g_y_0_yyz_xyyyyyz[k];

                g_y_0_yyzz_xyyyyz[k] = -g_y_0_yyz_xyyyyz[k] * ab_z + g_y_0_yyz_xyyyyzz[k];

                g_y_0_yyzz_xyyyzz[k] = -g_y_0_yyz_xyyyzz[k] * ab_z + g_y_0_yyz_xyyyzzz[k];

                g_y_0_yyzz_xyyzzz[k] = -g_y_0_yyz_xyyzzz[k] * ab_z + g_y_0_yyz_xyyzzzz[k];

                g_y_0_yyzz_xyzzzz[k] = -g_y_0_yyz_xyzzzz[k] * ab_z + g_y_0_yyz_xyzzzzz[k];

                g_y_0_yyzz_xzzzzz[k] = -g_y_0_yyz_xzzzzz[k] * ab_z + g_y_0_yyz_xzzzzzz[k];

                g_y_0_yyzz_yyyyyy[k] = -g_y_0_yyz_yyyyyy[k] * ab_z + g_y_0_yyz_yyyyyyz[k];

                g_y_0_yyzz_yyyyyz[k] = -g_y_0_yyz_yyyyyz[k] * ab_z + g_y_0_yyz_yyyyyzz[k];

                g_y_0_yyzz_yyyyzz[k] = -g_y_0_yyz_yyyyzz[k] * ab_z + g_y_0_yyz_yyyyzzz[k];

                g_y_0_yyzz_yyyzzz[k] = -g_y_0_yyz_yyyzzz[k] * ab_z + g_y_0_yyz_yyyzzzz[k];

                g_y_0_yyzz_yyzzzz[k] = -g_y_0_yyz_yyzzzz[k] * ab_z + g_y_0_yyz_yyzzzzz[k];

                g_y_0_yyzz_yzzzzz[k] = -g_y_0_yyz_yzzzzz[k] * ab_z + g_y_0_yyz_yzzzzzz[k];

                g_y_0_yyzz_zzzzzz[k] = -g_y_0_yyz_zzzzzz[k] * ab_z + g_y_0_yyz_zzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxxz, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxyz, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxxzz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyyz, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxyzz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxxzzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyyz, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyyzz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxyzzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxxzzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyyz, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyyzz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyyzzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxyzzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xxzzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyyz, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyyzz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyyzzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyyzzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xyzzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_xzzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyyz, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyyzz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyyzzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyyzzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yyzzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_yzzzzzz, g_y_0_yzz_zzzzzz, g_y_0_yzz_zzzzzzz, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_xxxxxx[k] = -g_y_0_yzz_xxxxxx[k] * ab_z + g_y_0_yzz_xxxxxxz[k];

                g_y_0_yzzz_xxxxxy[k] = -g_y_0_yzz_xxxxxy[k] * ab_z + g_y_0_yzz_xxxxxyz[k];

                g_y_0_yzzz_xxxxxz[k] = -g_y_0_yzz_xxxxxz[k] * ab_z + g_y_0_yzz_xxxxxzz[k];

                g_y_0_yzzz_xxxxyy[k] = -g_y_0_yzz_xxxxyy[k] * ab_z + g_y_0_yzz_xxxxyyz[k];

                g_y_0_yzzz_xxxxyz[k] = -g_y_0_yzz_xxxxyz[k] * ab_z + g_y_0_yzz_xxxxyzz[k];

                g_y_0_yzzz_xxxxzz[k] = -g_y_0_yzz_xxxxzz[k] * ab_z + g_y_0_yzz_xxxxzzz[k];

                g_y_0_yzzz_xxxyyy[k] = -g_y_0_yzz_xxxyyy[k] * ab_z + g_y_0_yzz_xxxyyyz[k];

                g_y_0_yzzz_xxxyyz[k] = -g_y_0_yzz_xxxyyz[k] * ab_z + g_y_0_yzz_xxxyyzz[k];

                g_y_0_yzzz_xxxyzz[k] = -g_y_0_yzz_xxxyzz[k] * ab_z + g_y_0_yzz_xxxyzzz[k];

                g_y_0_yzzz_xxxzzz[k] = -g_y_0_yzz_xxxzzz[k] * ab_z + g_y_0_yzz_xxxzzzz[k];

                g_y_0_yzzz_xxyyyy[k] = -g_y_0_yzz_xxyyyy[k] * ab_z + g_y_0_yzz_xxyyyyz[k];

                g_y_0_yzzz_xxyyyz[k] = -g_y_0_yzz_xxyyyz[k] * ab_z + g_y_0_yzz_xxyyyzz[k];

                g_y_0_yzzz_xxyyzz[k] = -g_y_0_yzz_xxyyzz[k] * ab_z + g_y_0_yzz_xxyyzzz[k];

                g_y_0_yzzz_xxyzzz[k] = -g_y_0_yzz_xxyzzz[k] * ab_z + g_y_0_yzz_xxyzzzz[k];

                g_y_0_yzzz_xxzzzz[k] = -g_y_0_yzz_xxzzzz[k] * ab_z + g_y_0_yzz_xxzzzzz[k];

                g_y_0_yzzz_xyyyyy[k] = -g_y_0_yzz_xyyyyy[k] * ab_z + g_y_0_yzz_xyyyyyz[k];

                g_y_0_yzzz_xyyyyz[k] = -g_y_0_yzz_xyyyyz[k] * ab_z + g_y_0_yzz_xyyyyzz[k];

                g_y_0_yzzz_xyyyzz[k] = -g_y_0_yzz_xyyyzz[k] * ab_z + g_y_0_yzz_xyyyzzz[k];

                g_y_0_yzzz_xyyzzz[k] = -g_y_0_yzz_xyyzzz[k] * ab_z + g_y_0_yzz_xyyzzzz[k];

                g_y_0_yzzz_xyzzzz[k] = -g_y_0_yzz_xyzzzz[k] * ab_z + g_y_0_yzz_xyzzzzz[k];

                g_y_0_yzzz_xzzzzz[k] = -g_y_0_yzz_xzzzzz[k] * ab_z + g_y_0_yzz_xzzzzzz[k];

                g_y_0_yzzz_yyyyyy[k] = -g_y_0_yzz_yyyyyy[k] * ab_z + g_y_0_yzz_yyyyyyz[k];

                g_y_0_yzzz_yyyyyz[k] = -g_y_0_yzz_yyyyyz[k] * ab_z + g_y_0_yzz_yyyyyzz[k];

                g_y_0_yzzz_yyyyzz[k] = -g_y_0_yzz_yyyyzz[k] * ab_z + g_y_0_yzz_yyyyzzz[k];

                g_y_0_yzzz_yyyzzz[k] = -g_y_0_yzz_yyyzzz[k] * ab_z + g_y_0_yzz_yyyzzzz[k];

                g_y_0_yzzz_yyzzzz[k] = -g_y_0_yzz_yyzzzz[k] * ab_z + g_y_0_yzz_yyzzzzz[k];

                g_y_0_yzzz_yzzzzz[k] = -g_y_0_yzz_yzzzzz[k] * ab_z + g_y_0_yzz_yzzzzzz[k];

                g_y_0_yzzz_zzzzzz[k] = -g_y_0_yzz_zzzzzz[k] * ab_z + g_y_0_yzz_zzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxxz, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxyz, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxxzz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyyz, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxyzz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxxzzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyyz, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyyzz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxyzzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxxzzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyyz, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyyzz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyyzzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxyzzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xxzzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyyz, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyyzz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyyzzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyyzzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xyzzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_xzzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyyz, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyyzz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyyzzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyyzzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yyzzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_yzzzzzz, g_y_0_zzz_zzzzzz, g_y_0_zzz_zzzzzzz, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_xxxxxx[k] = -g_y_0_zzz_xxxxxx[k] * ab_z + g_y_0_zzz_xxxxxxz[k];

                g_y_0_zzzz_xxxxxy[k] = -g_y_0_zzz_xxxxxy[k] * ab_z + g_y_0_zzz_xxxxxyz[k];

                g_y_0_zzzz_xxxxxz[k] = -g_y_0_zzz_xxxxxz[k] * ab_z + g_y_0_zzz_xxxxxzz[k];

                g_y_0_zzzz_xxxxyy[k] = -g_y_0_zzz_xxxxyy[k] * ab_z + g_y_0_zzz_xxxxyyz[k];

                g_y_0_zzzz_xxxxyz[k] = -g_y_0_zzz_xxxxyz[k] * ab_z + g_y_0_zzz_xxxxyzz[k];

                g_y_0_zzzz_xxxxzz[k] = -g_y_0_zzz_xxxxzz[k] * ab_z + g_y_0_zzz_xxxxzzz[k];

                g_y_0_zzzz_xxxyyy[k] = -g_y_0_zzz_xxxyyy[k] * ab_z + g_y_0_zzz_xxxyyyz[k];

                g_y_0_zzzz_xxxyyz[k] = -g_y_0_zzz_xxxyyz[k] * ab_z + g_y_0_zzz_xxxyyzz[k];

                g_y_0_zzzz_xxxyzz[k] = -g_y_0_zzz_xxxyzz[k] * ab_z + g_y_0_zzz_xxxyzzz[k];

                g_y_0_zzzz_xxxzzz[k] = -g_y_0_zzz_xxxzzz[k] * ab_z + g_y_0_zzz_xxxzzzz[k];

                g_y_0_zzzz_xxyyyy[k] = -g_y_0_zzz_xxyyyy[k] * ab_z + g_y_0_zzz_xxyyyyz[k];

                g_y_0_zzzz_xxyyyz[k] = -g_y_0_zzz_xxyyyz[k] * ab_z + g_y_0_zzz_xxyyyzz[k];

                g_y_0_zzzz_xxyyzz[k] = -g_y_0_zzz_xxyyzz[k] * ab_z + g_y_0_zzz_xxyyzzz[k];

                g_y_0_zzzz_xxyzzz[k] = -g_y_0_zzz_xxyzzz[k] * ab_z + g_y_0_zzz_xxyzzzz[k];

                g_y_0_zzzz_xxzzzz[k] = -g_y_0_zzz_xxzzzz[k] * ab_z + g_y_0_zzz_xxzzzzz[k];

                g_y_0_zzzz_xyyyyy[k] = -g_y_0_zzz_xyyyyy[k] * ab_z + g_y_0_zzz_xyyyyyz[k];

                g_y_0_zzzz_xyyyyz[k] = -g_y_0_zzz_xyyyyz[k] * ab_z + g_y_0_zzz_xyyyyzz[k];

                g_y_0_zzzz_xyyyzz[k] = -g_y_0_zzz_xyyyzz[k] * ab_z + g_y_0_zzz_xyyyzzz[k];

                g_y_0_zzzz_xyyzzz[k] = -g_y_0_zzz_xyyzzz[k] * ab_z + g_y_0_zzz_xyyzzzz[k];

                g_y_0_zzzz_xyzzzz[k] = -g_y_0_zzz_xyzzzz[k] * ab_z + g_y_0_zzz_xyzzzzz[k];

                g_y_0_zzzz_xzzzzz[k] = -g_y_0_zzz_xzzzzz[k] * ab_z + g_y_0_zzz_xzzzzzz[k];

                g_y_0_zzzz_yyyyyy[k] = -g_y_0_zzz_yyyyyy[k] * ab_z + g_y_0_zzz_yyyyyyz[k];

                g_y_0_zzzz_yyyyyz[k] = -g_y_0_zzz_yyyyyz[k] * ab_z + g_y_0_zzz_yyyyyzz[k];

                g_y_0_zzzz_yyyyzz[k] = -g_y_0_zzz_yyyyzz[k] * ab_z + g_y_0_zzz_yyyyzzz[k];

                g_y_0_zzzz_yyyzzz[k] = -g_y_0_zzz_yyyzzz[k] * ab_z + g_y_0_zzz_yyyzzzz[k];

                g_y_0_zzzz_yyzzzz[k] = -g_y_0_zzz_yyzzzz[k] * ab_z + g_y_0_zzz_yyzzzzz[k];

                g_y_0_zzzz_yzzzzz[k] = -g_y_0_zzz_yzzzzz[k] * ab_z + g_y_0_zzz_yzzzzzz[k];

                g_y_0_zzzz_zzzzzz[k] = -g_y_0_zzz_zzzzzz[k] * ab_z + g_y_0_zzz_zzzzzzz[k];
            }

            /// Set up 840-868 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 854 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 867 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_xxxxxx, g_z_0_xxx_xxxxxxx, g_z_0_xxx_xxxxxxy, g_z_0_xxx_xxxxxxz, g_z_0_xxx_xxxxxy, g_z_0_xxx_xxxxxyy, g_z_0_xxx_xxxxxyz, g_z_0_xxx_xxxxxz, g_z_0_xxx_xxxxxzz, g_z_0_xxx_xxxxyy, g_z_0_xxx_xxxxyyy, g_z_0_xxx_xxxxyyz, g_z_0_xxx_xxxxyz, g_z_0_xxx_xxxxyzz, g_z_0_xxx_xxxxzz, g_z_0_xxx_xxxxzzz, g_z_0_xxx_xxxyyy, g_z_0_xxx_xxxyyyy, g_z_0_xxx_xxxyyyz, g_z_0_xxx_xxxyyz, g_z_0_xxx_xxxyyzz, g_z_0_xxx_xxxyzz, g_z_0_xxx_xxxyzzz, g_z_0_xxx_xxxzzz, g_z_0_xxx_xxxzzzz, g_z_0_xxx_xxyyyy, g_z_0_xxx_xxyyyyy, g_z_0_xxx_xxyyyyz, g_z_0_xxx_xxyyyz, g_z_0_xxx_xxyyyzz, g_z_0_xxx_xxyyzz, g_z_0_xxx_xxyyzzz, g_z_0_xxx_xxyzzz, g_z_0_xxx_xxyzzzz, g_z_0_xxx_xxzzzz, g_z_0_xxx_xxzzzzz, g_z_0_xxx_xyyyyy, g_z_0_xxx_xyyyyyy, g_z_0_xxx_xyyyyyz, g_z_0_xxx_xyyyyz, g_z_0_xxx_xyyyyzz, g_z_0_xxx_xyyyzz, g_z_0_xxx_xyyyzzz, g_z_0_xxx_xyyzzz, g_z_0_xxx_xyyzzzz, g_z_0_xxx_xyzzzz, g_z_0_xxx_xyzzzzz, g_z_0_xxx_xzzzzz, g_z_0_xxx_xzzzzzz, g_z_0_xxx_yyyyyy, g_z_0_xxx_yyyyyz, g_z_0_xxx_yyyyzz, g_z_0_xxx_yyyzzz, g_z_0_xxx_yyzzzz, g_z_0_xxx_yzzzzz, g_z_0_xxx_zzzzzz, g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_yyyyyy, g_z_0_xxxx_yyyyyz, g_z_0_xxxx_yyyyzz, g_z_0_xxxx_yyyzzz, g_z_0_xxxx_yyzzzz, g_z_0_xxxx_yzzzzz, g_z_0_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_xxxxxx[k] = -g_z_0_xxx_xxxxxx[k] * ab_x + g_z_0_xxx_xxxxxxx[k];

                g_z_0_xxxx_xxxxxy[k] = -g_z_0_xxx_xxxxxy[k] * ab_x + g_z_0_xxx_xxxxxxy[k];

                g_z_0_xxxx_xxxxxz[k] = -g_z_0_xxx_xxxxxz[k] * ab_x + g_z_0_xxx_xxxxxxz[k];

                g_z_0_xxxx_xxxxyy[k] = -g_z_0_xxx_xxxxyy[k] * ab_x + g_z_0_xxx_xxxxxyy[k];

                g_z_0_xxxx_xxxxyz[k] = -g_z_0_xxx_xxxxyz[k] * ab_x + g_z_0_xxx_xxxxxyz[k];

                g_z_0_xxxx_xxxxzz[k] = -g_z_0_xxx_xxxxzz[k] * ab_x + g_z_0_xxx_xxxxxzz[k];

                g_z_0_xxxx_xxxyyy[k] = -g_z_0_xxx_xxxyyy[k] * ab_x + g_z_0_xxx_xxxxyyy[k];

                g_z_0_xxxx_xxxyyz[k] = -g_z_0_xxx_xxxyyz[k] * ab_x + g_z_0_xxx_xxxxyyz[k];

                g_z_0_xxxx_xxxyzz[k] = -g_z_0_xxx_xxxyzz[k] * ab_x + g_z_0_xxx_xxxxyzz[k];

                g_z_0_xxxx_xxxzzz[k] = -g_z_0_xxx_xxxzzz[k] * ab_x + g_z_0_xxx_xxxxzzz[k];

                g_z_0_xxxx_xxyyyy[k] = -g_z_0_xxx_xxyyyy[k] * ab_x + g_z_0_xxx_xxxyyyy[k];

                g_z_0_xxxx_xxyyyz[k] = -g_z_0_xxx_xxyyyz[k] * ab_x + g_z_0_xxx_xxxyyyz[k];

                g_z_0_xxxx_xxyyzz[k] = -g_z_0_xxx_xxyyzz[k] * ab_x + g_z_0_xxx_xxxyyzz[k];

                g_z_0_xxxx_xxyzzz[k] = -g_z_0_xxx_xxyzzz[k] * ab_x + g_z_0_xxx_xxxyzzz[k];

                g_z_0_xxxx_xxzzzz[k] = -g_z_0_xxx_xxzzzz[k] * ab_x + g_z_0_xxx_xxxzzzz[k];

                g_z_0_xxxx_xyyyyy[k] = -g_z_0_xxx_xyyyyy[k] * ab_x + g_z_0_xxx_xxyyyyy[k];

                g_z_0_xxxx_xyyyyz[k] = -g_z_0_xxx_xyyyyz[k] * ab_x + g_z_0_xxx_xxyyyyz[k];

                g_z_0_xxxx_xyyyzz[k] = -g_z_0_xxx_xyyyzz[k] * ab_x + g_z_0_xxx_xxyyyzz[k];

                g_z_0_xxxx_xyyzzz[k] = -g_z_0_xxx_xyyzzz[k] * ab_x + g_z_0_xxx_xxyyzzz[k];

                g_z_0_xxxx_xyzzzz[k] = -g_z_0_xxx_xyzzzz[k] * ab_x + g_z_0_xxx_xxyzzzz[k];

                g_z_0_xxxx_xzzzzz[k] = -g_z_0_xxx_xzzzzz[k] * ab_x + g_z_0_xxx_xxzzzzz[k];

                g_z_0_xxxx_yyyyyy[k] = -g_z_0_xxx_yyyyyy[k] * ab_x + g_z_0_xxx_xyyyyyy[k];

                g_z_0_xxxx_yyyyyz[k] = -g_z_0_xxx_yyyyyz[k] * ab_x + g_z_0_xxx_xyyyyyz[k];

                g_z_0_xxxx_yyyyzz[k] = -g_z_0_xxx_yyyyzz[k] * ab_x + g_z_0_xxx_xyyyyzz[k];

                g_z_0_xxxx_yyyzzz[k] = -g_z_0_xxx_yyyzzz[k] * ab_x + g_z_0_xxx_xyyyzzz[k];

                g_z_0_xxxx_yyzzzz[k] = -g_z_0_xxx_yyzzzz[k] * ab_x + g_z_0_xxx_xyyzzzz[k];

                g_z_0_xxxx_yzzzzz[k] = -g_z_0_xxx_yzzzzz[k] * ab_x + g_z_0_xxx_xyzzzzz[k];

                g_z_0_xxxx_zzzzzz[k] = -g_z_0_xxx_zzzzzz[k] * ab_x + g_z_0_xxx_xzzzzzz[k];
            }

            /// Set up 868-896 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 869 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 895 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_yyyyyy, g_z_0_xxxy_yyyyyz, g_z_0_xxxy_yyyyzz, g_z_0_xxxy_yyyzzz, g_z_0_xxxy_yyzzzz, g_z_0_xxxy_yzzzzz, g_z_0_xxxy_zzzzzz, g_z_0_xxy_xxxxxx, g_z_0_xxy_xxxxxxx, g_z_0_xxy_xxxxxxy, g_z_0_xxy_xxxxxxz, g_z_0_xxy_xxxxxy, g_z_0_xxy_xxxxxyy, g_z_0_xxy_xxxxxyz, g_z_0_xxy_xxxxxz, g_z_0_xxy_xxxxxzz, g_z_0_xxy_xxxxyy, g_z_0_xxy_xxxxyyy, g_z_0_xxy_xxxxyyz, g_z_0_xxy_xxxxyz, g_z_0_xxy_xxxxyzz, g_z_0_xxy_xxxxzz, g_z_0_xxy_xxxxzzz, g_z_0_xxy_xxxyyy, g_z_0_xxy_xxxyyyy, g_z_0_xxy_xxxyyyz, g_z_0_xxy_xxxyyz, g_z_0_xxy_xxxyyzz, g_z_0_xxy_xxxyzz, g_z_0_xxy_xxxyzzz, g_z_0_xxy_xxxzzz, g_z_0_xxy_xxxzzzz, g_z_0_xxy_xxyyyy, g_z_0_xxy_xxyyyyy, g_z_0_xxy_xxyyyyz, g_z_0_xxy_xxyyyz, g_z_0_xxy_xxyyyzz, g_z_0_xxy_xxyyzz, g_z_0_xxy_xxyyzzz, g_z_0_xxy_xxyzzz, g_z_0_xxy_xxyzzzz, g_z_0_xxy_xxzzzz, g_z_0_xxy_xxzzzzz, g_z_0_xxy_xyyyyy, g_z_0_xxy_xyyyyyy, g_z_0_xxy_xyyyyyz, g_z_0_xxy_xyyyyz, g_z_0_xxy_xyyyyzz, g_z_0_xxy_xyyyzz, g_z_0_xxy_xyyyzzz, g_z_0_xxy_xyyzzz, g_z_0_xxy_xyyzzzz, g_z_0_xxy_xyzzzz, g_z_0_xxy_xyzzzzz, g_z_0_xxy_xzzzzz, g_z_0_xxy_xzzzzzz, g_z_0_xxy_yyyyyy, g_z_0_xxy_yyyyyz, g_z_0_xxy_yyyyzz, g_z_0_xxy_yyyzzz, g_z_0_xxy_yyzzzz, g_z_0_xxy_yzzzzz, g_z_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_xxxxxx[k] = -g_z_0_xxy_xxxxxx[k] * ab_x + g_z_0_xxy_xxxxxxx[k];

                g_z_0_xxxy_xxxxxy[k] = -g_z_0_xxy_xxxxxy[k] * ab_x + g_z_0_xxy_xxxxxxy[k];

                g_z_0_xxxy_xxxxxz[k] = -g_z_0_xxy_xxxxxz[k] * ab_x + g_z_0_xxy_xxxxxxz[k];

                g_z_0_xxxy_xxxxyy[k] = -g_z_0_xxy_xxxxyy[k] * ab_x + g_z_0_xxy_xxxxxyy[k];

                g_z_0_xxxy_xxxxyz[k] = -g_z_0_xxy_xxxxyz[k] * ab_x + g_z_0_xxy_xxxxxyz[k];

                g_z_0_xxxy_xxxxzz[k] = -g_z_0_xxy_xxxxzz[k] * ab_x + g_z_0_xxy_xxxxxzz[k];

                g_z_0_xxxy_xxxyyy[k] = -g_z_0_xxy_xxxyyy[k] * ab_x + g_z_0_xxy_xxxxyyy[k];

                g_z_0_xxxy_xxxyyz[k] = -g_z_0_xxy_xxxyyz[k] * ab_x + g_z_0_xxy_xxxxyyz[k];

                g_z_0_xxxy_xxxyzz[k] = -g_z_0_xxy_xxxyzz[k] * ab_x + g_z_0_xxy_xxxxyzz[k];

                g_z_0_xxxy_xxxzzz[k] = -g_z_0_xxy_xxxzzz[k] * ab_x + g_z_0_xxy_xxxxzzz[k];

                g_z_0_xxxy_xxyyyy[k] = -g_z_0_xxy_xxyyyy[k] * ab_x + g_z_0_xxy_xxxyyyy[k];

                g_z_0_xxxy_xxyyyz[k] = -g_z_0_xxy_xxyyyz[k] * ab_x + g_z_0_xxy_xxxyyyz[k];

                g_z_0_xxxy_xxyyzz[k] = -g_z_0_xxy_xxyyzz[k] * ab_x + g_z_0_xxy_xxxyyzz[k];

                g_z_0_xxxy_xxyzzz[k] = -g_z_0_xxy_xxyzzz[k] * ab_x + g_z_0_xxy_xxxyzzz[k];

                g_z_0_xxxy_xxzzzz[k] = -g_z_0_xxy_xxzzzz[k] * ab_x + g_z_0_xxy_xxxzzzz[k];

                g_z_0_xxxy_xyyyyy[k] = -g_z_0_xxy_xyyyyy[k] * ab_x + g_z_0_xxy_xxyyyyy[k];

                g_z_0_xxxy_xyyyyz[k] = -g_z_0_xxy_xyyyyz[k] * ab_x + g_z_0_xxy_xxyyyyz[k];

                g_z_0_xxxy_xyyyzz[k] = -g_z_0_xxy_xyyyzz[k] * ab_x + g_z_0_xxy_xxyyyzz[k];

                g_z_0_xxxy_xyyzzz[k] = -g_z_0_xxy_xyyzzz[k] * ab_x + g_z_0_xxy_xxyyzzz[k];

                g_z_0_xxxy_xyzzzz[k] = -g_z_0_xxy_xyzzzz[k] * ab_x + g_z_0_xxy_xxyzzzz[k];

                g_z_0_xxxy_xzzzzz[k] = -g_z_0_xxy_xzzzzz[k] * ab_x + g_z_0_xxy_xxzzzzz[k];

                g_z_0_xxxy_yyyyyy[k] = -g_z_0_xxy_yyyyyy[k] * ab_x + g_z_0_xxy_xyyyyyy[k];

                g_z_0_xxxy_yyyyyz[k] = -g_z_0_xxy_yyyyyz[k] * ab_x + g_z_0_xxy_xyyyyyz[k];

                g_z_0_xxxy_yyyyzz[k] = -g_z_0_xxy_yyyyzz[k] * ab_x + g_z_0_xxy_xyyyyzz[k];

                g_z_0_xxxy_yyyzzz[k] = -g_z_0_xxy_yyyzzz[k] * ab_x + g_z_0_xxy_xyyyzzz[k];

                g_z_0_xxxy_yyzzzz[k] = -g_z_0_xxy_yyzzzz[k] * ab_x + g_z_0_xxy_xyyzzzz[k];

                g_z_0_xxxy_yzzzzz[k] = -g_z_0_xxy_yzzzzz[k] * ab_x + g_z_0_xxy_xyzzzzz[k];

                g_z_0_xxxy_zzzzzz[k] = -g_z_0_xxy_zzzzzz[k] * ab_x + g_z_0_xxy_xzzzzzz[k];
            }

            /// Set up 896-924 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_yyyyyy, g_z_0_xxxz_yyyyyz, g_z_0_xxxz_yyyyzz, g_z_0_xxxz_yyyzzz, g_z_0_xxxz_yyzzzz, g_z_0_xxxz_yzzzzz, g_z_0_xxxz_zzzzzz, g_z_0_xxz_xxxxxx, g_z_0_xxz_xxxxxxx, g_z_0_xxz_xxxxxxy, g_z_0_xxz_xxxxxxz, g_z_0_xxz_xxxxxy, g_z_0_xxz_xxxxxyy, g_z_0_xxz_xxxxxyz, g_z_0_xxz_xxxxxz, g_z_0_xxz_xxxxxzz, g_z_0_xxz_xxxxyy, g_z_0_xxz_xxxxyyy, g_z_0_xxz_xxxxyyz, g_z_0_xxz_xxxxyz, g_z_0_xxz_xxxxyzz, g_z_0_xxz_xxxxzz, g_z_0_xxz_xxxxzzz, g_z_0_xxz_xxxyyy, g_z_0_xxz_xxxyyyy, g_z_0_xxz_xxxyyyz, g_z_0_xxz_xxxyyz, g_z_0_xxz_xxxyyzz, g_z_0_xxz_xxxyzz, g_z_0_xxz_xxxyzzz, g_z_0_xxz_xxxzzz, g_z_0_xxz_xxxzzzz, g_z_0_xxz_xxyyyy, g_z_0_xxz_xxyyyyy, g_z_0_xxz_xxyyyyz, g_z_0_xxz_xxyyyz, g_z_0_xxz_xxyyyzz, g_z_0_xxz_xxyyzz, g_z_0_xxz_xxyyzzz, g_z_0_xxz_xxyzzz, g_z_0_xxz_xxyzzzz, g_z_0_xxz_xxzzzz, g_z_0_xxz_xxzzzzz, g_z_0_xxz_xyyyyy, g_z_0_xxz_xyyyyyy, g_z_0_xxz_xyyyyyz, g_z_0_xxz_xyyyyz, g_z_0_xxz_xyyyyzz, g_z_0_xxz_xyyyzz, g_z_0_xxz_xyyyzzz, g_z_0_xxz_xyyzzz, g_z_0_xxz_xyyzzzz, g_z_0_xxz_xyzzzz, g_z_0_xxz_xyzzzzz, g_z_0_xxz_xzzzzz, g_z_0_xxz_xzzzzzz, g_z_0_xxz_yyyyyy, g_z_0_xxz_yyyyyz, g_z_0_xxz_yyyyzz, g_z_0_xxz_yyyzzz, g_z_0_xxz_yyzzzz, g_z_0_xxz_yzzzzz, g_z_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_xxxxxx[k] = -g_z_0_xxz_xxxxxx[k] * ab_x + g_z_0_xxz_xxxxxxx[k];

                g_z_0_xxxz_xxxxxy[k] = -g_z_0_xxz_xxxxxy[k] * ab_x + g_z_0_xxz_xxxxxxy[k];

                g_z_0_xxxz_xxxxxz[k] = -g_z_0_xxz_xxxxxz[k] * ab_x + g_z_0_xxz_xxxxxxz[k];

                g_z_0_xxxz_xxxxyy[k] = -g_z_0_xxz_xxxxyy[k] * ab_x + g_z_0_xxz_xxxxxyy[k];

                g_z_0_xxxz_xxxxyz[k] = -g_z_0_xxz_xxxxyz[k] * ab_x + g_z_0_xxz_xxxxxyz[k];

                g_z_0_xxxz_xxxxzz[k] = -g_z_0_xxz_xxxxzz[k] * ab_x + g_z_0_xxz_xxxxxzz[k];

                g_z_0_xxxz_xxxyyy[k] = -g_z_0_xxz_xxxyyy[k] * ab_x + g_z_0_xxz_xxxxyyy[k];

                g_z_0_xxxz_xxxyyz[k] = -g_z_0_xxz_xxxyyz[k] * ab_x + g_z_0_xxz_xxxxyyz[k];

                g_z_0_xxxz_xxxyzz[k] = -g_z_0_xxz_xxxyzz[k] * ab_x + g_z_0_xxz_xxxxyzz[k];

                g_z_0_xxxz_xxxzzz[k] = -g_z_0_xxz_xxxzzz[k] * ab_x + g_z_0_xxz_xxxxzzz[k];

                g_z_0_xxxz_xxyyyy[k] = -g_z_0_xxz_xxyyyy[k] * ab_x + g_z_0_xxz_xxxyyyy[k];

                g_z_0_xxxz_xxyyyz[k] = -g_z_0_xxz_xxyyyz[k] * ab_x + g_z_0_xxz_xxxyyyz[k];

                g_z_0_xxxz_xxyyzz[k] = -g_z_0_xxz_xxyyzz[k] * ab_x + g_z_0_xxz_xxxyyzz[k];

                g_z_0_xxxz_xxyzzz[k] = -g_z_0_xxz_xxyzzz[k] * ab_x + g_z_0_xxz_xxxyzzz[k];

                g_z_0_xxxz_xxzzzz[k] = -g_z_0_xxz_xxzzzz[k] * ab_x + g_z_0_xxz_xxxzzzz[k];

                g_z_0_xxxz_xyyyyy[k] = -g_z_0_xxz_xyyyyy[k] * ab_x + g_z_0_xxz_xxyyyyy[k];

                g_z_0_xxxz_xyyyyz[k] = -g_z_0_xxz_xyyyyz[k] * ab_x + g_z_0_xxz_xxyyyyz[k];

                g_z_0_xxxz_xyyyzz[k] = -g_z_0_xxz_xyyyzz[k] * ab_x + g_z_0_xxz_xxyyyzz[k];

                g_z_0_xxxz_xyyzzz[k] = -g_z_0_xxz_xyyzzz[k] * ab_x + g_z_0_xxz_xxyyzzz[k];

                g_z_0_xxxz_xyzzzz[k] = -g_z_0_xxz_xyzzzz[k] * ab_x + g_z_0_xxz_xxyzzzz[k];

                g_z_0_xxxz_xzzzzz[k] = -g_z_0_xxz_xzzzzz[k] * ab_x + g_z_0_xxz_xxzzzzz[k];

                g_z_0_xxxz_yyyyyy[k] = -g_z_0_xxz_yyyyyy[k] * ab_x + g_z_0_xxz_xyyyyyy[k];

                g_z_0_xxxz_yyyyyz[k] = -g_z_0_xxz_yyyyyz[k] * ab_x + g_z_0_xxz_xyyyyyz[k];

                g_z_0_xxxz_yyyyzz[k] = -g_z_0_xxz_yyyyzz[k] * ab_x + g_z_0_xxz_xyyyyzz[k];

                g_z_0_xxxz_yyyzzz[k] = -g_z_0_xxz_yyyzzz[k] * ab_x + g_z_0_xxz_xyyyzzz[k];

                g_z_0_xxxz_yyzzzz[k] = -g_z_0_xxz_yyzzzz[k] * ab_x + g_z_0_xxz_xyyzzzz[k];

                g_z_0_xxxz_yzzzzz[k] = -g_z_0_xxz_yzzzzz[k] * ab_x + g_z_0_xxz_xyzzzzz[k];

                g_z_0_xxxz_zzzzzz[k] = -g_z_0_xxz_zzzzzz[k] * ab_x + g_z_0_xxz_xzzzzzz[k];
            }

            /// Set up 924-952 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 944 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 951 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_yyyyyy, g_z_0_xxyy_yyyyyz, g_z_0_xxyy_yyyyzz, g_z_0_xxyy_yyyzzz, g_z_0_xxyy_yyzzzz, g_z_0_xxyy_yzzzzz, g_z_0_xxyy_zzzzzz, g_z_0_xyy_xxxxxx, g_z_0_xyy_xxxxxxx, g_z_0_xyy_xxxxxxy, g_z_0_xyy_xxxxxxz, g_z_0_xyy_xxxxxy, g_z_0_xyy_xxxxxyy, g_z_0_xyy_xxxxxyz, g_z_0_xyy_xxxxxz, g_z_0_xyy_xxxxxzz, g_z_0_xyy_xxxxyy, g_z_0_xyy_xxxxyyy, g_z_0_xyy_xxxxyyz, g_z_0_xyy_xxxxyz, g_z_0_xyy_xxxxyzz, g_z_0_xyy_xxxxzz, g_z_0_xyy_xxxxzzz, g_z_0_xyy_xxxyyy, g_z_0_xyy_xxxyyyy, g_z_0_xyy_xxxyyyz, g_z_0_xyy_xxxyyz, g_z_0_xyy_xxxyyzz, g_z_0_xyy_xxxyzz, g_z_0_xyy_xxxyzzz, g_z_0_xyy_xxxzzz, g_z_0_xyy_xxxzzzz, g_z_0_xyy_xxyyyy, g_z_0_xyy_xxyyyyy, g_z_0_xyy_xxyyyyz, g_z_0_xyy_xxyyyz, g_z_0_xyy_xxyyyzz, g_z_0_xyy_xxyyzz, g_z_0_xyy_xxyyzzz, g_z_0_xyy_xxyzzz, g_z_0_xyy_xxyzzzz, g_z_0_xyy_xxzzzz, g_z_0_xyy_xxzzzzz, g_z_0_xyy_xyyyyy, g_z_0_xyy_xyyyyyy, g_z_0_xyy_xyyyyyz, g_z_0_xyy_xyyyyz, g_z_0_xyy_xyyyyzz, g_z_0_xyy_xyyyzz, g_z_0_xyy_xyyyzzz, g_z_0_xyy_xyyzzz, g_z_0_xyy_xyyzzzz, g_z_0_xyy_xyzzzz, g_z_0_xyy_xyzzzzz, g_z_0_xyy_xzzzzz, g_z_0_xyy_xzzzzzz, g_z_0_xyy_yyyyyy, g_z_0_xyy_yyyyyz, g_z_0_xyy_yyyyzz, g_z_0_xyy_yyyzzz, g_z_0_xyy_yyzzzz, g_z_0_xyy_yzzzzz, g_z_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_xxxxxx[k] = -g_z_0_xyy_xxxxxx[k] * ab_x + g_z_0_xyy_xxxxxxx[k];

                g_z_0_xxyy_xxxxxy[k] = -g_z_0_xyy_xxxxxy[k] * ab_x + g_z_0_xyy_xxxxxxy[k];

                g_z_0_xxyy_xxxxxz[k] = -g_z_0_xyy_xxxxxz[k] * ab_x + g_z_0_xyy_xxxxxxz[k];

                g_z_0_xxyy_xxxxyy[k] = -g_z_0_xyy_xxxxyy[k] * ab_x + g_z_0_xyy_xxxxxyy[k];

                g_z_0_xxyy_xxxxyz[k] = -g_z_0_xyy_xxxxyz[k] * ab_x + g_z_0_xyy_xxxxxyz[k];

                g_z_0_xxyy_xxxxzz[k] = -g_z_0_xyy_xxxxzz[k] * ab_x + g_z_0_xyy_xxxxxzz[k];

                g_z_0_xxyy_xxxyyy[k] = -g_z_0_xyy_xxxyyy[k] * ab_x + g_z_0_xyy_xxxxyyy[k];

                g_z_0_xxyy_xxxyyz[k] = -g_z_0_xyy_xxxyyz[k] * ab_x + g_z_0_xyy_xxxxyyz[k];

                g_z_0_xxyy_xxxyzz[k] = -g_z_0_xyy_xxxyzz[k] * ab_x + g_z_0_xyy_xxxxyzz[k];

                g_z_0_xxyy_xxxzzz[k] = -g_z_0_xyy_xxxzzz[k] * ab_x + g_z_0_xyy_xxxxzzz[k];

                g_z_0_xxyy_xxyyyy[k] = -g_z_0_xyy_xxyyyy[k] * ab_x + g_z_0_xyy_xxxyyyy[k];

                g_z_0_xxyy_xxyyyz[k] = -g_z_0_xyy_xxyyyz[k] * ab_x + g_z_0_xyy_xxxyyyz[k];

                g_z_0_xxyy_xxyyzz[k] = -g_z_0_xyy_xxyyzz[k] * ab_x + g_z_0_xyy_xxxyyzz[k];

                g_z_0_xxyy_xxyzzz[k] = -g_z_0_xyy_xxyzzz[k] * ab_x + g_z_0_xyy_xxxyzzz[k];

                g_z_0_xxyy_xxzzzz[k] = -g_z_0_xyy_xxzzzz[k] * ab_x + g_z_0_xyy_xxxzzzz[k];

                g_z_0_xxyy_xyyyyy[k] = -g_z_0_xyy_xyyyyy[k] * ab_x + g_z_0_xyy_xxyyyyy[k];

                g_z_0_xxyy_xyyyyz[k] = -g_z_0_xyy_xyyyyz[k] * ab_x + g_z_0_xyy_xxyyyyz[k];

                g_z_0_xxyy_xyyyzz[k] = -g_z_0_xyy_xyyyzz[k] * ab_x + g_z_0_xyy_xxyyyzz[k];

                g_z_0_xxyy_xyyzzz[k] = -g_z_0_xyy_xyyzzz[k] * ab_x + g_z_0_xyy_xxyyzzz[k];

                g_z_0_xxyy_xyzzzz[k] = -g_z_0_xyy_xyzzzz[k] * ab_x + g_z_0_xyy_xxyzzzz[k];

                g_z_0_xxyy_xzzzzz[k] = -g_z_0_xyy_xzzzzz[k] * ab_x + g_z_0_xyy_xxzzzzz[k];

                g_z_0_xxyy_yyyyyy[k] = -g_z_0_xyy_yyyyyy[k] * ab_x + g_z_0_xyy_xyyyyyy[k];

                g_z_0_xxyy_yyyyyz[k] = -g_z_0_xyy_yyyyyz[k] * ab_x + g_z_0_xyy_xyyyyyz[k];

                g_z_0_xxyy_yyyyzz[k] = -g_z_0_xyy_yyyyzz[k] * ab_x + g_z_0_xyy_xyyyyzz[k];

                g_z_0_xxyy_yyyzzz[k] = -g_z_0_xyy_yyyzzz[k] * ab_x + g_z_0_xyy_xyyyzzz[k];

                g_z_0_xxyy_yyzzzz[k] = -g_z_0_xyy_yyzzzz[k] * ab_x + g_z_0_xyy_xyyzzzz[k];

                g_z_0_xxyy_yzzzzz[k] = -g_z_0_xyy_yzzzzz[k] * ab_x + g_z_0_xyy_xyzzzzz[k];

                g_z_0_xxyy_zzzzzz[k] = -g_z_0_xyy_zzzzzz[k] * ab_x + g_z_0_xyy_xzzzzzz[k];
            }

            /// Set up 952-980 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 959 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 974 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_yyyyyy, g_z_0_xxyz_yyyyyz, g_z_0_xxyz_yyyyzz, g_z_0_xxyz_yyyzzz, g_z_0_xxyz_yyzzzz, g_z_0_xxyz_yzzzzz, g_z_0_xxyz_zzzzzz, g_z_0_xyz_xxxxxx, g_z_0_xyz_xxxxxxx, g_z_0_xyz_xxxxxxy, g_z_0_xyz_xxxxxxz, g_z_0_xyz_xxxxxy, g_z_0_xyz_xxxxxyy, g_z_0_xyz_xxxxxyz, g_z_0_xyz_xxxxxz, g_z_0_xyz_xxxxxzz, g_z_0_xyz_xxxxyy, g_z_0_xyz_xxxxyyy, g_z_0_xyz_xxxxyyz, g_z_0_xyz_xxxxyz, g_z_0_xyz_xxxxyzz, g_z_0_xyz_xxxxzz, g_z_0_xyz_xxxxzzz, g_z_0_xyz_xxxyyy, g_z_0_xyz_xxxyyyy, g_z_0_xyz_xxxyyyz, g_z_0_xyz_xxxyyz, g_z_0_xyz_xxxyyzz, g_z_0_xyz_xxxyzz, g_z_0_xyz_xxxyzzz, g_z_0_xyz_xxxzzz, g_z_0_xyz_xxxzzzz, g_z_0_xyz_xxyyyy, g_z_0_xyz_xxyyyyy, g_z_0_xyz_xxyyyyz, g_z_0_xyz_xxyyyz, g_z_0_xyz_xxyyyzz, g_z_0_xyz_xxyyzz, g_z_0_xyz_xxyyzzz, g_z_0_xyz_xxyzzz, g_z_0_xyz_xxyzzzz, g_z_0_xyz_xxzzzz, g_z_0_xyz_xxzzzzz, g_z_0_xyz_xyyyyy, g_z_0_xyz_xyyyyyy, g_z_0_xyz_xyyyyyz, g_z_0_xyz_xyyyyz, g_z_0_xyz_xyyyyzz, g_z_0_xyz_xyyyzz, g_z_0_xyz_xyyyzzz, g_z_0_xyz_xyyzzz, g_z_0_xyz_xyyzzzz, g_z_0_xyz_xyzzzz, g_z_0_xyz_xyzzzzz, g_z_0_xyz_xzzzzz, g_z_0_xyz_xzzzzzz, g_z_0_xyz_yyyyyy, g_z_0_xyz_yyyyyz, g_z_0_xyz_yyyyzz, g_z_0_xyz_yyyzzz, g_z_0_xyz_yyzzzz, g_z_0_xyz_yzzzzz, g_z_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_xxxxxx[k] = -g_z_0_xyz_xxxxxx[k] * ab_x + g_z_0_xyz_xxxxxxx[k];

                g_z_0_xxyz_xxxxxy[k] = -g_z_0_xyz_xxxxxy[k] * ab_x + g_z_0_xyz_xxxxxxy[k];

                g_z_0_xxyz_xxxxxz[k] = -g_z_0_xyz_xxxxxz[k] * ab_x + g_z_0_xyz_xxxxxxz[k];

                g_z_0_xxyz_xxxxyy[k] = -g_z_0_xyz_xxxxyy[k] * ab_x + g_z_0_xyz_xxxxxyy[k];

                g_z_0_xxyz_xxxxyz[k] = -g_z_0_xyz_xxxxyz[k] * ab_x + g_z_0_xyz_xxxxxyz[k];

                g_z_0_xxyz_xxxxzz[k] = -g_z_0_xyz_xxxxzz[k] * ab_x + g_z_0_xyz_xxxxxzz[k];

                g_z_0_xxyz_xxxyyy[k] = -g_z_0_xyz_xxxyyy[k] * ab_x + g_z_0_xyz_xxxxyyy[k];

                g_z_0_xxyz_xxxyyz[k] = -g_z_0_xyz_xxxyyz[k] * ab_x + g_z_0_xyz_xxxxyyz[k];

                g_z_0_xxyz_xxxyzz[k] = -g_z_0_xyz_xxxyzz[k] * ab_x + g_z_0_xyz_xxxxyzz[k];

                g_z_0_xxyz_xxxzzz[k] = -g_z_0_xyz_xxxzzz[k] * ab_x + g_z_0_xyz_xxxxzzz[k];

                g_z_0_xxyz_xxyyyy[k] = -g_z_0_xyz_xxyyyy[k] * ab_x + g_z_0_xyz_xxxyyyy[k];

                g_z_0_xxyz_xxyyyz[k] = -g_z_0_xyz_xxyyyz[k] * ab_x + g_z_0_xyz_xxxyyyz[k];

                g_z_0_xxyz_xxyyzz[k] = -g_z_0_xyz_xxyyzz[k] * ab_x + g_z_0_xyz_xxxyyzz[k];

                g_z_0_xxyz_xxyzzz[k] = -g_z_0_xyz_xxyzzz[k] * ab_x + g_z_0_xyz_xxxyzzz[k];

                g_z_0_xxyz_xxzzzz[k] = -g_z_0_xyz_xxzzzz[k] * ab_x + g_z_0_xyz_xxxzzzz[k];

                g_z_0_xxyz_xyyyyy[k] = -g_z_0_xyz_xyyyyy[k] * ab_x + g_z_0_xyz_xxyyyyy[k];

                g_z_0_xxyz_xyyyyz[k] = -g_z_0_xyz_xyyyyz[k] * ab_x + g_z_0_xyz_xxyyyyz[k];

                g_z_0_xxyz_xyyyzz[k] = -g_z_0_xyz_xyyyzz[k] * ab_x + g_z_0_xyz_xxyyyzz[k];

                g_z_0_xxyz_xyyzzz[k] = -g_z_0_xyz_xyyzzz[k] * ab_x + g_z_0_xyz_xxyyzzz[k];

                g_z_0_xxyz_xyzzzz[k] = -g_z_0_xyz_xyzzzz[k] * ab_x + g_z_0_xyz_xxyzzzz[k];

                g_z_0_xxyz_xzzzzz[k] = -g_z_0_xyz_xzzzzz[k] * ab_x + g_z_0_xyz_xxzzzzz[k];

                g_z_0_xxyz_yyyyyy[k] = -g_z_0_xyz_yyyyyy[k] * ab_x + g_z_0_xyz_xyyyyyy[k];

                g_z_0_xxyz_yyyyyz[k] = -g_z_0_xyz_yyyyyz[k] * ab_x + g_z_0_xyz_xyyyyyz[k];

                g_z_0_xxyz_yyyyzz[k] = -g_z_0_xyz_yyyyzz[k] * ab_x + g_z_0_xyz_xyyyyzz[k];

                g_z_0_xxyz_yyyzzz[k] = -g_z_0_xyz_yyyzzz[k] * ab_x + g_z_0_xyz_xyyyzzz[k];

                g_z_0_xxyz_yyzzzz[k] = -g_z_0_xyz_yyzzzz[k] * ab_x + g_z_0_xyz_xyyzzzz[k];

                g_z_0_xxyz_yzzzzz[k] = -g_z_0_xyz_yzzzzz[k] * ab_x + g_z_0_xyz_xyzzzzz[k];

                g_z_0_xxyz_zzzzzz[k] = -g_z_0_xyz_zzzzzz[k] * ab_x + g_z_0_xyz_xzzzzzz[k];
            }

            /// Set up 980-1008 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 989 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1004 * ccomps * dcomps);

            auto g_z_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_yyyyyy, g_z_0_xxzz_yyyyyz, g_z_0_xxzz_yyyyzz, g_z_0_xxzz_yyyzzz, g_z_0_xxzz_yyzzzz, g_z_0_xxzz_yzzzzz, g_z_0_xxzz_zzzzzz, g_z_0_xzz_xxxxxx, g_z_0_xzz_xxxxxxx, g_z_0_xzz_xxxxxxy, g_z_0_xzz_xxxxxxz, g_z_0_xzz_xxxxxy, g_z_0_xzz_xxxxxyy, g_z_0_xzz_xxxxxyz, g_z_0_xzz_xxxxxz, g_z_0_xzz_xxxxxzz, g_z_0_xzz_xxxxyy, g_z_0_xzz_xxxxyyy, g_z_0_xzz_xxxxyyz, g_z_0_xzz_xxxxyz, g_z_0_xzz_xxxxyzz, g_z_0_xzz_xxxxzz, g_z_0_xzz_xxxxzzz, g_z_0_xzz_xxxyyy, g_z_0_xzz_xxxyyyy, g_z_0_xzz_xxxyyyz, g_z_0_xzz_xxxyyz, g_z_0_xzz_xxxyyzz, g_z_0_xzz_xxxyzz, g_z_0_xzz_xxxyzzz, g_z_0_xzz_xxxzzz, g_z_0_xzz_xxxzzzz, g_z_0_xzz_xxyyyy, g_z_0_xzz_xxyyyyy, g_z_0_xzz_xxyyyyz, g_z_0_xzz_xxyyyz, g_z_0_xzz_xxyyyzz, g_z_0_xzz_xxyyzz, g_z_0_xzz_xxyyzzz, g_z_0_xzz_xxyzzz, g_z_0_xzz_xxyzzzz, g_z_0_xzz_xxzzzz, g_z_0_xzz_xxzzzzz, g_z_0_xzz_xyyyyy, g_z_0_xzz_xyyyyyy, g_z_0_xzz_xyyyyyz, g_z_0_xzz_xyyyyz, g_z_0_xzz_xyyyyzz, g_z_0_xzz_xyyyzz, g_z_0_xzz_xyyyzzz, g_z_0_xzz_xyyzzz, g_z_0_xzz_xyyzzzz, g_z_0_xzz_xyzzzz, g_z_0_xzz_xyzzzzz, g_z_0_xzz_xzzzzz, g_z_0_xzz_xzzzzzz, g_z_0_xzz_yyyyyy, g_z_0_xzz_yyyyyz, g_z_0_xzz_yyyyzz, g_z_0_xzz_yyyzzz, g_z_0_xzz_yyzzzz, g_z_0_xzz_yzzzzz, g_z_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_xxxxxx[k] = -g_z_0_xzz_xxxxxx[k] * ab_x + g_z_0_xzz_xxxxxxx[k];

                g_z_0_xxzz_xxxxxy[k] = -g_z_0_xzz_xxxxxy[k] * ab_x + g_z_0_xzz_xxxxxxy[k];

                g_z_0_xxzz_xxxxxz[k] = -g_z_0_xzz_xxxxxz[k] * ab_x + g_z_0_xzz_xxxxxxz[k];

                g_z_0_xxzz_xxxxyy[k] = -g_z_0_xzz_xxxxyy[k] * ab_x + g_z_0_xzz_xxxxxyy[k];

                g_z_0_xxzz_xxxxyz[k] = -g_z_0_xzz_xxxxyz[k] * ab_x + g_z_0_xzz_xxxxxyz[k];

                g_z_0_xxzz_xxxxzz[k] = -g_z_0_xzz_xxxxzz[k] * ab_x + g_z_0_xzz_xxxxxzz[k];

                g_z_0_xxzz_xxxyyy[k] = -g_z_0_xzz_xxxyyy[k] * ab_x + g_z_0_xzz_xxxxyyy[k];

                g_z_0_xxzz_xxxyyz[k] = -g_z_0_xzz_xxxyyz[k] * ab_x + g_z_0_xzz_xxxxyyz[k];

                g_z_0_xxzz_xxxyzz[k] = -g_z_0_xzz_xxxyzz[k] * ab_x + g_z_0_xzz_xxxxyzz[k];

                g_z_0_xxzz_xxxzzz[k] = -g_z_0_xzz_xxxzzz[k] * ab_x + g_z_0_xzz_xxxxzzz[k];

                g_z_0_xxzz_xxyyyy[k] = -g_z_0_xzz_xxyyyy[k] * ab_x + g_z_0_xzz_xxxyyyy[k];

                g_z_0_xxzz_xxyyyz[k] = -g_z_0_xzz_xxyyyz[k] * ab_x + g_z_0_xzz_xxxyyyz[k];

                g_z_0_xxzz_xxyyzz[k] = -g_z_0_xzz_xxyyzz[k] * ab_x + g_z_0_xzz_xxxyyzz[k];

                g_z_0_xxzz_xxyzzz[k] = -g_z_0_xzz_xxyzzz[k] * ab_x + g_z_0_xzz_xxxyzzz[k];

                g_z_0_xxzz_xxzzzz[k] = -g_z_0_xzz_xxzzzz[k] * ab_x + g_z_0_xzz_xxxzzzz[k];

                g_z_0_xxzz_xyyyyy[k] = -g_z_0_xzz_xyyyyy[k] * ab_x + g_z_0_xzz_xxyyyyy[k];

                g_z_0_xxzz_xyyyyz[k] = -g_z_0_xzz_xyyyyz[k] * ab_x + g_z_0_xzz_xxyyyyz[k];

                g_z_0_xxzz_xyyyzz[k] = -g_z_0_xzz_xyyyzz[k] * ab_x + g_z_0_xzz_xxyyyzz[k];

                g_z_0_xxzz_xyyzzz[k] = -g_z_0_xzz_xyyzzz[k] * ab_x + g_z_0_xzz_xxyyzzz[k];

                g_z_0_xxzz_xyzzzz[k] = -g_z_0_xzz_xyzzzz[k] * ab_x + g_z_0_xzz_xxyzzzz[k];

                g_z_0_xxzz_xzzzzz[k] = -g_z_0_xzz_xzzzzz[k] * ab_x + g_z_0_xzz_xxzzzzz[k];

                g_z_0_xxzz_yyyyyy[k] = -g_z_0_xzz_yyyyyy[k] * ab_x + g_z_0_xzz_xyyyyyy[k];

                g_z_0_xxzz_yyyyyz[k] = -g_z_0_xzz_yyyyyz[k] * ab_x + g_z_0_xzz_xyyyyyz[k];

                g_z_0_xxzz_yyyyzz[k] = -g_z_0_xzz_yyyyzz[k] * ab_x + g_z_0_xzz_xyyyyzz[k];

                g_z_0_xxzz_yyyzzz[k] = -g_z_0_xzz_yyyzzz[k] * ab_x + g_z_0_xzz_xyyyzzz[k];

                g_z_0_xxzz_yyzzzz[k] = -g_z_0_xzz_yyzzzz[k] * ab_x + g_z_0_xzz_xyyzzzz[k];

                g_z_0_xxzz_yzzzzz[k] = -g_z_0_xzz_yzzzzz[k] * ab_x + g_z_0_xzz_xyzzzzz[k];

                g_z_0_xxzz_zzzzzz[k] = -g_z_0_xzz_zzzzzz[k] * ab_x + g_z_0_xzz_xzzzzzz[k];
            }

            /// Set up 1008-1036 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 1019 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 1034 * ccomps * dcomps);

            auto g_z_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 1035 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_yyyyyy, g_z_0_xyyy_yyyyyz, g_z_0_xyyy_yyyyzz, g_z_0_xyyy_yyyzzz, g_z_0_xyyy_yyzzzz, g_z_0_xyyy_yzzzzz, g_z_0_xyyy_zzzzzz, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxxx, g_z_0_yyy_xxxxxxy, g_z_0_yyy_xxxxxxz, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxyy, g_z_0_yyy_xxxxxyz, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxxzz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyyy, g_z_0_yyy_xxxxyyz, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxyzz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxxzzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyyy, g_z_0_yyy_xxxyyyz, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyyzz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxyzzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxxzzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyyy, g_z_0_yyy_xxyyyyz, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyyzz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyyzzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxyzzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xxzzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyyy, g_z_0_yyy_xyyyyyz, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyyzz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyyzzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyyzzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xyzzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_xzzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_xxxxxx[k] = -g_z_0_yyy_xxxxxx[k] * ab_x + g_z_0_yyy_xxxxxxx[k];

                g_z_0_xyyy_xxxxxy[k] = -g_z_0_yyy_xxxxxy[k] * ab_x + g_z_0_yyy_xxxxxxy[k];

                g_z_0_xyyy_xxxxxz[k] = -g_z_0_yyy_xxxxxz[k] * ab_x + g_z_0_yyy_xxxxxxz[k];

                g_z_0_xyyy_xxxxyy[k] = -g_z_0_yyy_xxxxyy[k] * ab_x + g_z_0_yyy_xxxxxyy[k];

                g_z_0_xyyy_xxxxyz[k] = -g_z_0_yyy_xxxxyz[k] * ab_x + g_z_0_yyy_xxxxxyz[k];

                g_z_0_xyyy_xxxxzz[k] = -g_z_0_yyy_xxxxzz[k] * ab_x + g_z_0_yyy_xxxxxzz[k];

                g_z_0_xyyy_xxxyyy[k] = -g_z_0_yyy_xxxyyy[k] * ab_x + g_z_0_yyy_xxxxyyy[k];

                g_z_0_xyyy_xxxyyz[k] = -g_z_0_yyy_xxxyyz[k] * ab_x + g_z_0_yyy_xxxxyyz[k];

                g_z_0_xyyy_xxxyzz[k] = -g_z_0_yyy_xxxyzz[k] * ab_x + g_z_0_yyy_xxxxyzz[k];

                g_z_0_xyyy_xxxzzz[k] = -g_z_0_yyy_xxxzzz[k] * ab_x + g_z_0_yyy_xxxxzzz[k];

                g_z_0_xyyy_xxyyyy[k] = -g_z_0_yyy_xxyyyy[k] * ab_x + g_z_0_yyy_xxxyyyy[k];

                g_z_0_xyyy_xxyyyz[k] = -g_z_0_yyy_xxyyyz[k] * ab_x + g_z_0_yyy_xxxyyyz[k];

                g_z_0_xyyy_xxyyzz[k] = -g_z_0_yyy_xxyyzz[k] * ab_x + g_z_0_yyy_xxxyyzz[k];

                g_z_0_xyyy_xxyzzz[k] = -g_z_0_yyy_xxyzzz[k] * ab_x + g_z_0_yyy_xxxyzzz[k];

                g_z_0_xyyy_xxzzzz[k] = -g_z_0_yyy_xxzzzz[k] * ab_x + g_z_0_yyy_xxxzzzz[k];

                g_z_0_xyyy_xyyyyy[k] = -g_z_0_yyy_xyyyyy[k] * ab_x + g_z_0_yyy_xxyyyyy[k];

                g_z_0_xyyy_xyyyyz[k] = -g_z_0_yyy_xyyyyz[k] * ab_x + g_z_0_yyy_xxyyyyz[k];

                g_z_0_xyyy_xyyyzz[k] = -g_z_0_yyy_xyyyzz[k] * ab_x + g_z_0_yyy_xxyyyzz[k];

                g_z_0_xyyy_xyyzzz[k] = -g_z_0_yyy_xyyzzz[k] * ab_x + g_z_0_yyy_xxyyzzz[k];

                g_z_0_xyyy_xyzzzz[k] = -g_z_0_yyy_xyzzzz[k] * ab_x + g_z_0_yyy_xxyzzzz[k];

                g_z_0_xyyy_xzzzzz[k] = -g_z_0_yyy_xzzzzz[k] * ab_x + g_z_0_yyy_xxzzzzz[k];

                g_z_0_xyyy_yyyyyy[k] = -g_z_0_yyy_yyyyyy[k] * ab_x + g_z_0_yyy_xyyyyyy[k];

                g_z_0_xyyy_yyyyyz[k] = -g_z_0_yyy_yyyyyz[k] * ab_x + g_z_0_yyy_xyyyyyz[k];

                g_z_0_xyyy_yyyyzz[k] = -g_z_0_yyy_yyyyzz[k] * ab_x + g_z_0_yyy_xyyyyzz[k];

                g_z_0_xyyy_yyyzzz[k] = -g_z_0_yyy_yyyzzz[k] * ab_x + g_z_0_yyy_xyyyzzz[k];

                g_z_0_xyyy_yyzzzz[k] = -g_z_0_yyy_yyzzzz[k] * ab_x + g_z_0_yyy_xyyzzzz[k];

                g_z_0_xyyy_yzzzzz[k] = -g_z_0_yyy_yzzzzz[k] * ab_x + g_z_0_yyy_xyzzzzz[k];

                g_z_0_xyyy_zzzzzz[k] = -g_z_0_yyy_zzzzzz[k] * ab_x + g_z_0_yyy_xzzzzzz[k];
            }

            /// Set up 1036-1064 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 1049 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 1063 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_yyyyyy, g_z_0_xyyz_yyyyyz, g_z_0_xyyz_yyyyzz, g_z_0_xyyz_yyyzzz, g_z_0_xyyz_yyzzzz, g_z_0_xyyz_yzzzzz, g_z_0_xyyz_zzzzzz, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxxx, g_z_0_yyz_xxxxxxy, g_z_0_yyz_xxxxxxz, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxyy, g_z_0_yyz_xxxxxyz, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxxzz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyyy, g_z_0_yyz_xxxxyyz, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxyzz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxxzzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyyy, g_z_0_yyz_xxxyyyz, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyyzz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxyzzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxxzzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyyy, g_z_0_yyz_xxyyyyz, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyyzz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyyzzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxyzzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xxzzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyyy, g_z_0_yyz_xyyyyyz, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyyzz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyyzzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyyzzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xyzzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_xzzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_xxxxxx[k] = -g_z_0_yyz_xxxxxx[k] * ab_x + g_z_0_yyz_xxxxxxx[k];

                g_z_0_xyyz_xxxxxy[k] = -g_z_0_yyz_xxxxxy[k] * ab_x + g_z_0_yyz_xxxxxxy[k];

                g_z_0_xyyz_xxxxxz[k] = -g_z_0_yyz_xxxxxz[k] * ab_x + g_z_0_yyz_xxxxxxz[k];

                g_z_0_xyyz_xxxxyy[k] = -g_z_0_yyz_xxxxyy[k] * ab_x + g_z_0_yyz_xxxxxyy[k];

                g_z_0_xyyz_xxxxyz[k] = -g_z_0_yyz_xxxxyz[k] * ab_x + g_z_0_yyz_xxxxxyz[k];

                g_z_0_xyyz_xxxxzz[k] = -g_z_0_yyz_xxxxzz[k] * ab_x + g_z_0_yyz_xxxxxzz[k];

                g_z_0_xyyz_xxxyyy[k] = -g_z_0_yyz_xxxyyy[k] * ab_x + g_z_0_yyz_xxxxyyy[k];

                g_z_0_xyyz_xxxyyz[k] = -g_z_0_yyz_xxxyyz[k] * ab_x + g_z_0_yyz_xxxxyyz[k];

                g_z_0_xyyz_xxxyzz[k] = -g_z_0_yyz_xxxyzz[k] * ab_x + g_z_0_yyz_xxxxyzz[k];

                g_z_0_xyyz_xxxzzz[k] = -g_z_0_yyz_xxxzzz[k] * ab_x + g_z_0_yyz_xxxxzzz[k];

                g_z_0_xyyz_xxyyyy[k] = -g_z_0_yyz_xxyyyy[k] * ab_x + g_z_0_yyz_xxxyyyy[k];

                g_z_0_xyyz_xxyyyz[k] = -g_z_0_yyz_xxyyyz[k] * ab_x + g_z_0_yyz_xxxyyyz[k];

                g_z_0_xyyz_xxyyzz[k] = -g_z_0_yyz_xxyyzz[k] * ab_x + g_z_0_yyz_xxxyyzz[k];

                g_z_0_xyyz_xxyzzz[k] = -g_z_0_yyz_xxyzzz[k] * ab_x + g_z_0_yyz_xxxyzzz[k];

                g_z_0_xyyz_xxzzzz[k] = -g_z_0_yyz_xxzzzz[k] * ab_x + g_z_0_yyz_xxxzzzz[k];

                g_z_0_xyyz_xyyyyy[k] = -g_z_0_yyz_xyyyyy[k] * ab_x + g_z_0_yyz_xxyyyyy[k];

                g_z_0_xyyz_xyyyyz[k] = -g_z_0_yyz_xyyyyz[k] * ab_x + g_z_0_yyz_xxyyyyz[k];

                g_z_0_xyyz_xyyyzz[k] = -g_z_0_yyz_xyyyzz[k] * ab_x + g_z_0_yyz_xxyyyzz[k];

                g_z_0_xyyz_xyyzzz[k] = -g_z_0_yyz_xyyzzz[k] * ab_x + g_z_0_yyz_xxyyzzz[k];

                g_z_0_xyyz_xyzzzz[k] = -g_z_0_yyz_xyzzzz[k] * ab_x + g_z_0_yyz_xxyzzzz[k];

                g_z_0_xyyz_xzzzzz[k] = -g_z_0_yyz_xzzzzz[k] * ab_x + g_z_0_yyz_xxzzzzz[k];

                g_z_0_xyyz_yyyyyy[k] = -g_z_0_yyz_yyyyyy[k] * ab_x + g_z_0_yyz_xyyyyyy[k];

                g_z_0_xyyz_yyyyyz[k] = -g_z_0_yyz_yyyyyz[k] * ab_x + g_z_0_yyz_xyyyyyz[k];

                g_z_0_xyyz_yyyyzz[k] = -g_z_0_yyz_yyyyzz[k] * ab_x + g_z_0_yyz_xyyyyzz[k];

                g_z_0_xyyz_yyyzzz[k] = -g_z_0_yyz_yyyzzz[k] * ab_x + g_z_0_yyz_xyyyzzz[k];

                g_z_0_xyyz_yyzzzz[k] = -g_z_0_yyz_yyzzzz[k] * ab_x + g_z_0_yyz_xyyzzzz[k];

                g_z_0_xyyz_yzzzzz[k] = -g_z_0_yyz_yzzzzz[k] * ab_x + g_z_0_yyz_xyzzzzz[k];

                g_z_0_xyyz_zzzzzz[k] = -g_z_0_yyz_zzzzzz[k] * ab_x + g_z_0_yyz_xzzzzzz[k];
            }

            /// Set up 1064-1092 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1064 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1079 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_yyyyyy, g_z_0_xyzz_yyyyyz, g_z_0_xyzz_yyyyzz, g_z_0_xyzz_yyyzzz, g_z_0_xyzz_yyzzzz, g_z_0_xyzz_yzzzzz, g_z_0_xyzz_zzzzzz, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxxx, g_z_0_yzz_xxxxxxy, g_z_0_yzz_xxxxxxz, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxyy, g_z_0_yzz_xxxxxyz, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxxzz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyyy, g_z_0_yzz_xxxxyyz, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxyzz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxxzzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyyy, g_z_0_yzz_xxxyyyz, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyyzz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxyzzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxxzzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyyy, g_z_0_yzz_xxyyyyz, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyyzz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyyzzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxyzzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xxzzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyyy, g_z_0_yzz_xyyyyyz, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyyzz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyyzzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyyzzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xyzzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_xzzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_xxxxxx[k] = -g_z_0_yzz_xxxxxx[k] * ab_x + g_z_0_yzz_xxxxxxx[k];

                g_z_0_xyzz_xxxxxy[k] = -g_z_0_yzz_xxxxxy[k] * ab_x + g_z_0_yzz_xxxxxxy[k];

                g_z_0_xyzz_xxxxxz[k] = -g_z_0_yzz_xxxxxz[k] * ab_x + g_z_0_yzz_xxxxxxz[k];

                g_z_0_xyzz_xxxxyy[k] = -g_z_0_yzz_xxxxyy[k] * ab_x + g_z_0_yzz_xxxxxyy[k];

                g_z_0_xyzz_xxxxyz[k] = -g_z_0_yzz_xxxxyz[k] * ab_x + g_z_0_yzz_xxxxxyz[k];

                g_z_0_xyzz_xxxxzz[k] = -g_z_0_yzz_xxxxzz[k] * ab_x + g_z_0_yzz_xxxxxzz[k];

                g_z_0_xyzz_xxxyyy[k] = -g_z_0_yzz_xxxyyy[k] * ab_x + g_z_0_yzz_xxxxyyy[k];

                g_z_0_xyzz_xxxyyz[k] = -g_z_0_yzz_xxxyyz[k] * ab_x + g_z_0_yzz_xxxxyyz[k];

                g_z_0_xyzz_xxxyzz[k] = -g_z_0_yzz_xxxyzz[k] * ab_x + g_z_0_yzz_xxxxyzz[k];

                g_z_0_xyzz_xxxzzz[k] = -g_z_0_yzz_xxxzzz[k] * ab_x + g_z_0_yzz_xxxxzzz[k];

                g_z_0_xyzz_xxyyyy[k] = -g_z_0_yzz_xxyyyy[k] * ab_x + g_z_0_yzz_xxxyyyy[k];

                g_z_0_xyzz_xxyyyz[k] = -g_z_0_yzz_xxyyyz[k] * ab_x + g_z_0_yzz_xxxyyyz[k];

                g_z_0_xyzz_xxyyzz[k] = -g_z_0_yzz_xxyyzz[k] * ab_x + g_z_0_yzz_xxxyyzz[k];

                g_z_0_xyzz_xxyzzz[k] = -g_z_0_yzz_xxyzzz[k] * ab_x + g_z_0_yzz_xxxyzzz[k];

                g_z_0_xyzz_xxzzzz[k] = -g_z_0_yzz_xxzzzz[k] * ab_x + g_z_0_yzz_xxxzzzz[k];

                g_z_0_xyzz_xyyyyy[k] = -g_z_0_yzz_xyyyyy[k] * ab_x + g_z_0_yzz_xxyyyyy[k];

                g_z_0_xyzz_xyyyyz[k] = -g_z_0_yzz_xyyyyz[k] * ab_x + g_z_0_yzz_xxyyyyz[k];

                g_z_0_xyzz_xyyyzz[k] = -g_z_0_yzz_xyyyzz[k] * ab_x + g_z_0_yzz_xxyyyzz[k];

                g_z_0_xyzz_xyyzzz[k] = -g_z_0_yzz_xyyzzz[k] * ab_x + g_z_0_yzz_xxyyzzz[k];

                g_z_0_xyzz_xyzzzz[k] = -g_z_0_yzz_xyzzzz[k] * ab_x + g_z_0_yzz_xxyzzzz[k];

                g_z_0_xyzz_xzzzzz[k] = -g_z_0_yzz_xzzzzz[k] * ab_x + g_z_0_yzz_xxzzzzz[k];

                g_z_0_xyzz_yyyyyy[k] = -g_z_0_yzz_yyyyyy[k] * ab_x + g_z_0_yzz_xyyyyyy[k];

                g_z_0_xyzz_yyyyyz[k] = -g_z_0_yzz_yyyyyz[k] * ab_x + g_z_0_yzz_xyyyyyz[k];

                g_z_0_xyzz_yyyyzz[k] = -g_z_0_yzz_yyyyzz[k] * ab_x + g_z_0_yzz_xyyyyzz[k];

                g_z_0_xyzz_yyyzzz[k] = -g_z_0_yzz_yyyzzz[k] * ab_x + g_z_0_yzz_xyyyzzz[k];

                g_z_0_xyzz_yyzzzz[k] = -g_z_0_yzz_yyzzzz[k] * ab_x + g_z_0_yzz_xyyzzzz[k];

                g_z_0_xyzz_yzzzzz[k] = -g_z_0_yzz_yzzzzz[k] * ab_x + g_z_0_yzz_xyzzzzz[k];

                g_z_0_xyzz_zzzzzz[k] = -g_z_0_yzz_zzzzzz[k] * ab_x + g_z_0_yzz_xzzzzzz[k];
            }

            /// Set up 1092-1120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1094 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1108 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1109 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1110 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1111 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1112 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1113 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1114 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1115 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_yyyyyy, g_z_0_xzzz_yyyyyz, g_z_0_xzzz_yyyyzz, g_z_0_xzzz_yyyzzz, g_z_0_xzzz_yyzzzz, g_z_0_xzzz_yzzzzz, g_z_0_xzzz_zzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxx, g_z_0_zzz_xxxxxxy, g_z_0_zzz_xxxxxxz, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxxzz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxxzzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxxzzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xxzzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_xzzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_xxxxxx[k] = -g_z_0_zzz_xxxxxx[k] * ab_x + g_z_0_zzz_xxxxxxx[k];

                g_z_0_xzzz_xxxxxy[k] = -g_z_0_zzz_xxxxxy[k] * ab_x + g_z_0_zzz_xxxxxxy[k];

                g_z_0_xzzz_xxxxxz[k] = -g_z_0_zzz_xxxxxz[k] * ab_x + g_z_0_zzz_xxxxxxz[k];

                g_z_0_xzzz_xxxxyy[k] = -g_z_0_zzz_xxxxyy[k] * ab_x + g_z_0_zzz_xxxxxyy[k];

                g_z_0_xzzz_xxxxyz[k] = -g_z_0_zzz_xxxxyz[k] * ab_x + g_z_0_zzz_xxxxxyz[k];

                g_z_0_xzzz_xxxxzz[k] = -g_z_0_zzz_xxxxzz[k] * ab_x + g_z_0_zzz_xxxxxzz[k];

                g_z_0_xzzz_xxxyyy[k] = -g_z_0_zzz_xxxyyy[k] * ab_x + g_z_0_zzz_xxxxyyy[k];

                g_z_0_xzzz_xxxyyz[k] = -g_z_0_zzz_xxxyyz[k] * ab_x + g_z_0_zzz_xxxxyyz[k];

                g_z_0_xzzz_xxxyzz[k] = -g_z_0_zzz_xxxyzz[k] * ab_x + g_z_0_zzz_xxxxyzz[k];

                g_z_0_xzzz_xxxzzz[k] = -g_z_0_zzz_xxxzzz[k] * ab_x + g_z_0_zzz_xxxxzzz[k];

                g_z_0_xzzz_xxyyyy[k] = -g_z_0_zzz_xxyyyy[k] * ab_x + g_z_0_zzz_xxxyyyy[k];

                g_z_0_xzzz_xxyyyz[k] = -g_z_0_zzz_xxyyyz[k] * ab_x + g_z_0_zzz_xxxyyyz[k];

                g_z_0_xzzz_xxyyzz[k] = -g_z_0_zzz_xxyyzz[k] * ab_x + g_z_0_zzz_xxxyyzz[k];

                g_z_0_xzzz_xxyzzz[k] = -g_z_0_zzz_xxyzzz[k] * ab_x + g_z_0_zzz_xxxyzzz[k];

                g_z_0_xzzz_xxzzzz[k] = -g_z_0_zzz_xxzzzz[k] * ab_x + g_z_0_zzz_xxxzzzz[k];

                g_z_0_xzzz_xyyyyy[k] = -g_z_0_zzz_xyyyyy[k] * ab_x + g_z_0_zzz_xxyyyyy[k];

                g_z_0_xzzz_xyyyyz[k] = -g_z_0_zzz_xyyyyz[k] * ab_x + g_z_0_zzz_xxyyyyz[k];

                g_z_0_xzzz_xyyyzz[k] = -g_z_0_zzz_xyyyzz[k] * ab_x + g_z_0_zzz_xxyyyzz[k];

                g_z_0_xzzz_xyyzzz[k] = -g_z_0_zzz_xyyzzz[k] * ab_x + g_z_0_zzz_xxyyzzz[k];

                g_z_0_xzzz_xyzzzz[k] = -g_z_0_zzz_xyzzzz[k] * ab_x + g_z_0_zzz_xxyzzzz[k];

                g_z_0_xzzz_xzzzzz[k] = -g_z_0_zzz_xzzzzz[k] * ab_x + g_z_0_zzz_xxzzzzz[k];

                g_z_0_xzzz_yyyyyy[k] = -g_z_0_zzz_yyyyyy[k] * ab_x + g_z_0_zzz_xyyyyyy[k];

                g_z_0_xzzz_yyyyyz[k] = -g_z_0_zzz_yyyyyz[k] * ab_x + g_z_0_zzz_xyyyyyz[k];

                g_z_0_xzzz_yyyyzz[k] = -g_z_0_zzz_yyyyzz[k] * ab_x + g_z_0_zzz_xyyyyzz[k];

                g_z_0_xzzz_yyyzzz[k] = -g_z_0_zzz_yyyzzz[k] * ab_x + g_z_0_zzz_xyyyzzz[k];

                g_z_0_xzzz_yyzzzz[k] = -g_z_0_zzz_yyzzzz[k] * ab_x + g_z_0_zzz_xyyzzzz[k];

                g_z_0_xzzz_yzzzzz[k] = -g_z_0_zzz_yzzzzz[k] * ab_x + g_z_0_zzz_xyzzzzz[k];

                g_z_0_xzzz_zzzzzz[k] = -g_z_0_zzz_zzzzzz[k] * ab_x + g_z_0_zzz_xzzzzzz[k];
            }

            /// Set up 1120-1148 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 1124 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 1139 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 1144 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 1145 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 1146 * ccomps * dcomps);

            auto g_z_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 1147 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxxy, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxyy, g_z_0_yyy_xxxxxyz, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyyy, g_z_0_yyy_xxxxyyz, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxyzz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyyy, g_z_0_yyy_xxxyyyz, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyyzz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxyzzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyyy, g_z_0_yyy_xxyyyyz, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyyzz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyyzzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxyzzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyyy, g_z_0_yyy_xyyyyyz, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyyzz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyyzzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyyzzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xyzzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyyy, g_z_0_yyy_yyyyyyz, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyyzz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyyzzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyyzzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yyzzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_yzzzzzz, g_z_0_yyy_zzzzzz, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_xxxxxx[k] = -g_z_0_yyy_xxxxxx[k] * ab_y + g_z_0_yyy_xxxxxxy[k];

                g_z_0_yyyy_xxxxxy[k] = -g_z_0_yyy_xxxxxy[k] * ab_y + g_z_0_yyy_xxxxxyy[k];

                g_z_0_yyyy_xxxxxz[k] = -g_z_0_yyy_xxxxxz[k] * ab_y + g_z_0_yyy_xxxxxyz[k];

                g_z_0_yyyy_xxxxyy[k] = -g_z_0_yyy_xxxxyy[k] * ab_y + g_z_0_yyy_xxxxyyy[k];

                g_z_0_yyyy_xxxxyz[k] = -g_z_0_yyy_xxxxyz[k] * ab_y + g_z_0_yyy_xxxxyyz[k];

                g_z_0_yyyy_xxxxzz[k] = -g_z_0_yyy_xxxxzz[k] * ab_y + g_z_0_yyy_xxxxyzz[k];

                g_z_0_yyyy_xxxyyy[k] = -g_z_0_yyy_xxxyyy[k] * ab_y + g_z_0_yyy_xxxyyyy[k];

                g_z_0_yyyy_xxxyyz[k] = -g_z_0_yyy_xxxyyz[k] * ab_y + g_z_0_yyy_xxxyyyz[k];

                g_z_0_yyyy_xxxyzz[k] = -g_z_0_yyy_xxxyzz[k] * ab_y + g_z_0_yyy_xxxyyzz[k];

                g_z_0_yyyy_xxxzzz[k] = -g_z_0_yyy_xxxzzz[k] * ab_y + g_z_0_yyy_xxxyzzz[k];

                g_z_0_yyyy_xxyyyy[k] = -g_z_0_yyy_xxyyyy[k] * ab_y + g_z_0_yyy_xxyyyyy[k];

                g_z_0_yyyy_xxyyyz[k] = -g_z_0_yyy_xxyyyz[k] * ab_y + g_z_0_yyy_xxyyyyz[k];

                g_z_0_yyyy_xxyyzz[k] = -g_z_0_yyy_xxyyzz[k] * ab_y + g_z_0_yyy_xxyyyzz[k];

                g_z_0_yyyy_xxyzzz[k] = -g_z_0_yyy_xxyzzz[k] * ab_y + g_z_0_yyy_xxyyzzz[k];

                g_z_0_yyyy_xxzzzz[k] = -g_z_0_yyy_xxzzzz[k] * ab_y + g_z_0_yyy_xxyzzzz[k];

                g_z_0_yyyy_xyyyyy[k] = -g_z_0_yyy_xyyyyy[k] * ab_y + g_z_0_yyy_xyyyyyy[k];

                g_z_0_yyyy_xyyyyz[k] = -g_z_0_yyy_xyyyyz[k] * ab_y + g_z_0_yyy_xyyyyyz[k];

                g_z_0_yyyy_xyyyzz[k] = -g_z_0_yyy_xyyyzz[k] * ab_y + g_z_0_yyy_xyyyyzz[k];

                g_z_0_yyyy_xyyzzz[k] = -g_z_0_yyy_xyyzzz[k] * ab_y + g_z_0_yyy_xyyyzzz[k];

                g_z_0_yyyy_xyzzzz[k] = -g_z_0_yyy_xyzzzz[k] * ab_y + g_z_0_yyy_xyyzzzz[k];

                g_z_0_yyyy_xzzzzz[k] = -g_z_0_yyy_xzzzzz[k] * ab_y + g_z_0_yyy_xyzzzzz[k];

                g_z_0_yyyy_yyyyyy[k] = -g_z_0_yyy_yyyyyy[k] * ab_y + g_z_0_yyy_yyyyyyy[k];

                g_z_0_yyyy_yyyyyz[k] = -g_z_0_yyy_yyyyyz[k] * ab_y + g_z_0_yyy_yyyyyyz[k];

                g_z_0_yyyy_yyyyzz[k] = -g_z_0_yyy_yyyyzz[k] * ab_y + g_z_0_yyy_yyyyyzz[k];

                g_z_0_yyyy_yyyzzz[k] = -g_z_0_yyy_yyyzzz[k] * ab_y + g_z_0_yyy_yyyyzzz[k];

                g_z_0_yyyy_yyzzzz[k] = -g_z_0_yyy_yyzzzz[k] * ab_y + g_z_0_yyy_yyyzzzz[k];

                g_z_0_yyyy_yzzzzz[k] = -g_z_0_yyy_yzzzzz[k] * ab_y + g_z_0_yyy_yyzzzzz[k];

                g_z_0_yyyy_zzzzzz[k] = -g_z_0_yyy_zzzzzz[k] * ab_y + g_z_0_yyy_yzzzzzz[k];
            }

            /// Set up 1148-1176 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 1148 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 1149 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 1150 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 1151 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 1154 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 1169 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzzz, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxxy, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxyy, g_z_0_yyz_xxxxxyz, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyyy, g_z_0_yyz_xxxxyyz, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxyzz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyyy, g_z_0_yyz_xxxyyyz, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyyzz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxyzzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyyy, g_z_0_yyz_xxyyyyz, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyyzz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyyzzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxyzzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyyy, g_z_0_yyz_xyyyyyz, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyyzz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyyzzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyyzzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xyzzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyyy, g_z_0_yyz_yyyyyyz, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyyzz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyyzzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyyzzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yyzzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_yzzzzzz, g_z_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_xxxxxx[k] = -g_z_0_yyz_xxxxxx[k] * ab_y + g_z_0_yyz_xxxxxxy[k];

                g_z_0_yyyz_xxxxxy[k] = -g_z_0_yyz_xxxxxy[k] * ab_y + g_z_0_yyz_xxxxxyy[k];

                g_z_0_yyyz_xxxxxz[k] = -g_z_0_yyz_xxxxxz[k] * ab_y + g_z_0_yyz_xxxxxyz[k];

                g_z_0_yyyz_xxxxyy[k] = -g_z_0_yyz_xxxxyy[k] * ab_y + g_z_0_yyz_xxxxyyy[k];

                g_z_0_yyyz_xxxxyz[k] = -g_z_0_yyz_xxxxyz[k] * ab_y + g_z_0_yyz_xxxxyyz[k];

                g_z_0_yyyz_xxxxzz[k] = -g_z_0_yyz_xxxxzz[k] * ab_y + g_z_0_yyz_xxxxyzz[k];

                g_z_0_yyyz_xxxyyy[k] = -g_z_0_yyz_xxxyyy[k] * ab_y + g_z_0_yyz_xxxyyyy[k];

                g_z_0_yyyz_xxxyyz[k] = -g_z_0_yyz_xxxyyz[k] * ab_y + g_z_0_yyz_xxxyyyz[k];

                g_z_0_yyyz_xxxyzz[k] = -g_z_0_yyz_xxxyzz[k] * ab_y + g_z_0_yyz_xxxyyzz[k];

                g_z_0_yyyz_xxxzzz[k] = -g_z_0_yyz_xxxzzz[k] * ab_y + g_z_0_yyz_xxxyzzz[k];

                g_z_0_yyyz_xxyyyy[k] = -g_z_0_yyz_xxyyyy[k] * ab_y + g_z_0_yyz_xxyyyyy[k];

                g_z_0_yyyz_xxyyyz[k] = -g_z_0_yyz_xxyyyz[k] * ab_y + g_z_0_yyz_xxyyyyz[k];

                g_z_0_yyyz_xxyyzz[k] = -g_z_0_yyz_xxyyzz[k] * ab_y + g_z_0_yyz_xxyyyzz[k];

                g_z_0_yyyz_xxyzzz[k] = -g_z_0_yyz_xxyzzz[k] * ab_y + g_z_0_yyz_xxyyzzz[k];

                g_z_0_yyyz_xxzzzz[k] = -g_z_0_yyz_xxzzzz[k] * ab_y + g_z_0_yyz_xxyzzzz[k];

                g_z_0_yyyz_xyyyyy[k] = -g_z_0_yyz_xyyyyy[k] * ab_y + g_z_0_yyz_xyyyyyy[k];

                g_z_0_yyyz_xyyyyz[k] = -g_z_0_yyz_xyyyyz[k] * ab_y + g_z_0_yyz_xyyyyyz[k];

                g_z_0_yyyz_xyyyzz[k] = -g_z_0_yyz_xyyyzz[k] * ab_y + g_z_0_yyz_xyyyyzz[k];

                g_z_0_yyyz_xyyzzz[k] = -g_z_0_yyz_xyyzzz[k] * ab_y + g_z_0_yyz_xyyyzzz[k];

                g_z_0_yyyz_xyzzzz[k] = -g_z_0_yyz_xyzzzz[k] * ab_y + g_z_0_yyz_xyyzzzz[k];

                g_z_0_yyyz_xzzzzz[k] = -g_z_0_yyz_xzzzzz[k] * ab_y + g_z_0_yyz_xyzzzzz[k];

                g_z_0_yyyz_yyyyyy[k] = -g_z_0_yyz_yyyyyy[k] * ab_y + g_z_0_yyz_yyyyyyy[k];

                g_z_0_yyyz_yyyyyz[k] = -g_z_0_yyz_yyyyyz[k] * ab_y + g_z_0_yyz_yyyyyyz[k];

                g_z_0_yyyz_yyyyzz[k] = -g_z_0_yyz_yyyyzz[k] * ab_y + g_z_0_yyz_yyyyyzz[k];

                g_z_0_yyyz_yyyzzz[k] = -g_z_0_yyz_yyyzzz[k] * ab_y + g_z_0_yyz_yyyyzzz[k];

                g_z_0_yyyz_yyzzzz[k] = -g_z_0_yyz_yyzzzz[k] * ab_y + g_z_0_yyz_yyyzzzz[k];

                g_z_0_yyyz_yzzzzz[k] = -g_z_0_yyz_yzzzzz[k] * ab_y + g_z_0_yyz_yyzzzzz[k];

                g_z_0_yyyz_zzzzzz[k] = -g_z_0_yyz_zzzzzz[k] * ab_y + g_z_0_yyz_yzzzzzz[k];
            }

            /// Set up 1176-1204 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1203 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzzz, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxxy, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxyy, g_z_0_yzz_xxxxxyz, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyyy, g_z_0_yzz_xxxxyyz, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxyzz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyyy, g_z_0_yzz_xxxyyyz, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyyzz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxyzzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyyy, g_z_0_yzz_xxyyyyz, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyyzz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyyzzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxyzzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyyy, g_z_0_yzz_xyyyyyz, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyyzz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyyzzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyyzzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xyzzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyyy, g_z_0_yzz_yyyyyyz, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyyzz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyyzzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyyzzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yyzzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_yzzzzzz, g_z_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_xxxxxx[k] = -g_z_0_yzz_xxxxxx[k] * ab_y + g_z_0_yzz_xxxxxxy[k];

                g_z_0_yyzz_xxxxxy[k] = -g_z_0_yzz_xxxxxy[k] * ab_y + g_z_0_yzz_xxxxxyy[k];

                g_z_0_yyzz_xxxxxz[k] = -g_z_0_yzz_xxxxxz[k] * ab_y + g_z_0_yzz_xxxxxyz[k];

                g_z_0_yyzz_xxxxyy[k] = -g_z_0_yzz_xxxxyy[k] * ab_y + g_z_0_yzz_xxxxyyy[k];

                g_z_0_yyzz_xxxxyz[k] = -g_z_0_yzz_xxxxyz[k] * ab_y + g_z_0_yzz_xxxxyyz[k];

                g_z_0_yyzz_xxxxzz[k] = -g_z_0_yzz_xxxxzz[k] * ab_y + g_z_0_yzz_xxxxyzz[k];

                g_z_0_yyzz_xxxyyy[k] = -g_z_0_yzz_xxxyyy[k] * ab_y + g_z_0_yzz_xxxyyyy[k];

                g_z_0_yyzz_xxxyyz[k] = -g_z_0_yzz_xxxyyz[k] * ab_y + g_z_0_yzz_xxxyyyz[k];

                g_z_0_yyzz_xxxyzz[k] = -g_z_0_yzz_xxxyzz[k] * ab_y + g_z_0_yzz_xxxyyzz[k];

                g_z_0_yyzz_xxxzzz[k] = -g_z_0_yzz_xxxzzz[k] * ab_y + g_z_0_yzz_xxxyzzz[k];

                g_z_0_yyzz_xxyyyy[k] = -g_z_0_yzz_xxyyyy[k] * ab_y + g_z_0_yzz_xxyyyyy[k];

                g_z_0_yyzz_xxyyyz[k] = -g_z_0_yzz_xxyyyz[k] * ab_y + g_z_0_yzz_xxyyyyz[k];

                g_z_0_yyzz_xxyyzz[k] = -g_z_0_yzz_xxyyzz[k] * ab_y + g_z_0_yzz_xxyyyzz[k];

                g_z_0_yyzz_xxyzzz[k] = -g_z_0_yzz_xxyzzz[k] * ab_y + g_z_0_yzz_xxyyzzz[k];

                g_z_0_yyzz_xxzzzz[k] = -g_z_0_yzz_xxzzzz[k] * ab_y + g_z_0_yzz_xxyzzzz[k];

                g_z_0_yyzz_xyyyyy[k] = -g_z_0_yzz_xyyyyy[k] * ab_y + g_z_0_yzz_xyyyyyy[k];

                g_z_0_yyzz_xyyyyz[k] = -g_z_0_yzz_xyyyyz[k] * ab_y + g_z_0_yzz_xyyyyyz[k];

                g_z_0_yyzz_xyyyzz[k] = -g_z_0_yzz_xyyyzz[k] * ab_y + g_z_0_yzz_xyyyyzz[k];

                g_z_0_yyzz_xyyzzz[k] = -g_z_0_yzz_xyyzzz[k] * ab_y + g_z_0_yzz_xyyyzzz[k];

                g_z_0_yyzz_xyzzzz[k] = -g_z_0_yzz_xyzzzz[k] * ab_y + g_z_0_yzz_xyyzzzz[k];

                g_z_0_yyzz_xzzzzz[k] = -g_z_0_yzz_xzzzzz[k] * ab_y + g_z_0_yzz_xyzzzzz[k];

                g_z_0_yyzz_yyyyyy[k] = -g_z_0_yzz_yyyyyy[k] * ab_y + g_z_0_yzz_yyyyyyy[k];

                g_z_0_yyzz_yyyyyz[k] = -g_z_0_yzz_yyyyyz[k] * ab_y + g_z_0_yzz_yyyyyyz[k];

                g_z_0_yyzz_yyyyzz[k] = -g_z_0_yzz_yyyyzz[k] * ab_y + g_z_0_yzz_yyyyyzz[k];

                g_z_0_yyzz_yyyzzz[k] = -g_z_0_yzz_yyyzzz[k] * ab_y + g_z_0_yzz_yyyyzzz[k];

                g_z_0_yyzz_yyzzzz[k] = -g_z_0_yzz_yyzzzz[k] * ab_y + g_z_0_yzz_yyyzzzz[k];

                g_z_0_yyzz_yzzzzz[k] = -g_z_0_yzz_yzzzzz[k] * ab_y + g_z_0_yzz_yyzzzzz[k];

                g_z_0_yyzz_zzzzzz[k] = -g_z_0_yzz_zzzzzz[k] * ab_y + g_z_0_yzz_yzzzzzz[k];
            }

            /// Set up 1204-1232 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1231 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxy, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyyy, g_z_0_zzz_yyyyyyz, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyyzz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyyzzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyyzzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yyzzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_yzzzzzz, g_z_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_xxxxxx[k] = -g_z_0_zzz_xxxxxx[k] * ab_y + g_z_0_zzz_xxxxxxy[k];

                g_z_0_yzzz_xxxxxy[k] = -g_z_0_zzz_xxxxxy[k] * ab_y + g_z_0_zzz_xxxxxyy[k];

                g_z_0_yzzz_xxxxxz[k] = -g_z_0_zzz_xxxxxz[k] * ab_y + g_z_0_zzz_xxxxxyz[k];

                g_z_0_yzzz_xxxxyy[k] = -g_z_0_zzz_xxxxyy[k] * ab_y + g_z_0_zzz_xxxxyyy[k];

                g_z_0_yzzz_xxxxyz[k] = -g_z_0_zzz_xxxxyz[k] * ab_y + g_z_0_zzz_xxxxyyz[k];

                g_z_0_yzzz_xxxxzz[k] = -g_z_0_zzz_xxxxzz[k] * ab_y + g_z_0_zzz_xxxxyzz[k];

                g_z_0_yzzz_xxxyyy[k] = -g_z_0_zzz_xxxyyy[k] * ab_y + g_z_0_zzz_xxxyyyy[k];

                g_z_0_yzzz_xxxyyz[k] = -g_z_0_zzz_xxxyyz[k] * ab_y + g_z_0_zzz_xxxyyyz[k];

                g_z_0_yzzz_xxxyzz[k] = -g_z_0_zzz_xxxyzz[k] * ab_y + g_z_0_zzz_xxxyyzz[k];

                g_z_0_yzzz_xxxzzz[k] = -g_z_0_zzz_xxxzzz[k] * ab_y + g_z_0_zzz_xxxyzzz[k];

                g_z_0_yzzz_xxyyyy[k] = -g_z_0_zzz_xxyyyy[k] * ab_y + g_z_0_zzz_xxyyyyy[k];

                g_z_0_yzzz_xxyyyz[k] = -g_z_0_zzz_xxyyyz[k] * ab_y + g_z_0_zzz_xxyyyyz[k];

                g_z_0_yzzz_xxyyzz[k] = -g_z_0_zzz_xxyyzz[k] * ab_y + g_z_0_zzz_xxyyyzz[k];

                g_z_0_yzzz_xxyzzz[k] = -g_z_0_zzz_xxyzzz[k] * ab_y + g_z_0_zzz_xxyyzzz[k];

                g_z_0_yzzz_xxzzzz[k] = -g_z_0_zzz_xxzzzz[k] * ab_y + g_z_0_zzz_xxyzzzz[k];

                g_z_0_yzzz_xyyyyy[k] = -g_z_0_zzz_xyyyyy[k] * ab_y + g_z_0_zzz_xyyyyyy[k];

                g_z_0_yzzz_xyyyyz[k] = -g_z_0_zzz_xyyyyz[k] * ab_y + g_z_0_zzz_xyyyyyz[k];

                g_z_0_yzzz_xyyyzz[k] = -g_z_0_zzz_xyyyzz[k] * ab_y + g_z_0_zzz_xyyyyzz[k];

                g_z_0_yzzz_xyyzzz[k] = -g_z_0_zzz_xyyzzz[k] * ab_y + g_z_0_zzz_xyyyzzz[k];

                g_z_0_yzzz_xyzzzz[k] = -g_z_0_zzz_xyzzzz[k] * ab_y + g_z_0_zzz_xyyzzzz[k];

                g_z_0_yzzz_xzzzzz[k] = -g_z_0_zzz_xzzzzz[k] * ab_y + g_z_0_zzz_xyzzzzz[k];

                g_z_0_yzzz_yyyyyy[k] = -g_z_0_zzz_yyyyyy[k] * ab_y + g_z_0_zzz_yyyyyyy[k];

                g_z_0_yzzz_yyyyyz[k] = -g_z_0_zzz_yyyyyz[k] * ab_y + g_z_0_zzz_yyyyyyz[k];

                g_z_0_yzzz_yyyyzz[k] = -g_z_0_zzz_yyyyzz[k] * ab_y + g_z_0_zzz_yyyyyzz[k];

                g_z_0_yzzz_yyyzzz[k] = -g_z_0_zzz_yyyzzz[k] * ab_y + g_z_0_zzz_yyyyzzz[k];

                g_z_0_yzzz_yyzzzz[k] = -g_z_0_zzz_yyzzzz[k] * ab_y + g_z_0_zzz_yyyzzzz[k];

                g_z_0_yzzz_yzzzzz[k] = -g_z_0_zzz_yzzzzz[k] * ab_y + g_z_0_zzz_yyzzzzz[k];

                g_z_0_yzzz_zzzzzz[k] = -g_z_0_zzz_zzzzzz[k] * ab_y + g_z_0_zzz_yzzzzzz[k];
            }

            /// Set up 1232-1260 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxz, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxxzz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxxzzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxxzzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xxzzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_xzzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyyz, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyyzz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyyzzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyyzzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yyzzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_yzzzzzz, g_z_0_zzz_zzzzzz, g_z_0_zzz_zzzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzzz, g_zzz_xxxxxx, g_zzz_xxxxxy, g_zzz_xxxxxz, g_zzz_xxxxyy, g_zzz_xxxxyz, g_zzz_xxxxzz, g_zzz_xxxyyy, g_zzz_xxxyyz, g_zzz_xxxyzz, g_zzz_xxxzzz, g_zzz_xxyyyy, g_zzz_xxyyyz, g_zzz_xxyyzz, g_zzz_xxyzzz, g_zzz_xxzzzz, g_zzz_xyyyyy, g_zzz_xyyyyz, g_zzz_xyyyzz, g_zzz_xyyzzz, g_zzz_xyzzzz, g_zzz_xzzzzz, g_zzz_yyyyyy, g_zzz_yyyyyz, g_zzz_yyyyzz, g_zzz_yyyzzz, g_zzz_yyzzzz, g_zzz_yzzzzz, g_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_xxxxxx[k] = -g_zzz_xxxxxx[k] - g_z_0_zzz_xxxxxx[k] * ab_z + g_z_0_zzz_xxxxxxz[k];

                g_z_0_zzzz_xxxxxy[k] = -g_zzz_xxxxxy[k] - g_z_0_zzz_xxxxxy[k] * ab_z + g_z_0_zzz_xxxxxyz[k];

                g_z_0_zzzz_xxxxxz[k] = -g_zzz_xxxxxz[k] - g_z_0_zzz_xxxxxz[k] * ab_z + g_z_0_zzz_xxxxxzz[k];

                g_z_0_zzzz_xxxxyy[k] = -g_zzz_xxxxyy[k] - g_z_0_zzz_xxxxyy[k] * ab_z + g_z_0_zzz_xxxxyyz[k];

                g_z_0_zzzz_xxxxyz[k] = -g_zzz_xxxxyz[k] - g_z_0_zzz_xxxxyz[k] * ab_z + g_z_0_zzz_xxxxyzz[k];

                g_z_0_zzzz_xxxxzz[k] = -g_zzz_xxxxzz[k] - g_z_0_zzz_xxxxzz[k] * ab_z + g_z_0_zzz_xxxxzzz[k];

                g_z_0_zzzz_xxxyyy[k] = -g_zzz_xxxyyy[k] - g_z_0_zzz_xxxyyy[k] * ab_z + g_z_0_zzz_xxxyyyz[k];

                g_z_0_zzzz_xxxyyz[k] = -g_zzz_xxxyyz[k] - g_z_0_zzz_xxxyyz[k] * ab_z + g_z_0_zzz_xxxyyzz[k];

                g_z_0_zzzz_xxxyzz[k] = -g_zzz_xxxyzz[k] - g_z_0_zzz_xxxyzz[k] * ab_z + g_z_0_zzz_xxxyzzz[k];

                g_z_0_zzzz_xxxzzz[k] = -g_zzz_xxxzzz[k] - g_z_0_zzz_xxxzzz[k] * ab_z + g_z_0_zzz_xxxzzzz[k];

                g_z_0_zzzz_xxyyyy[k] = -g_zzz_xxyyyy[k] - g_z_0_zzz_xxyyyy[k] * ab_z + g_z_0_zzz_xxyyyyz[k];

                g_z_0_zzzz_xxyyyz[k] = -g_zzz_xxyyyz[k] - g_z_0_zzz_xxyyyz[k] * ab_z + g_z_0_zzz_xxyyyzz[k];

                g_z_0_zzzz_xxyyzz[k] = -g_zzz_xxyyzz[k] - g_z_0_zzz_xxyyzz[k] * ab_z + g_z_0_zzz_xxyyzzz[k];

                g_z_0_zzzz_xxyzzz[k] = -g_zzz_xxyzzz[k] - g_z_0_zzz_xxyzzz[k] * ab_z + g_z_0_zzz_xxyzzzz[k];

                g_z_0_zzzz_xxzzzz[k] = -g_zzz_xxzzzz[k] - g_z_0_zzz_xxzzzz[k] * ab_z + g_z_0_zzz_xxzzzzz[k];

                g_z_0_zzzz_xyyyyy[k] = -g_zzz_xyyyyy[k] - g_z_0_zzz_xyyyyy[k] * ab_z + g_z_0_zzz_xyyyyyz[k];

                g_z_0_zzzz_xyyyyz[k] = -g_zzz_xyyyyz[k] - g_z_0_zzz_xyyyyz[k] * ab_z + g_z_0_zzz_xyyyyzz[k];

                g_z_0_zzzz_xyyyzz[k] = -g_zzz_xyyyzz[k] - g_z_0_zzz_xyyyzz[k] * ab_z + g_z_0_zzz_xyyyzzz[k];

                g_z_0_zzzz_xyyzzz[k] = -g_zzz_xyyzzz[k] - g_z_0_zzz_xyyzzz[k] * ab_z + g_z_0_zzz_xyyzzzz[k];

                g_z_0_zzzz_xyzzzz[k] = -g_zzz_xyzzzz[k] - g_z_0_zzz_xyzzzz[k] * ab_z + g_z_0_zzz_xyzzzzz[k];

                g_z_0_zzzz_xzzzzz[k] = -g_zzz_xzzzzz[k] - g_z_0_zzz_xzzzzz[k] * ab_z + g_z_0_zzz_xzzzzzz[k];

                g_z_0_zzzz_yyyyyy[k] = -g_zzz_yyyyyy[k] - g_z_0_zzz_yyyyyy[k] * ab_z + g_z_0_zzz_yyyyyyz[k];

                g_z_0_zzzz_yyyyyz[k] = -g_zzz_yyyyyz[k] - g_z_0_zzz_yyyyyz[k] * ab_z + g_z_0_zzz_yyyyyzz[k];

                g_z_0_zzzz_yyyyzz[k] = -g_zzz_yyyyzz[k] - g_z_0_zzz_yyyyzz[k] * ab_z + g_z_0_zzz_yyyyzzz[k];

                g_z_0_zzzz_yyyzzz[k] = -g_zzz_yyyzzz[k] - g_z_0_zzz_yyyzzz[k] * ab_z + g_z_0_zzz_yyyzzzz[k];

                g_z_0_zzzz_yyzzzz[k] = -g_zzz_yyzzzz[k] - g_z_0_zzz_yyzzzz[k] * ab_z + g_z_0_zzz_yyzzzzz[k];

                g_z_0_zzzz_yzzzzz[k] = -g_zzz_yzzzzz[k] - g_z_0_zzz_yzzzzz[k] * ab_z + g_z_0_zzz_yzzzzzz[k];

                g_z_0_zzzz_zzzzzz[k] = -g_zzz_zzzzzz[k] - g_z_0_zzz_zzzzzz[k] * ab_z + g_z_0_zzz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

