#include "ElectronRepulsionGeom0100ContrRecKFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_kfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kfxx,
                                            const size_t idx_ifxx,
                                            const size_t idx_geom_01_ifxx,
                                            const size_t idx_geom_01_igxx,
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
            /// Set up components of auxilary buffer : IFSS

            const auto if_off = idx_ifxx + i * dcomps + j;

            auto g_xxxxxx_xxx = cbuffer.data(if_off + 0 * ccomps * dcomps);

            auto g_xxxxxx_xxy = cbuffer.data(if_off + 1 * ccomps * dcomps);

            auto g_xxxxxx_xxz = cbuffer.data(if_off + 2 * ccomps * dcomps);

            auto g_xxxxxx_xyy = cbuffer.data(if_off + 3 * ccomps * dcomps);

            auto g_xxxxxx_xyz = cbuffer.data(if_off + 4 * ccomps * dcomps);

            auto g_xxxxxx_xzz = cbuffer.data(if_off + 5 * ccomps * dcomps);

            auto g_xxxxxx_yyy = cbuffer.data(if_off + 6 * ccomps * dcomps);

            auto g_xxxxxx_yyz = cbuffer.data(if_off + 7 * ccomps * dcomps);

            auto g_xxxxxx_yzz = cbuffer.data(if_off + 8 * ccomps * dcomps);

            auto g_xxxxxx_zzz = cbuffer.data(if_off + 9 * ccomps * dcomps);

            auto g_xxxxxy_xxx = cbuffer.data(if_off + 10 * ccomps * dcomps);

            auto g_xxxxxy_xxy = cbuffer.data(if_off + 11 * ccomps * dcomps);

            auto g_xxxxxy_xxz = cbuffer.data(if_off + 12 * ccomps * dcomps);

            auto g_xxxxxy_xyy = cbuffer.data(if_off + 13 * ccomps * dcomps);

            auto g_xxxxxy_xyz = cbuffer.data(if_off + 14 * ccomps * dcomps);

            auto g_xxxxxy_xzz = cbuffer.data(if_off + 15 * ccomps * dcomps);

            auto g_xxxxxy_yyy = cbuffer.data(if_off + 16 * ccomps * dcomps);

            auto g_xxxxxy_yyz = cbuffer.data(if_off + 17 * ccomps * dcomps);

            auto g_xxxxxy_yzz = cbuffer.data(if_off + 18 * ccomps * dcomps);

            auto g_xxxxxy_zzz = cbuffer.data(if_off + 19 * ccomps * dcomps);

            auto g_xxxxxz_xxx = cbuffer.data(if_off + 20 * ccomps * dcomps);

            auto g_xxxxxz_xxy = cbuffer.data(if_off + 21 * ccomps * dcomps);

            auto g_xxxxxz_xxz = cbuffer.data(if_off + 22 * ccomps * dcomps);

            auto g_xxxxxz_xyy = cbuffer.data(if_off + 23 * ccomps * dcomps);

            auto g_xxxxxz_xyz = cbuffer.data(if_off + 24 * ccomps * dcomps);

            auto g_xxxxxz_xzz = cbuffer.data(if_off + 25 * ccomps * dcomps);

            auto g_xxxxxz_yyy = cbuffer.data(if_off + 26 * ccomps * dcomps);

            auto g_xxxxxz_yyz = cbuffer.data(if_off + 27 * ccomps * dcomps);

            auto g_xxxxxz_yzz = cbuffer.data(if_off + 28 * ccomps * dcomps);

            auto g_xxxxxz_zzz = cbuffer.data(if_off + 29 * ccomps * dcomps);

            auto g_xxxxyy_xxx = cbuffer.data(if_off + 30 * ccomps * dcomps);

            auto g_xxxxyy_xxy = cbuffer.data(if_off + 31 * ccomps * dcomps);

            auto g_xxxxyy_xxz = cbuffer.data(if_off + 32 * ccomps * dcomps);

            auto g_xxxxyy_xyy = cbuffer.data(if_off + 33 * ccomps * dcomps);

            auto g_xxxxyy_xyz = cbuffer.data(if_off + 34 * ccomps * dcomps);

            auto g_xxxxyy_xzz = cbuffer.data(if_off + 35 * ccomps * dcomps);

            auto g_xxxxyy_yyy = cbuffer.data(if_off + 36 * ccomps * dcomps);

            auto g_xxxxyy_yyz = cbuffer.data(if_off + 37 * ccomps * dcomps);

            auto g_xxxxyy_yzz = cbuffer.data(if_off + 38 * ccomps * dcomps);

            auto g_xxxxyy_zzz = cbuffer.data(if_off + 39 * ccomps * dcomps);

            auto g_xxxxyz_xxx = cbuffer.data(if_off + 40 * ccomps * dcomps);

            auto g_xxxxyz_xxy = cbuffer.data(if_off + 41 * ccomps * dcomps);

            auto g_xxxxyz_xxz = cbuffer.data(if_off + 42 * ccomps * dcomps);

            auto g_xxxxyz_xyy = cbuffer.data(if_off + 43 * ccomps * dcomps);

            auto g_xxxxyz_xyz = cbuffer.data(if_off + 44 * ccomps * dcomps);

            auto g_xxxxyz_xzz = cbuffer.data(if_off + 45 * ccomps * dcomps);

            auto g_xxxxyz_yyy = cbuffer.data(if_off + 46 * ccomps * dcomps);

            auto g_xxxxyz_yyz = cbuffer.data(if_off + 47 * ccomps * dcomps);

            auto g_xxxxyz_yzz = cbuffer.data(if_off + 48 * ccomps * dcomps);

            auto g_xxxxyz_zzz = cbuffer.data(if_off + 49 * ccomps * dcomps);

            auto g_xxxxzz_xxx = cbuffer.data(if_off + 50 * ccomps * dcomps);

            auto g_xxxxzz_xxy = cbuffer.data(if_off + 51 * ccomps * dcomps);

            auto g_xxxxzz_xxz = cbuffer.data(if_off + 52 * ccomps * dcomps);

            auto g_xxxxzz_xyy = cbuffer.data(if_off + 53 * ccomps * dcomps);

            auto g_xxxxzz_xyz = cbuffer.data(if_off + 54 * ccomps * dcomps);

            auto g_xxxxzz_xzz = cbuffer.data(if_off + 55 * ccomps * dcomps);

            auto g_xxxxzz_yyy = cbuffer.data(if_off + 56 * ccomps * dcomps);

            auto g_xxxxzz_yyz = cbuffer.data(if_off + 57 * ccomps * dcomps);

            auto g_xxxxzz_yzz = cbuffer.data(if_off + 58 * ccomps * dcomps);

            auto g_xxxxzz_zzz = cbuffer.data(if_off + 59 * ccomps * dcomps);

            auto g_xxxyyy_xxx = cbuffer.data(if_off + 60 * ccomps * dcomps);

            auto g_xxxyyy_xxy = cbuffer.data(if_off + 61 * ccomps * dcomps);

            auto g_xxxyyy_xxz = cbuffer.data(if_off + 62 * ccomps * dcomps);

            auto g_xxxyyy_xyy = cbuffer.data(if_off + 63 * ccomps * dcomps);

            auto g_xxxyyy_xyz = cbuffer.data(if_off + 64 * ccomps * dcomps);

            auto g_xxxyyy_xzz = cbuffer.data(if_off + 65 * ccomps * dcomps);

            auto g_xxxyyy_yyy = cbuffer.data(if_off + 66 * ccomps * dcomps);

            auto g_xxxyyy_yyz = cbuffer.data(if_off + 67 * ccomps * dcomps);

            auto g_xxxyyy_yzz = cbuffer.data(if_off + 68 * ccomps * dcomps);

            auto g_xxxyyy_zzz = cbuffer.data(if_off + 69 * ccomps * dcomps);

            auto g_xxxyyz_xxx = cbuffer.data(if_off + 70 * ccomps * dcomps);

            auto g_xxxyyz_xxy = cbuffer.data(if_off + 71 * ccomps * dcomps);

            auto g_xxxyyz_xxz = cbuffer.data(if_off + 72 * ccomps * dcomps);

            auto g_xxxyyz_xyy = cbuffer.data(if_off + 73 * ccomps * dcomps);

            auto g_xxxyyz_xyz = cbuffer.data(if_off + 74 * ccomps * dcomps);

            auto g_xxxyyz_xzz = cbuffer.data(if_off + 75 * ccomps * dcomps);

            auto g_xxxyyz_yyy = cbuffer.data(if_off + 76 * ccomps * dcomps);

            auto g_xxxyyz_yyz = cbuffer.data(if_off + 77 * ccomps * dcomps);

            auto g_xxxyyz_yzz = cbuffer.data(if_off + 78 * ccomps * dcomps);

            auto g_xxxyyz_zzz = cbuffer.data(if_off + 79 * ccomps * dcomps);

            auto g_xxxyzz_xxx = cbuffer.data(if_off + 80 * ccomps * dcomps);

            auto g_xxxyzz_xxy = cbuffer.data(if_off + 81 * ccomps * dcomps);

            auto g_xxxyzz_xxz = cbuffer.data(if_off + 82 * ccomps * dcomps);

            auto g_xxxyzz_xyy = cbuffer.data(if_off + 83 * ccomps * dcomps);

            auto g_xxxyzz_xyz = cbuffer.data(if_off + 84 * ccomps * dcomps);

            auto g_xxxyzz_xzz = cbuffer.data(if_off + 85 * ccomps * dcomps);

            auto g_xxxyzz_yyy = cbuffer.data(if_off + 86 * ccomps * dcomps);

            auto g_xxxyzz_yyz = cbuffer.data(if_off + 87 * ccomps * dcomps);

            auto g_xxxyzz_yzz = cbuffer.data(if_off + 88 * ccomps * dcomps);

            auto g_xxxyzz_zzz = cbuffer.data(if_off + 89 * ccomps * dcomps);

            auto g_xxxzzz_xxx = cbuffer.data(if_off + 90 * ccomps * dcomps);

            auto g_xxxzzz_xxy = cbuffer.data(if_off + 91 * ccomps * dcomps);

            auto g_xxxzzz_xxz = cbuffer.data(if_off + 92 * ccomps * dcomps);

            auto g_xxxzzz_xyy = cbuffer.data(if_off + 93 * ccomps * dcomps);

            auto g_xxxzzz_xyz = cbuffer.data(if_off + 94 * ccomps * dcomps);

            auto g_xxxzzz_xzz = cbuffer.data(if_off + 95 * ccomps * dcomps);

            auto g_xxxzzz_yyy = cbuffer.data(if_off + 96 * ccomps * dcomps);

            auto g_xxxzzz_yyz = cbuffer.data(if_off + 97 * ccomps * dcomps);

            auto g_xxxzzz_yzz = cbuffer.data(if_off + 98 * ccomps * dcomps);

            auto g_xxxzzz_zzz = cbuffer.data(if_off + 99 * ccomps * dcomps);

            auto g_xxyyyy_xxx = cbuffer.data(if_off + 100 * ccomps * dcomps);

            auto g_xxyyyy_xxy = cbuffer.data(if_off + 101 * ccomps * dcomps);

            auto g_xxyyyy_xxz = cbuffer.data(if_off + 102 * ccomps * dcomps);

            auto g_xxyyyy_xyy = cbuffer.data(if_off + 103 * ccomps * dcomps);

            auto g_xxyyyy_xyz = cbuffer.data(if_off + 104 * ccomps * dcomps);

            auto g_xxyyyy_xzz = cbuffer.data(if_off + 105 * ccomps * dcomps);

            auto g_xxyyyy_yyy = cbuffer.data(if_off + 106 * ccomps * dcomps);

            auto g_xxyyyy_yyz = cbuffer.data(if_off + 107 * ccomps * dcomps);

            auto g_xxyyyy_yzz = cbuffer.data(if_off + 108 * ccomps * dcomps);

            auto g_xxyyyy_zzz = cbuffer.data(if_off + 109 * ccomps * dcomps);

            auto g_xxyyyz_xxx = cbuffer.data(if_off + 110 * ccomps * dcomps);

            auto g_xxyyyz_xxy = cbuffer.data(if_off + 111 * ccomps * dcomps);

            auto g_xxyyyz_xxz = cbuffer.data(if_off + 112 * ccomps * dcomps);

            auto g_xxyyyz_xyy = cbuffer.data(if_off + 113 * ccomps * dcomps);

            auto g_xxyyyz_xyz = cbuffer.data(if_off + 114 * ccomps * dcomps);

            auto g_xxyyyz_xzz = cbuffer.data(if_off + 115 * ccomps * dcomps);

            auto g_xxyyyz_yyy = cbuffer.data(if_off + 116 * ccomps * dcomps);

            auto g_xxyyyz_yyz = cbuffer.data(if_off + 117 * ccomps * dcomps);

            auto g_xxyyyz_yzz = cbuffer.data(if_off + 118 * ccomps * dcomps);

            auto g_xxyyyz_zzz = cbuffer.data(if_off + 119 * ccomps * dcomps);

            auto g_xxyyzz_xxx = cbuffer.data(if_off + 120 * ccomps * dcomps);

            auto g_xxyyzz_xxy = cbuffer.data(if_off + 121 * ccomps * dcomps);

            auto g_xxyyzz_xxz = cbuffer.data(if_off + 122 * ccomps * dcomps);

            auto g_xxyyzz_xyy = cbuffer.data(if_off + 123 * ccomps * dcomps);

            auto g_xxyyzz_xyz = cbuffer.data(if_off + 124 * ccomps * dcomps);

            auto g_xxyyzz_xzz = cbuffer.data(if_off + 125 * ccomps * dcomps);

            auto g_xxyyzz_yyy = cbuffer.data(if_off + 126 * ccomps * dcomps);

            auto g_xxyyzz_yyz = cbuffer.data(if_off + 127 * ccomps * dcomps);

            auto g_xxyyzz_yzz = cbuffer.data(if_off + 128 * ccomps * dcomps);

            auto g_xxyyzz_zzz = cbuffer.data(if_off + 129 * ccomps * dcomps);

            auto g_xxyzzz_xxx = cbuffer.data(if_off + 130 * ccomps * dcomps);

            auto g_xxyzzz_xxy = cbuffer.data(if_off + 131 * ccomps * dcomps);

            auto g_xxyzzz_xxz = cbuffer.data(if_off + 132 * ccomps * dcomps);

            auto g_xxyzzz_xyy = cbuffer.data(if_off + 133 * ccomps * dcomps);

            auto g_xxyzzz_xyz = cbuffer.data(if_off + 134 * ccomps * dcomps);

            auto g_xxyzzz_xzz = cbuffer.data(if_off + 135 * ccomps * dcomps);

            auto g_xxyzzz_yyy = cbuffer.data(if_off + 136 * ccomps * dcomps);

            auto g_xxyzzz_yyz = cbuffer.data(if_off + 137 * ccomps * dcomps);

            auto g_xxyzzz_yzz = cbuffer.data(if_off + 138 * ccomps * dcomps);

            auto g_xxyzzz_zzz = cbuffer.data(if_off + 139 * ccomps * dcomps);

            auto g_xxzzzz_xxx = cbuffer.data(if_off + 140 * ccomps * dcomps);

            auto g_xxzzzz_xxy = cbuffer.data(if_off + 141 * ccomps * dcomps);

            auto g_xxzzzz_xxz = cbuffer.data(if_off + 142 * ccomps * dcomps);

            auto g_xxzzzz_xyy = cbuffer.data(if_off + 143 * ccomps * dcomps);

            auto g_xxzzzz_xyz = cbuffer.data(if_off + 144 * ccomps * dcomps);

            auto g_xxzzzz_xzz = cbuffer.data(if_off + 145 * ccomps * dcomps);

            auto g_xxzzzz_yyy = cbuffer.data(if_off + 146 * ccomps * dcomps);

            auto g_xxzzzz_yyz = cbuffer.data(if_off + 147 * ccomps * dcomps);

            auto g_xxzzzz_yzz = cbuffer.data(if_off + 148 * ccomps * dcomps);

            auto g_xxzzzz_zzz = cbuffer.data(if_off + 149 * ccomps * dcomps);

            auto g_xyyyyy_xxx = cbuffer.data(if_off + 150 * ccomps * dcomps);

            auto g_xyyyyy_xxy = cbuffer.data(if_off + 151 * ccomps * dcomps);

            auto g_xyyyyy_xxz = cbuffer.data(if_off + 152 * ccomps * dcomps);

            auto g_xyyyyy_xyy = cbuffer.data(if_off + 153 * ccomps * dcomps);

            auto g_xyyyyy_xyz = cbuffer.data(if_off + 154 * ccomps * dcomps);

            auto g_xyyyyy_xzz = cbuffer.data(if_off + 155 * ccomps * dcomps);

            auto g_xyyyyy_yyy = cbuffer.data(if_off + 156 * ccomps * dcomps);

            auto g_xyyyyy_yyz = cbuffer.data(if_off + 157 * ccomps * dcomps);

            auto g_xyyyyy_yzz = cbuffer.data(if_off + 158 * ccomps * dcomps);

            auto g_xyyyyy_zzz = cbuffer.data(if_off + 159 * ccomps * dcomps);

            auto g_xyyyyz_xxx = cbuffer.data(if_off + 160 * ccomps * dcomps);

            auto g_xyyyyz_xxy = cbuffer.data(if_off + 161 * ccomps * dcomps);

            auto g_xyyyyz_xxz = cbuffer.data(if_off + 162 * ccomps * dcomps);

            auto g_xyyyyz_xyy = cbuffer.data(if_off + 163 * ccomps * dcomps);

            auto g_xyyyyz_xyz = cbuffer.data(if_off + 164 * ccomps * dcomps);

            auto g_xyyyyz_xzz = cbuffer.data(if_off + 165 * ccomps * dcomps);

            auto g_xyyyyz_yyy = cbuffer.data(if_off + 166 * ccomps * dcomps);

            auto g_xyyyyz_yyz = cbuffer.data(if_off + 167 * ccomps * dcomps);

            auto g_xyyyyz_yzz = cbuffer.data(if_off + 168 * ccomps * dcomps);

            auto g_xyyyyz_zzz = cbuffer.data(if_off + 169 * ccomps * dcomps);

            auto g_xyyyzz_xxx = cbuffer.data(if_off + 170 * ccomps * dcomps);

            auto g_xyyyzz_xxy = cbuffer.data(if_off + 171 * ccomps * dcomps);

            auto g_xyyyzz_xxz = cbuffer.data(if_off + 172 * ccomps * dcomps);

            auto g_xyyyzz_xyy = cbuffer.data(if_off + 173 * ccomps * dcomps);

            auto g_xyyyzz_xyz = cbuffer.data(if_off + 174 * ccomps * dcomps);

            auto g_xyyyzz_xzz = cbuffer.data(if_off + 175 * ccomps * dcomps);

            auto g_xyyyzz_yyy = cbuffer.data(if_off + 176 * ccomps * dcomps);

            auto g_xyyyzz_yyz = cbuffer.data(if_off + 177 * ccomps * dcomps);

            auto g_xyyyzz_yzz = cbuffer.data(if_off + 178 * ccomps * dcomps);

            auto g_xyyyzz_zzz = cbuffer.data(if_off + 179 * ccomps * dcomps);

            auto g_xyyzzz_xxx = cbuffer.data(if_off + 180 * ccomps * dcomps);

            auto g_xyyzzz_xxy = cbuffer.data(if_off + 181 * ccomps * dcomps);

            auto g_xyyzzz_xxz = cbuffer.data(if_off + 182 * ccomps * dcomps);

            auto g_xyyzzz_xyy = cbuffer.data(if_off + 183 * ccomps * dcomps);

            auto g_xyyzzz_xyz = cbuffer.data(if_off + 184 * ccomps * dcomps);

            auto g_xyyzzz_xzz = cbuffer.data(if_off + 185 * ccomps * dcomps);

            auto g_xyyzzz_yyy = cbuffer.data(if_off + 186 * ccomps * dcomps);

            auto g_xyyzzz_yyz = cbuffer.data(if_off + 187 * ccomps * dcomps);

            auto g_xyyzzz_yzz = cbuffer.data(if_off + 188 * ccomps * dcomps);

            auto g_xyyzzz_zzz = cbuffer.data(if_off + 189 * ccomps * dcomps);

            auto g_xyzzzz_xxx = cbuffer.data(if_off + 190 * ccomps * dcomps);

            auto g_xyzzzz_xxy = cbuffer.data(if_off + 191 * ccomps * dcomps);

            auto g_xyzzzz_xxz = cbuffer.data(if_off + 192 * ccomps * dcomps);

            auto g_xyzzzz_xyy = cbuffer.data(if_off + 193 * ccomps * dcomps);

            auto g_xyzzzz_xyz = cbuffer.data(if_off + 194 * ccomps * dcomps);

            auto g_xyzzzz_xzz = cbuffer.data(if_off + 195 * ccomps * dcomps);

            auto g_xyzzzz_yyy = cbuffer.data(if_off + 196 * ccomps * dcomps);

            auto g_xyzzzz_yyz = cbuffer.data(if_off + 197 * ccomps * dcomps);

            auto g_xyzzzz_yzz = cbuffer.data(if_off + 198 * ccomps * dcomps);

            auto g_xyzzzz_zzz = cbuffer.data(if_off + 199 * ccomps * dcomps);

            auto g_xzzzzz_xxx = cbuffer.data(if_off + 200 * ccomps * dcomps);

            auto g_xzzzzz_xxy = cbuffer.data(if_off + 201 * ccomps * dcomps);

            auto g_xzzzzz_xxz = cbuffer.data(if_off + 202 * ccomps * dcomps);

            auto g_xzzzzz_xyy = cbuffer.data(if_off + 203 * ccomps * dcomps);

            auto g_xzzzzz_xyz = cbuffer.data(if_off + 204 * ccomps * dcomps);

            auto g_xzzzzz_xzz = cbuffer.data(if_off + 205 * ccomps * dcomps);

            auto g_xzzzzz_yyy = cbuffer.data(if_off + 206 * ccomps * dcomps);

            auto g_xzzzzz_yyz = cbuffer.data(if_off + 207 * ccomps * dcomps);

            auto g_xzzzzz_yzz = cbuffer.data(if_off + 208 * ccomps * dcomps);

            auto g_xzzzzz_zzz = cbuffer.data(if_off + 209 * ccomps * dcomps);

            auto g_yyyyyy_xxx = cbuffer.data(if_off + 210 * ccomps * dcomps);

            auto g_yyyyyy_xxy = cbuffer.data(if_off + 211 * ccomps * dcomps);

            auto g_yyyyyy_xxz = cbuffer.data(if_off + 212 * ccomps * dcomps);

            auto g_yyyyyy_xyy = cbuffer.data(if_off + 213 * ccomps * dcomps);

            auto g_yyyyyy_xyz = cbuffer.data(if_off + 214 * ccomps * dcomps);

            auto g_yyyyyy_xzz = cbuffer.data(if_off + 215 * ccomps * dcomps);

            auto g_yyyyyy_yyy = cbuffer.data(if_off + 216 * ccomps * dcomps);

            auto g_yyyyyy_yyz = cbuffer.data(if_off + 217 * ccomps * dcomps);

            auto g_yyyyyy_yzz = cbuffer.data(if_off + 218 * ccomps * dcomps);

            auto g_yyyyyy_zzz = cbuffer.data(if_off + 219 * ccomps * dcomps);

            auto g_yyyyyz_xxx = cbuffer.data(if_off + 220 * ccomps * dcomps);

            auto g_yyyyyz_xxy = cbuffer.data(if_off + 221 * ccomps * dcomps);

            auto g_yyyyyz_xxz = cbuffer.data(if_off + 222 * ccomps * dcomps);

            auto g_yyyyyz_xyy = cbuffer.data(if_off + 223 * ccomps * dcomps);

            auto g_yyyyyz_xyz = cbuffer.data(if_off + 224 * ccomps * dcomps);

            auto g_yyyyyz_xzz = cbuffer.data(if_off + 225 * ccomps * dcomps);

            auto g_yyyyyz_yyy = cbuffer.data(if_off + 226 * ccomps * dcomps);

            auto g_yyyyyz_yyz = cbuffer.data(if_off + 227 * ccomps * dcomps);

            auto g_yyyyyz_yzz = cbuffer.data(if_off + 228 * ccomps * dcomps);

            auto g_yyyyyz_zzz = cbuffer.data(if_off + 229 * ccomps * dcomps);

            auto g_yyyyzz_xxx = cbuffer.data(if_off + 230 * ccomps * dcomps);

            auto g_yyyyzz_xxy = cbuffer.data(if_off + 231 * ccomps * dcomps);

            auto g_yyyyzz_xxz = cbuffer.data(if_off + 232 * ccomps * dcomps);

            auto g_yyyyzz_xyy = cbuffer.data(if_off + 233 * ccomps * dcomps);

            auto g_yyyyzz_xyz = cbuffer.data(if_off + 234 * ccomps * dcomps);

            auto g_yyyyzz_xzz = cbuffer.data(if_off + 235 * ccomps * dcomps);

            auto g_yyyyzz_yyy = cbuffer.data(if_off + 236 * ccomps * dcomps);

            auto g_yyyyzz_yyz = cbuffer.data(if_off + 237 * ccomps * dcomps);

            auto g_yyyyzz_yzz = cbuffer.data(if_off + 238 * ccomps * dcomps);

            auto g_yyyyzz_zzz = cbuffer.data(if_off + 239 * ccomps * dcomps);

            auto g_yyyzzz_xxx = cbuffer.data(if_off + 240 * ccomps * dcomps);

            auto g_yyyzzz_xxy = cbuffer.data(if_off + 241 * ccomps * dcomps);

            auto g_yyyzzz_xxz = cbuffer.data(if_off + 242 * ccomps * dcomps);

            auto g_yyyzzz_xyy = cbuffer.data(if_off + 243 * ccomps * dcomps);

            auto g_yyyzzz_xyz = cbuffer.data(if_off + 244 * ccomps * dcomps);

            auto g_yyyzzz_xzz = cbuffer.data(if_off + 245 * ccomps * dcomps);

            auto g_yyyzzz_yyy = cbuffer.data(if_off + 246 * ccomps * dcomps);

            auto g_yyyzzz_yyz = cbuffer.data(if_off + 247 * ccomps * dcomps);

            auto g_yyyzzz_yzz = cbuffer.data(if_off + 248 * ccomps * dcomps);

            auto g_yyyzzz_zzz = cbuffer.data(if_off + 249 * ccomps * dcomps);

            auto g_yyzzzz_xxx = cbuffer.data(if_off + 250 * ccomps * dcomps);

            auto g_yyzzzz_xxy = cbuffer.data(if_off + 251 * ccomps * dcomps);

            auto g_yyzzzz_xxz = cbuffer.data(if_off + 252 * ccomps * dcomps);

            auto g_yyzzzz_xyy = cbuffer.data(if_off + 253 * ccomps * dcomps);

            auto g_yyzzzz_xyz = cbuffer.data(if_off + 254 * ccomps * dcomps);

            auto g_yyzzzz_xzz = cbuffer.data(if_off + 255 * ccomps * dcomps);

            auto g_yyzzzz_yyy = cbuffer.data(if_off + 256 * ccomps * dcomps);

            auto g_yyzzzz_yyz = cbuffer.data(if_off + 257 * ccomps * dcomps);

            auto g_yyzzzz_yzz = cbuffer.data(if_off + 258 * ccomps * dcomps);

            auto g_yyzzzz_zzz = cbuffer.data(if_off + 259 * ccomps * dcomps);

            auto g_yzzzzz_xxx = cbuffer.data(if_off + 260 * ccomps * dcomps);

            auto g_yzzzzz_xxy = cbuffer.data(if_off + 261 * ccomps * dcomps);

            auto g_yzzzzz_xxz = cbuffer.data(if_off + 262 * ccomps * dcomps);

            auto g_yzzzzz_xyy = cbuffer.data(if_off + 263 * ccomps * dcomps);

            auto g_yzzzzz_xyz = cbuffer.data(if_off + 264 * ccomps * dcomps);

            auto g_yzzzzz_xzz = cbuffer.data(if_off + 265 * ccomps * dcomps);

            auto g_yzzzzz_yyy = cbuffer.data(if_off + 266 * ccomps * dcomps);

            auto g_yzzzzz_yyz = cbuffer.data(if_off + 267 * ccomps * dcomps);

            auto g_yzzzzz_yzz = cbuffer.data(if_off + 268 * ccomps * dcomps);

            auto g_yzzzzz_zzz = cbuffer.data(if_off + 269 * ccomps * dcomps);

            auto g_zzzzzz_xxx = cbuffer.data(if_off + 270 * ccomps * dcomps);

            auto g_zzzzzz_xxy = cbuffer.data(if_off + 271 * ccomps * dcomps);

            auto g_zzzzzz_xxz = cbuffer.data(if_off + 272 * ccomps * dcomps);

            auto g_zzzzzz_xyy = cbuffer.data(if_off + 273 * ccomps * dcomps);

            auto g_zzzzzz_xyz = cbuffer.data(if_off + 274 * ccomps * dcomps);

            auto g_zzzzzz_xzz = cbuffer.data(if_off + 275 * ccomps * dcomps);

            auto g_zzzzzz_yyy = cbuffer.data(if_off + 276 * ccomps * dcomps);

            auto g_zzzzzz_yyz = cbuffer.data(if_off + 277 * ccomps * dcomps);

            auto g_zzzzzz_yzz = cbuffer.data(if_off + 278 * ccomps * dcomps);

            auto g_zzzzzz_zzz = cbuffer.data(if_off + 279 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IFSS

            const auto if_geom_01_off = idx_geom_01_ifxx + i * dcomps + j;

            auto g_0_x_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 839 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_kfxx

            const auto kf_geom_01_off = idx_geom_01_kfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxx_xxx = cbuffer.data(kf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxy = cbuffer.data(kf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxz = cbuffer.data(kf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xyy = cbuffer.data(kf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xyz = cbuffer.data(kf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xzz = cbuffer.data(kf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yyy = cbuffer.data(kf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yyz = cbuffer.data(kf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yzz = cbuffer.data(kf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_zzz = cbuffer.data(kf_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxx, g_0_x_xxxxxx_xxxx, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxy, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xyy, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_yyy, g_0_x_xxxxxx_yyz, g_0_x_xxxxxx_yzz, g_0_x_xxxxxx_zzz, g_0_x_xxxxxxx_xxx, g_0_x_xxxxxxx_xxy, g_0_x_xxxxxxx_xxz, g_0_x_xxxxxxx_xyy, g_0_x_xxxxxxx_xyz, g_0_x_xxxxxxx_xzz, g_0_x_xxxxxxx_yyy, g_0_x_xxxxxxx_yyz, g_0_x_xxxxxxx_yzz, g_0_x_xxxxxxx_zzz, g_xxxxxx_xxx, g_xxxxxx_xxy, g_xxxxxx_xxz, g_xxxxxx_xyy, g_xxxxxx_xyz, g_xxxxxx_xzz, g_xxxxxx_yyy, g_xxxxxx_yyz, g_xxxxxx_yzz, g_xxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxx_xxx[k] = g_xxxxxx_xxx[k] - g_0_x_xxxxxx_xxx[k] * ab_x + g_0_x_xxxxxx_xxxx[k];

                g_0_x_xxxxxxx_xxy[k] = g_xxxxxx_xxy[k] - g_0_x_xxxxxx_xxy[k] * ab_x + g_0_x_xxxxxx_xxxy[k];

                g_0_x_xxxxxxx_xxz[k] = g_xxxxxx_xxz[k] - g_0_x_xxxxxx_xxz[k] * ab_x + g_0_x_xxxxxx_xxxz[k];

                g_0_x_xxxxxxx_xyy[k] = g_xxxxxx_xyy[k] - g_0_x_xxxxxx_xyy[k] * ab_x + g_0_x_xxxxxx_xxyy[k];

                g_0_x_xxxxxxx_xyz[k] = g_xxxxxx_xyz[k] - g_0_x_xxxxxx_xyz[k] * ab_x + g_0_x_xxxxxx_xxyz[k];

                g_0_x_xxxxxxx_xzz[k] = g_xxxxxx_xzz[k] - g_0_x_xxxxxx_xzz[k] * ab_x + g_0_x_xxxxxx_xxzz[k];

                g_0_x_xxxxxxx_yyy[k] = g_xxxxxx_yyy[k] - g_0_x_xxxxxx_yyy[k] * ab_x + g_0_x_xxxxxx_xyyy[k];

                g_0_x_xxxxxxx_yyz[k] = g_xxxxxx_yyz[k] - g_0_x_xxxxxx_yyz[k] * ab_x + g_0_x_xxxxxx_xyyz[k];

                g_0_x_xxxxxxx_yzz[k] = g_xxxxxx_yzz[k] - g_0_x_xxxxxx_yzz[k] * ab_x + g_0_x_xxxxxx_xyzz[k];

                g_0_x_xxxxxxx_zzz[k] = g_xxxxxx_zzz[k] - g_0_x_xxxxxx_zzz[k] * ab_x + g_0_x_xxxxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxy_xxx = cbuffer.data(kf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxy = cbuffer.data(kf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxz = cbuffer.data(kf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xyy = cbuffer.data(kf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xyz = cbuffer.data(kf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xzz = cbuffer.data(kf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yyy = cbuffer.data(kf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yyz = cbuffer.data(kf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yzz = cbuffer.data(kf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_zzz = cbuffer.data(kf_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxx, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxy, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxz, g_0_x_xxxxxx_xyy, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xzz, g_0_x_xxxxxx_yyy, g_0_x_xxxxxx_yyyy, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_zzz, g_0_x_xxxxxxy_xxx, g_0_x_xxxxxxy_xxy, g_0_x_xxxxxxy_xxz, g_0_x_xxxxxxy_xyy, g_0_x_xxxxxxy_xyz, g_0_x_xxxxxxy_xzz, g_0_x_xxxxxxy_yyy, g_0_x_xxxxxxy_yyz, g_0_x_xxxxxxy_yzz, g_0_x_xxxxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxy_xxx[k] = -g_0_x_xxxxxx_xxx[k] * ab_y + g_0_x_xxxxxx_xxxy[k];

                g_0_x_xxxxxxy_xxy[k] = -g_0_x_xxxxxx_xxy[k] * ab_y + g_0_x_xxxxxx_xxyy[k];

                g_0_x_xxxxxxy_xxz[k] = -g_0_x_xxxxxx_xxz[k] * ab_y + g_0_x_xxxxxx_xxyz[k];

                g_0_x_xxxxxxy_xyy[k] = -g_0_x_xxxxxx_xyy[k] * ab_y + g_0_x_xxxxxx_xyyy[k];

                g_0_x_xxxxxxy_xyz[k] = -g_0_x_xxxxxx_xyz[k] * ab_y + g_0_x_xxxxxx_xyyz[k];

                g_0_x_xxxxxxy_xzz[k] = -g_0_x_xxxxxx_xzz[k] * ab_y + g_0_x_xxxxxx_xyzz[k];

                g_0_x_xxxxxxy_yyy[k] = -g_0_x_xxxxxx_yyy[k] * ab_y + g_0_x_xxxxxx_yyyy[k];

                g_0_x_xxxxxxy_yyz[k] = -g_0_x_xxxxxx_yyz[k] * ab_y + g_0_x_xxxxxx_yyyz[k];

                g_0_x_xxxxxxy_yzz[k] = -g_0_x_xxxxxx_yzz[k] * ab_y + g_0_x_xxxxxx_yyzz[k];

                g_0_x_xxxxxxy_zzz[k] = -g_0_x_xxxxxx_zzz[k] * ab_y + g_0_x_xxxxxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxz_xxx = cbuffer.data(kf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxy = cbuffer.data(kf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxz = cbuffer.data(kf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xyy = cbuffer.data(kf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xyz = cbuffer.data(kf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xzz = cbuffer.data(kf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yyy = cbuffer.data(kf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yyz = cbuffer.data(kf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yzz = cbuffer.data(kf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_zzz = cbuffer.data(kf_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxx, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxy, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xyy, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_yyy, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_zzz, g_0_x_xxxxxx_zzzz, g_0_x_xxxxxxz_xxx, g_0_x_xxxxxxz_xxy, g_0_x_xxxxxxz_xxz, g_0_x_xxxxxxz_xyy, g_0_x_xxxxxxz_xyz, g_0_x_xxxxxxz_xzz, g_0_x_xxxxxxz_yyy, g_0_x_xxxxxxz_yyz, g_0_x_xxxxxxz_yzz, g_0_x_xxxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxz_xxx[k] = -g_0_x_xxxxxx_xxx[k] * ab_z + g_0_x_xxxxxx_xxxz[k];

                g_0_x_xxxxxxz_xxy[k] = -g_0_x_xxxxxx_xxy[k] * ab_z + g_0_x_xxxxxx_xxyz[k];

                g_0_x_xxxxxxz_xxz[k] = -g_0_x_xxxxxx_xxz[k] * ab_z + g_0_x_xxxxxx_xxzz[k];

                g_0_x_xxxxxxz_xyy[k] = -g_0_x_xxxxxx_xyy[k] * ab_z + g_0_x_xxxxxx_xyyz[k];

                g_0_x_xxxxxxz_xyz[k] = -g_0_x_xxxxxx_xyz[k] * ab_z + g_0_x_xxxxxx_xyzz[k];

                g_0_x_xxxxxxz_xzz[k] = -g_0_x_xxxxxx_xzz[k] * ab_z + g_0_x_xxxxxx_xzzz[k];

                g_0_x_xxxxxxz_yyy[k] = -g_0_x_xxxxxx_yyy[k] * ab_z + g_0_x_xxxxxx_yyyz[k];

                g_0_x_xxxxxxz_yyz[k] = -g_0_x_xxxxxx_yyz[k] * ab_z + g_0_x_xxxxxx_yyzz[k];

                g_0_x_xxxxxxz_yzz[k] = -g_0_x_xxxxxx_yzz[k] * ab_z + g_0_x_xxxxxx_yzzz[k];

                g_0_x_xxxxxxz_zzz[k] = -g_0_x_xxxxxx_zzz[k] * ab_z + g_0_x_xxxxxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyy_xxx = cbuffer.data(kf_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxy = cbuffer.data(kf_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxz = cbuffer.data(kf_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xyy = cbuffer.data(kf_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xyz = cbuffer.data(kf_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xzz = cbuffer.data(kf_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yyy = cbuffer.data(kf_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yyz = cbuffer.data(kf_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yzz = cbuffer.data(kf_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_zzz = cbuffer.data(kf_geom_01_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxy_xxx, g_0_x_xxxxxy_xxxy, g_0_x_xxxxxy_xxy, g_0_x_xxxxxy_xxyy, g_0_x_xxxxxy_xxyz, g_0_x_xxxxxy_xxz, g_0_x_xxxxxy_xyy, g_0_x_xxxxxy_xyyy, g_0_x_xxxxxy_xyyz, g_0_x_xxxxxy_xyz, g_0_x_xxxxxy_xyzz, g_0_x_xxxxxy_xzz, g_0_x_xxxxxy_yyy, g_0_x_xxxxxy_yyyy, g_0_x_xxxxxy_yyyz, g_0_x_xxxxxy_yyz, g_0_x_xxxxxy_yyzz, g_0_x_xxxxxy_yzz, g_0_x_xxxxxy_yzzz, g_0_x_xxxxxy_zzz, g_0_x_xxxxxyy_xxx, g_0_x_xxxxxyy_xxy, g_0_x_xxxxxyy_xxz, g_0_x_xxxxxyy_xyy, g_0_x_xxxxxyy_xyz, g_0_x_xxxxxyy_xzz, g_0_x_xxxxxyy_yyy, g_0_x_xxxxxyy_yyz, g_0_x_xxxxxyy_yzz, g_0_x_xxxxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyy_xxx[k] = -g_0_x_xxxxxy_xxx[k] * ab_y + g_0_x_xxxxxy_xxxy[k];

                g_0_x_xxxxxyy_xxy[k] = -g_0_x_xxxxxy_xxy[k] * ab_y + g_0_x_xxxxxy_xxyy[k];

                g_0_x_xxxxxyy_xxz[k] = -g_0_x_xxxxxy_xxz[k] * ab_y + g_0_x_xxxxxy_xxyz[k];

                g_0_x_xxxxxyy_xyy[k] = -g_0_x_xxxxxy_xyy[k] * ab_y + g_0_x_xxxxxy_xyyy[k];

                g_0_x_xxxxxyy_xyz[k] = -g_0_x_xxxxxy_xyz[k] * ab_y + g_0_x_xxxxxy_xyyz[k];

                g_0_x_xxxxxyy_xzz[k] = -g_0_x_xxxxxy_xzz[k] * ab_y + g_0_x_xxxxxy_xyzz[k];

                g_0_x_xxxxxyy_yyy[k] = -g_0_x_xxxxxy_yyy[k] * ab_y + g_0_x_xxxxxy_yyyy[k];

                g_0_x_xxxxxyy_yyz[k] = -g_0_x_xxxxxy_yyz[k] * ab_y + g_0_x_xxxxxy_yyyz[k];

                g_0_x_xxxxxyy_yzz[k] = -g_0_x_xxxxxy_yzz[k] * ab_y + g_0_x_xxxxxy_yyzz[k];

                g_0_x_xxxxxyy_zzz[k] = -g_0_x_xxxxxy_zzz[k] * ab_y + g_0_x_xxxxxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyz_xxx = cbuffer.data(kf_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxy = cbuffer.data(kf_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxz = cbuffer.data(kf_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xyy = cbuffer.data(kf_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xyz = cbuffer.data(kf_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xzz = cbuffer.data(kf_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yyy = cbuffer.data(kf_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yyz = cbuffer.data(kf_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yzz = cbuffer.data(kf_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_zzz = cbuffer.data(kf_geom_01_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxyz_xxx, g_0_x_xxxxxyz_xxy, g_0_x_xxxxxyz_xxz, g_0_x_xxxxxyz_xyy, g_0_x_xxxxxyz_xyz, g_0_x_xxxxxyz_xzz, g_0_x_xxxxxyz_yyy, g_0_x_xxxxxyz_yyz, g_0_x_xxxxxyz_yzz, g_0_x_xxxxxyz_zzz, g_0_x_xxxxxz_xxx, g_0_x_xxxxxz_xxxy, g_0_x_xxxxxz_xxy, g_0_x_xxxxxz_xxyy, g_0_x_xxxxxz_xxyz, g_0_x_xxxxxz_xxz, g_0_x_xxxxxz_xyy, g_0_x_xxxxxz_xyyy, g_0_x_xxxxxz_xyyz, g_0_x_xxxxxz_xyz, g_0_x_xxxxxz_xyzz, g_0_x_xxxxxz_xzz, g_0_x_xxxxxz_yyy, g_0_x_xxxxxz_yyyy, g_0_x_xxxxxz_yyyz, g_0_x_xxxxxz_yyz, g_0_x_xxxxxz_yyzz, g_0_x_xxxxxz_yzz, g_0_x_xxxxxz_yzzz, g_0_x_xxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyz_xxx[k] = -g_0_x_xxxxxz_xxx[k] * ab_y + g_0_x_xxxxxz_xxxy[k];

                g_0_x_xxxxxyz_xxy[k] = -g_0_x_xxxxxz_xxy[k] * ab_y + g_0_x_xxxxxz_xxyy[k];

                g_0_x_xxxxxyz_xxz[k] = -g_0_x_xxxxxz_xxz[k] * ab_y + g_0_x_xxxxxz_xxyz[k];

                g_0_x_xxxxxyz_xyy[k] = -g_0_x_xxxxxz_xyy[k] * ab_y + g_0_x_xxxxxz_xyyy[k];

                g_0_x_xxxxxyz_xyz[k] = -g_0_x_xxxxxz_xyz[k] * ab_y + g_0_x_xxxxxz_xyyz[k];

                g_0_x_xxxxxyz_xzz[k] = -g_0_x_xxxxxz_xzz[k] * ab_y + g_0_x_xxxxxz_xyzz[k];

                g_0_x_xxxxxyz_yyy[k] = -g_0_x_xxxxxz_yyy[k] * ab_y + g_0_x_xxxxxz_yyyy[k];

                g_0_x_xxxxxyz_yyz[k] = -g_0_x_xxxxxz_yyz[k] * ab_y + g_0_x_xxxxxz_yyyz[k];

                g_0_x_xxxxxyz_yzz[k] = -g_0_x_xxxxxz_yzz[k] * ab_y + g_0_x_xxxxxz_yyzz[k];

                g_0_x_xxxxxyz_zzz[k] = -g_0_x_xxxxxz_zzz[k] * ab_y + g_0_x_xxxxxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxzz_xxx = cbuffer.data(kf_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxy = cbuffer.data(kf_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxz = cbuffer.data(kf_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xyy = cbuffer.data(kf_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xyz = cbuffer.data(kf_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xzz = cbuffer.data(kf_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yyy = cbuffer.data(kf_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yyz = cbuffer.data(kf_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yzz = cbuffer.data(kf_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_zzz = cbuffer.data(kf_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxz_xxx, g_0_x_xxxxxz_xxxz, g_0_x_xxxxxz_xxy, g_0_x_xxxxxz_xxyz, g_0_x_xxxxxz_xxz, g_0_x_xxxxxz_xxzz, g_0_x_xxxxxz_xyy, g_0_x_xxxxxz_xyyz, g_0_x_xxxxxz_xyz, g_0_x_xxxxxz_xyzz, g_0_x_xxxxxz_xzz, g_0_x_xxxxxz_xzzz, g_0_x_xxxxxz_yyy, g_0_x_xxxxxz_yyyz, g_0_x_xxxxxz_yyz, g_0_x_xxxxxz_yyzz, g_0_x_xxxxxz_yzz, g_0_x_xxxxxz_yzzz, g_0_x_xxxxxz_zzz, g_0_x_xxxxxz_zzzz, g_0_x_xxxxxzz_xxx, g_0_x_xxxxxzz_xxy, g_0_x_xxxxxzz_xxz, g_0_x_xxxxxzz_xyy, g_0_x_xxxxxzz_xyz, g_0_x_xxxxxzz_xzz, g_0_x_xxxxxzz_yyy, g_0_x_xxxxxzz_yyz, g_0_x_xxxxxzz_yzz, g_0_x_xxxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxzz_xxx[k] = -g_0_x_xxxxxz_xxx[k] * ab_z + g_0_x_xxxxxz_xxxz[k];

                g_0_x_xxxxxzz_xxy[k] = -g_0_x_xxxxxz_xxy[k] * ab_z + g_0_x_xxxxxz_xxyz[k];

                g_0_x_xxxxxzz_xxz[k] = -g_0_x_xxxxxz_xxz[k] * ab_z + g_0_x_xxxxxz_xxzz[k];

                g_0_x_xxxxxzz_xyy[k] = -g_0_x_xxxxxz_xyy[k] * ab_z + g_0_x_xxxxxz_xyyz[k];

                g_0_x_xxxxxzz_xyz[k] = -g_0_x_xxxxxz_xyz[k] * ab_z + g_0_x_xxxxxz_xyzz[k];

                g_0_x_xxxxxzz_xzz[k] = -g_0_x_xxxxxz_xzz[k] * ab_z + g_0_x_xxxxxz_xzzz[k];

                g_0_x_xxxxxzz_yyy[k] = -g_0_x_xxxxxz_yyy[k] * ab_z + g_0_x_xxxxxz_yyyz[k];

                g_0_x_xxxxxzz_yyz[k] = -g_0_x_xxxxxz_yyz[k] * ab_z + g_0_x_xxxxxz_yyzz[k];

                g_0_x_xxxxxzz_yzz[k] = -g_0_x_xxxxxz_yzz[k] * ab_z + g_0_x_xxxxxz_yzzz[k];

                g_0_x_xxxxxzz_zzz[k] = -g_0_x_xxxxxz_zzz[k] * ab_z + g_0_x_xxxxxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyy_xxx = cbuffer.data(kf_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxy = cbuffer.data(kf_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxz = cbuffer.data(kf_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xyy = cbuffer.data(kf_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xyz = cbuffer.data(kf_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xzz = cbuffer.data(kf_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yyy = cbuffer.data(kf_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yyz = cbuffer.data(kf_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yzz = cbuffer.data(kf_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_zzz = cbuffer.data(kf_geom_01_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyy_xxx, g_0_x_xxxxyy_xxxy, g_0_x_xxxxyy_xxy, g_0_x_xxxxyy_xxyy, g_0_x_xxxxyy_xxyz, g_0_x_xxxxyy_xxz, g_0_x_xxxxyy_xyy, g_0_x_xxxxyy_xyyy, g_0_x_xxxxyy_xyyz, g_0_x_xxxxyy_xyz, g_0_x_xxxxyy_xyzz, g_0_x_xxxxyy_xzz, g_0_x_xxxxyy_yyy, g_0_x_xxxxyy_yyyy, g_0_x_xxxxyy_yyyz, g_0_x_xxxxyy_yyz, g_0_x_xxxxyy_yyzz, g_0_x_xxxxyy_yzz, g_0_x_xxxxyy_yzzz, g_0_x_xxxxyy_zzz, g_0_x_xxxxyyy_xxx, g_0_x_xxxxyyy_xxy, g_0_x_xxxxyyy_xxz, g_0_x_xxxxyyy_xyy, g_0_x_xxxxyyy_xyz, g_0_x_xxxxyyy_xzz, g_0_x_xxxxyyy_yyy, g_0_x_xxxxyyy_yyz, g_0_x_xxxxyyy_yzz, g_0_x_xxxxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyy_xxx[k] = -g_0_x_xxxxyy_xxx[k] * ab_y + g_0_x_xxxxyy_xxxy[k];

                g_0_x_xxxxyyy_xxy[k] = -g_0_x_xxxxyy_xxy[k] * ab_y + g_0_x_xxxxyy_xxyy[k];

                g_0_x_xxxxyyy_xxz[k] = -g_0_x_xxxxyy_xxz[k] * ab_y + g_0_x_xxxxyy_xxyz[k];

                g_0_x_xxxxyyy_xyy[k] = -g_0_x_xxxxyy_xyy[k] * ab_y + g_0_x_xxxxyy_xyyy[k];

                g_0_x_xxxxyyy_xyz[k] = -g_0_x_xxxxyy_xyz[k] * ab_y + g_0_x_xxxxyy_xyyz[k];

                g_0_x_xxxxyyy_xzz[k] = -g_0_x_xxxxyy_xzz[k] * ab_y + g_0_x_xxxxyy_xyzz[k];

                g_0_x_xxxxyyy_yyy[k] = -g_0_x_xxxxyy_yyy[k] * ab_y + g_0_x_xxxxyy_yyyy[k];

                g_0_x_xxxxyyy_yyz[k] = -g_0_x_xxxxyy_yyz[k] * ab_y + g_0_x_xxxxyy_yyyz[k];

                g_0_x_xxxxyyy_yzz[k] = -g_0_x_xxxxyy_yzz[k] * ab_y + g_0_x_xxxxyy_yyzz[k];

                g_0_x_xxxxyyy_zzz[k] = -g_0_x_xxxxyy_zzz[k] * ab_y + g_0_x_xxxxyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyz_xxx = cbuffer.data(kf_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxy = cbuffer.data(kf_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxz = cbuffer.data(kf_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xyy = cbuffer.data(kf_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xyz = cbuffer.data(kf_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xzz = cbuffer.data(kf_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yyy = cbuffer.data(kf_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yyz = cbuffer.data(kf_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yzz = cbuffer.data(kf_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_zzz = cbuffer.data(kf_geom_01_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyyz_xxx, g_0_x_xxxxyyz_xxy, g_0_x_xxxxyyz_xxz, g_0_x_xxxxyyz_xyy, g_0_x_xxxxyyz_xyz, g_0_x_xxxxyyz_xzz, g_0_x_xxxxyyz_yyy, g_0_x_xxxxyyz_yyz, g_0_x_xxxxyyz_yzz, g_0_x_xxxxyyz_zzz, g_0_x_xxxxyz_xxx, g_0_x_xxxxyz_xxxy, g_0_x_xxxxyz_xxy, g_0_x_xxxxyz_xxyy, g_0_x_xxxxyz_xxyz, g_0_x_xxxxyz_xxz, g_0_x_xxxxyz_xyy, g_0_x_xxxxyz_xyyy, g_0_x_xxxxyz_xyyz, g_0_x_xxxxyz_xyz, g_0_x_xxxxyz_xyzz, g_0_x_xxxxyz_xzz, g_0_x_xxxxyz_yyy, g_0_x_xxxxyz_yyyy, g_0_x_xxxxyz_yyyz, g_0_x_xxxxyz_yyz, g_0_x_xxxxyz_yyzz, g_0_x_xxxxyz_yzz, g_0_x_xxxxyz_yzzz, g_0_x_xxxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyz_xxx[k] = -g_0_x_xxxxyz_xxx[k] * ab_y + g_0_x_xxxxyz_xxxy[k];

                g_0_x_xxxxyyz_xxy[k] = -g_0_x_xxxxyz_xxy[k] * ab_y + g_0_x_xxxxyz_xxyy[k];

                g_0_x_xxxxyyz_xxz[k] = -g_0_x_xxxxyz_xxz[k] * ab_y + g_0_x_xxxxyz_xxyz[k];

                g_0_x_xxxxyyz_xyy[k] = -g_0_x_xxxxyz_xyy[k] * ab_y + g_0_x_xxxxyz_xyyy[k];

                g_0_x_xxxxyyz_xyz[k] = -g_0_x_xxxxyz_xyz[k] * ab_y + g_0_x_xxxxyz_xyyz[k];

                g_0_x_xxxxyyz_xzz[k] = -g_0_x_xxxxyz_xzz[k] * ab_y + g_0_x_xxxxyz_xyzz[k];

                g_0_x_xxxxyyz_yyy[k] = -g_0_x_xxxxyz_yyy[k] * ab_y + g_0_x_xxxxyz_yyyy[k];

                g_0_x_xxxxyyz_yyz[k] = -g_0_x_xxxxyz_yyz[k] * ab_y + g_0_x_xxxxyz_yyyz[k];

                g_0_x_xxxxyyz_yzz[k] = -g_0_x_xxxxyz_yzz[k] * ab_y + g_0_x_xxxxyz_yyzz[k];

                g_0_x_xxxxyyz_zzz[k] = -g_0_x_xxxxyz_zzz[k] * ab_y + g_0_x_xxxxyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyzz_xxx = cbuffer.data(kf_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxy = cbuffer.data(kf_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxz = cbuffer.data(kf_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xyy = cbuffer.data(kf_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xyz = cbuffer.data(kf_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xzz = cbuffer.data(kf_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yyy = cbuffer.data(kf_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yyz = cbuffer.data(kf_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yzz = cbuffer.data(kf_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_zzz = cbuffer.data(kf_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyzz_xxx, g_0_x_xxxxyzz_xxy, g_0_x_xxxxyzz_xxz, g_0_x_xxxxyzz_xyy, g_0_x_xxxxyzz_xyz, g_0_x_xxxxyzz_xzz, g_0_x_xxxxyzz_yyy, g_0_x_xxxxyzz_yyz, g_0_x_xxxxyzz_yzz, g_0_x_xxxxyzz_zzz, g_0_x_xxxxzz_xxx, g_0_x_xxxxzz_xxxy, g_0_x_xxxxzz_xxy, g_0_x_xxxxzz_xxyy, g_0_x_xxxxzz_xxyz, g_0_x_xxxxzz_xxz, g_0_x_xxxxzz_xyy, g_0_x_xxxxzz_xyyy, g_0_x_xxxxzz_xyyz, g_0_x_xxxxzz_xyz, g_0_x_xxxxzz_xyzz, g_0_x_xxxxzz_xzz, g_0_x_xxxxzz_yyy, g_0_x_xxxxzz_yyyy, g_0_x_xxxxzz_yyyz, g_0_x_xxxxzz_yyz, g_0_x_xxxxzz_yyzz, g_0_x_xxxxzz_yzz, g_0_x_xxxxzz_yzzz, g_0_x_xxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyzz_xxx[k] = -g_0_x_xxxxzz_xxx[k] * ab_y + g_0_x_xxxxzz_xxxy[k];

                g_0_x_xxxxyzz_xxy[k] = -g_0_x_xxxxzz_xxy[k] * ab_y + g_0_x_xxxxzz_xxyy[k];

                g_0_x_xxxxyzz_xxz[k] = -g_0_x_xxxxzz_xxz[k] * ab_y + g_0_x_xxxxzz_xxyz[k];

                g_0_x_xxxxyzz_xyy[k] = -g_0_x_xxxxzz_xyy[k] * ab_y + g_0_x_xxxxzz_xyyy[k];

                g_0_x_xxxxyzz_xyz[k] = -g_0_x_xxxxzz_xyz[k] * ab_y + g_0_x_xxxxzz_xyyz[k];

                g_0_x_xxxxyzz_xzz[k] = -g_0_x_xxxxzz_xzz[k] * ab_y + g_0_x_xxxxzz_xyzz[k];

                g_0_x_xxxxyzz_yyy[k] = -g_0_x_xxxxzz_yyy[k] * ab_y + g_0_x_xxxxzz_yyyy[k];

                g_0_x_xxxxyzz_yyz[k] = -g_0_x_xxxxzz_yyz[k] * ab_y + g_0_x_xxxxzz_yyyz[k];

                g_0_x_xxxxyzz_yzz[k] = -g_0_x_xxxxzz_yzz[k] * ab_y + g_0_x_xxxxzz_yyzz[k];

                g_0_x_xxxxyzz_zzz[k] = -g_0_x_xxxxzz_zzz[k] * ab_y + g_0_x_xxxxzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzzz_xxx = cbuffer.data(kf_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxy = cbuffer.data(kf_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxz = cbuffer.data(kf_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xyy = cbuffer.data(kf_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xyz = cbuffer.data(kf_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xzz = cbuffer.data(kf_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yyy = cbuffer.data(kf_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yyz = cbuffer.data(kf_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yzz = cbuffer.data(kf_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_zzz = cbuffer.data(kf_geom_01_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxzz_xxx, g_0_x_xxxxzz_xxxz, g_0_x_xxxxzz_xxy, g_0_x_xxxxzz_xxyz, g_0_x_xxxxzz_xxz, g_0_x_xxxxzz_xxzz, g_0_x_xxxxzz_xyy, g_0_x_xxxxzz_xyyz, g_0_x_xxxxzz_xyz, g_0_x_xxxxzz_xyzz, g_0_x_xxxxzz_xzz, g_0_x_xxxxzz_xzzz, g_0_x_xxxxzz_yyy, g_0_x_xxxxzz_yyyz, g_0_x_xxxxzz_yyz, g_0_x_xxxxzz_yyzz, g_0_x_xxxxzz_yzz, g_0_x_xxxxzz_yzzz, g_0_x_xxxxzz_zzz, g_0_x_xxxxzz_zzzz, g_0_x_xxxxzzz_xxx, g_0_x_xxxxzzz_xxy, g_0_x_xxxxzzz_xxz, g_0_x_xxxxzzz_xyy, g_0_x_xxxxzzz_xyz, g_0_x_xxxxzzz_xzz, g_0_x_xxxxzzz_yyy, g_0_x_xxxxzzz_yyz, g_0_x_xxxxzzz_yzz, g_0_x_xxxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzzz_xxx[k] = -g_0_x_xxxxzz_xxx[k] * ab_z + g_0_x_xxxxzz_xxxz[k];

                g_0_x_xxxxzzz_xxy[k] = -g_0_x_xxxxzz_xxy[k] * ab_z + g_0_x_xxxxzz_xxyz[k];

                g_0_x_xxxxzzz_xxz[k] = -g_0_x_xxxxzz_xxz[k] * ab_z + g_0_x_xxxxzz_xxzz[k];

                g_0_x_xxxxzzz_xyy[k] = -g_0_x_xxxxzz_xyy[k] * ab_z + g_0_x_xxxxzz_xyyz[k];

                g_0_x_xxxxzzz_xyz[k] = -g_0_x_xxxxzz_xyz[k] * ab_z + g_0_x_xxxxzz_xyzz[k];

                g_0_x_xxxxzzz_xzz[k] = -g_0_x_xxxxzz_xzz[k] * ab_z + g_0_x_xxxxzz_xzzz[k];

                g_0_x_xxxxzzz_yyy[k] = -g_0_x_xxxxzz_yyy[k] * ab_z + g_0_x_xxxxzz_yyyz[k];

                g_0_x_xxxxzzz_yyz[k] = -g_0_x_xxxxzz_yyz[k] * ab_z + g_0_x_xxxxzz_yyzz[k];

                g_0_x_xxxxzzz_yzz[k] = -g_0_x_xxxxzz_yzz[k] * ab_z + g_0_x_xxxxzz_yzzz[k];

                g_0_x_xxxxzzz_zzz[k] = -g_0_x_xxxxzz_zzz[k] * ab_z + g_0_x_xxxxzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyy_xxx = cbuffer.data(kf_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxy = cbuffer.data(kf_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxz = cbuffer.data(kf_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xyy = cbuffer.data(kf_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xyz = cbuffer.data(kf_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xzz = cbuffer.data(kf_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yyy = cbuffer.data(kf_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yyz = cbuffer.data(kf_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yzz = cbuffer.data(kf_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_zzz = cbuffer.data(kf_geom_01_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyy_xxx, g_0_x_xxxyyy_xxxy, g_0_x_xxxyyy_xxy, g_0_x_xxxyyy_xxyy, g_0_x_xxxyyy_xxyz, g_0_x_xxxyyy_xxz, g_0_x_xxxyyy_xyy, g_0_x_xxxyyy_xyyy, g_0_x_xxxyyy_xyyz, g_0_x_xxxyyy_xyz, g_0_x_xxxyyy_xyzz, g_0_x_xxxyyy_xzz, g_0_x_xxxyyy_yyy, g_0_x_xxxyyy_yyyy, g_0_x_xxxyyy_yyyz, g_0_x_xxxyyy_yyz, g_0_x_xxxyyy_yyzz, g_0_x_xxxyyy_yzz, g_0_x_xxxyyy_yzzz, g_0_x_xxxyyy_zzz, g_0_x_xxxyyyy_xxx, g_0_x_xxxyyyy_xxy, g_0_x_xxxyyyy_xxz, g_0_x_xxxyyyy_xyy, g_0_x_xxxyyyy_xyz, g_0_x_xxxyyyy_xzz, g_0_x_xxxyyyy_yyy, g_0_x_xxxyyyy_yyz, g_0_x_xxxyyyy_yzz, g_0_x_xxxyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyy_xxx[k] = -g_0_x_xxxyyy_xxx[k] * ab_y + g_0_x_xxxyyy_xxxy[k];

                g_0_x_xxxyyyy_xxy[k] = -g_0_x_xxxyyy_xxy[k] * ab_y + g_0_x_xxxyyy_xxyy[k];

                g_0_x_xxxyyyy_xxz[k] = -g_0_x_xxxyyy_xxz[k] * ab_y + g_0_x_xxxyyy_xxyz[k];

                g_0_x_xxxyyyy_xyy[k] = -g_0_x_xxxyyy_xyy[k] * ab_y + g_0_x_xxxyyy_xyyy[k];

                g_0_x_xxxyyyy_xyz[k] = -g_0_x_xxxyyy_xyz[k] * ab_y + g_0_x_xxxyyy_xyyz[k];

                g_0_x_xxxyyyy_xzz[k] = -g_0_x_xxxyyy_xzz[k] * ab_y + g_0_x_xxxyyy_xyzz[k];

                g_0_x_xxxyyyy_yyy[k] = -g_0_x_xxxyyy_yyy[k] * ab_y + g_0_x_xxxyyy_yyyy[k];

                g_0_x_xxxyyyy_yyz[k] = -g_0_x_xxxyyy_yyz[k] * ab_y + g_0_x_xxxyyy_yyyz[k];

                g_0_x_xxxyyyy_yzz[k] = -g_0_x_xxxyyy_yzz[k] * ab_y + g_0_x_xxxyyy_yyzz[k];

                g_0_x_xxxyyyy_zzz[k] = -g_0_x_xxxyyy_zzz[k] * ab_y + g_0_x_xxxyyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyz_xxx = cbuffer.data(kf_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxy = cbuffer.data(kf_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxz = cbuffer.data(kf_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xyy = cbuffer.data(kf_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xyz = cbuffer.data(kf_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xzz = cbuffer.data(kf_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yyy = cbuffer.data(kf_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yyz = cbuffer.data(kf_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yzz = cbuffer.data(kf_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_zzz = cbuffer.data(kf_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyyz_xxx, g_0_x_xxxyyyz_xxy, g_0_x_xxxyyyz_xxz, g_0_x_xxxyyyz_xyy, g_0_x_xxxyyyz_xyz, g_0_x_xxxyyyz_xzz, g_0_x_xxxyyyz_yyy, g_0_x_xxxyyyz_yyz, g_0_x_xxxyyyz_yzz, g_0_x_xxxyyyz_zzz, g_0_x_xxxyyz_xxx, g_0_x_xxxyyz_xxxy, g_0_x_xxxyyz_xxy, g_0_x_xxxyyz_xxyy, g_0_x_xxxyyz_xxyz, g_0_x_xxxyyz_xxz, g_0_x_xxxyyz_xyy, g_0_x_xxxyyz_xyyy, g_0_x_xxxyyz_xyyz, g_0_x_xxxyyz_xyz, g_0_x_xxxyyz_xyzz, g_0_x_xxxyyz_xzz, g_0_x_xxxyyz_yyy, g_0_x_xxxyyz_yyyy, g_0_x_xxxyyz_yyyz, g_0_x_xxxyyz_yyz, g_0_x_xxxyyz_yyzz, g_0_x_xxxyyz_yzz, g_0_x_xxxyyz_yzzz, g_0_x_xxxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyz_xxx[k] = -g_0_x_xxxyyz_xxx[k] * ab_y + g_0_x_xxxyyz_xxxy[k];

                g_0_x_xxxyyyz_xxy[k] = -g_0_x_xxxyyz_xxy[k] * ab_y + g_0_x_xxxyyz_xxyy[k];

                g_0_x_xxxyyyz_xxz[k] = -g_0_x_xxxyyz_xxz[k] * ab_y + g_0_x_xxxyyz_xxyz[k];

                g_0_x_xxxyyyz_xyy[k] = -g_0_x_xxxyyz_xyy[k] * ab_y + g_0_x_xxxyyz_xyyy[k];

                g_0_x_xxxyyyz_xyz[k] = -g_0_x_xxxyyz_xyz[k] * ab_y + g_0_x_xxxyyz_xyyz[k];

                g_0_x_xxxyyyz_xzz[k] = -g_0_x_xxxyyz_xzz[k] * ab_y + g_0_x_xxxyyz_xyzz[k];

                g_0_x_xxxyyyz_yyy[k] = -g_0_x_xxxyyz_yyy[k] * ab_y + g_0_x_xxxyyz_yyyy[k];

                g_0_x_xxxyyyz_yyz[k] = -g_0_x_xxxyyz_yyz[k] * ab_y + g_0_x_xxxyyz_yyyz[k];

                g_0_x_xxxyyyz_yzz[k] = -g_0_x_xxxyyz_yzz[k] * ab_y + g_0_x_xxxyyz_yyzz[k];

                g_0_x_xxxyyyz_zzz[k] = -g_0_x_xxxyyz_zzz[k] * ab_y + g_0_x_xxxyyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyzz_xxx = cbuffer.data(kf_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxy = cbuffer.data(kf_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxz = cbuffer.data(kf_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xyy = cbuffer.data(kf_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xyz = cbuffer.data(kf_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xzz = cbuffer.data(kf_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yyy = cbuffer.data(kf_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yyz = cbuffer.data(kf_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yzz = cbuffer.data(kf_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_zzz = cbuffer.data(kf_geom_01_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyzz_xxx, g_0_x_xxxyyzz_xxy, g_0_x_xxxyyzz_xxz, g_0_x_xxxyyzz_xyy, g_0_x_xxxyyzz_xyz, g_0_x_xxxyyzz_xzz, g_0_x_xxxyyzz_yyy, g_0_x_xxxyyzz_yyz, g_0_x_xxxyyzz_yzz, g_0_x_xxxyyzz_zzz, g_0_x_xxxyzz_xxx, g_0_x_xxxyzz_xxxy, g_0_x_xxxyzz_xxy, g_0_x_xxxyzz_xxyy, g_0_x_xxxyzz_xxyz, g_0_x_xxxyzz_xxz, g_0_x_xxxyzz_xyy, g_0_x_xxxyzz_xyyy, g_0_x_xxxyzz_xyyz, g_0_x_xxxyzz_xyz, g_0_x_xxxyzz_xyzz, g_0_x_xxxyzz_xzz, g_0_x_xxxyzz_yyy, g_0_x_xxxyzz_yyyy, g_0_x_xxxyzz_yyyz, g_0_x_xxxyzz_yyz, g_0_x_xxxyzz_yyzz, g_0_x_xxxyzz_yzz, g_0_x_xxxyzz_yzzz, g_0_x_xxxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyzz_xxx[k] = -g_0_x_xxxyzz_xxx[k] * ab_y + g_0_x_xxxyzz_xxxy[k];

                g_0_x_xxxyyzz_xxy[k] = -g_0_x_xxxyzz_xxy[k] * ab_y + g_0_x_xxxyzz_xxyy[k];

                g_0_x_xxxyyzz_xxz[k] = -g_0_x_xxxyzz_xxz[k] * ab_y + g_0_x_xxxyzz_xxyz[k];

                g_0_x_xxxyyzz_xyy[k] = -g_0_x_xxxyzz_xyy[k] * ab_y + g_0_x_xxxyzz_xyyy[k];

                g_0_x_xxxyyzz_xyz[k] = -g_0_x_xxxyzz_xyz[k] * ab_y + g_0_x_xxxyzz_xyyz[k];

                g_0_x_xxxyyzz_xzz[k] = -g_0_x_xxxyzz_xzz[k] * ab_y + g_0_x_xxxyzz_xyzz[k];

                g_0_x_xxxyyzz_yyy[k] = -g_0_x_xxxyzz_yyy[k] * ab_y + g_0_x_xxxyzz_yyyy[k];

                g_0_x_xxxyyzz_yyz[k] = -g_0_x_xxxyzz_yyz[k] * ab_y + g_0_x_xxxyzz_yyyz[k];

                g_0_x_xxxyyzz_yzz[k] = -g_0_x_xxxyzz_yzz[k] * ab_y + g_0_x_xxxyzz_yyzz[k];

                g_0_x_xxxyyzz_zzz[k] = -g_0_x_xxxyzz_zzz[k] * ab_y + g_0_x_xxxyzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzzz_xxx = cbuffer.data(kf_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxy = cbuffer.data(kf_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxz = cbuffer.data(kf_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xyy = cbuffer.data(kf_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xyz = cbuffer.data(kf_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xzz = cbuffer.data(kf_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yyy = cbuffer.data(kf_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yyz = cbuffer.data(kf_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yzz = cbuffer.data(kf_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_zzz = cbuffer.data(kf_geom_01_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzzz_xxx, g_0_x_xxxyzzz_xxy, g_0_x_xxxyzzz_xxz, g_0_x_xxxyzzz_xyy, g_0_x_xxxyzzz_xyz, g_0_x_xxxyzzz_xzz, g_0_x_xxxyzzz_yyy, g_0_x_xxxyzzz_yyz, g_0_x_xxxyzzz_yzz, g_0_x_xxxyzzz_zzz, g_0_x_xxxzzz_xxx, g_0_x_xxxzzz_xxxy, g_0_x_xxxzzz_xxy, g_0_x_xxxzzz_xxyy, g_0_x_xxxzzz_xxyz, g_0_x_xxxzzz_xxz, g_0_x_xxxzzz_xyy, g_0_x_xxxzzz_xyyy, g_0_x_xxxzzz_xyyz, g_0_x_xxxzzz_xyz, g_0_x_xxxzzz_xyzz, g_0_x_xxxzzz_xzz, g_0_x_xxxzzz_yyy, g_0_x_xxxzzz_yyyy, g_0_x_xxxzzz_yyyz, g_0_x_xxxzzz_yyz, g_0_x_xxxzzz_yyzz, g_0_x_xxxzzz_yzz, g_0_x_xxxzzz_yzzz, g_0_x_xxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzzz_xxx[k] = -g_0_x_xxxzzz_xxx[k] * ab_y + g_0_x_xxxzzz_xxxy[k];

                g_0_x_xxxyzzz_xxy[k] = -g_0_x_xxxzzz_xxy[k] * ab_y + g_0_x_xxxzzz_xxyy[k];

                g_0_x_xxxyzzz_xxz[k] = -g_0_x_xxxzzz_xxz[k] * ab_y + g_0_x_xxxzzz_xxyz[k];

                g_0_x_xxxyzzz_xyy[k] = -g_0_x_xxxzzz_xyy[k] * ab_y + g_0_x_xxxzzz_xyyy[k];

                g_0_x_xxxyzzz_xyz[k] = -g_0_x_xxxzzz_xyz[k] * ab_y + g_0_x_xxxzzz_xyyz[k];

                g_0_x_xxxyzzz_xzz[k] = -g_0_x_xxxzzz_xzz[k] * ab_y + g_0_x_xxxzzz_xyzz[k];

                g_0_x_xxxyzzz_yyy[k] = -g_0_x_xxxzzz_yyy[k] * ab_y + g_0_x_xxxzzz_yyyy[k];

                g_0_x_xxxyzzz_yyz[k] = -g_0_x_xxxzzz_yyz[k] * ab_y + g_0_x_xxxzzz_yyyz[k];

                g_0_x_xxxyzzz_yzz[k] = -g_0_x_xxxzzz_yzz[k] * ab_y + g_0_x_xxxzzz_yyzz[k];

                g_0_x_xxxyzzz_zzz[k] = -g_0_x_xxxzzz_zzz[k] * ab_y + g_0_x_xxxzzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzzz_xxx = cbuffer.data(kf_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxy = cbuffer.data(kf_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxz = cbuffer.data(kf_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xyy = cbuffer.data(kf_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xyz = cbuffer.data(kf_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xzz = cbuffer.data(kf_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yyy = cbuffer.data(kf_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yyz = cbuffer.data(kf_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yzz = cbuffer.data(kf_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_zzz = cbuffer.data(kf_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzzz_xxx, g_0_x_xxxzzz_xxxz, g_0_x_xxxzzz_xxy, g_0_x_xxxzzz_xxyz, g_0_x_xxxzzz_xxz, g_0_x_xxxzzz_xxzz, g_0_x_xxxzzz_xyy, g_0_x_xxxzzz_xyyz, g_0_x_xxxzzz_xyz, g_0_x_xxxzzz_xyzz, g_0_x_xxxzzz_xzz, g_0_x_xxxzzz_xzzz, g_0_x_xxxzzz_yyy, g_0_x_xxxzzz_yyyz, g_0_x_xxxzzz_yyz, g_0_x_xxxzzz_yyzz, g_0_x_xxxzzz_yzz, g_0_x_xxxzzz_yzzz, g_0_x_xxxzzz_zzz, g_0_x_xxxzzz_zzzz, g_0_x_xxxzzzz_xxx, g_0_x_xxxzzzz_xxy, g_0_x_xxxzzzz_xxz, g_0_x_xxxzzzz_xyy, g_0_x_xxxzzzz_xyz, g_0_x_xxxzzzz_xzz, g_0_x_xxxzzzz_yyy, g_0_x_xxxzzzz_yyz, g_0_x_xxxzzzz_yzz, g_0_x_xxxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzzz_xxx[k] = -g_0_x_xxxzzz_xxx[k] * ab_z + g_0_x_xxxzzz_xxxz[k];

                g_0_x_xxxzzzz_xxy[k] = -g_0_x_xxxzzz_xxy[k] * ab_z + g_0_x_xxxzzz_xxyz[k];

                g_0_x_xxxzzzz_xxz[k] = -g_0_x_xxxzzz_xxz[k] * ab_z + g_0_x_xxxzzz_xxzz[k];

                g_0_x_xxxzzzz_xyy[k] = -g_0_x_xxxzzz_xyy[k] * ab_z + g_0_x_xxxzzz_xyyz[k];

                g_0_x_xxxzzzz_xyz[k] = -g_0_x_xxxzzz_xyz[k] * ab_z + g_0_x_xxxzzz_xyzz[k];

                g_0_x_xxxzzzz_xzz[k] = -g_0_x_xxxzzz_xzz[k] * ab_z + g_0_x_xxxzzz_xzzz[k];

                g_0_x_xxxzzzz_yyy[k] = -g_0_x_xxxzzz_yyy[k] * ab_z + g_0_x_xxxzzz_yyyz[k];

                g_0_x_xxxzzzz_yyz[k] = -g_0_x_xxxzzz_yyz[k] * ab_z + g_0_x_xxxzzz_yyzz[k];

                g_0_x_xxxzzzz_yzz[k] = -g_0_x_xxxzzz_yzz[k] * ab_z + g_0_x_xxxzzz_yzzz[k];

                g_0_x_xxxzzzz_zzz[k] = -g_0_x_xxxzzz_zzz[k] * ab_z + g_0_x_xxxzzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyy_xxx = cbuffer.data(kf_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxy = cbuffer.data(kf_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxz = cbuffer.data(kf_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xyy = cbuffer.data(kf_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xyz = cbuffer.data(kf_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xzz = cbuffer.data(kf_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yyy = cbuffer.data(kf_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yyz = cbuffer.data(kf_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yzz = cbuffer.data(kf_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_zzz = cbuffer.data(kf_geom_01_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyy_xxx, g_0_x_xxyyyy_xxxy, g_0_x_xxyyyy_xxy, g_0_x_xxyyyy_xxyy, g_0_x_xxyyyy_xxyz, g_0_x_xxyyyy_xxz, g_0_x_xxyyyy_xyy, g_0_x_xxyyyy_xyyy, g_0_x_xxyyyy_xyyz, g_0_x_xxyyyy_xyz, g_0_x_xxyyyy_xyzz, g_0_x_xxyyyy_xzz, g_0_x_xxyyyy_yyy, g_0_x_xxyyyy_yyyy, g_0_x_xxyyyy_yyyz, g_0_x_xxyyyy_yyz, g_0_x_xxyyyy_yyzz, g_0_x_xxyyyy_yzz, g_0_x_xxyyyy_yzzz, g_0_x_xxyyyy_zzz, g_0_x_xxyyyyy_xxx, g_0_x_xxyyyyy_xxy, g_0_x_xxyyyyy_xxz, g_0_x_xxyyyyy_xyy, g_0_x_xxyyyyy_xyz, g_0_x_xxyyyyy_xzz, g_0_x_xxyyyyy_yyy, g_0_x_xxyyyyy_yyz, g_0_x_xxyyyyy_yzz, g_0_x_xxyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyy_xxx[k] = -g_0_x_xxyyyy_xxx[k] * ab_y + g_0_x_xxyyyy_xxxy[k];

                g_0_x_xxyyyyy_xxy[k] = -g_0_x_xxyyyy_xxy[k] * ab_y + g_0_x_xxyyyy_xxyy[k];

                g_0_x_xxyyyyy_xxz[k] = -g_0_x_xxyyyy_xxz[k] * ab_y + g_0_x_xxyyyy_xxyz[k];

                g_0_x_xxyyyyy_xyy[k] = -g_0_x_xxyyyy_xyy[k] * ab_y + g_0_x_xxyyyy_xyyy[k];

                g_0_x_xxyyyyy_xyz[k] = -g_0_x_xxyyyy_xyz[k] * ab_y + g_0_x_xxyyyy_xyyz[k];

                g_0_x_xxyyyyy_xzz[k] = -g_0_x_xxyyyy_xzz[k] * ab_y + g_0_x_xxyyyy_xyzz[k];

                g_0_x_xxyyyyy_yyy[k] = -g_0_x_xxyyyy_yyy[k] * ab_y + g_0_x_xxyyyy_yyyy[k];

                g_0_x_xxyyyyy_yyz[k] = -g_0_x_xxyyyy_yyz[k] * ab_y + g_0_x_xxyyyy_yyyz[k];

                g_0_x_xxyyyyy_yzz[k] = -g_0_x_xxyyyy_yzz[k] * ab_y + g_0_x_xxyyyy_yyzz[k];

                g_0_x_xxyyyyy_zzz[k] = -g_0_x_xxyyyy_zzz[k] * ab_y + g_0_x_xxyyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyz_xxx = cbuffer.data(kf_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxy = cbuffer.data(kf_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxz = cbuffer.data(kf_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xyy = cbuffer.data(kf_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xyz = cbuffer.data(kf_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xzz = cbuffer.data(kf_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yyy = cbuffer.data(kf_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yyz = cbuffer.data(kf_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yzz = cbuffer.data(kf_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_zzz = cbuffer.data(kf_geom_01_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyyz_xxx, g_0_x_xxyyyyz_xxy, g_0_x_xxyyyyz_xxz, g_0_x_xxyyyyz_xyy, g_0_x_xxyyyyz_xyz, g_0_x_xxyyyyz_xzz, g_0_x_xxyyyyz_yyy, g_0_x_xxyyyyz_yyz, g_0_x_xxyyyyz_yzz, g_0_x_xxyyyyz_zzz, g_0_x_xxyyyz_xxx, g_0_x_xxyyyz_xxxy, g_0_x_xxyyyz_xxy, g_0_x_xxyyyz_xxyy, g_0_x_xxyyyz_xxyz, g_0_x_xxyyyz_xxz, g_0_x_xxyyyz_xyy, g_0_x_xxyyyz_xyyy, g_0_x_xxyyyz_xyyz, g_0_x_xxyyyz_xyz, g_0_x_xxyyyz_xyzz, g_0_x_xxyyyz_xzz, g_0_x_xxyyyz_yyy, g_0_x_xxyyyz_yyyy, g_0_x_xxyyyz_yyyz, g_0_x_xxyyyz_yyz, g_0_x_xxyyyz_yyzz, g_0_x_xxyyyz_yzz, g_0_x_xxyyyz_yzzz, g_0_x_xxyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyz_xxx[k] = -g_0_x_xxyyyz_xxx[k] * ab_y + g_0_x_xxyyyz_xxxy[k];

                g_0_x_xxyyyyz_xxy[k] = -g_0_x_xxyyyz_xxy[k] * ab_y + g_0_x_xxyyyz_xxyy[k];

                g_0_x_xxyyyyz_xxz[k] = -g_0_x_xxyyyz_xxz[k] * ab_y + g_0_x_xxyyyz_xxyz[k];

                g_0_x_xxyyyyz_xyy[k] = -g_0_x_xxyyyz_xyy[k] * ab_y + g_0_x_xxyyyz_xyyy[k];

                g_0_x_xxyyyyz_xyz[k] = -g_0_x_xxyyyz_xyz[k] * ab_y + g_0_x_xxyyyz_xyyz[k];

                g_0_x_xxyyyyz_xzz[k] = -g_0_x_xxyyyz_xzz[k] * ab_y + g_0_x_xxyyyz_xyzz[k];

                g_0_x_xxyyyyz_yyy[k] = -g_0_x_xxyyyz_yyy[k] * ab_y + g_0_x_xxyyyz_yyyy[k];

                g_0_x_xxyyyyz_yyz[k] = -g_0_x_xxyyyz_yyz[k] * ab_y + g_0_x_xxyyyz_yyyz[k];

                g_0_x_xxyyyyz_yzz[k] = -g_0_x_xxyyyz_yzz[k] * ab_y + g_0_x_xxyyyz_yyzz[k];

                g_0_x_xxyyyyz_zzz[k] = -g_0_x_xxyyyz_zzz[k] * ab_y + g_0_x_xxyyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyzz_xxx = cbuffer.data(kf_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxy = cbuffer.data(kf_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxz = cbuffer.data(kf_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xyy = cbuffer.data(kf_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xyz = cbuffer.data(kf_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xzz = cbuffer.data(kf_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yyy = cbuffer.data(kf_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yyz = cbuffer.data(kf_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yzz = cbuffer.data(kf_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_zzz = cbuffer.data(kf_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyzz_xxx, g_0_x_xxyyyzz_xxy, g_0_x_xxyyyzz_xxz, g_0_x_xxyyyzz_xyy, g_0_x_xxyyyzz_xyz, g_0_x_xxyyyzz_xzz, g_0_x_xxyyyzz_yyy, g_0_x_xxyyyzz_yyz, g_0_x_xxyyyzz_yzz, g_0_x_xxyyyzz_zzz, g_0_x_xxyyzz_xxx, g_0_x_xxyyzz_xxxy, g_0_x_xxyyzz_xxy, g_0_x_xxyyzz_xxyy, g_0_x_xxyyzz_xxyz, g_0_x_xxyyzz_xxz, g_0_x_xxyyzz_xyy, g_0_x_xxyyzz_xyyy, g_0_x_xxyyzz_xyyz, g_0_x_xxyyzz_xyz, g_0_x_xxyyzz_xyzz, g_0_x_xxyyzz_xzz, g_0_x_xxyyzz_yyy, g_0_x_xxyyzz_yyyy, g_0_x_xxyyzz_yyyz, g_0_x_xxyyzz_yyz, g_0_x_xxyyzz_yyzz, g_0_x_xxyyzz_yzz, g_0_x_xxyyzz_yzzz, g_0_x_xxyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyzz_xxx[k] = -g_0_x_xxyyzz_xxx[k] * ab_y + g_0_x_xxyyzz_xxxy[k];

                g_0_x_xxyyyzz_xxy[k] = -g_0_x_xxyyzz_xxy[k] * ab_y + g_0_x_xxyyzz_xxyy[k];

                g_0_x_xxyyyzz_xxz[k] = -g_0_x_xxyyzz_xxz[k] * ab_y + g_0_x_xxyyzz_xxyz[k];

                g_0_x_xxyyyzz_xyy[k] = -g_0_x_xxyyzz_xyy[k] * ab_y + g_0_x_xxyyzz_xyyy[k];

                g_0_x_xxyyyzz_xyz[k] = -g_0_x_xxyyzz_xyz[k] * ab_y + g_0_x_xxyyzz_xyyz[k];

                g_0_x_xxyyyzz_xzz[k] = -g_0_x_xxyyzz_xzz[k] * ab_y + g_0_x_xxyyzz_xyzz[k];

                g_0_x_xxyyyzz_yyy[k] = -g_0_x_xxyyzz_yyy[k] * ab_y + g_0_x_xxyyzz_yyyy[k];

                g_0_x_xxyyyzz_yyz[k] = -g_0_x_xxyyzz_yyz[k] * ab_y + g_0_x_xxyyzz_yyyz[k];

                g_0_x_xxyyyzz_yzz[k] = -g_0_x_xxyyzz_yzz[k] * ab_y + g_0_x_xxyyzz_yyzz[k];

                g_0_x_xxyyyzz_zzz[k] = -g_0_x_xxyyzz_zzz[k] * ab_y + g_0_x_xxyyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzzz_xxx = cbuffer.data(kf_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxy = cbuffer.data(kf_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxz = cbuffer.data(kf_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xyy = cbuffer.data(kf_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xyz = cbuffer.data(kf_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xzz = cbuffer.data(kf_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yyy = cbuffer.data(kf_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yyz = cbuffer.data(kf_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yzz = cbuffer.data(kf_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_zzz = cbuffer.data(kf_geom_01_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzzz_xxx, g_0_x_xxyyzzz_xxy, g_0_x_xxyyzzz_xxz, g_0_x_xxyyzzz_xyy, g_0_x_xxyyzzz_xyz, g_0_x_xxyyzzz_xzz, g_0_x_xxyyzzz_yyy, g_0_x_xxyyzzz_yyz, g_0_x_xxyyzzz_yzz, g_0_x_xxyyzzz_zzz, g_0_x_xxyzzz_xxx, g_0_x_xxyzzz_xxxy, g_0_x_xxyzzz_xxy, g_0_x_xxyzzz_xxyy, g_0_x_xxyzzz_xxyz, g_0_x_xxyzzz_xxz, g_0_x_xxyzzz_xyy, g_0_x_xxyzzz_xyyy, g_0_x_xxyzzz_xyyz, g_0_x_xxyzzz_xyz, g_0_x_xxyzzz_xyzz, g_0_x_xxyzzz_xzz, g_0_x_xxyzzz_yyy, g_0_x_xxyzzz_yyyy, g_0_x_xxyzzz_yyyz, g_0_x_xxyzzz_yyz, g_0_x_xxyzzz_yyzz, g_0_x_xxyzzz_yzz, g_0_x_xxyzzz_yzzz, g_0_x_xxyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzzz_xxx[k] = -g_0_x_xxyzzz_xxx[k] * ab_y + g_0_x_xxyzzz_xxxy[k];

                g_0_x_xxyyzzz_xxy[k] = -g_0_x_xxyzzz_xxy[k] * ab_y + g_0_x_xxyzzz_xxyy[k];

                g_0_x_xxyyzzz_xxz[k] = -g_0_x_xxyzzz_xxz[k] * ab_y + g_0_x_xxyzzz_xxyz[k];

                g_0_x_xxyyzzz_xyy[k] = -g_0_x_xxyzzz_xyy[k] * ab_y + g_0_x_xxyzzz_xyyy[k];

                g_0_x_xxyyzzz_xyz[k] = -g_0_x_xxyzzz_xyz[k] * ab_y + g_0_x_xxyzzz_xyyz[k];

                g_0_x_xxyyzzz_xzz[k] = -g_0_x_xxyzzz_xzz[k] * ab_y + g_0_x_xxyzzz_xyzz[k];

                g_0_x_xxyyzzz_yyy[k] = -g_0_x_xxyzzz_yyy[k] * ab_y + g_0_x_xxyzzz_yyyy[k];

                g_0_x_xxyyzzz_yyz[k] = -g_0_x_xxyzzz_yyz[k] * ab_y + g_0_x_xxyzzz_yyyz[k];

                g_0_x_xxyyzzz_yzz[k] = -g_0_x_xxyzzz_yzz[k] * ab_y + g_0_x_xxyzzz_yyzz[k];

                g_0_x_xxyyzzz_zzz[k] = -g_0_x_xxyzzz_zzz[k] * ab_y + g_0_x_xxyzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzzz_xxx = cbuffer.data(kf_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxy = cbuffer.data(kf_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxz = cbuffer.data(kf_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xyy = cbuffer.data(kf_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xyz = cbuffer.data(kf_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xzz = cbuffer.data(kf_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yyy = cbuffer.data(kf_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yyz = cbuffer.data(kf_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yzz = cbuffer.data(kf_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_zzz = cbuffer.data(kf_geom_01_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzzz_xxx, g_0_x_xxyzzzz_xxy, g_0_x_xxyzzzz_xxz, g_0_x_xxyzzzz_xyy, g_0_x_xxyzzzz_xyz, g_0_x_xxyzzzz_xzz, g_0_x_xxyzzzz_yyy, g_0_x_xxyzzzz_yyz, g_0_x_xxyzzzz_yzz, g_0_x_xxyzzzz_zzz, g_0_x_xxzzzz_xxx, g_0_x_xxzzzz_xxxy, g_0_x_xxzzzz_xxy, g_0_x_xxzzzz_xxyy, g_0_x_xxzzzz_xxyz, g_0_x_xxzzzz_xxz, g_0_x_xxzzzz_xyy, g_0_x_xxzzzz_xyyy, g_0_x_xxzzzz_xyyz, g_0_x_xxzzzz_xyz, g_0_x_xxzzzz_xyzz, g_0_x_xxzzzz_xzz, g_0_x_xxzzzz_yyy, g_0_x_xxzzzz_yyyy, g_0_x_xxzzzz_yyyz, g_0_x_xxzzzz_yyz, g_0_x_xxzzzz_yyzz, g_0_x_xxzzzz_yzz, g_0_x_xxzzzz_yzzz, g_0_x_xxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzzz_xxx[k] = -g_0_x_xxzzzz_xxx[k] * ab_y + g_0_x_xxzzzz_xxxy[k];

                g_0_x_xxyzzzz_xxy[k] = -g_0_x_xxzzzz_xxy[k] * ab_y + g_0_x_xxzzzz_xxyy[k];

                g_0_x_xxyzzzz_xxz[k] = -g_0_x_xxzzzz_xxz[k] * ab_y + g_0_x_xxzzzz_xxyz[k];

                g_0_x_xxyzzzz_xyy[k] = -g_0_x_xxzzzz_xyy[k] * ab_y + g_0_x_xxzzzz_xyyy[k];

                g_0_x_xxyzzzz_xyz[k] = -g_0_x_xxzzzz_xyz[k] * ab_y + g_0_x_xxzzzz_xyyz[k];

                g_0_x_xxyzzzz_xzz[k] = -g_0_x_xxzzzz_xzz[k] * ab_y + g_0_x_xxzzzz_xyzz[k];

                g_0_x_xxyzzzz_yyy[k] = -g_0_x_xxzzzz_yyy[k] * ab_y + g_0_x_xxzzzz_yyyy[k];

                g_0_x_xxyzzzz_yyz[k] = -g_0_x_xxzzzz_yyz[k] * ab_y + g_0_x_xxzzzz_yyyz[k];

                g_0_x_xxyzzzz_yzz[k] = -g_0_x_xxzzzz_yzz[k] * ab_y + g_0_x_xxzzzz_yyzz[k];

                g_0_x_xxyzzzz_zzz[k] = -g_0_x_xxzzzz_zzz[k] * ab_y + g_0_x_xxzzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzzz_xxx = cbuffer.data(kf_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxy = cbuffer.data(kf_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxz = cbuffer.data(kf_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xyy = cbuffer.data(kf_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xyz = cbuffer.data(kf_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xzz = cbuffer.data(kf_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yyy = cbuffer.data(kf_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yyz = cbuffer.data(kf_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yzz = cbuffer.data(kf_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_zzz = cbuffer.data(kf_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzzz_xxx, g_0_x_xxzzzz_xxxz, g_0_x_xxzzzz_xxy, g_0_x_xxzzzz_xxyz, g_0_x_xxzzzz_xxz, g_0_x_xxzzzz_xxzz, g_0_x_xxzzzz_xyy, g_0_x_xxzzzz_xyyz, g_0_x_xxzzzz_xyz, g_0_x_xxzzzz_xyzz, g_0_x_xxzzzz_xzz, g_0_x_xxzzzz_xzzz, g_0_x_xxzzzz_yyy, g_0_x_xxzzzz_yyyz, g_0_x_xxzzzz_yyz, g_0_x_xxzzzz_yyzz, g_0_x_xxzzzz_yzz, g_0_x_xxzzzz_yzzz, g_0_x_xxzzzz_zzz, g_0_x_xxzzzz_zzzz, g_0_x_xxzzzzz_xxx, g_0_x_xxzzzzz_xxy, g_0_x_xxzzzzz_xxz, g_0_x_xxzzzzz_xyy, g_0_x_xxzzzzz_xyz, g_0_x_xxzzzzz_xzz, g_0_x_xxzzzzz_yyy, g_0_x_xxzzzzz_yyz, g_0_x_xxzzzzz_yzz, g_0_x_xxzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzzz_xxx[k] = -g_0_x_xxzzzz_xxx[k] * ab_z + g_0_x_xxzzzz_xxxz[k];

                g_0_x_xxzzzzz_xxy[k] = -g_0_x_xxzzzz_xxy[k] * ab_z + g_0_x_xxzzzz_xxyz[k];

                g_0_x_xxzzzzz_xxz[k] = -g_0_x_xxzzzz_xxz[k] * ab_z + g_0_x_xxzzzz_xxzz[k];

                g_0_x_xxzzzzz_xyy[k] = -g_0_x_xxzzzz_xyy[k] * ab_z + g_0_x_xxzzzz_xyyz[k];

                g_0_x_xxzzzzz_xyz[k] = -g_0_x_xxzzzz_xyz[k] * ab_z + g_0_x_xxzzzz_xyzz[k];

                g_0_x_xxzzzzz_xzz[k] = -g_0_x_xxzzzz_xzz[k] * ab_z + g_0_x_xxzzzz_xzzz[k];

                g_0_x_xxzzzzz_yyy[k] = -g_0_x_xxzzzz_yyy[k] * ab_z + g_0_x_xxzzzz_yyyz[k];

                g_0_x_xxzzzzz_yyz[k] = -g_0_x_xxzzzz_yyz[k] * ab_z + g_0_x_xxzzzz_yyzz[k];

                g_0_x_xxzzzzz_yzz[k] = -g_0_x_xxzzzz_yzz[k] * ab_z + g_0_x_xxzzzz_yzzz[k];

                g_0_x_xxzzzzz_zzz[k] = -g_0_x_xxzzzz_zzz[k] * ab_z + g_0_x_xxzzzz_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyy_xxx, g_0_x_xyyyyy_xxxy, g_0_x_xyyyyy_xxy, g_0_x_xyyyyy_xxyy, g_0_x_xyyyyy_xxyz, g_0_x_xyyyyy_xxz, g_0_x_xyyyyy_xyy, g_0_x_xyyyyy_xyyy, g_0_x_xyyyyy_xyyz, g_0_x_xyyyyy_xyz, g_0_x_xyyyyy_xyzz, g_0_x_xyyyyy_xzz, g_0_x_xyyyyy_yyy, g_0_x_xyyyyy_yyyy, g_0_x_xyyyyy_yyyz, g_0_x_xyyyyy_yyz, g_0_x_xyyyyy_yyzz, g_0_x_xyyyyy_yzz, g_0_x_xyyyyy_yzzz, g_0_x_xyyyyy_zzz, g_0_x_xyyyyyy_xxx, g_0_x_xyyyyyy_xxy, g_0_x_xyyyyyy_xxz, g_0_x_xyyyyyy_xyy, g_0_x_xyyyyyy_xyz, g_0_x_xyyyyyy_xzz, g_0_x_xyyyyyy_yyy, g_0_x_xyyyyyy_yyz, g_0_x_xyyyyyy_yzz, g_0_x_xyyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyy_xxx[k] = -g_0_x_xyyyyy_xxx[k] * ab_y + g_0_x_xyyyyy_xxxy[k];

                g_0_x_xyyyyyy_xxy[k] = -g_0_x_xyyyyy_xxy[k] * ab_y + g_0_x_xyyyyy_xxyy[k];

                g_0_x_xyyyyyy_xxz[k] = -g_0_x_xyyyyy_xxz[k] * ab_y + g_0_x_xyyyyy_xxyz[k];

                g_0_x_xyyyyyy_xyy[k] = -g_0_x_xyyyyy_xyy[k] * ab_y + g_0_x_xyyyyy_xyyy[k];

                g_0_x_xyyyyyy_xyz[k] = -g_0_x_xyyyyy_xyz[k] * ab_y + g_0_x_xyyyyy_xyyz[k];

                g_0_x_xyyyyyy_xzz[k] = -g_0_x_xyyyyy_xzz[k] * ab_y + g_0_x_xyyyyy_xyzz[k];

                g_0_x_xyyyyyy_yyy[k] = -g_0_x_xyyyyy_yyy[k] * ab_y + g_0_x_xyyyyy_yyyy[k];

                g_0_x_xyyyyyy_yyz[k] = -g_0_x_xyyyyy_yyz[k] * ab_y + g_0_x_xyyyyy_yyyz[k];

                g_0_x_xyyyyyy_yzz[k] = -g_0_x_xyyyyy_yzz[k] * ab_y + g_0_x_xyyyyy_yyzz[k];

                g_0_x_xyyyyyy_zzz[k] = -g_0_x_xyyyyy_zzz[k] * ab_y + g_0_x_xyyyyy_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyyz_xxx, g_0_x_xyyyyyz_xxy, g_0_x_xyyyyyz_xxz, g_0_x_xyyyyyz_xyy, g_0_x_xyyyyyz_xyz, g_0_x_xyyyyyz_xzz, g_0_x_xyyyyyz_yyy, g_0_x_xyyyyyz_yyz, g_0_x_xyyyyyz_yzz, g_0_x_xyyyyyz_zzz, g_0_x_xyyyyz_xxx, g_0_x_xyyyyz_xxxy, g_0_x_xyyyyz_xxy, g_0_x_xyyyyz_xxyy, g_0_x_xyyyyz_xxyz, g_0_x_xyyyyz_xxz, g_0_x_xyyyyz_xyy, g_0_x_xyyyyz_xyyy, g_0_x_xyyyyz_xyyz, g_0_x_xyyyyz_xyz, g_0_x_xyyyyz_xyzz, g_0_x_xyyyyz_xzz, g_0_x_xyyyyz_yyy, g_0_x_xyyyyz_yyyy, g_0_x_xyyyyz_yyyz, g_0_x_xyyyyz_yyz, g_0_x_xyyyyz_yyzz, g_0_x_xyyyyz_yzz, g_0_x_xyyyyz_yzzz, g_0_x_xyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyz_xxx[k] = -g_0_x_xyyyyz_xxx[k] * ab_y + g_0_x_xyyyyz_xxxy[k];

                g_0_x_xyyyyyz_xxy[k] = -g_0_x_xyyyyz_xxy[k] * ab_y + g_0_x_xyyyyz_xxyy[k];

                g_0_x_xyyyyyz_xxz[k] = -g_0_x_xyyyyz_xxz[k] * ab_y + g_0_x_xyyyyz_xxyz[k];

                g_0_x_xyyyyyz_xyy[k] = -g_0_x_xyyyyz_xyy[k] * ab_y + g_0_x_xyyyyz_xyyy[k];

                g_0_x_xyyyyyz_xyz[k] = -g_0_x_xyyyyz_xyz[k] * ab_y + g_0_x_xyyyyz_xyyz[k];

                g_0_x_xyyyyyz_xzz[k] = -g_0_x_xyyyyz_xzz[k] * ab_y + g_0_x_xyyyyz_xyzz[k];

                g_0_x_xyyyyyz_yyy[k] = -g_0_x_xyyyyz_yyy[k] * ab_y + g_0_x_xyyyyz_yyyy[k];

                g_0_x_xyyyyyz_yyz[k] = -g_0_x_xyyyyz_yyz[k] * ab_y + g_0_x_xyyyyz_yyyz[k];

                g_0_x_xyyyyyz_yzz[k] = -g_0_x_xyyyyz_yzz[k] * ab_y + g_0_x_xyyyyz_yyzz[k];

                g_0_x_xyyyyyz_zzz[k] = -g_0_x_xyyyyz_zzz[k] * ab_y + g_0_x_xyyyyz_yzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyzz_xxx, g_0_x_xyyyyzz_xxy, g_0_x_xyyyyzz_xxz, g_0_x_xyyyyzz_xyy, g_0_x_xyyyyzz_xyz, g_0_x_xyyyyzz_xzz, g_0_x_xyyyyzz_yyy, g_0_x_xyyyyzz_yyz, g_0_x_xyyyyzz_yzz, g_0_x_xyyyyzz_zzz, g_0_x_xyyyzz_xxx, g_0_x_xyyyzz_xxxy, g_0_x_xyyyzz_xxy, g_0_x_xyyyzz_xxyy, g_0_x_xyyyzz_xxyz, g_0_x_xyyyzz_xxz, g_0_x_xyyyzz_xyy, g_0_x_xyyyzz_xyyy, g_0_x_xyyyzz_xyyz, g_0_x_xyyyzz_xyz, g_0_x_xyyyzz_xyzz, g_0_x_xyyyzz_xzz, g_0_x_xyyyzz_yyy, g_0_x_xyyyzz_yyyy, g_0_x_xyyyzz_yyyz, g_0_x_xyyyzz_yyz, g_0_x_xyyyzz_yyzz, g_0_x_xyyyzz_yzz, g_0_x_xyyyzz_yzzz, g_0_x_xyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyzz_xxx[k] = -g_0_x_xyyyzz_xxx[k] * ab_y + g_0_x_xyyyzz_xxxy[k];

                g_0_x_xyyyyzz_xxy[k] = -g_0_x_xyyyzz_xxy[k] * ab_y + g_0_x_xyyyzz_xxyy[k];

                g_0_x_xyyyyzz_xxz[k] = -g_0_x_xyyyzz_xxz[k] * ab_y + g_0_x_xyyyzz_xxyz[k];

                g_0_x_xyyyyzz_xyy[k] = -g_0_x_xyyyzz_xyy[k] * ab_y + g_0_x_xyyyzz_xyyy[k];

                g_0_x_xyyyyzz_xyz[k] = -g_0_x_xyyyzz_xyz[k] * ab_y + g_0_x_xyyyzz_xyyz[k];

                g_0_x_xyyyyzz_xzz[k] = -g_0_x_xyyyzz_xzz[k] * ab_y + g_0_x_xyyyzz_xyzz[k];

                g_0_x_xyyyyzz_yyy[k] = -g_0_x_xyyyzz_yyy[k] * ab_y + g_0_x_xyyyzz_yyyy[k];

                g_0_x_xyyyyzz_yyz[k] = -g_0_x_xyyyzz_yyz[k] * ab_y + g_0_x_xyyyzz_yyyz[k];

                g_0_x_xyyyyzz_yzz[k] = -g_0_x_xyyyzz_yzz[k] * ab_y + g_0_x_xyyyzz_yyzz[k];

                g_0_x_xyyyyzz_zzz[k] = -g_0_x_xyyyzz_zzz[k] * ab_y + g_0_x_xyyyzz_yzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzzz_xxx, g_0_x_xyyyzzz_xxy, g_0_x_xyyyzzz_xxz, g_0_x_xyyyzzz_xyy, g_0_x_xyyyzzz_xyz, g_0_x_xyyyzzz_xzz, g_0_x_xyyyzzz_yyy, g_0_x_xyyyzzz_yyz, g_0_x_xyyyzzz_yzz, g_0_x_xyyyzzz_zzz, g_0_x_xyyzzz_xxx, g_0_x_xyyzzz_xxxy, g_0_x_xyyzzz_xxy, g_0_x_xyyzzz_xxyy, g_0_x_xyyzzz_xxyz, g_0_x_xyyzzz_xxz, g_0_x_xyyzzz_xyy, g_0_x_xyyzzz_xyyy, g_0_x_xyyzzz_xyyz, g_0_x_xyyzzz_xyz, g_0_x_xyyzzz_xyzz, g_0_x_xyyzzz_xzz, g_0_x_xyyzzz_yyy, g_0_x_xyyzzz_yyyy, g_0_x_xyyzzz_yyyz, g_0_x_xyyzzz_yyz, g_0_x_xyyzzz_yyzz, g_0_x_xyyzzz_yzz, g_0_x_xyyzzz_yzzz, g_0_x_xyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzzz_xxx[k] = -g_0_x_xyyzzz_xxx[k] * ab_y + g_0_x_xyyzzz_xxxy[k];

                g_0_x_xyyyzzz_xxy[k] = -g_0_x_xyyzzz_xxy[k] * ab_y + g_0_x_xyyzzz_xxyy[k];

                g_0_x_xyyyzzz_xxz[k] = -g_0_x_xyyzzz_xxz[k] * ab_y + g_0_x_xyyzzz_xxyz[k];

                g_0_x_xyyyzzz_xyy[k] = -g_0_x_xyyzzz_xyy[k] * ab_y + g_0_x_xyyzzz_xyyy[k];

                g_0_x_xyyyzzz_xyz[k] = -g_0_x_xyyzzz_xyz[k] * ab_y + g_0_x_xyyzzz_xyyz[k];

                g_0_x_xyyyzzz_xzz[k] = -g_0_x_xyyzzz_xzz[k] * ab_y + g_0_x_xyyzzz_xyzz[k];

                g_0_x_xyyyzzz_yyy[k] = -g_0_x_xyyzzz_yyy[k] * ab_y + g_0_x_xyyzzz_yyyy[k];

                g_0_x_xyyyzzz_yyz[k] = -g_0_x_xyyzzz_yyz[k] * ab_y + g_0_x_xyyzzz_yyyz[k];

                g_0_x_xyyyzzz_yzz[k] = -g_0_x_xyyzzz_yzz[k] * ab_y + g_0_x_xyyzzz_yyzz[k];

                g_0_x_xyyyzzz_zzz[k] = -g_0_x_xyyzzz_zzz[k] * ab_y + g_0_x_xyyzzz_yzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzzz_xxx, g_0_x_xyyzzzz_xxy, g_0_x_xyyzzzz_xxz, g_0_x_xyyzzzz_xyy, g_0_x_xyyzzzz_xyz, g_0_x_xyyzzzz_xzz, g_0_x_xyyzzzz_yyy, g_0_x_xyyzzzz_yyz, g_0_x_xyyzzzz_yzz, g_0_x_xyyzzzz_zzz, g_0_x_xyzzzz_xxx, g_0_x_xyzzzz_xxxy, g_0_x_xyzzzz_xxy, g_0_x_xyzzzz_xxyy, g_0_x_xyzzzz_xxyz, g_0_x_xyzzzz_xxz, g_0_x_xyzzzz_xyy, g_0_x_xyzzzz_xyyy, g_0_x_xyzzzz_xyyz, g_0_x_xyzzzz_xyz, g_0_x_xyzzzz_xyzz, g_0_x_xyzzzz_xzz, g_0_x_xyzzzz_yyy, g_0_x_xyzzzz_yyyy, g_0_x_xyzzzz_yyyz, g_0_x_xyzzzz_yyz, g_0_x_xyzzzz_yyzz, g_0_x_xyzzzz_yzz, g_0_x_xyzzzz_yzzz, g_0_x_xyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzzz_xxx[k] = -g_0_x_xyzzzz_xxx[k] * ab_y + g_0_x_xyzzzz_xxxy[k];

                g_0_x_xyyzzzz_xxy[k] = -g_0_x_xyzzzz_xxy[k] * ab_y + g_0_x_xyzzzz_xxyy[k];

                g_0_x_xyyzzzz_xxz[k] = -g_0_x_xyzzzz_xxz[k] * ab_y + g_0_x_xyzzzz_xxyz[k];

                g_0_x_xyyzzzz_xyy[k] = -g_0_x_xyzzzz_xyy[k] * ab_y + g_0_x_xyzzzz_xyyy[k];

                g_0_x_xyyzzzz_xyz[k] = -g_0_x_xyzzzz_xyz[k] * ab_y + g_0_x_xyzzzz_xyyz[k];

                g_0_x_xyyzzzz_xzz[k] = -g_0_x_xyzzzz_xzz[k] * ab_y + g_0_x_xyzzzz_xyzz[k];

                g_0_x_xyyzzzz_yyy[k] = -g_0_x_xyzzzz_yyy[k] * ab_y + g_0_x_xyzzzz_yyyy[k];

                g_0_x_xyyzzzz_yyz[k] = -g_0_x_xyzzzz_yyz[k] * ab_y + g_0_x_xyzzzz_yyyz[k];

                g_0_x_xyyzzzz_yzz[k] = -g_0_x_xyzzzz_yzz[k] * ab_y + g_0_x_xyzzzz_yyzz[k];

                g_0_x_xyyzzzz_zzz[k] = -g_0_x_xyzzzz_zzz[k] * ab_y + g_0_x_xyzzzz_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzzz_xxx, g_0_x_xyzzzzz_xxy, g_0_x_xyzzzzz_xxz, g_0_x_xyzzzzz_xyy, g_0_x_xyzzzzz_xyz, g_0_x_xyzzzzz_xzz, g_0_x_xyzzzzz_yyy, g_0_x_xyzzzzz_yyz, g_0_x_xyzzzzz_yzz, g_0_x_xyzzzzz_zzz, g_0_x_xzzzzz_xxx, g_0_x_xzzzzz_xxxy, g_0_x_xzzzzz_xxy, g_0_x_xzzzzz_xxyy, g_0_x_xzzzzz_xxyz, g_0_x_xzzzzz_xxz, g_0_x_xzzzzz_xyy, g_0_x_xzzzzz_xyyy, g_0_x_xzzzzz_xyyz, g_0_x_xzzzzz_xyz, g_0_x_xzzzzz_xyzz, g_0_x_xzzzzz_xzz, g_0_x_xzzzzz_yyy, g_0_x_xzzzzz_yyyy, g_0_x_xzzzzz_yyyz, g_0_x_xzzzzz_yyz, g_0_x_xzzzzz_yyzz, g_0_x_xzzzzz_yzz, g_0_x_xzzzzz_yzzz, g_0_x_xzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzzz_xxx[k] = -g_0_x_xzzzzz_xxx[k] * ab_y + g_0_x_xzzzzz_xxxy[k];

                g_0_x_xyzzzzz_xxy[k] = -g_0_x_xzzzzz_xxy[k] * ab_y + g_0_x_xzzzzz_xxyy[k];

                g_0_x_xyzzzzz_xxz[k] = -g_0_x_xzzzzz_xxz[k] * ab_y + g_0_x_xzzzzz_xxyz[k];

                g_0_x_xyzzzzz_xyy[k] = -g_0_x_xzzzzz_xyy[k] * ab_y + g_0_x_xzzzzz_xyyy[k];

                g_0_x_xyzzzzz_xyz[k] = -g_0_x_xzzzzz_xyz[k] * ab_y + g_0_x_xzzzzz_xyyz[k];

                g_0_x_xyzzzzz_xzz[k] = -g_0_x_xzzzzz_xzz[k] * ab_y + g_0_x_xzzzzz_xyzz[k];

                g_0_x_xyzzzzz_yyy[k] = -g_0_x_xzzzzz_yyy[k] * ab_y + g_0_x_xzzzzz_yyyy[k];

                g_0_x_xyzzzzz_yyz[k] = -g_0_x_xzzzzz_yyz[k] * ab_y + g_0_x_xzzzzz_yyyz[k];

                g_0_x_xyzzzzz_yzz[k] = -g_0_x_xzzzzz_yzz[k] * ab_y + g_0_x_xzzzzz_yyzz[k];

                g_0_x_xyzzzzz_zzz[k] = -g_0_x_xzzzzz_zzz[k] * ab_y + g_0_x_xzzzzz_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzzz_xxx, g_0_x_xzzzzz_xxxz, g_0_x_xzzzzz_xxy, g_0_x_xzzzzz_xxyz, g_0_x_xzzzzz_xxz, g_0_x_xzzzzz_xxzz, g_0_x_xzzzzz_xyy, g_0_x_xzzzzz_xyyz, g_0_x_xzzzzz_xyz, g_0_x_xzzzzz_xyzz, g_0_x_xzzzzz_xzz, g_0_x_xzzzzz_xzzz, g_0_x_xzzzzz_yyy, g_0_x_xzzzzz_yyyz, g_0_x_xzzzzz_yyz, g_0_x_xzzzzz_yyzz, g_0_x_xzzzzz_yzz, g_0_x_xzzzzz_yzzz, g_0_x_xzzzzz_zzz, g_0_x_xzzzzz_zzzz, g_0_x_xzzzzzz_xxx, g_0_x_xzzzzzz_xxy, g_0_x_xzzzzzz_xxz, g_0_x_xzzzzzz_xyy, g_0_x_xzzzzzz_xyz, g_0_x_xzzzzzz_xzz, g_0_x_xzzzzzz_yyy, g_0_x_xzzzzzz_yyz, g_0_x_xzzzzzz_yzz, g_0_x_xzzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzzz_xxx[k] = -g_0_x_xzzzzz_xxx[k] * ab_z + g_0_x_xzzzzz_xxxz[k];

                g_0_x_xzzzzzz_xxy[k] = -g_0_x_xzzzzz_xxy[k] * ab_z + g_0_x_xzzzzz_xxyz[k];

                g_0_x_xzzzzzz_xxz[k] = -g_0_x_xzzzzz_xxz[k] * ab_z + g_0_x_xzzzzz_xxzz[k];

                g_0_x_xzzzzzz_xyy[k] = -g_0_x_xzzzzz_xyy[k] * ab_z + g_0_x_xzzzzz_xyyz[k];

                g_0_x_xzzzzzz_xyz[k] = -g_0_x_xzzzzz_xyz[k] * ab_z + g_0_x_xzzzzz_xyzz[k];

                g_0_x_xzzzzzz_xzz[k] = -g_0_x_xzzzzz_xzz[k] * ab_z + g_0_x_xzzzzz_xzzz[k];

                g_0_x_xzzzzzz_yyy[k] = -g_0_x_xzzzzz_yyy[k] * ab_z + g_0_x_xzzzzz_yyyz[k];

                g_0_x_xzzzzzz_yyz[k] = -g_0_x_xzzzzz_yyz[k] * ab_z + g_0_x_xzzzzz_yyzz[k];

                g_0_x_xzzzzzz_yzz[k] = -g_0_x_xzzzzz_yzz[k] * ab_z + g_0_x_xzzzzz_yzzz[k];

                g_0_x_xzzzzzz_zzz[k] = -g_0_x_xzzzzz_zzz[k] * ab_z + g_0_x_xzzzzz_zzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyy_xxx, g_0_x_yyyyyy_xxxy, g_0_x_yyyyyy_xxy, g_0_x_yyyyyy_xxyy, g_0_x_yyyyyy_xxyz, g_0_x_yyyyyy_xxz, g_0_x_yyyyyy_xyy, g_0_x_yyyyyy_xyyy, g_0_x_yyyyyy_xyyz, g_0_x_yyyyyy_xyz, g_0_x_yyyyyy_xyzz, g_0_x_yyyyyy_xzz, g_0_x_yyyyyy_yyy, g_0_x_yyyyyy_yyyy, g_0_x_yyyyyy_yyyz, g_0_x_yyyyyy_yyz, g_0_x_yyyyyy_yyzz, g_0_x_yyyyyy_yzz, g_0_x_yyyyyy_yzzz, g_0_x_yyyyyy_zzz, g_0_x_yyyyyyy_xxx, g_0_x_yyyyyyy_xxy, g_0_x_yyyyyyy_xxz, g_0_x_yyyyyyy_xyy, g_0_x_yyyyyyy_xyz, g_0_x_yyyyyyy_xzz, g_0_x_yyyyyyy_yyy, g_0_x_yyyyyyy_yyz, g_0_x_yyyyyyy_yzz, g_0_x_yyyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyy_xxx[k] = -g_0_x_yyyyyy_xxx[k] * ab_y + g_0_x_yyyyyy_xxxy[k];

                g_0_x_yyyyyyy_xxy[k] = -g_0_x_yyyyyy_xxy[k] * ab_y + g_0_x_yyyyyy_xxyy[k];

                g_0_x_yyyyyyy_xxz[k] = -g_0_x_yyyyyy_xxz[k] * ab_y + g_0_x_yyyyyy_xxyz[k];

                g_0_x_yyyyyyy_xyy[k] = -g_0_x_yyyyyy_xyy[k] * ab_y + g_0_x_yyyyyy_xyyy[k];

                g_0_x_yyyyyyy_xyz[k] = -g_0_x_yyyyyy_xyz[k] * ab_y + g_0_x_yyyyyy_xyyz[k];

                g_0_x_yyyyyyy_xzz[k] = -g_0_x_yyyyyy_xzz[k] * ab_y + g_0_x_yyyyyy_xyzz[k];

                g_0_x_yyyyyyy_yyy[k] = -g_0_x_yyyyyy_yyy[k] * ab_y + g_0_x_yyyyyy_yyyy[k];

                g_0_x_yyyyyyy_yyz[k] = -g_0_x_yyyyyy_yyz[k] * ab_y + g_0_x_yyyyyy_yyyz[k];

                g_0_x_yyyyyyy_yzz[k] = -g_0_x_yyyyyy_yzz[k] * ab_y + g_0_x_yyyyyy_yyzz[k];

                g_0_x_yyyyyyy_zzz[k] = -g_0_x_yyyyyy_zzz[k] * ab_y + g_0_x_yyyyyy_yzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyyz_xxx, g_0_x_yyyyyyz_xxy, g_0_x_yyyyyyz_xxz, g_0_x_yyyyyyz_xyy, g_0_x_yyyyyyz_xyz, g_0_x_yyyyyyz_xzz, g_0_x_yyyyyyz_yyy, g_0_x_yyyyyyz_yyz, g_0_x_yyyyyyz_yzz, g_0_x_yyyyyyz_zzz, g_0_x_yyyyyz_xxx, g_0_x_yyyyyz_xxxy, g_0_x_yyyyyz_xxy, g_0_x_yyyyyz_xxyy, g_0_x_yyyyyz_xxyz, g_0_x_yyyyyz_xxz, g_0_x_yyyyyz_xyy, g_0_x_yyyyyz_xyyy, g_0_x_yyyyyz_xyyz, g_0_x_yyyyyz_xyz, g_0_x_yyyyyz_xyzz, g_0_x_yyyyyz_xzz, g_0_x_yyyyyz_yyy, g_0_x_yyyyyz_yyyy, g_0_x_yyyyyz_yyyz, g_0_x_yyyyyz_yyz, g_0_x_yyyyyz_yyzz, g_0_x_yyyyyz_yzz, g_0_x_yyyyyz_yzzz, g_0_x_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyz_xxx[k] = -g_0_x_yyyyyz_xxx[k] * ab_y + g_0_x_yyyyyz_xxxy[k];

                g_0_x_yyyyyyz_xxy[k] = -g_0_x_yyyyyz_xxy[k] * ab_y + g_0_x_yyyyyz_xxyy[k];

                g_0_x_yyyyyyz_xxz[k] = -g_0_x_yyyyyz_xxz[k] * ab_y + g_0_x_yyyyyz_xxyz[k];

                g_0_x_yyyyyyz_xyy[k] = -g_0_x_yyyyyz_xyy[k] * ab_y + g_0_x_yyyyyz_xyyy[k];

                g_0_x_yyyyyyz_xyz[k] = -g_0_x_yyyyyz_xyz[k] * ab_y + g_0_x_yyyyyz_xyyz[k];

                g_0_x_yyyyyyz_xzz[k] = -g_0_x_yyyyyz_xzz[k] * ab_y + g_0_x_yyyyyz_xyzz[k];

                g_0_x_yyyyyyz_yyy[k] = -g_0_x_yyyyyz_yyy[k] * ab_y + g_0_x_yyyyyz_yyyy[k];

                g_0_x_yyyyyyz_yyz[k] = -g_0_x_yyyyyz_yyz[k] * ab_y + g_0_x_yyyyyz_yyyz[k];

                g_0_x_yyyyyyz_yzz[k] = -g_0_x_yyyyyz_yzz[k] * ab_y + g_0_x_yyyyyz_yyzz[k];

                g_0_x_yyyyyyz_zzz[k] = -g_0_x_yyyyyz_zzz[k] * ab_y + g_0_x_yyyyyz_yzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyzz_xxx, g_0_x_yyyyyzz_xxy, g_0_x_yyyyyzz_xxz, g_0_x_yyyyyzz_xyy, g_0_x_yyyyyzz_xyz, g_0_x_yyyyyzz_xzz, g_0_x_yyyyyzz_yyy, g_0_x_yyyyyzz_yyz, g_0_x_yyyyyzz_yzz, g_0_x_yyyyyzz_zzz, g_0_x_yyyyzz_xxx, g_0_x_yyyyzz_xxxy, g_0_x_yyyyzz_xxy, g_0_x_yyyyzz_xxyy, g_0_x_yyyyzz_xxyz, g_0_x_yyyyzz_xxz, g_0_x_yyyyzz_xyy, g_0_x_yyyyzz_xyyy, g_0_x_yyyyzz_xyyz, g_0_x_yyyyzz_xyz, g_0_x_yyyyzz_xyzz, g_0_x_yyyyzz_xzz, g_0_x_yyyyzz_yyy, g_0_x_yyyyzz_yyyy, g_0_x_yyyyzz_yyyz, g_0_x_yyyyzz_yyz, g_0_x_yyyyzz_yyzz, g_0_x_yyyyzz_yzz, g_0_x_yyyyzz_yzzz, g_0_x_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyzz_xxx[k] = -g_0_x_yyyyzz_xxx[k] * ab_y + g_0_x_yyyyzz_xxxy[k];

                g_0_x_yyyyyzz_xxy[k] = -g_0_x_yyyyzz_xxy[k] * ab_y + g_0_x_yyyyzz_xxyy[k];

                g_0_x_yyyyyzz_xxz[k] = -g_0_x_yyyyzz_xxz[k] * ab_y + g_0_x_yyyyzz_xxyz[k];

                g_0_x_yyyyyzz_xyy[k] = -g_0_x_yyyyzz_xyy[k] * ab_y + g_0_x_yyyyzz_xyyy[k];

                g_0_x_yyyyyzz_xyz[k] = -g_0_x_yyyyzz_xyz[k] * ab_y + g_0_x_yyyyzz_xyyz[k];

                g_0_x_yyyyyzz_xzz[k] = -g_0_x_yyyyzz_xzz[k] * ab_y + g_0_x_yyyyzz_xyzz[k];

                g_0_x_yyyyyzz_yyy[k] = -g_0_x_yyyyzz_yyy[k] * ab_y + g_0_x_yyyyzz_yyyy[k];

                g_0_x_yyyyyzz_yyz[k] = -g_0_x_yyyyzz_yyz[k] * ab_y + g_0_x_yyyyzz_yyyz[k];

                g_0_x_yyyyyzz_yzz[k] = -g_0_x_yyyyzz_yzz[k] * ab_y + g_0_x_yyyyzz_yyzz[k];

                g_0_x_yyyyyzz_zzz[k] = -g_0_x_yyyyzz_zzz[k] * ab_y + g_0_x_yyyyzz_yzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzzz_xxx, g_0_x_yyyyzzz_xxy, g_0_x_yyyyzzz_xxz, g_0_x_yyyyzzz_xyy, g_0_x_yyyyzzz_xyz, g_0_x_yyyyzzz_xzz, g_0_x_yyyyzzz_yyy, g_0_x_yyyyzzz_yyz, g_0_x_yyyyzzz_yzz, g_0_x_yyyyzzz_zzz, g_0_x_yyyzzz_xxx, g_0_x_yyyzzz_xxxy, g_0_x_yyyzzz_xxy, g_0_x_yyyzzz_xxyy, g_0_x_yyyzzz_xxyz, g_0_x_yyyzzz_xxz, g_0_x_yyyzzz_xyy, g_0_x_yyyzzz_xyyy, g_0_x_yyyzzz_xyyz, g_0_x_yyyzzz_xyz, g_0_x_yyyzzz_xyzz, g_0_x_yyyzzz_xzz, g_0_x_yyyzzz_yyy, g_0_x_yyyzzz_yyyy, g_0_x_yyyzzz_yyyz, g_0_x_yyyzzz_yyz, g_0_x_yyyzzz_yyzz, g_0_x_yyyzzz_yzz, g_0_x_yyyzzz_yzzz, g_0_x_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzzz_xxx[k] = -g_0_x_yyyzzz_xxx[k] * ab_y + g_0_x_yyyzzz_xxxy[k];

                g_0_x_yyyyzzz_xxy[k] = -g_0_x_yyyzzz_xxy[k] * ab_y + g_0_x_yyyzzz_xxyy[k];

                g_0_x_yyyyzzz_xxz[k] = -g_0_x_yyyzzz_xxz[k] * ab_y + g_0_x_yyyzzz_xxyz[k];

                g_0_x_yyyyzzz_xyy[k] = -g_0_x_yyyzzz_xyy[k] * ab_y + g_0_x_yyyzzz_xyyy[k];

                g_0_x_yyyyzzz_xyz[k] = -g_0_x_yyyzzz_xyz[k] * ab_y + g_0_x_yyyzzz_xyyz[k];

                g_0_x_yyyyzzz_xzz[k] = -g_0_x_yyyzzz_xzz[k] * ab_y + g_0_x_yyyzzz_xyzz[k];

                g_0_x_yyyyzzz_yyy[k] = -g_0_x_yyyzzz_yyy[k] * ab_y + g_0_x_yyyzzz_yyyy[k];

                g_0_x_yyyyzzz_yyz[k] = -g_0_x_yyyzzz_yyz[k] * ab_y + g_0_x_yyyzzz_yyyz[k];

                g_0_x_yyyyzzz_yzz[k] = -g_0_x_yyyzzz_yzz[k] * ab_y + g_0_x_yyyzzz_yyzz[k];

                g_0_x_yyyyzzz_zzz[k] = -g_0_x_yyyzzz_zzz[k] * ab_y + g_0_x_yyyzzz_yzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzzz_xxx, g_0_x_yyyzzzz_xxy, g_0_x_yyyzzzz_xxz, g_0_x_yyyzzzz_xyy, g_0_x_yyyzzzz_xyz, g_0_x_yyyzzzz_xzz, g_0_x_yyyzzzz_yyy, g_0_x_yyyzzzz_yyz, g_0_x_yyyzzzz_yzz, g_0_x_yyyzzzz_zzz, g_0_x_yyzzzz_xxx, g_0_x_yyzzzz_xxxy, g_0_x_yyzzzz_xxy, g_0_x_yyzzzz_xxyy, g_0_x_yyzzzz_xxyz, g_0_x_yyzzzz_xxz, g_0_x_yyzzzz_xyy, g_0_x_yyzzzz_xyyy, g_0_x_yyzzzz_xyyz, g_0_x_yyzzzz_xyz, g_0_x_yyzzzz_xyzz, g_0_x_yyzzzz_xzz, g_0_x_yyzzzz_yyy, g_0_x_yyzzzz_yyyy, g_0_x_yyzzzz_yyyz, g_0_x_yyzzzz_yyz, g_0_x_yyzzzz_yyzz, g_0_x_yyzzzz_yzz, g_0_x_yyzzzz_yzzz, g_0_x_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzzz_xxx[k] = -g_0_x_yyzzzz_xxx[k] * ab_y + g_0_x_yyzzzz_xxxy[k];

                g_0_x_yyyzzzz_xxy[k] = -g_0_x_yyzzzz_xxy[k] * ab_y + g_0_x_yyzzzz_xxyy[k];

                g_0_x_yyyzzzz_xxz[k] = -g_0_x_yyzzzz_xxz[k] * ab_y + g_0_x_yyzzzz_xxyz[k];

                g_0_x_yyyzzzz_xyy[k] = -g_0_x_yyzzzz_xyy[k] * ab_y + g_0_x_yyzzzz_xyyy[k];

                g_0_x_yyyzzzz_xyz[k] = -g_0_x_yyzzzz_xyz[k] * ab_y + g_0_x_yyzzzz_xyyz[k];

                g_0_x_yyyzzzz_xzz[k] = -g_0_x_yyzzzz_xzz[k] * ab_y + g_0_x_yyzzzz_xyzz[k];

                g_0_x_yyyzzzz_yyy[k] = -g_0_x_yyzzzz_yyy[k] * ab_y + g_0_x_yyzzzz_yyyy[k];

                g_0_x_yyyzzzz_yyz[k] = -g_0_x_yyzzzz_yyz[k] * ab_y + g_0_x_yyzzzz_yyyz[k];

                g_0_x_yyyzzzz_yzz[k] = -g_0_x_yyzzzz_yzz[k] * ab_y + g_0_x_yyzzzz_yyzz[k];

                g_0_x_yyyzzzz_zzz[k] = -g_0_x_yyzzzz_zzz[k] * ab_y + g_0_x_yyzzzz_yzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzzz_xxx, g_0_x_yyzzzzz_xxy, g_0_x_yyzzzzz_xxz, g_0_x_yyzzzzz_xyy, g_0_x_yyzzzzz_xyz, g_0_x_yyzzzzz_xzz, g_0_x_yyzzzzz_yyy, g_0_x_yyzzzzz_yyz, g_0_x_yyzzzzz_yzz, g_0_x_yyzzzzz_zzz, g_0_x_yzzzzz_xxx, g_0_x_yzzzzz_xxxy, g_0_x_yzzzzz_xxy, g_0_x_yzzzzz_xxyy, g_0_x_yzzzzz_xxyz, g_0_x_yzzzzz_xxz, g_0_x_yzzzzz_xyy, g_0_x_yzzzzz_xyyy, g_0_x_yzzzzz_xyyz, g_0_x_yzzzzz_xyz, g_0_x_yzzzzz_xyzz, g_0_x_yzzzzz_xzz, g_0_x_yzzzzz_yyy, g_0_x_yzzzzz_yyyy, g_0_x_yzzzzz_yyyz, g_0_x_yzzzzz_yyz, g_0_x_yzzzzz_yyzz, g_0_x_yzzzzz_yzz, g_0_x_yzzzzz_yzzz, g_0_x_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzzz_xxx[k] = -g_0_x_yzzzzz_xxx[k] * ab_y + g_0_x_yzzzzz_xxxy[k];

                g_0_x_yyzzzzz_xxy[k] = -g_0_x_yzzzzz_xxy[k] * ab_y + g_0_x_yzzzzz_xxyy[k];

                g_0_x_yyzzzzz_xxz[k] = -g_0_x_yzzzzz_xxz[k] * ab_y + g_0_x_yzzzzz_xxyz[k];

                g_0_x_yyzzzzz_xyy[k] = -g_0_x_yzzzzz_xyy[k] * ab_y + g_0_x_yzzzzz_xyyy[k];

                g_0_x_yyzzzzz_xyz[k] = -g_0_x_yzzzzz_xyz[k] * ab_y + g_0_x_yzzzzz_xyyz[k];

                g_0_x_yyzzzzz_xzz[k] = -g_0_x_yzzzzz_xzz[k] * ab_y + g_0_x_yzzzzz_xyzz[k];

                g_0_x_yyzzzzz_yyy[k] = -g_0_x_yzzzzz_yyy[k] * ab_y + g_0_x_yzzzzz_yyyy[k];

                g_0_x_yyzzzzz_yyz[k] = -g_0_x_yzzzzz_yyz[k] * ab_y + g_0_x_yzzzzz_yyyz[k];

                g_0_x_yyzzzzz_yzz[k] = -g_0_x_yzzzzz_yzz[k] * ab_y + g_0_x_yzzzzz_yyzz[k];

                g_0_x_yyzzzzz_zzz[k] = -g_0_x_yzzzzz_zzz[k] * ab_y + g_0_x_yzzzzz_yzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzzz_xxx, g_0_x_yzzzzzz_xxy, g_0_x_yzzzzzz_xxz, g_0_x_yzzzzzz_xyy, g_0_x_yzzzzzz_xyz, g_0_x_yzzzzzz_xzz, g_0_x_yzzzzzz_yyy, g_0_x_yzzzzzz_yyz, g_0_x_yzzzzzz_yzz, g_0_x_yzzzzzz_zzz, g_0_x_zzzzzz_xxx, g_0_x_zzzzzz_xxxy, g_0_x_zzzzzz_xxy, g_0_x_zzzzzz_xxyy, g_0_x_zzzzzz_xxyz, g_0_x_zzzzzz_xxz, g_0_x_zzzzzz_xyy, g_0_x_zzzzzz_xyyy, g_0_x_zzzzzz_xyyz, g_0_x_zzzzzz_xyz, g_0_x_zzzzzz_xyzz, g_0_x_zzzzzz_xzz, g_0_x_zzzzzz_yyy, g_0_x_zzzzzz_yyyy, g_0_x_zzzzzz_yyyz, g_0_x_zzzzzz_yyz, g_0_x_zzzzzz_yyzz, g_0_x_zzzzzz_yzz, g_0_x_zzzzzz_yzzz, g_0_x_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzzz_xxx[k] = -g_0_x_zzzzzz_xxx[k] * ab_y + g_0_x_zzzzzz_xxxy[k];

                g_0_x_yzzzzzz_xxy[k] = -g_0_x_zzzzzz_xxy[k] * ab_y + g_0_x_zzzzzz_xxyy[k];

                g_0_x_yzzzzzz_xxz[k] = -g_0_x_zzzzzz_xxz[k] * ab_y + g_0_x_zzzzzz_xxyz[k];

                g_0_x_yzzzzzz_xyy[k] = -g_0_x_zzzzzz_xyy[k] * ab_y + g_0_x_zzzzzz_xyyy[k];

                g_0_x_yzzzzzz_xyz[k] = -g_0_x_zzzzzz_xyz[k] * ab_y + g_0_x_zzzzzz_xyyz[k];

                g_0_x_yzzzzzz_xzz[k] = -g_0_x_zzzzzz_xzz[k] * ab_y + g_0_x_zzzzzz_xyzz[k];

                g_0_x_yzzzzzz_yyy[k] = -g_0_x_zzzzzz_yyy[k] * ab_y + g_0_x_zzzzzz_yyyy[k];

                g_0_x_yzzzzzz_yyz[k] = -g_0_x_zzzzzz_yyz[k] * ab_y + g_0_x_zzzzzz_yyyz[k];

                g_0_x_yzzzzzz_yzz[k] = -g_0_x_zzzzzz_yzz[k] * ab_y + g_0_x_zzzzzz_yyzz[k];

                g_0_x_yzzzzzz_zzz[k] = -g_0_x_zzzzzz_zzz[k] * ab_y + g_0_x_zzzzzz_yzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzzz_xxx, g_0_x_zzzzzz_xxxz, g_0_x_zzzzzz_xxy, g_0_x_zzzzzz_xxyz, g_0_x_zzzzzz_xxz, g_0_x_zzzzzz_xxzz, g_0_x_zzzzzz_xyy, g_0_x_zzzzzz_xyyz, g_0_x_zzzzzz_xyz, g_0_x_zzzzzz_xyzz, g_0_x_zzzzzz_xzz, g_0_x_zzzzzz_xzzz, g_0_x_zzzzzz_yyy, g_0_x_zzzzzz_yyyz, g_0_x_zzzzzz_yyz, g_0_x_zzzzzz_yyzz, g_0_x_zzzzzz_yzz, g_0_x_zzzzzz_yzzz, g_0_x_zzzzzz_zzz, g_0_x_zzzzzz_zzzz, g_0_x_zzzzzzz_xxx, g_0_x_zzzzzzz_xxy, g_0_x_zzzzzzz_xxz, g_0_x_zzzzzzz_xyy, g_0_x_zzzzzzz_xyz, g_0_x_zzzzzzz_xzz, g_0_x_zzzzzzz_yyy, g_0_x_zzzzzzz_yyz, g_0_x_zzzzzzz_yzz, g_0_x_zzzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzzz_xxx[k] = -g_0_x_zzzzzz_xxx[k] * ab_z + g_0_x_zzzzzz_xxxz[k];

                g_0_x_zzzzzzz_xxy[k] = -g_0_x_zzzzzz_xxy[k] * ab_z + g_0_x_zzzzzz_xxyz[k];

                g_0_x_zzzzzzz_xxz[k] = -g_0_x_zzzzzz_xxz[k] * ab_z + g_0_x_zzzzzz_xxzz[k];

                g_0_x_zzzzzzz_xyy[k] = -g_0_x_zzzzzz_xyy[k] * ab_z + g_0_x_zzzzzz_xyyz[k];

                g_0_x_zzzzzzz_xyz[k] = -g_0_x_zzzzzz_xyz[k] * ab_z + g_0_x_zzzzzz_xyzz[k];

                g_0_x_zzzzzzz_xzz[k] = -g_0_x_zzzzzz_xzz[k] * ab_z + g_0_x_zzzzzz_xzzz[k];

                g_0_x_zzzzzzz_yyy[k] = -g_0_x_zzzzzz_yyy[k] * ab_z + g_0_x_zzzzzz_yyyz[k];

                g_0_x_zzzzzzz_yyz[k] = -g_0_x_zzzzzz_yyz[k] * ab_z + g_0_x_zzzzzz_yyzz[k];

                g_0_x_zzzzzzz_yzz[k] = -g_0_x_zzzzzz_yzz[k] * ab_z + g_0_x_zzzzzz_yzzz[k];

                g_0_x_zzzzzzz_zzz[k] = -g_0_x_zzzzzz_zzz[k] * ab_z + g_0_x_zzzzzz_zzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxx_xxx = cbuffer.data(kf_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxy = cbuffer.data(kf_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxz = cbuffer.data(kf_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xyy = cbuffer.data(kf_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xyz = cbuffer.data(kf_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xzz = cbuffer.data(kf_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yyy = cbuffer.data(kf_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yyz = cbuffer.data(kf_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yzz = cbuffer.data(kf_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_zzz = cbuffer.data(kf_geom_01_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxx_xxx, g_0_y_xxxxxx_xxxx, g_0_y_xxxxxx_xxxy, g_0_y_xxxxxx_xxxz, g_0_y_xxxxxx_xxy, g_0_y_xxxxxx_xxyy, g_0_y_xxxxxx_xxyz, g_0_y_xxxxxx_xxz, g_0_y_xxxxxx_xxzz, g_0_y_xxxxxx_xyy, g_0_y_xxxxxx_xyyy, g_0_y_xxxxxx_xyyz, g_0_y_xxxxxx_xyz, g_0_y_xxxxxx_xyzz, g_0_y_xxxxxx_xzz, g_0_y_xxxxxx_xzzz, g_0_y_xxxxxx_yyy, g_0_y_xxxxxx_yyz, g_0_y_xxxxxx_yzz, g_0_y_xxxxxx_zzz, g_0_y_xxxxxxx_xxx, g_0_y_xxxxxxx_xxy, g_0_y_xxxxxxx_xxz, g_0_y_xxxxxxx_xyy, g_0_y_xxxxxxx_xyz, g_0_y_xxxxxxx_xzz, g_0_y_xxxxxxx_yyy, g_0_y_xxxxxxx_yyz, g_0_y_xxxxxxx_yzz, g_0_y_xxxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxx_xxx[k] = -g_0_y_xxxxxx_xxx[k] * ab_x + g_0_y_xxxxxx_xxxx[k];

                g_0_y_xxxxxxx_xxy[k] = -g_0_y_xxxxxx_xxy[k] * ab_x + g_0_y_xxxxxx_xxxy[k];

                g_0_y_xxxxxxx_xxz[k] = -g_0_y_xxxxxx_xxz[k] * ab_x + g_0_y_xxxxxx_xxxz[k];

                g_0_y_xxxxxxx_xyy[k] = -g_0_y_xxxxxx_xyy[k] * ab_x + g_0_y_xxxxxx_xxyy[k];

                g_0_y_xxxxxxx_xyz[k] = -g_0_y_xxxxxx_xyz[k] * ab_x + g_0_y_xxxxxx_xxyz[k];

                g_0_y_xxxxxxx_xzz[k] = -g_0_y_xxxxxx_xzz[k] * ab_x + g_0_y_xxxxxx_xxzz[k];

                g_0_y_xxxxxxx_yyy[k] = -g_0_y_xxxxxx_yyy[k] * ab_x + g_0_y_xxxxxx_xyyy[k];

                g_0_y_xxxxxxx_yyz[k] = -g_0_y_xxxxxx_yyz[k] * ab_x + g_0_y_xxxxxx_xyyz[k];

                g_0_y_xxxxxxx_yzz[k] = -g_0_y_xxxxxx_yzz[k] * ab_x + g_0_y_xxxxxx_xyzz[k];

                g_0_y_xxxxxxx_zzz[k] = -g_0_y_xxxxxx_zzz[k] * ab_x + g_0_y_xxxxxx_xzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxy_xxx = cbuffer.data(kf_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxy = cbuffer.data(kf_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxz = cbuffer.data(kf_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xyy = cbuffer.data(kf_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xyz = cbuffer.data(kf_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xzz = cbuffer.data(kf_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yyy = cbuffer.data(kf_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yyz = cbuffer.data(kf_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yzz = cbuffer.data(kf_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_zzz = cbuffer.data(kf_geom_01_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxy_xxx, g_0_y_xxxxxxy_xxy, g_0_y_xxxxxxy_xxz, g_0_y_xxxxxxy_xyy, g_0_y_xxxxxxy_xyz, g_0_y_xxxxxxy_xzz, g_0_y_xxxxxxy_yyy, g_0_y_xxxxxxy_yyz, g_0_y_xxxxxxy_yzz, g_0_y_xxxxxxy_zzz, g_0_y_xxxxxy_xxx, g_0_y_xxxxxy_xxxx, g_0_y_xxxxxy_xxxy, g_0_y_xxxxxy_xxxz, g_0_y_xxxxxy_xxy, g_0_y_xxxxxy_xxyy, g_0_y_xxxxxy_xxyz, g_0_y_xxxxxy_xxz, g_0_y_xxxxxy_xxzz, g_0_y_xxxxxy_xyy, g_0_y_xxxxxy_xyyy, g_0_y_xxxxxy_xyyz, g_0_y_xxxxxy_xyz, g_0_y_xxxxxy_xyzz, g_0_y_xxxxxy_xzz, g_0_y_xxxxxy_xzzz, g_0_y_xxxxxy_yyy, g_0_y_xxxxxy_yyz, g_0_y_xxxxxy_yzz, g_0_y_xxxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxy_xxx[k] = -g_0_y_xxxxxy_xxx[k] * ab_x + g_0_y_xxxxxy_xxxx[k];

                g_0_y_xxxxxxy_xxy[k] = -g_0_y_xxxxxy_xxy[k] * ab_x + g_0_y_xxxxxy_xxxy[k];

                g_0_y_xxxxxxy_xxz[k] = -g_0_y_xxxxxy_xxz[k] * ab_x + g_0_y_xxxxxy_xxxz[k];

                g_0_y_xxxxxxy_xyy[k] = -g_0_y_xxxxxy_xyy[k] * ab_x + g_0_y_xxxxxy_xxyy[k];

                g_0_y_xxxxxxy_xyz[k] = -g_0_y_xxxxxy_xyz[k] * ab_x + g_0_y_xxxxxy_xxyz[k];

                g_0_y_xxxxxxy_xzz[k] = -g_0_y_xxxxxy_xzz[k] * ab_x + g_0_y_xxxxxy_xxzz[k];

                g_0_y_xxxxxxy_yyy[k] = -g_0_y_xxxxxy_yyy[k] * ab_x + g_0_y_xxxxxy_xyyy[k];

                g_0_y_xxxxxxy_yyz[k] = -g_0_y_xxxxxy_yyz[k] * ab_x + g_0_y_xxxxxy_xyyz[k];

                g_0_y_xxxxxxy_yzz[k] = -g_0_y_xxxxxy_yzz[k] * ab_x + g_0_y_xxxxxy_xyzz[k];

                g_0_y_xxxxxxy_zzz[k] = -g_0_y_xxxxxy_zzz[k] * ab_x + g_0_y_xxxxxy_xzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxz_xxx = cbuffer.data(kf_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxy = cbuffer.data(kf_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxz = cbuffer.data(kf_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xyy = cbuffer.data(kf_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xyz = cbuffer.data(kf_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xzz = cbuffer.data(kf_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yyy = cbuffer.data(kf_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yyz = cbuffer.data(kf_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yzz = cbuffer.data(kf_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_zzz = cbuffer.data(kf_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxz_xxx, g_0_y_xxxxxxz_xxy, g_0_y_xxxxxxz_xxz, g_0_y_xxxxxxz_xyy, g_0_y_xxxxxxz_xyz, g_0_y_xxxxxxz_xzz, g_0_y_xxxxxxz_yyy, g_0_y_xxxxxxz_yyz, g_0_y_xxxxxxz_yzz, g_0_y_xxxxxxz_zzz, g_0_y_xxxxxz_xxx, g_0_y_xxxxxz_xxxx, g_0_y_xxxxxz_xxxy, g_0_y_xxxxxz_xxxz, g_0_y_xxxxxz_xxy, g_0_y_xxxxxz_xxyy, g_0_y_xxxxxz_xxyz, g_0_y_xxxxxz_xxz, g_0_y_xxxxxz_xxzz, g_0_y_xxxxxz_xyy, g_0_y_xxxxxz_xyyy, g_0_y_xxxxxz_xyyz, g_0_y_xxxxxz_xyz, g_0_y_xxxxxz_xyzz, g_0_y_xxxxxz_xzz, g_0_y_xxxxxz_xzzz, g_0_y_xxxxxz_yyy, g_0_y_xxxxxz_yyz, g_0_y_xxxxxz_yzz, g_0_y_xxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxz_xxx[k] = -g_0_y_xxxxxz_xxx[k] * ab_x + g_0_y_xxxxxz_xxxx[k];

                g_0_y_xxxxxxz_xxy[k] = -g_0_y_xxxxxz_xxy[k] * ab_x + g_0_y_xxxxxz_xxxy[k];

                g_0_y_xxxxxxz_xxz[k] = -g_0_y_xxxxxz_xxz[k] * ab_x + g_0_y_xxxxxz_xxxz[k];

                g_0_y_xxxxxxz_xyy[k] = -g_0_y_xxxxxz_xyy[k] * ab_x + g_0_y_xxxxxz_xxyy[k];

                g_0_y_xxxxxxz_xyz[k] = -g_0_y_xxxxxz_xyz[k] * ab_x + g_0_y_xxxxxz_xxyz[k];

                g_0_y_xxxxxxz_xzz[k] = -g_0_y_xxxxxz_xzz[k] * ab_x + g_0_y_xxxxxz_xxzz[k];

                g_0_y_xxxxxxz_yyy[k] = -g_0_y_xxxxxz_yyy[k] * ab_x + g_0_y_xxxxxz_xyyy[k];

                g_0_y_xxxxxxz_yyz[k] = -g_0_y_xxxxxz_yyz[k] * ab_x + g_0_y_xxxxxz_xyyz[k];

                g_0_y_xxxxxxz_yzz[k] = -g_0_y_xxxxxz_yzz[k] * ab_x + g_0_y_xxxxxz_xyzz[k];

                g_0_y_xxxxxxz_zzz[k] = -g_0_y_xxxxxz_zzz[k] * ab_x + g_0_y_xxxxxz_xzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyy_xxx = cbuffer.data(kf_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxy = cbuffer.data(kf_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxz = cbuffer.data(kf_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xyy = cbuffer.data(kf_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xyz = cbuffer.data(kf_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xzz = cbuffer.data(kf_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yyy = cbuffer.data(kf_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yyz = cbuffer.data(kf_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yzz = cbuffer.data(kf_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_zzz = cbuffer.data(kf_geom_01_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyy_xxx, g_0_y_xxxxxyy_xxy, g_0_y_xxxxxyy_xxz, g_0_y_xxxxxyy_xyy, g_0_y_xxxxxyy_xyz, g_0_y_xxxxxyy_xzz, g_0_y_xxxxxyy_yyy, g_0_y_xxxxxyy_yyz, g_0_y_xxxxxyy_yzz, g_0_y_xxxxxyy_zzz, g_0_y_xxxxyy_xxx, g_0_y_xxxxyy_xxxx, g_0_y_xxxxyy_xxxy, g_0_y_xxxxyy_xxxz, g_0_y_xxxxyy_xxy, g_0_y_xxxxyy_xxyy, g_0_y_xxxxyy_xxyz, g_0_y_xxxxyy_xxz, g_0_y_xxxxyy_xxzz, g_0_y_xxxxyy_xyy, g_0_y_xxxxyy_xyyy, g_0_y_xxxxyy_xyyz, g_0_y_xxxxyy_xyz, g_0_y_xxxxyy_xyzz, g_0_y_xxxxyy_xzz, g_0_y_xxxxyy_xzzz, g_0_y_xxxxyy_yyy, g_0_y_xxxxyy_yyz, g_0_y_xxxxyy_yzz, g_0_y_xxxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyy_xxx[k] = -g_0_y_xxxxyy_xxx[k] * ab_x + g_0_y_xxxxyy_xxxx[k];

                g_0_y_xxxxxyy_xxy[k] = -g_0_y_xxxxyy_xxy[k] * ab_x + g_0_y_xxxxyy_xxxy[k];

                g_0_y_xxxxxyy_xxz[k] = -g_0_y_xxxxyy_xxz[k] * ab_x + g_0_y_xxxxyy_xxxz[k];

                g_0_y_xxxxxyy_xyy[k] = -g_0_y_xxxxyy_xyy[k] * ab_x + g_0_y_xxxxyy_xxyy[k];

                g_0_y_xxxxxyy_xyz[k] = -g_0_y_xxxxyy_xyz[k] * ab_x + g_0_y_xxxxyy_xxyz[k];

                g_0_y_xxxxxyy_xzz[k] = -g_0_y_xxxxyy_xzz[k] * ab_x + g_0_y_xxxxyy_xxzz[k];

                g_0_y_xxxxxyy_yyy[k] = -g_0_y_xxxxyy_yyy[k] * ab_x + g_0_y_xxxxyy_xyyy[k];

                g_0_y_xxxxxyy_yyz[k] = -g_0_y_xxxxyy_yyz[k] * ab_x + g_0_y_xxxxyy_xyyz[k];

                g_0_y_xxxxxyy_yzz[k] = -g_0_y_xxxxyy_yzz[k] * ab_x + g_0_y_xxxxyy_xyzz[k];

                g_0_y_xxxxxyy_zzz[k] = -g_0_y_xxxxyy_zzz[k] * ab_x + g_0_y_xxxxyy_xzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyz_xxx = cbuffer.data(kf_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxy = cbuffer.data(kf_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxz = cbuffer.data(kf_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xyy = cbuffer.data(kf_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xyz = cbuffer.data(kf_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xzz = cbuffer.data(kf_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yyy = cbuffer.data(kf_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yyz = cbuffer.data(kf_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yzz = cbuffer.data(kf_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_zzz = cbuffer.data(kf_geom_01_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyz_xxx, g_0_y_xxxxxyz_xxy, g_0_y_xxxxxyz_xxz, g_0_y_xxxxxyz_xyy, g_0_y_xxxxxyz_xyz, g_0_y_xxxxxyz_xzz, g_0_y_xxxxxyz_yyy, g_0_y_xxxxxyz_yyz, g_0_y_xxxxxyz_yzz, g_0_y_xxxxxyz_zzz, g_0_y_xxxxyz_xxx, g_0_y_xxxxyz_xxxx, g_0_y_xxxxyz_xxxy, g_0_y_xxxxyz_xxxz, g_0_y_xxxxyz_xxy, g_0_y_xxxxyz_xxyy, g_0_y_xxxxyz_xxyz, g_0_y_xxxxyz_xxz, g_0_y_xxxxyz_xxzz, g_0_y_xxxxyz_xyy, g_0_y_xxxxyz_xyyy, g_0_y_xxxxyz_xyyz, g_0_y_xxxxyz_xyz, g_0_y_xxxxyz_xyzz, g_0_y_xxxxyz_xzz, g_0_y_xxxxyz_xzzz, g_0_y_xxxxyz_yyy, g_0_y_xxxxyz_yyz, g_0_y_xxxxyz_yzz, g_0_y_xxxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyz_xxx[k] = -g_0_y_xxxxyz_xxx[k] * ab_x + g_0_y_xxxxyz_xxxx[k];

                g_0_y_xxxxxyz_xxy[k] = -g_0_y_xxxxyz_xxy[k] * ab_x + g_0_y_xxxxyz_xxxy[k];

                g_0_y_xxxxxyz_xxz[k] = -g_0_y_xxxxyz_xxz[k] * ab_x + g_0_y_xxxxyz_xxxz[k];

                g_0_y_xxxxxyz_xyy[k] = -g_0_y_xxxxyz_xyy[k] * ab_x + g_0_y_xxxxyz_xxyy[k];

                g_0_y_xxxxxyz_xyz[k] = -g_0_y_xxxxyz_xyz[k] * ab_x + g_0_y_xxxxyz_xxyz[k];

                g_0_y_xxxxxyz_xzz[k] = -g_0_y_xxxxyz_xzz[k] * ab_x + g_0_y_xxxxyz_xxzz[k];

                g_0_y_xxxxxyz_yyy[k] = -g_0_y_xxxxyz_yyy[k] * ab_x + g_0_y_xxxxyz_xyyy[k];

                g_0_y_xxxxxyz_yyz[k] = -g_0_y_xxxxyz_yyz[k] * ab_x + g_0_y_xxxxyz_xyyz[k];

                g_0_y_xxxxxyz_yzz[k] = -g_0_y_xxxxyz_yzz[k] * ab_x + g_0_y_xxxxyz_xyzz[k];

                g_0_y_xxxxxyz_zzz[k] = -g_0_y_xxxxyz_zzz[k] * ab_x + g_0_y_xxxxyz_xzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxzz_xxx = cbuffer.data(kf_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxy = cbuffer.data(kf_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxz = cbuffer.data(kf_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xyy = cbuffer.data(kf_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xyz = cbuffer.data(kf_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xzz = cbuffer.data(kf_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yyy = cbuffer.data(kf_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yyz = cbuffer.data(kf_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yzz = cbuffer.data(kf_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_zzz = cbuffer.data(kf_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxzz_xxx, g_0_y_xxxxxzz_xxy, g_0_y_xxxxxzz_xxz, g_0_y_xxxxxzz_xyy, g_0_y_xxxxxzz_xyz, g_0_y_xxxxxzz_xzz, g_0_y_xxxxxzz_yyy, g_0_y_xxxxxzz_yyz, g_0_y_xxxxxzz_yzz, g_0_y_xxxxxzz_zzz, g_0_y_xxxxzz_xxx, g_0_y_xxxxzz_xxxx, g_0_y_xxxxzz_xxxy, g_0_y_xxxxzz_xxxz, g_0_y_xxxxzz_xxy, g_0_y_xxxxzz_xxyy, g_0_y_xxxxzz_xxyz, g_0_y_xxxxzz_xxz, g_0_y_xxxxzz_xxzz, g_0_y_xxxxzz_xyy, g_0_y_xxxxzz_xyyy, g_0_y_xxxxzz_xyyz, g_0_y_xxxxzz_xyz, g_0_y_xxxxzz_xyzz, g_0_y_xxxxzz_xzz, g_0_y_xxxxzz_xzzz, g_0_y_xxxxzz_yyy, g_0_y_xxxxzz_yyz, g_0_y_xxxxzz_yzz, g_0_y_xxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxzz_xxx[k] = -g_0_y_xxxxzz_xxx[k] * ab_x + g_0_y_xxxxzz_xxxx[k];

                g_0_y_xxxxxzz_xxy[k] = -g_0_y_xxxxzz_xxy[k] * ab_x + g_0_y_xxxxzz_xxxy[k];

                g_0_y_xxxxxzz_xxz[k] = -g_0_y_xxxxzz_xxz[k] * ab_x + g_0_y_xxxxzz_xxxz[k];

                g_0_y_xxxxxzz_xyy[k] = -g_0_y_xxxxzz_xyy[k] * ab_x + g_0_y_xxxxzz_xxyy[k];

                g_0_y_xxxxxzz_xyz[k] = -g_0_y_xxxxzz_xyz[k] * ab_x + g_0_y_xxxxzz_xxyz[k];

                g_0_y_xxxxxzz_xzz[k] = -g_0_y_xxxxzz_xzz[k] * ab_x + g_0_y_xxxxzz_xxzz[k];

                g_0_y_xxxxxzz_yyy[k] = -g_0_y_xxxxzz_yyy[k] * ab_x + g_0_y_xxxxzz_xyyy[k];

                g_0_y_xxxxxzz_yyz[k] = -g_0_y_xxxxzz_yyz[k] * ab_x + g_0_y_xxxxzz_xyyz[k];

                g_0_y_xxxxxzz_yzz[k] = -g_0_y_xxxxzz_yzz[k] * ab_x + g_0_y_xxxxzz_xyzz[k];

                g_0_y_xxxxxzz_zzz[k] = -g_0_y_xxxxzz_zzz[k] * ab_x + g_0_y_xxxxzz_xzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyy_xxx = cbuffer.data(kf_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxy = cbuffer.data(kf_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxz = cbuffer.data(kf_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xyy = cbuffer.data(kf_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xyz = cbuffer.data(kf_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xzz = cbuffer.data(kf_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yyy = cbuffer.data(kf_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yyz = cbuffer.data(kf_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yzz = cbuffer.data(kf_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_zzz = cbuffer.data(kf_geom_01_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyy_xxx, g_0_y_xxxxyyy_xxy, g_0_y_xxxxyyy_xxz, g_0_y_xxxxyyy_xyy, g_0_y_xxxxyyy_xyz, g_0_y_xxxxyyy_xzz, g_0_y_xxxxyyy_yyy, g_0_y_xxxxyyy_yyz, g_0_y_xxxxyyy_yzz, g_0_y_xxxxyyy_zzz, g_0_y_xxxyyy_xxx, g_0_y_xxxyyy_xxxx, g_0_y_xxxyyy_xxxy, g_0_y_xxxyyy_xxxz, g_0_y_xxxyyy_xxy, g_0_y_xxxyyy_xxyy, g_0_y_xxxyyy_xxyz, g_0_y_xxxyyy_xxz, g_0_y_xxxyyy_xxzz, g_0_y_xxxyyy_xyy, g_0_y_xxxyyy_xyyy, g_0_y_xxxyyy_xyyz, g_0_y_xxxyyy_xyz, g_0_y_xxxyyy_xyzz, g_0_y_xxxyyy_xzz, g_0_y_xxxyyy_xzzz, g_0_y_xxxyyy_yyy, g_0_y_xxxyyy_yyz, g_0_y_xxxyyy_yzz, g_0_y_xxxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyy_xxx[k] = -g_0_y_xxxyyy_xxx[k] * ab_x + g_0_y_xxxyyy_xxxx[k];

                g_0_y_xxxxyyy_xxy[k] = -g_0_y_xxxyyy_xxy[k] * ab_x + g_0_y_xxxyyy_xxxy[k];

                g_0_y_xxxxyyy_xxz[k] = -g_0_y_xxxyyy_xxz[k] * ab_x + g_0_y_xxxyyy_xxxz[k];

                g_0_y_xxxxyyy_xyy[k] = -g_0_y_xxxyyy_xyy[k] * ab_x + g_0_y_xxxyyy_xxyy[k];

                g_0_y_xxxxyyy_xyz[k] = -g_0_y_xxxyyy_xyz[k] * ab_x + g_0_y_xxxyyy_xxyz[k];

                g_0_y_xxxxyyy_xzz[k] = -g_0_y_xxxyyy_xzz[k] * ab_x + g_0_y_xxxyyy_xxzz[k];

                g_0_y_xxxxyyy_yyy[k] = -g_0_y_xxxyyy_yyy[k] * ab_x + g_0_y_xxxyyy_xyyy[k];

                g_0_y_xxxxyyy_yyz[k] = -g_0_y_xxxyyy_yyz[k] * ab_x + g_0_y_xxxyyy_xyyz[k];

                g_0_y_xxxxyyy_yzz[k] = -g_0_y_xxxyyy_yzz[k] * ab_x + g_0_y_xxxyyy_xyzz[k];

                g_0_y_xxxxyyy_zzz[k] = -g_0_y_xxxyyy_zzz[k] * ab_x + g_0_y_xxxyyy_xzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyz_xxx = cbuffer.data(kf_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxy = cbuffer.data(kf_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxz = cbuffer.data(kf_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xyy = cbuffer.data(kf_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xyz = cbuffer.data(kf_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xzz = cbuffer.data(kf_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yyy = cbuffer.data(kf_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yyz = cbuffer.data(kf_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yzz = cbuffer.data(kf_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_zzz = cbuffer.data(kf_geom_01_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyz_xxx, g_0_y_xxxxyyz_xxy, g_0_y_xxxxyyz_xxz, g_0_y_xxxxyyz_xyy, g_0_y_xxxxyyz_xyz, g_0_y_xxxxyyz_xzz, g_0_y_xxxxyyz_yyy, g_0_y_xxxxyyz_yyz, g_0_y_xxxxyyz_yzz, g_0_y_xxxxyyz_zzz, g_0_y_xxxyyz_xxx, g_0_y_xxxyyz_xxxx, g_0_y_xxxyyz_xxxy, g_0_y_xxxyyz_xxxz, g_0_y_xxxyyz_xxy, g_0_y_xxxyyz_xxyy, g_0_y_xxxyyz_xxyz, g_0_y_xxxyyz_xxz, g_0_y_xxxyyz_xxzz, g_0_y_xxxyyz_xyy, g_0_y_xxxyyz_xyyy, g_0_y_xxxyyz_xyyz, g_0_y_xxxyyz_xyz, g_0_y_xxxyyz_xyzz, g_0_y_xxxyyz_xzz, g_0_y_xxxyyz_xzzz, g_0_y_xxxyyz_yyy, g_0_y_xxxyyz_yyz, g_0_y_xxxyyz_yzz, g_0_y_xxxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyz_xxx[k] = -g_0_y_xxxyyz_xxx[k] * ab_x + g_0_y_xxxyyz_xxxx[k];

                g_0_y_xxxxyyz_xxy[k] = -g_0_y_xxxyyz_xxy[k] * ab_x + g_0_y_xxxyyz_xxxy[k];

                g_0_y_xxxxyyz_xxz[k] = -g_0_y_xxxyyz_xxz[k] * ab_x + g_0_y_xxxyyz_xxxz[k];

                g_0_y_xxxxyyz_xyy[k] = -g_0_y_xxxyyz_xyy[k] * ab_x + g_0_y_xxxyyz_xxyy[k];

                g_0_y_xxxxyyz_xyz[k] = -g_0_y_xxxyyz_xyz[k] * ab_x + g_0_y_xxxyyz_xxyz[k];

                g_0_y_xxxxyyz_xzz[k] = -g_0_y_xxxyyz_xzz[k] * ab_x + g_0_y_xxxyyz_xxzz[k];

                g_0_y_xxxxyyz_yyy[k] = -g_0_y_xxxyyz_yyy[k] * ab_x + g_0_y_xxxyyz_xyyy[k];

                g_0_y_xxxxyyz_yyz[k] = -g_0_y_xxxyyz_yyz[k] * ab_x + g_0_y_xxxyyz_xyyz[k];

                g_0_y_xxxxyyz_yzz[k] = -g_0_y_xxxyyz_yzz[k] * ab_x + g_0_y_xxxyyz_xyzz[k];

                g_0_y_xxxxyyz_zzz[k] = -g_0_y_xxxyyz_zzz[k] * ab_x + g_0_y_xxxyyz_xzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyzz_xxx = cbuffer.data(kf_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxy = cbuffer.data(kf_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxz = cbuffer.data(kf_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xyy = cbuffer.data(kf_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xyz = cbuffer.data(kf_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xzz = cbuffer.data(kf_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yyy = cbuffer.data(kf_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yyz = cbuffer.data(kf_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yzz = cbuffer.data(kf_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_zzz = cbuffer.data(kf_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyzz_xxx, g_0_y_xxxxyzz_xxy, g_0_y_xxxxyzz_xxz, g_0_y_xxxxyzz_xyy, g_0_y_xxxxyzz_xyz, g_0_y_xxxxyzz_xzz, g_0_y_xxxxyzz_yyy, g_0_y_xxxxyzz_yyz, g_0_y_xxxxyzz_yzz, g_0_y_xxxxyzz_zzz, g_0_y_xxxyzz_xxx, g_0_y_xxxyzz_xxxx, g_0_y_xxxyzz_xxxy, g_0_y_xxxyzz_xxxz, g_0_y_xxxyzz_xxy, g_0_y_xxxyzz_xxyy, g_0_y_xxxyzz_xxyz, g_0_y_xxxyzz_xxz, g_0_y_xxxyzz_xxzz, g_0_y_xxxyzz_xyy, g_0_y_xxxyzz_xyyy, g_0_y_xxxyzz_xyyz, g_0_y_xxxyzz_xyz, g_0_y_xxxyzz_xyzz, g_0_y_xxxyzz_xzz, g_0_y_xxxyzz_xzzz, g_0_y_xxxyzz_yyy, g_0_y_xxxyzz_yyz, g_0_y_xxxyzz_yzz, g_0_y_xxxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyzz_xxx[k] = -g_0_y_xxxyzz_xxx[k] * ab_x + g_0_y_xxxyzz_xxxx[k];

                g_0_y_xxxxyzz_xxy[k] = -g_0_y_xxxyzz_xxy[k] * ab_x + g_0_y_xxxyzz_xxxy[k];

                g_0_y_xxxxyzz_xxz[k] = -g_0_y_xxxyzz_xxz[k] * ab_x + g_0_y_xxxyzz_xxxz[k];

                g_0_y_xxxxyzz_xyy[k] = -g_0_y_xxxyzz_xyy[k] * ab_x + g_0_y_xxxyzz_xxyy[k];

                g_0_y_xxxxyzz_xyz[k] = -g_0_y_xxxyzz_xyz[k] * ab_x + g_0_y_xxxyzz_xxyz[k];

                g_0_y_xxxxyzz_xzz[k] = -g_0_y_xxxyzz_xzz[k] * ab_x + g_0_y_xxxyzz_xxzz[k];

                g_0_y_xxxxyzz_yyy[k] = -g_0_y_xxxyzz_yyy[k] * ab_x + g_0_y_xxxyzz_xyyy[k];

                g_0_y_xxxxyzz_yyz[k] = -g_0_y_xxxyzz_yyz[k] * ab_x + g_0_y_xxxyzz_xyyz[k];

                g_0_y_xxxxyzz_yzz[k] = -g_0_y_xxxyzz_yzz[k] * ab_x + g_0_y_xxxyzz_xyzz[k];

                g_0_y_xxxxyzz_zzz[k] = -g_0_y_xxxyzz_zzz[k] * ab_x + g_0_y_xxxyzz_xzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzzz_xxx = cbuffer.data(kf_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxy = cbuffer.data(kf_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxz = cbuffer.data(kf_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xyy = cbuffer.data(kf_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xyz = cbuffer.data(kf_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xzz = cbuffer.data(kf_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yyy = cbuffer.data(kf_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yyz = cbuffer.data(kf_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yzz = cbuffer.data(kf_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_zzz = cbuffer.data(kf_geom_01_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzzz_xxx, g_0_y_xxxxzzz_xxy, g_0_y_xxxxzzz_xxz, g_0_y_xxxxzzz_xyy, g_0_y_xxxxzzz_xyz, g_0_y_xxxxzzz_xzz, g_0_y_xxxxzzz_yyy, g_0_y_xxxxzzz_yyz, g_0_y_xxxxzzz_yzz, g_0_y_xxxxzzz_zzz, g_0_y_xxxzzz_xxx, g_0_y_xxxzzz_xxxx, g_0_y_xxxzzz_xxxy, g_0_y_xxxzzz_xxxz, g_0_y_xxxzzz_xxy, g_0_y_xxxzzz_xxyy, g_0_y_xxxzzz_xxyz, g_0_y_xxxzzz_xxz, g_0_y_xxxzzz_xxzz, g_0_y_xxxzzz_xyy, g_0_y_xxxzzz_xyyy, g_0_y_xxxzzz_xyyz, g_0_y_xxxzzz_xyz, g_0_y_xxxzzz_xyzz, g_0_y_xxxzzz_xzz, g_0_y_xxxzzz_xzzz, g_0_y_xxxzzz_yyy, g_0_y_xxxzzz_yyz, g_0_y_xxxzzz_yzz, g_0_y_xxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzzz_xxx[k] = -g_0_y_xxxzzz_xxx[k] * ab_x + g_0_y_xxxzzz_xxxx[k];

                g_0_y_xxxxzzz_xxy[k] = -g_0_y_xxxzzz_xxy[k] * ab_x + g_0_y_xxxzzz_xxxy[k];

                g_0_y_xxxxzzz_xxz[k] = -g_0_y_xxxzzz_xxz[k] * ab_x + g_0_y_xxxzzz_xxxz[k];

                g_0_y_xxxxzzz_xyy[k] = -g_0_y_xxxzzz_xyy[k] * ab_x + g_0_y_xxxzzz_xxyy[k];

                g_0_y_xxxxzzz_xyz[k] = -g_0_y_xxxzzz_xyz[k] * ab_x + g_0_y_xxxzzz_xxyz[k];

                g_0_y_xxxxzzz_xzz[k] = -g_0_y_xxxzzz_xzz[k] * ab_x + g_0_y_xxxzzz_xxzz[k];

                g_0_y_xxxxzzz_yyy[k] = -g_0_y_xxxzzz_yyy[k] * ab_x + g_0_y_xxxzzz_xyyy[k];

                g_0_y_xxxxzzz_yyz[k] = -g_0_y_xxxzzz_yyz[k] * ab_x + g_0_y_xxxzzz_xyyz[k];

                g_0_y_xxxxzzz_yzz[k] = -g_0_y_xxxzzz_yzz[k] * ab_x + g_0_y_xxxzzz_xyzz[k];

                g_0_y_xxxxzzz_zzz[k] = -g_0_y_xxxzzz_zzz[k] * ab_x + g_0_y_xxxzzz_xzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyy_xxx = cbuffer.data(kf_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxy = cbuffer.data(kf_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxz = cbuffer.data(kf_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xyy = cbuffer.data(kf_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xyz = cbuffer.data(kf_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xzz = cbuffer.data(kf_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yyy = cbuffer.data(kf_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yyz = cbuffer.data(kf_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yzz = cbuffer.data(kf_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_zzz = cbuffer.data(kf_geom_01_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyy_xxx, g_0_y_xxxyyyy_xxy, g_0_y_xxxyyyy_xxz, g_0_y_xxxyyyy_xyy, g_0_y_xxxyyyy_xyz, g_0_y_xxxyyyy_xzz, g_0_y_xxxyyyy_yyy, g_0_y_xxxyyyy_yyz, g_0_y_xxxyyyy_yzz, g_0_y_xxxyyyy_zzz, g_0_y_xxyyyy_xxx, g_0_y_xxyyyy_xxxx, g_0_y_xxyyyy_xxxy, g_0_y_xxyyyy_xxxz, g_0_y_xxyyyy_xxy, g_0_y_xxyyyy_xxyy, g_0_y_xxyyyy_xxyz, g_0_y_xxyyyy_xxz, g_0_y_xxyyyy_xxzz, g_0_y_xxyyyy_xyy, g_0_y_xxyyyy_xyyy, g_0_y_xxyyyy_xyyz, g_0_y_xxyyyy_xyz, g_0_y_xxyyyy_xyzz, g_0_y_xxyyyy_xzz, g_0_y_xxyyyy_xzzz, g_0_y_xxyyyy_yyy, g_0_y_xxyyyy_yyz, g_0_y_xxyyyy_yzz, g_0_y_xxyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyy_xxx[k] = -g_0_y_xxyyyy_xxx[k] * ab_x + g_0_y_xxyyyy_xxxx[k];

                g_0_y_xxxyyyy_xxy[k] = -g_0_y_xxyyyy_xxy[k] * ab_x + g_0_y_xxyyyy_xxxy[k];

                g_0_y_xxxyyyy_xxz[k] = -g_0_y_xxyyyy_xxz[k] * ab_x + g_0_y_xxyyyy_xxxz[k];

                g_0_y_xxxyyyy_xyy[k] = -g_0_y_xxyyyy_xyy[k] * ab_x + g_0_y_xxyyyy_xxyy[k];

                g_0_y_xxxyyyy_xyz[k] = -g_0_y_xxyyyy_xyz[k] * ab_x + g_0_y_xxyyyy_xxyz[k];

                g_0_y_xxxyyyy_xzz[k] = -g_0_y_xxyyyy_xzz[k] * ab_x + g_0_y_xxyyyy_xxzz[k];

                g_0_y_xxxyyyy_yyy[k] = -g_0_y_xxyyyy_yyy[k] * ab_x + g_0_y_xxyyyy_xyyy[k];

                g_0_y_xxxyyyy_yyz[k] = -g_0_y_xxyyyy_yyz[k] * ab_x + g_0_y_xxyyyy_xyyz[k];

                g_0_y_xxxyyyy_yzz[k] = -g_0_y_xxyyyy_yzz[k] * ab_x + g_0_y_xxyyyy_xyzz[k];

                g_0_y_xxxyyyy_zzz[k] = -g_0_y_xxyyyy_zzz[k] * ab_x + g_0_y_xxyyyy_xzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyz_xxx = cbuffer.data(kf_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxy = cbuffer.data(kf_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxz = cbuffer.data(kf_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xyy = cbuffer.data(kf_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xyz = cbuffer.data(kf_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xzz = cbuffer.data(kf_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yyy = cbuffer.data(kf_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yyz = cbuffer.data(kf_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yzz = cbuffer.data(kf_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_zzz = cbuffer.data(kf_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyz_xxx, g_0_y_xxxyyyz_xxy, g_0_y_xxxyyyz_xxz, g_0_y_xxxyyyz_xyy, g_0_y_xxxyyyz_xyz, g_0_y_xxxyyyz_xzz, g_0_y_xxxyyyz_yyy, g_0_y_xxxyyyz_yyz, g_0_y_xxxyyyz_yzz, g_0_y_xxxyyyz_zzz, g_0_y_xxyyyz_xxx, g_0_y_xxyyyz_xxxx, g_0_y_xxyyyz_xxxy, g_0_y_xxyyyz_xxxz, g_0_y_xxyyyz_xxy, g_0_y_xxyyyz_xxyy, g_0_y_xxyyyz_xxyz, g_0_y_xxyyyz_xxz, g_0_y_xxyyyz_xxzz, g_0_y_xxyyyz_xyy, g_0_y_xxyyyz_xyyy, g_0_y_xxyyyz_xyyz, g_0_y_xxyyyz_xyz, g_0_y_xxyyyz_xyzz, g_0_y_xxyyyz_xzz, g_0_y_xxyyyz_xzzz, g_0_y_xxyyyz_yyy, g_0_y_xxyyyz_yyz, g_0_y_xxyyyz_yzz, g_0_y_xxyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyz_xxx[k] = -g_0_y_xxyyyz_xxx[k] * ab_x + g_0_y_xxyyyz_xxxx[k];

                g_0_y_xxxyyyz_xxy[k] = -g_0_y_xxyyyz_xxy[k] * ab_x + g_0_y_xxyyyz_xxxy[k];

                g_0_y_xxxyyyz_xxz[k] = -g_0_y_xxyyyz_xxz[k] * ab_x + g_0_y_xxyyyz_xxxz[k];

                g_0_y_xxxyyyz_xyy[k] = -g_0_y_xxyyyz_xyy[k] * ab_x + g_0_y_xxyyyz_xxyy[k];

                g_0_y_xxxyyyz_xyz[k] = -g_0_y_xxyyyz_xyz[k] * ab_x + g_0_y_xxyyyz_xxyz[k];

                g_0_y_xxxyyyz_xzz[k] = -g_0_y_xxyyyz_xzz[k] * ab_x + g_0_y_xxyyyz_xxzz[k];

                g_0_y_xxxyyyz_yyy[k] = -g_0_y_xxyyyz_yyy[k] * ab_x + g_0_y_xxyyyz_xyyy[k];

                g_0_y_xxxyyyz_yyz[k] = -g_0_y_xxyyyz_yyz[k] * ab_x + g_0_y_xxyyyz_xyyz[k];

                g_0_y_xxxyyyz_yzz[k] = -g_0_y_xxyyyz_yzz[k] * ab_x + g_0_y_xxyyyz_xyzz[k];

                g_0_y_xxxyyyz_zzz[k] = -g_0_y_xxyyyz_zzz[k] * ab_x + g_0_y_xxyyyz_xzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyzz_xxx = cbuffer.data(kf_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxy = cbuffer.data(kf_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxz = cbuffer.data(kf_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xyy = cbuffer.data(kf_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xyz = cbuffer.data(kf_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xzz = cbuffer.data(kf_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yyy = cbuffer.data(kf_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yyz = cbuffer.data(kf_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yzz = cbuffer.data(kf_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_zzz = cbuffer.data(kf_geom_01_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyzz_xxx, g_0_y_xxxyyzz_xxy, g_0_y_xxxyyzz_xxz, g_0_y_xxxyyzz_xyy, g_0_y_xxxyyzz_xyz, g_0_y_xxxyyzz_xzz, g_0_y_xxxyyzz_yyy, g_0_y_xxxyyzz_yyz, g_0_y_xxxyyzz_yzz, g_0_y_xxxyyzz_zzz, g_0_y_xxyyzz_xxx, g_0_y_xxyyzz_xxxx, g_0_y_xxyyzz_xxxy, g_0_y_xxyyzz_xxxz, g_0_y_xxyyzz_xxy, g_0_y_xxyyzz_xxyy, g_0_y_xxyyzz_xxyz, g_0_y_xxyyzz_xxz, g_0_y_xxyyzz_xxzz, g_0_y_xxyyzz_xyy, g_0_y_xxyyzz_xyyy, g_0_y_xxyyzz_xyyz, g_0_y_xxyyzz_xyz, g_0_y_xxyyzz_xyzz, g_0_y_xxyyzz_xzz, g_0_y_xxyyzz_xzzz, g_0_y_xxyyzz_yyy, g_0_y_xxyyzz_yyz, g_0_y_xxyyzz_yzz, g_0_y_xxyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyzz_xxx[k] = -g_0_y_xxyyzz_xxx[k] * ab_x + g_0_y_xxyyzz_xxxx[k];

                g_0_y_xxxyyzz_xxy[k] = -g_0_y_xxyyzz_xxy[k] * ab_x + g_0_y_xxyyzz_xxxy[k];

                g_0_y_xxxyyzz_xxz[k] = -g_0_y_xxyyzz_xxz[k] * ab_x + g_0_y_xxyyzz_xxxz[k];

                g_0_y_xxxyyzz_xyy[k] = -g_0_y_xxyyzz_xyy[k] * ab_x + g_0_y_xxyyzz_xxyy[k];

                g_0_y_xxxyyzz_xyz[k] = -g_0_y_xxyyzz_xyz[k] * ab_x + g_0_y_xxyyzz_xxyz[k];

                g_0_y_xxxyyzz_xzz[k] = -g_0_y_xxyyzz_xzz[k] * ab_x + g_0_y_xxyyzz_xxzz[k];

                g_0_y_xxxyyzz_yyy[k] = -g_0_y_xxyyzz_yyy[k] * ab_x + g_0_y_xxyyzz_xyyy[k];

                g_0_y_xxxyyzz_yyz[k] = -g_0_y_xxyyzz_yyz[k] * ab_x + g_0_y_xxyyzz_xyyz[k];

                g_0_y_xxxyyzz_yzz[k] = -g_0_y_xxyyzz_yzz[k] * ab_x + g_0_y_xxyyzz_xyzz[k];

                g_0_y_xxxyyzz_zzz[k] = -g_0_y_xxyyzz_zzz[k] * ab_x + g_0_y_xxyyzz_xzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzzz_xxx = cbuffer.data(kf_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxy = cbuffer.data(kf_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxz = cbuffer.data(kf_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xyy = cbuffer.data(kf_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xyz = cbuffer.data(kf_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xzz = cbuffer.data(kf_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yyy = cbuffer.data(kf_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yyz = cbuffer.data(kf_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yzz = cbuffer.data(kf_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_zzz = cbuffer.data(kf_geom_01_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzzz_xxx, g_0_y_xxxyzzz_xxy, g_0_y_xxxyzzz_xxz, g_0_y_xxxyzzz_xyy, g_0_y_xxxyzzz_xyz, g_0_y_xxxyzzz_xzz, g_0_y_xxxyzzz_yyy, g_0_y_xxxyzzz_yyz, g_0_y_xxxyzzz_yzz, g_0_y_xxxyzzz_zzz, g_0_y_xxyzzz_xxx, g_0_y_xxyzzz_xxxx, g_0_y_xxyzzz_xxxy, g_0_y_xxyzzz_xxxz, g_0_y_xxyzzz_xxy, g_0_y_xxyzzz_xxyy, g_0_y_xxyzzz_xxyz, g_0_y_xxyzzz_xxz, g_0_y_xxyzzz_xxzz, g_0_y_xxyzzz_xyy, g_0_y_xxyzzz_xyyy, g_0_y_xxyzzz_xyyz, g_0_y_xxyzzz_xyz, g_0_y_xxyzzz_xyzz, g_0_y_xxyzzz_xzz, g_0_y_xxyzzz_xzzz, g_0_y_xxyzzz_yyy, g_0_y_xxyzzz_yyz, g_0_y_xxyzzz_yzz, g_0_y_xxyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzzz_xxx[k] = -g_0_y_xxyzzz_xxx[k] * ab_x + g_0_y_xxyzzz_xxxx[k];

                g_0_y_xxxyzzz_xxy[k] = -g_0_y_xxyzzz_xxy[k] * ab_x + g_0_y_xxyzzz_xxxy[k];

                g_0_y_xxxyzzz_xxz[k] = -g_0_y_xxyzzz_xxz[k] * ab_x + g_0_y_xxyzzz_xxxz[k];

                g_0_y_xxxyzzz_xyy[k] = -g_0_y_xxyzzz_xyy[k] * ab_x + g_0_y_xxyzzz_xxyy[k];

                g_0_y_xxxyzzz_xyz[k] = -g_0_y_xxyzzz_xyz[k] * ab_x + g_0_y_xxyzzz_xxyz[k];

                g_0_y_xxxyzzz_xzz[k] = -g_0_y_xxyzzz_xzz[k] * ab_x + g_0_y_xxyzzz_xxzz[k];

                g_0_y_xxxyzzz_yyy[k] = -g_0_y_xxyzzz_yyy[k] * ab_x + g_0_y_xxyzzz_xyyy[k];

                g_0_y_xxxyzzz_yyz[k] = -g_0_y_xxyzzz_yyz[k] * ab_x + g_0_y_xxyzzz_xyyz[k];

                g_0_y_xxxyzzz_yzz[k] = -g_0_y_xxyzzz_yzz[k] * ab_x + g_0_y_xxyzzz_xyzz[k];

                g_0_y_xxxyzzz_zzz[k] = -g_0_y_xxyzzz_zzz[k] * ab_x + g_0_y_xxyzzz_xzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzzz_xxx = cbuffer.data(kf_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxy = cbuffer.data(kf_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxz = cbuffer.data(kf_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xyy = cbuffer.data(kf_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xyz = cbuffer.data(kf_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xzz = cbuffer.data(kf_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yyy = cbuffer.data(kf_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yyz = cbuffer.data(kf_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yzz = cbuffer.data(kf_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_zzz = cbuffer.data(kf_geom_01_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzzz_xxx, g_0_y_xxxzzzz_xxy, g_0_y_xxxzzzz_xxz, g_0_y_xxxzzzz_xyy, g_0_y_xxxzzzz_xyz, g_0_y_xxxzzzz_xzz, g_0_y_xxxzzzz_yyy, g_0_y_xxxzzzz_yyz, g_0_y_xxxzzzz_yzz, g_0_y_xxxzzzz_zzz, g_0_y_xxzzzz_xxx, g_0_y_xxzzzz_xxxx, g_0_y_xxzzzz_xxxy, g_0_y_xxzzzz_xxxz, g_0_y_xxzzzz_xxy, g_0_y_xxzzzz_xxyy, g_0_y_xxzzzz_xxyz, g_0_y_xxzzzz_xxz, g_0_y_xxzzzz_xxzz, g_0_y_xxzzzz_xyy, g_0_y_xxzzzz_xyyy, g_0_y_xxzzzz_xyyz, g_0_y_xxzzzz_xyz, g_0_y_xxzzzz_xyzz, g_0_y_xxzzzz_xzz, g_0_y_xxzzzz_xzzz, g_0_y_xxzzzz_yyy, g_0_y_xxzzzz_yyz, g_0_y_xxzzzz_yzz, g_0_y_xxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzzz_xxx[k] = -g_0_y_xxzzzz_xxx[k] * ab_x + g_0_y_xxzzzz_xxxx[k];

                g_0_y_xxxzzzz_xxy[k] = -g_0_y_xxzzzz_xxy[k] * ab_x + g_0_y_xxzzzz_xxxy[k];

                g_0_y_xxxzzzz_xxz[k] = -g_0_y_xxzzzz_xxz[k] * ab_x + g_0_y_xxzzzz_xxxz[k];

                g_0_y_xxxzzzz_xyy[k] = -g_0_y_xxzzzz_xyy[k] * ab_x + g_0_y_xxzzzz_xxyy[k];

                g_0_y_xxxzzzz_xyz[k] = -g_0_y_xxzzzz_xyz[k] * ab_x + g_0_y_xxzzzz_xxyz[k];

                g_0_y_xxxzzzz_xzz[k] = -g_0_y_xxzzzz_xzz[k] * ab_x + g_0_y_xxzzzz_xxzz[k];

                g_0_y_xxxzzzz_yyy[k] = -g_0_y_xxzzzz_yyy[k] * ab_x + g_0_y_xxzzzz_xyyy[k];

                g_0_y_xxxzzzz_yyz[k] = -g_0_y_xxzzzz_yyz[k] * ab_x + g_0_y_xxzzzz_xyyz[k];

                g_0_y_xxxzzzz_yzz[k] = -g_0_y_xxzzzz_yzz[k] * ab_x + g_0_y_xxzzzz_xyzz[k];

                g_0_y_xxxzzzz_zzz[k] = -g_0_y_xxzzzz_zzz[k] * ab_x + g_0_y_xxzzzz_xzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyy_xxx = cbuffer.data(kf_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxy = cbuffer.data(kf_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxz = cbuffer.data(kf_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xyy = cbuffer.data(kf_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xyz = cbuffer.data(kf_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xzz = cbuffer.data(kf_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yyy = cbuffer.data(kf_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yyz = cbuffer.data(kf_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yzz = cbuffer.data(kf_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_zzz = cbuffer.data(kf_geom_01_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyy_xxx, g_0_y_xxyyyyy_xxy, g_0_y_xxyyyyy_xxz, g_0_y_xxyyyyy_xyy, g_0_y_xxyyyyy_xyz, g_0_y_xxyyyyy_xzz, g_0_y_xxyyyyy_yyy, g_0_y_xxyyyyy_yyz, g_0_y_xxyyyyy_yzz, g_0_y_xxyyyyy_zzz, g_0_y_xyyyyy_xxx, g_0_y_xyyyyy_xxxx, g_0_y_xyyyyy_xxxy, g_0_y_xyyyyy_xxxz, g_0_y_xyyyyy_xxy, g_0_y_xyyyyy_xxyy, g_0_y_xyyyyy_xxyz, g_0_y_xyyyyy_xxz, g_0_y_xyyyyy_xxzz, g_0_y_xyyyyy_xyy, g_0_y_xyyyyy_xyyy, g_0_y_xyyyyy_xyyz, g_0_y_xyyyyy_xyz, g_0_y_xyyyyy_xyzz, g_0_y_xyyyyy_xzz, g_0_y_xyyyyy_xzzz, g_0_y_xyyyyy_yyy, g_0_y_xyyyyy_yyz, g_0_y_xyyyyy_yzz, g_0_y_xyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyy_xxx[k] = -g_0_y_xyyyyy_xxx[k] * ab_x + g_0_y_xyyyyy_xxxx[k];

                g_0_y_xxyyyyy_xxy[k] = -g_0_y_xyyyyy_xxy[k] * ab_x + g_0_y_xyyyyy_xxxy[k];

                g_0_y_xxyyyyy_xxz[k] = -g_0_y_xyyyyy_xxz[k] * ab_x + g_0_y_xyyyyy_xxxz[k];

                g_0_y_xxyyyyy_xyy[k] = -g_0_y_xyyyyy_xyy[k] * ab_x + g_0_y_xyyyyy_xxyy[k];

                g_0_y_xxyyyyy_xyz[k] = -g_0_y_xyyyyy_xyz[k] * ab_x + g_0_y_xyyyyy_xxyz[k];

                g_0_y_xxyyyyy_xzz[k] = -g_0_y_xyyyyy_xzz[k] * ab_x + g_0_y_xyyyyy_xxzz[k];

                g_0_y_xxyyyyy_yyy[k] = -g_0_y_xyyyyy_yyy[k] * ab_x + g_0_y_xyyyyy_xyyy[k];

                g_0_y_xxyyyyy_yyz[k] = -g_0_y_xyyyyy_yyz[k] * ab_x + g_0_y_xyyyyy_xyyz[k];

                g_0_y_xxyyyyy_yzz[k] = -g_0_y_xyyyyy_yzz[k] * ab_x + g_0_y_xyyyyy_xyzz[k];

                g_0_y_xxyyyyy_zzz[k] = -g_0_y_xyyyyy_zzz[k] * ab_x + g_0_y_xyyyyy_xzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyz_xxx = cbuffer.data(kf_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxy = cbuffer.data(kf_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxz = cbuffer.data(kf_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xyy = cbuffer.data(kf_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xyz = cbuffer.data(kf_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xzz = cbuffer.data(kf_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yyy = cbuffer.data(kf_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yyz = cbuffer.data(kf_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yzz = cbuffer.data(kf_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_zzz = cbuffer.data(kf_geom_01_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyz_xxx, g_0_y_xxyyyyz_xxy, g_0_y_xxyyyyz_xxz, g_0_y_xxyyyyz_xyy, g_0_y_xxyyyyz_xyz, g_0_y_xxyyyyz_xzz, g_0_y_xxyyyyz_yyy, g_0_y_xxyyyyz_yyz, g_0_y_xxyyyyz_yzz, g_0_y_xxyyyyz_zzz, g_0_y_xyyyyz_xxx, g_0_y_xyyyyz_xxxx, g_0_y_xyyyyz_xxxy, g_0_y_xyyyyz_xxxz, g_0_y_xyyyyz_xxy, g_0_y_xyyyyz_xxyy, g_0_y_xyyyyz_xxyz, g_0_y_xyyyyz_xxz, g_0_y_xyyyyz_xxzz, g_0_y_xyyyyz_xyy, g_0_y_xyyyyz_xyyy, g_0_y_xyyyyz_xyyz, g_0_y_xyyyyz_xyz, g_0_y_xyyyyz_xyzz, g_0_y_xyyyyz_xzz, g_0_y_xyyyyz_xzzz, g_0_y_xyyyyz_yyy, g_0_y_xyyyyz_yyz, g_0_y_xyyyyz_yzz, g_0_y_xyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyz_xxx[k] = -g_0_y_xyyyyz_xxx[k] * ab_x + g_0_y_xyyyyz_xxxx[k];

                g_0_y_xxyyyyz_xxy[k] = -g_0_y_xyyyyz_xxy[k] * ab_x + g_0_y_xyyyyz_xxxy[k];

                g_0_y_xxyyyyz_xxz[k] = -g_0_y_xyyyyz_xxz[k] * ab_x + g_0_y_xyyyyz_xxxz[k];

                g_0_y_xxyyyyz_xyy[k] = -g_0_y_xyyyyz_xyy[k] * ab_x + g_0_y_xyyyyz_xxyy[k];

                g_0_y_xxyyyyz_xyz[k] = -g_0_y_xyyyyz_xyz[k] * ab_x + g_0_y_xyyyyz_xxyz[k];

                g_0_y_xxyyyyz_xzz[k] = -g_0_y_xyyyyz_xzz[k] * ab_x + g_0_y_xyyyyz_xxzz[k];

                g_0_y_xxyyyyz_yyy[k] = -g_0_y_xyyyyz_yyy[k] * ab_x + g_0_y_xyyyyz_xyyy[k];

                g_0_y_xxyyyyz_yyz[k] = -g_0_y_xyyyyz_yyz[k] * ab_x + g_0_y_xyyyyz_xyyz[k];

                g_0_y_xxyyyyz_yzz[k] = -g_0_y_xyyyyz_yzz[k] * ab_x + g_0_y_xyyyyz_xyzz[k];

                g_0_y_xxyyyyz_zzz[k] = -g_0_y_xyyyyz_zzz[k] * ab_x + g_0_y_xyyyyz_xzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyzz_xxx = cbuffer.data(kf_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxy = cbuffer.data(kf_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxz = cbuffer.data(kf_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xyy = cbuffer.data(kf_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xyz = cbuffer.data(kf_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xzz = cbuffer.data(kf_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yyy = cbuffer.data(kf_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yyz = cbuffer.data(kf_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yzz = cbuffer.data(kf_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_zzz = cbuffer.data(kf_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyzz_xxx, g_0_y_xxyyyzz_xxy, g_0_y_xxyyyzz_xxz, g_0_y_xxyyyzz_xyy, g_0_y_xxyyyzz_xyz, g_0_y_xxyyyzz_xzz, g_0_y_xxyyyzz_yyy, g_0_y_xxyyyzz_yyz, g_0_y_xxyyyzz_yzz, g_0_y_xxyyyzz_zzz, g_0_y_xyyyzz_xxx, g_0_y_xyyyzz_xxxx, g_0_y_xyyyzz_xxxy, g_0_y_xyyyzz_xxxz, g_0_y_xyyyzz_xxy, g_0_y_xyyyzz_xxyy, g_0_y_xyyyzz_xxyz, g_0_y_xyyyzz_xxz, g_0_y_xyyyzz_xxzz, g_0_y_xyyyzz_xyy, g_0_y_xyyyzz_xyyy, g_0_y_xyyyzz_xyyz, g_0_y_xyyyzz_xyz, g_0_y_xyyyzz_xyzz, g_0_y_xyyyzz_xzz, g_0_y_xyyyzz_xzzz, g_0_y_xyyyzz_yyy, g_0_y_xyyyzz_yyz, g_0_y_xyyyzz_yzz, g_0_y_xyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyzz_xxx[k] = -g_0_y_xyyyzz_xxx[k] * ab_x + g_0_y_xyyyzz_xxxx[k];

                g_0_y_xxyyyzz_xxy[k] = -g_0_y_xyyyzz_xxy[k] * ab_x + g_0_y_xyyyzz_xxxy[k];

                g_0_y_xxyyyzz_xxz[k] = -g_0_y_xyyyzz_xxz[k] * ab_x + g_0_y_xyyyzz_xxxz[k];

                g_0_y_xxyyyzz_xyy[k] = -g_0_y_xyyyzz_xyy[k] * ab_x + g_0_y_xyyyzz_xxyy[k];

                g_0_y_xxyyyzz_xyz[k] = -g_0_y_xyyyzz_xyz[k] * ab_x + g_0_y_xyyyzz_xxyz[k];

                g_0_y_xxyyyzz_xzz[k] = -g_0_y_xyyyzz_xzz[k] * ab_x + g_0_y_xyyyzz_xxzz[k];

                g_0_y_xxyyyzz_yyy[k] = -g_0_y_xyyyzz_yyy[k] * ab_x + g_0_y_xyyyzz_xyyy[k];

                g_0_y_xxyyyzz_yyz[k] = -g_0_y_xyyyzz_yyz[k] * ab_x + g_0_y_xyyyzz_xyyz[k];

                g_0_y_xxyyyzz_yzz[k] = -g_0_y_xyyyzz_yzz[k] * ab_x + g_0_y_xyyyzz_xyzz[k];

                g_0_y_xxyyyzz_zzz[k] = -g_0_y_xyyyzz_zzz[k] * ab_x + g_0_y_xyyyzz_xzzz[k];
            }

            /// Set up 540-550 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzzz_xxx = cbuffer.data(kf_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxy = cbuffer.data(kf_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxz = cbuffer.data(kf_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xyy = cbuffer.data(kf_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xyz = cbuffer.data(kf_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xzz = cbuffer.data(kf_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yyy = cbuffer.data(kf_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yyz = cbuffer.data(kf_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yzz = cbuffer.data(kf_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_zzz = cbuffer.data(kf_geom_01_off + 549 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzzz_xxx, g_0_y_xxyyzzz_xxy, g_0_y_xxyyzzz_xxz, g_0_y_xxyyzzz_xyy, g_0_y_xxyyzzz_xyz, g_0_y_xxyyzzz_xzz, g_0_y_xxyyzzz_yyy, g_0_y_xxyyzzz_yyz, g_0_y_xxyyzzz_yzz, g_0_y_xxyyzzz_zzz, g_0_y_xyyzzz_xxx, g_0_y_xyyzzz_xxxx, g_0_y_xyyzzz_xxxy, g_0_y_xyyzzz_xxxz, g_0_y_xyyzzz_xxy, g_0_y_xyyzzz_xxyy, g_0_y_xyyzzz_xxyz, g_0_y_xyyzzz_xxz, g_0_y_xyyzzz_xxzz, g_0_y_xyyzzz_xyy, g_0_y_xyyzzz_xyyy, g_0_y_xyyzzz_xyyz, g_0_y_xyyzzz_xyz, g_0_y_xyyzzz_xyzz, g_0_y_xyyzzz_xzz, g_0_y_xyyzzz_xzzz, g_0_y_xyyzzz_yyy, g_0_y_xyyzzz_yyz, g_0_y_xyyzzz_yzz, g_0_y_xyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzzz_xxx[k] = -g_0_y_xyyzzz_xxx[k] * ab_x + g_0_y_xyyzzz_xxxx[k];

                g_0_y_xxyyzzz_xxy[k] = -g_0_y_xyyzzz_xxy[k] * ab_x + g_0_y_xyyzzz_xxxy[k];

                g_0_y_xxyyzzz_xxz[k] = -g_0_y_xyyzzz_xxz[k] * ab_x + g_0_y_xyyzzz_xxxz[k];

                g_0_y_xxyyzzz_xyy[k] = -g_0_y_xyyzzz_xyy[k] * ab_x + g_0_y_xyyzzz_xxyy[k];

                g_0_y_xxyyzzz_xyz[k] = -g_0_y_xyyzzz_xyz[k] * ab_x + g_0_y_xyyzzz_xxyz[k];

                g_0_y_xxyyzzz_xzz[k] = -g_0_y_xyyzzz_xzz[k] * ab_x + g_0_y_xyyzzz_xxzz[k];

                g_0_y_xxyyzzz_yyy[k] = -g_0_y_xyyzzz_yyy[k] * ab_x + g_0_y_xyyzzz_xyyy[k];

                g_0_y_xxyyzzz_yyz[k] = -g_0_y_xyyzzz_yyz[k] * ab_x + g_0_y_xyyzzz_xyyz[k];

                g_0_y_xxyyzzz_yzz[k] = -g_0_y_xyyzzz_yzz[k] * ab_x + g_0_y_xyyzzz_xyzz[k];

                g_0_y_xxyyzzz_zzz[k] = -g_0_y_xyyzzz_zzz[k] * ab_x + g_0_y_xyyzzz_xzzz[k];
            }

            /// Set up 550-560 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzzz_xxx = cbuffer.data(kf_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxy = cbuffer.data(kf_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxz = cbuffer.data(kf_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xyy = cbuffer.data(kf_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xyz = cbuffer.data(kf_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xzz = cbuffer.data(kf_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yyy = cbuffer.data(kf_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yyz = cbuffer.data(kf_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yzz = cbuffer.data(kf_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_zzz = cbuffer.data(kf_geom_01_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzzz_xxx, g_0_y_xxyzzzz_xxy, g_0_y_xxyzzzz_xxz, g_0_y_xxyzzzz_xyy, g_0_y_xxyzzzz_xyz, g_0_y_xxyzzzz_xzz, g_0_y_xxyzzzz_yyy, g_0_y_xxyzzzz_yyz, g_0_y_xxyzzzz_yzz, g_0_y_xxyzzzz_zzz, g_0_y_xyzzzz_xxx, g_0_y_xyzzzz_xxxx, g_0_y_xyzzzz_xxxy, g_0_y_xyzzzz_xxxz, g_0_y_xyzzzz_xxy, g_0_y_xyzzzz_xxyy, g_0_y_xyzzzz_xxyz, g_0_y_xyzzzz_xxz, g_0_y_xyzzzz_xxzz, g_0_y_xyzzzz_xyy, g_0_y_xyzzzz_xyyy, g_0_y_xyzzzz_xyyz, g_0_y_xyzzzz_xyz, g_0_y_xyzzzz_xyzz, g_0_y_xyzzzz_xzz, g_0_y_xyzzzz_xzzz, g_0_y_xyzzzz_yyy, g_0_y_xyzzzz_yyz, g_0_y_xyzzzz_yzz, g_0_y_xyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzzz_xxx[k] = -g_0_y_xyzzzz_xxx[k] * ab_x + g_0_y_xyzzzz_xxxx[k];

                g_0_y_xxyzzzz_xxy[k] = -g_0_y_xyzzzz_xxy[k] * ab_x + g_0_y_xyzzzz_xxxy[k];

                g_0_y_xxyzzzz_xxz[k] = -g_0_y_xyzzzz_xxz[k] * ab_x + g_0_y_xyzzzz_xxxz[k];

                g_0_y_xxyzzzz_xyy[k] = -g_0_y_xyzzzz_xyy[k] * ab_x + g_0_y_xyzzzz_xxyy[k];

                g_0_y_xxyzzzz_xyz[k] = -g_0_y_xyzzzz_xyz[k] * ab_x + g_0_y_xyzzzz_xxyz[k];

                g_0_y_xxyzzzz_xzz[k] = -g_0_y_xyzzzz_xzz[k] * ab_x + g_0_y_xyzzzz_xxzz[k];

                g_0_y_xxyzzzz_yyy[k] = -g_0_y_xyzzzz_yyy[k] * ab_x + g_0_y_xyzzzz_xyyy[k];

                g_0_y_xxyzzzz_yyz[k] = -g_0_y_xyzzzz_yyz[k] * ab_x + g_0_y_xyzzzz_xyyz[k];

                g_0_y_xxyzzzz_yzz[k] = -g_0_y_xyzzzz_yzz[k] * ab_x + g_0_y_xyzzzz_xyzz[k];

                g_0_y_xxyzzzz_zzz[k] = -g_0_y_xyzzzz_zzz[k] * ab_x + g_0_y_xyzzzz_xzzz[k];
            }

            /// Set up 560-570 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzzz_xxx = cbuffer.data(kf_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxy = cbuffer.data(kf_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxz = cbuffer.data(kf_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xyy = cbuffer.data(kf_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xyz = cbuffer.data(kf_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xzz = cbuffer.data(kf_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yyy = cbuffer.data(kf_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yyz = cbuffer.data(kf_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yzz = cbuffer.data(kf_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_zzz = cbuffer.data(kf_geom_01_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzzz_xxx, g_0_y_xxzzzzz_xxy, g_0_y_xxzzzzz_xxz, g_0_y_xxzzzzz_xyy, g_0_y_xxzzzzz_xyz, g_0_y_xxzzzzz_xzz, g_0_y_xxzzzzz_yyy, g_0_y_xxzzzzz_yyz, g_0_y_xxzzzzz_yzz, g_0_y_xxzzzzz_zzz, g_0_y_xzzzzz_xxx, g_0_y_xzzzzz_xxxx, g_0_y_xzzzzz_xxxy, g_0_y_xzzzzz_xxxz, g_0_y_xzzzzz_xxy, g_0_y_xzzzzz_xxyy, g_0_y_xzzzzz_xxyz, g_0_y_xzzzzz_xxz, g_0_y_xzzzzz_xxzz, g_0_y_xzzzzz_xyy, g_0_y_xzzzzz_xyyy, g_0_y_xzzzzz_xyyz, g_0_y_xzzzzz_xyz, g_0_y_xzzzzz_xyzz, g_0_y_xzzzzz_xzz, g_0_y_xzzzzz_xzzz, g_0_y_xzzzzz_yyy, g_0_y_xzzzzz_yyz, g_0_y_xzzzzz_yzz, g_0_y_xzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzzz_xxx[k] = -g_0_y_xzzzzz_xxx[k] * ab_x + g_0_y_xzzzzz_xxxx[k];

                g_0_y_xxzzzzz_xxy[k] = -g_0_y_xzzzzz_xxy[k] * ab_x + g_0_y_xzzzzz_xxxy[k];

                g_0_y_xxzzzzz_xxz[k] = -g_0_y_xzzzzz_xxz[k] * ab_x + g_0_y_xzzzzz_xxxz[k];

                g_0_y_xxzzzzz_xyy[k] = -g_0_y_xzzzzz_xyy[k] * ab_x + g_0_y_xzzzzz_xxyy[k];

                g_0_y_xxzzzzz_xyz[k] = -g_0_y_xzzzzz_xyz[k] * ab_x + g_0_y_xzzzzz_xxyz[k];

                g_0_y_xxzzzzz_xzz[k] = -g_0_y_xzzzzz_xzz[k] * ab_x + g_0_y_xzzzzz_xxzz[k];

                g_0_y_xxzzzzz_yyy[k] = -g_0_y_xzzzzz_yyy[k] * ab_x + g_0_y_xzzzzz_xyyy[k];

                g_0_y_xxzzzzz_yyz[k] = -g_0_y_xzzzzz_yyz[k] * ab_x + g_0_y_xzzzzz_xyyz[k];

                g_0_y_xxzzzzz_yzz[k] = -g_0_y_xzzzzz_yzz[k] * ab_x + g_0_y_xzzzzz_xyzz[k];

                g_0_y_xxzzzzz_zzz[k] = -g_0_y_xzzzzz_zzz[k] * ab_x + g_0_y_xzzzzz_xzzz[k];
            }

            /// Set up 570-580 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 579 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyy_xxx, g_0_y_xyyyyyy_xxy, g_0_y_xyyyyyy_xxz, g_0_y_xyyyyyy_xyy, g_0_y_xyyyyyy_xyz, g_0_y_xyyyyyy_xzz, g_0_y_xyyyyyy_yyy, g_0_y_xyyyyyy_yyz, g_0_y_xyyyyyy_yzz, g_0_y_xyyyyyy_zzz, g_0_y_yyyyyy_xxx, g_0_y_yyyyyy_xxxx, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxy, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xyy, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_yyy, g_0_y_yyyyyy_yyz, g_0_y_yyyyyy_yzz, g_0_y_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyy_xxx[k] = -g_0_y_yyyyyy_xxx[k] * ab_x + g_0_y_yyyyyy_xxxx[k];

                g_0_y_xyyyyyy_xxy[k] = -g_0_y_yyyyyy_xxy[k] * ab_x + g_0_y_yyyyyy_xxxy[k];

                g_0_y_xyyyyyy_xxz[k] = -g_0_y_yyyyyy_xxz[k] * ab_x + g_0_y_yyyyyy_xxxz[k];

                g_0_y_xyyyyyy_xyy[k] = -g_0_y_yyyyyy_xyy[k] * ab_x + g_0_y_yyyyyy_xxyy[k];

                g_0_y_xyyyyyy_xyz[k] = -g_0_y_yyyyyy_xyz[k] * ab_x + g_0_y_yyyyyy_xxyz[k];

                g_0_y_xyyyyyy_xzz[k] = -g_0_y_yyyyyy_xzz[k] * ab_x + g_0_y_yyyyyy_xxzz[k];

                g_0_y_xyyyyyy_yyy[k] = -g_0_y_yyyyyy_yyy[k] * ab_x + g_0_y_yyyyyy_xyyy[k];

                g_0_y_xyyyyyy_yyz[k] = -g_0_y_yyyyyy_yyz[k] * ab_x + g_0_y_yyyyyy_xyyz[k];

                g_0_y_xyyyyyy_yzz[k] = -g_0_y_yyyyyy_yzz[k] * ab_x + g_0_y_yyyyyy_xyzz[k];

                g_0_y_xyyyyyy_zzz[k] = -g_0_y_yyyyyy_zzz[k] * ab_x + g_0_y_yyyyyy_xzzz[k];
            }

            /// Set up 580-590 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 589 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyz_xxx, g_0_y_xyyyyyz_xxy, g_0_y_xyyyyyz_xxz, g_0_y_xyyyyyz_xyy, g_0_y_xyyyyyz_xyz, g_0_y_xyyyyyz_xzz, g_0_y_xyyyyyz_yyy, g_0_y_xyyyyyz_yyz, g_0_y_xyyyyyz_yzz, g_0_y_xyyyyyz_zzz, g_0_y_yyyyyz_xxx, g_0_y_yyyyyz_xxxx, g_0_y_yyyyyz_xxxy, g_0_y_yyyyyz_xxxz, g_0_y_yyyyyz_xxy, g_0_y_yyyyyz_xxyy, g_0_y_yyyyyz_xxyz, g_0_y_yyyyyz_xxz, g_0_y_yyyyyz_xxzz, g_0_y_yyyyyz_xyy, g_0_y_yyyyyz_xyyy, g_0_y_yyyyyz_xyyz, g_0_y_yyyyyz_xyz, g_0_y_yyyyyz_xyzz, g_0_y_yyyyyz_xzz, g_0_y_yyyyyz_xzzz, g_0_y_yyyyyz_yyy, g_0_y_yyyyyz_yyz, g_0_y_yyyyyz_yzz, g_0_y_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyz_xxx[k] = -g_0_y_yyyyyz_xxx[k] * ab_x + g_0_y_yyyyyz_xxxx[k];

                g_0_y_xyyyyyz_xxy[k] = -g_0_y_yyyyyz_xxy[k] * ab_x + g_0_y_yyyyyz_xxxy[k];

                g_0_y_xyyyyyz_xxz[k] = -g_0_y_yyyyyz_xxz[k] * ab_x + g_0_y_yyyyyz_xxxz[k];

                g_0_y_xyyyyyz_xyy[k] = -g_0_y_yyyyyz_xyy[k] * ab_x + g_0_y_yyyyyz_xxyy[k];

                g_0_y_xyyyyyz_xyz[k] = -g_0_y_yyyyyz_xyz[k] * ab_x + g_0_y_yyyyyz_xxyz[k];

                g_0_y_xyyyyyz_xzz[k] = -g_0_y_yyyyyz_xzz[k] * ab_x + g_0_y_yyyyyz_xxzz[k];

                g_0_y_xyyyyyz_yyy[k] = -g_0_y_yyyyyz_yyy[k] * ab_x + g_0_y_yyyyyz_xyyy[k];

                g_0_y_xyyyyyz_yyz[k] = -g_0_y_yyyyyz_yyz[k] * ab_x + g_0_y_yyyyyz_xyyz[k];

                g_0_y_xyyyyyz_yzz[k] = -g_0_y_yyyyyz_yzz[k] * ab_x + g_0_y_yyyyyz_xyzz[k];

                g_0_y_xyyyyyz_zzz[k] = -g_0_y_yyyyyz_zzz[k] * ab_x + g_0_y_yyyyyz_xzzz[k];
            }

            /// Set up 590-600 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyzz_xxx, g_0_y_xyyyyzz_xxy, g_0_y_xyyyyzz_xxz, g_0_y_xyyyyzz_xyy, g_0_y_xyyyyzz_xyz, g_0_y_xyyyyzz_xzz, g_0_y_xyyyyzz_yyy, g_0_y_xyyyyzz_yyz, g_0_y_xyyyyzz_yzz, g_0_y_xyyyyzz_zzz, g_0_y_yyyyzz_xxx, g_0_y_yyyyzz_xxxx, g_0_y_yyyyzz_xxxy, g_0_y_yyyyzz_xxxz, g_0_y_yyyyzz_xxy, g_0_y_yyyyzz_xxyy, g_0_y_yyyyzz_xxyz, g_0_y_yyyyzz_xxz, g_0_y_yyyyzz_xxzz, g_0_y_yyyyzz_xyy, g_0_y_yyyyzz_xyyy, g_0_y_yyyyzz_xyyz, g_0_y_yyyyzz_xyz, g_0_y_yyyyzz_xyzz, g_0_y_yyyyzz_xzz, g_0_y_yyyyzz_xzzz, g_0_y_yyyyzz_yyy, g_0_y_yyyyzz_yyz, g_0_y_yyyyzz_yzz, g_0_y_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyzz_xxx[k] = -g_0_y_yyyyzz_xxx[k] * ab_x + g_0_y_yyyyzz_xxxx[k];

                g_0_y_xyyyyzz_xxy[k] = -g_0_y_yyyyzz_xxy[k] * ab_x + g_0_y_yyyyzz_xxxy[k];

                g_0_y_xyyyyzz_xxz[k] = -g_0_y_yyyyzz_xxz[k] * ab_x + g_0_y_yyyyzz_xxxz[k];

                g_0_y_xyyyyzz_xyy[k] = -g_0_y_yyyyzz_xyy[k] * ab_x + g_0_y_yyyyzz_xxyy[k];

                g_0_y_xyyyyzz_xyz[k] = -g_0_y_yyyyzz_xyz[k] * ab_x + g_0_y_yyyyzz_xxyz[k];

                g_0_y_xyyyyzz_xzz[k] = -g_0_y_yyyyzz_xzz[k] * ab_x + g_0_y_yyyyzz_xxzz[k];

                g_0_y_xyyyyzz_yyy[k] = -g_0_y_yyyyzz_yyy[k] * ab_x + g_0_y_yyyyzz_xyyy[k];

                g_0_y_xyyyyzz_yyz[k] = -g_0_y_yyyyzz_yyz[k] * ab_x + g_0_y_yyyyzz_xyyz[k];

                g_0_y_xyyyyzz_yzz[k] = -g_0_y_yyyyzz_yzz[k] * ab_x + g_0_y_yyyyzz_xyzz[k];

                g_0_y_xyyyyzz_zzz[k] = -g_0_y_yyyyzz_zzz[k] * ab_x + g_0_y_yyyyzz_xzzz[k];
            }

            /// Set up 600-610 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 609 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzzz_xxx, g_0_y_xyyyzzz_xxy, g_0_y_xyyyzzz_xxz, g_0_y_xyyyzzz_xyy, g_0_y_xyyyzzz_xyz, g_0_y_xyyyzzz_xzz, g_0_y_xyyyzzz_yyy, g_0_y_xyyyzzz_yyz, g_0_y_xyyyzzz_yzz, g_0_y_xyyyzzz_zzz, g_0_y_yyyzzz_xxx, g_0_y_yyyzzz_xxxx, g_0_y_yyyzzz_xxxy, g_0_y_yyyzzz_xxxz, g_0_y_yyyzzz_xxy, g_0_y_yyyzzz_xxyy, g_0_y_yyyzzz_xxyz, g_0_y_yyyzzz_xxz, g_0_y_yyyzzz_xxzz, g_0_y_yyyzzz_xyy, g_0_y_yyyzzz_xyyy, g_0_y_yyyzzz_xyyz, g_0_y_yyyzzz_xyz, g_0_y_yyyzzz_xyzz, g_0_y_yyyzzz_xzz, g_0_y_yyyzzz_xzzz, g_0_y_yyyzzz_yyy, g_0_y_yyyzzz_yyz, g_0_y_yyyzzz_yzz, g_0_y_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzzz_xxx[k] = -g_0_y_yyyzzz_xxx[k] * ab_x + g_0_y_yyyzzz_xxxx[k];

                g_0_y_xyyyzzz_xxy[k] = -g_0_y_yyyzzz_xxy[k] * ab_x + g_0_y_yyyzzz_xxxy[k];

                g_0_y_xyyyzzz_xxz[k] = -g_0_y_yyyzzz_xxz[k] * ab_x + g_0_y_yyyzzz_xxxz[k];

                g_0_y_xyyyzzz_xyy[k] = -g_0_y_yyyzzz_xyy[k] * ab_x + g_0_y_yyyzzz_xxyy[k];

                g_0_y_xyyyzzz_xyz[k] = -g_0_y_yyyzzz_xyz[k] * ab_x + g_0_y_yyyzzz_xxyz[k];

                g_0_y_xyyyzzz_xzz[k] = -g_0_y_yyyzzz_xzz[k] * ab_x + g_0_y_yyyzzz_xxzz[k];

                g_0_y_xyyyzzz_yyy[k] = -g_0_y_yyyzzz_yyy[k] * ab_x + g_0_y_yyyzzz_xyyy[k];

                g_0_y_xyyyzzz_yyz[k] = -g_0_y_yyyzzz_yyz[k] * ab_x + g_0_y_yyyzzz_xyyz[k];

                g_0_y_xyyyzzz_yzz[k] = -g_0_y_yyyzzz_yzz[k] * ab_x + g_0_y_yyyzzz_xyzz[k];

                g_0_y_xyyyzzz_zzz[k] = -g_0_y_yyyzzz_zzz[k] * ab_x + g_0_y_yyyzzz_xzzz[k];
            }

            /// Set up 610-620 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 619 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzzz_xxx, g_0_y_xyyzzzz_xxy, g_0_y_xyyzzzz_xxz, g_0_y_xyyzzzz_xyy, g_0_y_xyyzzzz_xyz, g_0_y_xyyzzzz_xzz, g_0_y_xyyzzzz_yyy, g_0_y_xyyzzzz_yyz, g_0_y_xyyzzzz_yzz, g_0_y_xyyzzzz_zzz, g_0_y_yyzzzz_xxx, g_0_y_yyzzzz_xxxx, g_0_y_yyzzzz_xxxy, g_0_y_yyzzzz_xxxz, g_0_y_yyzzzz_xxy, g_0_y_yyzzzz_xxyy, g_0_y_yyzzzz_xxyz, g_0_y_yyzzzz_xxz, g_0_y_yyzzzz_xxzz, g_0_y_yyzzzz_xyy, g_0_y_yyzzzz_xyyy, g_0_y_yyzzzz_xyyz, g_0_y_yyzzzz_xyz, g_0_y_yyzzzz_xyzz, g_0_y_yyzzzz_xzz, g_0_y_yyzzzz_xzzz, g_0_y_yyzzzz_yyy, g_0_y_yyzzzz_yyz, g_0_y_yyzzzz_yzz, g_0_y_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzzz_xxx[k] = -g_0_y_yyzzzz_xxx[k] * ab_x + g_0_y_yyzzzz_xxxx[k];

                g_0_y_xyyzzzz_xxy[k] = -g_0_y_yyzzzz_xxy[k] * ab_x + g_0_y_yyzzzz_xxxy[k];

                g_0_y_xyyzzzz_xxz[k] = -g_0_y_yyzzzz_xxz[k] * ab_x + g_0_y_yyzzzz_xxxz[k];

                g_0_y_xyyzzzz_xyy[k] = -g_0_y_yyzzzz_xyy[k] * ab_x + g_0_y_yyzzzz_xxyy[k];

                g_0_y_xyyzzzz_xyz[k] = -g_0_y_yyzzzz_xyz[k] * ab_x + g_0_y_yyzzzz_xxyz[k];

                g_0_y_xyyzzzz_xzz[k] = -g_0_y_yyzzzz_xzz[k] * ab_x + g_0_y_yyzzzz_xxzz[k];

                g_0_y_xyyzzzz_yyy[k] = -g_0_y_yyzzzz_yyy[k] * ab_x + g_0_y_yyzzzz_xyyy[k];

                g_0_y_xyyzzzz_yyz[k] = -g_0_y_yyzzzz_yyz[k] * ab_x + g_0_y_yyzzzz_xyyz[k];

                g_0_y_xyyzzzz_yzz[k] = -g_0_y_yyzzzz_yzz[k] * ab_x + g_0_y_yyzzzz_xyzz[k];

                g_0_y_xyyzzzz_zzz[k] = -g_0_y_yyzzzz_zzz[k] * ab_x + g_0_y_yyzzzz_xzzz[k];
            }

            /// Set up 620-630 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzzz_xxx, g_0_y_xyzzzzz_xxy, g_0_y_xyzzzzz_xxz, g_0_y_xyzzzzz_xyy, g_0_y_xyzzzzz_xyz, g_0_y_xyzzzzz_xzz, g_0_y_xyzzzzz_yyy, g_0_y_xyzzzzz_yyz, g_0_y_xyzzzzz_yzz, g_0_y_xyzzzzz_zzz, g_0_y_yzzzzz_xxx, g_0_y_yzzzzz_xxxx, g_0_y_yzzzzz_xxxy, g_0_y_yzzzzz_xxxz, g_0_y_yzzzzz_xxy, g_0_y_yzzzzz_xxyy, g_0_y_yzzzzz_xxyz, g_0_y_yzzzzz_xxz, g_0_y_yzzzzz_xxzz, g_0_y_yzzzzz_xyy, g_0_y_yzzzzz_xyyy, g_0_y_yzzzzz_xyyz, g_0_y_yzzzzz_xyz, g_0_y_yzzzzz_xyzz, g_0_y_yzzzzz_xzz, g_0_y_yzzzzz_xzzz, g_0_y_yzzzzz_yyy, g_0_y_yzzzzz_yyz, g_0_y_yzzzzz_yzz, g_0_y_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzzz_xxx[k] = -g_0_y_yzzzzz_xxx[k] * ab_x + g_0_y_yzzzzz_xxxx[k];

                g_0_y_xyzzzzz_xxy[k] = -g_0_y_yzzzzz_xxy[k] * ab_x + g_0_y_yzzzzz_xxxy[k];

                g_0_y_xyzzzzz_xxz[k] = -g_0_y_yzzzzz_xxz[k] * ab_x + g_0_y_yzzzzz_xxxz[k];

                g_0_y_xyzzzzz_xyy[k] = -g_0_y_yzzzzz_xyy[k] * ab_x + g_0_y_yzzzzz_xxyy[k];

                g_0_y_xyzzzzz_xyz[k] = -g_0_y_yzzzzz_xyz[k] * ab_x + g_0_y_yzzzzz_xxyz[k];

                g_0_y_xyzzzzz_xzz[k] = -g_0_y_yzzzzz_xzz[k] * ab_x + g_0_y_yzzzzz_xxzz[k];

                g_0_y_xyzzzzz_yyy[k] = -g_0_y_yzzzzz_yyy[k] * ab_x + g_0_y_yzzzzz_xyyy[k];

                g_0_y_xyzzzzz_yyz[k] = -g_0_y_yzzzzz_yyz[k] * ab_x + g_0_y_yzzzzz_xyyz[k];

                g_0_y_xyzzzzz_yzz[k] = -g_0_y_yzzzzz_yzz[k] * ab_x + g_0_y_yzzzzz_xyzz[k];

                g_0_y_xyzzzzz_zzz[k] = -g_0_y_yzzzzz_zzz[k] * ab_x + g_0_y_yzzzzz_xzzz[k];
            }

            /// Set up 630-640 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 639 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzzz_xxx, g_0_y_xzzzzzz_xxy, g_0_y_xzzzzzz_xxz, g_0_y_xzzzzzz_xyy, g_0_y_xzzzzzz_xyz, g_0_y_xzzzzzz_xzz, g_0_y_xzzzzzz_yyy, g_0_y_xzzzzzz_yyz, g_0_y_xzzzzzz_yzz, g_0_y_xzzzzzz_zzz, g_0_y_zzzzzz_xxx, g_0_y_zzzzzz_xxxx, g_0_y_zzzzzz_xxxy, g_0_y_zzzzzz_xxxz, g_0_y_zzzzzz_xxy, g_0_y_zzzzzz_xxyy, g_0_y_zzzzzz_xxyz, g_0_y_zzzzzz_xxz, g_0_y_zzzzzz_xxzz, g_0_y_zzzzzz_xyy, g_0_y_zzzzzz_xyyy, g_0_y_zzzzzz_xyyz, g_0_y_zzzzzz_xyz, g_0_y_zzzzzz_xyzz, g_0_y_zzzzzz_xzz, g_0_y_zzzzzz_xzzz, g_0_y_zzzzzz_yyy, g_0_y_zzzzzz_yyz, g_0_y_zzzzzz_yzz, g_0_y_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzzz_xxx[k] = -g_0_y_zzzzzz_xxx[k] * ab_x + g_0_y_zzzzzz_xxxx[k];

                g_0_y_xzzzzzz_xxy[k] = -g_0_y_zzzzzz_xxy[k] * ab_x + g_0_y_zzzzzz_xxxy[k];

                g_0_y_xzzzzzz_xxz[k] = -g_0_y_zzzzzz_xxz[k] * ab_x + g_0_y_zzzzzz_xxxz[k];

                g_0_y_xzzzzzz_xyy[k] = -g_0_y_zzzzzz_xyy[k] * ab_x + g_0_y_zzzzzz_xxyy[k];

                g_0_y_xzzzzzz_xyz[k] = -g_0_y_zzzzzz_xyz[k] * ab_x + g_0_y_zzzzzz_xxyz[k];

                g_0_y_xzzzzzz_xzz[k] = -g_0_y_zzzzzz_xzz[k] * ab_x + g_0_y_zzzzzz_xxzz[k];

                g_0_y_xzzzzzz_yyy[k] = -g_0_y_zzzzzz_yyy[k] * ab_x + g_0_y_zzzzzz_xyyy[k];

                g_0_y_xzzzzzz_yyz[k] = -g_0_y_zzzzzz_yyz[k] * ab_x + g_0_y_zzzzzz_xyyz[k];

                g_0_y_xzzzzzz_yzz[k] = -g_0_y_zzzzzz_yzz[k] * ab_x + g_0_y_zzzzzz_xyzz[k];

                g_0_y_xzzzzzz_zzz[k] = -g_0_y_zzzzzz_zzz[k] * ab_x + g_0_y_zzzzzz_xzzz[k];
            }

            /// Set up 640-650 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 649 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xxx, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxy, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxz, g_0_y_yyyyyy_xyy, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xzz, g_0_y_yyyyyy_yyy, g_0_y_yyyyyy_yyyy, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_zzz, g_0_y_yyyyyyy_xxx, g_0_y_yyyyyyy_xxy, g_0_y_yyyyyyy_xxz, g_0_y_yyyyyyy_xyy, g_0_y_yyyyyyy_xyz, g_0_y_yyyyyyy_xzz, g_0_y_yyyyyyy_yyy, g_0_y_yyyyyyy_yyz, g_0_y_yyyyyyy_yzz, g_0_y_yyyyyyy_zzz, g_yyyyyy_xxx, g_yyyyyy_xxy, g_yyyyyy_xxz, g_yyyyyy_xyy, g_yyyyyy_xyz, g_yyyyyy_xzz, g_yyyyyy_yyy, g_yyyyyy_yyz, g_yyyyyy_yzz, g_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyy_xxx[k] = g_yyyyyy_xxx[k] - g_0_y_yyyyyy_xxx[k] * ab_y + g_0_y_yyyyyy_xxxy[k];

                g_0_y_yyyyyyy_xxy[k] = g_yyyyyy_xxy[k] - g_0_y_yyyyyy_xxy[k] * ab_y + g_0_y_yyyyyy_xxyy[k];

                g_0_y_yyyyyyy_xxz[k] = g_yyyyyy_xxz[k] - g_0_y_yyyyyy_xxz[k] * ab_y + g_0_y_yyyyyy_xxyz[k];

                g_0_y_yyyyyyy_xyy[k] = g_yyyyyy_xyy[k] - g_0_y_yyyyyy_xyy[k] * ab_y + g_0_y_yyyyyy_xyyy[k];

                g_0_y_yyyyyyy_xyz[k] = g_yyyyyy_xyz[k] - g_0_y_yyyyyy_xyz[k] * ab_y + g_0_y_yyyyyy_xyyz[k];

                g_0_y_yyyyyyy_xzz[k] = g_yyyyyy_xzz[k] - g_0_y_yyyyyy_xzz[k] * ab_y + g_0_y_yyyyyy_xyzz[k];

                g_0_y_yyyyyyy_yyy[k] = g_yyyyyy_yyy[k] - g_0_y_yyyyyy_yyy[k] * ab_y + g_0_y_yyyyyy_yyyy[k];

                g_0_y_yyyyyyy_yyz[k] = g_yyyyyy_yyz[k] - g_0_y_yyyyyy_yyz[k] * ab_y + g_0_y_yyyyyy_yyyz[k];

                g_0_y_yyyyyyy_yzz[k] = g_yyyyyy_yzz[k] - g_0_y_yyyyyy_yzz[k] * ab_y + g_0_y_yyyyyy_yyzz[k];

                g_0_y_yyyyyyy_zzz[k] = g_yyyyyy_zzz[k] - g_0_y_yyyyyy_zzz[k] * ab_y + g_0_y_yyyyyy_yzzz[k];
            }

            /// Set up 650-660 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xxx, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxy, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xyy, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_yyy, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_zzz, g_0_y_yyyyyy_zzzz, g_0_y_yyyyyyz_xxx, g_0_y_yyyyyyz_xxy, g_0_y_yyyyyyz_xxz, g_0_y_yyyyyyz_xyy, g_0_y_yyyyyyz_xyz, g_0_y_yyyyyyz_xzz, g_0_y_yyyyyyz_yyy, g_0_y_yyyyyyz_yyz, g_0_y_yyyyyyz_yzz, g_0_y_yyyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyz_xxx[k] = -g_0_y_yyyyyy_xxx[k] * ab_z + g_0_y_yyyyyy_xxxz[k];

                g_0_y_yyyyyyz_xxy[k] = -g_0_y_yyyyyy_xxy[k] * ab_z + g_0_y_yyyyyy_xxyz[k];

                g_0_y_yyyyyyz_xxz[k] = -g_0_y_yyyyyy_xxz[k] * ab_z + g_0_y_yyyyyy_xxzz[k];

                g_0_y_yyyyyyz_xyy[k] = -g_0_y_yyyyyy_xyy[k] * ab_z + g_0_y_yyyyyy_xyyz[k];

                g_0_y_yyyyyyz_xyz[k] = -g_0_y_yyyyyy_xyz[k] * ab_z + g_0_y_yyyyyy_xyzz[k];

                g_0_y_yyyyyyz_xzz[k] = -g_0_y_yyyyyy_xzz[k] * ab_z + g_0_y_yyyyyy_xzzz[k];

                g_0_y_yyyyyyz_yyy[k] = -g_0_y_yyyyyy_yyy[k] * ab_z + g_0_y_yyyyyy_yyyz[k];

                g_0_y_yyyyyyz_yyz[k] = -g_0_y_yyyyyy_yyz[k] * ab_z + g_0_y_yyyyyy_yyzz[k];

                g_0_y_yyyyyyz_yzz[k] = -g_0_y_yyyyyy_yzz[k] * ab_z + g_0_y_yyyyyy_yzzz[k];

                g_0_y_yyyyyyz_zzz[k] = -g_0_y_yyyyyy_zzz[k] * ab_z + g_0_y_yyyyyy_zzzz[k];
            }

            /// Set up 660-670 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 669 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyz_xxx, g_0_y_yyyyyz_xxxz, g_0_y_yyyyyz_xxy, g_0_y_yyyyyz_xxyz, g_0_y_yyyyyz_xxz, g_0_y_yyyyyz_xxzz, g_0_y_yyyyyz_xyy, g_0_y_yyyyyz_xyyz, g_0_y_yyyyyz_xyz, g_0_y_yyyyyz_xyzz, g_0_y_yyyyyz_xzz, g_0_y_yyyyyz_xzzz, g_0_y_yyyyyz_yyy, g_0_y_yyyyyz_yyyz, g_0_y_yyyyyz_yyz, g_0_y_yyyyyz_yyzz, g_0_y_yyyyyz_yzz, g_0_y_yyyyyz_yzzz, g_0_y_yyyyyz_zzz, g_0_y_yyyyyz_zzzz, g_0_y_yyyyyzz_xxx, g_0_y_yyyyyzz_xxy, g_0_y_yyyyyzz_xxz, g_0_y_yyyyyzz_xyy, g_0_y_yyyyyzz_xyz, g_0_y_yyyyyzz_xzz, g_0_y_yyyyyzz_yyy, g_0_y_yyyyyzz_yyz, g_0_y_yyyyyzz_yzz, g_0_y_yyyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyzz_xxx[k] = -g_0_y_yyyyyz_xxx[k] * ab_z + g_0_y_yyyyyz_xxxz[k];

                g_0_y_yyyyyzz_xxy[k] = -g_0_y_yyyyyz_xxy[k] * ab_z + g_0_y_yyyyyz_xxyz[k];

                g_0_y_yyyyyzz_xxz[k] = -g_0_y_yyyyyz_xxz[k] * ab_z + g_0_y_yyyyyz_xxzz[k];

                g_0_y_yyyyyzz_xyy[k] = -g_0_y_yyyyyz_xyy[k] * ab_z + g_0_y_yyyyyz_xyyz[k];

                g_0_y_yyyyyzz_xyz[k] = -g_0_y_yyyyyz_xyz[k] * ab_z + g_0_y_yyyyyz_xyzz[k];

                g_0_y_yyyyyzz_xzz[k] = -g_0_y_yyyyyz_xzz[k] * ab_z + g_0_y_yyyyyz_xzzz[k];

                g_0_y_yyyyyzz_yyy[k] = -g_0_y_yyyyyz_yyy[k] * ab_z + g_0_y_yyyyyz_yyyz[k];

                g_0_y_yyyyyzz_yyz[k] = -g_0_y_yyyyyz_yyz[k] * ab_z + g_0_y_yyyyyz_yyzz[k];

                g_0_y_yyyyyzz_yzz[k] = -g_0_y_yyyyyz_yzz[k] * ab_z + g_0_y_yyyyyz_yzzz[k];

                g_0_y_yyyyyzz_zzz[k] = -g_0_y_yyyyyz_zzz[k] * ab_z + g_0_y_yyyyyz_zzzz[k];
            }

            /// Set up 670-680 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 679 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyzz_xxx, g_0_y_yyyyzz_xxxz, g_0_y_yyyyzz_xxy, g_0_y_yyyyzz_xxyz, g_0_y_yyyyzz_xxz, g_0_y_yyyyzz_xxzz, g_0_y_yyyyzz_xyy, g_0_y_yyyyzz_xyyz, g_0_y_yyyyzz_xyz, g_0_y_yyyyzz_xyzz, g_0_y_yyyyzz_xzz, g_0_y_yyyyzz_xzzz, g_0_y_yyyyzz_yyy, g_0_y_yyyyzz_yyyz, g_0_y_yyyyzz_yyz, g_0_y_yyyyzz_yyzz, g_0_y_yyyyzz_yzz, g_0_y_yyyyzz_yzzz, g_0_y_yyyyzz_zzz, g_0_y_yyyyzz_zzzz, g_0_y_yyyyzzz_xxx, g_0_y_yyyyzzz_xxy, g_0_y_yyyyzzz_xxz, g_0_y_yyyyzzz_xyy, g_0_y_yyyyzzz_xyz, g_0_y_yyyyzzz_xzz, g_0_y_yyyyzzz_yyy, g_0_y_yyyyzzz_yyz, g_0_y_yyyyzzz_yzz, g_0_y_yyyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzzz_xxx[k] = -g_0_y_yyyyzz_xxx[k] * ab_z + g_0_y_yyyyzz_xxxz[k];

                g_0_y_yyyyzzz_xxy[k] = -g_0_y_yyyyzz_xxy[k] * ab_z + g_0_y_yyyyzz_xxyz[k];

                g_0_y_yyyyzzz_xxz[k] = -g_0_y_yyyyzz_xxz[k] * ab_z + g_0_y_yyyyzz_xxzz[k];

                g_0_y_yyyyzzz_xyy[k] = -g_0_y_yyyyzz_xyy[k] * ab_z + g_0_y_yyyyzz_xyyz[k];

                g_0_y_yyyyzzz_xyz[k] = -g_0_y_yyyyzz_xyz[k] * ab_z + g_0_y_yyyyzz_xyzz[k];

                g_0_y_yyyyzzz_xzz[k] = -g_0_y_yyyyzz_xzz[k] * ab_z + g_0_y_yyyyzz_xzzz[k];

                g_0_y_yyyyzzz_yyy[k] = -g_0_y_yyyyzz_yyy[k] * ab_z + g_0_y_yyyyzz_yyyz[k];

                g_0_y_yyyyzzz_yyz[k] = -g_0_y_yyyyzz_yyz[k] * ab_z + g_0_y_yyyyzz_yyzz[k];

                g_0_y_yyyyzzz_yzz[k] = -g_0_y_yyyyzz_yzz[k] * ab_z + g_0_y_yyyyzz_yzzz[k];

                g_0_y_yyyyzzz_zzz[k] = -g_0_y_yyyyzz_zzz[k] * ab_z + g_0_y_yyyyzz_zzzz[k];
            }

            /// Set up 680-690 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzzz_xxx, g_0_y_yyyzzz_xxxz, g_0_y_yyyzzz_xxy, g_0_y_yyyzzz_xxyz, g_0_y_yyyzzz_xxz, g_0_y_yyyzzz_xxzz, g_0_y_yyyzzz_xyy, g_0_y_yyyzzz_xyyz, g_0_y_yyyzzz_xyz, g_0_y_yyyzzz_xyzz, g_0_y_yyyzzz_xzz, g_0_y_yyyzzz_xzzz, g_0_y_yyyzzz_yyy, g_0_y_yyyzzz_yyyz, g_0_y_yyyzzz_yyz, g_0_y_yyyzzz_yyzz, g_0_y_yyyzzz_yzz, g_0_y_yyyzzz_yzzz, g_0_y_yyyzzz_zzz, g_0_y_yyyzzz_zzzz, g_0_y_yyyzzzz_xxx, g_0_y_yyyzzzz_xxy, g_0_y_yyyzzzz_xxz, g_0_y_yyyzzzz_xyy, g_0_y_yyyzzzz_xyz, g_0_y_yyyzzzz_xzz, g_0_y_yyyzzzz_yyy, g_0_y_yyyzzzz_yyz, g_0_y_yyyzzzz_yzz, g_0_y_yyyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzzz_xxx[k] = -g_0_y_yyyzzz_xxx[k] * ab_z + g_0_y_yyyzzz_xxxz[k];

                g_0_y_yyyzzzz_xxy[k] = -g_0_y_yyyzzz_xxy[k] * ab_z + g_0_y_yyyzzz_xxyz[k];

                g_0_y_yyyzzzz_xxz[k] = -g_0_y_yyyzzz_xxz[k] * ab_z + g_0_y_yyyzzz_xxzz[k];

                g_0_y_yyyzzzz_xyy[k] = -g_0_y_yyyzzz_xyy[k] * ab_z + g_0_y_yyyzzz_xyyz[k];

                g_0_y_yyyzzzz_xyz[k] = -g_0_y_yyyzzz_xyz[k] * ab_z + g_0_y_yyyzzz_xyzz[k];

                g_0_y_yyyzzzz_xzz[k] = -g_0_y_yyyzzz_xzz[k] * ab_z + g_0_y_yyyzzz_xzzz[k];

                g_0_y_yyyzzzz_yyy[k] = -g_0_y_yyyzzz_yyy[k] * ab_z + g_0_y_yyyzzz_yyyz[k];

                g_0_y_yyyzzzz_yyz[k] = -g_0_y_yyyzzz_yyz[k] * ab_z + g_0_y_yyyzzz_yyzz[k];

                g_0_y_yyyzzzz_yzz[k] = -g_0_y_yyyzzz_yzz[k] * ab_z + g_0_y_yyyzzz_yzzz[k];

                g_0_y_yyyzzzz_zzz[k] = -g_0_y_yyyzzz_zzz[k] * ab_z + g_0_y_yyyzzz_zzzz[k];
            }

            /// Set up 690-700 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzzz_xxx, g_0_y_yyzzzz_xxxz, g_0_y_yyzzzz_xxy, g_0_y_yyzzzz_xxyz, g_0_y_yyzzzz_xxz, g_0_y_yyzzzz_xxzz, g_0_y_yyzzzz_xyy, g_0_y_yyzzzz_xyyz, g_0_y_yyzzzz_xyz, g_0_y_yyzzzz_xyzz, g_0_y_yyzzzz_xzz, g_0_y_yyzzzz_xzzz, g_0_y_yyzzzz_yyy, g_0_y_yyzzzz_yyyz, g_0_y_yyzzzz_yyz, g_0_y_yyzzzz_yyzz, g_0_y_yyzzzz_yzz, g_0_y_yyzzzz_yzzz, g_0_y_yyzzzz_zzz, g_0_y_yyzzzz_zzzz, g_0_y_yyzzzzz_xxx, g_0_y_yyzzzzz_xxy, g_0_y_yyzzzzz_xxz, g_0_y_yyzzzzz_xyy, g_0_y_yyzzzzz_xyz, g_0_y_yyzzzzz_xzz, g_0_y_yyzzzzz_yyy, g_0_y_yyzzzzz_yyz, g_0_y_yyzzzzz_yzz, g_0_y_yyzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzzz_xxx[k] = -g_0_y_yyzzzz_xxx[k] * ab_z + g_0_y_yyzzzz_xxxz[k];

                g_0_y_yyzzzzz_xxy[k] = -g_0_y_yyzzzz_xxy[k] * ab_z + g_0_y_yyzzzz_xxyz[k];

                g_0_y_yyzzzzz_xxz[k] = -g_0_y_yyzzzz_xxz[k] * ab_z + g_0_y_yyzzzz_xxzz[k];

                g_0_y_yyzzzzz_xyy[k] = -g_0_y_yyzzzz_xyy[k] * ab_z + g_0_y_yyzzzz_xyyz[k];

                g_0_y_yyzzzzz_xyz[k] = -g_0_y_yyzzzz_xyz[k] * ab_z + g_0_y_yyzzzz_xyzz[k];

                g_0_y_yyzzzzz_xzz[k] = -g_0_y_yyzzzz_xzz[k] * ab_z + g_0_y_yyzzzz_xzzz[k];

                g_0_y_yyzzzzz_yyy[k] = -g_0_y_yyzzzz_yyy[k] * ab_z + g_0_y_yyzzzz_yyyz[k];

                g_0_y_yyzzzzz_yyz[k] = -g_0_y_yyzzzz_yyz[k] * ab_z + g_0_y_yyzzzz_yyzz[k];

                g_0_y_yyzzzzz_yzz[k] = -g_0_y_yyzzzz_yzz[k] * ab_z + g_0_y_yyzzzz_yzzz[k];

                g_0_y_yyzzzzz_zzz[k] = -g_0_y_yyzzzz_zzz[k] * ab_z + g_0_y_yyzzzz_zzzz[k];
            }

            /// Set up 700-710 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 709 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzzz_xxx, g_0_y_yzzzzz_xxxz, g_0_y_yzzzzz_xxy, g_0_y_yzzzzz_xxyz, g_0_y_yzzzzz_xxz, g_0_y_yzzzzz_xxzz, g_0_y_yzzzzz_xyy, g_0_y_yzzzzz_xyyz, g_0_y_yzzzzz_xyz, g_0_y_yzzzzz_xyzz, g_0_y_yzzzzz_xzz, g_0_y_yzzzzz_xzzz, g_0_y_yzzzzz_yyy, g_0_y_yzzzzz_yyyz, g_0_y_yzzzzz_yyz, g_0_y_yzzzzz_yyzz, g_0_y_yzzzzz_yzz, g_0_y_yzzzzz_yzzz, g_0_y_yzzzzz_zzz, g_0_y_yzzzzz_zzzz, g_0_y_yzzzzzz_xxx, g_0_y_yzzzzzz_xxy, g_0_y_yzzzzzz_xxz, g_0_y_yzzzzzz_xyy, g_0_y_yzzzzzz_xyz, g_0_y_yzzzzzz_xzz, g_0_y_yzzzzzz_yyy, g_0_y_yzzzzzz_yyz, g_0_y_yzzzzzz_yzz, g_0_y_yzzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzzz_xxx[k] = -g_0_y_yzzzzz_xxx[k] * ab_z + g_0_y_yzzzzz_xxxz[k];

                g_0_y_yzzzzzz_xxy[k] = -g_0_y_yzzzzz_xxy[k] * ab_z + g_0_y_yzzzzz_xxyz[k];

                g_0_y_yzzzzzz_xxz[k] = -g_0_y_yzzzzz_xxz[k] * ab_z + g_0_y_yzzzzz_xxzz[k];

                g_0_y_yzzzzzz_xyy[k] = -g_0_y_yzzzzz_xyy[k] * ab_z + g_0_y_yzzzzz_xyyz[k];

                g_0_y_yzzzzzz_xyz[k] = -g_0_y_yzzzzz_xyz[k] * ab_z + g_0_y_yzzzzz_xyzz[k];

                g_0_y_yzzzzzz_xzz[k] = -g_0_y_yzzzzz_xzz[k] * ab_z + g_0_y_yzzzzz_xzzz[k];

                g_0_y_yzzzzzz_yyy[k] = -g_0_y_yzzzzz_yyy[k] * ab_z + g_0_y_yzzzzz_yyyz[k];

                g_0_y_yzzzzzz_yyz[k] = -g_0_y_yzzzzz_yyz[k] * ab_z + g_0_y_yzzzzz_yyzz[k];

                g_0_y_yzzzzzz_yzz[k] = -g_0_y_yzzzzz_yzz[k] * ab_z + g_0_y_yzzzzz_yzzz[k];

                g_0_y_yzzzzzz_zzz[k] = -g_0_y_yzzzzz_zzz[k] * ab_z + g_0_y_yzzzzz_zzzz[k];
            }

            /// Set up 710-720 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzzz_xxx, g_0_y_zzzzzz_xxxz, g_0_y_zzzzzz_xxy, g_0_y_zzzzzz_xxyz, g_0_y_zzzzzz_xxz, g_0_y_zzzzzz_xxzz, g_0_y_zzzzzz_xyy, g_0_y_zzzzzz_xyyz, g_0_y_zzzzzz_xyz, g_0_y_zzzzzz_xyzz, g_0_y_zzzzzz_xzz, g_0_y_zzzzzz_xzzz, g_0_y_zzzzzz_yyy, g_0_y_zzzzzz_yyyz, g_0_y_zzzzzz_yyz, g_0_y_zzzzzz_yyzz, g_0_y_zzzzzz_yzz, g_0_y_zzzzzz_yzzz, g_0_y_zzzzzz_zzz, g_0_y_zzzzzz_zzzz, g_0_y_zzzzzzz_xxx, g_0_y_zzzzzzz_xxy, g_0_y_zzzzzzz_xxz, g_0_y_zzzzzzz_xyy, g_0_y_zzzzzzz_xyz, g_0_y_zzzzzzz_xzz, g_0_y_zzzzzzz_yyy, g_0_y_zzzzzzz_yyz, g_0_y_zzzzzzz_yzz, g_0_y_zzzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzzz_xxx[k] = -g_0_y_zzzzzz_xxx[k] * ab_z + g_0_y_zzzzzz_xxxz[k];

                g_0_y_zzzzzzz_xxy[k] = -g_0_y_zzzzzz_xxy[k] * ab_z + g_0_y_zzzzzz_xxyz[k];

                g_0_y_zzzzzzz_xxz[k] = -g_0_y_zzzzzz_xxz[k] * ab_z + g_0_y_zzzzzz_xxzz[k];

                g_0_y_zzzzzzz_xyy[k] = -g_0_y_zzzzzz_xyy[k] * ab_z + g_0_y_zzzzzz_xyyz[k];

                g_0_y_zzzzzzz_xyz[k] = -g_0_y_zzzzzz_xyz[k] * ab_z + g_0_y_zzzzzz_xyzz[k];

                g_0_y_zzzzzzz_xzz[k] = -g_0_y_zzzzzz_xzz[k] * ab_z + g_0_y_zzzzzz_xzzz[k];

                g_0_y_zzzzzzz_yyy[k] = -g_0_y_zzzzzz_yyy[k] * ab_z + g_0_y_zzzzzz_yyyz[k];

                g_0_y_zzzzzzz_yyz[k] = -g_0_y_zzzzzz_yyz[k] * ab_z + g_0_y_zzzzzz_yyzz[k];

                g_0_y_zzzzzzz_yzz[k] = -g_0_y_zzzzzz_yzz[k] * ab_z + g_0_y_zzzzzz_yzzz[k];

                g_0_y_zzzzzzz_zzz[k] = -g_0_y_zzzzzz_zzz[k] * ab_z + g_0_y_zzzzzz_zzzz[k];
            }

            /// Set up 720-730 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxx_xxx = cbuffer.data(kf_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxy = cbuffer.data(kf_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxz = cbuffer.data(kf_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xyy = cbuffer.data(kf_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xyz = cbuffer.data(kf_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xzz = cbuffer.data(kf_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yyy = cbuffer.data(kf_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yyz = cbuffer.data(kf_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yzz = cbuffer.data(kf_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_zzz = cbuffer.data(kf_geom_01_off + 729 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxx_xxx, g_0_z_xxxxxx_xxxx, g_0_z_xxxxxx_xxxy, g_0_z_xxxxxx_xxxz, g_0_z_xxxxxx_xxy, g_0_z_xxxxxx_xxyy, g_0_z_xxxxxx_xxyz, g_0_z_xxxxxx_xxz, g_0_z_xxxxxx_xxzz, g_0_z_xxxxxx_xyy, g_0_z_xxxxxx_xyyy, g_0_z_xxxxxx_xyyz, g_0_z_xxxxxx_xyz, g_0_z_xxxxxx_xyzz, g_0_z_xxxxxx_xzz, g_0_z_xxxxxx_xzzz, g_0_z_xxxxxx_yyy, g_0_z_xxxxxx_yyz, g_0_z_xxxxxx_yzz, g_0_z_xxxxxx_zzz, g_0_z_xxxxxxx_xxx, g_0_z_xxxxxxx_xxy, g_0_z_xxxxxxx_xxz, g_0_z_xxxxxxx_xyy, g_0_z_xxxxxxx_xyz, g_0_z_xxxxxxx_xzz, g_0_z_xxxxxxx_yyy, g_0_z_xxxxxxx_yyz, g_0_z_xxxxxxx_yzz, g_0_z_xxxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxx_xxx[k] = -g_0_z_xxxxxx_xxx[k] * ab_x + g_0_z_xxxxxx_xxxx[k];

                g_0_z_xxxxxxx_xxy[k] = -g_0_z_xxxxxx_xxy[k] * ab_x + g_0_z_xxxxxx_xxxy[k];

                g_0_z_xxxxxxx_xxz[k] = -g_0_z_xxxxxx_xxz[k] * ab_x + g_0_z_xxxxxx_xxxz[k];

                g_0_z_xxxxxxx_xyy[k] = -g_0_z_xxxxxx_xyy[k] * ab_x + g_0_z_xxxxxx_xxyy[k];

                g_0_z_xxxxxxx_xyz[k] = -g_0_z_xxxxxx_xyz[k] * ab_x + g_0_z_xxxxxx_xxyz[k];

                g_0_z_xxxxxxx_xzz[k] = -g_0_z_xxxxxx_xzz[k] * ab_x + g_0_z_xxxxxx_xxzz[k];

                g_0_z_xxxxxxx_yyy[k] = -g_0_z_xxxxxx_yyy[k] * ab_x + g_0_z_xxxxxx_xyyy[k];

                g_0_z_xxxxxxx_yyz[k] = -g_0_z_xxxxxx_yyz[k] * ab_x + g_0_z_xxxxxx_xyyz[k];

                g_0_z_xxxxxxx_yzz[k] = -g_0_z_xxxxxx_yzz[k] * ab_x + g_0_z_xxxxxx_xyzz[k];

                g_0_z_xxxxxxx_zzz[k] = -g_0_z_xxxxxx_zzz[k] * ab_x + g_0_z_xxxxxx_xzzz[k];
            }

            /// Set up 730-740 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxy_xxx = cbuffer.data(kf_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxy = cbuffer.data(kf_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxz = cbuffer.data(kf_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xyy = cbuffer.data(kf_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xyz = cbuffer.data(kf_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xzz = cbuffer.data(kf_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yyy = cbuffer.data(kf_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yyz = cbuffer.data(kf_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yzz = cbuffer.data(kf_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_zzz = cbuffer.data(kf_geom_01_off + 739 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxy_xxx, g_0_z_xxxxxxy_xxy, g_0_z_xxxxxxy_xxz, g_0_z_xxxxxxy_xyy, g_0_z_xxxxxxy_xyz, g_0_z_xxxxxxy_xzz, g_0_z_xxxxxxy_yyy, g_0_z_xxxxxxy_yyz, g_0_z_xxxxxxy_yzz, g_0_z_xxxxxxy_zzz, g_0_z_xxxxxy_xxx, g_0_z_xxxxxy_xxxx, g_0_z_xxxxxy_xxxy, g_0_z_xxxxxy_xxxz, g_0_z_xxxxxy_xxy, g_0_z_xxxxxy_xxyy, g_0_z_xxxxxy_xxyz, g_0_z_xxxxxy_xxz, g_0_z_xxxxxy_xxzz, g_0_z_xxxxxy_xyy, g_0_z_xxxxxy_xyyy, g_0_z_xxxxxy_xyyz, g_0_z_xxxxxy_xyz, g_0_z_xxxxxy_xyzz, g_0_z_xxxxxy_xzz, g_0_z_xxxxxy_xzzz, g_0_z_xxxxxy_yyy, g_0_z_xxxxxy_yyz, g_0_z_xxxxxy_yzz, g_0_z_xxxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxy_xxx[k] = -g_0_z_xxxxxy_xxx[k] * ab_x + g_0_z_xxxxxy_xxxx[k];

                g_0_z_xxxxxxy_xxy[k] = -g_0_z_xxxxxy_xxy[k] * ab_x + g_0_z_xxxxxy_xxxy[k];

                g_0_z_xxxxxxy_xxz[k] = -g_0_z_xxxxxy_xxz[k] * ab_x + g_0_z_xxxxxy_xxxz[k];

                g_0_z_xxxxxxy_xyy[k] = -g_0_z_xxxxxy_xyy[k] * ab_x + g_0_z_xxxxxy_xxyy[k];

                g_0_z_xxxxxxy_xyz[k] = -g_0_z_xxxxxy_xyz[k] * ab_x + g_0_z_xxxxxy_xxyz[k];

                g_0_z_xxxxxxy_xzz[k] = -g_0_z_xxxxxy_xzz[k] * ab_x + g_0_z_xxxxxy_xxzz[k];

                g_0_z_xxxxxxy_yyy[k] = -g_0_z_xxxxxy_yyy[k] * ab_x + g_0_z_xxxxxy_xyyy[k];

                g_0_z_xxxxxxy_yyz[k] = -g_0_z_xxxxxy_yyz[k] * ab_x + g_0_z_xxxxxy_xyyz[k];

                g_0_z_xxxxxxy_yzz[k] = -g_0_z_xxxxxy_yzz[k] * ab_x + g_0_z_xxxxxy_xyzz[k];

                g_0_z_xxxxxxy_zzz[k] = -g_0_z_xxxxxy_zzz[k] * ab_x + g_0_z_xxxxxy_xzzz[k];
            }

            /// Set up 740-750 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxz_xxx = cbuffer.data(kf_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxy = cbuffer.data(kf_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxz = cbuffer.data(kf_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xyy = cbuffer.data(kf_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xyz = cbuffer.data(kf_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xzz = cbuffer.data(kf_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yyy = cbuffer.data(kf_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yyz = cbuffer.data(kf_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yzz = cbuffer.data(kf_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_zzz = cbuffer.data(kf_geom_01_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxz_xxx, g_0_z_xxxxxxz_xxy, g_0_z_xxxxxxz_xxz, g_0_z_xxxxxxz_xyy, g_0_z_xxxxxxz_xyz, g_0_z_xxxxxxz_xzz, g_0_z_xxxxxxz_yyy, g_0_z_xxxxxxz_yyz, g_0_z_xxxxxxz_yzz, g_0_z_xxxxxxz_zzz, g_0_z_xxxxxz_xxx, g_0_z_xxxxxz_xxxx, g_0_z_xxxxxz_xxxy, g_0_z_xxxxxz_xxxz, g_0_z_xxxxxz_xxy, g_0_z_xxxxxz_xxyy, g_0_z_xxxxxz_xxyz, g_0_z_xxxxxz_xxz, g_0_z_xxxxxz_xxzz, g_0_z_xxxxxz_xyy, g_0_z_xxxxxz_xyyy, g_0_z_xxxxxz_xyyz, g_0_z_xxxxxz_xyz, g_0_z_xxxxxz_xyzz, g_0_z_xxxxxz_xzz, g_0_z_xxxxxz_xzzz, g_0_z_xxxxxz_yyy, g_0_z_xxxxxz_yyz, g_0_z_xxxxxz_yzz, g_0_z_xxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxz_xxx[k] = -g_0_z_xxxxxz_xxx[k] * ab_x + g_0_z_xxxxxz_xxxx[k];

                g_0_z_xxxxxxz_xxy[k] = -g_0_z_xxxxxz_xxy[k] * ab_x + g_0_z_xxxxxz_xxxy[k];

                g_0_z_xxxxxxz_xxz[k] = -g_0_z_xxxxxz_xxz[k] * ab_x + g_0_z_xxxxxz_xxxz[k];

                g_0_z_xxxxxxz_xyy[k] = -g_0_z_xxxxxz_xyy[k] * ab_x + g_0_z_xxxxxz_xxyy[k];

                g_0_z_xxxxxxz_xyz[k] = -g_0_z_xxxxxz_xyz[k] * ab_x + g_0_z_xxxxxz_xxyz[k];

                g_0_z_xxxxxxz_xzz[k] = -g_0_z_xxxxxz_xzz[k] * ab_x + g_0_z_xxxxxz_xxzz[k];

                g_0_z_xxxxxxz_yyy[k] = -g_0_z_xxxxxz_yyy[k] * ab_x + g_0_z_xxxxxz_xyyy[k];

                g_0_z_xxxxxxz_yyz[k] = -g_0_z_xxxxxz_yyz[k] * ab_x + g_0_z_xxxxxz_xyyz[k];

                g_0_z_xxxxxxz_yzz[k] = -g_0_z_xxxxxz_yzz[k] * ab_x + g_0_z_xxxxxz_xyzz[k];

                g_0_z_xxxxxxz_zzz[k] = -g_0_z_xxxxxz_zzz[k] * ab_x + g_0_z_xxxxxz_xzzz[k];
            }

            /// Set up 750-760 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyy_xxx = cbuffer.data(kf_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxy = cbuffer.data(kf_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxz = cbuffer.data(kf_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xyy = cbuffer.data(kf_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xyz = cbuffer.data(kf_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xzz = cbuffer.data(kf_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yyy = cbuffer.data(kf_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yyz = cbuffer.data(kf_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yzz = cbuffer.data(kf_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_zzz = cbuffer.data(kf_geom_01_off + 759 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyy_xxx, g_0_z_xxxxxyy_xxy, g_0_z_xxxxxyy_xxz, g_0_z_xxxxxyy_xyy, g_0_z_xxxxxyy_xyz, g_0_z_xxxxxyy_xzz, g_0_z_xxxxxyy_yyy, g_0_z_xxxxxyy_yyz, g_0_z_xxxxxyy_yzz, g_0_z_xxxxxyy_zzz, g_0_z_xxxxyy_xxx, g_0_z_xxxxyy_xxxx, g_0_z_xxxxyy_xxxy, g_0_z_xxxxyy_xxxz, g_0_z_xxxxyy_xxy, g_0_z_xxxxyy_xxyy, g_0_z_xxxxyy_xxyz, g_0_z_xxxxyy_xxz, g_0_z_xxxxyy_xxzz, g_0_z_xxxxyy_xyy, g_0_z_xxxxyy_xyyy, g_0_z_xxxxyy_xyyz, g_0_z_xxxxyy_xyz, g_0_z_xxxxyy_xyzz, g_0_z_xxxxyy_xzz, g_0_z_xxxxyy_xzzz, g_0_z_xxxxyy_yyy, g_0_z_xxxxyy_yyz, g_0_z_xxxxyy_yzz, g_0_z_xxxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyy_xxx[k] = -g_0_z_xxxxyy_xxx[k] * ab_x + g_0_z_xxxxyy_xxxx[k];

                g_0_z_xxxxxyy_xxy[k] = -g_0_z_xxxxyy_xxy[k] * ab_x + g_0_z_xxxxyy_xxxy[k];

                g_0_z_xxxxxyy_xxz[k] = -g_0_z_xxxxyy_xxz[k] * ab_x + g_0_z_xxxxyy_xxxz[k];

                g_0_z_xxxxxyy_xyy[k] = -g_0_z_xxxxyy_xyy[k] * ab_x + g_0_z_xxxxyy_xxyy[k];

                g_0_z_xxxxxyy_xyz[k] = -g_0_z_xxxxyy_xyz[k] * ab_x + g_0_z_xxxxyy_xxyz[k];

                g_0_z_xxxxxyy_xzz[k] = -g_0_z_xxxxyy_xzz[k] * ab_x + g_0_z_xxxxyy_xxzz[k];

                g_0_z_xxxxxyy_yyy[k] = -g_0_z_xxxxyy_yyy[k] * ab_x + g_0_z_xxxxyy_xyyy[k];

                g_0_z_xxxxxyy_yyz[k] = -g_0_z_xxxxyy_yyz[k] * ab_x + g_0_z_xxxxyy_xyyz[k];

                g_0_z_xxxxxyy_yzz[k] = -g_0_z_xxxxyy_yzz[k] * ab_x + g_0_z_xxxxyy_xyzz[k];

                g_0_z_xxxxxyy_zzz[k] = -g_0_z_xxxxyy_zzz[k] * ab_x + g_0_z_xxxxyy_xzzz[k];
            }

            /// Set up 760-770 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyz_xxx = cbuffer.data(kf_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxy = cbuffer.data(kf_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxz = cbuffer.data(kf_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xyy = cbuffer.data(kf_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xyz = cbuffer.data(kf_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xzz = cbuffer.data(kf_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yyy = cbuffer.data(kf_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yyz = cbuffer.data(kf_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yzz = cbuffer.data(kf_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_zzz = cbuffer.data(kf_geom_01_off + 769 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyz_xxx, g_0_z_xxxxxyz_xxy, g_0_z_xxxxxyz_xxz, g_0_z_xxxxxyz_xyy, g_0_z_xxxxxyz_xyz, g_0_z_xxxxxyz_xzz, g_0_z_xxxxxyz_yyy, g_0_z_xxxxxyz_yyz, g_0_z_xxxxxyz_yzz, g_0_z_xxxxxyz_zzz, g_0_z_xxxxyz_xxx, g_0_z_xxxxyz_xxxx, g_0_z_xxxxyz_xxxy, g_0_z_xxxxyz_xxxz, g_0_z_xxxxyz_xxy, g_0_z_xxxxyz_xxyy, g_0_z_xxxxyz_xxyz, g_0_z_xxxxyz_xxz, g_0_z_xxxxyz_xxzz, g_0_z_xxxxyz_xyy, g_0_z_xxxxyz_xyyy, g_0_z_xxxxyz_xyyz, g_0_z_xxxxyz_xyz, g_0_z_xxxxyz_xyzz, g_0_z_xxxxyz_xzz, g_0_z_xxxxyz_xzzz, g_0_z_xxxxyz_yyy, g_0_z_xxxxyz_yyz, g_0_z_xxxxyz_yzz, g_0_z_xxxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyz_xxx[k] = -g_0_z_xxxxyz_xxx[k] * ab_x + g_0_z_xxxxyz_xxxx[k];

                g_0_z_xxxxxyz_xxy[k] = -g_0_z_xxxxyz_xxy[k] * ab_x + g_0_z_xxxxyz_xxxy[k];

                g_0_z_xxxxxyz_xxz[k] = -g_0_z_xxxxyz_xxz[k] * ab_x + g_0_z_xxxxyz_xxxz[k];

                g_0_z_xxxxxyz_xyy[k] = -g_0_z_xxxxyz_xyy[k] * ab_x + g_0_z_xxxxyz_xxyy[k];

                g_0_z_xxxxxyz_xyz[k] = -g_0_z_xxxxyz_xyz[k] * ab_x + g_0_z_xxxxyz_xxyz[k];

                g_0_z_xxxxxyz_xzz[k] = -g_0_z_xxxxyz_xzz[k] * ab_x + g_0_z_xxxxyz_xxzz[k];

                g_0_z_xxxxxyz_yyy[k] = -g_0_z_xxxxyz_yyy[k] * ab_x + g_0_z_xxxxyz_xyyy[k];

                g_0_z_xxxxxyz_yyz[k] = -g_0_z_xxxxyz_yyz[k] * ab_x + g_0_z_xxxxyz_xyyz[k];

                g_0_z_xxxxxyz_yzz[k] = -g_0_z_xxxxyz_yzz[k] * ab_x + g_0_z_xxxxyz_xyzz[k];

                g_0_z_xxxxxyz_zzz[k] = -g_0_z_xxxxyz_zzz[k] * ab_x + g_0_z_xxxxyz_xzzz[k];
            }

            /// Set up 770-780 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxzz_xxx = cbuffer.data(kf_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxy = cbuffer.data(kf_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxz = cbuffer.data(kf_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xyy = cbuffer.data(kf_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xyz = cbuffer.data(kf_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xzz = cbuffer.data(kf_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yyy = cbuffer.data(kf_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yyz = cbuffer.data(kf_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yzz = cbuffer.data(kf_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_zzz = cbuffer.data(kf_geom_01_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxzz_xxx, g_0_z_xxxxxzz_xxy, g_0_z_xxxxxzz_xxz, g_0_z_xxxxxzz_xyy, g_0_z_xxxxxzz_xyz, g_0_z_xxxxxzz_xzz, g_0_z_xxxxxzz_yyy, g_0_z_xxxxxzz_yyz, g_0_z_xxxxxzz_yzz, g_0_z_xxxxxzz_zzz, g_0_z_xxxxzz_xxx, g_0_z_xxxxzz_xxxx, g_0_z_xxxxzz_xxxy, g_0_z_xxxxzz_xxxz, g_0_z_xxxxzz_xxy, g_0_z_xxxxzz_xxyy, g_0_z_xxxxzz_xxyz, g_0_z_xxxxzz_xxz, g_0_z_xxxxzz_xxzz, g_0_z_xxxxzz_xyy, g_0_z_xxxxzz_xyyy, g_0_z_xxxxzz_xyyz, g_0_z_xxxxzz_xyz, g_0_z_xxxxzz_xyzz, g_0_z_xxxxzz_xzz, g_0_z_xxxxzz_xzzz, g_0_z_xxxxzz_yyy, g_0_z_xxxxzz_yyz, g_0_z_xxxxzz_yzz, g_0_z_xxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxzz_xxx[k] = -g_0_z_xxxxzz_xxx[k] * ab_x + g_0_z_xxxxzz_xxxx[k];

                g_0_z_xxxxxzz_xxy[k] = -g_0_z_xxxxzz_xxy[k] * ab_x + g_0_z_xxxxzz_xxxy[k];

                g_0_z_xxxxxzz_xxz[k] = -g_0_z_xxxxzz_xxz[k] * ab_x + g_0_z_xxxxzz_xxxz[k];

                g_0_z_xxxxxzz_xyy[k] = -g_0_z_xxxxzz_xyy[k] * ab_x + g_0_z_xxxxzz_xxyy[k];

                g_0_z_xxxxxzz_xyz[k] = -g_0_z_xxxxzz_xyz[k] * ab_x + g_0_z_xxxxzz_xxyz[k];

                g_0_z_xxxxxzz_xzz[k] = -g_0_z_xxxxzz_xzz[k] * ab_x + g_0_z_xxxxzz_xxzz[k];

                g_0_z_xxxxxzz_yyy[k] = -g_0_z_xxxxzz_yyy[k] * ab_x + g_0_z_xxxxzz_xyyy[k];

                g_0_z_xxxxxzz_yyz[k] = -g_0_z_xxxxzz_yyz[k] * ab_x + g_0_z_xxxxzz_xyyz[k];

                g_0_z_xxxxxzz_yzz[k] = -g_0_z_xxxxzz_yzz[k] * ab_x + g_0_z_xxxxzz_xyzz[k];

                g_0_z_xxxxxzz_zzz[k] = -g_0_z_xxxxzz_zzz[k] * ab_x + g_0_z_xxxxzz_xzzz[k];
            }

            /// Set up 780-790 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyy_xxx = cbuffer.data(kf_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxy = cbuffer.data(kf_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxz = cbuffer.data(kf_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xyy = cbuffer.data(kf_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xyz = cbuffer.data(kf_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xzz = cbuffer.data(kf_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yyy = cbuffer.data(kf_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yyz = cbuffer.data(kf_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yzz = cbuffer.data(kf_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_zzz = cbuffer.data(kf_geom_01_off + 789 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyy_xxx, g_0_z_xxxxyyy_xxy, g_0_z_xxxxyyy_xxz, g_0_z_xxxxyyy_xyy, g_0_z_xxxxyyy_xyz, g_0_z_xxxxyyy_xzz, g_0_z_xxxxyyy_yyy, g_0_z_xxxxyyy_yyz, g_0_z_xxxxyyy_yzz, g_0_z_xxxxyyy_zzz, g_0_z_xxxyyy_xxx, g_0_z_xxxyyy_xxxx, g_0_z_xxxyyy_xxxy, g_0_z_xxxyyy_xxxz, g_0_z_xxxyyy_xxy, g_0_z_xxxyyy_xxyy, g_0_z_xxxyyy_xxyz, g_0_z_xxxyyy_xxz, g_0_z_xxxyyy_xxzz, g_0_z_xxxyyy_xyy, g_0_z_xxxyyy_xyyy, g_0_z_xxxyyy_xyyz, g_0_z_xxxyyy_xyz, g_0_z_xxxyyy_xyzz, g_0_z_xxxyyy_xzz, g_0_z_xxxyyy_xzzz, g_0_z_xxxyyy_yyy, g_0_z_xxxyyy_yyz, g_0_z_xxxyyy_yzz, g_0_z_xxxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyy_xxx[k] = -g_0_z_xxxyyy_xxx[k] * ab_x + g_0_z_xxxyyy_xxxx[k];

                g_0_z_xxxxyyy_xxy[k] = -g_0_z_xxxyyy_xxy[k] * ab_x + g_0_z_xxxyyy_xxxy[k];

                g_0_z_xxxxyyy_xxz[k] = -g_0_z_xxxyyy_xxz[k] * ab_x + g_0_z_xxxyyy_xxxz[k];

                g_0_z_xxxxyyy_xyy[k] = -g_0_z_xxxyyy_xyy[k] * ab_x + g_0_z_xxxyyy_xxyy[k];

                g_0_z_xxxxyyy_xyz[k] = -g_0_z_xxxyyy_xyz[k] * ab_x + g_0_z_xxxyyy_xxyz[k];

                g_0_z_xxxxyyy_xzz[k] = -g_0_z_xxxyyy_xzz[k] * ab_x + g_0_z_xxxyyy_xxzz[k];

                g_0_z_xxxxyyy_yyy[k] = -g_0_z_xxxyyy_yyy[k] * ab_x + g_0_z_xxxyyy_xyyy[k];

                g_0_z_xxxxyyy_yyz[k] = -g_0_z_xxxyyy_yyz[k] * ab_x + g_0_z_xxxyyy_xyyz[k];

                g_0_z_xxxxyyy_yzz[k] = -g_0_z_xxxyyy_yzz[k] * ab_x + g_0_z_xxxyyy_xyzz[k];

                g_0_z_xxxxyyy_zzz[k] = -g_0_z_xxxyyy_zzz[k] * ab_x + g_0_z_xxxyyy_xzzz[k];
            }

            /// Set up 790-800 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyz_xxx = cbuffer.data(kf_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxy = cbuffer.data(kf_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxz = cbuffer.data(kf_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xyy = cbuffer.data(kf_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xyz = cbuffer.data(kf_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xzz = cbuffer.data(kf_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yyy = cbuffer.data(kf_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yyz = cbuffer.data(kf_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yzz = cbuffer.data(kf_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_zzz = cbuffer.data(kf_geom_01_off + 799 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyz_xxx, g_0_z_xxxxyyz_xxy, g_0_z_xxxxyyz_xxz, g_0_z_xxxxyyz_xyy, g_0_z_xxxxyyz_xyz, g_0_z_xxxxyyz_xzz, g_0_z_xxxxyyz_yyy, g_0_z_xxxxyyz_yyz, g_0_z_xxxxyyz_yzz, g_0_z_xxxxyyz_zzz, g_0_z_xxxyyz_xxx, g_0_z_xxxyyz_xxxx, g_0_z_xxxyyz_xxxy, g_0_z_xxxyyz_xxxz, g_0_z_xxxyyz_xxy, g_0_z_xxxyyz_xxyy, g_0_z_xxxyyz_xxyz, g_0_z_xxxyyz_xxz, g_0_z_xxxyyz_xxzz, g_0_z_xxxyyz_xyy, g_0_z_xxxyyz_xyyy, g_0_z_xxxyyz_xyyz, g_0_z_xxxyyz_xyz, g_0_z_xxxyyz_xyzz, g_0_z_xxxyyz_xzz, g_0_z_xxxyyz_xzzz, g_0_z_xxxyyz_yyy, g_0_z_xxxyyz_yyz, g_0_z_xxxyyz_yzz, g_0_z_xxxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyz_xxx[k] = -g_0_z_xxxyyz_xxx[k] * ab_x + g_0_z_xxxyyz_xxxx[k];

                g_0_z_xxxxyyz_xxy[k] = -g_0_z_xxxyyz_xxy[k] * ab_x + g_0_z_xxxyyz_xxxy[k];

                g_0_z_xxxxyyz_xxz[k] = -g_0_z_xxxyyz_xxz[k] * ab_x + g_0_z_xxxyyz_xxxz[k];

                g_0_z_xxxxyyz_xyy[k] = -g_0_z_xxxyyz_xyy[k] * ab_x + g_0_z_xxxyyz_xxyy[k];

                g_0_z_xxxxyyz_xyz[k] = -g_0_z_xxxyyz_xyz[k] * ab_x + g_0_z_xxxyyz_xxyz[k];

                g_0_z_xxxxyyz_xzz[k] = -g_0_z_xxxyyz_xzz[k] * ab_x + g_0_z_xxxyyz_xxzz[k];

                g_0_z_xxxxyyz_yyy[k] = -g_0_z_xxxyyz_yyy[k] * ab_x + g_0_z_xxxyyz_xyyy[k];

                g_0_z_xxxxyyz_yyz[k] = -g_0_z_xxxyyz_yyz[k] * ab_x + g_0_z_xxxyyz_xyyz[k];

                g_0_z_xxxxyyz_yzz[k] = -g_0_z_xxxyyz_yzz[k] * ab_x + g_0_z_xxxyyz_xyzz[k];

                g_0_z_xxxxyyz_zzz[k] = -g_0_z_xxxyyz_zzz[k] * ab_x + g_0_z_xxxyyz_xzzz[k];
            }

            /// Set up 800-810 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyzz_xxx = cbuffer.data(kf_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxy = cbuffer.data(kf_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxz = cbuffer.data(kf_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xyy = cbuffer.data(kf_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xyz = cbuffer.data(kf_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xzz = cbuffer.data(kf_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yyy = cbuffer.data(kf_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yyz = cbuffer.data(kf_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yzz = cbuffer.data(kf_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_zzz = cbuffer.data(kf_geom_01_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyzz_xxx, g_0_z_xxxxyzz_xxy, g_0_z_xxxxyzz_xxz, g_0_z_xxxxyzz_xyy, g_0_z_xxxxyzz_xyz, g_0_z_xxxxyzz_xzz, g_0_z_xxxxyzz_yyy, g_0_z_xxxxyzz_yyz, g_0_z_xxxxyzz_yzz, g_0_z_xxxxyzz_zzz, g_0_z_xxxyzz_xxx, g_0_z_xxxyzz_xxxx, g_0_z_xxxyzz_xxxy, g_0_z_xxxyzz_xxxz, g_0_z_xxxyzz_xxy, g_0_z_xxxyzz_xxyy, g_0_z_xxxyzz_xxyz, g_0_z_xxxyzz_xxz, g_0_z_xxxyzz_xxzz, g_0_z_xxxyzz_xyy, g_0_z_xxxyzz_xyyy, g_0_z_xxxyzz_xyyz, g_0_z_xxxyzz_xyz, g_0_z_xxxyzz_xyzz, g_0_z_xxxyzz_xzz, g_0_z_xxxyzz_xzzz, g_0_z_xxxyzz_yyy, g_0_z_xxxyzz_yyz, g_0_z_xxxyzz_yzz, g_0_z_xxxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyzz_xxx[k] = -g_0_z_xxxyzz_xxx[k] * ab_x + g_0_z_xxxyzz_xxxx[k];

                g_0_z_xxxxyzz_xxy[k] = -g_0_z_xxxyzz_xxy[k] * ab_x + g_0_z_xxxyzz_xxxy[k];

                g_0_z_xxxxyzz_xxz[k] = -g_0_z_xxxyzz_xxz[k] * ab_x + g_0_z_xxxyzz_xxxz[k];

                g_0_z_xxxxyzz_xyy[k] = -g_0_z_xxxyzz_xyy[k] * ab_x + g_0_z_xxxyzz_xxyy[k];

                g_0_z_xxxxyzz_xyz[k] = -g_0_z_xxxyzz_xyz[k] * ab_x + g_0_z_xxxyzz_xxyz[k];

                g_0_z_xxxxyzz_xzz[k] = -g_0_z_xxxyzz_xzz[k] * ab_x + g_0_z_xxxyzz_xxzz[k];

                g_0_z_xxxxyzz_yyy[k] = -g_0_z_xxxyzz_yyy[k] * ab_x + g_0_z_xxxyzz_xyyy[k];

                g_0_z_xxxxyzz_yyz[k] = -g_0_z_xxxyzz_yyz[k] * ab_x + g_0_z_xxxyzz_xyyz[k];

                g_0_z_xxxxyzz_yzz[k] = -g_0_z_xxxyzz_yzz[k] * ab_x + g_0_z_xxxyzz_xyzz[k];

                g_0_z_xxxxyzz_zzz[k] = -g_0_z_xxxyzz_zzz[k] * ab_x + g_0_z_xxxyzz_xzzz[k];
            }

            /// Set up 810-820 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzzz_xxx = cbuffer.data(kf_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxy = cbuffer.data(kf_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxz = cbuffer.data(kf_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xyy = cbuffer.data(kf_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xyz = cbuffer.data(kf_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xzz = cbuffer.data(kf_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yyy = cbuffer.data(kf_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yyz = cbuffer.data(kf_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yzz = cbuffer.data(kf_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_zzz = cbuffer.data(kf_geom_01_off + 819 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzzz_xxx, g_0_z_xxxxzzz_xxy, g_0_z_xxxxzzz_xxz, g_0_z_xxxxzzz_xyy, g_0_z_xxxxzzz_xyz, g_0_z_xxxxzzz_xzz, g_0_z_xxxxzzz_yyy, g_0_z_xxxxzzz_yyz, g_0_z_xxxxzzz_yzz, g_0_z_xxxxzzz_zzz, g_0_z_xxxzzz_xxx, g_0_z_xxxzzz_xxxx, g_0_z_xxxzzz_xxxy, g_0_z_xxxzzz_xxxz, g_0_z_xxxzzz_xxy, g_0_z_xxxzzz_xxyy, g_0_z_xxxzzz_xxyz, g_0_z_xxxzzz_xxz, g_0_z_xxxzzz_xxzz, g_0_z_xxxzzz_xyy, g_0_z_xxxzzz_xyyy, g_0_z_xxxzzz_xyyz, g_0_z_xxxzzz_xyz, g_0_z_xxxzzz_xyzz, g_0_z_xxxzzz_xzz, g_0_z_xxxzzz_xzzz, g_0_z_xxxzzz_yyy, g_0_z_xxxzzz_yyz, g_0_z_xxxzzz_yzz, g_0_z_xxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzzz_xxx[k] = -g_0_z_xxxzzz_xxx[k] * ab_x + g_0_z_xxxzzz_xxxx[k];

                g_0_z_xxxxzzz_xxy[k] = -g_0_z_xxxzzz_xxy[k] * ab_x + g_0_z_xxxzzz_xxxy[k];

                g_0_z_xxxxzzz_xxz[k] = -g_0_z_xxxzzz_xxz[k] * ab_x + g_0_z_xxxzzz_xxxz[k];

                g_0_z_xxxxzzz_xyy[k] = -g_0_z_xxxzzz_xyy[k] * ab_x + g_0_z_xxxzzz_xxyy[k];

                g_0_z_xxxxzzz_xyz[k] = -g_0_z_xxxzzz_xyz[k] * ab_x + g_0_z_xxxzzz_xxyz[k];

                g_0_z_xxxxzzz_xzz[k] = -g_0_z_xxxzzz_xzz[k] * ab_x + g_0_z_xxxzzz_xxzz[k];

                g_0_z_xxxxzzz_yyy[k] = -g_0_z_xxxzzz_yyy[k] * ab_x + g_0_z_xxxzzz_xyyy[k];

                g_0_z_xxxxzzz_yyz[k] = -g_0_z_xxxzzz_yyz[k] * ab_x + g_0_z_xxxzzz_xyyz[k];

                g_0_z_xxxxzzz_yzz[k] = -g_0_z_xxxzzz_yzz[k] * ab_x + g_0_z_xxxzzz_xyzz[k];

                g_0_z_xxxxzzz_zzz[k] = -g_0_z_xxxzzz_zzz[k] * ab_x + g_0_z_xxxzzz_xzzz[k];
            }

            /// Set up 820-830 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyy_xxx = cbuffer.data(kf_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxy = cbuffer.data(kf_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxz = cbuffer.data(kf_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xyy = cbuffer.data(kf_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xyz = cbuffer.data(kf_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xzz = cbuffer.data(kf_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yyy = cbuffer.data(kf_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yyz = cbuffer.data(kf_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yzz = cbuffer.data(kf_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_zzz = cbuffer.data(kf_geom_01_off + 829 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyy_xxx, g_0_z_xxxyyyy_xxy, g_0_z_xxxyyyy_xxz, g_0_z_xxxyyyy_xyy, g_0_z_xxxyyyy_xyz, g_0_z_xxxyyyy_xzz, g_0_z_xxxyyyy_yyy, g_0_z_xxxyyyy_yyz, g_0_z_xxxyyyy_yzz, g_0_z_xxxyyyy_zzz, g_0_z_xxyyyy_xxx, g_0_z_xxyyyy_xxxx, g_0_z_xxyyyy_xxxy, g_0_z_xxyyyy_xxxz, g_0_z_xxyyyy_xxy, g_0_z_xxyyyy_xxyy, g_0_z_xxyyyy_xxyz, g_0_z_xxyyyy_xxz, g_0_z_xxyyyy_xxzz, g_0_z_xxyyyy_xyy, g_0_z_xxyyyy_xyyy, g_0_z_xxyyyy_xyyz, g_0_z_xxyyyy_xyz, g_0_z_xxyyyy_xyzz, g_0_z_xxyyyy_xzz, g_0_z_xxyyyy_xzzz, g_0_z_xxyyyy_yyy, g_0_z_xxyyyy_yyz, g_0_z_xxyyyy_yzz, g_0_z_xxyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyy_xxx[k] = -g_0_z_xxyyyy_xxx[k] * ab_x + g_0_z_xxyyyy_xxxx[k];

                g_0_z_xxxyyyy_xxy[k] = -g_0_z_xxyyyy_xxy[k] * ab_x + g_0_z_xxyyyy_xxxy[k];

                g_0_z_xxxyyyy_xxz[k] = -g_0_z_xxyyyy_xxz[k] * ab_x + g_0_z_xxyyyy_xxxz[k];

                g_0_z_xxxyyyy_xyy[k] = -g_0_z_xxyyyy_xyy[k] * ab_x + g_0_z_xxyyyy_xxyy[k];

                g_0_z_xxxyyyy_xyz[k] = -g_0_z_xxyyyy_xyz[k] * ab_x + g_0_z_xxyyyy_xxyz[k];

                g_0_z_xxxyyyy_xzz[k] = -g_0_z_xxyyyy_xzz[k] * ab_x + g_0_z_xxyyyy_xxzz[k];

                g_0_z_xxxyyyy_yyy[k] = -g_0_z_xxyyyy_yyy[k] * ab_x + g_0_z_xxyyyy_xyyy[k];

                g_0_z_xxxyyyy_yyz[k] = -g_0_z_xxyyyy_yyz[k] * ab_x + g_0_z_xxyyyy_xyyz[k];

                g_0_z_xxxyyyy_yzz[k] = -g_0_z_xxyyyy_yzz[k] * ab_x + g_0_z_xxyyyy_xyzz[k];

                g_0_z_xxxyyyy_zzz[k] = -g_0_z_xxyyyy_zzz[k] * ab_x + g_0_z_xxyyyy_xzzz[k];
            }

            /// Set up 830-840 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyz_xxx = cbuffer.data(kf_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxy = cbuffer.data(kf_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxz = cbuffer.data(kf_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xyy = cbuffer.data(kf_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xyz = cbuffer.data(kf_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xzz = cbuffer.data(kf_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yyy = cbuffer.data(kf_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yyz = cbuffer.data(kf_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yzz = cbuffer.data(kf_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_zzz = cbuffer.data(kf_geom_01_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyz_xxx, g_0_z_xxxyyyz_xxy, g_0_z_xxxyyyz_xxz, g_0_z_xxxyyyz_xyy, g_0_z_xxxyyyz_xyz, g_0_z_xxxyyyz_xzz, g_0_z_xxxyyyz_yyy, g_0_z_xxxyyyz_yyz, g_0_z_xxxyyyz_yzz, g_0_z_xxxyyyz_zzz, g_0_z_xxyyyz_xxx, g_0_z_xxyyyz_xxxx, g_0_z_xxyyyz_xxxy, g_0_z_xxyyyz_xxxz, g_0_z_xxyyyz_xxy, g_0_z_xxyyyz_xxyy, g_0_z_xxyyyz_xxyz, g_0_z_xxyyyz_xxz, g_0_z_xxyyyz_xxzz, g_0_z_xxyyyz_xyy, g_0_z_xxyyyz_xyyy, g_0_z_xxyyyz_xyyz, g_0_z_xxyyyz_xyz, g_0_z_xxyyyz_xyzz, g_0_z_xxyyyz_xzz, g_0_z_xxyyyz_xzzz, g_0_z_xxyyyz_yyy, g_0_z_xxyyyz_yyz, g_0_z_xxyyyz_yzz, g_0_z_xxyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyz_xxx[k] = -g_0_z_xxyyyz_xxx[k] * ab_x + g_0_z_xxyyyz_xxxx[k];

                g_0_z_xxxyyyz_xxy[k] = -g_0_z_xxyyyz_xxy[k] * ab_x + g_0_z_xxyyyz_xxxy[k];

                g_0_z_xxxyyyz_xxz[k] = -g_0_z_xxyyyz_xxz[k] * ab_x + g_0_z_xxyyyz_xxxz[k];

                g_0_z_xxxyyyz_xyy[k] = -g_0_z_xxyyyz_xyy[k] * ab_x + g_0_z_xxyyyz_xxyy[k];

                g_0_z_xxxyyyz_xyz[k] = -g_0_z_xxyyyz_xyz[k] * ab_x + g_0_z_xxyyyz_xxyz[k];

                g_0_z_xxxyyyz_xzz[k] = -g_0_z_xxyyyz_xzz[k] * ab_x + g_0_z_xxyyyz_xxzz[k];

                g_0_z_xxxyyyz_yyy[k] = -g_0_z_xxyyyz_yyy[k] * ab_x + g_0_z_xxyyyz_xyyy[k];

                g_0_z_xxxyyyz_yyz[k] = -g_0_z_xxyyyz_yyz[k] * ab_x + g_0_z_xxyyyz_xyyz[k];

                g_0_z_xxxyyyz_yzz[k] = -g_0_z_xxyyyz_yzz[k] * ab_x + g_0_z_xxyyyz_xyzz[k];

                g_0_z_xxxyyyz_zzz[k] = -g_0_z_xxyyyz_zzz[k] * ab_x + g_0_z_xxyyyz_xzzz[k];
            }

            /// Set up 840-850 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyzz_xxx = cbuffer.data(kf_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxy = cbuffer.data(kf_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxz = cbuffer.data(kf_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xyy = cbuffer.data(kf_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xyz = cbuffer.data(kf_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xzz = cbuffer.data(kf_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yyy = cbuffer.data(kf_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yyz = cbuffer.data(kf_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yzz = cbuffer.data(kf_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_zzz = cbuffer.data(kf_geom_01_off + 849 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyzz_xxx, g_0_z_xxxyyzz_xxy, g_0_z_xxxyyzz_xxz, g_0_z_xxxyyzz_xyy, g_0_z_xxxyyzz_xyz, g_0_z_xxxyyzz_xzz, g_0_z_xxxyyzz_yyy, g_0_z_xxxyyzz_yyz, g_0_z_xxxyyzz_yzz, g_0_z_xxxyyzz_zzz, g_0_z_xxyyzz_xxx, g_0_z_xxyyzz_xxxx, g_0_z_xxyyzz_xxxy, g_0_z_xxyyzz_xxxz, g_0_z_xxyyzz_xxy, g_0_z_xxyyzz_xxyy, g_0_z_xxyyzz_xxyz, g_0_z_xxyyzz_xxz, g_0_z_xxyyzz_xxzz, g_0_z_xxyyzz_xyy, g_0_z_xxyyzz_xyyy, g_0_z_xxyyzz_xyyz, g_0_z_xxyyzz_xyz, g_0_z_xxyyzz_xyzz, g_0_z_xxyyzz_xzz, g_0_z_xxyyzz_xzzz, g_0_z_xxyyzz_yyy, g_0_z_xxyyzz_yyz, g_0_z_xxyyzz_yzz, g_0_z_xxyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyzz_xxx[k] = -g_0_z_xxyyzz_xxx[k] * ab_x + g_0_z_xxyyzz_xxxx[k];

                g_0_z_xxxyyzz_xxy[k] = -g_0_z_xxyyzz_xxy[k] * ab_x + g_0_z_xxyyzz_xxxy[k];

                g_0_z_xxxyyzz_xxz[k] = -g_0_z_xxyyzz_xxz[k] * ab_x + g_0_z_xxyyzz_xxxz[k];

                g_0_z_xxxyyzz_xyy[k] = -g_0_z_xxyyzz_xyy[k] * ab_x + g_0_z_xxyyzz_xxyy[k];

                g_0_z_xxxyyzz_xyz[k] = -g_0_z_xxyyzz_xyz[k] * ab_x + g_0_z_xxyyzz_xxyz[k];

                g_0_z_xxxyyzz_xzz[k] = -g_0_z_xxyyzz_xzz[k] * ab_x + g_0_z_xxyyzz_xxzz[k];

                g_0_z_xxxyyzz_yyy[k] = -g_0_z_xxyyzz_yyy[k] * ab_x + g_0_z_xxyyzz_xyyy[k];

                g_0_z_xxxyyzz_yyz[k] = -g_0_z_xxyyzz_yyz[k] * ab_x + g_0_z_xxyyzz_xyyz[k];

                g_0_z_xxxyyzz_yzz[k] = -g_0_z_xxyyzz_yzz[k] * ab_x + g_0_z_xxyyzz_xyzz[k];

                g_0_z_xxxyyzz_zzz[k] = -g_0_z_xxyyzz_zzz[k] * ab_x + g_0_z_xxyyzz_xzzz[k];
            }

            /// Set up 850-860 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzzz_xxx = cbuffer.data(kf_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxy = cbuffer.data(kf_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxz = cbuffer.data(kf_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xyy = cbuffer.data(kf_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xyz = cbuffer.data(kf_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xzz = cbuffer.data(kf_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yyy = cbuffer.data(kf_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yyz = cbuffer.data(kf_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yzz = cbuffer.data(kf_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_zzz = cbuffer.data(kf_geom_01_off + 859 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzzz_xxx, g_0_z_xxxyzzz_xxy, g_0_z_xxxyzzz_xxz, g_0_z_xxxyzzz_xyy, g_0_z_xxxyzzz_xyz, g_0_z_xxxyzzz_xzz, g_0_z_xxxyzzz_yyy, g_0_z_xxxyzzz_yyz, g_0_z_xxxyzzz_yzz, g_0_z_xxxyzzz_zzz, g_0_z_xxyzzz_xxx, g_0_z_xxyzzz_xxxx, g_0_z_xxyzzz_xxxy, g_0_z_xxyzzz_xxxz, g_0_z_xxyzzz_xxy, g_0_z_xxyzzz_xxyy, g_0_z_xxyzzz_xxyz, g_0_z_xxyzzz_xxz, g_0_z_xxyzzz_xxzz, g_0_z_xxyzzz_xyy, g_0_z_xxyzzz_xyyy, g_0_z_xxyzzz_xyyz, g_0_z_xxyzzz_xyz, g_0_z_xxyzzz_xyzz, g_0_z_xxyzzz_xzz, g_0_z_xxyzzz_xzzz, g_0_z_xxyzzz_yyy, g_0_z_xxyzzz_yyz, g_0_z_xxyzzz_yzz, g_0_z_xxyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzzz_xxx[k] = -g_0_z_xxyzzz_xxx[k] * ab_x + g_0_z_xxyzzz_xxxx[k];

                g_0_z_xxxyzzz_xxy[k] = -g_0_z_xxyzzz_xxy[k] * ab_x + g_0_z_xxyzzz_xxxy[k];

                g_0_z_xxxyzzz_xxz[k] = -g_0_z_xxyzzz_xxz[k] * ab_x + g_0_z_xxyzzz_xxxz[k];

                g_0_z_xxxyzzz_xyy[k] = -g_0_z_xxyzzz_xyy[k] * ab_x + g_0_z_xxyzzz_xxyy[k];

                g_0_z_xxxyzzz_xyz[k] = -g_0_z_xxyzzz_xyz[k] * ab_x + g_0_z_xxyzzz_xxyz[k];

                g_0_z_xxxyzzz_xzz[k] = -g_0_z_xxyzzz_xzz[k] * ab_x + g_0_z_xxyzzz_xxzz[k];

                g_0_z_xxxyzzz_yyy[k] = -g_0_z_xxyzzz_yyy[k] * ab_x + g_0_z_xxyzzz_xyyy[k];

                g_0_z_xxxyzzz_yyz[k] = -g_0_z_xxyzzz_yyz[k] * ab_x + g_0_z_xxyzzz_xyyz[k];

                g_0_z_xxxyzzz_yzz[k] = -g_0_z_xxyzzz_yzz[k] * ab_x + g_0_z_xxyzzz_xyzz[k];

                g_0_z_xxxyzzz_zzz[k] = -g_0_z_xxyzzz_zzz[k] * ab_x + g_0_z_xxyzzz_xzzz[k];
            }

            /// Set up 860-870 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzzz_xxx = cbuffer.data(kf_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxy = cbuffer.data(kf_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxz = cbuffer.data(kf_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xyy = cbuffer.data(kf_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xyz = cbuffer.data(kf_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xzz = cbuffer.data(kf_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yyy = cbuffer.data(kf_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yyz = cbuffer.data(kf_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yzz = cbuffer.data(kf_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_zzz = cbuffer.data(kf_geom_01_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzzz_xxx, g_0_z_xxxzzzz_xxy, g_0_z_xxxzzzz_xxz, g_0_z_xxxzzzz_xyy, g_0_z_xxxzzzz_xyz, g_0_z_xxxzzzz_xzz, g_0_z_xxxzzzz_yyy, g_0_z_xxxzzzz_yyz, g_0_z_xxxzzzz_yzz, g_0_z_xxxzzzz_zzz, g_0_z_xxzzzz_xxx, g_0_z_xxzzzz_xxxx, g_0_z_xxzzzz_xxxy, g_0_z_xxzzzz_xxxz, g_0_z_xxzzzz_xxy, g_0_z_xxzzzz_xxyy, g_0_z_xxzzzz_xxyz, g_0_z_xxzzzz_xxz, g_0_z_xxzzzz_xxzz, g_0_z_xxzzzz_xyy, g_0_z_xxzzzz_xyyy, g_0_z_xxzzzz_xyyz, g_0_z_xxzzzz_xyz, g_0_z_xxzzzz_xyzz, g_0_z_xxzzzz_xzz, g_0_z_xxzzzz_xzzz, g_0_z_xxzzzz_yyy, g_0_z_xxzzzz_yyz, g_0_z_xxzzzz_yzz, g_0_z_xxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzzz_xxx[k] = -g_0_z_xxzzzz_xxx[k] * ab_x + g_0_z_xxzzzz_xxxx[k];

                g_0_z_xxxzzzz_xxy[k] = -g_0_z_xxzzzz_xxy[k] * ab_x + g_0_z_xxzzzz_xxxy[k];

                g_0_z_xxxzzzz_xxz[k] = -g_0_z_xxzzzz_xxz[k] * ab_x + g_0_z_xxzzzz_xxxz[k];

                g_0_z_xxxzzzz_xyy[k] = -g_0_z_xxzzzz_xyy[k] * ab_x + g_0_z_xxzzzz_xxyy[k];

                g_0_z_xxxzzzz_xyz[k] = -g_0_z_xxzzzz_xyz[k] * ab_x + g_0_z_xxzzzz_xxyz[k];

                g_0_z_xxxzzzz_xzz[k] = -g_0_z_xxzzzz_xzz[k] * ab_x + g_0_z_xxzzzz_xxzz[k];

                g_0_z_xxxzzzz_yyy[k] = -g_0_z_xxzzzz_yyy[k] * ab_x + g_0_z_xxzzzz_xyyy[k];

                g_0_z_xxxzzzz_yyz[k] = -g_0_z_xxzzzz_yyz[k] * ab_x + g_0_z_xxzzzz_xyyz[k];

                g_0_z_xxxzzzz_yzz[k] = -g_0_z_xxzzzz_yzz[k] * ab_x + g_0_z_xxzzzz_xyzz[k];

                g_0_z_xxxzzzz_zzz[k] = -g_0_z_xxzzzz_zzz[k] * ab_x + g_0_z_xxzzzz_xzzz[k];
            }

            /// Set up 870-880 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyy_xxx = cbuffer.data(kf_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxy = cbuffer.data(kf_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxz = cbuffer.data(kf_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xyy = cbuffer.data(kf_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xyz = cbuffer.data(kf_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xzz = cbuffer.data(kf_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yyy = cbuffer.data(kf_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yyz = cbuffer.data(kf_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yzz = cbuffer.data(kf_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_zzz = cbuffer.data(kf_geom_01_off + 879 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyy_xxx, g_0_z_xxyyyyy_xxy, g_0_z_xxyyyyy_xxz, g_0_z_xxyyyyy_xyy, g_0_z_xxyyyyy_xyz, g_0_z_xxyyyyy_xzz, g_0_z_xxyyyyy_yyy, g_0_z_xxyyyyy_yyz, g_0_z_xxyyyyy_yzz, g_0_z_xxyyyyy_zzz, g_0_z_xyyyyy_xxx, g_0_z_xyyyyy_xxxx, g_0_z_xyyyyy_xxxy, g_0_z_xyyyyy_xxxz, g_0_z_xyyyyy_xxy, g_0_z_xyyyyy_xxyy, g_0_z_xyyyyy_xxyz, g_0_z_xyyyyy_xxz, g_0_z_xyyyyy_xxzz, g_0_z_xyyyyy_xyy, g_0_z_xyyyyy_xyyy, g_0_z_xyyyyy_xyyz, g_0_z_xyyyyy_xyz, g_0_z_xyyyyy_xyzz, g_0_z_xyyyyy_xzz, g_0_z_xyyyyy_xzzz, g_0_z_xyyyyy_yyy, g_0_z_xyyyyy_yyz, g_0_z_xyyyyy_yzz, g_0_z_xyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyy_xxx[k] = -g_0_z_xyyyyy_xxx[k] * ab_x + g_0_z_xyyyyy_xxxx[k];

                g_0_z_xxyyyyy_xxy[k] = -g_0_z_xyyyyy_xxy[k] * ab_x + g_0_z_xyyyyy_xxxy[k];

                g_0_z_xxyyyyy_xxz[k] = -g_0_z_xyyyyy_xxz[k] * ab_x + g_0_z_xyyyyy_xxxz[k];

                g_0_z_xxyyyyy_xyy[k] = -g_0_z_xyyyyy_xyy[k] * ab_x + g_0_z_xyyyyy_xxyy[k];

                g_0_z_xxyyyyy_xyz[k] = -g_0_z_xyyyyy_xyz[k] * ab_x + g_0_z_xyyyyy_xxyz[k];

                g_0_z_xxyyyyy_xzz[k] = -g_0_z_xyyyyy_xzz[k] * ab_x + g_0_z_xyyyyy_xxzz[k];

                g_0_z_xxyyyyy_yyy[k] = -g_0_z_xyyyyy_yyy[k] * ab_x + g_0_z_xyyyyy_xyyy[k];

                g_0_z_xxyyyyy_yyz[k] = -g_0_z_xyyyyy_yyz[k] * ab_x + g_0_z_xyyyyy_xyyz[k];

                g_0_z_xxyyyyy_yzz[k] = -g_0_z_xyyyyy_yzz[k] * ab_x + g_0_z_xyyyyy_xyzz[k];

                g_0_z_xxyyyyy_zzz[k] = -g_0_z_xyyyyy_zzz[k] * ab_x + g_0_z_xyyyyy_xzzz[k];
            }

            /// Set up 880-890 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyz_xxx = cbuffer.data(kf_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxy = cbuffer.data(kf_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxz = cbuffer.data(kf_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xyy = cbuffer.data(kf_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xyz = cbuffer.data(kf_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xzz = cbuffer.data(kf_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yyy = cbuffer.data(kf_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yyz = cbuffer.data(kf_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yzz = cbuffer.data(kf_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_zzz = cbuffer.data(kf_geom_01_off + 889 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyz_xxx, g_0_z_xxyyyyz_xxy, g_0_z_xxyyyyz_xxz, g_0_z_xxyyyyz_xyy, g_0_z_xxyyyyz_xyz, g_0_z_xxyyyyz_xzz, g_0_z_xxyyyyz_yyy, g_0_z_xxyyyyz_yyz, g_0_z_xxyyyyz_yzz, g_0_z_xxyyyyz_zzz, g_0_z_xyyyyz_xxx, g_0_z_xyyyyz_xxxx, g_0_z_xyyyyz_xxxy, g_0_z_xyyyyz_xxxz, g_0_z_xyyyyz_xxy, g_0_z_xyyyyz_xxyy, g_0_z_xyyyyz_xxyz, g_0_z_xyyyyz_xxz, g_0_z_xyyyyz_xxzz, g_0_z_xyyyyz_xyy, g_0_z_xyyyyz_xyyy, g_0_z_xyyyyz_xyyz, g_0_z_xyyyyz_xyz, g_0_z_xyyyyz_xyzz, g_0_z_xyyyyz_xzz, g_0_z_xyyyyz_xzzz, g_0_z_xyyyyz_yyy, g_0_z_xyyyyz_yyz, g_0_z_xyyyyz_yzz, g_0_z_xyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyz_xxx[k] = -g_0_z_xyyyyz_xxx[k] * ab_x + g_0_z_xyyyyz_xxxx[k];

                g_0_z_xxyyyyz_xxy[k] = -g_0_z_xyyyyz_xxy[k] * ab_x + g_0_z_xyyyyz_xxxy[k];

                g_0_z_xxyyyyz_xxz[k] = -g_0_z_xyyyyz_xxz[k] * ab_x + g_0_z_xyyyyz_xxxz[k];

                g_0_z_xxyyyyz_xyy[k] = -g_0_z_xyyyyz_xyy[k] * ab_x + g_0_z_xyyyyz_xxyy[k];

                g_0_z_xxyyyyz_xyz[k] = -g_0_z_xyyyyz_xyz[k] * ab_x + g_0_z_xyyyyz_xxyz[k];

                g_0_z_xxyyyyz_xzz[k] = -g_0_z_xyyyyz_xzz[k] * ab_x + g_0_z_xyyyyz_xxzz[k];

                g_0_z_xxyyyyz_yyy[k] = -g_0_z_xyyyyz_yyy[k] * ab_x + g_0_z_xyyyyz_xyyy[k];

                g_0_z_xxyyyyz_yyz[k] = -g_0_z_xyyyyz_yyz[k] * ab_x + g_0_z_xyyyyz_xyyz[k];

                g_0_z_xxyyyyz_yzz[k] = -g_0_z_xyyyyz_yzz[k] * ab_x + g_0_z_xyyyyz_xyzz[k];

                g_0_z_xxyyyyz_zzz[k] = -g_0_z_xyyyyz_zzz[k] * ab_x + g_0_z_xyyyyz_xzzz[k];
            }

            /// Set up 890-900 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyzz_xxx = cbuffer.data(kf_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxy = cbuffer.data(kf_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxz = cbuffer.data(kf_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xyy = cbuffer.data(kf_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xyz = cbuffer.data(kf_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xzz = cbuffer.data(kf_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yyy = cbuffer.data(kf_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yyz = cbuffer.data(kf_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yzz = cbuffer.data(kf_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_zzz = cbuffer.data(kf_geom_01_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyzz_xxx, g_0_z_xxyyyzz_xxy, g_0_z_xxyyyzz_xxz, g_0_z_xxyyyzz_xyy, g_0_z_xxyyyzz_xyz, g_0_z_xxyyyzz_xzz, g_0_z_xxyyyzz_yyy, g_0_z_xxyyyzz_yyz, g_0_z_xxyyyzz_yzz, g_0_z_xxyyyzz_zzz, g_0_z_xyyyzz_xxx, g_0_z_xyyyzz_xxxx, g_0_z_xyyyzz_xxxy, g_0_z_xyyyzz_xxxz, g_0_z_xyyyzz_xxy, g_0_z_xyyyzz_xxyy, g_0_z_xyyyzz_xxyz, g_0_z_xyyyzz_xxz, g_0_z_xyyyzz_xxzz, g_0_z_xyyyzz_xyy, g_0_z_xyyyzz_xyyy, g_0_z_xyyyzz_xyyz, g_0_z_xyyyzz_xyz, g_0_z_xyyyzz_xyzz, g_0_z_xyyyzz_xzz, g_0_z_xyyyzz_xzzz, g_0_z_xyyyzz_yyy, g_0_z_xyyyzz_yyz, g_0_z_xyyyzz_yzz, g_0_z_xyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyzz_xxx[k] = -g_0_z_xyyyzz_xxx[k] * ab_x + g_0_z_xyyyzz_xxxx[k];

                g_0_z_xxyyyzz_xxy[k] = -g_0_z_xyyyzz_xxy[k] * ab_x + g_0_z_xyyyzz_xxxy[k];

                g_0_z_xxyyyzz_xxz[k] = -g_0_z_xyyyzz_xxz[k] * ab_x + g_0_z_xyyyzz_xxxz[k];

                g_0_z_xxyyyzz_xyy[k] = -g_0_z_xyyyzz_xyy[k] * ab_x + g_0_z_xyyyzz_xxyy[k];

                g_0_z_xxyyyzz_xyz[k] = -g_0_z_xyyyzz_xyz[k] * ab_x + g_0_z_xyyyzz_xxyz[k];

                g_0_z_xxyyyzz_xzz[k] = -g_0_z_xyyyzz_xzz[k] * ab_x + g_0_z_xyyyzz_xxzz[k];

                g_0_z_xxyyyzz_yyy[k] = -g_0_z_xyyyzz_yyy[k] * ab_x + g_0_z_xyyyzz_xyyy[k];

                g_0_z_xxyyyzz_yyz[k] = -g_0_z_xyyyzz_yyz[k] * ab_x + g_0_z_xyyyzz_xyyz[k];

                g_0_z_xxyyyzz_yzz[k] = -g_0_z_xyyyzz_yzz[k] * ab_x + g_0_z_xyyyzz_xyzz[k];

                g_0_z_xxyyyzz_zzz[k] = -g_0_z_xyyyzz_zzz[k] * ab_x + g_0_z_xyyyzz_xzzz[k];
            }

            /// Set up 900-910 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzzz_xxx = cbuffer.data(kf_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxy = cbuffer.data(kf_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxz = cbuffer.data(kf_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xyy = cbuffer.data(kf_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xyz = cbuffer.data(kf_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xzz = cbuffer.data(kf_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yyy = cbuffer.data(kf_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yyz = cbuffer.data(kf_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yzz = cbuffer.data(kf_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_zzz = cbuffer.data(kf_geom_01_off + 909 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzzz_xxx, g_0_z_xxyyzzz_xxy, g_0_z_xxyyzzz_xxz, g_0_z_xxyyzzz_xyy, g_0_z_xxyyzzz_xyz, g_0_z_xxyyzzz_xzz, g_0_z_xxyyzzz_yyy, g_0_z_xxyyzzz_yyz, g_0_z_xxyyzzz_yzz, g_0_z_xxyyzzz_zzz, g_0_z_xyyzzz_xxx, g_0_z_xyyzzz_xxxx, g_0_z_xyyzzz_xxxy, g_0_z_xyyzzz_xxxz, g_0_z_xyyzzz_xxy, g_0_z_xyyzzz_xxyy, g_0_z_xyyzzz_xxyz, g_0_z_xyyzzz_xxz, g_0_z_xyyzzz_xxzz, g_0_z_xyyzzz_xyy, g_0_z_xyyzzz_xyyy, g_0_z_xyyzzz_xyyz, g_0_z_xyyzzz_xyz, g_0_z_xyyzzz_xyzz, g_0_z_xyyzzz_xzz, g_0_z_xyyzzz_xzzz, g_0_z_xyyzzz_yyy, g_0_z_xyyzzz_yyz, g_0_z_xyyzzz_yzz, g_0_z_xyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzzz_xxx[k] = -g_0_z_xyyzzz_xxx[k] * ab_x + g_0_z_xyyzzz_xxxx[k];

                g_0_z_xxyyzzz_xxy[k] = -g_0_z_xyyzzz_xxy[k] * ab_x + g_0_z_xyyzzz_xxxy[k];

                g_0_z_xxyyzzz_xxz[k] = -g_0_z_xyyzzz_xxz[k] * ab_x + g_0_z_xyyzzz_xxxz[k];

                g_0_z_xxyyzzz_xyy[k] = -g_0_z_xyyzzz_xyy[k] * ab_x + g_0_z_xyyzzz_xxyy[k];

                g_0_z_xxyyzzz_xyz[k] = -g_0_z_xyyzzz_xyz[k] * ab_x + g_0_z_xyyzzz_xxyz[k];

                g_0_z_xxyyzzz_xzz[k] = -g_0_z_xyyzzz_xzz[k] * ab_x + g_0_z_xyyzzz_xxzz[k];

                g_0_z_xxyyzzz_yyy[k] = -g_0_z_xyyzzz_yyy[k] * ab_x + g_0_z_xyyzzz_xyyy[k];

                g_0_z_xxyyzzz_yyz[k] = -g_0_z_xyyzzz_yyz[k] * ab_x + g_0_z_xyyzzz_xyyz[k];

                g_0_z_xxyyzzz_yzz[k] = -g_0_z_xyyzzz_yzz[k] * ab_x + g_0_z_xyyzzz_xyzz[k];

                g_0_z_xxyyzzz_zzz[k] = -g_0_z_xyyzzz_zzz[k] * ab_x + g_0_z_xyyzzz_xzzz[k];
            }

            /// Set up 910-920 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzzz_xxx = cbuffer.data(kf_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxy = cbuffer.data(kf_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxz = cbuffer.data(kf_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xyy = cbuffer.data(kf_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xyz = cbuffer.data(kf_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xzz = cbuffer.data(kf_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yyy = cbuffer.data(kf_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yyz = cbuffer.data(kf_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yzz = cbuffer.data(kf_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_zzz = cbuffer.data(kf_geom_01_off + 919 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzzz_xxx, g_0_z_xxyzzzz_xxy, g_0_z_xxyzzzz_xxz, g_0_z_xxyzzzz_xyy, g_0_z_xxyzzzz_xyz, g_0_z_xxyzzzz_xzz, g_0_z_xxyzzzz_yyy, g_0_z_xxyzzzz_yyz, g_0_z_xxyzzzz_yzz, g_0_z_xxyzzzz_zzz, g_0_z_xyzzzz_xxx, g_0_z_xyzzzz_xxxx, g_0_z_xyzzzz_xxxy, g_0_z_xyzzzz_xxxz, g_0_z_xyzzzz_xxy, g_0_z_xyzzzz_xxyy, g_0_z_xyzzzz_xxyz, g_0_z_xyzzzz_xxz, g_0_z_xyzzzz_xxzz, g_0_z_xyzzzz_xyy, g_0_z_xyzzzz_xyyy, g_0_z_xyzzzz_xyyz, g_0_z_xyzzzz_xyz, g_0_z_xyzzzz_xyzz, g_0_z_xyzzzz_xzz, g_0_z_xyzzzz_xzzz, g_0_z_xyzzzz_yyy, g_0_z_xyzzzz_yyz, g_0_z_xyzzzz_yzz, g_0_z_xyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzzz_xxx[k] = -g_0_z_xyzzzz_xxx[k] * ab_x + g_0_z_xyzzzz_xxxx[k];

                g_0_z_xxyzzzz_xxy[k] = -g_0_z_xyzzzz_xxy[k] * ab_x + g_0_z_xyzzzz_xxxy[k];

                g_0_z_xxyzzzz_xxz[k] = -g_0_z_xyzzzz_xxz[k] * ab_x + g_0_z_xyzzzz_xxxz[k];

                g_0_z_xxyzzzz_xyy[k] = -g_0_z_xyzzzz_xyy[k] * ab_x + g_0_z_xyzzzz_xxyy[k];

                g_0_z_xxyzzzz_xyz[k] = -g_0_z_xyzzzz_xyz[k] * ab_x + g_0_z_xyzzzz_xxyz[k];

                g_0_z_xxyzzzz_xzz[k] = -g_0_z_xyzzzz_xzz[k] * ab_x + g_0_z_xyzzzz_xxzz[k];

                g_0_z_xxyzzzz_yyy[k] = -g_0_z_xyzzzz_yyy[k] * ab_x + g_0_z_xyzzzz_xyyy[k];

                g_0_z_xxyzzzz_yyz[k] = -g_0_z_xyzzzz_yyz[k] * ab_x + g_0_z_xyzzzz_xyyz[k];

                g_0_z_xxyzzzz_yzz[k] = -g_0_z_xyzzzz_yzz[k] * ab_x + g_0_z_xyzzzz_xyzz[k];

                g_0_z_xxyzzzz_zzz[k] = -g_0_z_xyzzzz_zzz[k] * ab_x + g_0_z_xyzzzz_xzzz[k];
            }

            /// Set up 920-930 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzzz_xxx = cbuffer.data(kf_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxy = cbuffer.data(kf_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxz = cbuffer.data(kf_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xyy = cbuffer.data(kf_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xyz = cbuffer.data(kf_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xzz = cbuffer.data(kf_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yyy = cbuffer.data(kf_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yyz = cbuffer.data(kf_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yzz = cbuffer.data(kf_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_zzz = cbuffer.data(kf_geom_01_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzzz_xxx, g_0_z_xxzzzzz_xxy, g_0_z_xxzzzzz_xxz, g_0_z_xxzzzzz_xyy, g_0_z_xxzzzzz_xyz, g_0_z_xxzzzzz_xzz, g_0_z_xxzzzzz_yyy, g_0_z_xxzzzzz_yyz, g_0_z_xxzzzzz_yzz, g_0_z_xxzzzzz_zzz, g_0_z_xzzzzz_xxx, g_0_z_xzzzzz_xxxx, g_0_z_xzzzzz_xxxy, g_0_z_xzzzzz_xxxz, g_0_z_xzzzzz_xxy, g_0_z_xzzzzz_xxyy, g_0_z_xzzzzz_xxyz, g_0_z_xzzzzz_xxz, g_0_z_xzzzzz_xxzz, g_0_z_xzzzzz_xyy, g_0_z_xzzzzz_xyyy, g_0_z_xzzzzz_xyyz, g_0_z_xzzzzz_xyz, g_0_z_xzzzzz_xyzz, g_0_z_xzzzzz_xzz, g_0_z_xzzzzz_xzzz, g_0_z_xzzzzz_yyy, g_0_z_xzzzzz_yyz, g_0_z_xzzzzz_yzz, g_0_z_xzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzzz_xxx[k] = -g_0_z_xzzzzz_xxx[k] * ab_x + g_0_z_xzzzzz_xxxx[k];

                g_0_z_xxzzzzz_xxy[k] = -g_0_z_xzzzzz_xxy[k] * ab_x + g_0_z_xzzzzz_xxxy[k];

                g_0_z_xxzzzzz_xxz[k] = -g_0_z_xzzzzz_xxz[k] * ab_x + g_0_z_xzzzzz_xxxz[k];

                g_0_z_xxzzzzz_xyy[k] = -g_0_z_xzzzzz_xyy[k] * ab_x + g_0_z_xzzzzz_xxyy[k];

                g_0_z_xxzzzzz_xyz[k] = -g_0_z_xzzzzz_xyz[k] * ab_x + g_0_z_xzzzzz_xxyz[k];

                g_0_z_xxzzzzz_xzz[k] = -g_0_z_xzzzzz_xzz[k] * ab_x + g_0_z_xzzzzz_xxzz[k];

                g_0_z_xxzzzzz_yyy[k] = -g_0_z_xzzzzz_yyy[k] * ab_x + g_0_z_xzzzzz_xyyy[k];

                g_0_z_xxzzzzz_yyz[k] = -g_0_z_xzzzzz_yyz[k] * ab_x + g_0_z_xzzzzz_xyyz[k];

                g_0_z_xxzzzzz_yzz[k] = -g_0_z_xzzzzz_yzz[k] * ab_x + g_0_z_xzzzzz_xyzz[k];

                g_0_z_xxzzzzz_zzz[k] = -g_0_z_xzzzzz_zzz[k] * ab_x + g_0_z_xzzzzz_xzzz[k];
            }

            /// Set up 930-940 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 939 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyy_xxx, g_0_z_xyyyyyy_xxy, g_0_z_xyyyyyy_xxz, g_0_z_xyyyyyy_xyy, g_0_z_xyyyyyy_xyz, g_0_z_xyyyyyy_xzz, g_0_z_xyyyyyy_yyy, g_0_z_xyyyyyy_yyz, g_0_z_xyyyyyy_yzz, g_0_z_xyyyyyy_zzz, g_0_z_yyyyyy_xxx, g_0_z_yyyyyy_xxxx, g_0_z_yyyyyy_xxxy, g_0_z_yyyyyy_xxxz, g_0_z_yyyyyy_xxy, g_0_z_yyyyyy_xxyy, g_0_z_yyyyyy_xxyz, g_0_z_yyyyyy_xxz, g_0_z_yyyyyy_xxzz, g_0_z_yyyyyy_xyy, g_0_z_yyyyyy_xyyy, g_0_z_yyyyyy_xyyz, g_0_z_yyyyyy_xyz, g_0_z_yyyyyy_xyzz, g_0_z_yyyyyy_xzz, g_0_z_yyyyyy_xzzz, g_0_z_yyyyyy_yyy, g_0_z_yyyyyy_yyz, g_0_z_yyyyyy_yzz, g_0_z_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyy_xxx[k] = -g_0_z_yyyyyy_xxx[k] * ab_x + g_0_z_yyyyyy_xxxx[k];

                g_0_z_xyyyyyy_xxy[k] = -g_0_z_yyyyyy_xxy[k] * ab_x + g_0_z_yyyyyy_xxxy[k];

                g_0_z_xyyyyyy_xxz[k] = -g_0_z_yyyyyy_xxz[k] * ab_x + g_0_z_yyyyyy_xxxz[k];

                g_0_z_xyyyyyy_xyy[k] = -g_0_z_yyyyyy_xyy[k] * ab_x + g_0_z_yyyyyy_xxyy[k];

                g_0_z_xyyyyyy_xyz[k] = -g_0_z_yyyyyy_xyz[k] * ab_x + g_0_z_yyyyyy_xxyz[k];

                g_0_z_xyyyyyy_xzz[k] = -g_0_z_yyyyyy_xzz[k] * ab_x + g_0_z_yyyyyy_xxzz[k];

                g_0_z_xyyyyyy_yyy[k] = -g_0_z_yyyyyy_yyy[k] * ab_x + g_0_z_yyyyyy_xyyy[k];

                g_0_z_xyyyyyy_yyz[k] = -g_0_z_yyyyyy_yyz[k] * ab_x + g_0_z_yyyyyy_xyyz[k];

                g_0_z_xyyyyyy_yzz[k] = -g_0_z_yyyyyy_yzz[k] * ab_x + g_0_z_yyyyyy_xyzz[k];

                g_0_z_xyyyyyy_zzz[k] = -g_0_z_yyyyyy_zzz[k] * ab_x + g_0_z_yyyyyy_xzzz[k];
            }

            /// Set up 940-950 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 949 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyz_xxx, g_0_z_xyyyyyz_xxy, g_0_z_xyyyyyz_xxz, g_0_z_xyyyyyz_xyy, g_0_z_xyyyyyz_xyz, g_0_z_xyyyyyz_xzz, g_0_z_xyyyyyz_yyy, g_0_z_xyyyyyz_yyz, g_0_z_xyyyyyz_yzz, g_0_z_xyyyyyz_zzz, g_0_z_yyyyyz_xxx, g_0_z_yyyyyz_xxxx, g_0_z_yyyyyz_xxxy, g_0_z_yyyyyz_xxxz, g_0_z_yyyyyz_xxy, g_0_z_yyyyyz_xxyy, g_0_z_yyyyyz_xxyz, g_0_z_yyyyyz_xxz, g_0_z_yyyyyz_xxzz, g_0_z_yyyyyz_xyy, g_0_z_yyyyyz_xyyy, g_0_z_yyyyyz_xyyz, g_0_z_yyyyyz_xyz, g_0_z_yyyyyz_xyzz, g_0_z_yyyyyz_xzz, g_0_z_yyyyyz_xzzz, g_0_z_yyyyyz_yyy, g_0_z_yyyyyz_yyz, g_0_z_yyyyyz_yzz, g_0_z_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyz_xxx[k] = -g_0_z_yyyyyz_xxx[k] * ab_x + g_0_z_yyyyyz_xxxx[k];

                g_0_z_xyyyyyz_xxy[k] = -g_0_z_yyyyyz_xxy[k] * ab_x + g_0_z_yyyyyz_xxxy[k];

                g_0_z_xyyyyyz_xxz[k] = -g_0_z_yyyyyz_xxz[k] * ab_x + g_0_z_yyyyyz_xxxz[k];

                g_0_z_xyyyyyz_xyy[k] = -g_0_z_yyyyyz_xyy[k] * ab_x + g_0_z_yyyyyz_xxyy[k];

                g_0_z_xyyyyyz_xyz[k] = -g_0_z_yyyyyz_xyz[k] * ab_x + g_0_z_yyyyyz_xxyz[k];

                g_0_z_xyyyyyz_xzz[k] = -g_0_z_yyyyyz_xzz[k] * ab_x + g_0_z_yyyyyz_xxzz[k];

                g_0_z_xyyyyyz_yyy[k] = -g_0_z_yyyyyz_yyy[k] * ab_x + g_0_z_yyyyyz_xyyy[k];

                g_0_z_xyyyyyz_yyz[k] = -g_0_z_yyyyyz_yyz[k] * ab_x + g_0_z_yyyyyz_xyyz[k];

                g_0_z_xyyyyyz_yzz[k] = -g_0_z_yyyyyz_yzz[k] * ab_x + g_0_z_yyyyyz_xyzz[k];

                g_0_z_xyyyyyz_zzz[k] = -g_0_z_yyyyyz_zzz[k] * ab_x + g_0_z_yyyyyz_xzzz[k];
            }

            /// Set up 950-960 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyzz_xxx, g_0_z_xyyyyzz_xxy, g_0_z_xyyyyzz_xxz, g_0_z_xyyyyzz_xyy, g_0_z_xyyyyzz_xyz, g_0_z_xyyyyzz_xzz, g_0_z_xyyyyzz_yyy, g_0_z_xyyyyzz_yyz, g_0_z_xyyyyzz_yzz, g_0_z_xyyyyzz_zzz, g_0_z_yyyyzz_xxx, g_0_z_yyyyzz_xxxx, g_0_z_yyyyzz_xxxy, g_0_z_yyyyzz_xxxz, g_0_z_yyyyzz_xxy, g_0_z_yyyyzz_xxyy, g_0_z_yyyyzz_xxyz, g_0_z_yyyyzz_xxz, g_0_z_yyyyzz_xxzz, g_0_z_yyyyzz_xyy, g_0_z_yyyyzz_xyyy, g_0_z_yyyyzz_xyyz, g_0_z_yyyyzz_xyz, g_0_z_yyyyzz_xyzz, g_0_z_yyyyzz_xzz, g_0_z_yyyyzz_xzzz, g_0_z_yyyyzz_yyy, g_0_z_yyyyzz_yyz, g_0_z_yyyyzz_yzz, g_0_z_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyzz_xxx[k] = -g_0_z_yyyyzz_xxx[k] * ab_x + g_0_z_yyyyzz_xxxx[k];

                g_0_z_xyyyyzz_xxy[k] = -g_0_z_yyyyzz_xxy[k] * ab_x + g_0_z_yyyyzz_xxxy[k];

                g_0_z_xyyyyzz_xxz[k] = -g_0_z_yyyyzz_xxz[k] * ab_x + g_0_z_yyyyzz_xxxz[k];

                g_0_z_xyyyyzz_xyy[k] = -g_0_z_yyyyzz_xyy[k] * ab_x + g_0_z_yyyyzz_xxyy[k];

                g_0_z_xyyyyzz_xyz[k] = -g_0_z_yyyyzz_xyz[k] * ab_x + g_0_z_yyyyzz_xxyz[k];

                g_0_z_xyyyyzz_xzz[k] = -g_0_z_yyyyzz_xzz[k] * ab_x + g_0_z_yyyyzz_xxzz[k];

                g_0_z_xyyyyzz_yyy[k] = -g_0_z_yyyyzz_yyy[k] * ab_x + g_0_z_yyyyzz_xyyy[k];

                g_0_z_xyyyyzz_yyz[k] = -g_0_z_yyyyzz_yyz[k] * ab_x + g_0_z_yyyyzz_xyyz[k];

                g_0_z_xyyyyzz_yzz[k] = -g_0_z_yyyyzz_yzz[k] * ab_x + g_0_z_yyyyzz_xyzz[k];

                g_0_z_xyyyyzz_zzz[k] = -g_0_z_yyyyzz_zzz[k] * ab_x + g_0_z_yyyyzz_xzzz[k];
            }

            /// Set up 960-970 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 969 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzzz_xxx, g_0_z_xyyyzzz_xxy, g_0_z_xyyyzzz_xxz, g_0_z_xyyyzzz_xyy, g_0_z_xyyyzzz_xyz, g_0_z_xyyyzzz_xzz, g_0_z_xyyyzzz_yyy, g_0_z_xyyyzzz_yyz, g_0_z_xyyyzzz_yzz, g_0_z_xyyyzzz_zzz, g_0_z_yyyzzz_xxx, g_0_z_yyyzzz_xxxx, g_0_z_yyyzzz_xxxy, g_0_z_yyyzzz_xxxz, g_0_z_yyyzzz_xxy, g_0_z_yyyzzz_xxyy, g_0_z_yyyzzz_xxyz, g_0_z_yyyzzz_xxz, g_0_z_yyyzzz_xxzz, g_0_z_yyyzzz_xyy, g_0_z_yyyzzz_xyyy, g_0_z_yyyzzz_xyyz, g_0_z_yyyzzz_xyz, g_0_z_yyyzzz_xyzz, g_0_z_yyyzzz_xzz, g_0_z_yyyzzz_xzzz, g_0_z_yyyzzz_yyy, g_0_z_yyyzzz_yyz, g_0_z_yyyzzz_yzz, g_0_z_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzzz_xxx[k] = -g_0_z_yyyzzz_xxx[k] * ab_x + g_0_z_yyyzzz_xxxx[k];

                g_0_z_xyyyzzz_xxy[k] = -g_0_z_yyyzzz_xxy[k] * ab_x + g_0_z_yyyzzz_xxxy[k];

                g_0_z_xyyyzzz_xxz[k] = -g_0_z_yyyzzz_xxz[k] * ab_x + g_0_z_yyyzzz_xxxz[k];

                g_0_z_xyyyzzz_xyy[k] = -g_0_z_yyyzzz_xyy[k] * ab_x + g_0_z_yyyzzz_xxyy[k];

                g_0_z_xyyyzzz_xyz[k] = -g_0_z_yyyzzz_xyz[k] * ab_x + g_0_z_yyyzzz_xxyz[k];

                g_0_z_xyyyzzz_xzz[k] = -g_0_z_yyyzzz_xzz[k] * ab_x + g_0_z_yyyzzz_xxzz[k];

                g_0_z_xyyyzzz_yyy[k] = -g_0_z_yyyzzz_yyy[k] * ab_x + g_0_z_yyyzzz_xyyy[k];

                g_0_z_xyyyzzz_yyz[k] = -g_0_z_yyyzzz_yyz[k] * ab_x + g_0_z_yyyzzz_xyyz[k];

                g_0_z_xyyyzzz_yzz[k] = -g_0_z_yyyzzz_yzz[k] * ab_x + g_0_z_yyyzzz_xyzz[k];

                g_0_z_xyyyzzz_zzz[k] = -g_0_z_yyyzzz_zzz[k] * ab_x + g_0_z_yyyzzz_xzzz[k];
            }

            /// Set up 970-980 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzzz_xxx, g_0_z_xyyzzzz_xxy, g_0_z_xyyzzzz_xxz, g_0_z_xyyzzzz_xyy, g_0_z_xyyzzzz_xyz, g_0_z_xyyzzzz_xzz, g_0_z_xyyzzzz_yyy, g_0_z_xyyzzzz_yyz, g_0_z_xyyzzzz_yzz, g_0_z_xyyzzzz_zzz, g_0_z_yyzzzz_xxx, g_0_z_yyzzzz_xxxx, g_0_z_yyzzzz_xxxy, g_0_z_yyzzzz_xxxz, g_0_z_yyzzzz_xxy, g_0_z_yyzzzz_xxyy, g_0_z_yyzzzz_xxyz, g_0_z_yyzzzz_xxz, g_0_z_yyzzzz_xxzz, g_0_z_yyzzzz_xyy, g_0_z_yyzzzz_xyyy, g_0_z_yyzzzz_xyyz, g_0_z_yyzzzz_xyz, g_0_z_yyzzzz_xyzz, g_0_z_yyzzzz_xzz, g_0_z_yyzzzz_xzzz, g_0_z_yyzzzz_yyy, g_0_z_yyzzzz_yyz, g_0_z_yyzzzz_yzz, g_0_z_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzzz_xxx[k] = -g_0_z_yyzzzz_xxx[k] * ab_x + g_0_z_yyzzzz_xxxx[k];

                g_0_z_xyyzzzz_xxy[k] = -g_0_z_yyzzzz_xxy[k] * ab_x + g_0_z_yyzzzz_xxxy[k];

                g_0_z_xyyzzzz_xxz[k] = -g_0_z_yyzzzz_xxz[k] * ab_x + g_0_z_yyzzzz_xxxz[k];

                g_0_z_xyyzzzz_xyy[k] = -g_0_z_yyzzzz_xyy[k] * ab_x + g_0_z_yyzzzz_xxyy[k];

                g_0_z_xyyzzzz_xyz[k] = -g_0_z_yyzzzz_xyz[k] * ab_x + g_0_z_yyzzzz_xxyz[k];

                g_0_z_xyyzzzz_xzz[k] = -g_0_z_yyzzzz_xzz[k] * ab_x + g_0_z_yyzzzz_xxzz[k];

                g_0_z_xyyzzzz_yyy[k] = -g_0_z_yyzzzz_yyy[k] * ab_x + g_0_z_yyzzzz_xyyy[k];

                g_0_z_xyyzzzz_yyz[k] = -g_0_z_yyzzzz_yyz[k] * ab_x + g_0_z_yyzzzz_xyyz[k];

                g_0_z_xyyzzzz_yzz[k] = -g_0_z_yyzzzz_yzz[k] * ab_x + g_0_z_yyzzzz_xyzz[k];

                g_0_z_xyyzzzz_zzz[k] = -g_0_z_yyzzzz_zzz[k] * ab_x + g_0_z_yyzzzz_xzzz[k];
            }

            /// Set up 980-990 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzzz_xxx, g_0_z_xyzzzzz_xxy, g_0_z_xyzzzzz_xxz, g_0_z_xyzzzzz_xyy, g_0_z_xyzzzzz_xyz, g_0_z_xyzzzzz_xzz, g_0_z_xyzzzzz_yyy, g_0_z_xyzzzzz_yyz, g_0_z_xyzzzzz_yzz, g_0_z_xyzzzzz_zzz, g_0_z_yzzzzz_xxx, g_0_z_yzzzzz_xxxx, g_0_z_yzzzzz_xxxy, g_0_z_yzzzzz_xxxz, g_0_z_yzzzzz_xxy, g_0_z_yzzzzz_xxyy, g_0_z_yzzzzz_xxyz, g_0_z_yzzzzz_xxz, g_0_z_yzzzzz_xxzz, g_0_z_yzzzzz_xyy, g_0_z_yzzzzz_xyyy, g_0_z_yzzzzz_xyyz, g_0_z_yzzzzz_xyz, g_0_z_yzzzzz_xyzz, g_0_z_yzzzzz_xzz, g_0_z_yzzzzz_xzzz, g_0_z_yzzzzz_yyy, g_0_z_yzzzzz_yyz, g_0_z_yzzzzz_yzz, g_0_z_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzzz_xxx[k] = -g_0_z_yzzzzz_xxx[k] * ab_x + g_0_z_yzzzzz_xxxx[k];

                g_0_z_xyzzzzz_xxy[k] = -g_0_z_yzzzzz_xxy[k] * ab_x + g_0_z_yzzzzz_xxxy[k];

                g_0_z_xyzzzzz_xxz[k] = -g_0_z_yzzzzz_xxz[k] * ab_x + g_0_z_yzzzzz_xxxz[k];

                g_0_z_xyzzzzz_xyy[k] = -g_0_z_yzzzzz_xyy[k] * ab_x + g_0_z_yzzzzz_xxyy[k];

                g_0_z_xyzzzzz_xyz[k] = -g_0_z_yzzzzz_xyz[k] * ab_x + g_0_z_yzzzzz_xxyz[k];

                g_0_z_xyzzzzz_xzz[k] = -g_0_z_yzzzzz_xzz[k] * ab_x + g_0_z_yzzzzz_xxzz[k];

                g_0_z_xyzzzzz_yyy[k] = -g_0_z_yzzzzz_yyy[k] * ab_x + g_0_z_yzzzzz_xyyy[k];

                g_0_z_xyzzzzz_yyz[k] = -g_0_z_yzzzzz_yyz[k] * ab_x + g_0_z_yzzzzz_xyyz[k];

                g_0_z_xyzzzzz_yzz[k] = -g_0_z_yzzzzz_yzz[k] * ab_x + g_0_z_yzzzzz_xyzz[k];

                g_0_z_xyzzzzz_zzz[k] = -g_0_z_yzzzzz_zzz[k] * ab_x + g_0_z_yzzzzz_xzzz[k];
            }

            /// Set up 990-1000 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 999 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzzz_xxx, g_0_z_xzzzzzz_xxy, g_0_z_xzzzzzz_xxz, g_0_z_xzzzzzz_xyy, g_0_z_xzzzzzz_xyz, g_0_z_xzzzzzz_xzz, g_0_z_xzzzzzz_yyy, g_0_z_xzzzzzz_yyz, g_0_z_xzzzzzz_yzz, g_0_z_xzzzzzz_zzz, g_0_z_zzzzzz_xxx, g_0_z_zzzzzz_xxxx, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxy, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xyy, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_yyy, g_0_z_zzzzzz_yyz, g_0_z_zzzzzz_yzz, g_0_z_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzzz_xxx[k] = -g_0_z_zzzzzz_xxx[k] * ab_x + g_0_z_zzzzzz_xxxx[k];

                g_0_z_xzzzzzz_xxy[k] = -g_0_z_zzzzzz_xxy[k] * ab_x + g_0_z_zzzzzz_xxxy[k];

                g_0_z_xzzzzzz_xxz[k] = -g_0_z_zzzzzz_xxz[k] * ab_x + g_0_z_zzzzzz_xxxz[k];

                g_0_z_xzzzzzz_xyy[k] = -g_0_z_zzzzzz_xyy[k] * ab_x + g_0_z_zzzzzz_xxyy[k];

                g_0_z_xzzzzzz_xyz[k] = -g_0_z_zzzzzz_xyz[k] * ab_x + g_0_z_zzzzzz_xxyz[k];

                g_0_z_xzzzzzz_xzz[k] = -g_0_z_zzzzzz_xzz[k] * ab_x + g_0_z_zzzzzz_xxzz[k];

                g_0_z_xzzzzzz_yyy[k] = -g_0_z_zzzzzz_yyy[k] * ab_x + g_0_z_zzzzzz_xyyy[k];

                g_0_z_xzzzzzz_yyz[k] = -g_0_z_zzzzzz_yyz[k] * ab_x + g_0_z_zzzzzz_xyyz[k];

                g_0_z_xzzzzzz_yzz[k] = -g_0_z_zzzzzz_yzz[k] * ab_x + g_0_z_zzzzzz_xyzz[k];

                g_0_z_xzzzzzz_zzz[k] = -g_0_z_zzzzzz_zzz[k] * ab_x + g_0_z_zzzzzz_xzzz[k];
            }

            /// Set up 1000-1010 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyy_xxx = cbuffer.data(kf_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxy = cbuffer.data(kf_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxz = cbuffer.data(kf_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xyy = cbuffer.data(kf_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xyz = cbuffer.data(kf_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xzz = cbuffer.data(kf_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yyy = cbuffer.data(kf_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yyz = cbuffer.data(kf_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yzz = cbuffer.data(kf_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_zzz = cbuffer.data(kf_geom_01_off + 1009 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyy_xxx, g_0_z_yyyyyy_xxxy, g_0_z_yyyyyy_xxy, g_0_z_yyyyyy_xxyy, g_0_z_yyyyyy_xxyz, g_0_z_yyyyyy_xxz, g_0_z_yyyyyy_xyy, g_0_z_yyyyyy_xyyy, g_0_z_yyyyyy_xyyz, g_0_z_yyyyyy_xyz, g_0_z_yyyyyy_xyzz, g_0_z_yyyyyy_xzz, g_0_z_yyyyyy_yyy, g_0_z_yyyyyy_yyyy, g_0_z_yyyyyy_yyyz, g_0_z_yyyyyy_yyz, g_0_z_yyyyyy_yyzz, g_0_z_yyyyyy_yzz, g_0_z_yyyyyy_yzzz, g_0_z_yyyyyy_zzz, g_0_z_yyyyyyy_xxx, g_0_z_yyyyyyy_xxy, g_0_z_yyyyyyy_xxz, g_0_z_yyyyyyy_xyy, g_0_z_yyyyyyy_xyz, g_0_z_yyyyyyy_xzz, g_0_z_yyyyyyy_yyy, g_0_z_yyyyyyy_yyz, g_0_z_yyyyyyy_yzz, g_0_z_yyyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyy_xxx[k] = -g_0_z_yyyyyy_xxx[k] * ab_y + g_0_z_yyyyyy_xxxy[k];

                g_0_z_yyyyyyy_xxy[k] = -g_0_z_yyyyyy_xxy[k] * ab_y + g_0_z_yyyyyy_xxyy[k];

                g_0_z_yyyyyyy_xxz[k] = -g_0_z_yyyyyy_xxz[k] * ab_y + g_0_z_yyyyyy_xxyz[k];

                g_0_z_yyyyyyy_xyy[k] = -g_0_z_yyyyyy_xyy[k] * ab_y + g_0_z_yyyyyy_xyyy[k];

                g_0_z_yyyyyyy_xyz[k] = -g_0_z_yyyyyy_xyz[k] * ab_y + g_0_z_yyyyyy_xyyz[k];

                g_0_z_yyyyyyy_xzz[k] = -g_0_z_yyyyyy_xzz[k] * ab_y + g_0_z_yyyyyy_xyzz[k];

                g_0_z_yyyyyyy_yyy[k] = -g_0_z_yyyyyy_yyy[k] * ab_y + g_0_z_yyyyyy_yyyy[k];

                g_0_z_yyyyyyy_yyz[k] = -g_0_z_yyyyyy_yyz[k] * ab_y + g_0_z_yyyyyy_yyyz[k];

                g_0_z_yyyyyyy_yzz[k] = -g_0_z_yyyyyy_yzz[k] * ab_y + g_0_z_yyyyyy_yyzz[k];

                g_0_z_yyyyyyy_zzz[k] = -g_0_z_yyyyyy_zzz[k] * ab_y + g_0_z_yyyyyy_yzzz[k];
            }

            /// Set up 1010-1020 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyz_xxx = cbuffer.data(kf_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxy = cbuffer.data(kf_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxz = cbuffer.data(kf_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xyy = cbuffer.data(kf_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xyz = cbuffer.data(kf_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xzz = cbuffer.data(kf_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yyy = cbuffer.data(kf_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yyz = cbuffer.data(kf_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yzz = cbuffer.data(kf_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_zzz = cbuffer.data(kf_geom_01_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyyz_xxx, g_0_z_yyyyyyz_xxy, g_0_z_yyyyyyz_xxz, g_0_z_yyyyyyz_xyy, g_0_z_yyyyyyz_xyz, g_0_z_yyyyyyz_xzz, g_0_z_yyyyyyz_yyy, g_0_z_yyyyyyz_yyz, g_0_z_yyyyyyz_yzz, g_0_z_yyyyyyz_zzz, g_0_z_yyyyyz_xxx, g_0_z_yyyyyz_xxxy, g_0_z_yyyyyz_xxy, g_0_z_yyyyyz_xxyy, g_0_z_yyyyyz_xxyz, g_0_z_yyyyyz_xxz, g_0_z_yyyyyz_xyy, g_0_z_yyyyyz_xyyy, g_0_z_yyyyyz_xyyz, g_0_z_yyyyyz_xyz, g_0_z_yyyyyz_xyzz, g_0_z_yyyyyz_xzz, g_0_z_yyyyyz_yyy, g_0_z_yyyyyz_yyyy, g_0_z_yyyyyz_yyyz, g_0_z_yyyyyz_yyz, g_0_z_yyyyyz_yyzz, g_0_z_yyyyyz_yzz, g_0_z_yyyyyz_yzzz, g_0_z_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyz_xxx[k] = -g_0_z_yyyyyz_xxx[k] * ab_y + g_0_z_yyyyyz_xxxy[k];

                g_0_z_yyyyyyz_xxy[k] = -g_0_z_yyyyyz_xxy[k] * ab_y + g_0_z_yyyyyz_xxyy[k];

                g_0_z_yyyyyyz_xxz[k] = -g_0_z_yyyyyz_xxz[k] * ab_y + g_0_z_yyyyyz_xxyz[k];

                g_0_z_yyyyyyz_xyy[k] = -g_0_z_yyyyyz_xyy[k] * ab_y + g_0_z_yyyyyz_xyyy[k];

                g_0_z_yyyyyyz_xyz[k] = -g_0_z_yyyyyz_xyz[k] * ab_y + g_0_z_yyyyyz_xyyz[k];

                g_0_z_yyyyyyz_xzz[k] = -g_0_z_yyyyyz_xzz[k] * ab_y + g_0_z_yyyyyz_xyzz[k];

                g_0_z_yyyyyyz_yyy[k] = -g_0_z_yyyyyz_yyy[k] * ab_y + g_0_z_yyyyyz_yyyy[k];

                g_0_z_yyyyyyz_yyz[k] = -g_0_z_yyyyyz_yyz[k] * ab_y + g_0_z_yyyyyz_yyyz[k];

                g_0_z_yyyyyyz_yzz[k] = -g_0_z_yyyyyz_yzz[k] * ab_y + g_0_z_yyyyyz_yyzz[k];

                g_0_z_yyyyyyz_zzz[k] = -g_0_z_yyyyyz_zzz[k] * ab_y + g_0_z_yyyyyz_yzzz[k];
            }

            /// Set up 1020-1030 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyzz_xxx = cbuffer.data(kf_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxy = cbuffer.data(kf_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxz = cbuffer.data(kf_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xyy = cbuffer.data(kf_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xyz = cbuffer.data(kf_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xzz = cbuffer.data(kf_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yyy = cbuffer.data(kf_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yyz = cbuffer.data(kf_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yzz = cbuffer.data(kf_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_zzz = cbuffer.data(kf_geom_01_off + 1029 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyzz_xxx, g_0_z_yyyyyzz_xxy, g_0_z_yyyyyzz_xxz, g_0_z_yyyyyzz_xyy, g_0_z_yyyyyzz_xyz, g_0_z_yyyyyzz_xzz, g_0_z_yyyyyzz_yyy, g_0_z_yyyyyzz_yyz, g_0_z_yyyyyzz_yzz, g_0_z_yyyyyzz_zzz, g_0_z_yyyyzz_xxx, g_0_z_yyyyzz_xxxy, g_0_z_yyyyzz_xxy, g_0_z_yyyyzz_xxyy, g_0_z_yyyyzz_xxyz, g_0_z_yyyyzz_xxz, g_0_z_yyyyzz_xyy, g_0_z_yyyyzz_xyyy, g_0_z_yyyyzz_xyyz, g_0_z_yyyyzz_xyz, g_0_z_yyyyzz_xyzz, g_0_z_yyyyzz_xzz, g_0_z_yyyyzz_yyy, g_0_z_yyyyzz_yyyy, g_0_z_yyyyzz_yyyz, g_0_z_yyyyzz_yyz, g_0_z_yyyyzz_yyzz, g_0_z_yyyyzz_yzz, g_0_z_yyyyzz_yzzz, g_0_z_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyzz_xxx[k] = -g_0_z_yyyyzz_xxx[k] * ab_y + g_0_z_yyyyzz_xxxy[k];

                g_0_z_yyyyyzz_xxy[k] = -g_0_z_yyyyzz_xxy[k] * ab_y + g_0_z_yyyyzz_xxyy[k];

                g_0_z_yyyyyzz_xxz[k] = -g_0_z_yyyyzz_xxz[k] * ab_y + g_0_z_yyyyzz_xxyz[k];

                g_0_z_yyyyyzz_xyy[k] = -g_0_z_yyyyzz_xyy[k] * ab_y + g_0_z_yyyyzz_xyyy[k];

                g_0_z_yyyyyzz_xyz[k] = -g_0_z_yyyyzz_xyz[k] * ab_y + g_0_z_yyyyzz_xyyz[k];

                g_0_z_yyyyyzz_xzz[k] = -g_0_z_yyyyzz_xzz[k] * ab_y + g_0_z_yyyyzz_xyzz[k];

                g_0_z_yyyyyzz_yyy[k] = -g_0_z_yyyyzz_yyy[k] * ab_y + g_0_z_yyyyzz_yyyy[k];

                g_0_z_yyyyyzz_yyz[k] = -g_0_z_yyyyzz_yyz[k] * ab_y + g_0_z_yyyyzz_yyyz[k];

                g_0_z_yyyyyzz_yzz[k] = -g_0_z_yyyyzz_yzz[k] * ab_y + g_0_z_yyyyzz_yyzz[k];

                g_0_z_yyyyyzz_zzz[k] = -g_0_z_yyyyzz_zzz[k] * ab_y + g_0_z_yyyyzz_yzzz[k];
            }

            /// Set up 1030-1040 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzzz_xxx = cbuffer.data(kf_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxy = cbuffer.data(kf_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxz = cbuffer.data(kf_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xyy = cbuffer.data(kf_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xyz = cbuffer.data(kf_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xzz = cbuffer.data(kf_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yyy = cbuffer.data(kf_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yyz = cbuffer.data(kf_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yzz = cbuffer.data(kf_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_zzz = cbuffer.data(kf_geom_01_off + 1039 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzzz_xxx, g_0_z_yyyyzzz_xxy, g_0_z_yyyyzzz_xxz, g_0_z_yyyyzzz_xyy, g_0_z_yyyyzzz_xyz, g_0_z_yyyyzzz_xzz, g_0_z_yyyyzzz_yyy, g_0_z_yyyyzzz_yyz, g_0_z_yyyyzzz_yzz, g_0_z_yyyyzzz_zzz, g_0_z_yyyzzz_xxx, g_0_z_yyyzzz_xxxy, g_0_z_yyyzzz_xxy, g_0_z_yyyzzz_xxyy, g_0_z_yyyzzz_xxyz, g_0_z_yyyzzz_xxz, g_0_z_yyyzzz_xyy, g_0_z_yyyzzz_xyyy, g_0_z_yyyzzz_xyyz, g_0_z_yyyzzz_xyz, g_0_z_yyyzzz_xyzz, g_0_z_yyyzzz_xzz, g_0_z_yyyzzz_yyy, g_0_z_yyyzzz_yyyy, g_0_z_yyyzzz_yyyz, g_0_z_yyyzzz_yyz, g_0_z_yyyzzz_yyzz, g_0_z_yyyzzz_yzz, g_0_z_yyyzzz_yzzz, g_0_z_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzzz_xxx[k] = -g_0_z_yyyzzz_xxx[k] * ab_y + g_0_z_yyyzzz_xxxy[k];

                g_0_z_yyyyzzz_xxy[k] = -g_0_z_yyyzzz_xxy[k] * ab_y + g_0_z_yyyzzz_xxyy[k];

                g_0_z_yyyyzzz_xxz[k] = -g_0_z_yyyzzz_xxz[k] * ab_y + g_0_z_yyyzzz_xxyz[k];

                g_0_z_yyyyzzz_xyy[k] = -g_0_z_yyyzzz_xyy[k] * ab_y + g_0_z_yyyzzz_xyyy[k];

                g_0_z_yyyyzzz_xyz[k] = -g_0_z_yyyzzz_xyz[k] * ab_y + g_0_z_yyyzzz_xyyz[k];

                g_0_z_yyyyzzz_xzz[k] = -g_0_z_yyyzzz_xzz[k] * ab_y + g_0_z_yyyzzz_xyzz[k];

                g_0_z_yyyyzzz_yyy[k] = -g_0_z_yyyzzz_yyy[k] * ab_y + g_0_z_yyyzzz_yyyy[k];

                g_0_z_yyyyzzz_yyz[k] = -g_0_z_yyyzzz_yyz[k] * ab_y + g_0_z_yyyzzz_yyyz[k];

                g_0_z_yyyyzzz_yzz[k] = -g_0_z_yyyzzz_yzz[k] * ab_y + g_0_z_yyyzzz_yyzz[k];

                g_0_z_yyyyzzz_zzz[k] = -g_0_z_yyyzzz_zzz[k] * ab_y + g_0_z_yyyzzz_yzzz[k];
            }

            /// Set up 1040-1050 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzzz_xxx = cbuffer.data(kf_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxy = cbuffer.data(kf_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxz = cbuffer.data(kf_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xyy = cbuffer.data(kf_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xyz = cbuffer.data(kf_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xzz = cbuffer.data(kf_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yyy = cbuffer.data(kf_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yyz = cbuffer.data(kf_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yzz = cbuffer.data(kf_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_zzz = cbuffer.data(kf_geom_01_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzzz_xxx, g_0_z_yyyzzzz_xxy, g_0_z_yyyzzzz_xxz, g_0_z_yyyzzzz_xyy, g_0_z_yyyzzzz_xyz, g_0_z_yyyzzzz_xzz, g_0_z_yyyzzzz_yyy, g_0_z_yyyzzzz_yyz, g_0_z_yyyzzzz_yzz, g_0_z_yyyzzzz_zzz, g_0_z_yyzzzz_xxx, g_0_z_yyzzzz_xxxy, g_0_z_yyzzzz_xxy, g_0_z_yyzzzz_xxyy, g_0_z_yyzzzz_xxyz, g_0_z_yyzzzz_xxz, g_0_z_yyzzzz_xyy, g_0_z_yyzzzz_xyyy, g_0_z_yyzzzz_xyyz, g_0_z_yyzzzz_xyz, g_0_z_yyzzzz_xyzz, g_0_z_yyzzzz_xzz, g_0_z_yyzzzz_yyy, g_0_z_yyzzzz_yyyy, g_0_z_yyzzzz_yyyz, g_0_z_yyzzzz_yyz, g_0_z_yyzzzz_yyzz, g_0_z_yyzzzz_yzz, g_0_z_yyzzzz_yzzz, g_0_z_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzzz_xxx[k] = -g_0_z_yyzzzz_xxx[k] * ab_y + g_0_z_yyzzzz_xxxy[k];

                g_0_z_yyyzzzz_xxy[k] = -g_0_z_yyzzzz_xxy[k] * ab_y + g_0_z_yyzzzz_xxyy[k];

                g_0_z_yyyzzzz_xxz[k] = -g_0_z_yyzzzz_xxz[k] * ab_y + g_0_z_yyzzzz_xxyz[k];

                g_0_z_yyyzzzz_xyy[k] = -g_0_z_yyzzzz_xyy[k] * ab_y + g_0_z_yyzzzz_xyyy[k];

                g_0_z_yyyzzzz_xyz[k] = -g_0_z_yyzzzz_xyz[k] * ab_y + g_0_z_yyzzzz_xyyz[k];

                g_0_z_yyyzzzz_xzz[k] = -g_0_z_yyzzzz_xzz[k] * ab_y + g_0_z_yyzzzz_xyzz[k];

                g_0_z_yyyzzzz_yyy[k] = -g_0_z_yyzzzz_yyy[k] * ab_y + g_0_z_yyzzzz_yyyy[k];

                g_0_z_yyyzzzz_yyz[k] = -g_0_z_yyzzzz_yyz[k] * ab_y + g_0_z_yyzzzz_yyyz[k];

                g_0_z_yyyzzzz_yzz[k] = -g_0_z_yyzzzz_yzz[k] * ab_y + g_0_z_yyzzzz_yyzz[k];

                g_0_z_yyyzzzz_zzz[k] = -g_0_z_yyzzzz_zzz[k] * ab_y + g_0_z_yyzzzz_yzzz[k];
            }

            /// Set up 1050-1060 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzzz_xxx = cbuffer.data(kf_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxy = cbuffer.data(kf_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxz = cbuffer.data(kf_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xyy = cbuffer.data(kf_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xyz = cbuffer.data(kf_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xzz = cbuffer.data(kf_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yyy = cbuffer.data(kf_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yyz = cbuffer.data(kf_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yzz = cbuffer.data(kf_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_zzz = cbuffer.data(kf_geom_01_off + 1059 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzzz_xxx, g_0_z_yyzzzzz_xxy, g_0_z_yyzzzzz_xxz, g_0_z_yyzzzzz_xyy, g_0_z_yyzzzzz_xyz, g_0_z_yyzzzzz_xzz, g_0_z_yyzzzzz_yyy, g_0_z_yyzzzzz_yyz, g_0_z_yyzzzzz_yzz, g_0_z_yyzzzzz_zzz, g_0_z_yzzzzz_xxx, g_0_z_yzzzzz_xxxy, g_0_z_yzzzzz_xxy, g_0_z_yzzzzz_xxyy, g_0_z_yzzzzz_xxyz, g_0_z_yzzzzz_xxz, g_0_z_yzzzzz_xyy, g_0_z_yzzzzz_xyyy, g_0_z_yzzzzz_xyyz, g_0_z_yzzzzz_xyz, g_0_z_yzzzzz_xyzz, g_0_z_yzzzzz_xzz, g_0_z_yzzzzz_yyy, g_0_z_yzzzzz_yyyy, g_0_z_yzzzzz_yyyz, g_0_z_yzzzzz_yyz, g_0_z_yzzzzz_yyzz, g_0_z_yzzzzz_yzz, g_0_z_yzzzzz_yzzz, g_0_z_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzzz_xxx[k] = -g_0_z_yzzzzz_xxx[k] * ab_y + g_0_z_yzzzzz_xxxy[k];

                g_0_z_yyzzzzz_xxy[k] = -g_0_z_yzzzzz_xxy[k] * ab_y + g_0_z_yzzzzz_xxyy[k];

                g_0_z_yyzzzzz_xxz[k] = -g_0_z_yzzzzz_xxz[k] * ab_y + g_0_z_yzzzzz_xxyz[k];

                g_0_z_yyzzzzz_xyy[k] = -g_0_z_yzzzzz_xyy[k] * ab_y + g_0_z_yzzzzz_xyyy[k];

                g_0_z_yyzzzzz_xyz[k] = -g_0_z_yzzzzz_xyz[k] * ab_y + g_0_z_yzzzzz_xyyz[k];

                g_0_z_yyzzzzz_xzz[k] = -g_0_z_yzzzzz_xzz[k] * ab_y + g_0_z_yzzzzz_xyzz[k];

                g_0_z_yyzzzzz_yyy[k] = -g_0_z_yzzzzz_yyy[k] * ab_y + g_0_z_yzzzzz_yyyy[k];

                g_0_z_yyzzzzz_yyz[k] = -g_0_z_yzzzzz_yyz[k] * ab_y + g_0_z_yzzzzz_yyyz[k];

                g_0_z_yyzzzzz_yzz[k] = -g_0_z_yzzzzz_yzz[k] * ab_y + g_0_z_yzzzzz_yyzz[k];

                g_0_z_yyzzzzz_zzz[k] = -g_0_z_yzzzzz_zzz[k] * ab_y + g_0_z_yzzzzz_yzzz[k];
            }

            /// Set up 1060-1070 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 1069 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzzz_xxx, g_0_z_yzzzzzz_xxy, g_0_z_yzzzzzz_xxz, g_0_z_yzzzzzz_xyy, g_0_z_yzzzzzz_xyz, g_0_z_yzzzzzz_xzz, g_0_z_yzzzzzz_yyy, g_0_z_yzzzzzz_yyz, g_0_z_yzzzzzz_yzz, g_0_z_yzzzzzz_zzz, g_0_z_zzzzzz_xxx, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxy, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxz, g_0_z_zzzzzz_xyy, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xzz, g_0_z_zzzzzz_yyy, g_0_z_zzzzzz_yyyy, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzzz_xxx[k] = -g_0_z_zzzzzz_xxx[k] * ab_y + g_0_z_zzzzzz_xxxy[k];

                g_0_z_yzzzzzz_xxy[k] = -g_0_z_zzzzzz_xxy[k] * ab_y + g_0_z_zzzzzz_xxyy[k];

                g_0_z_yzzzzzz_xxz[k] = -g_0_z_zzzzzz_xxz[k] * ab_y + g_0_z_zzzzzz_xxyz[k];

                g_0_z_yzzzzzz_xyy[k] = -g_0_z_zzzzzz_xyy[k] * ab_y + g_0_z_zzzzzz_xyyy[k];

                g_0_z_yzzzzzz_xyz[k] = -g_0_z_zzzzzz_xyz[k] * ab_y + g_0_z_zzzzzz_xyyz[k];

                g_0_z_yzzzzzz_xzz[k] = -g_0_z_zzzzzz_xzz[k] * ab_y + g_0_z_zzzzzz_xyzz[k];

                g_0_z_yzzzzzz_yyy[k] = -g_0_z_zzzzzz_yyy[k] * ab_y + g_0_z_zzzzzz_yyyy[k];

                g_0_z_yzzzzzz_yyz[k] = -g_0_z_zzzzzz_yyz[k] * ab_y + g_0_z_zzzzzz_yyyz[k];

                g_0_z_yzzzzzz_yzz[k] = -g_0_z_zzzzzz_yzz[k] * ab_y + g_0_z_zzzzzz_yyzz[k];

                g_0_z_yzzzzzz_zzz[k] = -g_0_z_zzzzzz_zzz[k] * ab_y + g_0_z_zzzzzz_yzzz[k];
            }

            /// Set up 1070-1080 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzzz_xxx = cbuffer.data(kf_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxy = cbuffer.data(kf_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxz = cbuffer.data(kf_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xyy = cbuffer.data(kf_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xyz = cbuffer.data(kf_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xzz = cbuffer.data(kf_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yyy = cbuffer.data(kf_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yyz = cbuffer.data(kf_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yzz = cbuffer.data(kf_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_zzz = cbuffer.data(kf_geom_01_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzzz_xxx, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxy, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xyy, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_yyy, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_zzz, g_0_z_zzzzzz_zzzz, g_0_z_zzzzzzz_xxx, g_0_z_zzzzzzz_xxy, g_0_z_zzzzzzz_xxz, g_0_z_zzzzzzz_xyy, g_0_z_zzzzzzz_xyz, g_0_z_zzzzzzz_xzz, g_0_z_zzzzzzz_yyy, g_0_z_zzzzzzz_yyz, g_0_z_zzzzzzz_yzz, g_0_z_zzzzzzz_zzz, g_zzzzzz_xxx, g_zzzzzz_xxy, g_zzzzzz_xxz, g_zzzzzz_xyy, g_zzzzzz_xyz, g_zzzzzz_xzz, g_zzzzzz_yyy, g_zzzzzz_yyz, g_zzzzzz_yzz, g_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzzz_xxx[k] = g_zzzzzz_xxx[k] - g_0_z_zzzzzz_xxx[k] * ab_z + g_0_z_zzzzzz_xxxz[k];

                g_0_z_zzzzzzz_xxy[k] = g_zzzzzz_xxy[k] - g_0_z_zzzzzz_xxy[k] * ab_z + g_0_z_zzzzzz_xxyz[k];

                g_0_z_zzzzzzz_xxz[k] = g_zzzzzz_xxz[k] - g_0_z_zzzzzz_xxz[k] * ab_z + g_0_z_zzzzzz_xxzz[k];

                g_0_z_zzzzzzz_xyy[k] = g_zzzzzz_xyy[k] - g_0_z_zzzzzz_xyy[k] * ab_z + g_0_z_zzzzzz_xyyz[k];

                g_0_z_zzzzzzz_xyz[k] = g_zzzzzz_xyz[k] - g_0_z_zzzzzz_xyz[k] * ab_z + g_0_z_zzzzzz_xyzz[k];

                g_0_z_zzzzzzz_xzz[k] = g_zzzzzz_xzz[k] - g_0_z_zzzzzz_xzz[k] * ab_z + g_0_z_zzzzzz_xzzz[k];

                g_0_z_zzzzzzz_yyy[k] = g_zzzzzz_yyy[k] - g_0_z_zzzzzz_yyy[k] * ab_z + g_0_z_zzzzzz_yyyz[k];

                g_0_z_zzzzzzz_yyz[k] = g_zzzzzz_yyz[k] - g_0_z_zzzzzz_yyz[k] * ab_z + g_0_z_zzzzzz_yyzz[k];

                g_0_z_zzzzzzz_yzz[k] = g_zzzzzz_yzz[k] - g_0_z_zzzzzz_yzz[k] * ab_z + g_0_z_zzzzzz_yzzz[k];

                g_0_z_zzzzzzz_zzz[k] = g_zzzzzz_zzz[k] - g_0_z_zzzzzz_zzz[k] * ab_z + g_0_z_zzzzzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

