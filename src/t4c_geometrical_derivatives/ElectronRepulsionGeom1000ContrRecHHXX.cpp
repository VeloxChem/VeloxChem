#include "ElectronRepulsionGeom1000ContrRecHHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_hhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_hhxx,
                                            const size_t idx_ghxx,
                                            const size_t idx_geom_10_ghxx,
                                            const size_t idx_geom_10_gixx,
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

            const auto gh_geom_10_off = idx_geom_10_ghxx + i * dcomps + j;

            auto g_x_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 659 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 674 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 675 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 676 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 677 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 678 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 679 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 680 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 681 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 682 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 683 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 684 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 685 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 686 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 687 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 688 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 689 * ccomps * dcomps);

            auto g_z_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 690 * ccomps * dcomps);

            auto g_z_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 691 * ccomps * dcomps);

            auto g_z_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 692 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 693 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 694 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 695 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 696 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 697 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 698 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 699 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 700 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 701 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 702 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 703 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 704 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 705 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 706 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 707 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 708 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 709 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 710 * ccomps * dcomps);

            auto g_z_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 711 * ccomps * dcomps);

            auto g_z_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 712 * ccomps * dcomps);

            auto g_z_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 713 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 714 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 715 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 716 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 717 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 718 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 719 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 727 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 729 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 739 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 749 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 755 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 759 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 769 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 779 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 783 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 789 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 799 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 809 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 811 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 819 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 829 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 839 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 854 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 867 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 869 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 944 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GISS

            const auto gi_geom_10_off = idx_geom_10_gixx + i * dcomps + j;

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

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 54 * ccomps * dcomps);

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

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 138 * ccomps * dcomps);

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

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 250 * ccomps * dcomps);

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

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 390 * ccomps * dcomps);

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

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 755 * ccomps * dcomps);

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

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 783 * ccomps * dcomps);

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

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 811 * ccomps * dcomps);

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

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 839 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_hhxx

            const auto hh_geom_10_off = idx_geom_10_hhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzzz, g_xxxx_xxxxx, g_xxxx_xxxxy, g_xxxx_xxxxz, g_xxxx_xxxyy, g_xxxx_xxxyz, g_xxxx_xxxzz, g_xxxx_xxyyy, g_xxxx_xxyyz, g_xxxx_xxyzz, g_xxxx_xxzzz, g_xxxx_xyyyy, g_xxxx_xyyyz, g_xxxx_xyyzz, g_xxxx_xyzzz, g_xxxx_xzzzz, g_xxxx_yyyyy, g_xxxx_yyyyz, g_xxxx_yyyzz, g_xxxx_yyzzz, g_xxxx_yzzzz, g_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxxx[k] = -g_xxxx_xxxxx[k] - g_x_0_xxxx_xxxxx[k] * ab_x + g_x_0_xxxx_xxxxxx[k];

                g_x_0_xxxxx_xxxxy[k] = -g_xxxx_xxxxy[k] - g_x_0_xxxx_xxxxy[k] * ab_x + g_x_0_xxxx_xxxxxy[k];

                g_x_0_xxxxx_xxxxz[k] = -g_xxxx_xxxxz[k] - g_x_0_xxxx_xxxxz[k] * ab_x + g_x_0_xxxx_xxxxxz[k];

                g_x_0_xxxxx_xxxyy[k] = -g_xxxx_xxxyy[k] - g_x_0_xxxx_xxxyy[k] * ab_x + g_x_0_xxxx_xxxxyy[k];

                g_x_0_xxxxx_xxxyz[k] = -g_xxxx_xxxyz[k] - g_x_0_xxxx_xxxyz[k] * ab_x + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxx_xxxzz[k] = -g_xxxx_xxxzz[k] - g_x_0_xxxx_xxxzz[k] * ab_x + g_x_0_xxxx_xxxxzz[k];

                g_x_0_xxxxx_xxyyy[k] = -g_xxxx_xxyyy[k] - g_x_0_xxxx_xxyyy[k] * ab_x + g_x_0_xxxx_xxxyyy[k];

                g_x_0_xxxxx_xxyyz[k] = -g_xxxx_xxyyz[k] - g_x_0_xxxx_xxyyz[k] * ab_x + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxx_xxyzz[k] = -g_xxxx_xxyzz[k] - g_x_0_xxxx_xxyzz[k] * ab_x + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxx_xxzzz[k] = -g_xxxx_xxzzz[k] - g_x_0_xxxx_xxzzz[k] * ab_x + g_x_0_xxxx_xxxzzz[k];

                g_x_0_xxxxx_xyyyy[k] = -g_xxxx_xyyyy[k] - g_x_0_xxxx_xyyyy[k] * ab_x + g_x_0_xxxx_xxyyyy[k];

                g_x_0_xxxxx_xyyyz[k] = -g_xxxx_xyyyz[k] - g_x_0_xxxx_xyyyz[k] * ab_x + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxx_xyyzz[k] = -g_xxxx_xyyzz[k] - g_x_0_xxxx_xyyzz[k] * ab_x + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxx_xyzzz[k] = -g_xxxx_xyzzz[k] - g_x_0_xxxx_xyzzz[k] * ab_x + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxx_xzzzz[k] = -g_xxxx_xzzzz[k] - g_x_0_xxxx_xzzzz[k] * ab_x + g_x_0_xxxx_xxzzzz[k];

                g_x_0_xxxxx_yyyyy[k] = -g_xxxx_yyyyy[k] - g_x_0_xxxx_yyyyy[k] * ab_x + g_x_0_xxxx_xyyyyy[k];

                g_x_0_xxxxx_yyyyz[k] = -g_xxxx_yyyyz[k] - g_x_0_xxxx_yyyyz[k] * ab_x + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxx_yyyzz[k] = -g_xxxx_yyyzz[k] - g_x_0_xxxx_yyyzz[k] * ab_x + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxx_yyzzz[k] = -g_xxxx_yyzzz[k] - g_x_0_xxxx_yyzzz[k] * ab_x + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxx_yzzzz[k] = -g_xxxx_yzzzz[k] - g_x_0_xxxx_yzzzz[k] * ab_x + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxx_zzzzz[k] = -g_xxxx_zzzzz[k] - g_x_0_xxxx_zzzzz[k] * ab_x + g_x_0_xxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxy_xxxxx, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxxz, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxzz, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxzzz, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xzzzz, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxxx[k] = -g_x_0_xxxx_xxxxx[k] * ab_y + g_x_0_xxxx_xxxxxy[k];

                g_x_0_xxxxy_xxxxy[k] = -g_x_0_xxxx_xxxxy[k] * ab_y + g_x_0_xxxx_xxxxyy[k];

                g_x_0_xxxxy_xxxxz[k] = -g_x_0_xxxx_xxxxz[k] * ab_y + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxy_xxxyy[k] = -g_x_0_xxxx_xxxyy[k] * ab_y + g_x_0_xxxx_xxxyyy[k];

                g_x_0_xxxxy_xxxyz[k] = -g_x_0_xxxx_xxxyz[k] * ab_y + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxy_xxxzz[k] = -g_x_0_xxxx_xxxzz[k] * ab_y + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxy_xxyyy[k] = -g_x_0_xxxx_xxyyy[k] * ab_y + g_x_0_xxxx_xxyyyy[k];

                g_x_0_xxxxy_xxyyz[k] = -g_x_0_xxxx_xxyyz[k] * ab_y + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxy_xxyzz[k] = -g_x_0_xxxx_xxyzz[k] * ab_y + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxy_xxzzz[k] = -g_x_0_xxxx_xxzzz[k] * ab_y + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxy_xyyyy[k] = -g_x_0_xxxx_xyyyy[k] * ab_y + g_x_0_xxxx_xyyyyy[k];

                g_x_0_xxxxy_xyyyz[k] = -g_x_0_xxxx_xyyyz[k] * ab_y + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxy_xyyzz[k] = -g_x_0_xxxx_xyyzz[k] * ab_y + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxy_xyzzz[k] = -g_x_0_xxxx_xyzzz[k] * ab_y + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxy_xzzzz[k] = -g_x_0_xxxx_xzzzz[k] * ab_y + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxy_yyyyy[k] = -g_x_0_xxxx_yyyyy[k] * ab_y + g_x_0_xxxx_yyyyyy[k];

                g_x_0_xxxxy_yyyyz[k] = -g_x_0_xxxx_yyyyz[k] * ab_y + g_x_0_xxxx_yyyyyz[k];

                g_x_0_xxxxy_yyyzz[k] = -g_x_0_xxxx_yyyzz[k] * ab_y + g_x_0_xxxx_yyyyzz[k];

                g_x_0_xxxxy_yyzzz[k] = -g_x_0_xxxx_yyzzz[k] * ab_y + g_x_0_xxxx_yyyzzz[k];

                g_x_0_xxxxy_yzzzz[k] = -g_x_0_xxxx_yzzzz[k] * ab_y + g_x_0_xxxx_yyzzzz[k];

                g_x_0_xxxxy_zzzzz[k] = -g_x_0_xxxx_zzzzz[k] * ab_y + g_x_0_xxxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxxx[k] = -g_x_0_xxxx_xxxxx[k] * ab_z + g_x_0_xxxx_xxxxxz[k];

                g_x_0_xxxxz_xxxxy[k] = -g_x_0_xxxx_xxxxy[k] * ab_z + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxz_xxxxz[k] = -g_x_0_xxxx_xxxxz[k] * ab_z + g_x_0_xxxx_xxxxzz[k];

                g_x_0_xxxxz_xxxyy[k] = -g_x_0_xxxx_xxxyy[k] * ab_z + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxz_xxxyz[k] = -g_x_0_xxxx_xxxyz[k] * ab_z + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxz_xxxzz[k] = -g_x_0_xxxx_xxxzz[k] * ab_z + g_x_0_xxxx_xxxzzz[k];

                g_x_0_xxxxz_xxyyy[k] = -g_x_0_xxxx_xxyyy[k] * ab_z + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxz_xxyyz[k] = -g_x_0_xxxx_xxyyz[k] * ab_z + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxz_xxyzz[k] = -g_x_0_xxxx_xxyzz[k] * ab_z + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxz_xxzzz[k] = -g_x_0_xxxx_xxzzz[k] * ab_z + g_x_0_xxxx_xxzzzz[k];

                g_x_0_xxxxz_xyyyy[k] = -g_x_0_xxxx_xyyyy[k] * ab_z + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxz_xyyyz[k] = -g_x_0_xxxx_xyyyz[k] * ab_z + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxz_xyyzz[k] = -g_x_0_xxxx_xyyzz[k] * ab_z + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxz_xyzzz[k] = -g_x_0_xxxx_xyzzz[k] * ab_z + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxz_xzzzz[k] = -g_x_0_xxxx_xzzzz[k] * ab_z + g_x_0_xxxx_xzzzzz[k];

                g_x_0_xxxxz_yyyyy[k] = -g_x_0_xxxx_yyyyy[k] * ab_z + g_x_0_xxxx_yyyyyz[k];

                g_x_0_xxxxz_yyyyz[k] = -g_x_0_xxxx_yyyyz[k] * ab_z + g_x_0_xxxx_yyyyzz[k];

                g_x_0_xxxxz_yyyzz[k] = -g_x_0_xxxx_yyyzz[k] * ab_z + g_x_0_xxxx_yyyzzz[k];

                g_x_0_xxxxz_yyzzz[k] = -g_x_0_xxxx_yyzzz[k] * ab_z + g_x_0_xxxx_yyzzzz[k];

                g_x_0_xxxxz_yzzzz[k] = -g_x_0_xxxx_yzzzz[k] * ab_z + g_x_0_xxxx_yzzzzz[k];

                g_x_0_xxxxz_zzzzz[k] = -g_x_0_xxxx_zzzzz[k] * ab_z + g_x_0_xxxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxy_xxxxx, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxy, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxz, g_x_0_xxxy_xxxyy, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxzz, g_x_0_xxxy_xxyyy, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxzzz, g_x_0_xxxy_xyyyy, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xzzzz, g_x_0_xxxy_yyyyy, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_zzzzz, g_x_0_xxxyy_xxxxx, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxxz, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxzz, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxzzz, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xzzzz, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxxx[k] = -g_x_0_xxxy_xxxxx[k] * ab_y + g_x_0_xxxy_xxxxxy[k];

                g_x_0_xxxyy_xxxxy[k] = -g_x_0_xxxy_xxxxy[k] * ab_y + g_x_0_xxxy_xxxxyy[k];

                g_x_0_xxxyy_xxxxz[k] = -g_x_0_xxxy_xxxxz[k] * ab_y + g_x_0_xxxy_xxxxyz[k];

                g_x_0_xxxyy_xxxyy[k] = -g_x_0_xxxy_xxxyy[k] * ab_y + g_x_0_xxxy_xxxyyy[k];

                g_x_0_xxxyy_xxxyz[k] = -g_x_0_xxxy_xxxyz[k] * ab_y + g_x_0_xxxy_xxxyyz[k];

                g_x_0_xxxyy_xxxzz[k] = -g_x_0_xxxy_xxxzz[k] * ab_y + g_x_0_xxxy_xxxyzz[k];

                g_x_0_xxxyy_xxyyy[k] = -g_x_0_xxxy_xxyyy[k] * ab_y + g_x_0_xxxy_xxyyyy[k];

                g_x_0_xxxyy_xxyyz[k] = -g_x_0_xxxy_xxyyz[k] * ab_y + g_x_0_xxxy_xxyyyz[k];

                g_x_0_xxxyy_xxyzz[k] = -g_x_0_xxxy_xxyzz[k] * ab_y + g_x_0_xxxy_xxyyzz[k];

                g_x_0_xxxyy_xxzzz[k] = -g_x_0_xxxy_xxzzz[k] * ab_y + g_x_0_xxxy_xxyzzz[k];

                g_x_0_xxxyy_xyyyy[k] = -g_x_0_xxxy_xyyyy[k] * ab_y + g_x_0_xxxy_xyyyyy[k];

                g_x_0_xxxyy_xyyyz[k] = -g_x_0_xxxy_xyyyz[k] * ab_y + g_x_0_xxxy_xyyyyz[k];

                g_x_0_xxxyy_xyyzz[k] = -g_x_0_xxxy_xyyzz[k] * ab_y + g_x_0_xxxy_xyyyzz[k];

                g_x_0_xxxyy_xyzzz[k] = -g_x_0_xxxy_xyzzz[k] * ab_y + g_x_0_xxxy_xyyzzz[k];

                g_x_0_xxxyy_xzzzz[k] = -g_x_0_xxxy_xzzzz[k] * ab_y + g_x_0_xxxy_xyzzzz[k];

                g_x_0_xxxyy_yyyyy[k] = -g_x_0_xxxy_yyyyy[k] * ab_y + g_x_0_xxxy_yyyyyy[k];

                g_x_0_xxxyy_yyyyz[k] = -g_x_0_xxxy_yyyyz[k] * ab_y + g_x_0_xxxy_yyyyyz[k];

                g_x_0_xxxyy_yyyzz[k] = -g_x_0_xxxy_yyyzz[k] * ab_y + g_x_0_xxxy_yyyyzz[k];

                g_x_0_xxxyy_yyzzz[k] = -g_x_0_xxxy_yyzzz[k] * ab_y + g_x_0_xxxy_yyyzzz[k];

                g_x_0_xxxyy_yzzzz[k] = -g_x_0_xxxy_yzzzz[k] * ab_y + g_x_0_xxxy_yyzzzz[k];

                g_x_0_xxxyy_zzzzz[k] = -g_x_0_xxxy_zzzzz[k] * ab_y + g_x_0_xxxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxyz_xxxxx, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxxz, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxzz, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxzzz, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xzzzz, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_zzzzz, g_x_0_xxxz_xxxxx, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxxx[k] = -g_x_0_xxxz_xxxxx[k] * ab_y + g_x_0_xxxz_xxxxxy[k];

                g_x_0_xxxyz_xxxxy[k] = -g_x_0_xxxz_xxxxy[k] * ab_y + g_x_0_xxxz_xxxxyy[k];

                g_x_0_xxxyz_xxxxz[k] = -g_x_0_xxxz_xxxxz[k] * ab_y + g_x_0_xxxz_xxxxyz[k];

                g_x_0_xxxyz_xxxyy[k] = -g_x_0_xxxz_xxxyy[k] * ab_y + g_x_0_xxxz_xxxyyy[k];

                g_x_0_xxxyz_xxxyz[k] = -g_x_0_xxxz_xxxyz[k] * ab_y + g_x_0_xxxz_xxxyyz[k];

                g_x_0_xxxyz_xxxzz[k] = -g_x_0_xxxz_xxxzz[k] * ab_y + g_x_0_xxxz_xxxyzz[k];

                g_x_0_xxxyz_xxyyy[k] = -g_x_0_xxxz_xxyyy[k] * ab_y + g_x_0_xxxz_xxyyyy[k];

                g_x_0_xxxyz_xxyyz[k] = -g_x_0_xxxz_xxyyz[k] * ab_y + g_x_0_xxxz_xxyyyz[k];

                g_x_0_xxxyz_xxyzz[k] = -g_x_0_xxxz_xxyzz[k] * ab_y + g_x_0_xxxz_xxyyzz[k];

                g_x_0_xxxyz_xxzzz[k] = -g_x_0_xxxz_xxzzz[k] * ab_y + g_x_0_xxxz_xxyzzz[k];

                g_x_0_xxxyz_xyyyy[k] = -g_x_0_xxxz_xyyyy[k] * ab_y + g_x_0_xxxz_xyyyyy[k];

                g_x_0_xxxyz_xyyyz[k] = -g_x_0_xxxz_xyyyz[k] * ab_y + g_x_0_xxxz_xyyyyz[k];

                g_x_0_xxxyz_xyyzz[k] = -g_x_0_xxxz_xyyzz[k] * ab_y + g_x_0_xxxz_xyyyzz[k];

                g_x_0_xxxyz_xyzzz[k] = -g_x_0_xxxz_xyzzz[k] * ab_y + g_x_0_xxxz_xyyzzz[k];

                g_x_0_xxxyz_xzzzz[k] = -g_x_0_xxxz_xzzzz[k] * ab_y + g_x_0_xxxz_xyzzzz[k];

                g_x_0_xxxyz_yyyyy[k] = -g_x_0_xxxz_yyyyy[k] * ab_y + g_x_0_xxxz_yyyyyy[k];

                g_x_0_xxxyz_yyyyz[k] = -g_x_0_xxxz_yyyyz[k] * ab_y + g_x_0_xxxz_yyyyyz[k];

                g_x_0_xxxyz_yyyzz[k] = -g_x_0_xxxz_yyyzz[k] * ab_y + g_x_0_xxxz_yyyyzz[k];

                g_x_0_xxxyz_yyzzz[k] = -g_x_0_xxxz_yyzzz[k] * ab_y + g_x_0_xxxz_yyyzzz[k];

                g_x_0_xxxyz_yzzzz[k] = -g_x_0_xxxz_yzzzz[k] * ab_y + g_x_0_xxxz_yyzzzz[k];

                g_x_0_xxxyz_zzzzz[k] = -g_x_0_xxxz_zzzzz[k] * ab_y + g_x_0_xxxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxz_xxxxx, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzz, g_x_0_xxxz_zzzzzz, g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxxx[k] = -g_x_0_xxxz_xxxxx[k] * ab_z + g_x_0_xxxz_xxxxxz[k];

                g_x_0_xxxzz_xxxxy[k] = -g_x_0_xxxz_xxxxy[k] * ab_z + g_x_0_xxxz_xxxxyz[k];

                g_x_0_xxxzz_xxxxz[k] = -g_x_0_xxxz_xxxxz[k] * ab_z + g_x_0_xxxz_xxxxzz[k];

                g_x_0_xxxzz_xxxyy[k] = -g_x_0_xxxz_xxxyy[k] * ab_z + g_x_0_xxxz_xxxyyz[k];

                g_x_0_xxxzz_xxxyz[k] = -g_x_0_xxxz_xxxyz[k] * ab_z + g_x_0_xxxz_xxxyzz[k];

                g_x_0_xxxzz_xxxzz[k] = -g_x_0_xxxz_xxxzz[k] * ab_z + g_x_0_xxxz_xxxzzz[k];

                g_x_0_xxxzz_xxyyy[k] = -g_x_0_xxxz_xxyyy[k] * ab_z + g_x_0_xxxz_xxyyyz[k];

                g_x_0_xxxzz_xxyyz[k] = -g_x_0_xxxz_xxyyz[k] * ab_z + g_x_0_xxxz_xxyyzz[k];

                g_x_0_xxxzz_xxyzz[k] = -g_x_0_xxxz_xxyzz[k] * ab_z + g_x_0_xxxz_xxyzzz[k];

                g_x_0_xxxzz_xxzzz[k] = -g_x_0_xxxz_xxzzz[k] * ab_z + g_x_0_xxxz_xxzzzz[k];

                g_x_0_xxxzz_xyyyy[k] = -g_x_0_xxxz_xyyyy[k] * ab_z + g_x_0_xxxz_xyyyyz[k];

                g_x_0_xxxzz_xyyyz[k] = -g_x_0_xxxz_xyyyz[k] * ab_z + g_x_0_xxxz_xyyyzz[k];

                g_x_0_xxxzz_xyyzz[k] = -g_x_0_xxxz_xyyzz[k] * ab_z + g_x_0_xxxz_xyyzzz[k];

                g_x_0_xxxzz_xyzzz[k] = -g_x_0_xxxz_xyzzz[k] * ab_z + g_x_0_xxxz_xyzzzz[k];

                g_x_0_xxxzz_xzzzz[k] = -g_x_0_xxxz_xzzzz[k] * ab_z + g_x_0_xxxz_xzzzzz[k];

                g_x_0_xxxzz_yyyyy[k] = -g_x_0_xxxz_yyyyy[k] * ab_z + g_x_0_xxxz_yyyyyz[k];

                g_x_0_xxxzz_yyyyz[k] = -g_x_0_xxxz_yyyyz[k] * ab_z + g_x_0_xxxz_yyyyzz[k];

                g_x_0_xxxzz_yyyzz[k] = -g_x_0_xxxz_yyyzz[k] * ab_z + g_x_0_xxxz_yyyzzz[k];

                g_x_0_xxxzz_yyzzz[k] = -g_x_0_xxxz_yyzzz[k] * ab_z + g_x_0_xxxz_yyzzzz[k];

                g_x_0_xxxzz_yzzzz[k] = -g_x_0_xxxz_yzzzz[k] * ab_z + g_x_0_xxxz_yzzzzz[k];

                g_x_0_xxxzz_zzzzz[k] = -g_x_0_xxxz_zzzzz[k] * ab_z + g_x_0_xxxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyy_xxxxx, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxy, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxz, g_x_0_xxyy_xxxyy, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxzz, g_x_0_xxyy_xxyyy, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxzzz, g_x_0_xxyy_xyyyy, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xzzzz, g_x_0_xxyy_yyyyy, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_zzzzz, g_x_0_xxyyy_xxxxx, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxxz, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxzz, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxzzz, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xzzzz, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxxx[k] = -g_x_0_xxyy_xxxxx[k] * ab_y + g_x_0_xxyy_xxxxxy[k];

                g_x_0_xxyyy_xxxxy[k] = -g_x_0_xxyy_xxxxy[k] * ab_y + g_x_0_xxyy_xxxxyy[k];

                g_x_0_xxyyy_xxxxz[k] = -g_x_0_xxyy_xxxxz[k] * ab_y + g_x_0_xxyy_xxxxyz[k];

                g_x_0_xxyyy_xxxyy[k] = -g_x_0_xxyy_xxxyy[k] * ab_y + g_x_0_xxyy_xxxyyy[k];

                g_x_0_xxyyy_xxxyz[k] = -g_x_0_xxyy_xxxyz[k] * ab_y + g_x_0_xxyy_xxxyyz[k];

                g_x_0_xxyyy_xxxzz[k] = -g_x_0_xxyy_xxxzz[k] * ab_y + g_x_0_xxyy_xxxyzz[k];

                g_x_0_xxyyy_xxyyy[k] = -g_x_0_xxyy_xxyyy[k] * ab_y + g_x_0_xxyy_xxyyyy[k];

                g_x_0_xxyyy_xxyyz[k] = -g_x_0_xxyy_xxyyz[k] * ab_y + g_x_0_xxyy_xxyyyz[k];

                g_x_0_xxyyy_xxyzz[k] = -g_x_0_xxyy_xxyzz[k] * ab_y + g_x_0_xxyy_xxyyzz[k];

                g_x_0_xxyyy_xxzzz[k] = -g_x_0_xxyy_xxzzz[k] * ab_y + g_x_0_xxyy_xxyzzz[k];

                g_x_0_xxyyy_xyyyy[k] = -g_x_0_xxyy_xyyyy[k] * ab_y + g_x_0_xxyy_xyyyyy[k];

                g_x_0_xxyyy_xyyyz[k] = -g_x_0_xxyy_xyyyz[k] * ab_y + g_x_0_xxyy_xyyyyz[k];

                g_x_0_xxyyy_xyyzz[k] = -g_x_0_xxyy_xyyzz[k] * ab_y + g_x_0_xxyy_xyyyzz[k];

                g_x_0_xxyyy_xyzzz[k] = -g_x_0_xxyy_xyzzz[k] * ab_y + g_x_0_xxyy_xyyzzz[k];

                g_x_0_xxyyy_xzzzz[k] = -g_x_0_xxyy_xzzzz[k] * ab_y + g_x_0_xxyy_xyzzzz[k];

                g_x_0_xxyyy_yyyyy[k] = -g_x_0_xxyy_yyyyy[k] * ab_y + g_x_0_xxyy_yyyyyy[k];

                g_x_0_xxyyy_yyyyz[k] = -g_x_0_xxyy_yyyyz[k] * ab_y + g_x_0_xxyy_yyyyyz[k];

                g_x_0_xxyyy_yyyzz[k] = -g_x_0_xxyy_yyyzz[k] * ab_y + g_x_0_xxyy_yyyyzz[k];

                g_x_0_xxyyy_yyzzz[k] = -g_x_0_xxyy_yyzzz[k] * ab_y + g_x_0_xxyy_yyyzzz[k];

                g_x_0_xxyyy_yzzzz[k] = -g_x_0_xxyy_yzzzz[k] * ab_y + g_x_0_xxyy_yyzzzz[k];

                g_x_0_xxyyy_zzzzz[k] = -g_x_0_xxyy_zzzzz[k] * ab_y + g_x_0_xxyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyyz_xxxxx, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxxz, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxzz, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxzzz, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xzzzz, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_zzzzz, g_x_0_xxyz_xxxxx, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxy, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxz, g_x_0_xxyz_xxxyy, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxzz, g_x_0_xxyz_xxyyy, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxzzz, g_x_0_xxyz_xyyyy, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xzzzz, g_x_0_xxyz_yyyyy, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxxx[k] = -g_x_0_xxyz_xxxxx[k] * ab_y + g_x_0_xxyz_xxxxxy[k];

                g_x_0_xxyyz_xxxxy[k] = -g_x_0_xxyz_xxxxy[k] * ab_y + g_x_0_xxyz_xxxxyy[k];

                g_x_0_xxyyz_xxxxz[k] = -g_x_0_xxyz_xxxxz[k] * ab_y + g_x_0_xxyz_xxxxyz[k];

                g_x_0_xxyyz_xxxyy[k] = -g_x_0_xxyz_xxxyy[k] * ab_y + g_x_0_xxyz_xxxyyy[k];

                g_x_0_xxyyz_xxxyz[k] = -g_x_0_xxyz_xxxyz[k] * ab_y + g_x_0_xxyz_xxxyyz[k];

                g_x_0_xxyyz_xxxzz[k] = -g_x_0_xxyz_xxxzz[k] * ab_y + g_x_0_xxyz_xxxyzz[k];

                g_x_0_xxyyz_xxyyy[k] = -g_x_0_xxyz_xxyyy[k] * ab_y + g_x_0_xxyz_xxyyyy[k];

                g_x_0_xxyyz_xxyyz[k] = -g_x_0_xxyz_xxyyz[k] * ab_y + g_x_0_xxyz_xxyyyz[k];

                g_x_0_xxyyz_xxyzz[k] = -g_x_0_xxyz_xxyzz[k] * ab_y + g_x_0_xxyz_xxyyzz[k];

                g_x_0_xxyyz_xxzzz[k] = -g_x_0_xxyz_xxzzz[k] * ab_y + g_x_0_xxyz_xxyzzz[k];

                g_x_0_xxyyz_xyyyy[k] = -g_x_0_xxyz_xyyyy[k] * ab_y + g_x_0_xxyz_xyyyyy[k];

                g_x_0_xxyyz_xyyyz[k] = -g_x_0_xxyz_xyyyz[k] * ab_y + g_x_0_xxyz_xyyyyz[k];

                g_x_0_xxyyz_xyyzz[k] = -g_x_0_xxyz_xyyzz[k] * ab_y + g_x_0_xxyz_xyyyzz[k];

                g_x_0_xxyyz_xyzzz[k] = -g_x_0_xxyz_xyzzz[k] * ab_y + g_x_0_xxyz_xyyzzz[k];

                g_x_0_xxyyz_xzzzz[k] = -g_x_0_xxyz_xzzzz[k] * ab_y + g_x_0_xxyz_xyzzzz[k];

                g_x_0_xxyyz_yyyyy[k] = -g_x_0_xxyz_yyyyy[k] * ab_y + g_x_0_xxyz_yyyyyy[k];

                g_x_0_xxyyz_yyyyz[k] = -g_x_0_xxyz_yyyyz[k] * ab_y + g_x_0_xxyz_yyyyyz[k];

                g_x_0_xxyyz_yyyzz[k] = -g_x_0_xxyz_yyyzz[k] * ab_y + g_x_0_xxyz_yyyyzz[k];

                g_x_0_xxyyz_yyzzz[k] = -g_x_0_xxyz_yyzzz[k] * ab_y + g_x_0_xxyz_yyyzzz[k];

                g_x_0_xxyyz_yzzzz[k] = -g_x_0_xxyz_yzzzz[k] * ab_y + g_x_0_xxyz_yyzzzz[k];

                g_x_0_xxyyz_zzzzz[k] = -g_x_0_xxyz_zzzzz[k] * ab_y + g_x_0_xxyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyzz_xxxxx, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxxz, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxzz, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxzzz, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xzzzz, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_zzzzz, g_x_0_xxzz_xxxxx, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxxx[k] = -g_x_0_xxzz_xxxxx[k] * ab_y + g_x_0_xxzz_xxxxxy[k];

                g_x_0_xxyzz_xxxxy[k] = -g_x_0_xxzz_xxxxy[k] * ab_y + g_x_0_xxzz_xxxxyy[k];

                g_x_0_xxyzz_xxxxz[k] = -g_x_0_xxzz_xxxxz[k] * ab_y + g_x_0_xxzz_xxxxyz[k];

                g_x_0_xxyzz_xxxyy[k] = -g_x_0_xxzz_xxxyy[k] * ab_y + g_x_0_xxzz_xxxyyy[k];

                g_x_0_xxyzz_xxxyz[k] = -g_x_0_xxzz_xxxyz[k] * ab_y + g_x_0_xxzz_xxxyyz[k];

                g_x_0_xxyzz_xxxzz[k] = -g_x_0_xxzz_xxxzz[k] * ab_y + g_x_0_xxzz_xxxyzz[k];

                g_x_0_xxyzz_xxyyy[k] = -g_x_0_xxzz_xxyyy[k] * ab_y + g_x_0_xxzz_xxyyyy[k];

                g_x_0_xxyzz_xxyyz[k] = -g_x_0_xxzz_xxyyz[k] * ab_y + g_x_0_xxzz_xxyyyz[k];

                g_x_0_xxyzz_xxyzz[k] = -g_x_0_xxzz_xxyzz[k] * ab_y + g_x_0_xxzz_xxyyzz[k];

                g_x_0_xxyzz_xxzzz[k] = -g_x_0_xxzz_xxzzz[k] * ab_y + g_x_0_xxzz_xxyzzz[k];

                g_x_0_xxyzz_xyyyy[k] = -g_x_0_xxzz_xyyyy[k] * ab_y + g_x_0_xxzz_xyyyyy[k];

                g_x_0_xxyzz_xyyyz[k] = -g_x_0_xxzz_xyyyz[k] * ab_y + g_x_0_xxzz_xyyyyz[k];

                g_x_0_xxyzz_xyyzz[k] = -g_x_0_xxzz_xyyzz[k] * ab_y + g_x_0_xxzz_xyyyzz[k];

                g_x_0_xxyzz_xyzzz[k] = -g_x_0_xxzz_xyzzz[k] * ab_y + g_x_0_xxzz_xyyzzz[k];

                g_x_0_xxyzz_xzzzz[k] = -g_x_0_xxzz_xzzzz[k] * ab_y + g_x_0_xxzz_xyzzzz[k];

                g_x_0_xxyzz_yyyyy[k] = -g_x_0_xxzz_yyyyy[k] * ab_y + g_x_0_xxzz_yyyyyy[k];

                g_x_0_xxyzz_yyyyz[k] = -g_x_0_xxzz_yyyyz[k] * ab_y + g_x_0_xxzz_yyyyyz[k];

                g_x_0_xxyzz_yyyzz[k] = -g_x_0_xxzz_yyyzz[k] * ab_y + g_x_0_xxzz_yyyyzz[k];

                g_x_0_xxyzz_yyzzz[k] = -g_x_0_xxzz_yyzzz[k] * ab_y + g_x_0_xxzz_yyyzzz[k];

                g_x_0_xxyzz_yzzzz[k] = -g_x_0_xxzz_yzzzz[k] * ab_y + g_x_0_xxzz_yyzzzz[k];

                g_x_0_xxyzz_zzzzz[k] = -g_x_0_xxzz_zzzzz[k] * ab_y + g_x_0_xxzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxzz_xxxxx, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzz, g_x_0_xxzz_zzzzzz, g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxxx[k] = -g_x_0_xxzz_xxxxx[k] * ab_z + g_x_0_xxzz_xxxxxz[k];

                g_x_0_xxzzz_xxxxy[k] = -g_x_0_xxzz_xxxxy[k] * ab_z + g_x_0_xxzz_xxxxyz[k];

                g_x_0_xxzzz_xxxxz[k] = -g_x_0_xxzz_xxxxz[k] * ab_z + g_x_0_xxzz_xxxxzz[k];

                g_x_0_xxzzz_xxxyy[k] = -g_x_0_xxzz_xxxyy[k] * ab_z + g_x_0_xxzz_xxxyyz[k];

                g_x_0_xxzzz_xxxyz[k] = -g_x_0_xxzz_xxxyz[k] * ab_z + g_x_0_xxzz_xxxyzz[k];

                g_x_0_xxzzz_xxxzz[k] = -g_x_0_xxzz_xxxzz[k] * ab_z + g_x_0_xxzz_xxxzzz[k];

                g_x_0_xxzzz_xxyyy[k] = -g_x_0_xxzz_xxyyy[k] * ab_z + g_x_0_xxzz_xxyyyz[k];

                g_x_0_xxzzz_xxyyz[k] = -g_x_0_xxzz_xxyyz[k] * ab_z + g_x_0_xxzz_xxyyzz[k];

                g_x_0_xxzzz_xxyzz[k] = -g_x_0_xxzz_xxyzz[k] * ab_z + g_x_0_xxzz_xxyzzz[k];

                g_x_0_xxzzz_xxzzz[k] = -g_x_0_xxzz_xxzzz[k] * ab_z + g_x_0_xxzz_xxzzzz[k];

                g_x_0_xxzzz_xyyyy[k] = -g_x_0_xxzz_xyyyy[k] * ab_z + g_x_0_xxzz_xyyyyz[k];

                g_x_0_xxzzz_xyyyz[k] = -g_x_0_xxzz_xyyyz[k] * ab_z + g_x_0_xxzz_xyyyzz[k];

                g_x_0_xxzzz_xyyzz[k] = -g_x_0_xxzz_xyyzz[k] * ab_z + g_x_0_xxzz_xyyzzz[k];

                g_x_0_xxzzz_xyzzz[k] = -g_x_0_xxzz_xyzzz[k] * ab_z + g_x_0_xxzz_xyzzzz[k];

                g_x_0_xxzzz_xzzzz[k] = -g_x_0_xxzz_xzzzz[k] * ab_z + g_x_0_xxzz_xzzzzz[k];

                g_x_0_xxzzz_yyyyy[k] = -g_x_0_xxzz_yyyyy[k] * ab_z + g_x_0_xxzz_yyyyyz[k];

                g_x_0_xxzzz_yyyyz[k] = -g_x_0_xxzz_yyyyz[k] * ab_z + g_x_0_xxzz_yyyyzz[k];

                g_x_0_xxzzz_yyyzz[k] = -g_x_0_xxzz_yyyzz[k] * ab_z + g_x_0_xxzz_yyyzzz[k];

                g_x_0_xxzzz_yyzzz[k] = -g_x_0_xxzz_yyzzz[k] * ab_z + g_x_0_xxzz_yyzzzz[k];

                g_x_0_xxzzz_yzzzz[k] = -g_x_0_xxzz_yzzzz[k] * ab_z + g_x_0_xxzz_yzzzzz[k];

                g_x_0_xxzzz_zzzzz[k] = -g_x_0_xxzz_zzzzz[k] * ab_z + g_x_0_xxzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyy_xxxxx, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxy, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxz, g_x_0_xyyy_xxxyy, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxzz, g_x_0_xyyy_xxyyy, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxzzz, g_x_0_xyyy_xyyyy, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xzzzz, g_x_0_xyyy_yyyyy, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_zzzzz, g_x_0_xyyyy_xxxxx, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxxz, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxzz, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxzzz, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xzzzz, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxxx[k] = -g_x_0_xyyy_xxxxx[k] * ab_y + g_x_0_xyyy_xxxxxy[k];

                g_x_0_xyyyy_xxxxy[k] = -g_x_0_xyyy_xxxxy[k] * ab_y + g_x_0_xyyy_xxxxyy[k];

                g_x_0_xyyyy_xxxxz[k] = -g_x_0_xyyy_xxxxz[k] * ab_y + g_x_0_xyyy_xxxxyz[k];

                g_x_0_xyyyy_xxxyy[k] = -g_x_0_xyyy_xxxyy[k] * ab_y + g_x_0_xyyy_xxxyyy[k];

                g_x_0_xyyyy_xxxyz[k] = -g_x_0_xyyy_xxxyz[k] * ab_y + g_x_0_xyyy_xxxyyz[k];

                g_x_0_xyyyy_xxxzz[k] = -g_x_0_xyyy_xxxzz[k] * ab_y + g_x_0_xyyy_xxxyzz[k];

                g_x_0_xyyyy_xxyyy[k] = -g_x_0_xyyy_xxyyy[k] * ab_y + g_x_0_xyyy_xxyyyy[k];

                g_x_0_xyyyy_xxyyz[k] = -g_x_0_xyyy_xxyyz[k] * ab_y + g_x_0_xyyy_xxyyyz[k];

                g_x_0_xyyyy_xxyzz[k] = -g_x_0_xyyy_xxyzz[k] * ab_y + g_x_0_xyyy_xxyyzz[k];

                g_x_0_xyyyy_xxzzz[k] = -g_x_0_xyyy_xxzzz[k] * ab_y + g_x_0_xyyy_xxyzzz[k];

                g_x_0_xyyyy_xyyyy[k] = -g_x_0_xyyy_xyyyy[k] * ab_y + g_x_0_xyyy_xyyyyy[k];

                g_x_0_xyyyy_xyyyz[k] = -g_x_0_xyyy_xyyyz[k] * ab_y + g_x_0_xyyy_xyyyyz[k];

                g_x_0_xyyyy_xyyzz[k] = -g_x_0_xyyy_xyyzz[k] * ab_y + g_x_0_xyyy_xyyyzz[k];

                g_x_0_xyyyy_xyzzz[k] = -g_x_0_xyyy_xyzzz[k] * ab_y + g_x_0_xyyy_xyyzzz[k];

                g_x_0_xyyyy_xzzzz[k] = -g_x_0_xyyy_xzzzz[k] * ab_y + g_x_0_xyyy_xyzzzz[k];

                g_x_0_xyyyy_yyyyy[k] = -g_x_0_xyyy_yyyyy[k] * ab_y + g_x_0_xyyy_yyyyyy[k];

                g_x_0_xyyyy_yyyyz[k] = -g_x_0_xyyy_yyyyz[k] * ab_y + g_x_0_xyyy_yyyyyz[k];

                g_x_0_xyyyy_yyyzz[k] = -g_x_0_xyyy_yyyzz[k] * ab_y + g_x_0_xyyy_yyyyzz[k];

                g_x_0_xyyyy_yyzzz[k] = -g_x_0_xyyy_yyzzz[k] * ab_y + g_x_0_xyyy_yyyzzz[k];

                g_x_0_xyyyy_yzzzz[k] = -g_x_0_xyyy_yzzzz[k] * ab_y + g_x_0_xyyy_yyzzzz[k];

                g_x_0_xyyyy_zzzzz[k] = -g_x_0_xyyy_zzzzz[k] * ab_y + g_x_0_xyyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyyz_xxxxx, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxxz, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxzz, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxzzz, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xzzzz, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_zzzzz, g_x_0_xyyz_xxxxx, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxy, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxz, g_x_0_xyyz_xxxyy, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxzz, g_x_0_xyyz_xxyyy, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxzzz, g_x_0_xyyz_xyyyy, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xzzzz, g_x_0_xyyz_yyyyy, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxxx[k] = -g_x_0_xyyz_xxxxx[k] * ab_y + g_x_0_xyyz_xxxxxy[k];

                g_x_0_xyyyz_xxxxy[k] = -g_x_0_xyyz_xxxxy[k] * ab_y + g_x_0_xyyz_xxxxyy[k];

                g_x_0_xyyyz_xxxxz[k] = -g_x_0_xyyz_xxxxz[k] * ab_y + g_x_0_xyyz_xxxxyz[k];

                g_x_0_xyyyz_xxxyy[k] = -g_x_0_xyyz_xxxyy[k] * ab_y + g_x_0_xyyz_xxxyyy[k];

                g_x_0_xyyyz_xxxyz[k] = -g_x_0_xyyz_xxxyz[k] * ab_y + g_x_0_xyyz_xxxyyz[k];

                g_x_0_xyyyz_xxxzz[k] = -g_x_0_xyyz_xxxzz[k] * ab_y + g_x_0_xyyz_xxxyzz[k];

                g_x_0_xyyyz_xxyyy[k] = -g_x_0_xyyz_xxyyy[k] * ab_y + g_x_0_xyyz_xxyyyy[k];

                g_x_0_xyyyz_xxyyz[k] = -g_x_0_xyyz_xxyyz[k] * ab_y + g_x_0_xyyz_xxyyyz[k];

                g_x_0_xyyyz_xxyzz[k] = -g_x_0_xyyz_xxyzz[k] * ab_y + g_x_0_xyyz_xxyyzz[k];

                g_x_0_xyyyz_xxzzz[k] = -g_x_0_xyyz_xxzzz[k] * ab_y + g_x_0_xyyz_xxyzzz[k];

                g_x_0_xyyyz_xyyyy[k] = -g_x_0_xyyz_xyyyy[k] * ab_y + g_x_0_xyyz_xyyyyy[k];

                g_x_0_xyyyz_xyyyz[k] = -g_x_0_xyyz_xyyyz[k] * ab_y + g_x_0_xyyz_xyyyyz[k];

                g_x_0_xyyyz_xyyzz[k] = -g_x_0_xyyz_xyyzz[k] * ab_y + g_x_0_xyyz_xyyyzz[k];

                g_x_0_xyyyz_xyzzz[k] = -g_x_0_xyyz_xyzzz[k] * ab_y + g_x_0_xyyz_xyyzzz[k];

                g_x_0_xyyyz_xzzzz[k] = -g_x_0_xyyz_xzzzz[k] * ab_y + g_x_0_xyyz_xyzzzz[k];

                g_x_0_xyyyz_yyyyy[k] = -g_x_0_xyyz_yyyyy[k] * ab_y + g_x_0_xyyz_yyyyyy[k];

                g_x_0_xyyyz_yyyyz[k] = -g_x_0_xyyz_yyyyz[k] * ab_y + g_x_0_xyyz_yyyyyz[k];

                g_x_0_xyyyz_yyyzz[k] = -g_x_0_xyyz_yyyzz[k] * ab_y + g_x_0_xyyz_yyyyzz[k];

                g_x_0_xyyyz_yyzzz[k] = -g_x_0_xyyz_yyzzz[k] * ab_y + g_x_0_xyyz_yyyzzz[k];

                g_x_0_xyyyz_yzzzz[k] = -g_x_0_xyyz_yzzzz[k] * ab_y + g_x_0_xyyz_yyzzzz[k];

                g_x_0_xyyyz_zzzzz[k] = -g_x_0_xyyz_zzzzz[k] * ab_y + g_x_0_xyyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyzz_xxxxx, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxxz, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxzz, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxzzz, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xzzzz, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_zzzzz, g_x_0_xyzz_xxxxx, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxy, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxz, g_x_0_xyzz_xxxyy, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxzz, g_x_0_xyzz_xxyyy, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxzzz, g_x_0_xyzz_xyyyy, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xzzzz, g_x_0_xyzz_yyyyy, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxxx[k] = -g_x_0_xyzz_xxxxx[k] * ab_y + g_x_0_xyzz_xxxxxy[k];

                g_x_0_xyyzz_xxxxy[k] = -g_x_0_xyzz_xxxxy[k] * ab_y + g_x_0_xyzz_xxxxyy[k];

                g_x_0_xyyzz_xxxxz[k] = -g_x_0_xyzz_xxxxz[k] * ab_y + g_x_0_xyzz_xxxxyz[k];

                g_x_0_xyyzz_xxxyy[k] = -g_x_0_xyzz_xxxyy[k] * ab_y + g_x_0_xyzz_xxxyyy[k];

                g_x_0_xyyzz_xxxyz[k] = -g_x_0_xyzz_xxxyz[k] * ab_y + g_x_0_xyzz_xxxyyz[k];

                g_x_0_xyyzz_xxxzz[k] = -g_x_0_xyzz_xxxzz[k] * ab_y + g_x_0_xyzz_xxxyzz[k];

                g_x_0_xyyzz_xxyyy[k] = -g_x_0_xyzz_xxyyy[k] * ab_y + g_x_0_xyzz_xxyyyy[k];

                g_x_0_xyyzz_xxyyz[k] = -g_x_0_xyzz_xxyyz[k] * ab_y + g_x_0_xyzz_xxyyyz[k];

                g_x_0_xyyzz_xxyzz[k] = -g_x_0_xyzz_xxyzz[k] * ab_y + g_x_0_xyzz_xxyyzz[k];

                g_x_0_xyyzz_xxzzz[k] = -g_x_0_xyzz_xxzzz[k] * ab_y + g_x_0_xyzz_xxyzzz[k];

                g_x_0_xyyzz_xyyyy[k] = -g_x_0_xyzz_xyyyy[k] * ab_y + g_x_0_xyzz_xyyyyy[k];

                g_x_0_xyyzz_xyyyz[k] = -g_x_0_xyzz_xyyyz[k] * ab_y + g_x_0_xyzz_xyyyyz[k];

                g_x_0_xyyzz_xyyzz[k] = -g_x_0_xyzz_xyyzz[k] * ab_y + g_x_0_xyzz_xyyyzz[k];

                g_x_0_xyyzz_xyzzz[k] = -g_x_0_xyzz_xyzzz[k] * ab_y + g_x_0_xyzz_xyyzzz[k];

                g_x_0_xyyzz_xzzzz[k] = -g_x_0_xyzz_xzzzz[k] * ab_y + g_x_0_xyzz_xyzzzz[k];

                g_x_0_xyyzz_yyyyy[k] = -g_x_0_xyzz_yyyyy[k] * ab_y + g_x_0_xyzz_yyyyyy[k];

                g_x_0_xyyzz_yyyyz[k] = -g_x_0_xyzz_yyyyz[k] * ab_y + g_x_0_xyzz_yyyyyz[k];

                g_x_0_xyyzz_yyyzz[k] = -g_x_0_xyzz_yyyzz[k] * ab_y + g_x_0_xyzz_yyyyzz[k];

                g_x_0_xyyzz_yyzzz[k] = -g_x_0_xyzz_yyzzz[k] * ab_y + g_x_0_xyzz_yyyzzz[k];

                g_x_0_xyyzz_yzzzz[k] = -g_x_0_xyzz_yzzzz[k] * ab_y + g_x_0_xyzz_yyzzzz[k];

                g_x_0_xyyzz_zzzzz[k] = -g_x_0_xyzz_zzzzz[k] * ab_y + g_x_0_xyzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyzzz_xxxxx, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxxz, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxzz, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxzzz, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xzzzz, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_zzzzz, g_x_0_xzzz_xxxxx, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxxx[k] = -g_x_0_xzzz_xxxxx[k] * ab_y + g_x_0_xzzz_xxxxxy[k];

                g_x_0_xyzzz_xxxxy[k] = -g_x_0_xzzz_xxxxy[k] * ab_y + g_x_0_xzzz_xxxxyy[k];

                g_x_0_xyzzz_xxxxz[k] = -g_x_0_xzzz_xxxxz[k] * ab_y + g_x_0_xzzz_xxxxyz[k];

                g_x_0_xyzzz_xxxyy[k] = -g_x_0_xzzz_xxxyy[k] * ab_y + g_x_0_xzzz_xxxyyy[k];

                g_x_0_xyzzz_xxxyz[k] = -g_x_0_xzzz_xxxyz[k] * ab_y + g_x_0_xzzz_xxxyyz[k];

                g_x_0_xyzzz_xxxzz[k] = -g_x_0_xzzz_xxxzz[k] * ab_y + g_x_0_xzzz_xxxyzz[k];

                g_x_0_xyzzz_xxyyy[k] = -g_x_0_xzzz_xxyyy[k] * ab_y + g_x_0_xzzz_xxyyyy[k];

                g_x_0_xyzzz_xxyyz[k] = -g_x_0_xzzz_xxyyz[k] * ab_y + g_x_0_xzzz_xxyyyz[k];

                g_x_0_xyzzz_xxyzz[k] = -g_x_0_xzzz_xxyzz[k] * ab_y + g_x_0_xzzz_xxyyzz[k];

                g_x_0_xyzzz_xxzzz[k] = -g_x_0_xzzz_xxzzz[k] * ab_y + g_x_0_xzzz_xxyzzz[k];

                g_x_0_xyzzz_xyyyy[k] = -g_x_0_xzzz_xyyyy[k] * ab_y + g_x_0_xzzz_xyyyyy[k];

                g_x_0_xyzzz_xyyyz[k] = -g_x_0_xzzz_xyyyz[k] * ab_y + g_x_0_xzzz_xyyyyz[k];

                g_x_0_xyzzz_xyyzz[k] = -g_x_0_xzzz_xyyzz[k] * ab_y + g_x_0_xzzz_xyyyzz[k];

                g_x_0_xyzzz_xyzzz[k] = -g_x_0_xzzz_xyzzz[k] * ab_y + g_x_0_xzzz_xyyzzz[k];

                g_x_0_xyzzz_xzzzz[k] = -g_x_0_xzzz_xzzzz[k] * ab_y + g_x_0_xzzz_xyzzzz[k];

                g_x_0_xyzzz_yyyyy[k] = -g_x_0_xzzz_yyyyy[k] * ab_y + g_x_0_xzzz_yyyyyy[k];

                g_x_0_xyzzz_yyyyz[k] = -g_x_0_xzzz_yyyyz[k] * ab_y + g_x_0_xzzz_yyyyyz[k];

                g_x_0_xyzzz_yyyzz[k] = -g_x_0_xzzz_yyyzz[k] * ab_y + g_x_0_xzzz_yyyyzz[k];

                g_x_0_xyzzz_yyzzz[k] = -g_x_0_xzzz_yyzzz[k] * ab_y + g_x_0_xzzz_yyyzzz[k];

                g_x_0_xyzzz_yzzzz[k] = -g_x_0_xzzz_yzzzz[k] * ab_y + g_x_0_xzzz_yyzzzz[k];

                g_x_0_xyzzz_zzzzz[k] = -g_x_0_xzzz_zzzzz[k] * ab_y + g_x_0_xzzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xzzz_xxxxx, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzz, g_x_0_xzzz_zzzzzz, g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxxx[k] = -g_x_0_xzzz_xxxxx[k] * ab_z + g_x_0_xzzz_xxxxxz[k];

                g_x_0_xzzzz_xxxxy[k] = -g_x_0_xzzz_xxxxy[k] * ab_z + g_x_0_xzzz_xxxxyz[k];

                g_x_0_xzzzz_xxxxz[k] = -g_x_0_xzzz_xxxxz[k] * ab_z + g_x_0_xzzz_xxxxzz[k];

                g_x_0_xzzzz_xxxyy[k] = -g_x_0_xzzz_xxxyy[k] * ab_z + g_x_0_xzzz_xxxyyz[k];

                g_x_0_xzzzz_xxxyz[k] = -g_x_0_xzzz_xxxyz[k] * ab_z + g_x_0_xzzz_xxxyzz[k];

                g_x_0_xzzzz_xxxzz[k] = -g_x_0_xzzz_xxxzz[k] * ab_z + g_x_0_xzzz_xxxzzz[k];

                g_x_0_xzzzz_xxyyy[k] = -g_x_0_xzzz_xxyyy[k] * ab_z + g_x_0_xzzz_xxyyyz[k];

                g_x_0_xzzzz_xxyyz[k] = -g_x_0_xzzz_xxyyz[k] * ab_z + g_x_0_xzzz_xxyyzz[k];

                g_x_0_xzzzz_xxyzz[k] = -g_x_0_xzzz_xxyzz[k] * ab_z + g_x_0_xzzz_xxyzzz[k];

                g_x_0_xzzzz_xxzzz[k] = -g_x_0_xzzz_xxzzz[k] * ab_z + g_x_0_xzzz_xxzzzz[k];

                g_x_0_xzzzz_xyyyy[k] = -g_x_0_xzzz_xyyyy[k] * ab_z + g_x_0_xzzz_xyyyyz[k];

                g_x_0_xzzzz_xyyyz[k] = -g_x_0_xzzz_xyyyz[k] * ab_z + g_x_0_xzzz_xyyyzz[k];

                g_x_0_xzzzz_xyyzz[k] = -g_x_0_xzzz_xyyzz[k] * ab_z + g_x_0_xzzz_xyyzzz[k];

                g_x_0_xzzzz_xyzzz[k] = -g_x_0_xzzz_xyzzz[k] * ab_z + g_x_0_xzzz_xyzzzz[k];

                g_x_0_xzzzz_xzzzz[k] = -g_x_0_xzzz_xzzzz[k] * ab_z + g_x_0_xzzz_xzzzzz[k];

                g_x_0_xzzzz_yyyyy[k] = -g_x_0_xzzz_yyyyy[k] * ab_z + g_x_0_xzzz_yyyyyz[k];

                g_x_0_xzzzz_yyyyz[k] = -g_x_0_xzzz_yyyyz[k] * ab_z + g_x_0_xzzz_yyyyzz[k];

                g_x_0_xzzzz_yyyzz[k] = -g_x_0_xzzz_yyyzz[k] * ab_z + g_x_0_xzzz_yyyzzz[k];

                g_x_0_xzzzz_yyzzz[k] = -g_x_0_xzzz_yyzzz[k] * ab_z + g_x_0_xzzz_yyzzzz[k];

                g_x_0_xzzzz_yzzzz[k] = -g_x_0_xzzz_yzzzz[k] * ab_z + g_x_0_xzzz_yzzzzz[k];

                g_x_0_xzzzz_zzzzz[k] = -g_x_0_xzzz_zzzzz[k] * ab_z + g_x_0_xzzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyy_xxxxx, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxy, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxz, g_x_0_yyyy_xxxyy, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxzz, g_x_0_yyyy_xxyyy, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxzzz, g_x_0_yyyy_xyyyy, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xzzzz, g_x_0_yyyy_yyyyy, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_zzzzz, g_x_0_yyyyy_xxxxx, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxxz, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxzz, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxzzz, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xzzzz, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxxx[k] = -g_x_0_yyyy_xxxxx[k] * ab_y + g_x_0_yyyy_xxxxxy[k];

                g_x_0_yyyyy_xxxxy[k] = -g_x_0_yyyy_xxxxy[k] * ab_y + g_x_0_yyyy_xxxxyy[k];

                g_x_0_yyyyy_xxxxz[k] = -g_x_0_yyyy_xxxxz[k] * ab_y + g_x_0_yyyy_xxxxyz[k];

                g_x_0_yyyyy_xxxyy[k] = -g_x_0_yyyy_xxxyy[k] * ab_y + g_x_0_yyyy_xxxyyy[k];

                g_x_0_yyyyy_xxxyz[k] = -g_x_0_yyyy_xxxyz[k] * ab_y + g_x_0_yyyy_xxxyyz[k];

                g_x_0_yyyyy_xxxzz[k] = -g_x_0_yyyy_xxxzz[k] * ab_y + g_x_0_yyyy_xxxyzz[k];

                g_x_0_yyyyy_xxyyy[k] = -g_x_0_yyyy_xxyyy[k] * ab_y + g_x_0_yyyy_xxyyyy[k];

                g_x_0_yyyyy_xxyyz[k] = -g_x_0_yyyy_xxyyz[k] * ab_y + g_x_0_yyyy_xxyyyz[k];

                g_x_0_yyyyy_xxyzz[k] = -g_x_0_yyyy_xxyzz[k] * ab_y + g_x_0_yyyy_xxyyzz[k];

                g_x_0_yyyyy_xxzzz[k] = -g_x_0_yyyy_xxzzz[k] * ab_y + g_x_0_yyyy_xxyzzz[k];

                g_x_0_yyyyy_xyyyy[k] = -g_x_0_yyyy_xyyyy[k] * ab_y + g_x_0_yyyy_xyyyyy[k];

                g_x_0_yyyyy_xyyyz[k] = -g_x_0_yyyy_xyyyz[k] * ab_y + g_x_0_yyyy_xyyyyz[k];

                g_x_0_yyyyy_xyyzz[k] = -g_x_0_yyyy_xyyzz[k] * ab_y + g_x_0_yyyy_xyyyzz[k];

                g_x_0_yyyyy_xyzzz[k] = -g_x_0_yyyy_xyzzz[k] * ab_y + g_x_0_yyyy_xyyzzz[k];

                g_x_0_yyyyy_xzzzz[k] = -g_x_0_yyyy_xzzzz[k] * ab_y + g_x_0_yyyy_xyzzzz[k];

                g_x_0_yyyyy_yyyyy[k] = -g_x_0_yyyy_yyyyy[k] * ab_y + g_x_0_yyyy_yyyyyy[k];

                g_x_0_yyyyy_yyyyz[k] = -g_x_0_yyyy_yyyyz[k] * ab_y + g_x_0_yyyy_yyyyyz[k];

                g_x_0_yyyyy_yyyzz[k] = -g_x_0_yyyy_yyyzz[k] * ab_y + g_x_0_yyyy_yyyyzz[k];

                g_x_0_yyyyy_yyzzz[k] = -g_x_0_yyyy_yyzzz[k] * ab_y + g_x_0_yyyy_yyyzzz[k];

                g_x_0_yyyyy_yzzzz[k] = -g_x_0_yyyy_yzzzz[k] * ab_y + g_x_0_yyyy_yyzzzz[k];

                g_x_0_yyyyy_zzzzz[k] = -g_x_0_yyyy_zzzzz[k] * ab_y + g_x_0_yyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyyz_xxxxx, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxxz, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxzz, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxzzz, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xzzzz, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_zzzzz, g_x_0_yyyz_xxxxx, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxy, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxz, g_x_0_yyyz_xxxyy, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxzz, g_x_0_yyyz_xxyyy, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxzzz, g_x_0_yyyz_xyyyy, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xzzzz, g_x_0_yyyz_yyyyy, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxxx[k] = -g_x_0_yyyz_xxxxx[k] * ab_y + g_x_0_yyyz_xxxxxy[k];

                g_x_0_yyyyz_xxxxy[k] = -g_x_0_yyyz_xxxxy[k] * ab_y + g_x_0_yyyz_xxxxyy[k];

                g_x_0_yyyyz_xxxxz[k] = -g_x_0_yyyz_xxxxz[k] * ab_y + g_x_0_yyyz_xxxxyz[k];

                g_x_0_yyyyz_xxxyy[k] = -g_x_0_yyyz_xxxyy[k] * ab_y + g_x_0_yyyz_xxxyyy[k];

                g_x_0_yyyyz_xxxyz[k] = -g_x_0_yyyz_xxxyz[k] * ab_y + g_x_0_yyyz_xxxyyz[k];

                g_x_0_yyyyz_xxxzz[k] = -g_x_0_yyyz_xxxzz[k] * ab_y + g_x_0_yyyz_xxxyzz[k];

                g_x_0_yyyyz_xxyyy[k] = -g_x_0_yyyz_xxyyy[k] * ab_y + g_x_0_yyyz_xxyyyy[k];

                g_x_0_yyyyz_xxyyz[k] = -g_x_0_yyyz_xxyyz[k] * ab_y + g_x_0_yyyz_xxyyyz[k];

                g_x_0_yyyyz_xxyzz[k] = -g_x_0_yyyz_xxyzz[k] * ab_y + g_x_0_yyyz_xxyyzz[k];

                g_x_0_yyyyz_xxzzz[k] = -g_x_0_yyyz_xxzzz[k] * ab_y + g_x_0_yyyz_xxyzzz[k];

                g_x_0_yyyyz_xyyyy[k] = -g_x_0_yyyz_xyyyy[k] * ab_y + g_x_0_yyyz_xyyyyy[k];

                g_x_0_yyyyz_xyyyz[k] = -g_x_0_yyyz_xyyyz[k] * ab_y + g_x_0_yyyz_xyyyyz[k];

                g_x_0_yyyyz_xyyzz[k] = -g_x_0_yyyz_xyyzz[k] * ab_y + g_x_0_yyyz_xyyyzz[k];

                g_x_0_yyyyz_xyzzz[k] = -g_x_0_yyyz_xyzzz[k] * ab_y + g_x_0_yyyz_xyyzzz[k];

                g_x_0_yyyyz_xzzzz[k] = -g_x_0_yyyz_xzzzz[k] * ab_y + g_x_0_yyyz_xyzzzz[k];

                g_x_0_yyyyz_yyyyy[k] = -g_x_0_yyyz_yyyyy[k] * ab_y + g_x_0_yyyz_yyyyyy[k];

                g_x_0_yyyyz_yyyyz[k] = -g_x_0_yyyz_yyyyz[k] * ab_y + g_x_0_yyyz_yyyyyz[k];

                g_x_0_yyyyz_yyyzz[k] = -g_x_0_yyyz_yyyzz[k] * ab_y + g_x_0_yyyz_yyyyzz[k];

                g_x_0_yyyyz_yyzzz[k] = -g_x_0_yyyz_yyzzz[k] * ab_y + g_x_0_yyyz_yyyzzz[k];

                g_x_0_yyyyz_yzzzz[k] = -g_x_0_yyyz_yzzzz[k] * ab_y + g_x_0_yyyz_yyzzzz[k];

                g_x_0_yyyyz_zzzzz[k] = -g_x_0_yyyz_zzzzz[k] * ab_y + g_x_0_yyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyzz_xxxxx, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxxz, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxzz, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxzzz, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xzzzz, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_zzzzz, g_x_0_yyzz_xxxxx, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxy, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxz, g_x_0_yyzz_xxxyy, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxzz, g_x_0_yyzz_xxyyy, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxzzz, g_x_0_yyzz_xyyyy, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xzzzz, g_x_0_yyzz_yyyyy, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxxx[k] = -g_x_0_yyzz_xxxxx[k] * ab_y + g_x_0_yyzz_xxxxxy[k];

                g_x_0_yyyzz_xxxxy[k] = -g_x_0_yyzz_xxxxy[k] * ab_y + g_x_0_yyzz_xxxxyy[k];

                g_x_0_yyyzz_xxxxz[k] = -g_x_0_yyzz_xxxxz[k] * ab_y + g_x_0_yyzz_xxxxyz[k];

                g_x_0_yyyzz_xxxyy[k] = -g_x_0_yyzz_xxxyy[k] * ab_y + g_x_0_yyzz_xxxyyy[k];

                g_x_0_yyyzz_xxxyz[k] = -g_x_0_yyzz_xxxyz[k] * ab_y + g_x_0_yyzz_xxxyyz[k];

                g_x_0_yyyzz_xxxzz[k] = -g_x_0_yyzz_xxxzz[k] * ab_y + g_x_0_yyzz_xxxyzz[k];

                g_x_0_yyyzz_xxyyy[k] = -g_x_0_yyzz_xxyyy[k] * ab_y + g_x_0_yyzz_xxyyyy[k];

                g_x_0_yyyzz_xxyyz[k] = -g_x_0_yyzz_xxyyz[k] * ab_y + g_x_0_yyzz_xxyyyz[k];

                g_x_0_yyyzz_xxyzz[k] = -g_x_0_yyzz_xxyzz[k] * ab_y + g_x_0_yyzz_xxyyzz[k];

                g_x_0_yyyzz_xxzzz[k] = -g_x_0_yyzz_xxzzz[k] * ab_y + g_x_0_yyzz_xxyzzz[k];

                g_x_0_yyyzz_xyyyy[k] = -g_x_0_yyzz_xyyyy[k] * ab_y + g_x_0_yyzz_xyyyyy[k];

                g_x_0_yyyzz_xyyyz[k] = -g_x_0_yyzz_xyyyz[k] * ab_y + g_x_0_yyzz_xyyyyz[k];

                g_x_0_yyyzz_xyyzz[k] = -g_x_0_yyzz_xyyzz[k] * ab_y + g_x_0_yyzz_xyyyzz[k];

                g_x_0_yyyzz_xyzzz[k] = -g_x_0_yyzz_xyzzz[k] * ab_y + g_x_0_yyzz_xyyzzz[k];

                g_x_0_yyyzz_xzzzz[k] = -g_x_0_yyzz_xzzzz[k] * ab_y + g_x_0_yyzz_xyzzzz[k];

                g_x_0_yyyzz_yyyyy[k] = -g_x_0_yyzz_yyyyy[k] * ab_y + g_x_0_yyzz_yyyyyy[k];

                g_x_0_yyyzz_yyyyz[k] = -g_x_0_yyzz_yyyyz[k] * ab_y + g_x_0_yyzz_yyyyyz[k];

                g_x_0_yyyzz_yyyzz[k] = -g_x_0_yyzz_yyyzz[k] * ab_y + g_x_0_yyzz_yyyyzz[k];

                g_x_0_yyyzz_yyzzz[k] = -g_x_0_yyzz_yyzzz[k] * ab_y + g_x_0_yyzz_yyyzzz[k];

                g_x_0_yyyzz_yzzzz[k] = -g_x_0_yyzz_yzzzz[k] * ab_y + g_x_0_yyzz_yyzzzz[k];

                g_x_0_yyyzz_zzzzz[k] = -g_x_0_yyzz_zzzzz[k] * ab_y + g_x_0_yyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyzzz_xxxxx, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxxz, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxzz, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxzzz, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xzzzz, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_zzzzz, g_x_0_yzzz_xxxxx, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxy, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxz, g_x_0_yzzz_xxxyy, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxzz, g_x_0_yzzz_xxyyy, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxzzz, g_x_0_yzzz_xyyyy, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xzzzz, g_x_0_yzzz_yyyyy, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxxx[k] = -g_x_0_yzzz_xxxxx[k] * ab_y + g_x_0_yzzz_xxxxxy[k];

                g_x_0_yyzzz_xxxxy[k] = -g_x_0_yzzz_xxxxy[k] * ab_y + g_x_0_yzzz_xxxxyy[k];

                g_x_0_yyzzz_xxxxz[k] = -g_x_0_yzzz_xxxxz[k] * ab_y + g_x_0_yzzz_xxxxyz[k];

                g_x_0_yyzzz_xxxyy[k] = -g_x_0_yzzz_xxxyy[k] * ab_y + g_x_0_yzzz_xxxyyy[k];

                g_x_0_yyzzz_xxxyz[k] = -g_x_0_yzzz_xxxyz[k] * ab_y + g_x_0_yzzz_xxxyyz[k];

                g_x_0_yyzzz_xxxzz[k] = -g_x_0_yzzz_xxxzz[k] * ab_y + g_x_0_yzzz_xxxyzz[k];

                g_x_0_yyzzz_xxyyy[k] = -g_x_0_yzzz_xxyyy[k] * ab_y + g_x_0_yzzz_xxyyyy[k];

                g_x_0_yyzzz_xxyyz[k] = -g_x_0_yzzz_xxyyz[k] * ab_y + g_x_0_yzzz_xxyyyz[k];

                g_x_0_yyzzz_xxyzz[k] = -g_x_0_yzzz_xxyzz[k] * ab_y + g_x_0_yzzz_xxyyzz[k];

                g_x_0_yyzzz_xxzzz[k] = -g_x_0_yzzz_xxzzz[k] * ab_y + g_x_0_yzzz_xxyzzz[k];

                g_x_0_yyzzz_xyyyy[k] = -g_x_0_yzzz_xyyyy[k] * ab_y + g_x_0_yzzz_xyyyyy[k];

                g_x_0_yyzzz_xyyyz[k] = -g_x_0_yzzz_xyyyz[k] * ab_y + g_x_0_yzzz_xyyyyz[k];

                g_x_0_yyzzz_xyyzz[k] = -g_x_0_yzzz_xyyzz[k] * ab_y + g_x_0_yzzz_xyyyzz[k];

                g_x_0_yyzzz_xyzzz[k] = -g_x_0_yzzz_xyzzz[k] * ab_y + g_x_0_yzzz_xyyzzz[k];

                g_x_0_yyzzz_xzzzz[k] = -g_x_0_yzzz_xzzzz[k] * ab_y + g_x_0_yzzz_xyzzzz[k];

                g_x_0_yyzzz_yyyyy[k] = -g_x_0_yzzz_yyyyy[k] * ab_y + g_x_0_yzzz_yyyyyy[k];

                g_x_0_yyzzz_yyyyz[k] = -g_x_0_yzzz_yyyyz[k] * ab_y + g_x_0_yzzz_yyyyyz[k];

                g_x_0_yyzzz_yyyzz[k] = -g_x_0_yzzz_yyyzz[k] * ab_y + g_x_0_yzzz_yyyyzz[k];

                g_x_0_yyzzz_yyzzz[k] = -g_x_0_yzzz_yyzzz[k] * ab_y + g_x_0_yzzz_yyyzzz[k];

                g_x_0_yyzzz_yzzzz[k] = -g_x_0_yzzz_yzzzz[k] * ab_y + g_x_0_yzzz_yyzzzz[k];

                g_x_0_yyzzz_zzzzz[k] = -g_x_0_yzzz_zzzzz[k] * ab_y + g_x_0_yzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yzzzz_xxxxx, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxxz, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxzz, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxzzz, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xzzzz, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_zzzzz, g_x_0_zzzz_xxxxx, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxxx[k] = -g_x_0_zzzz_xxxxx[k] * ab_y + g_x_0_zzzz_xxxxxy[k];

                g_x_0_yzzzz_xxxxy[k] = -g_x_0_zzzz_xxxxy[k] * ab_y + g_x_0_zzzz_xxxxyy[k];

                g_x_0_yzzzz_xxxxz[k] = -g_x_0_zzzz_xxxxz[k] * ab_y + g_x_0_zzzz_xxxxyz[k];

                g_x_0_yzzzz_xxxyy[k] = -g_x_0_zzzz_xxxyy[k] * ab_y + g_x_0_zzzz_xxxyyy[k];

                g_x_0_yzzzz_xxxyz[k] = -g_x_0_zzzz_xxxyz[k] * ab_y + g_x_0_zzzz_xxxyyz[k];

                g_x_0_yzzzz_xxxzz[k] = -g_x_0_zzzz_xxxzz[k] * ab_y + g_x_0_zzzz_xxxyzz[k];

                g_x_0_yzzzz_xxyyy[k] = -g_x_0_zzzz_xxyyy[k] * ab_y + g_x_0_zzzz_xxyyyy[k];

                g_x_0_yzzzz_xxyyz[k] = -g_x_0_zzzz_xxyyz[k] * ab_y + g_x_0_zzzz_xxyyyz[k];

                g_x_0_yzzzz_xxyzz[k] = -g_x_0_zzzz_xxyzz[k] * ab_y + g_x_0_zzzz_xxyyzz[k];

                g_x_0_yzzzz_xxzzz[k] = -g_x_0_zzzz_xxzzz[k] * ab_y + g_x_0_zzzz_xxyzzz[k];

                g_x_0_yzzzz_xyyyy[k] = -g_x_0_zzzz_xyyyy[k] * ab_y + g_x_0_zzzz_xyyyyy[k];

                g_x_0_yzzzz_xyyyz[k] = -g_x_0_zzzz_xyyyz[k] * ab_y + g_x_0_zzzz_xyyyyz[k];

                g_x_0_yzzzz_xyyzz[k] = -g_x_0_zzzz_xyyzz[k] * ab_y + g_x_0_zzzz_xyyyzz[k];

                g_x_0_yzzzz_xyzzz[k] = -g_x_0_zzzz_xyzzz[k] * ab_y + g_x_0_zzzz_xyyzzz[k];

                g_x_0_yzzzz_xzzzz[k] = -g_x_0_zzzz_xzzzz[k] * ab_y + g_x_0_zzzz_xyzzzz[k];

                g_x_0_yzzzz_yyyyy[k] = -g_x_0_zzzz_yyyyy[k] * ab_y + g_x_0_zzzz_yyyyyy[k];

                g_x_0_yzzzz_yyyyz[k] = -g_x_0_zzzz_yyyyz[k] * ab_y + g_x_0_zzzz_yyyyyz[k];

                g_x_0_yzzzz_yyyzz[k] = -g_x_0_zzzz_yyyzz[k] * ab_y + g_x_0_zzzz_yyyyzz[k];

                g_x_0_yzzzz_yyzzz[k] = -g_x_0_zzzz_yyzzz[k] * ab_y + g_x_0_zzzz_yyyzzz[k];

                g_x_0_yzzzz_yzzzz[k] = -g_x_0_zzzz_yzzzz[k] * ab_y + g_x_0_zzzz_yyzzzz[k];

                g_x_0_yzzzz_zzzzz[k] = -g_x_0_zzzz_zzzzz[k] * ab_y + g_x_0_zzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_zzzz_xxxxx, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzz, g_x_0_zzzz_zzzzzz, g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxxx[k] = -g_x_0_zzzz_xxxxx[k] * ab_z + g_x_0_zzzz_xxxxxz[k];

                g_x_0_zzzzz_xxxxy[k] = -g_x_0_zzzz_xxxxy[k] * ab_z + g_x_0_zzzz_xxxxyz[k];

                g_x_0_zzzzz_xxxxz[k] = -g_x_0_zzzz_xxxxz[k] * ab_z + g_x_0_zzzz_xxxxzz[k];

                g_x_0_zzzzz_xxxyy[k] = -g_x_0_zzzz_xxxyy[k] * ab_z + g_x_0_zzzz_xxxyyz[k];

                g_x_0_zzzzz_xxxyz[k] = -g_x_0_zzzz_xxxyz[k] * ab_z + g_x_0_zzzz_xxxyzz[k];

                g_x_0_zzzzz_xxxzz[k] = -g_x_0_zzzz_xxxzz[k] * ab_z + g_x_0_zzzz_xxxzzz[k];

                g_x_0_zzzzz_xxyyy[k] = -g_x_0_zzzz_xxyyy[k] * ab_z + g_x_0_zzzz_xxyyyz[k];

                g_x_0_zzzzz_xxyyz[k] = -g_x_0_zzzz_xxyyz[k] * ab_z + g_x_0_zzzz_xxyyzz[k];

                g_x_0_zzzzz_xxyzz[k] = -g_x_0_zzzz_xxyzz[k] * ab_z + g_x_0_zzzz_xxyzzz[k];

                g_x_0_zzzzz_xxzzz[k] = -g_x_0_zzzz_xxzzz[k] * ab_z + g_x_0_zzzz_xxzzzz[k];

                g_x_0_zzzzz_xyyyy[k] = -g_x_0_zzzz_xyyyy[k] * ab_z + g_x_0_zzzz_xyyyyz[k];

                g_x_0_zzzzz_xyyyz[k] = -g_x_0_zzzz_xyyyz[k] * ab_z + g_x_0_zzzz_xyyyzz[k];

                g_x_0_zzzzz_xyyzz[k] = -g_x_0_zzzz_xyyzz[k] * ab_z + g_x_0_zzzz_xyyzzz[k];

                g_x_0_zzzzz_xyzzz[k] = -g_x_0_zzzz_xyzzz[k] * ab_z + g_x_0_zzzz_xyzzzz[k];

                g_x_0_zzzzz_xzzzz[k] = -g_x_0_zzzz_xzzzz[k] * ab_z + g_x_0_zzzz_xzzzzz[k];

                g_x_0_zzzzz_yyyyy[k] = -g_x_0_zzzz_yyyyy[k] * ab_z + g_x_0_zzzz_yyyyyz[k];

                g_x_0_zzzzz_yyyyz[k] = -g_x_0_zzzz_yyyyz[k] * ab_z + g_x_0_zzzz_yyyyzz[k];

                g_x_0_zzzzz_yyyzz[k] = -g_x_0_zzzz_yyyzz[k] * ab_z + g_x_0_zzzz_yyyzzz[k];

                g_x_0_zzzzz_yyzzz[k] = -g_x_0_zzzz_yyzzz[k] * ab_z + g_x_0_zzzz_yyzzzz[k];

                g_x_0_zzzzz_yzzzz[k] = -g_x_0_zzzz_yzzzz[k] * ab_z + g_x_0_zzzz_yzzzzz[k];

                g_x_0_zzzzz_zzzzz[k] = -g_x_0_zzzz_zzzzz[k] * ab_z + g_x_0_zzzz_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxx_xxxxx, g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxy, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxyy, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxyyy, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xyyyy, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_yyyyy, g_y_0_xxxx_yyyyz, g_y_0_xxxx_yyyzz, g_y_0_xxxx_yyzzz, g_y_0_xxxx_yzzzz, g_y_0_xxxx_zzzzz, g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_yyyyy, g_y_0_xxxxx_yyyyz, g_y_0_xxxxx_yyyzz, g_y_0_xxxxx_yyzzz, g_y_0_xxxxx_yzzzz, g_y_0_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxxx[k] = -g_y_0_xxxx_xxxxx[k] * ab_x + g_y_0_xxxx_xxxxxx[k];

                g_y_0_xxxxx_xxxxy[k] = -g_y_0_xxxx_xxxxy[k] * ab_x + g_y_0_xxxx_xxxxxy[k];

                g_y_0_xxxxx_xxxxz[k] = -g_y_0_xxxx_xxxxz[k] * ab_x + g_y_0_xxxx_xxxxxz[k];

                g_y_0_xxxxx_xxxyy[k] = -g_y_0_xxxx_xxxyy[k] * ab_x + g_y_0_xxxx_xxxxyy[k];

                g_y_0_xxxxx_xxxyz[k] = -g_y_0_xxxx_xxxyz[k] * ab_x + g_y_0_xxxx_xxxxyz[k];

                g_y_0_xxxxx_xxxzz[k] = -g_y_0_xxxx_xxxzz[k] * ab_x + g_y_0_xxxx_xxxxzz[k];

                g_y_0_xxxxx_xxyyy[k] = -g_y_0_xxxx_xxyyy[k] * ab_x + g_y_0_xxxx_xxxyyy[k];

                g_y_0_xxxxx_xxyyz[k] = -g_y_0_xxxx_xxyyz[k] * ab_x + g_y_0_xxxx_xxxyyz[k];

                g_y_0_xxxxx_xxyzz[k] = -g_y_0_xxxx_xxyzz[k] * ab_x + g_y_0_xxxx_xxxyzz[k];

                g_y_0_xxxxx_xxzzz[k] = -g_y_0_xxxx_xxzzz[k] * ab_x + g_y_0_xxxx_xxxzzz[k];

                g_y_0_xxxxx_xyyyy[k] = -g_y_0_xxxx_xyyyy[k] * ab_x + g_y_0_xxxx_xxyyyy[k];

                g_y_0_xxxxx_xyyyz[k] = -g_y_0_xxxx_xyyyz[k] * ab_x + g_y_0_xxxx_xxyyyz[k];

                g_y_0_xxxxx_xyyzz[k] = -g_y_0_xxxx_xyyzz[k] * ab_x + g_y_0_xxxx_xxyyzz[k];

                g_y_0_xxxxx_xyzzz[k] = -g_y_0_xxxx_xyzzz[k] * ab_x + g_y_0_xxxx_xxyzzz[k];

                g_y_0_xxxxx_xzzzz[k] = -g_y_0_xxxx_xzzzz[k] * ab_x + g_y_0_xxxx_xxzzzz[k];

                g_y_0_xxxxx_yyyyy[k] = -g_y_0_xxxx_yyyyy[k] * ab_x + g_y_0_xxxx_xyyyyy[k];

                g_y_0_xxxxx_yyyyz[k] = -g_y_0_xxxx_yyyyz[k] * ab_x + g_y_0_xxxx_xyyyyz[k];

                g_y_0_xxxxx_yyyzz[k] = -g_y_0_xxxx_yyyzz[k] * ab_x + g_y_0_xxxx_xyyyzz[k];

                g_y_0_xxxxx_yyzzz[k] = -g_y_0_xxxx_yyzzz[k] * ab_x + g_y_0_xxxx_xyyzzz[k];

                g_y_0_xxxxx_yzzzz[k] = -g_y_0_xxxx_yzzzz[k] * ab_x + g_y_0_xxxx_xyzzzz[k];

                g_y_0_xxxxx_zzzzz[k] = -g_y_0_xxxx_zzzzz[k] * ab_x + g_y_0_xxxx_xzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_yyyyy, g_y_0_xxxxy_yyyyz, g_y_0_xxxxy_yyyzz, g_y_0_xxxxy_yyzzz, g_y_0_xxxxy_yzzzz, g_y_0_xxxxy_zzzzz, g_y_0_xxxy_xxxxx, g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxy, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxyy, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxyyy, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xyyyy, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_yyyyy, g_y_0_xxxy_yyyyz, g_y_0_xxxy_yyyzz, g_y_0_xxxy_yyzzz, g_y_0_xxxy_yzzzz, g_y_0_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxxx[k] = -g_y_0_xxxy_xxxxx[k] * ab_x + g_y_0_xxxy_xxxxxx[k];

                g_y_0_xxxxy_xxxxy[k] = -g_y_0_xxxy_xxxxy[k] * ab_x + g_y_0_xxxy_xxxxxy[k];

                g_y_0_xxxxy_xxxxz[k] = -g_y_0_xxxy_xxxxz[k] * ab_x + g_y_0_xxxy_xxxxxz[k];

                g_y_0_xxxxy_xxxyy[k] = -g_y_0_xxxy_xxxyy[k] * ab_x + g_y_0_xxxy_xxxxyy[k];

                g_y_0_xxxxy_xxxyz[k] = -g_y_0_xxxy_xxxyz[k] * ab_x + g_y_0_xxxy_xxxxyz[k];

                g_y_0_xxxxy_xxxzz[k] = -g_y_0_xxxy_xxxzz[k] * ab_x + g_y_0_xxxy_xxxxzz[k];

                g_y_0_xxxxy_xxyyy[k] = -g_y_0_xxxy_xxyyy[k] * ab_x + g_y_0_xxxy_xxxyyy[k];

                g_y_0_xxxxy_xxyyz[k] = -g_y_0_xxxy_xxyyz[k] * ab_x + g_y_0_xxxy_xxxyyz[k];

                g_y_0_xxxxy_xxyzz[k] = -g_y_0_xxxy_xxyzz[k] * ab_x + g_y_0_xxxy_xxxyzz[k];

                g_y_0_xxxxy_xxzzz[k] = -g_y_0_xxxy_xxzzz[k] * ab_x + g_y_0_xxxy_xxxzzz[k];

                g_y_0_xxxxy_xyyyy[k] = -g_y_0_xxxy_xyyyy[k] * ab_x + g_y_0_xxxy_xxyyyy[k];

                g_y_0_xxxxy_xyyyz[k] = -g_y_0_xxxy_xyyyz[k] * ab_x + g_y_0_xxxy_xxyyyz[k];

                g_y_0_xxxxy_xyyzz[k] = -g_y_0_xxxy_xyyzz[k] * ab_x + g_y_0_xxxy_xxyyzz[k];

                g_y_0_xxxxy_xyzzz[k] = -g_y_0_xxxy_xyzzz[k] * ab_x + g_y_0_xxxy_xxyzzz[k];

                g_y_0_xxxxy_xzzzz[k] = -g_y_0_xxxy_xzzzz[k] * ab_x + g_y_0_xxxy_xxzzzz[k];

                g_y_0_xxxxy_yyyyy[k] = -g_y_0_xxxy_yyyyy[k] * ab_x + g_y_0_xxxy_xyyyyy[k];

                g_y_0_xxxxy_yyyyz[k] = -g_y_0_xxxy_yyyyz[k] * ab_x + g_y_0_xxxy_xyyyyz[k];

                g_y_0_xxxxy_yyyzz[k] = -g_y_0_xxxy_yyyzz[k] * ab_x + g_y_0_xxxy_xyyyzz[k];

                g_y_0_xxxxy_yyzzz[k] = -g_y_0_xxxy_yyzzz[k] * ab_x + g_y_0_xxxy_xyyzzz[k];

                g_y_0_xxxxy_yzzzz[k] = -g_y_0_xxxy_yzzzz[k] * ab_x + g_y_0_xxxy_xyzzzz[k];

                g_y_0_xxxxy_zzzzz[k] = -g_y_0_xxxy_zzzzz[k] * ab_x + g_y_0_xxxy_xzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_yyyyy, g_y_0_xxxxz_yyyyz, g_y_0_xxxxz_yyyzz, g_y_0_xxxxz_yyzzz, g_y_0_xxxxz_yzzzz, g_y_0_xxxxz_zzzzz, g_y_0_xxxz_xxxxx, g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxy, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxyy, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxyyy, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xyyyy, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_yyyyy, g_y_0_xxxz_yyyyz, g_y_0_xxxz_yyyzz, g_y_0_xxxz_yyzzz, g_y_0_xxxz_yzzzz, g_y_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxxx[k] = -g_y_0_xxxz_xxxxx[k] * ab_x + g_y_0_xxxz_xxxxxx[k];

                g_y_0_xxxxz_xxxxy[k] = -g_y_0_xxxz_xxxxy[k] * ab_x + g_y_0_xxxz_xxxxxy[k];

                g_y_0_xxxxz_xxxxz[k] = -g_y_0_xxxz_xxxxz[k] * ab_x + g_y_0_xxxz_xxxxxz[k];

                g_y_0_xxxxz_xxxyy[k] = -g_y_0_xxxz_xxxyy[k] * ab_x + g_y_0_xxxz_xxxxyy[k];

                g_y_0_xxxxz_xxxyz[k] = -g_y_0_xxxz_xxxyz[k] * ab_x + g_y_0_xxxz_xxxxyz[k];

                g_y_0_xxxxz_xxxzz[k] = -g_y_0_xxxz_xxxzz[k] * ab_x + g_y_0_xxxz_xxxxzz[k];

                g_y_0_xxxxz_xxyyy[k] = -g_y_0_xxxz_xxyyy[k] * ab_x + g_y_0_xxxz_xxxyyy[k];

                g_y_0_xxxxz_xxyyz[k] = -g_y_0_xxxz_xxyyz[k] * ab_x + g_y_0_xxxz_xxxyyz[k];

                g_y_0_xxxxz_xxyzz[k] = -g_y_0_xxxz_xxyzz[k] * ab_x + g_y_0_xxxz_xxxyzz[k];

                g_y_0_xxxxz_xxzzz[k] = -g_y_0_xxxz_xxzzz[k] * ab_x + g_y_0_xxxz_xxxzzz[k];

                g_y_0_xxxxz_xyyyy[k] = -g_y_0_xxxz_xyyyy[k] * ab_x + g_y_0_xxxz_xxyyyy[k];

                g_y_0_xxxxz_xyyyz[k] = -g_y_0_xxxz_xyyyz[k] * ab_x + g_y_0_xxxz_xxyyyz[k];

                g_y_0_xxxxz_xyyzz[k] = -g_y_0_xxxz_xyyzz[k] * ab_x + g_y_0_xxxz_xxyyzz[k];

                g_y_0_xxxxz_xyzzz[k] = -g_y_0_xxxz_xyzzz[k] * ab_x + g_y_0_xxxz_xxyzzz[k];

                g_y_0_xxxxz_xzzzz[k] = -g_y_0_xxxz_xzzzz[k] * ab_x + g_y_0_xxxz_xxzzzz[k];

                g_y_0_xxxxz_yyyyy[k] = -g_y_0_xxxz_yyyyy[k] * ab_x + g_y_0_xxxz_xyyyyy[k];

                g_y_0_xxxxz_yyyyz[k] = -g_y_0_xxxz_yyyyz[k] * ab_x + g_y_0_xxxz_xyyyyz[k];

                g_y_0_xxxxz_yyyzz[k] = -g_y_0_xxxz_yyyzz[k] * ab_x + g_y_0_xxxz_xyyyzz[k];

                g_y_0_xxxxz_yyzzz[k] = -g_y_0_xxxz_yyzzz[k] * ab_x + g_y_0_xxxz_xyyzzz[k];

                g_y_0_xxxxz_yzzzz[k] = -g_y_0_xxxz_yzzzz[k] * ab_x + g_y_0_xxxz_xyzzzz[k];

                g_y_0_xxxxz_zzzzz[k] = -g_y_0_xxxz_zzzzz[k] * ab_x + g_y_0_xxxz_xzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_yyyyy, g_y_0_xxxyy_yyyyz, g_y_0_xxxyy_yyyzz, g_y_0_xxxyy_yyzzz, g_y_0_xxxyy_yzzzz, g_y_0_xxxyy_zzzzz, g_y_0_xxyy_xxxxx, g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxy, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxyy, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxyyy, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xyyyy, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_yyyyy, g_y_0_xxyy_yyyyz, g_y_0_xxyy_yyyzz, g_y_0_xxyy_yyzzz, g_y_0_xxyy_yzzzz, g_y_0_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxxx[k] = -g_y_0_xxyy_xxxxx[k] * ab_x + g_y_0_xxyy_xxxxxx[k];

                g_y_0_xxxyy_xxxxy[k] = -g_y_0_xxyy_xxxxy[k] * ab_x + g_y_0_xxyy_xxxxxy[k];

                g_y_0_xxxyy_xxxxz[k] = -g_y_0_xxyy_xxxxz[k] * ab_x + g_y_0_xxyy_xxxxxz[k];

                g_y_0_xxxyy_xxxyy[k] = -g_y_0_xxyy_xxxyy[k] * ab_x + g_y_0_xxyy_xxxxyy[k];

                g_y_0_xxxyy_xxxyz[k] = -g_y_0_xxyy_xxxyz[k] * ab_x + g_y_0_xxyy_xxxxyz[k];

                g_y_0_xxxyy_xxxzz[k] = -g_y_0_xxyy_xxxzz[k] * ab_x + g_y_0_xxyy_xxxxzz[k];

                g_y_0_xxxyy_xxyyy[k] = -g_y_0_xxyy_xxyyy[k] * ab_x + g_y_0_xxyy_xxxyyy[k];

                g_y_0_xxxyy_xxyyz[k] = -g_y_0_xxyy_xxyyz[k] * ab_x + g_y_0_xxyy_xxxyyz[k];

                g_y_0_xxxyy_xxyzz[k] = -g_y_0_xxyy_xxyzz[k] * ab_x + g_y_0_xxyy_xxxyzz[k];

                g_y_0_xxxyy_xxzzz[k] = -g_y_0_xxyy_xxzzz[k] * ab_x + g_y_0_xxyy_xxxzzz[k];

                g_y_0_xxxyy_xyyyy[k] = -g_y_0_xxyy_xyyyy[k] * ab_x + g_y_0_xxyy_xxyyyy[k];

                g_y_0_xxxyy_xyyyz[k] = -g_y_0_xxyy_xyyyz[k] * ab_x + g_y_0_xxyy_xxyyyz[k];

                g_y_0_xxxyy_xyyzz[k] = -g_y_0_xxyy_xyyzz[k] * ab_x + g_y_0_xxyy_xxyyzz[k];

                g_y_0_xxxyy_xyzzz[k] = -g_y_0_xxyy_xyzzz[k] * ab_x + g_y_0_xxyy_xxyzzz[k];

                g_y_0_xxxyy_xzzzz[k] = -g_y_0_xxyy_xzzzz[k] * ab_x + g_y_0_xxyy_xxzzzz[k];

                g_y_0_xxxyy_yyyyy[k] = -g_y_0_xxyy_yyyyy[k] * ab_x + g_y_0_xxyy_xyyyyy[k];

                g_y_0_xxxyy_yyyyz[k] = -g_y_0_xxyy_yyyyz[k] * ab_x + g_y_0_xxyy_xyyyyz[k];

                g_y_0_xxxyy_yyyzz[k] = -g_y_0_xxyy_yyyzz[k] * ab_x + g_y_0_xxyy_xyyyzz[k];

                g_y_0_xxxyy_yyzzz[k] = -g_y_0_xxyy_yyzzz[k] * ab_x + g_y_0_xxyy_xyyzzz[k];

                g_y_0_xxxyy_yzzzz[k] = -g_y_0_xxyy_yzzzz[k] * ab_x + g_y_0_xxyy_xyzzzz[k];

                g_y_0_xxxyy_zzzzz[k] = -g_y_0_xxyy_zzzzz[k] * ab_x + g_y_0_xxyy_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_yyyyy, g_y_0_xxxyz_yyyyz, g_y_0_xxxyz_yyyzz, g_y_0_xxxyz_yyzzz, g_y_0_xxxyz_yzzzz, g_y_0_xxxyz_zzzzz, g_y_0_xxyz_xxxxx, g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxy, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxyy, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxyyy, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xyyyy, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_yyyyy, g_y_0_xxyz_yyyyz, g_y_0_xxyz_yyyzz, g_y_0_xxyz_yyzzz, g_y_0_xxyz_yzzzz, g_y_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxxx[k] = -g_y_0_xxyz_xxxxx[k] * ab_x + g_y_0_xxyz_xxxxxx[k];

                g_y_0_xxxyz_xxxxy[k] = -g_y_0_xxyz_xxxxy[k] * ab_x + g_y_0_xxyz_xxxxxy[k];

                g_y_0_xxxyz_xxxxz[k] = -g_y_0_xxyz_xxxxz[k] * ab_x + g_y_0_xxyz_xxxxxz[k];

                g_y_0_xxxyz_xxxyy[k] = -g_y_0_xxyz_xxxyy[k] * ab_x + g_y_0_xxyz_xxxxyy[k];

                g_y_0_xxxyz_xxxyz[k] = -g_y_0_xxyz_xxxyz[k] * ab_x + g_y_0_xxyz_xxxxyz[k];

                g_y_0_xxxyz_xxxzz[k] = -g_y_0_xxyz_xxxzz[k] * ab_x + g_y_0_xxyz_xxxxzz[k];

                g_y_0_xxxyz_xxyyy[k] = -g_y_0_xxyz_xxyyy[k] * ab_x + g_y_0_xxyz_xxxyyy[k];

                g_y_0_xxxyz_xxyyz[k] = -g_y_0_xxyz_xxyyz[k] * ab_x + g_y_0_xxyz_xxxyyz[k];

                g_y_0_xxxyz_xxyzz[k] = -g_y_0_xxyz_xxyzz[k] * ab_x + g_y_0_xxyz_xxxyzz[k];

                g_y_0_xxxyz_xxzzz[k] = -g_y_0_xxyz_xxzzz[k] * ab_x + g_y_0_xxyz_xxxzzz[k];

                g_y_0_xxxyz_xyyyy[k] = -g_y_0_xxyz_xyyyy[k] * ab_x + g_y_0_xxyz_xxyyyy[k];

                g_y_0_xxxyz_xyyyz[k] = -g_y_0_xxyz_xyyyz[k] * ab_x + g_y_0_xxyz_xxyyyz[k];

                g_y_0_xxxyz_xyyzz[k] = -g_y_0_xxyz_xyyzz[k] * ab_x + g_y_0_xxyz_xxyyzz[k];

                g_y_0_xxxyz_xyzzz[k] = -g_y_0_xxyz_xyzzz[k] * ab_x + g_y_0_xxyz_xxyzzz[k];

                g_y_0_xxxyz_xzzzz[k] = -g_y_0_xxyz_xzzzz[k] * ab_x + g_y_0_xxyz_xxzzzz[k];

                g_y_0_xxxyz_yyyyy[k] = -g_y_0_xxyz_yyyyy[k] * ab_x + g_y_0_xxyz_xyyyyy[k];

                g_y_0_xxxyz_yyyyz[k] = -g_y_0_xxyz_yyyyz[k] * ab_x + g_y_0_xxyz_xyyyyz[k];

                g_y_0_xxxyz_yyyzz[k] = -g_y_0_xxyz_yyyzz[k] * ab_x + g_y_0_xxyz_xyyyzz[k];

                g_y_0_xxxyz_yyzzz[k] = -g_y_0_xxyz_yyzzz[k] * ab_x + g_y_0_xxyz_xyyzzz[k];

                g_y_0_xxxyz_yzzzz[k] = -g_y_0_xxyz_yzzzz[k] * ab_x + g_y_0_xxyz_xyzzzz[k];

                g_y_0_xxxyz_zzzzz[k] = -g_y_0_xxyz_zzzzz[k] * ab_x + g_y_0_xxyz_xzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_yyyyy, g_y_0_xxxzz_yyyyz, g_y_0_xxxzz_yyyzz, g_y_0_xxxzz_yyzzz, g_y_0_xxxzz_yzzzz, g_y_0_xxxzz_zzzzz, g_y_0_xxzz_xxxxx, g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxy, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxyy, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxyyy, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xyyyy, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_yyyyy, g_y_0_xxzz_yyyyz, g_y_0_xxzz_yyyzz, g_y_0_xxzz_yyzzz, g_y_0_xxzz_yzzzz, g_y_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxxx[k] = -g_y_0_xxzz_xxxxx[k] * ab_x + g_y_0_xxzz_xxxxxx[k];

                g_y_0_xxxzz_xxxxy[k] = -g_y_0_xxzz_xxxxy[k] * ab_x + g_y_0_xxzz_xxxxxy[k];

                g_y_0_xxxzz_xxxxz[k] = -g_y_0_xxzz_xxxxz[k] * ab_x + g_y_0_xxzz_xxxxxz[k];

                g_y_0_xxxzz_xxxyy[k] = -g_y_0_xxzz_xxxyy[k] * ab_x + g_y_0_xxzz_xxxxyy[k];

                g_y_0_xxxzz_xxxyz[k] = -g_y_0_xxzz_xxxyz[k] * ab_x + g_y_0_xxzz_xxxxyz[k];

                g_y_0_xxxzz_xxxzz[k] = -g_y_0_xxzz_xxxzz[k] * ab_x + g_y_0_xxzz_xxxxzz[k];

                g_y_0_xxxzz_xxyyy[k] = -g_y_0_xxzz_xxyyy[k] * ab_x + g_y_0_xxzz_xxxyyy[k];

                g_y_0_xxxzz_xxyyz[k] = -g_y_0_xxzz_xxyyz[k] * ab_x + g_y_0_xxzz_xxxyyz[k];

                g_y_0_xxxzz_xxyzz[k] = -g_y_0_xxzz_xxyzz[k] * ab_x + g_y_0_xxzz_xxxyzz[k];

                g_y_0_xxxzz_xxzzz[k] = -g_y_0_xxzz_xxzzz[k] * ab_x + g_y_0_xxzz_xxxzzz[k];

                g_y_0_xxxzz_xyyyy[k] = -g_y_0_xxzz_xyyyy[k] * ab_x + g_y_0_xxzz_xxyyyy[k];

                g_y_0_xxxzz_xyyyz[k] = -g_y_0_xxzz_xyyyz[k] * ab_x + g_y_0_xxzz_xxyyyz[k];

                g_y_0_xxxzz_xyyzz[k] = -g_y_0_xxzz_xyyzz[k] * ab_x + g_y_0_xxzz_xxyyzz[k];

                g_y_0_xxxzz_xyzzz[k] = -g_y_0_xxzz_xyzzz[k] * ab_x + g_y_0_xxzz_xxyzzz[k];

                g_y_0_xxxzz_xzzzz[k] = -g_y_0_xxzz_xzzzz[k] * ab_x + g_y_0_xxzz_xxzzzz[k];

                g_y_0_xxxzz_yyyyy[k] = -g_y_0_xxzz_yyyyy[k] * ab_x + g_y_0_xxzz_xyyyyy[k];

                g_y_0_xxxzz_yyyyz[k] = -g_y_0_xxzz_yyyyz[k] * ab_x + g_y_0_xxzz_xyyyyz[k];

                g_y_0_xxxzz_yyyzz[k] = -g_y_0_xxzz_yyyzz[k] * ab_x + g_y_0_xxzz_xyyyzz[k];

                g_y_0_xxxzz_yyzzz[k] = -g_y_0_xxzz_yyzzz[k] * ab_x + g_y_0_xxzz_xyyzzz[k];

                g_y_0_xxxzz_yzzzz[k] = -g_y_0_xxzz_yzzzz[k] * ab_x + g_y_0_xxzz_xyzzzz[k];

                g_y_0_xxxzz_zzzzz[k] = -g_y_0_xxzz_zzzzz[k] * ab_x + g_y_0_xxzz_xzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_yyyyy, g_y_0_xxyyy_yyyyz, g_y_0_xxyyy_yyyzz, g_y_0_xxyyy_yyzzz, g_y_0_xxyyy_yzzzz, g_y_0_xxyyy_zzzzz, g_y_0_xyyy_xxxxx, g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxy, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxyy, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxyyy, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xyyyy, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_yyyyy, g_y_0_xyyy_yyyyz, g_y_0_xyyy_yyyzz, g_y_0_xyyy_yyzzz, g_y_0_xyyy_yzzzz, g_y_0_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxxx[k] = -g_y_0_xyyy_xxxxx[k] * ab_x + g_y_0_xyyy_xxxxxx[k];

                g_y_0_xxyyy_xxxxy[k] = -g_y_0_xyyy_xxxxy[k] * ab_x + g_y_0_xyyy_xxxxxy[k];

                g_y_0_xxyyy_xxxxz[k] = -g_y_0_xyyy_xxxxz[k] * ab_x + g_y_0_xyyy_xxxxxz[k];

                g_y_0_xxyyy_xxxyy[k] = -g_y_0_xyyy_xxxyy[k] * ab_x + g_y_0_xyyy_xxxxyy[k];

                g_y_0_xxyyy_xxxyz[k] = -g_y_0_xyyy_xxxyz[k] * ab_x + g_y_0_xyyy_xxxxyz[k];

                g_y_0_xxyyy_xxxzz[k] = -g_y_0_xyyy_xxxzz[k] * ab_x + g_y_0_xyyy_xxxxzz[k];

                g_y_0_xxyyy_xxyyy[k] = -g_y_0_xyyy_xxyyy[k] * ab_x + g_y_0_xyyy_xxxyyy[k];

                g_y_0_xxyyy_xxyyz[k] = -g_y_0_xyyy_xxyyz[k] * ab_x + g_y_0_xyyy_xxxyyz[k];

                g_y_0_xxyyy_xxyzz[k] = -g_y_0_xyyy_xxyzz[k] * ab_x + g_y_0_xyyy_xxxyzz[k];

                g_y_0_xxyyy_xxzzz[k] = -g_y_0_xyyy_xxzzz[k] * ab_x + g_y_0_xyyy_xxxzzz[k];

                g_y_0_xxyyy_xyyyy[k] = -g_y_0_xyyy_xyyyy[k] * ab_x + g_y_0_xyyy_xxyyyy[k];

                g_y_0_xxyyy_xyyyz[k] = -g_y_0_xyyy_xyyyz[k] * ab_x + g_y_0_xyyy_xxyyyz[k];

                g_y_0_xxyyy_xyyzz[k] = -g_y_0_xyyy_xyyzz[k] * ab_x + g_y_0_xyyy_xxyyzz[k];

                g_y_0_xxyyy_xyzzz[k] = -g_y_0_xyyy_xyzzz[k] * ab_x + g_y_0_xyyy_xxyzzz[k];

                g_y_0_xxyyy_xzzzz[k] = -g_y_0_xyyy_xzzzz[k] * ab_x + g_y_0_xyyy_xxzzzz[k];

                g_y_0_xxyyy_yyyyy[k] = -g_y_0_xyyy_yyyyy[k] * ab_x + g_y_0_xyyy_xyyyyy[k];

                g_y_0_xxyyy_yyyyz[k] = -g_y_0_xyyy_yyyyz[k] * ab_x + g_y_0_xyyy_xyyyyz[k];

                g_y_0_xxyyy_yyyzz[k] = -g_y_0_xyyy_yyyzz[k] * ab_x + g_y_0_xyyy_xyyyzz[k];

                g_y_0_xxyyy_yyzzz[k] = -g_y_0_xyyy_yyzzz[k] * ab_x + g_y_0_xyyy_xyyzzz[k];

                g_y_0_xxyyy_yzzzz[k] = -g_y_0_xyyy_yzzzz[k] * ab_x + g_y_0_xyyy_xyzzzz[k];

                g_y_0_xxyyy_zzzzz[k] = -g_y_0_xyyy_zzzzz[k] * ab_x + g_y_0_xyyy_xzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_yyyyy, g_y_0_xxyyz_yyyyz, g_y_0_xxyyz_yyyzz, g_y_0_xxyyz_yyzzz, g_y_0_xxyyz_yzzzz, g_y_0_xxyyz_zzzzz, g_y_0_xyyz_xxxxx, g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxy, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxyy, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxyyy, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xyyyy, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_yyyyy, g_y_0_xyyz_yyyyz, g_y_0_xyyz_yyyzz, g_y_0_xyyz_yyzzz, g_y_0_xyyz_yzzzz, g_y_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxxx[k] = -g_y_0_xyyz_xxxxx[k] * ab_x + g_y_0_xyyz_xxxxxx[k];

                g_y_0_xxyyz_xxxxy[k] = -g_y_0_xyyz_xxxxy[k] * ab_x + g_y_0_xyyz_xxxxxy[k];

                g_y_0_xxyyz_xxxxz[k] = -g_y_0_xyyz_xxxxz[k] * ab_x + g_y_0_xyyz_xxxxxz[k];

                g_y_0_xxyyz_xxxyy[k] = -g_y_0_xyyz_xxxyy[k] * ab_x + g_y_0_xyyz_xxxxyy[k];

                g_y_0_xxyyz_xxxyz[k] = -g_y_0_xyyz_xxxyz[k] * ab_x + g_y_0_xyyz_xxxxyz[k];

                g_y_0_xxyyz_xxxzz[k] = -g_y_0_xyyz_xxxzz[k] * ab_x + g_y_0_xyyz_xxxxzz[k];

                g_y_0_xxyyz_xxyyy[k] = -g_y_0_xyyz_xxyyy[k] * ab_x + g_y_0_xyyz_xxxyyy[k];

                g_y_0_xxyyz_xxyyz[k] = -g_y_0_xyyz_xxyyz[k] * ab_x + g_y_0_xyyz_xxxyyz[k];

                g_y_0_xxyyz_xxyzz[k] = -g_y_0_xyyz_xxyzz[k] * ab_x + g_y_0_xyyz_xxxyzz[k];

                g_y_0_xxyyz_xxzzz[k] = -g_y_0_xyyz_xxzzz[k] * ab_x + g_y_0_xyyz_xxxzzz[k];

                g_y_0_xxyyz_xyyyy[k] = -g_y_0_xyyz_xyyyy[k] * ab_x + g_y_0_xyyz_xxyyyy[k];

                g_y_0_xxyyz_xyyyz[k] = -g_y_0_xyyz_xyyyz[k] * ab_x + g_y_0_xyyz_xxyyyz[k];

                g_y_0_xxyyz_xyyzz[k] = -g_y_0_xyyz_xyyzz[k] * ab_x + g_y_0_xyyz_xxyyzz[k];

                g_y_0_xxyyz_xyzzz[k] = -g_y_0_xyyz_xyzzz[k] * ab_x + g_y_0_xyyz_xxyzzz[k];

                g_y_0_xxyyz_xzzzz[k] = -g_y_0_xyyz_xzzzz[k] * ab_x + g_y_0_xyyz_xxzzzz[k];

                g_y_0_xxyyz_yyyyy[k] = -g_y_0_xyyz_yyyyy[k] * ab_x + g_y_0_xyyz_xyyyyy[k];

                g_y_0_xxyyz_yyyyz[k] = -g_y_0_xyyz_yyyyz[k] * ab_x + g_y_0_xyyz_xyyyyz[k];

                g_y_0_xxyyz_yyyzz[k] = -g_y_0_xyyz_yyyzz[k] * ab_x + g_y_0_xyyz_xyyyzz[k];

                g_y_0_xxyyz_yyzzz[k] = -g_y_0_xyyz_yyzzz[k] * ab_x + g_y_0_xyyz_xyyzzz[k];

                g_y_0_xxyyz_yzzzz[k] = -g_y_0_xyyz_yzzzz[k] * ab_x + g_y_0_xyyz_xyzzzz[k];

                g_y_0_xxyyz_zzzzz[k] = -g_y_0_xyyz_zzzzz[k] * ab_x + g_y_0_xyyz_xzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_yyyyy, g_y_0_xxyzz_yyyyz, g_y_0_xxyzz_yyyzz, g_y_0_xxyzz_yyzzz, g_y_0_xxyzz_yzzzz, g_y_0_xxyzz_zzzzz, g_y_0_xyzz_xxxxx, g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxy, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxyy, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxyyy, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xyyyy, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_yyyyy, g_y_0_xyzz_yyyyz, g_y_0_xyzz_yyyzz, g_y_0_xyzz_yyzzz, g_y_0_xyzz_yzzzz, g_y_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxxx[k] = -g_y_0_xyzz_xxxxx[k] * ab_x + g_y_0_xyzz_xxxxxx[k];

                g_y_0_xxyzz_xxxxy[k] = -g_y_0_xyzz_xxxxy[k] * ab_x + g_y_0_xyzz_xxxxxy[k];

                g_y_0_xxyzz_xxxxz[k] = -g_y_0_xyzz_xxxxz[k] * ab_x + g_y_0_xyzz_xxxxxz[k];

                g_y_0_xxyzz_xxxyy[k] = -g_y_0_xyzz_xxxyy[k] * ab_x + g_y_0_xyzz_xxxxyy[k];

                g_y_0_xxyzz_xxxyz[k] = -g_y_0_xyzz_xxxyz[k] * ab_x + g_y_0_xyzz_xxxxyz[k];

                g_y_0_xxyzz_xxxzz[k] = -g_y_0_xyzz_xxxzz[k] * ab_x + g_y_0_xyzz_xxxxzz[k];

                g_y_0_xxyzz_xxyyy[k] = -g_y_0_xyzz_xxyyy[k] * ab_x + g_y_0_xyzz_xxxyyy[k];

                g_y_0_xxyzz_xxyyz[k] = -g_y_0_xyzz_xxyyz[k] * ab_x + g_y_0_xyzz_xxxyyz[k];

                g_y_0_xxyzz_xxyzz[k] = -g_y_0_xyzz_xxyzz[k] * ab_x + g_y_0_xyzz_xxxyzz[k];

                g_y_0_xxyzz_xxzzz[k] = -g_y_0_xyzz_xxzzz[k] * ab_x + g_y_0_xyzz_xxxzzz[k];

                g_y_0_xxyzz_xyyyy[k] = -g_y_0_xyzz_xyyyy[k] * ab_x + g_y_0_xyzz_xxyyyy[k];

                g_y_0_xxyzz_xyyyz[k] = -g_y_0_xyzz_xyyyz[k] * ab_x + g_y_0_xyzz_xxyyyz[k];

                g_y_0_xxyzz_xyyzz[k] = -g_y_0_xyzz_xyyzz[k] * ab_x + g_y_0_xyzz_xxyyzz[k];

                g_y_0_xxyzz_xyzzz[k] = -g_y_0_xyzz_xyzzz[k] * ab_x + g_y_0_xyzz_xxyzzz[k];

                g_y_0_xxyzz_xzzzz[k] = -g_y_0_xyzz_xzzzz[k] * ab_x + g_y_0_xyzz_xxzzzz[k];

                g_y_0_xxyzz_yyyyy[k] = -g_y_0_xyzz_yyyyy[k] * ab_x + g_y_0_xyzz_xyyyyy[k];

                g_y_0_xxyzz_yyyyz[k] = -g_y_0_xyzz_yyyyz[k] * ab_x + g_y_0_xyzz_xyyyyz[k];

                g_y_0_xxyzz_yyyzz[k] = -g_y_0_xyzz_yyyzz[k] * ab_x + g_y_0_xyzz_xyyyzz[k];

                g_y_0_xxyzz_yyzzz[k] = -g_y_0_xyzz_yyzzz[k] * ab_x + g_y_0_xyzz_xyyzzz[k];

                g_y_0_xxyzz_yzzzz[k] = -g_y_0_xyzz_yzzzz[k] * ab_x + g_y_0_xyzz_xyzzzz[k];

                g_y_0_xxyzz_zzzzz[k] = -g_y_0_xyzz_zzzzz[k] * ab_x + g_y_0_xyzz_xzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_yyyyy, g_y_0_xxzzz_yyyyz, g_y_0_xxzzz_yyyzz, g_y_0_xxzzz_yyzzz, g_y_0_xxzzz_yzzzz, g_y_0_xxzzz_zzzzz, g_y_0_xzzz_xxxxx, g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxy, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxyy, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxyyy, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xyyyy, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_yyyyy, g_y_0_xzzz_yyyyz, g_y_0_xzzz_yyyzz, g_y_0_xzzz_yyzzz, g_y_0_xzzz_yzzzz, g_y_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxxx[k] = -g_y_0_xzzz_xxxxx[k] * ab_x + g_y_0_xzzz_xxxxxx[k];

                g_y_0_xxzzz_xxxxy[k] = -g_y_0_xzzz_xxxxy[k] * ab_x + g_y_0_xzzz_xxxxxy[k];

                g_y_0_xxzzz_xxxxz[k] = -g_y_0_xzzz_xxxxz[k] * ab_x + g_y_0_xzzz_xxxxxz[k];

                g_y_0_xxzzz_xxxyy[k] = -g_y_0_xzzz_xxxyy[k] * ab_x + g_y_0_xzzz_xxxxyy[k];

                g_y_0_xxzzz_xxxyz[k] = -g_y_0_xzzz_xxxyz[k] * ab_x + g_y_0_xzzz_xxxxyz[k];

                g_y_0_xxzzz_xxxzz[k] = -g_y_0_xzzz_xxxzz[k] * ab_x + g_y_0_xzzz_xxxxzz[k];

                g_y_0_xxzzz_xxyyy[k] = -g_y_0_xzzz_xxyyy[k] * ab_x + g_y_0_xzzz_xxxyyy[k];

                g_y_0_xxzzz_xxyyz[k] = -g_y_0_xzzz_xxyyz[k] * ab_x + g_y_0_xzzz_xxxyyz[k];

                g_y_0_xxzzz_xxyzz[k] = -g_y_0_xzzz_xxyzz[k] * ab_x + g_y_0_xzzz_xxxyzz[k];

                g_y_0_xxzzz_xxzzz[k] = -g_y_0_xzzz_xxzzz[k] * ab_x + g_y_0_xzzz_xxxzzz[k];

                g_y_0_xxzzz_xyyyy[k] = -g_y_0_xzzz_xyyyy[k] * ab_x + g_y_0_xzzz_xxyyyy[k];

                g_y_0_xxzzz_xyyyz[k] = -g_y_0_xzzz_xyyyz[k] * ab_x + g_y_0_xzzz_xxyyyz[k];

                g_y_0_xxzzz_xyyzz[k] = -g_y_0_xzzz_xyyzz[k] * ab_x + g_y_0_xzzz_xxyyzz[k];

                g_y_0_xxzzz_xyzzz[k] = -g_y_0_xzzz_xyzzz[k] * ab_x + g_y_0_xzzz_xxyzzz[k];

                g_y_0_xxzzz_xzzzz[k] = -g_y_0_xzzz_xzzzz[k] * ab_x + g_y_0_xzzz_xxzzzz[k];

                g_y_0_xxzzz_yyyyy[k] = -g_y_0_xzzz_yyyyy[k] * ab_x + g_y_0_xzzz_xyyyyy[k];

                g_y_0_xxzzz_yyyyz[k] = -g_y_0_xzzz_yyyyz[k] * ab_x + g_y_0_xzzz_xyyyyz[k];

                g_y_0_xxzzz_yyyzz[k] = -g_y_0_xzzz_yyyzz[k] * ab_x + g_y_0_xzzz_xyyyzz[k];

                g_y_0_xxzzz_yyzzz[k] = -g_y_0_xzzz_yyzzz[k] * ab_x + g_y_0_xzzz_xyyzzz[k];

                g_y_0_xxzzz_yzzzz[k] = -g_y_0_xzzz_yzzzz[k] * ab_x + g_y_0_xzzz_xyzzzz[k];

                g_y_0_xxzzz_zzzzz[k] = -g_y_0_xzzz_zzzzz[k] * ab_x + g_y_0_xzzz_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_yyyyy, g_y_0_xyyyy_yyyyz, g_y_0_xyyyy_yyyzz, g_y_0_xyyyy_yyzzz, g_y_0_xyyyy_yzzzz, g_y_0_xyyyy_zzzzz, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxxx[k] = -g_y_0_yyyy_xxxxx[k] * ab_x + g_y_0_yyyy_xxxxxx[k];

                g_y_0_xyyyy_xxxxy[k] = -g_y_0_yyyy_xxxxy[k] * ab_x + g_y_0_yyyy_xxxxxy[k];

                g_y_0_xyyyy_xxxxz[k] = -g_y_0_yyyy_xxxxz[k] * ab_x + g_y_0_yyyy_xxxxxz[k];

                g_y_0_xyyyy_xxxyy[k] = -g_y_0_yyyy_xxxyy[k] * ab_x + g_y_0_yyyy_xxxxyy[k];

                g_y_0_xyyyy_xxxyz[k] = -g_y_0_yyyy_xxxyz[k] * ab_x + g_y_0_yyyy_xxxxyz[k];

                g_y_0_xyyyy_xxxzz[k] = -g_y_0_yyyy_xxxzz[k] * ab_x + g_y_0_yyyy_xxxxzz[k];

                g_y_0_xyyyy_xxyyy[k] = -g_y_0_yyyy_xxyyy[k] * ab_x + g_y_0_yyyy_xxxyyy[k];

                g_y_0_xyyyy_xxyyz[k] = -g_y_0_yyyy_xxyyz[k] * ab_x + g_y_0_yyyy_xxxyyz[k];

                g_y_0_xyyyy_xxyzz[k] = -g_y_0_yyyy_xxyzz[k] * ab_x + g_y_0_yyyy_xxxyzz[k];

                g_y_0_xyyyy_xxzzz[k] = -g_y_0_yyyy_xxzzz[k] * ab_x + g_y_0_yyyy_xxxzzz[k];

                g_y_0_xyyyy_xyyyy[k] = -g_y_0_yyyy_xyyyy[k] * ab_x + g_y_0_yyyy_xxyyyy[k];

                g_y_0_xyyyy_xyyyz[k] = -g_y_0_yyyy_xyyyz[k] * ab_x + g_y_0_yyyy_xxyyyz[k];

                g_y_0_xyyyy_xyyzz[k] = -g_y_0_yyyy_xyyzz[k] * ab_x + g_y_0_yyyy_xxyyzz[k];

                g_y_0_xyyyy_xyzzz[k] = -g_y_0_yyyy_xyzzz[k] * ab_x + g_y_0_yyyy_xxyzzz[k];

                g_y_0_xyyyy_xzzzz[k] = -g_y_0_yyyy_xzzzz[k] * ab_x + g_y_0_yyyy_xxzzzz[k];

                g_y_0_xyyyy_yyyyy[k] = -g_y_0_yyyy_yyyyy[k] * ab_x + g_y_0_yyyy_xyyyyy[k];

                g_y_0_xyyyy_yyyyz[k] = -g_y_0_yyyy_yyyyz[k] * ab_x + g_y_0_yyyy_xyyyyz[k];

                g_y_0_xyyyy_yyyzz[k] = -g_y_0_yyyy_yyyzz[k] * ab_x + g_y_0_yyyy_xyyyzz[k];

                g_y_0_xyyyy_yyzzz[k] = -g_y_0_yyyy_yyzzz[k] * ab_x + g_y_0_yyyy_xyyzzz[k];

                g_y_0_xyyyy_yzzzz[k] = -g_y_0_yyyy_yzzzz[k] * ab_x + g_y_0_yyyy_xyzzzz[k];

                g_y_0_xyyyy_zzzzz[k] = -g_y_0_yyyy_zzzzz[k] * ab_x + g_y_0_yyyy_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_yyyyy, g_y_0_xyyyz_yyyyz, g_y_0_xyyyz_yyyzz, g_y_0_xyyyz_yyzzz, g_y_0_xyyyz_yzzzz, g_y_0_xyyyz_zzzzz, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyy, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxxx[k] = -g_y_0_yyyz_xxxxx[k] * ab_x + g_y_0_yyyz_xxxxxx[k];

                g_y_0_xyyyz_xxxxy[k] = -g_y_0_yyyz_xxxxy[k] * ab_x + g_y_0_yyyz_xxxxxy[k];

                g_y_0_xyyyz_xxxxz[k] = -g_y_0_yyyz_xxxxz[k] * ab_x + g_y_0_yyyz_xxxxxz[k];

                g_y_0_xyyyz_xxxyy[k] = -g_y_0_yyyz_xxxyy[k] * ab_x + g_y_0_yyyz_xxxxyy[k];

                g_y_0_xyyyz_xxxyz[k] = -g_y_0_yyyz_xxxyz[k] * ab_x + g_y_0_yyyz_xxxxyz[k];

                g_y_0_xyyyz_xxxzz[k] = -g_y_0_yyyz_xxxzz[k] * ab_x + g_y_0_yyyz_xxxxzz[k];

                g_y_0_xyyyz_xxyyy[k] = -g_y_0_yyyz_xxyyy[k] * ab_x + g_y_0_yyyz_xxxyyy[k];

                g_y_0_xyyyz_xxyyz[k] = -g_y_0_yyyz_xxyyz[k] * ab_x + g_y_0_yyyz_xxxyyz[k];

                g_y_0_xyyyz_xxyzz[k] = -g_y_0_yyyz_xxyzz[k] * ab_x + g_y_0_yyyz_xxxyzz[k];

                g_y_0_xyyyz_xxzzz[k] = -g_y_0_yyyz_xxzzz[k] * ab_x + g_y_0_yyyz_xxxzzz[k];

                g_y_0_xyyyz_xyyyy[k] = -g_y_0_yyyz_xyyyy[k] * ab_x + g_y_0_yyyz_xxyyyy[k];

                g_y_0_xyyyz_xyyyz[k] = -g_y_0_yyyz_xyyyz[k] * ab_x + g_y_0_yyyz_xxyyyz[k];

                g_y_0_xyyyz_xyyzz[k] = -g_y_0_yyyz_xyyzz[k] * ab_x + g_y_0_yyyz_xxyyzz[k];

                g_y_0_xyyyz_xyzzz[k] = -g_y_0_yyyz_xyzzz[k] * ab_x + g_y_0_yyyz_xxyzzz[k];

                g_y_0_xyyyz_xzzzz[k] = -g_y_0_yyyz_xzzzz[k] * ab_x + g_y_0_yyyz_xxzzzz[k];

                g_y_0_xyyyz_yyyyy[k] = -g_y_0_yyyz_yyyyy[k] * ab_x + g_y_0_yyyz_xyyyyy[k];

                g_y_0_xyyyz_yyyyz[k] = -g_y_0_yyyz_yyyyz[k] * ab_x + g_y_0_yyyz_xyyyyz[k];

                g_y_0_xyyyz_yyyzz[k] = -g_y_0_yyyz_yyyzz[k] * ab_x + g_y_0_yyyz_xyyyzz[k];

                g_y_0_xyyyz_yyzzz[k] = -g_y_0_yyyz_yyzzz[k] * ab_x + g_y_0_yyyz_xyyzzz[k];

                g_y_0_xyyyz_yzzzz[k] = -g_y_0_yyyz_yzzzz[k] * ab_x + g_y_0_yyyz_xyzzzz[k];

                g_y_0_xyyyz_zzzzz[k] = -g_y_0_yyyz_zzzzz[k] * ab_x + g_y_0_yyyz_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_yyyyy, g_y_0_xyyzz_yyyyz, g_y_0_xyyzz_yyyzz, g_y_0_xyyzz_yyzzz, g_y_0_xyyzz_yzzzz, g_y_0_xyyzz_zzzzz, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyy, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxxx[k] = -g_y_0_yyzz_xxxxx[k] * ab_x + g_y_0_yyzz_xxxxxx[k];

                g_y_0_xyyzz_xxxxy[k] = -g_y_0_yyzz_xxxxy[k] * ab_x + g_y_0_yyzz_xxxxxy[k];

                g_y_0_xyyzz_xxxxz[k] = -g_y_0_yyzz_xxxxz[k] * ab_x + g_y_0_yyzz_xxxxxz[k];

                g_y_0_xyyzz_xxxyy[k] = -g_y_0_yyzz_xxxyy[k] * ab_x + g_y_0_yyzz_xxxxyy[k];

                g_y_0_xyyzz_xxxyz[k] = -g_y_0_yyzz_xxxyz[k] * ab_x + g_y_0_yyzz_xxxxyz[k];

                g_y_0_xyyzz_xxxzz[k] = -g_y_0_yyzz_xxxzz[k] * ab_x + g_y_0_yyzz_xxxxzz[k];

                g_y_0_xyyzz_xxyyy[k] = -g_y_0_yyzz_xxyyy[k] * ab_x + g_y_0_yyzz_xxxyyy[k];

                g_y_0_xyyzz_xxyyz[k] = -g_y_0_yyzz_xxyyz[k] * ab_x + g_y_0_yyzz_xxxyyz[k];

                g_y_0_xyyzz_xxyzz[k] = -g_y_0_yyzz_xxyzz[k] * ab_x + g_y_0_yyzz_xxxyzz[k];

                g_y_0_xyyzz_xxzzz[k] = -g_y_0_yyzz_xxzzz[k] * ab_x + g_y_0_yyzz_xxxzzz[k];

                g_y_0_xyyzz_xyyyy[k] = -g_y_0_yyzz_xyyyy[k] * ab_x + g_y_0_yyzz_xxyyyy[k];

                g_y_0_xyyzz_xyyyz[k] = -g_y_0_yyzz_xyyyz[k] * ab_x + g_y_0_yyzz_xxyyyz[k];

                g_y_0_xyyzz_xyyzz[k] = -g_y_0_yyzz_xyyzz[k] * ab_x + g_y_0_yyzz_xxyyzz[k];

                g_y_0_xyyzz_xyzzz[k] = -g_y_0_yyzz_xyzzz[k] * ab_x + g_y_0_yyzz_xxyzzz[k];

                g_y_0_xyyzz_xzzzz[k] = -g_y_0_yyzz_xzzzz[k] * ab_x + g_y_0_yyzz_xxzzzz[k];

                g_y_0_xyyzz_yyyyy[k] = -g_y_0_yyzz_yyyyy[k] * ab_x + g_y_0_yyzz_xyyyyy[k];

                g_y_0_xyyzz_yyyyz[k] = -g_y_0_yyzz_yyyyz[k] * ab_x + g_y_0_yyzz_xyyyyz[k];

                g_y_0_xyyzz_yyyzz[k] = -g_y_0_yyzz_yyyzz[k] * ab_x + g_y_0_yyzz_xyyyzz[k];

                g_y_0_xyyzz_yyzzz[k] = -g_y_0_yyzz_yyzzz[k] * ab_x + g_y_0_yyzz_xyyzzz[k];

                g_y_0_xyyzz_yzzzz[k] = -g_y_0_yyzz_yzzzz[k] * ab_x + g_y_0_yyzz_xyzzzz[k];

                g_y_0_xyyzz_zzzzz[k] = -g_y_0_yyzz_zzzzz[k] * ab_x + g_y_0_yyzz_xzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_yyyyy, g_y_0_xyzzz_yyyyz, g_y_0_xyzzz_yyyzz, g_y_0_xyzzz_yyzzz, g_y_0_xyzzz_yzzzz, g_y_0_xyzzz_zzzzz, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyy, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxxx[k] = -g_y_0_yzzz_xxxxx[k] * ab_x + g_y_0_yzzz_xxxxxx[k];

                g_y_0_xyzzz_xxxxy[k] = -g_y_0_yzzz_xxxxy[k] * ab_x + g_y_0_yzzz_xxxxxy[k];

                g_y_0_xyzzz_xxxxz[k] = -g_y_0_yzzz_xxxxz[k] * ab_x + g_y_0_yzzz_xxxxxz[k];

                g_y_0_xyzzz_xxxyy[k] = -g_y_0_yzzz_xxxyy[k] * ab_x + g_y_0_yzzz_xxxxyy[k];

                g_y_0_xyzzz_xxxyz[k] = -g_y_0_yzzz_xxxyz[k] * ab_x + g_y_0_yzzz_xxxxyz[k];

                g_y_0_xyzzz_xxxzz[k] = -g_y_0_yzzz_xxxzz[k] * ab_x + g_y_0_yzzz_xxxxzz[k];

                g_y_0_xyzzz_xxyyy[k] = -g_y_0_yzzz_xxyyy[k] * ab_x + g_y_0_yzzz_xxxyyy[k];

                g_y_0_xyzzz_xxyyz[k] = -g_y_0_yzzz_xxyyz[k] * ab_x + g_y_0_yzzz_xxxyyz[k];

                g_y_0_xyzzz_xxyzz[k] = -g_y_0_yzzz_xxyzz[k] * ab_x + g_y_0_yzzz_xxxyzz[k];

                g_y_0_xyzzz_xxzzz[k] = -g_y_0_yzzz_xxzzz[k] * ab_x + g_y_0_yzzz_xxxzzz[k];

                g_y_0_xyzzz_xyyyy[k] = -g_y_0_yzzz_xyyyy[k] * ab_x + g_y_0_yzzz_xxyyyy[k];

                g_y_0_xyzzz_xyyyz[k] = -g_y_0_yzzz_xyyyz[k] * ab_x + g_y_0_yzzz_xxyyyz[k];

                g_y_0_xyzzz_xyyzz[k] = -g_y_0_yzzz_xyyzz[k] * ab_x + g_y_0_yzzz_xxyyzz[k];

                g_y_0_xyzzz_xyzzz[k] = -g_y_0_yzzz_xyzzz[k] * ab_x + g_y_0_yzzz_xxyzzz[k];

                g_y_0_xyzzz_xzzzz[k] = -g_y_0_yzzz_xzzzz[k] * ab_x + g_y_0_yzzz_xxzzzz[k];

                g_y_0_xyzzz_yyyyy[k] = -g_y_0_yzzz_yyyyy[k] * ab_x + g_y_0_yzzz_xyyyyy[k];

                g_y_0_xyzzz_yyyyz[k] = -g_y_0_yzzz_yyyyz[k] * ab_x + g_y_0_yzzz_xyyyyz[k];

                g_y_0_xyzzz_yyyzz[k] = -g_y_0_yzzz_yyyzz[k] * ab_x + g_y_0_yzzz_xyyyzz[k];

                g_y_0_xyzzz_yyzzz[k] = -g_y_0_yzzz_yyzzz[k] * ab_x + g_y_0_yzzz_xyyzzz[k];

                g_y_0_xyzzz_yzzzz[k] = -g_y_0_yzzz_yzzzz[k] * ab_x + g_y_0_yzzz_xyzzzz[k];

                g_y_0_xyzzz_zzzzz[k] = -g_y_0_yzzz_zzzzz[k] * ab_x + g_y_0_yzzz_xzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_yyyyy, g_y_0_xzzzz_yyyyz, g_y_0_xzzzz_yyyzz, g_y_0_xzzzz_yyzzz, g_y_0_xzzzz_yzzzz, g_y_0_xzzzz_zzzzz, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyy, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxxx[k] = -g_y_0_zzzz_xxxxx[k] * ab_x + g_y_0_zzzz_xxxxxx[k];

                g_y_0_xzzzz_xxxxy[k] = -g_y_0_zzzz_xxxxy[k] * ab_x + g_y_0_zzzz_xxxxxy[k];

                g_y_0_xzzzz_xxxxz[k] = -g_y_0_zzzz_xxxxz[k] * ab_x + g_y_0_zzzz_xxxxxz[k];

                g_y_0_xzzzz_xxxyy[k] = -g_y_0_zzzz_xxxyy[k] * ab_x + g_y_0_zzzz_xxxxyy[k];

                g_y_0_xzzzz_xxxyz[k] = -g_y_0_zzzz_xxxyz[k] * ab_x + g_y_0_zzzz_xxxxyz[k];

                g_y_0_xzzzz_xxxzz[k] = -g_y_0_zzzz_xxxzz[k] * ab_x + g_y_0_zzzz_xxxxzz[k];

                g_y_0_xzzzz_xxyyy[k] = -g_y_0_zzzz_xxyyy[k] * ab_x + g_y_0_zzzz_xxxyyy[k];

                g_y_0_xzzzz_xxyyz[k] = -g_y_0_zzzz_xxyyz[k] * ab_x + g_y_0_zzzz_xxxyyz[k];

                g_y_0_xzzzz_xxyzz[k] = -g_y_0_zzzz_xxyzz[k] * ab_x + g_y_0_zzzz_xxxyzz[k];

                g_y_0_xzzzz_xxzzz[k] = -g_y_0_zzzz_xxzzz[k] * ab_x + g_y_0_zzzz_xxxzzz[k];

                g_y_0_xzzzz_xyyyy[k] = -g_y_0_zzzz_xyyyy[k] * ab_x + g_y_0_zzzz_xxyyyy[k];

                g_y_0_xzzzz_xyyyz[k] = -g_y_0_zzzz_xyyyz[k] * ab_x + g_y_0_zzzz_xxyyyz[k];

                g_y_0_xzzzz_xyyzz[k] = -g_y_0_zzzz_xyyzz[k] * ab_x + g_y_0_zzzz_xxyyzz[k];

                g_y_0_xzzzz_xyzzz[k] = -g_y_0_zzzz_xyzzz[k] * ab_x + g_y_0_zzzz_xxyzzz[k];

                g_y_0_xzzzz_xzzzz[k] = -g_y_0_zzzz_xzzzz[k] * ab_x + g_y_0_zzzz_xxzzzz[k];

                g_y_0_xzzzz_yyyyy[k] = -g_y_0_zzzz_yyyyy[k] * ab_x + g_y_0_zzzz_xyyyyy[k];

                g_y_0_xzzzz_yyyyz[k] = -g_y_0_zzzz_yyyyz[k] * ab_x + g_y_0_zzzz_xyyyyz[k];

                g_y_0_xzzzz_yyyzz[k] = -g_y_0_zzzz_yyyzz[k] * ab_x + g_y_0_zzzz_xyyyzz[k];

                g_y_0_xzzzz_yyzzz[k] = -g_y_0_zzzz_yyzzz[k] * ab_x + g_y_0_zzzz_xyyzzz[k];

                g_y_0_xzzzz_yzzzz[k] = -g_y_0_zzzz_yzzzz[k] * ab_x + g_y_0_zzzz_xyzzzz[k];

                g_y_0_xzzzz_zzzzz[k] = -g_y_0_zzzz_zzzzz[k] * ab_x + g_y_0_zzzz_xzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzzz, g_yyyy_xxxxx, g_yyyy_xxxxy, g_yyyy_xxxxz, g_yyyy_xxxyy, g_yyyy_xxxyz, g_yyyy_xxxzz, g_yyyy_xxyyy, g_yyyy_xxyyz, g_yyyy_xxyzz, g_yyyy_xxzzz, g_yyyy_xyyyy, g_yyyy_xyyyz, g_yyyy_xyyzz, g_yyyy_xyzzz, g_yyyy_xzzzz, g_yyyy_yyyyy, g_yyyy_yyyyz, g_yyyy_yyyzz, g_yyyy_yyzzz, g_yyyy_yzzzz, g_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxxx[k] = -g_yyyy_xxxxx[k] - g_y_0_yyyy_xxxxx[k] * ab_y + g_y_0_yyyy_xxxxxy[k];

                g_y_0_yyyyy_xxxxy[k] = -g_yyyy_xxxxy[k] - g_y_0_yyyy_xxxxy[k] * ab_y + g_y_0_yyyy_xxxxyy[k];

                g_y_0_yyyyy_xxxxz[k] = -g_yyyy_xxxxz[k] - g_y_0_yyyy_xxxxz[k] * ab_y + g_y_0_yyyy_xxxxyz[k];

                g_y_0_yyyyy_xxxyy[k] = -g_yyyy_xxxyy[k] - g_y_0_yyyy_xxxyy[k] * ab_y + g_y_0_yyyy_xxxyyy[k];

                g_y_0_yyyyy_xxxyz[k] = -g_yyyy_xxxyz[k] - g_y_0_yyyy_xxxyz[k] * ab_y + g_y_0_yyyy_xxxyyz[k];

                g_y_0_yyyyy_xxxzz[k] = -g_yyyy_xxxzz[k] - g_y_0_yyyy_xxxzz[k] * ab_y + g_y_0_yyyy_xxxyzz[k];

                g_y_0_yyyyy_xxyyy[k] = -g_yyyy_xxyyy[k] - g_y_0_yyyy_xxyyy[k] * ab_y + g_y_0_yyyy_xxyyyy[k];

                g_y_0_yyyyy_xxyyz[k] = -g_yyyy_xxyyz[k] - g_y_0_yyyy_xxyyz[k] * ab_y + g_y_0_yyyy_xxyyyz[k];

                g_y_0_yyyyy_xxyzz[k] = -g_yyyy_xxyzz[k] - g_y_0_yyyy_xxyzz[k] * ab_y + g_y_0_yyyy_xxyyzz[k];

                g_y_0_yyyyy_xxzzz[k] = -g_yyyy_xxzzz[k] - g_y_0_yyyy_xxzzz[k] * ab_y + g_y_0_yyyy_xxyzzz[k];

                g_y_0_yyyyy_xyyyy[k] = -g_yyyy_xyyyy[k] - g_y_0_yyyy_xyyyy[k] * ab_y + g_y_0_yyyy_xyyyyy[k];

                g_y_0_yyyyy_xyyyz[k] = -g_yyyy_xyyyz[k] - g_y_0_yyyy_xyyyz[k] * ab_y + g_y_0_yyyy_xyyyyz[k];

                g_y_0_yyyyy_xyyzz[k] = -g_yyyy_xyyzz[k] - g_y_0_yyyy_xyyzz[k] * ab_y + g_y_0_yyyy_xyyyzz[k];

                g_y_0_yyyyy_xyzzz[k] = -g_yyyy_xyzzz[k] - g_y_0_yyyy_xyzzz[k] * ab_y + g_y_0_yyyy_xyyzzz[k];

                g_y_0_yyyyy_xzzzz[k] = -g_yyyy_xzzzz[k] - g_y_0_yyyy_xzzzz[k] * ab_y + g_y_0_yyyy_xyzzzz[k];

                g_y_0_yyyyy_yyyyy[k] = -g_yyyy_yyyyy[k] - g_y_0_yyyy_yyyyy[k] * ab_y + g_y_0_yyyy_yyyyyy[k];

                g_y_0_yyyyy_yyyyz[k] = -g_yyyy_yyyyz[k] - g_y_0_yyyy_yyyyz[k] * ab_y + g_y_0_yyyy_yyyyyz[k];

                g_y_0_yyyyy_yyyzz[k] = -g_yyyy_yyyzz[k] - g_y_0_yyyy_yyyzz[k] * ab_y + g_y_0_yyyy_yyyyzz[k];

                g_y_0_yyyyy_yyzzz[k] = -g_yyyy_yyzzz[k] - g_y_0_yyyy_yyzzz[k] * ab_y + g_y_0_yyyy_yyyzzz[k];

                g_y_0_yyyyy_yzzzz[k] = -g_yyyy_yzzzz[k] - g_y_0_yyyy_yzzzz[k] * ab_y + g_y_0_yyyy_yyzzzz[k];

                g_y_0_yyyyy_zzzzz[k] = -g_yyyy_zzzzz[k] - g_y_0_yyyy_zzzzz[k] * ab_y + g_y_0_yyyy_yzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxxx[k] = -g_y_0_yyyy_xxxxx[k] * ab_z + g_y_0_yyyy_xxxxxz[k];

                g_y_0_yyyyz_xxxxy[k] = -g_y_0_yyyy_xxxxy[k] * ab_z + g_y_0_yyyy_xxxxyz[k];

                g_y_0_yyyyz_xxxxz[k] = -g_y_0_yyyy_xxxxz[k] * ab_z + g_y_0_yyyy_xxxxzz[k];

                g_y_0_yyyyz_xxxyy[k] = -g_y_0_yyyy_xxxyy[k] * ab_z + g_y_0_yyyy_xxxyyz[k];

                g_y_0_yyyyz_xxxyz[k] = -g_y_0_yyyy_xxxyz[k] * ab_z + g_y_0_yyyy_xxxyzz[k];

                g_y_0_yyyyz_xxxzz[k] = -g_y_0_yyyy_xxxzz[k] * ab_z + g_y_0_yyyy_xxxzzz[k];

                g_y_0_yyyyz_xxyyy[k] = -g_y_0_yyyy_xxyyy[k] * ab_z + g_y_0_yyyy_xxyyyz[k];

                g_y_0_yyyyz_xxyyz[k] = -g_y_0_yyyy_xxyyz[k] * ab_z + g_y_0_yyyy_xxyyzz[k];

                g_y_0_yyyyz_xxyzz[k] = -g_y_0_yyyy_xxyzz[k] * ab_z + g_y_0_yyyy_xxyzzz[k];

                g_y_0_yyyyz_xxzzz[k] = -g_y_0_yyyy_xxzzz[k] * ab_z + g_y_0_yyyy_xxzzzz[k];

                g_y_0_yyyyz_xyyyy[k] = -g_y_0_yyyy_xyyyy[k] * ab_z + g_y_0_yyyy_xyyyyz[k];

                g_y_0_yyyyz_xyyyz[k] = -g_y_0_yyyy_xyyyz[k] * ab_z + g_y_0_yyyy_xyyyzz[k];

                g_y_0_yyyyz_xyyzz[k] = -g_y_0_yyyy_xyyzz[k] * ab_z + g_y_0_yyyy_xyyzzz[k];

                g_y_0_yyyyz_xyzzz[k] = -g_y_0_yyyy_xyzzz[k] * ab_z + g_y_0_yyyy_xyzzzz[k];

                g_y_0_yyyyz_xzzzz[k] = -g_y_0_yyyy_xzzzz[k] * ab_z + g_y_0_yyyy_xzzzzz[k];

                g_y_0_yyyyz_yyyyy[k] = -g_y_0_yyyy_yyyyy[k] * ab_z + g_y_0_yyyy_yyyyyz[k];

                g_y_0_yyyyz_yyyyz[k] = -g_y_0_yyyy_yyyyz[k] * ab_z + g_y_0_yyyy_yyyyzz[k];

                g_y_0_yyyyz_yyyzz[k] = -g_y_0_yyyy_yyyzz[k] * ab_z + g_y_0_yyyy_yyyzzz[k];

                g_y_0_yyyyz_yyzzz[k] = -g_y_0_yyyy_yyzzz[k] * ab_z + g_y_0_yyyy_yyzzzz[k];

                g_y_0_yyyyz_yzzzz[k] = -g_y_0_yyyy_yzzzz[k] * ab_z + g_y_0_yyyy_yzzzzz[k];

                g_y_0_yyyyz_zzzzz[k] = -g_y_0_yyyy_zzzzz[k] * ab_z + g_y_0_yyyy_zzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzz, g_y_0_yyyz_zzzzzz, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxxx[k] = -g_y_0_yyyz_xxxxx[k] * ab_z + g_y_0_yyyz_xxxxxz[k];

                g_y_0_yyyzz_xxxxy[k] = -g_y_0_yyyz_xxxxy[k] * ab_z + g_y_0_yyyz_xxxxyz[k];

                g_y_0_yyyzz_xxxxz[k] = -g_y_0_yyyz_xxxxz[k] * ab_z + g_y_0_yyyz_xxxxzz[k];

                g_y_0_yyyzz_xxxyy[k] = -g_y_0_yyyz_xxxyy[k] * ab_z + g_y_0_yyyz_xxxyyz[k];

                g_y_0_yyyzz_xxxyz[k] = -g_y_0_yyyz_xxxyz[k] * ab_z + g_y_0_yyyz_xxxyzz[k];

                g_y_0_yyyzz_xxxzz[k] = -g_y_0_yyyz_xxxzz[k] * ab_z + g_y_0_yyyz_xxxzzz[k];

                g_y_0_yyyzz_xxyyy[k] = -g_y_0_yyyz_xxyyy[k] * ab_z + g_y_0_yyyz_xxyyyz[k];

                g_y_0_yyyzz_xxyyz[k] = -g_y_0_yyyz_xxyyz[k] * ab_z + g_y_0_yyyz_xxyyzz[k];

                g_y_0_yyyzz_xxyzz[k] = -g_y_0_yyyz_xxyzz[k] * ab_z + g_y_0_yyyz_xxyzzz[k];

                g_y_0_yyyzz_xxzzz[k] = -g_y_0_yyyz_xxzzz[k] * ab_z + g_y_0_yyyz_xxzzzz[k];

                g_y_0_yyyzz_xyyyy[k] = -g_y_0_yyyz_xyyyy[k] * ab_z + g_y_0_yyyz_xyyyyz[k];

                g_y_0_yyyzz_xyyyz[k] = -g_y_0_yyyz_xyyyz[k] * ab_z + g_y_0_yyyz_xyyyzz[k];

                g_y_0_yyyzz_xyyzz[k] = -g_y_0_yyyz_xyyzz[k] * ab_z + g_y_0_yyyz_xyyzzz[k];

                g_y_0_yyyzz_xyzzz[k] = -g_y_0_yyyz_xyzzz[k] * ab_z + g_y_0_yyyz_xyzzzz[k];

                g_y_0_yyyzz_xzzzz[k] = -g_y_0_yyyz_xzzzz[k] * ab_z + g_y_0_yyyz_xzzzzz[k];

                g_y_0_yyyzz_yyyyy[k] = -g_y_0_yyyz_yyyyy[k] * ab_z + g_y_0_yyyz_yyyyyz[k];

                g_y_0_yyyzz_yyyyz[k] = -g_y_0_yyyz_yyyyz[k] * ab_z + g_y_0_yyyz_yyyyzz[k];

                g_y_0_yyyzz_yyyzz[k] = -g_y_0_yyyz_yyyzz[k] * ab_z + g_y_0_yyyz_yyyzzz[k];

                g_y_0_yyyzz_yyzzz[k] = -g_y_0_yyyz_yyzzz[k] * ab_z + g_y_0_yyyz_yyzzzz[k];

                g_y_0_yyyzz_yzzzz[k] = -g_y_0_yyyz_yzzzz[k] * ab_z + g_y_0_yyyz_yzzzzz[k];

                g_y_0_yyyzz_zzzzz[k] = -g_y_0_yyyz_zzzzz[k] * ab_z + g_y_0_yyyz_zzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzz, g_y_0_yyzz_zzzzzz, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxxx[k] = -g_y_0_yyzz_xxxxx[k] * ab_z + g_y_0_yyzz_xxxxxz[k];

                g_y_0_yyzzz_xxxxy[k] = -g_y_0_yyzz_xxxxy[k] * ab_z + g_y_0_yyzz_xxxxyz[k];

                g_y_0_yyzzz_xxxxz[k] = -g_y_0_yyzz_xxxxz[k] * ab_z + g_y_0_yyzz_xxxxzz[k];

                g_y_0_yyzzz_xxxyy[k] = -g_y_0_yyzz_xxxyy[k] * ab_z + g_y_0_yyzz_xxxyyz[k];

                g_y_0_yyzzz_xxxyz[k] = -g_y_0_yyzz_xxxyz[k] * ab_z + g_y_0_yyzz_xxxyzz[k];

                g_y_0_yyzzz_xxxzz[k] = -g_y_0_yyzz_xxxzz[k] * ab_z + g_y_0_yyzz_xxxzzz[k];

                g_y_0_yyzzz_xxyyy[k] = -g_y_0_yyzz_xxyyy[k] * ab_z + g_y_0_yyzz_xxyyyz[k];

                g_y_0_yyzzz_xxyyz[k] = -g_y_0_yyzz_xxyyz[k] * ab_z + g_y_0_yyzz_xxyyzz[k];

                g_y_0_yyzzz_xxyzz[k] = -g_y_0_yyzz_xxyzz[k] * ab_z + g_y_0_yyzz_xxyzzz[k];

                g_y_0_yyzzz_xxzzz[k] = -g_y_0_yyzz_xxzzz[k] * ab_z + g_y_0_yyzz_xxzzzz[k];

                g_y_0_yyzzz_xyyyy[k] = -g_y_0_yyzz_xyyyy[k] * ab_z + g_y_0_yyzz_xyyyyz[k];

                g_y_0_yyzzz_xyyyz[k] = -g_y_0_yyzz_xyyyz[k] * ab_z + g_y_0_yyzz_xyyyzz[k];

                g_y_0_yyzzz_xyyzz[k] = -g_y_0_yyzz_xyyzz[k] * ab_z + g_y_0_yyzz_xyyzzz[k];

                g_y_0_yyzzz_xyzzz[k] = -g_y_0_yyzz_xyzzz[k] * ab_z + g_y_0_yyzz_xyzzzz[k];

                g_y_0_yyzzz_xzzzz[k] = -g_y_0_yyzz_xzzzz[k] * ab_z + g_y_0_yyzz_xzzzzz[k];

                g_y_0_yyzzz_yyyyy[k] = -g_y_0_yyzz_yyyyy[k] * ab_z + g_y_0_yyzz_yyyyyz[k];

                g_y_0_yyzzz_yyyyz[k] = -g_y_0_yyzz_yyyyz[k] * ab_z + g_y_0_yyzz_yyyyzz[k];

                g_y_0_yyzzz_yyyzz[k] = -g_y_0_yyzz_yyyzz[k] * ab_z + g_y_0_yyzz_yyyzzz[k];

                g_y_0_yyzzz_yyzzz[k] = -g_y_0_yyzz_yyzzz[k] * ab_z + g_y_0_yyzz_yyzzzz[k];

                g_y_0_yyzzz_yzzzz[k] = -g_y_0_yyzz_yzzzz[k] * ab_z + g_y_0_yyzz_yzzzzz[k];

                g_y_0_yyzzz_zzzzz[k] = -g_y_0_yyzz_zzzzz[k] * ab_z + g_y_0_yyzz_zzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzz, g_y_0_yzzz_zzzzzz, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxxx[k] = -g_y_0_yzzz_xxxxx[k] * ab_z + g_y_0_yzzz_xxxxxz[k];

                g_y_0_yzzzz_xxxxy[k] = -g_y_0_yzzz_xxxxy[k] * ab_z + g_y_0_yzzz_xxxxyz[k];

                g_y_0_yzzzz_xxxxz[k] = -g_y_0_yzzz_xxxxz[k] * ab_z + g_y_0_yzzz_xxxxzz[k];

                g_y_0_yzzzz_xxxyy[k] = -g_y_0_yzzz_xxxyy[k] * ab_z + g_y_0_yzzz_xxxyyz[k];

                g_y_0_yzzzz_xxxyz[k] = -g_y_0_yzzz_xxxyz[k] * ab_z + g_y_0_yzzz_xxxyzz[k];

                g_y_0_yzzzz_xxxzz[k] = -g_y_0_yzzz_xxxzz[k] * ab_z + g_y_0_yzzz_xxxzzz[k];

                g_y_0_yzzzz_xxyyy[k] = -g_y_0_yzzz_xxyyy[k] * ab_z + g_y_0_yzzz_xxyyyz[k];

                g_y_0_yzzzz_xxyyz[k] = -g_y_0_yzzz_xxyyz[k] * ab_z + g_y_0_yzzz_xxyyzz[k];

                g_y_0_yzzzz_xxyzz[k] = -g_y_0_yzzz_xxyzz[k] * ab_z + g_y_0_yzzz_xxyzzz[k];

                g_y_0_yzzzz_xxzzz[k] = -g_y_0_yzzz_xxzzz[k] * ab_z + g_y_0_yzzz_xxzzzz[k];

                g_y_0_yzzzz_xyyyy[k] = -g_y_0_yzzz_xyyyy[k] * ab_z + g_y_0_yzzz_xyyyyz[k];

                g_y_0_yzzzz_xyyyz[k] = -g_y_0_yzzz_xyyyz[k] * ab_z + g_y_0_yzzz_xyyyzz[k];

                g_y_0_yzzzz_xyyzz[k] = -g_y_0_yzzz_xyyzz[k] * ab_z + g_y_0_yzzz_xyyzzz[k];

                g_y_0_yzzzz_xyzzz[k] = -g_y_0_yzzz_xyzzz[k] * ab_z + g_y_0_yzzz_xyzzzz[k];

                g_y_0_yzzzz_xzzzz[k] = -g_y_0_yzzz_xzzzz[k] * ab_z + g_y_0_yzzz_xzzzzz[k];

                g_y_0_yzzzz_yyyyy[k] = -g_y_0_yzzz_yyyyy[k] * ab_z + g_y_0_yzzz_yyyyyz[k];

                g_y_0_yzzzz_yyyyz[k] = -g_y_0_yzzz_yyyyz[k] * ab_z + g_y_0_yzzz_yyyyzz[k];

                g_y_0_yzzzz_yyyzz[k] = -g_y_0_yzzz_yyyzz[k] * ab_z + g_y_0_yzzz_yyyzzz[k];

                g_y_0_yzzzz_yyzzz[k] = -g_y_0_yzzz_yyzzz[k] * ab_z + g_y_0_yzzz_yyzzzz[k];

                g_y_0_yzzzz_yzzzz[k] = -g_y_0_yzzz_yzzzz[k] * ab_z + g_y_0_yzzz_yzzzzz[k];

                g_y_0_yzzzz_zzzzz[k] = -g_y_0_yzzz_zzzzz[k] * ab_z + g_y_0_yzzz_zzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzz, g_y_0_zzzz_zzzzzz, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxxx[k] = -g_y_0_zzzz_xxxxx[k] * ab_z + g_y_0_zzzz_xxxxxz[k];

                g_y_0_zzzzz_xxxxy[k] = -g_y_0_zzzz_xxxxy[k] * ab_z + g_y_0_zzzz_xxxxyz[k];

                g_y_0_zzzzz_xxxxz[k] = -g_y_0_zzzz_xxxxz[k] * ab_z + g_y_0_zzzz_xxxxzz[k];

                g_y_0_zzzzz_xxxyy[k] = -g_y_0_zzzz_xxxyy[k] * ab_z + g_y_0_zzzz_xxxyyz[k];

                g_y_0_zzzzz_xxxyz[k] = -g_y_0_zzzz_xxxyz[k] * ab_z + g_y_0_zzzz_xxxyzz[k];

                g_y_0_zzzzz_xxxzz[k] = -g_y_0_zzzz_xxxzz[k] * ab_z + g_y_0_zzzz_xxxzzz[k];

                g_y_0_zzzzz_xxyyy[k] = -g_y_0_zzzz_xxyyy[k] * ab_z + g_y_0_zzzz_xxyyyz[k];

                g_y_0_zzzzz_xxyyz[k] = -g_y_0_zzzz_xxyyz[k] * ab_z + g_y_0_zzzz_xxyyzz[k];

                g_y_0_zzzzz_xxyzz[k] = -g_y_0_zzzz_xxyzz[k] * ab_z + g_y_0_zzzz_xxyzzz[k];

                g_y_0_zzzzz_xxzzz[k] = -g_y_0_zzzz_xxzzz[k] * ab_z + g_y_0_zzzz_xxzzzz[k];

                g_y_0_zzzzz_xyyyy[k] = -g_y_0_zzzz_xyyyy[k] * ab_z + g_y_0_zzzz_xyyyyz[k];

                g_y_0_zzzzz_xyyyz[k] = -g_y_0_zzzz_xyyyz[k] * ab_z + g_y_0_zzzz_xyyyzz[k];

                g_y_0_zzzzz_xyyzz[k] = -g_y_0_zzzz_xyyzz[k] * ab_z + g_y_0_zzzz_xyyzzz[k];

                g_y_0_zzzzz_xyzzz[k] = -g_y_0_zzzz_xyzzz[k] * ab_z + g_y_0_zzzz_xyzzzz[k];

                g_y_0_zzzzz_xzzzz[k] = -g_y_0_zzzz_xzzzz[k] * ab_z + g_y_0_zzzz_xzzzzz[k];

                g_y_0_zzzzz_yyyyy[k] = -g_y_0_zzzz_yyyyy[k] * ab_z + g_y_0_zzzz_yyyyyz[k];

                g_y_0_zzzzz_yyyyz[k] = -g_y_0_zzzz_yyyyz[k] * ab_z + g_y_0_zzzz_yyyyzz[k];

                g_y_0_zzzzz_yyyzz[k] = -g_y_0_zzzz_yyyzz[k] * ab_z + g_y_0_zzzz_yyyzzz[k];

                g_y_0_zzzzz_yyzzz[k] = -g_y_0_zzzz_yyzzz[k] * ab_z + g_y_0_zzzz_yyzzzz[k];

                g_y_0_zzzzz_yzzzz[k] = -g_y_0_zzzz_yzzzz[k] * ab_z + g_y_0_zzzz_yzzzzz[k];

                g_y_0_zzzzz_zzzzz[k] = -g_y_0_zzzz_zzzzz[k] * ab_z + g_y_0_zzzz_zzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxx_xxxxx, g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxy, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxyy, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxyyy, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xyyyy, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_yyyyy, g_z_0_xxxx_yyyyz, g_z_0_xxxx_yyyzz, g_z_0_xxxx_yyzzz, g_z_0_xxxx_yzzzz, g_z_0_xxxx_zzzzz, g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_yyyyy, g_z_0_xxxxx_yyyyz, g_z_0_xxxxx_yyyzz, g_z_0_xxxxx_yyzzz, g_z_0_xxxxx_yzzzz, g_z_0_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxxx[k] = -g_z_0_xxxx_xxxxx[k] * ab_x + g_z_0_xxxx_xxxxxx[k];

                g_z_0_xxxxx_xxxxy[k] = -g_z_0_xxxx_xxxxy[k] * ab_x + g_z_0_xxxx_xxxxxy[k];

                g_z_0_xxxxx_xxxxz[k] = -g_z_0_xxxx_xxxxz[k] * ab_x + g_z_0_xxxx_xxxxxz[k];

                g_z_0_xxxxx_xxxyy[k] = -g_z_0_xxxx_xxxyy[k] * ab_x + g_z_0_xxxx_xxxxyy[k];

                g_z_0_xxxxx_xxxyz[k] = -g_z_0_xxxx_xxxyz[k] * ab_x + g_z_0_xxxx_xxxxyz[k];

                g_z_0_xxxxx_xxxzz[k] = -g_z_0_xxxx_xxxzz[k] * ab_x + g_z_0_xxxx_xxxxzz[k];

                g_z_0_xxxxx_xxyyy[k] = -g_z_0_xxxx_xxyyy[k] * ab_x + g_z_0_xxxx_xxxyyy[k];

                g_z_0_xxxxx_xxyyz[k] = -g_z_0_xxxx_xxyyz[k] * ab_x + g_z_0_xxxx_xxxyyz[k];

                g_z_0_xxxxx_xxyzz[k] = -g_z_0_xxxx_xxyzz[k] * ab_x + g_z_0_xxxx_xxxyzz[k];

                g_z_0_xxxxx_xxzzz[k] = -g_z_0_xxxx_xxzzz[k] * ab_x + g_z_0_xxxx_xxxzzz[k];

                g_z_0_xxxxx_xyyyy[k] = -g_z_0_xxxx_xyyyy[k] * ab_x + g_z_0_xxxx_xxyyyy[k];

                g_z_0_xxxxx_xyyyz[k] = -g_z_0_xxxx_xyyyz[k] * ab_x + g_z_0_xxxx_xxyyyz[k];

                g_z_0_xxxxx_xyyzz[k] = -g_z_0_xxxx_xyyzz[k] * ab_x + g_z_0_xxxx_xxyyzz[k];

                g_z_0_xxxxx_xyzzz[k] = -g_z_0_xxxx_xyzzz[k] * ab_x + g_z_0_xxxx_xxyzzz[k];

                g_z_0_xxxxx_xzzzz[k] = -g_z_0_xxxx_xzzzz[k] * ab_x + g_z_0_xxxx_xxzzzz[k];

                g_z_0_xxxxx_yyyyy[k] = -g_z_0_xxxx_yyyyy[k] * ab_x + g_z_0_xxxx_xyyyyy[k];

                g_z_0_xxxxx_yyyyz[k] = -g_z_0_xxxx_yyyyz[k] * ab_x + g_z_0_xxxx_xyyyyz[k];

                g_z_0_xxxxx_yyyzz[k] = -g_z_0_xxxx_yyyzz[k] * ab_x + g_z_0_xxxx_xyyyzz[k];

                g_z_0_xxxxx_yyzzz[k] = -g_z_0_xxxx_yyzzz[k] * ab_x + g_z_0_xxxx_xyyzzz[k];

                g_z_0_xxxxx_yzzzz[k] = -g_z_0_xxxx_yzzzz[k] * ab_x + g_z_0_xxxx_xyzzzz[k];

                g_z_0_xxxxx_zzzzz[k] = -g_z_0_xxxx_zzzzz[k] * ab_x + g_z_0_xxxx_xzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_yyyyy, g_z_0_xxxxy_yyyyz, g_z_0_xxxxy_yyyzz, g_z_0_xxxxy_yyzzz, g_z_0_xxxxy_yzzzz, g_z_0_xxxxy_zzzzz, g_z_0_xxxy_xxxxx, g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxy, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxyy, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxyyy, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xyyyy, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_yyyyy, g_z_0_xxxy_yyyyz, g_z_0_xxxy_yyyzz, g_z_0_xxxy_yyzzz, g_z_0_xxxy_yzzzz, g_z_0_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxxx[k] = -g_z_0_xxxy_xxxxx[k] * ab_x + g_z_0_xxxy_xxxxxx[k];

                g_z_0_xxxxy_xxxxy[k] = -g_z_0_xxxy_xxxxy[k] * ab_x + g_z_0_xxxy_xxxxxy[k];

                g_z_0_xxxxy_xxxxz[k] = -g_z_0_xxxy_xxxxz[k] * ab_x + g_z_0_xxxy_xxxxxz[k];

                g_z_0_xxxxy_xxxyy[k] = -g_z_0_xxxy_xxxyy[k] * ab_x + g_z_0_xxxy_xxxxyy[k];

                g_z_0_xxxxy_xxxyz[k] = -g_z_0_xxxy_xxxyz[k] * ab_x + g_z_0_xxxy_xxxxyz[k];

                g_z_0_xxxxy_xxxzz[k] = -g_z_0_xxxy_xxxzz[k] * ab_x + g_z_0_xxxy_xxxxzz[k];

                g_z_0_xxxxy_xxyyy[k] = -g_z_0_xxxy_xxyyy[k] * ab_x + g_z_0_xxxy_xxxyyy[k];

                g_z_0_xxxxy_xxyyz[k] = -g_z_0_xxxy_xxyyz[k] * ab_x + g_z_0_xxxy_xxxyyz[k];

                g_z_0_xxxxy_xxyzz[k] = -g_z_0_xxxy_xxyzz[k] * ab_x + g_z_0_xxxy_xxxyzz[k];

                g_z_0_xxxxy_xxzzz[k] = -g_z_0_xxxy_xxzzz[k] * ab_x + g_z_0_xxxy_xxxzzz[k];

                g_z_0_xxxxy_xyyyy[k] = -g_z_0_xxxy_xyyyy[k] * ab_x + g_z_0_xxxy_xxyyyy[k];

                g_z_0_xxxxy_xyyyz[k] = -g_z_0_xxxy_xyyyz[k] * ab_x + g_z_0_xxxy_xxyyyz[k];

                g_z_0_xxxxy_xyyzz[k] = -g_z_0_xxxy_xyyzz[k] * ab_x + g_z_0_xxxy_xxyyzz[k];

                g_z_0_xxxxy_xyzzz[k] = -g_z_0_xxxy_xyzzz[k] * ab_x + g_z_0_xxxy_xxyzzz[k];

                g_z_0_xxxxy_xzzzz[k] = -g_z_0_xxxy_xzzzz[k] * ab_x + g_z_0_xxxy_xxzzzz[k];

                g_z_0_xxxxy_yyyyy[k] = -g_z_0_xxxy_yyyyy[k] * ab_x + g_z_0_xxxy_xyyyyy[k];

                g_z_0_xxxxy_yyyyz[k] = -g_z_0_xxxy_yyyyz[k] * ab_x + g_z_0_xxxy_xyyyyz[k];

                g_z_0_xxxxy_yyyzz[k] = -g_z_0_xxxy_yyyzz[k] * ab_x + g_z_0_xxxy_xyyyzz[k];

                g_z_0_xxxxy_yyzzz[k] = -g_z_0_xxxy_yyzzz[k] * ab_x + g_z_0_xxxy_xyyzzz[k];

                g_z_0_xxxxy_yzzzz[k] = -g_z_0_xxxy_yzzzz[k] * ab_x + g_z_0_xxxy_xyzzzz[k];

                g_z_0_xxxxy_zzzzz[k] = -g_z_0_xxxy_zzzzz[k] * ab_x + g_z_0_xxxy_xzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_yyyyy, g_z_0_xxxxz_yyyyz, g_z_0_xxxxz_yyyzz, g_z_0_xxxxz_yyzzz, g_z_0_xxxxz_yzzzz, g_z_0_xxxxz_zzzzz, g_z_0_xxxz_xxxxx, g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxy, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxyy, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxyyy, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xyyyy, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_yyyyy, g_z_0_xxxz_yyyyz, g_z_0_xxxz_yyyzz, g_z_0_xxxz_yyzzz, g_z_0_xxxz_yzzzz, g_z_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxxx[k] = -g_z_0_xxxz_xxxxx[k] * ab_x + g_z_0_xxxz_xxxxxx[k];

                g_z_0_xxxxz_xxxxy[k] = -g_z_0_xxxz_xxxxy[k] * ab_x + g_z_0_xxxz_xxxxxy[k];

                g_z_0_xxxxz_xxxxz[k] = -g_z_0_xxxz_xxxxz[k] * ab_x + g_z_0_xxxz_xxxxxz[k];

                g_z_0_xxxxz_xxxyy[k] = -g_z_0_xxxz_xxxyy[k] * ab_x + g_z_0_xxxz_xxxxyy[k];

                g_z_0_xxxxz_xxxyz[k] = -g_z_0_xxxz_xxxyz[k] * ab_x + g_z_0_xxxz_xxxxyz[k];

                g_z_0_xxxxz_xxxzz[k] = -g_z_0_xxxz_xxxzz[k] * ab_x + g_z_0_xxxz_xxxxzz[k];

                g_z_0_xxxxz_xxyyy[k] = -g_z_0_xxxz_xxyyy[k] * ab_x + g_z_0_xxxz_xxxyyy[k];

                g_z_0_xxxxz_xxyyz[k] = -g_z_0_xxxz_xxyyz[k] * ab_x + g_z_0_xxxz_xxxyyz[k];

                g_z_0_xxxxz_xxyzz[k] = -g_z_0_xxxz_xxyzz[k] * ab_x + g_z_0_xxxz_xxxyzz[k];

                g_z_0_xxxxz_xxzzz[k] = -g_z_0_xxxz_xxzzz[k] * ab_x + g_z_0_xxxz_xxxzzz[k];

                g_z_0_xxxxz_xyyyy[k] = -g_z_0_xxxz_xyyyy[k] * ab_x + g_z_0_xxxz_xxyyyy[k];

                g_z_0_xxxxz_xyyyz[k] = -g_z_0_xxxz_xyyyz[k] * ab_x + g_z_0_xxxz_xxyyyz[k];

                g_z_0_xxxxz_xyyzz[k] = -g_z_0_xxxz_xyyzz[k] * ab_x + g_z_0_xxxz_xxyyzz[k];

                g_z_0_xxxxz_xyzzz[k] = -g_z_0_xxxz_xyzzz[k] * ab_x + g_z_0_xxxz_xxyzzz[k];

                g_z_0_xxxxz_xzzzz[k] = -g_z_0_xxxz_xzzzz[k] * ab_x + g_z_0_xxxz_xxzzzz[k];

                g_z_0_xxxxz_yyyyy[k] = -g_z_0_xxxz_yyyyy[k] * ab_x + g_z_0_xxxz_xyyyyy[k];

                g_z_0_xxxxz_yyyyz[k] = -g_z_0_xxxz_yyyyz[k] * ab_x + g_z_0_xxxz_xyyyyz[k];

                g_z_0_xxxxz_yyyzz[k] = -g_z_0_xxxz_yyyzz[k] * ab_x + g_z_0_xxxz_xyyyzz[k];

                g_z_0_xxxxz_yyzzz[k] = -g_z_0_xxxz_yyzzz[k] * ab_x + g_z_0_xxxz_xyyzzz[k];

                g_z_0_xxxxz_yzzzz[k] = -g_z_0_xxxz_yzzzz[k] * ab_x + g_z_0_xxxz_xyzzzz[k];

                g_z_0_xxxxz_zzzzz[k] = -g_z_0_xxxz_zzzzz[k] * ab_x + g_z_0_xxxz_xzzzzz[k];
            }

            /// Set up 945-966 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_yyyyy, g_z_0_xxxyy_yyyyz, g_z_0_xxxyy_yyyzz, g_z_0_xxxyy_yyzzz, g_z_0_xxxyy_yzzzz, g_z_0_xxxyy_zzzzz, g_z_0_xxyy_xxxxx, g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxy, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxyy, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxyyy, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xyyyy, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_yyyyy, g_z_0_xxyy_yyyyz, g_z_0_xxyy_yyyzz, g_z_0_xxyy_yyzzz, g_z_0_xxyy_yzzzz, g_z_0_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxxx[k] = -g_z_0_xxyy_xxxxx[k] * ab_x + g_z_0_xxyy_xxxxxx[k];

                g_z_0_xxxyy_xxxxy[k] = -g_z_0_xxyy_xxxxy[k] * ab_x + g_z_0_xxyy_xxxxxy[k];

                g_z_0_xxxyy_xxxxz[k] = -g_z_0_xxyy_xxxxz[k] * ab_x + g_z_0_xxyy_xxxxxz[k];

                g_z_0_xxxyy_xxxyy[k] = -g_z_0_xxyy_xxxyy[k] * ab_x + g_z_0_xxyy_xxxxyy[k];

                g_z_0_xxxyy_xxxyz[k] = -g_z_0_xxyy_xxxyz[k] * ab_x + g_z_0_xxyy_xxxxyz[k];

                g_z_0_xxxyy_xxxzz[k] = -g_z_0_xxyy_xxxzz[k] * ab_x + g_z_0_xxyy_xxxxzz[k];

                g_z_0_xxxyy_xxyyy[k] = -g_z_0_xxyy_xxyyy[k] * ab_x + g_z_0_xxyy_xxxyyy[k];

                g_z_0_xxxyy_xxyyz[k] = -g_z_0_xxyy_xxyyz[k] * ab_x + g_z_0_xxyy_xxxyyz[k];

                g_z_0_xxxyy_xxyzz[k] = -g_z_0_xxyy_xxyzz[k] * ab_x + g_z_0_xxyy_xxxyzz[k];

                g_z_0_xxxyy_xxzzz[k] = -g_z_0_xxyy_xxzzz[k] * ab_x + g_z_0_xxyy_xxxzzz[k];

                g_z_0_xxxyy_xyyyy[k] = -g_z_0_xxyy_xyyyy[k] * ab_x + g_z_0_xxyy_xxyyyy[k];

                g_z_0_xxxyy_xyyyz[k] = -g_z_0_xxyy_xyyyz[k] * ab_x + g_z_0_xxyy_xxyyyz[k];

                g_z_0_xxxyy_xyyzz[k] = -g_z_0_xxyy_xyyzz[k] * ab_x + g_z_0_xxyy_xxyyzz[k];

                g_z_0_xxxyy_xyzzz[k] = -g_z_0_xxyy_xyzzz[k] * ab_x + g_z_0_xxyy_xxyzzz[k];

                g_z_0_xxxyy_xzzzz[k] = -g_z_0_xxyy_xzzzz[k] * ab_x + g_z_0_xxyy_xxzzzz[k];

                g_z_0_xxxyy_yyyyy[k] = -g_z_0_xxyy_yyyyy[k] * ab_x + g_z_0_xxyy_xyyyyy[k];

                g_z_0_xxxyy_yyyyz[k] = -g_z_0_xxyy_yyyyz[k] * ab_x + g_z_0_xxyy_xyyyyz[k];

                g_z_0_xxxyy_yyyzz[k] = -g_z_0_xxyy_yyyzz[k] * ab_x + g_z_0_xxyy_xyyyzz[k];

                g_z_0_xxxyy_yyzzz[k] = -g_z_0_xxyy_yyzzz[k] * ab_x + g_z_0_xxyy_xyyzzz[k];

                g_z_0_xxxyy_yzzzz[k] = -g_z_0_xxyy_yzzzz[k] * ab_x + g_z_0_xxyy_xyzzzz[k];

                g_z_0_xxxyy_zzzzz[k] = -g_z_0_xxyy_zzzzz[k] * ab_x + g_z_0_xxyy_xzzzzz[k];
            }

            /// Set up 966-987 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_yyyyy, g_z_0_xxxyz_yyyyz, g_z_0_xxxyz_yyyzz, g_z_0_xxxyz_yyzzz, g_z_0_xxxyz_yzzzz, g_z_0_xxxyz_zzzzz, g_z_0_xxyz_xxxxx, g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxy, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxyy, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxyyy, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xyyyy, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_yyyyy, g_z_0_xxyz_yyyyz, g_z_0_xxyz_yyyzz, g_z_0_xxyz_yyzzz, g_z_0_xxyz_yzzzz, g_z_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxxx[k] = -g_z_0_xxyz_xxxxx[k] * ab_x + g_z_0_xxyz_xxxxxx[k];

                g_z_0_xxxyz_xxxxy[k] = -g_z_0_xxyz_xxxxy[k] * ab_x + g_z_0_xxyz_xxxxxy[k];

                g_z_0_xxxyz_xxxxz[k] = -g_z_0_xxyz_xxxxz[k] * ab_x + g_z_0_xxyz_xxxxxz[k];

                g_z_0_xxxyz_xxxyy[k] = -g_z_0_xxyz_xxxyy[k] * ab_x + g_z_0_xxyz_xxxxyy[k];

                g_z_0_xxxyz_xxxyz[k] = -g_z_0_xxyz_xxxyz[k] * ab_x + g_z_0_xxyz_xxxxyz[k];

                g_z_0_xxxyz_xxxzz[k] = -g_z_0_xxyz_xxxzz[k] * ab_x + g_z_0_xxyz_xxxxzz[k];

                g_z_0_xxxyz_xxyyy[k] = -g_z_0_xxyz_xxyyy[k] * ab_x + g_z_0_xxyz_xxxyyy[k];

                g_z_0_xxxyz_xxyyz[k] = -g_z_0_xxyz_xxyyz[k] * ab_x + g_z_0_xxyz_xxxyyz[k];

                g_z_0_xxxyz_xxyzz[k] = -g_z_0_xxyz_xxyzz[k] * ab_x + g_z_0_xxyz_xxxyzz[k];

                g_z_0_xxxyz_xxzzz[k] = -g_z_0_xxyz_xxzzz[k] * ab_x + g_z_0_xxyz_xxxzzz[k];

                g_z_0_xxxyz_xyyyy[k] = -g_z_0_xxyz_xyyyy[k] * ab_x + g_z_0_xxyz_xxyyyy[k];

                g_z_0_xxxyz_xyyyz[k] = -g_z_0_xxyz_xyyyz[k] * ab_x + g_z_0_xxyz_xxyyyz[k];

                g_z_0_xxxyz_xyyzz[k] = -g_z_0_xxyz_xyyzz[k] * ab_x + g_z_0_xxyz_xxyyzz[k];

                g_z_0_xxxyz_xyzzz[k] = -g_z_0_xxyz_xyzzz[k] * ab_x + g_z_0_xxyz_xxyzzz[k];

                g_z_0_xxxyz_xzzzz[k] = -g_z_0_xxyz_xzzzz[k] * ab_x + g_z_0_xxyz_xxzzzz[k];

                g_z_0_xxxyz_yyyyy[k] = -g_z_0_xxyz_yyyyy[k] * ab_x + g_z_0_xxyz_xyyyyy[k];

                g_z_0_xxxyz_yyyyz[k] = -g_z_0_xxyz_yyyyz[k] * ab_x + g_z_0_xxyz_xyyyyz[k];

                g_z_0_xxxyz_yyyzz[k] = -g_z_0_xxyz_yyyzz[k] * ab_x + g_z_0_xxyz_xyyyzz[k];

                g_z_0_xxxyz_yyzzz[k] = -g_z_0_xxyz_yyzzz[k] * ab_x + g_z_0_xxyz_xyyzzz[k];

                g_z_0_xxxyz_yzzzz[k] = -g_z_0_xxyz_yzzzz[k] * ab_x + g_z_0_xxyz_xyzzzz[k];

                g_z_0_xxxyz_zzzzz[k] = -g_z_0_xxyz_zzzzz[k] * ab_x + g_z_0_xxyz_xzzzzz[k];
            }

            /// Set up 987-1008 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_yyyyy, g_z_0_xxxzz_yyyyz, g_z_0_xxxzz_yyyzz, g_z_0_xxxzz_yyzzz, g_z_0_xxxzz_yzzzz, g_z_0_xxxzz_zzzzz, g_z_0_xxzz_xxxxx, g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxy, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxyy, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxyyy, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xyyyy, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_yyyyy, g_z_0_xxzz_yyyyz, g_z_0_xxzz_yyyzz, g_z_0_xxzz_yyzzz, g_z_0_xxzz_yzzzz, g_z_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxxx[k] = -g_z_0_xxzz_xxxxx[k] * ab_x + g_z_0_xxzz_xxxxxx[k];

                g_z_0_xxxzz_xxxxy[k] = -g_z_0_xxzz_xxxxy[k] * ab_x + g_z_0_xxzz_xxxxxy[k];

                g_z_0_xxxzz_xxxxz[k] = -g_z_0_xxzz_xxxxz[k] * ab_x + g_z_0_xxzz_xxxxxz[k];

                g_z_0_xxxzz_xxxyy[k] = -g_z_0_xxzz_xxxyy[k] * ab_x + g_z_0_xxzz_xxxxyy[k];

                g_z_0_xxxzz_xxxyz[k] = -g_z_0_xxzz_xxxyz[k] * ab_x + g_z_0_xxzz_xxxxyz[k];

                g_z_0_xxxzz_xxxzz[k] = -g_z_0_xxzz_xxxzz[k] * ab_x + g_z_0_xxzz_xxxxzz[k];

                g_z_0_xxxzz_xxyyy[k] = -g_z_0_xxzz_xxyyy[k] * ab_x + g_z_0_xxzz_xxxyyy[k];

                g_z_0_xxxzz_xxyyz[k] = -g_z_0_xxzz_xxyyz[k] * ab_x + g_z_0_xxzz_xxxyyz[k];

                g_z_0_xxxzz_xxyzz[k] = -g_z_0_xxzz_xxyzz[k] * ab_x + g_z_0_xxzz_xxxyzz[k];

                g_z_0_xxxzz_xxzzz[k] = -g_z_0_xxzz_xxzzz[k] * ab_x + g_z_0_xxzz_xxxzzz[k];

                g_z_0_xxxzz_xyyyy[k] = -g_z_0_xxzz_xyyyy[k] * ab_x + g_z_0_xxzz_xxyyyy[k];

                g_z_0_xxxzz_xyyyz[k] = -g_z_0_xxzz_xyyyz[k] * ab_x + g_z_0_xxzz_xxyyyz[k];

                g_z_0_xxxzz_xyyzz[k] = -g_z_0_xxzz_xyyzz[k] * ab_x + g_z_0_xxzz_xxyyzz[k];

                g_z_0_xxxzz_xyzzz[k] = -g_z_0_xxzz_xyzzz[k] * ab_x + g_z_0_xxzz_xxyzzz[k];

                g_z_0_xxxzz_xzzzz[k] = -g_z_0_xxzz_xzzzz[k] * ab_x + g_z_0_xxzz_xxzzzz[k];

                g_z_0_xxxzz_yyyyy[k] = -g_z_0_xxzz_yyyyy[k] * ab_x + g_z_0_xxzz_xyyyyy[k];

                g_z_0_xxxzz_yyyyz[k] = -g_z_0_xxzz_yyyyz[k] * ab_x + g_z_0_xxzz_xyyyyz[k];

                g_z_0_xxxzz_yyyzz[k] = -g_z_0_xxzz_yyyzz[k] * ab_x + g_z_0_xxzz_xyyyzz[k];

                g_z_0_xxxzz_yyzzz[k] = -g_z_0_xxzz_yyzzz[k] * ab_x + g_z_0_xxzz_xyyzzz[k];

                g_z_0_xxxzz_yzzzz[k] = -g_z_0_xxzz_yzzzz[k] * ab_x + g_z_0_xxzz_xyzzzz[k];

                g_z_0_xxxzz_zzzzz[k] = -g_z_0_xxzz_zzzzz[k] * ab_x + g_z_0_xxzz_xzzzzz[k];
            }

            /// Set up 1008-1029 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_yyyyy, g_z_0_xxyyy_yyyyz, g_z_0_xxyyy_yyyzz, g_z_0_xxyyy_yyzzz, g_z_0_xxyyy_yzzzz, g_z_0_xxyyy_zzzzz, g_z_0_xyyy_xxxxx, g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxy, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxyy, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxyyy, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xyyyy, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_yyyyy, g_z_0_xyyy_yyyyz, g_z_0_xyyy_yyyzz, g_z_0_xyyy_yyzzz, g_z_0_xyyy_yzzzz, g_z_0_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxxx[k] = -g_z_0_xyyy_xxxxx[k] * ab_x + g_z_0_xyyy_xxxxxx[k];

                g_z_0_xxyyy_xxxxy[k] = -g_z_0_xyyy_xxxxy[k] * ab_x + g_z_0_xyyy_xxxxxy[k];

                g_z_0_xxyyy_xxxxz[k] = -g_z_0_xyyy_xxxxz[k] * ab_x + g_z_0_xyyy_xxxxxz[k];

                g_z_0_xxyyy_xxxyy[k] = -g_z_0_xyyy_xxxyy[k] * ab_x + g_z_0_xyyy_xxxxyy[k];

                g_z_0_xxyyy_xxxyz[k] = -g_z_0_xyyy_xxxyz[k] * ab_x + g_z_0_xyyy_xxxxyz[k];

                g_z_0_xxyyy_xxxzz[k] = -g_z_0_xyyy_xxxzz[k] * ab_x + g_z_0_xyyy_xxxxzz[k];

                g_z_0_xxyyy_xxyyy[k] = -g_z_0_xyyy_xxyyy[k] * ab_x + g_z_0_xyyy_xxxyyy[k];

                g_z_0_xxyyy_xxyyz[k] = -g_z_0_xyyy_xxyyz[k] * ab_x + g_z_0_xyyy_xxxyyz[k];

                g_z_0_xxyyy_xxyzz[k] = -g_z_0_xyyy_xxyzz[k] * ab_x + g_z_0_xyyy_xxxyzz[k];

                g_z_0_xxyyy_xxzzz[k] = -g_z_0_xyyy_xxzzz[k] * ab_x + g_z_0_xyyy_xxxzzz[k];

                g_z_0_xxyyy_xyyyy[k] = -g_z_0_xyyy_xyyyy[k] * ab_x + g_z_0_xyyy_xxyyyy[k];

                g_z_0_xxyyy_xyyyz[k] = -g_z_0_xyyy_xyyyz[k] * ab_x + g_z_0_xyyy_xxyyyz[k];

                g_z_0_xxyyy_xyyzz[k] = -g_z_0_xyyy_xyyzz[k] * ab_x + g_z_0_xyyy_xxyyzz[k];

                g_z_0_xxyyy_xyzzz[k] = -g_z_0_xyyy_xyzzz[k] * ab_x + g_z_0_xyyy_xxyzzz[k];

                g_z_0_xxyyy_xzzzz[k] = -g_z_0_xyyy_xzzzz[k] * ab_x + g_z_0_xyyy_xxzzzz[k];

                g_z_0_xxyyy_yyyyy[k] = -g_z_0_xyyy_yyyyy[k] * ab_x + g_z_0_xyyy_xyyyyy[k];

                g_z_0_xxyyy_yyyyz[k] = -g_z_0_xyyy_yyyyz[k] * ab_x + g_z_0_xyyy_xyyyyz[k];

                g_z_0_xxyyy_yyyzz[k] = -g_z_0_xyyy_yyyzz[k] * ab_x + g_z_0_xyyy_xyyyzz[k];

                g_z_0_xxyyy_yyzzz[k] = -g_z_0_xyyy_yyzzz[k] * ab_x + g_z_0_xyyy_xyyzzz[k];

                g_z_0_xxyyy_yzzzz[k] = -g_z_0_xyyy_yzzzz[k] * ab_x + g_z_0_xyyy_xyzzzz[k];

                g_z_0_xxyyy_zzzzz[k] = -g_z_0_xyyy_zzzzz[k] * ab_x + g_z_0_xyyy_xzzzzz[k];
            }

            /// Set up 1029-1050 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_yyyyy, g_z_0_xxyyz_yyyyz, g_z_0_xxyyz_yyyzz, g_z_0_xxyyz_yyzzz, g_z_0_xxyyz_yzzzz, g_z_0_xxyyz_zzzzz, g_z_0_xyyz_xxxxx, g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxy, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxyy, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxyyy, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xyyyy, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_yyyyy, g_z_0_xyyz_yyyyz, g_z_0_xyyz_yyyzz, g_z_0_xyyz_yyzzz, g_z_0_xyyz_yzzzz, g_z_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxxx[k] = -g_z_0_xyyz_xxxxx[k] * ab_x + g_z_0_xyyz_xxxxxx[k];

                g_z_0_xxyyz_xxxxy[k] = -g_z_0_xyyz_xxxxy[k] * ab_x + g_z_0_xyyz_xxxxxy[k];

                g_z_0_xxyyz_xxxxz[k] = -g_z_0_xyyz_xxxxz[k] * ab_x + g_z_0_xyyz_xxxxxz[k];

                g_z_0_xxyyz_xxxyy[k] = -g_z_0_xyyz_xxxyy[k] * ab_x + g_z_0_xyyz_xxxxyy[k];

                g_z_0_xxyyz_xxxyz[k] = -g_z_0_xyyz_xxxyz[k] * ab_x + g_z_0_xyyz_xxxxyz[k];

                g_z_0_xxyyz_xxxzz[k] = -g_z_0_xyyz_xxxzz[k] * ab_x + g_z_0_xyyz_xxxxzz[k];

                g_z_0_xxyyz_xxyyy[k] = -g_z_0_xyyz_xxyyy[k] * ab_x + g_z_0_xyyz_xxxyyy[k];

                g_z_0_xxyyz_xxyyz[k] = -g_z_0_xyyz_xxyyz[k] * ab_x + g_z_0_xyyz_xxxyyz[k];

                g_z_0_xxyyz_xxyzz[k] = -g_z_0_xyyz_xxyzz[k] * ab_x + g_z_0_xyyz_xxxyzz[k];

                g_z_0_xxyyz_xxzzz[k] = -g_z_0_xyyz_xxzzz[k] * ab_x + g_z_0_xyyz_xxxzzz[k];

                g_z_0_xxyyz_xyyyy[k] = -g_z_0_xyyz_xyyyy[k] * ab_x + g_z_0_xyyz_xxyyyy[k];

                g_z_0_xxyyz_xyyyz[k] = -g_z_0_xyyz_xyyyz[k] * ab_x + g_z_0_xyyz_xxyyyz[k];

                g_z_0_xxyyz_xyyzz[k] = -g_z_0_xyyz_xyyzz[k] * ab_x + g_z_0_xyyz_xxyyzz[k];

                g_z_0_xxyyz_xyzzz[k] = -g_z_0_xyyz_xyzzz[k] * ab_x + g_z_0_xyyz_xxyzzz[k];

                g_z_0_xxyyz_xzzzz[k] = -g_z_0_xyyz_xzzzz[k] * ab_x + g_z_0_xyyz_xxzzzz[k];

                g_z_0_xxyyz_yyyyy[k] = -g_z_0_xyyz_yyyyy[k] * ab_x + g_z_0_xyyz_xyyyyy[k];

                g_z_0_xxyyz_yyyyz[k] = -g_z_0_xyyz_yyyyz[k] * ab_x + g_z_0_xyyz_xyyyyz[k];

                g_z_0_xxyyz_yyyzz[k] = -g_z_0_xyyz_yyyzz[k] * ab_x + g_z_0_xyyz_xyyyzz[k];

                g_z_0_xxyyz_yyzzz[k] = -g_z_0_xyyz_yyzzz[k] * ab_x + g_z_0_xyyz_xyyzzz[k];

                g_z_0_xxyyz_yzzzz[k] = -g_z_0_xyyz_yzzzz[k] * ab_x + g_z_0_xyyz_xyzzzz[k];

                g_z_0_xxyyz_zzzzz[k] = -g_z_0_xyyz_zzzzz[k] * ab_x + g_z_0_xyyz_xzzzzz[k];
            }

            /// Set up 1050-1071 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_yyyyy, g_z_0_xxyzz_yyyyz, g_z_0_xxyzz_yyyzz, g_z_0_xxyzz_yyzzz, g_z_0_xxyzz_yzzzz, g_z_0_xxyzz_zzzzz, g_z_0_xyzz_xxxxx, g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxy, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxyy, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxyyy, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xyyyy, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_yyyyy, g_z_0_xyzz_yyyyz, g_z_0_xyzz_yyyzz, g_z_0_xyzz_yyzzz, g_z_0_xyzz_yzzzz, g_z_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxxx[k] = -g_z_0_xyzz_xxxxx[k] * ab_x + g_z_0_xyzz_xxxxxx[k];

                g_z_0_xxyzz_xxxxy[k] = -g_z_0_xyzz_xxxxy[k] * ab_x + g_z_0_xyzz_xxxxxy[k];

                g_z_0_xxyzz_xxxxz[k] = -g_z_0_xyzz_xxxxz[k] * ab_x + g_z_0_xyzz_xxxxxz[k];

                g_z_0_xxyzz_xxxyy[k] = -g_z_0_xyzz_xxxyy[k] * ab_x + g_z_0_xyzz_xxxxyy[k];

                g_z_0_xxyzz_xxxyz[k] = -g_z_0_xyzz_xxxyz[k] * ab_x + g_z_0_xyzz_xxxxyz[k];

                g_z_0_xxyzz_xxxzz[k] = -g_z_0_xyzz_xxxzz[k] * ab_x + g_z_0_xyzz_xxxxzz[k];

                g_z_0_xxyzz_xxyyy[k] = -g_z_0_xyzz_xxyyy[k] * ab_x + g_z_0_xyzz_xxxyyy[k];

                g_z_0_xxyzz_xxyyz[k] = -g_z_0_xyzz_xxyyz[k] * ab_x + g_z_0_xyzz_xxxyyz[k];

                g_z_0_xxyzz_xxyzz[k] = -g_z_0_xyzz_xxyzz[k] * ab_x + g_z_0_xyzz_xxxyzz[k];

                g_z_0_xxyzz_xxzzz[k] = -g_z_0_xyzz_xxzzz[k] * ab_x + g_z_0_xyzz_xxxzzz[k];

                g_z_0_xxyzz_xyyyy[k] = -g_z_0_xyzz_xyyyy[k] * ab_x + g_z_0_xyzz_xxyyyy[k];

                g_z_0_xxyzz_xyyyz[k] = -g_z_0_xyzz_xyyyz[k] * ab_x + g_z_0_xyzz_xxyyyz[k];

                g_z_0_xxyzz_xyyzz[k] = -g_z_0_xyzz_xyyzz[k] * ab_x + g_z_0_xyzz_xxyyzz[k];

                g_z_0_xxyzz_xyzzz[k] = -g_z_0_xyzz_xyzzz[k] * ab_x + g_z_0_xyzz_xxyzzz[k];

                g_z_0_xxyzz_xzzzz[k] = -g_z_0_xyzz_xzzzz[k] * ab_x + g_z_0_xyzz_xxzzzz[k];

                g_z_0_xxyzz_yyyyy[k] = -g_z_0_xyzz_yyyyy[k] * ab_x + g_z_0_xyzz_xyyyyy[k];

                g_z_0_xxyzz_yyyyz[k] = -g_z_0_xyzz_yyyyz[k] * ab_x + g_z_0_xyzz_xyyyyz[k];

                g_z_0_xxyzz_yyyzz[k] = -g_z_0_xyzz_yyyzz[k] * ab_x + g_z_0_xyzz_xyyyzz[k];

                g_z_0_xxyzz_yyzzz[k] = -g_z_0_xyzz_yyzzz[k] * ab_x + g_z_0_xyzz_xyyzzz[k];

                g_z_0_xxyzz_yzzzz[k] = -g_z_0_xyzz_yzzzz[k] * ab_x + g_z_0_xyzz_xyzzzz[k];

                g_z_0_xxyzz_zzzzz[k] = -g_z_0_xyzz_zzzzz[k] * ab_x + g_z_0_xyzz_xzzzzz[k];
            }

            /// Set up 1071-1092 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_yyyyy, g_z_0_xxzzz_yyyyz, g_z_0_xxzzz_yyyzz, g_z_0_xxzzz_yyzzz, g_z_0_xxzzz_yzzzz, g_z_0_xxzzz_zzzzz, g_z_0_xzzz_xxxxx, g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxy, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxyy, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxyyy, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xyyyy, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_yyyyy, g_z_0_xzzz_yyyyz, g_z_0_xzzz_yyyzz, g_z_0_xzzz_yyzzz, g_z_0_xzzz_yzzzz, g_z_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxxx[k] = -g_z_0_xzzz_xxxxx[k] * ab_x + g_z_0_xzzz_xxxxxx[k];

                g_z_0_xxzzz_xxxxy[k] = -g_z_0_xzzz_xxxxy[k] * ab_x + g_z_0_xzzz_xxxxxy[k];

                g_z_0_xxzzz_xxxxz[k] = -g_z_0_xzzz_xxxxz[k] * ab_x + g_z_0_xzzz_xxxxxz[k];

                g_z_0_xxzzz_xxxyy[k] = -g_z_0_xzzz_xxxyy[k] * ab_x + g_z_0_xzzz_xxxxyy[k];

                g_z_0_xxzzz_xxxyz[k] = -g_z_0_xzzz_xxxyz[k] * ab_x + g_z_0_xzzz_xxxxyz[k];

                g_z_0_xxzzz_xxxzz[k] = -g_z_0_xzzz_xxxzz[k] * ab_x + g_z_0_xzzz_xxxxzz[k];

                g_z_0_xxzzz_xxyyy[k] = -g_z_0_xzzz_xxyyy[k] * ab_x + g_z_0_xzzz_xxxyyy[k];

                g_z_0_xxzzz_xxyyz[k] = -g_z_0_xzzz_xxyyz[k] * ab_x + g_z_0_xzzz_xxxyyz[k];

                g_z_0_xxzzz_xxyzz[k] = -g_z_0_xzzz_xxyzz[k] * ab_x + g_z_0_xzzz_xxxyzz[k];

                g_z_0_xxzzz_xxzzz[k] = -g_z_0_xzzz_xxzzz[k] * ab_x + g_z_0_xzzz_xxxzzz[k];

                g_z_0_xxzzz_xyyyy[k] = -g_z_0_xzzz_xyyyy[k] * ab_x + g_z_0_xzzz_xxyyyy[k];

                g_z_0_xxzzz_xyyyz[k] = -g_z_0_xzzz_xyyyz[k] * ab_x + g_z_0_xzzz_xxyyyz[k];

                g_z_0_xxzzz_xyyzz[k] = -g_z_0_xzzz_xyyzz[k] * ab_x + g_z_0_xzzz_xxyyzz[k];

                g_z_0_xxzzz_xyzzz[k] = -g_z_0_xzzz_xyzzz[k] * ab_x + g_z_0_xzzz_xxyzzz[k];

                g_z_0_xxzzz_xzzzz[k] = -g_z_0_xzzz_xzzzz[k] * ab_x + g_z_0_xzzz_xxzzzz[k];

                g_z_0_xxzzz_yyyyy[k] = -g_z_0_xzzz_yyyyy[k] * ab_x + g_z_0_xzzz_xyyyyy[k];

                g_z_0_xxzzz_yyyyz[k] = -g_z_0_xzzz_yyyyz[k] * ab_x + g_z_0_xzzz_xyyyyz[k];

                g_z_0_xxzzz_yyyzz[k] = -g_z_0_xzzz_yyyzz[k] * ab_x + g_z_0_xzzz_xyyyzz[k];

                g_z_0_xxzzz_yyzzz[k] = -g_z_0_xzzz_yyzzz[k] * ab_x + g_z_0_xzzz_xyyzzz[k];

                g_z_0_xxzzz_yzzzz[k] = -g_z_0_xzzz_yzzzz[k] * ab_x + g_z_0_xzzz_xyzzzz[k];

                g_z_0_xxzzz_zzzzz[k] = -g_z_0_xzzz_zzzzz[k] * ab_x + g_z_0_xzzz_xzzzzz[k];
            }

            /// Set up 1092-1113 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_yyyyy, g_z_0_xyyyy_yyyyz, g_z_0_xyyyy_yyyzz, g_z_0_xyyyy_yyzzz, g_z_0_xyyyy_yzzzz, g_z_0_xyyyy_zzzzz, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxxx[k] = -g_z_0_yyyy_xxxxx[k] * ab_x + g_z_0_yyyy_xxxxxx[k];

                g_z_0_xyyyy_xxxxy[k] = -g_z_0_yyyy_xxxxy[k] * ab_x + g_z_0_yyyy_xxxxxy[k];

                g_z_0_xyyyy_xxxxz[k] = -g_z_0_yyyy_xxxxz[k] * ab_x + g_z_0_yyyy_xxxxxz[k];

                g_z_0_xyyyy_xxxyy[k] = -g_z_0_yyyy_xxxyy[k] * ab_x + g_z_0_yyyy_xxxxyy[k];

                g_z_0_xyyyy_xxxyz[k] = -g_z_0_yyyy_xxxyz[k] * ab_x + g_z_0_yyyy_xxxxyz[k];

                g_z_0_xyyyy_xxxzz[k] = -g_z_0_yyyy_xxxzz[k] * ab_x + g_z_0_yyyy_xxxxzz[k];

                g_z_0_xyyyy_xxyyy[k] = -g_z_0_yyyy_xxyyy[k] * ab_x + g_z_0_yyyy_xxxyyy[k];

                g_z_0_xyyyy_xxyyz[k] = -g_z_0_yyyy_xxyyz[k] * ab_x + g_z_0_yyyy_xxxyyz[k];

                g_z_0_xyyyy_xxyzz[k] = -g_z_0_yyyy_xxyzz[k] * ab_x + g_z_0_yyyy_xxxyzz[k];

                g_z_0_xyyyy_xxzzz[k] = -g_z_0_yyyy_xxzzz[k] * ab_x + g_z_0_yyyy_xxxzzz[k];

                g_z_0_xyyyy_xyyyy[k] = -g_z_0_yyyy_xyyyy[k] * ab_x + g_z_0_yyyy_xxyyyy[k];

                g_z_0_xyyyy_xyyyz[k] = -g_z_0_yyyy_xyyyz[k] * ab_x + g_z_0_yyyy_xxyyyz[k];

                g_z_0_xyyyy_xyyzz[k] = -g_z_0_yyyy_xyyzz[k] * ab_x + g_z_0_yyyy_xxyyzz[k];

                g_z_0_xyyyy_xyzzz[k] = -g_z_0_yyyy_xyzzz[k] * ab_x + g_z_0_yyyy_xxyzzz[k];

                g_z_0_xyyyy_xzzzz[k] = -g_z_0_yyyy_xzzzz[k] * ab_x + g_z_0_yyyy_xxzzzz[k];

                g_z_0_xyyyy_yyyyy[k] = -g_z_0_yyyy_yyyyy[k] * ab_x + g_z_0_yyyy_xyyyyy[k];

                g_z_0_xyyyy_yyyyz[k] = -g_z_0_yyyy_yyyyz[k] * ab_x + g_z_0_yyyy_xyyyyz[k];

                g_z_0_xyyyy_yyyzz[k] = -g_z_0_yyyy_yyyzz[k] * ab_x + g_z_0_yyyy_xyyyzz[k];

                g_z_0_xyyyy_yyzzz[k] = -g_z_0_yyyy_yyzzz[k] * ab_x + g_z_0_yyyy_xyyzzz[k];

                g_z_0_xyyyy_yzzzz[k] = -g_z_0_yyyy_yzzzz[k] * ab_x + g_z_0_yyyy_xyzzzz[k];

                g_z_0_xyyyy_zzzzz[k] = -g_z_0_yyyy_zzzzz[k] * ab_x + g_z_0_yyyy_xzzzzz[k];
            }

            /// Set up 1113-1134 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_yyyyy, g_z_0_xyyyz_yyyyz, g_z_0_xyyyz_yyyzz, g_z_0_xyyyz_yyzzz, g_z_0_xyyyz_yzzzz, g_z_0_xyyyz_zzzzz, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxxx[k] = -g_z_0_yyyz_xxxxx[k] * ab_x + g_z_0_yyyz_xxxxxx[k];

                g_z_0_xyyyz_xxxxy[k] = -g_z_0_yyyz_xxxxy[k] * ab_x + g_z_0_yyyz_xxxxxy[k];

                g_z_0_xyyyz_xxxxz[k] = -g_z_0_yyyz_xxxxz[k] * ab_x + g_z_0_yyyz_xxxxxz[k];

                g_z_0_xyyyz_xxxyy[k] = -g_z_0_yyyz_xxxyy[k] * ab_x + g_z_0_yyyz_xxxxyy[k];

                g_z_0_xyyyz_xxxyz[k] = -g_z_0_yyyz_xxxyz[k] * ab_x + g_z_0_yyyz_xxxxyz[k];

                g_z_0_xyyyz_xxxzz[k] = -g_z_0_yyyz_xxxzz[k] * ab_x + g_z_0_yyyz_xxxxzz[k];

                g_z_0_xyyyz_xxyyy[k] = -g_z_0_yyyz_xxyyy[k] * ab_x + g_z_0_yyyz_xxxyyy[k];

                g_z_0_xyyyz_xxyyz[k] = -g_z_0_yyyz_xxyyz[k] * ab_x + g_z_0_yyyz_xxxyyz[k];

                g_z_0_xyyyz_xxyzz[k] = -g_z_0_yyyz_xxyzz[k] * ab_x + g_z_0_yyyz_xxxyzz[k];

                g_z_0_xyyyz_xxzzz[k] = -g_z_0_yyyz_xxzzz[k] * ab_x + g_z_0_yyyz_xxxzzz[k];

                g_z_0_xyyyz_xyyyy[k] = -g_z_0_yyyz_xyyyy[k] * ab_x + g_z_0_yyyz_xxyyyy[k];

                g_z_0_xyyyz_xyyyz[k] = -g_z_0_yyyz_xyyyz[k] * ab_x + g_z_0_yyyz_xxyyyz[k];

                g_z_0_xyyyz_xyyzz[k] = -g_z_0_yyyz_xyyzz[k] * ab_x + g_z_0_yyyz_xxyyzz[k];

                g_z_0_xyyyz_xyzzz[k] = -g_z_0_yyyz_xyzzz[k] * ab_x + g_z_0_yyyz_xxyzzz[k];

                g_z_0_xyyyz_xzzzz[k] = -g_z_0_yyyz_xzzzz[k] * ab_x + g_z_0_yyyz_xxzzzz[k];

                g_z_0_xyyyz_yyyyy[k] = -g_z_0_yyyz_yyyyy[k] * ab_x + g_z_0_yyyz_xyyyyy[k];

                g_z_0_xyyyz_yyyyz[k] = -g_z_0_yyyz_yyyyz[k] * ab_x + g_z_0_yyyz_xyyyyz[k];

                g_z_0_xyyyz_yyyzz[k] = -g_z_0_yyyz_yyyzz[k] * ab_x + g_z_0_yyyz_xyyyzz[k];

                g_z_0_xyyyz_yyzzz[k] = -g_z_0_yyyz_yyzzz[k] * ab_x + g_z_0_yyyz_xyyzzz[k];

                g_z_0_xyyyz_yzzzz[k] = -g_z_0_yyyz_yzzzz[k] * ab_x + g_z_0_yyyz_xyzzzz[k];

                g_z_0_xyyyz_zzzzz[k] = -g_z_0_yyyz_zzzzz[k] * ab_x + g_z_0_yyyz_xzzzzz[k];
            }

            /// Set up 1134-1155 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_yyyyy, g_z_0_xyyzz_yyyyz, g_z_0_xyyzz_yyyzz, g_z_0_xyyzz_yyzzz, g_z_0_xyyzz_yzzzz, g_z_0_xyyzz_zzzzz, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxxx[k] = -g_z_0_yyzz_xxxxx[k] * ab_x + g_z_0_yyzz_xxxxxx[k];

                g_z_0_xyyzz_xxxxy[k] = -g_z_0_yyzz_xxxxy[k] * ab_x + g_z_0_yyzz_xxxxxy[k];

                g_z_0_xyyzz_xxxxz[k] = -g_z_0_yyzz_xxxxz[k] * ab_x + g_z_0_yyzz_xxxxxz[k];

                g_z_0_xyyzz_xxxyy[k] = -g_z_0_yyzz_xxxyy[k] * ab_x + g_z_0_yyzz_xxxxyy[k];

                g_z_0_xyyzz_xxxyz[k] = -g_z_0_yyzz_xxxyz[k] * ab_x + g_z_0_yyzz_xxxxyz[k];

                g_z_0_xyyzz_xxxzz[k] = -g_z_0_yyzz_xxxzz[k] * ab_x + g_z_0_yyzz_xxxxzz[k];

                g_z_0_xyyzz_xxyyy[k] = -g_z_0_yyzz_xxyyy[k] * ab_x + g_z_0_yyzz_xxxyyy[k];

                g_z_0_xyyzz_xxyyz[k] = -g_z_0_yyzz_xxyyz[k] * ab_x + g_z_0_yyzz_xxxyyz[k];

                g_z_0_xyyzz_xxyzz[k] = -g_z_0_yyzz_xxyzz[k] * ab_x + g_z_0_yyzz_xxxyzz[k];

                g_z_0_xyyzz_xxzzz[k] = -g_z_0_yyzz_xxzzz[k] * ab_x + g_z_0_yyzz_xxxzzz[k];

                g_z_0_xyyzz_xyyyy[k] = -g_z_0_yyzz_xyyyy[k] * ab_x + g_z_0_yyzz_xxyyyy[k];

                g_z_0_xyyzz_xyyyz[k] = -g_z_0_yyzz_xyyyz[k] * ab_x + g_z_0_yyzz_xxyyyz[k];

                g_z_0_xyyzz_xyyzz[k] = -g_z_0_yyzz_xyyzz[k] * ab_x + g_z_0_yyzz_xxyyzz[k];

                g_z_0_xyyzz_xyzzz[k] = -g_z_0_yyzz_xyzzz[k] * ab_x + g_z_0_yyzz_xxyzzz[k];

                g_z_0_xyyzz_xzzzz[k] = -g_z_0_yyzz_xzzzz[k] * ab_x + g_z_0_yyzz_xxzzzz[k];

                g_z_0_xyyzz_yyyyy[k] = -g_z_0_yyzz_yyyyy[k] * ab_x + g_z_0_yyzz_xyyyyy[k];

                g_z_0_xyyzz_yyyyz[k] = -g_z_0_yyzz_yyyyz[k] * ab_x + g_z_0_yyzz_xyyyyz[k];

                g_z_0_xyyzz_yyyzz[k] = -g_z_0_yyzz_yyyzz[k] * ab_x + g_z_0_yyzz_xyyyzz[k];

                g_z_0_xyyzz_yyzzz[k] = -g_z_0_yyzz_yyzzz[k] * ab_x + g_z_0_yyzz_xyyzzz[k];

                g_z_0_xyyzz_yzzzz[k] = -g_z_0_yyzz_yzzzz[k] * ab_x + g_z_0_yyzz_xyzzzz[k];

                g_z_0_xyyzz_zzzzz[k] = -g_z_0_yyzz_zzzzz[k] * ab_x + g_z_0_yyzz_xzzzzz[k];
            }

            /// Set up 1155-1176 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_yyyyy, g_z_0_xyzzz_yyyyz, g_z_0_xyzzz_yyyzz, g_z_0_xyzzz_yyzzz, g_z_0_xyzzz_yzzzz, g_z_0_xyzzz_zzzzz, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxxx[k] = -g_z_0_yzzz_xxxxx[k] * ab_x + g_z_0_yzzz_xxxxxx[k];

                g_z_0_xyzzz_xxxxy[k] = -g_z_0_yzzz_xxxxy[k] * ab_x + g_z_0_yzzz_xxxxxy[k];

                g_z_0_xyzzz_xxxxz[k] = -g_z_0_yzzz_xxxxz[k] * ab_x + g_z_0_yzzz_xxxxxz[k];

                g_z_0_xyzzz_xxxyy[k] = -g_z_0_yzzz_xxxyy[k] * ab_x + g_z_0_yzzz_xxxxyy[k];

                g_z_0_xyzzz_xxxyz[k] = -g_z_0_yzzz_xxxyz[k] * ab_x + g_z_0_yzzz_xxxxyz[k];

                g_z_0_xyzzz_xxxzz[k] = -g_z_0_yzzz_xxxzz[k] * ab_x + g_z_0_yzzz_xxxxzz[k];

                g_z_0_xyzzz_xxyyy[k] = -g_z_0_yzzz_xxyyy[k] * ab_x + g_z_0_yzzz_xxxyyy[k];

                g_z_0_xyzzz_xxyyz[k] = -g_z_0_yzzz_xxyyz[k] * ab_x + g_z_0_yzzz_xxxyyz[k];

                g_z_0_xyzzz_xxyzz[k] = -g_z_0_yzzz_xxyzz[k] * ab_x + g_z_0_yzzz_xxxyzz[k];

                g_z_0_xyzzz_xxzzz[k] = -g_z_0_yzzz_xxzzz[k] * ab_x + g_z_0_yzzz_xxxzzz[k];

                g_z_0_xyzzz_xyyyy[k] = -g_z_0_yzzz_xyyyy[k] * ab_x + g_z_0_yzzz_xxyyyy[k];

                g_z_0_xyzzz_xyyyz[k] = -g_z_0_yzzz_xyyyz[k] * ab_x + g_z_0_yzzz_xxyyyz[k];

                g_z_0_xyzzz_xyyzz[k] = -g_z_0_yzzz_xyyzz[k] * ab_x + g_z_0_yzzz_xxyyzz[k];

                g_z_0_xyzzz_xyzzz[k] = -g_z_0_yzzz_xyzzz[k] * ab_x + g_z_0_yzzz_xxyzzz[k];

                g_z_0_xyzzz_xzzzz[k] = -g_z_0_yzzz_xzzzz[k] * ab_x + g_z_0_yzzz_xxzzzz[k];

                g_z_0_xyzzz_yyyyy[k] = -g_z_0_yzzz_yyyyy[k] * ab_x + g_z_0_yzzz_xyyyyy[k];

                g_z_0_xyzzz_yyyyz[k] = -g_z_0_yzzz_yyyyz[k] * ab_x + g_z_0_yzzz_xyyyyz[k];

                g_z_0_xyzzz_yyyzz[k] = -g_z_0_yzzz_yyyzz[k] * ab_x + g_z_0_yzzz_xyyyzz[k];

                g_z_0_xyzzz_yyzzz[k] = -g_z_0_yzzz_yyzzz[k] * ab_x + g_z_0_yzzz_xyyzzz[k];

                g_z_0_xyzzz_yzzzz[k] = -g_z_0_yzzz_yzzzz[k] * ab_x + g_z_0_yzzz_xyzzzz[k];

                g_z_0_xyzzz_zzzzz[k] = -g_z_0_yzzz_zzzzz[k] * ab_x + g_z_0_yzzz_xzzzzz[k];
            }

            /// Set up 1176-1197 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_yyyyy, g_z_0_xzzzz_yyyyz, g_z_0_xzzzz_yyyzz, g_z_0_xzzzz_yyzzz, g_z_0_xzzzz_yzzzz, g_z_0_xzzzz_zzzzz, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxxx[k] = -g_z_0_zzzz_xxxxx[k] * ab_x + g_z_0_zzzz_xxxxxx[k];

                g_z_0_xzzzz_xxxxy[k] = -g_z_0_zzzz_xxxxy[k] * ab_x + g_z_0_zzzz_xxxxxy[k];

                g_z_0_xzzzz_xxxxz[k] = -g_z_0_zzzz_xxxxz[k] * ab_x + g_z_0_zzzz_xxxxxz[k];

                g_z_0_xzzzz_xxxyy[k] = -g_z_0_zzzz_xxxyy[k] * ab_x + g_z_0_zzzz_xxxxyy[k];

                g_z_0_xzzzz_xxxyz[k] = -g_z_0_zzzz_xxxyz[k] * ab_x + g_z_0_zzzz_xxxxyz[k];

                g_z_0_xzzzz_xxxzz[k] = -g_z_0_zzzz_xxxzz[k] * ab_x + g_z_0_zzzz_xxxxzz[k];

                g_z_0_xzzzz_xxyyy[k] = -g_z_0_zzzz_xxyyy[k] * ab_x + g_z_0_zzzz_xxxyyy[k];

                g_z_0_xzzzz_xxyyz[k] = -g_z_0_zzzz_xxyyz[k] * ab_x + g_z_0_zzzz_xxxyyz[k];

                g_z_0_xzzzz_xxyzz[k] = -g_z_0_zzzz_xxyzz[k] * ab_x + g_z_0_zzzz_xxxyzz[k];

                g_z_0_xzzzz_xxzzz[k] = -g_z_0_zzzz_xxzzz[k] * ab_x + g_z_0_zzzz_xxxzzz[k];

                g_z_0_xzzzz_xyyyy[k] = -g_z_0_zzzz_xyyyy[k] * ab_x + g_z_0_zzzz_xxyyyy[k];

                g_z_0_xzzzz_xyyyz[k] = -g_z_0_zzzz_xyyyz[k] * ab_x + g_z_0_zzzz_xxyyyz[k];

                g_z_0_xzzzz_xyyzz[k] = -g_z_0_zzzz_xyyzz[k] * ab_x + g_z_0_zzzz_xxyyzz[k];

                g_z_0_xzzzz_xyzzz[k] = -g_z_0_zzzz_xyzzz[k] * ab_x + g_z_0_zzzz_xxyzzz[k];

                g_z_0_xzzzz_xzzzz[k] = -g_z_0_zzzz_xzzzz[k] * ab_x + g_z_0_zzzz_xxzzzz[k];

                g_z_0_xzzzz_yyyyy[k] = -g_z_0_zzzz_yyyyy[k] * ab_x + g_z_0_zzzz_xyyyyy[k];

                g_z_0_xzzzz_yyyyz[k] = -g_z_0_zzzz_yyyyz[k] * ab_x + g_z_0_zzzz_xyyyyz[k];

                g_z_0_xzzzz_yyyzz[k] = -g_z_0_zzzz_yyyzz[k] * ab_x + g_z_0_zzzz_xyyyzz[k];

                g_z_0_xzzzz_yyzzz[k] = -g_z_0_zzzz_yyzzz[k] * ab_x + g_z_0_zzzz_xyyzzz[k];

                g_z_0_xzzzz_yzzzz[k] = -g_z_0_zzzz_yzzzz[k] * ab_x + g_z_0_zzzz_xyzzzz[k];

                g_z_0_xzzzz_zzzzz[k] = -g_z_0_zzzz_zzzzz[k] * ab_x + g_z_0_zzzz_xzzzzz[k];
            }

            /// Set up 1197-1218 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzz, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxxx[k] = -g_z_0_yyyy_xxxxx[k] * ab_y + g_z_0_yyyy_xxxxxy[k];

                g_z_0_yyyyy_xxxxy[k] = -g_z_0_yyyy_xxxxy[k] * ab_y + g_z_0_yyyy_xxxxyy[k];

                g_z_0_yyyyy_xxxxz[k] = -g_z_0_yyyy_xxxxz[k] * ab_y + g_z_0_yyyy_xxxxyz[k];

                g_z_0_yyyyy_xxxyy[k] = -g_z_0_yyyy_xxxyy[k] * ab_y + g_z_0_yyyy_xxxyyy[k];

                g_z_0_yyyyy_xxxyz[k] = -g_z_0_yyyy_xxxyz[k] * ab_y + g_z_0_yyyy_xxxyyz[k];

                g_z_0_yyyyy_xxxzz[k] = -g_z_0_yyyy_xxxzz[k] * ab_y + g_z_0_yyyy_xxxyzz[k];

                g_z_0_yyyyy_xxyyy[k] = -g_z_0_yyyy_xxyyy[k] * ab_y + g_z_0_yyyy_xxyyyy[k];

                g_z_0_yyyyy_xxyyz[k] = -g_z_0_yyyy_xxyyz[k] * ab_y + g_z_0_yyyy_xxyyyz[k];

                g_z_0_yyyyy_xxyzz[k] = -g_z_0_yyyy_xxyzz[k] * ab_y + g_z_0_yyyy_xxyyzz[k];

                g_z_0_yyyyy_xxzzz[k] = -g_z_0_yyyy_xxzzz[k] * ab_y + g_z_0_yyyy_xxyzzz[k];

                g_z_0_yyyyy_xyyyy[k] = -g_z_0_yyyy_xyyyy[k] * ab_y + g_z_0_yyyy_xyyyyy[k];

                g_z_0_yyyyy_xyyyz[k] = -g_z_0_yyyy_xyyyz[k] * ab_y + g_z_0_yyyy_xyyyyz[k];

                g_z_0_yyyyy_xyyzz[k] = -g_z_0_yyyy_xyyzz[k] * ab_y + g_z_0_yyyy_xyyyzz[k];

                g_z_0_yyyyy_xyzzz[k] = -g_z_0_yyyy_xyzzz[k] * ab_y + g_z_0_yyyy_xyyzzz[k];

                g_z_0_yyyyy_xzzzz[k] = -g_z_0_yyyy_xzzzz[k] * ab_y + g_z_0_yyyy_xyzzzz[k];

                g_z_0_yyyyy_yyyyy[k] = -g_z_0_yyyy_yyyyy[k] * ab_y + g_z_0_yyyy_yyyyyy[k];

                g_z_0_yyyyy_yyyyz[k] = -g_z_0_yyyy_yyyyz[k] * ab_y + g_z_0_yyyy_yyyyyz[k];

                g_z_0_yyyyy_yyyzz[k] = -g_z_0_yyyy_yyyzz[k] * ab_y + g_z_0_yyyy_yyyyzz[k];

                g_z_0_yyyyy_yyzzz[k] = -g_z_0_yyyy_yyzzz[k] * ab_y + g_z_0_yyyy_yyyzzz[k];

                g_z_0_yyyyy_yzzzz[k] = -g_z_0_yyyy_yzzzz[k] * ab_y + g_z_0_yyyy_yyzzzz[k];

                g_z_0_yyyyy_zzzzz[k] = -g_z_0_yyyy_zzzzz[k] * ab_y + g_z_0_yyyy_yzzzzz[k];
            }

            /// Set up 1218-1239 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzzz, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxxx[k] = -g_z_0_yyyz_xxxxx[k] * ab_y + g_z_0_yyyz_xxxxxy[k];

                g_z_0_yyyyz_xxxxy[k] = -g_z_0_yyyz_xxxxy[k] * ab_y + g_z_0_yyyz_xxxxyy[k];

                g_z_0_yyyyz_xxxxz[k] = -g_z_0_yyyz_xxxxz[k] * ab_y + g_z_0_yyyz_xxxxyz[k];

                g_z_0_yyyyz_xxxyy[k] = -g_z_0_yyyz_xxxyy[k] * ab_y + g_z_0_yyyz_xxxyyy[k];

                g_z_0_yyyyz_xxxyz[k] = -g_z_0_yyyz_xxxyz[k] * ab_y + g_z_0_yyyz_xxxyyz[k];

                g_z_0_yyyyz_xxxzz[k] = -g_z_0_yyyz_xxxzz[k] * ab_y + g_z_0_yyyz_xxxyzz[k];

                g_z_0_yyyyz_xxyyy[k] = -g_z_0_yyyz_xxyyy[k] * ab_y + g_z_0_yyyz_xxyyyy[k];

                g_z_0_yyyyz_xxyyz[k] = -g_z_0_yyyz_xxyyz[k] * ab_y + g_z_0_yyyz_xxyyyz[k];

                g_z_0_yyyyz_xxyzz[k] = -g_z_0_yyyz_xxyzz[k] * ab_y + g_z_0_yyyz_xxyyzz[k];

                g_z_0_yyyyz_xxzzz[k] = -g_z_0_yyyz_xxzzz[k] * ab_y + g_z_0_yyyz_xxyzzz[k];

                g_z_0_yyyyz_xyyyy[k] = -g_z_0_yyyz_xyyyy[k] * ab_y + g_z_0_yyyz_xyyyyy[k];

                g_z_0_yyyyz_xyyyz[k] = -g_z_0_yyyz_xyyyz[k] * ab_y + g_z_0_yyyz_xyyyyz[k];

                g_z_0_yyyyz_xyyzz[k] = -g_z_0_yyyz_xyyzz[k] * ab_y + g_z_0_yyyz_xyyyzz[k];

                g_z_0_yyyyz_xyzzz[k] = -g_z_0_yyyz_xyzzz[k] * ab_y + g_z_0_yyyz_xyyzzz[k];

                g_z_0_yyyyz_xzzzz[k] = -g_z_0_yyyz_xzzzz[k] * ab_y + g_z_0_yyyz_xyzzzz[k];

                g_z_0_yyyyz_yyyyy[k] = -g_z_0_yyyz_yyyyy[k] * ab_y + g_z_0_yyyz_yyyyyy[k];

                g_z_0_yyyyz_yyyyz[k] = -g_z_0_yyyz_yyyyz[k] * ab_y + g_z_0_yyyz_yyyyyz[k];

                g_z_0_yyyyz_yyyzz[k] = -g_z_0_yyyz_yyyzz[k] * ab_y + g_z_0_yyyz_yyyyzz[k];

                g_z_0_yyyyz_yyzzz[k] = -g_z_0_yyyz_yyzzz[k] * ab_y + g_z_0_yyyz_yyyzzz[k];

                g_z_0_yyyyz_yzzzz[k] = -g_z_0_yyyz_yzzzz[k] * ab_y + g_z_0_yyyz_yyzzzz[k];

                g_z_0_yyyyz_zzzzz[k] = -g_z_0_yyyz_zzzzz[k] * ab_y + g_z_0_yyyz_yzzzzz[k];
            }

            /// Set up 1239-1260 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzzz, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxxx[k] = -g_z_0_yyzz_xxxxx[k] * ab_y + g_z_0_yyzz_xxxxxy[k];

                g_z_0_yyyzz_xxxxy[k] = -g_z_0_yyzz_xxxxy[k] * ab_y + g_z_0_yyzz_xxxxyy[k];

                g_z_0_yyyzz_xxxxz[k] = -g_z_0_yyzz_xxxxz[k] * ab_y + g_z_0_yyzz_xxxxyz[k];

                g_z_0_yyyzz_xxxyy[k] = -g_z_0_yyzz_xxxyy[k] * ab_y + g_z_0_yyzz_xxxyyy[k];

                g_z_0_yyyzz_xxxyz[k] = -g_z_0_yyzz_xxxyz[k] * ab_y + g_z_0_yyzz_xxxyyz[k];

                g_z_0_yyyzz_xxxzz[k] = -g_z_0_yyzz_xxxzz[k] * ab_y + g_z_0_yyzz_xxxyzz[k];

                g_z_0_yyyzz_xxyyy[k] = -g_z_0_yyzz_xxyyy[k] * ab_y + g_z_0_yyzz_xxyyyy[k];

                g_z_0_yyyzz_xxyyz[k] = -g_z_0_yyzz_xxyyz[k] * ab_y + g_z_0_yyzz_xxyyyz[k];

                g_z_0_yyyzz_xxyzz[k] = -g_z_0_yyzz_xxyzz[k] * ab_y + g_z_0_yyzz_xxyyzz[k];

                g_z_0_yyyzz_xxzzz[k] = -g_z_0_yyzz_xxzzz[k] * ab_y + g_z_0_yyzz_xxyzzz[k];

                g_z_0_yyyzz_xyyyy[k] = -g_z_0_yyzz_xyyyy[k] * ab_y + g_z_0_yyzz_xyyyyy[k];

                g_z_0_yyyzz_xyyyz[k] = -g_z_0_yyzz_xyyyz[k] * ab_y + g_z_0_yyzz_xyyyyz[k];

                g_z_0_yyyzz_xyyzz[k] = -g_z_0_yyzz_xyyzz[k] * ab_y + g_z_0_yyzz_xyyyzz[k];

                g_z_0_yyyzz_xyzzz[k] = -g_z_0_yyzz_xyzzz[k] * ab_y + g_z_0_yyzz_xyyzzz[k];

                g_z_0_yyyzz_xzzzz[k] = -g_z_0_yyzz_xzzzz[k] * ab_y + g_z_0_yyzz_xyzzzz[k];

                g_z_0_yyyzz_yyyyy[k] = -g_z_0_yyzz_yyyyy[k] * ab_y + g_z_0_yyzz_yyyyyy[k];

                g_z_0_yyyzz_yyyyz[k] = -g_z_0_yyzz_yyyyz[k] * ab_y + g_z_0_yyzz_yyyyyz[k];

                g_z_0_yyyzz_yyyzz[k] = -g_z_0_yyzz_yyyzz[k] * ab_y + g_z_0_yyzz_yyyyzz[k];

                g_z_0_yyyzz_yyzzz[k] = -g_z_0_yyzz_yyzzz[k] * ab_y + g_z_0_yyzz_yyyzzz[k];

                g_z_0_yyyzz_yzzzz[k] = -g_z_0_yyzz_yzzzz[k] * ab_y + g_z_0_yyzz_yyzzzz[k];

                g_z_0_yyyzz_zzzzz[k] = -g_z_0_yyzz_zzzzz[k] * ab_y + g_z_0_yyzz_yzzzzz[k];
            }

            /// Set up 1260-1281 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzzz, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxxx[k] = -g_z_0_yzzz_xxxxx[k] * ab_y + g_z_0_yzzz_xxxxxy[k];

                g_z_0_yyzzz_xxxxy[k] = -g_z_0_yzzz_xxxxy[k] * ab_y + g_z_0_yzzz_xxxxyy[k];

                g_z_0_yyzzz_xxxxz[k] = -g_z_0_yzzz_xxxxz[k] * ab_y + g_z_0_yzzz_xxxxyz[k];

                g_z_0_yyzzz_xxxyy[k] = -g_z_0_yzzz_xxxyy[k] * ab_y + g_z_0_yzzz_xxxyyy[k];

                g_z_0_yyzzz_xxxyz[k] = -g_z_0_yzzz_xxxyz[k] * ab_y + g_z_0_yzzz_xxxyyz[k];

                g_z_0_yyzzz_xxxzz[k] = -g_z_0_yzzz_xxxzz[k] * ab_y + g_z_0_yzzz_xxxyzz[k];

                g_z_0_yyzzz_xxyyy[k] = -g_z_0_yzzz_xxyyy[k] * ab_y + g_z_0_yzzz_xxyyyy[k];

                g_z_0_yyzzz_xxyyz[k] = -g_z_0_yzzz_xxyyz[k] * ab_y + g_z_0_yzzz_xxyyyz[k];

                g_z_0_yyzzz_xxyzz[k] = -g_z_0_yzzz_xxyzz[k] * ab_y + g_z_0_yzzz_xxyyzz[k];

                g_z_0_yyzzz_xxzzz[k] = -g_z_0_yzzz_xxzzz[k] * ab_y + g_z_0_yzzz_xxyzzz[k];

                g_z_0_yyzzz_xyyyy[k] = -g_z_0_yzzz_xyyyy[k] * ab_y + g_z_0_yzzz_xyyyyy[k];

                g_z_0_yyzzz_xyyyz[k] = -g_z_0_yzzz_xyyyz[k] * ab_y + g_z_0_yzzz_xyyyyz[k];

                g_z_0_yyzzz_xyyzz[k] = -g_z_0_yzzz_xyyzz[k] * ab_y + g_z_0_yzzz_xyyyzz[k];

                g_z_0_yyzzz_xyzzz[k] = -g_z_0_yzzz_xyzzz[k] * ab_y + g_z_0_yzzz_xyyzzz[k];

                g_z_0_yyzzz_xzzzz[k] = -g_z_0_yzzz_xzzzz[k] * ab_y + g_z_0_yzzz_xyzzzz[k];

                g_z_0_yyzzz_yyyyy[k] = -g_z_0_yzzz_yyyyy[k] * ab_y + g_z_0_yzzz_yyyyyy[k];

                g_z_0_yyzzz_yyyyz[k] = -g_z_0_yzzz_yyyyz[k] * ab_y + g_z_0_yzzz_yyyyyz[k];

                g_z_0_yyzzz_yyyzz[k] = -g_z_0_yzzz_yyyzz[k] * ab_y + g_z_0_yzzz_yyyyzz[k];

                g_z_0_yyzzz_yyzzz[k] = -g_z_0_yzzz_yyzzz[k] * ab_y + g_z_0_yzzz_yyyzzz[k];

                g_z_0_yyzzz_yzzzz[k] = -g_z_0_yzzz_yzzzz[k] * ab_y + g_z_0_yzzz_yyzzzz[k];

                g_z_0_yyzzz_zzzzz[k] = -g_z_0_yzzz_zzzzz[k] * ab_y + g_z_0_yzzz_yzzzzz[k];
            }

            /// Set up 1281-1302 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzzz, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxxx[k] = -g_z_0_zzzz_xxxxx[k] * ab_y + g_z_0_zzzz_xxxxxy[k];

                g_z_0_yzzzz_xxxxy[k] = -g_z_0_zzzz_xxxxy[k] * ab_y + g_z_0_zzzz_xxxxyy[k];

                g_z_0_yzzzz_xxxxz[k] = -g_z_0_zzzz_xxxxz[k] * ab_y + g_z_0_zzzz_xxxxyz[k];

                g_z_0_yzzzz_xxxyy[k] = -g_z_0_zzzz_xxxyy[k] * ab_y + g_z_0_zzzz_xxxyyy[k];

                g_z_0_yzzzz_xxxyz[k] = -g_z_0_zzzz_xxxyz[k] * ab_y + g_z_0_zzzz_xxxyyz[k];

                g_z_0_yzzzz_xxxzz[k] = -g_z_0_zzzz_xxxzz[k] * ab_y + g_z_0_zzzz_xxxyzz[k];

                g_z_0_yzzzz_xxyyy[k] = -g_z_0_zzzz_xxyyy[k] * ab_y + g_z_0_zzzz_xxyyyy[k];

                g_z_0_yzzzz_xxyyz[k] = -g_z_0_zzzz_xxyyz[k] * ab_y + g_z_0_zzzz_xxyyyz[k];

                g_z_0_yzzzz_xxyzz[k] = -g_z_0_zzzz_xxyzz[k] * ab_y + g_z_0_zzzz_xxyyzz[k];

                g_z_0_yzzzz_xxzzz[k] = -g_z_0_zzzz_xxzzz[k] * ab_y + g_z_0_zzzz_xxyzzz[k];

                g_z_0_yzzzz_xyyyy[k] = -g_z_0_zzzz_xyyyy[k] * ab_y + g_z_0_zzzz_xyyyyy[k];

                g_z_0_yzzzz_xyyyz[k] = -g_z_0_zzzz_xyyyz[k] * ab_y + g_z_0_zzzz_xyyyyz[k];

                g_z_0_yzzzz_xyyzz[k] = -g_z_0_zzzz_xyyzz[k] * ab_y + g_z_0_zzzz_xyyyzz[k];

                g_z_0_yzzzz_xyzzz[k] = -g_z_0_zzzz_xyzzz[k] * ab_y + g_z_0_zzzz_xyyzzz[k];

                g_z_0_yzzzz_xzzzz[k] = -g_z_0_zzzz_xzzzz[k] * ab_y + g_z_0_zzzz_xyzzzz[k];

                g_z_0_yzzzz_yyyyy[k] = -g_z_0_zzzz_yyyyy[k] * ab_y + g_z_0_zzzz_yyyyyy[k];

                g_z_0_yzzzz_yyyyz[k] = -g_z_0_zzzz_yyyyz[k] * ab_y + g_z_0_zzzz_yyyyyz[k];

                g_z_0_yzzzz_yyyzz[k] = -g_z_0_zzzz_yyyzz[k] * ab_y + g_z_0_zzzz_yyyyzz[k];

                g_z_0_yzzzz_yyzzz[k] = -g_z_0_zzzz_yyzzz[k] * ab_y + g_z_0_zzzz_yyyzzz[k];

                g_z_0_yzzzz_yzzzz[k] = -g_z_0_zzzz_yzzzz[k] * ab_y + g_z_0_zzzz_yyzzzz[k];

                g_z_0_yzzzz_zzzzz[k] = -g_z_0_zzzz_zzzzz[k] * ab_y + g_z_0_zzzz_yzzzzz[k];
            }

            /// Set up 1302-1323 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzz, g_z_0_zzzz_zzzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzzz, g_zzzz_xxxxx, g_zzzz_xxxxy, g_zzzz_xxxxz, g_zzzz_xxxyy, g_zzzz_xxxyz, g_zzzz_xxxzz, g_zzzz_xxyyy, g_zzzz_xxyyz, g_zzzz_xxyzz, g_zzzz_xxzzz, g_zzzz_xyyyy, g_zzzz_xyyyz, g_zzzz_xyyzz, g_zzzz_xyzzz, g_zzzz_xzzzz, g_zzzz_yyyyy, g_zzzz_yyyyz, g_zzzz_yyyzz, g_zzzz_yyzzz, g_zzzz_yzzzz, g_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxxx[k] = -g_zzzz_xxxxx[k] - g_z_0_zzzz_xxxxx[k] * ab_z + g_z_0_zzzz_xxxxxz[k];

                g_z_0_zzzzz_xxxxy[k] = -g_zzzz_xxxxy[k] - g_z_0_zzzz_xxxxy[k] * ab_z + g_z_0_zzzz_xxxxyz[k];

                g_z_0_zzzzz_xxxxz[k] = -g_zzzz_xxxxz[k] - g_z_0_zzzz_xxxxz[k] * ab_z + g_z_0_zzzz_xxxxzz[k];

                g_z_0_zzzzz_xxxyy[k] = -g_zzzz_xxxyy[k] - g_z_0_zzzz_xxxyy[k] * ab_z + g_z_0_zzzz_xxxyyz[k];

                g_z_0_zzzzz_xxxyz[k] = -g_zzzz_xxxyz[k] - g_z_0_zzzz_xxxyz[k] * ab_z + g_z_0_zzzz_xxxyzz[k];

                g_z_0_zzzzz_xxxzz[k] = -g_zzzz_xxxzz[k] - g_z_0_zzzz_xxxzz[k] * ab_z + g_z_0_zzzz_xxxzzz[k];

                g_z_0_zzzzz_xxyyy[k] = -g_zzzz_xxyyy[k] - g_z_0_zzzz_xxyyy[k] * ab_z + g_z_0_zzzz_xxyyyz[k];

                g_z_0_zzzzz_xxyyz[k] = -g_zzzz_xxyyz[k] - g_z_0_zzzz_xxyyz[k] * ab_z + g_z_0_zzzz_xxyyzz[k];

                g_z_0_zzzzz_xxyzz[k] = -g_zzzz_xxyzz[k] - g_z_0_zzzz_xxyzz[k] * ab_z + g_z_0_zzzz_xxyzzz[k];

                g_z_0_zzzzz_xxzzz[k] = -g_zzzz_xxzzz[k] - g_z_0_zzzz_xxzzz[k] * ab_z + g_z_0_zzzz_xxzzzz[k];

                g_z_0_zzzzz_xyyyy[k] = -g_zzzz_xyyyy[k] - g_z_0_zzzz_xyyyy[k] * ab_z + g_z_0_zzzz_xyyyyz[k];

                g_z_0_zzzzz_xyyyz[k] = -g_zzzz_xyyyz[k] - g_z_0_zzzz_xyyyz[k] * ab_z + g_z_0_zzzz_xyyyzz[k];

                g_z_0_zzzzz_xyyzz[k] = -g_zzzz_xyyzz[k] - g_z_0_zzzz_xyyzz[k] * ab_z + g_z_0_zzzz_xyyzzz[k];

                g_z_0_zzzzz_xyzzz[k] = -g_zzzz_xyzzz[k] - g_z_0_zzzz_xyzzz[k] * ab_z + g_z_0_zzzz_xyzzzz[k];

                g_z_0_zzzzz_xzzzz[k] = -g_zzzz_xzzzz[k] - g_z_0_zzzz_xzzzz[k] * ab_z + g_z_0_zzzz_xzzzzz[k];

                g_z_0_zzzzz_yyyyy[k] = -g_zzzz_yyyyy[k] - g_z_0_zzzz_yyyyy[k] * ab_z + g_z_0_zzzz_yyyyyz[k];

                g_z_0_zzzzz_yyyyz[k] = -g_zzzz_yyyyz[k] - g_z_0_zzzz_yyyyz[k] * ab_z + g_z_0_zzzz_yyyyzz[k];

                g_z_0_zzzzz_yyyzz[k] = -g_zzzz_yyyzz[k] - g_z_0_zzzz_yyyzz[k] * ab_z + g_z_0_zzzz_yyyzzz[k];

                g_z_0_zzzzz_yyzzz[k] = -g_zzzz_yyzzz[k] - g_z_0_zzzz_yyzzz[k] * ab_z + g_z_0_zzzz_yyzzzz[k];

                g_z_0_zzzzz_yzzzz[k] = -g_zzzz_yzzzz[k] - g_z_0_zzzz_yzzzz[k] * ab_z + g_z_0_zzzz_yzzzzz[k];

                g_z_0_zzzzz_zzzzz[k] = -g_zzzz_zzzzz[k] - g_z_0_zzzz_zzzzz[k] * ab_z + g_z_0_zzzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

