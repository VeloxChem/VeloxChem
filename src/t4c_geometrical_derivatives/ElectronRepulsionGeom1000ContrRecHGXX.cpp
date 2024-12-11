#include "ElectronRepulsionGeom1000ContrRecHGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_hgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_hgxx,
                                            const size_t idx_ggxx,
                                            const size_t idx_geom_10_ggxx,
                                            const size_t idx_geom_10_ghxx,
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
            /// Set up components of auxilary buffer : GGSS

            const auto gg_off = idx_ggxx + i * dcomps + j;

            auto g_xxxx_xxxx = cbuffer.data(gg_off + 0 * ccomps * dcomps);

            auto g_xxxx_xxxy = cbuffer.data(gg_off + 1 * ccomps * dcomps);

            auto g_xxxx_xxxz = cbuffer.data(gg_off + 2 * ccomps * dcomps);

            auto g_xxxx_xxyy = cbuffer.data(gg_off + 3 * ccomps * dcomps);

            auto g_xxxx_xxyz = cbuffer.data(gg_off + 4 * ccomps * dcomps);

            auto g_xxxx_xxzz = cbuffer.data(gg_off + 5 * ccomps * dcomps);

            auto g_xxxx_xyyy = cbuffer.data(gg_off + 6 * ccomps * dcomps);

            auto g_xxxx_xyyz = cbuffer.data(gg_off + 7 * ccomps * dcomps);

            auto g_xxxx_xyzz = cbuffer.data(gg_off + 8 * ccomps * dcomps);

            auto g_xxxx_xzzz = cbuffer.data(gg_off + 9 * ccomps * dcomps);

            auto g_xxxx_yyyy = cbuffer.data(gg_off + 10 * ccomps * dcomps);

            auto g_xxxx_yyyz = cbuffer.data(gg_off + 11 * ccomps * dcomps);

            auto g_xxxx_yyzz = cbuffer.data(gg_off + 12 * ccomps * dcomps);

            auto g_xxxx_yzzz = cbuffer.data(gg_off + 13 * ccomps * dcomps);

            auto g_xxxx_zzzz = cbuffer.data(gg_off + 14 * ccomps * dcomps);

            auto g_yyyy_xxxx = cbuffer.data(gg_off + 150 * ccomps * dcomps);

            auto g_yyyy_xxxy = cbuffer.data(gg_off + 151 * ccomps * dcomps);

            auto g_yyyy_xxxz = cbuffer.data(gg_off + 152 * ccomps * dcomps);

            auto g_yyyy_xxyy = cbuffer.data(gg_off + 153 * ccomps * dcomps);

            auto g_yyyy_xxyz = cbuffer.data(gg_off + 154 * ccomps * dcomps);

            auto g_yyyy_xxzz = cbuffer.data(gg_off + 155 * ccomps * dcomps);

            auto g_yyyy_xyyy = cbuffer.data(gg_off + 156 * ccomps * dcomps);

            auto g_yyyy_xyyz = cbuffer.data(gg_off + 157 * ccomps * dcomps);

            auto g_yyyy_xyzz = cbuffer.data(gg_off + 158 * ccomps * dcomps);

            auto g_yyyy_xzzz = cbuffer.data(gg_off + 159 * ccomps * dcomps);

            auto g_yyyy_yyyy = cbuffer.data(gg_off + 160 * ccomps * dcomps);

            auto g_yyyy_yyyz = cbuffer.data(gg_off + 161 * ccomps * dcomps);

            auto g_yyyy_yyzz = cbuffer.data(gg_off + 162 * ccomps * dcomps);

            auto g_yyyy_yzzz = cbuffer.data(gg_off + 163 * ccomps * dcomps);

            auto g_yyyy_zzzz = cbuffer.data(gg_off + 164 * ccomps * dcomps);

            auto g_zzzz_xxxx = cbuffer.data(gg_off + 210 * ccomps * dcomps);

            auto g_zzzz_xxxy = cbuffer.data(gg_off + 211 * ccomps * dcomps);

            auto g_zzzz_xxxz = cbuffer.data(gg_off + 212 * ccomps * dcomps);

            auto g_zzzz_xxyy = cbuffer.data(gg_off + 213 * ccomps * dcomps);

            auto g_zzzz_xxyz = cbuffer.data(gg_off + 214 * ccomps * dcomps);

            auto g_zzzz_xxzz = cbuffer.data(gg_off + 215 * ccomps * dcomps);

            auto g_zzzz_xyyy = cbuffer.data(gg_off + 216 * ccomps * dcomps);

            auto g_zzzz_xyyz = cbuffer.data(gg_off + 217 * ccomps * dcomps);

            auto g_zzzz_xyzz = cbuffer.data(gg_off + 218 * ccomps * dcomps);

            auto g_zzzz_xzzz = cbuffer.data(gg_off + 219 * ccomps * dcomps);

            auto g_zzzz_yyyy = cbuffer.data(gg_off + 220 * ccomps * dcomps);

            auto g_zzzz_yyyz = cbuffer.data(gg_off + 221 * ccomps * dcomps);

            auto g_zzzz_yyzz = cbuffer.data(gg_off + 222 * ccomps * dcomps);

            auto g_zzzz_yzzz = cbuffer.data(gg_off + 223 * ccomps * dcomps);

            auto g_zzzz_zzzz = cbuffer.data(gg_off + 224 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GGSS

            const auto gg_geom_10_off = idx_geom_10_ggxx + i * dcomps + j;

            auto g_x_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 503 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 509 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 524 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 539 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 659 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 674 * ccomps * dcomps);

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

            auto g_x_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 40 * ccomps * dcomps);

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

            auto g_x_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 103 * ccomps * dcomps);

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

            auto g_x_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 187 * ccomps * dcomps);

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

            auto g_x_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 292 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_hgxx

            const auto hg_geom_10_off = idx_geom_10_hgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzzz, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzzz, g_xxxx_xxxx, g_xxxx_xxxy, g_xxxx_xxxz, g_xxxx_xxyy, g_xxxx_xxyz, g_xxxx_xxzz, g_xxxx_xyyy, g_xxxx_xyyz, g_xxxx_xyzz, g_xxxx_xzzz, g_xxxx_yyyy, g_xxxx_yyyz, g_xxxx_yyzz, g_xxxx_yzzz, g_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxx[k] = -g_xxxx_xxxx[k] - g_x_0_xxxx_xxxx[k] * ab_x + g_x_0_xxxx_xxxxx[k];

                g_x_0_xxxxx_xxxy[k] = -g_xxxx_xxxy[k] - g_x_0_xxxx_xxxy[k] * ab_x + g_x_0_xxxx_xxxxy[k];

                g_x_0_xxxxx_xxxz[k] = -g_xxxx_xxxz[k] - g_x_0_xxxx_xxxz[k] * ab_x + g_x_0_xxxx_xxxxz[k];

                g_x_0_xxxxx_xxyy[k] = -g_xxxx_xxyy[k] - g_x_0_xxxx_xxyy[k] * ab_x + g_x_0_xxxx_xxxyy[k];

                g_x_0_xxxxx_xxyz[k] = -g_xxxx_xxyz[k] - g_x_0_xxxx_xxyz[k] * ab_x + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxx_xxzz[k] = -g_xxxx_xxzz[k] - g_x_0_xxxx_xxzz[k] * ab_x + g_x_0_xxxx_xxxzz[k];

                g_x_0_xxxxx_xyyy[k] = -g_xxxx_xyyy[k] - g_x_0_xxxx_xyyy[k] * ab_x + g_x_0_xxxx_xxyyy[k];

                g_x_0_xxxxx_xyyz[k] = -g_xxxx_xyyz[k] - g_x_0_xxxx_xyyz[k] * ab_x + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxx_xyzz[k] = -g_xxxx_xyzz[k] - g_x_0_xxxx_xyzz[k] * ab_x + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxx_xzzz[k] = -g_xxxx_xzzz[k] - g_x_0_xxxx_xzzz[k] * ab_x + g_x_0_xxxx_xxzzz[k];

                g_x_0_xxxxx_yyyy[k] = -g_xxxx_yyyy[k] - g_x_0_xxxx_yyyy[k] * ab_x + g_x_0_xxxx_xyyyy[k];

                g_x_0_xxxxx_yyyz[k] = -g_xxxx_yyyz[k] - g_x_0_xxxx_yyyz[k] * ab_x + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxx_yyzz[k] = -g_xxxx_yyzz[k] - g_x_0_xxxx_yyzz[k] * ab_x + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxx_yzzz[k] = -g_xxxx_yzzz[k] - g_x_0_xxxx_yzzz[k] * ab_x + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxx_zzzz[k] = -g_xxxx_zzzz[k] - g_x_0_xxxx_zzzz[k] * ab_x + g_x_0_xxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzz, g_x_0_xxxxy_xxxx, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxxz, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxzz, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xzzz, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxx[k] = -g_x_0_xxxx_xxxx[k] * ab_y + g_x_0_xxxx_xxxxy[k];

                g_x_0_xxxxy_xxxy[k] = -g_x_0_xxxx_xxxy[k] * ab_y + g_x_0_xxxx_xxxyy[k];

                g_x_0_xxxxy_xxxz[k] = -g_x_0_xxxx_xxxz[k] * ab_y + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxy_xxyy[k] = -g_x_0_xxxx_xxyy[k] * ab_y + g_x_0_xxxx_xxyyy[k];

                g_x_0_xxxxy_xxyz[k] = -g_x_0_xxxx_xxyz[k] * ab_y + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxy_xxzz[k] = -g_x_0_xxxx_xxzz[k] * ab_y + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxy_xyyy[k] = -g_x_0_xxxx_xyyy[k] * ab_y + g_x_0_xxxx_xyyyy[k];

                g_x_0_xxxxy_xyyz[k] = -g_x_0_xxxx_xyyz[k] * ab_y + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxy_xyzz[k] = -g_x_0_xxxx_xyzz[k] * ab_y + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxy_xzzz[k] = -g_x_0_xxxx_xzzz[k] * ab_y + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxy_yyyy[k] = -g_x_0_xxxx_yyyy[k] * ab_y + g_x_0_xxxx_yyyyy[k];

                g_x_0_xxxxy_yyyz[k] = -g_x_0_xxxx_yyyz[k] * ab_y + g_x_0_xxxx_yyyyz[k];

                g_x_0_xxxxy_yyzz[k] = -g_x_0_xxxx_yyzz[k] * ab_y + g_x_0_xxxx_yyyzz[k];

                g_x_0_xxxxy_yzzz[k] = -g_x_0_xxxx_yzzz[k] * ab_y + g_x_0_xxxx_yyzzz[k];

                g_x_0_xxxxy_zzzz[k] = -g_x_0_xxxx_zzzz[k] * ab_y + g_x_0_xxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxx[k] = -g_x_0_xxxx_xxxx[k] * ab_z + g_x_0_xxxx_xxxxz[k];

                g_x_0_xxxxz_xxxy[k] = -g_x_0_xxxx_xxxy[k] * ab_z + g_x_0_xxxx_xxxyz[k];

                g_x_0_xxxxz_xxxz[k] = -g_x_0_xxxx_xxxz[k] * ab_z + g_x_0_xxxx_xxxzz[k];

                g_x_0_xxxxz_xxyy[k] = -g_x_0_xxxx_xxyy[k] * ab_z + g_x_0_xxxx_xxyyz[k];

                g_x_0_xxxxz_xxyz[k] = -g_x_0_xxxx_xxyz[k] * ab_z + g_x_0_xxxx_xxyzz[k];

                g_x_0_xxxxz_xxzz[k] = -g_x_0_xxxx_xxzz[k] * ab_z + g_x_0_xxxx_xxzzz[k];

                g_x_0_xxxxz_xyyy[k] = -g_x_0_xxxx_xyyy[k] * ab_z + g_x_0_xxxx_xyyyz[k];

                g_x_0_xxxxz_xyyz[k] = -g_x_0_xxxx_xyyz[k] * ab_z + g_x_0_xxxx_xyyzz[k];

                g_x_0_xxxxz_xyzz[k] = -g_x_0_xxxx_xyzz[k] * ab_z + g_x_0_xxxx_xyzzz[k];

                g_x_0_xxxxz_xzzz[k] = -g_x_0_xxxx_xzzz[k] * ab_z + g_x_0_xxxx_xzzzz[k];

                g_x_0_xxxxz_yyyy[k] = -g_x_0_xxxx_yyyy[k] * ab_z + g_x_0_xxxx_yyyyz[k];

                g_x_0_xxxxz_yyyz[k] = -g_x_0_xxxx_yyyz[k] * ab_z + g_x_0_xxxx_yyyzz[k];

                g_x_0_xxxxz_yyzz[k] = -g_x_0_xxxx_yyzz[k] * ab_z + g_x_0_xxxx_yyzzz[k];

                g_x_0_xxxxz_yzzz[k] = -g_x_0_xxxx_yzzz[k] * ab_z + g_x_0_xxxx_yzzzz[k];

                g_x_0_xxxxz_zzzz[k] = -g_x_0_xxxx_zzzz[k] * ab_z + g_x_0_xxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxy_xxxx, g_x_0_xxxy_xxxxy, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxxyy, g_x_0_xxxy_xxxyz, g_x_0_xxxy_xxxz, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyyy, g_x_0_xxxy_xxyyz, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxyzz, g_x_0_xxxy_xxzz, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyyy, g_x_0_xxxy_xyyyz, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyyzz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xyzzz, g_x_0_xxxy_xzzz, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyyy, g_x_0_xxxy_yyyyz, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyyzz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yyzzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_yzzzz, g_x_0_xxxy_zzzz, g_x_0_xxxyy_xxxx, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxxz, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxzz, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xzzz, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxx[k] = -g_x_0_xxxy_xxxx[k] * ab_y + g_x_0_xxxy_xxxxy[k];

                g_x_0_xxxyy_xxxy[k] = -g_x_0_xxxy_xxxy[k] * ab_y + g_x_0_xxxy_xxxyy[k];

                g_x_0_xxxyy_xxxz[k] = -g_x_0_xxxy_xxxz[k] * ab_y + g_x_0_xxxy_xxxyz[k];

                g_x_0_xxxyy_xxyy[k] = -g_x_0_xxxy_xxyy[k] * ab_y + g_x_0_xxxy_xxyyy[k];

                g_x_0_xxxyy_xxyz[k] = -g_x_0_xxxy_xxyz[k] * ab_y + g_x_0_xxxy_xxyyz[k];

                g_x_0_xxxyy_xxzz[k] = -g_x_0_xxxy_xxzz[k] * ab_y + g_x_0_xxxy_xxyzz[k];

                g_x_0_xxxyy_xyyy[k] = -g_x_0_xxxy_xyyy[k] * ab_y + g_x_0_xxxy_xyyyy[k];

                g_x_0_xxxyy_xyyz[k] = -g_x_0_xxxy_xyyz[k] * ab_y + g_x_0_xxxy_xyyyz[k];

                g_x_0_xxxyy_xyzz[k] = -g_x_0_xxxy_xyzz[k] * ab_y + g_x_0_xxxy_xyyzz[k];

                g_x_0_xxxyy_xzzz[k] = -g_x_0_xxxy_xzzz[k] * ab_y + g_x_0_xxxy_xyzzz[k];

                g_x_0_xxxyy_yyyy[k] = -g_x_0_xxxy_yyyy[k] * ab_y + g_x_0_xxxy_yyyyy[k];

                g_x_0_xxxyy_yyyz[k] = -g_x_0_xxxy_yyyz[k] * ab_y + g_x_0_xxxy_yyyyz[k];

                g_x_0_xxxyy_yyzz[k] = -g_x_0_xxxy_yyzz[k] * ab_y + g_x_0_xxxy_yyyzz[k];

                g_x_0_xxxyy_yzzz[k] = -g_x_0_xxxy_yzzz[k] * ab_y + g_x_0_xxxy_yyzzz[k];

                g_x_0_xxxyy_zzzz[k] = -g_x_0_xxxy_zzzz[k] * ab_y + g_x_0_xxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxyz_xxxx, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxxz, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxzz, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xzzz, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_zzzz, g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxx[k] = -g_x_0_xxxz_xxxx[k] * ab_y + g_x_0_xxxz_xxxxy[k];

                g_x_0_xxxyz_xxxy[k] = -g_x_0_xxxz_xxxy[k] * ab_y + g_x_0_xxxz_xxxyy[k];

                g_x_0_xxxyz_xxxz[k] = -g_x_0_xxxz_xxxz[k] * ab_y + g_x_0_xxxz_xxxyz[k];

                g_x_0_xxxyz_xxyy[k] = -g_x_0_xxxz_xxyy[k] * ab_y + g_x_0_xxxz_xxyyy[k];

                g_x_0_xxxyz_xxyz[k] = -g_x_0_xxxz_xxyz[k] * ab_y + g_x_0_xxxz_xxyyz[k];

                g_x_0_xxxyz_xxzz[k] = -g_x_0_xxxz_xxzz[k] * ab_y + g_x_0_xxxz_xxyzz[k];

                g_x_0_xxxyz_xyyy[k] = -g_x_0_xxxz_xyyy[k] * ab_y + g_x_0_xxxz_xyyyy[k];

                g_x_0_xxxyz_xyyz[k] = -g_x_0_xxxz_xyyz[k] * ab_y + g_x_0_xxxz_xyyyz[k];

                g_x_0_xxxyz_xyzz[k] = -g_x_0_xxxz_xyzz[k] * ab_y + g_x_0_xxxz_xyyzz[k];

                g_x_0_xxxyz_xzzz[k] = -g_x_0_xxxz_xzzz[k] * ab_y + g_x_0_xxxz_xyzzz[k];

                g_x_0_xxxyz_yyyy[k] = -g_x_0_xxxz_yyyy[k] * ab_y + g_x_0_xxxz_yyyyy[k];

                g_x_0_xxxyz_yyyz[k] = -g_x_0_xxxz_yyyz[k] * ab_y + g_x_0_xxxz_yyyyz[k];

                g_x_0_xxxyz_yyzz[k] = -g_x_0_xxxz_yyzz[k] * ab_y + g_x_0_xxxz_yyyzz[k];

                g_x_0_xxxyz_yzzz[k] = -g_x_0_xxxz_yzzz[k] * ab_y + g_x_0_xxxz_yyzzz[k];

                g_x_0_xxxyz_zzzz[k] = -g_x_0_xxxz_zzzz[k] * ab_y + g_x_0_xxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_zzzz, g_x_0_xxxz_zzzzz, g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxx[k] = -g_x_0_xxxz_xxxx[k] * ab_z + g_x_0_xxxz_xxxxz[k];

                g_x_0_xxxzz_xxxy[k] = -g_x_0_xxxz_xxxy[k] * ab_z + g_x_0_xxxz_xxxyz[k];

                g_x_0_xxxzz_xxxz[k] = -g_x_0_xxxz_xxxz[k] * ab_z + g_x_0_xxxz_xxxzz[k];

                g_x_0_xxxzz_xxyy[k] = -g_x_0_xxxz_xxyy[k] * ab_z + g_x_0_xxxz_xxyyz[k];

                g_x_0_xxxzz_xxyz[k] = -g_x_0_xxxz_xxyz[k] * ab_z + g_x_0_xxxz_xxyzz[k];

                g_x_0_xxxzz_xxzz[k] = -g_x_0_xxxz_xxzz[k] * ab_z + g_x_0_xxxz_xxzzz[k];

                g_x_0_xxxzz_xyyy[k] = -g_x_0_xxxz_xyyy[k] * ab_z + g_x_0_xxxz_xyyyz[k];

                g_x_0_xxxzz_xyyz[k] = -g_x_0_xxxz_xyyz[k] * ab_z + g_x_0_xxxz_xyyzz[k];

                g_x_0_xxxzz_xyzz[k] = -g_x_0_xxxz_xyzz[k] * ab_z + g_x_0_xxxz_xyzzz[k];

                g_x_0_xxxzz_xzzz[k] = -g_x_0_xxxz_xzzz[k] * ab_z + g_x_0_xxxz_xzzzz[k];

                g_x_0_xxxzz_yyyy[k] = -g_x_0_xxxz_yyyy[k] * ab_z + g_x_0_xxxz_yyyyz[k];

                g_x_0_xxxzz_yyyz[k] = -g_x_0_xxxz_yyyz[k] * ab_z + g_x_0_xxxz_yyyzz[k];

                g_x_0_xxxzz_yyzz[k] = -g_x_0_xxxz_yyzz[k] * ab_z + g_x_0_xxxz_yyzzz[k];

                g_x_0_xxxzz_yzzz[k] = -g_x_0_xxxz_yzzz[k] * ab_z + g_x_0_xxxz_yzzzz[k];

                g_x_0_xxxzz_zzzz[k] = -g_x_0_xxxz_zzzz[k] * ab_z + g_x_0_xxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyy_xxxx, g_x_0_xxyy_xxxxy, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxxyy, g_x_0_xxyy_xxxyz, g_x_0_xxyy_xxxz, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyyy, g_x_0_xxyy_xxyyz, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxyzz, g_x_0_xxyy_xxzz, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyyy, g_x_0_xxyy_xyyyz, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyyzz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xyzzz, g_x_0_xxyy_xzzz, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyyy, g_x_0_xxyy_yyyyz, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyyzz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yyzzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_yzzzz, g_x_0_xxyy_zzzz, g_x_0_xxyyy_xxxx, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxxz, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxzz, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xzzz, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxx[k] = -g_x_0_xxyy_xxxx[k] * ab_y + g_x_0_xxyy_xxxxy[k];

                g_x_0_xxyyy_xxxy[k] = -g_x_0_xxyy_xxxy[k] * ab_y + g_x_0_xxyy_xxxyy[k];

                g_x_0_xxyyy_xxxz[k] = -g_x_0_xxyy_xxxz[k] * ab_y + g_x_0_xxyy_xxxyz[k];

                g_x_0_xxyyy_xxyy[k] = -g_x_0_xxyy_xxyy[k] * ab_y + g_x_0_xxyy_xxyyy[k];

                g_x_0_xxyyy_xxyz[k] = -g_x_0_xxyy_xxyz[k] * ab_y + g_x_0_xxyy_xxyyz[k];

                g_x_0_xxyyy_xxzz[k] = -g_x_0_xxyy_xxzz[k] * ab_y + g_x_0_xxyy_xxyzz[k];

                g_x_0_xxyyy_xyyy[k] = -g_x_0_xxyy_xyyy[k] * ab_y + g_x_0_xxyy_xyyyy[k];

                g_x_0_xxyyy_xyyz[k] = -g_x_0_xxyy_xyyz[k] * ab_y + g_x_0_xxyy_xyyyz[k];

                g_x_0_xxyyy_xyzz[k] = -g_x_0_xxyy_xyzz[k] * ab_y + g_x_0_xxyy_xyyzz[k];

                g_x_0_xxyyy_xzzz[k] = -g_x_0_xxyy_xzzz[k] * ab_y + g_x_0_xxyy_xyzzz[k];

                g_x_0_xxyyy_yyyy[k] = -g_x_0_xxyy_yyyy[k] * ab_y + g_x_0_xxyy_yyyyy[k];

                g_x_0_xxyyy_yyyz[k] = -g_x_0_xxyy_yyyz[k] * ab_y + g_x_0_xxyy_yyyyz[k];

                g_x_0_xxyyy_yyzz[k] = -g_x_0_xxyy_yyzz[k] * ab_y + g_x_0_xxyy_yyyzz[k];

                g_x_0_xxyyy_yzzz[k] = -g_x_0_xxyy_yzzz[k] * ab_y + g_x_0_xxyy_yyzzz[k];

                g_x_0_xxyyy_zzzz[k] = -g_x_0_xxyy_zzzz[k] * ab_y + g_x_0_xxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyyz_xxxx, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxxz, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxzz, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xzzz, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_zzzz, g_x_0_xxyz_xxxx, g_x_0_xxyz_xxxxy, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxxyy, g_x_0_xxyz_xxxyz, g_x_0_xxyz_xxxz, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyyy, g_x_0_xxyz_xxyyz, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxyzz, g_x_0_xxyz_xxzz, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyyy, g_x_0_xxyz_xyyyz, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyyzz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xyzzz, g_x_0_xxyz_xzzz, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyyy, g_x_0_xxyz_yyyyz, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyyzz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yyzzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_yzzzz, g_x_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxx[k] = -g_x_0_xxyz_xxxx[k] * ab_y + g_x_0_xxyz_xxxxy[k];

                g_x_0_xxyyz_xxxy[k] = -g_x_0_xxyz_xxxy[k] * ab_y + g_x_0_xxyz_xxxyy[k];

                g_x_0_xxyyz_xxxz[k] = -g_x_0_xxyz_xxxz[k] * ab_y + g_x_0_xxyz_xxxyz[k];

                g_x_0_xxyyz_xxyy[k] = -g_x_0_xxyz_xxyy[k] * ab_y + g_x_0_xxyz_xxyyy[k];

                g_x_0_xxyyz_xxyz[k] = -g_x_0_xxyz_xxyz[k] * ab_y + g_x_0_xxyz_xxyyz[k];

                g_x_0_xxyyz_xxzz[k] = -g_x_0_xxyz_xxzz[k] * ab_y + g_x_0_xxyz_xxyzz[k];

                g_x_0_xxyyz_xyyy[k] = -g_x_0_xxyz_xyyy[k] * ab_y + g_x_0_xxyz_xyyyy[k];

                g_x_0_xxyyz_xyyz[k] = -g_x_0_xxyz_xyyz[k] * ab_y + g_x_0_xxyz_xyyyz[k];

                g_x_0_xxyyz_xyzz[k] = -g_x_0_xxyz_xyzz[k] * ab_y + g_x_0_xxyz_xyyzz[k];

                g_x_0_xxyyz_xzzz[k] = -g_x_0_xxyz_xzzz[k] * ab_y + g_x_0_xxyz_xyzzz[k];

                g_x_0_xxyyz_yyyy[k] = -g_x_0_xxyz_yyyy[k] * ab_y + g_x_0_xxyz_yyyyy[k];

                g_x_0_xxyyz_yyyz[k] = -g_x_0_xxyz_yyyz[k] * ab_y + g_x_0_xxyz_yyyyz[k];

                g_x_0_xxyyz_yyzz[k] = -g_x_0_xxyz_yyzz[k] * ab_y + g_x_0_xxyz_yyyzz[k];

                g_x_0_xxyyz_yzzz[k] = -g_x_0_xxyz_yzzz[k] * ab_y + g_x_0_xxyz_yyzzz[k];

                g_x_0_xxyyz_zzzz[k] = -g_x_0_xxyz_zzzz[k] * ab_y + g_x_0_xxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxyzz_xxxx, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxxz, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxzz, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xzzz, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_zzzz, g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxx[k] = -g_x_0_xxzz_xxxx[k] * ab_y + g_x_0_xxzz_xxxxy[k];

                g_x_0_xxyzz_xxxy[k] = -g_x_0_xxzz_xxxy[k] * ab_y + g_x_0_xxzz_xxxyy[k];

                g_x_0_xxyzz_xxxz[k] = -g_x_0_xxzz_xxxz[k] * ab_y + g_x_0_xxzz_xxxyz[k];

                g_x_0_xxyzz_xxyy[k] = -g_x_0_xxzz_xxyy[k] * ab_y + g_x_0_xxzz_xxyyy[k];

                g_x_0_xxyzz_xxyz[k] = -g_x_0_xxzz_xxyz[k] * ab_y + g_x_0_xxzz_xxyyz[k];

                g_x_0_xxyzz_xxzz[k] = -g_x_0_xxzz_xxzz[k] * ab_y + g_x_0_xxzz_xxyzz[k];

                g_x_0_xxyzz_xyyy[k] = -g_x_0_xxzz_xyyy[k] * ab_y + g_x_0_xxzz_xyyyy[k];

                g_x_0_xxyzz_xyyz[k] = -g_x_0_xxzz_xyyz[k] * ab_y + g_x_0_xxzz_xyyyz[k];

                g_x_0_xxyzz_xyzz[k] = -g_x_0_xxzz_xyzz[k] * ab_y + g_x_0_xxzz_xyyzz[k];

                g_x_0_xxyzz_xzzz[k] = -g_x_0_xxzz_xzzz[k] * ab_y + g_x_0_xxzz_xyzzz[k];

                g_x_0_xxyzz_yyyy[k] = -g_x_0_xxzz_yyyy[k] * ab_y + g_x_0_xxzz_yyyyy[k];

                g_x_0_xxyzz_yyyz[k] = -g_x_0_xxzz_yyyz[k] * ab_y + g_x_0_xxzz_yyyyz[k];

                g_x_0_xxyzz_yyzz[k] = -g_x_0_xxzz_yyzz[k] * ab_y + g_x_0_xxzz_yyyzz[k];

                g_x_0_xxyzz_yzzz[k] = -g_x_0_xxzz_yzzz[k] * ab_y + g_x_0_xxzz_yyzzz[k];

                g_x_0_xxyzz_zzzz[k] = -g_x_0_xxzz_zzzz[k] * ab_y + g_x_0_xxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_zzzz, g_x_0_xxzz_zzzzz, g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxx[k] = -g_x_0_xxzz_xxxx[k] * ab_z + g_x_0_xxzz_xxxxz[k];

                g_x_0_xxzzz_xxxy[k] = -g_x_0_xxzz_xxxy[k] * ab_z + g_x_0_xxzz_xxxyz[k];

                g_x_0_xxzzz_xxxz[k] = -g_x_0_xxzz_xxxz[k] * ab_z + g_x_0_xxzz_xxxzz[k];

                g_x_0_xxzzz_xxyy[k] = -g_x_0_xxzz_xxyy[k] * ab_z + g_x_0_xxzz_xxyyz[k];

                g_x_0_xxzzz_xxyz[k] = -g_x_0_xxzz_xxyz[k] * ab_z + g_x_0_xxzz_xxyzz[k];

                g_x_0_xxzzz_xxzz[k] = -g_x_0_xxzz_xxzz[k] * ab_z + g_x_0_xxzz_xxzzz[k];

                g_x_0_xxzzz_xyyy[k] = -g_x_0_xxzz_xyyy[k] * ab_z + g_x_0_xxzz_xyyyz[k];

                g_x_0_xxzzz_xyyz[k] = -g_x_0_xxzz_xyyz[k] * ab_z + g_x_0_xxzz_xyyzz[k];

                g_x_0_xxzzz_xyzz[k] = -g_x_0_xxzz_xyzz[k] * ab_z + g_x_0_xxzz_xyzzz[k];

                g_x_0_xxzzz_xzzz[k] = -g_x_0_xxzz_xzzz[k] * ab_z + g_x_0_xxzz_xzzzz[k];

                g_x_0_xxzzz_yyyy[k] = -g_x_0_xxzz_yyyy[k] * ab_z + g_x_0_xxzz_yyyyz[k];

                g_x_0_xxzzz_yyyz[k] = -g_x_0_xxzz_yyyz[k] * ab_z + g_x_0_xxzz_yyyzz[k];

                g_x_0_xxzzz_yyzz[k] = -g_x_0_xxzz_yyzz[k] * ab_z + g_x_0_xxzz_yyzzz[k];

                g_x_0_xxzzz_yzzz[k] = -g_x_0_xxzz_yzzz[k] * ab_z + g_x_0_xxzz_yzzzz[k];

                g_x_0_xxzzz_zzzz[k] = -g_x_0_xxzz_zzzz[k] * ab_z + g_x_0_xxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyy_xxxx, g_x_0_xyyy_xxxxy, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxxyy, g_x_0_xyyy_xxxyz, g_x_0_xyyy_xxxz, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyyy, g_x_0_xyyy_xxyyz, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxyzz, g_x_0_xyyy_xxzz, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyyy, g_x_0_xyyy_xyyyz, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyyzz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xyzzz, g_x_0_xyyy_xzzz, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyyy, g_x_0_xyyy_yyyyz, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyyzz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yyzzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_yzzzz, g_x_0_xyyy_zzzz, g_x_0_xyyyy_xxxx, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxxz, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxzz, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xzzz, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxx[k] = -g_x_0_xyyy_xxxx[k] * ab_y + g_x_0_xyyy_xxxxy[k];

                g_x_0_xyyyy_xxxy[k] = -g_x_0_xyyy_xxxy[k] * ab_y + g_x_0_xyyy_xxxyy[k];

                g_x_0_xyyyy_xxxz[k] = -g_x_0_xyyy_xxxz[k] * ab_y + g_x_0_xyyy_xxxyz[k];

                g_x_0_xyyyy_xxyy[k] = -g_x_0_xyyy_xxyy[k] * ab_y + g_x_0_xyyy_xxyyy[k];

                g_x_0_xyyyy_xxyz[k] = -g_x_0_xyyy_xxyz[k] * ab_y + g_x_0_xyyy_xxyyz[k];

                g_x_0_xyyyy_xxzz[k] = -g_x_0_xyyy_xxzz[k] * ab_y + g_x_0_xyyy_xxyzz[k];

                g_x_0_xyyyy_xyyy[k] = -g_x_0_xyyy_xyyy[k] * ab_y + g_x_0_xyyy_xyyyy[k];

                g_x_0_xyyyy_xyyz[k] = -g_x_0_xyyy_xyyz[k] * ab_y + g_x_0_xyyy_xyyyz[k];

                g_x_0_xyyyy_xyzz[k] = -g_x_0_xyyy_xyzz[k] * ab_y + g_x_0_xyyy_xyyzz[k];

                g_x_0_xyyyy_xzzz[k] = -g_x_0_xyyy_xzzz[k] * ab_y + g_x_0_xyyy_xyzzz[k];

                g_x_0_xyyyy_yyyy[k] = -g_x_0_xyyy_yyyy[k] * ab_y + g_x_0_xyyy_yyyyy[k];

                g_x_0_xyyyy_yyyz[k] = -g_x_0_xyyy_yyyz[k] * ab_y + g_x_0_xyyy_yyyyz[k];

                g_x_0_xyyyy_yyzz[k] = -g_x_0_xyyy_yyzz[k] * ab_y + g_x_0_xyyy_yyyzz[k];

                g_x_0_xyyyy_yzzz[k] = -g_x_0_xyyy_yzzz[k] * ab_y + g_x_0_xyyy_yyzzz[k];

                g_x_0_xyyyy_zzzz[k] = -g_x_0_xyyy_zzzz[k] * ab_y + g_x_0_xyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyyz_xxxx, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxxz, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxzz, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xzzz, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_zzzz, g_x_0_xyyz_xxxx, g_x_0_xyyz_xxxxy, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxxyy, g_x_0_xyyz_xxxyz, g_x_0_xyyz_xxxz, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyyy, g_x_0_xyyz_xxyyz, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxyzz, g_x_0_xyyz_xxzz, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyyy, g_x_0_xyyz_xyyyz, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyyzz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xyzzz, g_x_0_xyyz_xzzz, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyyy, g_x_0_xyyz_yyyyz, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyyzz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yyzzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_yzzzz, g_x_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxx[k] = -g_x_0_xyyz_xxxx[k] * ab_y + g_x_0_xyyz_xxxxy[k];

                g_x_0_xyyyz_xxxy[k] = -g_x_0_xyyz_xxxy[k] * ab_y + g_x_0_xyyz_xxxyy[k];

                g_x_0_xyyyz_xxxz[k] = -g_x_0_xyyz_xxxz[k] * ab_y + g_x_0_xyyz_xxxyz[k];

                g_x_0_xyyyz_xxyy[k] = -g_x_0_xyyz_xxyy[k] * ab_y + g_x_0_xyyz_xxyyy[k];

                g_x_0_xyyyz_xxyz[k] = -g_x_0_xyyz_xxyz[k] * ab_y + g_x_0_xyyz_xxyyz[k];

                g_x_0_xyyyz_xxzz[k] = -g_x_0_xyyz_xxzz[k] * ab_y + g_x_0_xyyz_xxyzz[k];

                g_x_0_xyyyz_xyyy[k] = -g_x_0_xyyz_xyyy[k] * ab_y + g_x_0_xyyz_xyyyy[k];

                g_x_0_xyyyz_xyyz[k] = -g_x_0_xyyz_xyyz[k] * ab_y + g_x_0_xyyz_xyyyz[k];

                g_x_0_xyyyz_xyzz[k] = -g_x_0_xyyz_xyzz[k] * ab_y + g_x_0_xyyz_xyyzz[k];

                g_x_0_xyyyz_xzzz[k] = -g_x_0_xyyz_xzzz[k] * ab_y + g_x_0_xyyz_xyzzz[k];

                g_x_0_xyyyz_yyyy[k] = -g_x_0_xyyz_yyyy[k] * ab_y + g_x_0_xyyz_yyyyy[k];

                g_x_0_xyyyz_yyyz[k] = -g_x_0_xyyz_yyyz[k] * ab_y + g_x_0_xyyz_yyyyz[k];

                g_x_0_xyyyz_yyzz[k] = -g_x_0_xyyz_yyzz[k] * ab_y + g_x_0_xyyz_yyyzz[k];

                g_x_0_xyyyz_yzzz[k] = -g_x_0_xyyz_yzzz[k] * ab_y + g_x_0_xyyz_yyzzz[k];

                g_x_0_xyyyz_zzzz[k] = -g_x_0_xyyz_zzzz[k] * ab_y + g_x_0_xyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyyzz_xxxx, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxxz, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxzz, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xzzz, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_zzzz, g_x_0_xyzz_xxxx, g_x_0_xyzz_xxxxy, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxxyy, g_x_0_xyzz_xxxyz, g_x_0_xyzz_xxxz, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyyy, g_x_0_xyzz_xxyyz, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxyzz, g_x_0_xyzz_xxzz, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyyy, g_x_0_xyzz_xyyyz, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyyzz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xyzzz, g_x_0_xyzz_xzzz, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyyy, g_x_0_xyzz_yyyyz, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyyzz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yyzzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_yzzzz, g_x_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxx[k] = -g_x_0_xyzz_xxxx[k] * ab_y + g_x_0_xyzz_xxxxy[k];

                g_x_0_xyyzz_xxxy[k] = -g_x_0_xyzz_xxxy[k] * ab_y + g_x_0_xyzz_xxxyy[k];

                g_x_0_xyyzz_xxxz[k] = -g_x_0_xyzz_xxxz[k] * ab_y + g_x_0_xyzz_xxxyz[k];

                g_x_0_xyyzz_xxyy[k] = -g_x_0_xyzz_xxyy[k] * ab_y + g_x_0_xyzz_xxyyy[k];

                g_x_0_xyyzz_xxyz[k] = -g_x_0_xyzz_xxyz[k] * ab_y + g_x_0_xyzz_xxyyz[k];

                g_x_0_xyyzz_xxzz[k] = -g_x_0_xyzz_xxzz[k] * ab_y + g_x_0_xyzz_xxyzz[k];

                g_x_0_xyyzz_xyyy[k] = -g_x_0_xyzz_xyyy[k] * ab_y + g_x_0_xyzz_xyyyy[k];

                g_x_0_xyyzz_xyyz[k] = -g_x_0_xyzz_xyyz[k] * ab_y + g_x_0_xyzz_xyyyz[k];

                g_x_0_xyyzz_xyzz[k] = -g_x_0_xyzz_xyzz[k] * ab_y + g_x_0_xyzz_xyyzz[k];

                g_x_0_xyyzz_xzzz[k] = -g_x_0_xyzz_xzzz[k] * ab_y + g_x_0_xyzz_xyzzz[k];

                g_x_0_xyyzz_yyyy[k] = -g_x_0_xyzz_yyyy[k] * ab_y + g_x_0_xyzz_yyyyy[k];

                g_x_0_xyyzz_yyyz[k] = -g_x_0_xyzz_yyyz[k] * ab_y + g_x_0_xyzz_yyyyz[k];

                g_x_0_xyyzz_yyzz[k] = -g_x_0_xyzz_yyzz[k] * ab_y + g_x_0_xyzz_yyyzz[k];

                g_x_0_xyyzz_yzzz[k] = -g_x_0_xyzz_yzzz[k] * ab_y + g_x_0_xyzz_yyzzz[k];

                g_x_0_xyyzz_zzzz[k] = -g_x_0_xyzz_zzzz[k] * ab_y + g_x_0_xyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xyzzz_xxxx, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxxz, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxzz, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xzzz, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_zzzz, g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxx[k] = -g_x_0_xzzz_xxxx[k] * ab_y + g_x_0_xzzz_xxxxy[k];

                g_x_0_xyzzz_xxxy[k] = -g_x_0_xzzz_xxxy[k] * ab_y + g_x_0_xzzz_xxxyy[k];

                g_x_0_xyzzz_xxxz[k] = -g_x_0_xzzz_xxxz[k] * ab_y + g_x_0_xzzz_xxxyz[k];

                g_x_0_xyzzz_xxyy[k] = -g_x_0_xzzz_xxyy[k] * ab_y + g_x_0_xzzz_xxyyy[k];

                g_x_0_xyzzz_xxyz[k] = -g_x_0_xzzz_xxyz[k] * ab_y + g_x_0_xzzz_xxyyz[k];

                g_x_0_xyzzz_xxzz[k] = -g_x_0_xzzz_xxzz[k] * ab_y + g_x_0_xzzz_xxyzz[k];

                g_x_0_xyzzz_xyyy[k] = -g_x_0_xzzz_xyyy[k] * ab_y + g_x_0_xzzz_xyyyy[k];

                g_x_0_xyzzz_xyyz[k] = -g_x_0_xzzz_xyyz[k] * ab_y + g_x_0_xzzz_xyyyz[k];

                g_x_0_xyzzz_xyzz[k] = -g_x_0_xzzz_xyzz[k] * ab_y + g_x_0_xzzz_xyyzz[k];

                g_x_0_xyzzz_xzzz[k] = -g_x_0_xzzz_xzzz[k] * ab_y + g_x_0_xzzz_xyzzz[k];

                g_x_0_xyzzz_yyyy[k] = -g_x_0_xzzz_yyyy[k] * ab_y + g_x_0_xzzz_yyyyy[k];

                g_x_0_xyzzz_yyyz[k] = -g_x_0_xzzz_yyyz[k] * ab_y + g_x_0_xzzz_yyyyz[k];

                g_x_0_xyzzz_yyzz[k] = -g_x_0_xzzz_yyzz[k] * ab_y + g_x_0_xzzz_yyyzz[k];

                g_x_0_xyzzz_yzzz[k] = -g_x_0_xzzz_yzzz[k] * ab_y + g_x_0_xzzz_yyzzz[k];

                g_x_0_xyzzz_zzzz[k] = -g_x_0_xzzz_zzzz[k] * ab_y + g_x_0_xzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_zzzz, g_x_0_xzzz_zzzzz, g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxx[k] = -g_x_0_xzzz_xxxx[k] * ab_z + g_x_0_xzzz_xxxxz[k];

                g_x_0_xzzzz_xxxy[k] = -g_x_0_xzzz_xxxy[k] * ab_z + g_x_0_xzzz_xxxyz[k];

                g_x_0_xzzzz_xxxz[k] = -g_x_0_xzzz_xxxz[k] * ab_z + g_x_0_xzzz_xxxzz[k];

                g_x_0_xzzzz_xxyy[k] = -g_x_0_xzzz_xxyy[k] * ab_z + g_x_0_xzzz_xxyyz[k];

                g_x_0_xzzzz_xxyz[k] = -g_x_0_xzzz_xxyz[k] * ab_z + g_x_0_xzzz_xxyzz[k];

                g_x_0_xzzzz_xxzz[k] = -g_x_0_xzzz_xxzz[k] * ab_z + g_x_0_xzzz_xxzzz[k];

                g_x_0_xzzzz_xyyy[k] = -g_x_0_xzzz_xyyy[k] * ab_z + g_x_0_xzzz_xyyyz[k];

                g_x_0_xzzzz_xyyz[k] = -g_x_0_xzzz_xyyz[k] * ab_z + g_x_0_xzzz_xyyzz[k];

                g_x_0_xzzzz_xyzz[k] = -g_x_0_xzzz_xyzz[k] * ab_z + g_x_0_xzzz_xyzzz[k];

                g_x_0_xzzzz_xzzz[k] = -g_x_0_xzzz_xzzz[k] * ab_z + g_x_0_xzzz_xzzzz[k];

                g_x_0_xzzzz_yyyy[k] = -g_x_0_xzzz_yyyy[k] * ab_z + g_x_0_xzzz_yyyyz[k];

                g_x_0_xzzzz_yyyz[k] = -g_x_0_xzzz_yyyz[k] * ab_z + g_x_0_xzzz_yyyzz[k];

                g_x_0_xzzzz_yyzz[k] = -g_x_0_xzzz_yyzz[k] * ab_z + g_x_0_xzzz_yyzzz[k];

                g_x_0_xzzzz_yzzz[k] = -g_x_0_xzzz_yzzz[k] * ab_z + g_x_0_xzzz_yzzzz[k];

                g_x_0_xzzzz_zzzz[k] = -g_x_0_xzzz_zzzz[k] * ab_z + g_x_0_xzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyy_xxxx, g_x_0_yyyy_xxxxy, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxxyy, g_x_0_yyyy_xxxyz, g_x_0_yyyy_xxxz, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyyy, g_x_0_yyyy_xxyyz, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxyzz, g_x_0_yyyy_xxzz, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyyy, g_x_0_yyyy_xyyyz, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyyzz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xyzzz, g_x_0_yyyy_xzzz, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyyy, g_x_0_yyyy_yyyyz, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyyzz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yyzzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_yzzzz, g_x_0_yyyy_zzzz, g_x_0_yyyyy_xxxx, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxxz, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxzz, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xzzz, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxx[k] = -g_x_0_yyyy_xxxx[k] * ab_y + g_x_0_yyyy_xxxxy[k];

                g_x_0_yyyyy_xxxy[k] = -g_x_0_yyyy_xxxy[k] * ab_y + g_x_0_yyyy_xxxyy[k];

                g_x_0_yyyyy_xxxz[k] = -g_x_0_yyyy_xxxz[k] * ab_y + g_x_0_yyyy_xxxyz[k];

                g_x_0_yyyyy_xxyy[k] = -g_x_0_yyyy_xxyy[k] * ab_y + g_x_0_yyyy_xxyyy[k];

                g_x_0_yyyyy_xxyz[k] = -g_x_0_yyyy_xxyz[k] * ab_y + g_x_0_yyyy_xxyyz[k];

                g_x_0_yyyyy_xxzz[k] = -g_x_0_yyyy_xxzz[k] * ab_y + g_x_0_yyyy_xxyzz[k];

                g_x_0_yyyyy_xyyy[k] = -g_x_0_yyyy_xyyy[k] * ab_y + g_x_0_yyyy_xyyyy[k];

                g_x_0_yyyyy_xyyz[k] = -g_x_0_yyyy_xyyz[k] * ab_y + g_x_0_yyyy_xyyyz[k];

                g_x_0_yyyyy_xyzz[k] = -g_x_0_yyyy_xyzz[k] * ab_y + g_x_0_yyyy_xyyzz[k];

                g_x_0_yyyyy_xzzz[k] = -g_x_0_yyyy_xzzz[k] * ab_y + g_x_0_yyyy_xyzzz[k];

                g_x_0_yyyyy_yyyy[k] = -g_x_0_yyyy_yyyy[k] * ab_y + g_x_0_yyyy_yyyyy[k];

                g_x_0_yyyyy_yyyz[k] = -g_x_0_yyyy_yyyz[k] * ab_y + g_x_0_yyyy_yyyyz[k];

                g_x_0_yyyyy_yyzz[k] = -g_x_0_yyyy_yyzz[k] * ab_y + g_x_0_yyyy_yyyzz[k];

                g_x_0_yyyyy_yzzz[k] = -g_x_0_yyyy_yzzz[k] * ab_y + g_x_0_yyyy_yyzzz[k];

                g_x_0_yyyyy_zzzz[k] = -g_x_0_yyyy_zzzz[k] * ab_y + g_x_0_yyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyyz_xxxx, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxxz, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxzz, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xzzz, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_zzzz, g_x_0_yyyz_xxxx, g_x_0_yyyz_xxxxy, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxxyy, g_x_0_yyyz_xxxyz, g_x_0_yyyz_xxxz, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyyy, g_x_0_yyyz_xxyyz, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxyzz, g_x_0_yyyz_xxzz, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyyy, g_x_0_yyyz_xyyyz, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyyzz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xyzzz, g_x_0_yyyz_xzzz, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyyy, g_x_0_yyyz_yyyyz, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyyzz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yyzzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_yzzzz, g_x_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxx[k] = -g_x_0_yyyz_xxxx[k] * ab_y + g_x_0_yyyz_xxxxy[k];

                g_x_0_yyyyz_xxxy[k] = -g_x_0_yyyz_xxxy[k] * ab_y + g_x_0_yyyz_xxxyy[k];

                g_x_0_yyyyz_xxxz[k] = -g_x_0_yyyz_xxxz[k] * ab_y + g_x_0_yyyz_xxxyz[k];

                g_x_0_yyyyz_xxyy[k] = -g_x_0_yyyz_xxyy[k] * ab_y + g_x_0_yyyz_xxyyy[k];

                g_x_0_yyyyz_xxyz[k] = -g_x_0_yyyz_xxyz[k] * ab_y + g_x_0_yyyz_xxyyz[k];

                g_x_0_yyyyz_xxzz[k] = -g_x_0_yyyz_xxzz[k] * ab_y + g_x_0_yyyz_xxyzz[k];

                g_x_0_yyyyz_xyyy[k] = -g_x_0_yyyz_xyyy[k] * ab_y + g_x_0_yyyz_xyyyy[k];

                g_x_0_yyyyz_xyyz[k] = -g_x_0_yyyz_xyyz[k] * ab_y + g_x_0_yyyz_xyyyz[k];

                g_x_0_yyyyz_xyzz[k] = -g_x_0_yyyz_xyzz[k] * ab_y + g_x_0_yyyz_xyyzz[k];

                g_x_0_yyyyz_xzzz[k] = -g_x_0_yyyz_xzzz[k] * ab_y + g_x_0_yyyz_xyzzz[k];

                g_x_0_yyyyz_yyyy[k] = -g_x_0_yyyz_yyyy[k] * ab_y + g_x_0_yyyz_yyyyy[k];

                g_x_0_yyyyz_yyyz[k] = -g_x_0_yyyz_yyyz[k] * ab_y + g_x_0_yyyz_yyyyz[k];

                g_x_0_yyyyz_yyzz[k] = -g_x_0_yyyz_yyzz[k] * ab_y + g_x_0_yyyz_yyyzz[k];

                g_x_0_yyyyz_yzzz[k] = -g_x_0_yyyz_yzzz[k] * ab_y + g_x_0_yyyz_yyzzz[k];

                g_x_0_yyyyz_zzzz[k] = -g_x_0_yyyz_zzzz[k] * ab_y + g_x_0_yyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyyzz_xxxx, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxxz, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxzz, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xzzz, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_zzzz, g_x_0_yyzz_xxxx, g_x_0_yyzz_xxxxy, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxxyy, g_x_0_yyzz_xxxyz, g_x_0_yyzz_xxxz, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyyy, g_x_0_yyzz_xxyyz, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxyzz, g_x_0_yyzz_xxzz, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyyy, g_x_0_yyzz_xyyyz, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyyzz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xyzzz, g_x_0_yyzz_xzzz, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyyy, g_x_0_yyzz_yyyyz, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyyzz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yyzzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_yzzzz, g_x_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxx[k] = -g_x_0_yyzz_xxxx[k] * ab_y + g_x_0_yyzz_xxxxy[k];

                g_x_0_yyyzz_xxxy[k] = -g_x_0_yyzz_xxxy[k] * ab_y + g_x_0_yyzz_xxxyy[k];

                g_x_0_yyyzz_xxxz[k] = -g_x_0_yyzz_xxxz[k] * ab_y + g_x_0_yyzz_xxxyz[k];

                g_x_0_yyyzz_xxyy[k] = -g_x_0_yyzz_xxyy[k] * ab_y + g_x_0_yyzz_xxyyy[k];

                g_x_0_yyyzz_xxyz[k] = -g_x_0_yyzz_xxyz[k] * ab_y + g_x_0_yyzz_xxyyz[k];

                g_x_0_yyyzz_xxzz[k] = -g_x_0_yyzz_xxzz[k] * ab_y + g_x_0_yyzz_xxyzz[k];

                g_x_0_yyyzz_xyyy[k] = -g_x_0_yyzz_xyyy[k] * ab_y + g_x_0_yyzz_xyyyy[k];

                g_x_0_yyyzz_xyyz[k] = -g_x_0_yyzz_xyyz[k] * ab_y + g_x_0_yyzz_xyyyz[k];

                g_x_0_yyyzz_xyzz[k] = -g_x_0_yyzz_xyzz[k] * ab_y + g_x_0_yyzz_xyyzz[k];

                g_x_0_yyyzz_xzzz[k] = -g_x_0_yyzz_xzzz[k] * ab_y + g_x_0_yyzz_xyzzz[k];

                g_x_0_yyyzz_yyyy[k] = -g_x_0_yyzz_yyyy[k] * ab_y + g_x_0_yyzz_yyyyy[k];

                g_x_0_yyyzz_yyyz[k] = -g_x_0_yyzz_yyyz[k] * ab_y + g_x_0_yyzz_yyyyz[k];

                g_x_0_yyyzz_yyzz[k] = -g_x_0_yyzz_yyzz[k] * ab_y + g_x_0_yyzz_yyyzz[k];

                g_x_0_yyyzz_yzzz[k] = -g_x_0_yyzz_yzzz[k] * ab_y + g_x_0_yyzz_yyzzz[k];

                g_x_0_yyyzz_zzzz[k] = -g_x_0_yyzz_zzzz[k] * ab_y + g_x_0_yyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yyzzz_xxxx, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxxz, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxzz, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xzzz, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_zzzz, g_x_0_yzzz_xxxx, g_x_0_yzzz_xxxxy, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxxyy, g_x_0_yzzz_xxxyz, g_x_0_yzzz_xxxz, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyyy, g_x_0_yzzz_xxyyz, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxyzz, g_x_0_yzzz_xxzz, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyyy, g_x_0_yzzz_xyyyz, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyyzz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xyzzz, g_x_0_yzzz_xzzz, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyyy, g_x_0_yzzz_yyyyz, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyyzz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yyzzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_yzzzz, g_x_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxx[k] = -g_x_0_yzzz_xxxx[k] * ab_y + g_x_0_yzzz_xxxxy[k];

                g_x_0_yyzzz_xxxy[k] = -g_x_0_yzzz_xxxy[k] * ab_y + g_x_0_yzzz_xxxyy[k];

                g_x_0_yyzzz_xxxz[k] = -g_x_0_yzzz_xxxz[k] * ab_y + g_x_0_yzzz_xxxyz[k];

                g_x_0_yyzzz_xxyy[k] = -g_x_0_yzzz_xxyy[k] * ab_y + g_x_0_yzzz_xxyyy[k];

                g_x_0_yyzzz_xxyz[k] = -g_x_0_yzzz_xxyz[k] * ab_y + g_x_0_yzzz_xxyyz[k];

                g_x_0_yyzzz_xxzz[k] = -g_x_0_yzzz_xxzz[k] * ab_y + g_x_0_yzzz_xxyzz[k];

                g_x_0_yyzzz_xyyy[k] = -g_x_0_yzzz_xyyy[k] * ab_y + g_x_0_yzzz_xyyyy[k];

                g_x_0_yyzzz_xyyz[k] = -g_x_0_yzzz_xyyz[k] * ab_y + g_x_0_yzzz_xyyyz[k];

                g_x_0_yyzzz_xyzz[k] = -g_x_0_yzzz_xyzz[k] * ab_y + g_x_0_yzzz_xyyzz[k];

                g_x_0_yyzzz_xzzz[k] = -g_x_0_yzzz_xzzz[k] * ab_y + g_x_0_yzzz_xyzzz[k];

                g_x_0_yyzzz_yyyy[k] = -g_x_0_yzzz_yyyy[k] * ab_y + g_x_0_yzzz_yyyyy[k];

                g_x_0_yyzzz_yyyz[k] = -g_x_0_yzzz_yyyz[k] * ab_y + g_x_0_yzzz_yyyyz[k];

                g_x_0_yyzzz_yyzz[k] = -g_x_0_yzzz_yyzz[k] * ab_y + g_x_0_yzzz_yyyzz[k];

                g_x_0_yyzzz_yzzz[k] = -g_x_0_yzzz_yzzz[k] * ab_y + g_x_0_yzzz_yyzzz[k];

                g_x_0_yyzzz_zzzz[k] = -g_x_0_yzzz_zzzz[k] * ab_y + g_x_0_yzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_yzzzz_xxxx, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxxz, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxzz, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xzzz, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_zzzz, g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxx[k] = -g_x_0_zzzz_xxxx[k] * ab_y + g_x_0_zzzz_xxxxy[k];

                g_x_0_yzzzz_xxxy[k] = -g_x_0_zzzz_xxxy[k] * ab_y + g_x_0_zzzz_xxxyy[k];

                g_x_0_yzzzz_xxxz[k] = -g_x_0_zzzz_xxxz[k] * ab_y + g_x_0_zzzz_xxxyz[k];

                g_x_0_yzzzz_xxyy[k] = -g_x_0_zzzz_xxyy[k] * ab_y + g_x_0_zzzz_xxyyy[k];

                g_x_0_yzzzz_xxyz[k] = -g_x_0_zzzz_xxyz[k] * ab_y + g_x_0_zzzz_xxyyz[k];

                g_x_0_yzzzz_xxzz[k] = -g_x_0_zzzz_xxzz[k] * ab_y + g_x_0_zzzz_xxyzz[k];

                g_x_0_yzzzz_xyyy[k] = -g_x_0_zzzz_xyyy[k] * ab_y + g_x_0_zzzz_xyyyy[k];

                g_x_0_yzzzz_xyyz[k] = -g_x_0_zzzz_xyyz[k] * ab_y + g_x_0_zzzz_xyyyz[k];

                g_x_0_yzzzz_xyzz[k] = -g_x_0_zzzz_xyzz[k] * ab_y + g_x_0_zzzz_xyyzz[k];

                g_x_0_yzzzz_xzzz[k] = -g_x_0_zzzz_xzzz[k] * ab_y + g_x_0_zzzz_xyzzz[k];

                g_x_0_yzzzz_yyyy[k] = -g_x_0_zzzz_yyyy[k] * ab_y + g_x_0_zzzz_yyyyy[k];

                g_x_0_yzzzz_yyyz[k] = -g_x_0_zzzz_yyyz[k] * ab_y + g_x_0_zzzz_yyyyz[k];

                g_x_0_yzzzz_yyzz[k] = -g_x_0_zzzz_yyzz[k] * ab_y + g_x_0_zzzz_yyyzz[k];

                g_x_0_yzzzz_yzzz[k] = -g_x_0_zzzz_yzzz[k] * ab_y + g_x_0_zzzz_yyzzz[k];

                g_x_0_yzzzz_zzzz[k] = -g_x_0_zzzz_zzzz[k] * ab_y + g_x_0_zzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_zzzz, g_x_0_zzzz_zzzzz, g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxx[k] = -g_x_0_zzzz_xxxx[k] * ab_z + g_x_0_zzzz_xxxxz[k];

                g_x_0_zzzzz_xxxy[k] = -g_x_0_zzzz_xxxy[k] * ab_z + g_x_0_zzzz_xxxyz[k];

                g_x_0_zzzzz_xxxz[k] = -g_x_0_zzzz_xxxz[k] * ab_z + g_x_0_zzzz_xxxzz[k];

                g_x_0_zzzzz_xxyy[k] = -g_x_0_zzzz_xxyy[k] * ab_z + g_x_0_zzzz_xxyyz[k];

                g_x_0_zzzzz_xxyz[k] = -g_x_0_zzzz_xxyz[k] * ab_z + g_x_0_zzzz_xxyzz[k];

                g_x_0_zzzzz_xxzz[k] = -g_x_0_zzzz_xxzz[k] * ab_z + g_x_0_zzzz_xxzzz[k];

                g_x_0_zzzzz_xyyy[k] = -g_x_0_zzzz_xyyy[k] * ab_z + g_x_0_zzzz_xyyyz[k];

                g_x_0_zzzzz_xyyz[k] = -g_x_0_zzzz_xyyz[k] * ab_z + g_x_0_zzzz_xyyzz[k];

                g_x_0_zzzzz_xyzz[k] = -g_x_0_zzzz_xyzz[k] * ab_z + g_x_0_zzzz_xyzzz[k];

                g_x_0_zzzzz_xzzz[k] = -g_x_0_zzzz_xzzz[k] * ab_z + g_x_0_zzzz_xzzzz[k];

                g_x_0_zzzzz_yyyy[k] = -g_x_0_zzzz_yyyy[k] * ab_z + g_x_0_zzzz_yyyyz[k];

                g_x_0_zzzzz_yyyz[k] = -g_x_0_zzzz_yyyz[k] * ab_z + g_x_0_zzzz_yyyzz[k];

                g_x_0_zzzzz_yyzz[k] = -g_x_0_zzzz_yyzz[k] * ab_z + g_x_0_zzzz_yyzzz[k];

                g_x_0_zzzzz_yzzz[k] = -g_x_0_zzzz_yzzz[k] * ab_z + g_x_0_zzzz_yzzzz[k];

                g_x_0_zzzzz_zzzz[k] = -g_x_0_zzzz_zzzz[k] * ab_z + g_x_0_zzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxxx, g_y_0_xxxx_xxxxy, g_y_0_xxxx_xxxxz, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxyy, g_y_0_xxxx_xxxyz, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxxzz, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyyy, g_y_0_xxxx_xxyyz, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxyzz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xxzzz, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyyy, g_y_0_xxxx_xyyyz, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyyzz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xyzzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_xzzzz, g_y_0_xxxx_yyyy, g_y_0_xxxx_yyyz, g_y_0_xxxx_yyzz, g_y_0_xxxx_yzzz, g_y_0_xxxx_zzzz, g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_yyyy, g_y_0_xxxxx_yyyz, g_y_0_xxxxx_yyzz, g_y_0_xxxxx_yzzz, g_y_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxx[k] = -g_y_0_xxxx_xxxx[k] * ab_x + g_y_0_xxxx_xxxxx[k];

                g_y_0_xxxxx_xxxy[k] = -g_y_0_xxxx_xxxy[k] * ab_x + g_y_0_xxxx_xxxxy[k];

                g_y_0_xxxxx_xxxz[k] = -g_y_0_xxxx_xxxz[k] * ab_x + g_y_0_xxxx_xxxxz[k];

                g_y_0_xxxxx_xxyy[k] = -g_y_0_xxxx_xxyy[k] * ab_x + g_y_0_xxxx_xxxyy[k];

                g_y_0_xxxxx_xxyz[k] = -g_y_0_xxxx_xxyz[k] * ab_x + g_y_0_xxxx_xxxyz[k];

                g_y_0_xxxxx_xxzz[k] = -g_y_0_xxxx_xxzz[k] * ab_x + g_y_0_xxxx_xxxzz[k];

                g_y_0_xxxxx_xyyy[k] = -g_y_0_xxxx_xyyy[k] * ab_x + g_y_0_xxxx_xxyyy[k];

                g_y_0_xxxxx_xyyz[k] = -g_y_0_xxxx_xyyz[k] * ab_x + g_y_0_xxxx_xxyyz[k];

                g_y_0_xxxxx_xyzz[k] = -g_y_0_xxxx_xyzz[k] * ab_x + g_y_0_xxxx_xxyzz[k];

                g_y_0_xxxxx_xzzz[k] = -g_y_0_xxxx_xzzz[k] * ab_x + g_y_0_xxxx_xxzzz[k];

                g_y_0_xxxxx_yyyy[k] = -g_y_0_xxxx_yyyy[k] * ab_x + g_y_0_xxxx_xyyyy[k];

                g_y_0_xxxxx_yyyz[k] = -g_y_0_xxxx_yyyz[k] * ab_x + g_y_0_xxxx_xyyyz[k];

                g_y_0_xxxxx_yyzz[k] = -g_y_0_xxxx_yyzz[k] * ab_x + g_y_0_xxxx_xyyzz[k];

                g_y_0_xxxxx_yzzz[k] = -g_y_0_xxxx_yzzz[k] * ab_x + g_y_0_xxxx_xyzzz[k];

                g_y_0_xxxxx_zzzz[k] = -g_y_0_xxxx_zzzz[k] * ab_x + g_y_0_xxxx_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_yyyy, g_y_0_xxxxy_yyyz, g_y_0_xxxxy_yyzz, g_y_0_xxxxy_yzzz, g_y_0_xxxxy_zzzz, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxxx, g_y_0_xxxy_xxxxy, g_y_0_xxxy_xxxxz, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxyy, g_y_0_xxxy_xxxyz, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxxzz, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyyy, g_y_0_xxxy_xxyyz, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxyzz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xxzzz, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyyy, g_y_0_xxxy_xyyyz, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyyzz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xyzzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_xzzzz, g_y_0_xxxy_yyyy, g_y_0_xxxy_yyyz, g_y_0_xxxy_yyzz, g_y_0_xxxy_yzzz, g_y_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxx[k] = -g_y_0_xxxy_xxxx[k] * ab_x + g_y_0_xxxy_xxxxx[k];

                g_y_0_xxxxy_xxxy[k] = -g_y_0_xxxy_xxxy[k] * ab_x + g_y_0_xxxy_xxxxy[k];

                g_y_0_xxxxy_xxxz[k] = -g_y_0_xxxy_xxxz[k] * ab_x + g_y_0_xxxy_xxxxz[k];

                g_y_0_xxxxy_xxyy[k] = -g_y_0_xxxy_xxyy[k] * ab_x + g_y_0_xxxy_xxxyy[k];

                g_y_0_xxxxy_xxyz[k] = -g_y_0_xxxy_xxyz[k] * ab_x + g_y_0_xxxy_xxxyz[k];

                g_y_0_xxxxy_xxzz[k] = -g_y_0_xxxy_xxzz[k] * ab_x + g_y_0_xxxy_xxxzz[k];

                g_y_0_xxxxy_xyyy[k] = -g_y_0_xxxy_xyyy[k] * ab_x + g_y_0_xxxy_xxyyy[k];

                g_y_0_xxxxy_xyyz[k] = -g_y_0_xxxy_xyyz[k] * ab_x + g_y_0_xxxy_xxyyz[k];

                g_y_0_xxxxy_xyzz[k] = -g_y_0_xxxy_xyzz[k] * ab_x + g_y_0_xxxy_xxyzz[k];

                g_y_0_xxxxy_xzzz[k] = -g_y_0_xxxy_xzzz[k] * ab_x + g_y_0_xxxy_xxzzz[k];

                g_y_0_xxxxy_yyyy[k] = -g_y_0_xxxy_yyyy[k] * ab_x + g_y_0_xxxy_xyyyy[k];

                g_y_0_xxxxy_yyyz[k] = -g_y_0_xxxy_yyyz[k] * ab_x + g_y_0_xxxy_xyyyz[k];

                g_y_0_xxxxy_yyzz[k] = -g_y_0_xxxy_yyzz[k] * ab_x + g_y_0_xxxy_xyyzz[k];

                g_y_0_xxxxy_yzzz[k] = -g_y_0_xxxy_yzzz[k] * ab_x + g_y_0_xxxy_xyzzz[k];

                g_y_0_xxxxy_zzzz[k] = -g_y_0_xxxy_zzzz[k] * ab_x + g_y_0_xxxy_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_yyyy, g_y_0_xxxxz_yyyz, g_y_0_xxxxz_yyzz, g_y_0_xxxxz_yzzz, g_y_0_xxxxz_zzzz, g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxxx, g_y_0_xxxz_xxxxy, g_y_0_xxxz_xxxxz, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxyy, g_y_0_xxxz_xxxyz, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxxzz, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyyy, g_y_0_xxxz_xxyyz, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxyzz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xxzzz, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyyy, g_y_0_xxxz_xyyyz, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyyzz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xyzzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_xzzzz, g_y_0_xxxz_yyyy, g_y_0_xxxz_yyyz, g_y_0_xxxz_yyzz, g_y_0_xxxz_yzzz, g_y_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxx[k] = -g_y_0_xxxz_xxxx[k] * ab_x + g_y_0_xxxz_xxxxx[k];

                g_y_0_xxxxz_xxxy[k] = -g_y_0_xxxz_xxxy[k] * ab_x + g_y_0_xxxz_xxxxy[k];

                g_y_0_xxxxz_xxxz[k] = -g_y_0_xxxz_xxxz[k] * ab_x + g_y_0_xxxz_xxxxz[k];

                g_y_0_xxxxz_xxyy[k] = -g_y_0_xxxz_xxyy[k] * ab_x + g_y_0_xxxz_xxxyy[k];

                g_y_0_xxxxz_xxyz[k] = -g_y_0_xxxz_xxyz[k] * ab_x + g_y_0_xxxz_xxxyz[k];

                g_y_0_xxxxz_xxzz[k] = -g_y_0_xxxz_xxzz[k] * ab_x + g_y_0_xxxz_xxxzz[k];

                g_y_0_xxxxz_xyyy[k] = -g_y_0_xxxz_xyyy[k] * ab_x + g_y_0_xxxz_xxyyy[k];

                g_y_0_xxxxz_xyyz[k] = -g_y_0_xxxz_xyyz[k] * ab_x + g_y_0_xxxz_xxyyz[k];

                g_y_0_xxxxz_xyzz[k] = -g_y_0_xxxz_xyzz[k] * ab_x + g_y_0_xxxz_xxyzz[k];

                g_y_0_xxxxz_xzzz[k] = -g_y_0_xxxz_xzzz[k] * ab_x + g_y_0_xxxz_xxzzz[k];

                g_y_0_xxxxz_yyyy[k] = -g_y_0_xxxz_yyyy[k] * ab_x + g_y_0_xxxz_xyyyy[k];

                g_y_0_xxxxz_yyyz[k] = -g_y_0_xxxz_yyyz[k] * ab_x + g_y_0_xxxz_xyyyz[k];

                g_y_0_xxxxz_yyzz[k] = -g_y_0_xxxz_yyzz[k] * ab_x + g_y_0_xxxz_xyyzz[k];

                g_y_0_xxxxz_yzzz[k] = -g_y_0_xxxz_yzzz[k] * ab_x + g_y_0_xxxz_xyzzz[k];

                g_y_0_xxxxz_zzzz[k] = -g_y_0_xxxz_zzzz[k] * ab_x + g_y_0_xxxz_xzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_yyyy, g_y_0_xxxyy_yyyz, g_y_0_xxxyy_yyzz, g_y_0_xxxyy_yzzz, g_y_0_xxxyy_zzzz, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxxx, g_y_0_xxyy_xxxxy, g_y_0_xxyy_xxxxz, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxyy, g_y_0_xxyy_xxxyz, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxxzz, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyyy, g_y_0_xxyy_xxyyz, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxyzz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xxzzz, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyyy, g_y_0_xxyy_xyyyz, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyyzz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xyzzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_xzzzz, g_y_0_xxyy_yyyy, g_y_0_xxyy_yyyz, g_y_0_xxyy_yyzz, g_y_0_xxyy_yzzz, g_y_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxx[k] = -g_y_0_xxyy_xxxx[k] * ab_x + g_y_0_xxyy_xxxxx[k];

                g_y_0_xxxyy_xxxy[k] = -g_y_0_xxyy_xxxy[k] * ab_x + g_y_0_xxyy_xxxxy[k];

                g_y_0_xxxyy_xxxz[k] = -g_y_0_xxyy_xxxz[k] * ab_x + g_y_0_xxyy_xxxxz[k];

                g_y_0_xxxyy_xxyy[k] = -g_y_0_xxyy_xxyy[k] * ab_x + g_y_0_xxyy_xxxyy[k];

                g_y_0_xxxyy_xxyz[k] = -g_y_0_xxyy_xxyz[k] * ab_x + g_y_0_xxyy_xxxyz[k];

                g_y_0_xxxyy_xxzz[k] = -g_y_0_xxyy_xxzz[k] * ab_x + g_y_0_xxyy_xxxzz[k];

                g_y_0_xxxyy_xyyy[k] = -g_y_0_xxyy_xyyy[k] * ab_x + g_y_0_xxyy_xxyyy[k];

                g_y_0_xxxyy_xyyz[k] = -g_y_0_xxyy_xyyz[k] * ab_x + g_y_0_xxyy_xxyyz[k];

                g_y_0_xxxyy_xyzz[k] = -g_y_0_xxyy_xyzz[k] * ab_x + g_y_0_xxyy_xxyzz[k];

                g_y_0_xxxyy_xzzz[k] = -g_y_0_xxyy_xzzz[k] * ab_x + g_y_0_xxyy_xxzzz[k];

                g_y_0_xxxyy_yyyy[k] = -g_y_0_xxyy_yyyy[k] * ab_x + g_y_0_xxyy_xyyyy[k];

                g_y_0_xxxyy_yyyz[k] = -g_y_0_xxyy_yyyz[k] * ab_x + g_y_0_xxyy_xyyyz[k];

                g_y_0_xxxyy_yyzz[k] = -g_y_0_xxyy_yyzz[k] * ab_x + g_y_0_xxyy_xyyzz[k];

                g_y_0_xxxyy_yzzz[k] = -g_y_0_xxyy_yzzz[k] * ab_x + g_y_0_xxyy_xyzzz[k];

                g_y_0_xxxyy_zzzz[k] = -g_y_0_xxyy_zzzz[k] * ab_x + g_y_0_xxyy_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_yyyy, g_y_0_xxxyz_yyyz, g_y_0_xxxyz_yyzz, g_y_0_xxxyz_yzzz, g_y_0_xxxyz_zzzz, g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxxx, g_y_0_xxyz_xxxxy, g_y_0_xxyz_xxxxz, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxyy, g_y_0_xxyz_xxxyz, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxxzz, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyyy, g_y_0_xxyz_xxyyz, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxyzz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xxzzz, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyyy, g_y_0_xxyz_xyyyz, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyyzz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xyzzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_xzzzz, g_y_0_xxyz_yyyy, g_y_0_xxyz_yyyz, g_y_0_xxyz_yyzz, g_y_0_xxyz_yzzz, g_y_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxx[k] = -g_y_0_xxyz_xxxx[k] * ab_x + g_y_0_xxyz_xxxxx[k];

                g_y_0_xxxyz_xxxy[k] = -g_y_0_xxyz_xxxy[k] * ab_x + g_y_0_xxyz_xxxxy[k];

                g_y_0_xxxyz_xxxz[k] = -g_y_0_xxyz_xxxz[k] * ab_x + g_y_0_xxyz_xxxxz[k];

                g_y_0_xxxyz_xxyy[k] = -g_y_0_xxyz_xxyy[k] * ab_x + g_y_0_xxyz_xxxyy[k];

                g_y_0_xxxyz_xxyz[k] = -g_y_0_xxyz_xxyz[k] * ab_x + g_y_0_xxyz_xxxyz[k];

                g_y_0_xxxyz_xxzz[k] = -g_y_0_xxyz_xxzz[k] * ab_x + g_y_0_xxyz_xxxzz[k];

                g_y_0_xxxyz_xyyy[k] = -g_y_0_xxyz_xyyy[k] * ab_x + g_y_0_xxyz_xxyyy[k];

                g_y_0_xxxyz_xyyz[k] = -g_y_0_xxyz_xyyz[k] * ab_x + g_y_0_xxyz_xxyyz[k];

                g_y_0_xxxyz_xyzz[k] = -g_y_0_xxyz_xyzz[k] * ab_x + g_y_0_xxyz_xxyzz[k];

                g_y_0_xxxyz_xzzz[k] = -g_y_0_xxyz_xzzz[k] * ab_x + g_y_0_xxyz_xxzzz[k];

                g_y_0_xxxyz_yyyy[k] = -g_y_0_xxyz_yyyy[k] * ab_x + g_y_0_xxyz_xyyyy[k];

                g_y_0_xxxyz_yyyz[k] = -g_y_0_xxyz_yyyz[k] * ab_x + g_y_0_xxyz_xyyyz[k];

                g_y_0_xxxyz_yyzz[k] = -g_y_0_xxyz_yyzz[k] * ab_x + g_y_0_xxyz_xyyzz[k];

                g_y_0_xxxyz_yzzz[k] = -g_y_0_xxyz_yzzz[k] * ab_x + g_y_0_xxyz_xyzzz[k];

                g_y_0_xxxyz_zzzz[k] = -g_y_0_xxyz_zzzz[k] * ab_x + g_y_0_xxyz_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_yyyy, g_y_0_xxxzz_yyyz, g_y_0_xxxzz_yyzz, g_y_0_xxxzz_yzzz, g_y_0_xxxzz_zzzz, g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxxx, g_y_0_xxzz_xxxxy, g_y_0_xxzz_xxxxz, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxyy, g_y_0_xxzz_xxxyz, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxxzz, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyyy, g_y_0_xxzz_xxyyz, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxyzz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xxzzz, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyyy, g_y_0_xxzz_xyyyz, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyyzz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xyzzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_xzzzz, g_y_0_xxzz_yyyy, g_y_0_xxzz_yyyz, g_y_0_xxzz_yyzz, g_y_0_xxzz_yzzz, g_y_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxx[k] = -g_y_0_xxzz_xxxx[k] * ab_x + g_y_0_xxzz_xxxxx[k];

                g_y_0_xxxzz_xxxy[k] = -g_y_0_xxzz_xxxy[k] * ab_x + g_y_0_xxzz_xxxxy[k];

                g_y_0_xxxzz_xxxz[k] = -g_y_0_xxzz_xxxz[k] * ab_x + g_y_0_xxzz_xxxxz[k];

                g_y_0_xxxzz_xxyy[k] = -g_y_0_xxzz_xxyy[k] * ab_x + g_y_0_xxzz_xxxyy[k];

                g_y_0_xxxzz_xxyz[k] = -g_y_0_xxzz_xxyz[k] * ab_x + g_y_0_xxzz_xxxyz[k];

                g_y_0_xxxzz_xxzz[k] = -g_y_0_xxzz_xxzz[k] * ab_x + g_y_0_xxzz_xxxzz[k];

                g_y_0_xxxzz_xyyy[k] = -g_y_0_xxzz_xyyy[k] * ab_x + g_y_0_xxzz_xxyyy[k];

                g_y_0_xxxzz_xyyz[k] = -g_y_0_xxzz_xyyz[k] * ab_x + g_y_0_xxzz_xxyyz[k];

                g_y_0_xxxzz_xyzz[k] = -g_y_0_xxzz_xyzz[k] * ab_x + g_y_0_xxzz_xxyzz[k];

                g_y_0_xxxzz_xzzz[k] = -g_y_0_xxzz_xzzz[k] * ab_x + g_y_0_xxzz_xxzzz[k];

                g_y_0_xxxzz_yyyy[k] = -g_y_0_xxzz_yyyy[k] * ab_x + g_y_0_xxzz_xyyyy[k];

                g_y_0_xxxzz_yyyz[k] = -g_y_0_xxzz_yyyz[k] * ab_x + g_y_0_xxzz_xyyyz[k];

                g_y_0_xxxzz_yyzz[k] = -g_y_0_xxzz_yyzz[k] * ab_x + g_y_0_xxzz_xyyzz[k];

                g_y_0_xxxzz_yzzz[k] = -g_y_0_xxzz_yzzz[k] * ab_x + g_y_0_xxzz_xyzzz[k];

                g_y_0_xxxzz_zzzz[k] = -g_y_0_xxzz_zzzz[k] * ab_x + g_y_0_xxzz_xzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_yyyy, g_y_0_xxyyy_yyyz, g_y_0_xxyyy_yyzz, g_y_0_xxyyy_yzzz, g_y_0_xxyyy_zzzz, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxxx, g_y_0_xyyy_xxxxy, g_y_0_xyyy_xxxxz, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxyy, g_y_0_xyyy_xxxyz, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxxzz, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyyy, g_y_0_xyyy_xxyyz, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxyzz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xxzzz, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyyy, g_y_0_xyyy_xyyyz, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyyzz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xyzzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_xzzzz, g_y_0_xyyy_yyyy, g_y_0_xyyy_yyyz, g_y_0_xyyy_yyzz, g_y_0_xyyy_yzzz, g_y_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxx[k] = -g_y_0_xyyy_xxxx[k] * ab_x + g_y_0_xyyy_xxxxx[k];

                g_y_0_xxyyy_xxxy[k] = -g_y_0_xyyy_xxxy[k] * ab_x + g_y_0_xyyy_xxxxy[k];

                g_y_0_xxyyy_xxxz[k] = -g_y_0_xyyy_xxxz[k] * ab_x + g_y_0_xyyy_xxxxz[k];

                g_y_0_xxyyy_xxyy[k] = -g_y_0_xyyy_xxyy[k] * ab_x + g_y_0_xyyy_xxxyy[k];

                g_y_0_xxyyy_xxyz[k] = -g_y_0_xyyy_xxyz[k] * ab_x + g_y_0_xyyy_xxxyz[k];

                g_y_0_xxyyy_xxzz[k] = -g_y_0_xyyy_xxzz[k] * ab_x + g_y_0_xyyy_xxxzz[k];

                g_y_0_xxyyy_xyyy[k] = -g_y_0_xyyy_xyyy[k] * ab_x + g_y_0_xyyy_xxyyy[k];

                g_y_0_xxyyy_xyyz[k] = -g_y_0_xyyy_xyyz[k] * ab_x + g_y_0_xyyy_xxyyz[k];

                g_y_0_xxyyy_xyzz[k] = -g_y_0_xyyy_xyzz[k] * ab_x + g_y_0_xyyy_xxyzz[k];

                g_y_0_xxyyy_xzzz[k] = -g_y_0_xyyy_xzzz[k] * ab_x + g_y_0_xyyy_xxzzz[k];

                g_y_0_xxyyy_yyyy[k] = -g_y_0_xyyy_yyyy[k] * ab_x + g_y_0_xyyy_xyyyy[k];

                g_y_0_xxyyy_yyyz[k] = -g_y_0_xyyy_yyyz[k] * ab_x + g_y_0_xyyy_xyyyz[k];

                g_y_0_xxyyy_yyzz[k] = -g_y_0_xyyy_yyzz[k] * ab_x + g_y_0_xyyy_xyyzz[k];

                g_y_0_xxyyy_yzzz[k] = -g_y_0_xyyy_yzzz[k] * ab_x + g_y_0_xyyy_xyzzz[k];

                g_y_0_xxyyy_zzzz[k] = -g_y_0_xyyy_zzzz[k] * ab_x + g_y_0_xyyy_xzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_yyyy, g_y_0_xxyyz_yyyz, g_y_0_xxyyz_yyzz, g_y_0_xxyyz_yzzz, g_y_0_xxyyz_zzzz, g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxxx, g_y_0_xyyz_xxxxy, g_y_0_xyyz_xxxxz, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxyy, g_y_0_xyyz_xxxyz, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxxzz, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyyy, g_y_0_xyyz_xxyyz, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxyzz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xxzzz, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyyy, g_y_0_xyyz_xyyyz, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyyzz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xyzzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_xzzzz, g_y_0_xyyz_yyyy, g_y_0_xyyz_yyyz, g_y_0_xyyz_yyzz, g_y_0_xyyz_yzzz, g_y_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxx[k] = -g_y_0_xyyz_xxxx[k] * ab_x + g_y_0_xyyz_xxxxx[k];

                g_y_0_xxyyz_xxxy[k] = -g_y_0_xyyz_xxxy[k] * ab_x + g_y_0_xyyz_xxxxy[k];

                g_y_0_xxyyz_xxxz[k] = -g_y_0_xyyz_xxxz[k] * ab_x + g_y_0_xyyz_xxxxz[k];

                g_y_0_xxyyz_xxyy[k] = -g_y_0_xyyz_xxyy[k] * ab_x + g_y_0_xyyz_xxxyy[k];

                g_y_0_xxyyz_xxyz[k] = -g_y_0_xyyz_xxyz[k] * ab_x + g_y_0_xyyz_xxxyz[k];

                g_y_0_xxyyz_xxzz[k] = -g_y_0_xyyz_xxzz[k] * ab_x + g_y_0_xyyz_xxxzz[k];

                g_y_0_xxyyz_xyyy[k] = -g_y_0_xyyz_xyyy[k] * ab_x + g_y_0_xyyz_xxyyy[k];

                g_y_0_xxyyz_xyyz[k] = -g_y_0_xyyz_xyyz[k] * ab_x + g_y_0_xyyz_xxyyz[k];

                g_y_0_xxyyz_xyzz[k] = -g_y_0_xyyz_xyzz[k] * ab_x + g_y_0_xyyz_xxyzz[k];

                g_y_0_xxyyz_xzzz[k] = -g_y_0_xyyz_xzzz[k] * ab_x + g_y_0_xyyz_xxzzz[k];

                g_y_0_xxyyz_yyyy[k] = -g_y_0_xyyz_yyyy[k] * ab_x + g_y_0_xyyz_xyyyy[k];

                g_y_0_xxyyz_yyyz[k] = -g_y_0_xyyz_yyyz[k] * ab_x + g_y_0_xyyz_xyyyz[k];

                g_y_0_xxyyz_yyzz[k] = -g_y_0_xyyz_yyzz[k] * ab_x + g_y_0_xyyz_xyyzz[k];

                g_y_0_xxyyz_yzzz[k] = -g_y_0_xyyz_yzzz[k] * ab_x + g_y_0_xyyz_xyzzz[k];

                g_y_0_xxyyz_zzzz[k] = -g_y_0_xyyz_zzzz[k] * ab_x + g_y_0_xyyz_xzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_yyyy, g_y_0_xxyzz_yyyz, g_y_0_xxyzz_yyzz, g_y_0_xxyzz_yzzz, g_y_0_xxyzz_zzzz, g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxxx, g_y_0_xyzz_xxxxy, g_y_0_xyzz_xxxxz, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxyy, g_y_0_xyzz_xxxyz, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxxzz, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyyy, g_y_0_xyzz_xxyyz, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxyzz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xxzzz, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyyy, g_y_0_xyzz_xyyyz, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyyzz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xyzzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_xzzzz, g_y_0_xyzz_yyyy, g_y_0_xyzz_yyyz, g_y_0_xyzz_yyzz, g_y_0_xyzz_yzzz, g_y_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxx[k] = -g_y_0_xyzz_xxxx[k] * ab_x + g_y_0_xyzz_xxxxx[k];

                g_y_0_xxyzz_xxxy[k] = -g_y_0_xyzz_xxxy[k] * ab_x + g_y_0_xyzz_xxxxy[k];

                g_y_0_xxyzz_xxxz[k] = -g_y_0_xyzz_xxxz[k] * ab_x + g_y_0_xyzz_xxxxz[k];

                g_y_0_xxyzz_xxyy[k] = -g_y_0_xyzz_xxyy[k] * ab_x + g_y_0_xyzz_xxxyy[k];

                g_y_0_xxyzz_xxyz[k] = -g_y_0_xyzz_xxyz[k] * ab_x + g_y_0_xyzz_xxxyz[k];

                g_y_0_xxyzz_xxzz[k] = -g_y_0_xyzz_xxzz[k] * ab_x + g_y_0_xyzz_xxxzz[k];

                g_y_0_xxyzz_xyyy[k] = -g_y_0_xyzz_xyyy[k] * ab_x + g_y_0_xyzz_xxyyy[k];

                g_y_0_xxyzz_xyyz[k] = -g_y_0_xyzz_xyyz[k] * ab_x + g_y_0_xyzz_xxyyz[k];

                g_y_0_xxyzz_xyzz[k] = -g_y_0_xyzz_xyzz[k] * ab_x + g_y_0_xyzz_xxyzz[k];

                g_y_0_xxyzz_xzzz[k] = -g_y_0_xyzz_xzzz[k] * ab_x + g_y_0_xyzz_xxzzz[k];

                g_y_0_xxyzz_yyyy[k] = -g_y_0_xyzz_yyyy[k] * ab_x + g_y_0_xyzz_xyyyy[k];

                g_y_0_xxyzz_yyyz[k] = -g_y_0_xyzz_yyyz[k] * ab_x + g_y_0_xyzz_xyyyz[k];

                g_y_0_xxyzz_yyzz[k] = -g_y_0_xyzz_yyzz[k] * ab_x + g_y_0_xyzz_xyyzz[k];

                g_y_0_xxyzz_yzzz[k] = -g_y_0_xyzz_yzzz[k] * ab_x + g_y_0_xyzz_xyzzz[k];

                g_y_0_xxyzz_zzzz[k] = -g_y_0_xyzz_zzzz[k] * ab_x + g_y_0_xyzz_xzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_yyyy, g_y_0_xxzzz_yyyz, g_y_0_xxzzz_yyzz, g_y_0_xxzzz_yzzz, g_y_0_xxzzz_zzzz, g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxxx, g_y_0_xzzz_xxxxy, g_y_0_xzzz_xxxxz, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxyy, g_y_0_xzzz_xxxyz, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxxzz, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyyy, g_y_0_xzzz_xxyyz, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxyzz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xxzzz, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyyy, g_y_0_xzzz_xyyyz, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyyzz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xyzzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_xzzzz, g_y_0_xzzz_yyyy, g_y_0_xzzz_yyyz, g_y_0_xzzz_yyzz, g_y_0_xzzz_yzzz, g_y_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxx[k] = -g_y_0_xzzz_xxxx[k] * ab_x + g_y_0_xzzz_xxxxx[k];

                g_y_0_xxzzz_xxxy[k] = -g_y_0_xzzz_xxxy[k] * ab_x + g_y_0_xzzz_xxxxy[k];

                g_y_0_xxzzz_xxxz[k] = -g_y_0_xzzz_xxxz[k] * ab_x + g_y_0_xzzz_xxxxz[k];

                g_y_0_xxzzz_xxyy[k] = -g_y_0_xzzz_xxyy[k] * ab_x + g_y_0_xzzz_xxxyy[k];

                g_y_0_xxzzz_xxyz[k] = -g_y_0_xzzz_xxyz[k] * ab_x + g_y_0_xzzz_xxxyz[k];

                g_y_0_xxzzz_xxzz[k] = -g_y_0_xzzz_xxzz[k] * ab_x + g_y_0_xzzz_xxxzz[k];

                g_y_0_xxzzz_xyyy[k] = -g_y_0_xzzz_xyyy[k] * ab_x + g_y_0_xzzz_xxyyy[k];

                g_y_0_xxzzz_xyyz[k] = -g_y_0_xzzz_xyyz[k] * ab_x + g_y_0_xzzz_xxyyz[k];

                g_y_0_xxzzz_xyzz[k] = -g_y_0_xzzz_xyzz[k] * ab_x + g_y_0_xzzz_xxyzz[k];

                g_y_0_xxzzz_xzzz[k] = -g_y_0_xzzz_xzzz[k] * ab_x + g_y_0_xzzz_xxzzz[k];

                g_y_0_xxzzz_yyyy[k] = -g_y_0_xzzz_yyyy[k] * ab_x + g_y_0_xzzz_xyyyy[k];

                g_y_0_xxzzz_yyyz[k] = -g_y_0_xzzz_yyyz[k] * ab_x + g_y_0_xzzz_xyyyz[k];

                g_y_0_xxzzz_yyzz[k] = -g_y_0_xzzz_yyzz[k] * ab_x + g_y_0_xzzz_xyyzz[k];

                g_y_0_xxzzz_yzzz[k] = -g_y_0_xzzz_yzzz[k] * ab_x + g_y_0_xzzz_xyzzz[k];

                g_y_0_xxzzz_zzzz[k] = -g_y_0_xzzz_zzzz[k] * ab_x + g_y_0_xzzz_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_yyyy, g_y_0_xyyyy_yyyz, g_y_0_xyyyy_yyzz, g_y_0_xyyyy_yzzz, g_y_0_xyyyy_zzzz, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxx[k] = -g_y_0_yyyy_xxxx[k] * ab_x + g_y_0_yyyy_xxxxx[k];

                g_y_0_xyyyy_xxxy[k] = -g_y_0_yyyy_xxxy[k] * ab_x + g_y_0_yyyy_xxxxy[k];

                g_y_0_xyyyy_xxxz[k] = -g_y_0_yyyy_xxxz[k] * ab_x + g_y_0_yyyy_xxxxz[k];

                g_y_0_xyyyy_xxyy[k] = -g_y_0_yyyy_xxyy[k] * ab_x + g_y_0_yyyy_xxxyy[k];

                g_y_0_xyyyy_xxyz[k] = -g_y_0_yyyy_xxyz[k] * ab_x + g_y_0_yyyy_xxxyz[k];

                g_y_0_xyyyy_xxzz[k] = -g_y_0_yyyy_xxzz[k] * ab_x + g_y_0_yyyy_xxxzz[k];

                g_y_0_xyyyy_xyyy[k] = -g_y_0_yyyy_xyyy[k] * ab_x + g_y_0_yyyy_xxyyy[k];

                g_y_0_xyyyy_xyyz[k] = -g_y_0_yyyy_xyyz[k] * ab_x + g_y_0_yyyy_xxyyz[k];

                g_y_0_xyyyy_xyzz[k] = -g_y_0_yyyy_xyzz[k] * ab_x + g_y_0_yyyy_xxyzz[k];

                g_y_0_xyyyy_xzzz[k] = -g_y_0_yyyy_xzzz[k] * ab_x + g_y_0_yyyy_xxzzz[k];

                g_y_0_xyyyy_yyyy[k] = -g_y_0_yyyy_yyyy[k] * ab_x + g_y_0_yyyy_xyyyy[k];

                g_y_0_xyyyy_yyyz[k] = -g_y_0_yyyy_yyyz[k] * ab_x + g_y_0_yyyy_xyyyz[k];

                g_y_0_xyyyy_yyzz[k] = -g_y_0_yyyy_yyzz[k] * ab_x + g_y_0_yyyy_xyyzz[k];

                g_y_0_xyyyy_yzzz[k] = -g_y_0_yyyy_yzzz[k] * ab_x + g_y_0_yyyy_xyzzz[k];

                g_y_0_xyyyy_zzzz[k] = -g_y_0_yyyy_zzzz[k] * ab_x + g_y_0_yyyy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_yyyy, g_y_0_xyyyz_yyyz, g_y_0_xyyyz_yyzz, g_y_0_xyyyz_yzzz, g_y_0_xyyyz_zzzz, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxx[k] = -g_y_0_yyyz_xxxx[k] * ab_x + g_y_0_yyyz_xxxxx[k];

                g_y_0_xyyyz_xxxy[k] = -g_y_0_yyyz_xxxy[k] * ab_x + g_y_0_yyyz_xxxxy[k];

                g_y_0_xyyyz_xxxz[k] = -g_y_0_yyyz_xxxz[k] * ab_x + g_y_0_yyyz_xxxxz[k];

                g_y_0_xyyyz_xxyy[k] = -g_y_0_yyyz_xxyy[k] * ab_x + g_y_0_yyyz_xxxyy[k];

                g_y_0_xyyyz_xxyz[k] = -g_y_0_yyyz_xxyz[k] * ab_x + g_y_0_yyyz_xxxyz[k];

                g_y_0_xyyyz_xxzz[k] = -g_y_0_yyyz_xxzz[k] * ab_x + g_y_0_yyyz_xxxzz[k];

                g_y_0_xyyyz_xyyy[k] = -g_y_0_yyyz_xyyy[k] * ab_x + g_y_0_yyyz_xxyyy[k];

                g_y_0_xyyyz_xyyz[k] = -g_y_0_yyyz_xyyz[k] * ab_x + g_y_0_yyyz_xxyyz[k];

                g_y_0_xyyyz_xyzz[k] = -g_y_0_yyyz_xyzz[k] * ab_x + g_y_0_yyyz_xxyzz[k];

                g_y_0_xyyyz_xzzz[k] = -g_y_0_yyyz_xzzz[k] * ab_x + g_y_0_yyyz_xxzzz[k];

                g_y_0_xyyyz_yyyy[k] = -g_y_0_yyyz_yyyy[k] * ab_x + g_y_0_yyyz_xyyyy[k];

                g_y_0_xyyyz_yyyz[k] = -g_y_0_yyyz_yyyz[k] * ab_x + g_y_0_yyyz_xyyyz[k];

                g_y_0_xyyyz_yyzz[k] = -g_y_0_yyyz_yyzz[k] * ab_x + g_y_0_yyyz_xyyzz[k];

                g_y_0_xyyyz_yzzz[k] = -g_y_0_yyyz_yzzz[k] * ab_x + g_y_0_yyyz_xyzzz[k];

                g_y_0_xyyyz_zzzz[k] = -g_y_0_yyyz_zzzz[k] * ab_x + g_y_0_yyyz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_yyyy, g_y_0_xyyzz_yyyz, g_y_0_xyyzz_yyzz, g_y_0_xyyzz_yzzz, g_y_0_xyyzz_zzzz, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxx[k] = -g_y_0_yyzz_xxxx[k] * ab_x + g_y_0_yyzz_xxxxx[k];

                g_y_0_xyyzz_xxxy[k] = -g_y_0_yyzz_xxxy[k] * ab_x + g_y_0_yyzz_xxxxy[k];

                g_y_0_xyyzz_xxxz[k] = -g_y_0_yyzz_xxxz[k] * ab_x + g_y_0_yyzz_xxxxz[k];

                g_y_0_xyyzz_xxyy[k] = -g_y_0_yyzz_xxyy[k] * ab_x + g_y_0_yyzz_xxxyy[k];

                g_y_0_xyyzz_xxyz[k] = -g_y_0_yyzz_xxyz[k] * ab_x + g_y_0_yyzz_xxxyz[k];

                g_y_0_xyyzz_xxzz[k] = -g_y_0_yyzz_xxzz[k] * ab_x + g_y_0_yyzz_xxxzz[k];

                g_y_0_xyyzz_xyyy[k] = -g_y_0_yyzz_xyyy[k] * ab_x + g_y_0_yyzz_xxyyy[k];

                g_y_0_xyyzz_xyyz[k] = -g_y_0_yyzz_xyyz[k] * ab_x + g_y_0_yyzz_xxyyz[k];

                g_y_0_xyyzz_xyzz[k] = -g_y_0_yyzz_xyzz[k] * ab_x + g_y_0_yyzz_xxyzz[k];

                g_y_0_xyyzz_xzzz[k] = -g_y_0_yyzz_xzzz[k] * ab_x + g_y_0_yyzz_xxzzz[k];

                g_y_0_xyyzz_yyyy[k] = -g_y_0_yyzz_yyyy[k] * ab_x + g_y_0_yyzz_xyyyy[k];

                g_y_0_xyyzz_yyyz[k] = -g_y_0_yyzz_yyyz[k] * ab_x + g_y_0_yyzz_xyyyz[k];

                g_y_0_xyyzz_yyzz[k] = -g_y_0_yyzz_yyzz[k] * ab_x + g_y_0_yyzz_xyyzz[k];

                g_y_0_xyyzz_yzzz[k] = -g_y_0_yyzz_yzzz[k] * ab_x + g_y_0_yyzz_xyzzz[k];

                g_y_0_xyyzz_zzzz[k] = -g_y_0_yyzz_zzzz[k] * ab_x + g_y_0_yyzz_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_yyyy, g_y_0_xyzzz_yyyz, g_y_0_xyzzz_yyzz, g_y_0_xyzzz_yzzz, g_y_0_xyzzz_zzzz, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxx[k] = -g_y_0_yzzz_xxxx[k] * ab_x + g_y_0_yzzz_xxxxx[k];

                g_y_0_xyzzz_xxxy[k] = -g_y_0_yzzz_xxxy[k] * ab_x + g_y_0_yzzz_xxxxy[k];

                g_y_0_xyzzz_xxxz[k] = -g_y_0_yzzz_xxxz[k] * ab_x + g_y_0_yzzz_xxxxz[k];

                g_y_0_xyzzz_xxyy[k] = -g_y_0_yzzz_xxyy[k] * ab_x + g_y_0_yzzz_xxxyy[k];

                g_y_0_xyzzz_xxyz[k] = -g_y_0_yzzz_xxyz[k] * ab_x + g_y_0_yzzz_xxxyz[k];

                g_y_0_xyzzz_xxzz[k] = -g_y_0_yzzz_xxzz[k] * ab_x + g_y_0_yzzz_xxxzz[k];

                g_y_0_xyzzz_xyyy[k] = -g_y_0_yzzz_xyyy[k] * ab_x + g_y_0_yzzz_xxyyy[k];

                g_y_0_xyzzz_xyyz[k] = -g_y_0_yzzz_xyyz[k] * ab_x + g_y_0_yzzz_xxyyz[k];

                g_y_0_xyzzz_xyzz[k] = -g_y_0_yzzz_xyzz[k] * ab_x + g_y_0_yzzz_xxyzz[k];

                g_y_0_xyzzz_xzzz[k] = -g_y_0_yzzz_xzzz[k] * ab_x + g_y_0_yzzz_xxzzz[k];

                g_y_0_xyzzz_yyyy[k] = -g_y_0_yzzz_yyyy[k] * ab_x + g_y_0_yzzz_xyyyy[k];

                g_y_0_xyzzz_yyyz[k] = -g_y_0_yzzz_yyyz[k] * ab_x + g_y_0_yzzz_xyyyz[k];

                g_y_0_xyzzz_yyzz[k] = -g_y_0_yzzz_yyzz[k] * ab_x + g_y_0_yzzz_xyyzz[k];

                g_y_0_xyzzz_yzzz[k] = -g_y_0_yzzz_yzzz[k] * ab_x + g_y_0_yzzz_xyzzz[k];

                g_y_0_xyzzz_zzzz[k] = -g_y_0_yzzz_zzzz[k] * ab_x + g_y_0_yzzz_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_yyyy, g_y_0_xzzzz_yyyz, g_y_0_xzzzz_yyzz, g_y_0_xzzzz_yzzz, g_y_0_xzzzz_zzzz, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxx[k] = -g_y_0_zzzz_xxxx[k] * ab_x + g_y_0_zzzz_xxxxx[k];

                g_y_0_xzzzz_xxxy[k] = -g_y_0_zzzz_xxxy[k] * ab_x + g_y_0_zzzz_xxxxy[k];

                g_y_0_xzzzz_xxxz[k] = -g_y_0_zzzz_xxxz[k] * ab_x + g_y_0_zzzz_xxxxz[k];

                g_y_0_xzzzz_xxyy[k] = -g_y_0_zzzz_xxyy[k] * ab_x + g_y_0_zzzz_xxxyy[k];

                g_y_0_xzzzz_xxyz[k] = -g_y_0_zzzz_xxyz[k] * ab_x + g_y_0_zzzz_xxxyz[k];

                g_y_0_xzzzz_xxzz[k] = -g_y_0_zzzz_xxzz[k] * ab_x + g_y_0_zzzz_xxxzz[k];

                g_y_0_xzzzz_xyyy[k] = -g_y_0_zzzz_xyyy[k] * ab_x + g_y_0_zzzz_xxyyy[k];

                g_y_0_xzzzz_xyyz[k] = -g_y_0_zzzz_xyyz[k] * ab_x + g_y_0_zzzz_xxyyz[k];

                g_y_0_xzzzz_xyzz[k] = -g_y_0_zzzz_xyzz[k] * ab_x + g_y_0_zzzz_xxyzz[k];

                g_y_0_xzzzz_xzzz[k] = -g_y_0_zzzz_xzzz[k] * ab_x + g_y_0_zzzz_xxzzz[k];

                g_y_0_xzzzz_yyyy[k] = -g_y_0_zzzz_yyyy[k] * ab_x + g_y_0_zzzz_xyyyy[k];

                g_y_0_xzzzz_yyyz[k] = -g_y_0_zzzz_yyyz[k] * ab_x + g_y_0_zzzz_xyyyz[k];

                g_y_0_xzzzz_yyzz[k] = -g_y_0_zzzz_yyzz[k] * ab_x + g_y_0_zzzz_xyyzz[k];

                g_y_0_xzzzz_yzzz[k] = -g_y_0_zzzz_yzzz[k] * ab_x + g_y_0_zzzz_xyzzz[k];

                g_y_0_xzzzz_zzzz[k] = -g_y_0_zzzz_zzzz[k] * ab_x + g_y_0_zzzz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzz, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzzz, g_yyyy_xxxx, g_yyyy_xxxy, g_yyyy_xxxz, g_yyyy_xxyy, g_yyyy_xxyz, g_yyyy_xxzz, g_yyyy_xyyy, g_yyyy_xyyz, g_yyyy_xyzz, g_yyyy_xzzz, g_yyyy_yyyy, g_yyyy_yyyz, g_yyyy_yyzz, g_yyyy_yzzz, g_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxx[k] = -g_yyyy_xxxx[k] - g_y_0_yyyy_xxxx[k] * ab_y + g_y_0_yyyy_xxxxy[k];

                g_y_0_yyyyy_xxxy[k] = -g_yyyy_xxxy[k] - g_y_0_yyyy_xxxy[k] * ab_y + g_y_0_yyyy_xxxyy[k];

                g_y_0_yyyyy_xxxz[k] = -g_yyyy_xxxz[k] - g_y_0_yyyy_xxxz[k] * ab_y + g_y_0_yyyy_xxxyz[k];

                g_y_0_yyyyy_xxyy[k] = -g_yyyy_xxyy[k] - g_y_0_yyyy_xxyy[k] * ab_y + g_y_0_yyyy_xxyyy[k];

                g_y_0_yyyyy_xxyz[k] = -g_yyyy_xxyz[k] - g_y_0_yyyy_xxyz[k] * ab_y + g_y_0_yyyy_xxyyz[k];

                g_y_0_yyyyy_xxzz[k] = -g_yyyy_xxzz[k] - g_y_0_yyyy_xxzz[k] * ab_y + g_y_0_yyyy_xxyzz[k];

                g_y_0_yyyyy_xyyy[k] = -g_yyyy_xyyy[k] - g_y_0_yyyy_xyyy[k] * ab_y + g_y_0_yyyy_xyyyy[k];

                g_y_0_yyyyy_xyyz[k] = -g_yyyy_xyyz[k] - g_y_0_yyyy_xyyz[k] * ab_y + g_y_0_yyyy_xyyyz[k];

                g_y_0_yyyyy_xyzz[k] = -g_yyyy_xyzz[k] - g_y_0_yyyy_xyzz[k] * ab_y + g_y_0_yyyy_xyyzz[k];

                g_y_0_yyyyy_xzzz[k] = -g_yyyy_xzzz[k] - g_y_0_yyyy_xzzz[k] * ab_y + g_y_0_yyyy_xyzzz[k];

                g_y_0_yyyyy_yyyy[k] = -g_yyyy_yyyy[k] - g_y_0_yyyy_yyyy[k] * ab_y + g_y_0_yyyy_yyyyy[k];

                g_y_0_yyyyy_yyyz[k] = -g_yyyy_yyyz[k] - g_y_0_yyyy_yyyz[k] * ab_y + g_y_0_yyyy_yyyyz[k];

                g_y_0_yyyyy_yyzz[k] = -g_yyyy_yyzz[k] - g_y_0_yyyy_yyzz[k] * ab_y + g_y_0_yyyy_yyyzz[k];

                g_y_0_yyyyy_yzzz[k] = -g_yyyy_yzzz[k] - g_y_0_yyyy_yzzz[k] * ab_y + g_y_0_yyyy_yyzzz[k];

                g_y_0_yyyyy_zzzz[k] = -g_yyyy_zzzz[k] - g_y_0_yyyy_zzzz[k] * ab_y + g_y_0_yyyy_yzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxx[k] = -g_y_0_yyyy_xxxx[k] * ab_z + g_y_0_yyyy_xxxxz[k];

                g_y_0_yyyyz_xxxy[k] = -g_y_0_yyyy_xxxy[k] * ab_z + g_y_0_yyyy_xxxyz[k];

                g_y_0_yyyyz_xxxz[k] = -g_y_0_yyyy_xxxz[k] * ab_z + g_y_0_yyyy_xxxzz[k];

                g_y_0_yyyyz_xxyy[k] = -g_y_0_yyyy_xxyy[k] * ab_z + g_y_0_yyyy_xxyyz[k];

                g_y_0_yyyyz_xxyz[k] = -g_y_0_yyyy_xxyz[k] * ab_z + g_y_0_yyyy_xxyzz[k];

                g_y_0_yyyyz_xxzz[k] = -g_y_0_yyyy_xxzz[k] * ab_z + g_y_0_yyyy_xxzzz[k];

                g_y_0_yyyyz_xyyy[k] = -g_y_0_yyyy_xyyy[k] * ab_z + g_y_0_yyyy_xyyyz[k];

                g_y_0_yyyyz_xyyz[k] = -g_y_0_yyyy_xyyz[k] * ab_z + g_y_0_yyyy_xyyzz[k];

                g_y_0_yyyyz_xyzz[k] = -g_y_0_yyyy_xyzz[k] * ab_z + g_y_0_yyyy_xyzzz[k];

                g_y_0_yyyyz_xzzz[k] = -g_y_0_yyyy_xzzz[k] * ab_z + g_y_0_yyyy_xzzzz[k];

                g_y_0_yyyyz_yyyy[k] = -g_y_0_yyyy_yyyy[k] * ab_z + g_y_0_yyyy_yyyyz[k];

                g_y_0_yyyyz_yyyz[k] = -g_y_0_yyyy_yyyz[k] * ab_z + g_y_0_yyyy_yyyzz[k];

                g_y_0_yyyyz_yyzz[k] = -g_y_0_yyyy_yyzz[k] * ab_z + g_y_0_yyyy_yyzzz[k];

                g_y_0_yyyyz_yzzz[k] = -g_y_0_yyyy_yzzz[k] * ab_z + g_y_0_yyyy_yzzzz[k];

                g_y_0_yyyyz_zzzz[k] = -g_y_0_yyyy_zzzz[k] * ab_z + g_y_0_yyyy_zzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_zzzz, g_y_0_yyyz_zzzzz, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxx[k] = -g_y_0_yyyz_xxxx[k] * ab_z + g_y_0_yyyz_xxxxz[k];

                g_y_0_yyyzz_xxxy[k] = -g_y_0_yyyz_xxxy[k] * ab_z + g_y_0_yyyz_xxxyz[k];

                g_y_0_yyyzz_xxxz[k] = -g_y_0_yyyz_xxxz[k] * ab_z + g_y_0_yyyz_xxxzz[k];

                g_y_0_yyyzz_xxyy[k] = -g_y_0_yyyz_xxyy[k] * ab_z + g_y_0_yyyz_xxyyz[k];

                g_y_0_yyyzz_xxyz[k] = -g_y_0_yyyz_xxyz[k] * ab_z + g_y_0_yyyz_xxyzz[k];

                g_y_0_yyyzz_xxzz[k] = -g_y_0_yyyz_xxzz[k] * ab_z + g_y_0_yyyz_xxzzz[k];

                g_y_0_yyyzz_xyyy[k] = -g_y_0_yyyz_xyyy[k] * ab_z + g_y_0_yyyz_xyyyz[k];

                g_y_0_yyyzz_xyyz[k] = -g_y_0_yyyz_xyyz[k] * ab_z + g_y_0_yyyz_xyyzz[k];

                g_y_0_yyyzz_xyzz[k] = -g_y_0_yyyz_xyzz[k] * ab_z + g_y_0_yyyz_xyzzz[k];

                g_y_0_yyyzz_xzzz[k] = -g_y_0_yyyz_xzzz[k] * ab_z + g_y_0_yyyz_xzzzz[k];

                g_y_0_yyyzz_yyyy[k] = -g_y_0_yyyz_yyyy[k] * ab_z + g_y_0_yyyz_yyyyz[k];

                g_y_0_yyyzz_yyyz[k] = -g_y_0_yyyz_yyyz[k] * ab_z + g_y_0_yyyz_yyyzz[k];

                g_y_0_yyyzz_yyzz[k] = -g_y_0_yyyz_yyzz[k] * ab_z + g_y_0_yyyz_yyzzz[k];

                g_y_0_yyyzz_yzzz[k] = -g_y_0_yyyz_yzzz[k] * ab_z + g_y_0_yyyz_yzzzz[k];

                g_y_0_yyyzz_zzzz[k] = -g_y_0_yyyz_zzzz[k] * ab_z + g_y_0_yyyz_zzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_zzzz, g_y_0_yyzz_zzzzz, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxx[k] = -g_y_0_yyzz_xxxx[k] * ab_z + g_y_0_yyzz_xxxxz[k];

                g_y_0_yyzzz_xxxy[k] = -g_y_0_yyzz_xxxy[k] * ab_z + g_y_0_yyzz_xxxyz[k];

                g_y_0_yyzzz_xxxz[k] = -g_y_0_yyzz_xxxz[k] * ab_z + g_y_0_yyzz_xxxzz[k];

                g_y_0_yyzzz_xxyy[k] = -g_y_0_yyzz_xxyy[k] * ab_z + g_y_0_yyzz_xxyyz[k];

                g_y_0_yyzzz_xxyz[k] = -g_y_0_yyzz_xxyz[k] * ab_z + g_y_0_yyzz_xxyzz[k];

                g_y_0_yyzzz_xxzz[k] = -g_y_0_yyzz_xxzz[k] * ab_z + g_y_0_yyzz_xxzzz[k];

                g_y_0_yyzzz_xyyy[k] = -g_y_0_yyzz_xyyy[k] * ab_z + g_y_0_yyzz_xyyyz[k];

                g_y_0_yyzzz_xyyz[k] = -g_y_0_yyzz_xyyz[k] * ab_z + g_y_0_yyzz_xyyzz[k];

                g_y_0_yyzzz_xyzz[k] = -g_y_0_yyzz_xyzz[k] * ab_z + g_y_0_yyzz_xyzzz[k];

                g_y_0_yyzzz_xzzz[k] = -g_y_0_yyzz_xzzz[k] * ab_z + g_y_0_yyzz_xzzzz[k];

                g_y_0_yyzzz_yyyy[k] = -g_y_0_yyzz_yyyy[k] * ab_z + g_y_0_yyzz_yyyyz[k];

                g_y_0_yyzzz_yyyz[k] = -g_y_0_yyzz_yyyz[k] * ab_z + g_y_0_yyzz_yyyzz[k];

                g_y_0_yyzzz_yyzz[k] = -g_y_0_yyzz_yyzz[k] * ab_z + g_y_0_yyzz_yyzzz[k];

                g_y_0_yyzzz_yzzz[k] = -g_y_0_yyzz_yzzz[k] * ab_z + g_y_0_yyzz_yzzzz[k];

                g_y_0_yyzzz_zzzz[k] = -g_y_0_yyzz_zzzz[k] * ab_z + g_y_0_yyzz_zzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_zzzz, g_y_0_yzzz_zzzzz, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxx[k] = -g_y_0_yzzz_xxxx[k] * ab_z + g_y_0_yzzz_xxxxz[k];

                g_y_0_yzzzz_xxxy[k] = -g_y_0_yzzz_xxxy[k] * ab_z + g_y_0_yzzz_xxxyz[k];

                g_y_0_yzzzz_xxxz[k] = -g_y_0_yzzz_xxxz[k] * ab_z + g_y_0_yzzz_xxxzz[k];

                g_y_0_yzzzz_xxyy[k] = -g_y_0_yzzz_xxyy[k] * ab_z + g_y_0_yzzz_xxyyz[k];

                g_y_0_yzzzz_xxyz[k] = -g_y_0_yzzz_xxyz[k] * ab_z + g_y_0_yzzz_xxyzz[k];

                g_y_0_yzzzz_xxzz[k] = -g_y_0_yzzz_xxzz[k] * ab_z + g_y_0_yzzz_xxzzz[k];

                g_y_0_yzzzz_xyyy[k] = -g_y_0_yzzz_xyyy[k] * ab_z + g_y_0_yzzz_xyyyz[k];

                g_y_0_yzzzz_xyyz[k] = -g_y_0_yzzz_xyyz[k] * ab_z + g_y_0_yzzz_xyyzz[k];

                g_y_0_yzzzz_xyzz[k] = -g_y_0_yzzz_xyzz[k] * ab_z + g_y_0_yzzz_xyzzz[k];

                g_y_0_yzzzz_xzzz[k] = -g_y_0_yzzz_xzzz[k] * ab_z + g_y_0_yzzz_xzzzz[k];

                g_y_0_yzzzz_yyyy[k] = -g_y_0_yzzz_yyyy[k] * ab_z + g_y_0_yzzz_yyyyz[k];

                g_y_0_yzzzz_yyyz[k] = -g_y_0_yzzz_yyyz[k] * ab_z + g_y_0_yzzz_yyyzz[k];

                g_y_0_yzzzz_yyzz[k] = -g_y_0_yzzz_yyzz[k] * ab_z + g_y_0_yzzz_yyzzz[k];

                g_y_0_yzzzz_yzzz[k] = -g_y_0_yzzz_yzzz[k] * ab_z + g_y_0_yzzz_yzzzz[k];

                g_y_0_yzzzz_zzzz[k] = -g_y_0_yzzz_zzzz[k] * ab_z + g_y_0_yzzz_zzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_zzzz, g_y_0_zzzz_zzzzz, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxx[k] = -g_y_0_zzzz_xxxx[k] * ab_z + g_y_0_zzzz_xxxxz[k];

                g_y_0_zzzzz_xxxy[k] = -g_y_0_zzzz_xxxy[k] * ab_z + g_y_0_zzzz_xxxyz[k];

                g_y_0_zzzzz_xxxz[k] = -g_y_0_zzzz_xxxz[k] * ab_z + g_y_0_zzzz_xxxzz[k];

                g_y_0_zzzzz_xxyy[k] = -g_y_0_zzzz_xxyy[k] * ab_z + g_y_0_zzzz_xxyyz[k];

                g_y_0_zzzzz_xxyz[k] = -g_y_0_zzzz_xxyz[k] * ab_z + g_y_0_zzzz_xxyzz[k];

                g_y_0_zzzzz_xxzz[k] = -g_y_0_zzzz_xxzz[k] * ab_z + g_y_0_zzzz_xxzzz[k];

                g_y_0_zzzzz_xyyy[k] = -g_y_0_zzzz_xyyy[k] * ab_z + g_y_0_zzzz_xyyyz[k];

                g_y_0_zzzzz_xyyz[k] = -g_y_0_zzzz_xyyz[k] * ab_z + g_y_0_zzzz_xyyzz[k];

                g_y_0_zzzzz_xyzz[k] = -g_y_0_zzzz_xyzz[k] * ab_z + g_y_0_zzzz_xyzzz[k];

                g_y_0_zzzzz_xzzz[k] = -g_y_0_zzzz_xzzz[k] * ab_z + g_y_0_zzzz_xzzzz[k];

                g_y_0_zzzzz_yyyy[k] = -g_y_0_zzzz_yyyy[k] * ab_z + g_y_0_zzzz_yyyyz[k];

                g_y_0_zzzzz_yyyz[k] = -g_y_0_zzzz_yyyz[k] * ab_z + g_y_0_zzzz_yyyzz[k];

                g_y_0_zzzzz_yyzz[k] = -g_y_0_zzzz_yyzz[k] * ab_z + g_y_0_zzzz_yyzzz[k];

                g_y_0_zzzzz_yzzz[k] = -g_y_0_zzzz_yzzz[k] * ab_z + g_y_0_zzzz_yzzzz[k];

                g_y_0_zzzzz_zzzz[k] = -g_y_0_zzzz_zzzz[k] * ab_z + g_y_0_zzzz_zzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxxx, g_z_0_xxxx_xxxxy, g_z_0_xxxx_xxxxz, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxyy, g_z_0_xxxx_xxxyz, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxxzz, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyyy, g_z_0_xxxx_xxyyz, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxyzz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xxzzz, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyyy, g_z_0_xxxx_xyyyz, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyyzz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xyzzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_xzzzz, g_z_0_xxxx_yyyy, g_z_0_xxxx_yyyz, g_z_0_xxxx_yyzz, g_z_0_xxxx_yzzz, g_z_0_xxxx_zzzz, g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_yyyy, g_z_0_xxxxx_yyyz, g_z_0_xxxxx_yyzz, g_z_0_xxxxx_yzzz, g_z_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxx[k] = -g_z_0_xxxx_xxxx[k] * ab_x + g_z_0_xxxx_xxxxx[k];

                g_z_0_xxxxx_xxxy[k] = -g_z_0_xxxx_xxxy[k] * ab_x + g_z_0_xxxx_xxxxy[k];

                g_z_0_xxxxx_xxxz[k] = -g_z_0_xxxx_xxxz[k] * ab_x + g_z_0_xxxx_xxxxz[k];

                g_z_0_xxxxx_xxyy[k] = -g_z_0_xxxx_xxyy[k] * ab_x + g_z_0_xxxx_xxxyy[k];

                g_z_0_xxxxx_xxyz[k] = -g_z_0_xxxx_xxyz[k] * ab_x + g_z_0_xxxx_xxxyz[k];

                g_z_0_xxxxx_xxzz[k] = -g_z_0_xxxx_xxzz[k] * ab_x + g_z_0_xxxx_xxxzz[k];

                g_z_0_xxxxx_xyyy[k] = -g_z_0_xxxx_xyyy[k] * ab_x + g_z_0_xxxx_xxyyy[k];

                g_z_0_xxxxx_xyyz[k] = -g_z_0_xxxx_xyyz[k] * ab_x + g_z_0_xxxx_xxyyz[k];

                g_z_0_xxxxx_xyzz[k] = -g_z_0_xxxx_xyzz[k] * ab_x + g_z_0_xxxx_xxyzz[k];

                g_z_0_xxxxx_xzzz[k] = -g_z_0_xxxx_xzzz[k] * ab_x + g_z_0_xxxx_xxzzz[k];

                g_z_0_xxxxx_yyyy[k] = -g_z_0_xxxx_yyyy[k] * ab_x + g_z_0_xxxx_xyyyy[k];

                g_z_0_xxxxx_yyyz[k] = -g_z_0_xxxx_yyyz[k] * ab_x + g_z_0_xxxx_xyyyz[k];

                g_z_0_xxxxx_yyzz[k] = -g_z_0_xxxx_yyzz[k] * ab_x + g_z_0_xxxx_xyyzz[k];

                g_z_0_xxxxx_yzzz[k] = -g_z_0_xxxx_yzzz[k] * ab_x + g_z_0_xxxx_xyzzz[k];

                g_z_0_xxxxx_zzzz[k] = -g_z_0_xxxx_zzzz[k] * ab_x + g_z_0_xxxx_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_yyyy, g_z_0_xxxxy_yyyz, g_z_0_xxxxy_yyzz, g_z_0_xxxxy_yzzz, g_z_0_xxxxy_zzzz, g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxxx, g_z_0_xxxy_xxxxy, g_z_0_xxxy_xxxxz, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxyy, g_z_0_xxxy_xxxyz, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxxzz, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyyy, g_z_0_xxxy_xxyyz, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxyzz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xxzzz, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyyy, g_z_0_xxxy_xyyyz, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyyzz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xyzzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_xzzzz, g_z_0_xxxy_yyyy, g_z_0_xxxy_yyyz, g_z_0_xxxy_yyzz, g_z_0_xxxy_yzzz, g_z_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxx[k] = -g_z_0_xxxy_xxxx[k] * ab_x + g_z_0_xxxy_xxxxx[k];

                g_z_0_xxxxy_xxxy[k] = -g_z_0_xxxy_xxxy[k] * ab_x + g_z_0_xxxy_xxxxy[k];

                g_z_0_xxxxy_xxxz[k] = -g_z_0_xxxy_xxxz[k] * ab_x + g_z_0_xxxy_xxxxz[k];

                g_z_0_xxxxy_xxyy[k] = -g_z_0_xxxy_xxyy[k] * ab_x + g_z_0_xxxy_xxxyy[k];

                g_z_0_xxxxy_xxyz[k] = -g_z_0_xxxy_xxyz[k] * ab_x + g_z_0_xxxy_xxxyz[k];

                g_z_0_xxxxy_xxzz[k] = -g_z_0_xxxy_xxzz[k] * ab_x + g_z_0_xxxy_xxxzz[k];

                g_z_0_xxxxy_xyyy[k] = -g_z_0_xxxy_xyyy[k] * ab_x + g_z_0_xxxy_xxyyy[k];

                g_z_0_xxxxy_xyyz[k] = -g_z_0_xxxy_xyyz[k] * ab_x + g_z_0_xxxy_xxyyz[k];

                g_z_0_xxxxy_xyzz[k] = -g_z_0_xxxy_xyzz[k] * ab_x + g_z_0_xxxy_xxyzz[k];

                g_z_0_xxxxy_xzzz[k] = -g_z_0_xxxy_xzzz[k] * ab_x + g_z_0_xxxy_xxzzz[k];

                g_z_0_xxxxy_yyyy[k] = -g_z_0_xxxy_yyyy[k] * ab_x + g_z_0_xxxy_xyyyy[k];

                g_z_0_xxxxy_yyyz[k] = -g_z_0_xxxy_yyyz[k] * ab_x + g_z_0_xxxy_xyyyz[k];

                g_z_0_xxxxy_yyzz[k] = -g_z_0_xxxy_yyzz[k] * ab_x + g_z_0_xxxy_xyyzz[k];

                g_z_0_xxxxy_yzzz[k] = -g_z_0_xxxy_yzzz[k] * ab_x + g_z_0_xxxy_xyzzz[k];

                g_z_0_xxxxy_zzzz[k] = -g_z_0_xxxy_zzzz[k] * ab_x + g_z_0_xxxy_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_yyyy, g_z_0_xxxxz_yyyz, g_z_0_xxxxz_yyzz, g_z_0_xxxxz_yzzz, g_z_0_xxxxz_zzzz, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxxx, g_z_0_xxxz_xxxxy, g_z_0_xxxz_xxxxz, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxyy, g_z_0_xxxz_xxxyz, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxxzz, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyyy, g_z_0_xxxz_xxyyz, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxyzz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xxzzz, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyyy, g_z_0_xxxz_xyyyz, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyyzz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xyzzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_xzzzz, g_z_0_xxxz_yyyy, g_z_0_xxxz_yyyz, g_z_0_xxxz_yyzz, g_z_0_xxxz_yzzz, g_z_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxx[k] = -g_z_0_xxxz_xxxx[k] * ab_x + g_z_0_xxxz_xxxxx[k];

                g_z_0_xxxxz_xxxy[k] = -g_z_0_xxxz_xxxy[k] * ab_x + g_z_0_xxxz_xxxxy[k];

                g_z_0_xxxxz_xxxz[k] = -g_z_0_xxxz_xxxz[k] * ab_x + g_z_0_xxxz_xxxxz[k];

                g_z_0_xxxxz_xxyy[k] = -g_z_0_xxxz_xxyy[k] * ab_x + g_z_0_xxxz_xxxyy[k];

                g_z_0_xxxxz_xxyz[k] = -g_z_0_xxxz_xxyz[k] * ab_x + g_z_0_xxxz_xxxyz[k];

                g_z_0_xxxxz_xxzz[k] = -g_z_0_xxxz_xxzz[k] * ab_x + g_z_0_xxxz_xxxzz[k];

                g_z_0_xxxxz_xyyy[k] = -g_z_0_xxxz_xyyy[k] * ab_x + g_z_0_xxxz_xxyyy[k];

                g_z_0_xxxxz_xyyz[k] = -g_z_0_xxxz_xyyz[k] * ab_x + g_z_0_xxxz_xxyyz[k];

                g_z_0_xxxxz_xyzz[k] = -g_z_0_xxxz_xyzz[k] * ab_x + g_z_0_xxxz_xxyzz[k];

                g_z_0_xxxxz_xzzz[k] = -g_z_0_xxxz_xzzz[k] * ab_x + g_z_0_xxxz_xxzzz[k];

                g_z_0_xxxxz_yyyy[k] = -g_z_0_xxxz_yyyy[k] * ab_x + g_z_0_xxxz_xyyyy[k];

                g_z_0_xxxxz_yyyz[k] = -g_z_0_xxxz_yyyz[k] * ab_x + g_z_0_xxxz_xyyyz[k];

                g_z_0_xxxxz_yyzz[k] = -g_z_0_xxxz_yyzz[k] * ab_x + g_z_0_xxxz_xyyzz[k];

                g_z_0_xxxxz_yzzz[k] = -g_z_0_xxxz_yzzz[k] * ab_x + g_z_0_xxxz_xyzzz[k];

                g_z_0_xxxxz_zzzz[k] = -g_z_0_xxxz_zzzz[k] * ab_x + g_z_0_xxxz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_yyyy, g_z_0_xxxyy_yyyz, g_z_0_xxxyy_yyzz, g_z_0_xxxyy_yzzz, g_z_0_xxxyy_zzzz, g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxxx, g_z_0_xxyy_xxxxy, g_z_0_xxyy_xxxxz, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxyy, g_z_0_xxyy_xxxyz, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxxzz, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyyy, g_z_0_xxyy_xxyyz, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxyzz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xxzzz, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyyy, g_z_0_xxyy_xyyyz, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyyzz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xyzzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_xzzzz, g_z_0_xxyy_yyyy, g_z_0_xxyy_yyyz, g_z_0_xxyy_yyzz, g_z_0_xxyy_yzzz, g_z_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxx[k] = -g_z_0_xxyy_xxxx[k] * ab_x + g_z_0_xxyy_xxxxx[k];

                g_z_0_xxxyy_xxxy[k] = -g_z_0_xxyy_xxxy[k] * ab_x + g_z_0_xxyy_xxxxy[k];

                g_z_0_xxxyy_xxxz[k] = -g_z_0_xxyy_xxxz[k] * ab_x + g_z_0_xxyy_xxxxz[k];

                g_z_0_xxxyy_xxyy[k] = -g_z_0_xxyy_xxyy[k] * ab_x + g_z_0_xxyy_xxxyy[k];

                g_z_0_xxxyy_xxyz[k] = -g_z_0_xxyy_xxyz[k] * ab_x + g_z_0_xxyy_xxxyz[k];

                g_z_0_xxxyy_xxzz[k] = -g_z_0_xxyy_xxzz[k] * ab_x + g_z_0_xxyy_xxxzz[k];

                g_z_0_xxxyy_xyyy[k] = -g_z_0_xxyy_xyyy[k] * ab_x + g_z_0_xxyy_xxyyy[k];

                g_z_0_xxxyy_xyyz[k] = -g_z_0_xxyy_xyyz[k] * ab_x + g_z_0_xxyy_xxyyz[k];

                g_z_0_xxxyy_xyzz[k] = -g_z_0_xxyy_xyzz[k] * ab_x + g_z_0_xxyy_xxyzz[k];

                g_z_0_xxxyy_xzzz[k] = -g_z_0_xxyy_xzzz[k] * ab_x + g_z_0_xxyy_xxzzz[k];

                g_z_0_xxxyy_yyyy[k] = -g_z_0_xxyy_yyyy[k] * ab_x + g_z_0_xxyy_xyyyy[k];

                g_z_0_xxxyy_yyyz[k] = -g_z_0_xxyy_yyyz[k] * ab_x + g_z_0_xxyy_xyyyz[k];

                g_z_0_xxxyy_yyzz[k] = -g_z_0_xxyy_yyzz[k] * ab_x + g_z_0_xxyy_xyyzz[k];

                g_z_0_xxxyy_yzzz[k] = -g_z_0_xxyy_yzzz[k] * ab_x + g_z_0_xxyy_xyzzz[k];

                g_z_0_xxxyy_zzzz[k] = -g_z_0_xxyy_zzzz[k] * ab_x + g_z_0_xxyy_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_yyyy, g_z_0_xxxyz_yyyz, g_z_0_xxxyz_yyzz, g_z_0_xxxyz_yzzz, g_z_0_xxxyz_zzzz, g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxxx, g_z_0_xxyz_xxxxy, g_z_0_xxyz_xxxxz, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxyy, g_z_0_xxyz_xxxyz, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxxzz, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyyy, g_z_0_xxyz_xxyyz, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxyzz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xxzzz, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyyy, g_z_0_xxyz_xyyyz, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyyzz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xyzzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_xzzzz, g_z_0_xxyz_yyyy, g_z_0_xxyz_yyyz, g_z_0_xxyz_yyzz, g_z_0_xxyz_yzzz, g_z_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxx[k] = -g_z_0_xxyz_xxxx[k] * ab_x + g_z_0_xxyz_xxxxx[k];

                g_z_0_xxxyz_xxxy[k] = -g_z_0_xxyz_xxxy[k] * ab_x + g_z_0_xxyz_xxxxy[k];

                g_z_0_xxxyz_xxxz[k] = -g_z_0_xxyz_xxxz[k] * ab_x + g_z_0_xxyz_xxxxz[k];

                g_z_0_xxxyz_xxyy[k] = -g_z_0_xxyz_xxyy[k] * ab_x + g_z_0_xxyz_xxxyy[k];

                g_z_0_xxxyz_xxyz[k] = -g_z_0_xxyz_xxyz[k] * ab_x + g_z_0_xxyz_xxxyz[k];

                g_z_0_xxxyz_xxzz[k] = -g_z_0_xxyz_xxzz[k] * ab_x + g_z_0_xxyz_xxxzz[k];

                g_z_0_xxxyz_xyyy[k] = -g_z_0_xxyz_xyyy[k] * ab_x + g_z_0_xxyz_xxyyy[k];

                g_z_0_xxxyz_xyyz[k] = -g_z_0_xxyz_xyyz[k] * ab_x + g_z_0_xxyz_xxyyz[k];

                g_z_0_xxxyz_xyzz[k] = -g_z_0_xxyz_xyzz[k] * ab_x + g_z_0_xxyz_xxyzz[k];

                g_z_0_xxxyz_xzzz[k] = -g_z_0_xxyz_xzzz[k] * ab_x + g_z_0_xxyz_xxzzz[k];

                g_z_0_xxxyz_yyyy[k] = -g_z_0_xxyz_yyyy[k] * ab_x + g_z_0_xxyz_xyyyy[k];

                g_z_0_xxxyz_yyyz[k] = -g_z_0_xxyz_yyyz[k] * ab_x + g_z_0_xxyz_xyyyz[k];

                g_z_0_xxxyz_yyzz[k] = -g_z_0_xxyz_yyzz[k] * ab_x + g_z_0_xxyz_xyyzz[k];

                g_z_0_xxxyz_yzzz[k] = -g_z_0_xxyz_yzzz[k] * ab_x + g_z_0_xxyz_xyzzz[k];

                g_z_0_xxxyz_zzzz[k] = -g_z_0_xxyz_zzzz[k] * ab_x + g_z_0_xxyz_xzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_yyyy, g_z_0_xxxzz_yyyz, g_z_0_xxxzz_yyzz, g_z_0_xxxzz_yzzz, g_z_0_xxxzz_zzzz, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxxx, g_z_0_xxzz_xxxxy, g_z_0_xxzz_xxxxz, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxyy, g_z_0_xxzz_xxxyz, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxxzz, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyyy, g_z_0_xxzz_xxyyz, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxyzz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xxzzz, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyyy, g_z_0_xxzz_xyyyz, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyyzz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xyzzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_xzzzz, g_z_0_xxzz_yyyy, g_z_0_xxzz_yyyz, g_z_0_xxzz_yyzz, g_z_0_xxzz_yzzz, g_z_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxx[k] = -g_z_0_xxzz_xxxx[k] * ab_x + g_z_0_xxzz_xxxxx[k];

                g_z_0_xxxzz_xxxy[k] = -g_z_0_xxzz_xxxy[k] * ab_x + g_z_0_xxzz_xxxxy[k];

                g_z_0_xxxzz_xxxz[k] = -g_z_0_xxzz_xxxz[k] * ab_x + g_z_0_xxzz_xxxxz[k];

                g_z_0_xxxzz_xxyy[k] = -g_z_0_xxzz_xxyy[k] * ab_x + g_z_0_xxzz_xxxyy[k];

                g_z_0_xxxzz_xxyz[k] = -g_z_0_xxzz_xxyz[k] * ab_x + g_z_0_xxzz_xxxyz[k];

                g_z_0_xxxzz_xxzz[k] = -g_z_0_xxzz_xxzz[k] * ab_x + g_z_0_xxzz_xxxzz[k];

                g_z_0_xxxzz_xyyy[k] = -g_z_0_xxzz_xyyy[k] * ab_x + g_z_0_xxzz_xxyyy[k];

                g_z_0_xxxzz_xyyz[k] = -g_z_0_xxzz_xyyz[k] * ab_x + g_z_0_xxzz_xxyyz[k];

                g_z_0_xxxzz_xyzz[k] = -g_z_0_xxzz_xyzz[k] * ab_x + g_z_0_xxzz_xxyzz[k];

                g_z_0_xxxzz_xzzz[k] = -g_z_0_xxzz_xzzz[k] * ab_x + g_z_0_xxzz_xxzzz[k];

                g_z_0_xxxzz_yyyy[k] = -g_z_0_xxzz_yyyy[k] * ab_x + g_z_0_xxzz_xyyyy[k];

                g_z_0_xxxzz_yyyz[k] = -g_z_0_xxzz_yyyz[k] * ab_x + g_z_0_xxzz_xyyyz[k];

                g_z_0_xxxzz_yyzz[k] = -g_z_0_xxzz_yyzz[k] * ab_x + g_z_0_xxzz_xyyzz[k];

                g_z_0_xxxzz_yzzz[k] = -g_z_0_xxzz_yzzz[k] * ab_x + g_z_0_xxzz_xyzzz[k];

                g_z_0_xxxzz_zzzz[k] = -g_z_0_xxzz_zzzz[k] * ab_x + g_z_0_xxzz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_yyyy, g_z_0_xxyyy_yyyz, g_z_0_xxyyy_yyzz, g_z_0_xxyyy_yzzz, g_z_0_xxyyy_zzzz, g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxxx, g_z_0_xyyy_xxxxy, g_z_0_xyyy_xxxxz, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxyy, g_z_0_xyyy_xxxyz, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxxzz, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyyy, g_z_0_xyyy_xxyyz, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxyzz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xxzzz, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyyy, g_z_0_xyyy_xyyyz, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyyzz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xyzzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_xzzzz, g_z_0_xyyy_yyyy, g_z_0_xyyy_yyyz, g_z_0_xyyy_yyzz, g_z_0_xyyy_yzzz, g_z_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxx[k] = -g_z_0_xyyy_xxxx[k] * ab_x + g_z_0_xyyy_xxxxx[k];

                g_z_0_xxyyy_xxxy[k] = -g_z_0_xyyy_xxxy[k] * ab_x + g_z_0_xyyy_xxxxy[k];

                g_z_0_xxyyy_xxxz[k] = -g_z_0_xyyy_xxxz[k] * ab_x + g_z_0_xyyy_xxxxz[k];

                g_z_0_xxyyy_xxyy[k] = -g_z_0_xyyy_xxyy[k] * ab_x + g_z_0_xyyy_xxxyy[k];

                g_z_0_xxyyy_xxyz[k] = -g_z_0_xyyy_xxyz[k] * ab_x + g_z_0_xyyy_xxxyz[k];

                g_z_0_xxyyy_xxzz[k] = -g_z_0_xyyy_xxzz[k] * ab_x + g_z_0_xyyy_xxxzz[k];

                g_z_0_xxyyy_xyyy[k] = -g_z_0_xyyy_xyyy[k] * ab_x + g_z_0_xyyy_xxyyy[k];

                g_z_0_xxyyy_xyyz[k] = -g_z_0_xyyy_xyyz[k] * ab_x + g_z_0_xyyy_xxyyz[k];

                g_z_0_xxyyy_xyzz[k] = -g_z_0_xyyy_xyzz[k] * ab_x + g_z_0_xyyy_xxyzz[k];

                g_z_0_xxyyy_xzzz[k] = -g_z_0_xyyy_xzzz[k] * ab_x + g_z_0_xyyy_xxzzz[k];

                g_z_0_xxyyy_yyyy[k] = -g_z_0_xyyy_yyyy[k] * ab_x + g_z_0_xyyy_xyyyy[k];

                g_z_0_xxyyy_yyyz[k] = -g_z_0_xyyy_yyyz[k] * ab_x + g_z_0_xyyy_xyyyz[k];

                g_z_0_xxyyy_yyzz[k] = -g_z_0_xyyy_yyzz[k] * ab_x + g_z_0_xyyy_xyyzz[k];

                g_z_0_xxyyy_yzzz[k] = -g_z_0_xyyy_yzzz[k] * ab_x + g_z_0_xyyy_xyzzz[k];

                g_z_0_xxyyy_zzzz[k] = -g_z_0_xyyy_zzzz[k] * ab_x + g_z_0_xyyy_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_yyyy, g_z_0_xxyyz_yyyz, g_z_0_xxyyz_yyzz, g_z_0_xxyyz_yzzz, g_z_0_xxyyz_zzzz, g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxxx, g_z_0_xyyz_xxxxy, g_z_0_xyyz_xxxxz, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxyy, g_z_0_xyyz_xxxyz, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxxzz, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyyy, g_z_0_xyyz_xxyyz, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxyzz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xxzzz, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyyy, g_z_0_xyyz_xyyyz, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyyzz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xyzzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_xzzzz, g_z_0_xyyz_yyyy, g_z_0_xyyz_yyyz, g_z_0_xyyz_yyzz, g_z_0_xyyz_yzzz, g_z_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxx[k] = -g_z_0_xyyz_xxxx[k] * ab_x + g_z_0_xyyz_xxxxx[k];

                g_z_0_xxyyz_xxxy[k] = -g_z_0_xyyz_xxxy[k] * ab_x + g_z_0_xyyz_xxxxy[k];

                g_z_0_xxyyz_xxxz[k] = -g_z_0_xyyz_xxxz[k] * ab_x + g_z_0_xyyz_xxxxz[k];

                g_z_0_xxyyz_xxyy[k] = -g_z_0_xyyz_xxyy[k] * ab_x + g_z_0_xyyz_xxxyy[k];

                g_z_0_xxyyz_xxyz[k] = -g_z_0_xyyz_xxyz[k] * ab_x + g_z_0_xyyz_xxxyz[k];

                g_z_0_xxyyz_xxzz[k] = -g_z_0_xyyz_xxzz[k] * ab_x + g_z_0_xyyz_xxxzz[k];

                g_z_0_xxyyz_xyyy[k] = -g_z_0_xyyz_xyyy[k] * ab_x + g_z_0_xyyz_xxyyy[k];

                g_z_0_xxyyz_xyyz[k] = -g_z_0_xyyz_xyyz[k] * ab_x + g_z_0_xyyz_xxyyz[k];

                g_z_0_xxyyz_xyzz[k] = -g_z_0_xyyz_xyzz[k] * ab_x + g_z_0_xyyz_xxyzz[k];

                g_z_0_xxyyz_xzzz[k] = -g_z_0_xyyz_xzzz[k] * ab_x + g_z_0_xyyz_xxzzz[k];

                g_z_0_xxyyz_yyyy[k] = -g_z_0_xyyz_yyyy[k] * ab_x + g_z_0_xyyz_xyyyy[k];

                g_z_0_xxyyz_yyyz[k] = -g_z_0_xyyz_yyyz[k] * ab_x + g_z_0_xyyz_xyyyz[k];

                g_z_0_xxyyz_yyzz[k] = -g_z_0_xyyz_yyzz[k] * ab_x + g_z_0_xyyz_xyyzz[k];

                g_z_0_xxyyz_yzzz[k] = -g_z_0_xyyz_yzzz[k] * ab_x + g_z_0_xyyz_xyzzz[k];

                g_z_0_xxyyz_zzzz[k] = -g_z_0_xyyz_zzzz[k] * ab_x + g_z_0_xyyz_xzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_yyyy, g_z_0_xxyzz_yyyz, g_z_0_xxyzz_yyzz, g_z_0_xxyzz_yzzz, g_z_0_xxyzz_zzzz, g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxxx, g_z_0_xyzz_xxxxy, g_z_0_xyzz_xxxxz, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxyy, g_z_0_xyzz_xxxyz, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxxzz, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyyy, g_z_0_xyzz_xxyyz, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxyzz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xxzzz, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyyy, g_z_0_xyzz_xyyyz, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyyzz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xyzzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_xzzzz, g_z_0_xyzz_yyyy, g_z_0_xyzz_yyyz, g_z_0_xyzz_yyzz, g_z_0_xyzz_yzzz, g_z_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxx[k] = -g_z_0_xyzz_xxxx[k] * ab_x + g_z_0_xyzz_xxxxx[k];

                g_z_0_xxyzz_xxxy[k] = -g_z_0_xyzz_xxxy[k] * ab_x + g_z_0_xyzz_xxxxy[k];

                g_z_0_xxyzz_xxxz[k] = -g_z_0_xyzz_xxxz[k] * ab_x + g_z_0_xyzz_xxxxz[k];

                g_z_0_xxyzz_xxyy[k] = -g_z_0_xyzz_xxyy[k] * ab_x + g_z_0_xyzz_xxxyy[k];

                g_z_0_xxyzz_xxyz[k] = -g_z_0_xyzz_xxyz[k] * ab_x + g_z_0_xyzz_xxxyz[k];

                g_z_0_xxyzz_xxzz[k] = -g_z_0_xyzz_xxzz[k] * ab_x + g_z_0_xyzz_xxxzz[k];

                g_z_0_xxyzz_xyyy[k] = -g_z_0_xyzz_xyyy[k] * ab_x + g_z_0_xyzz_xxyyy[k];

                g_z_0_xxyzz_xyyz[k] = -g_z_0_xyzz_xyyz[k] * ab_x + g_z_0_xyzz_xxyyz[k];

                g_z_0_xxyzz_xyzz[k] = -g_z_0_xyzz_xyzz[k] * ab_x + g_z_0_xyzz_xxyzz[k];

                g_z_0_xxyzz_xzzz[k] = -g_z_0_xyzz_xzzz[k] * ab_x + g_z_0_xyzz_xxzzz[k];

                g_z_0_xxyzz_yyyy[k] = -g_z_0_xyzz_yyyy[k] * ab_x + g_z_0_xyzz_xyyyy[k];

                g_z_0_xxyzz_yyyz[k] = -g_z_0_xyzz_yyyz[k] * ab_x + g_z_0_xyzz_xyyyz[k];

                g_z_0_xxyzz_yyzz[k] = -g_z_0_xyzz_yyzz[k] * ab_x + g_z_0_xyzz_xyyzz[k];

                g_z_0_xxyzz_yzzz[k] = -g_z_0_xyzz_yzzz[k] * ab_x + g_z_0_xyzz_xyzzz[k];

                g_z_0_xxyzz_zzzz[k] = -g_z_0_xyzz_zzzz[k] * ab_x + g_z_0_xyzz_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_yyyy, g_z_0_xxzzz_yyyz, g_z_0_xxzzz_yyzz, g_z_0_xxzzz_yzzz, g_z_0_xxzzz_zzzz, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxxx, g_z_0_xzzz_xxxxy, g_z_0_xzzz_xxxxz, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxyy, g_z_0_xzzz_xxxyz, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxxzz, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyyy, g_z_0_xzzz_xxyyz, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxyzz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xxzzz, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyyy, g_z_0_xzzz_xyyyz, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyyzz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xyzzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_xzzzz, g_z_0_xzzz_yyyy, g_z_0_xzzz_yyyz, g_z_0_xzzz_yyzz, g_z_0_xzzz_yzzz, g_z_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxx[k] = -g_z_0_xzzz_xxxx[k] * ab_x + g_z_0_xzzz_xxxxx[k];

                g_z_0_xxzzz_xxxy[k] = -g_z_0_xzzz_xxxy[k] * ab_x + g_z_0_xzzz_xxxxy[k];

                g_z_0_xxzzz_xxxz[k] = -g_z_0_xzzz_xxxz[k] * ab_x + g_z_0_xzzz_xxxxz[k];

                g_z_0_xxzzz_xxyy[k] = -g_z_0_xzzz_xxyy[k] * ab_x + g_z_0_xzzz_xxxyy[k];

                g_z_0_xxzzz_xxyz[k] = -g_z_0_xzzz_xxyz[k] * ab_x + g_z_0_xzzz_xxxyz[k];

                g_z_0_xxzzz_xxzz[k] = -g_z_0_xzzz_xxzz[k] * ab_x + g_z_0_xzzz_xxxzz[k];

                g_z_0_xxzzz_xyyy[k] = -g_z_0_xzzz_xyyy[k] * ab_x + g_z_0_xzzz_xxyyy[k];

                g_z_0_xxzzz_xyyz[k] = -g_z_0_xzzz_xyyz[k] * ab_x + g_z_0_xzzz_xxyyz[k];

                g_z_0_xxzzz_xyzz[k] = -g_z_0_xzzz_xyzz[k] * ab_x + g_z_0_xzzz_xxyzz[k];

                g_z_0_xxzzz_xzzz[k] = -g_z_0_xzzz_xzzz[k] * ab_x + g_z_0_xzzz_xxzzz[k];

                g_z_0_xxzzz_yyyy[k] = -g_z_0_xzzz_yyyy[k] * ab_x + g_z_0_xzzz_xyyyy[k];

                g_z_0_xxzzz_yyyz[k] = -g_z_0_xzzz_yyyz[k] * ab_x + g_z_0_xzzz_xyyyz[k];

                g_z_0_xxzzz_yyzz[k] = -g_z_0_xzzz_yyzz[k] * ab_x + g_z_0_xzzz_xyyzz[k];

                g_z_0_xxzzz_yzzz[k] = -g_z_0_xzzz_yzzz[k] * ab_x + g_z_0_xzzz_xyzzz[k];

                g_z_0_xxzzz_zzzz[k] = -g_z_0_xzzz_zzzz[k] * ab_x + g_z_0_xzzz_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_yyyy, g_z_0_xyyyy_yyyz, g_z_0_xyyyy_yyzz, g_z_0_xyyyy_yzzz, g_z_0_xyyyy_zzzz, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxx[k] = -g_z_0_yyyy_xxxx[k] * ab_x + g_z_0_yyyy_xxxxx[k];

                g_z_0_xyyyy_xxxy[k] = -g_z_0_yyyy_xxxy[k] * ab_x + g_z_0_yyyy_xxxxy[k];

                g_z_0_xyyyy_xxxz[k] = -g_z_0_yyyy_xxxz[k] * ab_x + g_z_0_yyyy_xxxxz[k];

                g_z_0_xyyyy_xxyy[k] = -g_z_0_yyyy_xxyy[k] * ab_x + g_z_0_yyyy_xxxyy[k];

                g_z_0_xyyyy_xxyz[k] = -g_z_0_yyyy_xxyz[k] * ab_x + g_z_0_yyyy_xxxyz[k];

                g_z_0_xyyyy_xxzz[k] = -g_z_0_yyyy_xxzz[k] * ab_x + g_z_0_yyyy_xxxzz[k];

                g_z_0_xyyyy_xyyy[k] = -g_z_0_yyyy_xyyy[k] * ab_x + g_z_0_yyyy_xxyyy[k];

                g_z_0_xyyyy_xyyz[k] = -g_z_0_yyyy_xyyz[k] * ab_x + g_z_0_yyyy_xxyyz[k];

                g_z_0_xyyyy_xyzz[k] = -g_z_0_yyyy_xyzz[k] * ab_x + g_z_0_yyyy_xxyzz[k];

                g_z_0_xyyyy_xzzz[k] = -g_z_0_yyyy_xzzz[k] * ab_x + g_z_0_yyyy_xxzzz[k];

                g_z_0_xyyyy_yyyy[k] = -g_z_0_yyyy_yyyy[k] * ab_x + g_z_0_yyyy_xyyyy[k];

                g_z_0_xyyyy_yyyz[k] = -g_z_0_yyyy_yyyz[k] * ab_x + g_z_0_yyyy_xyyyz[k];

                g_z_0_xyyyy_yyzz[k] = -g_z_0_yyyy_yyzz[k] * ab_x + g_z_0_yyyy_xyyzz[k];

                g_z_0_xyyyy_yzzz[k] = -g_z_0_yyyy_yzzz[k] * ab_x + g_z_0_yyyy_xyzzz[k];

                g_z_0_xyyyy_zzzz[k] = -g_z_0_yyyy_zzzz[k] * ab_x + g_z_0_yyyy_xzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_yyyy, g_z_0_xyyyz_yyyz, g_z_0_xyyyz_yyzz, g_z_0_xyyyz_yzzz, g_z_0_xyyyz_zzzz, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxx[k] = -g_z_0_yyyz_xxxx[k] * ab_x + g_z_0_yyyz_xxxxx[k];

                g_z_0_xyyyz_xxxy[k] = -g_z_0_yyyz_xxxy[k] * ab_x + g_z_0_yyyz_xxxxy[k];

                g_z_0_xyyyz_xxxz[k] = -g_z_0_yyyz_xxxz[k] * ab_x + g_z_0_yyyz_xxxxz[k];

                g_z_0_xyyyz_xxyy[k] = -g_z_0_yyyz_xxyy[k] * ab_x + g_z_0_yyyz_xxxyy[k];

                g_z_0_xyyyz_xxyz[k] = -g_z_0_yyyz_xxyz[k] * ab_x + g_z_0_yyyz_xxxyz[k];

                g_z_0_xyyyz_xxzz[k] = -g_z_0_yyyz_xxzz[k] * ab_x + g_z_0_yyyz_xxxzz[k];

                g_z_0_xyyyz_xyyy[k] = -g_z_0_yyyz_xyyy[k] * ab_x + g_z_0_yyyz_xxyyy[k];

                g_z_0_xyyyz_xyyz[k] = -g_z_0_yyyz_xyyz[k] * ab_x + g_z_0_yyyz_xxyyz[k];

                g_z_0_xyyyz_xyzz[k] = -g_z_0_yyyz_xyzz[k] * ab_x + g_z_0_yyyz_xxyzz[k];

                g_z_0_xyyyz_xzzz[k] = -g_z_0_yyyz_xzzz[k] * ab_x + g_z_0_yyyz_xxzzz[k];

                g_z_0_xyyyz_yyyy[k] = -g_z_0_yyyz_yyyy[k] * ab_x + g_z_0_yyyz_xyyyy[k];

                g_z_0_xyyyz_yyyz[k] = -g_z_0_yyyz_yyyz[k] * ab_x + g_z_0_yyyz_xyyyz[k];

                g_z_0_xyyyz_yyzz[k] = -g_z_0_yyyz_yyzz[k] * ab_x + g_z_0_yyyz_xyyzz[k];

                g_z_0_xyyyz_yzzz[k] = -g_z_0_yyyz_yzzz[k] * ab_x + g_z_0_yyyz_xyzzz[k];

                g_z_0_xyyyz_zzzz[k] = -g_z_0_yyyz_zzzz[k] * ab_x + g_z_0_yyyz_xzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_yyyy, g_z_0_xyyzz_yyyz, g_z_0_xyyzz_yyzz, g_z_0_xyyzz_yzzz, g_z_0_xyyzz_zzzz, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxx[k] = -g_z_0_yyzz_xxxx[k] * ab_x + g_z_0_yyzz_xxxxx[k];

                g_z_0_xyyzz_xxxy[k] = -g_z_0_yyzz_xxxy[k] * ab_x + g_z_0_yyzz_xxxxy[k];

                g_z_0_xyyzz_xxxz[k] = -g_z_0_yyzz_xxxz[k] * ab_x + g_z_0_yyzz_xxxxz[k];

                g_z_0_xyyzz_xxyy[k] = -g_z_0_yyzz_xxyy[k] * ab_x + g_z_0_yyzz_xxxyy[k];

                g_z_0_xyyzz_xxyz[k] = -g_z_0_yyzz_xxyz[k] * ab_x + g_z_0_yyzz_xxxyz[k];

                g_z_0_xyyzz_xxzz[k] = -g_z_0_yyzz_xxzz[k] * ab_x + g_z_0_yyzz_xxxzz[k];

                g_z_0_xyyzz_xyyy[k] = -g_z_0_yyzz_xyyy[k] * ab_x + g_z_0_yyzz_xxyyy[k];

                g_z_0_xyyzz_xyyz[k] = -g_z_0_yyzz_xyyz[k] * ab_x + g_z_0_yyzz_xxyyz[k];

                g_z_0_xyyzz_xyzz[k] = -g_z_0_yyzz_xyzz[k] * ab_x + g_z_0_yyzz_xxyzz[k];

                g_z_0_xyyzz_xzzz[k] = -g_z_0_yyzz_xzzz[k] * ab_x + g_z_0_yyzz_xxzzz[k];

                g_z_0_xyyzz_yyyy[k] = -g_z_0_yyzz_yyyy[k] * ab_x + g_z_0_yyzz_xyyyy[k];

                g_z_0_xyyzz_yyyz[k] = -g_z_0_yyzz_yyyz[k] * ab_x + g_z_0_yyzz_xyyyz[k];

                g_z_0_xyyzz_yyzz[k] = -g_z_0_yyzz_yyzz[k] * ab_x + g_z_0_yyzz_xyyzz[k];

                g_z_0_xyyzz_yzzz[k] = -g_z_0_yyzz_yzzz[k] * ab_x + g_z_0_yyzz_xyzzz[k];

                g_z_0_xyyzz_zzzz[k] = -g_z_0_yyzz_zzzz[k] * ab_x + g_z_0_yyzz_xzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_yyyy, g_z_0_xyzzz_yyyz, g_z_0_xyzzz_yyzz, g_z_0_xyzzz_yzzz, g_z_0_xyzzz_zzzz, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxx[k] = -g_z_0_yzzz_xxxx[k] * ab_x + g_z_0_yzzz_xxxxx[k];

                g_z_0_xyzzz_xxxy[k] = -g_z_0_yzzz_xxxy[k] * ab_x + g_z_0_yzzz_xxxxy[k];

                g_z_0_xyzzz_xxxz[k] = -g_z_0_yzzz_xxxz[k] * ab_x + g_z_0_yzzz_xxxxz[k];

                g_z_0_xyzzz_xxyy[k] = -g_z_0_yzzz_xxyy[k] * ab_x + g_z_0_yzzz_xxxyy[k];

                g_z_0_xyzzz_xxyz[k] = -g_z_0_yzzz_xxyz[k] * ab_x + g_z_0_yzzz_xxxyz[k];

                g_z_0_xyzzz_xxzz[k] = -g_z_0_yzzz_xxzz[k] * ab_x + g_z_0_yzzz_xxxzz[k];

                g_z_0_xyzzz_xyyy[k] = -g_z_0_yzzz_xyyy[k] * ab_x + g_z_0_yzzz_xxyyy[k];

                g_z_0_xyzzz_xyyz[k] = -g_z_0_yzzz_xyyz[k] * ab_x + g_z_0_yzzz_xxyyz[k];

                g_z_0_xyzzz_xyzz[k] = -g_z_0_yzzz_xyzz[k] * ab_x + g_z_0_yzzz_xxyzz[k];

                g_z_0_xyzzz_xzzz[k] = -g_z_0_yzzz_xzzz[k] * ab_x + g_z_0_yzzz_xxzzz[k];

                g_z_0_xyzzz_yyyy[k] = -g_z_0_yzzz_yyyy[k] * ab_x + g_z_0_yzzz_xyyyy[k];

                g_z_0_xyzzz_yyyz[k] = -g_z_0_yzzz_yyyz[k] * ab_x + g_z_0_yzzz_xyyyz[k];

                g_z_0_xyzzz_yyzz[k] = -g_z_0_yzzz_yyzz[k] * ab_x + g_z_0_yzzz_xyyzz[k];

                g_z_0_xyzzz_yzzz[k] = -g_z_0_yzzz_yzzz[k] * ab_x + g_z_0_yzzz_xyzzz[k];

                g_z_0_xyzzz_zzzz[k] = -g_z_0_yzzz_zzzz[k] * ab_x + g_z_0_yzzz_xzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_yyyy, g_z_0_xzzzz_yyyz, g_z_0_xzzzz_yyzz, g_z_0_xzzzz_yzzz, g_z_0_xzzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] * ab_x + g_z_0_zzzz_xxxxx[k];

                g_z_0_xzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] * ab_x + g_z_0_zzzz_xxxxy[k];

                g_z_0_xzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] * ab_x + g_z_0_zzzz_xxxxz[k];

                g_z_0_xzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] * ab_x + g_z_0_zzzz_xxxyy[k];

                g_z_0_xzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] * ab_x + g_z_0_zzzz_xxxyz[k];

                g_z_0_xzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] * ab_x + g_z_0_zzzz_xxxzz[k];

                g_z_0_xzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] * ab_x + g_z_0_zzzz_xxyyy[k];

                g_z_0_xzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] * ab_x + g_z_0_zzzz_xxyyz[k];

                g_z_0_xzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] * ab_x + g_z_0_zzzz_xxyzz[k];

                g_z_0_xzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] * ab_x + g_z_0_zzzz_xxzzz[k];

                g_z_0_xzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] * ab_x + g_z_0_zzzz_xyyyy[k];

                g_z_0_xzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] * ab_x + g_z_0_zzzz_xyyyz[k];

                g_z_0_xzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] * ab_x + g_z_0_zzzz_xyyzz[k];

                g_z_0_xzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] * ab_x + g_z_0_zzzz_xyzzz[k];

                g_z_0_xzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] * ab_x + g_z_0_zzzz_xzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_zzzz, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxx[k] = -g_z_0_yyyy_xxxx[k] * ab_y + g_z_0_yyyy_xxxxy[k];

                g_z_0_yyyyy_xxxy[k] = -g_z_0_yyyy_xxxy[k] * ab_y + g_z_0_yyyy_xxxyy[k];

                g_z_0_yyyyy_xxxz[k] = -g_z_0_yyyy_xxxz[k] * ab_y + g_z_0_yyyy_xxxyz[k];

                g_z_0_yyyyy_xxyy[k] = -g_z_0_yyyy_xxyy[k] * ab_y + g_z_0_yyyy_xxyyy[k];

                g_z_0_yyyyy_xxyz[k] = -g_z_0_yyyy_xxyz[k] * ab_y + g_z_0_yyyy_xxyyz[k];

                g_z_0_yyyyy_xxzz[k] = -g_z_0_yyyy_xxzz[k] * ab_y + g_z_0_yyyy_xxyzz[k];

                g_z_0_yyyyy_xyyy[k] = -g_z_0_yyyy_xyyy[k] * ab_y + g_z_0_yyyy_xyyyy[k];

                g_z_0_yyyyy_xyyz[k] = -g_z_0_yyyy_xyyz[k] * ab_y + g_z_0_yyyy_xyyyz[k];

                g_z_0_yyyyy_xyzz[k] = -g_z_0_yyyy_xyzz[k] * ab_y + g_z_0_yyyy_xyyzz[k];

                g_z_0_yyyyy_xzzz[k] = -g_z_0_yyyy_xzzz[k] * ab_y + g_z_0_yyyy_xyzzz[k];

                g_z_0_yyyyy_yyyy[k] = -g_z_0_yyyy_yyyy[k] * ab_y + g_z_0_yyyy_yyyyy[k];

                g_z_0_yyyyy_yyyz[k] = -g_z_0_yyyy_yyyz[k] * ab_y + g_z_0_yyyy_yyyyz[k];

                g_z_0_yyyyy_yyzz[k] = -g_z_0_yyyy_yyzz[k] * ab_y + g_z_0_yyyy_yyyzz[k];

                g_z_0_yyyyy_yzzz[k] = -g_z_0_yyyy_yzzz[k] * ab_y + g_z_0_yyyy_yyzzz[k];

                g_z_0_yyyyy_zzzz[k] = -g_z_0_yyyy_zzzz[k] * ab_y + g_z_0_yyyy_yzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzzz, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxx[k] = -g_z_0_yyyz_xxxx[k] * ab_y + g_z_0_yyyz_xxxxy[k];

                g_z_0_yyyyz_xxxy[k] = -g_z_0_yyyz_xxxy[k] * ab_y + g_z_0_yyyz_xxxyy[k];

                g_z_0_yyyyz_xxxz[k] = -g_z_0_yyyz_xxxz[k] * ab_y + g_z_0_yyyz_xxxyz[k];

                g_z_0_yyyyz_xxyy[k] = -g_z_0_yyyz_xxyy[k] * ab_y + g_z_0_yyyz_xxyyy[k];

                g_z_0_yyyyz_xxyz[k] = -g_z_0_yyyz_xxyz[k] * ab_y + g_z_0_yyyz_xxyyz[k];

                g_z_0_yyyyz_xxzz[k] = -g_z_0_yyyz_xxzz[k] * ab_y + g_z_0_yyyz_xxyzz[k];

                g_z_0_yyyyz_xyyy[k] = -g_z_0_yyyz_xyyy[k] * ab_y + g_z_0_yyyz_xyyyy[k];

                g_z_0_yyyyz_xyyz[k] = -g_z_0_yyyz_xyyz[k] * ab_y + g_z_0_yyyz_xyyyz[k];

                g_z_0_yyyyz_xyzz[k] = -g_z_0_yyyz_xyzz[k] * ab_y + g_z_0_yyyz_xyyzz[k];

                g_z_0_yyyyz_xzzz[k] = -g_z_0_yyyz_xzzz[k] * ab_y + g_z_0_yyyz_xyzzz[k];

                g_z_0_yyyyz_yyyy[k] = -g_z_0_yyyz_yyyy[k] * ab_y + g_z_0_yyyz_yyyyy[k];

                g_z_0_yyyyz_yyyz[k] = -g_z_0_yyyz_yyyz[k] * ab_y + g_z_0_yyyz_yyyyz[k];

                g_z_0_yyyyz_yyzz[k] = -g_z_0_yyyz_yyzz[k] * ab_y + g_z_0_yyyz_yyyzz[k];

                g_z_0_yyyyz_yzzz[k] = -g_z_0_yyyz_yzzz[k] * ab_y + g_z_0_yyyz_yyzzz[k];

                g_z_0_yyyyz_zzzz[k] = -g_z_0_yyyz_zzzz[k] * ab_y + g_z_0_yyyz_yzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzzz, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxx[k] = -g_z_0_yyzz_xxxx[k] * ab_y + g_z_0_yyzz_xxxxy[k];

                g_z_0_yyyzz_xxxy[k] = -g_z_0_yyzz_xxxy[k] * ab_y + g_z_0_yyzz_xxxyy[k];

                g_z_0_yyyzz_xxxz[k] = -g_z_0_yyzz_xxxz[k] * ab_y + g_z_0_yyzz_xxxyz[k];

                g_z_0_yyyzz_xxyy[k] = -g_z_0_yyzz_xxyy[k] * ab_y + g_z_0_yyzz_xxyyy[k];

                g_z_0_yyyzz_xxyz[k] = -g_z_0_yyzz_xxyz[k] * ab_y + g_z_0_yyzz_xxyyz[k];

                g_z_0_yyyzz_xxzz[k] = -g_z_0_yyzz_xxzz[k] * ab_y + g_z_0_yyzz_xxyzz[k];

                g_z_0_yyyzz_xyyy[k] = -g_z_0_yyzz_xyyy[k] * ab_y + g_z_0_yyzz_xyyyy[k];

                g_z_0_yyyzz_xyyz[k] = -g_z_0_yyzz_xyyz[k] * ab_y + g_z_0_yyzz_xyyyz[k];

                g_z_0_yyyzz_xyzz[k] = -g_z_0_yyzz_xyzz[k] * ab_y + g_z_0_yyzz_xyyzz[k];

                g_z_0_yyyzz_xzzz[k] = -g_z_0_yyzz_xzzz[k] * ab_y + g_z_0_yyzz_xyzzz[k];

                g_z_0_yyyzz_yyyy[k] = -g_z_0_yyzz_yyyy[k] * ab_y + g_z_0_yyzz_yyyyy[k];

                g_z_0_yyyzz_yyyz[k] = -g_z_0_yyzz_yyyz[k] * ab_y + g_z_0_yyzz_yyyyz[k];

                g_z_0_yyyzz_yyzz[k] = -g_z_0_yyzz_yyzz[k] * ab_y + g_z_0_yyzz_yyyzz[k];

                g_z_0_yyyzz_yzzz[k] = -g_z_0_yyzz_yzzz[k] * ab_y + g_z_0_yyzz_yyzzz[k];

                g_z_0_yyyzz_zzzz[k] = -g_z_0_yyzz_zzzz[k] * ab_y + g_z_0_yyzz_yzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzzz, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxx[k] = -g_z_0_yzzz_xxxx[k] * ab_y + g_z_0_yzzz_xxxxy[k];

                g_z_0_yyzzz_xxxy[k] = -g_z_0_yzzz_xxxy[k] * ab_y + g_z_0_yzzz_xxxyy[k];

                g_z_0_yyzzz_xxxz[k] = -g_z_0_yzzz_xxxz[k] * ab_y + g_z_0_yzzz_xxxyz[k];

                g_z_0_yyzzz_xxyy[k] = -g_z_0_yzzz_xxyy[k] * ab_y + g_z_0_yzzz_xxyyy[k];

                g_z_0_yyzzz_xxyz[k] = -g_z_0_yzzz_xxyz[k] * ab_y + g_z_0_yzzz_xxyyz[k];

                g_z_0_yyzzz_xxzz[k] = -g_z_0_yzzz_xxzz[k] * ab_y + g_z_0_yzzz_xxyzz[k];

                g_z_0_yyzzz_xyyy[k] = -g_z_0_yzzz_xyyy[k] * ab_y + g_z_0_yzzz_xyyyy[k];

                g_z_0_yyzzz_xyyz[k] = -g_z_0_yzzz_xyyz[k] * ab_y + g_z_0_yzzz_xyyyz[k];

                g_z_0_yyzzz_xyzz[k] = -g_z_0_yzzz_xyzz[k] * ab_y + g_z_0_yzzz_xyyzz[k];

                g_z_0_yyzzz_xzzz[k] = -g_z_0_yzzz_xzzz[k] * ab_y + g_z_0_yzzz_xyzzz[k];

                g_z_0_yyzzz_yyyy[k] = -g_z_0_yzzz_yyyy[k] * ab_y + g_z_0_yzzz_yyyyy[k];

                g_z_0_yyzzz_yyyz[k] = -g_z_0_yzzz_yyyz[k] * ab_y + g_z_0_yzzz_yyyyz[k];

                g_z_0_yyzzz_yyzz[k] = -g_z_0_yzzz_yyzz[k] * ab_y + g_z_0_yzzz_yyyzz[k];

                g_z_0_yyzzz_yzzz[k] = -g_z_0_yzzz_yzzz[k] * ab_y + g_z_0_yzzz_yyzzz[k];

                g_z_0_yyzzz_zzzz[k] = -g_z_0_yzzz_zzzz[k] * ab_y + g_z_0_yzzz_yzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] * ab_y + g_z_0_zzzz_xxxxy[k];

                g_z_0_yzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] * ab_y + g_z_0_zzzz_xxxyy[k];

                g_z_0_yzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] * ab_y + g_z_0_zzzz_xxxyz[k];

                g_z_0_yzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] * ab_y + g_z_0_zzzz_xxyyy[k];

                g_z_0_yzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] * ab_y + g_z_0_zzzz_xxyyz[k];

                g_z_0_yzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] * ab_y + g_z_0_zzzz_xxyzz[k];

                g_z_0_yzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] * ab_y + g_z_0_zzzz_xyyyy[k];

                g_z_0_yzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] * ab_y + g_z_0_zzzz_xyyyz[k];

                g_z_0_yzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] * ab_y + g_z_0_zzzz_xyyzz[k];

                g_z_0_yzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] * ab_y + g_z_0_zzzz_xyzzz[k];

                g_z_0_yzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] * ab_y + g_z_0_zzzz_yyyyy[k];

                g_z_0_yzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] * ab_y + g_z_0_zzzz_yyyyz[k];

                g_z_0_yzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] * ab_y + g_z_0_zzzz_yyyzz[k];

                g_z_0_yzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] * ab_y + g_z_0_zzzz_yyzzz[k];

                g_z_0_yzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] * ab_y + g_z_0_zzzz_yzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzz, g_z_0_zzzz_zzzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzzz, g_zzzz_xxxx, g_zzzz_xxxy, g_zzzz_xxxz, g_zzzz_xxyy, g_zzzz_xxyz, g_zzzz_xxzz, g_zzzz_xyyy, g_zzzz_xyyz, g_zzzz_xyzz, g_zzzz_xzzz, g_zzzz_yyyy, g_zzzz_yyyz, g_zzzz_yyzz, g_zzzz_yzzz, g_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxx[k] = -g_zzzz_xxxx[k] - g_z_0_zzzz_xxxx[k] * ab_z + g_z_0_zzzz_xxxxz[k];

                g_z_0_zzzzz_xxxy[k] = -g_zzzz_xxxy[k] - g_z_0_zzzz_xxxy[k] * ab_z + g_z_0_zzzz_xxxyz[k];

                g_z_0_zzzzz_xxxz[k] = -g_zzzz_xxxz[k] - g_z_0_zzzz_xxxz[k] * ab_z + g_z_0_zzzz_xxxzz[k];

                g_z_0_zzzzz_xxyy[k] = -g_zzzz_xxyy[k] - g_z_0_zzzz_xxyy[k] * ab_z + g_z_0_zzzz_xxyyz[k];

                g_z_0_zzzzz_xxyz[k] = -g_zzzz_xxyz[k] - g_z_0_zzzz_xxyz[k] * ab_z + g_z_0_zzzz_xxyzz[k];

                g_z_0_zzzzz_xxzz[k] = -g_zzzz_xxzz[k] - g_z_0_zzzz_xxzz[k] * ab_z + g_z_0_zzzz_xxzzz[k];

                g_z_0_zzzzz_xyyy[k] = -g_zzzz_xyyy[k] - g_z_0_zzzz_xyyy[k] * ab_z + g_z_0_zzzz_xyyyz[k];

                g_z_0_zzzzz_xyyz[k] = -g_zzzz_xyyz[k] - g_z_0_zzzz_xyyz[k] * ab_z + g_z_0_zzzz_xyyzz[k];

                g_z_0_zzzzz_xyzz[k] = -g_zzzz_xyzz[k] - g_z_0_zzzz_xyzz[k] * ab_z + g_z_0_zzzz_xyzzz[k];

                g_z_0_zzzzz_xzzz[k] = -g_zzzz_xzzz[k] - g_z_0_zzzz_xzzz[k] * ab_z + g_z_0_zzzz_xzzzz[k];

                g_z_0_zzzzz_yyyy[k] = -g_zzzz_yyyy[k] - g_z_0_zzzz_yyyy[k] * ab_z + g_z_0_zzzz_yyyyz[k];

                g_z_0_zzzzz_yyyz[k] = -g_zzzz_yyyz[k] - g_z_0_zzzz_yyyz[k] * ab_z + g_z_0_zzzz_yyyzz[k];

                g_z_0_zzzzz_yyzz[k] = -g_zzzz_yyzz[k] - g_z_0_zzzz_yyzz[k] * ab_z + g_z_0_zzzz_yyzzz[k];

                g_z_0_zzzzz_yzzz[k] = -g_zzzz_yzzz[k] - g_z_0_zzzz_yzzz[k] * ab_z + g_z_0_zzzz_yzzzz[k];

                g_z_0_zzzzz_zzzz[k] = -g_zzzz_zzzz[k] - g_z_0_zzzz_zzzz[k] * ab_z + g_z_0_zzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

