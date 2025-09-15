#include "ElectronRepulsionGeom2000ContrRecHGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_hgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_hgxx,
                                            const size_t idx_geom_10_ggxx,
                                            const size_t idx_geom_20_ggxx,
                                            const size_t idx_geom_20_ghxx,
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

            /// Set up components of auxilary buffer : GGSS

            const auto gg_geom_20_off = idx_geom_20_ggxx + i * dcomps + j;

            auto g_xx_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 149 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 164 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xx_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xx_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xx_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 179 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xx_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xx_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xx_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 194 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xx_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xx_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xx_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 209 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xx_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xx_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xx_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 269 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 299 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 300 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 301 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 302 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 303 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 304 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 305 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 306 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 307 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 308 * ccomps * dcomps);

            auto g_xy_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 309 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 310 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 311 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 312 * ccomps * dcomps);

            auto g_xy_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 313 * ccomps * dcomps);

            auto g_xy_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 314 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 317 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 323 * ccomps * dcomps);

            auto g_xy_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 329 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 335 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 336 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 337 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 338 * ccomps * dcomps);

            auto g_xy_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 339 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 340 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 341 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 342 * ccomps * dcomps);

            auto g_xy_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 343 * ccomps * dcomps);

            auto g_xy_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 344 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 345 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 346 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 347 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 348 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 349 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 350 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 351 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 352 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 353 * ccomps * dcomps);

            auto g_xy_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 354 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 355 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 356 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 357 * ccomps * dcomps);

            auto g_xy_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 358 * ccomps * dcomps);

            auto g_xy_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 359 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 360 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 361 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 362 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 363 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 364 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 365 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 366 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 367 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 368 * ccomps * dcomps);

            auto g_xy_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 369 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 370 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 371 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 372 * ccomps * dcomps);

            auto g_xy_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 373 * ccomps * dcomps);

            auto g_xy_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 374 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 375 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 376 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 377 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 378 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 379 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 380 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 381 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 382 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 383 * ccomps * dcomps);

            auto g_xy_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 384 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 385 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 386 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 387 * ccomps * dcomps);

            auto g_xy_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 388 * ccomps * dcomps);

            auto g_xy_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 389 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 390 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 391 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 392 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 393 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 394 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 395 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 396 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 397 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 398 * ccomps * dcomps);

            auto g_xy_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 399 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 400 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 401 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 402 * ccomps * dcomps);

            auto g_xy_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 403 * ccomps * dcomps);

            auto g_xy_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 404 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 405 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 406 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 407 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 408 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 409 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 410 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 411 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 412 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 413 * ccomps * dcomps);

            auto g_xy_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 414 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 415 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 416 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 417 * ccomps * dcomps);

            auto g_xy_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 418 * ccomps * dcomps);

            auto g_xy_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 419 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 420 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 421 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 422 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 423 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 424 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 425 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 426 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 427 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 428 * ccomps * dcomps);

            auto g_xy_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 429 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 430 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 431 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 432 * ccomps * dcomps);

            auto g_xy_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 433 * ccomps * dcomps);

            auto g_xy_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 434 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 435 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 436 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 437 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 438 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 439 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 440 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 441 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 442 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 443 * ccomps * dcomps);

            auto g_xy_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 444 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 445 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 446 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 447 * ccomps * dcomps);

            auto g_xy_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 448 * ccomps * dcomps);

            auto g_xy_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 449 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 450 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 451 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 452 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 453 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 454 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 455 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 456 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 457 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 458 * ccomps * dcomps);

            auto g_xz_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 459 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 460 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 461 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 462 * ccomps * dcomps);

            auto g_xz_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 463 * ccomps * dcomps);

            auto g_xz_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 464 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 465 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 466 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 467 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 468 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 469 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 470 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 471 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 472 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 473 * ccomps * dcomps);

            auto g_xz_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 474 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 475 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 476 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 477 * ccomps * dcomps);

            auto g_xz_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 478 * ccomps * dcomps);

            auto g_xz_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 479 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 480 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 481 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 482 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 483 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 484 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 485 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 486 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 487 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 488 * ccomps * dcomps);

            auto g_xz_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 489 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 490 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 491 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 492 * ccomps * dcomps);

            auto g_xz_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 493 * ccomps * dcomps);

            auto g_xz_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 494 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 495 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 496 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 497 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 498 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 499 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 500 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 501 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 502 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 503 * ccomps * dcomps);

            auto g_xz_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 504 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 505 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 506 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 507 * ccomps * dcomps);

            auto g_xz_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 508 * ccomps * dcomps);

            auto g_xz_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 509 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 510 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 511 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 512 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 513 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 514 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 515 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 516 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 517 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 518 * ccomps * dcomps);

            auto g_xz_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 519 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 520 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 521 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 522 * ccomps * dcomps);

            auto g_xz_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 523 * ccomps * dcomps);

            auto g_xz_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 524 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 525 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 526 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 527 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 528 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 529 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 530 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 531 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 532 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 533 * ccomps * dcomps);

            auto g_xz_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 534 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 535 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 536 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 537 * ccomps * dcomps);

            auto g_xz_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 538 * ccomps * dcomps);

            auto g_xz_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 539 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 540 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 541 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 542 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 543 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 544 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 545 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 546 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 547 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 548 * ccomps * dcomps);

            auto g_xz_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 549 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 550 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 551 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 552 * ccomps * dcomps);

            auto g_xz_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 553 * ccomps * dcomps);

            auto g_xz_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 554 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 555 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 556 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 557 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 558 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 559 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 560 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 561 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 562 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 563 * ccomps * dcomps);

            auto g_xz_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 564 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 565 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 566 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 567 * ccomps * dcomps);

            auto g_xz_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 568 * ccomps * dcomps);

            auto g_xz_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 569 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 570 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 571 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 572 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 573 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 574 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 575 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 576 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 577 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 578 * ccomps * dcomps);

            auto g_xz_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 579 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 580 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 581 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 582 * ccomps * dcomps);

            auto g_xz_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 583 * ccomps * dcomps);

            auto g_xz_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 584 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 585 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 586 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 587 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 588 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 589 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 590 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 591 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 592 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 593 * ccomps * dcomps);

            auto g_xz_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 594 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 595 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 596 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 597 * ccomps * dcomps);

            auto g_xz_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 598 * ccomps * dcomps);

            auto g_xz_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 599 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 600 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 601 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 602 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 603 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 604 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 605 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 606 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 607 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 608 * ccomps * dcomps);

            auto g_xz_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 609 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 610 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 611 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 612 * ccomps * dcomps);

            auto g_xz_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 613 * ccomps * dcomps);

            auto g_xz_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 614 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 615 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 616 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 617 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 618 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 619 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 620 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 621 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 622 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 623 * ccomps * dcomps);

            auto g_xz_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 624 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 625 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 626 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 627 * ccomps * dcomps);

            auto g_xz_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 628 * ccomps * dcomps);

            auto g_xz_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 629 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 630 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 631 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 632 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 633 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 634 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 635 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 636 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 637 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 638 * ccomps * dcomps);

            auto g_xz_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 639 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 640 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 641 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 642 * ccomps * dcomps);

            auto g_xz_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 643 * ccomps * dcomps);

            auto g_xz_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 644 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 645 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 646 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 647 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 648 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 649 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 650 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 651 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 652 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 653 * ccomps * dcomps);

            auto g_xz_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 654 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 655 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 656 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 657 * ccomps * dcomps);

            auto g_xz_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 658 * ccomps * dcomps);

            auto g_xz_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 659 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 660 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 661 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 662 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 663 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 664 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 665 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 666 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 667 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 668 * ccomps * dcomps);

            auto g_xz_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 669 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 670 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 671 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 672 * ccomps * dcomps);

            auto g_xz_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 673 * ccomps * dcomps);

            auto g_xz_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 674 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 675 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 676 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 677 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 678 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 679 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 680 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 681 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 682 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 683 * ccomps * dcomps);

            auto g_yy_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 684 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 685 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 686 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 687 * ccomps * dcomps);

            auto g_yy_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 688 * ccomps * dcomps);

            auto g_yy_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 689 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 690 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 691 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 692 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 693 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 694 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 695 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 696 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 697 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 698 * ccomps * dcomps);

            auto g_yy_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 699 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 700 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 701 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 702 * ccomps * dcomps);

            auto g_yy_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 703 * ccomps * dcomps);

            auto g_yy_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 704 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 705 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 706 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 707 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 708 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 709 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 710 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 711 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 712 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 713 * ccomps * dcomps);

            auto g_yy_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 714 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 715 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 716 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 717 * ccomps * dcomps);

            auto g_yy_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 718 * ccomps * dcomps);

            auto g_yy_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 719 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 720 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 721 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 722 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 723 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 724 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 725 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 726 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 727 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 728 * ccomps * dcomps);

            auto g_yy_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 729 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 730 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 731 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 732 * ccomps * dcomps);

            auto g_yy_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 733 * ccomps * dcomps);

            auto g_yy_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 734 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 735 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 736 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 737 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 738 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 739 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 740 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 741 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 742 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 743 * ccomps * dcomps);

            auto g_yy_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 744 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 745 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 746 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 747 * ccomps * dcomps);

            auto g_yy_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 748 * ccomps * dcomps);

            auto g_yy_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 749 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 750 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 751 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 752 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 753 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 754 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 755 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 756 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 757 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 758 * ccomps * dcomps);

            auto g_yy_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 759 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 760 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 761 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 762 * ccomps * dcomps);

            auto g_yy_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 763 * ccomps * dcomps);

            auto g_yy_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 764 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 765 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 766 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 767 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 768 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 769 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 770 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 771 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 772 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 773 * ccomps * dcomps);

            auto g_yy_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 774 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 775 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 776 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 777 * ccomps * dcomps);

            auto g_yy_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 778 * ccomps * dcomps);

            auto g_yy_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 779 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 780 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 781 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 782 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 783 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 784 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 785 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 786 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 787 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 788 * ccomps * dcomps);

            auto g_yy_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 789 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 790 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 791 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 792 * ccomps * dcomps);

            auto g_yy_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 793 * ccomps * dcomps);

            auto g_yy_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 794 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 795 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 796 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 797 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 798 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 799 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 800 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 801 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 802 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 803 * ccomps * dcomps);

            auto g_yy_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 804 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 805 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 806 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 807 * ccomps * dcomps);

            auto g_yy_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 808 * ccomps * dcomps);

            auto g_yy_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 809 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 810 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 811 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 812 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 813 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 814 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 815 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 816 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 817 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 818 * ccomps * dcomps);

            auto g_yy_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 819 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 820 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 821 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 822 * ccomps * dcomps);

            auto g_yy_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 823 * ccomps * dcomps);

            auto g_yy_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 824 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 825 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 826 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 827 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 828 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 829 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 830 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 831 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 832 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 833 * ccomps * dcomps);

            auto g_yy_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 834 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 835 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 836 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 837 * ccomps * dcomps);

            auto g_yy_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 838 * ccomps * dcomps);

            auto g_yy_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 839 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 840 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 841 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 842 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 843 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 844 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 845 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 846 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 847 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 848 * ccomps * dcomps);

            auto g_yy_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 849 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 850 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 851 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 852 * ccomps * dcomps);

            auto g_yy_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 853 * ccomps * dcomps);

            auto g_yy_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 854 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 855 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 856 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 857 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 858 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 859 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 860 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 861 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 862 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 863 * ccomps * dcomps);

            auto g_yy_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 864 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 865 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 866 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 867 * ccomps * dcomps);

            auto g_yy_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 868 * ccomps * dcomps);

            auto g_yy_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 869 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 870 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 871 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 872 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 873 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 874 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 875 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 876 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 877 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 878 * ccomps * dcomps);

            auto g_yy_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 879 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 880 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 881 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 882 * ccomps * dcomps);

            auto g_yy_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 883 * ccomps * dcomps);

            auto g_yy_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 884 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 885 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 886 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 887 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 888 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 889 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 890 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 891 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 892 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 893 * ccomps * dcomps);

            auto g_yy_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 894 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 895 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 896 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 897 * ccomps * dcomps);

            auto g_yy_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 898 * ccomps * dcomps);

            auto g_yy_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 899 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 900 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 901 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 902 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 903 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 904 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 905 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 906 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 907 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 908 * ccomps * dcomps);

            auto g_yz_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 909 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 910 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 911 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 912 * ccomps * dcomps);

            auto g_yz_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 913 * ccomps * dcomps);

            auto g_yz_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 914 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 915 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 916 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 917 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 918 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 919 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 920 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 921 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 922 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 923 * ccomps * dcomps);

            auto g_yz_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 924 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 925 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 926 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 927 * ccomps * dcomps);

            auto g_yz_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 928 * ccomps * dcomps);

            auto g_yz_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 929 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 930 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 931 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 932 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 933 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 934 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 935 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 936 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 937 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 938 * ccomps * dcomps);

            auto g_yz_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 939 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 940 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 941 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 942 * ccomps * dcomps);

            auto g_yz_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 943 * ccomps * dcomps);

            auto g_yz_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 944 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 945 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 946 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 947 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 948 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 949 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 950 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 951 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 952 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 953 * ccomps * dcomps);

            auto g_yz_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 954 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 955 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 956 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 957 * ccomps * dcomps);

            auto g_yz_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 958 * ccomps * dcomps);

            auto g_yz_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 959 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 960 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 961 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 962 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 963 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 964 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 965 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 966 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 967 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 968 * ccomps * dcomps);

            auto g_yz_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 969 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 970 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 971 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 972 * ccomps * dcomps);

            auto g_yz_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 973 * ccomps * dcomps);

            auto g_yz_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 974 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 975 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 976 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 977 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 978 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 979 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 980 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 981 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 982 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 983 * ccomps * dcomps);

            auto g_yz_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 984 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 985 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 986 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 987 * ccomps * dcomps);

            auto g_yz_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 988 * ccomps * dcomps);

            auto g_yz_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 989 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 990 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 991 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 992 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 993 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 994 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 995 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 996 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 997 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 998 * ccomps * dcomps);

            auto g_yz_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 999 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 1000 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 1001 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 1002 * ccomps * dcomps);

            auto g_yz_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 1003 * ccomps * dcomps);

            auto g_yz_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 1004 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 1005 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 1006 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 1007 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 1008 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 1009 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 1010 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 1011 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 1012 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 1013 * ccomps * dcomps);

            auto g_yz_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 1014 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 1015 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 1016 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 1017 * ccomps * dcomps);

            auto g_yz_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 1018 * ccomps * dcomps);

            auto g_yz_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 1019 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 1020 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 1021 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 1022 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 1023 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 1024 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 1025 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 1026 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 1027 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 1028 * ccomps * dcomps);

            auto g_yz_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 1029 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 1030 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 1031 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 1032 * ccomps * dcomps);

            auto g_yz_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 1033 * ccomps * dcomps);

            auto g_yz_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 1034 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 1035 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 1036 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 1037 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 1038 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 1039 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 1040 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 1041 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 1042 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 1043 * ccomps * dcomps);

            auto g_yz_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 1044 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 1045 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 1046 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 1047 * ccomps * dcomps);

            auto g_yz_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 1048 * ccomps * dcomps);

            auto g_yz_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 1049 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 1050 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 1051 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 1052 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 1053 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 1054 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 1055 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 1056 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 1057 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 1058 * ccomps * dcomps);

            auto g_yz_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 1059 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 1060 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 1061 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 1062 * ccomps * dcomps);

            auto g_yz_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 1063 * ccomps * dcomps);

            auto g_yz_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 1064 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 1065 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 1066 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 1067 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 1068 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 1069 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 1070 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 1071 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 1072 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 1073 * ccomps * dcomps);

            auto g_yz_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 1074 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 1075 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 1076 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 1077 * ccomps * dcomps);

            auto g_yz_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 1078 * ccomps * dcomps);

            auto g_yz_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 1079 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 1080 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 1081 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 1082 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 1083 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 1084 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 1085 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 1086 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 1087 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 1088 * ccomps * dcomps);

            auto g_yz_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 1089 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 1090 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 1091 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 1092 * ccomps * dcomps);

            auto g_yz_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 1093 * ccomps * dcomps);

            auto g_yz_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 1094 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 1095 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 1096 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 1097 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 1098 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 1099 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 1100 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 1101 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 1102 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 1103 * ccomps * dcomps);

            auto g_yz_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 1104 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 1105 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 1106 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 1107 * ccomps * dcomps);

            auto g_yz_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 1108 * ccomps * dcomps);

            auto g_yz_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 1109 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 1110 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 1111 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 1112 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 1113 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 1114 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 1115 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 1116 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 1117 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 1118 * ccomps * dcomps);

            auto g_yz_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 1119 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 1120 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 1121 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 1122 * ccomps * dcomps);

            auto g_yz_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 1123 * ccomps * dcomps);

            auto g_yz_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 1124 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxx = cbuffer.data(gg_geom_20_off + 1125 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxy = cbuffer.data(gg_geom_20_off + 1126 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxz = cbuffer.data(gg_geom_20_off + 1127 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxyy = cbuffer.data(gg_geom_20_off + 1128 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxyz = cbuffer.data(gg_geom_20_off + 1129 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxzz = cbuffer.data(gg_geom_20_off + 1130 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyyy = cbuffer.data(gg_geom_20_off + 1131 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyyz = cbuffer.data(gg_geom_20_off + 1132 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyzz = cbuffer.data(gg_geom_20_off + 1133 * ccomps * dcomps);

            auto g_zz_0_xxxx_xzzz = cbuffer.data(gg_geom_20_off + 1134 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyyy = cbuffer.data(gg_geom_20_off + 1135 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyyz = cbuffer.data(gg_geom_20_off + 1136 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyzz = cbuffer.data(gg_geom_20_off + 1137 * ccomps * dcomps);

            auto g_zz_0_xxxx_yzzz = cbuffer.data(gg_geom_20_off + 1138 * ccomps * dcomps);

            auto g_zz_0_xxxx_zzzz = cbuffer.data(gg_geom_20_off + 1139 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxx = cbuffer.data(gg_geom_20_off + 1140 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxy = cbuffer.data(gg_geom_20_off + 1141 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxz = cbuffer.data(gg_geom_20_off + 1142 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxyy = cbuffer.data(gg_geom_20_off + 1143 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxyz = cbuffer.data(gg_geom_20_off + 1144 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxzz = cbuffer.data(gg_geom_20_off + 1145 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyyy = cbuffer.data(gg_geom_20_off + 1146 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyyz = cbuffer.data(gg_geom_20_off + 1147 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyzz = cbuffer.data(gg_geom_20_off + 1148 * ccomps * dcomps);

            auto g_zz_0_xxxy_xzzz = cbuffer.data(gg_geom_20_off + 1149 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyyy = cbuffer.data(gg_geom_20_off + 1150 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyyz = cbuffer.data(gg_geom_20_off + 1151 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyzz = cbuffer.data(gg_geom_20_off + 1152 * ccomps * dcomps);

            auto g_zz_0_xxxy_yzzz = cbuffer.data(gg_geom_20_off + 1153 * ccomps * dcomps);

            auto g_zz_0_xxxy_zzzz = cbuffer.data(gg_geom_20_off + 1154 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxx = cbuffer.data(gg_geom_20_off + 1155 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxy = cbuffer.data(gg_geom_20_off + 1156 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxz = cbuffer.data(gg_geom_20_off + 1157 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxyy = cbuffer.data(gg_geom_20_off + 1158 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxyz = cbuffer.data(gg_geom_20_off + 1159 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxzz = cbuffer.data(gg_geom_20_off + 1160 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyyy = cbuffer.data(gg_geom_20_off + 1161 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyyz = cbuffer.data(gg_geom_20_off + 1162 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyzz = cbuffer.data(gg_geom_20_off + 1163 * ccomps * dcomps);

            auto g_zz_0_xxxz_xzzz = cbuffer.data(gg_geom_20_off + 1164 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyyy = cbuffer.data(gg_geom_20_off + 1165 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyyz = cbuffer.data(gg_geom_20_off + 1166 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyzz = cbuffer.data(gg_geom_20_off + 1167 * ccomps * dcomps);

            auto g_zz_0_xxxz_yzzz = cbuffer.data(gg_geom_20_off + 1168 * ccomps * dcomps);

            auto g_zz_0_xxxz_zzzz = cbuffer.data(gg_geom_20_off + 1169 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxx = cbuffer.data(gg_geom_20_off + 1170 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxy = cbuffer.data(gg_geom_20_off + 1171 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxz = cbuffer.data(gg_geom_20_off + 1172 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxyy = cbuffer.data(gg_geom_20_off + 1173 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxyz = cbuffer.data(gg_geom_20_off + 1174 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxzz = cbuffer.data(gg_geom_20_off + 1175 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyyy = cbuffer.data(gg_geom_20_off + 1176 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyyz = cbuffer.data(gg_geom_20_off + 1177 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyzz = cbuffer.data(gg_geom_20_off + 1178 * ccomps * dcomps);

            auto g_zz_0_xxyy_xzzz = cbuffer.data(gg_geom_20_off + 1179 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyyy = cbuffer.data(gg_geom_20_off + 1180 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyyz = cbuffer.data(gg_geom_20_off + 1181 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyzz = cbuffer.data(gg_geom_20_off + 1182 * ccomps * dcomps);

            auto g_zz_0_xxyy_yzzz = cbuffer.data(gg_geom_20_off + 1183 * ccomps * dcomps);

            auto g_zz_0_xxyy_zzzz = cbuffer.data(gg_geom_20_off + 1184 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxx = cbuffer.data(gg_geom_20_off + 1185 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxy = cbuffer.data(gg_geom_20_off + 1186 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxz = cbuffer.data(gg_geom_20_off + 1187 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxyy = cbuffer.data(gg_geom_20_off + 1188 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxyz = cbuffer.data(gg_geom_20_off + 1189 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxzz = cbuffer.data(gg_geom_20_off + 1190 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyyy = cbuffer.data(gg_geom_20_off + 1191 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyyz = cbuffer.data(gg_geom_20_off + 1192 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyzz = cbuffer.data(gg_geom_20_off + 1193 * ccomps * dcomps);

            auto g_zz_0_xxyz_xzzz = cbuffer.data(gg_geom_20_off + 1194 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyyy = cbuffer.data(gg_geom_20_off + 1195 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyyz = cbuffer.data(gg_geom_20_off + 1196 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyzz = cbuffer.data(gg_geom_20_off + 1197 * ccomps * dcomps);

            auto g_zz_0_xxyz_yzzz = cbuffer.data(gg_geom_20_off + 1198 * ccomps * dcomps);

            auto g_zz_0_xxyz_zzzz = cbuffer.data(gg_geom_20_off + 1199 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxx = cbuffer.data(gg_geom_20_off + 1200 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxy = cbuffer.data(gg_geom_20_off + 1201 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxz = cbuffer.data(gg_geom_20_off + 1202 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxyy = cbuffer.data(gg_geom_20_off + 1203 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxyz = cbuffer.data(gg_geom_20_off + 1204 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxzz = cbuffer.data(gg_geom_20_off + 1205 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyyy = cbuffer.data(gg_geom_20_off + 1206 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyyz = cbuffer.data(gg_geom_20_off + 1207 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyzz = cbuffer.data(gg_geom_20_off + 1208 * ccomps * dcomps);

            auto g_zz_0_xxzz_xzzz = cbuffer.data(gg_geom_20_off + 1209 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyyy = cbuffer.data(gg_geom_20_off + 1210 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyyz = cbuffer.data(gg_geom_20_off + 1211 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyzz = cbuffer.data(gg_geom_20_off + 1212 * ccomps * dcomps);

            auto g_zz_0_xxzz_yzzz = cbuffer.data(gg_geom_20_off + 1213 * ccomps * dcomps);

            auto g_zz_0_xxzz_zzzz = cbuffer.data(gg_geom_20_off + 1214 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxx = cbuffer.data(gg_geom_20_off + 1215 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxy = cbuffer.data(gg_geom_20_off + 1216 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxz = cbuffer.data(gg_geom_20_off + 1217 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxyy = cbuffer.data(gg_geom_20_off + 1218 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxyz = cbuffer.data(gg_geom_20_off + 1219 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxzz = cbuffer.data(gg_geom_20_off + 1220 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyyy = cbuffer.data(gg_geom_20_off + 1221 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyyz = cbuffer.data(gg_geom_20_off + 1222 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyzz = cbuffer.data(gg_geom_20_off + 1223 * ccomps * dcomps);

            auto g_zz_0_xyyy_xzzz = cbuffer.data(gg_geom_20_off + 1224 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyyy = cbuffer.data(gg_geom_20_off + 1225 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyyz = cbuffer.data(gg_geom_20_off + 1226 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyzz = cbuffer.data(gg_geom_20_off + 1227 * ccomps * dcomps);

            auto g_zz_0_xyyy_yzzz = cbuffer.data(gg_geom_20_off + 1228 * ccomps * dcomps);

            auto g_zz_0_xyyy_zzzz = cbuffer.data(gg_geom_20_off + 1229 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxx = cbuffer.data(gg_geom_20_off + 1230 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxy = cbuffer.data(gg_geom_20_off + 1231 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxz = cbuffer.data(gg_geom_20_off + 1232 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxyy = cbuffer.data(gg_geom_20_off + 1233 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxyz = cbuffer.data(gg_geom_20_off + 1234 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxzz = cbuffer.data(gg_geom_20_off + 1235 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyyy = cbuffer.data(gg_geom_20_off + 1236 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyyz = cbuffer.data(gg_geom_20_off + 1237 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyzz = cbuffer.data(gg_geom_20_off + 1238 * ccomps * dcomps);

            auto g_zz_0_xyyz_xzzz = cbuffer.data(gg_geom_20_off + 1239 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyyy = cbuffer.data(gg_geom_20_off + 1240 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyyz = cbuffer.data(gg_geom_20_off + 1241 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyzz = cbuffer.data(gg_geom_20_off + 1242 * ccomps * dcomps);

            auto g_zz_0_xyyz_yzzz = cbuffer.data(gg_geom_20_off + 1243 * ccomps * dcomps);

            auto g_zz_0_xyyz_zzzz = cbuffer.data(gg_geom_20_off + 1244 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxx = cbuffer.data(gg_geom_20_off + 1245 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxy = cbuffer.data(gg_geom_20_off + 1246 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxz = cbuffer.data(gg_geom_20_off + 1247 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxyy = cbuffer.data(gg_geom_20_off + 1248 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxyz = cbuffer.data(gg_geom_20_off + 1249 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxzz = cbuffer.data(gg_geom_20_off + 1250 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyyy = cbuffer.data(gg_geom_20_off + 1251 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyyz = cbuffer.data(gg_geom_20_off + 1252 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyzz = cbuffer.data(gg_geom_20_off + 1253 * ccomps * dcomps);

            auto g_zz_0_xyzz_xzzz = cbuffer.data(gg_geom_20_off + 1254 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyyy = cbuffer.data(gg_geom_20_off + 1255 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyyz = cbuffer.data(gg_geom_20_off + 1256 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyzz = cbuffer.data(gg_geom_20_off + 1257 * ccomps * dcomps);

            auto g_zz_0_xyzz_yzzz = cbuffer.data(gg_geom_20_off + 1258 * ccomps * dcomps);

            auto g_zz_0_xyzz_zzzz = cbuffer.data(gg_geom_20_off + 1259 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxx = cbuffer.data(gg_geom_20_off + 1260 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxy = cbuffer.data(gg_geom_20_off + 1261 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxz = cbuffer.data(gg_geom_20_off + 1262 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxyy = cbuffer.data(gg_geom_20_off + 1263 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxyz = cbuffer.data(gg_geom_20_off + 1264 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxzz = cbuffer.data(gg_geom_20_off + 1265 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyyy = cbuffer.data(gg_geom_20_off + 1266 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyyz = cbuffer.data(gg_geom_20_off + 1267 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyzz = cbuffer.data(gg_geom_20_off + 1268 * ccomps * dcomps);

            auto g_zz_0_xzzz_xzzz = cbuffer.data(gg_geom_20_off + 1269 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyyy = cbuffer.data(gg_geom_20_off + 1270 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyyz = cbuffer.data(gg_geom_20_off + 1271 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyzz = cbuffer.data(gg_geom_20_off + 1272 * ccomps * dcomps);

            auto g_zz_0_xzzz_yzzz = cbuffer.data(gg_geom_20_off + 1273 * ccomps * dcomps);

            auto g_zz_0_xzzz_zzzz = cbuffer.data(gg_geom_20_off + 1274 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxx = cbuffer.data(gg_geom_20_off + 1275 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxy = cbuffer.data(gg_geom_20_off + 1276 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxz = cbuffer.data(gg_geom_20_off + 1277 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxyy = cbuffer.data(gg_geom_20_off + 1278 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxyz = cbuffer.data(gg_geom_20_off + 1279 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxzz = cbuffer.data(gg_geom_20_off + 1280 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyyy = cbuffer.data(gg_geom_20_off + 1281 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyyz = cbuffer.data(gg_geom_20_off + 1282 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyzz = cbuffer.data(gg_geom_20_off + 1283 * ccomps * dcomps);

            auto g_zz_0_yyyy_xzzz = cbuffer.data(gg_geom_20_off + 1284 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyyy = cbuffer.data(gg_geom_20_off + 1285 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyyz = cbuffer.data(gg_geom_20_off + 1286 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyzz = cbuffer.data(gg_geom_20_off + 1287 * ccomps * dcomps);

            auto g_zz_0_yyyy_yzzz = cbuffer.data(gg_geom_20_off + 1288 * ccomps * dcomps);

            auto g_zz_0_yyyy_zzzz = cbuffer.data(gg_geom_20_off + 1289 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxx = cbuffer.data(gg_geom_20_off + 1290 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxy = cbuffer.data(gg_geom_20_off + 1291 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxz = cbuffer.data(gg_geom_20_off + 1292 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxyy = cbuffer.data(gg_geom_20_off + 1293 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxyz = cbuffer.data(gg_geom_20_off + 1294 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxzz = cbuffer.data(gg_geom_20_off + 1295 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyyy = cbuffer.data(gg_geom_20_off + 1296 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyyz = cbuffer.data(gg_geom_20_off + 1297 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyzz = cbuffer.data(gg_geom_20_off + 1298 * ccomps * dcomps);

            auto g_zz_0_yyyz_xzzz = cbuffer.data(gg_geom_20_off + 1299 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyyy = cbuffer.data(gg_geom_20_off + 1300 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyyz = cbuffer.data(gg_geom_20_off + 1301 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyzz = cbuffer.data(gg_geom_20_off + 1302 * ccomps * dcomps);

            auto g_zz_0_yyyz_yzzz = cbuffer.data(gg_geom_20_off + 1303 * ccomps * dcomps);

            auto g_zz_0_yyyz_zzzz = cbuffer.data(gg_geom_20_off + 1304 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxx = cbuffer.data(gg_geom_20_off + 1305 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxy = cbuffer.data(gg_geom_20_off + 1306 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxz = cbuffer.data(gg_geom_20_off + 1307 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxyy = cbuffer.data(gg_geom_20_off + 1308 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxyz = cbuffer.data(gg_geom_20_off + 1309 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxzz = cbuffer.data(gg_geom_20_off + 1310 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyyy = cbuffer.data(gg_geom_20_off + 1311 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyyz = cbuffer.data(gg_geom_20_off + 1312 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyzz = cbuffer.data(gg_geom_20_off + 1313 * ccomps * dcomps);

            auto g_zz_0_yyzz_xzzz = cbuffer.data(gg_geom_20_off + 1314 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyyy = cbuffer.data(gg_geom_20_off + 1315 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyyz = cbuffer.data(gg_geom_20_off + 1316 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyzz = cbuffer.data(gg_geom_20_off + 1317 * ccomps * dcomps);

            auto g_zz_0_yyzz_yzzz = cbuffer.data(gg_geom_20_off + 1318 * ccomps * dcomps);

            auto g_zz_0_yyzz_zzzz = cbuffer.data(gg_geom_20_off + 1319 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxx = cbuffer.data(gg_geom_20_off + 1320 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxy = cbuffer.data(gg_geom_20_off + 1321 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxz = cbuffer.data(gg_geom_20_off + 1322 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxyy = cbuffer.data(gg_geom_20_off + 1323 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxyz = cbuffer.data(gg_geom_20_off + 1324 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxzz = cbuffer.data(gg_geom_20_off + 1325 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyyy = cbuffer.data(gg_geom_20_off + 1326 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyyz = cbuffer.data(gg_geom_20_off + 1327 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyzz = cbuffer.data(gg_geom_20_off + 1328 * ccomps * dcomps);

            auto g_zz_0_yzzz_xzzz = cbuffer.data(gg_geom_20_off + 1329 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyyy = cbuffer.data(gg_geom_20_off + 1330 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyyz = cbuffer.data(gg_geom_20_off + 1331 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyzz = cbuffer.data(gg_geom_20_off + 1332 * ccomps * dcomps);

            auto g_zz_0_yzzz_yzzz = cbuffer.data(gg_geom_20_off + 1333 * ccomps * dcomps);

            auto g_zz_0_yzzz_zzzz = cbuffer.data(gg_geom_20_off + 1334 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxx = cbuffer.data(gg_geom_20_off + 1335 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxy = cbuffer.data(gg_geom_20_off + 1336 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxz = cbuffer.data(gg_geom_20_off + 1337 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxyy = cbuffer.data(gg_geom_20_off + 1338 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxyz = cbuffer.data(gg_geom_20_off + 1339 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxzz = cbuffer.data(gg_geom_20_off + 1340 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyyy = cbuffer.data(gg_geom_20_off + 1341 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyyz = cbuffer.data(gg_geom_20_off + 1342 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyzz = cbuffer.data(gg_geom_20_off + 1343 * ccomps * dcomps);

            auto g_zz_0_zzzz_xzzz = cbuffer.data(gg_geom_20_off + 1344 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyyy = cbuffer.data(gg_geom_20_off + 1345 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyyz = cbuffer.data(gg_geom_20_off + 1346 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyzz = cbuffer.data(gg_geom_20_off + 1347 * ccomps * dcomps);

            auto g_zz_0_zzzz_yzzz = cbuffer.data(gg_geom_20_off + 1348 * ccomps * dcomps);

            auto g_zz_0_zzzz_zzzz = cbuffer.data(gg_geom_20_off + 1349 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GHSS

            const auto gh_geom_20_off = idx_geom_20_ghxx + i * dcomps + j;

            auto g_xx_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 149 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 155 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 161 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 164 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 167 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 168 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 169 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 170 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 171 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 172 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 173 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 174 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 175 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 176 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 177 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 178 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 179 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 180 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 181 * ccomps * dcomps);

            auto g_xx_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 182 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 183 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 184 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 185 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 186 * ccomps * dcomps);

            auto g_xx_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 187 * ccomps * dcomps);

            auto g_xx_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 188 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 189 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 190 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 191 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 192 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 193 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 194 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 195 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 196 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 197 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 198 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 199 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 200 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 201 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 202 * ccomps * dcomps);

            auto g_xx_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 203 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 204 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 205 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 206 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 207 * ccomps * dcomps);

            auto g_xx_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 208 * ccomps * dcomps);

            auto g_xx_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 209 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 210 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 211 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 212 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 213 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 214 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 215 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 216 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 217 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 218 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 219 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 220 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 221 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 222 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 223 * ccomps * dcomps);

            auto g_xx_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 224 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 225 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 226 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 227 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 228 * ccomps * dcomps);

            auto g_xx_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 229 * ccomps * dcomps);

            auto g_xx_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 230 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 231 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 232 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 233 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 234 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 235 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 236 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 237 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 238 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 239 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 240 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 241 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 242 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 243 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 244 * ccomps * dcomps);

            auto g_xx_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 245 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 246 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 247 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 248 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 249 * ccomps * dcomps);

            auto g_xx_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 250 * ccomps * dcomps);

            auto g_xx_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 251 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 252 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 253 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 254 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 255 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 256 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 257 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 258 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 259 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 260 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 261 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 262 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 263 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 264 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 265 * ccomps * dcomps);

            auto g_xx_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 266 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 267 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 268 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 269 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 270 * ccomps * dcomps);

            auto g_xx_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 271 * ccomps * dcomps);

            auto g_xx_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 272 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 273 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 274 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 275 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 276 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 277 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 278 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 279 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 280 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 281 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 282 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 283 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 284 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 285 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 286 * ccomps * dcomps);

            auto g_xx_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 287 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 288 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 289 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 290 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 291 * ccomps * dcomps);

            auto g_xx_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 292 * ccomps * dcomps);

            auto g_xx_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 293 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 294 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 295 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 296 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 297 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 298 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 299 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 300 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 301 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 302 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 303 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 304 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 305 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 306 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 307 * ccomps * dcomps);

            auto g_xx_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 308 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 309 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 310 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 311 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 312 * ccomps * dcomps);

            auto g_xx_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 313 * ccomps * dcomps);

            auto g_xx_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 314 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 317 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 323 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 329 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 335 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 336 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 337 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 338 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 339 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 340 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 341 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 342 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 343 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 344 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 345 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 346 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 347 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 348 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 349 * ccomps * dcomps);

            auto g_xy_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 350 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 351 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 352 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 353 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 354 * ccomps * dcomps);

            auto g_xy_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 355 * ccomps * dcomps);

            auto g_xy_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 356 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 357 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 358 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 359 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 360 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 361 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 362 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 363 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 364 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 365 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 366 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 367 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 368 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 369 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 370 * ccomps * dcomps);

            auto g_xy_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 371 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 372 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 373 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 374 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 375 * ccomps * dcomps);

            auto g_xy_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 376 * ccomps * dcomps);

            auto g_xy_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 377 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 378 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 379 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 380 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 381 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 382 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 383 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 384 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 385 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 386 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 387 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 388 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 389 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 390 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 391 * ccomps * dcomps);

            auto g_xy_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 392 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 393 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 394 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 395 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 396 * ccomps * dcomps);

            auto g_xy_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 397 * ccomps * dcomps);

            auto g_xy_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 398 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 399 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 400 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 401 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 402 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 403 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 404 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 405 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 406 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 407 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 408 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 409 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 410 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 411 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 412 * ccomps * dcomps);

            auto g_xy_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 413 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 414 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 415 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 416 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 417 * ccomps * dcomps);

            auto g_xy_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 418 * ccomps * dcomps);

            auto g_xy_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 419 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 420 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 421 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 422 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 423 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 424 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 425 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 426 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 427 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 428 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 429 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 430 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 431 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 432 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 433 * ccomps * dcomps);

            auto g_xy_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 434 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 435 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 436 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 437 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 438 * ccomps * dcomps);

            auto g_xy_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 439 * ccomps * dcomps);

            auto g_xy_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 440 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 441 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 442 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 443 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 444 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 445 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 446 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 447 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 448 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 449 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 450 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 451 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 452 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 453 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 454 * ccomps * dcomps);

            auto g_xy_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 455 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 456 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 457 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 458 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 459 * ccomps * dcomps);

            auto g_xy_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 460 * ccomps * dcomps);

            auto g_xy_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 461 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 462 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 463 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 464 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 465 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 466 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 467 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 468 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 469 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 470 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 471 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 472 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 473 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 474 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 475 * ccomps * dcomps);

            auto g_xy_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 476 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 477 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 478 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 479 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 480 * ccomps * dcomps);

            auto g_xy_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 481 * ccomps * dcomps);

            auto g_xy_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 482 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 483 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 484 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 485 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 486 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 487 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 488 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 489 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 490 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 491 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 492 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 493 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 494 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 495 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 496 * ccomps * dcomps);

            auto g_xy_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 497 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 498 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 499 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 500 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 501 * ccomps * dcomps);

            auto g_xy_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 502 * ccomps * dcomps);

            auto g_xy_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 503 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 504 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 505 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 506 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 507 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 508 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 509 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 510 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 511 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 512 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 513 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 514 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 515 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 516 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 517 * ccomps * dcomps);

            auto g_xy_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 518 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 519 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 520 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 521 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 522 * ccomps * dcomps);

            auto g_xy_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 523 * ccomps * dcomps);

            auto g_xy_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 524 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 525 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 526 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 527 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 528 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 529 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 530 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 531 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 532 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 533 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 534 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 535 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 536 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 537 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 538 * ccomps * dcomps);

            auto g_xy_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 539 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 540 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 541 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 542 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 543 * ccomps * dcomps);

            auto g_xy_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 544 * ccomps * dcomps);

            auto g_xy_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 545 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 546 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 547 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 548 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 549 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 550 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 551 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 552 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 553 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 554 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 555 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 556 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 557 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 558 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 559 * ccomps * dcomps);

            auto g_xy_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 560 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 561 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 562 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 563 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 564 * ccomps * dcomps);

            auto g_xy_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 565 * ccomps * dcomps);

            auto g_xy_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 566 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 567 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 568 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 569 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 570 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 571 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 572 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 573 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 574 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 575 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 576 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 577 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 578 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 579 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 580 * ccomps * dcomps);

            auto g_xy_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 581 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 582 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 583 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 584 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 585 * ccomps * dcomps);

            auto g_xy_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 586 * ccomps * dcomps);

            auto g_xy_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 587 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 588 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 589 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 590 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 591 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 592 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 593 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 594 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 595 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 596 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 597 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 598 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 599 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 600 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 601 * ccomps * dcomps);

            auto g_xy_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 602 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 603 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 604 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 605 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 606 * ccomps * dcomps);

            auto g_xy_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 607 * ccomps * dcomps);

            auto g_xy_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 608 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 609 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 610 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 611 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 612 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 613 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 614 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 615 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 616 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 617 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 618 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 619 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 620 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 621 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 622 * ccomps * dcomps);

            auto g_xy_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 623 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 624 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 625 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 626 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 627 * ccomps * dcomps);

            auto g_xy_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 628 * ccomps * dcomps);

            auto g_xy_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 629 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 630 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 631 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 632 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 633 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 634 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 635 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 636 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 637 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 638 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 639 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 640 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 641 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 642 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 643 * ccomps * dcomps);

            auto g_xz_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 644 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 645 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 646 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 647 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 648 * ccomps * dcomps);

            auto g_xz_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 649 * ccomps * dcomps);

            auto g_xz_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 650 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 651 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 652 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 653 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 654 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 655 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 656 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 657 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 658 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 659 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 660 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 661 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 662 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 663 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 664 * ccomps * dcomps);

            auto g_xz_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 665 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 666 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 667 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 668 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 669 * ccomps * dcomps);

            auto g_xz_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 670 * ccomps * dcomps);

            auto g_xz_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 671 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 672 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 673 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 674 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 675 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 676 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 677 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 678 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 679 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 680 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 681 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 682 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 683 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 684 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 685 * ccomps * dcomps);

            auto g_xz_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 686 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 687 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 688 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 689 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 690 * ccomps * dcomps);

            auto g_xz_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 691 * ccomps * dcomps);

            auto g_xz_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 692 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 693 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 694 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 695 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 696 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 697 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 698 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 699 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 700 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 701 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 702 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 703 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 704 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 705 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 706 * ccomps * dcomps);

            auto g_xz_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 707 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 708 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 709 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 710 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 711 * ccomps * dcomps);

            auto g_xz_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 712 * ccomps * dcomps);

            auto g_xz_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 713 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 714 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 715 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 716 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 717 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 718 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 719 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 720 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 721 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 722 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 723 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 724 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 725 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 726 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 727 * ccomps * dcomps);

            auto g_xz_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 728 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 729 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 730 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 731 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 732 * ccomps * dcomps);

            auto g_xz_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 733 * ccomps * dcomps);

            auto g_xz_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 734 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 735 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 736 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 737 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 738 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 739 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 740 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 741 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 742 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 743 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 744 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 745 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 746 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 747 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 748 * ccomps * dcomps);

            auto g_xz_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 749 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 750 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 751 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 752 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 753 * ccomps * dcomps);

            auto g_xz_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 754 * ccomps * dcomps);

            auto g_xz_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 755 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 756 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 757 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 758 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 759 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 760 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 761 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 762 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 763 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 764 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 765 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 766 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 767 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 768 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 769 * ccomps * dcomps);

            auto g_xz_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 770 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 771 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 772 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 773 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 774 * ccomps * dcomps);

            auto g_xz_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 775 * ccomps * dcomps);

            auto g_xz_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 776 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 777 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 778 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 779 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 780 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 781 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 782 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 783 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 784 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 785 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 786 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 787 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 788 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 789 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 790 * ccomps * dcomps);

            auto g_xz_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 791 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 792 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 793 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 794 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 795 * ccomps * dcomps);

            auto g_xz_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 796 * ccomps * dcomps);

            auto g_xz_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 797 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 798 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 799 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 800 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 801 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 802 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 803 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 804 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 805 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 806 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 807 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 808 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 809 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 810 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 811 * ccomps * dcomps);

            auto g_xz_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 812 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 813 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 814 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 815 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 816 * ccomps * dcomps);

            auto g_xz_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 817 * ccomps * dcomps);

            auto g_xz_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 818 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 819 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 820 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 821 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 822 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 823 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 824 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 825 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 826 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 827 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 828 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 829 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 830 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 831 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 832 * ccomps * dcomps);

            auto g_xz_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 833 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 834 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 835 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 836 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 837 * ccomps * dcomps);

            auto g_xz_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 838 * ccomps * dcomps);

            auto g_xz_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 839 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 840 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 841 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 842 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 843 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 844 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 845 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 846 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 847 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 848 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 849 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 850 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 851 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 852 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 853 * ccomps * dcomps);

            auto g_xz_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 854 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 855 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 856 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 857 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 858 * ccomps * dcomps);

            auto g_xz_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 859 * ccomps * dcomps);

            auto g_xz_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 860 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 861 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 862 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 863 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 864 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 865 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 866 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 867 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 868 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 869 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 870 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 871 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 872 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 873 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 874 * ccomps * dcomps);

            auto g_xz_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 875 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 876 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 877 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 878 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 879 * ccomps * dcomps);

            auto g_xz_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 880 * ccomps * dcomps);

            auto g_xz_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 881 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 882 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 883 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 884 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 885 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 886 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 887 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 888 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 889 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 890 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 891 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 892 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 893 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 894 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 895 * ccomps * dcomps);

            auto g_xz_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 896 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 897 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 898 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 899 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 900 * ccomps * dcomps);

            auto g_xz_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 901 * ccomps * dcomps);

            auto g_xz_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 902 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 903 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 904 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 905 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 906 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 907 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 908 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 909 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 910 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 911 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 912 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 913 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 914 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 915 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 916 * ccomps * dcomps);

            auto g_xz_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 917 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 918 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 919 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 920 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 921 * ccomps * dcomps);

            auto g_xz_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 922 * ccomps * dcomps);

            auto g_xz_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 923 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 924 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 925 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 926 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 927 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 928 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 929 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 930 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 931 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 932 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 933 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 934 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 935 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 936 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 937 * ccomps * dcomps);

            auto g_xz_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 938 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 939 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 940 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 941 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 942 * ccomps * dcomps);

            auto g_xz_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 943 * ccomps * dcomps);

            auto g_xz_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 944 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 945 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 946 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 947 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 948 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 949 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 950 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 951 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 952 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 953 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 954 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 955 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 956 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 957 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 958 * ccomps * dcomps);

            auto g_yy_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 959 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 960 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 961 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 962 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 963 * ccomps * dcomps);

            auto g_yy_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 964 * ccomps * dcomps);

            auto g_yy_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 965 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 966 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 967 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 968 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 969 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 970 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 971 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 972 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 973 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 974 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 975 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 976 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 977 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 978 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 979 * ccomps * dcomps);

            auto g_yy_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 980 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 981 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 982 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 983 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 984 * ccomps * dcomps);

            auto g_yy_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 985 * ccomps * dcomps);

            auto g_yy_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 986 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 987 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 988 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 989 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 990 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 991 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 992 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 993 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 994 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 995 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 996 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 997 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 998 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 999 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 1000 * ccomps * dcomps);

            auto g_yy_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 1001 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 1002 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 1003 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 1004 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 1005 * ccomps * dcomps);

            auto g_yy_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 1006 * ccomps * dcomps);

            auto g_yy_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 1007 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 1008 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 1009 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 1010 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 1011 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 1012 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 1013 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 1014 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 1015 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 1016 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 1017 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 1018 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 1019 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 1020 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 1021 * ccomps * dcomps);

            auto g_yy_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 1022 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 1023 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 1024 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 1025 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 1026 * ccomps * dcomps);

            auto g_yy_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 1027 * ccomps * dcomps);

            auto g_yy_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 1028 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 1029 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 1030 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 1031 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 1032 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 1033 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 1034 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 1035 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 1036 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 1037 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 1038 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 1039 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 1040 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 1041 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 1042 * ccomps * dcomps);

            auto g_yy_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 1043 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 1044 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 1045 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 1046 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 1047 * ccomps * dcomps);

            auto g_yy_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 1048 * ccomps * dcomps);

            auto g_yy_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 1049 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 1050 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 1051 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 1052 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 1053 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 1054 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 1055 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 1056 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 1057 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 1058 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 1059 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 1060 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 1061 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 1062 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 1063 * ccomps * dcomps);

            auto g_yy_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 1064 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 1065 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 1066 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 1067 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 1068 * ccomps * dcomps);

            auto g_yy_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 1069 * ccomps * dcomps);

            auto g_yy_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 1070 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1071 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1072 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1073 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1074 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1075 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1076 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1077 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1078 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1079 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1080 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1081 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1082 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1083 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1084 * ccomps * dcomps);

            auto g_yy_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1085 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1086 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1087 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1088 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1089 * ccomps * dcomps);

            auto g_yy_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1090 * ccomps * dcomps);

            auto g_yy_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1091 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1092 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1093 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1094 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1095 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1096 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1097 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1098 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1099 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1100 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1101 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1102 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1103 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1104 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1105 * ccomps * dcomps);

            auto g_yy_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1106 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1107 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1108 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1109 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1110 * ccomps * dcomps);

            auto g_yy_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1111 * ccomps * dcomps);

            auto g_yy_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1112 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1113 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1114 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1115 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1116 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1117 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1118 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1119 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1120 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1121 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1122 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1123 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1124 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1125 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1126 * ccomps * dcomps);

            auto g_yy_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1127 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1128 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1129 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1130 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1131 * ccomps * dcomps);

            auto g_yy_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1132 * ccomps * dcomps);

            auto g_yy_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1133 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1134 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1135 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1136 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1137 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1138 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1139 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1140 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1141 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1142 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1143 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1144 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1145 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1146 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1147 * ccomps * dcomps);

            auto g_yy_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1148 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1149 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1150 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1151 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1152 * ccomps * dcomps);

            auto g_yy_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1153 * ccomps * dcomps);

            auto g_yy_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1154 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1155 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1156 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1157 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1158 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1159 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1160 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1161 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1162 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1163 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1164 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1165 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1166 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1167 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1168 * ccomps * dcomps);

            auto g_yy_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1169 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1170 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1171 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1172 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1173 * ccomps * dcomps);

            auto g_yy_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1174 * ccomps * dcomps);

            auto g_yy_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1175 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1176 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1177 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1178 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1179 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1180 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1181 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1182 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1183 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1184 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1185 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1186 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1187 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1188 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1189 * ccomps * dcomps);

            auto g_yy_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1190 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1191 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1192 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1193 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1194 * ccomps * dcomps);

            auto g_yy_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1195 * ccomps * dcomps);

            auto g_yy_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1196 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1197 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1198 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1199 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1200 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1201 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1202 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1203 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1204 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1205 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1206 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1207 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1208 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1209 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1210 * ccomps * dcomps);

            auto g_yy_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1211 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1212 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1213 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1214 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1215 * ccomps * dcomps);

            auto g_yy_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1216 * ccomps * dcomps);

            auto g_yy_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1217 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1218 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1219 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1220 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1221 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1222 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1223 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1224 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1225 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1226 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1227 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1228 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1229 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1230 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1231 * ccomps * dcomps);

            auto g_yy_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1232 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1233 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1234 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1235 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1236 * ccomps * dcomps);

            auto g_yy_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1237 * ccomps * dcomps);

            auto g_yy_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1238 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1239 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1240 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1241 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1242 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1243 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1244 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1245 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1246 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1247 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1248 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1249 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1250 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1251 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1252 * ccomps * dcomps);

            auto g_yy_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1253 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1254 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1255 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1256 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1257 * ccomps * dcomps);

            auto g_yy_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1258 * ccomps * dcomps);

            auto g_yy_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1259 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 1260 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 1261 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 1262 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 1263 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 1264 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 1265 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 1266 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 1267 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 1268 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 1269 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 1270 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 1271 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 1272 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 1273 * ccomps * dcomps);

            auto g_yz_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 1274 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 1275 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 1276 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 1277 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 1278 * ccomps * dcomps);

            auto g_yz_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 1279 * ccomps * dcomps);

            auto g_yz_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 1280 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 1281 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 1282 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 1283 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 1284 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 1285 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 1286 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 1287 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 1288 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 1289 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 1290 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 1291 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 1292 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 1293 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 1294 * ccomps * dcomps);

            auto g_yz_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 1295 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 1296 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 1297 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 1298 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 1299 * ccomps * dcomps);

            auto g_yz_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 1300 * ccomps * dcomps);

            auto g_yz_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 1301 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 1302 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 1303 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 1304 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 1305 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 1306 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 1307 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 1308 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 1309 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 1310 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 1311 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 1312 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 1313 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 1314 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 1315 * ccomps * dcomps);

            auto g_yz_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 1316 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 1317 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 1318 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 1319 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 1320 * ccomps * dcomps);

            auto g_yz_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 1321 * ccomps * dcomps);

            auto g_yz_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 1322 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 1323 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 1324 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 1325 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 1326 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 1327 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 1328 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 1329 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 1330 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 1331 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 1332 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 1333 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 1334 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 1335 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 1336 * ccomps * dcomps);

            auto g_yz_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 1337 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 1338 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 1339 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 1340 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 1341 * ccomps * dcomps);

            auto g_yz_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 1342 * ccomps * dcomps);

            auto g_yz_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 1343 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 1344 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 1345 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 1346 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 1347 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 1348 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 1349 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 1350 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 1351 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 1352 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 1353 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 1354 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 1355 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 1356 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 1357 * ccomps * dcomps);

            auto g_yz_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 1358 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 1359 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 1360 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 1361 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 1362 * ccomps * dcomps);

            auto g_yz_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 1363 * ccomps * dcomps);

            auto g_yz_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 1364 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 1365 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 1366 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 1367 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 1368 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 1369 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 1370 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 1371 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 1372 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 1373 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 1374 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 1375 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 1376 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 1377 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 1378 * ccomps * dcomps);

            auto g_yz_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 1379 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 1380 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 1381 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 1382 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 1383 * ccomps * dcomps);

            auto g_yz_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 1384 * ccomps * dcomps);

            auto g_yz_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 1385 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1386 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1387 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1388 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1389 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1390 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1391 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1392 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1393 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1394 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1395 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1396 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1397 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1398 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1399 * ccomps * dcomps);

            auto g_yz_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1400 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1401 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1402 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1403 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1404 * ccomps * dcomps);

            auto g_yz_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1405 * ccomps * dcomps);

            auto g_yz_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1406 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1407 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1408 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1409 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1410 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1411 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1412 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1413 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1414 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1415 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1416 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1417 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1418 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1419 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1420 * ccomps * dcomps);

            auto g_yz_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1421 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1422 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1423 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1424 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1425 * ccomps * dcomps);

            auto g_yz_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1426 * ccomps * dcomps);

            auto g_yz_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1427 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1428 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1429 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1430 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1431 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1432 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1433 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1434 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1435 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1436 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1437 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1438 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1439 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1440 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1441 * ccomps * dcomps);

            auto g_yz_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1442 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1443 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1444 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1445 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1446 * ccomps * dcomps);

            auto g_yz_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1447 * ccomps * dcomps);

            auto g_yz_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1448 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1449 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1450 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1451 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1452 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1453 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1454 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1455 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1456 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1457 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1458 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1459 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1460 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1461 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1462 * ccomps * dcomps);

            auto g_yz_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1463 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1464 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1465 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1466 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1467 * ccomps * dcomps);

            auto g_yz_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1468 * ccomps * dcomps);

            auto g_yz_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1469 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1470 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1471 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1472 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1473 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1474 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1475 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1476 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1477 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1478 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1479 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1480 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1481 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1482 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1483 * ccomps * dcomps);

            auto g_yz_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1484 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1485 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1486 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1487 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1488 * ccomps * dcomps);

            auto g_yz_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1489 * ccomps * dcomps);

            auto g_yz_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1490 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1491 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1492 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1493 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1494 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1495 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1496 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1497 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1498 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1499 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1500 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1501 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1502 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1503 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1504 * ccomps * dcomps);

            auto g_yz_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1505 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1506 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1507 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1508 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1509 * ccomps * dcomps);

            auto g_yz_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1510 * ccomps * dcomps);

            auto g_yz_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1511 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1512 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1513 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1514 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1515 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1516 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1517 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1518 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1519 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1520 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1521 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1522 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1523 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1524 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1525 * ccomps * dcomps);

            auto g_yz_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1526 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1527 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1528 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1529 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1530 * ccomps * dcomps);

            auto g_yz_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1531 * ccomps * dcomps);

            auto g_yz_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1532 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1533 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1534 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1535 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1536 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1537 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1538 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1539 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1540 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1541 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1542 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1543 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1544 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1545 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1546 * ccomps * dcomps);

            auto g_yz_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1547 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1548 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1549 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1550 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1551 * ccomps * dcomps);

            auto g_yz_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1552 * ccomps * dcomps);

            auto g_yz_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1553 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1554 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1555 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1556 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1557 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1558 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1559 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1560 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1561 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1562 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1563 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1564 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1565 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1566 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1567 * ccomps * dcomps);

            auto g_yz_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1568 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1569 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1570 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1571 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1572 * ccomps * dcomps);

            auto g_yz_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1573 * ccomps * dcomps);

            auto g_yz_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1574 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxxx = cbuffer.data(gh_geom_20_off + 1575 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxxy = cbuffer.data(gh_geom_20_off + 1576 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxxz = cbuffer.data(gh_geom_20_off + 1577 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxyy = cbuffer.data(gh_geom_20_off + 1578 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxyz = cbuffer.data(gh_geom_20_off + 1579 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxxzz = cbuffer.data(gh_geom_20_off + 1580 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxyyy = cbuffer.data(gh_geom_20_off + 1581 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxyyz = cbuffer.data(gh_geom_20_off + 1582 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxyzz = cbuffer.data(gh_geom_20_off + 1583 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxzzz = cbuffer.data(gh_geom_20_off + 1584 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyyyy = cbuffer.data(gh_geom_20_off + 1585 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyyyz = cbuffer.data(gh_geom_20_off + 1586 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyyzz = cbuffer.data(gh_geom_20_off + 1587 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyzzz = cbuffer.data(gh_geom_20_off + 1588 * ccomps * dcomps);

            auto g_zz_0_xxxx_xzzzz = cbuffer.data(gh_geom_20_off + 1589 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyyyy = cbuffer.data(gh_geom_20_off + 1590 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyyyz = cbuffer.data(gh_geom_20_off + 1591 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyyzz = cbuffer.data(gh_geom_20_off + 1592 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyzzz = cbuffer.data(gh_geom_20_off + 1593 * ccomps * dcomps);

            auto g_zz_0_xxxx_yzzzz = cbuffer.data(gh_geom_20_off + 1594 * ccomps * dcomps);

            auto g_zz_0_xxxx_zzzzz = cbuffer.data(gh_geom_20_off + 1595 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxxx = cbuffer.data(gh_geom_20_off + 1596 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxxy = cbuffer.data(gh_geom_20_off + 1597 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxxz = cbuffer.data(gh_geom_20_off + 1598 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxyy = cbuffer.data(gh_geom_20_off + 1599 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxyz = cbuffer.data(gh_geom_20_off + 1600 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxxzz = cbuffer.data(gh_geom_20_off + 1601 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxyyy = cbuffer.data(gh_geom_20_off + 1602 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxyyz = cbuffer.data(gh_geom_20_off + 1603 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxyzz = cbuffer.data(gh_geom_20_off + 1604 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxzzz = cbuffer.data(gh_geom_20_off + 1605 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyyyy = cbuffer.data(gh_geom_20_off + 1606 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyyyz = cbuffer.data(gh_geom_20_off + 1607 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyyzz = cbuffer.data(gh_geom_20_off + 1608 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyzzz = cbuffer.data(gh_geom_20_off + 1609 * ccomps * dcomps);

            auto g_zz_0_xxxy_xzzzz = cbuffer.data(gh_geom_20_off + 1610 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyyyy = cbuffer.data(gh_geom_20_off + 1611 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyyyz = cbuffer.data(gh_geom_20_off + 1612 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyyzz = cbuffer.data(gh_geom_20_off + 1613 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyzzz = cbuffer.data(gh_geom_20_off + 1614 * ccomps * dcomps);

            auto g_zz_0_xxxy_yzzzz = cbuffer.data(gh_geom_20_off + 1615 * ccomps * dcomps);

            auto g_zz_0_xxxy_zzzzz = cbuffer.data(gh_geom_20_off + 1616 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxxx = cbuffer.data(gh_geom_20_off + 1617 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxxy = cbuffer.data(gh_geom_20_off + 1618 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxxz = cbuffer.data(gh_geom_20_off + 1619 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxyy = cbuffer.data(gh_geom_20_off + 1620 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxyz = cbuffer.data(gh_geom_20_off + 1621 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxxzz = cbuffer.data(gh_geom_20_off + 1622 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxyyy = cbuffer.data(gh_geom_20_off + 1623 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxyyz = cbuffer.data(gh_geom_20_off + 1624 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxyzz = cbuffer.data(gh_geom_20_off + 1625 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxzzz = cbuffer.data(gh_geom_20_off + 1626 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyyyy = cbuffer.data(gh_geom_20_off + 1627 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyyyz = cbuffer.data(gh_geom_20_off + 1628 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyyzz = cbuffer.data(gh_geom_20_off + 1629 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyzzz = cbuffer.data(gh_geom_20_off + 1630 * ccomps * dcomps);

            auto g_zz_0_xxxz_xzzzz = cbuffer.data(gh_geom_20_off + 1631 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyyyy = cbuffer.data(gh_geom_20_off + 1632 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyyyz = cbuffer.data(gh_geom_20_off + 1633 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyyzz = cbuffer.data(gh_geom_20_off + 1634 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyzzz = cbuffer.data(gh_geom_20_off + 1635 * ccomps * dcomps);

            auto g_zz_0_xxxz_yzzzz = cbuffer.data(gh_geom_20_off + 1636 * ccomps * dcomps);

            auto g_zz_0_xxxz_zzzzz = cbuffer.data(gh_geom_20_off + 1637 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxxx = cbuffer.data(gh_geom_20_off + 1638 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxxy = cbuffer.data(gh_geom_20_off + 1639 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxxz = cbuffer.data(gh_geom_20_off + 1640 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxyy = cbuffer.data(gh_geom_20_off + 1641 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxyz = cbuffer.data(gh_geom_20_off + 1642 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxxzz = cbuffer.data(gh_geom_20_off + 1643 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxyyy = cbuffer.data(gh_geom_20_off + 1644 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxyyz = cbuffer.data(gh_geom_20_off + 1645 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxyzz = cbuffer.data(gh_geom_20_off + 1646 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxzzz = cbuffer.data(gh_geom_20_off + 1647 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyyyy = cbuffer.data(gh_geom_20_off + 1648 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyyyz = cbuffer.data(gh_geom_20_off + 1649 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyyzz = cbuffer.data(gh_geom_20_off + 1650 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyzzz = cbuffer.data(gh_geom_20_off + 1651 * ccomps * dcomps);

            auto g_zz_0_xxyy_xzzzz = cbuffer.data(gh_geom_20_off + 1652 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyyyy = cbuffer.data(gh_geom_20_off + 1653 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyyyz = cbuffer.data(gh_geom_20_off + 1654 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyyzz = cbuffer.data(gh_geom_20_off + 1655 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyzzz = cbuffer.data(gh_geom_20_off + 1656 * ccomps * dcomps);

            auto g_zz_0_xxyy_yzzzz = cbuffer.data(gh_geom_20_off + 1657 * ccomps * dcomps);

            auto g_zz_0_xxyy_zzzzz = cbuffer.data(gh_geom_20_off + 1658 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxxx = cbuffer.data(gh_geom_20_off + 1659 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxxy = cbuffer.data(gh_geom_20_off + 1660 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxxz = cbuffer.data(gh_geom_20_off + 1661 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxyy = cbuffer.data(gh_geom_20_off + 1662 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxyz = cbuffer.data(gh_geom_20_off + 1663 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxxzz = cbuffer.data(gh_geom_20_off + 1664 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxyyy = cbuffer.data(gh_geom_20_off + 1665 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxyyz = cbuffer.data(gh_geom_20_off + 1666 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxyzz = cbuffer.data(gh_geom_20_off + 1667 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxzzz = cbuffer.data(gh_geom_20_off + 1668 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyyyy = cbuffer.data(gh_geom_20_off + 1669 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyyyz = cbuffer.data(gh_geom_20_off + 1670 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyyzz = cbuffer.data(gh_geom_20_off + 1671 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyzzz = cbuffer.data(gh_geom_20_off + 1672 * ccomps * dcomps);

            auto g_zz_0_xxyz_xzzzz = cbuffer.data(gh_geom_20_off + 1673 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyyyy = cbuffer.data(gh_geom_20_off + 1674 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyyyz = cbuffer.data(gh_geom_20_off + 1675 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyyzz = cbuffer.data(gh_geom_20_off + 1676 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyzzz = cbuffer.data(gh_geom_20_off + 1677 * ccomps * dcomps);

            auto g_zz_0_xxyz_yzzzz = cbuffer.data(gh_geom_20_off + 1678 * ccomps * dcomps);

            auto g_zz_0_xxyz_zzzzz = cbuffer.data(gh_geom_20_off + 1679 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxxx = cbuffer.data(gh_geom_20_off + 1680 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxxy = cbuffer.data(gh_geom_20_off + 1681 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxxz = cbuffer.data(gh_geom_20_off + 1682 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxyy = cbuffer.data(gh_geom_20_off + 1683 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxyz = cbuffer.data(gh_geom_20_off + 1684 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxxzz = cbuffer.data(gh_geom_20_off + 1685 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxyyy = cbuffer.data(gh_geom_20_off + 1686 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxyyz = cbuffer.data(gh_geom_20_off + 1687 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxyzz = cbuffer.data(gh_geom_20_off + 1688 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxzzz = cbuffer.data(gh_geom_20_off + 1689 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyyyy = cbuffer.data(gh_geom_20_off + 1690 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyyyz = cbuffer.data(gh_geom_20_off + 1691 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyyzz = cbuffer.data(gh_geom_20_off + 1692 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyzzz = cbuffer.data(gh_geom_20_off + 1693 * ccomps * dcomps);

            auto g_zz_0_xxzz_xzzzz = cbuffer.data(gh_geom_20_off + 1694 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyyyy = cbuffer.data(gh_geom_20_off + 1695 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyyyz = cbuffer.data(gh_geom_20_off + 1696 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyyzz = cbuffer.data(gh_geom_20_off + 1697 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyzzz = cbuffer.data(gh_geom_20_off + 1698 * ccomps * dcomps);

            auto g_zz_0_xxzz_yzzzz = cbuffer.data(gh_geom_20_off + 1699 * ccomps * dcomps);

            auto g_zz_0_xxzz_zzzzz = cbuffer.data(gh_geom_20_off + 1700 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1701 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1702 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1703 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1704 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1705 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1706 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1707 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1708 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1709 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1710 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1711 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1712 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1713 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1714 * ccomps * dcomps);

            auto g_zz_0_xyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1715 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1716 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1717 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1718 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1719 * ccomps * dcomps);

            auto g_zz_0_xyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1720 * ccomps * dcomps);

            auto g_zz_0_xyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1721 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1722 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1723 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1724 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1725 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1726 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1727 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1728 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1729 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1730 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1731 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1732 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1733 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1734 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1735 * ccomps * dcomps);

            auto g_zz_0_xyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1736 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1737 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1738 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1739 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1740 * ccomps * dcomps);

            auto g_zz_0_xyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1741 * ccomps * dcomps);

            auto g_zz_0_xyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1742 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1743 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1744 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1745 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1746 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1747 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1748 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1749 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1750 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1751 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1752 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1753 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1754 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1755 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1756 * ccomps * dcomps);

            auto g_zz_0_xyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1757 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1758 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1759 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1760 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1761 * ccomps * dcomps);

            auto g_zz_0_xyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1762 * ccomps * dcomps);

            auto g_zz_0_xyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1763 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1764 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1765 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1766 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1767 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1768 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1769 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1770 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1771 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1772 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1773 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1774 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1775 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1776 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1777 * ccomps * dcomps);

            auto g_zz_0_xzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1778 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1779 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1780 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1781 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1782 * ccomps * dcomps);

            auto g_zz_0_xzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1783 * ccomps * dcomps);

            auto g_zz_0_xzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1784 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxxx = cbuffer.data(gh_geom_20_off + 1785 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxxy = cbuffer.data(gh_geom_20_off + 1786 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxxz = cbuffer.data(gh_geom_20_off + 1787 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxyy = cbuffer.data(gh_geom_20_off + 1788 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxyz = cbuffer.data(gh_geom_20_off + 1789 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxxzz = cbuffer.data(gh_geom_20_off + 1790 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxyyy = cbuffer.data(gh_geom_20_off + 1791 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxyyz = cbuffer.data(gh_geom_20_off + 1792 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxyzz = cbuffer.data(gh_geom_20_off + 1793 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxzzz = cbuffer.data(gh_geom_20_off + 1794 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyyyy = cbuffer.data(gh_geom_20_off + 1795 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyyyz = cbuffer.data(gh_geom_20_off + 1796 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyyzz = cbuffer.data(gh_geom_20_off + 1797 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyzzz = cbuffer.data(gh_geom_20_off + 1798 * ccomps * dcomps);

            auto g_zz_0_yyyy_xzzzz = cbuffer.data(gh_geom_20_off + 1799 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyyyy = cbuffer.data(gh_geom_20_off + 1800 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyyyz = cbuffer.data(gh_geom_20_off + 1801 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyyzz = cbuffer.data(gh_geom_20_off + 1802 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyzzz = cbuffer.data(gh_geom_20_off + 1803 * ccomps * dcomps);

            auto g_zz_0_yyyy_yzzzz = cbuffer.data(gh_geom_20_off + 1804 * ccomps * dcomps);

            auto g_zz_0_yyyy_zzzzz = cbuffer.data(gh_geom_20_off + 1805 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxxx = cbuffer.data(gh_geom_20_off + 1806 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxxy = cbuffer.data(gh_geom_20_off + 1807 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxxz = cbuffer.data(gh_geom_20_off + 1808 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxyy = cbuffer.data(gh_geom_20_off + 1809 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxyz = cbuffer.data(gh_geom_20_off + 1810 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxxzz = cbuffer.data(gh_geom_20_off + 1811 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxyyy = cbuffer.data(gh_geom_20_off + 1812 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxyyz = cbuffer.data(gh_geom_20_off + 1813 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxyzz = cbuffer.data(gh_geom_20_off + 1814 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxzzz = cbuffer.data(gh_geom_20_off + 1815 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyyyy = cbuffer.data(gh_geom_20_off + 1816 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyyyz = cbuffer.data(gh_geom_20_off + 1817 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyyzz = cbuffer.data(gh_geom_20_off + 1818 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyzzz = cbuffer.data(gh_geom_20_off + 1819 * ccomps * dcomps);

            auto g_zz_0_yyyz_xzzzz = cbuffer.data(gh_geom_20_off + 1820 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyyyy = cbuffer.data(gh_geom_20_off + 1821 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyyyz = cbuffer.data(gh_geom_20_off + 1822 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyyzz = cbuffer.data(gh_geom_20_off + 1823 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyzzz = cbuffer.data(gh_geom_20_off + 1824 * ccomps * dcomps);

            auto g_zz_0_yyyz_yzzzz = cbuffer.data(gh_geom_20_off + 1825 * ccomps * dcomps);

            auto g_zz_0_yyyz_zzzzz = cbuffer.data(gh_geom_20_off + 1826 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxxx = cbuffer.data(gh_geom_20_off + 1827 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxxy = cbuffer.data(gh_geom_20_off + 1828 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxxz = cbuffer.data(gh_geom_20_off + 1829 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxyy = cbuffer.data(gh_geom_20_off + 1830 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxyz = cbuffer.data(gh_geom_20_off + 1831 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxxzz = cbuffer.data(gh_geom_20_off + 1832 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxyyy = cbuffer.data(gh_geom_20_off + 1833 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxyyz = cbuffer.data(gh_geom_20_off + 1834 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxyzz = cbuffer.data(gh_geom_20_off + 1835 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxzzz = cbuffer.data(gh_geom_20_off + 1836 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyyyy = cbuffer.data(gh_geom_20_off + 1837 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyyyz = cbuffer.data(gh_geom_20_off + 1838 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyyzz = cbuffer.data(gh_geom_20_off + 1839 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyzzz = cbuffer.data(gh_geom_20_off + 1840 * ccomps * dcomps);

            auto g_zz_0_yyzz_xzzzz = cbuffer.data(gh_geom_20_off + 1841 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyyyy = cbuffer.data(gh_geom_20_off + 1842 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyyyz = cbuffer.data(gh_geom_20_off + 1843 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyyzz = cbuffer.data(gh_geom_20_off + 1844 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyzzz = cbuffer.data(gh_geom_20_off + 1845 * ccomps * dcomps);

            auto g_zz_0_yyzz_yzzzz = cbuffer.data(gh_geom_20_off + 1846 * ccomps * dcomps);

            auto g_zz_0_yyzz_zzzzz = cbuffer.data(gh_geom_20_off + 1847 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1848 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1849 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1850 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1851 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1852 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1853 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1854 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1855 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1856 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1857 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1858 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1859 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1860 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1861 * ccomps * dcomps);

            auto g_zz_0_yzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1862 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1863 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1864 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1865 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1866 * ccomps * dcomps);

            auto g_zz_0_yzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1867 * ccomps * dcomps);

            auto g_zz_0_yzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1868 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxxx = cbuffer.data(gh_geom_20_off + 1869 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxxy = cbuffer.data(gh_geom_20_off + 1870 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxxz = cbuffer.data(gh_geom_20_off + 1871 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxyy = cbuffer.data(gh_geom_20_off + 1872 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxyz = cbuffer.data(gh_geom_20_off + 1873 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxxzz = cbuffer.data(gh_geom_20_off + 1874 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxyyy = cbuffer.data(gh_geom_20_off + 1875 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxyyz = cbuffer.data(gh_geom_20_off + 1876 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxyzz = cbuffer.data(gh_geom_20_off + 1877 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxzzz = cbuffer.data(gh_geom_20_off + 1878 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyyyy = cbuffer.data(gh_geom_20_off + 1879 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyyyz = cbuffer.data(gh_geom_20_off + 1880 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyyzz = cbuffer.data(gh_geom_20_off + 1881 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyzzz = cbuffer.data(gh_geom_20_off + 1882 * ccomps * dcomps);

            auto g_zz_0_zzzz_xzzzz = cbuffer.data(gh_geom_20_off + 1883 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyyyy = cbuffer.data(gh_geom_20_off + 1884 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyyyz = cbuffer.data(gh_geom_20_off + 1885 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyyzz = cbuffer.data(gh_geom_20_off + 1886 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyzzz = cbuffer.data(gh_geom_20_off + 1887 * ccomps * dcomps);

            auto g_zz_0_zzzz_yzzzz = cbuffer.data(gh_geom_20_off + 1888 * ccomps * dcomps);

            auto g_zz_0_zzzz_zzzzz = cbuffer.data(gh_geom_20_off + 1889 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hgxx

            const auto hg_geom_20_off = idx_geom_20_hgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzzz, g_xx_0_xxxx_xxxx, g_xx_0_xxxx_xxxxx, g_xx_0_xxxx_xxxxy, g_xx_0_xxxx_xxxxz, g_xx_0_xxxx_xxxy, g_xx_0_xxxx_xxxyy, g_xx_0_xxxx_xxxyz, g_xx_0_xxxx_xxxz, g_xx_0_xxxx_xxxzz, g_xx_0_xxxx_xxyy, g_xx_0_xxxx_xxyyy, g_xx_0_xxxx_xxyyz, g_xx_0_xxxx_xxyz, g_xx_0_xxxx_xxyzz, g_xx_0_xxxx_xxzz, g_xx_0_xxxx_xxzzz, g_xx_0_xxxx_xyyy, g_xx_0_xxxx_xyyyy, g_xx_0_xxxx_xyyyz, g_xx_0_xxxx_xyyz, g_xx_0_xxxx_xyyzz, g_xx_0_xxxx_xyzz, g_xx_0_xxxx_xyzzz, g_xx_0_xxxx_xzzz, g_xx_0_xxxx_xzzzz, g_xx_0_xxxx_yyyy, g_xx_0_xxxx_yyyz, g_xx_0_xxxx_yyzz, g_xx_0_xxxx_yzzz, g_xx_0_xxxx_zzzz, g_xx_0_xxxxx_xxxx, g_xx_0_xxxxx_xxxy, g_xx_0_xxxxx_xxxz, g_xx_0_xxxxx_xxyy, g_xx_0_xxxxx_xxyz, g_xx_0_xxxxx_xxzz, g_xx_0_xxxxx_xyyy, g_xx_0_xxxxx_xyyz, g_xx_0_xxxxx_xyzz, g_xx_0_xxxxx_xzzz, g_xx_0_xxxxx_yyyy, g_xx_0_xxxxx_yyyz, g_xx_0_xxxxx_yyzz, g_xx_0_xxxxx_yzzz, g_xx_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxx_xxxx[k] = -2.0 * g_x_0_xxxx_xxxx[k] - g_xx_0_xxxx_xxxx[k] * ab_x + g_xx_0_xxxx_xxxxx[k];

                g_xx_0_xxxxx_xxxy[k] = -2.0 * g_x_0_xxxx_xxxy[k] - g_xx_0_xxxx_xxxy[k] * ab_x + g_xx_0_xxxx_xxxxy[k];

                g_xx_0_xxxxx_xxxz[k] = -2.0 * g_x_0_xxxx_xxxz[k] - g_xx_0_xxxx_xxxz[k] * ab_x + g_xx_0_xxxx_xxxxz[k];

                g_xx_0_xxxxx_xxyy[k] = -2.0 * g_x_0_xxxx_xxyy[k] - g_xx_0_xxxx_xxyy[k] * ab_x + g_xx_0_xxxx_xxxyy[k];

                g_xx_0_xxxxx_xxyz[k] = -2.0 * g_x_0_xxxx_xxyz[k] - g_xx_0_xxxx_xxyz[k] * ab_x + g_xx_0_xxxx_xxxyz[k];

                g_xx_0_xxxxx_xxzz[k] = -2.0 * g_x_0_xxxx_xxzz[k] - g_xx_0_xxxx_xxzz[k] * ab_x + g_xx_0_xxxx_xxxzz[k];

                g_xx_0_xxxxx_xyyy[k] = -2.0 * g_x_0_xxxx_xyyy[k] - g_xx_0_xxxx_xyyy[k] * ab_x + g_xx_0_xxxx_xxyyy[k];

                g_xx_0_xxxxx_xyyz[k] = -2.0 * g_x_0_xxxx_xyyz[k] - g_xx_0_xxxx_xyyz[k] * ab_x + g_xx_0_xxxx_xxyyz[k];

                g_xx_0_xxxxx_xyzz[k] = -2.0 * g_x_0_xxxx_xyzz[k] - g_xx_0_xxxx_xyzz[k] * ab_x + g_xx_0_xxxx_xxyzz[k];

                g_xx_0_xxxxx_xzzz[k] = -2.0 * g_x_0_xxxx_xzzz[k] - g_xx_0_xxxx_xzzz[k] * ab_x + g_xx_0_xxxx_xxzzz[k];

                g_xx_0_xxxxx_yyyy[k] = -2.0 * g_x_0_xxxx_yyyy[k] - g_xx_0_xxxx_yyyy[k] * ab_x + g_xx_0_xxxx_xyyyy[k];

                g_xx_0_xxxxx_yyyz[k] = -2.0 * g_x_0_xxxx_yyyz[k] - g_xx_0_xxxx_yyyz[k] * ab_x + g_xx_0_xxxx_xyyyz[k];

                g_xx_0_xxxxx_yyzz[k] = -2.0 * g_x_0_xxxx_yyzz[k] - g_xx_0_xxxx_yyzz[k] * ab_x + g_xx_0_xxxx_xyyzz[k];

                g_xx_0_xxxxx_yzzz[k] = -2.0 * g_x_0_xxxx_yzzz[k] - g_xx_0_xxxx_yzzz[k] * ab_x + g_xx_0_xxxx_xyzzz[k];

                g_xx_0_xxxxx_zzzz[k] = -2.0 * g_x_0_xxxx_zzzz[k] - g_xx_0_xxxx_zzzz[k] * ab_x + g_xx_0_xxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_xxxx, g_xx_0_xxxx_xxxxy, g_xx_0_xxxx_xxxy, g_xx_0_xxxx_xxxyy, g_xx_0_xxxx_xxxyz, g_xx_0_xxxx_xxxz, g_xx_0_xxxx_xxyy, g_xx_0_xxxx_xxyyy, g_xx_0_xxxx_xxyyz, g_xx_0_xxxx_xxyz, g_xx_0_xxxx_xxyzz, g_xx_0_xxxx_xxzz, g_xx_0_xxxx_xyyy, g_xx_0_xxxx_xyyyy, g_xx_0_xxxx_xyyyz, g_xx_0_xxxx_xyyz, g_xx_0_xxxx_xyyzz, g_xx_0_xxxx_xyzz, g_xx_0_xxxx_xyzzz, g_xx_0_xxxx_xzzz, g_xx_0_xxxx_yyyy, g_xx_0_xxxx_yyyyy, g_xx_0_xxxx_yyyyz, g_xx_0_xxxx_yyyz, g_xx_0_xxxx_yyyzz, g_xx_0_xxxx_yyzz, g_xx_0_xxxx_yyzzz, g_xx_0_xxxx_yzzz, g_xx_0_xxxx_yzzzz, g_xx_0_xxxx_zzzz, g_xx_0_xxxxy_xxxx, g_xx_0_xxxxy_xxxy, g_xx_0_xxxxy_xxxz, g_xx_0_xxxxy_xxyy, g_xx_0_xxxxy_xxyz, g_xx_0_xxxxy_xxzz, g_xx_0_xxxxy_xyyy, g_xx_0_xxxxy_xyyz, g_xx_0_xxxxy_xyzz, g_xx_0_xxxxy_xzzz, g_xx_0_xxxxy_yyyy, g_xx_0_xxxxy_yyyz, g_xx_0_xxxxy_yyzz, g_xx_0_xxxxy_yzzz, g_xx_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxy_xxxx[k] = -g_xx_0_xxxx_xxxx[k] * ab_y + g_xx_0_xxxx_xxxxy[k];

                g_xx_0_xxxxy_xxxy[k] = -g_xx_0_xxxx_xxxy[k] * ab_y + g_xx_0_xxxx_xxxyy[k];

                g_xx_0_xxxxy_xxxz[k] = -g_xx_0_xxxx_xxxz[k] * ab_y + g_xx_0_xxxx_xxxyz[k];

                g_xx_0_xxxxy_xxyy[k] = -g_xx_0_xxxx_xxyy[k] * ab_y + g_xx_0_xxxx_xxyyy[k];

                g_xx_0_xxxxy_xxyz[k] = -g_xx_0_xxxx_xxyz[k] * ab_y + g_xx_0_xxxx_xxyyz[k];

                g_xx_0_xxxxy_xxzz[k] = -g_xx_0_xxxx_xxzz[k] * ab_y + g_xx_0_xxxx_xxyzz[k];

                g_xx_0_xxxxy_xyyy[k] = -g_xx_0_xxxx_xyyy[k] * ab_y + g_xx_0_xxxx_xyyyy[k];

                g_xx_0_xxxxy_xyyz[k] = -g_xx_0_xxxx_xyyz[k] * ab_y + g_xx_0_xxxx_xyyyz[k];

                g_xx_0_xxxxy_xyzz[k] = -g_xx_0_xxxx_xyzz[k] * ab_y + g_xx_0_xxxx_xyyzz[k];

                g_xx_0_xxxxy_xzzz[k] = -g_xx_0_xxxx_xzzz[k] * ab_y + g_xx_0_xxxx_xyzzz[k];

                g_xx_0_xxxxy_yyyy[k] = -g_xx_0_xxxx_yyyy[k] * ab_y + g_xx_0_xxxx_yyyyy[k];

                g_xx_0_xxxxy_yyyz[k] = -g_xx_0_xxxx_yyyz[k] * ab_y + g_xx_0_xxxx_yyyyz[k];

                g_xx_0_xxxxy_yyzz[k] = -g_xx_0_xxxx_yyzz[k] * ab_y + g_xx_0_xxxx_yyyzz[k];

                g_xx_0_xxxxy_yzzz[k] = -g_xx_0_xxxx_yzzz[k] * ab_y + g_xx_0_xxxx_yyzzz[k];

                g_xx_0_xxxxy_zzzz[k] = -g_xx_0_xxxx_zzzz[k] * ab_y + g_xx_0_xxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_xxxx, g_xx_0_xxxx_xxxxz, g_xx_0_xxxx_xxxy, g_xx_0_xxxx_xxxyz, g_xx_0_xxxx_xxxz, g_xx_0_xxxx_xxxzz, g_xx_0_xxxx_xxyy, g_xx_0_xxxx_xxyyz, g_xx_0_xxxx_xxyz, g_xx_0_xxxx_xxyzz, g_xx_0_xxxx_xxzz, g_xx_0_xxxx_xxzzz, g_xx_0_xxxx_xyyy, g_xx_0_xxxx_xyyyz, g_xx_0_xxxx_xyyz, g_xx_0_xxxx_xyyzz, g_xx_0_xxxx_xyzz, g_xx_0_xxxx_xyzzz, g_xx_0_xxxx_xzzz, g_xx_0_xxxx_xzzzz, g_xx_0_xxxx_yyyy, g_xx_0_xxxx_yyyyz, g_xx_0_xxxx_yyyz, g_xx_0_xxxx_yyyzz, g_xx_0_xxxx_yyzz, g_xx_0_xxxx_yyzzz, g_xx_0_xxxx_yzzz, g_xx_0_xxxx_yzzzz, g_xx_0_xxxx_zzzz, g_xx_0_xxxx_zzzzz, g_xx_0_xxxxz_xxxx, g_xx_0_xxxxz_xxxy, g_xx_0_xxxxz_xxxz, g_xx_0_xxxxz_xxyy, g_xx_0_xxxxz_xxyz, g_xx_0_xxxxz_xxzz, g_xx_0_xxxxz_xyyy, g_xx_0_xxxxz_xyyz, g_xx_0_xxxxz_xyzz, g_xx_0_xxxxz_xzzz, g_xx_0_xxxxz_yyyy, g_xx_0_xxxxz_yyyz, g_xx_0_xxxxz_yyzz, g_xx_0_xxxxz_yzzz, g_xx_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxz_xxxx[k] = -g_xx_0_xxxx_xxxx[k] * ab_z + g_xx_0_xxxx_xxxxz[k];

                g_xx_0_xxxxz_xxxy[k] = -g_xx_0_xxxx_xxxy[k] * ab_z + g_xx_0_xxxx_xxxyz[k];

                g_xx_0_xxxxz_xxxz[k] = -g_xx_0_xxxx_xxxz[k] * ab_z + g_xx_0_xxxx_xxxzz[k];

                g_xx_0_xxxxz_xxyy[k] = -g_xx_0_xxxx_xxyy[k] * ab_z + g_xx_0_xxxx_xxyyz[k];

                g_xx_0_xxxxz_xxyz[k] = -g_xx_0_xxxx_xxyz[k] * ab_z + g_xx_0_xxxx_xxyzz[k];

                g_xx_0_xxxxz_xxzz[k] = -g_xx_0_xxxx_xxzz[k] * ab_z + g_xx_0_xxxx_xxzzz[k];

                g_xx_0_xxxxz_xyyy[k] = -g_xx_0_xxxx_xyyy[k] * ab_z + g_xx_0_xxxx_xyyyz[k];

                g_xx_0_xxxxz_xyyz[k] = -g_xx_0_xxxx_xyyz[k] * ab_z + g_xx_0_xxxx_xyyzz[k];

                g_xx_0_xxxxz_xyzz[k] = -g_xx_0_xxxx_xyzz[k] * ab_z + g_xx_0_xxxx_xyzzz[k];

                g_xx_0_xxxxz_xzzz[k] = -g_xx_0_xxxx_xzzz[k] * ab_z + g_xx_0_xxxx_xzzzz[k];

                g_xx_0_xxxxz_yyyy[k] = -g_xx_0_xxxx_yyyy[k] * ab_z + g_xx_0_xxxx_yyyyz[k];

                g_xx_0_xxxxz_yyyz[k] = -g_xx_0_xxxx_yyyz[k] * ab_z + g_xx_0_xxxx_yyyzz[k];

                g_xx_0_xxxxz_yyzz[k] = -g_xx_0_xxxx_yyzz[k] * ab_z + g_xx_0_xxxx_yyzzz[k];

                g_xx_0_xxxxz_yzzz[k] = -g_xx_0_xxxx_yzzz[k] * ab_z + g_xx_0_xxxx_yzzzz[k];

                g_xx_0_xxxxz_zzzz[k] = -g_xx_0_xxxx_zzzz[k] * ab_z + g_xx_0_xxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxy_xxxx, g_xx_0_xxxy_xxxxy, g_xx_0_xxxy_xxxy, g_xx_0_xxxy_xxxyy, g_xx_0_xxxy_xxxyz, g_xx_0_xxxy_xxxz, g_xx_0_xxxy_xxyy, g_xx_0_xxxy_xxyyy, g_xx_0_xxxy_xxyyz, g_xx_0_xxxy_xxyz, g_xx_0_xxxy_xxyzz, g_xx_0_xxxy_xxzz, g_xx_0_xxxy_xyyy, g_xx_0_xxxy_xyyyy, g_xx_0_xxxy_xyyyz, g_xx_0_xxxy_xyyz, g_xx_0_xxxy_xyyzz, g_xx_0_xxxy_xyzz, g_xx_0_xxxy_xyzzz, g_xx_0_xxxy_xzzz, g_xx_0_xxxy_yyyy, g_xx_0_xxxy_yyyyy, g_xx_0_xxxy_yyyyz, g_xx_0_xxxy_yyyz, g_xx_0_xxxy_yyyzz, g_xx_0_xxxy_yyzz, g_xx_0_xxxy_yyzzz, g_xx_0_xxxy_yzzz, g_xx_0_xxxy_yzzzz, g_xx_0_xxxy_zzzz, g_xx_0_xxxyy_xxxx, g_xx_0_xxxyy_xxxy, g_xx_0_xxxyy_xxxz, g_xx_0_xxxyy_xxyy, g_xx_0_xxxyy_xxyz, g_xx_0_xxxyy_xxzz, g_xx_0_xxxyy_xyyy, g_xx_0_xxxyy_xyyz, g_xx_0_xxxyy_xyzz, g_xx_0_xxxyy_xzzz, g_xx_0_xxxyy_yyyy, g_xx_0_xxxyy_yyyz, g_xx_0_xxxyy_yyzz, g_xx_0_xxxyy_yzzz, g_xx_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyy_xxxx[k] = -g_xx_0_xxxy_xxxx[k] * ab_y + g_xx_0_xxxy_xxxxy[k];

                g_xx_0_xxxyy_xxxy[k] = -g_xx_0_xxxy_xxxy[k] * ab_y + g_xx_0_xxxy_xxxyy[k];

                g_xx_0_xxxyy_xxxz[k] = -g_xx_0_xxxy_xxxz[k] * ab_y + g_xx_0_xxxy_xxxyz[k];

                g_xx_0_xxxyy_xxyy[k] = -g_xx_0_xxxy_xxyy[k] * ab_y + g_xx_0_xxxy_xxyyy[k];

                g_xx_0_xxxyy_xxyz[k] = -g_xx_0_xxxy_xxyz[k] * ab_y + g_xx_0_xxxy_xxyyz[k];

                g_xx_0_xxxyy_xxzz[k] = -g_xx_0_xxxy_xxzz[k] * ab_y + g_xx_0_xxxy_xxyzz[k];

                g_xx_0_xxxyy_xyyy[k] = -g_xx_0_xxxy_xyyy[k] * ab_y + g_xx_0_xxxy_xyyyy[k];

                g_xx_0_xxxyy_xyyz[k] = -g_xx_0_xxxy_xyyz[k] * ab_y + g_xx_0_xxxy_xyyyz[k];

                g_xx_0_xxxyy_xyzz[k] = -g_xx_0_xxxy_xyzz[k] * ab_y + g_xx_0_xxxy_xyyzz[k];

                g_xx_0_xxxyy_xzzz[k] = -g_xx_0_xxxy_xzzz[k] * ab_y + g_xx_0_xxxy_xyzzz[k];

                g_xx_0_xxxyy_yyyy[k] = -g_xx_0_xxxy_yyyy[k] * ab_y + g_xx_0_xxxy_yyyyy[k];

                g_xx_0_xxxyy_yyyz[k] = -g_xx_0_xxxy_yyyz[k] * ab_y + g_xx_0_xxxy_yyyyz[k];

                g_xx_0_xxxyy_yyzz[k] = -g_xx_0_xxxy_yyzz[k] * ab_y + g_xx_0_xxxy_yyyzz[k];

                g_xx_0_xxxyy_yzzz[k] = -g_xx_0_xxxy_yzzz[k] * ab_y + g_xx_0_xxxy_yyzzz[k];

                g_xx_0_xxxyy_zzzz[k] = -g_xx_0_xxxy_zzzz[k] * ab_y + g_xx_0_xxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyz_xxxx, g_xx_0_xxxyz_xxxy, g_xx_0_xxxyz_xxxz, g_xx_0_xxxyz_xxyy, g_xx_0_xxxyz_xxyz, g_xx_0_xxxyz_xxzz, g_xx_0_xxxyz_xyyy, g_xx_0_xxxyz_xyyz, g_xx_0_xxxyz_xyzz, g_xx_0_xxxyz_xzzz, g_xx_0_xxxyz_yyyy, g_xx_0_xxxyz_yyyz, g_xx_0_xxxyz_yyzz, g_xx_0_xxxyz_yzzz, g_xx_0_xxxyz_zzzz, g_xx_0_xxxz_xxxx, g_xx_0_xxxz_xxxxy, g_xx_0_xxxz_xxxy, g_xx_0_xxxz_xxxyy, g_xx_0_xxxz_xxxyz, g_xx_0_xxxz_xxxz, g_xx_0_xxxz_xxyy, g_xx_0_xxxz_xxyyy, g_xx_0_xxxz_xxyyz, g_xx_0_xxxz_xxyz, g_xx_0_xxxz_xxyzz, g_xx_0_xxxz_xxzz, g_xx_0_xxxz_xyyy, g_xx_0_xxxz_xyyyy, g_xx_0_xxxz_xyyyz, g_xx_0_xxxz_xyyz, g_xx_0_xxxz_xyyzz, g_xx_0_xxxz_xyzz, g_xx_0_xxxz_xyzzz, g_xx_0_xxxz_xzzz, g_xx_0_xxxz_yyyy, g_xx_0_xxxz_yyyyy, g_xx_0_xxxz_yyyyz, g_xx_0_xxxz_yyyz, g_xx_0_xxxz_yyyzz, g_xx_0_xxxz_yyzz, g_xx_0_xxxz_yyzzz, g_xx_0_xxxz_yzzz, g_xx_0_xxxz_yzzzz, g_xx_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyz_xxxx[k] = -g_xx_0_xxxz_xxxx[k] * ab_y + g_xx_0_xxxz_xxxxy[k];

                g_xx_0_xxxyz_xxxy[k] = -g_xx_0_xxxz_xxxy[k] * ab_y + g_xx_0_xxxz_xxxyy[k];

                g_xx_0_xxxyz_xxxz[k] = -g_xx_0_xxxz_xxxz[k] * ab_y + g_xx_0_xxxz_xxxyz[k];

                g_xx_0_xxxyz_xxyy[k] = -g_xx_0_xxxz_xxyy[k] * ab_y + g_xx_0_xxxz_xxyyy[k];

                g_xx_0_xxxyz_xxyz[k] = -g_xx_0_xxxz_xxyz[k] * ab_y + g_xx_0_xxxz_xxyyz[k];

                g_xx_0_xxxyz_xxzz[k] = -g_xx_0_xxxz_xxzz[k] * ab_y + g_xx_0_xxxz_xxyzz[k];

                g_xx_0_xxxyz_xyyy[k] = -g_xx_0_xxxz_xyyy[k] * ab_y + g_xx_0_xxxz_xyyyy[k];

                g_xx_0_xxxyz_xyyz[k] = -g_xx_0_xxxz_xyyz[k] * ab_y + g_xx_0_xxxz_xyyyz[k];

                g_xx_0_xxxyz_xyzz[k] = -g_xx_0_xxxz_xyzz[k] * ab_y + g_xx_0_xxxz_xyyzz[k];

                g_xx_0_xxxyz_xzzz[k] = -g_xx_0_xxxz_xzzz[k] * ab_y + g_xx_0_xxxz_xyzzz[k];

                g_xx_0_xxxyz_yyyy[k] = -g_xx_0_xxxz_yyyy[k] * ab_y + g_xx_0_xxxz_yyyyy[k];

                g_xx_0_xxxyz_yyyz[k] = -g_xx_0_xxxz_yyyz[k] * ab_y + g_xx_0_xxxz_yyyyz[k];

                g_xx_0_xxxyz_yyzz[k] = -g_xx_0_xxxz_yyzz[k] * ab_y + g_xx_0_xxxz_yyyzz[k];

                g_xx_0_xxxyz_yzzz[k] = -g_xx_0_xxxz_yzzz[k] * ab_y + g_xx_0_xxxz_yyzzz[k];

                g_xx_0_xxxyz_zzzz[k] = -g_xx_0_xxxz_zzzz[k] * ab_y + g_xx_0_xxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxz_xxxx, g_xx_0_xxxz_xxxxz, g_xx_0_xxxz_xxxy, g_xx_0_xxxz_xxxyz, g_xx_0_xxxz_xxxz, g_xx_0_xxxz_xxxzz, g_xx_0_xxxz_xxyy, g_xx_0_xxxz_xxyyz, g_xx_0_xxxz_xxyz, g_xx_0_xxxz_xxyzz, g_xx_0_xxxz_xxzz, g_xx_0_xxxz_xxzzz, g_xx_0_xxxz_xyyy, g_xx_0_xxxz_xyyyz, g_xx_0_xxxz_xyyz, g_xx_0_xxxz_xyyzz, g_xx_0_xxxz_xyzz, g_xx_0_xxxz_xyzzz, g_xx_0_xxxz_xzzz, g_xx_0_xxxz_xzzzz, g_xx_0_xxxz_yyyy, g_xx_0_xxxz_yyyyz, g_xx_0_xxxz_yyyz, g_xx_0_xxxz_yyyzz, g_xx_0_xxxz_yyzz, g_xx_0_xxxz_yyzzz, g_xx_0_xxxz_yzzz, g_xx_0_xxxz_yzzzz, g_xx_0_xxxz_zzzz, g_xx_0_xxxz_zzzzz, g_xx_0_xxxzz_xxxx, g_xx_0_xxxzz_xxxy, g_xx_0_xxxzz_xxxz, g_xx_0_xxxzz_xxyy, g_xx_0_xxxzz_xxyz, g_xx_0_xxxzz_xxzz, g_xx_0_xxxzz_xyyy, g_xx_0_xxxzz_xyyz, g_xx_0_xxxzz_xyzz, g_xx_0_xxxzz_xzzz, g_xx_0_xxxzz_yyyy, g_xx_0_xxxzz_yyyz, g_xx_0_xxxzz_yyzz, g_xx_0_xxxzz_yzzz, g_xx_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzz_xxxx[k] = -g_xx_0_xxxz_xxxx[k] * ab_z + g_xx_0_xxxz_xxxxz[k];

                g_xx_0_xxxzz_xxxy[k] = -g_xx_0_xxxz_xxxy[k] * ab_z + g_xx_0_xxxz_xxxyz[k];

                g_xx_0_xxxzz_xxxz[k] = -g_xx_0_xxxz_xxxz[k] * ab_z + g_xx_0_xxxz_xxxzz[k];

                g_xx_0_xxxzz_xxyy[k] = -g_xx_0_xxxz_xxyy[k] * ab_z + g_xx_0_xxxz_xxyyz[k];

                g_xx_0_xxxzz_xxyz[k] = -g_xx_0_xxxz_xxyz[k] * ab_z + g_xx_0_xxxz_xxyzz[k];

                g_xx_0_xxxzz_xxzz[k] = -g_xx_0_xxxz_xxzz[k] * ab_z + g_xx_0_xxxz_xxzzz[k];

                g_xx_0_xxxzz_xyyy[k] = -g_xx_0_xxxz_xyyy[k] * ab_z + g_xx_0_xxxz_xyyyz[k];

                g_xx_0_xxxzz_xyyz[k] = -g_xx_0_xxxz_xyyz[k] * ab_z + g_xx_0_xxxz_xyyzz[k];

                g_xx_0_xxxzz_xyzz[k] = -g_xx_0_xxxz_xyzz[k] * ab_z + g_xx_0_xxxz_xyzzz[k];

                g_xx_0_xxxzz_xzzz[k] = -g_xx_0_xxxz_xzzz[k] * ab_z + g_xx_0_xxxz_xzzzz[k];

                g_xx_0_xxxzz_yyyy[k] = -g_xx_0_xxxz_yyyy[k] * ab_z + g_xx_0_xxxz_yyyyz[k];

                g_xx_0_xxxzz_yyyz[k] = -g_xx_0_xxxz_yyyz[k] * ab_z + g_xx_0_xxxz_yyyzz[k];

                g_xx_0_xxxzz_yyzz[k] = -g_xx_0_xxxz_yyzz[k] * ab_z + g_xx_0_xxxz_yyzzz[k];

                g_xx_0_xxxzz_yzzz[k] = -g_xx_0_xxxz_yzzz[k] * ab_z + g_xx_0_xxxz_yzzzz[k];

                g_xx_0_xxxzz_zzzz[k] = -g_xx_0_xxxz_zzzz[k] * ab_z + g_xx_0_xxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyy_xxxx, g_xx_0_xxyy_xxxxy, g_xx_0_xxyy_xxxy, g_xx_0_xxyy_xxxyy, g_xx_0_xxyy_xxxyz, g_xx_0_xxyy_xxxz, g_xx_0_xxyy_xxyy, g_xx_0_xxyy_xxyyy, g_xx_0_xxyy_xxyyz, g_xx_0_xxyy_xxyz, g_xx_0_xxyy_xxyzz, g_xx_0_xxyy_xxzz, g_xx_0_xxyy_xyyy, g_xx_0_xxyy_xyyyy, g_xx_0_xxyy_xyyyz, g_xx_0_xxyy_xyyz, g_xx_0_xxyy_xyyzz, g_xx_0_xxyy_xyzz, g_xx_0_xxyy_xyzzz, g_xx_0_xxyy_xzzz, g_xx_0_xxyy_yyyy, g_xx_0_xxyy_yyyyy, g_xx_0_xxyy_yyyyz, g_xx_0_xxyy_yyyz, g_xx_0_xxyy_yyyzz, g_xx_0_xxyy_yyzz, g_xx_0_xxyy_yyzzz, g_xx_0_xxyy_yzzz, g_xx_0_xxyy_yzzzz, g_xx_0_xxyy_zzzz, g_xx_0_xxyyy_xxxx, g_xx_0_xxyyy_xxxy, g_xx_0_xxyyy_xxxz, g_xx_0_xxyyy_xxyy, g_xx_0_xxyyy_xxyz, g_xx_0_xxyyy_xxzz, g_xx_0_xxyyy_xyyy, g_xx_0_xxyyy_xyyz, g_xx_0_xxyyy_xyzz, g_xx_0_xxyyy_xzzz, g_xx_0_xxyyy_yyyy, g_xx_0_xxyyy_yyyz, g_xx_0_xxyyy_yyzz, g_xx_0_xxyyy_yzzz, g_xx_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyy_xxxx[k] = -g_xx_0_xxyy_xxxx[k] * ab_y + g_xx_0_xxyy_xxxxy[k];

                g_xx_0_xxyyy_xxxy[k] = -g_xx_0_xxyy_xxxy[k] * ab_y + g_xx_0_xxyy_xxxyy[k];

                g_xx_0_xxyyy_xxxz[k] = -g_xx_0_xxyy_xxxz[k] * ab_y + g_xx_0_xxyy_xxxyz[k];

                g_xx_0_xxyyy_xxyy[k] = -g_xx_0_xxyy_xxyy[k] * ab_y + g_xx_0_xxyy_xxyyy[k];

                g_xx_0_xxyyy_xxyz[k] = -g_xx_0_xxyy_xxyz[k] * ab_y + g_xx_0_xxyy_xxyyz[k];

                g_xx_0_xxyyy_xxzz[k] = -g_xx_0_xxyy_xxzz[k] * ab_y + g_xx_0_xxyy_xxyzz[k];

                g_xx_0_xxyyy_xyyy[k] = -g_xx_0_xxyy_xyyy[k] * ab_y + g_xx_0_xxyy_xyyyy[k];

                g_xx_0_xxyyy_xyyz[k] = -g_xx_0_xxyy_xyyz[k] * ab_y + g_xx_0_xxyy_xyyyz[k];

                g_xx_0_xxyyy_xyzz[k] = -g_xx_0_xxyy_xyzz[k] * ab_y + g_xx_0_xxyy_xyyzz[k];

                g_xx_0_xxyyy_xzzz[k] = -g_xx_0_xxyy_xzzz[k] * ab_y + g_xx_0_xxyy_xyzzz[k];

                g_xx_0_xxyyy_yyyy[k] = -g_xx_0_xxyy_yyyy[k] * ab_y + g_xx_0_xxyy_yyyyy[k];

                g_xx_0_xxyyy_yyyz[k] = -g_xx_0_xxyy_yyyz[k] * ab_y + g_xx_0_xxyy_yyyyz[k];

                g_xx_0_xxyyy_yyzz[k] = -g_xx_0_xxyy_yyzz[k] * ab_y + g_xx_0_xxyy_yyyzz[k];

                g_xx_0_xxyyy_yzzz[k] = -g_xx_0_xxyy_yzzz[k] * ab_y + g_xx_0_xxyy_yyzzz[k];

                g_xx_0_xxyyy_zzzz[k] = -g_xx_0_xxyy_zzzz[k] * ab_y + g_xx_0_xxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyz_xxxx, g_xx_0_xxyyz_xxxy, g_xx_0_xxyyz_xxxz, g_xx_0_xxyyz_xxyy, g_xx_0_xxyyz_xxyz, g_xx_0_xxyyz_xxzz, g_xx_0_xxyyz_xyyy, g_xx_0_xxyyz_xyyz, g_xx_0_xxyyz_xyzz, g_xx_0_xxyyz_xzzz, g_xx_0_xxyyz_yyyy, g_xx_0_xxyyz_yyyz, g_xx_0_xxyyz_yyzz, g_xx_0_xxyyz_yzzz, g_xx_0_xxyyz_zzzz, g_xx_0_xxyz_xxxx, g_xx_0_xxyz_xxxxy, g_xx_0_xxyz_xxxy, g_xx_0_xxyz_xxxyy, g_xx_0_xxyz_xxxyz, g_xx_0_xxyz_xxxz, g_xx_0_xxyz_xxyy, g_xx_0_xxyz_xxyyy, g_xx_0_xxyz_xxyyz, g_xx_0_xxyz_xxyz, g_xx_0_xxyz_xxyzz, g_xx_0_xxyz_xxzz, g_xx_0_xxyz_xyyy, g_xx_0_xxyz_xyyyy, g_xx_0_xxyz_xyyyz, g_xx_0_xxyz_xyyz, g_xx_0_xxyz_xyyzz, g_xx_0_xxyz_xyzz, g_xx_0_xxyz_xyzzz, g_xx_0_xxyz_xzzz, g_xx_0_xxyz_yyyy, g_xx_0_xxyz_yyyyy, g_xx_0_xxyz_yyyyz, g_xx_0_xxyz_yyyz, g_xx_0_xxyz_yyyzz, g_xx_0_xxyz_yyzz, g_xx_0_xxyz_yyzzz, g_xx_0_xxyz_yzzz, g_xx_0_xxyz_yzzzz, g_xx_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyz_xxxx[k] = -g_xx_0_xxyz_xxxx[k] * ab_y + g_xx_0_xxyz_xxxxy[k];

                g_xx_0_xxyyz_xxxy[k] = -g_xx_0_xxyz_xxxy[k] * ab_y + g_xx_0_xxyz_xxxyy[k];

                g_xx_0_xxyyz_xxxz[k] = -g_xx_0_xxyz_xxxz[k] * ab_y + g_xx_0_xxyz_xxxyz[k];

                g_xx_0_xxyyz_xxyy[k] = -g_xx_0_xxyz_xxyy[k] * ab_y + g_xx_0_xxyz_xxyyy[k];

                g_xx_0_xxyyz_xxyz[k] = -g_xx_0_xxyz_xxyz[k] * ab_y + g_xx_0_xxyz_xxyyz[k];

                g_xx_0_xxyyz_xxzz[k] = -g_xx_0_xxyz_xxzz[k] * ab_y + g_xx_0_xxyz_xxyzz[k];

                g_xx_0_xxyyz_xyyy[k] = -g_xx_0_xxyz_xyyy[k] * ab_y + g_xx_0_xxyz_xyyyy[k];

                g_xx_0_xxyyz_xyyz[k] = -g_xx_0_xxyz_xyyz[k] * ab_y + g_xx_0_xxyz_xyyyz[k];

                g_xx_0_xxyyz_xyzz[k] = -g_xx_0_xxyz_xyzz[k] * ab_y + g_xx_0_xxyz_xyyzz[k];

                g_xx_0_xxyyz_xzzz[k] = -g_xx_0_xxyz_xzzz[k] * ab_y + g_xx_0_xxyz_xyzzz[k];

                g_xx_0_xxyyz_yyyy[k] = -g_xx_0_xxyz_yyyy[k] * ab_y + g_xx_0_xxyz_yyyyy[k];

                g_xx_0_xxyyz_yyyz[k] = -g_xx_0_xxyz_yyyz[k] * ab_y + g_xx_0_xxyz_yyyyz[k];

                g_xx_0_xxyyz_yyzz[k] = -g_xx_0_xxyz_yyzz[k] * ab_y + g_xx_0_xxyz_yyyzz[k];

                g_xx_0_xxyyz_yzzz[k] = -g_xx_0_xxyz_yzzz[k] * ab_y + g_xx_0_xxyz_yyzzz[k];

                g_xx_0_xxyyz_zzzz[k] = -g_xx_0_xxyz_zzzz[k] * ab_y + g_xx_0_xxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzz_xxxx, g_xx_0_xxyzz_xxxy, g_xx_0_xxyzz_xxxz, g_xx_0_xxyzz_xxyy, g_xx_0_xxyzz_xxyz, g_xx_0_xxyzz_xxzz, g_xx_0_xxyzz_xyyy, g_xx_0_xxyzz_xyyz, g_xx_0_xxyzz_xyzz, g_xx_0_xxyzz_xzzz, g_xx_0_xxyzz_yyyy, g_xx_0_xxyzz_yyyz, g_xx_0_xxyzz_yyzz, g_xx_0_xxyzz_yzzz, g_xx_0_xxyzz_zzzz, g_xx_0_xxzz_xxxx, g_xx_0_xxzz_xxxxy, g_xx_0_xxzz_xxxy, g_xx_0_xxzz_xxxyy, g_xx_0_xxzz_xxxyz, g_xx_0_xxzz_xxxz, g_xx_0_xxzz_xxyy, g_xx_0_xxzz_xxyyy, g_xx_0_xxzz_xxyyz, g_xx_0_xxzz_xxyz, g_xx_0_xxzz_xxyzz, g_xx_0_xxzz_xxzz, g_xx_0_xxzz_xyyy, g_xx_0_xxzz_xyyyy, g_xx_0_xxzz_xyyyz, g_xx_0_xxzz_xyyz, g_xx_0_xxzz_xyyzz, g_xx_0_xxzz_xyzz, g_xx_0_xxzz_xyzzz, g_xx_0_xxzz_xzzz, g_xx_0_xxzz_yyyy, g_xx_0_xxzz_yyyyy, g_xx_0_xxzz_yyyyz, g_xx_0_xxzz_yyyz, g_xx_0_xxzz_yyyzz, g_xx_0_xxzz_yyzz, g_xx_0_xxzz_yyzzz, g_xx_0_xxzz_yzzz, g_xx_0_xxzz_yzzzz, g_xx_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzz_xxxx[k] = -g_xx_0_xxzz_xxxx[k] * ab_y + g_xx_0_xxzz_xxxxy[k];

                g_xx_0_xxyzz_xxxy[k] = -g_xx_0_xxzz_xxxy[k] * ab_y + g_xx_0_xxzz_xxxyy[k];

                g_xx_0_xxyzz_xxxz[k] = -g_xx_0_xxzz_xxxz[k] * ab_y + g_xx_0_xxzz_xxxyz[k];

                g_xx_0_xxyzz_xxyy[k] = -g_xx_0_xxzz_xxyy[k] * ab_y + g_xx_0_xxzz_xxyyy[k];

                g_xx_0_xxyzz_xxyz[k] = -g_xx_0_xxzz_xxyz[k] * ab_y + g_xx_0_xxzz_xxyyz[k];

                g_xx_0_xxyzz_xxzz[k] = -g_xx_0_xxzz_xxzz[k] * ab_y + g_xx_0_xxzz_xxyzz[k];

                g_xx_0_xxyzz_xyyy[k] = -g_xx_0_xxzz_xyyy[k] * ab_y + g_xx_0_xxzz_xyyyy[k];

                g_xx_0_xxyzz_xyyz[k] = -g_xx_0_xxzz_xyyz[k] * ab_y + g_xx_0_xxzz_xyyyz[k];

                g_xx_0_xxyzz_xyzz[k] = -g_xx_0_xxzz_xyzz[k] * ab_y + g_xx_0_xxzz_xyyzz[k];

                g_xx_0_xxyzz_xzzz[k] = -g_xx_0_xxzz_xzzz[k] * ab_y + g_xx_0_xxzz_xyzzz[k];

                g_xx_0_xxyzz_yyyy[k] = -g_xx_0_xxzz_yyyy[k] * ab_y + g_xx_0_xxzz_yyyyy[k];

                g_xx_0_xxyzz_yyyz[k] = -g_xx_0_xxzz_yyyz[k] * ab_y + g_xx_0_xxzz_yyyyz[k];

                g_xx_0_xxyzz_yyzz[k] = -g_xx_0_xxzz_yyzz[k] * ab_y + g_xx_0_xxzz_yyyzz[k];

                g_xx_0_xxyzz_yzzz[k] = -g_xx_0_xxzz_yzzz[k] * ab_y + g_xx_0_xxzz_yyzzz[k];

                g_xx_0_xxyzz_zzzz[k] = -g_xx_0_xxzz_zzzz[k] * ab_y + g_xx_0_xxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzz_xxxx, g_xx_0_xxzz_xxxxz, g_xx_0_xxzz_xxxy, g_xx_0_xxzz_xxxyz, g_xx_0_xxzz_xxxz, g_xx_0_xxzz_xxxzz, g_xx_0_xxzz_xxyy, g_xx_0_xxzz_xxyyz, g_xx_0_xxzz_xxyz, g_xx_0_xxzz_xxyzz, g_xx_0_xxzz_xxzz, g_xx_0_xxzz_xxzzz, g_xx_0_xxzz_xyyy, g_xx_0_xxzz_xyyyz, g_xx_0_xxzz_xyyz, g_xx_0_xxzz_xyyzz, g_xx_0_xxzz_xyzz, g_xx_0_xxzz_xyzzz, g_xx_0_xxzz_xzzz, g_xx_0_xxzz_xzzzz, g_xx_0_xxzz_yyyy, g_xx_0_xxzz_yyyyz, g_xx_0_xxzz_yyyz, g_xx_0_xxzz_yyyzz, g_xx_0_xxzz_yyzz, g_xx_0_xxzz_yyzzz, g_xx_0_xxzz_yzzz, g_xx_0_xxzz_yzzzz, g_xx_0_xxzz_zzzz, g_xx_0_xxzz_zzzzz, g_xx_0_xxzzz_xxxx, g_xx_0_xxzzz_xxxy, g_xx_0_xxzzz_xxxz, g_xx_0_xxzzz_xxyy, g_xx_0_xxzzz_xxyz, g_xx_0_xxzzz_xxzz, g_xx_0_xxzzz_xyyy, g_xx_0_xxzzz_xyyz, g_xx_0_xxzzz_xyzz, g_xx_0_xxzzz_xzzz, g_xx_0_xxzzz_yyyy, g_xx_0_xxzzz_yyyz, g_xx_0_xxzzz_yyzz, g_xx_0_xxzzz_yzzz, g_xx_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzz_xxxx[k] = -g_xx_0_xxzz_xxxx[k] * ab_z + g_xx_0_xxzz_xxxxz[k];

                g_xx_0_xxzzz_xxxy[k] = -g_xx_0_xxzz_xxxy[k] * ab_z + g_xx_0_xxzz_xxxyz[k];

                g_xx_0_xxzzz_xxxz[k] = -g_xx_0_xxzz_xxxz[k] * ab_z + g_xx_0_xxzz_xxxzz[k];

                g_xx_0_xxzzz_xxyy[k] = -g_xx_0_xxzz_xxyy[k] * ab_z + g_xx_0_xxzz_xxyyz[k];

                g_xx_0_xxzzz_xxyz[k] = -g_xx_0_xxzz_xxyz[k] * ab_z + g_xx_0_xxzz_xxyzz[k];

                g_xx_0_xxzzz_xxzz[k] = -g_xx_0_xxzz_xxzz[k] * ab_z + g_xx_0_xxzz_xxzzz[k];

                g_xx_0_xxzzz_xyyy[k] = -g_xx_0_xxzz_xyyy[k] * ab_z + g_xx_0_xxzz_xyyyz[k];

                g_xx_0_xxzzz_xyyz[k] = -g_xx_0_xxzz_xyyz[k] * ab_z + g_xx_0_xxzz_xyyzz[k];

                g_xx_0_xxzzz_xyzz[k] = -g_xx_0_xxzz_xyzz[k] * ab_z + g_xx_0_xxzz_xyzzz[k];

                g_xx_0_xxzzz_xzzz[k] = -g_xx_0_xxzz_xzzz[k] * ab_z + g_xx_0_xxzz_xzzzz[k];

                g_xx_0_xxzzz_yyyy[k] = -g_xx_0_xxzz_yyyy[k] * ab_z + g_xx_0_xxzz_yyyyz[k];

                g_xx_0_xxzzz_yyyz[k] = -g_xx_0_xxzz_yyyz[k] * ab_z + g_xx_0_xxzz_yyyzz[k];

                g_xx_0_xxzzz_yyzz[k] = -g_xx_0_xxzz_yyzz[k] * ab_z + g_xx_0_xxzz_yyzzz[k];

                g_xx_0_xxzzz_yzzz[k] = -g_xx_0_xxzz_yzzz[k] * ab_z + g_xx_0_xxzz_yzzzz[k];

                g_xx_0_xxzzz_zzzz[k] = -g_xx_0_xxzz_zzzz[k] * ab_z + g_xx_0_xxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyy_xxxx, g_xx_0_xyyy_xxxxy, g_xx_0_xyyy_xxxy, g_xx_0_xyyy_xxxyy, g_xx_0_xyyy_xxxyz, g_xx_0_xyyy_xxxz, g_xx_0_xyyy_xxyy, g_xx_0_xyyy_xxyyy, g_xx_0_xyyy_xxyyz, g_xx_0_xyyy_xxyz, g_xx_0_xyyy_xxyzz, g_xx_0_xyyy_xxzz, g_xx_0_xyyy_xyyy, g_xx_0_xyyy_xyyyy, g_xx_0_xyyy_xyyyz, g_xx_0_xyyy_xyyz, g_xx_0_xyyy_xyyzz, g_xx_0_xyyy_xyzz, g_xx_0_xyyy_xyzzz, g_xx_0_xyyy_xzzz, g_xx_0_xyyy_yyyy, g_xx_0_xyyy_yyyyy, g_xx_0_xyyy_yyyyz, g_xx_0_xyyy_yyyz, g_xx_0_xyyy_yyyzz, g_xx_0_xyyy_yyzz, g_xx_0_xyyy_yyzzz, g_xx_0_xyyy_yzzz, g_xx_0_xyyy_yzzzz, g_xx_0_xyyy_zzzz, g_xx_0_xyyyy_xxxx, g_xx_0_xyyyy_xxxy, g_xx_0_xyyyy_xxxz, g_xx_0_xyyyy_xxyy, g_xx_0_xyyyy_xxyz, g_xx_0_xyyyy_xxzz, g_xx_0_xyyyy_xyyy, g_xx_0_xyyyy_xyyz, g_xx_0_xyyyy_xyzz, g_xx_0_xyyyy_xzzz, g_xx_0_xyyyy_yyyy, g_xx_0_xyyyy_yyyz, g_xx_0_xyyyy_yyzz, g_xx_0_xyyyy_yzzz, g_xx_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyy_xxxx[k] = -g_xx_0_xyyy_xxxx[k] * ab_y + g_xx_0_xyyy_xxxxy[k];

                g_xx_0_xyyyy_xxxy[k] = -g_xx_0_xyyy_xxxy[k] * ab_y + g_xx_0_xyyy_xxxyy[k];

                g_xx_0_xyyyy_xxxz[k] = -g_xx_0_xyyy_xxxz[k] * ab_y + g_xx_0_xyyy_xxxyz[k];

                g_xx_0_xyyyy_xxyy[k] = -g_xx_0_xyyy_xxyy[k] * ab_y + g_xx_0_xyyy_xxyyy[k];

                g_xx_0_xyyyy_xxyz[k] = -g_xx_0_xyyy_xxyz[k] * ab_y + g_xx_0_xyyy_xxyyz[k];

                g_xx_0_xyyyy_xxzz[k] = -g_xx_0_xyyy_xxzz[k] * ab_y + g_xx_0_xyyy_xxyzz[k];

                g_xx_0_xyyyy_xyyy[k] = -g_xx_0_xyyy_xyyy[k] * ab_y + g_xx_0_xyyy_xyyyy[k];

                g_xx_0_xyyyy_xyyz[k] = -g_xx_0_xyyy_xyyz[k] * ab_y + g_xx_0_xyyy_xyyyz[k];

                g_xx_0_xyyyy_xyzz[k] = -g_xx_0_xyyy_xyzz[k] * ab_y + g_xx_0_xyyy_xyyzz[k];

                g_xx_0_xyyyy_xzzz[k] = -g_xx_0_xyyy_xzzz[k] * ab_y + g_xx_0_xyyy_xyzzz[k];

                g_xx_0_xyyyy_yyyy[k] = -g_xx_0_xyyy_yyyy[k] * ab_y + g_xx_0_xyyy_yyyyy[k];

                g_xx_0_xyyyy_yyyz[k] = -g_xx_0_xyyy_yyyz[k] * ab_y + g_xx_0_xyyy_yyyyz[k];

                g_xx_0_xyyyy_yyzz[k] = -g_xx_0_xyyy_yyzz[k] * ab_y + g_xx_0_xyyy_yyyzz[k];

                g_xx_0_xyyyy_yzzz[k] = -g_xx_0_xyyy_yzzz[k] * ab_y + g_xx_0_xyyy_yyzzz[k];

                g_xx_0_xyyyy_zzzz[k] = -g_xx_0_xyyy_zzzz[k] * ab_y + g_xx_0_xyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xx_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyz_xxxx, g_xx_0_xyyyz_xxxy, g_xx_0_xyyyz_xxxz, g_xx_0_xyyyz_xxyy, g_xx_0_xyyyz_xxyz, g_xx_0_xyyyz_xxzz, g_xx_0_xyyyz_xyyy, g_xx_0_xyyyz_xyyz, g_xx_0_xyyyz_xyzz, g_xx_0_xyyyz_xzzz, g_xx_0_xyyyz_yyyy, g_xx_0_xyyyz_yyyz, g_xx_0_xyyyz_yyzz, g_xx_0_xyyyz_yzzz, g_xx_0_xyyyz_zzzz, g_xx_0_xyyz_xxxx, g_xx_0_xyyz_xxxxy, g_xx_0_xyyz_xxxy, g_xx_0_xyyz_xxxyy, g_xx_0_xyyz_xxxyz, g_xx_0_xyyz_xxxz, g_xx_0_xyyz_xxyy, g_xx_0_xyyz_xxyyy, g_xx_0_xyyz_xxyyz, g_xx_0_xyyz_xxyz, g_xx_0_xyyz_xxyzz, g_xx_0_xyyz_xxzz, g_xx_0_xyyz_xyyy, g_xx_0_xyyz_xyyyy, g_xx_0_xyyz_xyyyz, g_xx_0_xyyz_xyyz, g_xx_0_xyyz_xyyzz, g_xx_0_xyyz_xyzz, g_xx_0_xyyz_xyzzz, g_xx_0_xyyz_xzzz, g_xx_0_xyyz_yyyy, g_xx_0_xyyz_yyyyy, g_xx_0_xyyz_yyyyz, g_xx_0_xyyz_yyyz, g_xx_0_xyyz_yyyzz, g_xx_0_xyyz_yyzz, g_xx_0_xyyz_yyzzz, g_xx_0_xyyz_yzzz, g_xx_0_xyyz_yzzzz, g_xx_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyz_xxxx[k] = -g_xx_0_xyyz_xxxx[k] * ab_y + g_xx_0_xyyz_xxxxy[k];

                g_xx_0_xyyyz_xxxy[k] = -g_xx_0_xyyz_xxxy[k] * ab_y + g_xx_0_xyyz_xxxyy[k];

                g_xx_0_xyyyz_xxxz[k] = -g_xx_0_xyyz_xxxz[k] * ab_y + g_xx_0_xyyz_xxxyz[k];

                g_xx_0_xyyyz_xxyy[k] = -g_xx_0_xyyz_xxyy[k] * ab_y + g_xx_0_xyyz_xxyyy[k];

                g_xx_0_xyyyz_xxyz[k] = -g_xx_0_xyyz_xxyz[k] * ab_y + g_xx_0_xyyz_xxyyz[k];

                g_xx_0_xyyyz_xxzz[k] = -g_xx_0_xyyz_xxzz[k] * ab_y + g_xx_0_xyyz_xxyzz[k];

                g_xx_0_xyyyz_xyyy[k] = -g_xx_0_xyyz_xyyy[k] * ab_y + g_xx_0_xyyz_xyyyy[k];

                g_xx_0_xyyyz_xyyz[k] = -g_xx_0_xyyz_xyyz[k] * ab_y + g_xx_0_xyyz_xyyyz[k];

                g_xx_0_xyyyz_xyzz[k] = -g_xx_0_xyyz_xyzz[k] * ab_y + g_xx_0_xyyz_xyyzz[k];

                g_xx_0_xyyyz_xzzz[k] = -g_xx_0_xyyz_xzzz[k] * ab_y + g_xx_0_xyyz_xyzzz[k];

                g_xx_0_xyyyz_yyyy[k] = -g_xx_0_xyyz_yyyy[k] * ab_y + g_xx_0_xyyz_yyyyy[k];

                g_xx_0_xyyyz_yyyz[k] = -g_xx_0_xyyz_yyyz[k] * ab_y + g_xx_0_xyyz_yyyyz[k];

                g_xx_0_xyyyz_yyzz[k] = -g_xx_0_xyyz_yyzz[k] * ab_y + g_xx_0_xyyz_yyyzz[k];

                g_xx_0_xyyyz_yzzz[k] = -g_xx_0_xyyz_yzzz[k] * ab_y + g_xx_0_xyyz_yyzzz[k];

                g_xx_0_xyyyz_zzzz[k] = -g_xx_0_xyyz_zzzz[k] * ab_y + g_xx_0_xyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xx_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzz_xxxx, g_xx_0_xyyzz_xxxy, g_xx_0_xyyzz_xxxz, g_xx_0_xyyzz_xxyy, g_xx_0_xyyzz_xxyz, g_xx_0_xyyzz_xxzz, g_xx_0_xyyzz_xyyy, g_xx_0_xyyzz_xyyz, g_xx_0_xyyzz_xyzz, g_xx_0_xyyzz_xzzz, g_xx_0_xyyzz_yyyy, g_xx_0_xyyzz_yyyz, g_xx_0_xyyzz_yyzz, g_xx_0_xyyzz_yzzz, g_xx_0_xyyzz_zzzz, g_xx_0_xyzz_xxxx, g_xx_0_xyzz_xxxxy, g_xx_0_xyzz_xxxy, g_xx_0_xyzz_xxxyy, g_xx_0_xyzz_xxxyz, g_xx_0_xyzz_xxxz, g_xx_0_xyzz_xxyy, g_xx_0_xyzz_xxyyy, g_xx_0_xyzz_xxyyz, g_xx_0_xyzz_xxyz, g_xx_0_xyzz_xxyzz, g_xx_0_xyzz_xxzz, g_xx_0_xyzz_xyyy, g_xx_0_xyzz_xyyyy, g_xx_0_xyzz_xyyyz, g_xx_0_xyzz_xyyz, g_xx_0_xyzz_xyyzz, g_xx_0_xyzz_xyzz, g_xx_0_xyzz_xyzzz, g_xx_0_xyzz_xzzz, g_xx_0_xyzz_yyyy, g_xx_0_xyzz_yyyyy, g_xx_0_xyzz_yyyyz, g_xx_0_xyzz_yyyz, g_xx_0_xyzz_yyyzz, g_xx_0_xyzz_yyzz, g_xx_0_xyzz_yyzzz, g_xx_0_xyzz_yzzz, g_xx_0_xyzz_yzzzz, g_xx_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzz_xxxx[k] = -g_xx_0_xyzz_xxxx[k] * ab_y + g_xx_0_xyzz_xxxxy[k];

                g_xx_0_xyyzz_xxxy[k] = -g_xx_0_xyzz_xxxy[k] * ab_y + g_xx_0_xyzz_xxxyy[k];

                g_xx_0_xyyzz_xxxz[k] = -g_xx_0_xyzz_xxxz[k] * ab_y + g_xx_0_xyzz_xxxyz[k];

                g_xx_0_xyyzz_xxyy[k] = -g_xx_0_xyzz_xxyy[k] * ab_y + g_xx_0_xyzz_xxyyy[k];

                g_xx_0_xyyzz_xxyz[k] = -g_xx_0_xyzz_xxyz[k] * ab_y + g_xx_0_xyzz_xxyyz[k];

                g_xx_0_xyyzz_xxzz[k] = -g_xx_0_xyzz_xxzz[k] * ab_y + g_xx_0_xyzz_xxyzz[k];

                g_xx_0_xyyzz_xyyy[k] = -g_xx_0_xyzz_xyyy[k] * ab_y + g_xx_0_xyzz_xyyyy[k];

                g_xx_0_xyyzz_xyyz[k] = -g_xx_0_xyzz_xyyz[k] * ab_y + g_xx_0_xyzz_xyyyz[k];

                g_xx_0_xyyzz_xyzz[k] = -g_xx_0_xyzz_xyzz[k] * ab_y + g_xx_0_xyzz_xyyzz[k];

                g_xx_0_xyyzz_xzzz[k] = -g_xx_0_xyzz_xzzz[k] * ab_y + g_xx_0_xyzz_xyzzz[k];

                g_xx_0_xyyzz_yyyy[k] = -g_xx_0_xyzz_yyyy[k] * ab_y + g_xx_0_xyzz_yyyyy[k];

                g_xx_0_xyyzz_yyyz[k] = -g_xx_0_xyzz_yyyz[k] * ab_y + g_xx_0_xyzz_yyyyz[k];

                g_xx_0_xyyzz_yyzz[k] = -g_xx_0_xyzz_yyzz[k] * ab_y + g_xx_0_xyzz_yyyzz[k];

                g_xx_0_xyyzz_yzzz[k] = -g_xx_0_xyzz_yzzz[k] * ab_y + g_xx_0_xyzz_yyzzz[k];

                g_xx_0_xyyzz_zzzz[k] = -g_xx_0_xyzz_zzzz[k] * ab_y + g_xx_0_xyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xx_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzz_xxxx, g_xx_0_xyzzz_xxxy, g_xx_0_xyzzz_xxxz, g_xx_0_xyzzz_xxyy, g_xx_0_xyzzz_xxyz, g_xx_0_xyzzz_xxzz, g_xx_0_xyzzz_xyyy, g_xx_0_xyzzz_xyyz, g_xx_0_xyzzz_xyzz, g_xx_0_xyzzz_xzzz, g_xx_0_xyzzz_yyyy, g_xx_0_xyzzz_yyyz, g_xx_0_xyzzz_yyzz, g_xx_0_xyzzz_yzzz, g_xx_0_xyzzz_zzzz, g_xx_0_xzzz_xxxx, g_xx_0_xzzz_xxxxy, g_xx_0_xzzz_xxxy, g_xx_0_xzzz_xxxyy, g_xx_0_xzzz_xxxyz, g_xx_0_xzzz_xxxz, g_xx_0_xzzz_xxyy, g_xx_0_xzzz_xxyyy, g_xx_0_xzzz_xxyyz, g_xx_0_xzzz_xxyz, g_xx_0_xzzz_xxyzz, g_xx_0_xzzz_xxzz, g_xx_0_xzzz_xyyy, g_xx_0_xzzz_xyyyy, g_xx_0_xzzz_xyyyz, g_xx_0_xzzz_xyyz, g_xx_0_xzzz_xyyzz, g_xx_0_xzzz_xyzz, g_xx_0_xzzz_xyzzz, g_xx_0_xzzz_xzzz, g_xx_0_xzzz_yyyy, g_xx_0_xzzz_yyyyy, g_xx_0_xzzz_yyyyz, g_xx_0_xzzz_yyyz, g_xx_0_xzzz_yyyzz, g_xx_0_xzzz_yyzz, g_xx_0_xzzz_yyzzz, g_xx_0_xzzz_yzzz, g_xx_0_xzzz_yzzzz, g_xx_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzz_xxxx[k] = -g_xx_0_xzzz_xxxx[k] * ab_y + g_xx_0_xzzz_xxxxy[k];

                g_xx_0_xyzzz_xxxy[k] = -g_xx_0_xzzz_xxxy[k] * ab_y + g_xx_0_xzzz_xxxyy[k];

                g_xx_0_xyzzz_xxxz[k] = -g_xx_0_xzzz_xxxz[k] * ab_y + g_xx_0_xzzz_xxxyz[k];

                g_xx_0_xyzzz_xxyy[k] = -g_xx_0_xzzz_xxyy[k] * ab_y + g_xx_0_xzzz_xxyyy[k];

                g_xx_0_xyzzz_xxyz[k] = -g_xx_0_xzzz_xxyz[k] * ab_y + g_xx_0_xzzz_xxyyz[k];

                g_xx_0_xyzzz_xxzz[k] = -g_xx_0_xzzz_xxzz[k] * ab_y + g_xx_0_xzzz_xxyzz[k];

                g_xx_0_xyzzz_xyyy[k] = -g_xx_0_xzzz_xyyy[k] * ab_y + g_xx_0_xzzz_xyyyy[k];

                g_xx_0_xyzzz_xyyz[k] = -g_xx_0_xzzz_xyyz[k] * ab_y + g_xx_0_xzzz_xyyyz[k];

                g_xx_0_xyzzz_xyzz[k] = -g_xx_0_xzzz_xyzz[k] * ab_y + g_xx_0_xzzz_xyyzz[k];

                g_xx_0_xyzzz_xzzz[k] = -g_xx_0_xzzz_xzzz[k] * ab_y + g_xx_0_xzzz_xyzzz[k];

                g_xx_0_xyzzz_yyyy[k] = -g_xx_0_xzzz_yyyy[k] * ab_y + g_xx_0_xzzz_yyyyy[k];

                g_xx_0_xyzzz_yyyz[k] = -g_xx_0_xzzz_yyyz[k] * ab_y + g_xx_0_xzzz_yyyyz[k];

                g_xx_0_xyzzz_yyzz[k] = -g_xx_0_xzzz_yyzz[k] * ab_y + g_xx_0_xzzz_yyyzz[k];

                g_xx_0_xyzzz_yzzz[k] = -g_xx_0_xzzz_yzzz[k] * ab_y + g_xx_0_xzzz_yyzzz[k];

                g_xx_0_xyzzz_zzzz[k] = -g_xx_0_xzzz_zzzz[k] * ab_y + g_xx_0_xzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xx_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzz_xxxx, g_xx_0_xzzz_xxxxz, g_xx_0_xzzz_xxxy, g_xx_0_xzzz_xxxyz, g_xx_0_xzzz_xxxz, g_xx_0_xzzz_xxxzz, g_xx_0_xzzz_xxyy, g_xx_0_xzzz_xxyyz, g_xx_0_xzzz_xxyz, g_xx_0_xzzz_xxyzz, g_xx_0_xzzz_xxzz, g_xx_0_xzzz_xxzzz, g_xx_0_xzzz_xyyy, g_xx_0_xzzz_xyyyz, g_xx_0_xzzz_xyyz, g_xx_0_xzzz_xyyzz, g_xx_0_xzzz_xyzz, g_xx_0_xzzz_xyzzz, g_xx_0_xzzz_xzzz, g_xx_0_xzzz_xzzzz, g_xx_0_xzzz_yyyy, g_xx_0_xzzz_yyyyz, g_xx_0_xzzz_yyyz, g_xx_0_xzzz_yyyzz, g_xx_0_xzzz_yyzz, g_xx_0_xzzz_yyzzz, g_xx_0_xzzz_yzzz, g_xx_0_xzzz_yzzzz, g_xx_0_xzzz_zzzz, g_xx_0_xzzz_zzzzz, g_xx_0_xzzzz_xxxx, g_xx_0_xzzzz_xxxy, g_xx_0_xzzzz_xxxz, g_xx_0_xzzzz_xxyy, g_xx_0_xzzzz_xxyz, g_xx_0_xzzzz_xxzz, g_xx_0_xzzzz_xyyy, g_xx_0_xzzzz_xyyz, g_xx_0_xzzzz_xyzz, g_xx_0_xzzzz_xzzz, g_xx_0_xzzzz_yyyy, g_xx_0_xzzzz_yyyz, g_xx_0_xzzzz_yyzz, g_xx_0_xzzzz_yzzz, g_xx_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzz_xxxx[k] = -g_xx_0_xzzz_xxxx[k] * ab_z + g_xx_0_xzzz_xxxxz[k];

                g_xx_0_xzzzz_xxxy[k] = -g_xx_0_xzzz_xxxy[k] * ab_z + g_xx_0_xzzz_xxxyz[k];

                g_xx_0_xzzzz_xxxz[k] = -g_xx_0_xzzz_xxxz[k] * ab_z + g_xx_0_xzzz_xxxzz[k];

                g_xx_0_xzzzz_xxyy[k] = -g_xx_0_xzzz_xxyy[k] * ab_z + g_xx_0_xzzz_xxyyz[k];

                g_xx_0_xzzzz_xxyz[k] = -g_xx_0_xzzz_xxyz[k] * ab_z + g_xx_0_xzzz_xxyzz[k];

                g_xx_0_xzzzz_xxzz[k] = -g_xx_0_xzzz_xxzz[k] * ab_z + g_xx_0_xzzz_xxzzz[k];

                g_xx_0_xzzzz_xyyy[k] = -g_xx_0_xzzz_xyyy[k] * ab_z + g_xx_0_xzzz_xyyyz[k];

                g_xx_0_xzzzz_xyyz[k] = -g_xx_0_xzzz_xyyz[k] * ab_z + g_xx_0_xzzz_xyyzz[k];

                g_xx_0_xzzzz_xyzz[k] = -g_xx_0_xzzz_xyzz[k] * ab_z + g_xx_0_xzzz_xyzzz[k];

                g_xx_0_xzzzz_xzzz[k] = -g_xx_0_xzzz_xzzz[k] * ab_z + g_xx_0_xzzz_xzzzz[k];

                g_xx_0_xzzzz_yyyy[k] = -g_xx_0_xzzz_yyyy[k] * ab_z + g_xx_0_xzzz_yyyyz[k];

                g_xx_0_xzzzz_yyyz[k] = -g_xx_0_xzzz_yyyz[k] * ab_z + g_xx_0_xzzz_yyyzz[k];

                g_xx_0_xzzzz_yyzz[k] = -g_xx_0_xzzz_yyzz[k] * ab_z + g_xx_0_xzzz_yyzzz[k];

                g_xx_0_xzzzz_yzzz[k] = -g_xx_0_xzzz_yzzz[k] * ab_z + g_xx_0_xzzz_yzzzz[k];

                g_xx_0_xzzzz_zzzz[k] = -g_xx_0_xzzz_zzzz[k] * ab_z + g_xx_0_xzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xx_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyy_xxxx, g_xx_0_yyyy_xxxxy, g_xx_0_yyyy_xxxy, g_xx_0_yyyy_xxxyy, g_xx_0_yyyy_xxxyz, g_xx_0_yyyy_xxxz, g_xx_0_yyyy_xxyy, g_xx_0_yyyy_xxyyy, g_xx_0_yyyy_xxyyz, g_xx_0_yyyy_xxyz, g_xx_0_yyyy_xxyzz, g_xx_0_yyyy_xxzz, g_xx_0_yyyy_xyyy, g_xx_0_yyyy_xyyyy, g_xx_0_yyyy_xyyyz, g_xx_0_yyyy_xyyz, g_xx_0_yyyy_xyyzz, g_xx_0_yyyy_xyzz, g_xx_0_yyyy_xyzzz, g_xx_0_yyyy_xzzz, g_xx_0_yyyy_yyyy, g_xx_0_yyyy_yyyyy, g_xx_0_yyyy_yyyyz, g_xx_0_yyyy_yyyz, g_xx_0_yyyy_yyyzz, g_xx_0_yyyy_yyzz, g_xx_0_yyyy_yyzzz, g_xx_0_yyyy_yzzz, g_xx_0_yyyy_yzzzz, g_xx_0_yyyy_zzzz, g_xx_0_yyyyy_xxxx, g_xx_0_yyyyy_xxxy, g_xx_0_yyyyy_xxxz, g_xx_0_yyyyy_xxyy, g_xx_0_yyyyy_xxyz, g_xx_0_yyyyy_xxzz, g_xx_0_yyyyy_xyyy, g_xx_0_yyyyy_xyyz, g_xx_0_yyyyy_xyzz, g_xx_0_yyyyy_xzzz, g_xx_0_yyyyy_yyyy, g_xx_0_yyyyy_yyyz, g_xx_0_yyyyy_yyzz, g_xx_0_yyyyy_yzzz, g_xx_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyy_xxxx[k] = -g_xx_0_yyyy_xxxx[k] * ab_y + g_xx_0_yyyy_xxxxy[k];

                g_xx_0_yyyyy_xxxy[k] = -g_xx_0_yyyy_xxxy[k] * ab_y + g_xx_0_yyyy_xxxyy[k];

                g_xx_0_yyyyy_xxxz[k] = -g_xx_0_yyyy_xxxz[k] * ab_y + g_xx_0_yyyy_xxxyz[k];

                g_xx_0_yyyyy_xxyy[k] = -g_xx_0_yyyy_xxyy[k] * ab_y + g_xx_0_yyyy_xxyyy[k];

                g_xx_0_yyyyy_xxyz[k] = -g_xx_0_yyyy_xxyz[k] * ab_y + g_xx_0_yyyy_xxyyz[k];

                g_xx_0_yyyyy_xxzz[k] = -g_xx_0_yyyy_xxzz[k] * ab_y + g_xx_0_yyyy_xxyzz[k];

                g_xx_0_yyyyy_xyyy[k] = -g_xx_0_yyyy_xyyy[k] * ab_y + g_xx_0_yyyy_xyyyy[k];

                g_xx_0_yyyyy_xyyz[k] = -g_xx_0_yyyy_xyyz[k] * ab_y + g_xx_0_yyyy_xyyyz[k];

                g_xx_0_yyyyy_xyzz[k] = -g_xx_0_yyyy_xyzz[k] * ab_y + g_xx_0_yyyy_xyyzz[k];

                g_xx_0_yyyyy_xzzz[k] = -g_xx_0_yyyy_xzzz[k] * ab_y + g_xx_0_yyyy_xyzzz[k];

                g_xx_0_yyyyy_yyyy[k] = -g_xx_0_yyyy_yyyy[k] * ab_y + g_xx_0_yyyy_yyyyy[k];

                g_xx_0_yyyyy_yyyz[k] = -g_xx_0_yyyy_yyyz[k] * ab_y + g_xx_0_yyyy_yyyyz[k];

                g_xx_0_yyyyy_yyzz[k] = -g_xx_0_yyyy_yyzz[k] * ab_y + g_xx_0_yyyy_yyyzz[k];

                g_xx_0_yyyyy_yzzz[k] = -g_xx_0_yyyy_yzzz[k] * ab_y + g_xx_0_yyyy_yyzzz[k];

                g_xx_0_yyyyy_zzzz[k] = -g_xx_0_yyyy_zzzz[k] * ab_y + g_xx_0_yyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xx_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyz_xxxx, g_xx_0_yyyyz_xxxy, g_xx_0_yyyyz_xxxz, g_xx_0_yyyyz_xxyy, g_xx_0_yyyyz_xxyz, g_xx_0_yyyyz_xxzz, g_xx_0_yyyyz_xyyy, g_xx_0_yyyyz_xyyz, g_xx_0_yyyyz_xyzz, g_xx_0_yyyyz_xzzz, g_xx_0_yyyyz_yyyy, g_xx_0_yyyyz_yyyz, g_xx_0_yyyyz_yyzz, g_xx_0_yyyyz_yzzz, g_xx_0_yyyyz_zzzz, g_xx_0_yyyz_xxxx, g_xx_0_yyyz_xxxxy, g_xx_0_yyyz_xxxy, g_xx_0_yyyz_xxxyy, g_xx_0_yyyz_xxxyz, g_xx_0_yyyz_xxxz, g_xx_0_yyyz_xxyy, g_xx_0_yyyz_xxyyy, g_xx_0_yyyz_xxyyz, g_xx_0_yyyz_xxyz, g_xx_0_yyyz_xxyzz, g_xx_0_yyyz_xxzz, g_xx_0_yyyz_xyyy, g_xx_0_yyyz_xyyyy, g_xx_0_yyyz_xyyyz, g_xx_0_yyyz_xyyz, g_xx_0_yyyz_xyyzz, g_xx_0_yyyz_xyzz, g_xx_0_yyyz_xyzzz, g_xx_0_yyyz_xzzz, g_xx_0_yyyz_yyyy, g_xx_0_yyyz_yyyyy, g_xx_0_yyyz_yyyyz, g_xx_0_yyyz_yyyz, g_xx_0_yyyz_yyyzz, g_xx_0_yyyz_yyzz, g_xx_0_yyyz_yyzzz, g_xx_0_yyyz_yzzz, g_xx_0_yyyz_yzzzz, g_xx_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyz_xxxx[k] = -g_xx_0_yyyz_xxxx[k] * ab_y + g_xx_0_yyyz_xxxxy[k];

                g_xx_0_yyyyz_xxxy[k] = -g_xx_0_yyyz_xxxy[k] * ab_y + g_xx_0_yyyz_xxxyy[k];

                g_xx_0_yyyyz_xxxz[k] = -g_xx_0_yyyz_xxxz[k] * ab_y + g_xx_0_yyyz_xxxyz[k];

                g_xx_0_yyyyz_xxyy[k] = -g_xx_0_yyyz_xxyy[k] * ab_y + g_xx_0_yyyz_xxyyy[k];

                g_xx_0_yyyyz_xxyz[k] = -g_xx_0_yyyz_xxyz[k] * ab_y + g_xx_0_yyyz_xxyyz[k];

                g_xx_0_yyyyz_xxzz[k] = -g_xx_0_yyyz_xxzz[k] * ab_y + g_xx_0_yyyz_xxyzz[k];

                g_xx_0_yyyyz_xyyy[k] = -g_xx_0_yyyz_xyyy[k] * ab_y + g_xx_0_yyyz_xyyyy[k];

                g_xx_0_yyyyz_xyyz[k] = -g_xx_0_yyyz_xyyz[k] * ab_y + g_xx_0_yyyz_xyyyz[k];

                g_xx_0_yyyyz_xyzz[k] = -g_xx_0_yyyz_xyzz[k] * ab_y + g_xx_0_yyyz_xyyzz[k];

                g_xx_0_yyyyz_xzzz[k] = -g_xx_0_yyyz_xzzz[k] * ab_y + g_xx_0_yyyz_xyzzz[k];

                g_xx_0_yyyyz_yyyy[k] = -g_xx_0_yyyz_yyyy[k] * ab_y + g_xx_0_yyyz_yyyyy[k];

                g_xx_0_yyyyz_yyyz[k] = -g_xx_0_yyyz_yyyz[k] * ab_y + g_xx_0_yyyz_yyyyz[k];

                g_xx_0_yyyyz_yyzz[k] = -g_xx_0_yyyz_yyzz[k] * ab_y + g_xx_0_yyyz_yyyzz[k];

                g_xx_0_yyyyz_yzzz[k] = -g_xx_0_yyyz_yzzz[k] * ab_y + g_xx_0_yyyz_yyzzz[k];

                g_xx_0_yyyyz_zzzz[k] = -g_xx_0_yyyz_zzzz[k] * ab_y + g_xx_0_yyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xx_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzz_xxxx, g_xx_0_yyyzz_xxxy, g_xx_0_yyyzz_xxxz, g_xx_0_yyyzz_xxyy, g_xx_0_yyyzz_xxyz, g_xx_0_yyyzz_xxzz, g_xx_0_yyyzz_xyyy, g_xx_0_yyyzz_xyyz, g_xx_0_yyyzz_xyzz, g_xx_0_yyyzz_xzzz, g_xx_0_yyyzz_yyyy, g_xx_0_yyyzz_yyyz, g_xx_0_yyyzz_yyzz, g_xx_0_yyyzz_yzzz, g_xx_0_yyyzz_zzzz, g_xx_0_yyzz_xxxx, g_xx_0_yyzz_xxxxy, g_xx_0_yyzz_xxxy, g_xx_0_yyzz_xxxyy, g_xx_0_yyzz_xxxyz, g_xx_0_yyzz_xxxz, g_xx_0_yyzz_xxyy, g_xx_0_yyzz_xxyyy, g_xx_0_yyzz_xxyyz, g_xx_0_yyzz_xxyz, g_xx_0_yyzz_xxyzz, g_xx_0_yyzz_xxzz, g_xx_0_yyzz_xyyy, g_xx_0_yyzz_xyyyy, g_xx_0_yyzz_xyyyz, g_xx_0_yyzz_xyyz, g_xx_0_yyzz_xyyzz, g_xx_0_yyzz_xyzz, g_xx_0_yyzz_xyzzz, g_xx_0_yyzz_xzzz, g_xx_0_yyzz_yyyy, g_xx_0_yyzz_yyyyy, g_xx_0_yyzz_yyyyz, g_xx_0_yyzz_yyyz, g_xx_0_yyzz_yyyzz, g_xx_0_yyzz_yyzz, g_xx_0_yyzz_yyzzz, g_xx_0_yyzz_yzzz, g_xx_0_yyzz_yzzzz, g_xx_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzz_xxxx[k] = -g_xx_0_yyzz_xxxx[k] * ab_y + g_xx_0_yyzz_xxxxy[k];

                g_xx_0_yyyzz_xxxy[k] = -g_xx_0_yyzz_xxxy[k] * ab_y + g_xx_0_yyzz_xxxyy[k];

                g_xx_0_yyyzz_xxxz[k] = -g_xx_0_yyzz_xxxz[k] * ab_y + g_xx_0_yyzz_xxxyz[k];

                g_xx_0_yyyzz_xxyy[k] = -g_xx_0_yyzz_xxyy[k] * ab_y + g_xx_0_yyzz_xxyyy[k];

                g_xx_0_yyyzz_xxyz[k] = -g_xx_0_yyzz_xxyz[k] * ab_y + g_xx_0_yyzz_xxyyz[k];

                g_xx_0_yyyzz_xxzz[k] = -g_xx_0_yyzz_xxzz[k] * ab_y + g_xx_0_yyzz_xxyzz[k];

                g_xx_0_yyyzz_xyyy[k] = -g_xx_0_yyzz_xyyy[k] * ab_y + g_xx_0_yyzz_xyyyy[k];

                g_xx_0_yyyzz_xyyz[k] = -g_xx_0_yyzz_xyyz[k] * ab_y + g_xx_0_yyzz_xyyyz[k];

                g_xx_0_yyyzz_xyzz[k] = -g_xx_0_yyzz_xyzz[k] * ab_y + g_xx_0_yyzz_xyyzz[k];

                g_xx_0_yyyzz_xzzz[k] = -g_xx_0_yyzz_xzzz[k] * ab_y + g_xx_0_yyzz_xyzzz[k];

                g_xx_0_yyyzz_yyyy[k] = -g_xx_0_yyzz_yyyy[k] * ab_y + g_xx_0_yyzz_yyyyy[k];

                g_xx_0_yyyzz_yyyz[k] = -g_xx_0_yyzz_yyyz[k] * ab_y + g_xx_0_yyzz_yyyyz[k];

                g_xx_0_yyyzz_yyzz[k] = -g_xx_0_yyzz_yyzz[k] * ab_y + g_xx_0_yyzz_yyyzz[k];

                g_xx_0_yyyzz_yzzz[k] = -g_xx_0_yyzz_yzzz[k] * ab_y + g_xx_0_yyzz_yyzzz[k];

                g_xx_0_yyyzz_zzzz[k] = -g_xx_0_yyzz_zzzz[k] * ab_y + g_xx_0_yyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 270 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 271 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 272 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 273 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 274 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 275 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 276 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 277 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 278 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 279 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 280 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 281 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 282 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 283 * ccomps * dcomps);

            auto g_xx_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzz_xxxx, g_xx_0_yyzzz_xxxy, g_xx_0_yyzzz_xxxz, g_xx_0_yyzzz_xxyy, g_xx_0_yyzzz_xxyz, g_xx_0_yyzzz_xxzz, g_xx_0_yyzzz_xyyy, g_xx_0_yyzzz_xyyz, g_xx_0_yyzzz_xyzz, g_xx_0_yyzzz_xzzz, g_xx_0_yyzzz_yyyy, g_xx_0_yyzzz_yyyz, g_xx_0_yyzzz_yyzz, g_xx_0_yyzzz_yzzz, g_xx_0_yyzzz_zzzz, g_xx_0_yzzz_xxxx, g_xx_0_yzzz_xxxxy, g_xx_0_yzzz_xxxy, g_xx_0_yzzz_xxxyy, g_xx_0_yzzz_xxxyz, g_xx_0_yzzz_xxxz, g_xx_0_yzzz_xxyy, g_xx_0_yzzz_xxyyy, g_xx_0_yzzz_xxyyz, g_xx_0_yzzz_xxyz, g_xx_0_yzzz_xxyzz, g_xx_0_yzzz_xxzz, g_xx_0_yzzz_xyyy, g_xx_0_yzzz_xyyyy, g_xx_0_yzzz_xyyyz, g_xx_0_yzzz_xyyz, g_xx_0_yzzz_xyyzz, g_xx_0_yzzz_xyzz, g_xx_0_yzzz_xyzzz, g_xx_0_yzzz_xzzz, g_xx_0_yzzz_yyyy, g_xx_0_yzzz_yyyyy, g_xx_0_yzzz_yyyyz, g_xx_0_yzzz_yyyz, g_xx_0_yzzz_yyyzz, g_xx_0_yzzz_yyzz, g_xx_0_yzzz_yyzzz, g_xx_0_yzzz_yzzz, g_xx_0_yzzz_yzzzz, g_xx_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzz_xxxx[k] = -g_xx_0_yzzz_xxxx[k] * ab_y + g_xx_0_yzzz_xxxxy[k];

                g_xx_0_yyzzz_xxxy[k] = -g_xx_0_yzzz_xxxy[k] * ab_y + g_xx_0_yzzz_xxxyy[k];

                g_xx_0_yyzzz_xxxz[k] = -g_xx_0_yzzz_xxxz[k] * ab_y + g_xx_0_yzzz_xxxyz[k];

                g_xx_0_yyzzz_xxyy[k] = -g_xx_0_yzzz_xxyy[k] * ab_y + g_xx_0_yzzz_xxyyy[k];

                g_xx_0_yyzzz_xxyz[k] = -g_xx_0_yzzz_xxyz[k] * ab_y + g_xx_0_yzzz_xxyyz[k];

                g_xx_0_yyzzz_xxzz[k] = -g_xx_0_yzzz_xxzz[k] * ab_y + g_xx_0_yzzz_xxyzz[k];

                g_xx_0_yyzzz_xyyy[k] = -g_xx_0_yzzz_xyyy[k] * ab_y + g_xx_0_yzzz_xyyyy[k];

                g_xx_0_yyzzz_xyyz[k] = -g_xx_0_yzzz_xyyz[k] * ab_y + g_xx_0_yzzz_xyyyz[k];

                g_xx_0_yyzzz_xyzz[k] = -g_xx_0_yzzz_xyzz[k] * ab_y + g_xx_0_yzzz_xyyzz[k];

                g_xx_0_yyzzz_xzzz[k] = -g_xx_0_yzzz_xzzz[k] * ab_y + g_xx_0_yzzz_xyzzz[k];

                g_xx_0_yyzzz_yyyy[k] = -g_xx_0_yzzz_yyyy[k] * ab_y + g_xx_0_yzzz_yyyyy[k];

                g_xx_0_yyzzz_yyyz[k] = -g_xx_0_yzzz_yyyz[k] * ab_y + g_xx_0_yzzz_yyyyz[k];

                g_xx_0_yyzzz_yyzz[k] = -g_xx_0_yzzz_yyzz[k] * ab_y + g_xx_0_yzzz_yyyzz[k];

                g_xx_0_yyzzz_yzzz[k] = -g_xx_0_yzzz_yzzz[k] * ab_y + g_xx_0_yzzz_yyzzz[k];

                g_xx_0_yyzzz_zzzz[k] = -g_xx_0_yzzz_zzzz[k] * ab_y + g_xx_0_yzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 285 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 286 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 287 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 288 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 289 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 290 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 291 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 292 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 293 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 294 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 295 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 296 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 297 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 298 * ccomps * dcomps);

            auto g_xx_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzz_xxxx, g_xx_0_yzzzz_xxxy, g_xx_0_yzzzz_xxxz, g_xx_0_yzzzz_xxyy, g_xx_0_yzzzz_xxyz, g_xx_0_yzzzz_xxzz, g_xx_0_yzzzz_xyyy, g_xx_0_yzzzz_xyyz, g_xx_0_yzzzz_xyzz, g_xx_0_yzzzz_xzzz, g_xx_0_yzzzz_yyyy, g_xx_0_yzzzz_yyyz, g_xx_0_yzzzz_yyzz, g_xx_0_yzzzz_yzzz, g_xx_0_yzzzz_zzzz, g_xx_0_zzzz_xxxx, g_xx_0_zzzz_xxxxy, g_xx_0_zzzz_xxxy, g_xx_0_zzzz_xxxyy, g_xx_0_zzzz_xxxyz, g_xx_0_zzzz_xxxz, g_xx_0_zzzz_xxyy, g_xx_0_zzzz_xxyyy, g_xx_0_zzzz_xxyyz, g_xx_0_zzzz_xxyz, g_xx_0_zzzz_xxyzz, g_xx_0_zzzz_xxzz, g_xx_0_zzzz_xyyy, g_xx_0_zzzz_xyyyy, g_xx_0_zzzz_xyyyz, g_xx_0_zzzz_xyyz, g_xx_0_zzzz_xyyzz, g_xx_0_zzzz_xyzz, g_xx_0_zzzz_xyzzz, g_xx_0_zzzz_xzzz, g_xx_0_zzzz_yyyy, g_xx_0_zzzz_yyyyy, g_xx_0_zzzz_yyyyz, g_xx_0_zzzz_yyyz, g_xx_0_zzzz_yyyzz, g_xx_0_zzzz_yyzz, g_xx_0_zzzz_yyzzz, g_xx_0_zzzz_yzzz, g_xx_0_zzzz_yzzzz, g_xx_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzz_xxxx[k] = -g_xx_0_zzzz_xxxx[k] * ab_y + g_xx_0_zzzz_xxxxy[k];

                g_xx_0_yzzzz_xxxy[k] = -g_xx_0_zzzz_xxxy[k] * ab_y + g_xx_0_zzzz_xxxyy[k];

                g_xx_0_yzzzz_xxxz[k] = -g_xx_0_zzzz_xxxz[k] * ab_y + g_xx_0_zzzz_xxxyz[k];

                g_xx_0_yzzzz_xxyy[k] = -g_xx_0_zzzz_xxyy[k] * ab_y + g_xx_0_zzzz_xxyyy[k];

                g_xx_0_yzzzz_xxyz[k] = -g_xx_0_zzzz_xxyz[k] * ab_y + g_xx_0_zzzz_xxyyz[k];

                g_xx_0_yzzzz_xxzz[k] = -g_xx_0_zzzz_xxzz[k] * ab_y + g_xx_0_zzzz_xxyzz[k];

                g_xx_0_yzzzz_xyyy[k] = -g_xx_0_zzzz_xyyy[k] * ab_y + g_xx_0_zzzz_xyyyy[k];

                g_xx_0_yzzzz_xyyz[k] = -g_xx_0_zzzz_xyyz[k] * ab_y + g_xx_0_zzzz_xyyyz[k];

                g_xx_0_yzzzz_xyzz[k] = -g_xx_0_zzzz_xyzz[k] * ab_y + g_xx_0_zzzz_xyyzz[k];

                g_xx_0_yzzzz_xzzz[k] = -g_xx_0_zzzz_xzzz[k] * ab_y + g_xx_0_zzzz_xyzzz[k];

                g_xx_0_yzzzz_yyyy[k] = -g_xx_0_zzzz_yyyy[k] * ab_y + g_xx_0_zzzz_yyyyy[k];

                g_xx_0_yzzzz_yyyz[k] = -g_xx_0_zzzz_yyyz[k] * ab_y + g_xx_0_zzzz_yyyyz[k];

                g_xx_0_yzzzz_yyzz[k] = -g_xx_0_zzzz_yyzz[k] * ab_y + g_xx_0_zzzz_yyyzz[k];

                g_xx_0_yzzzz_yzzz[k] = -g_xx_0_zzzz_yzzz[k] * ab_y + g_xx_0_zzzz_yyzzz[k];

                g_xx_0_yzzzz_zzzz[k] = -g_xx_0_zzzz_zzzz[k] * ab_y + g_xx_0_zzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 300 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 301 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 302 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 303 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 304 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 305 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 306 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 307 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 308 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 309 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 310 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 311 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 312 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 313 * ccomps * dcomps);

            auto g_xx_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzz_xxxx, g_xx_0_zzzz_xxxxz, g_xx_0_zzzz_xxxy, g_xx_0_zzzz_xxxyz, g_xx_0_zzzz_xxxz, g_xx_0_zzzz_xxxzz, g_xx_0_zzzz_xxyy, g_xx_0_zzzz_xxyyz, g_xx_0_zzzz_xxyz, g_xx_0_zzzz_xxyzz, g_xx_0_zzzz_xxzz, g_xx_0_zzzz_xxzzz, g_xx_0_zzzz_xyyy, g_xx_0_zzzz_xyyyz, g_xx_0_zzzz_xyyz, g_xx_0_zzzz_xyyzz, g_xx_0_zzzz_xyzz, g_xx_0_zzzz_xyzzz, g_xx_0_zzzz_xzzz, g_xx_0_zzzz_xzzzz, g_xx_0_zzzz_yyyy, g_xx_0_zzzz_yyyyz, g_xx_0_zzzz_yyyz, g_xx_0_zzzz_yyyzz, g_xx_0_zzzz_yyzz, g_xx_0_zzzz_yyzzz, g_xx_0_zzzz_yzzz, g_xx_0_zzzz_yzzzz, g_xx_0_zzzz_zzzz, g_xx_0_zzzz_zzzzz, g_xx_0_zzzzz_xxxx, g_xx_0_zzzzz_xxxy, g_xx_0_zzzzz_xxxz, g_xx_0_zzzzz_xxyy, g_xx_0_zzzzz_xxyz, g_xx_0_zzzzz_xxzz, g_xx_0_zzzzz_xyyy, g_xx_0_zzzzz_xyyz, g_xx_0_zzzzz_xyzz, g_xx_0_zzzzz_xzzz, g_xx_0_zzzzz_yyyy, g_xx_0_zzzzz_yyyz, g_xx_0_zzzzz_yyzz, g_xx_0_zzzzz_yzzz, g_xx_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzz_xxxx[k] = -g_xx_0_zzzz_xxxx[k] * ab_z + g_xx_0_zzzz_xxxxz[k];

                g_xx_0_zzzzz_xxxy[k] = -g_xx_0_zzzz_xxxy[k] * ab_z + g_xx_0_zzzz_xxxyz[k];

                g_xx_0_zzzzz_xxxz[k] = -g_xx_0_zzzz_xxxz[k] * ab_z + g_xx_0_zzzz_xxxzz[k];

                g_xx_0_zzzzz_xxyy[k] = -g_xx_0_zzzz_xxyy[k] * ab_z + g_xx_0_zzzz_xxyyz[k];

                g_xx_0_zzzzz_xxyz[k] = -g_xx_0_zzzz_xxyz[k] * ab_z + g_xx_0_zzzz_xxyzz[k];

                g_xx_0_zzzzz_xxzz[k] = -g_xx_0_zzzz_xxzz[k] * ab_z + g_xx_0_zzzz_xxzzz[k];

                g_xx_0_zzzzz_xyyy[k] = -g_xx_0_zzzz_xyyy[k] * ab_z + g_xx_0_zzzz_xyyyz[k];

                g_xx_0_zzzzz_xyyz[k] = -g_xx_0_zzzz_xyyz[k] * ab_z + g_xx_0_zzzz_xyyzz[k];

                g_xx_0_zzzzz_xyzz[k] = -g_xx_0_zzzz_xyzz[k] * ab_z + g_xx_0_zzzz_xyzzz[k];

                g_xx_0_zzzzz_xzzz[k] = -g_xx_0_zzzz_xzzz[k] * ab_z + g_xx_0_zzzz_xzzzz[k];

                g_xx_0_zzzzz_yyyy[k] = -g_xx_0_zzzz_yyyy[k] * ab_z + g_xx_0_zzzz_yyyyz[k];

                g_xx_0_zzzzz_yyyz[k] = -g_xx_0_zzzz_yyyz[k] * ab_z + g_xx_0_zzzz_yyyzz[k];

                g_xx_0_zzzzz_yyzz[k] = -g_xx_0_zzzz_yyzz[k] * ab_z + g_xx_0_zzzz_yyzzz[k];

                g_xx_0_zzzzz_yzzz[k] = -g_xx_0_zzzz_yzzz[k] * ab_z + g_xx_0_zzzz_yzzzz[k];

                g_xx_0_zzzzz_zzzz[k] = -g_xx_0_zzzz_zzzz[k] * ab_z + g_xx_0_zzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 317 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 323 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_xxxx, g_xy_0_xxxx_xxxxx, g_xy_0_xxxx_xxxxy, g_xy_0_xxxx_xxxxz, g_xy_0_xxxx_xxxy, g_xy_0_xxxx_xxxyy, g_xy_0_xxxx_xxxyz, g_xy_0_xxxx_xxxz, g_xy_0_xxxx_xxxzz, g_xy_0_xxxx_xxyy, g_xy_0_xxxx_xxyyy, g_xy_0_xxxx_xxyyz, g_xy_0_xxxx_xxyz, g_xy_0_xxxx_xxyzz, g_xy_0_xxxx_xxzz, g_xy_0_xxxx_xxzzz, g_xy_0_xxxx_xyyy, g_xy_0_xxxx_xyyyy, g_xy_0_xxxx_xyyyz, g_xy_0_xxxx_xyyz, g_xy_0_xxxx_xyyzz, g_xy_0_xxxx_xyzz, g_xy_0_xxxx_xyzzz, g_xy_0_xxxx_xzzz, g_xy_0_xxxx_xzzzz, g_xy_0_xxxx_yyyy, g_xy_0_xxxx_yyyz, g_xy_0_xxxx_yyzz, g_xy_0_xxxx_yzzz, g_xy_0_xxxx_zzzz, g_xy_0_xxxxx_xxxx, g_xy_0_xxxxx_xxxy, g_xy_0_xxxxx_xxxz, g_xy_0_xxxxx_xxyy, g_xy_0_xxxxx_xxyz, g_xy_0_xxxxx_xxzz, g_xy_0_xxxxx_xyyy, g_xy_0_xxxxx_xyyz, g_xy_0_xxxxx_xyzz, g_xy_0_xxxxx_xzzz, g_xy_0_xxxxx_yyyy, g_xy_0_xxxxx_yyyz, g_xy_0_xxxxx_yyzz, g_xy_0_xxxxx_yzzz, g_xy_0_xxxxx_zzzz, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_yyyy, g_y_0_xxxx_yyyz, g_y_0_xxxx_yyzz, g_y_0_xxxx_yzzz, g_y_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxx_xxxx[k] = -g_y_0_xxxx_xxxx[k] - g_xy_0_xxxx_xxxx[k] * ab_x + g_xy_0_xxxx_xxxxx[k];

                g_xy_0_xxxxx_xxxy[k] = -g_y_0_xxxx_xxxy[k] - g_xy_0_xxxx_xxxy[k] * ab_x + g_xy_0_xxxx_xxxxy[k];

                g_xy_0_xxxxx_xxxz[k] = -g_y_0_xxxx_xxxz[k] - g_xy_0_xxxx_xxxz[k] * ab_x + g_xy_0_xxxx_xxxxz[k];

                g_xy_0_xxxxx_xxyy[k] = -g_y_0_xxxx_xxyy[k] - g_xy_0_xxxx_xxyy[k] * ab_x + g_xy_0_xxxx_xxxyy[k];

                g_xy_0_xxxxx_xxyz[k] = -g_y_0_xxxx_xxyz[k] - g_xy_0_xxxx_xxyz[k] * ab_x + g_xy_0_xxxx_xxxyz[k];

                g_xy_0_xxxxx_xxzz[k] = -g_y_0_xxxx_xxzz[k] - g_xy_0_xxxx_xxzz[k] * ab_x + g_xy_0_xxxx_xxxzz[k];

                g_xy_0_xxxxx_xyyy[k] = -g_y_0_xxxx_xyyy[k] - g_xy_0_xxxx_xyyy[k] * ab_x + g_xy_0_xxxx_xxyyy[k];

                g_xy_0_xxxxx_xyyz[k] = -g_y_0_xxxx_xyyz[k] - g_xy_0_xxxx_xyyz[k] * ab_x + g_xy_0_xxxx_xxyyz[k];

                g_xy_0_xxxxx_xyzz[k] = -g_y_0_xxxx_xyzz[k] - g_xy_0_xxxx_xyzz[k] * ab_x + g_xy_0_xxxx_xxyzz[k];

                g_xy_0_xxxxx_xzzz[k] = -g_y_0_xxxx_xzzz[k] - g_xy_0_xxxx_xzzz[k] * ab_x + g_xy_0_xxxx_xxzzz[k];

                g_xy_0_xxxxx_yyyy[k] = -g_y_0_xxxx_yyyy[k] - g_xy_0_xxxx_yyyy[k] * ab_x + g_xy_0_xxxx_xyyyy[k];

                g_xy_0_xxxxx_yyyz[k] = -g_y_0_xxxx_yyyz[k] - g_xy_0_xxxx_yyyz[k] * ab_x + g_xy_0_xxxx_xyyyz[k];

                g_xy_0_xxxxx_yyzz[k] = -g_y_0_xxxx_yyzz[k] - g_xy_0_xxxx_yyzz[k] * ab_x + g_xy_0_xxxx_xyyzz[k];

                g_xy_0_xxxxx_yzzz[k] = -g_y_0_xxxx_yzzz[k] - g_xy_0_xxxx_yzzz[k] * ab_x + g_xy_0_xxxx_xyzzz[k];

                g_xy_0_xxxxx_zzzz[k] = -g_y_0_xxxx_zzzz[k] - g_xy_0_xxxx_zzzz[k] * ab_x + g_xy_0_xxxx_xzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 335 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 336 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 337 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 338 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 339 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 340 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 341 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 342 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 343 * ccomps * dcomps);

            auto g_xy_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_xxxx, g_xy_0_xxxxy_xxxy, g_xy_0_xxxxy_xxxz, g_xy_0_xxxxy_xxyy, g_xy_0_xxxxy_xxyz, g_xy_0_xxxxy_xxzz, g_xy_0_xxxxy_xyyy, g_xy_0_xxxxy_xyyz, g_xy_0_xxxxy_xyzz, g_xy_0_xxxxy_xzzz, g_xy_0_xxxxy_yyyy, g_xy_0_xxxxy_yyyz, g_xy_0_xxxxy_yyzz, g_xy_0_xxxxy_yzzz, g_xy_0_xxxxy_zzzz, g_xy_0_xxxy_xxxx, g_xy_0_xxxy_xxxxx, g_xy_0_xxxy_xxxxy, g_xy_0_xxxy_xxxxz, g_xy_0_xxxy_xxxy, g_xy_0_xxxy_xxxyy, g_xy_0_xxxy_xxxyz, g_xy_0_xxxy_xxxz, g_xy_0_xxxy_xxxzz, g_xy_0_xxxy_xxyy, g_xy_0_xxxy_xxyyy, g_xy_0_xxxy_xxyyz, g_xy_0_xxxy_xxyz, g_xy_0_xxxy_xxyzz, g_xy_0_xxxy_xxzz, g_xy_0_xxxy_xxzzz, g_xy_0_xxxy_xyyy, g_xy_0_xxxy_xyyyy, g_xy_0_xxxy_xyyyz, g_xy_0_xxxy_xyyz, g_xy_0_xxxy_xyyzz, g_xy_0_xxxy_xyzz, g_xy_0_xxxy_xyzzz, g_xy_0_xxxy_xzzz, g_xy_0_xxxy_xzzzz, g_xy_0_xxxy_yyyy, g_xy_0_xxxy_yyyz, g_xy_0_xxxy_yyzz, g_xy_0_xxxy_yzzz, g_xy_0_xxxy_zzzz, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_yyyy, g_y_0_xxxy_yyyz, g_y_0_xxxy_yyzz, g_y_0_xxxy_yzzz, g_y_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxy_xxxx[k] = -g_y_0_xxxy_xxxx[k] - g_xy_0_xxxy_xxxx[k] * ab_x + g_xy_0_xxxy_xxxxx[k];

                g_xy_0_xxxxy_xxxy[k] = -g_y_0_xxxy_xxxy[k] - g_xy_0_xxxy_xxxy[k] * ab_x + g_xy_0_xxxy_xxxxy[k];

                g_xy_0_xxxxy_xxxz[k] = -g_y_0_xxxy_xxxz[k] - g_xy_0_xxxy_xxxz[k] * ab_x + g_xy_0_xxxy_xxxxz[k];

                g_xy_0_xxxxy_xxyy[k] = -g_y_0_xxxy_xxyy[k] - g_xy_0_xxxy_xxyy[k] * ab_x + g_xy_0_xxxy_xxxyy[k];

                g_xy_0_xxxxy_xxyz[k] = -g_y_0_xxxy_xxyz[k] - g_xy_0_xxxy_xxyz[k] * ab_x + g_xy_0_xxxy_xxxyz[k];

                g_xy_0_xxxxy_xxzz[k] = -g_y_0_xxxy_xxzz[k] - g_xy_0_xxxy_xxzz[k] * ab_x + g_xy_0_xxxy_xxxzz[k];

                g_xy_0_xxxxy_xyyy[k] = -g_y_0_xxxy_xyyy[k] - g_xy_0_xxxy_xyyy[k] * ab_x + g_xy_0_xxxy_xxyyy[k];

                g_xy_0_xxxxy_xyyz[k] = -g_y_0_xxxy_xyyz[k] - g_xy_0_xxxy_xyyz[k] * ab_x + g_xy_0_xxxy_xxyyz[k];

                g_xy_0_xxxxy_xyzz[k] = -g_y_0_xxxy_xyzz[k] - g_xy_0_xxxy_xyzz[k] * ab_x + g_xy_0_xxxy_xxyzz[k];

                g_xy_0_xxxxy_xzzz[k] = -g_y_0_xxxy_xzzz[k] - g_xy_0_xxxy_xzzz[k] * ab_x + g_xy_0_xxxy_xxzzz[k];

                g_xy_0_xxxxy_yyyy[k] = -g_y_0_xxxy_yyyy[k] - g_xy_0_xxxy_yyyy[k] * ab_x + g_xy_0_xxxy_xyyyy[k];

                g_xy_0_xxxxy_yyyz[k] = -g_y_0_xxxy_yyyz[k] - g_xy_0_xxxy_yyyz[k] * ab_x + g_xy_0_xxxy_xyyyz[k];

                g_xy_0_xxxxy_yyzz[k] = -g_y_0_xxxy_yyzz[k] - g_xy_0_xxxy_yyzz[k] * ab_x + g_xy_0_xxxy_xyyzz[k];

                g_xy_0_xxxxy_yzzz[k] = -g_y_0_xxxy_yzzz[k] - g_xy_0_xxxy_yzzz[k] * ab_x + g_xy_0_xxxy_xyzzz[k];

                g_xy_0_xxxxy_zzzz[k] = -g_y_0_xxxy_zzzz[k] - g_xy_0_xxxy_zzzz[k] * ab_x + g_xy_0_xxxy_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 345 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 346 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 347 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 348 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 349 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 350 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 351 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 352 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 353 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 354 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 355 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 356 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 357 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 358 * ccomps * dcomps);

            auto g_xy_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_xxxx, g_xy_0_xxxx_xxxxz, g_xy_0_xxxx_xxxy, g_xy_0_xxxx_xxxyz, g_xy_0_xxxx_xxxz, g_xy_0_xxxx_xxxzz, g_xy_0_xxxx_xxyy, g_xy_0_xxxx_xxyyz, g_xy_0_xxxx_xxyz, g_xy_0_xxxx_xxyzz, g_xy_0_xxxx_xxzz, g_xy_0_xxxx_xxzzz, g_xy_0_xxxx_xyyy, g_xy_0_xxxx_xyyyz, g_xy_0_xxxx_xyyz, g_xy_0_xxxx_xyyzz, g_xy_0_xxxx_xyzz, g_xy_0_xxxx_xyzzz, g_xy_0_xxxx_xzzz, g_xy_0_xxxx_xzzzz, g_xy_0_xxxx_yyyy, g_xy_0_xxxx_yyyyz, g_xy_0_xxxx_yyyz, g_xy_0_xxxx_yyyzz, g_xy_0_xxxx_yyzz, g_xy_0_xxxx_yyzzz, g_xy_0_xxxx_yzzz, g_xy_0_xxxx_yzzzz, g_xy_0_xxxx_zzzz, g_xy_0_xxxx_zzzzz, g_xy_0_xxxxz_xxxx, g_xy_0_xxxxz_xxxy, g_xy_0_xxxxz_xxxz, g_xy_0_xxxxz_xxyy, g_xy_0_xxxxz_xxyz, g_xy_0_xxxxz_xxzz, g_xy_0_xxxxz_xyyy, g_xy_0_xxxxz_xyyz, g_xy_0_xxxxz_xyzz, g_xy_0_xxxxz_xzzz, g_xy_0_xxxxz_yyyy, g_xy_0_xxxxz_yyyz, g_xy_0_xxxxz_yyzz, g_xy_0_xxxxz_yzzz, g_xy_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxz_xxxx[k] = -g_xy_0_xxxx_xxxx[k] * ab_z + g_xy_0_xxxx_xxxxz[k];

                g_xy_0_xxxxz_xxxy[k] = -g_xy_0_xxxx_xxxy[k] * ab_z + g_xy_0_xxxx_xxxyz[k];

                g_xy_0_xxxxz_xxxz[k] = -g_xy_0_xxxx_xxxz[k] * ab_z + g_xy_0_xxxx_xxxzz[k];

                g_xy_0_xxxxz_xxyy[k] = -g_xy_0_xxxx_xxyy[k] * ab_z + g_xy_0_xxxx_xxyyz[k];

                g_xy_0_xxxxz_xxyz[k] = -g_xy_0_xxxx_xxyz[k] * ab_z + g_xy_0_xxxx_xxyzz[k];

                g_xy_0_xxxxz_xxzz[k] = -g_xy_0_xxxx_xxzz[k] * ab_z + g_xy_0_xxxx_xxzzz[k];

                g_xy_0_xxxxz_xyyy[k] = -g_xy_0_xxxx_xyyy[k] * ab_z + g_xy_0_xxxx_xyyyz[k];

                g_xy_0_xxxxz_xyyz[k] = -g_xy_0_xxxx_xyyz[k] * ab_z + g_xy_0_xxxx_xyyzz[k];

                g_xy_0_xxxxz_xyzz[k] = -g_xy_0_xxxx_xyzz[k] * ab_z + g_xy_0_xxxx_xyzzz[k];

                g_xy_0_xxxxz_xzzz[k] = -g_xy_0_xxxx_xzzz[k] * ab_z + g_xy_0_xxxx_xzzzz[k];

                g_xy_0_xxxxz_yyyy[k] = -g_xy_0_xxxx_yyyy[k] * ab_z + g_xy_0_xxxx_yyyyz[k];

                g_xy_0_xxxxz_yyyz[k] = -g_xy_0_xxxx_yyyz[k] * ab_z + g_xy_0_xxxx_yyyzz[k];

                g_xy_0_xxxxz_yyzz[k] = -g_xy_0_xxxx_yyzz[k] * ab_z + g_xy_0_xxxx_yyzzz[k];

                g_xy_0_xxxxz_yzzz[k] = -g_xy_0_xxxx_yzzz[k] * ab_z + g_xy_0_xxxx_yzzzz[k];

                g_xy_0_xxxxz_zzzz[k] = -g_xy_0_xxxx_zzzz[k] * ab_z + g_xy_0_xxxx_zzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 360 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 361 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 362 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 363 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 364 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 365 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 366 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 367 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 368 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 369 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 370 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 371 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 372 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 373 * ccomps * dcomps);

            auto g_xy_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_xxxx, g_xy_0_xxxyy_xxxy, g_xy_0_xxxyy_xxxz, g_xy_0_xxxyy_xxyy, g_xy_0_xxxyy_xxyz, g_xy_0_xxxyy_xxzz, g_xy_0_xxxyy_xyyy, g_xy_0_xxxyy_xyyz, g_xy_0_xxxyy_xyzz, g_xy_0_xxxyy_xzzz, g_xy_0_xxxyy_yyyy, g_xy_0_xxxyy_yyyz, g_xy_0_xxxyy_yyzz, g_xy_0_xxxyy_yzzz, g_xy_0_xxxyy_zzzz, g_xy_0_xxyy_xxxx, g_xy_0_xxyy_xxxxx, g_xy_0_xxyy_xxxxy, g_xy_0_xxyy_xxxxz, g_xy_0_xxyy_xxxy, g_xy_0_xxyy_xxxyy, g_xy_0_xxyy_xxxyz, g_xy_0_xxyy_xxxz, g_xy_0_xxyy_xxxzz, g_xy_0_xxyy_xxyy, g_xy_0_xxyy_xxyyy, g_xy_0_xxyy_xxyyz, g_xy_0_xxyy_xxyz, g_xy_0_xxyy_xxyzz, g_xy_0_xxyy_xxzz, g_xy_0_xxyy_xxzzz, g_xy_0_xxyy_xyyy, g_xy_0_xxyy_xyyyy, g_xy_0_xxyy_xyyyz, g_xy_0_xxyy_xyyz, g_xy_0_xxyy_xyyzz, g_xy_0_xxyy_xyzz, g_xy_0_xxyy_xyzzz, g_xy_0_xxyy_xzzz, g_xy_0_xxyy_xzzzz, g_xy_0_xxyy_yyyy, g_xy_0_xxyy_yyyz, g_xy_0_xxyy_yyzz, g_xy_0_xxyy_yzzz, g_xy_0_xxyy_zzzz, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_yyyy, g_y_0_xxyy_yyyz, g_y_0_xxyy_yyzz, g_y_0_xxyy_yzzz, g_y_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyy_xxxx[k] = -g_y_0_xxyy_xxxx[k] - g_xy_0_xxyy_xxxx[k] * ab_x + g_xy_0_xxyy_xxxxx[k];

                g_xy_0_xxxyy_xxxy[k] = -g_y_0_xxyy_xxxy[k] - g_xy_0_xxyy_xxxy[k] * ab_x + g_xy_0_xxyy_xxxxy[k];

                g_xy_0_xxxyy_xxxz[k] = -g_y_0_xxyy_xxxz[k] - g_xy_0_xxyy_xxxz[k] * ab_x + g_xy_0_xxyy_xxxxz[k];

                g_xy_0_xxxyy_xxyy[k] = -g_y_0_xxyy_xxyy[k] - g_xy_0_xxyy_xxyy[k] * ab_x + g_xy_0_xxyy_xxxyy[k];

                g_xy_0_xxxyy_xxyz[k] = -g_y_0_xxyy_xxyz[k] - g_xy_0_xxyy_xxyz[k] * ab_x + g_xy_0_xxyy_xxxyz[k];

                g_xy_0_xxxyy_xxzz[k] = -g_y_0_xxyy_xxzz[k] - g_xy_0_xxyy_xxzz[k] * ab_x + g_xy_0_xxyy_xxxzz[k];

                g_xy_0_xxxyy_xyyy[k] = -g_y_0_xxyy_xyyy[k] - g_xy_0_xxyy_xyyy[k] * ab_x + g_xy_0_xxyy_xxyyy[k];

                g_xy_0_xxxyy_xyyz[k] = -g_y_0_xxyy_xyyz[k] - g_xy_0_xxyy_xyyz[k] * ab_x + g_xy_0_xxyy_xxyyz[k];

                g_xy_0_xxxyy_xyzz[k] = -g_y_0_xxyy_xyzz[k] - g_xy_0_xxyy_xyzz[k] * ab_x + g_xy_0_xxyy_xxyzz[k];

                g_xy_0_xxxyy_xzzz[k] = -g_y_0_xxyy_xzzz[k] - g_xy_0_xxyy_xzzz[k] * ab_x + g_xy_0_xxyy_xxzzz[k];

                g_xy_0_xxxyy_yyyy[k] = -g_y_0_xxyy_yyyy[k] - g_xy_0_xxyy_yyyy[k] * ab_x + g_xy_0_xxyy_xyyyy[k];

                g_xy_0_xxxyy_yyyz[k] = -g_y_0_xxyy_yyyz[k] - g_xy_0_xxyy_yyyz[k] * ab_x + g_xy_0_xxyy_xyyyz[k];

                g_xy_0_xxxyy_yyzz[k] = -g_y_0_xxyy_yyzz[k] - g_xy_0_xxyy_yyzz[k] * ab_x + g_xy_0_xxyy_xyyzz[k];

                g_xy_0_xxxyy_yzzz[k] = -g_y_0_xxyy_yzzz[k] - g_xy_0_xxyy_yzzz[k] * ab_x + g_xy_0_xxyy_xyzzz[k];

                g_xy_0_xxxyy_zzzz[k] = -g_y_0_xxyy_zzzz[k] - g_xy_0_xxyy_zzzz[k] * ab_x + g_xy_0_xxyy_xzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 375 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 376 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 377 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 378 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 379 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 380 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 381 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 382 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 383 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 384 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 385 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 386 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 387 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 388 * ccomps * dcomps);

            auto g_xy_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_xxxx, g_xy_0_xxxy_xxxxz, g_xy_0_xxxy_xxxy, g_xy_0_xxxy_xxxyz, g_xy_0_xxxy_xxxz, g_xy_0_xxxy_xxxzz, g_xy_0_xxxy_xxyy, g_xy_0_xxxy_xxyyz, g_xy_0_xxxy_xxyz, g_xy_0_xxxy_xxyzz, g_xy_0_xxxy_xxzz, g_xy_0_xxxy_xxzzz, g_xy_0_xxxy_xyyy, g_xy_0_xxxy_xyyyz, g_xy_0_xxxy_xyyz, g_xy_0_xxxy_xyyzz, g_xy_0_xxxy_xyzz, g_xy_0_xxxy_xyzzz, g_xy_0_xxxy_xzzz, g_xy_0_xxxy_xzzzz, g_xy_0_xxxy_yyyy, g_xy_0_xxxy_yyyyz, g_xy_0_xxxy_yyyz, g_xy_0_xxxy_yyyzz, g_xy_0_xxxy_yyzz, g_xy_0_xxxy_yyzzz, g_xy_0_xxxy_yzzz, g_xy_0_xxxy_yzzzz, g_xy_0_xxxy_zzzz, g_xy_0_xxxy_zzzzz, g_xy_0_xxxyz_xxxx, g_xy_0_xxxyz_xxxy, g_xy_0_xxxyz_xxxz, g_xy_0_xxxyz_xxyy, g_xy_0_xxxyz_xxyz, g_xy_0_xxxyz_xxzz, g_xy_0_xxxyz_xyyy, g_xy_0_xxxyz_xyyz, g_xy_0_xxxyz_xyzz, g_xy_0_xxxyz_xzzz, g_xy_0_xxxyz_yyyy, g_xy_0_xxxyz_yyyz, g_xy_0_xxxyz_yyzz, g_xy_0_xxxyz_yzzz, g_xy_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyz_xxxx[k] = -g_xy_0_xxxy_xxxx[k] * ab_z + g_xy_0_xxxy_xxxxz[k];

                g_xy_0_xxxyz_xxxy[k] = -g_xy_0_xxxy_xxxy[k] * ab_z + g_xy_0_xxxy_xxxyz[k];

                g_xy_0_xxxyz_xxxz[k] = -g_xy_0_xxxy_xxxz[k] * ab_z + g_xy_0_xxxy_xxxzz[k];

                g_xy_0_xxxyz_xxyy[k] = -g_xy_0_xxxy_xxyy[k] * ab_z + g_xy_0_xxxy_xxyyz[k];

                g_xy_0_xxxyz_xxyz[k] = -g_xy_0_xxxy_xxyz[k] * ab_z + g_xy_0_xxxy_xxyzz[k];

                g_xy_0_xxxyz_xxzz[k] = -g_xy_0_xxxy_xxzz[k] * ab_z + g_xy_0_xxxy_xxzzz[k];

                g_xy_0_xxxyz_xyyy[k] = -g_xy_0_xxxy_xyyy[k] * ab_z + g_xy_0_xxxy_xyyyz[k];

                g_xy_0_xxxyz_xyyz[k] = -g_xy_0_xxxy_xyyz[k] * ab_z + g_xy_0_xxxy_xyyzz[k];

                g_xy_0_xxxyz_xyzz[k] = -g_xy_0_xxxy_xyzz[k] * ab_z + g_xy_0_xxxy_xyzzz[k];

                g_xy_0_xxxyz_xzzz[k] = -g_xy_0_xxxy_xzzz[k] * ab_z + g_xy_0_xxxy_xzzzz[k];

                g_xy_0_xxxyz_yyyy[k] = -g_xy_0_xxxy_yyyy[k] * ab_z + g_xy_0_xxxy_yyyyz[k];

                g_xy_0_xxxyz_yyyz[k] = -g_xy_0_xxxy_yyyz[k] * ab_z + g_xy_0_xxxy_yyyzz[k];

                g_xy_0_xxxyz_yyzz[k] = -g_xy_0_xxxy_yyzz[k] * ab_z + g_xy_0_xxxy_yyzzz[k];

                g_xy_0_xxxyz_yzzz[k] = -g_xy_0_xxxy_yzzz[k] * ab_z + g_xy_0_xxxy_yzzzz[k];

                g_xy_0_xxxyz_zzzz[k] = -g_xy_0_xxxy_zzzz[k] * ab_z + g_xy_0_xxxy_zzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 390 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 391 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 392 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 393 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 394 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 395 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 396 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 397 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 398 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 399 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 400 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 401 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 402 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 403 * ccomps * dcomps);

            auto g_xy_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxz_xxxx, g_xy_0_xxxz_xxxxz, g_xy_0_xxxz_xxxy, g_xy_0_xxxz_xxxyz, g_xy_0_xxxz_xxxz, g_xy_0_xxxz_xxxzz, g_xy_0_xxxz_xxyy, g_xy_0_xxxz_xxyyz, g_xy_0_xxxz_xxyz, g_xy_0_xxxz_xxyzz, g_xy_0_xxxz_xxzz, g_xy_0_xxxz_xxzzz, g_xy_0_xxxz_xyyy, g_xy_0_xxxz_xyyyz, g_xy_0_xxxz_xyyz, g_xy_0_xxxz_xyyzz, g_xy_0_xxxz_xyzz, g_xy_0_xxxz_xyzzz, g_xy_0_xxxz_xzzz, g_xy_0_xxxz_xzzzz, g_xy_0_xxxz_yyyy, g_xy_0_xxxz_yyyyz, g_xy_0_xxxz_yyyz, g_xy_0_xxxz_yyyzz, g_xy_0_xxxz_yyzz, g_xy_0_xxxz_yyzzz, g_xy_0_xxxz_yzzz, g_xy_0_xxxz_yzzzz, g_xy_0_xxxz_zzzz, g_xy_0_xxxz_zzzzz, g_xy_0_xxxzz_xxxx, g_xy_0_xxxzz_xxxy, g_xy_0_xxxzz_xxxz, g_xy_0_xxxzz_xxyy, g_xy_0_xxxzz_xxyz, g_xy_0_xxxzz_xxzz, g_xy_0_xxxzz_xyyy, g_xy_0_xxxzz_xyyz, g_xy_0_xxxzz_xyzz, g_xy_0_xxxzz_xzzz, g_xy_0_xxxzz_yyyy, g_xy_0_xxxzz_yyyz, g_xy_0_xxxzz_yyzz, g_xy_0_xxxzz_yzzz, g_xy_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzz_xxxx[k] = -g_xy_0_xxxz_xxxx[k] * ab_z + g_xy_0_xxxz_xxxxz[k];

                g_xy_0_xxxzz_xxxy[k] = -g_xy_0_xxxz_xxxy[k] * ab_z + g_xy_0_xxxz_xxxyz[k];

                g_xy_0_xxxzz_xxxz[k] = -g_xy_0_xxxz_xxxz[k] * ab_z + g_xy_0_xxxz_xxxzz[k];

                g_xy_0_xxxzz_xxyy[k] = -g_xy_0_xxxz_xxyy[k] * ab_z + g_xy_0_xxxz_xxyyz[k];

                g_xy_0_xxxzz_xxyz[k] = -g_xy_0_xxxz_xxyz[k] * ab_z + g_xy_0_xxxz_xxyzz[k];

                g_xy_0_xxxzz_xxzz[k] = -g_xy_0_xxxz_xxzz[k] * ab_z + g_xy_0_xxxz_xxzzz[k];

                g_xy_0_xxxzz_xyyy[k] = -g_xy_0_xxxz_xyyy[k] * ab_z + g_xy_0_xxxz_xyyyz[k];

                g_xy_0_xxxzz_xyyz[k] = -g_xy_0_xxxz_xyyz[k] * ab_z + g_xy_0_xxxz_xyyzz[k];

                g_xy_0_xxxzz_xyzz[k] = -g_xy_0_xxxz_xyzz[k] * ab_z + g_xy_0_xxxz_xyzzz[k];

                g_xy_0_xxxzz_xzzz[k] = -g_xy_0_xxxz_xzzz[k] * ab_z + g_xy_0_xxxz_xzzzz[k];

                g_xy_0_xxxzz_yyyy[k] = -g_xy_0_xxxz_yyyy[k] * ab_z + g_xy_0_xxxz_yyyyz[k];

                g_xy_0_xxxzz_yyyz[k] = -g_xy_0_xxxz_yyyz[k] * ab_z + g_xy_0_xxxz_yyyzz[k];

                g_xy_0_xxxzz_yyzz[k] = -g_xy_0_xxxz_yyzz[k] * ab_z + g_xy_0_xxxz_yyzzz[k];

                g_xy_0_xxxzz_yzzz[k] = -g_xy_0_xxxz_yzzz[k] * ab_z + g_xy_0_xxxz_yzzzz[k];

                g_xy_0_xxxzz_zzzz[k] = -g_xy_0_xxxz_zzzz[k] * ab_z + g_xy_0_xxxz_zzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 405 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 406 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 407 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 408 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 409 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 410 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 411 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 412 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 413 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 414 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 415 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 416 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 417 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 418 * ccomps * dcomps);

            auto g_xy_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_xxxx, g_xy_0_xxyyy_xxxy, g_xy_0_xxyyy_xxxz, g_xy_0_xxyyy_xxyy, g_xy_0_xxyyy_xxyz, g_xy_0_xxyyy_xxzz, g_xy_0_xxyyy_xyyy, g_xy_0_xxyyy_xyyz, g_xy_0_xxyyy_xyzz, g_xy_0_xxyyy_xzzz, g_xy_0_xxyyy_yyyy, g_xy_0_xxyyy_yyyz, g_xy_0_xxyyy_yyzz, g_xy_0_xxyyy_yzzz, g_xy_0_xxyyy_zzzz, g_xy_0_xyyy_xxxx, g_xy_0_xyyy_xxxxx, g_xy_0_xyyy_xxxxy, g_xy_0_xyyy_xxxxz, g_xy_0_xyyy_xxxy, g_xy_0_xyyy_xxxyy, g_xy_0_xyyy_xxxyz, g_xy_0_xyyy_xxxz, g_xy_0_xyyy_xxxzz, g_xy_0_xyyy_xxyy, g_xy_0_xyyy_xxyyy, g_xy_0_xyyy_xxyyz, g_xy_0_xyyy_xxyz, g_xy_0_xyyy_xxyzz, g_xy_0_xyyy_xxzz, g_xy_0_xyyy_xxzzz, g_xy_0_xyyy_xyyy, g_xy_0_xyyy_xyyyy, g_xy_0_xyyy_xyyyz, g_xy_0_xyyy_xyyz, g_xy_0_xyyy_xyyzz, g_xy_0_xyyy_xyzz, g_xy_0_xyyy_xyzzz, g_xy_0_xyyy_xzzz, g_xy_0_xyyy_xzzzz, g_xy_0_xyyy_yyyy, g_xy_0_xyyy_yyyz, g_xy_0_xyyy_yyzz, g_xy_0_xyyy_yzzz, g_xy_0_xyyy_zzzz, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_yyyy, g_y_0_xyyy_yyyz, g_y_0_xyyy_yyzz, g_y_0_xyyy_yzzz, g_y_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyy_xxxx[k] = -g_y_0_xyyy_xxxx[k] - g_xy_0_xyyy_xxxx[k] * ab_x + g_xy_0_xyyy_xxxxx[k];

                g_xy_0_xxyyy_xxxy[k] = -g_y_0_xyyy_xxxy[k] - g_xy_0_xyyy_xxxy[k] * ab_x + g_xy_0_xyyy_xxxxy[k];

                g_xy_0_xxyyy_xxxz[k] = -g_y_0_xyyy_xxxz[k] - g_xy_0_xyyy_xxxz[k] * ab_x + g_xy_0_xyyy_xxxxz[k];

                g_xy_0_xxyyy_xxyy[k] = -g_y_0_xyyy_xxyy[k] - g_xy_0_xyyy_xxyy[k] * ab_x + g_xy_0_xyyy_xxxyy[k];

                g_xy_0_xxyyy_xxyz[k] = -g_y_0_xyyy_xxyz[k] - g_xy_0_xyyy_xxyz[k] * ab_x + g_xy_0_xyyy_xxxyz[k];

                g_xy_0_xxyyy_xxzz[k] = -g_y_0_xyyy_xxzz[k] - g_xy_0_xyyy_xxzz[k] * ab_x + g_xy_0_xyyy_xxxzz[k];

                g_xy_0_xxyyy_xyyy[k] = -g_y_0_xyyy_xyyy[k] - g_xy_0_xyyy_xyyy[k] * ab_x + g_xy_0_xyyy_xxyyy[k];

                g_xy_0_xxyyy_xyyz[k] = -g_y_0_xyyy_xyyz[k] - g_xy_0_xyyy_xyyz[k] * ab_x + g_xy_0_xyyy_xxyyz[k];

                g_xy_0_xxyyy_xyzz[k] = -g_y_0_xyyy_xyzz[k] - g_xy_0_xyyy_xyzz[k] * ab_x + g_xy_0_xyyy_xxyzz[k];

                g_xy_0_xxyyy_xzzz[k] = -g_y_0_xyyy_xzzz[k] - g_xy_0_xyyy_xzzz[k] * ab_x + g_xy_0_xyyy_xxzzz[k];

                g_xy_0_xxyyy_yyyy[k] = -g_y_0_xyyy_yyyy[k] - g_xy_0_xyyy_yyyy[k] * ab_x + g_xy_0_xyyy_xyyyy[k];

                g_xy_0_xxyyy_yyyz[k] = -g_y_0_xyyy_yyyz[k] - g_xy_0_xyyy_yyyz[k] * ab_x + g_xy_0_xyyy_xyyyz[k];

                g_xy_0_xxyyy_yyzz[k] = -g_y_0_xyyy_yyzz[k] - g_xy_0_xyyy_yyzz[k] * ab_x + g_xy_0_xyyy_xyyzz[k];

                g_xy_0_xxyyy_yzzz[k] = -g_y_0_xyyy_yzzz[k] - g_xy_0_xyyy_yzzz[k] * ab_x + g_xy_0_xyyy_xyzzz[k];

                g_xy_0_xxyyy_zzzz[k] = -g_y_0_xyyy_zzzz[k] - g_xy_0_xyyy_zzzz[k] * ab_x + g_xy_0_xyyy_xzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 420 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 421 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 422 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 423 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 424 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 425 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 426 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 427 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 428 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 429 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 430 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 431 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 432 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 433 * ccomps * dcomps);

            auto g_xy_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_xxxx, g_xy_0_xxyy_xxxxz, g_xy_0_xxyy_xxxy, g_xy_0_xxyy_xxxyz, g_xy_0_xxyy_xxxz, g_xy_0_xxyy_xxxzz, g_xy_0_xxyy_xxyy, g_xy_0_xxyy_xxyyz, g_xy_0_xxyy_xxyz, g_xy_0_xxyy_xxyzz, g_xy_0_xxyy_xxzz, g_xy_0_xxyy_xxzzz, g_xy_0_xxyy_xyyy, g_xy_0_xxyy_xyyyz, g_xy_0_xxyy_xyyz, g_xy_0_xxyy_xyyzz, g_xy_0_xxyy_xyzz, g_xy_0_xxyy_xyzzz, g_xy_0_xxyy_xzzz, g_xy_0_xxyy_xzzzz, g_xy_0_xxyy_yyyy, g_xy_0_xxyy_yyyyz, g_xy_0_xxyy_yyyz, g_xy_0_xxyy_yyyzz, g_xy_0_xxyy_yyzz, g_xy_0_xxyy_yyzzz, g_xy_0_xxyy_yzzz, g_xy_0_xxyy_yzzzz, g_xy_0_xxyy_zzzz, g_xy_0_xxyy_zzzzz, g_xy_0_xxyyz_xxxx, g_xy_0_xxyyz_xxxy, g_xy_0_xxyyz_xxxz, g_xy_0_xxyyz_xxyy, g_xy_0_xxyyz_xxyz, g_xy_0_xxyyz_xxzz, g_xy_0_xxyyz_xyyy, g_xy_0_xxyyz_xyyz, g_xy_0_xxyyz_xyzz, g_xy_0_xxyyz_xzzz, g_xy_0_xxyyz_yyyy, g_xy_0_xxyyz_yyyz, g_xy_0_xxyyz_yyzz, g_xy_0_xxyyz_yzzz, g_xy_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyz_xxxx[k] = -g_xy_0_xxyy_xxxx[k] * ab_z + g_xy_0_xxyy_xxxxz[k];

                g_xy_0_xxyyz_xxxy[k] = -g_xy_0_xxyy_xxxy[k] * ab_z + g_xy_0_xxyy_xxxyz[k];

                g_xy_0_xxyyz_xxxz[k] = -g_xy_0_xxyy_xxxz[k] * ab_z + g_xy_0_xxyy_xxxzz[k];

                g_xy_0_xxyyz_xxyy[k] = -g_xy_0_xxyy_xxyy[k] * ab_z + g_xy_0_xxyy_xxyyz[k];

                g_xy_0_xxyyz_xxyz[k] = -g_xy_0_xxyy_xxyz[k] * ab_z + g_xy_0_xxyy_xxyzz[k];

                g_xy_0_xxyyz_xxzz[k] = -g_xy_0_xxyy_xxzz[k] * ab_z + g_xy_0_xxyy_xxzzz[k];

                g_xy_0_xxyyz_xyyy[k] = -g_xy_0_xxyy_xyyy[k] * ab_z + g_xy_0_xxyy_xyyyz[k];

                g_xy_0_xxyyz_xyyz[k] = -g_xy_0_xxyy_xyyz[k] * ab_z + g_xy_0_xxyy_xyyzz[k];

                g_xy_0_xxyyz_xyzz[k] = -g_xy_0_xxyy_xyzz[k] * ab_z + g_xy_0_xxyy_xyzzz[k];

                g_xy_0_xxyyz_xzzz[k] = -g_xy_0_xxyy_xzzz[k] * ab_z + g_xy_0_xxyy_xzzzz[k];

                g_xy_0_xxyyz_yyyy[k] = -g_xy_0_xxyy_yyyy[k] * ab_z + g_xy_0_xxyy_yyyyz[k];

                g_xy_0_xxyyz_yyyz[k] = -g_xy_0_xxyy_yyyz[k] * ab_z + g_xy_0_xxyy_yyyzz[k];

                g_xy_0_xxyyz_yyzz[k] = -g_xy_0_xxyy_yyzz[k] * ab_z + g_xy_0_xxyy_yyzzz[k];

                g_xy_0_xxyyz_yzzz[k] = -g_xy_0_xxyy_yzzz[k] * ab_z + g_xy_0_xxyy_yzzzz[k];

                g_xy_0_xxyyz_zzzz[k] = -g_xy_0_xxyy_zzzz[k] * ab_z + g_xy_0_xxyy_zzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 435 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 436 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 437 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 438 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 439 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 440 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 441 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 442 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 443 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 444 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 445 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 446 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 447 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 448 * ccomps * dcomps);

            auto g_xy_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyz_xxxx, g_xy_0_xxyz_xxxxz, g_xy_0_xxyz_xxxy, g_xy_0_xxyz_xxxyz, g_xy_0_xxyz_xxxz, g_xy_0_xxyz_xxxzz, g_xy_0_xxyz_xxyy, g_xy_0_xxyz_xxyyz, g_xy_0_xxyz_xxyz, g_xy_0_xxyz_xxyzz, g_xy_0_xxyz_xxzz, g_xy_0_xxyz_xxzzz, g_xy_0_xxyz_xyyy, g_xy_0_xxyz_xyyyz, g_xy_0_xxyz_xyyz, g_xy_0_xxyz_xyyzz, g_xy_0_xxyz_xyzz, g_xy_0_xxyz_xyzzz, g_xy_0_xxyz_xzzz, g_xy_0_xxyz_xzzzz, g_xy_0_xxyz_yyyy, g_xy_0_xxyz_yyyyz, g_xy_0_xxyz_yyyz, g_xy_0_xxyz_yyyzz, g_xy_0_xxyz_yyzz, g_xy_0_xxyz_yyzzz, g_xy_0_xxyz_yzzz, g_xy_0_xxyz_yzzzz, g_xy_0_xxyz_zzzz, g_xy_0_xxyz_zzzzz, g_xy_0_xxyzz_xxxx, g_xy_0_xxyzz_xxxy, g_xy_0_xxyzz_xxxz, g_xy_0_xxyzz_xxyy, g_xy_0_xxyzz_xxyz, g_xy_0_xxyzz_xxzz, g_xy_0_xxyzz_xyyy, g_xy_0_xxyzz_xyyz, g_xy_0_xxyzz_xyzz, g_xy_0_xxyzz_xzzz, g_xy_0_xxyzz_yyyy, g_xy_0_xxyzz_yyyz, g_xy_0_xxyzz_yyzz, g_xy_0_xxyzz_yzzz, g_xy_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzz_xxxx[k] = -g_xy_0_xxyz_xxxx[k] * ab_z + g_xy_0_xxyz_xxxxz[k];

                g_xy_0_xxyzz_xxxy[k] = -g_xy_0_xxyz_xxxy[k] * ab_z + g_xy_0_xxyz_xxxyz[k];

                g_xy_0_xxyzz_xxxz[k] = -g_xy_0_xxyz_xxxz[k] * ab_z + g_xy_0_xxyz_xxxzz[k];

                g_xy_0_xxyzz_xxyy[k] = -g_xy_0_xxyz_xxyy[k] * ab_z + g_xy_0_xxyz_xxyyz[k];

                g_xy_0_xxyzz_xxyz[k] = -g_xy_0_xxyz_xxyz[k] * ab_z + g_xy_0_xxyz_xxyzz[k];

                g_xy_0_xxyzz_xxzz[k] = -g_xy_0_xxyz_xxzz[k] * ab_z + g_xy_0_xxyz_xxzzz[k];

                g_xy_0_xxyzz_xyyy[k] = -g_xy_0_xxyz_xyyy[k] * ab_z + g_xy_0_xxyz_xyyyz[k];

                g_xy_0_xxyzz_xyyz[k] = -g_xy_0_xxyz_xyyz[k] * ab_z + g_xy_0_xxyz_xyyzz[k];

                g_xy_0_xxyzz_xyzz[k] = -g_xy_0_xxyz_xyzz[k] * ab_z + g_xy_0_xxyz_xyzzz[k];

                g_xy_0_xxyzz_xzzz[k] = -g_xy_0_xxyz_xzzz[k] * ab_z + g_xy_0_xxyz_xzzzz[k];

                g_xy_0_xxyzz_yyyy[k] = -g_xy_0_xxyz_yyyy[k] * ab_z + g_xy_0_xxyz_yyyyz[k];

                g_xy_0_xxyzz_yyyz[k] = -g_xy_0_xxyz_yyyz[k] * ab_z + g_xy_0_xxyz_yyyzz[k];

                g_xy_0_xxyzz_yyzz[k] = -g_xy_0_xxyz_yyzz[k] * ab_z + g_xy_0_xxyz_yyzzz[k];

                g_xy_0_xxyzz_yzzz[k] = -g_xy_0_xxyz_yzzz[k] * ab_z + g_xy_0_xxyz_yzzzz[k];

                g_xy_0_xxyzz_zzzz[k] = -g_xy_0_xxyz_zzzz[k] * ab_z + g_xy_0_xxyz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 450 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 451 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 452 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 453 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 454 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 455 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 456 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 457 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 458 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 459 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 460 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 461 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 462 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 463 * ccomps * dcomps);

            auto g_xy_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzz_xxxx, g_xy_0_xxzz_xxxxz, g_xy_0_xxzz_xxxy, g_xy_0_xxzz_xxxyz, g_xy_0_xxzz_xxxz, g_xy_0_xxzz_xxxzz, g_xy_0_xxzz_xxyy, g_xy_0_xxzz_xxyyz, g_xy_0_xxzz_xxyz, g_xy_0_xxzz_xxyzz, g_xy_0_xxzz_xxzz, g_xy_0_xxzz_xxzzz, g_xy_0_xxzz_xyyy, g_xy_0_xxzz_xyyyz, g_xy_0_xxzz_xyyz, g_xy_0_xxzz_xyyzz, g_xy_0_xxzz_xyzz, g_xy_0_xxzz_xyzzz, g_xy_0_xxzz_xzzz, g_xy_0_xxzz_xzzzz, g_xy_0_xxzz_yyyy, g_xy_0_xxzz_yyyyz, g_xy_0_xxzz_yyyz, g_xy_0_xxzz_yyyzz, g_xy_0_xxzz_yyzz, g_xy_0_xxzz_yyzzz, g_xy_0_xxzz_yzzz, g_xy_0_xxzz_yzzzz, g_xy_0_xxzz_zzzz, g_xy_0_xxzz_zzzzz, g_xy_0_xxzzz_xxxx, g_xy_0_xxzzz_xxxy, g_xy_0_xxzzz_xxxz, g_xy_0_xxzzz_xxyy, g_xy_0_xxzzz_xxyz, g_xy_0_xxzzz_xxzz, g_xy_0_xxzzz_xyyy, g_xy_0_xxzzz_xyyz, g_xy_0_xxzzz_xyzz, g_xy_0_xxzzz_xzzz, g_xy_0_xxzzz_yyyy, g_xy_0_xxzzz_yyyz, g_xy_0_xxzzz_yyzz, g_xy_0_xxzzz_yzzz, g_xy_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzz_xxxx[k] = -g_xy_0_xxzz_xxxx[k] * ab_z + g_xy_0_xxzz_xxxxz[k];

                g_xy_0_xxzzz_xxxy[k] = -g_xy_0_xxzz_xxxy[k] * ab_z + g_xy_0_xxzz_xxxyz[k];

                g_xy_0_xxzzz_xxxz[k] = -g_xy_0_xxzz_xxxz[k] * ab_z + g_xy_0_xxzz_xxxzz[k];

                g_xy_0_xxzzz_xxyy[k] = -g_xy_0_xxzz_xxyy[k] * ab_z + g_xy_0_xxzz_xxyyz[k];

                g_xy_0_xxzzz_xxyz[k] = -g_xy_0_xxzz_xxyz[k] * ab_z + g_xy_0_xxzz_xxyzz[k];

                g_xy_0_xxzzz_xxzz[k] = -g_xy_0_xxzz_xxzz[k] * ab_z + g_xy_0_xxzz_xxzzz[k];

                g_xy_0_xxzzz_xyyy[k] = -g_xy_0_xxzz_xyyy[k] * ab_z + g_xy_0_xxzz_xyyyz[k];

                g_xy_0_xxzzz_xyyz[k] = -g_xy_0_xxzz_xyyz[k] * ab_z + g_xy_0_xxzz_xyyzz[k];

                g_xy_0_xxzzz_xyzz[k] = -g_xy_0_xxzz_xyzz[k] * ab_z + g_xy_0_xxzz_xyzzz[k];

                g_xy_0_xxzzz_xzzz[k] = -g_xy_0_xxzz_xzzz[k] * ab_z + g_xy_0_xxzz_xzzzz[k];

                g_xy_0_xxzzz_yyyy[k] = -g_xy_0_xxzz_yyyy[k] * ab_z + g_xy_0_xxzz_yyyyz[k];

                g_xy_0_xxzzz_yyyz[k] = -g_xy_0_xxzz_yyyz[k] * ab_z + g_xy_0_xxzz_yyyzz[k];

                g_xy_0_xxzzz_yyzz[k] = -g_xy_0_xxzz_yyzz[k] * ab_z + g_xy_0_xxzz_yyzzz[k];

                g_xy_0_xxzzz_yzzz[k] = -g_xy_0_xxzz_yzzz[k] * ab_z + g_xy_0_xxzz_yzzzz[k];

                g_xy_0_xxzzz_zzzz[k] = -g_xy_0_xxzz_zzzz[k] * ab_z + g_xy_0_xxzz_zzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 465 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 466 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 467 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 468 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 469 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 470 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 471 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 472 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 473 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 474 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 475 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 476 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 477 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 478 * ccomps * dcomps);

            auto g_xy_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_xxxx, g_xy_0_xyyyy_xxxy, g_xy_0_xyyyy_xxxz, g_xy_0_xyyyy_xxyy, g_xy_0_xyyyy_xxyz, g_xy_0_xyyyy_xxzz, g_xy_0_xyyyy_xyyy, g_xy_0_xyyyy_xyyz, g_xy_0_xyyyy_xyzz, g_xy_0_xyyyy_xzzz, g_xy_0_xyyyy_yyyy, g_xy_0_xyyyy_yyyz, g_xy_0_xyyyy_yyzz, g_xy_0_xyyyy_yzzz, g_xy_0_xyyyy_zzzz, g_xy_0_yyyy_xxxx, g_xy_0_yyyy_xxxxx, g_xy_0_yyyy_xxxxy, g_xy_0_yyyy_xxxxz, g_xy_0_yyyy_xxxy, g_xy_0_yyyy_xxxyy, g_xy_0_yyyy_xxxyz, g_xy_0_yyyy_xxxz, g_xy_0_yyyy_xxxzz, g_xy_0_yyyy_xxyy, g_xy_0_yyyy_xxyyy, g_xy_0_yyyy_xxyyz, g_xy_0_yyyy_xxyz, g_xy_0_yyyy_xxyzz, g_xy_0_yyyy_xxzz, g_xy_0_yyyy_xxzzz, g_xy_0_yyyy_xyyy, g_xy_0_yyyy_xyyyy, g_xy_0_yyyy_xyyyz, g_xy_0_yyyy_xyyz, g_xy_0_yyyy_xyyzz, g_xy_0_yyyy_xyzz, g_xy_0_yyyy_xyzzz, g_xy_0_yyyy_xzzz, g_xy_0_yyyy_xzzzz, g_xy_0_yyyy_yyyy, g_xy_0_yyyy_yyyz, g_xy_0_yyyy_yyzz, g_xy_0_yyyy_yzzz, g_xy_0_yyyy_zzzz, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyy_xxxx[k] = -g_y_0_yyyy_xxxx[k] - g_xy_0_yyyy_xxxx[k] * ab_x + g_xy_0_yyyy_xxxxx[k];

                g_xy_0_xyyyy_xxxy[k] = -g_y_0_yyyy_xxxy[k] - g_xy_0_yyyy_xxxy[k] * ab_x + g_xy_0_yyyy_xxxxy[k];

                g_xy_0_xyyyy_xxxz[k] = -g_y_0_yyyy_xxxz[k] - g_xy_0_yyyy_xxxz[k] * ab_x + g_xy_0_yyyy_xxxxz[k];

                g_xy_0_xyyyy_xxyy[k] = -g_y_0_yyyy_xxyy[k] - g_xy_0_yyyy_xxyy[k] * ab_x + g_xy_0_yyyy_xxxyy[k];

                g_xy_0_xyyyy_xxyz[k] = -g_y_0_yyyy_xxyz[k] - g_xy_0_yyyy_xxyz[k] * ab_x + g_xy_0_yyyy_xxxyz[k];

                g_xy_0_xyyyy_xxzz[k] = -g_y_0_yyyy_xxzz[k] - g_xy_0_yyyy_xxzz[k] * ab_x + g_xy_0_yyyy_xxxzz[k];

                g_xy_0_xyyyy_xyyy[k] = -g_y_0_yyyy_xyyy[k] - g_xy_0_yyyy_xyyy[k] * ab_x + g_xy_0_yyyy_xxyyy[k];

                g_xy_0_xyyyy_xyyz[k] = -g_y_0_yyyy_xyyz[k] - g_xy_0_yyyy_xyyz[k] * ab_x + g_xy_0_yyyy_xxyyz[k];

                g_xy_0_xyyyy_xyzz[k] = -g_y_0_yyyy_xyzz[k] - g_xy_0_yyyy_xyzz[k] * ab_x + g_xy_0_yyyy_xxyzz[k];

                g_xy_0_xyyyy_xzzz[k] = -g_y_0_yyyy_xzzz[k] - g_xy_0_yyyy_xzzz[k] * ab_x + g_xy_0_yyyy_xxzzz[k];

                g_xy_0_xyyyy_yyyy[k] = -g_y_0_yyyy_yyyy[k] - g_xy_0_yyyy_yyyy[k] * ab_x + g_xy_0_yyyy_xyyyy[k];

                g_xy_0_xyyyy_yyyz[k] = -g_y_0_yyyy_yyyz[k] - g_xy_0_yyyy_yyyz[k] * ab_x + g_xy_0_yyyy_xyyyz[k];

                g_xy_0_xyyyy_yyzz[k] = -g_y_0_yyyy_yyzz[k] - g_xy_0_yyyy_yyzz[k] * ab_x + g_xy_0_yyyy_xyyzz[k];

                g_xy_0_xyyyy_yzzz[k] = -g_y_0_yyyy_yzzz[k] - g_xy_0_yyyy_yzzz[k] * ab_x + g_xy_0_yyyy_xyzzz[k];

                g_xy_0_xyyyy_zzzz[k] = -g_y_0_yyyy_zzzz[k] - g_xy_0_yyyy_zzzz[k] * ab_x + g_xy_0_yyyy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 480 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 481 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 482 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 483 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 484 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 485 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 486 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 487 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 488 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 489 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 490 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 491 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 492 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 493 * ccomps * dcomps);

            auto g_xy_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_xxxx, g_xy_0_xyyy_xxxxz, g_xy_0_xyyy_xxxy, g_xy_0_xyyy_xxxyz, g_xy_0_xyyy_xxxz, g_xy_0_xyyy_xxxzz, g_xy_0_xyyy_xxyy, g_xy_0_xyyy_xxyyz, g_xy_0_xyyy_xxyz, g_xy_0_xyyy_xxyzz, g_xy_0_xyyy_xxzz, g_xy_0_xyyy_xxzzz, g_xy_0_xyyy_xyyy, g_xy_0_xyyy_xyyyz, g_xy_0_xyyy_xyyz, g_xy_0_xyyy_xyyzz, g_xy_0_xyyy_xyzz, g_xy_0_xyyy_xyzzz, g_xy_0_xyyy_xzzz, g_xy_0_xyyy_xzzzz, g_xy_0_xyyy_yyyy, g_xy_0_xyyy_yyyyz, g_xy_0_xyyy_yyyz, g_xy_0_xyyy_yyyzz, g_xy_0_xyyy_yyzz, g_xy_0_xyyy_yyzzz, g_xy_0_xyyy_yzzz, g_xy_0_xyyy_yzzzz, g_xy_0_xyyy_zzzz, g_xy_0_xyyy_zzzzz, g_xy_0_xyyyz_xxxx, g_xy_0_xyyyz_xxxy, g_xy_0_xyyyz_xxxz, g_xy_0_xyyyz_xxyy, g_xy_0_xyyyz_xxyz, g_xy_0_xyyyz_xxzz, g_xy_0_xyyyz_xyyy, g_xy_0_xyyyz_xyyz, g_xy_0_xyyyz_xyzz, g_xy_0_xyyyz_xzzz, g_xy_0_xyyyz_yyyy, g_xy_0_xyyyz_yyyz, g_xy_0_xyyyz_yyzz, g_xy_0_xyyyz_yzzz, g_xy_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyz_xxxx[k] = -g_xy_0_xyyy_xxxx[k] * ab_z + g_xy_0_xyyy_xxxxz[k];

                g_xy_0_xyyyz_xxxy[k] = -g_xy_0_xyyy_xxxy[k] * ab_z + g_xy_0_xyyy_xxxyz[k];

                g_xy_0_xyyyz_xxxz[k] = -g_xy_0_xyyy_xxxz[k] * ab_z + g_xy_0_xyyy_xxxzz[k];

                g_xy_0_xyyyz_xxyy[k] = -g_xy_0_xyyy_xxyy[k] * ab_z + g_xy_0_xyyy_xxyyz[k];

                g_xy_0_xyyyz_xxyz[k] = -g_xy_0_xyyy_xxyz[k] * ab_z + g_xy_0_xyyy_xxyzz[k];

                g_xy_0_xyyyz_xxzz[k] = -g_xy_0_xyyy_xxzz[k] * ab_z + g_xy_0_xyyy_xxzzz[k];

                g_xy_0_xyyyz_xyyy[k] = -g_xy_0_xyyy_xyyy[k] * ab_z + g_xy_0_xyyy_xyyyz[k];

                g_xy_0_xyyyz_xyyz[k] = -g_xy_0_xyyy_xyyz[k] * ab_z + g_xy_0_xyyy_xyyzz[k];

                g_xy_0_xyyyz_xyzz[k] = -g_xy_0_xyyy_xyzz[k] * ab_z + g_xy_0_xyyy_xyzzz[k];

                g_xy_0_xyyyz_xzzz[k] = -g_xy_0_xyyy_xzzz[k] * ab_z + g_xy_0_xyyy_xzzzz[k];

                g_xy_0_xyyyz_yyyy[k] = -g_xy_0_xyyy_yyyy[k] * ab_z + g_xy_0_xyyy_yyyyz[k];

                g_xy_0_xyyyz_yyyz[k] = -g_xy_0_xyyy_yyyz[k] * ab_z + g_xy_0_xyyy_yyyzz[k];

                g_xy_0_xyyyz_yyzz[k] = -g_xy_0_xyyy_yyzz[k] * ab_z + g_xy_0_xyyy_yyzzz[k];

                g_xy_0_xyyyz_yzzz[k] = -g_xy_0_xyyy_yzzz[k] * ab_z + g_xy_0_xyyy_yzzzz[k];

                g_xy_0_xyyyz_zzzz[k] = -g_xy_0_xyyy_zzzz[k] * ab_z + g_xy_0_xyyy_zzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 495 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 496 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 497 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 498 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 499 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 500 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 501 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 502 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 503 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 504 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 505 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 506 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 507 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 508 * ccomps * dcomps);

            auto g_xy_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyz_xxxx, g_xy_0_xyyz_xxxxz, g_xy_0_xyyz_xxxy, g_xy_0_xyyz_xxxyz, g_xy_0_xyyz_xxxz, g_xy_0_xyyz_xxxzz, g_xy_0_xyyz_xxyy, g_xy_0_xyyz_xxyyz, g_xy_0_xyyz_xxyz, g_xy_0_xyyz_xxyzz, g_xy_0_xyyz_xxzz, g_xy_0_xyyz_xxzzz, g_xy_0_xyyz_xyyy, g_xy_0_xyyz_xyyyz, g_xy_0_xyyz_xyyz, g_xy_0_xyyz_xyyzz, g_xy_0_xyyz_xyzz, g_xy_0_xyyz_xyzzz, g_xy_0_xyyz_xzzz, g_xy_0_xyyz_xzzzz, g_xy_0_xyyz_yyyy, g_xy_0_xyyz_yyyyz, g_xy_0_xyyz_yyyz, g_xy_0_xyyz_yyyzz, g_xy_0_xyyz_yyzz, g_xy_0_xyyz_yyzzz, g_xy_0_xyyz_yzzz, g_xy_0_xyyz_yzzzz, g_xy_0_xyyz_zzzz, g_xy_0_xyyz_zzzzz, g_xy_0_xyyzz_xxxx, g_xy_0_xyyzz_xxxy, g_xy_0_xyyzz_xxxz, g_xy_0_xyyzz_xxyy, g_xy_0_xyyzz_xxyz, g_xy_0_xyyzz_xxzz, g_xy_0_xyyzz_xyyy, g_xy_0_xyyzz_xyyz, g_xy_0_xyyzz_xyzz, g_xy_0_xyyzz_xzzz, g_xy_0_xyyzz_yyyy, g_xy_0_xyyzz_yyyz, g_xy_0_xyyzz_yyzz, g_xy_0_xyyzz_yzzz, g_xy_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzz_xxxx[k] = -g_xy_0_xyyz_xxxx[k] * ab_z + g_xy_0_xyyz_xxxxz[k];

                g_xy_0_xyyzz_xxxy[k] = -g_xy_0_xyyz_xxxy[k] * ab_z + g_xy_0_xyyz_xxxyz[k];

                g_xy_0_xyyzz_xxxz[k] = -g_xy_0_xyyz_xxxz[k] * ab_z + g_xy_0_xyyz_xxxzz[k];

                g_xy_0_xyyzz_xxyy[k] = -g_xy_0_xyyz_xxyy[k] * ab_z + g_xy_0_xyyz_xxyyz[k];

                g_xy_0_xyyzz_xxyz[k] = -g_xy_0_xyyz_xxyz[k] * ab_z + g_xy_0_xyyz_xxyzz[k];

                g_xy_0_xyyzz_xxzz[k] = -g_xy_0_xyyz_xxzz[k] * ab_z + g_xy_0_xyyz_xxzzz[k];

                g_xy_0_xyyzz_xyyy[k] = -g_xy_0_xyyz_xyyy[k] * ab_z + g_xy_0_xyyz_xyyyz[k];

                g_xy_0_xyyzz_xyyz[k] = -g_xy_0_xyyz_xyyz[k] * ab_z + g_xy_0_xyyz_xyyzz[k];

                g_xy_0_xyyzz_xyzz[k] = -g_xy_0_xyyz_xyzz[k] * ab_z + g_xy_0_xyyz_xyzzz[k];

                g_xy_0_xyyzz_xzzz[k] = -g_xy_0_xyyz_xzzz[k] * ab_z + g_xy_0_xyyz_xzzzz[k];

                g_xy_0_xyyzz_yyyy[k] = -g_xy_0_xyyz_yyyy[k] * ab_z + g_xy_0_xyyz_yyyyz[k];

                g_xy_0_xyyzz_yyyz[k] = -g_xy_0_xyyz_yyyz[k] * ab_z + g_xy_0_xyyz_yyyzz[k];

                g_xy_0_xyyzz_yyzz[k] = -g_xy_0_xyyz_yyzz[k] * ab_z + g_xy_0_xyyz_yyzzz[k];

                g_xy_0_xyyzz_yzzz[k] = -g_xy_0_xyyz_yzzz[k] * ab_z + g_xy_0_xyyz_yzzzz[k];

                g_xy_0_xyyzz_zzzz[k] = -g_xy_0_xyyz_zzzz[k] * ab_z + g_xy_0_xyyz_zzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 510 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 511 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 512 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 513 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 514 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 515 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 516 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 517 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 518 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 519 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 520 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 521 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 522 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 523 * ccomps * dcomps);

            auto g_xy_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzz_xxxx, g_xy_0_xyzz_xxxxz, g_xy_0_xyzz_xxxy, g_xy_0_xyzz_xxxyz, g_xy_0_xyzz_xxxz, g_xy_0_xyzz_xxxzz, g_xy_0_xyzz_xxyy, g_xy_0_xyzz_xxyyz, g_xy_0_xyzz_xxyz, g_xy_0_xyzz_xxyzz, g_xy_0_xyzz_xxzz, g_xy_0_xyzz_xxzzz, g_xy_0_xyzz_xyyy, g_xy_0_xyzz_xyyyz, g_xy_0_xyzz_xyyz, g_xy_0_xyzz_xyyzz, g_xy_0_xyzz_xyzz, g_xy_0_xyzz_xyzzz, g_xy_0_xyzz_xzzz, g_xy_0_xyzz_xzzzz, g_xy_0_xyzz_yyyy, g_xy_0_xyzz_yyyyz, g_xy_0_xyzz_yyyz, g_xy_0_xyzz_yyyzz, g_xy_0_xyzz_yyzz, g_xy_0_xyzz_yyzzz, g_xy_0_xyzz_yzzz, g_xy_0_xyzz_yzzzz, g_xy_0_xyzz_zzzz, g_xy_0_xyzz_zzzzz, g_xy_0_xyzzz_xxxx, g_xy_0_xyzzz_xxxy, g_xy_0_xyzzz_xxxz, g_xy_0_xyzzz_xxyy, g_xy_0_xyzzz_xxyz, g_xy_0_xyzzz_xxzz, g_xy_0_xyzzz_xyyy, g_xy_0_xyzzz_xyyz, g_xy_0_xyzzz_xyzz, g_xy_0_xyzzz_xzzz, g_xy_0_xyzzz_yyyy, g_xy_0_xyzzz_yyyz, g_xy_0_xyzzz_yyzz, g_xy_0_xyzzz_yzzz, g_xy_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzz_xxxx[k] = -g_xy_0_xyzz_xxxx[k] * ab_z + g_xy_0_xyzz_xxxxz[k];

                g_xy_0_xyzzz_xxxy[k] = -g_xy_0_xyzz_xxxy[k] * ab_z + g_xy_0_xyzz_xxxyz[k];

                g_xy_0_xyzzz_xxxz[k] = -g_xy_0_xyzz_xxxz[k] * ab_z + g_xy_0_xyzz_xxxzz[k];

                g_xy_0_xyzzz_xxyy[k] = -g_xy_0_xyzz_xxyy[k] * ab_z + g_xy_0_xyzz_xxyyz[k];

                g_xy_0_xyzzz_xxyz[k] = -g_xy_0_xyzz_xxyz[k] * ab_z + g_xy_0_xyzz_xxyzz[k];

                g_xy_0_xyzzz_xxzz[k] = -g_xy_0_xyzz_xxzz[k] * ab_z + g_xy_0_xyzz_xxzzz[k];

                g_xy_0_xyzzz_xyyy[k] = -g_xy_0_xyzz_xyyy[k] * ab_z + g_xy_0_xyzz_xyyyz[k];

                g_xy_0_xyzzz_xyyz[k] = -g_xy_0_xyzz_xyyz[k] * ab_z + g_xy_0_xyzz_xyyzz[k];

                g_xy_0_xyzzz_xyzz[k] = -g_xy_0_xyzz_xyzz[k] * ab_z + g_xy_0_xyzz_xyzzz[k];

                g_xy_0_xyzzz_xzzz[k] = -g_xy_0_xyzz_xzzz[k] * ab_z + g_xy_0_xyzz_xzzzz[k];

                g_xy_0_xyzzz_yyyy[k] = -g_xy_0_xyzz_yyyy[k] * ab_z + g_xy_0_xyzz_yyyyz[k];

                g_xy_0_xyzzz_yyyz[k] = -g_xy_0_xyzz_yyyz[k] * ab_z + g_xy_0_xyzz_yyyzz[k];

                g_xy_0_xyzzz_yyzz[k] = -g_xy_0_xyzz_yyzz[k] * ab_z + g_xy_0_xyzz_yyzzz[k];

                g_xy_0_xyzzz_yzzz[k] = -g_xy_0_xyzz_yzzz[k] * ab_z + g_xy_0_xyzz_yzzzz[k];

                g_xy_0_xyzzz_zzzz[k] = -g_xy_0_xyzz_zzzz[k] * ab_z + g_xy_0_xyzz_zzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 525 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 526 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 527 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 528 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 529 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 530 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 531 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 532 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 533 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 534 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 535 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 536 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 537 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 538 * ccomps * dcomps);

            auto g_xy_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzz_xxxx, g_xy_0_xzzz_xxxxz, g_xy_0_xzzz_xxxy, g_xy_0_xzzz_xxxyz, g_xy_0_xzzz_xxxz, g_xy_0_xzzz_xxxzz, g_xy_0_xzzz_xxyy, g_xy_0_xzzz_xxyyz, g_xy_0_xzzz_xxyz, g_xy_0_xzzz_xxyzz, g_xy_0_xzzz_xxzz, g_xy_0_xzzz_xxzzz, g_xy_0_xzzz_xyyy, g_xy_0_xzzz_xyyyz, g_xy_0_xzzz_xyyz, g_xy_0_xzzz_xyyzz, g_xy_0_xzzz_xyzz, g_xy_0_xzzz_xyzzz, g_xy_0_xzzz_xzzz, g_xy_0_xzzz_xzzzz, g_xy_0_xzzz_yyyy, g_xy_0_xzzz_yyyyz, g_xy_0_xzzz_yyyz, g_xy_0_xzzz_yyyzz, g_xy_0_xzzz_yyzz, g_xy_0_xzzz_yyzzz, g_xy_0_xzzz_yzzz, g_xy_0_xzzz_yzzzz, g_xy_0_xzzz_zzzz, g_xy_0_xzzz_zzzzz, g_xy_0_xzzzz_xxxx, g_xy_0_xzzzz_xxxy, g_xy_0_xzzzz_xxxz, g_xy_0_xzzzz_xxyy, g_xy_0_xzzzz_xxyz, g_xy_0_xzzzz_xxzz, g_xy_0_xzzzz_xyyy, g_xy_0_xzzzz_xyyz, g_xy_0_xzzzz_xyzz, g_xy_0_xzzzz_xzzz, g_xy_0_xzzzz_yyyy, g_xy_0_xzzzz_yyyz, g_xy_0_xzzzz_yyzz, g_xy_0_xzzzz_yzzz, g_xy_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzz_xxxx[k] = -g_xy_0_xzzz_xxxx[k] * ab_z + g_xy_0_xzzz_xxxxz[k];

                g_xy_0_xzzzz_xxxy[k] = -g_xy_0_xzzz_xxxy[k] * ab_z + g_xy_0_xzzz_xxxyz[k];

                g_xy_0_xzzzz_xxxz[k] = -g_xy_0_xzzz_xxxz[k] * ab_z + g_xy_0_xzzz_xxxzz[k];

                g_xy_0_xzzzz_xxyy[k] = -g_xy_0_xzzz_xxyy[k] * ab_z + g_xy_0_xzzz_xxyyz[k];

                g_xy_0_xzzzz_xxyz[k] = -g_xy_0_xzzz_xxyz[k] * ab_z + g_xy_0_xzzz_xxyzz[k];

                g_xy_0_xzzzz_xxzz[k] = -g_xy_0_xzzz_xxzz[k] * ab_z + g_xy_0_xzzz_xxzzz[k];

                g_xy_0_xzzzz_xyyy[k] = -g_xy_0_xzzz_xyyy[k] * ab_z + g_xy_0_xzzz_xyyyz[k];

                g_xy_0_xzzzz_xyyz[k] = -g_xy_0_xzzz_xyyz[k] * ab_z + g_xy_0_xzzz_xyyzz[k];

                g_xy_0_xzzzz_xyzz[k] = -g_xy_0_xzzz_xyzz[k] * ab_z + g_xy_0_xzzz_xyzzz[k];

                g_xy_0_xzzzz_xzzz[k] = -g_xy_0_xzzz_xzzz[k] * ab_z + g_xy_0_xzzz_xzzzz[k];

                g_xy_0_xzzzz_yyyy[k] = -g_xy_0_xzzz_yyyy[k] * ab_z + g_xy_0_xzzz_yyyyz[k];

                g_xy_0_xzzzz_yyyz[k] = -g_xy_0_xzzz_yyyz[k] * ab_z + g_xy_0_xzzz_yyyzz[k];

                g_xy_0_xzzzz_yyzz[k] = -g_xy_0_xzzz_yyzz[k] * ab_z + g_xy_0_xzzz_yyzzz[k];

                g_xy_0_xzzzz_yzzz[k] = -g_xy_0_xzzz_yzzz[k] * ab_z + g_xy_0_xzzz_yzzzz[k];

                g_xy_0_xzzzz_zzzz[k] = -g_xy_0_xzzz_zzzz[k] * ab_z + g_xy_0_xzzz_zzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 540 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 541 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 542 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 543 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 544 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 545 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 546 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 547 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 548 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 549 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 550 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 551 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 552 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 553 * ccomps * dcomps);

            auto g_xy_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_xxxx, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxxz, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxzz, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xzzz, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_zzzz, g_xy_0_yyyy_xxxx, g_xy_0_yyyy_xxxxy, g_xy_0_yyyy_xxxy, g_xy_0_yyyy_xxxyy, g_xy_0_yyyy_xxxyz, g_xy_0_yyyy_xxxz, g_xy_0_yyyy_xxyy, g_xy_0_yyyy_xxyyy, g_xy_0_yyyy_xxyyz, g_xy_0_yyyy_xxyz, g_xy_0_yyyy_xxyzz, g_xy_0_yyyy_xxzz, g_xy_0_yyyy_xyyy, g_xy_0_yyyy_xyyyy, g_xy_0_yyyy_xyyyz, g_xy_0_yyyy_xyyz, g_xy_0_yyyy_xyyzz, g_xy_0_yyyy_xyzz, g_xy_0_yyyy_xyzzz, g_xy_0_yyyy_xzzz, g_xy_0_yyyy_yyyy, g_xy_0_yyyy_yyyyy, g_xy_0_yyyy_yyyyz, g_xy_0_yyyy_yyyz, g_xy_0_yyyy_yyyzz, g_xy_0_yyyy_yyzz, g_xy_0_yyyy_yyzzz, g_xy_0_yyyy_yzzz, g_xy_0_yyyy_yzzzz, g_xy_0_yyyy_zzzz, g_xy_0_yyyyy_xxxx, g_xy_0_yyyyy_xxxy, g_xy_0_yyyyy_xxxz, g_xy_0_yyyyy_xxyy, g_xy_0_yyyyy_xxyz, g_xy_0_yyyyy_xxzz, g_xy_0_yyyyy_xyyy, g_xy_0_yyyyy_xyyz, g_xy_0_yyyyy_xyzz, g_xy_0_yyyyy_xzzz, g_xy_0_yyyyy_yyyy, g_xy_0_yyyyy_yyyz, g_xy_0_yyyyy_yyzz, g_xy_0_yyyyy_yzzz, g_xy_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyy_xxxx[k] = -g_x_0_yyyy_xxxx[k] - g_xy_0_yyyy_xxxx[k] * ab_y + g_xy_0_yyyy_xxxxy[k];

                g_xy_0_yyyyy_xxxy[k] = -g_x_0_yyyy_xxxy[k] - g_xy_0_yyyy_xxxy[k] * ab_y + g_xy_0_yyyy_xxxyy[k];

                g_xy_0_yyyyy_xxxz[k] = -g_x_0_yyyy_xxxz[k] - g_xy_0_yyyy_xxxz[k] * ab_y + g_xy_0_yyyy_xxxyz[k];

                g_xy_0_yyyyy_xxyy[k] = -g_x_0_yyyy_xxyy[k] - g_xy_0_yyyy_xxyy[k] * ab_y + g_xy_0_yyyy_xxyyy[k];

                g_xy_0_yyyyy_xxyz[k] = -g_x_0_yyyy_xxyz[k] - g_xy_0_yyyy_xxyz[k] * ab_y + g_xy_0_yyyy_xxyyz[k];

                g_xy_0_yyyyy_xxzz[k] = -g_x_0_yyyy_xxzz[k] - g_xy_0_yyyy_xxzz[k] * ab_y + g_xy_0_yyyy_xxyzz[k];

                g_xy_0_yyyyy_xyyy[k] = -g_x_0_yyyy_xyyy[k] - g_xy_0_yyyy_xyyy[k] * ab_y + g_xy_0_yyyy_xyyyy[k];

                g_xy_0_yyyyy_xyyz[k] = -g_x_0_yyyy_xyyz[k] - g_xy_0_yyyy_xyyz[k] * ab_y + g_xy_0_yyyy_xyyyz[k];

                g_xy_0_yyyyy_xyzz[k] = -g_x_0_yyyy_xyzz[k] - g_xy_0_yyyy_xyzz[k] * ab_y + g_xy_0_yyyy_xyyzz[k];

                g_xy_0_yyyyy_xzzz[k] = -g_x_0_yyyy_xzzz[k] - g_xy_0_yyyy_xzzz[k] * ab_y + g_xy_0_yyyy_xyzzz[k];

                g_xy_0_yyyyy_yyyy[k] = -g_x_0_yyyy_yyyy[k] - g_xy_0_yyyy_yyyy[k] * ab_y + g_xy_0_yyyy_yyyyy[k];

                g_xy_0_yyyyy_yyyz[k] = -g_x_0_yyyy_yyyz[k] - g_xy_0_yyyy_yyyz[k] * ab_y + g_xy_0_yyyy_yyyyz[k];

                g_xy_0_yyyyy_yyzz[k] = -g_x_0_yyyy_yyzz[k] - g_xy_0_yyyy_yyzz[k] * ab_y + g_xy_0_yyyy_yyyzz[k];

                g_xy_0_yyyyy_yzzz[k] = -g_x_0_yyyy_yzzz[k] - g_xy_0_yyyy_yzzz[k] * ab_y + g_xy_0_yyyy_yyzzz[k];

                g_xy_0_yyyyy_zzzz[k] = -g_x_0_yyyy_zzzz[k] - g_xy_0_yyyy_zzzz[k] * ab_y + g_xy_0_yyyy_yzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 555 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 556 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 557 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 558 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 559 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 560 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 561 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 562 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 563 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 564 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 565 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 566 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 567 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 568 * ccomps * dcomps);

            auto g_xy_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyy_xxxx, g_xy_0_yyyy_xxxxz, g_xy_0_yyyy_xxxy, g_xy_0_yyyy_xxxyz, g_xy_0_yyyy_xxxz, g_xy_0_yyyy_xxxzz, g_xy_0_yyyy_xxyy, g_xy_0_yyyy_xxyyz, g_xy_0_yyyy_xxyz, g_xy_0_yyyy_xxyzz, g_xy_0_yyyy_xxzz, g_xy_0_yyyy_xxzzz, g_xy_0_yyyy_xyyy, g_xy_0_yyyy_xyyyz, g_xy_0_yyyy_xyyz, g_xy_0_yyyy_xyyzz, g_xy_0_yyyy_xyzz, g_xy_0_yyyy_xyzzz, g_xy_0_yyyy_xzzz, g_xy_0_yyyy_xzzzz, g_xy_0_yyyy_yyyy, g_xy_0_yyyy_yyyyz, g_xy_0_yyyy_yyyz, g_xy_0_yyyy_yyyzz, g_xy_0_yyyy_yyzz, g_xy_0_yyyy_yyzzz, g_xy_0_yyyy_yzzz, g_xy_0_yyyy_yzzzz, g_xy_0_yyyy_zzzz, g_xy_0_yyyy_zzzzz, g_xy_0_yyyyz_xxxx, g_xy_0_yyyyz_xxxy, g_xy_0_yyyyz_xxxz, g_xy_0_yyyyz_xxyy, g_xy_0_yyyyz_xxyz, g_xy_0_yyyyz_xxzz, g_xy_0_yyyyz_xyyy, g_xy_0_yyyyz_xyyz, g_xy_0_yyyyz_xyzz, g_xy_0_yyyyz_xzzz, g_xy_0_yyyyz_yyyy, g_xy_0_yyyyz_yyyz, g_xy_0_yyyyz_yyzz, g_xy_0_yyyyz_yzzz, g_xy_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyz_xxxx[k] = -g_xy_0_yyyy_xxxx[k] * ab_z + g_xy_0_yyyy_xxxxz[k];

                g_xy_0_yyyyz_xxxy[k] = -g_xy_0_yyyy_xxxy[k] * ab_z + g_xy_0_yyyy_xxxyz[k];

                g_xy_0_yyyyz_xxxz[k] = -g_xy_0_yyyy_xxxz[k] * ab_z + g_xy_0_yyyy_xxxzz[k];

                g_xy_0_yyyyz_xxyy[k] = -g_xy_0_yyyy_xxyy[k] * ab_z + g_xy_0_yyyy_xxyyz[k];

                g_xy_0_yyyyz_xxyz[k] = -g_xy_0_yyyy_xxyz[k] * ab_z + g_xy_0_yyyy_xxyzz[k];

                g_xy_0_yyyyz_xxzz[k] = -g_xy_0_yyyy_xxzz[k] * ab_z + g_xy_0_yyyy_xxzzz[k];

                g_xy_0_yyyyz_xyyy[k] = -g_xy_0_yyyy_xyyy[k] * ab_z + g_xy_0_yyyy_xyyyz[k];

                g_xy_0_yyyyz_xyyz[k] = -g_xy_0_yyyy_xyyz[k] * ab_z + g_xy_0_yyyy_xyyzz[k];

                g_xy_0_yyyyz_xyzz[k] = -g_xy_0_yyyy_xyzz[k] * ab_z + g_xy_0_yyyy_xyzzz[k];

                g_xy_0_yyyyz_xzzz[k] = -g_xy_0_yyyy_xzzz[k] * ab_z + g_xy_0_yyyy_xzzzz[k];

                g_xy_0_yyyyz_yyyy[k] = -g_xy_0_yyyy_yyyy[k] * ab_z + g_xy_0_yyyy_yyyyz[k];

                g_xy_0_yyyyz_yyyz[k] = -g_xy_0_yyyy_yyyz[k] * ab_z + g_xy_0_yyyy_yyyzz[k];

                g_xy_0_yyyyz_yyzz[k] = -g_xy_0_yyyy_yyzz[k] * ab_z + g_xy_0_yyyy_yyzzz[k];

                g_xy_0_yyyyz_yzzz[k] = -g_xy_0_yyyy_yzzz[k] * ab_z + g_xy_0_yyyy_yzzzz[k];

                g_xy_0_yyyyz_zzzz[k] = -g_xy_0_yyyy_zzzz[k] * ab_z + g_xy_0_yyyy_zzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 570 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 571 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 572 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 573 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 574 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 575 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 576 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 577 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 578 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 579 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 580 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 581 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 582 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 583 * ccomps * dcomps);

            auto g_xy_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyz_xxxx, g_xy_0_yyyz_xxxxz, g_xy_0_yyyz_xxxy, g_xy_0_yyyz_xxxyz, g_xy_0_yyyz_xxxz, g_xy_0_yyyz_xxxzz, g_xy_0_yyyz_xxyy, g_xy_0_yyyz_xxyyz, g_xy_0_yyyz_xxyz, g_xy_0_yyyz_xxyzz, g_xy_0_yyyz_xxzz, g_xy_0_yyyz_xxzzz, g_xy_0_yyyz_xyyy, g_xy_0_yyyz_xyyyz, g_xy_0_yyyz_xyyz, g_xy_0_yyyz_xyyzz, g_xy_0_yyyz_xyzz, g_xy_0_yyyz_xyzzz, g_xy_0_yyyz_xzzz, g_xy_0_yyyz_xzzzz, g_xy_0_yyyz_yyyy, g_xy_0_yyyz_yyyyz, g_xy_0_yyyz_yyyz, g_xy_0_yyyz_yyyzz, g_xy_0_yyyz_yyzz, g_xy_0_yyyz_yyzzz, g_xy_0_yyyz_yzzz, g_xy_0_yyyz_yzzzz, g_xy_0_yyyz_zzzz, g_xy_0_yyyz_zzzzz, g_xy_0_yyyzz_xxxx, g_xy_0_yyyzz_xxxy, g_xy_0_yyyzz_xxxz, g_xy_0_yyyzz_xxyy, g_xy_0_yyyzz_xxyz, g_xy_0_yyyzz_xxzz, g_xy_0_yyyzz_xyyy, g_xy_0_yyyzz_xyyz, g_xy_0_yyyzz_xyzz, g_xy_0_yyyzz_xzzz, g_xy_0_yyyzz_yyyy, g_xy_0_yyyzz_yyyz, g_xy_0_yyyzz_yyzz, g_xy_0_yyyzz_yzzz, g_xy_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzz_xxxx[k] = -g_xy_0_yyyz_xxxx[k] * ab_z + g_xy_0_yyyz_xxxxz[k];

                g_xy_0_yyyzz_xxxy[k] = -g_xy_0_yyyz_xxxy[k] * ab_z + g_xy_0_yyyz_xxxyz[k];

                g_xy_0_yyyzz_xxxz[k] = -g_xy_0_yyyz_xxxz[k] * ab_z + g_xy_0_yyyz_xxxzz[k];

                g_xy_0_yyyzz_xxyy[k] = -g_xy_0_yyyz_xxyy[k] * ab_z + g_xy_0_yyyz_xxyyz[k];

                g_xy_0_yyyzz_xxyz[k] = -g_xy_0_yyyz_xxyz[k] * ab_z + g_xy_0_yyyz_xxyzz[k];

                g_xy_0_yyyzz_xxzz[k] = -g_xy_0_yyyz_xxzz[k] * ab_z + g_xy_0_yyyz_xxzzz[k];

                g_xy_0_yyyzz_xyyy[k] = -g_xy_0_yyyz_xyyy[k] * ab_z + g_xy_0_yyyz_xyyyz[k];

                g_xy_0_yyyzz_xyyz[k] = -g_xy_0_yyyz_xyyz[k] * ab_z + g_xy_0_yyyz_xyyzz[k];

                g_xy_0_yyyzz_xyzz[k] = -g_xy_0_yyyz_xyzz[k] * ab_z + g_xy_0_yyyz_xyzzz[k];

                g_xy_0_yyyzz_xzzz[k] = -g_xy_0_yyyz_xzzz[k] * ab_z + g_xy_0_yyyz_xzzzz[k];

                g_xy_0_yyyzz_yyyy[k] = -g_xy_0_yyyz_yyyy[k] * ab_z + g_xy_0_yyyz_yyyyz[k];

                g_xy_0_yyyzz_yyyz[k] = -g_xy_0_yyyz_yyyz[k] * ab_z + g_xy_0_yyyz_yyyzz[k];

                g_xy_0_yyyzz_yyzz[k] = -g_xy_0_yyyz_yyzz[k] * ab_z + g_xy_0_yyyz_yyzzz[k];

                g_xy_0_yyyzz_yzzz[k] = -g_xy_0_yyyz_yzzz[k] * ab_z + g_xy_0_yyyz_yzzzz[k];

                g_xy_0_yyyzz_zzzz[k] = -g_xy_0_yyyz_zzzz[k] * ab_z + g_xy_0_yyyz_zzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 585 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 586 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 587 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 588 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 589 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 590 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 591 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 592 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 593 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 594 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 595 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 596 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 597 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 598 * ccomps * dcomps);

            auto g_xy_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzz_xxxx, g_xy_0_yyzz_xxxxz, g_xy_0_yyzz_xxxy, g_xy_0_yyzz_xxxyz, g_xy_0_yyzz_xxxz, g_xy_0_yyzz_xxxzz, g_xy_0_yyzz_xxyy, g_xy_0_yyzz_xxyyz, g_xy_0_yyzz_xxyz, g_xy_0_yyzz_xxyzz, g_xy_0_yyzz_xxzz, g_xy_0_yyzz_xxzzz, g_xy_0_yyzz_xyyy, g_xy_0_yyzz_xyyyz, g_xy_0_yyzz_xyyz, g_xy_0_yyzz_xyyzz, g_xy_0_yyzz_xyzz, g_xy_0_yyzz_xyzzz, g_xy_0_yyzz_xzzz, g_xy_0_yyzz_xzzzz, g_xy_0_yyzz_yyyy, g_xy_0_yyzz_yyyyz, g_xy_0_yyzz_yyyz, g_xy_0_yyzz_yyyzz, g_xy_0_yyzz_yyzz, g_xy_0_yyzz_yyzzz, g_xy_0_yyzz_yzzz, g_xy_0_yyzz_yzzzz, g_xy_0_yyzz_zzzz, g_xy_0_yyzz_zzzzz, g_xy_0_yyzzz_xxxx, g_xy_0_yyzzz_xxxy, g_xy_0_yyzzz_xxxz, g_xy_0_yyzzz_xxyy, g_xy_0_yyzzz_xxyz, g_xy_0_yyzzz_xxzz, g_xy_0_yyzzz_xyyy, g_xy_0_yyzzz_xyyz, g_xy_0_yyzzz_xyzz, g_xy_0_yyzzz_xzzz, g_xy_0_yyzzz_yyyy, g_xy_0_yyzzz_yyyz, g_xy_0_yyzzz_yyzz, g_xy_0_yyzzz_yzzz, g_xy_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzz_xxxx[k] = -g_xy_0_yyzz_xxxx[k] * ab_z + g_xy_0_yyzz_xxxxz[k];

                g_xy_0_yyzzz_xxxy[k] = -g_xy_0_yyzz_xxxy[k] * ab_z + g_xy_0_yyzz_xxxyz[k];

                g_xy_0_yyzzz_xxxz[k] = -g_xy_0_yyzz_xxxz[k] * ab_z + g_xy_0_yyzz_xxxzz[k];

                g_xy_0_yyzzz_xxyy[k] = -g_xy_0_yyzz_xxyy[k] * ab_z + g_xy_0_yyzz_xxyyz[k];

                g_xy_0_yyzzz_xxyz[k] = -g_xy_0_yyzz_xxyz[k] * ab_z + g_xy_0_yyzz_xxyzz[k];

                g_xy_0_yyzzz_xxzz[k] = -g_xy_0_yyzz_xxzz[k] * ab_z + g_xy_0_yyzz_xxzzz[k];

                g_xy_0_yyzzz_xyyy[k] = -g_xy_0_yyzz_xyyy[k] * ab_z + g_xy_0_yyzz_xyyyz[k];

                g_xy_0_yyzzz_xyyz[k] = -g_xy_0_yyzz_xyyz[k] * ab_z + g_xy_0_yyzz_xyyzz[k];

                g_xy_0_yyzzz_xyzz[k] = -g_xy_0_yyzz_xyzz[k] * ab_z + g_xy_0_yyzz_xyzzz[k];

                g_xy_0_yyzzz_xzzz[k] = -g_xy_0_yyzz_xzzz[k] * ab_z + g_xy_0_yyzz_xzzzz[k];

                g_xy_0_yyzzz_yyyy[k] = -g_xy_0_yyzz_yyyy[k] * ab_z + g_xy_0_yyzz_yyyyz[k];

                g_xy_0_yyzzz_yyyz[k] = -g_xy_0_yyzz_yyyz[k] * ab_z + g_xy_0_yyzz_yyyzz[k];

                g_xy_0_yyzzz_yyzz[k] = -g_xy_0_yyzz_yyzz[k] * ab_z + g_xy_0_yyzz_yyzzz[k];

                g_xy_0_yyzzz_yzzz[k] = -g_xy_0_yyzz_yzzz[k] * ab_z + g_xy_0_yyzz_yzzzz[k];

                g_xy_0_yyzzz_zzzz[k] = -g_xy_0_yyzz_zzzz[k] * ab_z + g_xy_0_yyzz_zzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 600 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 601 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 602 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 603 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 604 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 605 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 606 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 607 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 608 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 609 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 610 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 611 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 612 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 613 * ccomps * dcomps);

            auto g_xy_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzz_xxxx, g_xy_0_yzzz_xxxxz, g_xy_0_yzzz_xxxy, g_xy_0_yzzz_xxxyz, g_xy_0_yzzz_xxxz, g_xy_0_yzzz_xxxzz, g_xy_0_yzzz_xxyy, g_xy_0_yzzz_xxyyz, g_xy_0_yzzz_xxyz, g_xy_0_yzzz_xxyzz, g_xy_0_yzzz_xxzz, g_xy_0_yzzz_xxzzz, g_xy_0_yzzz_xyyy, g_xy_0_yzzz_xyyyz, g_xy_0_yzzz_xyyz, g_xy_0_yzzz_xyyzz, g_xy_0_yzzz_xyzz, g_xy_0_yzzz_xyzzz, g_xy_0_yzzz_xzzz, g_xy_0_yzzz_xzzzz, g_xy_0_yzzz_yyyy, g_xy_0_yzzz_yyyyz, g_xy_0_yzzz_yyyz, g_xy_0_yzzz_yyyzz, g_xy_0_yzzz_yyzz, g_xy_0_yzzz_yyzzz, g_xy_0_yzzz_yzzz, g_xy_0_yzzz_yzzzz, g_xy_0_yzzz_zzzz, g_xy_0_yzzz_zzzzz, g_xy_0_yzzzz_xxxx, g_xy_0_yzzzz_xxxy, g_xy_0_yzzzz_xxxz, g_xy_0_yzzzz_xxyy, g_xy_0_yzzzz_xxyz, g_xy_0_yzzzz_xxzz, g_xy_0_yzzzz_xyyy, g_xy_0_yzzzz_xyyz, g_xy_0_yzzzz_xyzz, g_xy_0_yzzzz_xzzz, g_xy_0_yzzzz_yyyy, g_xy_0_yzzzz_yyyz, g_xy_0_yzzzz_yyzz, g_xy_0_yzzzz_yzzz, g_xy_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzz_xxxx[k] = -g_xy_0_yzzz_xxxx[k] * ab_z + g_xy_0_yzzz_xxxxz[k];

                g_xy_0_yzzzz_xxxy[k] = -g_xy_0_yzzz_xxxy[k] * ab_z + g_xy_0_yzzz_xxxyz[k];

                g_xy_0_yzzzz_xxxz[k] = -g_xy_0_yzzz_xxxz[k] * ab_z + g_xy_0_yzzz_xxxzz[k];

                g_xy_0_yzzzz_xxyy[k] = -g_xy_0_yzzz_xxyy[k] * ab_z + g_xy_0_yzzz_xxyyz[k];

                g_xy_0_yzzzz_xxyz[k] = -g_xy_0_yzzz_xxyz[k] * ab_z + g_xy_0_yzzz_xxyzz[k];

                g_xy_0_yzzzz_xxzz[k] = -g_xy_0_yzzz_xxzz[k] * ab_z + g_xy_0_yzzz_xxzzz[k];

                g_xy_0_yzzzz_xyyy[k] = -g_xy_0_yzzz_xyyy[k] * ab_z + g_xy_0_yzzz_xyyyz[k];

                g_xy_0_yzzzz_xyyz[k] = -g_xy_0_yzzz_xyyz[k] * ab_z + g_xy_0_yzzz_xyyzz[k];

                g_xy_0_yzzzz_xyzz[k] = -g_xy_0_yzzz_xyzz[k] * ab_z + g_xy_0_yzzz_xyzzz[k];

                g_xy_0_yzzzz_xzzz[k] = -g_xy_0_yzzz_xzzz[k] * ab_z + g_xy_0_yzzz_xzzzz[k];

                g_xy_0_yzzzz_yyyy[k] = -g_xy_0_yzzz_yyyy[k] * ab_z + g_xy_0_yzzz_yyyyz[k];

                g_xy_0_yzzzz_yyyz[k] = -g_xy_0_yzzz_yyyz[k] * ab_z + g_xy_0_yzzz_yyyzz[k];

                g_xy_0_yzzzz_yyzz[k] = -g_xy_0_yzzz_yyzz[k] * ab_z + g_xy_0_yzzz_yyzzz[k];

                g_xy_0_yzzzz_yzzz[k] = -g_xy_0_yzzz_yzzz[k] * ab_z + g_xy_0_yzzz_yzzzz[k];

                g_xy_0_yzzzz_zzzz[k] = -g_xy_0_yzzz_zzzz[k] * ab_z + g_xy_0_yzzz_zzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 615 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 616 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 617 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 618 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 619 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 620 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 621 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 622 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 623 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 624 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 625 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 626 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 627 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 628 * ccomps * dcomps);

            auto g_xy_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzz_xxxx, g_xy_0_zzzz_xxxxz, g_xy_0_zzzz_xxxy, g_xy_0_zzzz_xxxyz, g_xy_0_zzzz_xxxz, g_xy_0_zzzz_xxxzz, g_xy_0_zzzz_xxyy, g_xy_0_zzzz_xxyyz, g_xy_0_zzzz_xxyz, g_xy_0_zzzz_xxyzz, g_xy_0_zzzz_xxzz, g_xy_0_zzzz_xxzzz, g_xy_0_zzzz_xyyy, g_xy_0_zzzz_xyyyz, g_xy_0_zzzz_xyyz, g_xy_0_zzzz_xyyzz, g_xy_0_zzzz_xyzz, g_xy_0_zzzz_xyzzz, g_xy_0_zzzz_xzzz, g_xy_0_zzzz_xzzzz, g_xy_0_zzzz_yyyy, g_xy_0_zzzz_yyyyz, g_xy_0_zzzz_yyyz, g_xy_0_zzzz_yyyzz, g_xy_0_zzzz_yyzz, g_xy_0_zzzz_yyzzz, g_xy_0_zzzz_yzzz, g_xy_0_zzzz_yzzzz, g_xy_0_zzzz_zzzz, g_xy_0_zzzz_zzzzz, g_xy_0_zzzzz_xxxx, g_xy_0_zzzzz_xxxy, g_xy_0_zzzzz_xxxz, g_xy_0_zzzzz_xxyy, g_xy_0_zzzzz_xxyz, g_xy_0_zzzzz_xxzz, g_xy_0_zzzzz_xyyy, g_xy_0_zzzzz_xyyz, g_xy_0_zzzzz_xyzz, g_xy_0_zzzzz_xzzz, g_xy_0_zzzzz_yyyy, g_xy_0_zzzzz_yyyz, g_xy_0_zzzzz_yyzz, g_xy_0_zzzzz_yzzz, g_xy_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzz_xxxx[k] = -g_xy_0_zzzz_xxxx[k] * ab_z + g_xy_0_zzzz_xxxxz[k];

                g_xy_0_zzzzz_xxxy[k] = -g_xy_0_zzzz_xxxy[k] * ab_z + g_xy_0_zzzz_xxxyz[k];

                g_xy_0_zzzzz_xxxz[k] = -g_xy_0_zzzz_xxxz[k] * ab_z + g_xy_0_zzzz_xxxzz[k];

                g_xy_0_zzzzz_xxyy[k] = -g_xy_0_zzzz_xxyy[k] * ab_z + g_xy_0_zzzz_xxyyz[k];

                g_xy_0_zzzzz_xxyz[k] = -g_xy_0_zzzz_xxyz[k] * ab_z + g_xy_0_zzzz_xxyzz[k];

                g_xy_0_zzzzz_xxzz[k] = -g_xy_0_zzzz_xxzz[k] * ab_z + g_xy_0_zzzz_xxzzz[k];

                g_xy_0_zzzzz_xyyy[k] = -g_xy_0_zzzz_xyyy[k] * ab_z + g_xy_0_zzzz_xyyyz[k];

                g_xy_0_zzzzz_xyyz[k] = -g_xy_0_zzzz_xyyz[k] * ab_z + g_xy_0_zzzz_xyyzz[k];

                g_xy_0_zzzzz_xyzz[k] = -g_xy_0_zzzz_xyzz[k] * ab_z + g_xy_0_zzzz_xyzzz[k];

                g_xy_0_zzzzz_xzzz[k] = -g_xy_0_zzzz_xzzz[k] * ab_z + g_xy_0_zzzz_xzzzz[k];

                g_xy_0_zzzzz_yyyy[k] = -g_xy_0_zzzz_yyyy[k] * ab_z + g_xy_0_zzzz_yyyyz[k];

                g_xy_0_zzzzz_yyyz[k] = -g_xy_0_zzzz_yyyz[k] * ab_z + g_xy_0_zzzz_yyyzz[k];

                g_xy_0_zzzzz_yyzz[k] = -g_xy_0_zzzz_yyzz[k] * ab_z + g_xy_0_zzzz_yyzzz[k];

                g_xy_0_zzzzz_yzzz[k] = -g_xy_0_zzzz_yzzz[k] * ab_z + g_xy_0_zzzz_yzzzz[k];

                g_xy_0_zzzzz_zzzz[k] = -g_xy_0_zzzz_zzzz[k] * ab_z + g_xy_0_zzzz_zzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 630 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 631 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 632 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 633 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 634 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 635 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 636 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 637 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 638 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 639 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 640 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 641 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 642 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 643 * ccomps * dcomps);

            auto g_xz_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_xxxx, g_xz_0_xxxx_xxxxx, g_xz_0_xxxx_xxxxy, g_xz_0_xxxx_xxxxz, g_xz_0_xxxx_xxxy, g_xz_0_xxxx_xxxyy, g_xz_0_xxxx_xxxyz, g_xz_0_xxxx_xxxz, g_xz_0_xxxx_xxxzz, g_xz_0_xxxx_xxyy, g_xz_0_xxxx_xxyyy, g_xz_0_xxxx_xxyyz, g_xz_0_xxxx_xxyz, g_xz_0_xxxx_xxyzz, g_xz_0_xxxx_xxzz, g_xz_0_xxxx_xxzzz, g_xz_0_xxxx_xyyy, g_xz_0_xxxx_xyyyy, g_xz_0_xxxx_xyyyz, g_xz_0_xxxx_xyyz, g_xz_0_xxxx_xyyzz, g_xz_0_xxxx_xyzz, g_xz_0_xxxx_xyzzz, g_xz_0_xxxx_xzzz, g_xz_0_xxxx_xzzzz, g_xz_0_xxxx_yyyy, g_xz_0_xxxx_yyyz, g_xz_0_xxxx_yyzz, g_xz_0_xxxx_yzzz, g_xz_0_xxxx_zzzz, g_xz_0_xxxxx_xxxx, g_xz_0_xxxxx_xxxy, g_xz_0_xxxxx_xxxz, g_xz_0_xxxxx_xxyy, g_xz_0_xxxxx_xxyz, g_xz_0_xxxxx_xxzz, g_xz_0_xxxxx_xyyy, g_xz_0_xxxxx_xyyz, g_xz_0_xxxxx_xyzz, g_xz_0_xxxxx_xzzz, g_xz_0_xxxxx_yyyy, g_xz_0_xxxxx_yyyz, g_xz_0_xxxxx_yyzz, g_xz_0_xxxxx_yzzz, g_xz_0_xxxxx_zzzz, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_yyyy, g_z_0_xxxx_yyyz, g_z_0_xxxx_yyzz, g_z_0_xxxx_yzzz, g_z_0_xxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxx_xxxx[k] = -g_z_0_xxxx_xxxx[k] - g_xz_0_xxxx_xxxx[k] * ab_x + g_xz_0_xxxx_xxxxx[k];

                g_xz_0_xxxxx_xxxy[k] = -g_z_0_xxxx_xxxy[k] - g_xz_0_xxxx_xxxy[k] * ab_x + g_xz_0_xxxx_xxxxy[k];

                g_xz_0_xxxxx_xxxz[k] = -g_z_0_xxxx_xxxz[k] - g_xz_0_xxxx_xxxz[k] * ab_x + g_xz_0_xxxx_xxxxz[k];

                g_xz_0_xxxxx_xxyy[k] = -g_z_0_xxxx_xxyy[k] - g_xz_0_xxxx_xxyy[k] * ab_x + g_xz_0_xxxx_xxxyy[k];

                g_xz_0_xxxxx_xxyz[k] = -g_z_0_xxxx_xxyz[k] - g_xz_0_xxxx_xxyz[k] * ab_x + g_xz_0_xxxx_xxxyz[k];

                g_xz_0_xxxxx_xxzz[k] = -g_z_0_xxxx_xxzz[k] - g_xz_0_xxxx_xxzz[k] * ab_x + g_xz_0_xxxx_xxxzz[k];

                g_xz_0_xxxxx_xyyy[k] = -g_z_0_xxxx_xyyy[k] - g_xz_0_xxxx_xyyy[k] * ab_x + g_xz_0_xxxx_xxyyy[k];

                g_xz_0_xxxxx_xyyz[k] = -g_z_0_xxxx_xyyz[k] - g_xz_0_xxxx_xyyz[k] * ab_x + g_xz_0_xxxx_xxyyz[k];

                g_xz_0_xxxxx_xyzz[k] = -g_z_0_xxxx_xyzz[k] - g_xz_0_xxxx_xyzz[k] * ab_x + g_xz_0_xxxx_xxyzz[k];

                g_xz_0_xxxxx_xzzz[k] = -g_z_0_xxxx_xzzz[k] - g_xz_0_xxxx_xzzz[k] * ab_x + g_xz_0_xxxx_xxzzz[k];

                g_xz_0_xxxxx_yyyy[k] = -g_z_0_xxxx_yyyy[k] - g_xz_0_xxxx_yyyy[k] * ab_x + g_xz_0_xxxx_xyyyy[k];

                g_xz_0_xxxxx_yyyz[k] = -g_z_0_xxxx_yyyz[k] - g_xz_0_xxxx_yyyz[k] * ab_x + g_xz_0_xxxx_xyyyz[k];

                g_xz_0_xxxxx_yyzz[k] = -g_z_0_xxxx_yyzz[k] - g_xz_0_xxxx_yyzz[k] * ab_x + g_xz_0_xxxx_xyyzz[k];

                g_xz_0_xxxxx_yzzz[k] = -g_z_0_xxxx_yzzz[k] - g_xz_0_xxxx_yzzz[k] * ab_x + g_xz_0_xxxx_xyzzz[k];

                g_xz_0_xxxxx_zzzz[k] = -g_z_0_xxxx_zzzz[k] - g_xz_0_xxxx_zzzz[k] * ab_x + g_xz_0_xxxx_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 645 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 646 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 647 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 648 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 649 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 650 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 651 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 652 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 653 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 654 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 655 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 656 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 657 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 658 * ccomps * dcomps);

            auto g_xz_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_xxxx, g_xz_0_xxxx_xxxxy, g_xz_0_xxxx_xxxy, g_xz_0_xxxx_xxxyy, g_xz_0_xxxx_xxxyz, g_xz_0_xxxx_xxxz, g_xz_0_xxxx_xxyy, g_xz_0_xxxx_xxyyy, g_xz_0_xxxx_xxyyz, g_xz_0_xxxx_xxyz, g_xz_0_xxxx_xxyzz, g_xz_0_xxxx_xxzz, g_xz_0_xxxx_xyyy, g_xz_0_xxxx_xyyyy, g_xz_0_xxxx_xyyyz, g_xz_0_xxxx_xyyz, g_xz_0_xxxx_xyyzz, g_xz_0_xxxx_xyzz, g_xz_0_xxxx_xyzzz, g_xz_0_xxxx_xzzz, g_xz_0_xxxx_yyyy, g_xz_0_xxxx_yyyyy, g_xz_0_xxxx_yyyyz, g_xz_0_xxxx_yyyz, g_xz_0_xxxx_yyyzz, g_xz_0_xxxx_yyzz, g_xz_0_xxxx_yyzzz, g_xz_0_xxxx_yzzz, g_xz_0_xxxx_yzzzz, g_xz_0_xxxx_zzzz, g_xz_0_xxxxy_xxxx, g_xz_0_xxxxy_xxxy, g_xz_0_xxxxy_xxxz, g_xz_0_xxxxy_xxyy, g_xz_0_xxxxy_xxyz, g_xz_0_xxxxy_xxzz, g_xz_0_xxxxy_xyyy, g_xz_0_xxxxy_xyyz, g_xz_0_xxxxy_xyzz, g_xz_0_xxxxy_xzzz, g_xz_0_xxxxy_yyyy, g_xz_0_xxxxy_yyyz, g_xz_0_xxxxy_yyzz, g_xz_0_xxxxy_yzzz, g_xz_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxy_xxxx[k] = -g_xz_0_xxxx_xxxx[k] * ab_y + g_xz_0_xxxx_xxxxy[k];

                g_xz_0_xxxxy_xxxy[k] = -g_xz_0_xxxx_xxxy[k] * ab_y + g_xz_0_xxxx_xxxyy[k];

                g_xz_0_xxxxy_xxxz[k] = -g_xz_0_xxxx_xxxz[k] * ab_y + g_xz_0_xxxx_xxxyz[k];

                g_xz_0_xxxxy_xxyy[k] = -g_xz_0_xxxx_xxyy[k] * ab_y + g_xz_0_xxxx_xxyyy[k];

                g_xz_0_xxxxy_xxyz[k] = -g_xz_0_xxxx_xxyz[k] * ab_y + g_xz_0_xxxx_xxyyz[k];

                g_xz_0_xxxxy_xxzz[k] = -g_xz_0_xxxx_xxzz[k] * ab_y + g_xz_0_xxxx_xxyzz[k];

                g_xz_0_xxxxy_xyyy[k] = -g_xz_0_xxxx_xyyy[k] * ab_y + g_xz_0_xxxx_xyyyy[k];

                g_xz_0_xxxxy_xyyz[k] = -g_xz_0_xxxx_xyyz[k] * ab_y + g_xz_0_xxxx_xyyyz[k];

                g_xz_0_xxxxy_xyzz[k] = -g_xz_0_xxxx_xyzz[k] * ab_y + g_xz_0_xxxx_xyyzz[k];

                g_xz_0_xxxxy_xzzz[k] = -g_xz_0_xxxx_xzzz[k] * ab_y + g_xz_0_xxxx_xyzzz[k];

                g_xz_0_xxxxy_yyyy[k] = -g_xz_0_xxxx_yyyy[k] * ab_y + g_xz_0_xxxx_yyyyy[k];

                g_xz_0_xxxxy_yyyz[k] = -g_xz_0_xxxx_yyyz[k] * ab_y + g_xz_0_xxxx_yyyyz[k];

                g_xz_0_xxxxy_yyzz[k] = -g_xz_0_xxxx_yyzz[k] * ab_y + g_xz_0_xxxx_yyyzz[k];

                g_xz_0_xxxxy_yzzz[k] = -g_xz_0_xxxx_yzzz[k] * ab_y + g_xz_0_xxxx_yyzzz[k];

                g_xz_0_xxxxy_zzzz[k] = -g_xz_0_xxxx_zzzz[k] * ab_y + g_xz_0_xxxx_yzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 660 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 661 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 662 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 663 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 664 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 665 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 666 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 667 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 668 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 669 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 670 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 671 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 672 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 673 * ccomps * dcomps);

            auto g_xz_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxz_xxxx, g_xz_0_xxxxz_xxxy, g_xz_0_xxxxz_xxxz, g_xz_0_xxxxz_xxyy, g_xz_0_xxxxz_xxyz, g_xz_0_xxxxz_xxzz, g_xz_0_xxxxz_xyyy, g_xz_0_xxxxz_xyyz, g_xz_0_xxxxz_xyzz, g_xz_0_xxxxz_xzzz, g_xz_0_xxxxz_yyyy, g_xz_0_xxxxz_yyyz, g_xz_0_xxxxz_yyzz, g_xz_0_xxxxz_yzzz, g_xz_0_xxxxz_zzzz, g_xz_0_xxxz_xxxx, g_xz_0_xxxz_xxxxx, g_xz_0_xxxz_xxxxy, g_xz_0_xxxz_xxxxz, g_xz_0_xxxz_xxxy, g_xz_0_xxxz_xxxyy, g_xz_0_xxxz_xxxyz, g_xz_0_xxxz_xxxz, g_xz_0_xxxz_xxxzz, g_xz_0_xxxz_xxyy, g_xz_0_xxxz_xxyyy, g_xz_0_xxxz_xxyyz, g_xz_0_xxxz_xxyz, g_xz_0_xxxz_xxyzz, g_xz_0_xxxz_xxzz, g_xz_0_xxxz_xxzzz, g_xz_0_xxxz_xyyy, g_xz_0_xxxz_xyyyy, g_xz_0_xxxz_xyyyz, g_xz_0_xxxz_xyyz, g_xz_0_xxxz_xyyzz, g_xz_0_xxxz_xyzz, g_xz_0_xxxz_xyzzz, g_xz_0_xxxz_xzzz, g_xz_0_xxxz_xzzzz, g_xz_0_xxxz_yyyy, g_xz_0_xxxz_yyyz, g_xz_0_xxxz_yyzz, g_xz_0_xxxz_yzzz, g_xz_0_xxxz_zzzz, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_yyyy, g_z_0_xxxz_yyyz, g_z_0_xxxz_yyzz, g_z_0_xxxz_yzzz, g_z_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxz_xxxx[k] = -g_z_0_xxxz_xxxx[k] - g_xz_0_xxxz_xxxx[k] * ab_x + g_xz_0_xxxz_xxxxx[k];

                g_xz_0_xxxxz_xxxy[k] = -g_z_0_xxxz_xxxy[k] - g_xz_0_xxxz_xxxy[k] * ab_x + g_xz_0_xxxz_xxxxy[k];

                g_xz_0_xxxxz_xxxz[k] = -g_z_0_xxxz_xxxz[k] - g_xz_0_xxxz_xxxz[k] * ab_x + g_xz_0_xxxz_xxxxz[k];

                g_xz_0_xxxxz_xxyy[k] = -g_z_0_xxxz_xxyy[k] - g_xz_0_xxxz_xxyy[k] * ab_x + g_xz_0_xxxz_xxxyy[k];

                g_xz_0_xxxxz_xxyz[k] = -g_z_0_xxxz_xxyz[k] - g_xz_0_xxxz_xxyz[k] * ab_x + g_xz_0_xxxz_xxxyz[k];

                g_xz_0_xxxxz_xxzz[k] = -g_z_0_xxxz_xxzz[k] - g_xz_0_xxxz_xxzz[k] * ab_x + g_xz_0_xxxz_xxxzz[k];

                g_xz_0_xxxxz_xyyy[k] = -g_z_0_xxxz_xyyy[k] - g_xz_0_xxxz_xyyy[k] * ab_x + g_xz_0_xxxz_xxyyy[k];

                g_xz_0_xxxxz_xyyz[k] = -g_z_0_xxxz_xyyz[k] - g_xz_0_xxxz_xyyz[k] * ab_x + g_xz_0_xxxz_xxyyz[k];

                g_xz_0_xxxxz_xyzz[k] = -g_z_0_xxxz_xyzz[k] - g_xz_0_xxxz_xyzz[k] * ab_x + g_xz_0_xxxz_xxyzz[k];

                g_xz_0_xxxxz_xzzz[k] = -g_z_0_xxxz_xzzz[k] - g_xz_0_xxxz_xzzz[k] * ab_x + g_xz_0_xxxz_xxzzz[k];

                g_xz_0_xxxxz_yyyy[k] = -g_z_0_xxxz_yyyy[k] - g_xz_0_xxxz_yyyy[k] * ab_x + g_xz_0_xxxz_xyyyy[k];

                g_xz_0_xxxxz_yyyz[k] = -g_z_0_xxxz_yyyz[k] - g_xz_0_xxxz_yyyz[k] * ab_x + g_xz_0_xxxz_xyyyz[k];

                g_xz_0_xxxxz_yyzz[k] = -g_z_0_xxxz_yyzz[k] - g_xz_0_xxxz_yyzz[k] * ab_x + g_xz_0_xxxz_xyyzz[k];

                g_xz_0_xxxxz_yzzz[k] = -g_z_0_xxxz_yzzz[k] - g_xz_0_xxxz_yzzz[k] * ab_x + g_xz_0_xxxz_xyzzz[k];

                g_xz_0_xxxxz_zzzz[k] = -g_z_0_xxxz_zzzz[k] - g_xz_0_xxxz_zzzz[k] * ab_x + g_xz_0_xxxz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 675 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 676 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 677 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 678 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 679 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 680 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 681 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 682 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 683 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 684 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 685 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 686 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 687 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 688 * ccomps * dcomps);

            auto g_xz_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxy_xxxx, g_xz_0_xxxy_xxxxy, g_xz_0_xxxy_xxxy, g_xz_0_xxxy_xxxyy, g_xz_0_xxxy_xxxyz, g_xz_0_xxxy_xxxz, g_xz_0_xxxy_xxyy, g_xz_0_xxxy_xxyyy, g_xz_0_xxxy_xxyyz, g_xz_0_xxxy_xxyz, g_xz_0_xxxy_xxyzz, g_xz_0_xxxy_xxzz, g_xz_0_xxxy_xyyy, g_xz_0_xxxy_xyyyy, g_xz_0_xxxy_xyyyz, g_xz_0_xxxy_xyyz, g_xz_0_xxxy_xyyzz, g_xz_0_xxxy_xyzz, g_xz_0_xxxy_xyzzz, g_xz_0_xxxy_xzzz, g_xz_0_xxxy_yyyy, g_xz_0_xxxy_yyyyy, g_xz_0_xxxy_yyyyz, g_xz_0_xxxy_yyyz, g_xz_0_xxxy_yyyzz, g_xz_0_xxxy_yyzz, g_xz_0_xxxy_yyzzz, g_xz_0_xxxy_yzzz, g_xz_0_xxxy_yzzzz, g_xz_0_xxxy_zzzz, g_xz_0_xxxyy_xxxx, g_xz_0_xxxyy_xxxy, g_xz_0_xxxyy_xxxz, g_xz_0_xxxyy_xxyy, g_xz_0_xxxyy_xxyz, g_xz_0_xxxyy_xxzz, g_xz_0_xxxyy_xyyy, g_xz_0_xxxyy_xyyz, g_xz_0_xxxyy_xyzz, g_xz_0_xxxyy_xzzz, g_xz_0_xxxyy_yyyy, g_xz_0_xxxyy_yyyz, g_xz_0_xxxyy_yyzz, g_xz_0_xxxyy_yzzz, g_xz_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyy_xxxx[k] = -g_xz_0_xxxy_xxxx[k] * ab_y + g_xz_0_xxxy_xxxxy[k];

                g_xz_0_xxxyy_xxxy[k] = -g_xz_0_xxxy_xxxy[k] * ab_y + g_xz_0_xxxy_xxxyy[k];

                g_xz_0_xxxyy_xxxz[k] = -g_xz_0_xxxy_xxxz[k] * ab_y + g_xz_0_xxxy_xxxyz[k];

                g_xz_0_xxxyy_xxyy[k] = -g_xz_0_xxxy_xxyy[k] * ab_y + g_xz_0_xxxy_xxyyy[k];

                g_xz_0_xxxyy_xxyz[k] = -g_xz_0_xxxy_xxyz[k] * ab_y + g_xz_0_xxxy_xxyyz[k];

                g_xz_0_xxxyy_xxzz[k] = -g_xz_0_xxxy_xxzz[k] * ab_y + g_xz_0_xxxy_xxyzz[k];

                g_xz_0_xxxyy_xyyy[k] = -g_xz_0_xxxy_xyyy[k] * ab_y + g_xz_0_xxxy_xyyyy[k];

                g_xz_0_xxxyy_xyyz[k] = -g_xz_0_xxxy_xyyz[k] * ab_y + g_xz_0_xxxy_xyyyz[k];

                g_xz_0_xxxyy_xyzz[k] = -g_xz_0_xxxy_xyzz[k] * ab_y + g_xz_0_xxxy_xyyzz[k];

                g_xz_0_xxxyy_xzzz[k] = -g_xz_0_xxxy_xzzz[k] * ab_y + g_xz_0_xxxy_xyzzz[k];

                g_xz_0_xxxyy_yyyy[k] = -g_xz_0_xxxy_yyyy[k] * ab_y + g_xz_0_xxxy_yyyyy[k];

                g_xz_0_xxxyy_yyyz[k] = -g_xz_0_xxxy_yyyz[k] * ab_y + g_xz_0_xxxy_yyyyz[k];

                g_xz_0_xxxyy_yyzz[k] = -g_xz_0_xxxy_yyzz[k] * ab_y + g_xz_0_xxxy_yyyzz[k];

                g_xz_0_xxxyy_yzzz[k] = -g_xz_0_xxxy_yzzz[k] * ab_y + g_xz_0_xxxy_yyzzz[k];

                g_xz_0_xxxyy_zzzz[k] = -g_xz_0_xxxy_zzzz[k] * ab_y + g_xz_0_xxxy_yzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 690 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 691 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 692 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 693 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 694 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 695 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 696 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 697 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 698 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 699 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 700 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 701 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 702 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 703 * ccomps * dcomps);

            auto g_xz_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyz_xxxx, g_xz_0_xxxyz_xxxy, g_xz_0_xxxyz_xxxz, g_xz_0_xxxyz_xxyy, g_xz_0_xxxyz_xxyz, g_xz_0_xxxyz_xxzz, g_xz_0_xxxyz_xyyy, g_xz_0_xxxyz_xyyz, g_xz_0_xxxyz_xyzz, g_xz_0_xxxyz_xzzz, g_xz_0_xxxyz_yyyy, g_xz_0_xxxyz_yyyz, g_xz_0_xxxyz_yyzz, g_xz_0_xxxyz_yzzz, g_xz_0_xxxyz_zzzz, g_xz_0_xxxz_xxxx, g_xz_0_xxxz_xxxxy, g_xz_0_xxxz_xxxy, g_xz_0_xxxz_xxxyy, g_xz_0_xxxz_xxxyz, g_xz_0_xxxz_xxxz, g_xz_0_xxxz_xxyy, g_xz_0_xxxz_xxyyy, g_xz_0_xxxz_xxyyz, g_xz_0_xxxz_xxyz, g_xz_0_xxxz_xxyzz, g_xz_0_xxxz_xxzz, g_xz_0_xxxz_xyyy, g_xz_0_xxxz_xyyyy, g_xz_0_xxxz_xyyyz, g_xz_0_xxxz_xyyz, g_xz_0_xxxz_xyyzz, g_xz_0_xxxz_xyzz, g_xz_0_xxxz_xyzzz, g_xz_0_xxxz_xzzz, g_xz_0_xxxz_yyyy, g_xz_0_xxxz_yyyyy, g_xz_0_xxxz_yyyyz, g_xz_0_xxxz_yyyz, g_xz_0_xxxz_yyyzz, g_xz_0_xxxz_yyzz, g_xz_0_xxxz_yyzzz, g_xz_0_xxxz_yzzz, g_xz_0_xxxz_yzzzz, g_xz_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyz_xxxx[k] = -g_xz_0_xxxz_xxxx[k] * ab_y + g_xz_0_xxxz_xxxxy[k];

                g_xz_0_xxxyz_xxxy[k] = -g_xz_0_xxxz_xxxy[k] * ab_y + g_xz_0_xxxz_xxxyy[k];

                g_xz_0_xxxyz_xxxz[k] = -g_xz_0_xxxz_xxxz[k] * ab_y + g_xz_0_xxxz_xxxyz[k];

                g_xz_0_xxxyz_xxyy[k] = -g_xz_0_xxxz_xxyy[k] * ab_y + g_xz_0_xxxz_xxyyy[k];

                g_xz_0_xxxyz_xxyz[k] = -g_xz_0_xxxz_xxyz[k] * ab_y + g_xz_0_xxxz_xxyyz[k];

                g_xz_0_xxxyz_xxzz[k] = -g_xz_0_xxxz_xxzz[k] * ab_y + g_xz_0_xxxz_xxyzz[k];

                g_xz_0_xxxyz_xyyy[k] = -g_xz_0_xxxz_xyyy[k] * ab_y + g_xz_0_xxxz_xyyyy[k];

                g_xz_0_xxxyz_xyyz[k] = -g_xz_0_xxxz_xyyz[k] * ab_y + g_xz_0_xxxz_xyyyz[k];

                g_xz_0_xxxyz_xyzz[k] = -g_xz_0_xxxz_xyzz[k] * ab_y + g_xz_0_xxxz_xyyzz[k];

                g_xz_0_xxxyz_xzzz[k] = -g_xz_0_xxxz_xzzz[k] * ab_y + g_xz_0_xxxz_xyzzz[k];

                g_xz_0_xxxyz_yyyy[k] = -g_xz_0_xxxz_yyyy[k] * ab_y + g_xz_0_xxxz_yyyyy[k];

                g_xz_0_xxxyz_yyyz[k] = -g_xz_0_xxxz_yyyz[k] * ab_y + g_xz_0_xxxz_yyyyz[k];

                g_xz_0_xxxyz_yyzz[k] = -g_xz_0_xxxz_yyzz[k] * ab_y + g_xz_0_xxxz_yyyzz[k];

                g_xz_0_xxxyz_yzzz[k] = -g_xz_0_xxxz_yzzz[k] * ab_y + g_xz_0_xxxz_yyzzz[k];

                g_xz_0_xxxyz_zzzz[k] = -g_xz_0_xxxz_zzzz[k] * ab_y + g_xz_0_xxxz_yzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 705 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 706 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 707 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 708 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 709 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 710 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 711 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 712 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 713 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 714 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 715 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 716 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 717 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 718 * ccomps * dcomps);

            auto g_xz_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzz_xxxx, g_xz_0_xxxzz_xxxy, g_xz_0_xxxzz_xxxz, g_xz_0_xxxzz_xxyy, g_xz_0_xxxzz_xxyz, g_xz_0_xxxzz_xxzz, g_xz_0_xxxzz_xyyy, g_xz_0_xxxzz_xyyz, g_xz_0_xxxzz_xyzz, g_xz_0_xxxzz_xzzz, g_xz_0_xxxzz_yyyy, g_xz_0_xxxzz_yyyz, g_xz_0_xxxzz_yyzz, g_xz_0_xxxzz_yzzz, g_xz_0_xxxzz_zzzz, g_xz_0_xxzz_xxxx, g_xz_0_xxzz_xxxxx, g_xz_0_xxzz_xxxxy, g_xz_0_xxzz_xxxxz, g_xz_0_xxzz_xxxy, g_xz_0_xxzz_xxxyy, g_xz_0_xxzz_xxxyz, g_xz_0_xxzz_xxxz, g_xz_0_xxzz_xxxzz, g_xz_0_xxzz_xxyy, g_xz_0_xxzz_xxyyy, g_xz_0_xxzz_xxyyz, g_xz_0_xxzz_xxyz, g_xz_0_xxzz_xxyzz, g_xz_0_xxzz_xxzz, g_xz_0_xxzz_xxzzz, g_xz_0_xxzz_xyyy, g_xz_0_xxzz_xyyyy, g_xz_0_xxzz_xyyyz, g_xz_0_xxzz_xyyz, g_xz_0_xxzz_xyyzz, g_xz_0_xxzz_xyzz, g_xz_0_xxzz_xyzzz, g_xz_0_xxzz_xzzz, g_xz_0_xxzz_xzzzz, g_xz_0_xxzz_yyyy, g_xz_0_xxzz_yyyz, g_xz_0_xxzz_yyzz, g_xz_0_xxzz_yzzz, g_xz_0_xxzz_zzzz, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_yyyy, g_z_0_xxzz_yyyz, g_z_0_xxzz_yyzz, g_z_0_xxzz_yzzz, g_z_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzz_xxxx[k] = -g_z_0_xxzz_xxxx[k] - g_xz_0_xxzz_xxxx[k] * ab_x + g_xz_0_xxzz_xxxxx[k];

                g_xz_0_xxxzz_xxxy[k] = -g_z_0_xxzz_xxxy[k] - g_xz_0_xxzz_xxxy[k] * ab_x + g_xz_0_xxzz_xxxxy[k];

                g_xz_0_xxxzz_xxxz[k] = -g_z_0_xxzz_xxxz[k] - g_xz_0_xxzz_xxxz[k] * ab_x + g_xz_0_xxzz_xxxxz[k];

                g_xz_0_xxxzz_xxyy[k] = -g_z_0_xxzz_xxyy[k] - g_xz_0_xxzz_xxyy[k] * ab_x + g_xz_0_xxzz_xxxyy[k];

                g_xz_0_xxxzz_xxyz[k] = -g_z_0_xxzz_xxyz[k] - g_xz_0_xxzz_xxyz[k] * ab_x + g_xz_0_xxzz_xxxyz[k];

                g_xz_0_xxxzz_xxzz[k] = -g_z_0_xxzz_xxzz[k] - g_xz_0_xxzz_xxzz[k] * ab_x + g_xz_0_xxzz_xxxzz[k];

                g_xz_0_xxxzz_xyyy[k] = -g_z_0_xxzz_xyyy[k] - g_xz_0_xxzz_xyyy[k] * ab_x + g_xz_0_xxzz_xxyyy[k];

                g_xz_0_xxxzz_xyyz[k] = -g_z_0_xxzz_xyyz[k] - g_xz_0_xxzz_xyyz[k] * ab_x + g_xz_0_xxzz_xxyyz[k];

                g_xz_0_xxxzz_xyzz[k] = -g_z_0_xxzz_xyzz[k] - g_xz_0_xxzz_xyzz[k] * ab_x + g_xz_0_xxzz_xxyzz[k];

                g_xz_0_xxxzz_xzzz[k] = -g_z_0_xxzz_xzzz[k] - g_xz_0_xxzz_xzzz[k] * ab_x + g_xz_0_xxzz_xxzzz[k];

                g_xz_0_xxxzz_yyyy[k] = -g_z_0_xxzz_yyyy[k] - g_xz_0_xxzz_yyyy[k] * ab_x + g_xz_0_xxzz_xyyyy[k];

                g_xz_0_xxxzz_yyyz[k] = -g_z_0_xxzz_yyyz[k] - g_xz_0_xxzz_yyyz[k] * ab_x + g_xz_0_xxzz_xyyyz[k];

                g_xz_0_xxxzz_yyzz[k] = -g_z_0_xxzz_yyzz[k] - g_xz_0_xxzz_yyzz[k] * ab_x + g_xz_0_xxzz_xyyzz[k];

                g_xz_0_xxxzz_yzzz[k] = -g_z_0_xxzz_yzzz[k] - g_xz_0_xxzz_yzzz[k] * ab_x + g_xz_0_xxzz_xyzzz[k];

                g_xz_0_xxxzz_zzzz[k] = -g_z_0_xxzz_zzzz[k] - g_xz_0_xxzz_zzzz[k] * ab_x + g_xz_0_xxzz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 720 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 721 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 722 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 723 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 724 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 725 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 726 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 727 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 728 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 729 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 730 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 731 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 732 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 733 * ccomps * dcomps);

            auto g_xz_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyy_xxxx, g_xz_0_xxyy_xxxxy, g_xz_0_xxyy_xxxy, g_xz_0_xxyy_xxxyy, g_xz_0_xxyy_xxxyz, g_xz_0_xxyy_xxxz, g_xz_0_xxyy_xxyy, g_xz_0_xxyy_xxyyy, g_xz_0_xxyy_xxyyz, g_xz_0_xxyy_xxyz, g_xz_0_xxyy_xxyzz, g_xz_0_xxyy_xxzz, g_xz_0_xxyy_xyyy, g_xz_0_xxyy_xyyyy, g_xz_0_xxyy_xyyyz, g_xz_0_xxyy_xyyz, g_xz_0_xxyy_xyyzz, g_xz_0_xxyy_xyzz, g_xz_0_xxyy_xyzzz, g_xz_0_xxyy_xzzz, g_xz_0_xxyy_yyyy, g_xz_0_xxyy_yyyyy, g_xz_0_xxyy_yyyyz, g_xz_0_xxyy_yyyz, g_xz_0_xxyy_yyyzz, g_xz_0_xxyy_yyzz, g_xz_0_xxyy_yyzzz, g_xz_0_xxyy_yzzz, g_xz_0_xxyy_yzzzz, g_xz_0_xxyy_zzzz, g_xz_0_xxyyy_xxxx, g_xz_0_xxyyy_xxxy, g_xz_0_xxyyy_xxxz, g_xz_0_xxyyy_xxyy, g_xz_0_xxyyy_xxyz, g_xz_0_xxyyy_xxzz, g_xz_0_xxyyy_xyyy, g_xz_0_xxyyy_xyyz, g_xz_0_xxyyy_xyzz, g_xz_0_xxyyy_xzzz, g_xz_0_xxyyy_yyyy, g_xz_0_xxyyy_yyyz, g_xz_0_xxyyy_yyzz, g_xz_0_xxyyy_yzzz, g_xz_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyy_xxxx[k] = -g_xz_0_xxyy_xxxx[k] * ab_y + g_xz_0_xxyy_xxxxy[k];

                g_xz_0_xxyyy_xxxy[k] = -g_xz_0_xxyy_xxxy[k] * ab_y + g_xz_0_xxyy_xxxyy[k];

                g_xz_0_xxyyy_xxxz[k] = -g_xz_0_xxyy_xxxz[k] * ab_y + g_xz_0_xxyy_xxxyz[k];

                g_xz_0_xxyyy_xxyy[k] = -g_xz_0_xxyy_xxyy[k] * ab_y + g_xz_0_xxyy_xxyyy[k];

                g_xz_0_xxyyy_xxyz[k] = -g_xz_0_xxyy_xxyz[k] * ab_y + g_xz_0_xxyy_xxyyz[k];

                g_xz_0_xxyyy_xxzz[k] = -g_xz_0_xxyy_xxzz[k] * ab_y + g_xz_0_xxyy_xxyzz[k];

                g_xz_0_xxyyy_xyyy[k] = -g_xz_0_xxyy_xyyy[k] * ab_y + g_xz_0_xxyy_xyyyy[k];

                g_xz_0_xxyyy_xyyz[k] = -g_xz_0_xxyy_xyyz[k] * ab_y + g_xz_0_xxyy_xyyyz[k];

                g_xz_0_xxyyy_xyzz[k] = -g_xz_0_xxyy_xyzz[k] * ab_y + g_xz_0_xxyy_xyyzz[k];

                g_xz_0_xxyyy_xzzz[k] = -g_xz_0_xxyy_xzzz[k] * ab_y + g_xz_0_xxyy_xyzzz[k];

                g_xz_0_xxyyy_yyyy[k] = -g_xz_0_xxyy_yyyy[k] * ab_y + g_xz_0_xxyy_yyyyy[k];

                g_xz_0_xxyyy_yyyz[k] = -g_xz_0_xxyy_yyyz[k] * ab_y + g_xz_0_xxyy_yyyyz[k];

                g_xz_0_xxyyy_yyzz[k] = -g_xz_0_xxyy_yyzz[k] * ab_y + g_xz_0_xxyy_yyyzz[k];

                g_xz_0_xxyyy_yzzz[k] = -g_xz_0_xxyy_yzzz[k] * ab_y + g_xz_0_xxyy_yyzzz[k];

                g_xz_0_xxyyy_zzzz[k] = -g_xz_0_xxyy_zzzz[k] * ab_y + g_xz_0_xxyy_yzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 735 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 736 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 737 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 738 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 739 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 740 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 741 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 742 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 743 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 744 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 745 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 746 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 747 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 748 * ccomps * dcomps);

            auto g_xz_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyz_xxxx, g_xz_0_xxyyz_xxxy, g_xz_0_xxyyz_xxxz, g_xz_0_xxyyz_xxyy, g_xz_0_xxyyz_xxyz, g_xz_0_xxyyz_xxzz, g_xz_0_xxyyz_xyyy, g_xz_0_xxyyz_xyyz, g_xz_0_xxyyz_xyzz, g_xz_0_xxyyz_xzzz, g_xz_0_xxyyz_yyyy, g_xz_0_xxyyz_yyyz, g_xz_0_xxyyz_yyzz, g_xz_0_xxyyz_yzzz, g_xz_0_xxyyz_zzzz, g_xz_0_xxyz_xxxx, g_xz_0_xxyz_xxxxy, g_xz_0_xxyz_xxxy, g_xz_0_xxyz_xxxyy, g_xz_0_xxyz_xxxyz, g_xz_0_xxyz_xxxz, g_xz_0_xxyz_xxyy, g_xz_0_xxyz_xxyyy, g_xz_0_xxyz_xxyyz, g_xz_0_xxyz_xxyz, g_xz_0_xxyz_xxyzz, g_xz_0_xxyz_xxzz, g_xz_0_xxyz_xyyy, g_xz_0_xxyz_xyyyy, g_xz_0_xxyz_xyyyz, g_xz_0_xxyz_xyyz, g_xz_0_xxyz_xyyzz, g_xz_0_xxyz_xyzz, g_xz_0_xxyz_xyzzz, g_xz_0_xxyz_xzzz, g_xz_0_xxyz_yyyy, g_xz_0_xxyz_yyyyy, g_xz_0_xxyz_yyyyz, g_xz_0_xxyz_yyyz, g_xz_0_xxyz_yyyzz, g_xz_0_xxyz_yyzz, g_xz_0_xxyz_yyzzz, g_xz_0_xxyz_yzzz, g_xz_0_xxyz_yzzzz, g_xz_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyz_xxxx[k] = -g_xz_0_xxyz_xxxx[k] * ab_y + g_xz_0_xxyz_xxxxy[k];

                g_xz_0_xxyyz_xxxy[k] = -g_xz_0_xxyz_xxxy[k] * ab_y + g_xz_0_xxyz_xxxyy[k];

                g_xz_0_xxyyz_xxxz[k] = -g_xz_0_xxyz_xxxz[k] * ab_y + g_xz_0_xxyz_xxxyz[k];

                g_xz_0_xxyyz_xxyy[k] = -g_xz_0_xxyz_xxyy[k] * ab_y + g_xz_0_xxyz_xxyyy[k];

                g_xz_0_xxyyz_xxyz[k] = -g_xz_0_xxyz_xxyz[k] * ab_y + g_xz_0_xxyz_xxyyz[k];

                g_xz_0_xxyyz_xxzz[k] = -g_xz_0_xxyz_xxzz[k] * ab_y + g_xz_0_xxyz_xxyzz[k];

                g_xz_0_xxyyz_xyyy[k] = -g_xz_0_xxyz_xyyy[k] * ab_y + g_xz_0_xxyz_xyyyy[k];

                g_xz_0_xxyyz_xyyz[k] = -g_xz_0_xxyz_xyyz[k] * ab_y + g_xz_0_xxyz_xyyyz[k];

                g_xz_0_xxyyz_xyzz[k] = -g_xz_0_xxyz_xyzz[k] * ab_y + g_xz_0_xxyz_xyyzz[k];

                g_xz_0_xxyyz_xzzz[k] = -g_xz_0_xxyz_xzzz[k] * ab_y + g_xz_0_xxyz_xyzzz[k];

                g_xz_0_xxyyz_yyyy[k] = -g_xz_0_xxyz_yyyy[k] * ab_y + g_xz_0_xxyz_yyyyy[k];

                g_xz_0_xxyyz_yyyz[k] = -g_xz_0_xxyz_yyyz[k] * ab_y + g_xz_0_xxyz_yyyyz[k];

                g_xz_0_xxyyz_yyzz[k] = -g_xz_0_xxyz_yyzz[k] * ab_y + g_xz_0_xxyz_yyyzz[k];

                g_xz_0_xxyyz_yzzz[k] = -g_xz_0_xxyz_yzzz[k] * ab_y + g_xz_0_xxyz_yyzzz[k];

                g_xz_0_xxyyz_zzzz[k] = -g_xz_0_xxyz_zzzz[k] * ab_y + g_xz_0_xxyz_yzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 750 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 751 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 752 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 753 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 754 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 755 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 756 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 757 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 758 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 759 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 760 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 761 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 762 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 763 * ccomps * dcomps);

            auto g_xz_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzz_xxxx, g_xz_0_xxyzz_xxxy, g_xz_0_xxyzz_xxxz, g_xz_0_xxyzz_xxyy, g_xz_0_xxyzz_xxyz, g_xz_0_xxyzz_xxzz, g_xz_0_xxyzz_xyyy, g_xz_0_xxyzz_xyyz, g_xz_0_xxyzz_xyzz, g_xz_0_xxyzz_xzzz, g_xz_0_xxyzz_yyyy, g_xz_0_xxyzz_yyyz, g_xz_0_xxyzz_yyzz, g_xz_0_xxyzz_yzzz, g_xz_0_xxyzz_zzzz, g_xz_0_xxzz_xxxx, g_xz_0_xxzz_xxxxy, g_xz_0_xxzz_xxxy, g_xz_0_xxzz_xxxyy, g_xz_0_xxzz_xxxyz, g_xz_0_xxzz_xxxz, g_xz_0_xxzz_xxyy, g_xz_0_xxzz_xxyyy, g_xz_0_xxzz_xxyyz, g_xz_0_xxzz_xxyz, g_xz_0_xxzz_xxyzz, g_xz_0_xxzz_xxzz, g_xz_0_xxzz_xyyy, g_xz_0_xxzz_xyyyy, g_xz_0_xxzz_xyyyz, g_xz_0_xxzz_xyyz, g_xz_0_xxzz_xyyzz, g_xz_0_xxzz_xyzz, g_xz_0_xxzz_xyzzz, g_xz_0_xxzz_xzzz, g_xz_0_xxzz_yyyy, g_xz_0_xxzz_yyyyy, g_xz_0_xxzz_yyyyz, g_xz_0_xxzz_yyyz, g_xz_0_xxzz_yyyzz, g_xz_0_xxzz_yyzz, g_xz_0_xxzz_yyzzz, g_xz_0_xxzz_yzzz, g_xz_0_xxzz_yzzzz, g_xz_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzz_xxxx[k] = -g_xz_0_xxzz_xxxx[k] * ab_y + g_xz_0_xxzz_xxxxy[k];

                g_xz_0_xxyzz_xxxy[k] = -g_xz_0_xxzz_xxxy[k] * ab_y + g_xz_0_xxzz_xxxyy[k];

                g_xz_0_xxyzz_xxxz[k] = -g_xz_0_xxzz_xxxz[k] * ab_y + g_xz_0_xxzz_xxxyz[k];

                g_xz_0_xxyzz_xxyy[k] = -g_xz_0_xxzz_xxyy[k] * ab_y + g_xz_0_xxzz_xxyyy[k];

                g_xz_0_xxyzz_xxyz[k] = -g_xz_0_xxzz_xxyz[k] * ab_y + g_xz_0_xxzz_xxyyz[k];

                g_xz_0_xxyzz_xxzz[k] = -g_xz_0_xxzz_xxzz[k] * ab_y + g_xz_0_xxzz_xxyzz[k];

                g_xz_0_xxyzz_xyyy[k] = -g_xz_0_xxzz_xyyy[k] * ab_y + g_xz_0_xxzz_xyyyy[k];

                g_xz_0_xxyzz_xyyz[k] = -g_xz_0_xxzz_xyyz[k] * ab_y + g_xz_0_xxzz_xyyyz[k];

                g_xz_0_xxyzz_xyzz[k] = -g_xz_0_xxzz_xyzz[k] * ab_y + g_xz_0_xxzz_xyyzz[k];

                g_xz_0_xxyzz_xzzz[k] = -g_xz_0_xxzz_xzzz[k] * ab_y + g_xz_0_xxzz_xyzzz[k];

                g_xz_0_xxyzz_yyyy[k] = -g_xz_0_xxzz_yyyy[k] * ab_y + g_xz_0_xxzz_yyyyy[k];

                g_xz_0_xxyzz_yyyz[k] = -g_xz_0_xxzz_yyyz[k] * ab_y + g_xz_0_xxzz_yyyyz[k];

                g_xz_0_xxyzz_yyzz[k] = -g_xz_0_xxzz_yyzz[k] * ab_y + g_xz_0_xxzz_yyyzz[k];

                g_xz_0_xxyzz_yzzz[k] = -g_xz_0_xxzz_yzzz[k] * ab_y + g_xz_0_xxzz_yyzzz[k];

                g_xz_0_xxyzz_zzzz[k] = -g_xz_0_xxzz_zzzz[k] * ab_y + g_xz_0_xxzz_yzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 765 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 766 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 767 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 768 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 769 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 770 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 771 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 772 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 773 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 774 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 775 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 776 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 777 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 778 * ccomps * dcomps);

            auto g_xz_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzz_xxxx, g_xz_0_xxzzz_xxxy, g_xz_0_xxzzz_xxxz, g_xz_0_xxzzz_xxyy, g_xz_0_xxzzz_xxyz, g_xz_0_xxzzz_xxzz, g_xz_0_xxzzz_xyyy, g_xz_0_xxzzz_xyyz, g_xz_0_xxzzz_xyzz, g_xz_0_xxzzz_xzzz, g_xz_0_xxzzz_yyyy, g_xz_0_xxzzz_yyyz, g_xz_0_xxzzz_yyzz, g_xz_0_xxzzz_yzzz, g_xz_0_xxzzz_zzzz, g_xz_0_xzzz_xxxx, g_xz_0_xzzz_xxxxx, g_xz_0_xzzz_xxxxy, g_xz_0_xzzz_xxxxz, g_xz_0_xzzz_xxxy, g_xz_0_xzzz_xxxyy, g_xz_0_xzzz_xxxyz, g_xz_0_xzzz_xxxz, g_xz_0_xzzz_xxxzz, g_xz_0_xzzz_xxyy, g_xz_0_xzzz_xxyyy, g_xz_0_xzzz_xxyyz, g_xz_0_xzzz_xxyz, g_xz_0_xzzz_xxyzz, g_xz_0_xzzz_xxzz, g_xz_0_xzzz_xxzzz, g_xz_0_xzzz_xyyy, g_xz_0_xzzz_xyyyy, g_xz_0_xzzz_xyyyz, g_xz_0_xzzz_xyyz, g_xz_0_xzzz_xyyzz, g_xz_0_xzzz_xyzz, g_xz_0_xzzz_xyzzz, g_xz_0_xzzz_xzzz, g_xz_0_xzzz_xzzzz, g_xz_0_xzzz_yyyy, g_xz_0_xzzz_yyyz, g_xz_0_xzzz_yyzz, g_xz_0_xzzz_yzzz, g_xz_0_xzzz_zzzz, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_yyyy, g_z_0_xzzz_yyyz, g_z_0_xzzz_yyzz, g_z_0_xzzz_yzzz, g_z_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzz_xxxx[k] = -g_z_0_xzzz_xxxx[k] - g_xz_0_xzzz_xxxx[k] * ab_x + g_xz_0_xzzz_xxxxx[k];

                g_xz_0_xxzzz_xxxy[k] = -g_z_0_xzzz_xxxy[k] - g_xz_0_xzzz_xxxy[k] * ab_x + g_xz_0_xzzz_xxxxy[k];

                g_xz_0_xxzzz_xxxz[k] = -g_z_0_xzzz_xxxz[k] - g_xz_0_xzzz_xxxz[k] * ab_x + g_xz_0_xzzz_xxxxz[k];

                g_xz_0_xxzzz_xxyy[k] = -g_z_0_xzzz_xxyy[k] - g_xz_0_xzzz_xxyy[k] * ab_x + g_xz_0_xzzz_xxxyy[k];

                g_xz_0_xxzzz_xxyz[k] = -g_z_0_xzzz_xxyz[k] - g_xz_0_xzzz_xxyz[k] * ab_x + g_xz_0_xzzz_xxxyz[k];

                g_xz_0_xxzzz_xxzz[k] = -g_z_0_xzzz_xxzz[k] - g_xz_0_xzzz_xxzz[k] * ab_x + g_xz_0_xzzz_xxxzz[k];

                g_xz_0_xxzzz_xyyy[k] = -g_z_0_xzzz_xyyy[k] - g_xz_0_xzzz_xyyy[k] * ab_x + g_xz_0_xzzz_xxyyy[k];

                g_xz_0_xxzzz_xyyz[k] = -g_z_0_xzzz_xyyz[k] - g_xz_0_xzzz_xyyz[k] * ab_x + g_xz_0_xzzz_xxyyz[k];

                g_xz_0_xxzzz_xyzz[k] = -g_z_0_xzzz_xyzz[k] - g_xz_0_xzzz_xyzz[k] * ab_x + g_xz_0_xzzz_xxyzz[k];

                g_xz_0_xxzzz_xzzz[k] = -g_z_0_xzzz_xzzz[k] - g_xz_0_xzzz_xzzz[k] * ab_x + g_xz_0_xzzz_xxzzz[k];

                g_xz_0_xxzzz_yyyy[k] = -g_z_0_xzzz_yyyy[k] - g_xz_0_xzzz_yyyy[k] * ab_x + g_xz_0_xzzz_xyyyy[k];

                g_xz_0_xxzzz_yyyz[k] = -g_z_0_xzzz_yyyz[k] - g_xz_0_xzzz_yyyz[k] * ab_x + g_xz_0_xzzz_xyyyz[k];

                g_xz_0_xxzzz_yyzz[k] = -g_z_0_xzzz_yyzz[k] - g_xz_0_xzzz_yyzz[k] * ab_x + g_xz_0_xzzz_xyyzz[k];

                g_xz_0_xxzzz_yzzz[k] = -g_z_0_xzzz_yzzz[k] - g_xz_0_xzzz_yzzz[k] * ab_x + g_xz_0_xzzz_xyzzz[k];

                g_xz_0_xxzzz_zzzz[k] = -g_z_0_xzzz_zzzz[k] - g_xz_0_xzzz_zzzz[k] * ab_x + g_xz_0_xzzz_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 780 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 781 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 782 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 783 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 784 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 785 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 786 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 787 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 788 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 789 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 790 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 791 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 792 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 793 * ccomps * dcomps);

            auto g_xz_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyy_xxxx, g_xz_0_xyyy_xxxxy, g_xz_0_xyyy_xxxy, g_xz_0_xyyy_xxxyy, g_xz_0_xyyy_xxxyz, g_xz_0_xyyy_xxxz, g_xz_0_xyyy_xxyy, g_xz_0_xyyy_xxyyy, g_xz_0_xyyy_xxyyz, g_xz_0_xyyy_xxyz, g_xz_0_xyyy_xxyzz, g_xz_0_xyyy_xxzz, g_xz_0_xyyy_xyyy, g_xz_0_xyyy_xyyyy, g_xz_0_xyyy_xyyyz, g_xz_0_xyyy_xyyz, g_xz_0_xyyy_xyyzz, g_xz_0_xyyy_xyzz, g_xz_0_xyyy_xyzzz, g_xz_0_xyyy_xzzz, g_xz_0_xyyy_yyyy, g_xz_0_xyyy_yyyyy, g_xz_0_xyyy_yyyyz, g_xz_0_xyyy_yyyz, g_xz_0_xyyy_yyyzz, g_xz_0_xyyy_yyzz, g_xz_0_xyyy_yyzzz, g_xz_0_xyyy_yzzz, g_xz_0_xyyy_yzzzz, g_xz_0_xyyy_zzzz, g_xz_0_xyyyy_xxxx, g_xz_0_xyyyy_xxxy, g_xz_0_xyyyy_xxxz, g_xz_0_xyyyy_xxyy, g_xz_0_xyyyy_xxyz, g_xz_0_xyyyy_xxzz, g_xz_0_xyyyy_xyyy, g_xz_0_xyyyy_xyyz, g_xz_0_xyyyy_xyzz, g_xz_0_xyyyy_xzzz, g_xz_0_xyyyy_yyyy, g_xz_0_xyyyy_yyyz, g_xz_0_xyyyy_yyzz, g_xz_0_xyyyy_yzzz, g_xz_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyy_xxxx[k] = -g_xz_0_xyyy_xxxx[k] * ab_y + g_xz_0_xyyy_xxxxy[k];

                g_xz_0_xyyyy_xxxy[k] = -g_xz_0_xyyy_xxxy[k] * ab_y + g_xz_0_xyyy_xxxyy[k];

                g_xz_0_xyyyy_xxxz[k] = -g_xz_0_xyyy_xxxz[k] * ab_y + g_xz_0_xyyy_xxxyz[k];

                g_xz_0_xyyyy_xxyy[k] = -g_xz_0_xyyy_xxyy[k] * ab_y + g_xz_0_xyyy_xxyyy[k];

                g_xz_0_xyyyy_xxyz[k] = -g_xz_0_xyyy_xxyz[k] * ab_y + g_xz_0_xyyy_xxyyz[k];

                g_xz_0_xyyyy_xxzz[k] = -g_xz_0_xyyy_xxzz[k] * ab_y + g_xz_0_xyyy_xxyzz[k];

                g_xz_0_xyyyy_xyyy[k] = -g_xz_0_xyyy_xyyy[k] * ab_y + g_xz_0_xyyy_xyyyy[k];

                g_xz_0_xyyyy_xyyz[k] = -g_xz_0_xyyy_xyyz[k] * ab_y + g_xz_0_xyyy_xyyyz[k];

                g_xz_0_xyyyy_xyzz[k] = -g_xz_0_xyyy_xyzz[k] * ab_y + g_xz_0_xyyy_xyyzz[k];

                g_xz_0_xyyyy_xzzz[k] = -g_xz_0_xyyy_xzzz[k] * ab_y + g_xz_0_xyyy_xyzzz[k];

                g_xz_0_xyyyy_yyyy[k] = -g_xz_0_xyyy_yyyy[k] * ab_y + g_xz_0_xyyy_yyyyy[k];

                g_xz_0_xyyyy_yyyz[k] = -g_xz_0_xyyy_yyyz[k] * ab_y + g_xz_0_xyyy_yyyyz[k];

                g_xz_0_xyyyy_yyzz[k] = -g_xz_0_xyyy_yyzz[k] * ab_y + g_xz_0_xyyy_yyyzz[k];

                g_xz_0_xyyyy_yzzz[k] = -g_xz_0_xyyy_yzzz[k] * ab_y + g_xz_0_xyyy_yyzzz[k];

                g_xz_0_xyyyy_zzzz[k] = -g_xz_0_xyyy_zzzz[k] * ab_y + g_xz_0_xyyy_yzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 795 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 796 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 797 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 798 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 799 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 800 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 801 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 802 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 803 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 804 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 805 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 806 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 807 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 808 * ccomps * dcomps);

            auto g_xz_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyz_xxxx, g_xz_0_xyyyz_xxxy, g_xz_0_xyyyz_xxxz, g_xz_0_xyyyz_xxyy, g_xz_0_xyyyz_xxyz, g_xz_0_xyyyz_xxzz, g_xz_0_xyyyz_xyyy, g_xz_0_xyyyz_xyyz, g_xz_0_xyyyz_xyzz, g_xz_0_xyyyz_xzzz, g_xz_0_xyyyz_yyyy, g_xz_0_xyyyz_yyyz, g_xz_0_xyyyz_yyzz, g_xz_0_xyyyz_yzzz, g_xz_0_xyyyz_zzzz, g_xz_0_xyyz_xxxx, g_xz_0_xyyz_xxxxy, g_xz_0_xyyz_xxxy, g_xz_0_xyyz_xxxyy, g_xz_0_xyyz_xxxyz, g_xz_0_xyyz_xxxz, g_xz_0_xyyz_xxyy, g_xz_0_xyyz_xxyyy, g_xz_0_xyyz_xxyyz, g_xz_0_xyyz_xxyz, g_xz_0_xyyz_xxyzz, g_xz_0_xyyz_xxzz, g_xz_0_xyyz_xyyy, g_xz_0_xyyz_xyyyy, g_xz_0_xyyz_xyyyz, g_xz_0_xyyz_xyyz, g_xz_0_xyyz_xyyzz, g_xz_0_xyyz_xyzz, g_xz_0_xyyz_xyzzz, g_xz_0_xyyz_xzzz, g_xz_0_xyyz_yyyy, g_xz_0_xyyz_yyyyy, g_xz_0_xyyz_yyyyz, g_xz_0_xyyz_yyyz, g_xz_0_xyyz_yyyzz, g_xz_0_xyyz_yyzz, g_xz_0_xyyz_yyzzz, g_xz_0_xyyz_yzzz, g_xz_0_xyyz_yzzzz, g_xz_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyz_xxxx[k] = -g_xz_0_xyyz_xxxx[k] * ab_y + g_xz_0_xyyz_xxxxy[k];

                g_xz_0_xyyyz_xxxy[k] = -g_xz_0_xyyz_xxxy[k] * ab_y + g_xz_0_xyyz_xxxyy[k];

                g_xz_0_xyyyz_xxxz[k] = -g_xz_0_xyyz_xxxz[k] * ab_y + g_xz_0_xyyz_xxxyz[k];

                g_xz_0_xyyyz_xxyy[k] = -g_xz_0_xyyz_xxyy[k] * ab_y + g_xz_0_xyyz_xxyyy[k];

                g_xz_0_xyyyz_xxyz[k] = -g_xz_0_xyyz_xxyz[k] * ab_y + g_xz_0_xyyz_xxyyz[k];

                g_xz_0_xyyyz_xxzz[k] = -g_xz_0_xyyz_xxzz[k] * ab_y + g_xz_0_xyyz_xxyzz[k];

                g_xz_0_xyyyz_xyyy[k] = -g_xz_0_xyyz_xyyy[k] * ab_y + g_xz_0_xyyz_xyyyy[k];

                g_xz_0_xyyyz_xyyz[k] = -g_xz_0_xyyz_xyyz[k] * ab_y + g_xz_0_xyyz_xyyyz[k];

                g_xz_0_xyyyz_xyzz[k] = -g_xz_0_xyyz_xyzz[k] * ab_y + g_xz_0_xyyz_xyyzz[k];

                g_xz_0_xyyyz_xzzz[k] = -g_xz_0_xyyz_xzzz[k] * ab_y + g_xz_0_xyyz_xyzzz[k];

                g_xz_0_xyyyz_yyyy[k] = -g_xz_0_xyyz_yyyy[k] * ab_y + g_xz_0_xyyz_yyyyy[k];

                g_xz_0_xyyyz_yyyz[k] = -g_xz_0_xyyz_yyyz[k] * ab_y + g_xz_0_xyyz_yyyyz[k];

                g_xz_0_xyyyz_yyzz[k] = -g_xz_0_xyyz_yyzz[k] * ab_y + g_xz_0_xyyz_yyyzz[k];

                g_xz_0_xyyyz_yzzz[k] = -g_xz_0_xyyz_yzzz[k] * ab_y + g_xz_0_xyyz_yyzzz[k];

                g_xz_0_xyyyz_zzzz[k] = -g_xz_0_xyyz_zzzz[k] * ab_y + g_xz_0_xyyz_yzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 810 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 811 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 812 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 813 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 814 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 815 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 816 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 817 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 818 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 819 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 820 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 821 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 822 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 823 * ccomps * dcomps);

            auto g_xz_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzz_xxxx, g_xz_0_xyyzz_xxxy, g_xz_0_xyyzz_xxxz, g_xz_0_xyyzz_xxyy, g_xz_0_xyyzz_xxyz, g_xz_0_xyyzz_xxzz, g_xz_0_xyyzz_xyyy, g_xz_0_xyyzz_xyyz, g_xz_0_xyyzz_xyzz, g_xz_0_xyyzz_xzzz, g_xz_0_xyyzz_yyyy, g_xz_0_xyyzz_yyyz, g_xz_0_xyyzz_yyzz, g_xz_0_xyyzz_yzzz, g_xz_0_xyyzz_zzzz, g_xz_0_xyzz_xxxx, g_xz_0_xyzz_xxxxy, g_xz_0_xyzz_xxxy, g_xz_0_xyzz_xxxyy, g_xz_0_xyzz_xxxyz, g_xz_0_xyzz_xxxz, g_xz_0_xyzz_xxyy, g_xz_0_xyzz_xxyyy, g_xz_0_xyzz_xxyyz, g_xz_0_xyzz_xxyz, g_xz_0_xyzz_xxyzz, g_xz_0_xyzz_xxzz, g_xz_0_xyzz_xyyy, g_xz_0_xyzz_xyyyy, g_xz_0_xyzz_xyyyz, g_xz_0_xyzz_xyyz, g_xz_0_xyzz_xyyzz, g_xz_0_xyzz_xyzz, g_xz_0_xyzz_xyzzz, g_xz_0_xyzz_xzzz, g_xz_0_xyzz_yyyy, g_xz_0_xyzz_yyyyy, g_xz_0_xyzz_yyyyz, g_xz_0_xyzz_yyyz, g_xz_0_xyzz_yyyzz, g_xz_0_xyzz_yyzz, g_xz_0_xyzz_yyzzz, g_xz_0_xyzz_yzzz, g_xz_0_xyzz_yzzzz, g_xz_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzz_xxxx[k] = -g_xz_0_xyzz_xxxx[k] * ab_y + g_xz_0_xyzz_xxxxy[k];

                g_xz_0_xyyzz_xxxy[k] = -g_xz_0_xyzz_xxxy[k] * ab_y + g_xz_0_xyzz_xxxyy[k];

                g_xz_0_xyyzz_xxxz[k] = -g_xz_0_xyzz_xxxz[k] * ab_y + g_xz_0_xyzz_xxxyz[k];

                g_xz_0_xyyzz_xxyy[k] = -g_xz_0_xyzz_xxyy[k] * ab_y + g_xz_0_xyzz_xxyyy[k];

                g_xz_0_xyyzz_xxyz[k] = -g_xz_0_xyzz_xxyz[k] * ab_y + g_xz_0_xyzz_xxyyz[k];

                g_xz_0_xyyzz_xxzz[k] = -g_xz_0_xyzz_xxzz[k] * ab_y + g_xz_0_xyzz_xxyzz[k];

                g_xz_0_xyyzz_xyyy[k] = -g_xz_0_xyzz_xyyy[k] * ab_y + g_xz_0_xyzz_xyyyy[k];

                g_xz_0_xyyzz_xyyz[k] = -g_xz_0_xyzz_xyyz[k] * ab_y + g_xz_0_xyzz_xyyyz[k];

                g_xz_0_xyyzz_xyzz[k] = -g_xz_0_xyzz_xyzz[k] * ab_y + g_xz_0_xyzz_xyyzz[k];

                g_xz_0_xyyzz_xzzz[k] = -g_xz_0_xyzz_xzzz[k] * ab_y + g_xz_0_xyzz_xyzzz[k];

                g_xz_0_xyyzz_yyyy[k] = -g_xz_0_xyzz_yyyy[k] * ab_y + g_xz_0_xyzz_yyyyy[k];

                g_xz_0_xyyzz_yyyz[k] = -g_xz_0_xyzz_yyyz[k] * ab_y + g_xz_0_xyzz_yyyyz[k];

                g_xz_0_xyyzz_yyzz[k] = -g_xz_0_xyzz_yyzz[k] * ab_y + g_xz_0_xyzz_yyyzz[k];

                g_xz_0_xyyzz_yzzz[k] = -g_xz_0_xyzz_yzzz[k] * ab_y + g_xz_0_xyzz_yyzzz[k];

                g_xz_0_xyyzz_zzzz[k] = -g_xz_0_xyzz_zzzz[k] * ab_y + g_xz_0_xyzz_yzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 825 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 826 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 827 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 828 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 829 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 830 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 831 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 832 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 833 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 834 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 835 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 836 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 837 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 838 * ccomps * dcomps);

            auto g_xz_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzz_xxxx, g_xz_0_xyzzz_xxxy, g_xz_0_xyzzz_xxxz, g_xz_0_xyzzz_xxyy, g_xz_0_xyzzz_xxyz, g_xz_0_xyzzz_xxzz, g_xz_0_xyzzz_xyyy, g_xz_0_xyzzz_xyyz, g_xz_0_xyzzz_xyzz, g_xz_0_xyzzz_xzzz, g_xz_0_xyzzz_yyyy, g_xz_0_xyzzz_yyyz, g_xz_0_xyzzz_yyzz, g_xz_0_xyzzz_yzzz, g_xz_0_xyzzz_zzzz, g_xz_0_xzzz_xxxx, g_xz_0_xzzz_xxxxy, g_xz_0_xzzz_xxxy, g_xz_0_xzzz_xxxyy, g_xz_0_xzzz_xxxyz, g_xz_0_xzzz_xxxz, g_xz_0_xzzz_xxyy, g_xz_0_xzzz_xxyyy, g_xz_0_xzzz_xxyyz, g_xz_0_xzzz_xxyz, g_xz_0_xzzz_xxyzz, g_xz_0_xzzz_xxzz, g_xz_0_xzzz_xyyy, g_xz_0_xzzz_xyyyy, g_xz_0_xzzz_xyyyz, g_xz_0_xzzz_xyyz, g_xz_0_xzzz_xyyzz, g_xz_0_xzzz_xyzz, g_xz_0_xzzz_xyzzz, g_xz_0_xzzz_xzzz, g_xz_0_xzzz_yyyy, g_xz_0_xzzz_yyyyy, g_xz_0_xzzz_yyyyz, g_xz_0_xzzz_yyyz, g_xz_0_xzzz_yyyzz, g_xz_0_xzzz_yyzz, g_xz_0_xzzz_yyzzz, g_xz_0_xzzz_yzzz, g_xz_0_xzzz_yzzzz, g_xz_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzz_xxxx[k] = -g_xz_0_xzzz_xxxx[k] * ab_y + g_xz_0_xzzz_xxxxy[k];

                g_xz_0_xyzzz_xxxy[k] = -g_xz_0_xzzz_xxxy[k] * ab_y + g_xz_0_xzzz_xxxyy[k];

                g_xz_0_xyzzz_xxxz[k] = -g_xz_0_xzzz_xxxz[k] * ab_y + g_xz_0_xzzz_xxxyz[k];

                g_xz_0_xyzzz_xxyy[k] = -g_xz_0_xzzz_xxyy[k] * ab_y + g_xz_0_xzzz_xxyyy[k];

                g_xz_0_xyzzz_xxyz[k] = -g_xz_0_xzzz_xxyz[k] * ab_y + g_xz_0_xzzz_xxyyz[k];

                g_xz_0_xyzzz_xxzz[k] = -g_xz_0_xzzz_xxzz[k] * ab_y + g_xz_0_xzzz_xxyzz[k];

                g_xz_0_xyzzz_xyyy[k] = -g_xz_0_xzzz_xyyy[k] * ab_y + g_xz_0_xzzz_xyyyy[k];

                g_xz_0_xyzzz_xyyz[k] = -g_xz_0_xzzz_xyyz[k] * ab_y + g_xz_0_xzzz_xyyyz[k];

                g_xz_0_xyzzz_xyzz[k] = -g_xz_0_xzzz_xyzz[k] * ab_y + g_xz_0_xzzz_xyyzz[k];

                g_xz_0_xyzzz_xzzz[k] = -g_xz_0_xzzz_xzzz[k] * ab_y + g_xz_0_xzzz_xyzzz[k];

                g_xz_0_xyzzz_yyyy[k] = -g_xz_0_xzzz_yyyy[k] * ab_y + g_xz_0_xzzz_yyyyy[k];

                g_xz_0_xyzzz_yyyz[k] = -g_xz_0_xzzz_yyyz[k] * ab_y + g_xz_0_xzzz_yyyyz[k];

                g_xz_0_xyzzz_yyzz[k] = -g_xz_0_xzzz_yyzz[k] * ab_y + g_xz_0_xzzz_yyyzz[k];

                g_xz_0_xyzzz_yzzz[k] = -g_xz_0_xzzz_yzzz[k] * ab_y + g_xz_0_xzzz_yyzzz[k];

                g_xz_0_xyzzz_zzzz[k] = -g_xz_0_xzzz_zzzz[k] * ab_y + g_xz_0_xzzz_yzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 840 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 841 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 842 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 843 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 844 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 845 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 846 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 847 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 848 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 849 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 850 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 851 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 852 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 853 * ccomps * dcomps);

            auto g_xz_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzz_xxxx, g_xz_0_xzzzz_xxxy, g_xz_0_xzzzz_xxxz, g_xz_0_xzzzz_xxyy, g_xz_0_xzzzz_xxyz, g_xz_0_xzzzz_xxzz, g_xz_0_xzzzz_xyyy, g_xz_0_xzzzz_xyyz, g_xz_0_xzzzz_xyzz, g_xz_0_xzzzz_xzzz, g_xz_0_xzzzz_yyyy, g_xz_0_xzzzz_yyyz, g_xz_0_xzzzz_yyzz, g_xz_0_xzzzz_yzzz, g_xz_0_xzzzz_zzzz, g_xz_0_zzzz_xxxx, g_xz_0_zzzz_xxxxx, g_xz_0_zzzz_xxxxy, g_xz_0_zzzz_xxxxz, g_xz_0_zzzz_xxxy, g_xz_0_zzzz_xxxyy, g_xz_0_zzzz_xxxyz, g_xz_0_zzzz_xxxz, g_xz_0_zzzz_xxxzz, g_xz_0_zzzz_xxyy, g_xz_0_zzzz_xxyyy, g_xz_0_zzzz_xxyyz, g_xz_0_zzzz_xxyz, g_xz_0_zzzz_xxyzz, g_xz_0_zzzz_xxzz, g_xz_0_zzzz_xxzzz, g_xz_0_zzzz_xyyy, g_xz_0_zzzz_xyyyy, g_xz_0_zzzz_xyyyz, g_xz_0_zzzz_xyyz, g_xz_0_zzzz_xyyzz, g_xz_0_zzzz_xyzz, g_xz_0_zzzz_xyzzz, g_xz_0_zzzz_xzzz, g_xz_0_zzzz_xzzzz, g_xz_0_zzzz_yyyy, g_xz_0_zzzz_yyyz, g_xz_0_zzzz_yyzz, g_xz_0_zzzz_yzzz, g_xz_0_zzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] - g_xz_0_zzzz_xxxx[k] * ab_x + g_xz_0_zzzz_xxxxx[k];

                g_xz_0_xzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] - g_xz_0_zzzz_xxxy[k] * ab_x + g_xz_0_zzzz_xxxxy[k];

                g_xz_0_xzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] - g_xz_0_zzzz_xxxz[k] * ab_x + g_xz_0_zzzz_xxxxz[k];

                g_xz_0_xzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] - g_xz_0_zzzz_xxyy[k] * ab_x + g_xz_0_zzzz_xxxyy[k];

                g_xz_0_xzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] - g_xz_0_zzzz_xxyz[k] * ab_x + g_xz_0_zzzz_xxxyz[k];

                g_xz_0_xzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] - g_xz_0_zzzz_xxzz[k] * ab_x + g_xz_0_zzzz_xxxzz[k];

                g_xz_0_xzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] - g_xz_0_zzzz_xyyy[k] * ab_x + g_xz_0_zzzz_xxyyy[k];

                g_xz_0_xzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] - g_xz_0_zzzz_xyyz[k] * ab_x + g_xz_0_zzzz_xxyyz[k];

                g_xz_0_xzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] - g_xz_0_zzzz_xyzz[k] * ab_x + g_xz_0_zzzz_xxyzz[k];

                g_xz_0_xzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] - g_xz_0_zzzz_xzzz[k] * ab_x + g_xz_0_zzzz_xxzzz[k];

                g_xz_0_xzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] - g_xz_0_zzzz_yyyy[k] * ab_x + g_xz_0_zzzz_xyyyy[k];

                g_xz_0_xzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] - g_xz_0_zzzz_yyyz[k] * ab_x + g_xz_0_zzzz_xyyyz[k];

                g_xz_0_xzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] - g_xz_0_zzzz_yyzz[k] * ab_x + g_xz_0_zzzz_xyyzz[k];

                g_xz_0_xzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] - g_xz_0_zzzz_yzzz[k] * ab_x + g_xz_0_zzzz_xyzzz[k];

                g_xz_0_xzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] - g_xz_0_zzzz_zzzz[k] * ab_x + g_xz_0_zzzz_xzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 855 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 856 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 857 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 858 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 859 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 860 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 861 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 862 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 863 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 864 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 865 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 866 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 867 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 868 * ccomps * dcomps);

            auto g_xz_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyy_xxxx, g_xz_0_yyyy_xxxxy, g_xz_0_yyyy_xxxy, g_xz_0_yyyy_xxxyy, g_xz_0_yyyy_xxxyz, g_xz_0_yyyy_xxxz, g_xz_0_yyyy_xxyy, g_xz_0_yyyy_xxyyy, g_xz_0_yyyy_xxyyz, g_xz_0_yyyy_xxyz, g_xz_0_yyyy_xxyzz, g_xz_0_yyyy_xxzz, g_xz_0_yyyy_xyyy, g_xz_0_yyyy_xyyyy, g_xz_0_yyyy_xyyyz, g_xz_0_yyyy_xyyz, g_xz_0_yyyy_xyyzz, g_xz_0_yyyy_xyzz, g_xz_0_yyyy_xyzzz, g_xz_0_yyyy_xzzz, g_xz_0_yyyy_yyyy, g_xz_0_yyyy_yyyyy, g_xz_0_yyyy_yyyyz, g_xz_0_yyyy_yyyz, g_xz_0_yyyy_yyyzz, g_xz_0_yyyy_yyzz, g_xz_0_yyyy_yyzzz, g_xz_0_yyyy_yzzz, g_xz_0_yyyy_yzzzz, g_xz_0_yyyy_zzzz, g_xz_0_yyyyy_xxxx, g_xz_0_yyyyy_xxxy, g_xz_0_yyyyy_xxxz, g_xz_0_yyyyy_xxyy, g_xz_0_yyyyy_xxyz, g_xz_0_yyyyy_xxzz, g_xz_0_yyyyy_xyyy, g_xz_0_yyyyy_xyyz, g_xz_0_yyyyy_xyzz, g_xz_0_yyyyy_xzzz, g_xz_0_yyyyy_yyyy, g_xz_0_yyyyy_yyyz, g_xz_0_yyyyy_yyzz, g_xz_0_yyyyy_yzzz, g_xz_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyy_xxxx[k] = -g_xz_0_yyyy_xxxx[k] * ab_y + g_xz_0_yyyy_xxxxy[k];

                g_xz_0_yyyyy_xxxy[k] = -g_xz_0_yyyy_xxxy[k] * ab_y + g_xz_0_yyyy_xxxyy[k];

                g_xz_0_yyyyy_xxxz[k] = -g_xz_0_yyyy_xxxz[k] * ab_y + g_xz_0_yyyy_xxxyz[k];

                g_xz_0_yyyyy_xxyy[k] = -g_xz_0_yyyy_xxyy[k] * ab_y + g_xz_0_yyyy_xxyyy[k];

                g_xz_0_yyyyy_xxyz[k] = -g_xz_0_yyyy_xxyz[k] * ab_y + g_xz_0_yyyy_xxyyz[k];

                g_xz_0_yyyyy_xxzz[k] = -g_xz_0_yyyy_xxzz[k] * ab_y + g_xz_0_yyyy_xxyzz[k];

                g_xz_0_yyyyy_xyyy[k] = -g_xz_0_yyyy_xyyy[k] * ab_y + g_xz_0_yyyy_xyyyy[k];

                g_xz_0_yyyyy_xyyz[k] = -g_xz_0_yyyy_xyyz[k] * ab_y + g_xz_0_yyyy_xyyyz[k];

                g_xz_0_yyyyy_xyzz[k] = -g_xz_0_yyyy_xyzz[k] * ab_y + g_xz_0_yyyy_xyyzz[k];

                g_xz_0_yyyyy_xzzz[k] = -g_xz_0_yyyy_xzzz[k] * ab_y + g_xz_0_yyyy_xyzzz[k];

                g_xz_0_yyyyy_yyyy[k] = -g_xz_0_yyyy_yyyy[k] * ab_y + g_xz_0_yyyy_yyyyy[k];

                g_xz_0_yyyyy_yyyz[k] = -g_xz_0_yyyy_yyyz[k] * ab_y + g_xz_0_yyyy_yyyyz[k];

                g_xz_0_yyyyy_yyzz[k] = -g_xz_0_yyyy_yyzz[k] * ab_y + g_xz_0_yyyy_yyyzz[k];

                g_xz_0_yyyyy_yzzz[k] = -g_xz_0_yyyy_yzzz[k] * ab_y + g_xz_0_yyyy_yyzzz[k];

                g_xz_0_yyyyy_zzzz[k] = -g_xz_0_yyyy_zzzz[k] * ab_y + g_xz_0_yyyy_yzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 870 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 871 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 872 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 873 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 874 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 875 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 876 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 877 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 878 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 879 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 880 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 881 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 882 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 883 * ccomps * dcomps);

            auto g_xz_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyz_xxxx, g_xz_0_yyyyz_xxxy, g_xz_0_yyyyz_xxxz, g_xz_0_yyyyz_xxyy, g_xz_0_yyyyz_xxyz, g_xz_0_yyyyz_xxzz, g_xz_0_yyyyz_xyyy, g_xz_0_yyyyz_xyyz, g_xz_0_yyyyz_xyzz, g_xz_0_yyyyz_xzzz, g_xz_0_yyyyz_yyyy, g_xz_0_yyyyz_yyyz, g_xz_0_yyyyz_yyzz, g_xz_0_yyyyz_yzzz, g_xz_0_yyyyz_zzzz, g_xz_0_yyyz_xxxx, g_xz_0_yyyz_xxxxy, g_xz_0_yyyz_xxxy, g_xz_0_yyyz_xxxyy, g_xz_0_yyyz_xxxyz, g_xz_0_yyyz_xxxz, g_xz_0_yyyz_xxyy, g_xz_0_yyyz_xxyyy, g_xz_0_yyyz_xxyyz, g_xz_0_yyyz_xxyz, g_xz_0_yyyz_xxyzz, g_xz_0_yyyz_xxzz, g_xz_0_yyyz_xyyy, g_xz_0_yyyz_xyyyy, g_xz_0_yyyz_xyyyz, g_xz_0_yyyz_xyyz, g_xz_0_yyyz_xyyzz, g_xz_0_yyyz_xyzz, g_xz_0_yyyz_xyzzz, g_xz_0_yyyz_xzzz, g_xz_0_yyyz_yyyy, g_xz_0_yyyz_yyyyy, g_xz_0_yyyz_yyyyz, g_xz_0_yyyz_yyyz, g_xz_0_yyyz_yyyzz, g_xz_0_yyyz_yyzz, g_xz_0_yyyz_yyzzz, g_xz_0_yyyz_yzzz, g_xz_0_yyyz_yzzzz, g_xz_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyz_xxxx[k] = -g_xz_0_yyyz_xxxx[k] * ab_y + g_xz_0_yyyz_xxxxy[k];

                g_xz_0_yyyyz_xxxy[k] = -g_xz_0_yyyz_xxxy[k] * ab_y + g_xz_0_yyyz_xxxyy[k];

                g_xz_0_yyyyz_xxxz[k] = -g_xz_0_yyyz_xxxz[k] * ab_y + g_xz_0_yyyz_xxxyz[k];

                g_xz_0_yyyyz_xxyy[k] = -g_xz_0_yyyz_xxyy[k] * ab_y + g_xz_0_yyyz_xxyyy[k];

                g_xz_0_yyyyz_xxyz[k] = -g_xz_0_yyyz_xxyz[k] * ab_y + g_xz_0_yyyz_xxyyz[k];

                g_xz_0_yyyyz_xxzz[k] = -g_xz_0_yyyz_xxzz[k] * ab_y + g_xz_0_yyyz_xxyzz[k];

                g_xz_0_yyyyz_xyyy[k] = -g_xz_0_yyyz_xyyy[k] * ab_y + g_xz_0_yyyz_xyyyy[k];

                g_xz_0_yyyyz_xyyz[k] = -g_xz_0_yyyz_xyyz[k] * ab_y + g_xz_0_yyyz_xyyyz[k];

                g_xz_0_yyyyz_xyzz[k] = -g_xz_0_yyyz_xyzz[k] * ab_y + g_xz_0_yyyz_xyyzz[k];

                g_xz_0_yyyyz_xzzz[k] = -g_xz_0_yyyz_xzzz[k] * ab_y + g_xz_0_yyyz_xyzzz[k];

                g_xz_0_yyyyz_yyyy[k] = -g_xz_0_yyyz_yyyy[k] * ab_y + g_xz_0_yyyz_yyyyy[k];

                g_xz_0_yyyyz_yyyz[k] = -g_xz_0_yyyz_yyyz[k] * ab_y + g_xz_0_yyyz_yyyyz[k];

                g_xz_0_yyyyz_yyzz[k] = -g_xz_0_yyyz_yyzz[k] * ab_y + g_xz_0_yyyz_yyyzz[k];

                g_xz_0_yyyyz_yzzz[k] = -g_xz_0_yyyz_yzzz[k] * ab_y + g_xz_0_yyyz_yyzzz[k];

                g_xz_0_yyyyz_zzzz[k] = -g_xz_0_yyyz_zzzz[k] * ab_y + g_xz_0_yyyz_yzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 885 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 886 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 887 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 888 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 889 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 890 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 891 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 892 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 893 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 894 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 895 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 896 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 897 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 898 * ccomps * dcomps);

            auto g_xz_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzz_xxxx, g_xz_0_yyyzz_xxxy, g_xz_0_yyyzz_xxxz, g_xz_0_yyyzz_xxyy, g_xz_0_yyyzz_xxyz, g_xz_0_yyyzz_xxzz, g_xz_0_yyyzz_xyyy, g_xz_0_yyyzz_xyyz, g_xz_0_yyyzz_xyzz, g_xz_0_yyyzz_xzzz, g_xz_0_yyyzz_yyyy, g_xz_0_yyyzz_yyyz, g_xz_0_yyyzz_yyzz, g_xz_0_yyyzz_yzzz, g_xz_0_yyyzz_zzzz, g_xz_0_yyzz_xxxx, g_xz_0_yyzz_xxxxy, g_xz_0_yyzz_xxxy, g_xz_0_yyzz_xxxyy, g_xz_0_yyzz_xxxyz, g_xz_0_yyzz_xxxz, g_xz_0_yyzz_xxyy, g_xz_0_yyzz_xxyyy, g_xz_0_yyzz_xxyyz, g_xz_0_yyzz_xxyz, g_xz_0_yyzz_xxyzz, g_xz_0_yyzz_xxzz, g_xz_0_yyzz_xyyy, g_xz_0_yyzz_xyyyy, g_xz_0_yyzz_xyyyz, g_xz_0_yyzz_xyyz, g_xz_0_yyzz_xyyzz, g_xz_0_yyzz_xyzz, g_xz_0_yyzz_xyzzz, g_xz_0_yyzz_xzzz, g_xz_0_yyzz_yyyy, g_xz_0_yyzz_yyyyy, g_xz_0_yyzz_yyyyz, g_xz_0_yyzz_yyyz, g_xz_0_yyzz_yyyzz, g_xz_0_yyzz_yyzz, g_xz_0_yyzz_yyzzz, g_xz_0_yyzz_yzzz, g_xz_0_yyzz_yzzzz, g_xz_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzz_xxxx[k] = -g_xz_0_yyzz_xxxx[k] * ab_y + g_xz_0_yyzz_xxxxy[k];

                g_xz_0_yyyzz_xxxy[k] = -g_xz_0_yyzz_xxxy[k] * ab_y + g_xz_0_yyzz_xxxyy[k];

                g_xz_0_yyyzz_xxxz[k] = -g_xz_0_yyzz_xxxz[k] * ab_y + g_xz_0_yyzz_xxxyz[k];

                g_xz_0_yyyzz_xxyy[k] = -g_xz_0_yyzz_xxyy[k] * ab_y + g_xz_0_yyzz_xxyyy[k];

                g_xz_0_yyyzz_xxyz[k] = -g_xz_0_yyzz_xxyz[k] * ab_y + g_xz_0_yyzz_xxyyz[k];

                g_xz_0_yyyzz_xxzz[k] = -g_xz_0_yyzz_xxzz[k] * ab_y + g_xz_0_yyzz_xxyzz[k];

                g_xz_0_yyyzz_xyyy[k] = -g_xz_0_yyzz_xyyy[k] * ab_y + g_xz_0_yyzz_xyyyy[k];

                g_xz_0_yyyzz_xyyz[k] = -g_xz_0_yyzz_xyyz[k] * ab_y + g_xz_0_yyzz_xyyyz[k];

                g_xz_0_yyyzz_xyzz[k] = -g_xz_0_yyzz_xyzz[k] * ab_y + g_xz_0_yyzz_xyyzz[k];

                g_xz_0_yyyzz_xzzz[k] = -g_xz_0_yyzz_xzzz[k] * ab_y + g_xz_0_yyzz_xyzzz[k];

                g_xz_0_yyyzz_yyyy[k] = -g_xz_0_yyzz_yyyy[k] * ab_y + g_xz_0_yyzz_yyyyy[k];

                g_xz_0_yyyzz_yyyz[k] = -g_xz_0_yyzz_yyyz[k] * ab_y + g_xz_0_yyzz_yyyyz[k];

                g_xz_0_yyyzz_yyzz[k] = -g_xz_0_yyzz_yyzz[k] * ab_y + g_xz_0_yyzz_yyyzz[k];

                g_xz_0_yyyzz_yzzz[k] = -g_xz_0_yyzz_yzzz[k] * ab_y + g_xz_0_yyzz_yyzzz[k];

                g_xz_0_yyyzz_zzzz[k] = -g_xz_0_yyzz_zzzz[k] * ab_y + g_xz_0_yyzz_yzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 900 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 901 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 902 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 903 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 904 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 905 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 906 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 907 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 908 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 909 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 910 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 911 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 912 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 913 * ccomps * dcomps);

            auto g_xz_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzz_xxxx, g_xz_0_yyzzz_xxxy, g_xz_0_yyzzz_xxxz, g_xz_0_yyzzz_xxyy, g_xz_0_yyzzz_xxyz, g_xz_0_yyzzz_xxzz, g_xz_0_yyzzz_xyyy, g_xz_0_yyzzz_xyyz, g_xz_0_yyzzz_xyzz, g_xz_0_yyzzz_xzzz, g_xz_0_yyzzz_yyyy, g_xz_0_yyzzz_yyyz, g_xz_0_yyzzz_yyzz, g_xz_0_yyzzz_yzzz, g_xz_0_yyzzz_zzzz, g_xz_0_yzzz_xxxx, g_xz_0_yzzz_xxxxy, g_xz_0_yzzz_xxxy, g_xz_0_yzzz_xxxyy, g_xz_0_yzzz_xxxyz, g_xz_0_yzzz_xxxz, g_xz_0_yzzz_xxyy, g_xz_0_yzzz_xxyyy, g_xz_0_yzzz_xxyyz, g_xz_0_yzzz_xxyz, g_xz_0_yzzz_xxyzz, g_xz_0_yzzz_xxzz, g_xz_0_yzzz_xyyy, g_xz_0_yzzz_xyyyy, g_xz_0_yzzz_xyyyz, g_xz_0_yzzz_xyyz, g_xz_0_yzzz_xyyzz, g_xz_0_yzzz_xyzz, g_xz_0_yzzz_xyzzz, g_xz_0_yzzz_xzzz, g_xz_0_yzzz_yyyy, g_xz_0_yzzz_yyyyy, g_xz_0_yzzz_yyyyz, g_xz_0_yzzz_yyyz, g_xz_0_yzzz_yyyzz, g_xz_0_yzzz_yyzz, g_xz_0_yzzz_yyzzz, g_xz_0_yzzz_yzzz, g_xz_0_yzzz_yzzzz, g_xz_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzz_xxxx[k] = -g_xz_0_yzzz_xxxx[k] * ab_y + g_xz_0_yzzz_xxxxy[k];

                g_xz_0_yyzzz_xxxy[k] = -g_xz_0_yzzz_xxxy[k] * ab_y + g_xz_0_yzzz_xxxyy[k];

                g_xz_0_yyzzz_xxxz[k] = -g_xz_0_yzzz_xxxz[k] * ab_y + g_xz_0_yzzz_xxxyz[k];

                g_xz_0_yyzzz_xxyy[k] = -g_xz_0_yzzz_xxyy[k] * ab_y + g_xz_0_yzzz_xxyyy[k];

                g_xz_0_yyzzz_xxyz[k] = -g_xz_0_yzzz_xxyz[k] * ab_y + g_xz_0_yzzz_xxyyz[k];

                g_xz_0_yyzzz_xxzz[k] = -g_xz_0_yzzz_xxzz[k] * ab_y + g_xz_0_yzzz_xxyzz[k];

                g_xz_0_yyzzz_xyyy[k] = -g_xz_0_yzzz_xyyy[k] * ab_y + g_xz_0_yzzz_xyyyy[k];

                g_xz_0_yyzzz_xyyz[k] = -g_xz_0_yzzz_xyyz[k] * ab_y + g_xz_0_yzzz_xyyyz[k];

                g_xz_0_yyzzz_xyzz[k] = -g_xz_0_yzzz_xyzz[k] * ab_y + g_xz_0_yzzz_xyyzz[k];

                g_xz_0_yyzzz_xzzz[k] = -g_xz_0_yzzz_xzzz[k] * ab_y + g_xz_0_yzzz_xyzzz[k];

                g_xz_0_yyzzz_yyyy[k] = -g_xz_0_yzzz_yyyy[k] * ab_y + g_xz_0_yzzz_yyyyy[k];

                g_xz_0_yyzzz_yyyz[k] = -g_xz_0_yzzz_yyyz[k] * ab_y + g_xz_0_yzzz_yyyyz[k];

                g_xz_0_yyzzz_yyzz[k] = -g_xz_0_yzzz_yyzz[k] * ab_y + g_xz_0_yzzz_yyyzz[k];

                g_xz_0_yyzzz_yzzz[k] = -g_xz_0_yzzz_yzzz[k] * ab_y + g_xz_0_yzzz_yyzzz[k];

                g_xz_0_yyzzz_zzzz[k] = -g_xz_0_yzzz_zzzz[k] * ab_y + g_xz_0_yzzz_yzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 915 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 916 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 917 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 918 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 919 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 920 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 921 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 922 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 923 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 924 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 925 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 926 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 927 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 928 * ccomps * dcomps);

            auto g_xz_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzz_xxxx, g_xz_0_yzzzz_xxxy, g_xz_0_yzzzz_xxxz, g_xz_0_yzzzz_xxyy, g_xz_0_yzzzz_xxyz, g_xz_0_yzzzz_xxzz, g_xz_0_yzzzz_xyyy, g_xz_0_yzzzz_xyyz, g_xz_0_yzzzz_xyzz, g_xz_0_yzzzz_xzzz, g_xz_0_yzzzz_yyyy, g_xz_0_yzzzz_yyyz, g_xz_0_yzzzz_yyzz, g_xz_0_yzzzz_yzzz, g_xz_0_yzzzz_zzzz, g_xz_0_zzzz_xxxx, g_xz_0_zzzz_xxxxy, g_xz_0_zzzz_xxxy, g_xz_0_zzzz_xxxyy, g_xz_0_zzzz_xxxyz, g_xz_0_zzzz_xxxz, g_xz_0_zzzz_xxyy, g_xz_0_zzzz_xxyyy, g_xz_0_zzzz_xxyyz, g_xz_0_zzzz_xxyz, g_xz_0_zzzz_xxyzz, g_xz_0_zzzz_xxzz, g_xz_0_zzzz_xyyy, g_xz_0_zzzz_xyyyy, g_xz_0_zzzz_xyyyz, g_xz_0_zzzz_xyyz, g_xz_0_zzzz_xyyzz, g_xz_0_zzzz_xyzz, g_xz_0_zzzz_xyzzz, g_xz_0_zzzz_xzzz, g_xz_0_zzzz_yyyy, g_xz_0_zzzz_yyyyy, g_xz_0_zzzz_yyyyz, g_xz_0_zzzz_yyyz, g_xz_0_zzzz_yyyzz, g_xz_0_zzzz_yyzz, g_xz_0_zzzz_yyzzz, g_xz_0_zzzz_yzzz, g_xz_0_zzzz_yzzzz, g_xz_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzz_xxxx[k] = -g_xz_0_zzzz_xxxx[k] * ab_y + g_xz_0_zzzz_xxxxy[k];

                g_xz_0_yzzzz_xxxy[k] = -g_xz_0_zzzz_xxxy[k] * ab_y + g_xz_0_zzzz_xxxyy[k];

                g_xz_0_yzzzz_xxxz[k] = -g_xz_0_zzzz_xxxz[k] * ab_y + g_xz_0_zzzz_xxxyz[k];

                g_xz_0_yzzzz_xxyy[k] = -g_xz_0_zzzz_xxyy[k] * ab_y + g_xz_0_zzzz_xxyyy[k];

                g_xz_0_yzzzz_xxyz[k] = -g_xz_0_zzzz_xxyz[k] * ab_y + g_xz_0_zzzz_xxyyz[k];

                g_xz_0_yzzzz_xxzz[k] = -g_xz_0_zzzz_xxzz[k] * ab_y + g_xz_0_zzzz_xxyzz[k];

                g_xz_0_yzzzz_xyyy[k] = -g_xz_0_zzzz_xyyy[k] * ab_y + g_xz_0_zzzz_xyyyy[k];

                g_xz_0_yzzzz_xyyz[k] = -g_xz_0_zzzz_xyyz[k] * ab_y + g_xz_0_zzzz_xyyyz[k];

                g_xz_0_yzzzz_xyzz[k] = -g_xz_0_zzzz_xyzz[k] * ab_y + g_xz_0_zzzz_xyyzz[k];

                g_xz_0_yzzzz_xzzz[k] = -g_xz_0_zzzz_xzzz[k] * ab_y + g_xz_0_zzzz_xyzzz[k];

                g_xz_0_yzzzz_yyyy[k] = -g_xz_0_zzzz_yyyy[k] * ab_y + g_xz_0_zzzz_yyyyy[k];

                g_xz_0_yzzzz_yyyz[k] = -g_xz_0_zzzz_yyyz[k] * ab_y + g_xz_0_zzzz_yyyyz[k];

                g_xz_0_yzzzz_yyzz[k] = -g_xz_0_zzzz_yyzz[k] * ab_y + g_xz_0_zzzz_yyyzz[k];

                g_xz_0_yzzzz_yzzz[k] = -g_xz_0_zzzz_yzzz[k] * ab_y + g_xz_0_zzzz_yyzzz[k];

                g_xz_0_yzzzz_zzzz[k] = -g_xz_0_zzzz_zzzz[k] * ab_y + g_xz_0_zzzz_yzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 930 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 931 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 932 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 933 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 934 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 935 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 936 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 937 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 938 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 939 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 940 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 941 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 942 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 943 * ccomps * dcomps);

            auto g_xz_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzzz, g_xz_0_zzzz_xxxx, g_xz_0_zzzz_xxxxz, g_xz_0_zzzz_xxxy, g_xz_0_zzzz_xxxyz, g_xz_0_zzzz_xxxz, g_xz_0_zzzz_xxxzz, g_xz_0_zzzz_xxyy, g_xz_0_zzzz_xxyyz, g_xz_0_zzzz_xxyz, g_xz_0_zzzz_xxyzz, g_xz_0_zzzz_xxzz, g_xz_0_zzzz_xxzzz, g_xz_0_zzzz_xyyy, g_xz_0_zzzz_xyyyz, g_xz_0_zzzz_xyyz, g_xz_0_zzzz_xyyzz, g_xz_0_zzzz_xyzz, g_xz_0_zzzz_xyzzz, g_xz_0_zzzz_xzzz, g_xz_0_zzzz_xzzzz, g_xz_0_zzzz_yyyy, g_xz_0_zzzz_yyyyz, g_xz_0_zzzz_yyyz, g_xz_0_zzzz_yyyzz, g_xz_0_zzzz_yyzz, g_xz_0_zzzz_yyzzz, g_xz_0_zzzz_yzzz, g_xz_0_zzzz_yzzzz, g_xz_0_zzzz_zzzz, g_xz_0_zzzz_zzzzz, g_xz_0_zzzzz_xxxx, g_xz_0_zzzzz_xxxy, g_xz_0_zzzzz_xxxz, g_xz_0_zzzzz_xxyy, g_xz_0_zzzzz_xxyz, g_xz_0_zzzzz_xxzz, g_xz_0_zzzzz_xyyy, g_xz_0_zzzzz_xyyz, g_xz_0_zzzzz_xyzz, g_xz_0_zzzzz_xzzz, g_xz_0_zzzzz_yyyy, g_xz_0_zzzzz_yyyz, g_xz_0_zzzzz_yyzz, g_xz_0_zzzzz_yzzz, g_xz_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzz_xxxx[k] = -g_x_0_zzzz_xxxx[k] - g_xz_0_zzzz_xxxx[k] * ab_z + g_xz_0_zzzz_xxxxz[k];

                g_xz_0_zzzzz_xxxy[k] = -g_x_0_zzzz_xxxy[k] - g_xz_0_zzzz_xxxy[k] * ab_z + g_xz_0_zzzz_xxxyz[k];

                g_xz_0_zzzzz_xxxz[k] = -g_x_0_zzzz_xxxz[k] - g_xz_0_zzzz_xxxz[k] * ab_z + g_xz_0_zzzz_xxxzz[k];

                g_xz_0_zzzzz_xxyy[k] = -g_x_0_zzzz_xxyy[k] - g_xz_0_zzzz_xxyy[k] * ab_z + g_xz_0_zzzz_xxyyz[k];

                g_xz_0_zzzzz_xxyz[k] = -g_x_0_zzzz_xxyz[k] - g_xz_0_zzzz_xxyz[k] * ab_z + g_xz_0_zzzz_xxyzz[k];

                g_xz_0_zzzzz_xxzz[k] = -g_x_0_zzzz_xxzz[k] - g_xz_0_zzzz_xxzz[k] * ab_z + g_xz_0_zzzz_xxzzz[k];

                g_xz_0_zzzzz_xyyy[k] = -g_x_0_zzzz_xyyy[k] - g_xz_0_zzzz_xyyy[k] * ab_z + g_xz_0_zzzz_xyyyz[k];

                g_xz_0_zzzzz_xyyz[k] = -g_x_0_zzzz_xyyz[k] - g_xz_0_zzzz_xyyz[k] * ab_z + g_xz_0_zzzz_xyyzz[k];

                g_xz_0_zzzzz_xyzz[k] = -g_x_0_zzzz_xyzz[k] - g_xz_0_zzzz_xyzz[k] * ab_z + g_xz_0_zzzz_xyzzz[k];

                g_xz_0_zzzzz_xzzz[k] = -g_x_0_zzzz_xzzz[k] - g_xz_0_zzzz_xzzz[k] * ab_z + g_xz_0_zzzz_xzzzz[k];

                g_xz_0_zzzzz_yyyy[k] = -g_x_0_zzzz_yyyy[k] - g_xz_0_zzzz_yyyy[k] * ab_z + g_xz_0_zzzz_yyyyz[k];

                g_xz_0_zzzzz_yyyz[k] = -g_x_0_zzzz_yyyz[k] - g_xz_0_zzzz_yyyz[k] * ab_z + g_xz_0_zzzz_yyyzz[k];

                g_xz_0_zzzzz_yyzz[k] = -g_x_0_zzzz_yyzz[k] - g_xz_0_zzzz_yyzz[k] * ab_z + g_xz_0_zzzz_yyzzz[k];

                g_xz_0_zzzzz_yzzz[k] = -g_x_0_zzzz_yzzz[k] - g_xz_0_zzzz_yzzz[k] * ab_z + g_xz_0_zzzz_yzzzz[k];

                g_xz_0_zzzzz_zzzz[k] = -g_x_0_zzzz_zzzz[k] - g_xz_0_zzzz_zzzz[k] * ab_z + g_xz_0_zzzz_zzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 945 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 946 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 947 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 948 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 949 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 950 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 951 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 952 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 953 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 954 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 955 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 956 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 957 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 958 * ccomps * dcomps);

            auto g_yy_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxx_xxxx, g_yy_0_xxxx_xxxxx, g_yy_0_xxxx_xxxxy, g_yy_0_xxxx_xxxxz, g_yy_0_xxxx_xxxy, g_yy_0_xxxx_xxxyy, g_yy_0_xxxx_xxxyz, g_yy_0_xxxx_xxxz, g_yy_0_xxxx_xxxzz, g_yy_0_xxxx_xxyy, g_yy_0_xxxx_xxyyy, g_yy_0_xxxx_xxyyz, g_yy_0_xxxx_xxyz, g_yy_0_xxxx_xxyzz, g_yy_0_xxxx_xxzz, g_yy_0_xxxx_xxzzz, g_yy_0_xxxx_xyyy, g_yy_0_xxxx_xyyyy, g_yy_0_xxxx_xyyyz, g_yy_0_xxxx_xyyz, g_yy_0_xxxx_xyyzz, g_yy_0_xxxx_xyzz, g_yy_0_xxxx_xyzzz, g_yy_0_xxxx_xzzz, g_yy_0_xxxx_xzzzz, g_yy_0_xxxx_yyyy, g_yy_0_xxxx_yyyz, g_yy_0_xxxx_yyzz, g_yy_0_xxxx_yzzz, g_yy_0_xxxx_zzzz, g_yy_0_xxxxx_xxxx, g_yy_0_xxxxx_xxxy, g_yy_0_xxxxx_xxxz, g_yy_0_xxxxx_xxyy, g_yy_0_xxxxx_xxyz, g_yy_0_xxxxx_xxzz, g_yy_0_xxxxx_xyyy, g_yy_0_xxxxx_xyyz, g_yy_0_xxxxx_xyzz, g_yy_0_xxxxx_xzzz, g_yy_0_xxxxx_yyyy, g_yy_0_xxxxx_yyyz, g_yy_0_xxxxx_yyzz, g_yy_0_xxxxx_yzzz, g_yy_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxx_xxxx[k] = -g_yy_0_xxxx_xxxx[k] * ab_x + g_yy_0_xxxx_xxxxx[k];

                g_yy_0_xxxxx_xxxy[k] = -g_yy_0_xxxx_xxxy[k] * ab_x + g_yy_0_xxxx_xxxxy[k];

                g_yy_0_xxxxx_xxxz[k] = -g_yy_0_xxxx_xxxz[k] * ab_x + g_yy_0_xxxx_xxxxz[k];

                g_yy_0_xxxxx_xxyy[k] = -g_yy_0_xxxx_xxyy[k] * ab_x + g_yy_0_xxxx_xxxyy[k];

                g_yy_0_xxxxx_xxyz[k] = -g_yy_0_xxxx_xxyz[k] * ab_x + g_yy_0_xxxx_xxxyz[k];

                g_yy_0_xxxxx_xxzz[k] = -g_yy_0_xxxx_xxzz[k] * ab_x + g_yy_0_xxxx_xxxzz[k];

                g_yy_0_xxxxx_xyyy[k] = -g_yy_0_xxxx_xyyy[k] * ab_x + g_yy_0_xxxx_xxyyy[k];

                g_yy_0_xxxxx_xyyz[k] = -g_yy_0_xxxx_xyyz[k] * ab_x + g_yy_0_xxxx_xxyyz[k];

                g_yy_0_xxxxx_xyzz[k] = -g_yy_0_xxxx_xyzz[k] * ab_x + g_yy_0_xxxx_xxyzz[k];

                g_yy_0_xxxxx_xzzz[k] = -g_yy_0_xxxx_xzzz[k] * ab_x + g_yy_0_xxxx_xxzzz[k];

                g_yy_0_xxxxx_yyyy[k] = -g_yy_0_xxxx_yyyy[k] * ab_x + g_yy_0_xxxx_xyyyy[k];

                g_yy_0_xxxxx_yyyz[k] = -g_yy_0_xxxx_yyyz[k] * ab_x + g_yy_0_xxxx_xyyyz[k];

                g_yy_0_xxxxx_yyzz[k] = -g_yy_0_xxxx_yyzz[k] * ab_x + g_yy_0_xxxx_xyyzz[k];

                g_yy_0_xxxxx_yzzz[k] = -g_yy_0_xxxx_yzzz[k] * ab_x + g_yy_0_xxxx_xyzzz[k];

                g_yy_0_xxxxx_zzzz[k] = -g_yy_0_xxxx_zzzz[k] * ab_x + g_yy_0_xxxx_xzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 960 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 961 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 962 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 963 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 964 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 965 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 966 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 967 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 968 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 969 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 970 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 971 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 972 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 973 * ccomps * dcomps);

            auto g_yy_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxy_xxxx, g_yy_0_xxxxy_xxxy, g_yy_0_xxxxy_xxxz, g_yy_0_xxxxy_xxyy, g_yy_0_xxxxy_xxyz, g_yy_0_xxxxy_xxzz, g_yy_0_xxxxy_xyyy, g_yy_0_xxxxy_xyyz, g_yy_0_xxxxy_xyzz, g_yy_0_xxxxy_xzzz, g_yy_0_xxxxy_yyyy, g_yy_0_xxxxy_yyyz, g_yy_0_xxxxy_yyzz, g_yy_0_xxxxy_yzzz, g_yy_0_xxxxy_zzzz, g_yy_0_xxxy_xxxx, g_yy_0_xxxy_xxxxx, g_yy_0_xxxy_xxxxy, g_yy_0_xxxy_xxxxz, g_yy_0_xxxy_xxxy, g_yy_0_xxxy_xxxyy, g_yy_0_xxxy_xxxyz, g_yy_0_xxxy_xxxz, g_yy_0_xxxy_xxxzz, g_yy_0_xxxy_xxyy, g_yy_0_xxxy_xxyyy, g_yy_0_xxxy_xxyyz, g_yy_0_xxxy_xxyz, g_yy_0_xxxy_xxyzz, g_yy_0_xxxy_xxzz, g_yy_0_xxxy_xxzzz, g_yy_0_xxxy_xyyy, g_yy_0_xxxy_xyyyy, g_yy_0_xxxy_xyyyz, g_yy_0_xxxy_xyyz, g_yy_0_xxxy_xyyzz, g_yy_0_xxxy_xyzz, g_yy_0_xxxy_xyzzz, g_yy_0_xxxy_xzzz, g_yy_0_xxxy_xzzzz, g_yy_0_xxxy_yyyy, g_yy_0_xxxy_yyyz, g_yy_0_xxxy_yyzz, g_yy_0_xxxy_yzzz, g_yy_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxy_xxxx[k] = -g_yy_0_xxxy_xxxx[k] * ab_x + g_yy_0_xxxy_xxxxx[k];

                g_yy_0_xxxxy_xxxy[k] = -g_yy_0_xxxy_xxxy[k] * ab_x + g_yy_0_xxxy_xxxxy[k];

                g_yy_0_xxxxy_xxxz[k] = -g_yy_0_xxxy_xxxz[k] * ab_x + g_yy_0_xxxy_xxxxz[k];

                g_yy_0_xxxxy_xxyy[k] = -g_yy_0_xxxy_xxyy[k] * ab_x + g_yy_0_xxxy_xxxyy[k];

                g_yy_0_xxxxy_xxyz[k] = -g_yy_0_xxxy_xxyz[k] * ab_x + g_yy_0_xxxy_xxxyz[k];

                g_yy_0_xxxxy_xxzz[k] = -g_yy_0_xxxy_xxzz[k] * ab_x + g_yy_0_xxxy_xxxzz[k];

                g_yy_0_xxxxy_xyyy[k] = -g_yy_0_xxxy_xyyy[k] * ab_x + g_yy_0_xxxy_xxyyy[k];

                g_yy_0_xxxxy_xyyz[k] = -g_yy_0_xxxy_xyyz[k] * ab_x + g_yy_0_xxxy_xxyyz[k];

                g_yy_0_xxxxy_xyzz[k] = -g_yy_0_xxxy_xyzz[k] * ab_x + g_yy_0_xxxy_xxyzz[k];

                g_yy_0_xxxxy_xzzz[k] = -g_yy_0_xxxy_xzzz[k] * ab_x + g_yy_0_xxxy_xxzzz[k];

                g_yy_0_xxxxy_yyyy[k] = -g_yy_0_xxxy_yyyy[k] * ab_x + g_yy_0_xxxy_xyyyy[k];

                g_yy_0_xxxxy_yyyz[k] = -g_yy_0_xxxy_yyyz[k] * ab_x + g_yy_0_xxxy_xyyyz[k];

                g_yy_0_xxxxy_yyzz[k] = -g_yy_0_xxxy_yyzz[k] * ab_x + g_yy_0_xxxy_xyyzz[k];

                g_yy_0_xxxxy_yzzz[k] = -g_yy_0_xxxy_yzzz[k] * ab_x + g_yy_0_xxxy_xyzzz[k];

                g_yy_0_xxxxy_zzzz[k] = -g_yy_0_xxxy_zzzz[k] * ab_x + g_yy_0_xxxy_xzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 975 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 976 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 977 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 978 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 979 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 980 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 981 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 982 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 983 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 984 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 985 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 986 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 987 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 988 * ccomps * dcomps);

            auto g_yy_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxz_xxxx, g_yy_0_xxxxz_xxxy, g_yy_0_xxxxz_xxxz, g_yy_0_xxxxz_xxyy, g_yy_0_xxxxz_xxyz, g_yy_0_xxxxz_xxzz, g_yy_0_xxxxz_xyyy, g_yy_0_xxxxz_xyyz, g_yy_0_xxxxz_xyzz, g_yy_0_xxxxz_xzzz, g_yy_0_xxxxz_yyyy, g_yy_0_xxxxz_yyyz, g_yy_0_xxxxz_yyzz, g_yy_0_xxxxz_yzzz, g_yy_0_xxxxz_zzzz, g_yy_0_xxxz_xxxx, g_yy_0_xxxz_xxxxx, g_yy_0_xxxz_xxxxy, g_yy_0_xxxz_xxxxz, g_yy_0_xxxz_xxxy, g_yy_0_xxxz_xxxyy, g_yy_0_xxxz_xxxyz, g_yy_0_xxxz_xxxz, g_yy_0_xxxz_xxxzz, g_yy_0_xxxz_xxyy, g_yy_0_xxxz_xxyyy, g_yy_0_xxxz_xxyyz, g_yy_0_xxxz_xxyz, g_yy_0_xxxz_xxyzz, g_yy_0_xxxz_xxzz, g_yy_0_xxxz_xxzzz, g_yy_0_xxxz_xyyy, g_yy_0_xxxz_xyyyy, g_yy_0_xxxz_xyyyz, g_yy_0_xxxz_xyyz, g_yy_0_xxxz_xyyzz, g_yy_0_xxxz_xyzz, g_yy_0_xxxz_xyzzz, g_yy_0_xxxz_xzzz, g_yy_0_xxxz_xzzzz, g_yy_0_xxxz_yyyy, g_yy_0_xxxz_yyyz, g_yy_0_xxxz_yyzz, g_yy_0_xxxz_yzzz, g_yy_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxz_xxxx[k] = -g_yy_0_xxxz_xxxx[k] * ab_x + g_yy_0_xxxz_xxxxx[k];

                g_yy_0_xxxxz_xxxy[k] = -g_yy_0_xxxz_xxxy[k] * ab_x + g_yy_0_xxxz_xxxxy[k];

                g_yy_0_xxxxz_xxxz[k] = -g_yy_0_xxxz_xxxz[k] * ab_x + g_yy_0_xxxz_xxxxz[k];

                g_yy_0_xxxxz_xxyy[k] = -g_yy_0_xxxz_xxyy[k] * ab_x + g_yy_0_xxxz_xxxyy[k];

                g_yy_0_xxxxz_xxyz[k] = -g_yy_0_xxxz_xxyz[k] * ab_x + g_yy_0_xxxz_xxxyz[k];

                g_yy_0_xxxxz_xxzz[k] = -g_yy_0_xxxz_xxzz[k] * ab_x + g_yy_0_xxxz_xxxzz[k];

                g_yy_0_xxxxz_xyyy[k] = -g_yy_0_xxxz_xyyy[k] * ab_x + g_yy_0_xxxz_xxyyy[k];

                g_yy_0_xxxxz_xyyz[k] = -g_yy_0_xxxz_xyyz[k] * ab_x + g_yy_0_xxxz_xxyyz[k];

                g_yy_0_xxxxz_xyzz[k] = -g_yy_0_xxxz_xyzz[k] * ab_x + g_yy_0_xxxz_xxyzz[k];

                g_yy_0_xxxxz_xzzz[k] = -g_yy_0_xxxz_xzzz[k] * ab_x + g_yy_0_xxxz_xxzzz[k];

                g_yy_0_xxxxz_yyyy[k] = -g_yy_0_xxxz_yyyy[k] * ab_x + g_yy_0_xxxz_xyyyy[k];

                g_yy_0_xxxxz_yyyz[k] = -g_yy_0_xxxz_yyyz[k] * ab_x + g_yy_0_xxxz_xyyyz[k];

                g_yy_0_xxxxz_yyzz[k] = -g_yy_0_xxxz_yyzz[k] * ab_x + g_yy_0_xxxz_xyyzz[k];

                g_yy_0_xxxxz_yzzz[k] = -g_yy_0_xxxz_yzzz[k] * ab_x + g_yy_0_xxxz_xyzzz[k];

                g_yy_0_xxxxz_zzzz[k] = -g_yy_0_xxxz_zzzz[k] * ab_x + g_yy_0_xxxz_xzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 990 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 991 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 992 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 993 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 994 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 995 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 996 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 997 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 998 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 999 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 1000 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 1001 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 1002 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 1003 * ccomps * dcomps);

            auto g_yy_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyy_xxxx, g_yy_0_xxxyy_xxxy, g_yy_0_xxxyy_xxxz, g_yy_0_xxxyy_xxyy, g_yy_0_xxxyy_xxyz, g_yy_0_xxxyy_xxzz, g_yy_0_xxxyy_xyyy, g_yy_0_xxxyy_xyyz, g_yy_0_xxxyy_xyzz, g_yy_0_xxxyy_xzzz, g_yy_0_xxxyy_yyyy, g_yy_0_xxxyy_yyyz, g_yy_0_xxxyy_yyzz, g_yy_0_xxxyy_yzzz, g_yy_0_xxxyy_zzzz, g_yy_0_xxyy_xxxx, g_yy_0_xxyy_xxxxx, g_yy_0_xxyy_xxxxy, g_yy_0_xxyy_xxxxz, g_yy_0_xxyy_xxxy, g_yy_0_xxyy_xxxyy, g_yy_0_xxyy_xxxyz, g_yy_0_xxyy_xxxz, g_yy_0_xxyy_xxxzz, g_yy_0_xxyy_xxyy, g_yy_0_xxyy_xxyyy, g_yy_0_xxyy_xxyyz, g_yy_0_xxyy_xxyz, g_yy_0_xxyy_xxyzz, g_yy_0_xxyy_xxzz, g_yy_0_xxyy_xxzzz, g_yy_0_xxyy_xyyy, g_yy_0_xxyy_xyyyy, g_yy_0_xxyy_xyyyz, g_yy_0_xxyy_xyyz, g_yy_0_xxyy_xyyzz, g_yy_0_xxyy_xyzz, g_yy_0_xxyy_xyzzz, g_yy_0_xxyy_xzzz, g_yy_0_xxyy_xzzzz, g_yy_0_xxyy_yyyy, g_yy_0_xxyy_yyyz, g_yy_0_xxyy_yyzz, g_yy_0_xxyy_yzzz, g_yy_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyy_xxxx[k] = -g_yy_0_xxyy_xxxx[k] * ab_x + g_yy_0_xxyy_xxxxx[k];

                g_yy_0_xxxyy_xxxy[k] = -g_yy_0_xxyy_xxxy[k] * ab_x + g_yy_0_xxyy_xxxxy[k];

                g_yy_0_xxxyy_xxxz[k] = -g_yy_0_xxyy_xxxz[k] * ab_x + g_yy_0_xxyy_xxxxz[k];

                g_yy_0_xxxyy_xxyy[k] = -g_yy_0_xxyy_xxyy[k] * ab_x + g_yy_0_xxyy_xxxyy[k];

                g_yy_0_xxxyy_xxyz[k] = -g_yy_0_xxyy_xxyz[k] * ab_x + g_yy_0_xxyy_xxxyz[k];

                g_yy_0_xxxyy_xxzz[k] = -g_yy_0_xxyy_xxzz[k] * ab_x + g_yy_0_xxyy_xxxzz[k];

                g_yy_0_xxxyy_xyyy[k] = -g_yy_0_xxyy_xyyy[k] * ab_x + g_yy_0_xxyy_xxyyy[k];

                g_yy_0_xxxyy_xyyz[k] = -g_yy_0_xxyy_xyyz[k] * ab_x + g_yy_0_xxyy_xxyyz[k];

                g_yy_0_xxxyy_xyzz[k] = -g_yy_0_xxyy_xyzz[k] * ab_x + g_yy_0_xxyy_xxyzz[k];

                g_yy_0_xxxyy_xzzz[k] = -g_yy_0_xxyy_xzzz[k] * ab_x + g_yy_0_xxyy_xxzzz[k];

                g_yy_0_xxxyy_yyyy[k] = -g_yy_0_xxyy_yyyy[k] * ab_x + g_yy_0_xxyy_xyyyy[k];

                g_yy_0_xxxyy_yyyz[k] = -g_yy_0_xxyy_yyyz[k] * ab_x + g_yy_0_xxyy_xyyyz[k];

                g_yy_0_xxxyy_yyzz[k] = -g_yy_0_xxyy_yyzz[k] * ab_x + g_yy_0_xxyy_xyyzz[k];

                g_yy_0_xxxyy_yzzz[k] = -g_yy_0_xxyy_yzzz[k] * ab_x + g_yy_0_xxyy_xyzzz[k];

                g_yy_0_xxxyy_zzzz[k] = -g_yy_0_xxyy_zzzz[k] * ab_x + g_yy_0_xxyy_xzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 1005 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 1006 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 1007 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 1008 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 1009 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 1010 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 1011 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 1012 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 1013 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 1014 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 1015 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 1016 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 1017 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 1018 * ccomps * dcomps);

            auto g_yy_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyz_xxxx, g_yy_0_xxxyz_xxxy, g_yy_0_xxxyz_xxxz, g_yy_0_xxxyz_xxyy, g_yy_0_xxxyz_xxyz, g_yy_0_xxxyz_xxzz, g_yy_0_xxxyz_xyyy, g_yy_0_xxxyz_xyyz, g_yy_0_xxxyz_xyzz, g_yy_0_xxxyz_xzzz, g_yy_0_xxxyz_yyyy, g_yy_0_xxxyz_yyyz, g_yy_0_xxxyz_yyzz, g_yy_0_xxxyz_yzzz, g_yy_0_xxxyz_zzzz, g_yy_0_xxyz_xxxx, g_yy_0_xxyz_xxxxx, g_yy_0_xxyz_xxxxy, g_yy_0_xxyz_xxxxz, g_yy_0_xxyz_xxxy, g_yy_0_xxyz_xxxyy, g_yy_0_xxyz_xxxyz, g_yy_0_xxyz_xxxz, g_yy_0_xxyz_xxxzz, g_yy_0_xxyz_xxyy, g_yy_0_xxyz_xxyyy, g_yy_0_xxyz_xxyyz, g_yy_0_xxyz_xxyz, g_yy_0_xxyz_xxyzz, g_yy_0_xxyz_xxzz, g_yy_0_xxyz_xxzzz, g_yy_0_xxyz_xyyy, g_yy_0_xxyz_xyyyy, g_yy_0_xxyz_xyyyz, g_yy_0_xxyz_xyyz, g_yy_0_xxyz_xyyzz, g_yy_0_xxyz_xyzz, g_yy_0_xxyz_xyzzz, g_yy_0_xxyz_xzzz, g_yy_0_xxyz_xzzzz, g_yy_0_xxyz_yyyy, g_yy_0_xxyz_yyyz, g_yy_0_xxyz_yyzz, g_yy_0_xxyz_yzzz, g_yy_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyz_xxxx[k] = -g_yy_0_xxyz_xxxx[k] * ab_x + g_yy_0_xxyz_xxxxx[k];

                g_yy_0_xxxyz_xxxy[k] = -g_yy_0_xxyz_xxxy[k] * ab_x + g_yy_0_xxyz_xxxxy[k];

                g_yy_0_xxxyz_xxxz[k] = -g_yy_0_xxyz_xxxz[k] * ab_x + g_yy_0_xxyz_xxxxz[k];

                g_yy_0_xxxyz_xxyy[k] = -g_yy_0_xxyz_xxyy[k] * ab_x + g_yy_0_xxyz_xxxyy[k];

                g_yy_0_xxxyz_xxyz[k] = -g_yy_0_xxyz_xxyz[k] * ab_x + g_yy_0_xxyz_xxxyz[k];

                g_yy_0_xxxyz_xxzz[k] = -g_yy_0_xxyz_xxzz[k] * ab_x + g_yy_0_xxyz_xxxzz[k];

                g_yy_0_xxxyz_xyyy[k] = -g_yy_0_xxyz_xyyy[k] * ab_x + g_yy_0_xxyz_xxyyy[k];

                g_yy_0_xxxyz_xyyz[k] = -g_yy_0_xxyz_xyyz[k] * ab_x + g_yy_0_xxyz_xxyyz[k];

                g_yy_0_xxxyz_xyzz[k] = -g_yy_0_xxyz_xyzz[k] * ab_x + g_yy_0_xxyz_xxyzz[k];

                g_yy_0_xxxyz_xzzz[k] = -g_yy_0_xxyz_xzzz[k] * ab_x + g_yy_0_xxyz_xxzzz[k];

                g_yy_0_xxxyz_yyyy[k] = -g_yy_0_xxyz_yyyy[k] * ab_x + g_yy_0_xxyz_xyyyy[k];

                g_yy_0_xxxyz_yyyz[k] = -g_yy_0_xxyz_yyyz[k] * ab_x + g_yy_0_xxyz_xyyyz[k];

                g_yy_0_xxxyz_yyzz[k] = -g_yy_0_xxyz_yyzz[k] * ab_x + g_yy_0_xxyz_xyyzz[k];

                g_yy_0_xxxyz_yzzz[k] = -g_yy_0_xxyz_yzzz[k] * ab_x + g_yy_0_xxyz_xyzzz[k];

                g_yy_0_xxxyz_zzzz[k] = -g_yy_0_xxyz_zzzz[k] * ab_x + g_yy_0_xxyz_xzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 1020 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 1021 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 1022 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 1023 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 1024 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 1025 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 1026 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 1027 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 1028 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 1029 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 1030 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 1031 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 1032 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 1033 * ccomps * dcomps);

            auto g_yy_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzz_xxxx, g_yy_0_xxxzz_xxxy, g_yy_0_xxxzz_xxxz, g_yy_0_xxxzz_xxyy, g_yy_0_xxxzz_xxyz, g_yy_0_xxxzz_xxzz, g_yy_0_xxxzz_xyyy, g_yy_0_xxxzz_xyyz, g_yy_0_xxxzz_xyzz, g_yy_0_xxxzz_xzzz, g_yy_0_xxxzz_yyyy, g_yy_0_xxxzz_yyyz, g_yy_0_xxxzz_yyzz, g_yy_0_xxxzz_yzzz, g_yy_0_xxxzz_zzzz, g_yy_0_xxzz_xxxx, g_yy_0_xxzz_xxxxx, g_yy_0_xxzz_xxxxy, g_yy_0_xxzz_xxxxz, g_yy_0_xxzz_xxxy, g_yy_0_xxzz_xxxyy, g_yy_0_xxzz_xxxyz, g_yy_0_xxzz_xxxz, g_yy_0_xxzz_xxxzz, g_yy_0_xxzz_xxyy, g_yy_0_xxzz_xxyyy, g_yy_0_xxzz_xxyyz, g_yy_0_xxzz_xxyz, g_yy_0_xxzz_xxyzz, g_yy_0_xxzz_xxzz, g_yy_0_xxzz_xxzzz, g_yy_0_xxzz_xyyy, g_yy_0_xxzz_xyyyy, g_yy_0_xxzz_xyyyz, g_yy_0_xxzz_xyyz, g_yy_0_xxzz_xyyzz, g_yy_0_xxzz_xyzz, g_yy_0_xxzz_xyzzz, g_yy_0_xxzz_xzzz, g_yy_0_xxzz_xzzzz, g_yy_0_xxzz_yyyy, g_yy_0_xxzz_yyyz, g_yy_0_xxzz_yyzz, g_yy_0_xxzz_yzzz, g_yy_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzz_xxxx[k] = -g_yy_0_xxzz_xxxx[k] * ab_x + g_yy_0_xxzz_xxxxx[k];

                g_yy_0_xxxzz_xxxy[k] = -g_yy_0_xxzz_xxxy[k] * ab_x + g_yy_0_xxzz_xxxxy[k];

                g_yy_0_xxxzz_xxxz[k] = -g_yy_0_xxzz_xxxz[k] * ab_x + g_yy_0_xxzz_xxxxz[k];

                g_yy_0_xxxzz_xxyy[k] = -g_yy_0_xxzz_xxyy[k] * ab_x + g_yy_0_xxzz_xxxyy[k];

                g_yy_0_xxxzz_xxyz[k] = -g_yy_0_xxzz_xxyz[k] * ab_x + g_yy_0_xxzz_xxxyz[k];

                g_yy_0_xxxzz_xxzz[k] = -g_yy_0_xxzz_xxzz[k] * ab_x + g_yy_0_xxzz_xxxzz[k];

                g_yy_0_xxxzz_xyyy[k] = -g_yy_0_xxzz_xyyy[k] * ab_x + g_yy_0_xxzz_xxyyy[k];

                g_yy_0_xxxzz_xyyz[k] = -g_yy_0_xxzz_xyyz[k] * ab_x + g_yy_0_xxzz_xxyyz[k];

                g_yy_0_xxxzz_xyzz[k] = -g_yy_0_xxzz_xyzz[k] * ab_x + g_yy_0_xxzz_xxyzz[k];

                g_yy_0_xxxzz_xzzz[k] = -g_yy_0_xxzz_xzzz[k] * ab_x + g_yy_0_xxzz_xxzzz[k];

                g_yy_0_xxxzz_yyyy[k] = -g_yy_0_xxzz_yyyy[k] * ab_x + g_yy_0_xxzz_xyyyy[k];

                g_yy_0_xxxzz_yyyz[k] = -g_yy_0_xxzz_yyyz[k] * ab_x + g_yy_0_xxzz_xyyyz[k];

                g_yy_0_xxxzz_yyzz[k] = -g_yy_0_xxzz_yyzz[k] * ab_x + g_yy_0_xxzz_xyyzz[k];

                g_yy_0_xxxzz_yzzz[k] = -g_yy_0_xxzz_yzzz[k] * ab_x + g_yy_0_xxzz_xyzzz[k];

                g_yy_0_xxxzz_zzzz[k] = -g_yy_0_xxzz_zzzz[k] * ab_x + g_yy_0_xxzz_xzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 1035 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 1036 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 1037 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 1038 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 1039 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 1040 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 1041 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 1042 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 1043 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 1044 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 1045 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 1046 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 1047 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 1048 * ccomps * dcomps);

            auto g_yy_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyy_xxxx, g_yy_0_xxyyy_xxxy, g_yy_0_xxyyy_xxxz, g_yy_0_xxyyy_xxyy, g_yy_0_xxyyy_xxyz, g_yy_0_xxyyy_xxzz, g_yy_0_xxyyy_xyyy, g_yy_0_xxyyy_xyyz, g_yy_0_xxyyy_xyzz, g_yy_0_xxyyy_xzzz, g_yy_0_xxyyy_yyyy, g_yy_0_xxyyy_yyyz, g_yy_0_xxyyy_yyzz, g_yy_0_xxyyy_yzzz, g_yy_0_xxyyy_zzzz, g_yy_0_xyyy_xxxx, g_yy_0_xyyy_xxxxx, g_yy_0_xyyy_xxxxy, g_yy_0_xyyy_xxxxz, g_yy_0_xyyy_xxxy, g_yy_0_xyyy_xxxyy, g_yy_0_xyyy_xxxyz, g_yy_0_xyyy_xxxz, g_yy_0_xyyy_xxxzz, g_yy_0_xyyy_xxyy, g_yy_0_xyyy_xxyyy, g_yy_0_xyyy_xxyyz, g_yy_0_xyyy_xxyz, g_yy_0_xyyy_xxyzz, g_yy_0_xyyy_xxzz, g_yy_0_xyyy_xxzzz, g_yy_0_xyyy_xyyy, g_yy_0_xyyy_xyyyy, g_yy_0_xyyy_xyyyz, g_yy_0_xyyy_xyyz, g_yy_0_xyyy_xyyzz, g_yy_0_xyyy_xyzz, g_yy_0_xyyy_xyzzz, g_yy_0_xyyy_xzzz, g_yy_0_xyyy_xzzzz, g_yy_0_xyyy_yyyy, g_yy_0_xyyy_yyyz, g_yy_0_xyyy_yyzz, g_yy_0_xyyy_yzzz, g_yy_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyy_xxxx[k] = -g_yy_0_xyyy_xxxx[k] * ab_x + g_yy_0_xyyy_xxxxx[k];

                g_yy_0_xxyyy_xxxy[k] = -g_yy_0_xyyy_xxxy[k] * ab_x + g_yy_0_xyyy_xxxxy[k];

                g_yy_0_xxyyy_xxxz[k] = -g_yy_0_xyyy_xxxz[k] * ab_x + g_yy_0_xyyy_xxxxz[k];

                g_yy_0_xxyyy_xxyy[k] = -g_yy_0_xyyy_xxyy[k] * ab_x + g_yy_0_xyyy_xxxyy[k];

                g_yy_0_xxyyy_xxyz[k] = -g_yy_0_xyyy_xxyz[k] * ab_x + g_yy_0_xyyy_xxxyz[k];

                g_yy_0_xxyyy_xxzz[k] = -g_yy_0_xyyy_xxzz[k] * ab_x + g_yy_0_xyyy_xxxzz[k];

                g_yy_0_xxyyy_xyyy[k] = -g_yy_0_xyyy_xyyy[k] * ab_x + g_yy_0_xyyy_xxyyy[k];

                g_yy_0_xxyyy_xyyz[k] = -g_yy_0_xyyy_xyyz[k] * ab_x + g_yy_0_xyyy_xxyyz[k];

                g_yy_0_xxyyy_xyzz[k] = -g_yy_0_xyyy_xyzz[k] * ab_x + g_yy_0_xyyy_xxyzz[k];

                g_yy_0_xxyyy_xzzz[k] = -g_yy_0_xyyy_xzzz[k] * ab_x + g_yy_0_xyyy_xxzzz[k];

                g_yy_0_xxyyy_yyyy[k] = -g_yy_0_xyyy_yyyy[k] * ab_x + g_yy_0_xyyy_xyyyy[k];

                g_yy_0_xxyyy_yyyz[k] = -g_yy_0_xyyy_yyyz[k] * ab_x + g_yy_0_xyyy_xyyyz[k];

                g_yy_0_xxyyy_yyzz[k] = -g_yy_0_xyyy_yyzz[k] * ab_x + g_yy_0_xyyy_xyyzz[k];

                g_yy_0_xxyyy_yzzz[k] = -g_yy_0_xyyy_yzzz[k] * ab_x + g_yy_0_xyyy_xyzzz[k];

                g_yy_0_xxyyy_zzzz[k] = -g_yy_0_xyyy_zzzz[k] * ab_x + g_yy_0_xyyy_xzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 1050 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 1051 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 1052 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 1053 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 1054 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 1055 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 1056 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 1057 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 1058 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 1059 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 1060 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 1061 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 1062 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 1063 * ccomps * dcomps);

            auto g_yy_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyz_xxxx, g_yy_0_xxyyz_xxxy, g_yy_0_xxyyz_xxxz, g_yy_0_xxyyz_xxyy, g_yy_0_xxyyz_xxyz, g_yy_0_xxyyz_xxzz, g_yy_0_xxyyz_xyyy, g_yy_0_xxyyz_xyyz, g_yy_0_xxyyz_xyzz, g_yy_0_xxyyz_xzzz, g_yy_0_xxyyz_yyyy, g_yy_0_xxyyz_yyyz, g_yy_0_xxyyz_yyzz, g_yy_0_xxyyz_yzzz, g_yy_0_xxyyz_zzzz, g_yy_0_xyyz_xxxx, g_yy_0_xyyz_xxxxx, g_yy_0_xyyz_xxxxy, g_yy_0_xyyz_xxxxz, g_yy_0_xyyz_xxxy, g_yy_0_xyyz_xxxyy, g_yy_0_xyyz_xxxyz, g_yy_0_xyyz_xxxz, g_yy_0_xyyz_xxxzz, g_yy_0_xyyz_xxyy, g_yy_0_xyyz_xxyyy, g_yy_0_xyyz_xxyyz, g_yy_0_xyyz_xxyz, g_yy_0_xyyz_xxyzz, g_yy_0_xyyz_xxzz, g_yy_0_xyyz_xxzzz, g_yy_0_xyyz_xyyy, g_yy_0_xyyz_xyyyy, g_yy_0_xyyz_xyyyz, g_yy_0_xyyz_xyyz, g_yy_0_xyyz_xyyzz, g_yy_0_xyyz_xyzz, g_yy_0_xyyz_xyzzz, g_yy_0_xyyz_xzzz, g_yy_0_xyyz_xzzzz, g_yy_0_xyyz_yyyy, g_yy_0_xyyz_yyyz, g_yy_0_xyyz_yyzz, g_yy_0_xyyz_yzzz, g_yy_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyz_xxxx[k] = -g_yy_0_xyyz_xxxx[k] * ab_x + g_yy_0_xyyz_xxxxx[k];

                g_yy_0_xxyyz_xxxy[k] = -g_yy_0_xyyz_xxxy[k] * ab_x + g_yy_0_xyyz_xxxxy[k];

                g_yy_0_xxyyz_xxxz[k] = -g_yy_0_xyyz_xxxz[k] * ab_x + g_yy_0_xyyz_xxxxz[k];

                g_yy_0_xxyyz_xxyy[k] = -g_yy_0_xyyz_xxyy[k] * ab_x + g_yy_0_xyyz_xxxyy[k];

                g_yy_0_xxyyz_xxyz[k] = -g_yy_0_xyyz_xxyz[k] * ab_x + g_yy_0_xyyz_xxxyz[k];

                g_yy_0_xxyyz_xxzz[k] = -g_yy_0_xyyz_xxzz[k] * ab_x + g_yy_0_xyyz_xxxzz[k];

                g_yy_0_xxyyz_xyyy[k] = -g_yy_0_xyyz_xyyy[k] * ab_x + g_yy_0_xyyz_xxyyy[k];

                g_yy_0_xxyyz_xyyz[k] = -g_yy_0_xyyz_xyyz[k] * ab_x + g_yy_0_xyyz_xxyyz[k];

                g_yy_0_xxyyz_xyzz[k] = -g_yy_0_xyyz_xyzz[k] * ab_x + g_yy_0_xyyz_xxyzz[k];

                g_yy_0_xxyyz_xzzz[k] = -g_yy_0_xyyz_xzzz[k] * ab_x + g_yy_0_xyyz_xxzzz[k];

                g_yy_0_xxyyz_yyyy[k] = -g_yy_0_xyyz_yyyy[k] * ab_x + g_yy_0_xyyz_xyyyy[k];

                g_yy_0_xxyyz_yyyz[k] = -g_yy_0_xyyz_yyyz[k] * ab_x + g_yy_0_xyyz_xyyyz[k];

                g_yy_0_xxyyz_yyzz[k] = -g_yy_0_xyyz_yyzz[k] * ab_x + g_yy_0_xyyz_xyyzz[k];

                g_yy_0_xxyyz_yzzz[k] = -g_yy_0_xyyz_yzzz[k] * ab_x + g_yy_0_xyyz_xyzzz[k];

                g_yy_0_xxyyz_zzzz[k] = -g_yy_0_xyyz_zzzz[k] * ab_x + g_yy_0_xyyz_xzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 1065 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 1066 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 1067 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 1068 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 1069 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 1070 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 1071 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 1072 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 1073 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 1074 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 1075 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 1076 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 1077 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 1078 * ccomps * dcomps);

            auto g_yy_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzz_xxxx, g_yy_0_xxyzz_xxxy, g_yy_0_xxyzz_xxxz, g_yy_0_xxyzz_xxyy, g_yy_0_xxyzz_xxyz, g_yy_0_xxyzz_xxzz, g_yy_0_xxyzz_xyyy, g_yy_0_xxyzz_xyyz, g_yy_0_xxyzz_xyzz, g_yy_0_xxyzz_xzzz, g_yy_0_xxyzz_yyyy, g_yy_0_xxyzz_yyyz, g_yy_0_xxyzz_yyzz, g_yy_0_xxyzz_yzzz, g_yy_0_xxyzz_zzzz, g_yy_0_xyzz_xxxx, g_yy_0_xyzz_xxxxx, g_yy_0_xyzz_xxxxy, g_yy_0_xyzz_xxxxz, g_yy_0_xyzz_xxxy, g_yy_0_xyzz_xxxyy, g_yy_0_xyzz_xxxyz, g_yy_0_xyzz_xxxz, g_yy_0_xyzz_xxxzz, g_yy_0_xyzz_xxyy, g_yy_0_xyzz_xxyyy, g_yy_0_xyzz_xxyyz, g_yy_0_xyzz_xxyz, g_yy_0_xyzz_xxyzz, g_yy_0_xyzz_xxzz, g_yy_0_xyzz_xxzzz, g_yy_0_xyzz_xyyy, g_yy_0_xyzz_xyyyy, g_yy_0_xyzz_xyyyz, g_yy_0_xyzz_xyyz, g_yy_0_xyzz_xyyzz, g_yy_0_xyzz_xyzz, g_yy_0_xyzz_xyzzz, g_yy_0_xyzz_xzzz, g_yy_0_xyzz_xzzzz, g_yy_0_xyzz_yyyy, g_yy_0_xyzz_yyyz, g_yy_0_xyzz_yyzz, g_yy_0_xyzz_yzzz, g_yy_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzz_xxxx[k] = -g_yy_0_xyzz_xxxx[k] * ab_x + g_yy_0_xyzz_xxxxx[k];

                g_yy_0_xxyzz_xxxy[k] = -g_yy_0_xyzz_xxxy[k] * ab_x + g_yy_0_xyzz_xxxxy[k];

                g_yy_0_xxyzz_xxxz[k] = -g_yy_0_xyzz_xxxz[k] * ab_x + g_yy_0_xyzz_xxxxz[k];

                g_yy_0_xxyzz_xxyy[k] = -g_yy_0_xyzz_xxyy[k] * ab_x + g_yy_0_xyzz_xxxyy[k];

                g_yy_0_xxyzz_xxyz[k] = -g_yy_0_xyzz_xxyz[k] * ab_x + g_yy_0_xyzz_xxxyz[k];

                g_yy_0_xxyzz_xxzz[k] = -g_yy_0_xyzz_xxzz[k] * ab_x + g_yy_0_xyzz_xxxzz[k];

                g_yy_0_xxyzz_xyyy[k] = -g_yy_0_xyzz_xyyy[k] * ab_x + g_yy_0_xyzz_xxyyy[k];

                g_yy_0_xxyzz_xyyz[k] = -g_yy_0_xyzz_xyyz[k] * ab_x + g_yy_0_xyzz_xxyyz[k];

                g_yy_0_xxyzz_xyzz[k] = -g_yy_0_xyzz_xyzz[k] * ab_x + g_yy_0_xyzz_xxyzz[k];

                g_yy_0_xxyzz_xzzz[k] = -g_yy_0_xyzz_xzzz[k] * ab_x + g_yy_0_xyzz_xxzzz[k];

                g_yy_0_xxyzz_yyyy[k] = -g_yy_0_xyzz_yyyy[k] * ab_x + g_yy_0_xyzz_xyyyy[k];

                g_yy_0_xxyzz_yyyz[k] = -g_yy_0_xyzz_yyyz[k] * ab_x + g_yy_0_xyzz_xyyyz[k];

                g_yy_0_xxyzz_yyzz[k] = -g_yy_0_xyzz_yyzz[k] * ab_x + g_yy_0_xyzz_xyyzz[k];

                g_yy_0_xxyzz_yzzz[k] = -g_yy_0_xyzz_yzzz[k] * ab_x + g_yy_0_xyzz_xyzzz[k];

                g_yy_0_xxyzz_zzzz[k] = -g_yy_0_xyzz_zzzz[k] * ab_x + g_yy_0_xyzz_xzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 1080 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 1081 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 1082 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 1083 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 1084 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 1085 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 1086 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 1087 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 1088 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 1089 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 1090 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 1091 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 1092 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 1093 * ccomps * dcomps);

            auto g_yy_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzz_xxxx, g_yy_0_xxzzz_xxxy, g_yy_0_xxzzz_xxxz, g_yy_0_xxzzz_xxyy, g_yy_0_xxzzz_xxyz, g_yy_0_xxzzz_xxzz, g_yy_0_xxzzz_xyyy, g_yy_0_xxzzz_xyyz, g_yy_0_xxzzz_xyzz, g_yy_0_xxzzz_xzzz, g_yy_0_xxzzz_yyyy, g_yy_0_xxzzz_yyyz, g_yy_0_xxzzz_yyzz, g_yy_0_xxzzz_yzzz, g_yy_0_xxzzz_zzzz, g_yy_0_xzzz_xxxx, g_yy_0_xzzz_xxxxx, g_yy_0_xzzz_xxxxy, g_yy_0_xzzz_xxxxz, g_yy_0_xzzz_xxxy, g_yy_0_xzzz_xxxyy, g_yy_0_xzzz_xxxyz, g_yy_0_xzzz_xxxz, g_yy_0_xzzz_xxxzz, g_yy_0_xzzz_xxyy, g_yy_0_xzzz_xxyyy, g_yy_0_xzzz_xxyyz, g_yy_0_xzzz_xxyz, g_yy_0_xzzz_xxyzz, g_yy_0_xzzz_xxzz, g_yy_0_xzzz_xxzzz, g_yy_0_xzzz_xyyy, g_yy_0_xzzz_xyyyy, g_yy_0_xzzz_xyyyz, g_yy_0_xzzz_xyyz, g_yy_0_xzzz_xyyzz, g_yy_0_xzzz_xyzz, g_yy_0_xzzz_xyzzz, g_yy_0_xzzz_xzzz, g_yy_0_xzzz_xzzzz, g_yy_0_xzzz_yyyy, g_yy_0_xzzz_yyyz, g_yy_0_xzzz_yyzz, g_yy_0_xzzz_yzzz, g_yy_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzz_xxxx[k] = -g_yy_0_xzzz_xxxx[k] * ab_x + g_yy_0_xzzz_xxxxx[k];

                g_yy_0_xxzzz_xxxy[k] = -g_yy_0_xzzz_xxxy[k] * ab_x + g_yy_0_xzzz_xxxxy[k];

                g_yy_0_xxzzz_xxxz[k] = -g_yy_0_xzzz_xxxz[k] * ab_x + g_yy_0_xzzz_xxxxz[k];

                g_yy_0_xxzzz_xxyy[k] = -g_yy_0_xzzz_xxyy[k] * ab_x + g_yy_0_xzzz_xxxyy[k];

                g_yy_0_xxzzz_xxyz[k] = -g_yy_0_xzzz_xxyz[k] * ab_x + g_yy_0_xzzz_xxxyz[k];

                g_yy_0_xxzzz_xxzz[k] = -g_yy_0_xzzz_xxzz[k] * ab_x + g_yy_0_xzzz_xxxzz[k];

                g_yy_0_xxzzz_xyyy[k] = -g_yy_0_xzzz_xyyy[k] * ab_x + g_yy_0_xzzz_xxyyy[k];

                g_yy_0_xxzzz_xyyz[k] = -g_yy_0_xzzz_xyyz[k] * ab_x + g_yy_0_xzzz_xxyyz[k];

                g_yy_0_xxzzz_xyzz[k] = -g_yy_0_xzzz_xyzz[k] * ab_x + g_yy_0_xzzz_xxyzz[k];

                g_yy_0_xxzzz_xzzz[k] = -g_yy_0_xzzz_xzzz[k] * ab_x + g_yy_0_xzzz_xxzzz[k];

                g_yy_0_xxzzz_yyyy[k] = -g_yy_0_xzzz_yyyy[k] * ab_x + g_yy_0_xzzz_xyyyy[k];

                g_yy_0_xxzzz_yyyz[k] = -g_yy_0_xzzz_yyyz[k] * ab_x + g_yy_0_xzzz_xyyyz[k];

                g_yy_0_xxzzz_yyzz[k] = -g_yy_0_xzzz_yyzz[k] * ab_x + g_yy_0_xzzz_xyyzz[k];

                g_yy_0_xxzzz_yzzz[k] = -g_yy_0_xzzz_yzzz[k] * ab_x + g_yy_0_xzzz_xyzzz[k];

                g_yy_0_xxzzz_zzzz[k] = -g_yy_0_xzzz_zzzz[k] * ab_x + g_yy_0_xzzz_xzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1095 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1096 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1097 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1098 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1099 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1100 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1101 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1102 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1103 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1104 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1105 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1106 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1107 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1108 * ccomps * dcomps);

            auto g_yy_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyy_xxxx, g_yy_0_xyyyy_xxxy, g_yy_0_xyyyy_xxxz, g_yy_0_xyyyy_xxyy, g_yy_0_xyyyy_xxyz, g_yy_0_xyyyy_xxzz, g_yy_0_xyyyy_xyyy, g_yy_0_xyyyy_xyyz, g_yy_0_xyyyy_xyzz, g_yy_0_xyyyy_xzzz, g_yy_0_xyyyy_yyyy, g_yy_0_xyyyy_yyyz, g_yy_0_xyyyy_yyzz, g_yy_0_xyyyy_yzzz, g_yy_0_xyyyy_zzzz, g_yy_0_yyyy_xxxx, g_yy_0_yyyy_xxxxx, g_yy_0_yyyy_xxxxy, g_yy_0_yyyy_xxxxz, g_yy_0_yyyy_xxxy, g_yy_0_yyyy_xxxyy, g_yy_0_yyyy_xxxyz, g_yy_0_yyyy_xxxz, g_yy_0_yyyy_xxxzz, g_yy_0_yyyy_xxyy, g_yy_0_yyyy_xxyyy, g_yy_0_yyyy_xxyyz, g_yy_0_yyyy_xxyz, g_yy_0_yyyy_xxyzz, g_yy_0_yyyy_xxzz, g_yy_0_yyyy_xxzzz, g_yy_0_yyyy_xyyy, g_yy_0_yyyy_xyyyy, g_yy_0_yyyy_xyyyz, g_yy_0_yyyy_xyyz, g_yy_0_yyyy_xyyzz, g_yy_0_yyyy_xyzz, g_yy_0_yyyy_xyzzz, g_yy_0_yyyy_xzzz, g_yy_0_yyyy_xzzzz, g_yy_0_yyyy_yyyy, g_yy_0_yyyy_yyyz, g_yy_0_yyyy_yyzz, g_yy_0_yyyy_yzzz, g_yy_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyy_xxxx[k] = -g_yy_0_yyyy_xxxx[k] * ab_x + g_yy_0_yyyy_xxxxx[k];

                g_yy_0_xyyyy_xxxy[k] = -g_yy_0_yyyy_xxxy[k] * ab_x + g_yy_0_yyyy_xxxxy[k];

                g_yy_0_xyyyy_xxxz[k] = -g_yy_0_yyyy_xxxz[k] * ab_x + g_yy_0_yyyy_xxxxz[k];

                g_yy_0_xyyyy_xxyy[k] = -g_yy_0_yyyy_xxyy[k] * ab_x + g_yy_0_yyyy_xxxyy[k];

                g_yy_0_xyyyy_xxyz[k] = -g_yy_0_yyyy_xxyz[k] * ab_x + g_yy_0_yyyy_xxxyz[k];

                g_yy_0_xyyyy_xxzz[k] = -g_yy_0_yyyy_xxzz[k] * ab_x + g_yy_0_yyyy_xxxzz[k];

                g_yy_0_xyyyy_xyyy[k] = -g_yy_0_yyyy_xyyy[k] * ab_x + g_yy_0_yyyy_xxyyy[k];

                g_yy_0_xyyyy_xyyz[k] = -g_yy_0_yyyy_xyyz[k] * ab_x + g_yy_0_yyyy_xxyyz[k];

                g_yy_0_xyyyy_xyzz[k] = -g_yy_0_yyyy_xyzz[k] * ab_x + g_yy_0_yyyy_xxyzz[k];

                g_yy_0_xyyyy_xzzz[k] = -g_yy_0_yyyy_xzzz[k] * ab_x + g_yy_0_yyyy_xxzzz[k];

                g_yy_0_xyyyy_yyyy[k] = -g_yy_0_yyyy_yyyy[k] * ab_x + g_yy_0_yyyy_xyyyy[k];

                g_yy_0_xyyyy_yyyz[k] = -g_yy_0_yyyy_yyyz[k] * ab_x + g_yy_0_yyyy_xyyyz[k];

                g_yy_0_xyyyy_yyzz[k] = -g_yy_0_yyyy_yyzz[k] * ab_x + g_yy_0_yyyy_xyyzz[k];

                g_yy_0_xyyyy_yzzz[k] = -g_yy_0_yyyy_yzzz[k] * ab_x + g_yy_0_yyyy_xyzzz[k];

                g_yy_0_xyyyy_zzzz[k] = -g_yy_0_yyyy_zzzz[k] * ab_x + g_yy_0_yyyy_xzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1110 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1111 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1112 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1113 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1114 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1115 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1116 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1117 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1118 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1119 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1120 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1121 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1122 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1123 * ccomps * dcomps);

            auto g_yy_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyz_xxxx, g_yy_0_xyyyz_xxxy, g_yy_0_xyyyz_xxxz, g_yy_0_xyyyz_xxyy, g_yy_0_xyyyz_xxyz, g_yy_0_xyyyz_xxzz, g_yy_0_xyyyz_xyyy, g_yy_0_xyyyz_xyyz, g_yy_0_xyyyz_xyzz, g_yy_0_xyyyz_xzzz, g_yy_0_xyyyz_yyyy, g_yy_0_xyyyz_yyyz, g_yy_0_xyyyz_yyzz, g_yy_0_xyyyz_yzzz, g_yy_0_xyyyz_zzzz, g_yy_0_yyyz_xxxx, g_yy_0_yyyz_xxxxx, g_yy_0_yyyz_xxxxy, g_yy_0_yyyz_xxxxz, g_yy_0_yyyz_xxxy, g_yy_0_yyyz_xxxyy, g_yy_0_yyyz_xxxyz, g_yy_0_yyyz_xxxz, g_yy_0_yyyz_xxxzz, g_yy_0_yyyz_xxyy, g_yy_0_yyyz_xxyyy, g_yy_0_yyyz_xxyyz, g_yy_0_yyyz_xxyz, g_yy_0_yyyz_xxyzz, g_yy_0_yyyz_xxzz, g_yy_0_yyyz_xxzzz, g_yy_0_yyyz_xyyy, g_yy_0_yyyz_xyyyy, g_yy_0_yyyz_xyyyz, g_yy_0_yyyz_xyyz, g_yy_0_yyyz_xyyzz, g_yy_0_yyyz_xyzz, g_yy_0_yyyz_xyzzz, g_yy_0_yyyz_xzzz, g_yy_0_yyyz_xzzzz, g_yy_0_yyyz_yyyy, g_yy_0_yyyz_yyyz, g_yy_0_yyyz_yyzz, g_yy_0_yyyz_yzzz, g_yy_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyz_xxxx[k] = -g_yy_0_yyyz_xxxx[k] * ab_x + g_yy_0_yyyz_xxxxx[k];

                g_yy_0_xyyyz_xxxy[k] = -g_yy_0_yyyz_xxxy[k] * ab_x + g_yy_0_yyyz_xxxxy[k];

                g_yy_0_xyyyz_xxxz[k] = -g_yy_0_yyyz_xxxz[k] * ab_x + g_yy_0_yyyz_xxxxz[k];

                g_yy_0_xyyyz_xxyy[k] = -g_yy_0_yyyz_xxyy[k] * ab_x + g_yy_0_yyyz_xxxyy[k];

                g_yy_0_xyyyz_xxyz[k] = -g_yy_0_yyyz_xxyz[k] * ab_x + g_yy_0_yyyz_xxxyz[k];

                g_yy_0_xyyyz_xxzz[k] = -g_yy_0_yyyz_xxzz[k] * ab_x + g_yy_0_yyyz_xxxzz[k];

                g_yy_0_xyyyz_xyyy[k] = -g_yy_0_yyyz_xyyy[k] * ab_x + g_yy_0_yyyz_xxyyy[k];

                g_yy_0_xyyyz_xyyz[k] = -g_yy_0_yyyz_xyyz[k] * ab_x + g_yy_0_yyyz_xxyyz[k];

                g_yy_0_xyyyz_xyzz[k] = -g_yy_0_yyyz_xyzz[k] * ab_x + g_yy_0_yyyz_xxyzz[k];

                g_yy_0_xyyyz_xzzz[k] = -g_yy_0_yyyz_xzzz[k] * ab_x + g_yy_0_yyyz_xxzzz[k];

                g_yy_0_xyyyz_yyyy[k] = -g_yy_0_yyyz_yyyy[k] * ab_x + g_yy_0_yyyz_xyyyy[k];

                g_yy_0_xyyyz_yyyz[k] = -g_yy_0_yyyz_yyyz[k] * ab_x + g_yy_0_yyyz_xyyyz[k];

                g_yy_0_xyyyz_yyzz[k] = -g_yy_0_yyyz_yyzz[k] * ab_x + g_yy_0_yyyz_xyyzz[k];

                g_yy_0_xyyyz_yzzz[k] = -g_yy_0_yyyz_yzzz[k] * ab_x + g_yy_0_yyyz_xyzzz[k];

                g_yy_0_xyyyz_zzzz[k] = -g_yy_0_yyyz_zzzz[k] * ab_x + g_yy_0_yyyz_xzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1125 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1126 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1127 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1128 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1129 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1130 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1131 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1132 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1133 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1134 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1135 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1136 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1137 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1138 * ccomps * dcomps);

            auto g_yy_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzz_xxxx, g_yy_0_xyyzz_xxxy, g_yy_0_xyyzz_xxxz, g_yy_0_xyyzz_xxyy, g_yy_0_xyyzz_xxyz, g_yy_0_xyyzz_xxzz, g_yy_0_xyyzz_xyyy, g_yy_0_xyyzz_xyyz, g_yy_0_xyyzz_xyzz, g_yy_0_xyyzz_xzzz, g_yy_0_xyyzz_yyyy, g_yy_0_xyyzz_yyyz, g_yy_0_xyyzz_yyzz, g_yy_0_xyyzz_yzzz, g_yy_0_xyyzz_zzzz, g_yy_0_yyzz_xxxx, g_yy_0_yyzz_xxxxx, g_yy_0_yyzz_xxxxy, g_yy_0_yyzz_xxxxz, g_yy_0_yyzz_xxxy, g_yy_0_yyzz_xxxyy, g_yy_0_yyzz_xxxyz, g_yy_0_yyzz_xxxz, g_yy_0_yyzz_xxxzz, g_yy_0_yyzz_xxyy, g_yy_0_yyzz_xxyyy, g_yy_0_yyzz_xxyyz, g_yy_0_yyzz_xxyz, g_yy_0_yyzz_xxyzz, g_yy_0_yyzz_xxzz, g_yy_0_yyzz_xxzzz, g_yy_0_yyzz_xyyy, g_yy_0_yyzz_xyyyy, g_yy_0_yyzz_xyyyz, g_yy_0_yyzz_xyyz, g_yy_0_yyzz_xyyzz, g_yy_0_yyzz_xyzz, g_yy_0_yyzz_xyzzz, g_yy_0_yyzz_xzzz, g_yy_0_yyzz_xzzzz, g_yy_0_yyzz_yyyy, g_yy_0_yyzz_yyyz, g_yy_0_yyzz_yyzz, g_yy_0_yyzz_yzzz, g_yy_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzz_xxxx[k] = -g_yy_0_yyzz_xxxx[k] * ab_x + g_yy_0_yyzz_xxxxx[k];

                g_yy_0_xyyzz_xxxy[k] = -g_yy_0_yyzz_xxxy[k] * ab_x + g_yy_0_yyzz_xxxxy[k];

                g_yy_0_xyyzz_xxxz[k] = -g_yy_0_yyzz_xxxz[k] * ab_x + g_yy_0_yyzz_xxxxz[k];

                g_yy_0_xyyzz_xxyy[k] = -g_yy_0_yyzz_xxyy[k] * ab_x + g_yy_0_yyzz_xxxyy[k];

                g_yy_0_xyyzz_xxyz[k] = -g_yy_0_yyzz_xxyz[k] * ab_x + g_yy_0_yyzz_xxxyz[k];

                g_yy_0_xyyzz_xxzz[k] = -g_yy_0_yyzz_xxzz[k] * ab_x + g_yy_0_yyzz_xxxzz[k];

                g_yy_0_xyyzz_xyyy[k] = -g_yy_0_yyzz_xyyy[k] * ab_x + g_yy_0_yyzz_xxyyy[k];

                g_yy_0_xyyzz_xyyz[k] = -g_yy_0_yyzz_xyyz[k] * ab_x + g_yy_0_yyzz_xxyyz[k];

                g_yy_0_xyyzz_xyzz[k] = -g_yy_0_yyzz_xyzz[k] * ab_x + g_yy_0_yyzz_xxyzz[k];

                g_yy_0_xyyzz_xzzz[k] = -g_yy_0_yyzz_xzzz[k] * ab_x + g_yy_0_yyzz_xxzzz[k];

                g_yy_0_xyyzz_yyyy[k] = -g_yy_0_yyzz_yyyy[k] * ab_x + g_yy_0_yyzz_xyyyy[k];

                g_yy_0_xyyzz_yyyz[k] = -g_yy_0_yyzz_yyyz[k] * ab_x + g_yy_0_yyzz_xyyyz[k];

                g_yy_0_xyyzz_yyzz[k] = -g_yy_0_yyzz_yyzz[k] * ab_x + g_yy_0_yyzz_xyyzz[k];

                g_yy_0_xyyzz_yzzz[k] = -g_yy_0_yyzz_yzzz[k] * ab_x + g_yy_0_yyzz_xyzzz[k];

                g_yy_0_xyyzz_zzzz[k] = -g_yy_0_yyzz_zzzz[k] * ab_x + g_yy_0_yyzz_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1140 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1141 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1142 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1143 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1144 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1145 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1146 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1147 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1148 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1149 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1150 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1151 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1152 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1153 * ccomps * dcomps);

            auto g_yy_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzz_xxxx, g_yy_0_xyzzz_xxxy, g_yy_0_xyzzz_xxxz, g_yy_0_xyzzz_xxyy, g_yy_0_xyzzz_xxyz, g_yy_0_xyzzz_xxzz, g_yy_0_xyzzz_xyyy, g_yy_0_xyzzz_xyyz, g_yy_0_xyzzz_xyzz, g_yy_0_xyzzz_xzzz, g_yy_0_xyzzz_yyyy, g_yy_0_xyzzz_yyyz, g_yy_0_xyzzz_yyzz, g_yy_0_xyzzz_yzzz, g_yy_0_xyzzz_zzzz, g_yy_0_yzzz_xxxx, g_yy_0_yzzz_xxxxx, g_yy_0_yzzz_xxxxy, g_yy_0_yzzz_xxxxz, g_yy_0_yzzz_xxxy, g_yy_0_yzzz_xxxyy, g_yy_0_yzzz_xxxyz, g_yy_0_yzzz_xxxz, g_yy_0_yzzz_xxxzz, g_yy_0_yzzz_xxyy, g_yy_0_yzzz_xxyyy, g_yy_0_yzzz_xxyyz, g_yy_0_yzzz_xxyz, g_yy_0_yzzz_xxyzz, g_yy_0_yzzz_xxzz, g_yy_0_yzzz_xxzzz, g_yy_0_yzzz_xyyy, g_yy_0_yzzz_xyyyy, g_yy_0_yzzz_xyyyz, g_yy_0_yzzz_xyyz, g_yy_0_yzzz_xyyzz, g_yy_0_yzzz_xyzz, g_yy_0_yzzz_xyzzz, g_yy_0_yzzz_xzzz, g_yy_0_yzzz_xzzzz, g_yy_0_yzzz_yyyy, g_yy_0_yzzz_yyyz, g_yy_0_yzzz_yyzz, g_yy_0_yzzz_yzzz, g_yy_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzz_xxxx[k] = -g_yy_0_yzzz_xxxx[k] * ab_x + g_yy_0_yzzz_xxxxx[k];

                g_yy_0_xyzzz_xxxy[k] = -g_yy_0_yzzz_xxxy[k] * ab_x + g_yy_0_yzzz_xxxxy[k];

                g_yy_0_xyzzz_xxxz[k] = -g_yy_0_yzzz_xxxz[k] * ab_x + g_yy_0_yzzz_xxxxz[k];

                g_yy_0_xyzzz_xxyy[k] = -g_yy_0_yzzz_xxyy[k] * ab_x + g_yy_0_yzzz_xxxyy[k];

                g_yy_0_xyzzz_xxyz[k] = -g_yy_0_yzzz_xxyz[k] * ab_x + g_yy_0_yzzz_xxxyz[k];

                g_yy_0_xyzzz_xxzz[k] = -g_yy_0_yzzz_xxzz[k] * ab_x + g_yy_0_yzzz_xxxzz[k];

                g_yy_0_xyzzz_xyyy[k] = -g_yy_0_yzzz_xyyy[k] * ab_x + g_yy_0_yzzz_xxyyy[k];

                g_yy_0_xyzzz_xyyz[k] = -g_yy_0_yzzz_xyyz[k] * ab_x + g_yy_0_yzzz_xxyyz[k];

                g_yy_0_xyzzz_xyzz[k] = -g_yy_0_yzzz_xyzz[k] * ab_x + g_yy_0_yzzz_xxyzz[k];

                g_yy_0_xyzzz_xzzz[k] = -g_yy_0_yzzz_xzzz[k] * ab_x + g_yy_0_yzzz_xxzzz[k];

                g_yy_0_xyzzz_yyyy[k] = -g_yy_0_yzzz_yyyy[k] * ab_x + g_yy_0_yzzz_xyyyy[k];

                g_yy_0_xyzzz_yyyz[k] = -g_yy_0_yzzz_yyyz[k] * ab_x + g_yy_0_yzzz_xyyyz[k];

                g_yy_0_xyzzz_yyzz[k] = -g_yy_0_yzzz_yyzz[k] * ab_x + g_yy_0_yzzz_xyyzz[k];

                g_yy_0_xyzzz_yzzz[k] = -g_yy_0_yzzz_yzzz[k] * ab_x + g_yy_0_yzzz_xyzzz[k];

                g_yy_0_xyzzz_zzzz[k] = -g_yy_0_yzzz_zzzz[k] * ab_x + g_yy_0_yzzz_xzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1155 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1156 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1157 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1158 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1159 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1160 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1161 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1162 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1163 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1164 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1165 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1166 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1167 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1168 * ccomps * dcomps);

            auto g_yy_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzz_xxxx, g_yy_0_xzzzz_xxxy, g_yy_0_xzzzz_xxxz, g_yy_0_xzzzz_xxyy, g_yy_0_xzzzz_xxyz, g_yy_0_xzzzz_xxzz, g_yy_0_xzzzz_xyyy, g_yy_0_xzzzz_xyyz, g_yy_0_xzzzz_xyzz, g_yy_0_xzzzz_xzzz, g_yy_0_xzzzz_yyyy, g_yy_0_xzzzz_yyyz, g_yy_0_xzzzz_yyzz, g_yy_0_xzzzz_yzzz, g_yy_0_xzzzz_zzzz, g_yy_0_zzzz_xxxx, g_yy_0_zzzz_xxxxx, g_yy_0_zzzz_xxxxy, g_yy_0_zzzz_xxxxz, g_yy_0_zzzz_xxxy, g_yy_0_zzzz_xxxyy, g_yy_0_zzzz_xxxyz, g_yy_0_zzzz_xxxz, g_yy_0_zzzz_xxxzz, g_yy_0_zzzz_xxyy, g_yy_0_zzzz_xxyyy, g_yy_0_zzzz_xxyyz, g_yy_0_zzzz_xxyz, g_yy_0_zzzz_xxyzz, g_yy_0_zzzz_xxzz, g_yy_0_zzzz_xxzzz, g_yy_0_zzzz_xyyy, g_yy_0_zzzz_xyyyy, g_yy_0_zzzz_xyyyz, g_yy_0_zzzz_xyyz, g_yy_0_zzzz_xyyzz, g_yy_0_zzzz_xyzz, g_yy_0_zzzz_xyzzz, g_yy_0_zzzz_xzzz, g_yy_0_zzzz_xzzzz, g_yy_0_zzzz_yyyy, g_yy_0_zzzz_yyyz, g_yy_0_zzzz_yyzz, g_yy_0_zzzz_yzzz, g_yy_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzz_xxxx[k] = -g_yy_0_zzzz_xxxx[k] * ab_x + g_yy_0_zzzz_xxxxx[k];

                g_yy_0_xzzzz_xxxy[k] = -g_yy_0_zzzz_xxxy[k] * ab_x + g_yy_0_zzzz_xxxxy[k];

                g_yy_0_xzzzz_xxxz[k] = -g_yy_0_zzzz_xxxz[k] * ab_x + g_yy_0_zzzz_xxxxz[k];

                g_yy_0_xzzzz_xxyy[k] = -g_yy_0_zzzz_xxyy[k] * ab_x + g_yy_0_zzzz_xxxyy[k];

                g_yy_0_xzzzz_xxyz[k] = -g_yy_0_zzzz_xxyz[k] * ab_x + g_yy_0_zzzz_xxxyz[k];

                g_yy_0_xzzzz_xxzz[k] = -g_yy_0_zzzz_xxzz[k] * ab_x + g_yy_0_zzzz_xxxzz[k];

                g_yy_0_xzzzz_xyyy[k] = -g_yy_0_zzzz_xyyy[k] * ab_x + g_yy_0_zzzz_xxyyy[k];

                g_yy_0_xzzzz_xyyz[k] = -g_yy_0_zzzz_xyyz[k] * ab_x + g_yy_0_zzzz_xxyyz[k];

                g_yy_0_xzzzz_xyzz[k] = -g_yy_0_zzzz_xyzz[k] * ab_x + g_yy_0_zzzz_xxyzz[k];

                g_yy_0_xzzzz_xzzz[k] = -g_yy_0_zzzz_xzzz[k] * ab_x + g_yy_0_zzzz_xxzzz[k];

                g_yy_0_xzzzz_yyyy[k] = -g_yy_0_zzzz_yyyy[k] * ab_x + g_yy_0_zzzz_xyyyy[k];

                g_yy_0_xzzzz_yyyz[k] = -g_yy_0_zzzz_yyyz[k] * ab_x + g_yy_0_zzzz_xyyyz[k];

                g_yy_0_xzzzz_yyzz[k] = -g_yy_0_zzzz_yyzz[k] * ab_x + g_yy_0_zzzz_xyyzz[k];

                g_yy_0_xzzzz_yzzz[k] = -g_yy_0_zzzz_yzzz[k] * ab_x + g_yy_0_zzzz_xyzzz[k];

                g_yy_0_xzzzz_zzzz[k] = -g_yy_0_zzzz_zzzz[k] * ab_x + g_yy_0_zzzz_xzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1170 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1171 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1172 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1173 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1174 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1175 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1176 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1177 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1178 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1179 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1180 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1181 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1182 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1183 * ccomps * dcomps);

            auto g_yy_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz, g_yy_0_yyyy_xxxx, g_yy_0_yyyy_xxxxy, g_yy_0_yyyy_xxxy, g_yy_0_yyyy_xxxyy, g_yy_0_yyyy_xxxyz, g_yy_0_yyyy_xxxz, g_yy_0_yyyy_xxyy, g_yy_0_yyyy_xxyyy, g_yy_0_yyyy_xxyyz, g_yy_0_yyyy_xxyz, g_yy_0_yyyy_xxyzz, g_yy_0_yyyy_xxzz, g_yy_0_yyyy_xyyy, g_yy_0_yyyy_xyyyy, g_yy_0_yyyy_xyyyz, g_yy_0_yyyy_xyyz, g_yy_0_yyyy_xyyzz, g_yy_0_yyyy_xyzz, g_yy_0_yyyy_xyzzz, g_yy_0_yyyy_xzzz, g_yy_0_yyyy_yyyy, g_yy_0_yyyy_yyyyy, g_yy_0_yyyy_yyyyz, g_yy_0_yyyy_yyyz, g_yy_0_yyyy_yyyzz, g_yy_0_yyyy_yyzz, g_yy_0_yyyy_yyzzz, g_yy_0_yyyy_yzzz, g_yy_0_yyyy_yzzzz, g_yy_0_yyyy_zzzz, g_yy_0_yyyyy_xxxx, g_yy_0_yyyyy_xxxy, g_yy_0_yyyyy_xxxz, g_yy_0_yyyyy_xxyy, g_yy_0_yyyyy_xxyz, g_yy_0_yyyyy_xxzz, g_yy_0_yyyyy_xyyy, g_yy_0_yyyyy_xyyz, g_yy_0_yyyyy_xyzz, g_yy_0_yyyyy_xzzz, g_yy_0_yyyyy_yyyy, g_yy_0_yyyyy_yyyz, g_yy_0_yyyyy_yyzz, g_yy_0_yyyyy_yzzz, g_yy_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyy_xxxx[k] = -2.0 * g_y_0_yyyy_xxxx[k] - g_yy_0_yyyy_xxxx[k] * ab_y + g_yy_0_yyyy_xxxxy[k];

                g_yy_0_yyyyy_xxxy[k] = -2.0 * g_y_0_yyyy_xxxy[k] - g_yy_0_yyyy_xxxy[k] * ab_y + g_yy_0_yyyy_xxxyy[k];

                g_yy_0_yyyyy_xxxz[k] = -2.0 * g_y_0_yyyy_xxxz[k] - g_yy_0_yyyy_xxxz[k] * ab_y + g_yy_0_yyyy_xxxyz[k];

                g_yy_0_yyyyy_xxyy[k] = -2.0 * g_y_0_yyyy_xxyy[k] - g_yy_0_yyyy_xxyy[k] * ab_y + g_yy_0_yyyy_xxyyy[k];

                g_yy_0_yyyyy_xxyz[k] = -2.0 * g_y_0_yyyy_xxyz[k] - g_yy_0_yyyy_xxyz[k] * ab_y + g_yy_0_yyyy_xxyyz[k];

                g_yy_0_yyyyy_xxzz[k] = -2.0 * g_y_0_yyyy_xxzz[k] - g_yy_0_yyyy_xxzz[k] * ab_y + g_yy_0_yyyy_xxyzz[k];

                g_yy_0_yyyyy_xyyy[k] = -2.0 * g_y_0_yyyy_xyyy[k] - g_yy_0_yyyy_xyyy[k] * ab_y + g_yy_0_yyyy_xyyyy[k];

                g_yy_0_yyyyy_xyyz[k] = -2.0 * g_y_0_yyyy_xyyz[k] - g_yy_0_yyyy_xyyz[k] * ab_y + g_yy_0_yyyy_xyyyz[k];

                g_yy_0_yyyyy_xyzz[k] = -2.0 * g_y_0_yyyy_xyzz[k] - g_yy_0_yyyy_xyzz[k] * ab_y + g_yy_0_yyyy_xyyzz[k];

                g_yy_0_yyyyy_xzzz[k] = -2.0 * g_y_0_yyyy_xzzz[k] - g_yy_0_yyyy_xzzz[k] * ab_y + g_yy_0_yyyy_xyzzz[k];

                g_yy_0_yyyyy_yyyy[k] = -2.0 * g_y_0_yyyy_yyyy[k] - g_yy_0_yyyy_yyyy[k] * ab_y + g_yy_0_yyyy_yyyyy[k];

                g_yy_0_yyyyy_yyyz[k] = -2.0 * g_y_0_yyyy_yyyz[k] - g_yy_0_yyyy_yyyz[k] * ab_y + g_yy_0_yyyy_yyyyz[k];

                g_yy_0_yyyyy_yyzz[k] = -2.0 * g_y_0_yyyy_yyzz[k] - g_yy_0_yyyy_yyzz[k] * ab_y + g_yy_0_yyyy_yyyzz[k];

                g_yy_0_yyyyy_yzzz[k] = -2.0 * g_y_0_yyyy_yzzz[k] - g_yy_0_yyyy_yzzz[k] * ab_y + g_yy_0_yyyy_yyzzz[k];

                g_yy_0_yyyyy_zzzz[k] = -2.0 * g_y_0_yyyy_zzzz[k] - g_yy_0_yyyy_zzzz[k] * ab_y + g_yy_0_yyyy_yzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1185 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1186 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1187 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1188 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1189 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1190 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1191 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1192 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1193 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1194 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1195 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1196 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1197 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1198 * ccomps * dcomps);

            auto g_yy_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyy_xxxx, g_yy_0_yyyy_xxxxz, g_yy_0_yyyy_xxxy, g_yy_0_yyyy_xxxyz, g_yy_0_yyyy_xxxz, g_yy_0_yyyy_xxxzz, g_yy_0_yyyy_xxyy, g_yy_0_yyyy_xxyyz, g_yy_0_yyyy_xxyz, g_yy_0_yyyy_xxyzz, g_yy_0_yyyy_xxzz, g_yy_0_yyyy_xxzzz, g_yy_0_yyyy_xyyy, g_yy_0_yyyy_xyyyz, g_yy_0_yyyy_xyyz, g_yy_0_yyyy_xyyzz, g_yy_0_yyyy_xyzz, g_yy_0_yyyy_xyzzz, g_yy_0_yyyy_xzzz, g_yy_0_yyyy_xzzzz, g_yy_0_yyyy_yyyy, g_yy_0_yyyy_yyyyz, g_yy_0_yyyy_yyyz, g_yy_0_yyyy_yyyzz, g_yy_0_yyyy_yyzz, g_yy_0_yyyy_yyzzz, g_yy_0_yyyy_yzzz, g_yy_0_yyyy_yzzzz, g_yy_0_yyyy_zzzz, g_yy_0_yyyy_zzzzz, g_yy_0_yyyyz_xxxx, g_yy_0_yyyyz_xxxy, g_yy_0_yyyyz_xxxz, g_yy_0_yyyyz_xxyy, g_yy_0_yyyyz_xxyz, g_yy_0_yyyyz_xxzz, g_yy_0_yyyyz_xyyy, g_yy_0_yyyyz_xyyz, g_yy_0_yyyyz_xyzz, g_yy_0_yyyyz_xzzz, g_yy_0_yyyyz_yyyy, g_yy_0_yyyyz_yyyz, g_yy_0_yyyyz_yyzz, g_yy_0_yyyyz_yzzz, g_yy_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyz_xxxx[k] = -g_yy_0_yyyy_xxxx[k] * ab_z + g_yy_0_yyyy_xxxxz[k];

                g_yy_0_yyyyz_xxxy[k] = -g_yy_0_yyyy_xxxy[k] * ab_z + g_yy_0_yyyy_xxxyz[k];

                g_yy_0_yyyyz_xxxz[k] = -g_yy_0_yyyy_xxxz[k] * ab_z + g_yy_0_yyyy_xxxzz[k];

                g_yy_0_yyyyz_xxyy[k] = -g_yy_0_yyyy_xxyy[k] * ab_z + g_yy_0_yyyy_xxyyz[k];

                g_yy_0_yyyyz_xxyz[k] = -g_yy_0_yyyy_xxyz[k] * ab_z + g_yy_0_yyyy_xxyzz[k];

                g_yy_0_yyyyz_xxzz[k] = -g_yy_0_yyyy_xxzz[k] * ab_z + g_yy_0_yyyy_xxzzz[k];

                g_yy_0_yyyyz_xyyy[k] = -g_yy_0_yyyy_xyyy[k] * ab_z + g_yy_0_yyyy_xyyyz[k];

                g_yy_0_yyyyz_xyyz[k] = -g_yy_0_yyyy_xyyz[k] * ab_z + g_yy_0_yyyy_xyyzz[k];

                g_yy_0_yyyyz_xyzz[k] = -g_yy_0_yyyy_xyzz[k] * ab_z + g_yy_0_yyyy_xyzzz[k];

                g_yy_0_yyyyz_xzzz[k] = -g_yy_0_yyyy_xzzz[k] * ab_z + g_yy_0_yyyy_xzzzz[k];

                g_yy_0_yyyyz_yyyy[k] = -g_yy_0_yyyy_yyyy[k] * ab_z + g_yy_0_yyyy_yyyyz[k];

                g_yy_0_yyyyz_yyyz[k] = -g_yy_0_yyyy_yyyz[k] * ab_z + g_yy_0_yyyy_yyyzz[k];

                g_yy_0_yyyyz_yyzz[k] = -g_yy_0_yyyy_yyzz[k] * ab_z + g_yy_0_yyyy_yyzzz[k];

                g_yy_0_yyyyz_yzzz[k] = -g_yy_0_yyyy_yzzz[k] * ab_z + g_yy_0_yyyy_yzzzz[k];

                g_yy_0_yyyyz_zzzz[k] = -g_yy_0_yyyy_zzzz[k] * ab_z + g_yy_0_yyyy_zzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1200 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1201 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1202 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1203 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1204 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1205 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1206 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1207 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1208 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1209 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1210 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1211 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1212 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1213 * ccomps * dcomps);

            auto g_yy_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyz_xxxx, g_yy_0_yyyz_xxxxz, g_yy_0_yyyz_xxxy, g_yy_0_yyyz_xxxyz, g_yy_0_yyyz_xxxz, g_yy_0_yyyz_xxxzz, g_yy_0_yyyz_xxyy, g_yy_0_yyyz_xxyyz, g_yy_0_yyyz_xxyz, g_yy_0_yyyz_xxyzz, g_yy_0_yyyz_xxzz, g_yy_0_yyyz_xxzzz, g_yy_0_yyyz_xyyy, g_yy_0_yyyz_xyyyz, g_yy_0_yyyz_xyyz, g_yy_0_yyyz_xyyzz, g_yy_0_yyyz_xyzz, g_yy_0_yyyz_xyzzz, g_yy_0_yyyz_xzzz, g_yy_0_yyyz_xzzzz, g_yy_0_yyyz_yyyy, g_yy_0_yyyz_yyyyz, g_yy_0_yyyz_yyyz, g_yy_0_yyyz_yyyzz, g_yy_0_yyyz_yyzz, g_yy_0_yyyz_yyzzz, g_yy_0_yyyz_yzzz, g_yy_0_yyyz_yzzzz, g_yy_0_yyyz_zzzz, g_yy_0_yyyz_zzzzz, g_yy_0_yyyzz_xxxx, g_yy_0_yyyzz_xxxy, g_yy_0_yyyzz_xxxz, g_yy_0_yyyzz_xxyy, g_yy_0_yyyzz_xxyz, g_yy_0_yyyzz_xxzz, g_yy_0_yyyzz_xyyy, g_yy_0_yyyzz_xyyz, g_yy_0_yyyzz_xyzz, g_yy_0_yyyzz_xzzz, g_yy_0_yyyzz_yyyy, g_yy_0_yyyzz_yyyz, g_yy_0_yyyzz_yyzz, g_yy_0_yyyzz_yzzz, g_yy_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzz_xxxx[k] = -g_yy_0_yyyz_xxxx[k] * ab_z + g_yy_0_yyyz_xxxxz[k];

                g_yy_0_yyyzz_xxxy[k] = -g_yy_0_yyyz_xxxy[k] * ab_z + g_yy_0_yyyz_xxxyz[k];

                g_yy_0_yyyzz_xxxz[k] = -g_yy_0_yyyz_xxxz[k] * ab_z + g_yy_0_yyyz_xxxzz[k];

                g_yy_0_yyyzz_xxyy[k] = -g_yy_0_yyyz_xxyy[k] * ab_z + g_yy_0_yyyz_xxyyz[k];

                g_yy_0_yyyzz_xxyz[k] = -g_yy_0_yyyz_xxyz[k] * ab_z + g_yy_0_yyyz_xxyzz[k];

                g_yy_0_yyyzz_xxzz[k] = -g_yy_0_yyyz_xxzz[k] * ab_z + g_yy_0_yyyz_xxzzz[k];

                g_yy_0_yyyzz_xyyy[k] = -g_yy_0_yyyz_xyyy[k] * ab_z + g_yy_0_yyyz_xyyyz[k];

                g_yy_0_yyyzz_xyyz[k] = -g_yy_0_yyyz_xyyz[k] * ab_z + g_yy_0_yyyz_xyyzz[k];

                g_yy_0_yyyzz_xyzz[k] = -g_yy_0_yyyz_xyzz[k] * ab_z + g_yy_0_yyyz_xyzzz[k];

                g_yy_0_yyyzz_xzzz[k] = -g_yy_0_yyyz_xzzz[k] * ab_z + g_yy_0_yyyz_xzzzz[k];

                g_yy_0_yyyzz_yyyy[k] = -g_yy_0_yyyz_yyyy[k] * ab_z + g_yy_0_yyyz_yyyyz[k];

                g_yy_0_yyyzz_yyyz[k] = -g_yy_0_yyyz_yyyz[k] * ab_z + g_yy_0_yyyz_yyyzz[k];

                g_yy_0_yyyzz_yyzz[k] = -g_yy_0_yyyz_yyzz[k] * ab_z + g_yy_0_yyyz_yyzzz[k];

                g_yy_0_yyyzz_yzzz[k] = -g_yy_0_yyyz_yzzz[k] * ab_z + g_yy_0_yyyz_yzzzz[k];

                g_yy_0_yyyzz_zzzz[k] = -g_yy_0_yyyz_zzzz[k] * ab_z + g_yy_0_yyyz_zzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1215 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1216 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1217 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1218 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1219 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1220 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1221 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1222 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1223 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1224 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1225 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1226 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1227 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1228 * ccomps * dcomps);

            auto g_yy_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzz_xxxx, g_yy_0_yyzz_xxxxz, g_yy_0_yyzz_xxxy, g_yy_0_yyzz_xxxyz, g_yy_0_yyzz_xxxz, g_yy_0_yyzz_xxxzz, g_yy_0_yyzz_xxyy, g_yy_0_yyzz_xxyyz, g_yy_0_yyzz_xxyz, g_yy_0_yyzz_xxyzz, g_yy_0_yyzz_xxzz, g_yy_0_yyzz_xxzzz, g_yy_0_yyzz_xyyy, g_yy_0_yyzz_xyyyz, g_yy_0_yyzz_xyyz, g_yy_0_yyzz_xyyzz, g_yy_0_yyzz_xyzz, g_yy_0_yyzz_xyzzz, g_yy_0_yyzz_xzzz, g_yy_0_yyzz_xzzzz, g_yy_0_yyzz_yyyy, g_yy_0_yyzz_yyyyz, g_yy_0_yyzz_yyyz, g_yy_0_yyzz_yyyzz, g_yy_0_yyzz_yyzz, g_yy_0_yyzz_yyzzz, g_yy_0_yyzz_yzzz, g_yy_0_yyzz_yzzzz, g_yy_0_yyzz_zzzz, g_yy_0_yyzz_zzzzz, g_yy_0_yyzzz_xxxx, g_yy_0_yyzzz_xxxy, g_yy_0_yyzzz_xxxz, g_yy_0_yyzzz_xxyy, g_yy_0_yyzzz_xxyz, g_yy_0_yyzzz_xxzz, g_yy_0_yyzzz_xyyy, g_yy_0_yyzzz_xyyz, g_yy_0_yyzzz_xyzz, g_yy_0_yyzzz_xzzz, g_yy_0_yyzzz_yyyy, g_yy_0_yyzzz_yyyz, g_yy_0_yyzzz_yyzz, g_yy_0_yyzzz_yzzz, g_yy_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzz_xxxx[k] = -g_yy_0_yyzz_xxxx[k] * ab_z + g_yy_0_yyzz_xxxxz[k];

                g_yy_0_yyzzz_xxxy[k] = -g_yy_0_yyzz_xxxy[k] * ab_z + g_yy_0_yyzz_xxxyz[k];

                g_yy_0_yyzzz_xxxz[k] = -g_yy_0_yyzz_xxxz[k] * ab_z + g_yy_0_yyzz_xxxzz[k];

                g_yy_0_yyzzz_xxyy[k] = -g_yy_0_yyzz_xxyy[k] * ab_z + g_yy_0_yyzz_xxyyz[k];

                g_yy_0_yyzzz_xxyz[k] = -g_yy_0_yyzz_xxyz[k] * ab_z + g_yy_0_yyzz_xxyzz[k];

                g_yy_0_yyzzz_xxzz[k] = -g_yy_0_yyzz_xxzz[k] * ab_z + g_yy_0_yyzz_xxzzz[k];

                g_yy_0_yyzzz_xyyy[k] = -g_yy_0_yyzz_xyyy[k] * ab_z + g_yy_0_yyzz_xyyyz[k];

                g_yy_0_yyzzz_xyyz[k] = -g_yy_0_yyzz_xyyz[k] * ab_z + g_yy_0_yyzz_xyyzz[k];

                g_yy_0_yyzzz_xyzz[k] = -g_yy_0_yyzz_xyzz[k] * ab_z + g_yy_0_yyzz_xyzzz[k];

                g_yy_0_yyzzz_xzzz[k] = -g_yy_0_yyzz_xzzz[k] * ab_z + g_yy_0_yyzz_xzzzz[k];

                g_yy_0_yyzzz_yyyy[k] = -g_yy_0_yyzz_yyyy[k] * ab_z + g_yy_0_yyzz_yyyyz[k];

                g_yy_0_yyzzz_yyyz[k] = -g_yy_0_yyzz_yyyz[k] * ab_z + g_yy_0_yyzz_yyyzz[k];

                g_yy_0_yyzzz_yyzz[k] = -g_yy_0_yyzz_yyzz[k] * ab_z + g_yy_0_yyzz_yyzzz[k];

                g_yy_0_yyzzz_yzzz[k] = -g_yy_0_yyzz_yzzz[k] * ab_z + g_yy_0_yyzz_yzzzz[k];

                g_yy_0_yyzzz_zzzz[k] = -g_yy_0_yyzz_zzzz[k] * ab_z + g_yy_0_yyzz_zzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1230 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1231 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1232 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1233 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1234 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1235 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1236 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1237 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1238 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1239 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1240 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1241 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1242 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1243 * ccomps * dcomps);

            auto g_yy_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzz_xxxx, g_yy_0_yzzz_xxxxz, g_yy_0_yzzz_xxxy, g_yy_0_yzzz_xxxyz, g_yy_0_yzzz_xxxz, g_yy_0_yzzz_xxxzz, g_yy_0_yzzz_xxyy, g_yy_0_yzzz_xxyyz, g_yy_0_yzzz_xxyz, g_yy_0_yzzz_xxyzz, g_yy_0_yzzz_xxzz, g_yy_0_yzzz_xxzzz, g_yy_0_yzzz_xyyy, g_yy_0_yzzz_xyyyz, g_yy_0_yzzz_xyyz, g_yy_0_yzzz_xyyzz, g_yy_0_yzzz_xyzz, g_yy_0_yzzz_xyzzz, g_yy_0_yzzz_xzzz, g_yy_0_yzzz_xzzzz, g_yy_0_yzzz_yyyy, g_yy_0_yzzz_yyyyz, g_yy_0_yzzz_yyyz, g_yy_0_yzzz_yyyzz, g_yy_0_yzzz_yyzz, g_yy_0_yzzz_yyzzz, g_yy_0_yzzz_yzzz, g_yy_0_yzzz_yzzzz, g_yy_0_yzzz_zzzz, g_yy_0_yzzz_zzzzz, g_yy_0_yzzzz_xxxx, g_yy_0_yzzzz_xxxy, g_yy_0_yzzzz_xxxz, g_yy_0_yzzzz_xxyy, g_yy_0_yzzzz_xxyz, g_yy_0_yzzzz_xxzz, g_yy_0_yzzzz_xyyy, g_yy_0_yzzzz_xyyz, g_yy_0_yzzzz_xyzz, g_yy_0_yzzzz_xzzz, g_yy_0_yzzzz_yyyy, g_yy_0_yzzzz_yyyz, g_yy_0_yzzzz_yyzz, g_yy_0_yzzzz_yzzz, g_yy_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzz_xxxx[k] = -g_yy_0_yzzz_xxxx[k] * ab_z + g_yy_0_yzzz_xxxxz[k];

                g_yy_0_yzzzz_xxxy[k] = -g_yy_0_yzzz_xxxy[k] * ab_z + g_yy_0_yzzz_xxxyz[k];

                g_yy_0_yzzzz_xxxz[k] = -g_yy_0_yzzz_xxxz[k] * ab_z + g_yy_0_yzzz_xxxzz[k];

                g_yy_0_yzzzz_xxyy[k] = -g_yy_0_yzzz_xxyy[k] * ab_z + g_yy_0_yzzz_xxyyz[k];

                g_yy_0_yzzzz_xxyz[k] = -g_yy_0_yzzz_xxyz[k] * ab_z + g_yy_0_yzzz_xxyzz[k];

                g_yy_0_yzzzz_xxzz[k] = -g_yy_0_yzzz_xxzz[k] * ab_z + g_yy_0_yzzz_xxzzz[k];

                g_yy_0_yzzzz_xyyy[k] = -g_yy_0_yzzz_xyyy[k] * ab_z + g_yy_0_yzzz_xyyyz[k];

                g_yy_0_yzzzz_xyyz[k] = -g_yy_0_yzzz_xyyz[k] * ab_z + g_yy_0_yzzz_xyyzz[k];

                g_yy_0_yzzzz_xyzz[k] = -g_yy_0_yzzz_xyzz[k] * ab_z + g_yy_0_yzzz_xyzzz[k];

                g_yy_0_yzzzz_xzzz[k] = -g_yy_0_yzzz_xzzz[k] * ab_z + g_yy_0_yzzz_xzzzz[k];

                g_yy_0_yzzzz_yyyy[k] = -g_yy_0_yzzz_yyyy[k] * ab_z + g_yy_0_yzzz_yyyyz[k];

                g_yy_0_yzzzz_yyyz[k] = -g_yy_0_yzzz_yyyz[k] * ab_z + g_yy_0_yzzz_yyyzz[k];

                g_yy_0_yzzzz_yyzz[k] = -g_yy_0_yzzz_yyzz[k] * ab_z + g_yy_0_yzzz_yyzzz[k];

                g_yy_0_yzzzz_yzzz[k] = -g_yy_0_yzzz_yzzz[k] * ab_z + g_yy_0_yzzz_yzzzz[k];

                g_yy_0_yzzzz_zzzz[k] = -g_yy_0_yzzz_zzzz[k] * ab_z + g_yy_0_yzzz_zzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1245 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1246 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1247 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1248 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1249 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1250 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1251 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1252 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1253 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1254 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1255 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1256 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1257 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1258 * ccomps * dcomps);

            auto g_yy_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzz_xxxx, g_yy_0_zzzz_xxxxz, g_yy_0_zzzz_xxxy, g_yy_0_zzzz_xxxyz, g_yy_0_zzzz_xxxz, g_yy_0_zzzz_xxxzz, g_yy_0_zzzz_xxyy, g_yy_0_zzzz_xxyyz, g_yy_0_zzzz_xxyz, g_yy_0_zzzz_xxyzz, g_yy_0_zzzz_xxzz, g_yy_0_zzzz_xxzzz, g_yy_0_zzzz_xyyy, g_yy_0_zzzz_xyyyz, g_yy_0_zzzz_xyyz, g_yy_0_zzzz_xyyzz, g_yy_0_zzzz_xyzz, g_yy_0_zzzz_xyzzz, g_yy_0_zzzz_xzzz, g_yy_0_zzzz_xzzzz, g_yy_0_zzzz_yyyy, g_yy_0_zzzz_yyyyz, g_yy_0_zzzz_yyyz, g_yy_0_zzzz_yyyzz, g_yy_0_zzzz_yyzz, g_yy_0_zzzz_yyzzz, g_yy_0_zzzz_yzzz, g_yy_0_zzzz_yzzzz, g_yy_0_zzzz_zzzz, g_yy_0_zzzz_zzzzz, g_yy_0_zzzzz_xxxx, g_yy_0_zzzzz_xxxy, g_yy_0_zzzzz_xxxz, g_yy_0_zzzzz_xxyy, g_yy_0_zzzzz_xxyz, g_yy_0_zzzzz_xxzz, g_yy_0_zzzzz_xyyy, g_yy_0_zzzzz_xyyz, g_yy_0_zzzzz_xyzz, g_yy_0_zzzzz_xzzz, g_yy_0_zzzzz_yyyy, g_yy_0_zzzzz_yyyz, g_yy_0_zzzzz_yyzz, g_yy_0_zzzzz_yzzz, g_yy_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzz_xxxx[k] = -g_yy_0_zzzz_xxxx[k] * ab_z + g_yy_0_zzzz_xxxxz[k];

                g_yy_0_zzzzz_xxxy[k] = -g_yy_0_zzzz_xxxy[k] * ab_z + g_yy_0_zzzz_xxxyz[k];

                g_yy_0_zzzzz_xxxz[k] = -g_yy_0_zzzz_xxxz[k] * ab_z + g_yy_0_zzzz_xxxzz[k];

                g_yy_0_zzzzz_xxyy[k] = -g_yy_0_zzzz_xxyy[k] * ab_z + g_yy_0_zzzz_xxyyz[k];

                g_yy_0_zzzzz_xxyz[k] = -g_yy_0_zzzz_xxyz[k] * ab_z + g_yy_0_zzzz_xxyzz[k];

                g_yy_0_zzzzz_xxzz[k] = -g_yy_0_zzzz_xxzz[k] * ab_z + g_yy_0_zzzz_xxzzz[k];

                g_yy_0_zzzzz_xyyy[k] = -g_yy_0_zzzz_xyyy[k] * ab_z + g_yy_0_zzzz_xyyyz[k];

                g_yy_0_zzzzz_xyyz[k] = -g_yy_0_zzzz_xyyz[k] * ab_z + g_yy_0_zzzz_xyyzz[k];

                g_yy_0_zzzzz_xyzz[k] = -g_yy_0_zzzz_xyzz[k] * ab_z + g_yy_0_zzzz_xyzzz[k];

                g_yy_0_zzzzz_xzzz[k] = -g_yy_0_zzzz_xzzz[k] * ab_z + g_yy_0_zzzz_xzzzz[k];

                g_yy_0_zzzzz_yyyy[k] = -g_yy_0_zzzz_yyyy[k] * ab_z + g_yy_0_zzzz_yyyyz[k];

                g_yy_0_zzzzz_yyyz[k] = -g_yy_0_zzzz_yyyz[k] * ab_z + g_yy_0_zzzz_yyyzz[k];

                g_yy_0_zzzzz_yyzz[k] = -g_yy_0_zzzz_yyzz[k] * ab_z + g_yy_0_zzzz_yyzzz[k];

                g_yy_0_zzzzz_yzzz[k] = -g_yy_0_zzzz_yzzz[k] * ab_z + g_yy_0_zzzz_yzzzz[k];

                g_yy_0_zzzzz_zzzz[k] = -g_yy_0_zzzz_zzzz[k] * ab_z + g_yy_0_zzzz_zzzzz[k];
            }

            /// Set up 1260-1275 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 1260 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 1261 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 1262 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 1263 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 1264 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 1265 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 1266 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 1267 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 1268 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 1269 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 1270 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 1271 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 1272 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 1273 * ccomps * dcomps);

            auto g_yz_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 1274 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxx_xxxx, g_yz_0_xxxx_xxxxx, g_yz_0_xxxx_xxxxy, g_yz_0_xxxx_xxxxz, g_yz_0_xxxx_xxxy, g_yz_0_xxxx_xxxyy, g_yz_0_xxxx_xxxyz, g_yz_0_xxxx_xxxz, g_yz_0_xxxx_xxxzz, g_yz_0_xxxx_xxyy, g_yz_0_xxxx_xxyyy, g_yz_0_xxxx_xxyyz, g_yz_0_xxxx_xxyz, g_yz_0_xxxx_xxyzz, g_yz_0_xxxx_xxzz, g_yz_0_xxxx_xxzzz, g_yz_0_xxxx_xyyy, g_yz_0_xxxx_xyyyy, g_yz_0_xxxx_xyyyz, g_yz_0_xxxx_xyyz, g_yz_0_xxxx_xyyzz, g_yz_0_xxxx_xyzz, g_yz_0_xxxx_xyzzz, g_yz_0_xxxx_xzzz, g_yz_0_xxxx_xzzzz, g_yz_0_xxxx_yyyy, g_yz_0_xxxx_yyyz, g_yz_0_xxxx_yyzz, g_yz_0_xxxx_yzzz, g_yz_0_xxxx_zzzz, g_yz_0_xxxxx_xxxx, g_yz_0_xxxxx_xxxy, g_yz_0_xxxxx_xxxz, g_yz_0_xxxxx_xxyy, g_yz_0_xxxxx_xxyz, g_yz_0_xxxxx_xxzz, g_yz_0_xxxxx_xyyy, g_yz_0_xxxxx_xyyz, g_yz_0_xxxxx_xyzz, g_yz_0_xxxxx_xzzz, g_yz_0_xxxxx_yyyy, g_yz_0_xxxxx_yyyz, g_yz_0_xxxxx_yyzz, g_yz_0_xxxxx_yzzz, g_yz_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxx_xxxx[k] = -g_yz_0_xxxx_xxxx[k] * ab_x + g_yz_0_xxxx_xxxxx[k];

                g_yz_0_xxxxx_xxxy[k] = -g_yz_0_xxxx_xxxy[k] * ab_x + g_yz_0_xxxx_xxxxy[k];

                g_yz_0_xxxxx_xxxz[k] = -g_yz_0_xxxx_xxxz[k] * ab_x + g_yz_0_xxxx_xxxxz[k];

                g_yz_0_xxxxx_xxyy[k] = -g_yz_0_xxxx_xxyy[k] * ab_x + g_yz_0_xxxx_xxxyy[k];

                g_yz_0_xxxxx_xxyz[k] = -g_yz_0_xxxx_xxyz[k] * ab_x + g_yz_0_xxxx_xxxyz[k];

                g_yz_0_xxxxx_xxzz[k] = -g_yz_0_xxxx_xxzz[k] * ab_x + g_yz_0_xxxx_xxxzz[k];

                g_yz_0_xxxxx_xyyy[k] = -g_yz_0_xxxx_xyyy[k] * ab_x + g_yz_0_xxxx_xxyyy[k];

                g_yz_0_xxxxx_xyyz[k] = -g_yz_0_xxxx_xyyz[k] * ab_x + g_yz_0_xxxx_xxyyz[k];

                g_yz_0_xxxxx_xyzz[k] = -g_yz_0_xxxx_xyzz[k] * ab_x + g_yz_0_xxxx_xxyzz[k];

                g_yz_0_xxxxx_xzzz[k] = -g_yz_0_xxxx_xzzz[k] * ab_x + g_yz_0_xxxx_xxzzz[k];

                g_yz_0_xxxxx_yyyy[k] = -g_yz_0_xxxx_yyyy[k] * ab_x + g_yz_0_xxxx_xyyyy[k];

                g_yz_0_xxxxx_yyyz[k] = -g_yz_0_xxxx_yyyz[k] * ab_x + g_yz_0_xxxx_xyyyz[k];

                g_yz_0_xxxxx_yyzz[k] = -g_yz_0_xxxx_yyzz[k] * ab_x + g_yz_0_xxxx_xyyzz[k];

                g_yz_0_xxxxx_yzzz[k] = -g_yz_0_xxxx_yzzz[k] * ab_x + g_yz_0_xxxx_xyzzz[k];

                g_yz_0_xxxxx_zzzz[k] = -g_yz_0_xxxx_zzzz[k] * ab_x + g_yz_0_xxxx_xzzzz[k];
            }

            /// Set up 1275-1290 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 1275 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 1276 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 1277 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 1278 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 1279 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 1280 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 1281 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 1282 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 1283 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 1284 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 1285 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 1286 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 1287 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 1288 * ccomps * dcomps);

            auto g_yz_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 1289 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxy_xxxx, g_yz_0_xxxxy_xxxy, g_yz_0_xxxxy_xxxz, g_yz_0_xxxxy_xxyy, g_yz_0_xxxxy_xxyz, g_yz_0_xxxxy_xxzz, g_yz_0_xxxxy_xyyy, g_yz_0_xxxxy_xyyz, g_yz_0_xxxxy_xyzz, g_yz_0_xxxxy_xzzz, g_yz_0_xxxxy_yyyy, g_yz_0_xxxxy_yyyz, g_yz_0_xxxxy_yyzz, g_yz_0_xxxxy_yzzz, g_yz_0_xxxxy_zzzz, g_yz_0_xxxy_xxxx, g_yz_0_xxxy_xxxxx, g_yz_0_xxxy_xxxxy, g_yz_0_xxxy_xxxxz, g_yz_0_xxxy_xxxy, g_yz_0_xxxy_xxxyy, g_yz_0_xxxy_xxxyz, g_yz_0_xxxy_xxxz, g_yz_0_xxxy_xxxzz, g_yz_0_xxxy_xxyy, g_yz_0_xxxy_xxyyy, g_yz_0_xxxy_xxyyz, g_yz_0_xxxy_xxyz, g_yz_0_xxxy_xxyzz, g_yz_0_xxxy_xxzz, g_yz_0_xxxy_xxzzz, g_yz_0_xxxy_xyyy, g_yz_0_xxxy_xyyyy, g_yz_0_xxxy_xyyyz, g_yz_0_xxxy_xyyz, g_yz_0_xxxy_xyyzz, g_yz_0_xxxy_xyzz, g_yz_0_xxxy_xyzzz, g_yz_0_xxxy_xzzz, g_yz_0_xxxy_xzzzz, g_yz_0_xxxy_yyyy, g_yz_0_xxxy_yyyz, g_yz_0_xxxy_yyzz, g_yz_0_xxxy_yzzz, g_yz_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxy_xxxx[k] = -g_yz_0_xxxy_xxxx[k] * ab_x + g_yz_0_xxxy_xxxxx[k];

                g_yz_0_xxxxy_xxxy[k] = -g_yz_0_xxxy_xxxy[k] * ab_x + g_yz_0_xxxy_xxxxy[k];

                g_yz_0_xxxxy_xxxz[k] = -g_yz_0_xxxy_xxxz[k] * ab_x + g_yz_0_xxxy_xxxxz[k];

                g_yz_0_xxxxy_xxyy[k] = -g_yz_0_xxxy_xxyy[k] * ab_x + g_yz_0_xxxy_xxxyy[k];

                g_yz_0_xxxxy_xxyz[k] = -g_yz_0_xxxy_xxyz[k] * ab_x + g_yz_0_xxxy_xxxyz[k];

                g_yz_0_xxxxy_xxzz[k] = -g_yz_0_xxxy_xxzz[k] * ab_x + g_yz_0_xxxy_xxxzz[k];

                g_yz_0_xxxxy_xyyy[k] = -g_yz_0_xxxy_xyyy[k] * ab_x + g_yz_0_xxxy_xxyyy[k];

                g_yz_0_xxxxy_xyyz[k] = -g_yz_0_xxxy_xyyz[k] * ab_x + g_yz_0_xxxy_xxyyz[k];

                g_yz_0_xxxxy_xyzz[k] = -g_yz_0_xxxy_xyzz[k] * ab_x + g_yz_0_xxxy_xxyzz[k];

                g_yz_0_xxxxy_xzzz[k] = -g_yz_0_xxxy_xzzz[k] * ab_x + g_yz_0_xxxy_xxzzz[k];

                g_yz_0_xxxxy_yyyy[k] = -g_yz_0_xxxy_yyyy[k] * ab_x + g_yz_0_xxxy_xyyyy[k];

                g_yz_0_xxxxy_yyyz[k] = -g_yz_0_xxxy_yyyz[k] * ab_x + g_yz_0_xxxy_xyyyz[k];

                g_yz_0_xxxxy_yyzz[k] = -g_yz_0_xxxy_yyzz[k] * ab_x + g_yz_0_xxxy_xyyzz[k];

                g_yz_0_xxxxy_yzzz[k] = -g_yz_0_xxxy_yzzz[k] * ab_x + g_yz_0_xxxy_xyzzz[k];

                g_yz_0_xxxxy_zzzz[k] = -g_yz_0_xxxy_zzzz[k] * ab_x + g_yz_0_xxxy_xzzzz[k];
            }

            /// Set up 1290-1305 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 1290 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 1291 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 1292 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 1293 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 1294 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 1295 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 1296 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 1297 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 1298 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 1299 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 1300 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 1301 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 1302 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 1303 * ccomps * dcomps);

            auto g_yz_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 1304 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxz_xxxx, g_yz_0_xxxxz_xxxy, g_yz_0_xxxxz_xxxz, g_yz_0_xxxxz_xxyy, g_yz_0_xxxxz_xxyz, g_yz_0_xxxxz_xxzz, g_yz_0_xxxxz_xyyy, g_yz_0_xxxxz_xyyz, g_yz_0_xxxxz_xyzz, g_yz_0_xxxxz_xzzz, g_yz_0_xxxxz_yyyy, g_yz_0_xxxxz_yyyz, g_yz_0_xxxxz_yyzz, g_yz_0_xxxxz_yzzz, g_yz_0_xxxxz_zzzz, g_yz_0_xxxz_xxxx, g_yz_0_xxxz_xxxxx, g_yz_0_xxxz_xxxxy, g_yz_0_xxxz_xxxxz, g_yz_0_xxxz_xxxy, g_yz_0_xxxz_xxxyy, g_yz_0_xxxz_xxxyz, g_yz_0_xxxz_xxxz, g_yz_0_xxxz_xxxzz, g_yz_0_xxxz_xxyy, g_yz_0_xxxz_xxyyy, g_yz_0_xxxz_xxyyz, g_yz_0_xxxz_xxyz, g_yz_0_xxxz_xxyzz, g_yz_0_xxxz_xxzz, g_yz_0_xxxz_xxzzz, g_yz_0_xxxz_xyyy, g_yz_0_xxxz_xyyyy, g_yz_0_xxxz_xyyyz, g_yz_0_xxxz_xyyz, g_yz_0_xxxz_xyyzz, g_yz_0_xxxz_xyzz, g_yz_0_xxxz_xyzzz, g_yz_0_xxxz_xzzz, g_yz_0_xxxz_xzzzz, g_yz_0_xxxz_yyyy, g_yz_0_xxxz_yyyz, g_yz_0_xxxz_yyzz, g_yz_0_xxxz_yzzz, g_yz_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxz_xxxx[k] = -g_yz_0_xxxz_xxxx[k] * ab_x + g_yz_0_xxxz_xxxxx[k];

                g_yz_0_xxxxz_xxxy[k] = -g_yz_0_xxxz_xxxy[k] * ab_x + g_yz_0_xxxz_xxxxy[k];

                g_yz_0_xxxxz_xxxz[k] = -g_yz_0_xxxz_xxxz[k] * ab_x + g_yz_0_xxxz_xxxxz[k];

                g_yz_0_xxxxz_xxyy[k] = -g_yz_0_xxxz_xxyy[k] * ab_x + g_yz_0_xxxz_xxxyy[k];

                g_yz_0_xxxxz_xxyz[k] = -g_yz_0_xxxz_xxyz[k] * ab_x + g_yz_0_xxxz_xxxyz[k];

                g_yz_0_xxxxz_xxzz[k] = -g_yz_0_xxxz_xxzz[k] * ab_x + g_yz_0_xxxz_xxxzz[k];

                g_yz_0_xxxxz_xyyy[k] = -g_yz_0_xxxz_xyyy[k] * ab_x + g_yz_0_xxxz_xxyyy[k];

                g_yz_0_xxxxz_xyyz[k] = -g_yz_0_xxxz_xyyz[k] * ab_x + g_yz_0_xxxz_xxyyz[k];

                g_yz_0_xxxxz_xyzz[k] = -g_yz_0_xxxz_xyzz[k] * ab_x + g_yz_0_xxxz_xxyzz[k];

                g_yz_0_xxxxz_xzzz[k] = -g_yz_0_xxxz_xzzz[k] * ab_x + g_yz_0_xxxz_xxzzz[k];

                g_yz_0_xxxxz_yyyy[k] = -g_yz_0_xxxz_yyyy[k] * ab_x + g_yz_0_xxxz_xyyyy[k];

                g_yz_0_xxxxz_yyyz[k] = -g_yz_0_xxxz_yyyz[k] * ab_x + g_yz_0_xxxz_xyyyz[k];

                g_yz_0_xxxxz_yyzz[k] = -g_yz_0_xxxz_yyzz[k] * ab_x + g_yz_0_xxxz_xyyzz[k];

                g_yz_0_xxxxz_yzzz[k] = -g_yz_0_xxxz_yzzz[k] * ab_x + g_yz_0_xxxz_xyzzz[k];

                g_yz_0_xxxxz_zzzz[k] = -g_yz_0_xxxz_zzzz[k] * ab_x + g_yz_0_xxxz_xzzzz[k];
            }

            /// Set up 1305-1320 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 1305 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 1306 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 1307 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 1308 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 1309 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 1310 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 1311 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 1312 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 1313 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 1314 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 1315 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 1316 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 1317 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 1318 * ccomps * dcomps);

            auto g_yz_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 1319 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyy_xxxx, g_yz_0_xxxyy_xxxy, g_yz_0_xxxyy_xxxz, g_yz_0_xxxyy_xxyy, g_yz_0_xxxyy_xxyz, g_yz_0_xxxyy_xxzz, g_yz_0_xxxyy_xyyy, g_yz_0_xxxyy_xyyz, g_yz_0_xxxyy_xyzz, g_yz_0_xxxyy_xzzz, g_yz_0_xxxyy_yyyy, g_yz_0_xxxyy_yyyz, g_yz_0_xxxyy_yyzz, g_yz_0_xxxyy_yzzz, g_yz_0_xxxyy_zzzz, g_yz_0_xxyy_xxxx, g_yz_0_xxyy_xxxxx, g_yz_0_xxyy_xxxxy, g_yz_0_xxyy_xxxxz, g_yz_0_xxyy_xxxy, g_yz_0_xxyy_xxxyy, g_yz_0_xxyy_xxxyz, g_yz_0_xxyy_xxxz, g_yz_0_xxyy_xxxzz, g_yz_0_xxyy_xxyy, g_yz_0_xxyy_xxyyy, g_yz_0_xxyy_xxyyz, g_yz_0_xxyy_xxyz, g_yz_0_xxyy_xxyzz, g_yz_0_xxyy_xxzz, g_yz_0_xxyy_xxzzz, g_yz_0_xxyy_xyyy, g_yz_0_xxyy_xyyyy, g_yz_0_xxyy_xyyyz, g_yz_0_xxyy_xyyz, g_yz_0_xxyy_xyyzz, g_yz_0_xxyy_xyzz, g_yz_0_xxyy_xyzzz, g_yz_0_xxyy_xzzz, g_yz_0_xxyy_xzzzz, g_yz_0_xxyy_yyyy, g_yz_0_xxyy_yyyz, g_yz_0_xxyy_yyzz, g_yz_0_xxyy_yzzz, g_yz_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyy_xxxx[k] = -g_yz_0_xxyy_xxxx[k] * ab_x + g_yz_0_xxyy_xxxxx[k];

                g_yz_0_xxxyy_xxxy[k] = -g_yz_0_xxyy_xxxy[k] * ab_x + g_yz_0_xxyy_xxxxy[k];

                g_yz_0_xxxyy_xxxz[k] = -g_yz_0_xxyy_xxxz[k] * ab_x + g_yz_0_xxyy_xxxxz[k];

                g_yz_0_xxxyy_xxyy[k] = -g_yz_0_xxyy_xxyy[k] * ab_x + g_yz_0_xxyy_xxxyy[k];

                g_yz_0_xxxyy_xxyz[k] = -g_yz_0_xxyy_xxyz[k] * ab_x + g_yz_0_xxyy_xxxyz[k];

                g_yz_0_xxxyy_xxzz[k] = -g_yz_0_xxyy_xxzz[k] * ab_x + g_yz_0_xxyy_xxxzz[k];

                g_yz_0_xxxyy_xyyy[k] = -g_yz_0_xxyy_xyyy[k] * ab_x + g_yz_0_xxyy_xxyyy[k];

                g_yz_0_xxxyy_xyyz[k] = -g_yz_0_xxyy_xyyz[k] * ab_x + g_yz_0_xxyy_xxyyz[k];

                g_yz_0_xxxyy_xyzz[k] = -g_yz_0_xxyy_xyzz[k] * ab_x + g_yz_0_xxyy_xxyzz[k];

                g_yz_0_xxxyy_xzzz[k] = -g_yz_0_xxyy_xzzz[k] * ab_x + g_yz_0_xxyy_xxzzz[k];

                g_yz_0_xxxyy_yyyy[k] = -g_yz_0_xxyy_yyyy[k] * ab_x + g_yz_0_xxyy_xyyyy[k];

                g_yz_0_xxxyy_yyyz[k] = -g_yz_0_xxyy_yyyz[k] * ab_x + g_yz_0_xxyy_xyyyz[k];

                g_yz_0_xxxyy_yyzz[k] = -g_yz_0_xxyy_yyzz[k] * ab_x + g_yz_0_xxyy_xyyzz[k];

                g_yz_0_xxxyy_yzzz[k] = -g_yz_0_xxyy_yzzz[k] * ab_x + g_yz_0_xxyy_xyzzz[k];

                g_yz_0_xxxyy_zzzz[k] = -g_yz_0_xxyy_zzzz[k] * ab_x + g_yz_0_xxyy_xzzzz[k];
            }

            /// Set up 1320-1335 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 1320 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 1321 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 1322 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 1323 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 1324 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 1325 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 1326 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 1327 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 1328 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 1329 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 1330 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 1331 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 1332 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 1333 * ccomps * dcomps);

            auto g_yz_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 1334 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyz_xxxx, g_yz_0_xxxyz_xxxy, g_yz_0_xxxyz_xxxz, g_yz_0_xxxyz_xxyy, g_yz_0_xxxyz_xxyz, g_yz_0_xxxyz_xxzz, g_yz_0_xxxyz_xyyy, g_yz_0_xxxyz_xyyz, g_yz_0_xxxyz_xyzz, g_yz_0_xxxyz_xzzz, g_yz_0_xxxyz_yyyy, g_yz_0_xxxyz_yyyz, g_yz_0_xxxyz_yyzz, g_yz_0_xxxyz_yzzz, g_yz_0_xxxyz_zzzz, g_yz_0_xxyz_xxxx, g_yz_0_xxyz_xxxxx, g_yz_0_xxyz_xxxxy, g_yz_0_xxyz_xxxxz, g_yz_0_xxyz_xxxy, g_yz_0_xxyz_xxxyy, g_yz_0_xxyz_xxxyz, g_yz_0_xxyz_xxxz, g_yz_0_xxyz_xxxzz, g_yz_0_xxyz_xxyy, g_yz_0_xxyz_xxyyy, g_yz_0_xxyz_xxyyz, g_yz_0_xxyz_xxyz, g_yz_0_xxyz_xxyzz, g_yz_0_xxyz_xxzz, g_yz_0_xxyz_xxzzz, g_yz_0_xxyz_xyyy, g_yz_0_xxyz_xyyyy, g_yz_0_xxyz_xyyyz, g_yz_0_xxyz_xyyz, g_yz_0_xxyz_xyyzz, g_yz_0_xxyz_xyzz, g_yz_0_xxyz_xyzzz, g_yz_0_xxyz_xzzz, g_yz_0_xxyz_xzzzz, g_yz_0_xxyz_yyyy, g_yz_0_xxyz_yyyz, g_yz_0_xxyz_yyzz, g_yz_0_xxyz_yzzz, g_yz_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyz_xxxx[k] = -g_yz_0_xxyz_xxxx[k] * ab_x + g_yz_0_xxyz_xxxxx[k];

                g_yz_0_xxxyz_xxxy[k] = -g_yz_0_xxyz_xxxy[k] * ab_x + g_yz_0_xxyz_xxxxy[k];

                g_yz_0_xxxyz_xxxz[k] = -g_yz_0_xxyz_xxxz[k] * ab_x + g_yz_0_xxyz_xxxxz[k];

                g_yz_0_xxxyz_xxyy[k] = -g_yz_0_xxyz_xxyy[k] * ab_x + g_yz_0_xxyz_xxxyy[k];

                g_yz_0_xxxyz_xxyz[k] = -g_yz_0_xxyz_xxyz[k] * ab_x + g_yz_0_xxyz_xxxyz[k];

                g_yz_0_xxxyz_xxzz[k] = -g_yz_0_xxyz_xxzz[k] * ab_x + g_yz_0_xxyz_xxxzz[k];

                g_yz_0_xxxyz_xyyy[k] = -g_yz_0_xxyz_xyyy[k] * ab_x + g_yz_0_xxyz_xxyyy[k];

                g_yz_0_xxxyz_xyyz[k] = -g_yz_0_xxyz_xyyz[k] * ab_x + g_yz_0_xxyz_xxyyz[k];

                g_yz_0_xxxyz_xyzz[k] = -g_yz_0_xxyz_xyzz[k] * ab_x + g_yz_0_xxyz_xxyzz[k];

                g_yz_0_xxxyz_xzzz[k] = -g_yz_0_xxyz_xzzz[k] * ab_x + g_yz_0_xxyz_xxzzz[k];

                g_yz_0_xxxyz_yyyy[k] = -g_yz_0_xxyz_yyyy[k] * ab_x + g_yz_0_xxyz_xyyyy[k];

                g_yz_0_xxxyz_yyyz[k] = -g_yz_0_xxyz_yyyz[k] * ab_x + g_yz_0_xxyz_xyyyz[k];

                g_yz_0_xxxyz_yyzz[k] = -g_yz_0_xxyz_yyzz[k] * ab_x + g_yz_0_xxyz_xyyzz[k];

                g_yz_0_xxxyz_yzzz[k] = -g_yz_0_xxyz_yzzz[k] * ab_x + g_yz_0_xxyz_xyzzz[k];

                g_yz_0_xxxyz_zzzz[k] = -g_yz_0_xxyz_zzzz[k] * ab_x + g_yz_0_xxyz_xzzzz[k];
            }

            /// Set up 1335-1350 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 1335 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 1336 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 1337 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 1338 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 1339 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 1340 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 1341 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 1342 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 1343 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 1344 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 1345 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 1346 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 1347 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 1348 * ccomps * dcomps);

            auto g_yz_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 1349 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzz_xxxx, g_yz_0_xxxzz_xxxy, g_yz_0_xxxzz_xxxz, g_yz_0_xxxzz_xxyy, g_yz_0_xxxzz_xxyz, g_yz_0_xxxzz_xxzz, g_yz_0_xxxzz_xyyy, g_yz_0_xxxzz_xyyz, g_yz_0_xxxzz_xyzz, g_yz_0_xxxzz_xzzz, g_yz_0_xxxzz_yyyy, g_yz_0_xxxzz_yyyz, g_yz_0_xxxzz_yyzz, g_yz_0_xxxzz_yzzz, g_yz_0_xxxzz_zzzz, g_yz_0_xxzz_xxxx, g_yz_0_xxzz_xxxxx, g_yz_0_xxzz_xxxxy, g_yz_0_xxzz_xxxxz, g_yz_0_xxzz_xxxy, g_yz_0_xxzz_xxxyy, g_yz_0_xxzz_xxxyz, g_yz_0_xxzz_xxxz, g_yz_0_xxzz_xxxzz, g_yz_0_xxzz_xxyy, g_yz_0_xxzz_xxyyy, g_yz_0_xxzz_xxyyz, g_yz_0_xxzz_xxyz, g_yz_0_xxzz_xxyzz, g_yz_0_xxzz_xxzz, g_yz_0_xxzz_xxzzz, g_yz_0_xxzz_xyyy, g_yz_0_xxzz_xyyyy, g_yz_0_xxzz_xyyyz, g_yz_0_xxzz_xyyz, g_yz_0_xxzz_xyyzz, g_yz_0_xxzz_xyzz, g_yz_0_xxzz_xyzzz, g_yz_0_xxzz_xzzz, g_yz_0_xxzz_xzzzz, g_yz_0_xxzz_yyyy, g_yz_0_xxzz_yyyz, g_yz_0_xxzz_yyzz, g_yz_0_xxzz_yzzz, g_yz_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzz_xxxx[k] = -g_yz_0_xxzz_xxxx[k] * ab_x + g_yz_0_xxzz_xxxxx[k];

                g_yz_0_xxxzz_xxxy[k] = -g_yz_0_xxzz_xxxy[k] * ab_x + g_yz_0_xxzz_xxxxy[k];

                g_yz_0_xxxzz_xxxz[k] = -g_yz_0_xxzz_xxxz[k] * ab_x + g_yz_0_xxzz_xxxxz[k];

                g_yz_0_xxxzz_xxyy[k] = -g_yz_0_xxzz_xxyy[k] * ab_x + g_yz_0_xxzz_xxxyy[k];

                g_yz_0_xxxzz_xxyz[k] = -g_yz_0_xxzz_xxyz[k] * ab_x + g_yz_0_xxzz_xxxyz[k];

                g_yz_0_xxxzz_xxzz[k] = -g_yz_0_xxzz_xxzz[k] * ab_x + g_yz_0_xxzz_xxxzz[k];

                g_yz_0_xxxzz_xyyy[k] = -g_yz_0_xxzz_xyyy[k] * ab_x + g_yz_0_xxzz_xxyyy[k];

                g_yz_0_xxxzz_xyyz[k] = -g_yz_0_xxzz_xyyz[k] * ab_x + g_yz_0_xxzz_xxyyz[k];

                g_yz_0_xxxzz_xyzz[k] = -g_yz_0_xxzz_xyzz[k] * ab_x + g_yz_0_xxzz_xxyzz[k];

                g_yz_0_xxxzz_xzzz[k] = -g_yz_0_xxzz_xzzz[k] * ab_x + g_yz_0_xxzz_xxzzz[k];

                g_yz_0_xxxzz_yyyy[k] = -g_yz_0_xxzz_yyyy[k] * ab_x + g_yz_0_xxzz_xyyyy[k];

                g_yz_0_xxxzz_yyyz[k] = -g_yz_0_xxzz_yyyz[k] * ab_x + g_yz_0_xxzz_xyyyz[k];

                g_yz_0_xxxzz_yyzz[k] = -g_yz_0_xxzz_yyzz[k] * ab_x + g_yz_0_xxzz_xyyzz[k];

                g_yz_0_xxxzz_yzzz[k] = -g_yz_0_xxzz_yzzz[k] * ab_x + g_yz_0_xxzz_xyzzz[k];

                g_yz_0_xxxzz_zzzz[k] = -g_yz_0_xxzz_zzzz[k] * ab_x + g_yz_0_xxzz_xzzzz[k];
            }

            /// Set up 1350-1365 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 1350 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 1351 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 1352 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 1353 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 1354 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 1355 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 1356 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 1357 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 1358 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 1359 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 1360 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 1361 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 1362 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 1363 * ccomps * dcomps);

            auto g_yz_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 1364 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyy_xxxx, g_yz_0_xxyyy_xxxy, g_yz_0_xxyyy_xxxz, g_yz_0_xxyyy_xxyy, g_yz_0_xxyyy_xxyz, g_yz_0_xxyyy_xxzz, g_yz_0_xxyyy_xyyy, g_yz_0_xxyyy_xyyz, g_yz_0_xxyyy_xyzz, g_yz_0_xxyyy_xzzz, g_yz_0_xxyyy_yyyy, g_yz_0_xxyyy_yyyz, g_yz_0_xxyyy_yyzz, g_yz_0_xxyyy_yzzz, g_yz_0_xxyyy_zzzz, g_yz_0_xyyy_xxxx, g_yz_0_xyyy_xxxxx, g_yz_0_xyyy_xxxxy, g_yz_0_xyyy_xxxxz, g_yz_0_xyyy_xxxy, g_yz_0_xyyy_xxxyy, g_yz_0_xyyy_xxxyz, g_yz_0_xyyy_xxxz, g_yz_0_xyyy_xxxzz, g_yz_0_xyyy_xxyy, g_yz_0_xyyy_xxyyy, g_yz_0_xyyy_xxyyz, g_yz_0_xyyy_xxyz, g_yz_0_xyyy_xxyzz, g_yz_0_xyyy_xxzz, g_yz_0_xyyy_xxzzz, g_yz_0_xyyy_xyyy, g_yz_0_xyyy_xyyyy, g_yz_0_xyyy_xyyyz, g_yz_0_xyyy_xyyz, g_yz_0_xyyy_xyyzz, g_yz_0_xyyy_xyzz, g_yz_0_xyyy_xyzzz, g_yz_0_xyyy_xzzz, g_yz_0_xyyy_xzzzz, g_yz_0_xyyy_yyyy, g_yz_0_xyyy_yyyz, g_yz_0_xyyy_yyzz, g_yz_0_xyyy_yzzz, g_yz_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyy_xxxx[k] = -g_yz_0_xyyy_xxxx[k] * ab_x + g_yz_0_xyyy_xxxxx[k];

                g_yz_0_xxyyy_xxxy[k] = -g_yz_0_xyyy_xxxy[k] * ab_x + g_yz_0_xyyy_xxxxy[k];

                g_yz_0_xxyyy_xxxz[k] = -g_yz_0_xyyy_xxxz[k] * ab_x + g_yz_0_xyyy_xxxxz[k];

                g_yz_0_xxyyy_xxyy[k] = -g_yz_0_xyyy_xxyy[k] * ab_x + g_yz_0_xyyy_xxxyy[k];

                g_yz_0_xxyyy_xxyz[k] = -g_yz_0_xyyy_xxyz[k] * ab_x + g_yz_0_xyyy_xxxyz[k];

                g_yz_0_xxyyy_xxzz[k] = -g_yz_0_xyyy_xxzz[k] * ab_x + g_yz_0_xyyy_xxxzz[k];

                g_yz_0_xxyyy_xyyy[k] = -g_yz_0_xyyy_xyyy[k] * ab_x + g_yz_0_xyyy_xxyyy[k];

                g_yz_0_xxyyy_xyyz[k] = -g_yz_0_xyyy_xyyz[k] * ab_x + g_yz_0_xyyy_xxyyz[k];

                g_yz_0_xxyyy_xyzz[k] = -g_yz_0_xyyy_xyzz[k] * ab_x + g_yz_0_xyyy_xxyzz[k];

                g_yz_0_xxyyy_xzzz[k] = -g_yz_0_xyyy_xzzz[k] * ab_x + g_yz_0_xyyy_xxzzz[k];

                g_yz_0_xxyyy_yyyy[k] = -g_yz_0_xyyy_yyyy[k] * ab_x + g_yz_0_xyyy_xyyyy[k];

                g_yz_0_xxyyy_yyyz[k] = -g_yz_0_xyyy_yyyz[k] * ab_x + g_yz_0_xyyy_xyyyz[k];

                g_yz_0_xxyyy_yyzz[k] = -g_yz_0_xyyy_yyzz[k] * ab_x + g_yz_0_xyyy_xyyzz[k];

                g_yz_0_xxyyy_yzzz[k] = -g_yz_0_xyyy_yzzz[k] * ab_x + g_yz_0_xyyy_xyzzz[k];

                g_yz_0_xxyyy_zzzz[k] = -g_yz_0_xyyy_zzzz[k] * ab_x + g_yz_0_xyyy_xzzzz[k];
            }

            /// Set up 1365-1380 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 1365 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 1366 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 1367 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 1368 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 1369 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 1370 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 1371 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 1372 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 1373 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 1374 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 1375 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 1376 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 1377 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 1378 * ccomps * dcomps);

            auto g_yz_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 1379 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyz_xxxx, g_yz_0_xxyyz_xxxy, g_yz_0_xxyyz_xxxz, g_yz_0_xxyyz_xxyy, g_yz_0_xxyyz_xxyz, g_yz_0_xxyyz_xxzz, g_yz_0_xxyyz_xyyy, g_yz_0_xxyyz_xyyz, g_yz_0_xxyyz_xyzz, g_yz_0_xxyyz_xzzz, g_yz_0_xxyyz_yyyy, g_yz_0_xxyyz_yyyz, g_yz_0_xxyyz_yyzz, g_yz_0_xxyyz_yzzz, g_yz_0_xxyyz_zzzz, g_yz_0_xyyz_xxxx, g_yz_0_xyyz_xxxxx, g_yz_0_xyyz_xxxxy, g_yz_0_xyyz_xxxxz, g_yz_0_xyyz_xxxy, g_yz_0_xyyz_xxxyy, g_yz_0_xyyz_xxxyz, g_yz_0_xyyz_xxxz, g_yz_0_xyyz_xxxzz, g_yz_0_xyyz_xxyy, g_yz_0_xyyz_xxyyy, g_yz_0_xyyz_xxyyz, g_yz_0_xyyz_xxyz, g_yz_0_xyyz_xxyzz, g_yz_0_xyyz_xxzz, g_yz_0_xyyz_xxzzz, g_yz_0_xyyz_xyyy, g_yz_0_xyyz_xyyyy, g_yz_0_xyyz_xyyyz, g_yz_0_xyyz_xyyz, g_yz_0_xyyz_xyyzz, g_yz_0_xyyz_xyzz, g_yz_0_xyyz_xyzzz, g_yz_0_xyyz_xzzz, g_yz_0_xyyz_xzzzz, g_yz_0_xyyz_yyyy, g_yz_0_xyyz_yyyz, g_yz_0_xyyz_yyzz, g_yz_0_xyyz_yzzz, g_yz_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyz_xxxx[k] = -g_yz_0_xyyz_xxxx[k] * ab_x + g_yz_0_xyyz_xxxxx[k];

                g_yz_0_xxyyz_xxxy[k] = -g_yz_0_xyyz_xxxy[k] * ab_x + g_yz_0_xyyz_xxxxy[k];

                g_yz_0_xxyyz_xxxz[k] = -g_yz_0_xyyz_xxxz[k] * ab_x + g_yz_0_xyyz_xxxxz[k];

                g_yz_0_xxyyz_xxyy[k] = -g_yz_0_xyyz_xxyy[k] * ab_x + g_yz_0_xyyz_xxxyy[k];

                g_yz_0_xxyyz_xxyz[k] = -g_yz_0_xyyz_xxyz[k] * ab_x + g_yz_0_xyyz_xxxyz[k];

                g_yz_0_xxyyz_xxzz[k] = -g_yz_0_xyyz_xxzz[k] * ab_x + g_yz_0_xyyz_xxxzz[k];

                g_yz_0_xxyyz_xyyy[k] = -g_yz_0_xyyz_xyyy[k] * ab_x + g_yz_0_xyyz_xxyyy[k];

                g_yz_0_xxyyz_xyyz[k] = -g_yz_0_xyyz_xyyz[k] * ab_x + g_yz_0_xyyz_xxyyz[k];

                g_yz_0_xxyyz_xyzz[k] = -g_yz_0_xyyz_xyzz[k] * ab_x + g_yz_0_xyyz_xxyzz[k];

                g_yz_0_xxyyz_xzzz[k] = -g_yz_0_xyyz_xzzz[k] * ab_x + g_yz_0_xyyz_xxzzz[k];

                g_yz_0_xxyyz_yyyy[k] = -g_yz_0_xyyz_yyyy[k] * ab_x + g_yz_0_xyyz_xyyyy[k];

                g_yz_0_xxyyz_yyyz[k] = -g_yz_0_xyyz_yyyz[k] * ab_x + g_yz_0_xyyz_xyyyz[k];

                g_yz_0_xxyyz_yyzz[k] = -g_yz_0_xyyz_yyzz[k] * ab_x + g_yz_0_xyyz_xyyzz[k];

                g_yz_0_xxyyz_yzzz[k] = -g_yz_0_xyyz_yzzz[k] * ab_x + g_yz_0_xyyz_xyzzz[k];

                g_yz_0_xxyyz_zzzz[k] = -g_yz_0_xyyz_zzzz[k] * ab_x + g_yz_0_xyyz_xzzzz[k];
            }

            /// Set up 1380-1395 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 1380 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 1381 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 1382 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 1383 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 1384 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 1385 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 1386 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 1387 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 1388 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 1389 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 1390 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 1391 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 1392 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 1393 * ccomps * dcomps);

            auto g_yz_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 1394 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzz_xxxx, g_yz_0_xxyzz_xxxy, g_yz_0_xxyzz_xxxz, g_yz_0_xxyzz_xxyy, g_yz_0_xxyzz_xxyz, g_yz_0_xxyzz_xxzz, g_yz_0_xxyzz_xyyy, g_yz_0_xxyzz_xyyz, g_yz_0_xxyzz_xyzz, g_yz_0_xxyzz_xzzz, g_yz_0_xxyzz_yyyy, g_yz_0_xxyzz_yyyz, g_yz_0_xxyzz_yyzz, g_yz_0_xxyzz_yzzz, g_yz_0_xxyzz_zzzz, g_yz_0_xyzz_xxxx, g_yz_0_xyzz_xxxxx, g_yz_0_xyzz_xxxxy, g_yz_0_xyzz_xxxxz, g_yz_0_xyzz_xxxy, g_yz_0_xyzz_xxxyy, g_yz_0_xyzz_xxxyz, g_yz_0_xyzz_xxxz, g_yz_0_xyzz_xxxzz, g_yz_0_xyzz_xxyy, g_yz_0_xyzz_xxyyy, g_yz_0_xyzz_xxyyz, g_yz_0_xyzz_xxyz, g_yz_0_xyzz_xxyzz, g_yz_0_xyzz_xxzz, g_yz_0_xyzz_xxzzz, g_yz_0_xyzz_xyyy, g_yz_0_xyzz_xyyyy, g_yz_0_xyzz_xyyyz, g_yz_0_xyzz_xyyz, g_yz_0_xyzz_xyyzz, g_yz_0_xyzz_xyzz, g_yz_0_xyzz_xyzzz, g_yz_0_xyzz_xzzz, g_yz_0_xyzz_xzzzz, g_yz_0_xyzz_yyyy, g_yz_0_xyzz_yyyz, g_yz_0_xyzz_yyzz, g_yz_0_xyzz_yzzz, g_yz_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzz_xxxx[k] = -g_yz_0_xyzz_xxxx[k] * ab_x + g_yz_0_xyzz_xxxxx[k];

                g_yz_0_xxyzz_xxxy[k] = -g_yz_0_xyzz_xxxy[k] * ab_x + g_yz_0_xyzz_xxxxy[k];

                g_yz_0_xxyzz_xxxz[k] = -g_yz_0_xyzz_xxxz[k] * ab_x + g_yz_0_xyzz_xxxxz[k];

                g_yz_0_xxyzz_xxyy[k] = -g_yz_0_xyzz_xxyy[k] * ab_x + g_yz_0_xyzz_xxxyy[k];

                g_yz_0_xxyzz_xxyz[k] = -g_yz_0_xyzz_xxyz[k] * ab_x + g_yz_0_xyzz_xxxyz[k];

                g_yz_0_xxyzz_xxzz[k] = -g_yz_0_xyzz_xxzz[k] * ab_x + g_yz_0_xyzz_xxxzz[k];

                g_yz_0_xxyzz_xyyy[k] = -g_yz_0_xyzz_xyyy[k] * ab_x + g_yz_0_xyzz_xxyyy[k];

                g_yz_0_xxyzz_xyyz[k] = -g_yz_0_xyzz_xyyz[k] * ab_x + g_yz_0_xyzz_xxyyz[k];

                g_yz_0_xxyzz_xyzz[k] = -g_yz_0_xyzz_xyzz[k] * ab_x + g_yz_0_xyzz_xxyzz[k];

                g_yz_0_xxyzz_xzzz[k] = -g_yz_0_xyzz_xzzz[k] * ab_x + g_yz_0_xyzz_xxzzz[k];

                g_yz_0_xxyzz_yyyy[k] = -g_yz_0_xyzz_yyyy[k] * ab_x + g_yz_0_xyzz_xyyyy[k];

                g_yz_0_xxyzz_yyyz[k] = -g_yz_0_xyzz_yyyz[k] * ab_x + g_yz_0_xyzz_xyyyz[k];

                g_yz_0_xxyzz_yyzz[k] = -g_yz_0_xyzz_yyzz[k] * ab_x + g_yz_0_xyzz_xyyzz[k];

                g_yz_0_xxyzz_yzzz[k] = -g_yz_0_xyzz_yzzz[k] * ab_x + g_yz_0_xyzz_xyzzz[k];

                g_yz_0_xxyzz_zzzz[k] = -g_yz_0_xyzz_zzzz[k] * ab_x + g_yz_0_xyzz_xzzzz[k];
            }

            /// Set up 1395-1410 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 1395 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 1396 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 1397 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 1398 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 1399 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 1400 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 1401 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 1402 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 1403 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 1404 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 1405 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 1406 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 1407 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 1408 * ccomps * dcomps);

            auto g_yz_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 1409 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzz_xxxx, g_yz_0_xxzzz_xxxy, g_yz_0_xxzzz_xxxz, g_yz_0_xxzzz_xxyy, g_yz_0_xxzzz_xxyz, g_yz_0_xxzzz_xxzz, g_yz_0_xxzzz_xyyy, g_yz_0_xxzzz_xyyz, g_yz_0_xxzzz_xyzz, g_yz_0_xxzzz_xzzz, g_yz_0_xxzzz_yyyy, g_yz_0_xxzzz_yyyz, g_yz_0_xxzzz_yyzz, g_yz_0_xxzzz_yzzz, g_yz_0_xxzzz_zzzz, g_yz_0_xzzz_xxxx, g_yz_0_xzzz_xxxxx, g_yz_0_xzzz_xxxxy, g_yz_0_xzzz_xxxxz, g_yz_0_xzzz_xxxy, g_yz_0_xzzz_xxxyy, g_yz_0_xzzz_xxxyz, g_yz_0_xzzz_xxxz, g_yz_0_xzzz_xxxzz, g_yz_0_xzzz_xxyy, g_yz_0_xzzz_xxyyy, g_yz_0_xzzz_xxyyz, g_yz_0_xzzz_xxyz, g_yz_0_xzzz_xxyzz, g_yz_0_xzzz_xxzz, g_yz_0_xzzz_xxzzz, g_yz_0_xzzz_xyyy, g_yz_0_xzzz_xyyyy, g_yz_0_xzzz_xyyyz, g_yz_0_xzzz_xyyz, g_yz_0_xzzz_xyyzz, g_yz_0_xzzz_xyzz, g_yz_0_xzzz_xyzzz, g_yz_0_xzzz_xzzz, g_yz_0_xzzz_xzzzz, g_yz_0_xzzz_yyyy, g_yz_0_xzzz_yyyz, g_yz_0_xzzz_yyzz, g_yz_0_xzzz_yzzz, g_yz_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzz_xxxx[k] = -g_yz_0_xzzz_xxxx[k] * ab_x + g_yz_0_xzzz_xxxxx[k];

                g_yz_0_xxzzz_xxxy[k] = -g_yz_0_xzzz_xxxy[k] * ab_x + g_yz_0_xzzz_xxxxy[k];

                g_yz_0_xxzzz_xxxz[k] = -g_yz_0_xzzz_xxxz[k] * ab_x + g_yz_0_xzzz_xxxxz[k];

                g_yz_0_xxzzz_xxyy[k] = -g_yz_0_xzzz_xxyy[k] * ab_x + g_yz_0_xzzz_xxxyy[k];

                g_yz_0_xxzzz_xxyz[k] = -g_yz_0_xzzz_xxyz[k] * ab_x + g_yz_0_xzzz_xxxyz[k];

                g_yz_0_xxzzz_xxzz[k] = -g_yz_0_xzzz_xxzz[k] * ab_x + g_yz_0_xzzz_xxxzz[k];

                g_yz_0_xxzzz_xyyy[k] = -g_yz_0_xzzz_xyyy[k] * ab_x + g_yz_0_xzzz_xxyyy[k];

                g_yz_0_xxzzz_xyyz[k] = -g_yz_0_xzzz_xyyz[k] * ab_x + g_yz_0_xzzz_xxyyz[k];

                g_yz_0_xxzzz_xyzz[k] = -g_yz_0_xzzz_xyzz[k] * ab_x + g_yz_0_xzzz_xxyzz[k];

                g_yz_0_xxzzz_xzzz[k] = -g_yz_0_xzzz_xzzz[k] * ab_x + g_yz_0_xzzz_xxzzz[k];

                g_yz_0_xxzzz_yyyy[k] = -g_yz_0_xzzz_yyyy[k] * ab_x + g_yz_0_xzzz_xyyyy[k];

                g_yz_0_xxzzz_yyyz[k] = -g_yz_0_xzzz_yyyz[k] * ab_x + g_yz_0_xzzz_xyyyz[k];

                g_yz_0_xxzzz_yyzz[k] = -g_yz_0_xzzz_yyzz[k] * ab_x + g_yz_0_xzzz_xyyzz[k];

                g_yz_0_xxzzz_yzzz[k] = -g_yz_0_xzzz_yzzz[k] * ab_x + g_yz_0_xzzz_xyzzz[k];

                g_yz_0_xxzzz_zzzz[k] = -g_yz_0_xzzz_zzzz[k] * ab_x + g_yz_0_xzzz_xzzzz[k];
            }

            /// Set up 1410-1425 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1410 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1411 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1412 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1413 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1414 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1415 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1416 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1417 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1418 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1419 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1420 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1421 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1422 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1423 * ccomps * dcomps);

            auto g_yz_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1424 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyy_xxxx, g_yz_0_xyyyy_xxxy, g_yz_0_xyyyy_xxxz, g_yz_0_xyyyy_xxyy, g_yz_0_xyyyy_xxyz, g_yz_0_xyyyy_xxzz, g_yz_0_xyyyy_xyyy, g_yz_0_xyyyy_xyyz, g_yz_0_xyyyy_xyzz, g_yz_0_xyyyy_xzzz, g_yz_0_xyyyy_yyyy, g_yz_0_xyyyy_yyyz, g_yz_0_xyyyy_yyzz, g_yz_0_xyyyy_yzzz, g_yz_0_xyyyy_zzzz, g_yz_0_yyyy_xxxx, g_yz_0_yyyy_xxxxx, g_yz_0_yyyy_xxxxy, g_yz_0_yyyy_xxxxz, g_yz_0_yyyy_xxxy, g_yz_0_yyyy_xxxyy, g_yz_0_yyyy_xxxyz, g_yz_0_yyyy_xxxz, g_yz_0_yyyy_xxxzz, g_yz_0_yyyy_xxyy, g_yz_0_yyyy_xxyyy, g_yz_0_yyyy_xxyyz, g_yz_0_yyyy_xxyz, g_yz_0_yyyy_xxyzz, g_yz_0_yyyy_xxzz, g_yz_0_yyyy_xxzzz, g_yz_0_yyyy_xyyy, g_yz_0_yyyy_xyyyy, g_yz_0_yyyy_xyyyz, g_yz_0_yyyy_xyyz, g_yz_0_yyyy_xyyzz, g_yz_0_yyyy_xyzz, g_yz_0_yyyy_xyzzz, g_yz_0_yyyy_xzzz, g_yz_0_yyyy_xzzzz, g_yz_0_yyyy_yyyy, g_yz_0_yyyy_yyyz, g_yz_0_yyyy_yyzz, g_yz_0_yyyy_yzzz, g_yz_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyy_xxxx[k] = -g_yz_0_yyyy_xxxx[k] * ab_x + g_yz_0_yyyy_xxxxx[k];

                g_yz_0_xyyyy_xxxy[k] = -g_yz_0_yyyy_xxxy[k] * ab_x + g_yz_0_yyyy_xxxxy[k];

                g_yz_0_xyyyy_xxxz[k] = -g_yz_0_yyyy_xxxz[k] * ab_x + g_yz_0_yyyy_xxxxz[k];

                g_yz_0_xyyyy_xxyy[k] = -g_yz_0_yyyy_xxyy[k] * ab_x + g_yz_0_yyyy_xxxyy[k];

                g_yz_0_xyyyy_xxyz[k] = -g_yz_0_yyyy_xxyz[k] * ab_x + g_yz_0_yyyy_xxxyz[k];

                g_yz_0_xyyyy_xxzz[k] = -g_yz_0_yyyy_xxzz[k] * ab_x + g_yz_0_yyyy_xxxzz[k];

                g_yz_0_xyyyy_xyyy[k] = -g_yz_0_yyyy_xyyy[k] * ab_x + g_yz_0_yyyy_xxyyy[k];

                g_yz_0_xyyyy_xyyz[k] = -g_yz_0_yyyy_xyyz[k] * ab_x + g_yz_0_yyyy_xxyyz[k];

                g_yz_0_xyyyy_xyzz[k] = -g_yz_0_yyyy_xyzz[k] * ab_x + g_yz_0_yyyy_xxyzz[k];

                g_yz_0_xyyyy_xzzz[k] = -g_yz_0_yyyy_xzzz[k] * ab_x + g_yz_0_yyyy_xxzzz[k];

                g_yz_0_xyyyy_yyyy[k] = -g_yz_0_yyyy_yyyy[k] * ab_x + g_yz_0_yyyy_xyyyy[k];

                g_yz_0_xyyyy_yyyz[k] = -g_yz_0_yyyy_yyyz[k] * ab_x + g_yz_0_yyyy_xyyyz[k];

                g_yz_0_xyyyy_yyzz[k] = -g_yz_0_yyyy_yyzz[k] * ab_x + g_yz_0_yyyy_xyyzz[k];

                g_yz_0_xyyyy_yzzz[k] = -g_yz_0_yyyy_yzzz[k] * ab_x + g_yz_0_yyyy_xyzzz[k];

                g_yz_0_xyyyy_zzzz[k] = -g_yz_0_yyyy_zzzz[k] * ab_x + g_yz_0_yyyy_xzzzz[k];
            }

            /// Set up 1425-1440 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1425 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1426 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1427 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1428 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1429 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1430 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1431 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1432 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1433 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1434 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1435 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1436 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1437 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1438 * ccomps * dcomps);

            auto g_yz_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1439 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyz_xxxx, g_yz_0_xyyyz_xxxy, g_yz_0_xyyyz_xxxz, g_yz_0_xyyyz_xxyy, g_yz_0_xyyyz_xxyz, g_yz_0_xyyyz_xxzz, g_yz_0_xyyyz_xyyy, g_yz_0_xyyyz_xyyz, g_yz_0_xyyyz_xyzz, g_yz_0_xyyyz_xzzz, g_yz_0_xyyyz_yyyy, g_yz_0_xyyyz_yyyz, g_yz_0_xyyyz_yyzz, g_yz_0_xyyyz_yzzz, g_yz_0_xyyyz_zzzz, g_yz_0_yyyz_xxxx, g_yz_0_yyyz_xxxxx, g_yz_0_yyyz_xxxxy, g_yz_0_yyyz_xxxxz, g_yz_0_yyyz_xxxy, g_yz_0_yyyz_xxxyy, g_yz_0_yyyz_xxxyz, g_yz_0_yyyz_xxxz, g_yz_0_yyyz_xxxzz, g_yz_0_yyyz_xxyy, g_yz_0_yyyz_xxyyy, g_yz_0_yyyz_xxyyz, g_yz_0_yyyz_xxyz, g_yz_0_yyyz_xxyzz, g_yz_0_yyyz_xxzz, g_yz_0_yyyz_xxzzz, g_yz_0_yyyz_xyyy, g_yz_0_yyyz_xyyyy, g_yz_0_yyyz_xyyyz, g_yz_0_yyyz_xyyz, g_yz_0_yyyz_xyyzz, g_yz_0_yyyz_xyzz, g_yz_0_yyyz_xyzzz, g_yz_0_yyyz_xzzz, g_yz_0_yyyz_xzzzz, g_yz_0_yyyz_yyyy, g_yz_0_yyyz_yyyz, g_yz_0_yyyz_yyzz, g_yz_0_yyyz_yzzz, g_yz_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyz_xxxx[k] = -g_yz_0_yyyz_xxxx[k] * ab_x + g_yz_0_yyyz_xxxxx[k];

                g_yz_0_xyyyz_xxxy[k] = -g_yz_0_yyyz_xxxy[k] * ab_x + g_yz_0_yyyz_xxxxy[k];

                g_yz_0_xyyyz_xxxz[k] = -g_yz_0_yyyz_xxxz[k] * ab_x + g_yz_0_yyyz_xxxxz[k];

                g_yz_0_xyyyz_xxyy[k] = -g_yz_0_yyyz_xxyy[k] * ab_x + g_yz_0_yyyz_xxxyy[k];

                g_yz_0_xyyyz_xxyz[k] = -g_yz_0_yyyz_xxyz[k] * ab_x + g_yz_0_yyyz_xxxyz[k];

                g_yz_0_xyyyz_xxzz[k] = -g_yz_0_yyyz_xxzz[k] * ab_x + g_yz_0_yyyz_xxxzz[k];

                g_yz_0_xyyyz_xyyy[k] = -g_yz_0_yyyz_xyyy[k] * ab_x + g_yz_0_yyyz_xxyyy[k];

                g_yz_0_xyyyz_xyyz[k] = -g_yz_0_yyyz_xyyz[k] * ab_x + g_yz_0_yyyz_xxyyz[k];

                g_yz_0_xyyyz_xyzz[k] = -g_yz_0_yyyz_xyzz[k] * ab_x + g_yz_0_yyyz_xxyzz[k];

                g_yz_0_xyyyz_xzzz[k] = -g_yz_0_yyyz_xzzz[k] * ab_x + g_yz_0_yyyz_xxzzz[k];

                g_yz_0_xyyyz_yyyy[k] = -g_yz_0_yyyz_yyyy[k] * ab_x + g_yz_0_yyyz_xyyyy[k];

                g_yz_0_xyyyz_yyyz[k] = -g_yz_0_yyyz_yyyz[k] * ab_x + g_yz_0_yyyz_xyyyz[k];

                g_yz_0_xyyyz_yyzz[k] = -g_yz_0_yyyz_yyzz[k] * ab_x + g_yz_0_yyyz_xyyzz[k];

                g_yz_0_xyyyz_yzzz[k] = -g_yz_0_yyyz_yzzz[k] * ab_x + g_yz_0_yyyz_xyzzz[k];

                g_yz_0_xyyyz_zzzz[k] = -g_yz_0_yyyz_zzzz[k] * ab_x + g_yz_0_yyyz_xzzzz[k];
            }

            /// Set up 1440-1455 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1440 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1441 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1442 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1443 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1444 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1445 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1446 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1447 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1448 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1449 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1450 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1451 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1452 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1453 * ccomps * dcomps);

            auto g_yz_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1454 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzz_xxxx, g_yz_0_xyyzz_xxxy, g_yz_0_xyyzz_xxxz, g_yz_0_xyyzz_xxyy, g_yz_0_xyyzz_xxyz, g_yz_0_xyyzz_xxzz, g_yz_0_xyyzz_xyyy, g_yz_0_xyyzz_xyyz, g_yz_0_xyyzz_xyzz, g_yz_0_xyyzz_xzzz, g_yz_0_xyyzz_yyyy, g_yz_0_xyyzz_yyyz, g_yz_0_xyyzz_yyzz, g_yz_0_xyyzz_yzzz, g_yz_0_xyyzz_zzzz, g_yz_0_yyzz_xxxx, g_yz_0_yyzz_xxxxx, g_yz_0_yyzz_xxxxy, g_yz_0_yyzz_xxxxz, g_yz_0_yyzz_xxxy, g_yz_0_yyzz_xxxyy, g_yz_0_yyzz_xxxyz, g_yz_0_yyzz_xxxz, g_yz_0_yyzz_xxxzz, g_yz_0_yyzz_xxyy, g_yz_0_yyzz_xxyyy, g_yz_0_yyzz_xxyyz, g_yz_0_yyzz_xxyz, g_yz_0_yyzz_xxyzz, g_yz_0_yyzz_xxzz, g_yz_0_yyzz_xxzzz, g_yz_0_yyzz_xyyy, g_yz_0_yyzz_xyyyy, g_yz_0_yyzz_xyyyz, g_yz_0_yyzz_xyyz, g_yz_0_yyzz_xyyzz, g_yz_0_yyzz_xyzz, g_yz_0_yyzz_xyzzz, g_yz_0_yyzz_xzzz, g_yz_0_yyzz_xzzzz, g_yz_0_yyzz_yyyy, g_yz_0_yyzz_yyyz, g_yz_0_yyzz_yyzz, g_yz_0_yyzz_yzzz, g_yz_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzz_xxxx[k] = -g_yz_0_yyzz_xxxx[k] * ab_x + g_yz_0_yyzz_xxxxx[k];

                g_yz_0_xyyzz_xxxy[k] = -g_yz_0_yyzz_xxxy[k] * ab_x + g_yz_0_yyzz_xxxxy[k];

                g_yz_0_xyyzz_xxxz[k] = -g_yz_0_yyzz_xxxz[k] * ab_x + g_yz_0_yyzz_xxxxz[k];

                g_yz_0_xyyzz_xxyy[k] = -g_yz_0_yyzz_xxyy[k] * ab_x + g_yz_0_yyzz_xxxyy[k];

                g_yz_0_xyyzz_xxyz[k] = -g_yz_0_yyzz_xxyz[k] * ab_x + g_yz_0_yyzz_xxxyz[k];

                g_yz_0_xyyzz_xxzz[k] = -g_yz_0_yyzz_xxzz[k] * ab_x + g_yz_0_yyzz_xxxzz[k];

                g_yz_0_xyyzz_xyyy[k] = -g_yz_0_yyzz_xyyy[k] * ab_x + g_yz_0_yyzz_xxyyy[k];

                g_yz_0_xyyzz_xyyz[k] = -g_yz_0_yyzz_xyyz[k] * ab_x + g_yz_0_yyzz_xxyyz[k];

                g_yz_0_xyyzz_xyzz[k] = -g_yz_0_yyzz_xyzz[k] * ab_x + g_yz_0_yyzz_xxyzz[k];

                g_yz_0_xyyzz_xzzz[k] = -g_yz_0_yyzz_xzzz[k] * ab_x + g_yz_0_yyzz_xxzzz[k];

                g_yz_0_xyyzz_yyyy[k] = -g_yz_0_yyzz_yyyy[k] * ab_x + g_yz_0_yyzz_xyyyy[k];

                g_yz_0_xyyzz_yyyz[k] = -g_yz_0_yyzz_yyyz[k] * ab_x + g_yz_0_yyzz_xyyyz[k];

                g_yz_0_xyyzz_yyzz[k] = -g_yz_0_yyzz_yyzz[k] * ab_x + g_yz_0_yyzz_xyyzz[k];

                g_yz_0_xyyzz_yzzz[k] = -g_yz_0_yyzz_yzzz[k] * ab_x + g_yz_0_yyzz_xyzzz[k];

                g_yz_0_xyyzz_zzzz[k] = -g_yz_0_yyzz_zzzz[k] * ab_x + g_yz_0_yyzz_xzzzz[k];
            }

            /// Set up 1455-1470 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1455 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1456 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1457 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1458 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1459 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1460 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1461 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1462 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1463 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1464 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1465 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1466 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1467 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1468 * ccomps * dcomps);

            auto g_yz_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzz_xxxx, g_yz_0_xyzzz_xxxy, g_yz_0_xyzzz_xxxz, g_yz_0_xyzzz_xxyy, g_yz_0_xyzzz_xxyz, g_yz_0_xyzzz_xxzz, g_yz_0_xyzzz_xyyy, g_yz_0_xyzzz_xyyz, g_yz_0_xyzzz_xyzz, g_yz_0_xyzzz_xzzz, g_yz_0_xyzzz_yyyy, g_yz_0_xyzzz_yyyz, g_yz_0_xyzzz_yyzz, g_yz_0_xyzzz_yzzz, g_yz_0_xyzzz_zzzz, g_yz_0_yzzz_xxxx, g_yz_0_yzzz_xxxxx, g_yz_0_yzzz_xxxxy, g_yz_0_yzzz_xxxxz, g_yz_0_yzzz_xxxy, g_yz_0_yzzz_xxxyy, g_yz_0_yzzz_xxxyz, g_yz_0_yzzz_xxxz, g_yz_0_yzzz_xxxzz, g_yz_0_yzzz_xxyy, g_yz_0_yzzz_xxyyy, g_yz_0_yzzz_xxyyz, g_yz_0_yzzz_xxyz, g_yz_0_yzzz_xxyzz, g_yz_0_yzzz_xxzz, g_yz_0_yzzz_xxzzz, g_yz_0_yzzz_xyyy, g_yz_0_yzzz_xyyyy, g_yz_0_yzzz_xyyyz, g_yz_0_yzzz_xyyz, g_yz_0_yzzz_xyyzz, g_yz_0_yzzz_xyzz, g_yz_0_yzzz_xyzzz, g_yz_0_yzzz_xzzz, g_yz_0_yzzz_xzzzz, g_yz_0_yzzz_yyyy, g_yz_0_yzzz_yyyz, g_yz_0_yzzz_yyzz, g_yz_0_yzzz_yzzz, g_yz_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzz_xxxx[k] = -g_yz_0_yzzz_xxxx[k] * ab_x + g_yz_0_yzzz_xxxxx[k];

                g_yz_0_xyzzz_xxxy[k] = -g_yz_0_yzzz_xxxy[k] * ab_x + g_yz_0_yzzz_xxxxy[k];

                g_yz_0_xyzzz_xxxz[k] = -g_yz_0_yzzz_xxxz[k] * ab_x + g_yz_0_yzzz_xxxxz[k];

                g_yz_0_xyzzz_xxyy[k] = -g_yz_0_yzzz_xxyy[k] * ab_x + g_yz_0_yzzz_xxxyy[k];

                g_yz_0_xyzzz_xxyz[k] = -g_yz_0_yzzz_xxyz[k] * ab_x + g_yz_0_yzzz_xxxyz[k];

                g_yz_0_xyzzz_xxzz[k] = -g_yz_0_yzzz_xxzz[k] * ab_x + g_yz_0_yzzz_xxxzz[k];

                g_yz_0_xyzzz_xyyy[k] = -g_yz_0_yzzz_xyyy[k] * ab_x + g_yz_0_yzzz_xxyyy[k];

                g_yz_0_xyzzz_xyyz[k] = -g_yz_0_yzzz_xyyz[k] * ab_x + g_yz_0_yzzz_xxyyz[k];

                g_yz_0_xyzzz_xyzz[k] = -g_yz_0_yzzz_xyzz[k] * ab_x + g_yz_0_yzzz_xxyzz[k];

                g_yz_0_xyzzz_xzzz[k] = -g_yz_0_yzzz_xzzz[k] * ab_x + g_yz_0_yzzz_xxzzz[k];

                g_yz_0_xyzzz_yyyy[k] = -g_yz_0_yzzz_yyyy[k] * ab_x + g_yz_0_yzzz_xyyyy[k];

                g_yz_0_xyzzz_yyyz[k] = -g_yz_0_yzzz_yyyz[k] * ab_x + g_yz_0_yzzz_xyyyz[k];

                g_yz_0_xyzzz_yyzz[k] = -g_yz_0_yzzz_yyzz[k] * ab_x + g_yz_0_yzzz_xyyzz[k];

                g_yz_0_xyzzz_yzzz[k] = -g_yz_0_yzzz_yzzz[k] * ab_x + g_yz_0_yzzz_xyzzz[k];

                g_yz_0_xyzzz_zzzz[k] = -g_yz_0_yzzz_zzzz[k] * ab_x + g_yz_0_yzzz_xzzzz[k];
            }

            /// Set up 1470-1485 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1470 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1471 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1472 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1473 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1474 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1475 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1476 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1477 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1478 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1479 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1480 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1481 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1482 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1483 * ccomps * dcomps);

            auto g_yz_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1484 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzz_xxxx, g_yz_0_xzzzz_xxxy, g_yz_0_xzzzz_xxxz, g_yz_0_xzzzz_xxyy, g_yz_0_xzzzz_xxyz, g_yz_0_xzzzz_xxzz, g_yz_0_xzzzz_xyyy, g_yz_0_xzzzz_xyyz, g_yz_0_xzzzz_xyzz, g_yz_0_xzzzz_xzzz, g_yz_0_xzzzz_yyyy, g_yz_0_xzzzz_yyyz, g_yz_0_xzzzz_yyzz, g_yz_0_xzzzz_yzzz, g_yz_0_xzzzz_zzzz, g_yz_0_zzzz_xxxx, g_yz_0_zzzz_xxxxx, g_yz_0_zzzz_xxxxy, g_yz_0_zzzz_xxxxz, g_yz_0_zzzz_xxxy, g_yz_0_zzzz_xxxyy, g_yz_0_zzzz_xxxyz, g_yz_0_zzzz_xxxz, g_yz_0_zzzz_xxxzz, g_yz_0_zzzz_xxyy, g_yz_0_zzzz_xxyyy, g_yz_0_zzzz_xxyyz, g_yz_0_zzzz_xxyz, g_yz_0_zzzz_xxyzz, g_yz_0_zzzz_xxzz, g_yz_0_zzzz_xxzzz, g_yz_0_zzzz_xyyy, g_yz_0_zzzz_xyyyy, g_yz_0_zzzz_xyyyz, g_yz_0_zzzz_xyyz, g_yz_0_zzzz_xyyzz, g_yz_0_zzzz_xyzz, g_yz_0_zzzz_xyzzz, g_yz_0_zzzz_xzzz, g_yz_0_zzzz_xzzzz, g_yz_0_zzzz_yyyy, g_yz_0_zzzz_yyyz, g_yz_0_zzzz_yyzz, g_yz_0_zzzz_yzzz, g_yz_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzz_xxxx[k] = -g_yz_0_zzzz_xxxx[k] * ab_x + g_yz_0_zzzz_xxxxx[k];

                g_yz_0_xzzzz_xxxy[k] = -g_yz_0_zzzz_xxxy[k] * ab_x + g_yz_0_zzzz_xxxxy[k];

                g_yz_0_xzzzz_xxxz[k] = -g_yz_0_zzzz_xxxz[k] * ab_x + g_yz_0_zzzz_xxxxz[k];

                g_yz_0_xzzzz_xxyy[k] = -g_yz_0_zzzz_xxyy[k] * ab_x + g_yz_0_zzzz_xxxyy[k];

                g_yz_0_xzzzz_xxyz[k] = -g_yz_0_zzzz_xxyz[k] * ab_x + g_yz_0_zzzz_xxxyz[k];

                g_yz_0_xzzzz_xxzz[k] = -g_yz_0_zzzz_xxzz[k] * ab_x + g_yz_0_zzzz_xxxzz[k];

                g_yz_0_xzzzz_xyyy[k] = -g_yz_0_zzzz_xyyy[k] * ab_x + g_yz_0_zzzz_xxyyy[k];

                g_yz_0_xzzzz_xyyz[k] = -g_yz_0_zzzz_xyyz[k] * ab_x + g_yz_0_zzzz_xxyyz[k];

                g_yz_0_xzzzz_xyzz[k] = -g_yz_0_zzzz_xyzz[k] * ab_x + g_yz_0_zzzz_xxyzz[k];

                g_yz_0_xzzzz_xzzz[k] = -g_yz_0_zzzz_xzzz[k] * ab_x + g_yz_0_zzzz_xxzzz[k];

                g_yz_0_xzzzz_yyyy[k] = -g_yz_0_zzzz_yyyy[k] * ab_x + g_yz_0_zzzz_xyyyy[k];

                g_yz_0_xzzzz_yyyz[k] = -g_yz_0_zzzz_yyyz[k] * ab_x + g_yz_0_zzzz_xyyyz[k];

                g_yz_0_xzzzz_yyzz[k] = -g_yz_0_zzzz_yyzz[k] * ab_x + g_yz_0_zzzz_xyyzz[k];

                g_yz_0_xzzzz_yzzz[k] = -g_yz_0_zzzz_yzzz[k] * ab_x + g_yz_0_zzzz_xyzzz[k];

                g_yz_0_xzzzz_zzzz[k] = -g_yz_0_zzzz_zzzz[k] * ab_x + g_yz_0_zzzz_xzzzz[k];
            }

            /// Set up 1485-1500 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1485 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1486 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1487 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1488 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1489 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1490 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1491 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1492 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1493 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1494 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1495 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1496 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1497 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1498 * ccomps * dcomps);

            auto g_yz_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1499 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyy_xxxx, g_yz_0_yyyy_xxxxy, g_yz_0_yyyy_xxxy, g_yz_0_yyyy_xxxyy, g_yz_0_yyyy_xxxyz, g_yz_0_yyyy_xxxz, g_yz_0_yyyy_xxyy, g_yz_0_yyyy_xxyyy, g_yz_0_yyyy_xxyyz, g_yz_0_yyyy_xxyz, g_yz_0_yyyy_xxyzz, g_yz_0_yyyy_xxzz, g_yz_0_yyyy_xyyy, g_yz_0_yyyy_xyyyy, g_yz_0_yyyy_xyyyz, g_yz_0_yyyy_xyyz, g_yz_0_yyyy_xyyzz, g_yz_0_yyyy_xyzz, g_yz_0_yyyy_xyzzz, g_yz_0_yyyy_xzzz, g_yz_0_yyyy_yyyy, g_yz_0_yyyy_yyyyy, g_yz_0_yyyy_yyyyz, g_yz_0_yyyy_yyyz, g_yz_0_yyyy_yyyzz, g_yz_0_yyyy_yyzz, g_yz_0_yyyy_yyzzz, g_yz_0_yyyy_yzzz, g_yz_0_yyyy_yzzzz, g_yz_0_yyyy_zzzz, g_yz_0_yyyyy_xxxx, g_yz_0_yyyyy_xxxy, g_yz_0_yyyyy_xxxz, g_yz_0_yyyyy_xxyy, g_yz_0_yyyyy_xxyz, g_yz_0_yyyyy_xxzz, g_yz_0_yyyyy_xyyy, g_yz_0_yyyyy_xyyz, g_yz_0_yyyyy_xyzz, g_yz_0_yyyyy_xzzz, g_yz_0_yyyyy_yyyy, g_yz_0_yyyyy_yyyz, g_yz_0_yyyyy_yyzz, g_yz_0_yyyyy_yzzz, g_yz_0_yyyyy_zzzz, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyy_xxxx[k] = -g_z_0_yyyy_xxxx[k] - g_yz_0_yyyy_xxxx[k] * ab_y + g_yz_0_yyyy_xxxxy[k];

                g_yz_0_yyyyy_xxxy[k] = -g_z_0_yyyy_xxxy[k] - g_yz_0_yyyy_xxxy[k] * ab_y + g_yz_0_yyyy_xxxyy[k];

                g_yz_0_yyyyy_xxxz[k] = -g_z_0_yyyy_xxxz[k] - g_yz_0_yyyy_xxxz[k] * ab_y + g_yz_0_yyyy_xxxyz[k];

                g_yz_0_yyyyy_xxyy[k] = -g_z_0_yyyy_xxyy[k] - g_yz_0_yyyy_xxyy[k] * ab_y + g_yz_0_yyyy_xxyyy[k];

                g_yz_0_yyyyy_xxyz[k] = -g_z_0_yyyy_xxyz[k] - g_yz_0_yyyy_xxyz[k] * ab_y + g_yz_0_yyyy_xxyyz[k];

                g_yz_0_yyyyy_xxzz[k] = -g_z_0_yyyy_xxzz[k] - g_yz_0_yyyy_xxzz[k] * ab_y + g_yz_0_yyyy_xxyzz[k];

                g_yz_0_yyyyy_xyyy[k] = -g_z_0_yyyy_xyyy[k] - g_yz_0_yyyy_xyyy[k] * ab_y + g_yz_0_yyyy_xyyyy[k];

                g_yz_0_yyyyy_xyyz[k] = -g_z_0_yyyy_xyyz[k] - g_yz_0_yyyy_xyyz[k] * ab_y + g_yz_0_yyyy_xyyyz[k];

                g_yz_0_yyyyy_xyzz[k] = -g_z_0_yyyy_xyzz[k] - g_yz_0_yyyy_xyzz[k] * ab_y + g_yz_0_yyyy_xyyzz[k];

                g_yz_0_yyyyy_xzzz[k] = -g_z_0_yyyy_xzzz[k] - g_yz_0_yyyy_xzzz[k] * ab_y + g_yz_0_yyyy_xyzzz[k];

                g_yz_0_yyyyy_yyyy[k] = -g_z_0_yyyy_yyyy[k] - g_yz_0_yyyy_yyyy[k] * ab_y + g_yz_0_yyyy_yyyyy[k];

                g_yz_0_yyyyy_yyyz[k] = -g_z_0_yyyy_yyyz[k] - g_yz_0_yyyy_yyyz[k] * ab_y + g_yz_0_yyyy_yyyyz[k];

                g_yz_0_yyyyy_yyzz[k] = -g_z_0_yyyy_yyzz[k] - g_yz_0_yyyy_yyzz[k] * ab_y + g_yz_0_yyyy_yyyzz[k];

                g_yz_0_yyyyy_yzzz[k] = -g_z_0_yyyy_yzzz[k] - g_yz_0_yyyy_yzzz[k] * ab_y + g_yz_0_yyyy_yyzzz[k];

                g_yz_0_yyyyy_zzzz[k] = -g_z_0_yyyy_zzzz[k] - g_yz_0_yyyy_zzzz[k] * ab_y + g_yz_0_yyyy_yzzzz[k];
            }

            /// Set up 1500-1515 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1500 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1501 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1502 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1503 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1504 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1505 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1506 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1507 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1508 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1509 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1510 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1511 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1512 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1513 * ccomps * dcomps);

            auto g_yz_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1514 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyz_xxxx, g_yz_0_yyyyz_xxxy, g_yz_0_yyyyz_xxxz, g_yz_0_yyyyz_xxyy, g_yz_0_yyyyz_xxyz, g_yz_0_yyyyz_xxzz, g_yz_0_yyyyz_xyyy, g_yz_0_yyyyz_xyyz, g_yz_0_yyyyz_xyzz, g_yz_0_yyyyz_xzzz, g_yz_0_yyyyz_yyyy, g_yz_0_yyyyz_yyyz, g_yz_0_yyyyz_yyzz, g_yz_0_yyyyz_yzzz, g_yz_0_yyyyz_zzzz, g_yz_0_yyyz_xxxx, g_yz_0_yyyz_xxxxy, g_yz_0_yyyz_xxxy, g_yz_0_yyyz_xxxyy, g_yz_0_yyyz_xxxyz, g_yz_0_yyyz_xxxz, g_yz_0_yyyz_xxyy, g_yz_0_yyyz_xxyyy, g_yz_0_yyyz_xxyyz, g_yz_0_yyyz_xxyz, g_yz_0_yyyz_xxyzz, g_yz_0_yyyz_xxzz, g_yz_0_yyyz_xyyy, g_yz_0_yyyz_xyyyy, g_yz_0_yyyz_xyyyz, g_yz_0_yyyz_xyyz, g_yz_0_yyyz_xyyzz, g_yz_0_yyyz_xyzz, g_yz_0_yyyz_xyzzz, g_yz_0_yyyz_xzzz, g_yz_0_yyyz_yyyy, g_yz_0_yyyz_yyyyy, g_yz_0_yyyz_yyyyz, g_yz_0_yyyz_yyyz, g_yz_0_yyyz_yyyzz, g_yz_0_yyyz_yyzz, g_yz_0_yyyz_yyzzz, g_yz_0_yyyz_yzzz, g_yz_0_yyyz_yzzzz, g_yz_0_yyyz_zzzz, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyz_xxxx[k] = -g_z_0_yyyz_xxxx[k] - g_yz_0_yyyz_xxxx[k] * ab_y + g_yz_0_yyyz_xxxxy[k];

                g_yz_0_yyyyz_xxxy[k] = -g_z_0_yyyz_xxxy[k] - g_yz_0_yyyz_xxxy[k] * ab_y + g_yz_0_yyyz_xxxyy[k];

                g_yz_0_yyyyz_xxxz[k] = -g_z_0_yyyz_xxxz[k] - g_yz_0_yyyz_xxxz[k] * ab_y + g_yz_0_yyyz_xxxyz[k];

                g_yz_0_yyyyz_xxyy[k] = -g_z_0_yyyz_xxyy[k] - g_yz_0_yyyz_xxyy[k] * ab_y + g_yz_0_yyyz_xxyyy[k];

                g_yz_0_yyyyz_xxyz[k] = -g_z_0_yyyz_xxyz[k] - g_yz_0_yyyz_xxyz[k] * ab_y + g_yz_0_yyyz_xxyyz[k];

                g_yz_0_yyyyz_xxzz[k] = -g_z_0_yyyz_xxzz[k] - g_yz_0_yyyz_xxzz[k] * ab_y + g_yz_0_yyyz_xxyzz[k];

                g_yz_0_yyyyz_xyyy[k] = -g_z_0_yyyz_xyyy[k] - g_yz_0_yyyz_xyyy[k] * ab_y + g_yz_0_yyyz_xyyyy[k];

                g_yz_0_yyyyz_xyyz[k] = -g_z_0_yyyz_xyyz[k] - g_yz_0_yyyz_xyyz[k] * ab_y + g_yz_0_yyyz_xyyyz[k];

                g_yz_0_yyyyz_xyzz[k] = -g_z_0_yyyz_xyzz[k] - g_yz_0_yyyz_xyzz[k] * ab_y + g_yz_0_yyyz_xyyzz[k];

                g_yz_0_yyyyz_xzzz[k] = -g_z_0_yyyz_xzzz[k] - g_yz_0_yyyz_xzzz[k] * ab_y + g_yz_0_yyyz_xyzzz[k];

                g_yz_0_yyyyz_yyyy[k] = -g_z_0_yyyz_yyyy[k] - g_yz_0_yyyz_yyyy[k] * ab_y + g_yz_0_yyyz_yyyyy[k];

                g_yz_0_yyyyz_yyyz[k] = -g_z_0_yyyz_yyyz[k] - g_yz_0_yyyz_yyyz[k] * ab_y + g_yz_0_yyyz_yyyyz[k];

                g_yz_0_yyyyz_yyzz[k] = -g_z_0_yyyz_yyzz[k] - g_yz_0_yyyz_yyzz[k] * ab_y + g_yz_0_yyyz_yyyzz[k];

                g_yz_0_yyyyz_yzzz[k] = -g_z_0_yyyz_yzzz[k] - g_yz_0_yyyz_yzzz[k] * ab_y + g_yz_0_yyyz_yyzzz[k];

                g_yz_0_yyyyz_zzzz[k] = -g_z_0_yyyz_zzzz[k] - g_yz_0_yyyz_zzzz[k] * ab_y + g_yz_0_yyyz_yzzzz[k];
            }

            /// Set up 1515-1530 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1515 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1516 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1517 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1518 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1519 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1520 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1521 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1522 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1523 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1524 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1525 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1526 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1527 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1528 * ccomps * dcomps);

            auto g_yz_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1529 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzz_xxxx, g_yz_0_yyyzz_xxxy, g_yz_0_yyyzz_xxxz, g_yz_0_yyyzz_xxyy, g_yz_0_yyyzz_xxyz, g_yz_0_yyyzz_xxzz, g_yz_0_yyyzz_xyyy, g_yz_0_yyyzz_xyyz, g_yz_0_yyyzz_xyzz, g_yz_0_yyyzz_xzzz, g_yz_0_yyyzz_yyyy, g_yz_0_yyyzz_yyyz, g_yz_0_yyyzz_yyzz, g_yz_0_yyyzz_yzzz, g_yz_0_yyyzz_zzzz, g_yz_0_yyzz_xxxx, g_yz_0_yyzz_xxxxy, g_yz_0_yyzz_xxxy, g_yz_0_yyzz_xxxyy, g_yz_0_yyzz_xxxyz, g_yz_0_yyzz_xxxz, g_yz_0_yyzz_xxyy, g_yz_0_yyzz_xxyyy, g_yz_0_yyzz_xxyyz, g_yz_0_yyzz_xxyz, g_yz_0_yyzz_xxyzz, g_yz_0_yyzz_xxzz, g_yz_0_yyzz_xyyy, g_yz_0_yyzz_xyyyy, g_yz_0_yyzz_xyyyz, g_yz_0_yyzz_xyyz, g_yz_0_yyzz_xyyzz, g_yz_0_yyzz_xyzz, g_yz_0_yyzz_xyzzz, g_yz_0_yyzz_xzzz, g_yz_0_yyzz_yyyy, g_yz_0_yyzz_yyyyy, g_yz_0_yyzz_yyyyz, g_yz_0_yyzz_yyyz, g_yz_0_yyzz_yyyzz, g_yz_0_yyzz_yyzz, g_yz_0_yyzz_yyzzz, g_yz_0_yyzz_yzzz, g_yz_0_yyzz_yzzzz, g_yz_0_yyzz_zzzz, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzz_xxxx[k] = -g_z_0_yyzz_xxxx[k] - g_yz_0_yyzz_xxxx[k] * ab_y + g_yz_0_yyzz_xxxxy[k];

                g_yz_0_yyyzz_xxxy[k] = -g_z_0_yyzz_xxxy[k] - g_yz_0_yyzz_xxxy[k] * ab_y + g_yz_0_yyzz_xxxyy[k];

                g_yz_0_yyyzz_xxxz[k] = -g_z_0_yyzz_xxxz[k] - g_yz_0_yyzz_xxxz[k] * ab_y + g_yz_0_yyzz_xxxyz[k];

                g_yz_0_yyyzz_xxyy[k] = -g_z_0_yyzz_xxyy[k] - g_yz_0_yyzz_xxyy[k] * ab_y + g_yz_0_yyzz_xxyyy[k];

                g_yz_0_yyyzz_xxyz[k] = -g_z_0_yyzz_xxyz[k] - g_yz_0_yyzz_xxyz[k] * ab_y + g_yz_0_yyzz_xxyyz[k];

                g_yz_0_yyyzz_xxzz[k] = -g_z_0_yyzz_xxzz[k] - g_yz_0_yyzz_xxzz[k] * ab_y + g_yz_0_yyzz_xxyzz[k];

                g_yz_0_yyyzz_xyyy[k] = -g_z_0_yyzz_xyyy[k] - g_yz_0_yyzz_xyyy[k] * ab_y + g_yz_0_yyzz_xyyyy[k];

                g_yz_0_yyyzz_xyyz[k] = -g_z_0_yyzz_xyyz[k] - g_yz_0_yyzz_xyyz[k] * ab_y + g_yz_0_yyzz_xyyyz[k];

                g_yz_0_yyyzz_xyzz[k] = -g_z_0_yyzz_xyzz[k] - g_yz_0_yyzz_xyzz[k] * ab_y + g_yz_0_yyzz_xyyzz[k];

                g_yz_0_yyyzz_xzzz[k] = -g_z_0_yyzz_xzzz[k] - g_yz_0_yyzz_xzzz[k] * ab_y + g_yz_0_yyzz_xyzzz[k];

                g_yz_0_yyyzz_yyyy[k] = -g_z_0_yyzz_yyyy[k] - g_yz_0_yyzz_yyyy[k] * ab_y + g_yz_0_yyzz_yyyyy[k];

                g_yz_0_yyyzz_yyyz[k] = -g_z_0_yyzz_yyyz[k] - g_yz_0_yyzz_yyyz[k] * ab_y + g_yz_0_yyzz_yyyyz[k];

                g_yz_0_yyyzz_yyzz[k] = -g_z_0_yyzz_yyzz[k] - g_yz_0_yyzz_yyzz[k] * ab_y + g_yz_0_yyzz_yyyzz[k];

                g_yz_0_yyyzz_yzzz[k] = -g_z_0_yyzz_yzzz[k] - g_yz_0_yyzz_yzzz[k] * ab_y + g_yz_0_yyzz_yyzzz[k];

                g_yz_0_yyyzz_zzzz[k] = -g_z_0_yyzz_zzzz[k] - g_yz_0_yyzz_zzzz[k] * ab_y + g_yz_0_yyzz_yzzzz[k];
            }

            /// Set up 1530-1545 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1530 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1531 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1532 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1533 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1534 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1535 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1536 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1537 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1538 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1539 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1540 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1541 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1542 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1543 * ccomps * dcomps);

            auto g_yz_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1544 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzz_xxxx, g_yz_0_yyzzz_xxxy, g_yz_0_yyzzz_xxxz, g_yz_0_yyzzz_xxyy, g_yz_0_yyzzz_xxyz, g_yz_0_yyzzz_xxzz, g_yz_0_yyzzz_xyyy, g_yz_0_yyzzz_xyyz, g_yz_0_yyzzz_xyzz, g_yz_0_yyzzz_xzzz, g_yz_0_yyzzz_yyyy, g_yz_0_yyzzz_yyyz, g_yz_0_yyzzz_yyzz, g_yz_0_yyzzz_yzzz, g_yz_0_yyzzz_zzzz, g_yz_0_yzzz_xxxx, g_yz_0_yzzz_xxxxy, g_yz_0_yzzz_xxxy, g_yz_0_yzzz_xxxyy, g_yz_0_yzzz_xxxyz, g_yz_0_yzzz_xxxz, g_yz_0_yzzz_xxyy, g_yz_0_yzzz_xxyyy, g_yz_0_yzzz_xxyyz, g_yz_0_yzzz_xxyz, g_yz_0_yzzz_xxyzz, g_yz_0_yzzz_xxzz, g_yz_0_yzzz_xyyy, g_yz_0_yzzz_xyyyy, g_yz_0_yzzz_xyyyz, g_yz_0_yzzz_xyyz, g_yz_0_yzzz_xyyzz, g_yz_0_yzzz_xyzz, g_yz_0_yzzz_xyzzz, g_yz_0_yzzz_xzzz, g_yz_0_yzzz_yyyy, g_yz_0_yzzz_yyyyy, g_yz_0_yzzz_yyyyz, g_yz_0_yzzz_yyyz, g_yz_0_yzzz_yyyzz, g_yz_0_yzzz_yyzz, g_yz_0_yzzz_yyzzz, g_yz_0_yzzz_yzzz, g_yz_0_yzzz_yzzzz, g_yz_0_yzzz_zzzz, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzz_xxxx[k] = -g_z_0_yzzz_xxxx[k] - g_yz_0_yzzz_xxxx[k] * ab_y + g_yz_0_yzzz_xxxxy[k];

                g_yz_0_yyzzz_xxxy[k] = -g_z_0_yzzz_xxxy[k] - g_yz_0_yzzz_xxxy[k] * ab_y + g_yz_0_yzzz_xxxyy[k];

                g_yz_0_yyzzz_xxxz[k] = -g_z_0_yzzz_xxxz[k] - g_yz_0_yzzz_xxxz[k] * ab_y + g_yz_0_yzzz_xxxyz[k];

                g_yz_0_yyzzz_xxyy[k] = -g_z_0_yzzz_xxyy[k] - g_yz_0_yzzz_xxyy[k] * ab_y + g_yz_0_yzzz_xxyyy[k];

                g_yz_0_yyzzz_xxyz[k] = -g_z_0_yzzz_xxyz[k] - g_yz_0_yzzz_xxyz[k] * ab_y + g_yz_0_yzzz_xxyyz[k];

                g_yz_0_yyzzz_xxzz[k] = -g_z_0_yzzz_xxzz[k] - g_yz_0_yzzz_xxzz[k] * ab_y + g_yz_0_yzzz_xxyzz[k];

                g_yz_0_yyzzz_xyyy[k] = -g_z_0_yzzz_xyyy[k] - g_yz_0_yzzz_xyyy[k] * ab_y + g_yz_0_yzzz_xyyyy[k];

                g_yz_0_yyzzz_xyyz[k] = -g_z_0_yzzz_xyyz[k] - g_yz_0_yzzz_xyyz[k] * ab_y + g_yz_0_yzzz_xyyyz[k];

                g_yz_0_yyzzz_xyzz[k] = -g_z_0_yzzz_xyzz[k] - g_yz_0_yzzz_xyzz[k] * ab_y + g_yz_0_yzzz_xyyzz[k];

                g_yz_0_yyzzz_xzzz[k] = -g_z_0_yzzz_xzzz[k] - g_yz_0_yzzz_xzzz[k] * ab_y + g_yz_0_yzzz_xyzzz[k];

                g_yz_0_yyzzz_yyyy[k] = -g_z_0_yzzz_yyyy[k] - g_yz_0_yzzz_yyyy[k] * ab_y + g_yz_0_yzzz_yyyyy[k];

                g_yz_0_yyzzz_yyyz[k] = -g_z_0_yzzz_yyyz[k] - g_yz_0_yzzz_yyyz[k] * ab_y + g_yz_0_yzzz_yyyyz[k];

                g_yz_0_yyzzz_yyzz[k] = -g_z_0_yzzz_yyzz[k] - g_yz_0_yzzz_yyzz[k] * ab_y + g_yz_0_yzzz_yyyzz[k];

                g_yz_0_yyzzz_yzzz[k] = -g_z_0_yzzz_yzzz[k] - g_yz_0_yzzz_yzzz[k] * ab_y + g_yz_0_yzzz_yyzzz[k];

                g_yz_0_yyzzz_zzzz[k] = -g_z_0_yzzz_zzzz[k] - g_yz_0_yzzz_zzzz[k] * ab_y + g_yz_0_yzzz_yzzzz[k];
            }

            /// Set up 1545-1560 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1545 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1546 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1547 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1548 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1549 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1550 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1551 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1552 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1553 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1554 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1555 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1556 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1557 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1558 * ccomps * dcomps);

            auto g_yz_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1559 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzz_xxxx, g_yz_0_yzzzz_xxxy, g_yz_0_yzzzz_xxxz, g_yz_0_yzzzz_xxyy, g_yz_0_yzzzz_xxyz, g_yz_0_yzzzz_xxzz, g_yz_0_yzzzz_xyyy, g_yz_0_yzzzz_xyyz, g_yz_0_yzzzz_xyzz, g_yz_0_yzzzz_xzzz, g_yz_0_yzzzz_yyyy, g_yz_0_yzzzz_yyyz, g_yz_0_yzzzz_yyzz, g_yz_0_yzzzz_yzzz, g_yz_0_yzzzz_zzzz, g_yz_0_zzzz_xxxx, g_yz_0_zzzz_xxxxy, g_yz_0_zzzz_xxxy, g_yz_0_zzzz_xxxyy, g_yz_0_zzzz_xxxyz, g_yz_0_zzzz_xxxz, g_yz_0_zzzz_xxyy, g_yz_0_zzzz_xxyyy, g_yz_0_zzzz_xxyyz, g_yz_0_zzzz_xxyz, g_yz_0_zzzz_xxyzz, g_yz_0_zzzz_xxzz, g_yz_0_zzzz_xyyy, g_yz_0_zzzz_xyyyy, g_yz_0_zzzz_xyyyz, g_yz_0_zzzz_xyyz, g_yz_0_zzzz_xyyzz, g_yz_0_zzzz_xyzz, g_yz_0_zzzz_xyzzz, g_yz_0_zzzz_xzzz, g_yz_0_zzzz_yyyy, g_yz_0_zzzz_yyyyy, g_yz_0_zzzz_yyyyz, g_yz_0_zzzz_yyyz, g_yz_0_zzzz_yyyzz, g_yz_0_zzzz_yyzz, g_yz_0_zzzz_yyzzz, g_yz_0_zzzz_yzzz, g_yz_0_zzzz_yzzzz, g_yz_0_zzzz_zzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzz_xxxx[k] = -g_z_0_zzzz_xxxx[k] - g_yz_0_zzzz_xxxx[k] * ab_y + g_yz_0_zzzz_xxxxy[k];

                g_yz_0_yzzzz_xxxy[k] = -g_z_0_zzzz_xxxy[k] - g_yz_0_zzzz_xxxy[k] * ab_y + g_yz_0_zzzz_xxxyy[k];

                g_yz_0_yzzzz_xxxz[k] = -g_z_0_zzzz_xxxz[k] - g_yz_0_zzzz_xxxz[k] * ab_y + g_yz_0_zzzz_xxxyz[k];

                g_yz_0_yzzzz_xxyy[k] = -g_z_0_zzzz_xxyy[k] - g_yz_0_zzzz_xxyy[k] * ab_y + g_yz_0_zzzz_xxyyy[k];

                g_yz_0_yzzzz_xxyz[k] = -g_z_0_zzzz_xxyz[k] - g_yz_0_zzzz_xxyz[k] * ab_y + g_yz_0_zzzz_xxyyz[k];

                g_yz_0_yzzzz_xxzz[k] = -g_z_0_zzzz_xxzz[k] - g_yz_0_zzzz_xxzz[k] * ab_y + g_yz_0_zzzz_xxyzz[k];

                g_yz_0_yzzzz_xyyy[k] = -g_z_0_zzzz_xyyy[k] - g_yz_0_zzzz_xyyy[k] * ab_y + g_yz_0_zzzz_xyyyy[k];

                g_yz_0_yzzzz_xyyz[k] = -g_z_0_zzzz_xyyz[k] - g_yz_0_zzzz_xyyz[k] * ab_y + g_yz_0_zzzz_xyyyz[k];

                g_yz_0_yzzzz_xyzz[k] = -g_z_0_zzzz_xyzz[k] - g_yz_0_zzzz_xyzz[k] * ab_y + g_yz_0_zzzz_xyyzz[k];

                g_yz_0_yzzzz_xzzz[k] = -g_z_0_zzzz_xzzz[k] - g_yz_0_zzzz_xzzz[k] * ab_y + g_yz_0_zzzz_xyzzz[k];

                g_yz_0_yzzzz_yyyy[k] = -g_z_0_zzzz_yyyy[k] - g_yz_0_zzzz_yyyy[k] * ab_y + g_yz_0_zzzz_yyyyy[k];

                g_yz_0_yzzzz_yyyz[k] = -g_z_0_zzzz_yyyz[k] - g_yz_0_zzzz_yyyz[k] * ab_y + g_yz_0_zzzz_yyyyz[k];

                g_yz_0_yzzzz_yyzz[k] = -g_z_0_zzzz_yyzz[k] - g_yz_0_zzzz_yyzz[k] * ab_y + g_yz_0_zzzz_yyyzz[k];

                g_yz_0_yzzzz_yzzz[k] = -g_z_0_zzzz_yzzz[k] - g_yz_0_zzzz_yzzz[k] * ab_y + g_yz_0_zzzz_yyzzz[k];

                g_yz_0_yzzzz_zzzz[k] = -g_z_0_zzzz_zzzz[k] - g_yz_0_zzzz_zzzz[k] * ab_y + g_yz_0_zzzz_yzzzz[k];
            }

            /// Set up 1560-1575 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1560 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1561 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1562 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1563 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1564 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1565 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1566 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1567 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1568 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1569 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1570 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1571 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1572 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1573 * ccomps * dcomps);

            auto g_yz_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1574 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzzz, g_yz_0_zzzz_xxxx, g_yz_0_zzzz_xxxxz, g_yz_0_zzzz_xxxy, g_yz_0_zzzz_xxxyz, g_yz_0_zzzz_xxxz, g_yz_0_zzzz_xxxzz, g_yz_0_zzzz_xxyy, g_yz_0_zzzz_xxyyz, g_yz_0_zzzz_xxyz, g_yz_0_zzzz_xxyzz, g_yz_0_zzzz_xxzz, g_yz_0_zzzz_xxzzz, g_yz_0_zzzz_xyyy, g_yz_0_zzzz_xyyyz, g_yz_0_zzzz_xyyz, g_yz_0_zzzz_xyyzz, g_yz_0_zzzz_xyzz, g_yz_0_zzzz_xyzzz, g_yz_0_zzzz_xzzz, g_yz_0_zzzz_xzzzz, g_yz_0_zzzz_yyyy, g_yz_0_zzzz_yyyyz, g_yz_0_zzzz_yyyz, g_yz_0_zzzz_yyyzz, g_yz_0_zzzz_yyzz, g_yz_0_zzzz_yyzzz, g_yz_0_zzzz_yzzz, g_yz_0_zzzz_yzzzz, g_yz_0_zzzz_zzzz, g_yz_0_zzzz_zzzzz, g_yz_0_zzzzz_xxxx, g_yz_0_zzzzz_xxxy, g_yz_0_zzzzz_xxxz, g_yz_0_zzzzz_xxyy, g_yz_0_zzzzz_xxyz, g_yz_0_zzzzz_xxzz, g_yz_0_zzzzz_xyyy, g_yz_0_zzzzz_xyyz, g_yz_0_zzzzz_xyzz, g_yz_0_zzzzz_xzzz, g_yz_0_zzzzz_yyyy, g_yz_0_zzzzz_yyyz, g_yz_0_zzzzz_yyzz, g_yz_0_zzzzz_yzzz, g_yz_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzz_xxxx[k] = -g_y_0_zzzz_xxxx[k] - g_yz_0_zzzz_xxxx[k] * ab_z + g_yz_0_zzzz_xxxxz[k];

                g_yz_0_zzzzz_xxxy[k] = -g_y_0_zzzz_xxxy[k] - g_yz_0_zzzz_xxxy[k] * ab_z + g_yz_0_zzzz_xxxyz[k];

                g_yz_0_zzzzz_xxxz[k] = -g_y_0_zzzz_xxxz[k] - g_yz_0_zzzz_xxxz[k] * ab_z + g_yz_0_zzzz_xxxzz[k];

                g_yz_0_zzzzz_xxyy[k] = -g_y_0_zzzz_xxyy[k] - g_yz_0_zzzz_xxyy[k] * ab_z + g_yz_0_zzzz_xxyyz[k];

                g_yz_0_zzzzz_xxyz[k] = -g_y_0_zzzz_xxyz[k] - g_yz_0_zzzz_xxyz[k] * ab_z + g_yz_0_zzzz_xxyzz[k];

                g_yz_0_zzzzz_xxzz[k] = -g_y_0_zzzz_xxzz[k] - g_yz_0_zzzz_xxzz[k] * ab_z + g_yz_0_zzzz_xxzzz[k];

                g_yz_0_zzzzz_xyyy[k] = -g_y_0_zzzz_xyyy[k] - g_yz_0_zzzz_xyyy[k] * ab_z + g_yz_0_zzzz_xyyyz[k];

                g_yz_0_zzzzz_xyyz[k] = -g_y_0_zzzz_xyyz[k] - g_yz_0_zzzz_xyyz[k] * ab_z + g_yz_0_zzzz_xyyzz[k];

                g_yz_0_zzzzz_xyzz[k] = -g_y_0_zzzz_xyzz[k] - g_yz_0_zzzz_xyzz[k] * ab_z + g_yz_0_zzzz_xyzzz[k];

                g_yz_0_zzzzz_xzzz[k] = -g_y_0_zzzz_xzzz[k] - g_yz_0_zzzz_xzzz[k] * ab_z + g_yz_0_zzzz_xzzzz[k];

                g_yz_0_zzzzz_yyyy[k] = -g_y_0_zzzz_yyyy[k] - g_yz_0_zzzz_yyyy[k] * ab_z + g_yz_0_zzzz_yyyyz[k];

                g_yz_0_zzzzz_yyyz[k] = -g_y_0_zzzz_yyyz[k] - g_yz_0_zzzz_yyyz[k] * ab_z + g_yz_0_zzzz_yyyzz[k];

                g_yz_0_zzzzz_yyzz[k] = -g_y_0_zzzz_yyzz[k] - g_yz_0_zzzz_yyzz[k] * ab_z + g_yz_0_zzzz_yyzzz[k];

                g_yz_0_zzzzz_yzzz[k] = -g_y_0_zzzz_yzzz[k] - g_yz_0_zzzz_yzzz[k] * ab_z + g_yz_0_zzzz_yzzzz[k];

                g_yz_0_zzzzz_zzzz[k] = -g_y_0_zzzz_zzzz[k] - g_yz_0_zzzz_zzzz[k] * ab_z + g_yz_0_zzzz_zzzzz[k];
            }

            /// Set up 1575-1590 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxx_xxxx = cbuffer.data(hg_geom_20_off + 1575 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxxy = cbuffer.data(hg_geom_20_off + 1576 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxxz = cbuffer.data(hg_geom_20_off + 1577 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxyy = cbuffer.data(hg_geom_20_off + 1578 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxyz = cbuffer.data(hg_geom_20_off + 1579 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xxzz = cbuffer.data(hg_geom_20_off + 1580 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xyyy = cbuffer.data(hg_geom_20_off + 1581 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xyyz = cbuffer.data(hg_geom_20_off + 1582 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xyzz = cbuffer.data(hg_geom_20_off + 1583 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xzzz = cbuffer.data(hg_geom_20_off + 1584 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yyyy = cbuffer.data(hg_geom_20_off + 1585 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yyyz = cbuffer.data(hg_geom_20_off + 1586 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yyzz = cbuffer.data(hg_geom_20_off + 1587 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yzzz = cbuffer.data(hg_geom_20_off + 1588 * ccomps * dcomps);

            auto g_zz_0_xxxxx_zzzz = cbuffer.data(hg_geom_20_off + 1589 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxx_xxxx, g_zz_0_xxxx_xxxxx, g_zz_0_xxxx_xxxxy, g_zz_0_xxxx_xxxxz, g_zz_0_xxxx_xxxy, g_zz_0_xxxx_xxxyy, g_zz_0_xxxx_xxxyz, g_zz_0_xxxx_xxxz, g_zz_0_xxxx_xxxzz, g_zz_0_xxxx_xxyy, g_zz_0_xxxx_xxyyy, g_zz_0_xxxx_xxyyz, g_zz_0_xxxx_xxyz, g_zz_0_xxxx_xxyzz, g_zz_0_xxxx_xxzz, g_zz_0_xxxx_xxzzz, g_zz_0_xxxx_xyyy, g_zz_0_xxxx_xyyyy, g_zz_0_xxxx_xyyyz, g_zz_0_xxxx_xyyz, g_zz_0_xxxx_xyyzz, g_zz_0_xxxx_xyzz, g_zz_0_xxxx_xyzzz, g_zz_0_xxxx_xzzz, g_zz_0_xxxx_xzzzz, g_zz_0_xxxx_yyyy, g_zz_0_xxxx_yyyz, g_zz_0_xxxx_yyzz, g_zz_0_xxxx_yzzz, g_zz_0_xxxx_zzzz, g_zz_0_xxxxx_xxxx, g_zz_0_xxxxx_xxxy, g_zz_0_xxxxx_xxxz, g_zz_0_xxxxx_xxyy, g_zz_0_xxxxx_xxyz, g_zz_0_xxxxx_xxzz, g_zz_0_xxxxx_xyyy, g_zz_0_xxxxx_xyyz, g_zz_0_xxxxx_xyzz, g_zz_0_xxxxx_xzzz, g_zz_0_xxxxx_yyyy, g_zz_0_xxxxx_yyyz, g_zz_0_xxxxx_yyzz, g_zz_0_xxxxx_yzzz, g_zz_0_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxx_xxxx[k] = -g_zz_0_xxxx_xxxx[k] * ab_x + g_zz_0_xxxx_xxxxx[k];

                g_zz_0_xxxxx_xxxy[k] = -g_zz_0_xxxx_xxxy[k] * ab_x + g_zz_0_xxxx_xxxxy[k];

                g_zz_0_xxxxx_xxxz[k] = -g_zz_0_xxxx_xxxz[k] * ab_x + g_zz_0_xxxx_xxxxz[k];

                g_zz_0_xxxxx_xxyy[k] = -g_zz_0_xxxx_xxyy[k] * ab_x + g_zz_0_xxxx_xxxyy[k];

                g_zz_0_xxxxx_xxyz[k] = -g_zz_0_xxxx_xxyz[k] * ab_x + g_zz_0_xxxx_xxxyz[k];

                g_zz_0_xxxxx_xxzz[k] = -g_zz_0_xxxx_xxzz[k] * ab_x + g_zz_0_xxxx_xxxzz[k];

                g_zz_0_xxxxx_xyyy[k] = -g_zz_0_xxxx_xyyy[k] * ab_x + g_zz_0_xxxx_xxyyy[k];

                g_zz_0_xxxxx_xyyz[k] = -g_zz_0_xxxx_xyyz[k] * ab_x + g_zz_0_xxxx_xxyyz[k];

                g_zz_0_xxxxx_xyzz[k] = -g_zz_0_xxxx_xyzz[k] * ab_x + g_zz_0_xxxx_xxyzz[k];

                g_zz_0_xxxxx_xzzz[k] = -g_zz_0_xxxx_xzzz[k] * ab_x + g_zz_0_xxxx_xxzzz[k];

                g_zz_0_xxxxx_yyyy[k] = -g_zz_0_xxxx_yyyy[k] * ab_x + g_zz_0_xxxx_xyyyy[k];

                g_zz_0_xxxxx_yyyz[k] = -g_zz_0_xxxx_yyyz[k] * ab_x + g_zz_0_xxxx_xyyyz[k];

                g_zz_0_xxxxx_yyzz[k] = -g_zz_0_xxxx_yyzz[k] * ab_x + g_zz_0_xxxx_xyyzz[k];

                g_zz_0_xxxxx_yzzz[k] = -g_zz_0_xxxx_yzzz[k] * ab_x + g_zz_0_xxxx_xyzzz[k];

                g_zz_0_xxxxx_zzzz[k] = -g_zz_0_xxxx_zzzz[k] * ab_x + g_zz_0_xxxx_xzzzz[k];
            }

            /// Set up 1590-1605 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxy_xxxx = cbuffer.data(hg_geom_20_off + 1590 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxxy = cbuffer.data(hg_geom_20_off + 1591 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxxz = cbuffer.data(hg_geom_20_off + 1592 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxyy = cbuffer.data(hg_geom_20_off + 1593 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxyz = cbuffer.data(hg_geom_20_off + 1594 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xxzz = cbuffer.data(hg_geom_20_off + 1595 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xyyy = cbuffer.data(hg_geom_20_off + 1596 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xyyz = cbuffer.data(hg_geom_20_off + 1597 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xyzz = cbuffer.data(hg_geom_20_off + 1598 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xzzz = cbuffer.data(hg_geom_20_off + 1599 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yyyy = cbuffer.data(hg_geom_20_off + 1600 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yyyz = cbuffer.data(hg_geom_20_off + 1601 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yyzz = cbuffer.data(hg_geom_20_off + 1602 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yzzz = cbuffer.data(hg_geom_20_off + 1603 * ccomps * dcomps);

            auto g_zz_0_xxxxy_zzzz = cbuffer.data(hg_geom_20_off + 1604 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxy_xxxx, g_zz_0_xxxxy_xxxy, g_zz_0_xxxxy_xxxz, g_zz_0_xxxxy_xxyy, g_zz_0_xxxxy_xxyz, g_zz_0_xxxxy_xxzz, g_zz_0_xxxxy_xyyy, g_zz_0_xxxxy_xyyz, g_zz_0_xxxxy_xyzz, g_zz_0_xxxxy_xzzz, g_zz_0_xxxxy_yyyy, g_zz_0_xxxxy_yyyz, g_zz_0_xxxxy_yyzz, g_zz_0_xxxxy_yzzz, g_zz_0_xxxxy_zzzz, g_zz_0_xxxy_xxxx, g_zz_0_xxxy_xxxxx, g_zz_0_xxxy_xxxxy, g_zz_0_xxxy_xxxxz, g_zz_0_xxxy_xxxy, g_zz_0_xxxy_xxxyy, g_zz_0_xxxy_xxxyz, g_zz_0_xxxy_xxxz, g_zz_0_xxxy_xxxzz, g_zz_0_xxxy_xxyy, g_zz_0_xxxy_xxyyy, g_zz_0_xxxy_xxyyz, g_zz_0_xxxy_xxyz, g_zz_0_xxxy_xxyzz, g_zz_0_xxxy_xxzz, g_zz_0_xxxy_xxzzz, g_zz_0_xxxy_xyyy, g_zz_0_xxxy_xyyyy, g_zz_0_xxxy_xyyyz, g_zz_0_xxxy_xyyz, g_zz_0_xxxy_xyyzz, g_zz_0_xxxy_xyzz, g_zz_0_xxxy_xyzzz, g_zz_0_xxxy_xzzz, g_zz_0_xxxy_xzzzz, g_zz_0_xxxy_yyyy, g_zz_0_xxxy_yyyz, g_zz_0_xxxy_yyzz, g_zz_0_xxxy_yzzz, g_zz_0_xxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxy_xxxx[k] = -g_zz_0_xxxy_xxxx[k] * ab_x + g_zz_0_xxxy_xxxxx[k];

                g_zz_0_xxxxy_xxxy[k] = -g_zz_0_xxxy_xxxy[k] * ab_x + g_zz_0_xxxy_xxxxy[k];

                g_zz_0_xxxxy_xxxz[k] = -g_zz_0_xxxy_xxxz[k] * ab_x + g_zz_0_xxxy_xxxxz[k];

                g_zz_0_xxxxy_xxyy[k] = -g_zz_0_xxxy_xxyy[k] * ab_x + g_zz_0_xxxy_xxxyy[k];

                g_zz_0_xxxxy_xxyz[k] = -g_zz_0_xxxy_xxyz[k] * ab_x + g_zz_0_xxxy_xxxyz[k];

                g_zz_0_xxxxy_xxzz[k] = -g_zz_0_xxxy_xxzz[k] * ab_x + g_zz_0_xxxy_xxxzz[k];

                g_zz_0_xxxxy_xyyy[k] = -g_zz_0_xxxy_xyyy[k] * ab_x + g_zz_0_xxxy_xxyyy[k];

                g_zz_0_xxxxy_xyyz[k] = -g_zz_0_xxxy_xyyz[k] * ab_x + g_zz_0_xxxy_xxyyz[k];

                g_zz_0_xxxxy_xyzz[k] = -g_zz_0_xxxy_xyzz[k] * ab_x + g_zz_0_xxxy_xxyzz[k];

                g_zz_0_xxxxy_xzzz[k] = -g_zz_0_xxxy_xzzz[k] * ab_x + g_zz_0_xxxy_xxzzz[k];

                g_zz_0_xxxxy_yyyy[k] = -g_zz_0_xxxy_yyyy[k] * ab_x + g_zz_0_xxxy_xyyyy[k];

                g_zz_0_xxxxy_yyyz[k] = -g_zz_0_xxxy_yyyz[k] * ab_x + g_zz_0_xxxy_xyyyz[k];

                g_zz_0_xxxxy_yyzz[k] = -g_zz_0_xxxy_yyzz[k] * ab_x + g_zz_0_xxxy_xyyzz[k];

                g_zz_0_xxxxy_yzzz[k] = -g_zz_0_xxxy_yzzz[k] * ab_x + g_zz_0_xxxy_xyzzz[k];

                g_zz_0_xxxxy_zzzz[k] = -g_zz_0_xxxy_zzzz[k] * ab_x + g_zz_0_xxxy_xzzzz[k];
            }

            /// Set up 1605-1620 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxz_xxxx = cbuffer.data(hg_geom_20_off + 1605 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxxy = cbuffer.data(hg_geom_20_off + 1606 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxxz = cbuffer.data(hg_geom_20_off + 1607 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxyy = cbuffer.data(hg_geom_20_off + 1608 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxyz = cbuffer.data(hg_geom_20_off + 1609 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xxzz = cbuffer.data(hg_geom_20_off + 1610 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xyyy = cbuffer.data(hg_geom_20_off + 1611 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xyyz = cbuffer.data(hg_geom_20_off + 1612 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xyzz = cbuffer.data(hg_geom_20_off + 1613 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xzzz = cbuffer.data(hg_geom_20_off + 1614 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yyyy = cbuffer.data(hg_geom_20_off + 1615 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yyyz = cbuffer.data(hg_geom_20_off + 1616 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yyzz = cbuffer.data(hg_geom_20_off + 1617 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yzzz = cbuffer.data(hg_geom_20_off + 1618 * ccomps * dcomps);

            auto g_zz_0_xxxxz_zzzz = cbuffer.data(hg_geom_20_off + 1619 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxz_xxxx, g_zz_0_xxxxz_xxxy, g_zz_0_xxxxz_xxxz, g_zz_0_xxxxz_xxyy, g_zz_0_xxxxz_xxyz, g_zz_0_xxxxz_xxzz, g_zz_0_xxxxz_xyyy, g_zz_0_xxxxz_xyyz, g_zz_0_xxxxz_xyzz, g_zz_0_xxxxz_xzzz, g_zz_0_xxxxz_yyyy, g_zz_0_xxxxz_yyyz, g_zz_0_xxxxz_yyzz, g_zz_0_xxxxz_yzzz, g_zz_0_xxxxz_zzzz, g_zz_0_xxxz_xxxx, g_zz_0_xxxz_xxxxx, g_zz_0_xxxz_xxxxy, g_zz_0_xxxz_xxxxz, g_zz_0_xxxz_xxxy, g_zz_0_xxxz_xxxyy, g_zz_0_xxxz_xxxyz, g_zz_0_xxxz_xxxz, g_zz_0_xxxz_xxxzz, g_zz_0_xxxz_xxyy, g_zz_0_xxxz_xxyyy, g_zz_0_xxxz_xxyyz, g_zz_0_xxxz_xxyz, g_zz_0_xxxz_xxyzz, g_zz_0_xxxz_xxzz, g_zz_0_xxxz_xxzzz, g_zz_0_xxxz_xyyy, g_zz_0_xxxz_xyyyy, g_zz_0_xxxz_xyyyz, g_zz_0_xxxz_xyyz, g_zz_0_xxxz_xyyzz, g_zz_0_xxxz_xyzz, g_zz_0_xxxz_xyzzz, g_zz_0_xxxz_xzzz, g_zz_0_xxxz_xzzzz, g_zz_0_xxxz_yyyy, g_zz_0_xxxz_yyyz, g_zz_0_xxxz_yyzz, g_zz_0_xxxz_yzzz, g_zz_0_xxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxz_xxxx[k] = -g_zz_0_xxxz_xxxx[k] * ab_x + g_zz_0_xxxz_xxxxx[k];

                g_zz_0_xxxxz_xxxy[k] = -g_zz_0_xxxz_xxxy[k] * ab_x + g_zz_0_xxxz_xxxxy[k];

                g_zz_0_xxxxz_xxxz[k] = -g_zz_0_xxxz_xxxz[k] * ab_x + g_zz_0_xxxz_xxxxz[k];

                g_zz_0_xxxxz_xxyy[k] = -g_zz_0_xxxz_xxyy[k] * ab_x + g_zz_0_xxxz_xxxyy[k];

                g_zz_0_xxxxz_xxyz[k] = -g_zz_0_xxxz_xxyz[k] * ab_x + g_zz_0_xxxz_xxxyz[k];

                g_zz_0_xxxxz_xxzz[k] = -g_zz_0_xxxz_xxzz[k] * ab_x + g_zz_0_xxxz_xxxzz[k];

                g_zz_0_xxxxz_xyyy[k] = -g_zz_0_xxxz_xyyy[k] * ab_x + g_zz_0_xxxz_xxyyy[k];

                g_zz_0_xxxxz_xyyz[k] = -g_zz_0_xxxz_xyyz[k] * ab_x + g_zz_0_xxxz_xxyyz[k];

                g_zz_0_xxxxz_xyzz[k] = -g_zz_0_xxxz_xyzz[k] * ab_x + g_zz_0_xxxz_xxyzz[k];

                g_zz_0_xxxxz_xzzz[k] = -g_zz_0_xxxz_xzzz[k] * ab_x + g_zz_0_xxxz_xxzzz[k];

                g_zz_0_xxxxz_yyyy[k] = -g_zz_0_xxxz_yyyy[k] * ab_x + g_zz_0_xxxz_xyyyy[k];

                g_zz_0_xxxxz_yyyz[k] = -g_zz_0_xxxz_yyyz[k] * ab_x + g_zz_0_xxxz_xyyyz[k];

                g_zz_0_xxxxz_yyzz[k] = -g_zz_0_xxxz_yyzz[k] * ab_x + g_zz_0_xxxz_xyyzz[k];

                g_zz_0_xxxxz_yzzz[k] = -g_zz_0_xxxz_yzzz[k] * ab_x + g_zz_0_xxxz_xyzzz[k];

                g_zz_0_xxxxz_zzzz[k] = -g_zz_0_xxxz_zzzz[k] * ab_x + g_zz_0_xxxz_xzzzz[k];
            }

            /// Set up 1620-1635 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyy_xxxx = cbuffer.data(hg_geom_20_off + 1620 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxxy = cbuffer.data(hg_geom_20_off + 1621 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxxz = cbuffer.data(hg_geom_20_off + 1622 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxyy = cbuffer.data(hg_geom_20_off + 1623 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxyz = cbuffer.data(hg_geom_20_off + 1624 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xxzz = cbuffer.data(hg_geom_20_off + 1625 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xyyy = cbuffer.data(hg_geom_20_off + 1626 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xyyz = cbuffer.data(hg_geom_20_off + 1627 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xyzz = cbuffer.data(hg_geom_20_off + 1628 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xzzz = cbuffer.data(hg_geom_20_off + 1629 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yyyy = cbuffer.data(hg_geom_20_off + 1630 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yyyz = cbuffer.data(hg_geom_20_off + 1631 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yyzz = cbuffer.data(hg_geom_20_off + 1632 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yzzz = cbuffer.data(hg_geom_20_off + 1633 * ccomps * dcomps);

            auto g_zz_0_xxxyy_zzzz = cbuffer.data(hg_geom_20_off + 1634 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyy_xxxx, g_zz_0_xxxyy_xxxy, g_zz_0_xxxyy_xxxz, g_zz_0_xxxyy_xxyy, g_zz_0_xxxyy_xxyz, g_zz_0_xxxyy_xxzz, g_zz_0_xxxyy_xyyy, g_zz_0_xxxyy_xyyz, g_zz_0_xxxyy_xyzz, g_zz_0_xxxyy_xzzz, g_zz_0_xxxyy_yyyy, g_zz_0_xxxyy_yyyz, g_zz_0_xxxyy_yyzz, g_zz_0_xxxyy_yzzz, g_zz_0_xxxyy_zzzz, g_zz_0_xxyy_xxxx, g_zz_0_xxyy_xxxxx, g_zz_0_xxyy_xxxxy, g_zz_0_xxyy_xxxxz, g_zz_0_xxyy_xxxy, g_zz_0_xxyy_xxxyy, g_zz_0_xxyy_xxxyz, g_zz_0_xxyy_xxxz, g_zz_0_xxyy_xxxzz, g_zz_0_xxyy_xxyy, g_zz_0_xxyy_xxyyy, g_zz_0_xxyy_xxyyz, g_zz_0_xxyy_xxyz, g_zz_0_xxyy_xxyzz, g_zz_0_xxyy_xxzz, g_zz_0_xxyy_xxzzz, g_zz_0_xxyy_xyyy, g_zz_0_xxyy_xyyyy, g_zz_0_xxyy_xyyyz, g_zz_0_xxyy_xyyz, g_zz_0_xxyy_xyyzz, g_zz_0_xxyy_xyzz, g_zz_0_xxyy_xyzzz, g_zz_0_xxyy_xzzz, g_zz_0_xxyy_xzzzz, g_zz_0_xxyy_yyyy, g_zz_0_xxyy_yyyz, g_zz_0_xxyy_yyzz, g_zz_0_xxyy_yzzz, g_zz_0_xxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyy_xxxx[k] = -g_zz_0_xxyy_xxxx[k] * ab_x + g_zz_0_xxyy_xxxxx[k];

                g_zz_0_xxxyy_xxxy[k] = -g_zz_0_xxyy_xxxy[k] * ab_x + g_zz_0_xxyy_xxxxy[k];

                g_zz_0_xxxyy_xxxz[k] = -g_zz_0_xxyy_xxxz[k] * ab_x + g_zz_0_xxyy_xxxxz[k];

                g_zz_0_xxxyy_xxyy[k] = -g_zz_0_xxyy_xxyy[k] * ab_x + g_zz_0_xxyy_xxxyy[k];

                g_zz_0_xxxyy_xxyz[k] = -g_zz_0_xxyy_xxyz[k] * ab_x + g_zz_0_xxyy_xxxyz[k];

                g_zz_0_xxxyy_xxzz[k] = -g_zz_0_xxyy_xxzz[k] * ab_x + g_zz_0_xxyy_xxxzz[k];

                g_zz_0_xxxyy_xyyy[k] = -g_zz_0_xxyy_xyyy[k] * ab_x + g_zz_0_xxyy_xxyyy[k];

                g_zz_0_xxxyy_xyyz[k] = -g_zz_0_xxyy_xyyz[k] * ab_x + g_zz_0_xxyy_xxyyz[k];

                g_zz_0_xxxyy_xyzz[k] = -g_zz_0_xxyy_xyzz[k] * ab_x + g_zz_0_xxyy_xxyzz[k];

                g_zz_0_xxxyy_xzzz[k] = -g_zz_0_xxyy_xzzz[k] * ab_x + g_zz_0_xxyy_xxzzz[k];

                g_zz_0_xxxyy_yyyy[k] = -g_zz_0_xxyy_yyyy[k] * ab_x + g_zz_0_xxyy_xyyyy[k];

                g_zz_0_xxxyy_yyyz[k] = -g_zz_0_xxyy_yyyz[k] * ab_x + g_zz_0_xxyy_xyyyz[k];

                g_zz_0_xxxyy_yyzz[k] = -g_zz_0_xxyy_yyzz[k] * ab_x + g_zz_0_xxyy_xyyzz[k];

                g_zz_0_xxxyy_yzzz[k] = -g_zz_0_xxyy_yzzz[k] * ab_x + g_zz_0_xxyy_xyzzz[k];

                g_zz_0_xxxyy_zzzz[k] = -g_zz_0_xxyy_zzzz[k] * ab_x + g_zz_0_xxyy_xzzzz[k];
            }

            /// Set up 1635-1650 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyz_xxxx = cbuffer.data(hg_geom_20_off + 1635 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxxy = cbuffer.data(hg_geom_20_off + 1636 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxxz = cbuffer.data(hg_geom_20_off + 1637 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxyy = cbuffer.data(hg_geom_20_off + 1638 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxyz = cbuffer.data(hg_geom_20_off + 1639 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xxzz = cbuffer.data(hg_geom_20_off + 1640 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xyyy = cbuffer.data(hg_geom_20_off + 1641 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xyyz = cbuffer.data(hg_geom_20_off + 1642 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xyzz = cbuffer.data(hg_geom_20_off + 1643 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xzzz = cbuffer.data(hg_geom_20_off + 1644 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yyyy = cbuffer.data(hg_geom_20_off + 1645 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yyyz = cbuffer.data(hg_geom_20_off + 1646 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yyzz = cbuffer.data(hg_geom_20_off + 1647 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yzzz = cbuffer.data(hg_geom_20_off + 1648 * ccomps * dcomps);

            auto g_zz_0_xxxyz_zzzz = cbuffer.data(hg_geom_20_off + 1649 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyz_xxxx, g_zz_0_xxxyz_xxxy, g_zz_0_xxxyz_xxxz, g_zz_0_xxxyz_xxyy, g_zz_0_xxxyz_xxyz, g_zz_0_xxxyz_xxzz, g_zz_0_xxxyz_xyyy, g_zz_0_xxxyz_xyyz, g_zz_0_xxxyz_xyzz, g_zz_0_xxxyz_xzzz, g_zz_0_xxxyz_yyyy, g_zz_0_xxxyz_yyyz, g_zz_0_xxxyz_yyzz, g_zz_0_xxxyz_yzzz, g_zz_0_xxxyz_zzzz, g_zz_0_xxyz_xxxx, g_zz_0_xxyz_xxxxx, g_zz_0_xxyz_xxxxy, g_zz_0_xxyz_xxxxz, g_zz_0_xxyz_xxxy, g_zz_0_xxyz_xxxyy, g_zz_0_xxyz_xxxyz, g_zz_0_xxyz_xxxz, g_zz_0_xxyz_xxxzz, g_zz_0_xxyz_xxyy, g_zz_0_xxyz_xxyyy, g_zz_0_xxyz_xxyyz, g_zz_0_xxyz_xxyz, g_zz_0_xxyz_xxyzz, g_zz_0_xxyz_xxzz, g_zz_0_xxyz_xxzzz, g_zz_0_xxyz_xyyy, g_zz_0_xxyz_xyyyy, g_zz_0_xxyz_xyyyz, g_zz_0_xxyz_xyyz, g_zz_0_xxyz_xyyzz, g_zz_0_xxyz_xyzz, g_zz_0_xxyz_xyzzz, g_zz_0_xxyz_xzzz, g_zz_0_xxyz_xzzzz, g_zz_0_xxyz_yyyy, g_zz_0_xxyz_yyyz, g_zz_0_xxyz_yyzz, g_zz_0_xxyz_yzzz, g_zz_0_xxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyz_xxxx[k] = -g_zz_0_xxyz_xxxx[k] * ab_x + g_zz_0_xxyz_xxxxx[k];

                g_zz_0_xxxyz_xxxy[k] = -g_zz_0_xxyz_xxxy[k] * ab_x + g_zz_0_xxyz_xxxxy[k];

                g_zz_0_xxxyz_xxxz[k] = -g_zz_0_xxyz_xxxz[k] * ab_x + g_zz_0_xxyz_xxxxz[k];

                g_zz_0_xxxyz_xxyy[k] = -g_zz_0_xxyz_xxyy[k] * ab_x + g_zz_0_xxyz_xxxyy[k];

                g_zz_0_xxxyz_xxyz[k] = -g_zz_0_xxyz_xxyz[k] * ab_x + g_zz_0_xxyz_xxxyz[k];

                g_zz_0_xxxyz_xxzz[k] = -g_zz_0_xxyz_xxzz[k] * ab_x + g_zz_0_xxyz_xxxzz[k];

                g_zz_0_xxxyz_xyyy[k] = -g_zz_0_xxyz_xyyy[k] * ab_x + g_zz_0_xxyz_xxyyy[k];

                g_zz_0_xxxyz_xyyz[k] = -g_zz_0_xxyz_xyyz[k] * ab_x + g_zz_0_xxyz_xxyyz[k];

                g_zz_0_xxxyz_xyzz[k] = -g_zz_0_xxyz_xyzz[k] * ab_x + g_zz_0_xxyz_xxyzz[k];

                g_zz_0_xxxyz_xzzz[k] = -g_zz_0_xxyz_xzzz[k] * ab_x + g_zz_0_xxyz_xxzzz[k];

                g_zz_0_xxxyz_yyyy[k] = -g_zz_0_xxyz_yyyy[k] * ab_x + g_zz_0_xxyz_xyyyy[k];

                g_zz_0_xxxyz_yyyz[k] = -g_zz_0_xxyz_yyyz[k] * ab_x + g_zz_0_xxyz_xyyyz[k];

                g_zz_0_xxxyz_yyzz[k] = -g_zz_0_xxyz_yyzz[k] * ab_x + g_zz_0_xxyz_xyyzz[k];

                g_zz_0_xxxyz_yzzz[k] = -g_zz_0_xxyz_yzzz[k] * ab_x + g_zz_0_xxyz_xyzzz[k];

                g_zz_0_xxxyz_zzzz[k] = -g_zz_0_xxyz_zzzz[k] * ab_x + g_zz_0_xxyz_xzzzz[k];
            }

            /// Set up 1650-1665 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzz_xxxx = cbuffer.data(hg_geom_20_off + 1650 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxxy = cbuffer.data(hg_geom_20_off + 1651 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxxz = cbuffer.data(hg_geom_20_off + 1652 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxyy = cbuffer.data(hg_geom_20_off + 1653 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxyz = cbuffer.data(hg_geom_20_off + 1654 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xxzz = cbuffer.data(hg_geom_20_off + 1655 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xyyy = cbuffer.data(hg_geom_20_off + 1656 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xyyz = cbuffer.data(hg_geom_20_off + 1657 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xyzz = cbuffer.data(hg_geom_20_off + 1658 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xzzz = cbuffer.data(hg_geom_20_off + 1659 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yyyy = cbuffer.data(hg_geom_20_off + 1660 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yyyz = cbuffer.data(hg_geom_20_off + 1661 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yyzz = cbuffer.data(hg_geom_20_off + 1662 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yzzz = cbuffer.data(hg_geom_20_off + 1663 * ccomps * dcomps);

            auto g_zz_0_xxxzz_zzzz = cbuffer.data(hg_geom_20_off + 1664 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzz_xxxx, g_zz_0_xxxzz_xxxy, g_zz_0_xxxzz_xxxz, g_zz_0_xxxzz_xxyy, g_zz_0_xxxzz_xxyz, g_zz_0_xxxzz_xxzz, g_zz_0_xxxzz_xyyy, g_zz_0_xxxzz_xyyz, g_zz_0_xxxzz_xyzz, g_zz_0_xxxzz_xzzz, g_zz_0_xxxzz_yyyy, g_zz_0_xxxzz_yyyz, g_zz_0_xxxzz_yyzz, g_zz_0_xxxzz_yzzz, g_zz_0_xxxzz_zzzz, g_zz_0_xxzz_xxxx, g_zz_0_xxzz_xxxxx, g_zz_0_xxzz_xxxxy, g_zz_0_xxzz_xxxxz, g_zz_0_xxzz_xxxy, g_zz_0_xxzz_xxxyy, g_zz_0_xxzz_xxxyz, g_zz_0_xxzz_xxxz, g_zz_0_xxzz_xxxzz, g_zz_0_xxzz_xxyy, g_zz_0_xxzz_xxyyy, g_zz_0_xxzz_xxyyz, g_zz_0_xxzz_xxyz, g_zz_0_xxzz_xxyzz, g_zz_0_xxzz_xxzz, g_zz_0_xxzz_xxzzz, g_zz_0_xxzz_xyyy, g_zz_0_xxzz_xyyyy, g_zz_0_xxzz_xyyyz, g_zz_0_xxzz_xyyz, g_zz_0_xxzz_xyyzz, g_zz_0_xxzz_xyzz, g_zz_0_xxzz_xyzzz, g_zz_0_xxzz_xzzz, g_zz_0_xxzz_xzzzz, g_zz_0_xxzz_yyyy, g_zz_0_xxzz_yyyz, g_zz_0_xxzz_yyzz, g_zz_0_xxzz_yzzz, g_zz_0_xxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzz_xxxx[k] = -g_zz_0_xxzz_xxxx[k] * ab_x + g_zz_0_xxzz_xxxxx[k];

                g_zz_0_xxxzz_xxxy[k] = -g_zz_0_xxzz_xxxy[k] * ab_x + g_zz_0_xxzz_xxxxy[k];

                g_zz_0_xxxzz_xxxz[k] = -g_zz_0_xxzz_xxxz[k] * ab_x + g_zz_0_xxzz_xxxxz[k];

                g_zz_0_xxxzz_xxyy[k] = -g_zz_0_xxzz_xxyy[k] * ab_x + g_zz_0_xxzz_xxxyy[k];

                g_zz_0_xxxzz_xxyz[k] = -g_zz_0_xxzz_xxyz[k] * ab_x + g_zz_0_xxzz_xxxyz[k];

                g_zz_0_xxxzz_xxzz[k] = -g_zz_0_xxzz_xxzz[k] * ab_x + g_zz_0_xxzz_xxxzz[k];

                g_zz_0_xxxzz_xyyy[k] = -g_zz_0_xxzz_xyyy[k] * ab_x + g_zz_0_xxzz_xxyyy[k];

                g_zz_0_xxxzz_xyyz[k] = -g_zz_0_xxzz_xyyz[k] * ab_x + g_zz_0_xxzz_xxyyz[k];

                g_zz_0_xxxzz_xyzz[k] = -g_zz_0_xxzz_xyzz[k] * ab_x + g_zz_0_xxzz_xxyzz[k];

                g_zz_0_xxxzz_xzzz[k] = -g_zz_0_xxzz_xzzz[k] * ab_x + g_zz_0_xxzz_xxzzz[k];

                g_zz_0_xxxzz_yyyy[k] = -g_zz_0_xxzz_yyyy[k] * ab_x + g_zz_0_xxzz_xyyyy[k];

                g_zz_0_xxxzz_yyyz[k] = -g_zz_0_xxzz_yyyz[k] * ab_x + g_zz_0_xxzz_xyyyz[k];

                g_zz_0_xxxzz_yyzz[k] = -g_zz_0_xxzz_yyzz[k] * ab_x + g_zz_0_xxzz_xyyzz[k];

                g_zz_0_xxxzz_yzzz[k] = -g_zz_0_xxzz_yzzz[k] * ab_x + g_zz_0_xxzz_xyzzz[k];

                g_zz_0_xxxzz_zzzz[k] = -g_zz_0_xxzz_zzzz[k] * ab_x + g_zz_0_xxzz_xzzzz[k];
            }

            /// Set up 1665-1680 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyy_xxxx = cbuffer.data(hg_geom_20_off + 1665 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxxy = cbuffer.data(hg_geom_20_off + 1666 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxxz = cbuffer.data(hg_geom_20_off + 1667 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxyy = cbuffer.data(hg_geom_20_off + 1668 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxyz = cbuffer.data(hg_geom_20_off + 1669 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xxzz = cbuffer.data(hg_geom_20_off + 1670 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xyyy = cbuffer.data(hg_geom_20_off + 1671 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xyyz = cbuffer.data(hg_geom_20_off + 1672 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xyzz = cbuffer.data(hg_geom_20_off + 1673 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xzzz = cbuffer.data(hg_geom_20_off + 1674 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yyyy = cbuffer.data(hg_geom_20_off + 1675 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yyyz = cbuffer.data(hg_geom_20_off + 1676 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yyzz = cbuffer.data(hg_geom_20_off + 1677 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yzzz = cbuffer.data(hg_geom_20_off + 1678 * ccomps * dcomps);

            auto g_zz_0_xxyyy_zzzz = cbuffer.data(hg_geom_20_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyy_xxxx, g_zz_0_xxyyy_xxxy, g_zz_0_xxyyy_xxxz, g_zz_0_xxyyy_xxyy, g_zz_0_xxyyy_xxyz, g_zz_0_xxyyy_xxzz, g_zz_0_xxyyy_xyyy, g_zz_0_xxyyy_xyyz, g_zz_0_xxyyy_xyzz, g_zz_0_xxyyy_xzzz, g_zz_0_xxyyy_yyyy, g_zz_0_xxyyy_yyyz, g_zz_0_xxyyy_yyzz, g_zz_0_xxyyy_yzzz, g_zz_0_xxyyy_zzzz, g_zz_0_xyyy_xxxx, g_zz_0_xyyy_xxxxx, g_zz_0_xyyy_xxxxy, g_zz_0_xyyy_xxxxz, g_zz_0_xyyy_xxxy, g_zz_0_xyyy_xxxyy, g_zz_0_xyyy_xxxyz, g_zz_0_xyyy_xxxz, g_zz_0_xyyy_xxxzz, g_zz_0_xyyy_xxyy, g_zz_0_xyyy_xxyyy, g_zz_0_xyyy_xxyyz, g_zz_0_xyyy_xxyz, g_zz_0_xyyy_xxyzz, g_zz_0_xyyy_xxzz, g_zz_0_xyyy_xxzzz, g_zz_0_xyyy_xyyy, g_zz_0_xyyy_xyyyy, g_zz_0_xyyy_xyyyz, g_zz_0_xyyy_xyyz, g_zz_0_xyyy_xyyzz, g_zz_0_xyyy_xyzz, g_zz_0_xyyy_xyzzz, g_zz_0_xyyy_xzzz, g_zz_0_xyyy_xzzzz, g_zz_0_xyyy_yyyy, g_zz_0_xyyy_yyyz, g_zz_0_xyyy_yyzz, g_zz_0_xyyy_yzzz, g_zz_0_xyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyy_xxxx[k] = -g_zz_0_xyyy_xxxx[k] * ab_x + g_zz_0_xyyy_xxxxx[k];

                g_zz_0_xxyyy_xxxy[k] = -g_zz_0_xyyy_xxxy[k] * ab_x + g_zz_0_xyyy_xxxxy[k];

                g_zz_0_xxyyy_xxxz[k] = -g_zz_0_xyyy_xxxz[k] * ab_x + g_zz_0_xyyy_xxxxz[k];

                g_zz_0_xxyyy_xxyy[k] = -g_zz_0_xyyy_xxyy[k] * ab_x + g_zz_0_xyyy_xxxyy[k];

                g_zz_0_xxyyy_xxyz[k] = -g_zz_0_xyyy_xxyz[k] * ab_x + g_zz_0_xyyy_xxxyz[k];

                g_zz_0_xxyyy_xxzz[k] = -g_zz_0_xyyy_xxzz[k] * ab_x + g_zz_0_xyyy_xxxzz[k];

                g_zz_0_xxyyy_xyyy[k] = -g_zz_0_xyyy_xyyy[k] * ab_x + g_zz_0_xyyy_xxyyy[k];

                g_zz_0_xxyyy_xyyz[k] = -g_zz_0_xyyy_xyyz[k] * ab_x + g_zz_0_xyyy_xxyyz[k];

                g_zz_0_xxyyy_xyzz[k] = -g_zz_0_xyyy_xyzz[k] * ab_x + g_zz_0_xyyy_xxyzz[k];

                g_zz_0_xxyyy_xzzz[k] = -g_zz_0_xyyy_xzzz[k] * ab_x + g_zz_0_xyyy_xxzzz[k];

                g_zz_0_xxyyy_yyyy[k] = -g_zz_0_xyyy_yyyy[k] * ab_x + g_zz_0_xyyy_xyyyy[k];

                g_zz_0_xxyyy_yyyz[k] = -g_zz_0_xyyy_yyyz[k] * ab_x + g_zz_0_xyyy_xyyyz[k];

                g_zz_0_xxyyy_yyzz[k] = -g_zz_0_xyyy_yyzz[k] * ab_x + g_zz_0_xyyy_xyyzz[k];

                g_zz_0_xxyyy_yzzz[k] = -g_zz_0_xyyy_yzzz[k] * ab_x + g_zz_0_xyyy_xyzzz[k];

                g_zz_0_xxyyy_zzzz[k] = -g_zz_0_xyyy_zzzz[k] * ab_x + g_zz_0_xyyy_xzzzz[k];
            }

            /// Set up 1680-1695 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyz_xxxx = cbuffer.data(hg_geom_20_off + 1680 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxxy = cbuffer.data(hg_geom_20_off + 1681 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxxz = cbuffer.data(hg_geom_20_off + 1682 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxyy = cbuffer.data(hg_geom_20_off + 1683 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxyz = cbuffer.data(hg_geom_20_off + 1684 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xxzz = cbuffer.data(hg_geom_20_off + 1685 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xyyy = cbuffer.data(hg_geom_20_off + 1686 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xyyz = cbuffer.data(hg_geom_20_off + 1687 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xyzz = cbuffer.data(hg_geom_20_off + 1688 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xzzz = cbuffer.data(hg_geom_20_off + 1689 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yyyy = cbuffer.data(hg_geom_20_off + 1690 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yyyz = cbuffer.data(hg_geom_20_off + 1691 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yyzz = cbuffer.data(hg_geom_20_off + 1692 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yzzz = cbuffer.data(hg_geom_20_off + 1693 * ccomps * dcomps);

            auto g_zz_0_xxyyz_zzzz = cbuffer.data(hg_geom_20_off + 1694 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyz_xxxx, g_zz_0_xxyyz_xxxy, g_zz_0_xxyyz_xxxz, g_zz_0_xxyyz_xxyy, g_zz_0_xxyyz_xxyz, g_zz_0_xxyyz_xxzz, g_zz_0_xxyyz_xyyy, g_zz_0_xxyyz_xyyz, g_zz_0_xxyyz_xyzz, g_zz_0_xxyyz_xzzz, g_zz_0_xxyyz_yyyy, g_zz_0_xxyyz_yyyz, g_zz_0_xxyyz_yyzz, g_zz_0_xxyyz_yzzz, g_zz_0_xxyyz_zzzz, g_zz_0_xyyz_xxxx, g_zz_0_xyyz_xxxxx, g_zz_0_xyyz_xxxxy, g_zz_0_xyyz_xxxxz, g_zz_0_xyyz_xxxy, g_zz_0_xyyz_xxxyy, g_zz_0_xyyz_xxxyz, g_zz_0_xyyz_xxxz, g_zz_0_xyyz_xxxzz, g_zz_0_xyyz_xxyy, g_zz_0_xyyz_xxyyy, g_zz_0_xyyz_xxyyz, g_zz_0_xyyz_xxyz, g_zz_0_xyyz_xxyzz, g_zz_0_xyyz_xxzz, g_zz_0_xyyz_xxzzz, g_zz_0_xyyz_xyyy, g_zz_0_xyyz_xyyyy, g_zz_0_xyyz_xyyyz, g_zz_0_xyyz_xyyz, g_zz_0_xyyz_xyyzz, g_zz_0_xyyz_xyzz, g_zz_0_xyyz_xyzzz, g_zz_0_xyyz_xzzz, g_zz_0_xyyz_xzzzz, g_zz_0_xyyz_yyyy, g_zz_0_xyyz_yyyz, g_zz_0_xyyz_yyzz, g_zz_0_xyyz_yzzz, g_zz_0_xyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyz_xxxx[k] = -g_zz_0_xyyz_xxxx[k] * ab_x + g_zz_0_xyyz_xxxxx[k];

                g_zz_0_xxyyz_xxxy[k] = -g_zz_0_xyyz_xxxy[k] * ab_x + g_zz_0_xyyz_xxxxy[k];

                g_zz_0_xxyyz_xxxz[k] = -g_zz_0_xyyz_xxxz[k] * ab_x + g_zz_0_xyyz_xxxxz[k];

                g_zz_0_xxyyz_xxyy[k] = -g_zz_0_xyyz_xxyy[k] * ab_x + g_zz_0_xyyz_xxxyy[k];

                g_zz_0_xxyyz_xxyz[k] = -g_zz_0_xyyz_xxyz[k] * ab_x + g_zz_0_xyyz_xxxyz[k];

                g_zz_0_xxyyz_xxzz[k] = -g_zz_0_xyyz_xxzz[k] * ab_x + g_zz_0_xyyz_xxxzz[k];

                g_zz_0_xxyyz_xyyy[k] = -g_zz_0_xyyz_xyyy[k] * ab_x + g_zz_0_xyyz_xxyyy[k];

                g_zz_0_xxyyz_xyyz[k] = -g_zz_0_xyyz_xyyz[k] * ab_x + g_zz_0_xyyz_xxyyz[k];

                g_zz_0_xxyyz_xyzz[k] = -g_zz_0_xyyz_xyzz[k] * ab_x + g_zz_0_xyyz_xxyzz[k];

                g_zz_0_xxyyz_xzzz[k] = -g_zz_0_xyyz_xzzz[k] * ab_x + g_zz_0_xyyz_xxzzz[k];

                g_zz_0_xxyyz_yyyy[k] = -g_zz_0_xyyz_yyyy[k] * ab_x + g_zz_0_xyyz_xyyyy[k];

                g_zz_0_xxyyz_yyyz[k] = -g_zz_0_xyyz_yyyz[k] * ab_x + g_zz_0_xyyz_xyyyz[k];

                g_zz_0_xxyyz_yyzz[k] = -g_zz_0_xyyz_yyzz[k] * ab_x + g_zz_0_xyyz_xyyzz[k];

                g_zz_0_xxyyz_yzzz[k] = -g_zz_0_xyyz_yzzz[k] * ab_x + g_zz_0_xyyz_xyzzz[k];

                g_zz_0_xxyyz_zzzz[k] = -g_zz_0_xyyz_zzzz[k] * ab_x + g_zz_0_xyyz_xzzzz[k];
            }

            /// Set up 1695-1710 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzz_xxxx = cbuffer.data(hg_geom_20_off + 1695 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxxy = cbuffer.data(hg_geom_20_off + 1696 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxxz = cbuffer.data(hg_geom_20_off + 1697 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxyy = cbuffer.data(hg_geom_20_off + 1698 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxyz = cbuffer.data(hg_geom_20_off + 1699 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xxzz = cbuffer.data(hg_geom_20_off + 1700 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xyyy = cbuffer.data(hg_geom_20_off + 1701 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xyyz = cbuffer.data(hg_geom_20_off + 1702 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xyzz = cbuffer.data(hg_geom_20_off + 1703 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xzzz = cbuffer.data(hg_geom_20_off + 1704 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yyyy = cbuffer.data(hg_geom_20_off + 1705 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yyyz = cbuffer.data(hg_geom_20_off + 1706 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yyzz = cbuffer.data(hg_geom_20_off + 1707 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yzzz = cbuffer.data(hg_geom_20_off + 1708 * ccomps * dcomps);

            auto g_zz_0_xxyzz_zzzz = cbuffer.data(hg_geom_20_off + 1709 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzz_xxxx, g_zz_0_xxyzz_xxxy, g_zz_0_xxyzz_xxxz, g_zz_0_xxyzz_xxyy, g_zz_0_xxyzz_xxyz, g_zz_0_xxyzz_xxzz, g_zz_0_xxyzz_xyyy, g_zz_0_xxyzz_xyyz, g_zz_0_xxyzz_xyzz, g_zz_0_xxyzz_xzzz, g_zz_0_xxyzz_yyyy, g_zz_0_xxyzz_yyyz, g_zz_0_xxyzz_yyzz, g_zz_0_xxyzz_yzzz, g_zz_0_xxyzz_zzzz, g_zz_0_xyzz_xxxx, g_zz_0_xyzz_xxxxx, g_zz_0_xyzz_xxxxy, g_zz_0_xyzz_xxxxz, g_zz_0_xyzz_xxxy, g_zz_0_xyzz_xxxyy, g_zz_0_xyzz_xxxyz, g_zz_0_xyzz_xxxz, g_zz_0_xyzz_xxxzz, g_zz_0_xyzz_xxyy, g_zz_0_xyzz_xxyyy, g_zz_0_xyzz_xxyyz, g_zz_0_xyzz_xxyz, g_zz_0_xyzz_xxyzz, g_zz_0_xyzz_xxzz, g_zz_0_xyzz_xxzzz, g_zz_0_xyzz_xyyy, g_zz_0_xyzz_xyyyy, g_zz_0_xyzz_xyyyz, g_zz_0_xyzz_xyyz, g_zz_0_xyzz_xyyzz, g_zz_0_xyzz_xyzz, g_zz_0_xyzz_xyzzz, g_zz_0_xyzz_xzzz, g_zz_0_xyzz_xzzzz, g_zz_0_xyzz_yyyy, g_zz_0_xyzz_yyyz, g_zz_0_xyzz_yyzz, g_zz_0_xyzz_yzzz, g_zz_0_xyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzz_xxxx[k] = -g_zz_0_xyzz_xxxx[k] * ab_x + g_zz_0_xyzz_xxxxx[k];

                g_zz_0_xxyzz_xxxy[k] = -g_zz_0_xyzz_xxxy[k] * ab_x + g_zz_0_xyzz_xxxxy[k];

                g_zz_0_xxyzz_xxxz[k] = -g_zz_0_xyzz_xxxz[k] * ab_x + g_zz_0_xyzz_xxxxz[k];

                g_zz_0_xxyzz_xxyy[k] = -g_zz_0_xyzz_xxyy[k] * ab_x + g_zz_0_xyzz_xxxyy[k];

                g_zz_0_xxyzz_xxyz[k] = -g_zz_0_xyzz_xxyz[k] * ab_x + g_zz_0_xyzz_xxxyz[k];

                g_zz_0_xxyzz_xxzz[k] = -g_zz_0_xyzz_xxzz[k] * ab_x + g_zz_0_xyzz_xxxzz[k];

                g_zz_0_xxyzz_xyyy[k] = -g_zz_0_xyzz_xyyy[k] * ab_x + g_zz_0_xyzz_xxyyy[k];

                g_zz_0_xxyzz_xyyz[k] = -g_zz_0_xyzz_xyyz[k] * ab_x + g_zz_0_xyzz_xxyyz[k];

                g_zz_0_xxyzz_xyzz[k] = -g_zz_0_xyzz_xyzz[k] * ab_x + g_zz_0_xyzz_xxyzz[k];

                g_zz_0_xxyzz_xzzz[k] = -g_zz_0_xyzz_xzzz[k] * ab_x + g_zz_0_xyzz_xxzzz[k];

                g_zz_0_xxyzz_yyyy[k] = -g_zz_0_xyzz_yyyy[k] * ab_x + g_zz_0_xyzz_xyyyy[k];

                g_zz_0_xxyzz_yyyz[k] = -g_zz_0_xyzz_yyyz[k] * ab_x + g_zz_0_xyzz_xyyyz[k];

                g_zz_0_xxyzz_yyzz[k] = -g_zz_0_xyzz_yyzz[k] * ab_x + g_zz_0_xyzz_xyyzz[k];

                g_zz_0_xxyzz_yzzz[k] = -g_zz_0_xyzz_yzzz[k] * ab_x + g_zz_0_xyzz_xyzzz[k];

                g_zz_0_xxyzz_zzzz[k] = -g_zz_0_xyzz_zzzz[k] * ab_x + g_zz_0_xyzz_xzzzz[k];
            }

            /// Set up 1710-1725 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzz_xxxx = cbuffer.data(hg_geom_20_off + 1710 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxxy = cbuffer.data(hg_geom_20_off + 1711 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxxz = cbuffer.data(hg_geom_20_off + 1712 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxyy = cbuffer.data(hg_geom_20_off + 1713 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxyz = cbuffer.data(hg_geom_20_off + 1714 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xxzz = cbuffer.data(hg_geom_20_off + 1715 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xyyy = cbuffer.data(hg_geom_20_off + 1716 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xyyz = cbuffer.data(hg_geom_20_off + 1717 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xyzz = cbuffer.data(hg_geom_20_off + 1718 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xzzz = cbuffer.data(hg_geom_20_off + 1719 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yyyy = cbuffer.data(hg_geom_20_off + 1720 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yyyz = cbuffer.data(hg_geom_20_off + 1721 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yyzz = cbuffer.data(hg_geom_20_off + 1722 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yzzz = cbuffer.data(hg_geom_20_off + 1723 * ccomps * dcomps);

            auto g_zz_0_xxzzz_zzzz = cbuffer.data(hg_geom_20_off + 1724 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzz_xxxx, g_zz_0_xxzzz_xxxy, g_zz_0_xxzzz_xxxz, g_zz_0_xxzzz_xxyy, g_zz_0_xxzzz_xxyz, g_zz_0_xxzzz_xxzz, g_zz_0_xxzzz_xyyy, g_zz_0_xxzzz_xyyz, g_zz_0_xxzzz_xyzz, g_zz_0_xxzzz_xzzz, g_zz_0_xxzzz_yyyy, g_zz_0_xxzzz_yyyz, g_zz_0_xxzzz_yyzz, g_zz_0_xxzzz_yzzz, g_zz_0_xxzzz_zzzz, g_zz_0_xzzz_xxxx, g_zz_0_xzzz_xxxxx, g_zz_0_xzzz_xxxxy, g_zz_0_xzzz_xxxxz, g_zz_0_xzzz_xxxy, g_zz_0_xzzz_xxxyy, g_zz_0_xzzz_xxxyz, g_zz_0_xzzz_xxxz, g_zz_0_xzzz_xxxzz, g_zz_0_xzzz_xxyy, g_zz_0_xzzz_xxyyy, g_zz_0_xzzz_xxyyz, g_zz_0_xzzz_xxyz, g_zz_0_xzzz_xxyzz, g_zz_0_xzzz_xxzz, g_zz_0_xzzz_xxzzz, g_zz_0_xzzz_xyyy, g_zz_0_xzzz_xyyyy, g_zz_0_xzzz_xyyyz, g_zz_0_xzzz_xyyz, g_zz_0_xzzz_xyyzz, g_zz_0_xzzz_xyzz, g_zz_0_xzzz_xyzzz, g_zz_0_xzzz_xzzz, g_zz_0_xzzz_xzzzz, g_zz_0_xzzz_yyyy, g_zz_0_xzzz_yyyz, g_zz_0_xzzz_yyzz, g_zz_0_xzzz_yzzz, g_zz_0_xzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzz_xxxx[k] = -g_zz_0_xzzz_xxxx[k] * ab_x + g_zz_0_xzzz_xxxxx[k];

                g_zz_0_xxzzz_xxxy[k] = -g_zz_0_xzzz_xxxy[k] * ab_x + g_zz_0_xzzz_xxxxy[k];

                g_zz_0_xxzzz_xxxz[k] = -g_zz_0_xzzz_xxxz[k] * ab_x + g_zz_0_xzzz_xxxxz[k];

                g_zz_0_xxzzz_xxyy[k] = -g_zz_0_xzzz_xxyy[k] * ab_x + g_zz_0_xzzz_xxxyy[k];

                g_zz_0_xxzzz_xxyz[k] = -g_zz_0_xzzz_xxyz[k] * ab_x + g_zz_0_xzzz_xxxyz[k];

                g_zz_0_xxzzz_xxzz[k] = -g_zz_0_xzzz_xxzz[k] * ab_x + g_zz_0_xzzz_xxxzz[k];

                g_zz_0_xxzzz_xyyy[k] = -g_zz_0_xzzz_xyyy[k] * ab_x + g_zz_0_xzzz_xxyyy[k];

                g_zz_0_xxzzz_xyyz[k] = -g_zz_0_xzzz_xyyz[k] * ab_x + g_zz_0_xzzz_xxyyz[k];

                g_zz_0_xxzzz_xyzz[k] = -g_zz_0_xzzz_xyzz[k] * ab_x + g_zz_0_xzzz_xxyzz[k];

                g_zz_0_xxzzz_xzzz[k] = -g_zz_0_xzzz_xzzz[k] * ab_x + g_zz_0_xzzz_xxzzz[k];

                g_zz_0_xxzzz_yyyy[k] = -g_zz_0_xzzz_yyyy[k] * ab_x + g_zz_0_xzzz_xyyyy[k];

                g_zz_0_xxzzz_yyyz[k] = -g_zz_0_xzzz_yyyz[k] * ab_x + g_zz_0_xzzz_xyyyz[k];

                g_zz_0_xxzzz_yyzz[k] = -g_zz_0_xzzz_yyzz[k] * ab_x + g_zz_0_xzzz_xyyzz[k];

                g_zz_0_xxzzz_yzzz[k] = -g_zz_0_xzzz_yzzz[k] * ab_x + g_zz_0_xzzz_xyzzz[k];

                g_zz_0_xxzzz_zzzz[k] = -g_zz_0_xzzz_zzzz[k] * ab_x + g_zz_0_xzzz_xzzzz[k];
            }

            /// Set up 1725-1740 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1725 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1726 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1727 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1728 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1729 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1730 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1731 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1732 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1733 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1734 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1735 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1736 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1737 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1738 * ccomps * dcomps);

            auto g_zz_0_xyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1739 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyy_xxxx, g_zz_0_xyyyy_xxxy, g_zz_0_xyyyy_xxxz, g_zz_0_xyyyy_xxyy, g_zz_0_xyyyy_xxyz, g_zz_0_xyyyy_xxzz, g_zz_0_xyyyy_xyyy, g_zz_0_xyyyy_xyyz, g_zz_0_xyyyy_xyzz, g_zz_0_xyyyy_xzzz, g_zz_0_xyyyy_yyyy, g_zz_0_xyyyy_yyyz, g_zz_0_xyyyy_yyzz, g_zz_0_xyyyy_yzzz, g_zz_0_xyyyy_zzzz, g_zz_0_yyyy_xxxx, g_zz_0_yyyy_xxxxx, g_zz_0_yyyy_xxxxy, g_zz_0_yyyy_xxxxz, g_zz_0_yyyy_xxxy, g_zz_0_yyyy_xxxyy, g_zz_0_yyyy_xxxyz, g_zz_0_yyyy_xxxz, g_zz_0_yyyy_xxxzz, g_zz_0_yyyy_xxyy, g_zz_0_yyyy_xxyyy, g_zz_0_yyyy_xxyyz, g_zz_0_yyyy_xxyz, g_zz_0_yyyy_xxyzz, g_zz_0_yyyy_xxzz, g_zz_0_yyyy_xxzzz, g_zz_0_yyyy_xyyy, g_zz_0_yyyy_xyyyy, g_zz_0_yyyy_xyyyz, g_zz_0_yyyy_xyyz, g_zz_0_yyyy_xyyzz, g_zz_0_yyyy_xyzz, g_zz_0_yyyy_xyzzz, g_zz_0_yyyy_xzzz, g_zz_0_yyyy_xzzzz, g_zz_0_yyyy_yyyy, g_zz_0_yyyy_yyyz, g_zz_0_yyyy_yyzz, g_zz_0_yyyy_yzzz, g_zz_0_yyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyy_xxxx[k] = -g_zz_0_yyyy_xxxx[k] * ab_x + g_zz_0_yyyy_xxxxx[k];

                g_zz_0_xyyyy_xxxy[k] = -g_zz_0_yyyy_xxxy[k] * ab_x + g_zz_0_yyyy_xxxxy[k];

                g_zz_0_xyyyy_xxxz[k] = -g_zz_0_yyyy_xxxz[k] * ab_x + g_zz_0_yyyy_xxxxz[k];

                g_zz_0_xyyyy_xxyy[k] = -g_zz_0_yyyy_xxyy[k] * ab_x + g_zz_0_yyyy_xxxyy[k];

                g_zz_0_xyyyy_xxyz[k] = -g_zz_0_yyyy_xxyz[k] * ab_x + g_zz_0_yyyy_xxxyz[k];

                g_zz_0_xyyyy_xxzz[k] = -g_zz_0_yyyy_xxzz[k] * ab_x + g_zz_0_yyyy_xxxzz[k];

                g_zz_0_xyyyy_xyyy[k] = -g_zz_0_yyyy_xyyy[k] * ab_x + g_zz_0_yyyy_xxyyy[k];

                g_zz_0_xyyyy_xyyz[k] = -g_zz_0_yyyy_xyyz[k] * ab_x + g_zz_0_yyyy_xxyyz[k];

                g_zz_0_xyyyy_xyzz[k] = -g_zz_0_yyyy_xyzz[k] * ab_x + g_zz_0_yyyy_xxyzz[k];

                g_zz_0_xyyyy_xzzz[k] = -g_zz_0_yyyy_xzzz[k] * ab_x + g_zz_0_yyyy_xxzzz[k];

                g_zz_0_xyyyy_yyyy[k] = -g_zz_0_yyyy_yyyy[k] * ab_x + g_zz_0_yyyy_xyyyy[k];

                g_zz_0_xyyyy_yyyz[k] = -g_zz_0_yyyy_yyyz[k] * ab_x + g_zz_0_yyyy_xyyyz[k];

                g_zz_0_xyyyy_yyzz[k] = -g_zz_0_yyyy_yyzz[k] * ab_x + g_zz_0_yyyy_xyyzz[k];

                g_zz_0_xyyyy_yzzz[k] = -g_zz_0_yyyy_yzzz[k] * ab_x + g_zz_0_yyyy_xyzzz[k];

                g_zz_0_xyyyy_zzzz[k] = -g_zz_0_yyyy_zzzz[k] * ab_x + g_zz_0_yyyy_xzzzz[k];
            }

            /// Set up 1740-1755 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1740 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1741 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1742 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1743 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1744 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1745 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1746 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1747 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1748 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1749 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1750 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1751 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1752 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1753 * ccomps * dcomps);

            auto g_zz_0_xyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1754 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyz_xxxx, g_zz_0_xyyyz_xxxy, g_zz_0_xyyyz_xxxz, g_zz_0_xyyyz_xxyy, g_zz_0_xyyyz_xxyz, g_zz_0_xyyyz_xxzz, g_zz_0_xyyyz_xyyy, g_zz_0_xyyyz_xyyz, g_zz_0_xyyyz_xyzz, g_zz_0_xyyyz_xzzz, g_zz_0_xyyyz_yyyy, g_zz_0_xyyyz_yyyz, g_zz_0_xyyyz_yyzz, g_zz_0_xyyyz_yzzz, g_zz_0_xyyyz_zzzz, g_zz_0_yyyz_xxxx, g_zz_0_yyyz_xxxxx, g_zz_0_yyyz_xxxxy, g_zz_0_yyyz_xxxxz, g_zz_0_yyyz_xxxy, g_zz_0_yyyz_xxxyy, g_zz_0_yyyz_xxxyz, g_zz_0_yyyz_xxxz, g_zz_0_yyyz_xxxzz, g_zz_0_yyyz_xxyy, g_zz_0_yyyz_xxyyy, g_zz_0_yyyz_xxyyz, g_zz_0_yyyz_xxyz, g_zz_0_yyyz_xxyzz, g_zz_0_yyyz_xxzz, g_zz_0_yyyz_xxzzz, g_zz_0_yyyz_xyyy, g_zz_0_yyyz_xyyyy, g_zz_0_yyyz_xyyyz, g_zz_0_yyyz_xyyz, g_zz_0_yyyz_xyyzz, g_zz_0_yyyz_xyzz, g_zz_0_yyyz_xyzzz, g_zz_0_yyyz_xzzz, g_zz_0_yyyz_xzzzz, g_zz_0_yyyz_yyyy, g_zz_0_yyyz_yyyz, g_zz_0_yyyz_yyzz, g_zz_0_yyyz_yzzz, g_zz_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyz_xxxx[k] = -g_zz_0_yyyz_xxxx[k] * ab_x + g_zz_0_yyyz_xxxxx[k];

                g_zz_0_xyyyz_xxxy[k] = -g_zz_0_yyyz_xxxy[k] * ab_x + g_zz_0_yyyz_xxxxy[k];

                g_zz_0_xyyyz_xxxz[k] = -g_zz_0_yyyz_xxxz[k] * ab_x + g_zz_0_yyyz_xxxxz[k];

                g_zz_0_xyyyz_xxyy[k] = -g_zz_0_yyyz_xxyy[k] * ab_x + g_zz_0_yyyz_xxxyy[k];

                g_zz_0_xyyyz_xxyz[k] = -g_zz_0_yyyz_xxyz[k] * ab_x + g_zz_0_yyyz_xxxyz[k];

                g_zz_0_xyyyz_xxzz[k] = -g_zz_0_yyyz_xxzz[k] * ab_x + g_zz_0_yyyz_xxxzz[k];

                g_zz_0_xyyyz_xyyy[k] = -g_zz_0_yyyz_xyyy[k] * ab_x + g_zz_0_yyyz_xxyyy[k];

                g_zz_0_xyyyz_xyyz[k] = -g_zz_0_yyyz_xyyz[k] * ab_x + g_zz_0_yyyz_xxyyz[k];

                g_zz_0_xyyyz_xyzz[k] = -g_zz_0_yyyz_xyzz[k] * ab_x + g_zz_0_yyyz_xxyzz[k];

                g_zz_0_xyyyz_xzzz[k] = -g_zz_0_yyyz_xzzz[k] * ab_x + g_zz_0_yyyz_xxzzz[k];

                g_zz_0_xyyyz_yyyy[k] = -g_zz_0_yyyz_yyyy[k] * ab_x + g_zz_0_yyyz_xyyyy[k];

                g_zz_0_xyyyz_yyyz[k] = -g_zz_0_yyyz_yyyz[k] * ab_x + g_zz_0_yyyz_xyyyz[k];

                g_zz_0_xyyyz_yyzz[k] = -g_zz_0_yyyz_yyzz[k] * ab_x + g_zz_0_yyyz_xyyzz[k];

                g_zz_0_xyyyz_yzzz[k] = -g_zz_0_yyyz_yzzz[k] * ab_x + g_zz_0_yyyz_xyzzz[k];

                g_zz_0_xyyyz_zzzz[k] = -g_zz_0_yyyz_zzzz[k] * ab_x + g_zz_0_yyyz_xzzzz[k];
            }

            /// Set up 1755-1770 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1755 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1756 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1757 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1758 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1759 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1760 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1761 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1762 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1763 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1764 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1765 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1766 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1767 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1768 * ccomps * dcomps);

            auto g_zz_0_xyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1769 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzz_xxxx, g_zz_0_xyyzz_xxxy, g_zz_0_xyyzz_xxxz, g_zz_0_xyyzz_xxyy, g_zz_0_xyyzz_xxyz, g_zz_0_xyyzz_xxzz, g_zz_0_xyyzz_xyyy, g_zz_0_xyyzz_xyyz, g_zz_0_xyyzz_xyzz, g_zz_0_xyyzz_xzzz, g_zz_0_xyyzz_yyyy, g_zz_0_xyyzz_yyyz, g_zz_0_xyyzz_yyzz, g_zz_0_xyyzz_yzzz, g_zz_0_xyyzz_zzzz, g_zz_0_yyzz_xxxx, g_zz_0_yyzz_xxxxx, g_zz_0_yyzz_xxxxy, g_zz_0_yyzz_xxxxz, g_zz_0_yyzz_xxxy, g_zz_0_yyzz_xxxyy, g_zz_0_yyzz_xxxyz, g_zz_0_yyzz_xxxz, g_zz_0_yyzz_xxxzz, g_zz_0_yyzz_xxyy, g_zz_0_yyzz_xxyyy, g_zz_0_yyzz_xxyyz, g_zz_0_yyzz_xxyz, g_zz_0_yyzz_xxyzz, g_zz_0_yyzz_xxzz, g_zz_0_yyzz_xxzzz, g_zz_0_yyzz_xyyy, g_zz_0_yyzz_xyyyy, g_zz_0_yyzz_xyyyz, g_zz_0_yyzz_xyyz, g_zz_0_yyzz_xyyzz, g_zz_0_yyzz_xyzz, g_zz_0_yyzz_xyzzz, g_zz_0_yyzz_xzzz, g_zz_0_yyzz_xzzzz, g_zz_0_yyzz_yyyy, g_zz_0_yyzz_yyyz, g_zz_0_yyzz_yyzz, g_zz_0_yyzz_yzzz, g_zz_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzz_xxxx[k] = -g_zz_0_yyzz_xxxx[k] * ab_x + g_zz_0_yyzz_xxxxx[k];

                g_zz_0_xyyzz_xxxy[k] = -g_zz_0_yyzz_xxxy[k] * ab_x + g_zz_0_yyzz_xxxxy[k];

                g_zz_0_xyyzz_xxxz[k] = -g_zz_0_yyzz_xxxz[k] * ab_x + g_zz_0_yyzz_xxxxz[k];

                g_zz_0_xyyzz_xxyy[k] = -g_zz_0_yyzz_xxyy[k] * ab_x + g_zz_0_yyzz_xxxyy[k];

                g_zz_0_xyyzz_xxyz[k] = -g_zz_0_yyzz_xxyz[k] * ab_x + g_zz_0_yyzz_xxxyz[k];

                g_zz_0_xyyzz_xxzz[k] = -g_zz_0_yyzz_xxzz[k] * ab_x + g_zz_0_yyzz_xxxzz[k];

                g_zz_0_xyyzz_xyyy[k] = -g_zz_0_yyzz_xyyy[k] * ab_x + g_zz_0_yyzz_xxyyy[k];

                g_zz_0_xyyzz_xyyz[k] = -g_zz_0_yyzz_xyyz[k] * ab_x + g_zz_0_yyzz_xxyyz[k];

                g_zz_0_xyyzz_xyzz[k] = -g_zz_0_yyzz_xyzz[k] * ab_x + g_zz_0_yyzz_xxyzz[k];

                g_zz_0_xyyzz_xzzz[k] = -g_zz_0_yyzz_xzzz[k] * ab_x + g_zz_0_yyzz_xxzzz[k];

                g_zz_0_xyyzz_yyyy[k] = -g_zz_0_yyzz_yyyy[k] * ab_x + g_zz_0_yyzz_xyyyy[k];

                g_zz_0_xyyzz_yyyz[k] = -g_zz_0_yyzz_yyyz[k] * ab_x + g_zz_0_yyzz_xyyyz[k];

                g_zz_0_xyyzz_yyzz[k] = -g_zz_0_yyzz_yyzz[k] * ab_x + g_zz_0_yyzz_xyyzz[k];

                g_zz_0_xyyzz_yzzz[k] = -g_zz_0_yyzz_yzzz[k] * ab_x + g_zz_0_yyzz_xyzzz[k];

                g_zz_0_xyyzz_zzzz[k] = -g_zz_0_yyzz_zzzz[k] * ab_x + g_zz_0_yyzz_xzzzz[k];
            }

            /// Set up 1770-1785 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1770 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1771 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1772 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1773 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1774 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1775 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1776 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1777 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1778 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1779 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1780 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1781 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1782 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1783 * ccomps * dcomps);

            auto g_zz_0_xyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1784 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzz_xxxx, g_zz_0_xyzzz_xxxy, g_zz_0_xyzzz_xxxz, g_zz_0_xyzzz_xxyy, g_zz_0_xyzzz_xxyz, g_zz_0_xyzzz_xxzz, g_zz_0_xyzzz_xyyy, g_zz_0_xyzzz_xyyz, g_zz_0_xyzzz_xyzz, g_zz_0_xyzzz_xzzz, g_zz_0_xyzzz_yyyy, g_zz_0_xyzzz_yyyz, g_zz_0_xyzzz_yyzz, g_zz_0_xyzzz_yzzz, g_zz_0_xyzzz_zzzz, g_zz_0_yzzz_xxxx, g_zz_0_yzzz_xxxxx, g_zz_0_yzzz_xxxxy, g_zz_0_yzzz_xxxxz, g_zz_0_yzzz_xxxy, g_zz_0_yzzz_xxxyy, g_zz_0_yzzz_xxxyz, g_zz_0_yzzz_xxxz, g_zz_0_yzzz_xxxzz, g_zz_0_yzzz_xxyy, g_zz_0_yzzz_xxyyy, g_zz_0_yzzz_xxyyz, g_zz_0_yzzz_xxyz, g_zz_0_yzzz_xxyzz, g_zz_0_yzzz_xxzz, g_zz_0_yzzz_xxzzz, g_zz_0_yzzz_xyyy, g_zz_0_yzzz_xyyyy, g_zz_0_yzzz_xyyyz, g_zz_0_yzzz_xyyz, g_zz_0_yzzz_xyyzz, g_zz_0_yzzz_xyzz, g_zz_0_yzzz_xyzzz, g_zz_0_yzzz_xzzz, g_zz_0_yzzz_xzzzz, g_zz_0_yzzz_yyyy, g_zz_0_yzzz_yyyz, g_zz_0_yzzz_yyzz, g_zz_0_yzzz_yzzz, g_zz_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzz_xxxx[k] = -g_zz_0_yzzz_xxxx[k] * ab_x + g_zz_0_yzzz_xxxxx[k];

                g_zz_0_xyzzz_xxxy[k] = -g_zz_0_yzzz_xxxy[k] * ab_x + g_zz_0_yzzz_xxxxy[k];

                g_zz_0_xyzzz_xxxz[k] = -g_zz_0_yzzz_xxxz[k] * ab_x + g_zz_0_yzzz_xxxxz[k];

                g_zz_0_xyzzz_xxyy[k] = -g_zz_0_yzzz_xxyy[k] * ab_x + g_zz_0_yzzz_xxxyy[k];

                g_zz_0_xyzzz_xxyz[k] = -g_zz_0_yzzz_xxyz[k] * ab_x + g_zz_0_yzzz_xxxyz[k];

                g_zz_0_xyzzz_xxzz[k] = -g_zz_0_yzzz_xxzz[k] * ab_x + g_zz_0_yzzz_xxxzz[k];

                g_zz_0_xyzzz_xyyy[k] = -g_zz_0_yzzz_xyyy[k] * ab_x + g_zz_0_yzzz_xxyyy[k];

                g_zz_0_xyzzz_xyyz[k] = -g_zz_0_yzzz_xyyz[k] * ab_x + g_zz_0_yzzz_xxyyz[k];

                g_zz_0_xyzzz_xyzz[k] = -g_zz_0_yzzz_xyzz[k] * ab_x + g_zz_0_yzzz_xxyzz[k];

                g_zz_0_xyzzz_xzzz[k] = -g_zz_0_yzzz_xzzz[k] * ab_x + g_zz_0_yzzz_xxzzz[k];

                g_zz_0_xyzzz_yyyy[k] = -g_zz_0_yzzz_yyyy[k] * ab_x + g_zz_0_yzzz_xyyyy[k];

                g_zz_0_xyzzz_yyyz[k] = -g_zz_0_yzzz_yyyz[k] * ab_x + g_zz_0_yzzz_xyyyz[k];

                g_zz_0_xyzzz_yyzz[k] = -g_zz_0_yzzz_yyzz[k] * ab_x + g_zz_0_yzzz_xyyzz[k];

                g_zz_0_xyzzz_yzzz[k] = -g_zz_0_yzzz_yzzz[k] * ab_x + g_zz_0_yzzz_xyzzz[k];

                g_zz_0_xyzzz_zzzz[k] = -g_zz_0_yzzz_zzzz[k] * ab_x + g_zz_0_yzzz_xzzzz[k];
            }

            /// Set up 1785-1800 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1785 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1786 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1787 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1788 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1789 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1790 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1791 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1792 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1793 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1794 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1795 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1796 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1797 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1798 * ccomps * dcomps);

            auto g_zz_0_xzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1799 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzz_xxxx, g_zz_0_xzzzz_xxxy, g_zz_0_xzzzz_xxxz, g_zz_0_xzzzz_xxyy, g_zz_0_xzzzz_xxyz, g_zz_0_xzzzz_xxzz, g_zz_0_xzzzz_xyyy, g_zz_0_xzzzz_xyyz, g_zz_0_xzzzz_xyzz, g_zz_0_xzzzz_xzzz, g_zz_0_xzzzz_yyyy, g_zz_0_xzzzz_yyyz, g_zz_0_xzzzz_yyzz, g_zz_0_xzzzz_yzzz, g_zz_0_xzzzz_zzzz, g_zz_0_zzzz_xxxx, g_zz_0_zzzz_xxxxx, g_zz_0_zzzz_xxxxy, g_zz_0_zzzz_xxxxz, g_zz_0_zzzz_xxxy, g_zz_0_zzzz_xxxyy, g_zz_0_zzzz_xxxyz, g_zz_0_zzzz_xxxz, g_zz_0_zzzz_xxxzz, g_zz_0_zzzz_xxyy, g_zz_0_zzzz_xxyyy, g_zz_0_zzzz_xxyyz, g_zz_0_zzzz_xxyz, g_zz_0_zzzz_xxyzz, g_zz_0_zzzz_xxzz, g_zz_0_zzzz_xxzzz, g_zz_0_zzzz_xyyy, g_zz_0_zzzz_xyyyy, g_zz_0_zzzz_xyyyz, g_zz_0_zzzz_xyyz, g_zz_0_zzzz_xyyzz, g_zz_0_zzzz_xyzz, g_zz_0_zzzz_xyzzz, g_zz_0_zzzz_xzzz, g_zz_0_zzzz_xzzzz, g_zz_0_zzzz_yyyy, g_zz_0_zzzz_yyyz, g_zz_0_zzzz_yyzz, g_zz_0_zzzz_yzzz, g_zz_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzz_xxxx[k] = -g_zz_0_zzzz_xxxx[k] * ab_x + g_zz_0_zzzz_xxxxx[k];

                g_zz_0_xzzzz_xxxy[k] = -g_zz_0_zzzz_xxxy[k] * ab_x + g_zz_0_zzzz_xxxxy[k];

                g_zz_0_xzzzz_xxxz[k] = -g_zz_0_zzzz_xxxz[k] * ab_x + g_zz_0_zzzz_xxxxz[k];

                g_zz_0_xzzzz_xxyy[k] = -g_zz_0_zzzz_xxyy[k] * ab_x + g_zz_0_zzzz_xxxyy[k];

                g_zz_0_xzzzz_xxyz[k] = -g_zz_0_zzzz_xxyz[k] * ab_x + g_zz_0_zzzz_xxxyz[k];

                g_zz_0_xzzzz_xxzz[k] = -g_zz_0_zzzz_xxzz[k] * ab_x + g_zz_0_zzzz_xxxzz[k];

                g_zz_0_xzzzz_xyyy[k] = -g_zz_0_zzzz_xyyy[k] * ab_x + g_zz_0_zzzz_xxyyy[k];

                g_zz_0_xzzzz_xyyz[k] = -g_zz_0_zzzz_xyyz[k] * ab_x + g_zz_0_zzzz_xxyyz[k];

                g_zz_0_xzzzz_xyzz[k] = -g_zz_0_zzzz_xyzz[k] * ab_x + g_zz_0_zzzz_xxyzz[k];

                g_zz_0_xzzzz_xzzz[k] = -g_zz_0_zzzz_xzzz[k] * ab_x + g_zz_0_zzzz_xxzzz[k];

                g_zz_0_xzzzz_yyyy[k] = -g_zz_0_zzzz_yyyy[k] * ab_x + g_zz_0_zzzz_xyyyy[k];

                g_zz_0_xzzzz_yyyz[k] = -g_zz_0_zzzz_yyyz[k] * ab_x + g_zz_0_zzzz_xyyyz[k];

                g_zz_0_xzzzz_yyzz[k] = -g_zz_0_zzzz_yyzz[k] * ab_x + g_zz_0_zzzz_xyyzz[k];

                g_zz_0_xzzzz_yzzz[k] = -g_zz_0_zzzz_yzzz[k] * ab_x + g_zz_0_zzzz_xyzzz[k];

                g_zz_0_xzzzz_zzzz[k] = -g_zz_0_zzzz_zzzz[k] * ab_x + g_zz_0_zzzz_xzzzz[k];
            }

            /// Set up 1800-1815 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyy_xxxx = cbuffer.data(hg_geom_20_off + 1800 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxxy = cbuffer.data(hg_geom_20_off + 1801 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxxz = cbuffer.data(hg_geom_20_off + 1802 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxyy = cbuffer.data(hg_geom_20_off + 1803 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxyz = cbuffer.data(hg_geom_20_off + 1804 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xxzz = cbuffer.data(hg_geom_20_off + 1805 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xyyy = cbuffer.data(hg_geom_20_off + 1806 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xyyz = cbuffer.data(hg_geom_20_off + 1807 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xyzz = cbuffer.data(hg_geom_20_off + 1808 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xzzz = cbuffer.data(hg_geom_20_off + 1809 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yyyy = cbuffer.data(hg_geom_20_off + 1810 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yyyz = cbuffer.data(hg_geom_20_off + 1811 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yyzz = cbuffer.data(hg_geom_20_off + 1812 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yzzz = cbuffer.data(hg_geom_20_off + 1813 * ccomps * dcomps);

            auto g_zz_0_yyyyy_zzzz = cbuffer.data(hg_geom_20_off + 1814 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyy_xxxx, g_zz_0_yyyy_xxxxy, g_zz_0_yyyy_xxxy, g_zz_0_yyyy_xxxyy, g_zz_0_yyyy_xxxyz, g_zz_0_yyyy_xxxz, g_zz_0_yyyy_xxyy, g_zz_0_yyyy_xxyyy, g_zz_0_yyyy_xxyyz, g_zz_0_yyyy_xxyz, g_zz_0_yyyy_xxyzz, g_zz_0_yyyy_xxzz, g_zz_0_yyyy_xyyy, g_zz_0_yyyy_xyyyy, g_zz_0_yyyy_xyyyz, g_zz_0_yyyy_xyyz, g_zz_0_yyyy_xyyzz, g_zz_0_yyyy_xyzz, g_zz_0_yyyy_xyzzz, g_zz_0_yyyy_xzzz, g_zz_0_yyyy_yyyy, g_zz_0_yyyy_yyyyy, g_zz_0_yyyy_yyyyz, g_zz_0_yyyy_yyyz, g_zz_0_yyyy_yyyzz, g_zz_0_yyyy_yyzz, g_zz_0_yyyy_yyzzz, g_zz_0_yyyy_yzzz, g_zz_0_yyyy_yzzzz, g_zz_0_yyyy_zzzz, g_zz_0_yyyyy_xxxx, g_zz_0_yyyyy_xxxy, g_zz_0_yyyyy_xxxz, g_zz_0_yyyyy_xxyy, g_zz_0_yyyyy_xxyz, g_zz_0_yyyyy_xxzz, g_zz_0_yyyyy_xyyy, g_zz_0_yyyyy_xyyz, g_zz_0_yyyyy_xyzz, g_zz_0_yyyyy_xzzz, g_zz_0_yyyyy_yyyy, g_zz_0_yyyyy_yyyz, g_zz_0_yyyyy_yyzz, g_zz_0_yyyyy_yzzz, g_zz_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyy_xxxx[k] = -g_zz_0_yyyy_xxxx[k] * ab_y + g_zz_0_yyyy_xxxxy[k];

                g_zz_0_yyyyy_xxxy[k] = -g_zz_0_yyyy_xxxy[k] * ab_y + g_zz_0_yyyy_xxxyy[k];

                g_zz_0_yyyyy_xxxz[k] = -g_zz_0_yyyy_xxxz[k] * ab_y + g_zz_0_yyyy_xxxyz[k];

                g_zz_0_yyyyy_xxyy[k] = -g_zz_0_yyyy_xxyy[k] * ab_y + g_zz_0_yyyy_xxyyy[k];

                g_zz_0_yyyyy_xxyz[k] = -g_zz_0_yyyy_xxyz[k] * ab_y + g_zz_0_yyyy_xxyyz[k];

                g_zz_0_yyyyy_xxzz[k] = -g_zz_0_yyyy_xxzz[k] * ab_y + g_zz_0_yyyy_xxyzz[k];

                g_zz_0_yyyyy_xyyy[k] = -g_zz_0_yyyy_xyyy[k] * ab_y + g_zz_0_yyyy_xyyyy[k];

                g_zz_0_yyyyy_xyyz[k] = -g_zz_0_yyyy_xyyz[k] * ab_y + g_zz_0_yyyy_xyyyz[k];

                g_zz_0_yyyyy_xyzz[k] = -g_zz_0_yyyy_xyzz[k] * ab_y + g_zz_0_yyyy_xyyzz[k];

                g_zz_0_yyyyy_xzzz[k] = -g_zz_0_yyyy_xzzz[k] * ab_y + g_zz_0_yyyy_xyzzz[k];

                g_zz_0_yyyyy_yyyy[k] = -g_zz_0_yyyy_yyyy[k] * ab_y + g_zz_0_yyyy_yyyyy[k];

                g_zz_0_yyyyy_yyyz[k] = -g_zz_0_yyyy_yyyz[k] * ab_y + g_zz_0_yyyy_yyyyz[k];

                g_zz_0_yyyyy_yyzz[k] = -g_zz_0_yyyy_yyzz[k] * ab_y + g_zz_0_yyyy_yyyzz[k];

                g_zz_0_yyyyy_yzzz[k] = -g_zz_0_yyyy_yzzz[k] * ab_y + g_zz_0_yyyy_yyzzz[k];

                g_zz_0_yyyyy_zzzz[k] = -g_zz_0_yyyy_zzzz[k] * ab_y + g_zz_0_yyyy_yzzzz[k];
            }

            /// Set up 1815-1830 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyz_xxxx = cbuffer.data(hg_geom_20_off + 1815 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxxy = cbuffer.data(hg_geom_20_off + 1816 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxxz = cbuffer.data(hg_geom_20_off + 1817 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxyy = cbuffer.data(hg_geom_20_off + 1818 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxyz = cbuffer.data(hg_geom_20_off + 1819 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xxzz = cbuffer.data(hg_geom_20_off + 1820 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xyyy = cbuffer.data(hg_geom_20_off + 1821 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xyyz = cbuffer.data(hg_geom_20_off + 1822 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xyzz = cbuffer.data(hg_geom_20_off + 1823 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xzzz = cbuffer.data(hg_geom_20_off + 1824 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yyyy = cbuffer.data(hg_geom_20_off + 1825 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yyyz = cbuffer.data(hg_geom_20_off + 1826 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yyzz = cbuffer.data(hg_geom_20_off + 1827 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yzzz = cbuffer.data(hg_geom_20_off + 1828 * ccomps * dcomps);

            auto g_zz_0_yyyyz_zzzz = cbuffer.data(hg_geom_20_off + 1829 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyz_xxxx, g_zz_0_yyyyz_xxxy, g_zz_0_yyyyz_xxxz, g_zz_0_yyyyz_xxyy, g_zz_0_yyyyz_xxyz, g_zz_0_yyyyz_xxzz, g_zz_0_yyyyz_xyyy, g_zz_0_yyyyz_xyyz, g_zz_0_yyyyz_xyzz, g_zz_0_yyyyz_xzzz, g_zz_0_yyyyz_yyyy, g_zz_0_yyyyz_yyyz, g_zz_0_yyyyz_yyzz, g_zz_0_yyyyz_yzzz, g_zz_0_yyyyz_zzzz, g_zz_0_yyyz_xxxx, g_zz_0_yyyz_xxxxy, g_zz_0_yyyz_xxxy, g_zz_0_yyyz_xxxyy, g_zz_0_yyyz_xxxyz, g_zz_0_yyyz_xxxz, g_zz_0_yyyz_xxyy, g_zz_0_yyyz_xxyyy, g_zz_0_yyyz_xxyyz, g_zz_0_yyyz_xxyz, g_zz_0_yyyz_xxyzz, g_zz_0_yyyz_xxzz, g_zz_0_yyyz_xyyy, g_zz_0_yyyz_xyyyy, g_zz_0_yyyz_xyyyz, g_zz_0_yyyz_xyyz, g_zz_0_yyyz_xyyzz, g_zz_0_yyyz_xyzz, g_zz_0_yyyz_xyzzz, g_zz_0_yyyz_xzzz, g_zz_0_yyyz_yyyy, g_zz_0_yyyz_yyyyy, g_zz_0_yyyz_yyyyz, g_zz_0_yyyz_yyyz, g_zz_0_yyyz_yyyzz, g_zz_0_yyyz_yyzz, g_zz_0_yyyz_yyzzz, g_zz_0_yyyz_yzzz, g_zz_0_yyyz_yzzzz, g_zz_0_yyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyz_xxxx[k] = -g_zz_0_yyyz_xxxx[k] * ab_y + g_zz_0_yyyz_xxxxy[k];

                g_zz_0_yyyyz_xxxy[k] = -g_zz_0_yyyz_xxxy[k] * ab_y + g_zz_0_yyyz_xxxyy[k];

                g_zz_0_yyyyz_xxxz[k] = -g_zz_0_yyyz_xxxz[k] * ab_y + g_zz_0_yyyz_xxxyz[k];

                g_zz_0_yyyyz_xxyy[k] = -g_zz_0_yyyz_xxyy[k] * ab_y + g_zz_0_yyyz_xxyyy[k];

                g_zz_0_yyyyz_xxyz[k] = -g_zz_0_yyyz_xxyz[k] * ab_y + g_zz_0_yyyz_xxyyz[k];

                g_zz_0_yyyyz_xxzz[k] = -g_zz_0_yyyz_xxzz[k] * ab_y + g_zz_0_yyyz_xxyzz[k];

                g_zz_0_yyyyz_xyyy[k] = -g_zz_0_yyyz_xyyy[k] * ab_y + g_zz_0_yyyz_xyyyy[k];

                g_zz_0_yyyyz_xyyz[k] = -g_zz_0_yyyz_xyyz[k] * ab_y + g_zz_0_yyyz_xyyyz[k];

                g_zz_0_yyyyz_xyzz[k] = -g_zz_0_yyyz_xyzz[k] * ab_y + g_zz_0_yyyz_xyyzz[k];

                g_zz_0_yyyyz_xzzz[k] = -g_zz_0_yyyz_xzzz[k] * ab_y + g_zz_0_yyyz_xyzzz[k];

                g_zz_0_yyyyz_yyyy[k] = -g_zz_0_yyyz_yyyy[k] * ab_y + g_zz_0_yyyz_yyyyy[k];

                g_zz_0_yyyyz_yyyz[k] = -g_zz_0_yyyz_yyyz[k] * ab_y + g_zz_0_yyyz_yyyyz[k];

                g_zz_0_yyyyz_yyzz[k] = -g_zz_0_yyyz_yyzz[k] * ab_y + g_zz_0_yyyz_yyyzz[k];

                g_zz_0_yyyyz_yzzz[k] = -g_zz_0_yyyz_yzzz[k] * ab_y + g_zz_0_yyyz_yyzzz[k];

                g_zz_0_yyyyz_zzzz[k] = -g_zz_0_yyyz_zzzz[k] * ab_y + g_zz_0_yyyz_yzzzz[k];
            }

            /// Set up 1830-1845 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzz_xxxx = cbuffer.data(hg_geom_20_off + 1830 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxxy = cbuffer.data(hg_geom_20_off + 1831 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxxz = cbuffer.data(hg_geom_20_off + 1832 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxyy = cbuffer.data(hg_geom_20_off + 1833 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxyz = cbuffer.data(hg_geom_20_off + 1834 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xxzz = cbuffer.data(hg_geom_20_off + 1835 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xyyy = cbuffer.data(hg_geom_20_off + 1836 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xyyz = cbuffer.data(hg_geom_20_off + 1837 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xyzz = cbuffer.data(hg_geom_20_off + 1838 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xzzz = cbuffer.data(hg_geom_20_off + 1839 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yyyy = cbuffer.data(hg_geom_20_off + 1840 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yyyz = cbuffer.data(hg_geom_20_off + 1841 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yyzz = cbuffer.data(hg_geom_20_off + 1842 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yzzz = cbuffer.data(hg_geom_20_off + 1843 * ccomps * dcomps);

            auto g_zz_0_yyyzz_zzzz = cbuffer.data(hg_geom_20_off + 1844 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzz_xxxx, g_zz_0_yyyzz_xxxy, g_zz_0_yyyzz_xxxz, g_zz_0_yyyzz_xxyy, g_zz_0_yyyzz_xxyz, g_zz_0_yyyzz_xxzz, g_zz_0_yyyzz_xyyy, g_zz_0_yyyzz_xyyz, g_zz_0_yyyzz_xyzz, g_zz_0_yyyzz_xzzz, g_zz_0_yyyzz_yyyy, g_zz_0_yyyzz_yyyz, g_zz_0_yyyzz_yyzz, g_zz_0_yyyzz_yzzz, g_zz_0_yyyzz_zzzz, g_zz_0_yyzz_xxxx, g_zz_0_yyzz_xxxxy, g_zz_0_yyzz_xxxy, g_zz_0_yyzz_xxxyy, g_zz_0_yyzz_xxxyz, g_zz_0_yyzz_xxxz, g_zz_0_yyzz_xxyy, g_zz_0_yyzz_xxyyy, g_zz_0_yyzz_xxyyz, g_zz_0_yyzz_xxyz, g_zz_0_yyzz_xxyzz, g_zz_0_yyzz_xxzz, g_zz_0_yyzz_xyyy, g_zz_0_yyzz_xyyyy, g_zz_0_yyzz_xyyyz, g_zz_0_yyzz_xyyz, g_zz_0_yyzz_xyyzz, g_zz_0_yyzz_xyzz, g_zz_0_yyzz_xyzzz, g_zz_0_yyzz_xzzz, g_zz_0_yyzz_yyyy, g_zz_0_yyzz_yyyyy, g_zz_0_yyzz_yyyyz, g_zz_0_yyzz_yyyz, g_zz_0_yyzz_yyyzz, g_zz_0_yyzz_yyzz, g_zz_0_yyzz_yyzzz, g_zz_0_yyzz_yzzz, g_zz_0_yyzz_yzzzz, g_zz_0_yyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzz_xxxx[k] = -g_zz_0_yyzz_xxxx[k] * ab_y + g_zz_0_yyzz_xxxxy[k];

                g_zz_0_yyyzz_xxxy[k] = -g_zz_0_yyzz_xxxy[k] * ab_y + g_zz_0_yyzz_xxxyy[k];

                g_zz_0_yyyzz_xxxz[k] = -g_zz_0_yyzz_xxxz[k] * ab_y + g_zz_0_yyzz_xxxyz[k];

                g_zz_0_yyyzz_xxyy[k] = -g_zz_0_yyzz_xxyy[k] * ab_y + g_zz_0_yyzz_xxyyy[k];

                g_zz_0_yyyzz_xxyz[k] = -g_zz_0_yyzz_xxyz[k] * ab_y + g_zz_0_yyzz_xxyyz[k];

                g_zz_0_yyyzz_xxzz[k] = -g_zz_0_yyzz_xxzz[k] * ab_y + g_zz_0_yyzz_xxyzz[k];

                g_zz_0_yyyzz_xyyy[k] = -g_zz_0_yyzz_xyyy[k] * ab_y + g_zz_0_yyzz_xyyyy[k];

                g_zz_0_yyyzz_xyyz[k] = -g_zz_0_yyzz_xyyz[k] * ab_y + g_zz_0_yyzz_xyyyz[k];

                g_zz_0_yyyzz_xyzz[k] = -g_zz_0_yyzz_xyzz[k] * ab_y + g_zz_0_yyzz_xyyzz[k];

                g_zz_0_yyyzz_xzzz[k] = -g_zz_0_yyzz_xzzz[k] * ab_y + g_zz_0_yyzz_xyzzz[k];

                g_zz_0_yyyzz_yyyy[k] = -g_zz_0_yyzz_yyyy[k] * ab_y + g_zz_0_yyzz_yyyyy[k];

                g_zz_0_yyyzz_yyyz[k] = -g_zz_0_yyzz_yyyz[k] * ab_y + g_zz_0_yyzz_yyyyz[k];

                g_zz_0_yyyzz_yyzz[k] = -g_zz_0_yyzz_yyzz[k] * ab_y + g_zz_0_yyzz_yyyzz[k];

                g_zz_0_yyyzz_yzzz[k] = -g_zz_0_yyzz_yzzz[k] * ab_y + g_zz_0_yyzz_yyzzz[k];

                g_zz_0_yyyzz_zzzz[k] = -g_zz_0_yyzz_zzzz[k] * ab_y + g_zz_0_yyzz_yzzzz[k];
            }

            /// Set up 1845-1860 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzz_xxxx = cbuffer.data(hg_geom_20_off + 1845 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxxy = cbuffer.data(hg_geom_20_off + 1846 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxxz = cbuffer.data(hg_geom_20_off + 1847 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxyy = cbuffer.data(hg_geom_20_off + 1848 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxyz = cbuffer.data(hg_geom_20_off + 1849 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xxzz = cbuffer.data(hg_geom_20_off + 1850 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xyyy = cbuffer.data(hg_geom_20_off + 1851 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xyyz = cbuffer.data(hg_geom_20_off + 1852 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xyzz = cbuffer.data(hg_geom_20_off + 1853 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xzzz = cbuffer.data(hg_geom_20_off + 1854 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yyyy = cbuffer.data(hg_geom_20_off + 1855 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yyyz = cbuffer.data(hg_geom_20_off + 1856 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yyzz = cbuffer.data(hg_geom_20_off + 1857 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yzzz = cbuffer.data(hg_geom_20_off + 1858 * ccomps * dcomps);

            auto g_zz_0_yyzzz_zzzz = cbuffer.data(hg_geom_20_off + 1859 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzz_xxxx, g_zz_0_yyzzz_xxxy, g_zz_0_yyzzz_xxxz, g_zz_0_yyzzz_xxyy, g_zz_0_yyzzz_xxyz, g_zz_0_yyzzz_xxzz, g_zz_0_yyzzz_xyyy, g_zz_0_yyzzz_xyyz, g_zz_0_yyzzz_xyzz, g_zz_0_yyzzz_xzzz, g_zz_0_yyzzz_yyyy, g_zz_0_yyzzz_yyyz, g_zz_0_yyzzz_yyzz, g_zz_0_yyzzz_yzzz, g_zz_0_yyzzz_zzzz, g_zz_0_yzzz_xxxx, g_zz_0_yzzz_xxxxy, g_zz_0_yzzz_xxxy, g_zz_0_yzzz_xxxyy, g_zz_0_yzzz_xxxyz, g_zz_0_yzzz_xxxz, g_zz_0_yzzz_xxyy, g_zz_0_yzzz_xxyyy, g_zz_0_yzzz_xxyyz, g_zz_0_yzzz_xxyz, g_zz_0_yzzz_xxyzz, g_zz_0_yzzz_xxzz, g_zz_0_yzzz_xyyy, g_zz_0_yzzz_xyyyy, g_zz_0_yzzz_xyyyz, g_zz_0_yzzz_xyyz, g_zz_0_yzzz_xyyzz, g_zz_0_yzzz_xyzz, g_zz_0_yzzz_xyzzz, g_zz_0_yzzz_xzzz, g_zz_0_yzzz_yyyy, g_zz_0_yzzz_yyyyy, g_zz_0_yzzz_yyyyz, g_zz_0_yzzz_yyyz, g_zz_0_yzzz_yyyzz, g_zz_0_yzzz_yyzz, g_zz_0_yzzz_yyzzz, g_zz_0_yzzz_yzzz, g_zz_0_yzzz_yzzzz, g_zz_0_yzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzz_xxxx[k] = -g_zz_0_yzzz_xxxx[k] * ab_y + g_zz_0_yzzz_xxxxy[k];

                g_zz_0_yyzzz_xxxy[k] = -g_zz_0_yzzz_xxxy[k] * ab_y + g_zz_0_yzzz_xxxyy[k];

                g_zz_0_yyzzz_xxxz[k] = -g_zz_0_yzzz_xxxz[k] * ab_y + g_zz_0_yzzz_xxxyz[k];

                g_zz_0_yyzzz_xxyy[k] = -g_zz_0_yzzz_xxyy[k] * ab_y + g_zz_0_yzzz_xxyyy[k];

                g_zz_0_yyzzz_xxyz[k] = -g_zz_0_yzzz_xxyz[k] * ab_y + g_zz_0_yzzz_xxyyz[k];

                g_zz_0_yyzzz_xxzz[k] = -g_zz_0_yzzz_xxzz[k] * ab_y + g_zz_0_yzzz_xxyzz[k];

                g_zz_0_yyzzz_xyyy[k] = -g_zz_0_yzzz_xyyy[k] * ab_y + g_zz_0_yzzz_xyyyy[k];

                g_zz_0_yyzzz_xyyz[k] = -g_zz_0_yzzz_xyyz[k] * ab_y + g_zz_0_yzzz_xyyyz[k];

                g_zz_0_yyzzz_xyzz[k] = -g_zz_0_yzzz_xyzz[k] * ab_y + g_zz_0_yzzz_xyyzz[k];

                g_zz_0_yyzzz_xzzz[k] = -g_zz_0_yzzz_xzzz[k] * ab_y + g_zz_0_yzzz_xyzzz[k];

                g_zz_0_yyzzz_yyyy[k] = -g_zz_0_yzzz_yyyy[k] * ab_y + g_zz_0_yzzz_yyyyy[k];

                g_zz_0_yyzzz_yyyz[k] = -g_zz_0_yzzz_yyyz[k] * ab_y + g_zz_0_yzzz_yyyyz[k];

                g_zz_0_yyzzz_yyzz[k] = -g_zz_0_yzzz_yyzz[k] * ab_y + g_zz_0_yzzz_yyyzz[k];

                g_zz_0_yyzzz_yzzz[k] = -g_zz_0_yzzz_yzzz[k] * ab_y + g_zz_0_yzzz_yyzzz[k];

                g_zz_0_yyzzz_zzzz[k] = -g_zz_0_yzzz_zzzz[k] * ab_y + g_zz_0_yzzz_yzzzz[k];
            }

            /// Set up 1860-1875 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1860 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1861 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1862 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1863 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1864 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1865 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1866 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1867 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1868 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1869 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1870 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1871 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1872 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1873 * ccomps * dcomps);

            auto g_zz_0_yzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1874 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzz_xxxx, g_zz_0_yzzzz_xxxy, g_zz_0_yzzzz_xxxz, g_zz_0_yzzzz_xxyy, g_zz_0_yzzzz_xxyz, g_zz_0_yzzzz_xxzz, g_zz_0_yzzzz_xyyy, g_zz_0_yzzzz_xyyz, g_zz_0_yzzzz_xyzz, g_zz_0_yzzzz_xzzz, g_zz_0_yzzzz_yyyy, g_zz_0_yzzzz_yyyz, g_zz_0_yzzzz_yyzz, g_zz_0_yzzzz_yzzz, g_zz_0_yzzzz_zzzz, g_zz_0_zzzz_xxxx, g_zz_0_zzzz_xxxxy, g_zz_0_zzzz_xxxy, g_zz_0_zzzz_xxxyy, g_zz_0_zzzz_xxxyz, g_zz_0_zzzz_xxxz, g_zz_0_zzzz_xxyy, g_zz_0_zzzz_xxyyy, g_zz_0_zzzz_xxyyz, g_zz_0_zzzz_xxyz, g_zz_0_zzzz_xxyzz, g_zz_0_zzzz_xxzz, g_zz_0_zzzz_xyyy, g_zz_0_zzzz_xyyyy, g_zz_0_zzzz_xyyyz, g_zz_0_zzzz_xyyz, g_zz_0_zzzz_xyyzz, g_zz_0_zzzz_xyzz, g_zz_0_zzzz_xyzzz, g_zz_0_zzzz_xzzz, g_zz_0_zzzz_yyyy, g_zz_0_zzzz_yyyyy, g_zz_0_zzzz_yyyyz, g_zz_0_zzzz_yyyz, g_zz_0_zzzz_yyyzz, g_zz_0_zzzz_yyzz, g_zz_0_zzzz_yyzzz, g_zz_0_zzzz_yzzz, g_zz_0_zzzz_yzzzz, g_zz_0_zzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzz_xxxx[k] = -g_zz_0_zzzz_xxxx[k] * ab_y + g_zz_0_zzzz_xxxxy[k];

                g_zz_0_yzzzz_xxxy[k] = -g_zz_0_zzzz_xxxy[k] * ab_y + g_zz_0_zzzz_xxxyy[k];

                g_zz_0_yzzzz_xxxz[k] = -g_zz_0_zzzz_xxxz[k] * ab_y + g_zz_0_zzzz_xxxyz[k];

                g_zz_0_yzzzz_xxyy[k] = -g_zz_0_zzzz_xxyy[k] * ab_y + g_zz_0_zzzz_xxyyy[k];

                g_zz_0_yzzzz_xxyz[k] = -g_zz_0_zzzz_xxyz[k] * ab_y + g_zz_0_zzzz_xxyyz[k];

                g_zz_0_yzzzz_xxzz[k] = -g_zz_0_zzzz_xxzz[k] * ab_y + g_zz_0_zzzz_xxyzz[k];

                g_zz_0_yzzzz_xyyy[k] = -g_zz_0_zzzz_xyyy[k] * ab_y + g_zz_0_zzzz_xyyyy[k];

                g_zz_0_yzzzz_xyyz[k] = -g_zz_0_zzzz_xyyz[k] * ab_y + g_zz_0_zzzz_xyyyz[k];

                g_zz_0_yzzzz_xyzz[k] = -g_zz_0_zzzz_xyzz[k] * ab_y + g_zz_0_zzzz_xyyzz[k];

                g_zz_0_yzzzz_xzzz[k] = -g_zz_0_zzzz_xzzz[k] * ab_y + g_zz_0_zzzz_xyzzz[k];

                g_zz_0_yzzzz_yyyy[k] = -g_zz_0_zzzz_yyyy[k] * ab_y + g_zz_0_zzzz_yyyyy[k];

                g_zz_0_yzzzz_yyyz[k] = -g_zz_0_zzzz_yyyz[k] * ab_y + g_zz_0_zzzz_yyyyz[k];

                g_zz_0_yzzzz_yyzz[k] = -g_zz_0_zzzz_yyzz[k] * ab_y + g_zz_0_zzzz_yyyzz[k];

                g_zz_0_yzzzz_yzzz[k] = -g_zz_0_zzzz_yzzz[k] * ab_y + g_zz_0_zzzz_yyzzz[k];

                g_zz_0_yzzzz_zzzz[k] = -g_zz_0_zzzz_zzzz[k] * ab_y + g_zz_0_zzzz_yzzzz[k];
            }

            /// Set up 1875-1890 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzz_xxxx = cbuffer.data(hg_geom_20_off + 1875 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxxy = cbuffer.data(hg_geom_20_off + 1876 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxxz = cbuffer.data(hg_geom_20_off + 1877 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxyy = cbuffer.data(hg_geom_20_off + 1878 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxyz = cbuffer.data(hg_geom_20_off + 1879 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xxzz = cbuffer.data(hg_geom_20_off + 1880 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xyyy = cbuffer.data(hg_geom_20_off + 1881 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xyyz = cbuffer.data(hg_geom_20_off + 1882 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xyzz = cbuffer.data(hg_geom_20_off + 1883 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xzzz = cbuffer.data(hg_geom_20_off + 1884 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yyyy = cbuffer.data(hg_geom_20_off + 1885 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yyyz = cbuffer.data(hg_geom_20_off + 1886 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yyzz = cbuffer.data(hg_geom_20_off + 1887 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yzzz = cbuffer.data(hg_geom_20_off + 1888 * ccomps * dcomps);

            auto g_zz_0_zzzzz_zzzz = cbuffer.data(hg_geom_20_off + 1889 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz, g_zz_0_zzzz_xxxx, g_zz_0_zzzz_xxxxz, g_zz_0_zzzz_xxxy, g_zz_0_zzzz_xxxyz, g_zz_0_zzzz_xxxz, g_zz_0_zzzz_xxxzz, g_zz_0_zzzz_xxyy, g_zz_0_zzzz_xxyyz, g_zz_0_zzzz_xxyz, g_zz_0_zzzz_xxyzz, g_zz_0_zzzz_xxzz, g_zz_0_zzzz_xxzzz, g_zz_0_zzzz_xyyy, g_zz_0_zzzz_xyyyz, g_zz_0_zzzz_xyyz, g_zz_0_zzzz_xyyzz, g_zz_0_zzzz_xyzz, g_zz_0_zzzz_xyzzz, g_zz_0_zzzz_xzzz, g_zz_0_zzzz_xzzzz, g_zz_0_zzzz_yyyy, g_zz_0_zzzz_yyyyz, g_zz_0_zzzz_yyyz, g_zz_0_zzzz_yyyzz, g_zz_0_zzzz_yyzz, g_zz_0_zzzz_yyzzz, g_zz_0_zzzz_yzzz, g_zz_0_zzzz_yzzzz, g_zz_0_zzzz_zzzz, g_zz_0_zzzz_zzzzz, g_zz_0_zzzzz_xxxx, g_zz_0_zzzzz_xxxy, g_zz_0_zzzzz_xxxz, g_zz_0_zzzzz_xxyy, g_zz_0_zzzzz_xxyz, g_zz_0_zzzzz_xxzz, g_zz_0_zzzzz_xyyy, g_zz_0_zzzzz_xyyz, g_zz_0_zzzzz_xyzz, g_zz_0_zzzzz_xzzz, g_zz_0_zzzzz_yyyy, g_zz_0_zzzzz_yyyz, g_zz_0_zzzzz_yyzz, g_zz_0_zzzzz_yzzz, g_zz_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzz_xxxx[k] = -2.0 * g_z_0_zzzz_xxxx[k] - g_zz_0_zzzz_xxxx[k] * ab_z + g_zz_0_zzzz_xxxxz[k];

                g_zz_0_zzzzz_xxxy[k] = -2.0 * g_z_0_zzzz_xxxy[k] - g_zz_0_zzzz_xxxy[k] * ab_z + g_zz_0_zzzz_xxxyz[k];

                g_zz_0_zzzzz_xxxz[k] = -2.0 * g_z_0_zzzz_xxxz[k] - g_zz_0_zzzz_xxxz[k] * ab_z + g_zz_0_zzzz_xxxzz[k];

                g_zz_0_zzzzz_xxyy[k] = -2.0 * g_z_0_zzzz_xxyy[k] - g_zz_0_zzzz_xxyy[k] * ab_z + g_zz_0_zzzz_xxyyz[k];

                g_zz_0_zzzzz_xxyz[k] = -2.0 * g_z_0_zzzz_xxyz[k] - g_zz_0_zzzz_xxyz[k] * ab_z + g_zz_0_zzzz_xxyzz[k];

                g_zz_0_zzzzz_xxzz[k] = -2.0 * g_z_0_zzzz_xxzz[k] - g_zz_0_zzzz_xxzz[k] * ab_z + g_zz_0_zzzz_xxzzz[k];

                g_zz_0_zzzzz_xyyy[k] = -2.0 * g_z_0_zzzz_xyyy[k] - g_zz_0_zzzz_xyyy[k] * ab_z + g_zz_0_zzzz_xyyyz[k];

                g_zz_0_zzzzz_xyyz[k] = -2.0 * g_z_0_zzzz_xyyz[k] - g_zz_0_zzzz_xyyz[k] * ab_z + g_zz_0_zzzz_xyyzz[k];

                g_zz_0_zzzzz_xyzz[k] = -2.0 * g_z_0_zzzz_xyzz[k] - g_zz_0_zzzz_xyzz[k] * ab_z + g_zz_0_zzzz_xyzzz[k];

                g_zz_0_zzzzz_xzzz[k] = -2.0 * g_z_0_zzzz_xzzz[k] - g_zz_0_zzzz_xzzz[k] * ab_z + g_zz_0_zzzz_xzzzz[k];

                g_zz_0_zzzzz_yyyy[k] = -2.0 * g_z_0_zzzz_yyyy[k] - g_zz_0_zzzz_yyyy[k] * ab_z + g_zz_0_zzzz_yyyyz[k];

                g_zz_0_zzzzz_yyyz[k] = -2.0 * g_z_0_zzzz_yyyz[k] - g_zz_0_zzzz_yyyz[k] * ab_z + g_zz_0_zzzz_yyyzz[k];

                g_zz_0_zzzzz_yyzz[k] = -2.0 * g_z_0_zzzz_yyzz[k] - g_zz_0_zzzz_yyzz[k] * ab_z + g_zz_0_zzzz_yyzzz[k];

                g_zz_0_zzzzz_yzzz[k] = -2.0 * g_z_0_zzzz_yzzz[k] - g_zz_0_zzzz_yzzz[k] * ab_z + g_zz_0_zzzz_yzzzz[k];

                g_zz_0_zzzzz_zzzz[k] = -2.0 * g_z_0_zzzz_zzzz[k] - g_zz_0_zzzz_zzzz[k] * ab_z + g_zz_0_zzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

