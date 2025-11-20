#include "ElectronRepulsionGeom2000ContrRecFGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_fgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_fgxx,
                                            const size_t idx_geom_10_dgxx,
                                            const size_t idx_geom_20_dgxx,
                                            const size_t idx_geom_20_dhxx,
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
            /// Set up components of auxilary buffer : DGSS

            const auto dg_geom_10_off = idx_geom_10_dgxx + i * dcomps + j;

            auto g_x_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 188 * ccomps * dcomps);

            auto g_z_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 189 * ccomps * dcomps);

            auto g_z_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 190 * ccomps * dcomps);

            auto g_z_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 191 * ccomps * dcomps);

            auto g_z_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 192 * ccomps * dcomps);

            auto g_z_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 193 * ccomps * dcomps);

            auto g_z_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 194 * ccomps * dcomps);

            auto g_z_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 195 * ccomps * dcomps);

            auto g_z_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 196 * ccomps * dcomps);

            auto g_z_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 197 * ccomps * dcomps);

            auto g_z_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 198 * ccomps * dcomps);

            auto g_z_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 199 * ccomps * dcomps);

            auto g_z_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 200 * ccomps * dcomps);

            auto g_z_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 201 * ccomps * dcomps);

            auto g_z_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 202 * ccomps * dcomps);

            auto g_z_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 203 * ccomps * dcomps);

            auto g_z_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 204 * ccomps * dcomps);

            auto g_z_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 205 * ccomps * dcomps);

            auto g_z_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 206 * ccomps * dcomps);

            auto g_z_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 207 * ccomps * dcomps);

            auto g_z_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 208 * ccomps * dcomps);

            auto g_z_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 209 * ccomps * dcomps);

            auto g_z_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 210 * ccomps * dcomps);

            auto g_z_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 211 * ccomps * dcomps);

            auto g_z_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 212 * ccomps * dcomps);

            auto g_z_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 213 * ccomps * dcomps);

            auto g_z_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 214 * ccomps * dcomps);

            auto g_z_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 224 * ccomps * dcomps);

            auto g_z_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 239 * ccomps * dcomps);

            auto g_z_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 244 * ccomps * dcomps);

            auto g_z_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 245 * ccomps * dcomps);

            auto g_z_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 246 * ccomps * dcomps);

            auto g_z_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 247 * ccomps * dcomps);

            auto g_z_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 248 * ccomps * dcomps);

            auto g_z_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 249 * ccomps * dcomps);

            auto g_z_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 250 * ccomps * dcomps);

            auto g_z_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DGSS

            const auto dg_geom_20_off = idx_geom_20_dgxx + i * dcomps + j;

            auto g_xx_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 194 * ccomps * dcomps);

            auto g_xz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 269 * ccomps * dcomps);

            auto g_yy_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 275 * ccomps * dcomps);

            auto g_yy_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 281 * ccomps * dcomps);

            auto g_yy_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 284 * ccomps * dcomps);

            auto g_yy_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 287 * ccomps * dcomps);

            auto g_yy_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 293 * ccomps * dcomps);

            auto g_yy_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 359 * ccomps * dcomps);

            auto g_yz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 365 * ccomps * dcomps);

            auto g_yz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 371 * ccomps * dcomps);

            auto g_yz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 374 * ccomps * dcomps);

            auto g_yz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 377 * ccomps * dcomps);

            auto g_yz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 383 * ccomps * dcomps);

            auto g_yz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 389 * ccomps * dcomps);

            auto g_yz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 395 * ccomps * dcomps);

            auto g_yz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 419 * ccomps * dcomps);

            auto g_yz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 449 * ccomps * dcomps);

            auto g_zz_0_xx_xxxx = cbuffer.data(dg_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xx_xxxy = cbuffer.data(dg_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xx_xxxz = cbuffer.data(dg_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_xx_xxyy = cbuffer.data(dg_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xx_xxyz = cbuffer.data(dg_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xx_xxzz = cbuffer.data(dg_geom_20_off + 455 * ccomps * dcomps);

            auto g_zz_0_xx_xyyy = cbuffer.data(dg_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xx_xyyz = cbuffer.data(dg_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xx_xyzz = cbuffer.data(dg_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_xx_xzzz = cbuffer.data(dg_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xx_yyyy = cbuffer.data(dg_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xx_yyyz = cbuffer.data(dg_geom_20_off + 461 * ccomps * dcomps);

            auto g_zz_0_xx_yyzz = cbuffer.data(dg_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xx_yzzz = cbuffer.data(dg_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xx_zzzz = cbuffer.data(dg_geom_20_off + 464 * ccomps * dcomps);

            auto g_zz_0_xy_xxxx = cbuffer.data(dg_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xy_xxxy = cbuffer.data(dg_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xy_xxxz = cbuffer.data(dg_geom_20_off + 467 * ccomps * dcomps);

            auto g_zz_0_xy_xxyy = cbuffer.data(dg_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xy_xxyz = cbuffer.data(dg_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xy_xxzz = cbuffer.data(dg_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_xy_xyyy = cbuffer.data(dg_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xy_xyyz = cbuffer.data(dg_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xy_xyzz = cbuffer.data(dg_geom_20_off + 473 * ccomps * dcomps);

            auto g_zz_0_xy_xzzz = cbuffer.data(dg_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xy_yyyy = cbuffer.data(dg_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xy_yyyz = cbuffer.data(dg_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_xy_yyzz = cbuffer.data(dg_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xy_yzzz = cbuffer.data(dg_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xy_zzzz = cbuffer.data(dg_geom_20_off + 479 * ccomps * dcomps);

            auto g_zz_0_xz_xxxx = cbuffer.data(dg_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xz_xxxy = cbuffer.data(dg_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xz_xxxz = cbuffer.data(dg_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_xz_xxyy = cbuffer.data(dg_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_xz_xxyz = cbuffer.data(dg_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_xz_xxzz = cbuffer.data(dg_geom_20_off + 485 * ccomps * dcomps);

            auto g_zz_0_xz_xyyy = cbuffer.data(dg_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_xz_xyyz = cbuffer.data(dg_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_xz_xyzz = cbuffer.data(dg_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_xz_xzzz = cbuffer.data(dg_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_xz_yyyy = cbuffer.data(dg_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_xz_yyyz = cbuffer.data(dg_geom_20_off + 491 * ccomps * dcomps);

            auto g_zz_0_xz_yyzz = cbuffer.data(dg_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_xz_yzzz = cbuffer.data(dg_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_xz_zzzz = cbuffer.data(dg_geom_20_off + 494 * ccomps * dcomps);

            auto g_zz_0_yy_xxxx = cbuffer.data(dg_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_yy_xxxy = cbuffer.data(dg_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_yy_xxxz = cbuffer.data(dg_geom_20_off + 497 * ccomps * dcomps);

            auto g_zz_0_yy_xxyy = cbuffer.data(dg_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_yy_xxyz = cbuffer.data(dg_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_yy_xxzz = cbuffer.data(dg_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_yy_xyyy = cbuffer.data(dg_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_yy_xyyz = cbuffer.data(dg_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_yy_xyzz = cbuffer.data(dg_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_yy_xzzz = cbuffer.data(dg_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_yy_yyyy = cbuffer.data(dg_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_yy_yyyz = cbuffer.data(dg_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_yy_yyzz = cbuffer.data(dg_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_yy_yzzz = cbuffer.data(dg_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_yy_zzzz = cbuffer.data(dg_geom_20_off + 509 * ccomps * dcomps);

            auto g_zz_0_yz_xxxx = cbuffer.data(dg_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_yz_xxxy = cbuffer.data(dg_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_yz_xxxz = cbuffer.data(dg_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_yz_xxyy = cbuffer.data(dg_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_yz_xxyz = cbuffer.data(dg_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_yz_xxzz = cbuffer.data(dg_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_yz_xyyy = cbuffer.data(dg_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_yz_xyyz = cbuffer.data(dg_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_yz_xyzz = cbuffer.data(dg_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_yz_xzzz = cbuffer.data(dg_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_yz_yyyy = cbuffer.data(dg_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_yz_yyyz = cbuffer.data(dg_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_yz_yyzz = cbuffer.data(dg_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_yz_yzzz = cbuffer.data(dg_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_yz_zzzz = cbuffer.data(dg_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_zz_xxxx = cbuffer.data(dg_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_zz_xxxy = cbuffer.data(dg_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_zz_xxxz = cbuffer.data(dg_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_zz_xxyy = cbuffer.data(dg_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_zz_xxyz = cbuffer.data(dg_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_zz_xxzz = cbuffer.data(dg_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_zz_xyyy = cbuffer.data(dg_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_zz_xyyz = cbuffer.data(dg_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_zz_xyzz = cbuffer.data(dg_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_zz_xzzz = cbuffer.data(dg_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_zz_yyyy = cbuffer.data(dg_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_zz_yyyz = cbuffer.data(dg_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_zz_yyzz = cbuffer.data(dg_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_zz_yzzz = cbuffer.data(dg_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_zz_zzzz = cbuffer.data(dg_geom_20_off + 539 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_20_off = idx_geom_20_dhxx + i * dcomps + j;

            auto g_xx_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 401 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 407 * ccomps * dcomps);

            auto g_yy_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 413 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 419 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 425 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 431 * ccomps * dcomps);

            auto g_yy_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 432 * ccomps * dcomps);

            auto g_yy_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 433 * ccomps * dcomps);

            auto g_yy_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 434 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 435 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 436 * ccomps * dcomps);

            auto g_yy_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 437 * ccomps * dcomps);

            auto g_yy_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 438 * ccomps * dcomps);

            auto g_yy_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 439 * ccomps * dcomps);

            auto g_yy_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 440 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 441 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 442 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 443 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 444 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 445 * ccomps * dcomps);

            auto g_yy_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 446 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 447 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 448 * ccomps * dcomps);

            auto g_yy_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 503 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 509 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 515 * ccomps * dcomps);

            auto g_yz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 521 * ccomps * dcomps);

            auto g_yz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 524 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 527 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 533 * ccomps * dcomps);

            auto g_yz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 539 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 540 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 541 * ccomps * dcomps);

            auto g_yz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 542 * ccomps * dcomps);

            auto g_yz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 543 * ccomps * dcomps);

            auto g_yz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 544 * ccomps * dcomps);

            auto g_yz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 545 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 546 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 547 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 548 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 549 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 550 * ccomps * dcomps);

            auto g_yz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 551 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 552 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 553 * ccomps * dcomps);

            auto g_yz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 554 * ccomps * dcomps);

            auto g_yz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 555 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 556 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 557 * ccomps * dcomps);

            auto g_yz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 558 * ccomps * dcomps);

            auto g_yz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 559 * ccomps * dcomps);

            auto g_yz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 560 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 561 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 562 * ccomps * dcomps);

            auto g_yz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 563 * ccomps * dcomps);

            auto g_yz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 564 * ccomps * dcomps);

            auto g_yz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 565 * ccomps * dcomps);

            auto g_yz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 566 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 567 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 568 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 569 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 570 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 571 * ccomps * dcomps);

            auto g_yz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 572 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 573 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 574 * ccomps * dcomps);

            auto g_yz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 575 * ccomps * dcomps);

            auto g_yz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 576 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 577 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 578 * ccomps * dcomps);

            auto g_yz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 579 * ccomps * dcomps);

            auto g_yz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 580 * ccomps * dcomps);

            auto g_yz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 581 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 582 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 583 * ccomps * dcomps);

            auto g_yz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 584 * ccomps * dcomps);

            auto g_yz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 585 * ccomps * dcomps);

            auto g_yz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 586 * ccomps * dcomps);

            auto g_yz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 587 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 588 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 589 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 590 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 591 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 592 * ccomps * dcomps);

            auto g_yz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 593 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 594 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 595 * ccomps * dcomps);

            auto g_yz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 596 * ccomps * dcomps);

            auto g_yz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 597 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 598 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 629 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 635 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 641 * ccomps * dcomps);

            auto g_zz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 647 * ccomps * dcomps);

            auto g_zz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 648 * ccomps * dcomps);

            auto g_zz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 649 * ccomps * dcomps);

            auto g_zz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 650 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 651 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 652 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 653 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 654 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 655 * ccomps * dcomps);

            auto g_zz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 656 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 657 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 658 * ccomps * dcomps);

            auto g_zz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 659 * ccomps * dcomps);

            auto g_zz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 660 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 661 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 662 * ccomps * dcomps);

            auto g_zz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 663 * ccomps * dcomps);

            auto g_zz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 664 * ccomps * dcomps);

            auto g_zz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 665 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 666 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 667 * ccomps * dcomps);

            auto g_zz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 668 * ccomps * dcomps);

            auto g_zz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 669 * ccomps * dcomps);

            auto g_zz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 670 * ccomps * dcomps);

            auto g_zz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 671 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 672 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 673 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 674 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 675 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 676 * ccomps * dcomps);

            auto g_zz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 677 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 678 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 679 * ccomps * dcomps);

            auto g_zz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 680 * ccomps * dcomps);

            auto g_zz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 681 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 682 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 683 * ccomps * dcomps);

            auto g_zz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 684 * ccomps * dcomps);

            auto g_zz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 685 * ccomps * dcomps);

            auto g_zz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 686 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 687 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 688 * ccomps * dcomps);

            auto g_zz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 689 * ccomps * dcomps);

            auto g_zz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 690 * ccomps * dcomps);

            auto g_zz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 691 * ccomps * dcomps);

            auto g_zz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 692 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 693 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 694 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 695 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 696 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 697 * ccomps * dcomps);

            auto g_zz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 698 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 699 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 700 * ccomps * dcomps);

            auto g_zz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 701 * ccomps * dcomps);

            auto g_zz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 702 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 703 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 704 * ccomps * dcomps);

            auto g_zz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 705 * ccomps * dcomps);

            auto g_zz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 706 * ccomps * dcomps);

            auto g_zz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 707 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 708 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 709 * ccomps * dcomps);

            auto g_zz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 710 * ccomps * dcomps);

            auto g_zz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 711 * ccomps * dcomps);

            auto g_zz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 712 * ccomps * dcomps);

            auto g_zz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 713 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 714 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 715 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 716 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 717 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 718 * ccomps * dcomps);

            auto g_zz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 719 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 720 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 721 * ccomps * dcomps);

            auto g_zz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 722 * ccomps * dcomps);

            auto g_zz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 723 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 724 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 725 * ccomps * dcomps);

            auto g_zz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 726 * ccomps * dcomps);

            auto g_zz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 727 * ccomps * dcomps);

            auto g_zz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 728 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 729 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 730 * ccomps * dcomps);

            auto g_zz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 731 * ccomps * dcomps);

            auto g_zz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 732 * ccomps * dcomps);

            auto g_zz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 733 * ccomps * dcomps);

            auto g_zz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 734 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 735 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 736 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 737 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 738 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 739 * ccomps * dcomps);

            auto g_zz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 740 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 741 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 742 * ccomps * dcomps);

            auto g_zz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 743 * ccomps * dcomps);

            auto g_zz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 744 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 745 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 746 * ccomps * dcomps);

            auto g_zz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 747 * ccomps * dcomps);

            auto g_zz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 748 * ccomps * dcomps);

            auto g_zz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 755 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fgxx

            const auto fg_geom_20_off = idx_geom_20_fgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxx, g_x_0_xx_xxxy, g_x_0_xx_xxxz, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxzz, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyzz, g_x_0_xx_xzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyzz, g_x_0_xx_yzzz, g_x_0_xx_zzzz, g_xx_0_xx_xxxx, g_xx_0_xx_xxxxx, g_xx_0_xx_xxxxy, g_xx_0_xx_xxxxz, g_xx_0_xx_xxxy, g_xx_0_xx_xxxyy, g_xx_0_xx_xxxyz, g_xx_0_xx_xxxz, g_xx_0_xx_xxxzz, g_xx_0_xx_xxyy, g_xx_0_xx_xxyyy, g_xx_0_xx_xxyyz, g_xx_0_xx_xxyz, g_xx_0_xx_xxyzz, g_xx_0_xx_xxzz, g_xx_0_xx_xxzzz, g_xx_0_xx_xyyy, g_xx_0_xx_xyyyy, g_xx_0_xx_xyyyz, g_xx_0_xx_xyyz, g_xx_0_xx_xyyzz, g_xx_0_xx_xyzz, g_xx_0_xx_xyzzz, g_xx_0_xx_xzzz, g_xx_0_xx_xzzzz, g_xx_0_xx_yyyy, g_xx_0_xx_yyyz, g_xx_0_xx_yyzz, g_xx_0_xx_yzzz, g_xx_0_xx_zzzz, g_xx_0_xxx_xxxx, g_xx_0_xxx_xxxy, g_xx_0_xxx_xxxz, g_xx_0_xxx_xxyy, g_xx_0_xxx_xxyz, g_xx_0_xxx_xxzz, g_xx_0_xxx_xyyy, g_xx_0_xxx_xyyz, g_xx_0_xxx_xyzz, g_xx_0_xxx_xzzz, g_xx_0_xxx_yyyy, g_xx_0_xxx_yyyz, g_xx_0_xxx_yyzz, g_xx_0_xxx_yzzz, g_xx_0_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxx_xxxx[k] = -2.0 * g_x_0_xx_xxxx[k] - g_xx_0_xx_xxxx[k] * ab_x + g_xx_0_xx_xxxxx[k];

                g_xx_0_xxx_xxxy[k] = -2.0 * g_x_0_xx_xxxy[k] - g_xx_0_xx_xxxy[k] * ab_x + g_xx_0_xx_xxxxy[k];

                g_xx_0_xxx_xxxz[k] = -2.0 * g_x_0_xx_xxxz[k] - g_xx_0_xx_xxxz[k] * ab_x + g_xx_0_xx_xxxxz[k];

                g_xx_0_xxx_xxyy[k] = -2.0 * g_x_0_xx_xxyy[k] - g_xx_0_xx_xxyy[k] * ab_x + g_xx_0_xx_xxxyy[k];

                g_xx_0_xxx_xxyz[k] = -2.0 * g_x_0_xx_xxyz[k] - g_xx_0_xx_xxyz[k] * ab_x + g_xx_0_xx_xxxyz[k];

                g_xx_0_xxx_xxzz[k] = -2.0 * g_x_0_xx_xxzz[k] - g_xx_0_xx_xxzz[k] * ab_x + g_xx_0_xx_xxxzz[k];

                g_xx_0_xxx_xyyy[k] = -2.0 * g_x_0_xx_xyyy[k] - g_xx_0_xx_xyyy[k] * ab_x + g_xx_0_xx_xxyyy[k];

                g_xx_0_xxx_xyyz[k] = -2.0 * g_x_0_xx_xyyz[k] - g_xx_0_xx_xyyz[k] * ab_x + g_xx_0_xx_xxyyz[k];

                g_xx_0_xxx_xyzz[k] = -2.0 * g_x_0_xx_xyzz[k] - g_xx_0_xx_xyzz[k] * ab_x + g_xx_0_xx_xxyzz[k];

                g_xx_0_xxx_xzzz[k] = -2.0 * g_x_0_xx_xzzz[k] - g_xx_0_xx_xzzz[k] * ab_x + g_xx_0_xx_xxzzz[k];

                g_xx_0_xxx_yyyy[k] = -2.0 * g_x_0_xx_yyyy[k] - g_xx_0_xx_yyyy[k] * ab_x + g_xx_0_xx_xyyyy[k];

                g_xx_0_xxx_yyyz[k] = -2.0 * g_x_0_xx_yyyz[k] - g_xx_0_xx_yyyz[k] * ab_x + g_xx_0_xx_xyyyz[k];

                g_xx_0_xxx_yyzz[k] = -2.0 * g_x_0_xx_yyzz[k] - g_xx_0_xx_yyzz[k] * ab_x + g_xx_0_xx_xyyzz[k];

                g_xx_0_xxx_yzzz[k] = -2.0 * g_x_0_xx_yzzz[k] - g_xx_0_xx_yzzz[k] * ab_x + g_xx_0_xx_xyzzz[k];

                g_xx_0_xxx_zzzz[k] = -2.0 * g_x_0_xx_zzzz[k] - g_xx_0_xx_zzzz[k] * ab_x + g_xx_0_xx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_xxxx, g_xx_0_xx_xxxxy, g_xx_0_xx_xxxy, g_xx_0_xx_xxxyy, g_xx_0_xx_xxxyz, g_xx_0_xx_xxxz, g_xx_0_xx_xxyy, g_xx_0_xx_xxyyy, g_xx_0_xx_xxyyz, g_xx_0_xx_xxyz, g_xx_0_xx_xxyzz, g_xx_0_xx_xxzz, g_xx_0_xx_xyyy, g_xx_0_xx_xyyyy, g_xx_0_xx_xyyyz, g_xx_0_xx_xyyz, g_xx_0_xx_xyyzz, g_xx_0_xx_xyzz, g_xx_0_xx_xyzzz, g_xx_0_xx_xzzz, g_xx_0_xx_yyyy, g_xx_0_xx_yyyyy, g_xx_0_xx_yyyyz, g_xx_0_xx_yyyz, g_xx_0_xx_yyyzz, g_xx_0_xx_yyzz, g_xx_0_xx_yyzzz, g_xx_0_xx_yzzz, g_xx_0_xx_yzzzz, g_xx_0_xx_zzzz, g_xx_0_xxy_xxxx, g_xx_0_xxy_xxxy, g_xx_0_xxy_xxxz, g_xx_0_xxy_xxyy, g_xx_0_xxy_xxyz, g_xx_0_xxy_xxzz, g_xx_0_xxy_xyyy, g_xx_0_xxy_xyyz, g_xx_0_xxy_xyzz, g_xx_0_xxy_xzzz, g_xx_0_xxy_yyyy, g_xx_0_xxy_yyyz, g_xx_0_xxy_yyzz, g_xx_0_xxy_yzzz, g_xx_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxy_xxxx[k] = -g_xx_0_xx_xxxx[k] * ab_y + g_xx_0_xx_xxxxy[k];

                g_xx_0_xxy_xxxy[k] = -g_xx_0_xx_xxxy[k] * ab_y + g_xx_0_xx_xxxyy[k];

                g_xx_0_xxy_xxxz[k] = -g_xx_0_xx_xxxz[k] * ab_y + g_xx_0_xx_xxxyz[k];

                g_xx_0_xxy_xxyy[k] = -g_xx_0_xx_xxyy[k] * ab_y + g_xx_0_xx_xxyyy[k];

                g_xx_0_xxy_xxyz[k] = -g_xx_0_xx_xxyz[k] * ab_y + g_xx_0_xx_xxyyz[k];

                g_xx_0_xxy_xxzz[k] = -g_xx_0_xx_xxzz[k] * ab_y + g_xx_0_xx_xxyzz[k];

                g_xx_0_xxy_xyyy[k] = -g_xx_0_xx_xyyy[k] * ab_y + g_xx_0_xx_xyyyy[k];

                g_xx_0_xxy_xyyz[k] = -g_xx_0_xx_xyyz[k] * ab_y + g_xx_0_xx_xyyyz[k];

                g_xx_0_xxy_xyzz[k] = -g_xx_0_xx_xyzz[k] * ab_y + g_xx_0_xx_xyyzz[k];

                g_xx_0_xxy_xzzz[k] = -g_xx_0_xx_xzzz[k] * ab_y + g_xx_0_xx_xyzzz[k];

                g_xx_0_xxy_yyyy[k] = -g_xx_0_xx_yyyy[k] * ab_y + g_xx_0_xx_yyyyy[k];

                g_xx_0_xxy_yyyz[k] = -g_xx_0_xx_yyyz[k] * ab_y + g_xx_0_xx_yyyyz[k];

                g_xx_0_xxy_yyzz[k] = -g_xx_0_xx_yyzz[k] * ab_y + g_xx_0_xx_yyyzz[k];

                g_xx_0_xxy_yzzz[k] = -g_xx_0_xx_yzzz[k] * ab_y + g_xx_0_xx_yyzzz[k];

                g_xx_0_xxy_zzzz[k] = -g_xx_0_xx_zzzz[k] * ab_y + g_xx_0_xx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_xxxx, g_xx_0_xx_xxxxz, g_xx_0_xx_xxxy, g_xx_0_xx_xxxyz, g_xx_0_xx_xxxz, g_xx_0_xx_xxxzz, g_xx_0_xx_xxyy, g_xx_0_xx_xxyyz, g_xx_0_xx_xxyz, g_xx_0_xx_xxyzz, g_xx_0_xx_xxzz, g_xx_0_xx_xxzzz, g_xx_0_xx_xyyy, g_xx_0_xx_xyyyz, g_xx_0_xx_xyyz, g_xx_0_xx_xyyzz, g_xx_0_xx_xyzz, g_xx_0_xx_xyzzz, g_xx_0_xx_xzzz, g_xx_0_xx_xzzzz, g_xx_0_xx_yyyy, g_xx_0_xx_yyyyz, g_xx_0_xx_yyyz, g_xx_0_xx_yyyzz, g_xx_0_xx_yyzz, g_xx_0_xx_yyzzz, g_xx_0_xx_yzzz, g_xx_0_xx_yzzzz, g_xx_0_xx_zzzz, g_xx_0_xx_zzzzz, g_xx_0_xxz_xxxx, g_xx_0_xxz_xxxy, g_xx_0_xxz_xxxz, g_xx_0_xxz_xxyy, g_xx_0_xxz_xxyz, g_xx_0_xxz_xxzz, g_xx_0_xxz_xyyy, g_xx_0_xxz_xyyz, g_xx_0_xxz_xyzz, g_xx_0_xxz_xzzz, g_xx_0_xxz_yyyy, g_xx_0_xxz_yyyz, g_xx_0_xxz_yyzz, g_xx_0_xxz_yzzz, g_xx_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxz_xxxx[k] = -g_xx_0_xx_xxxx[k] * ab_z + g_xx_0_xx_xxxxz[k];

                g_xx_0_xxz_xxxy[k] = -g_xx_0_xx_xxxy[k] * ab_z + g_xx_0_xx_xxxyz[k];

                g_xx_0_xxz_xxxz[k] = -g_xx_0_xx_xxxz[k] * ab_z + g_xx_0_xx_xxxzz[k];

                g_xx_0_xxz_xxyy[k] = -g_xx_0_xx_xxyy[k] * ab_z + g_xx_0_xx_xxyyz[k];

                g_xx_0_xxz_xxyz[k] = -g_xx_0_xx_xxyz[k] * ab_z + g_xx_0_xx_xxyzz[k];

                g_xx_0_xxz_xxzz[k] = -g_xx_0_xx_xxzz[k] * ab_z + g_xx_0_xx_xxzzz[k];

                g_xx_0_xxz_xyyy[k] = -g_xx_0_xx_xyyy[k] * ab_z + g_xx_0_xx_xyyyz[k];

                g_xx_0_xxz_xyyz[k] = -g_xx_0_xx_xyyz[k] * ab_z + g_xx_0_xx_xyyzz[k];

                g_xx_0_xxz_xyzz[k] = -g_xx_0_xx_xyzz[k] * ab_z + g_xx_0_xx_xyzzz[k];

                g_xx_0_xxz_xzzz[k] = -g_xx_0_xx_xzzz[k] * ab_z + g_xx_0_xx_xzzzz[k];

                g_xx_0_xxz_yyyy[k] = -g_xx_0_xx_yyyy[k] * ab_z + g_xx_0_xx_yyyyz[k];

                g_xx_0_xxz_yyyz[k] = -g_xx_0_xx_yyyz[k] * ab_z + g_xx_0_xx_yyyzz[k];

                g_xx_0_xxz_yyzz[k] = -g_xx_0_xx_yyzz[k] * ab_z + g_xx_0_xx_yyzzz[k];

                g_xx_0_xxz_yzzz[k] = -g_xx_0_xx_yzzz[k] * ab_z + g_xx_0_xx_yzzzz[k];

                g_xx_0_xxz_zzzz[k] = -g_xx_0_xx_zzzz[k] * ab_z + g_xx_0_xx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xy_xxxx, g_xx_0_xy_xxxxy, g_xx_0_xy_xxxy, g_xx_0_xy_xxxyy, g_xx_0_xy_xxxyz, g_xx_0_xy_xxxz, g_xx_0_xy_xxyy, g_xx_0_xy_xxyyy, g_xx_0_xy_xxyyz, g_xx_0_xy_xxyz, g_xx_0_xy_xxyzz, g_xx_0_xy_xxzz, g_xx_0_xy_xyyy, g_xx_0_xy_xyyyy, g_xx_0_xy_xyyyz, g_xx_0_xy_xyyz, g_xx_0_xy_xyyzz, g_xx_0_xy_xyzz, g_xx_0_xy_xyzzz, g_xx_0_xy_xzzz, g_xx_0_xy_yyyy, g_xx_0_xy_yyyyy, g_xx_0_xy_yyyyz, g_xx_0_xy_yyyz, g_xx_0_xy_yyyzz, g_xx_0_xy_yyzz, g_xx_0_xy_yyzzz, g_xx_0_xy_yzzz, g_xx_0_xy_yzzzz, g_xx_0_xy_zzzz, g_xx_0_xyy_xxxx, g_xx_0_xyy_xxxy, g_xx_0_xyy_xxxz, g_xx_0_xyy_xxyy, g_xx_0_xyy_xxyz, g_xx_0_xyy_xxzz, g_xx_0_xyy_xyyy, g_xx_0_xyy_xyyz, g_xx_0_xyy_xyzz, g_xx_0_xyy_xzzz, g_xx_0_xyy_yyyy, g_xx_0_xyy_yyyz, g_xx_0_xyy_yyzz, g_xx_0_xyy_yzzz, g_xx_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyy_xxxx[k] = -g_xx_0_xy_xxxx[k] * ab_y + g_xx_0_xy_xxxxy[k];

                g_xx_0_xyy_xxxy[k] = -g_xx_0_xy_xxxy[k] * ab_y + g_xx_0_xy_xxxyy[k];

                g_xx_0_xyy_xxxz[k] = -g_xx_0_xy_xxxz[k] * ab_y + g_xx_0_xy_xxxyz[k];

                g_xx_0_xyy_xxyy[k] = -g_xx_0_xy_xxyy[k] * ab_y + g_xx_0_xy_xxyyy[k];

                g_xx_0_xyy_xxyz[k] = -g_xx_0_xy_xxyz[k] * ab_y + g_xx_0_xy_xxyyz[k];

                g_xx_0_xyy_xxzz[k] = -g_xx_0_xy_xxzz[k] * ab_y + g_xx_0_xy_xxyzz[k];

                g_xx_0_xyy_xyyy[k] = -g_xx_0_xy_xyyy[k] * ab_y + g_xx_0_xy_xyyyy[k];

                g_xx_0_xyy_xyyz[k] = -g_xx_0_xy_xyyz[k] * ab_y + g_xx_0_xy_xyyyz[k];

                g_xx_0_xyy_xyzz[k] = -g_xx_0_xy_xyzz[k] * ab_y + g_xx_0_xy_xyyzz[k];

                g_xx_0_xyy_xzzz[k] = -g_xx_0_xy_xzzz[k] * ab_y + g_xx_0_xy_xyzzz[k];

                g_xx_0_xyy_yyyy[k] = -g_xx_0_xy_yyyy[k] * ab_y + g_xx_0_xy_yyyyy[k];

                g_xx_0_xyy_yyyz[k] = -g_xx_0_xy_yyyz[k] * ab_y + g_xx_0_xy_yyyyz[k];

                g_xx_0_xyy_yyzz[k] = -g_xx_0_xy_yyzz[k] * ab_y + g_xx_0_xy_yyyzz[k];

                g_xx_0_xyy_yzzz[k] = -g_xx_0_xy_yzzz[k] * ab_y + g_xx_0_xy_yyzzz[k];

                g_xx_0_xyy_zzzz[k] = -g_xx_0_xy_zzzz[k] * ab_y + g_xx_0_xy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyz_xxxx, g_xx_0_xyz_xxxy, g_xx_0_xyz_xxxz, g_xx_0_xyz_xxyy, g_xx_0_xyz_xxyz, g_xx_0_xyz_xxzz, g_xx_0_xyz_xyyy, g_xx_0_xyz_xyyz, g_xx_0_xyz_xyzz, g_xx_0_xyz_xzzz, g_xx_0_xyz_yyyy, g_xx_0_xyz_yyyz, g_xx_0_xyz_yyzz, g_xx_0_xyz_yzzz, g_xx_0_xyz_zzzz, g_xx_0_xz_xxxx, g_xx_0_xz_xxxxy, g_xx_0_xz_xxxy, g_xx_0_xz_xxxyy, g_xx_0_xz_xxxyz, g_xx_0_xz_xxxz, g_xx_0_xz_xxyy, g_xx_0_xz_xxyyy, g_xx_0_xz_xxyyz, g_xx_0_xz_xxyz, g_xx_0_xz_xxyzz, g_xx_0_xz_xxzz, g_xx_0_xz_xyyy, g_xx_0_xz_xyyyy, g_xx_0_xz_xyyyz, g_xx_0_xz_xyyz, g_xx_0_xz_xyyzz, g_xx_0_xz_xyzz, g_xx_0_xz_xyzzz, g_xx_0_xz_xzzz, g_xx_0_xz_yyyy, g_xx_0_xz_yyyyy, g_xx_0_xz_yyyyz, g_xx_0_xz_yyyz, g_xx_0_xz_yyyzz, g_xx_0_xz_yyzz, g_xx_0_xz_yyzzz, g_xx_0_xz_yzzz, g_xx_0_xz_yzzzz, g_xx_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyz_xxxx[k] = -g_xx_0_xz_xxxx[k] * ab_y + g_xx_0_xz_xxxxy[k];

                g_xx_0_xyz_xxxy[k] = -g_xx_0_xz_xxxy[k] * ab_y + g_xx_0_xz_xxxyy[k];

                g_xx_0_xyz_xxxz[k] = -g_xx_0_xz_xxxz[k] * ab_y + g_xx_0_xz_xxxyz[k];

                g_xx_0_xyz_xxyy[k] = -g_xx_0_xz_xxyy[k] * ab_y + g_xx_0_xz_xxyyy[k];

                g_xx_0_xyz_xxyz[k] = -g_xx_0_xz_xxyz[k] * ab_y + g_xx_0_xz_xxyyz[k];

                g_xx_0_xyz_xxzz[k] = -g_xx_0_xz_xxzz[k] * ab_y + g_xx_0_xz_xxyzz[k];

                g_xx_0_xyz_xyyy[k] = -g_xx_0_xz_xyyy[k] * ab_y + g_xx_0_xz_xyyyy[k];

                g_xx_0_xyz_xyyz[k] = -g_xx_0_xz_xyyz[k] * ab_y + g_xx_0_xz_xyyyz[k];

                g_xx_0_xyz_xyzz[k] = -g_xx_0_xz_xyzz[k] * ab_y + g_xx_0_xz_xyyzz[k];

                g_xx_0_xyz_xzzz[k] = -g_xx_0_xz_xzzz[k] * ab_y + g_xx_0_xz_xyzzz[k];

                g_xx_0_xyz_yyyy[k] = -g_xx_0_xz_yyyy[k] * ab_y + g_xx_0_xz_yyyyy[k];

                g_xx_0_xyz_yyyz[k] = -g_xx_0_xz_yyyz[k] * ab_y + g_xx_0_xz_yyyyz[k];

                g_xx_0_xyz_yyzz[k] = -g_xx_0_xz_yyzz[k] * ab_y + g_xx_0_xz_yyyzz[k];

                g_xx_0_xyz_yzzz[k] = -g_xx_0_xz_yzzz[k] * ab_y + g_xx_0_xz_yyzzz[k];

                g_xx_0_xyz_zzzz[k] = -g_xx_0_xz_zzzz[k] * ab_y + g_xx_0_xz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xz_xxxx, g_xx_0_xz_xxxxz, g_xx_0_xz_xxxy, g_xx_0_xz_xxxyz, g_xx_0_xz_xxxz, g_xx_0_xz_xxxzz, g_xx_0_xz_xxyy, g_xx_0_xz_xxyyz, g_xx_0_xz_xxyz, g_xx_0_xz_xxyzz, g_xx_0_xz_xxzz, g_xx_0_xz_xxzzz, g_xx_0_xz_xyyy, g_xx_0_xz_xyyyz, g_xx_0_xz_xyyz, g_xx_0_xz_xyyzz, g_xx_0_xz_xyzz, g_xx_0_xz_xyzzz, g_xx_0_xz_xzzz, g_xx_0_xz_xzzzz, g_xx_0_xz_yyyy, g_xx_0_xz_yyyyz, g_xx_0_xz_yyyz, g_xx_0_xz_yyyzz, g_xx_0_xz_yyzz, g_xx_0_xz_yyzzz, g_xx_0_xz_yzzz, g_xx_0_xz_yzzzz, g_xx_0_xz_zzzz, g_xx_0_xz_zzzzz, g_xx_0_xzz_xxxx, g_xx_0_xzz_xxxy, g_xx_0_xzz_xxxz, g_xx_0_xzz_xxyy, g_xx_0_xzz_xxyz, g_xx_0_xzz_xxzz, g_xx_0_xzz_xyyy, g_xx_0_xzz_xyyz, g_xx_0_xzz_xyzz, g_xx_0_xzz_xzzz, g_xx_0_xzz_yyyy, g_xx_0_xzz_yyyz, g_xx_0_xzz_yyzz, g_xx_0_xzz_yzzz, g_xx_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzz_xxxx[k] = -g_xx_0_xz_xxxx[k] * ab_z + g_xx_0_xz_xxxxz[k];

                g_xx_0_xzz_xxxy[k] = -g_xx_0_xz_xxxy[k] * ab_z + g_xx_0_xz_xxxyz[k];

                g_xx_0_xzz_xxxz[k] = -g_xx_0_xz_xxxz[k] * ab_z + g_xx_0_xz_xxxzz[k];

                g_xx_0_xzz_xxyy[k] = -g_xx_0_xz_xxyy[k] * ab_z + g_xx_0_xz_xxyyz[k];

                g_xx_0_xzz_xxyz[k] = -g_xx_0_xz_xxyz[k] * ab_z + g_xx_0_xz_xxyzz[k];

                g_xx_0_xzz_xxzz[k] = -g_xx_0_xz_xxzz[k] * ab_z + g_xx_0_xz_xxzzz[k];

                g_xx_0_xzz_xyyy[k] = -g_xx_0_xz_xyyy[k] * ab_z + g_xx_0_xz_xyyyz[k];

                g_xx_0_xzz_xyyz[k] = -g_xx_0_xz_xyyz[k] * ab_z + g_xx_0_xz_xyyzz[k];

                g_xx_0_xzz_xyzz[k] = -g_xx_0_xz_xyzz[k] * ab_z + g_xx_0_xz_xyzzz[k];

                g_xx_0_xzz_xzzz[k] = -g_xx_0_xz_xzzz[k] * ab_z + g_xx_0_xz_xzzzz[k];

                g_xx_0_xzz_yyyy[k] = -g_xx_0_xz_yyyy[k] * ab_z + g_xx_0_xz_yyyyz[k];

                g_xx_0_xzz_yyyz[k] = -g_xx_0_xz_yyyz[k] * ab_z + g_xx_0_xz_yyyzz[k];

                g_xx_0_xzz_yyzz[k] = -g_xx_0_xz_yyzz[k] * ab_z + g_xx_0_xz_yyzzz[k];

                g_xx_0_xzz_yzzz[k] = -g_xx_0_xz_yzzz[k] * ab_z + g_xx_0_xz_yzzzz[k];

                g_xx_0_xzz_zzzz[k] = -g_xx_0_xz_zzzz[k] * ab_z + g_xx_0_xz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yy_xxxx, g_xx_0_yy_xxxxy, g_xx_0_yy_xxxy, g_xx_0_yy_xxxyy, g_xx_0_yy_xxxyz, g_xx_0_yy_xxxz, g_xx_0_yy_xxyy, g_xx_0_yy_xxyyy, g_xx_0_yy_xxyyz, g_xx_0_yy_xxyz, g_xx_0_yy_xxyzz, g_xx_0_yy_xxzz, g_xx_0_yy_xyyy, g_xx_0_yy_xyyyy, g_xx_0_yy_xyyyz, g_xx_0_yy_xyyz, g_xx_0_yy_xyyzz, g_xx_0_yy_xyzz, g_xx_0_yy_xyzzz, g_xx_0_yy_xzzz, g_xx_0_yy_yyyy, g_xx_0_yy_yyyyy, g_xx_0_yy_yyyyz, g_xx_0_yy_yyyz, g_xx_0_yy_yyyzz, g_xx_0_yy_yyzz, g_xx_0_yy_yyzzz, g_xx_0_yy_yzzz, g_xx_0_yy_yzzzz, g_xx_0_yy_zzzz, g_xx_0_yyy_xxxx, g_xx_0_yyy_xxxy, g_xx_0_yyy_xxxz, g_xx_0_yyy_xxyy, g_xx_0_yyy_xxyz, g_xx_0_yyy_xxzz, g_xx_0_yyy_xyyy, g_xx_0_yyy_xyyz, g_xx_0_yyy_xyzz, g_xx_0_yyy_xzzz, g_xx_0_yyy_yyyy, g_xx_0_yyy_yyyz, g_xx_0_yyy_yyzz, g_xx_0_yyy_yzzz, g_xx_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyy_xxxx[k] = -g_xx_0_yy_xxxx[k] * ab_y + g_xx_0_yy_xxxxy[k];

                g_xx_0_yyy_xxxy[k] = -g_xx_0_yy_xxxy[k] * ab_y + g_xx_0_yy_xxxyy[k];

                g_xx_0_yyy_xxxz[k] = -g_xx_0_yy_xxxz[k] * ab_y + g_xx_0_yy_xxxyz[k];

                g_xx_0_yyy_xxyy[k] = -g_xx_0_yy_xxyy[k] * ab_y + g_xx_0_yy_xxyyy[k];

                g_xx_0_yyy_xxyz[k] = -g_xx_0_yy_xxyz[k] * ab_y + g_xx_0_yy_xxyyz[k];

                g_xx_0_yyy_xxzz[k] = -g_xx_0_yy_xxzz[k] * ab_y + g_xx_0_yy_xxyzz[k];

                g_xx_0_yyy_xyyy[k] = -g_xx_0_yy_xyyy[k] * ab_y + g_xx_0_yy_xyyyy[k];

                g_xx_0_yyy_xyyz[k] = -g_xx_0_yy_xyyz[k] * ab_y + g_xx_0_yy_xyyyz[k];

                g_xx_0_yyy_xyzz[k] = -g_xx_0_yy_xyzz[k] * ab_y + g_xx_0_yy_xyyzz[k];

                g_xx_0_yyy_xzzz[k] = -g_xx_0_yy_xzzz[k] * ab_y + g_xx_0_yy_xyzzz[k];

                g_xx_0_yyy_yyyy[k] = -g_xx_0_yy_yyyy[k] * ab_y + g_xx_0_yy_yyyyy[k];

                g_xx_0_yyy_yyyz[k] = -g_xx_0_yy_yyyz[k] * ab_y + g_xx_0_yy_yyyyz[k];

                g_xx_0_yyy_yyzz[k] = -g_xx_0_yy_yyzz[k] * ab_y + g_xx_0_yy_yyyzz[k];

                g_xx_0_yyy_yzzz[k] = -g_xx_0_yy_yzzz[k] * ab_y + g_xx_0_yy_yyzzz[k];

                g_xx_0_yyy_zzzz[k] = -g_xx_0_yy_zzzz[k] * ab_y + g_xx_0_yy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyz_xxxx, g_xx_0_yyz_xxxy, g_xx_0_yyz_xxxz, g_xx_0_yyz_xxyy, g_xx_0_yyz_xxyz, g_xx_0_yyz_xxzz, g_xx_0_yyz_xyyy, g_xx_0_yyz_xyyz, g_xx_0_yyz_xyzz, g_xx_0_yyz_xzzz, g_xx_0_yyz_yyyy, g_xx_0_yyz_yyyz, g_xx_0_yyz_yyzz, g_xx_0_yyz_yzzz, g_xx_0_yyz_zzzz, g_xx_0_yz_xxxx, g_xx_0_yz_xxxxy, g_xx_0_yz_xxxy, g_xx_0_yz_xxxyy, g_xx_0_yz_xxxyz, g_xx_0_yz_xxxz, g_xx_0_yz_xxyy, g_xx_0_yz_xxyyy, g_xx_0_yz_xxyyz, g_xx_0_yz_xxyz, g_xx_0_yz_xxyzz, g_xx_0_yz_xxzz, g_xx_0_yz_xyyy, g_xx_0_yz_xyyyy, g_xx_0_yz_xyyyz, g_xx_0_yz_xyyz, g_xx_0_yz_xyyzz, g_xx_0_yz_xyzz, g_xx_0_yz_xyzzz, g_xx_0_yz_xzzz, g_xx_0_yz_yyyy, g_xx_0_yz_yyyyy, g_xx_0_yz_yyyyz, g_xx_0_yz_yyyz, g_xx_0_yz_yyyzz, g_xx_0_yz_yyzz, g_xx_0_yz_yyzzz, g_xx_0_yz_yzzz, g_xx_0_yz_yzzzz, g_xx_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyz_xxxx[k] = -g_xx_0_yz_xxxx[k] * ab_y + g_xx_0_yz_xxxxy[k];

                g_xx_0_yyz_xxxy[k] = -g_xx_0_yz_xxxy[k] * ab_y + g_xx_0_yz_xxxyy[k];

                g_xx_0_yyz_xxxz[k] = -g_xx_0_yz_xxxz[k] * ab_y + g_xx_0_yz_xxxyz[k];

                g_xx_0_yyz_xxyy[k] = -g_xx_0_yz_xxyy[k] * ab_y + g_xx_0_yz_xxyyy[k];

                g_xx_0_yyz_xxyz[k] = -g_xx_0_yz_xxyz[k] * ab_y + g_xx_0_yz_xxyyz[k];

                g_xx_0_yyz_xxzz[k] = -g_xx_0_yz_xxzz[k] * ab_y + g_xx_0_yz_xxyzz[k];

                g_xx_0_yyz_xyyy[k] = -g_xx_0_yz_xyyy[k] * ab_y + g_xx_0_yz_xyyyy[k];

                g_xx_0_yyz_xyyz[k] = -g_xx_0_yz_xyyz[k] * ab_y + g_xx_0_yz_xyyyz[k];

                g_xx_0_yyz_xyzz[k] = -g_xx_0_yz_xyzz[k] * ab_y + g_xx_0_yz_xyyzz[k];

                g_xx_0_yyz_xzzz[k] = -g_xx_0_yz_xzzz[k] * ab_y + g_xx_0_yz_xyzzz[k];

                g_xx_0_yyz_yyyy[k] = -g_xx_0_yz_yyyy[k] * ab_y + g_xx_0_yz_yyyyy[k];

                g_xx_0_yyz_yyyz[k] = -g_xx_0_yz_yyyz[k] * ab_y + g_xx_0_yz_yyyyz[k];

                g_xx_0_yyz_yyzz[k] = -g_xx_0_yz_yyzz[k] * ab_y + g_xx_0_yz_yyyzz[k];

                g_xx_0_yyz_yzzz[k] = -g_xx_0_yz_yzzz[k] * ab_y + g_xx_0_yz_yyzzz[k];

                g_xx_0_yyz_zzzz[k] = -g_xx_0_yz_zzzz[k] * ab_y + g_xx_0_yz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzz_xxxx, g_xx_0_yzz_xxxy, g_xx_0_yzz_xxxz, g_xx_0_yzz_xxyy, g_xx_0_yzz_xxyz, g_xx_0_yzz_xxzz, g_xx_0_yzz_xyyy, g_xx_0_yzz_xyyz, g_xx_0_yzz_xyzz, g_xx_0_yzz_xzzz, g_xx_0_yzz_yyyy, g_xx_0_yzz_yyyz, g_xx_0_yzz_yyzz, g_xx_0_yzz_yzzz, g_xx_0_yzz_zzzz, g_xx_0_zz_xxxx, g_xx_0_zz_xxxxy, g_xx_0_zz_xxxy, g_xx_0_zz_xxxyy, g_xx_0_zz_xxxyz, g_xx_0_zz_xxxz, g_xx_0_zz_xxyy, g_xx_0_zz_xxyyy, g_xx_0_zz_xxyyz, g_xx_0_zz_xxyz, g_xx_0_zz_xxyzz, g_xx_0_zz_xxzz, g_xx_0_zz_xyyy, g_xx_0_zz_xyyyy, g_xx_0_zz_xyyyz, g_xx_0_zz_xyyz, g_xx_0_zz_xyyzz, g_xx_0_zz_xyzz, g_xx_0_zz_xyzzz, g_xx_0_zz_xzzz, g_xx_0_zz_yyyy, g_xx_0_zz_yyyyy, g_xx_0_zz_yyyyz, g_xx_0_zz_yyyz, g_xx_0_zz_yyyzz, g_xx_0_zz_yyzz, g_xx_0_zz_yyzzz, g_xx_0_zz_yzzz, g_xx_0_zz_yzzzz, g_xx_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzz_xxxx[k] = -g_xx_0_zz_xxxx[k] * ab_y + g_xx_0_zz_xxxxy[k];

                g_xx_0_yzz_xxxy[k] = -g_xx_0_zz_xxxy[k] * ab_y + g_xx_0_zz_xxxyy[k];

                g_xx_0_yzz_xxxz[k] = -g_xx_0_zz_xxxz[k] * ab_y + g_xx_0_zz_xxxyz[k];

                g_xx_0_yzz_xxyy[k] = -g_xx_0_zz_xxyy[k] * ab_y + g_xx_0_zz_xxyyy[k];

                g_xx_0_yzz_xxyz[k] = -g_xx_0_zz_xxyz[k] * ab_y + g_xx_0_zz_xxyyz[k];

                g_xx_0_yzz_xxzz[k] = -g_xx_0_zz_xxzz[k] * ab_y + g_xx_0_zz_xxyzz[k];

                g_xx_0_yzz_xyyy[k] = -g_xx_0_zz_xyyy[k] * ab_y + g_xx_0_zz_xyyyy[k];

                g_xx_0_yzz_xyyz[k] = -g_xx_0_zz_xyyz[k] * ab_y + g_xx_0_zz_xyyyz[k];

                g_xx_0_yzz_xyzz[k] = -g_xx_0_zz_xyzz[k] * ab_y + g_xx_0_zz_xyyzz[k];

                g_xx_0_yzz_xzzz[k] = -g_xx_0_zz_xzzz[k] * ab_y + g_xx_0_zz_xyzzz[k];

                g_xx_0_yzz_yyyy[k] = -g_xx_0_zz_yyyy[k] * ab_y + g_xx_0_zz_yyyyy[k];

                g_xx_0_yzz_yyyz[k] = -g_xx_0_zz_yyyz[k] * ab_y + g_xx_0_zz_yyyyz[k];

                g_xx_0_yzz_yyzz[k] = -g_xx_0_zz_yyzz[k] * ab_y + g_xx_0_zz_yyyzz[k];

                g_xx_0_yzz_yzzz[k] = -g_xx_0_zz_yzzz[k] * ab_y + g_xx_0_zz_yyzzz[k];

                g_xx_0_yzz_zzzz[k] = -g_xx_0_zz_zzzz[k] * ab_y + g_xx_0_zz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zz_xxxx, g_xx_0_zz_xxxxz, g_xx_0_zz_xxxy, g_xx_0_zz_xxxyz, g_xx_0_zz_xxxz, g_xx_0_zz_xxxzz, g_xx_0_zz_xxyy, g_xx_0_zz_xxyyz, g_xx_0_zz_xxyz, g_xx_0_zz_xxyzz, g_xx_0_zz_xxzz, g_xx_0_zz_xxzzz, g_xx_0_zz_xyyy, g_xx_0_zz_xyyyz, g_xx_0_zz_xyyz, g_xx_0_zz_xyyzz, g_xx_0_zz_xyzz, g_xx_0_zz_xyzzz, g_xx_0_zz_xzzz, g_xx_0_zz_xzzzz, g_xx_0_zz_yyyy, g_xx_0_zz_yyyyz, g_xx_0_zz_yyyz, g_xx_0_zz_yyyzz, g_xx_0_zz_yyzz, g_xx_0_zz_yyzzz, g_xx_0_zz_yzzz, g_xx_0_zz_yzzzz, g_xx_0_zz_zzzz, g_xx_0_zz_zzzzz, g_xx_0_zzz_xxxx, g_xx_0_zzz_xxxy, g_xx_0_zzz_xxxz, g_xx_0_zzz_xxyy, g_xx_0_zzz_xxyz, g_xx_0_zzz_xxzz, g_xx_0_zzz_xyyy, g_xx_0_zzz_xyyz, g_xx_0_zzz_xyzz, g_xx_0_zzz_xzzz, g_xx_0_zzz_yyyy, g_xx_0_zzz_yyyz, g_xx_0_zzz_yyzz, g_xx_0_zzz_yzzz, g_xx_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzz_xxxx[k] = -g_xx_0_zz_xxxx[k] * ab_z + g_xx_0_zz_xxxxz[k];

                g_xx_0_zzz_xxxy[k] = -g_xx_0_zz_xxxy[k] * ab_z + g_xx_0_zz_xxxyz[k];

                g_xx_0_zzz_xxxz[k] = -g_xx_0_zz_xxxz[k] * ab_z + g_xx_0_zz_xxxzz[k];

                g_xx_0_zzz_xxyy[k] = -g_xx_0_zz_xxyy[k] * ab_z + g_xx_0_zz_xxyyz[k];

                g_xx_0_zzz_xxyz[k] = -g_xx_0_zz_xxyz[k] * ab_z + g_xx_0_zz_xxyzz[k];

                g_xx_0_zzz_xxzz[k] = -g_xx_0_zz_xxzz[k] * ab_z + g_xx_0_zz_xxzzz[k];

                g_xx_0_zzz_xyyy[k] = -g_xx_0_zz_xyyy[k] * ab_z + g_xx_0_zz_xyyyz[k];

                g_xx_0_zzz_xyyz[k] = -g_xx_0_zz_xyyz[k] * ab_z + g_xx_0_zz_xyyzz[k];

                g_xx_0_zzz_xyzz[k] = -g_xx_0_zz_xyzz[k] * ab_z + g_xx_0_zz_xyzzz[k];

                g_xx_0_zzz_xzzz[k] = -g_xx_0_zz_xzzz[k] * ab_z + g_xx_0_zz_xzzzz[k];

                g_xx_0_zzz_yyyy[k] = -g_xx_0_zz_yyyy[k] * ab_z + g_xx_0_zz_yyyyz[k];

                g_xx_0_zzz_yyyz[k] = -g_xx_0_zz_yyyz[k] * ab_z + g_xx_0_zz_yyyzz[k];

                g_xx_0_zzz_yyzz[k] = -g_xx_0_zz_yyzz[k] * ab_z + g_xx_0_zz_yyzzz[k];

                g_xx_0_zzz_yzzz[k] = -g_xx_0_zz_yzzz[k] * ab_z + g_xx_0_zz_yzzzz[k];

                g_xx_0_zzz_zzzz[k] = -g_xx_0_zz_zzzz[k] * ab_z + g_xx_0_zz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_xxxx, g_xy_0_xx_xxxxx, g_xy_0_xx_xxxxy, g_xy_0_xx_xxxxz, g_xy_0_xx_xxxy, g_xy_0_xx_xxxyy, g_xy_0_xx_xxxyz, g_xy_0_xx_xxxz, g_xy_0_xx_xxxzz, g_xy_0_xx_xxyy, g_xy_0_xx_xxyyy, g_xy_0_xx_xxyyz, g_xy_0_xx_xxyz, g_xy_0_xx_xxyzz, g_xy_0_xx_xxzz, g_xy_0_xx_xxzzz, g_xy_0_xx_xyyy, g_xy_0_xx_xyyyy, g_xy_0_xx_xyyyz, g_xy_0_xx_xyyz, g_xy_0_xx_xyyzz, g_xy_0_xx_xyzz, g_xy_0_xx_xyzzz, g_xy_0_xx_xzzz, g_xy_0_xx_xzzzz, g_xy_0_xx_yyyy, g_xy_0_xx_yyyz, g_xy_0_xx_yyzz, g_xy_0_xx_yzzz, g_xy_0_xx_zzzz, g_xy_0_xxx_xxxx, g_xy_0_xxx_xxxy, g_xy_0_xxx_xxxz, g_xy_0_xxx_xxyy, g_xy_0_xxx_xxyz, g_xy_0_xxx_xxzz, g_xy_0_xxx_xyyy, g_xy_0_xxx_xyyz, g_xy_0_xxx_xyzz, g_xy_0_xxx_xzzz, g_xy_0_xxx_yyyy, g_xy_0_xxx_yyyz, g_xy_0_xxx_yyzz, g_xy_0_xxx_yzzz, g_xy_0_xxx_zzzz, g_y_0_xx_xxxx, g_y_0_xx_xxxy, g_y_0_xx_xxxz, g_y_0_xx_xxyy, g_y_0_xx_xxyz, g_y_0_xx_xxzz, g_y_0_xx_xyyy, g_y_0_xx_xyyz, g_y_0_xx_xyzz, g_y_0_xx_xzzz, g_y_0_xx_yyyy, g_y_0_xx_yyyz, g_y_0_xx_yyzz, g_y_0_xx_yzzz, g_y_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxx_xxxx[k] = -g_y_0_xx_xxxx[k] - g_xy_0_xx_xxxx[k] * ab_x + g_xy_0_xx_xxxxx[k];

                g_xy_0_xxx_xxxy[k] = -g_y_0_xx_xxxy[k] - g_xy_0_xx_xxxy[k] * ab_x + g_xy_0_xx_xxxxy[k];

                g_xy_0_xxx_xxxz[k] = -g_y_0_xx_xxxz[k] - g_xy_0_xx_xxxz[k] * ab_x + g_xy_0_xx_xxxxz[k];

                g_xy_0_xxx_xxyy[k] = -g_y_0_xx_xxyy[k] - g_xy_0_xx_xxyy[k] * ab_x + g_xy_0_xx_xxxyy[k];

                g_xy_0_xxx_xxyz[k] = -g_y_0_xx_xxyz[k] - g_xy_0_xx_xxyz[k] * ab_x + g_xy_0_xx_xxxyz[k];

                g_xy_0_xxx_xxzz[k] = -g_y_0_xx_xxzz[k] - g_xy_0_xx_xxzz[k] * ab_x + g_xy_0_xx_xxxzz[k];

                g_xy_0_xxx_xyyy[k] = -g_y_0_xx_xyyy[k] - g_xy_0_xx_xyyy[k] * ab_x + g_xy_0_xx_xxyyy[k];

                g_xy_0_xxx_xyyz[k] = -g_y_0_xx_xyyz[k] - g_xy_0_xx_xyyz[k] * ab_x + g_xy_0_xx_xxyyz[k];

                g_xy_0_xxx_xyzz[k] = -g_y_0_xx_xyzz[k] - g_xy_0_xx_xyzz[k] * ab_x + g_xy_0_xx_xxyzz[k];

                g_xy_0_xxx_xzzz[k] = -g_y_0_xx_xzzz[k] - g_xy_0_xx_xzzz[k] * ab_x + g_xy_0_xx_xxzzz[k];

                g_xy_0_xxx_yyyy[k] = -g_y_0_xx_yyyy[k] - g_xy_0_xx_yyyy[k] * ab_x + g_xy_0_xx_xyyyy[k];

                g_xy_0_xxx_yyyz[k] = -g_y_0_xx_yyyz[k] - g_xy_0_xx_yyyz[k] * ab_x + g_xy_0_xx_xyyyz[k];

                g_xy_0_xxx_yyzz[k] = -g_y_0_xx_yyzz[k] - g_xy_0_xx_yyzz[k] * ab_x + g_xy_0_xx_xyyzz[k];

                g_xy_0_xxx_yzzz[k] = -g_y_0_xx_yzzz[k] - g_xy_0_xx_yzzz[k] * ab_x + g_xy_0_xx_xyzzz[k];

                g_xy_0_xxx_zzzz[k] = -g_y_0_xx_zzzz[k] - g_xy_0_xx_zzzz[k] * ab_x + g_xy_0_xx_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_xxxx, g_xy_0_xxy_xxxy, g_xy_0_xxy_xxxz, g_xy_0_xxy_xxyy, g_xy_0_xxy_xxyz, g_xy_0_xxy_xxzz, g_xy_0_xxy_xyyy, g_xy_0_xxy_xyyz, g_xy_0_xxy_xyzz, g_xy_0_xxy_xzzz, g_xy_0_xxy_yyyy, g_xy_0_xxy_yyyz, g_xy_0_xxy_yyzz, g_xy_0_xxy_yzzz, g_xy_0_xxy_zzzz, g_xy_0_xy_xxxx, g_xy_0_xy_xxxxx, g_xy_0_xy_xxxxy, g_xy_0_xy_xxxxz, g_xy_0_xy_xxxy, g_xy_0_xy_xxxyy, g_xy_0_xy_xxxyz, g_xy_0_xy_xxxz, g_xy_0_xy_xxxzz, g_xy_0_xy_xxyy, g_xy_0_xy_xxyyy, g_xy_0_xy_xxyyz, g_xy_0_xy_xxyz, g_xy_0_xy_xxyzz, g_xy_0_xy_xxzz, g_xy_0_xy_xxzzz, g_xy_0_xy_xyyy, g_xy_0_xy_xyyyy, g_xy_0_xy_xyyyz, g_xy_0_xy_xyyz, g_xy_0_xy_xyyzz, g_xy_0_xy_xyzz, g_xy_0_xy_xyzzz, g_xy_0_xy_xzzz, g_xy_0_xy_xzzzz, g_xy_0_xy_yyyy, g_xy_0_xy_yyyz, g_xy_0_xy_yyzz, g_xy_0_xy_yzzz, g_xy_0_xy_zzzz, g_y_0_xy_xxxx, g_y_0_xy_xxxy, g_y_0_xy_xxxz, g_y_0_xy_xxyy, g_y_0_xy_xxyz, g_y_0_xy_xxzz, g_y_0_xy_xyyy, g_y_0_xy_xyyz, g_y_0_xy_xyzz, g_y_0_xy_xzzz, g_y_0_xy_yyyy, g_y_0_xy_yyyz, g_y_0_xy_yyzz, g_y_0_xy_yzzz, g_y_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxy_xxxx[k] = -g_y_0_xy_xxxx[k] - g_xy_0_xy_xxxx[k] * ab_x + g_xy_0_xy_xxxxx[k];

                g_xy_0_xxy_xxxy[k] = -g_y_0_xy_xxxy[k] - g_xy_0_xy_xxxy[k] * ab_x + g_xy_0_xy_xxxxy[k];

                g_xy_0_xxy_xxxz[k] = -g_y_0_xy_xxxz[k] - g_xy_0_xy_xxxz[k] * ab_x + g_xy_0_xy_xxxxz[k];

                g_xy_0_xxy_xxyy[k] = -g_y_0_xy_xxyy[k] - g_xy_0_xy_xxyy[k] * ab_x + g_xy_0_xy_xxxyy[k];

                g_xy_0_xxy_xxyz[k] = -g_y_0_xy_xxyz[k] - g_xy_0_xy_xxyz[k] * ab_x + g_xy_0_xy_xxxyz[k];

                g_xy_0_xxy_xxzz[k] = -g_y_0_xy_xxzz[k] - g_xy_0_xy_xxzz[k] * ab_x + g_xy_0_xy_xxxzz[k];

                g_xy_0_xxy_xyyy[k] = -g_y_0_xy_xyyy[k] - g_xy_0_xy_xyyy[k] * ab_x + g_xy_0_xy_xxyyy[k];

                g_xy_0_xxy_xyyz[k] = -g_y_0_xy_xyyz[k] - g_xy_0_xy_xyyz[k] * ab_x + g_xy_0_xy_xxyyz[k];

                g_xy_0_xxy_xyzz[k] = -g_y_0_xy_xyzz[k] - g_xy_0_xy_xyzz[k] * ab_x + g_xy_0_xy_xxyzz[k];

                g_xy_0_xxy_xzzz[k] = -g_y_0_xy_xzzz[k] - g_xy_0_xy_xzzz[k] * ab_x + g_xy_0_xy_xxzzz[k];

                g_xy_0_xxy_yyyy[k] = -g_y_0_xy_yyyy[k] - g_xy_0_xy_yyyy[k] * ab_x + g_xy_0_xy_xyyyy[k];

                g_xy_0_xxy_yyyz[k] = -g_y_0_xy_yyyz[k] - g_xy_0_xy_yyyz[k] * ab_x + g_xy_0_xy_xyyyz[k];

                g_xy_0_xxy_yyzz[k] = -g_y_0_xy_yyzz[k] - g_xy_0_xy_yyzz[k] * ab_x + g_xy_0_xy_xyyzz[k];

                g_xy_0_xxy_yzzz[k] = -g_y_0_xy_yzzz[k] - g_xy_0_xy_yzzz[k] * ab_x + g_xy_0_xy_xyzzz[k];

                g_xy_0_xxy_zzzz[k] = -g_y_0_xy_zzzz[k] - g_xy_0_xy_zzzz[k] * ab_x + g_xy_0_xy_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_xxxx, g_xy_0_xx_xxxxz, g_xy_0_xx_xxxy, g_xy_0_xx_xxxyz, g_xy_0_xx_xxxz, g_xy_0_xx_xxxzz, g_xy_0_xx_xxyy, g_xy_0_xx_xxyyz, g_xy_0_xx_xxyz, g_xy_0_xx_xxyzz, g_xy_0_xx_xxzz, g_xy_0_xx_xxzzz, g_xy_0_xx_xyyy, g_xy_0_xx_xyyyz, g_xy_0_xx_xyyz, g_xy_0_xx_xyyzz, g_xy_0_xx_xyzz, g_xy_0_xx_xyzzz, g_xy_0_xx_xzzz, g_xy_0_xx_xzzzz, g_xy_0_xx_yyyy, g_xy_0_xx_yyyyz, g_xy_0_xx_yyyz, g_xy_0_xx_yyyzz, g_xy_0_xx_yyzz, g_xy_0_xx_yyzzz, g_xy_0_xx_yzzz, g_xy_0_xx_yzzzz, g_xy_0_xx_zzzz, g_xy_0_xx_zzzzz, g_xy_0_xxz_xxxx, g_xy_0_xxz_xxxy, g_xy_0_xxz_xxxz, g_xy_0_xxz_xxyy, g_xy_0_xxz_xxyz, g_xy_0_xxz_xxzz, g_xy_0_xxz_xyyy, g_xy_0_xxz_xyyz, g_xy_0_xxz_xyzz, g_xy_0_xxz_xzzz, g_xy_0_xxz_yyyy, g_xy_0_xxz_yyyz, g_xy_0_xxz_yyzz, g_xy_0_xxz_yzzz, g_xy_0_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxz_xxxx[k] = -g_xy_0_xx_xxxx[k] * ab_z + g_xy_0_xx_xxxxz[k];

                g_xy_0_xxz_xxxy[k] = -g_xy_0_xx_xxxy[k] * ab_z + g_xy_0_xx_xxxyz[k];

                g_xy_0_xxz_xxxz[k] = -g_xy_0_xx_xxxz[k] * ab_z + g_xy_0_xx_xxxzz[k];

                g_xy_0_xxz_xxyy[k] = -g_xy_0_xx_xxyy[k] * ab_z + g_xy_0_xx_xxyyz[k];

                g_xy_0_xxz_xxyz[k] = -g_xy_0_xx_xxyz[k] * ab_z + g_xy_0_xx_xxyzz[k];

                g_xy_0_xxz_xxzz[k] = -g_xy_0_xx_xxzz[k] * ab_z + g_xy_0_xx_xxzzz[k];

                g_xy_0_xxz_xyyy[k] = -g_xy_0_xx_xyyy[k] * ab_z + g_xy_0_xx_xyyyz[k];

                g_xy_0_xxz_xyyz[k] = -g_xy_0_xx_xyyz[k] * ab_z + g_xy_0_xx_xyyzz[k];

                g_xy_0_xxz_xyzz[k] = -g_xy_0_xx_xyzz[k] * ab_z + g_xy_0_xx_xyzzz[k];

                g_xy_0_xxz_xzzz[k] = -g_xy_0_xx_xzzz[k] * ab_z + g_xy_0_xx_xzzzz[k];

                g_xy_0_xxz_yyyy[k] = -g_xy_0_xx_yyyy[k] * ab_z + g_xy_0_xx_yyyyz[k];

                g_xy_0_xxz_yyyz[k] = -g_xy_0_xx_yyyz[k] * ab_z + g_xy_0_xx_yyyzz[k];

                g_xy_0_xxz_yyzz[k] = -g_xy_0_xx_yyzz[k] * ab_z + g_xy_0_xx_yyzzz[k];

                g_xy_0_xxz_yzzz[k] = -g_xy_0_xx_yzzz[k] * ab_z + g_xy_0_xx_yzzzz[k];

                g_xy_0_xxz_zzzz[k] = -g_xy_0_xx_zzzz[k] * ab_z + g_xy_0_xx_zzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_xxxx, g_xy_0_xyy_xxxy, g_xy_0_xyy_xxxz, g_xy_0_xyy_xxyy, g_xy_0_xyy_xxyz, g_xy_0_xyy_xxzz, g_xy_0_xyy_xyyy, g_xy_0_xyy_xyyz, g_xy_0_xyy_xyzz, g_xy_0_xyy_xzzz, g_xy_0_xyy_yyyy, g_xy_0_xyy_yyyz, g_xy_0_xyy_yyzz, g_xy_0_xyy_yzzz, g_xy_0_xyy_zzzz, g_xy_0_yy_xxxx, g_xy_0_yy_xxxxx, g_xy_0_yy_xxxxy, g_xy_0_yy_xxxxz, g_xy_0_yy_xxxy, g_xy_0_yy_xxxyy, g_xy_0_yy_xxxyz, g_xy_0_yy_xxxz, g_xy_0_yy_xxxzz, g_xy_0_yy_xxyy, g_xy_0_yy_xxyyy, g_xy_0_yy_xxyyz, g_xy_0_yy_xxyz, g_xy_0_yy_xxyzz, g_xy_0_yy_xxzz, g_xy_0_yy_xxzzz, g_xy_0_yy_xyyy, g_xy_0_yy_xyyyy, g_xy_0_yy_xyyyz, g_xy_0_yy_xyyz, g_xy_0_yy_xyyzz, g_xy_0_yy_xyzz, g_xy_0_yy_xyzzz, g_xy_0_yy_xzzz, g_xy_0_yy_xzzzz, g_xy_0_yy_yyyy, g_xy_0_yy_yyyz, g_xy_0_yy_yyzz, g_xy_0_yy_yzzz, g_xy_0_yy_zzzz, g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyy_xxxx[k] = -g_y_0_yy_xxxx[k] - g_xy_0_yy_xxxx[k] * ab_x + g_xy_0_yy_xxxxx[k];

                g_xy_0_xyy_xxxy[k] = -g_y_0_yy_xxxy[k] - g_xy_0_yy_xxxy[k] * ab_x + g_xy_0_yy_xxxxy[k];

                g_xy_0_xyy_xxxz[k] = -g_y_0_yy_xxxz[k] - g_xy_0_yy_xxxz[k] * ab_x + g_xy_0_yy_xxxxz[k];

                g_xy_0_xyy_xxyy[k] = -g_y_0_yy_xxyy[k] - g_xy_0_yy_xxyy[k] * ab_x + g_xy_0_yy_xxxyy[k];

                g_xy_0_xyy_xxyz[k] = -g_y_0_yy_xxyz[k] - g_xy_0_yy_xxyz[k] * ab_x + g_xy_0_yy_xxxyz[k];

                g_xy_0_xyy_xxzz[k] = -g_y_0_yy_xxzz[k] - g_xy_0_yy_xxzz[k] * ab_x + g_xy_0_yy_xxxzz[k];

                g_xy_0_xyy_xyyy[k] = -g_y_0_yy_xyyy[k] - g_xy_0_yy_xyyy[k] * ab_x + g_xy_0_yy_xxyyy[k];

                g_xy_0_xyy_xyyz[k] = -g_y_0_yy_xyyz[k] - g_xy_0_yy_xyyz[k] * ab_x + g_xy_0_yy_xxyyz[k];

                g_xy_0_xyy_xyzz[k] = -g_y_0_yy_xyzz[k] - g_xy_0_yy_xyzz[k] * ab_x + g_xy_0_yy_xxyzz[k];

                g_xy_0_xyy_xzzz[k] = -g_y_0_yy_xzzz[k] - g_xy_0_yy_xzzz[k] * ab_x + g_xy_0_yy_xxzzz[k];

                g_xy_0_xyy_yyyy[k] = -g_y_0_yy_yyyy[k] - g_xy_0_yy_yyyy[k] * ab_x + g_xy_0_yy_xyyyy[k];

                g_xy_0_xyy_yyyz[k] = -g_y_0_yy_yyyz[k] - g_xy_0_yy_yyyz[k] * ab_x + g_xy_0_yy_xyyyz[k];

                g_xy_0_xyy_yyzz[k] = -g_y_0_yy_yyzz[k] - g_xy_0_yy_yyzz[k] * ab_x + g_xy_0_yy_xyyzz[k];

                g_xy_0_xyy_yzzz[k] = -g_y_0_yy_yzzz[k] - g_xy_0_yy_yzzz[k] * ab_x + g_xy_0_yy_xyzzz[k];

                g_xy_0_xyy_zzzz[k] = -g_y_0_yy_zzzz[k] - g_xy_0_yy_zzzz[k] * ab_x + g_xy_0_yy_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_xxxx, g_xy_0_xy_xxxxz, g_xy_0_xy_xxxy, g_xy_0_xy_xxxyz, g_xy_0_xy_xxxz, g_xy_0_xy_xxxzz, g_xy_0_xy_xxyy, g_xy_0_xy_xxyyz, g_xy_0_xy_xxyz, g_xy_0_xy_xxyzz, g_xy_0_xy_xxzz, g_xy_0_xy_xxzzz, g_xy_0_xy_xyyy, g_xy_0_xy_xyyyz, g_xy_0_xy_xyyz, g_xy_0_xy_xyyzz, g_xy_0_xy_xyzz, g_xy_0_xy_xyzzz, g_xy_0_xy_xzzz, g_xy_0_xy_xzzzz, g_xy_0_xy_yyyy, g_xy_0_xy_yyyyz, g_xy_0_xy_yyyz, g_xy_0_xy_yyyzz, g_xy_0_xy_yyzz, g_xy_0_xy_yyzzz, g_xy_0_xy_yzzz, g_xy_0_xy_yzzzz, g_xy_0_xy_zzzz, g_xy_0_xy_zzzzz, g_xy_0_xyz_xxxx, g_xy_0_xyz_xxxy, g_xy_0_xyz_xxxz, g_xy_0_xyz_xxyy, g_xy_0_xyz_xxyz, g_xy_0_xyz_xxzz, g_xy_0_xyz_xyyy, g_xy_0_xyz_xyyz, g_xy_0_xyz_xyzz, g_xy_0_xyz_xzzz, g_xy_0_xyz_yyyy, g_xy_0_xyz_yyyz, g_xy_0_xyz_yyzz, g_xy_0_xyz_yzzz, g_xy_0_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyz_xxxx[k] = -g_xy_0_xy_xxxx[k] * ab_z + g_xy_0_xy_xxxxz[k];

                g_xy_0_xyz_xxxy[k] = -g_xy_0_xy_xxxy[k] * ab_z + g_xy_0_xy_xxxyz[k];

                g_xy_0_xyz_xxxz[k] = -g_xy_0_xy_xxxz[k] * ab_z + g_xy_0_xy_xxxzz[k];

                g_xy_0_xyz_xxyy[k] = -g_xy_0_xy_xxyy[k] * ab_z + g_xy_0_xy_xxyyz[k];

                g_xy_0_xyz_xxyz[k] = -g_xy_0_xy_xxyz[k] * ab_z + g_xy_0_xy_xxyzz[k];

                g_xy_0_xyz_xxzz[k] = -g_xy_0_xy_xxzz[k] * ab_z + g_xy_0_xy_xxzzz[k];

                g_xy_0_xyz_xyyy[k] = -g_xy_0_xy_xyyy[k] * ab_z + g_xy_0_xy_xyyyz[k];

                g_xy_0_xyz_xyyz[k] = -g_xy_0_xy_xyyz[k] * ab_z + g_xy_0_xy_xyyzz[k];

                g_xy_0_xyz_xyzz[k] = -g_xy_0_xy_xyzz[k] * ab_z + g_xy_0_xy_xyzzz[k];

                g_xy_0_xyz_xzzz[k] = -g_xy_0_xy_xzzz[k] * ab_z + g_xy_0_xy_xzzzz[k];

                g_xy_0_xyz_yyyy[k] = -g_xy_0_xy_yyyy[k] * ab_z + g_xy_0_xy_yyyyz[k];

                g_xy_0_xyz_yyyz[k] = -g_xy_0_xy_yyyz[k] * ab_z + g_xy_0_xy_yyyzz[k];

                g_xy_0_xyz_yyzz[k] = -g_xy_0_xy_yyzz[k] * ab_z + g_xy_0_xy_yyzzz[k];

                g_xy_0_xyz_yzzz[k] = -g_xy_0_xy_yzzz[k] * ab_z + g_xy_0_xy_yzzzz[k];

                g_xy_0_xyz_zzzz[k] = -g_xy_0_xy_zzzz[k] * ab_z + g_xy_0_xy_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xz_xxxx, g_xy_0_xz_xxxxz, g_xy_0_xz_xxxy, g_xy_0_xz_xxxyz, g_xy_0_xz_xxxz, g_xy_0_xz_xxxzz, g_xy_0_xz_xxyy, g_xy_0_xz_xxyyz, g_xy_0_xz_xxyz, g_xy_0_xz_xxyzz, g_xy_0_xz_xxzz, g_xy_0_xz_xxzzz, g_xy_0_xz_xyyy, g_xy_0_xz_xyyyz, g_xy_0_xz_xyyz, g_xy_0_xz_xyyzz, g_xy_0_xz_xyzz, g_xy_0_xz_xyzzz, g_xy_0_xz_xzzz, g_xy_0_xz_xzzzz, g_xy_0_xz_yyyy, g_xy_0_xz_yyyyz, g_xy_0_xz_yyyz, g_xy_0_xz_yyyzz, g_xy_0_xz_yyzz, g_xy_0_xz_yyzzz, g_xy_0_xz_yzzz, g_xy_0_xz_yzzzz, g_xy_0_xz_zzzz, g_xy_0_xz_zzzzz, g_xy_0_xzz_xxxx, g_xy_0_xzz_xxxy, g_xy_0_xzz_xxxz, g_xy_0_xzz_xxyy, g_xy_0_xzz_xxyz, g_xy_0_xzz_xxzz, g_xy_0_xzz_xyyy, g_xy_0_xzz_xyyz, g_xy_0_xzz_xyzz, g_xy_0_xzz_xzzz, g_xy_0_xzz_yyyy, g_xy_0_xzz_yyyz, g_xy_0_xzz_yyzz, g_xy_0_xzz_yzzz, g_xy_0_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzz_xxxx[k] = -g_xy_0_xz_xxxx[k] * ab_z + g_xy_0_xz_xxxxz[k];

                g_xy_0_xzz_xxxy[k] = -g_xy_0_xz_xxxy[k] * ab_z + g_xy_0_xz_xxxyz[k];

                g_xy_0_xzz_xxxz[k] = -g_xy_0_xz_xxxz[k] * ab_z + g_xy_0_xz_xxxzz[k];

                g_xy_0_xzz_xxyy[k] = -g_xy_0_xz_xxyy[k] * ab_z + g_xy_0_xz_xxyyz[k];

                g_xy_0_xzz_xxyz[k] = -g_xy_0_xz_xxyz[k] * ab_z + g_xy_0_xz_xxyzz[k];

                g_xy_0_xzz_xxzz[k] = -g_xy_0_xz_xxzz[k] * ab_z + g_xy_0_xz_xxzzz[k];

                g_xy_0_xzz_xyyy[k] = -g_xy_0_xz_xyyy[k] * ab_z + g_xy_0_xz_xyyyz[k];

                g_xy_0_xzz_xyyz[k] = -g_xy_0_xz_xyyz[k] * ab_z + g_xy_0_xz_xyyzz[k];

                g_xy_0_xzz_xyzz[k] = -g_xy_0_xz_xyzz[k] * ab_z + g_xy_0_xz_xyzzz[k];

                g_xy_0_xzz_xzzz[k] = -g_xy_0_xz_xzzz[k] * ab_z + g_xy_0_xz_xzzzz[k];

                g_xy_0_xzz_yyyy[k] = -g_xy_0_xz_yyyy[k] * ab_z + g_xy_0_xz_yyyyz[k];

                g_xy_0_xzz_yyyz[k] = -g_xy_0_xz_yyyz[k] * ab_z + g_xy_0_xz_yyyzz[k];

                g_xy_0_xzz_yyzz[k] = -g_xy_0_xz_yyzz[k] * ab_z + g_xy_0_xz_yyzzz[k];

                g_xy_0_xzz_yzzz[k] = -g_xy_0_xz_yzzz[k] * ab_z + g_xy_0_xz_yzzzz[k];

                g_xy_0_xzz_zzzz[k] = -g_xy_0_xz_zzzz[k] * ab_z + g_xy_0_xz_zzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xxxx, g_x_0_yy_xxxy, g_x_0_yy_xxxz, g_x_0_yy_xxyy, g_x_0_yy_xxyz, g_x_0_yy_xxzz, g_x_0_yy_xyyy, g_x_0_yy_xyyz, g_x_0_yy_xyzz, g_x_0_yy_xzzz, g_x_0_yy_yyyy, g_x_0_yy_yyyz, g_x_0_yy_yyzz, g_x_0_yy_yzzz, g_x_0_yy_zzzz, g_xy_0_yy_xxxx, g_xy_0_yy_xxxxy, g_xy_0_yy_xxxy, g_xy_0_yy_xxxyy, g_xy_0_yy_xxxyz, g_xy_0_yy_xxxz, g_xy_0_yy_xxyy, g_xy_0_yy_xxyyy, g_xy_0_yy_xxyyz, g_xy_0_yy_xxyz, g_xy_0_yy_xxyzz, g_xy_0_yy_xxzz, g_xy_0_yy_xyyy, g_xy_0_yy_xyyyy, g_xy_0_yy_xyyyz, g_xy_0_yy_xyyz, g_xy_0_yy_xyyzz, g_xy_0_yy_xyzz, g_xy_0_yy_xyzzz, g_xy_0_yy_xzzz, g_xy_0_yy_yyyy, g_xy_0_yy_yyyyy, g_xy_0_yy_yyyyz, g_xy_0_yy_yyyz, g_xy_0_yy_yyyzz, g_xy_0_yy_yyzz, g_xy_0_yy_yyzzz, g_xy_0_yy_yzzz, g_xy_0_yy_yzzzz, g_xy_0_yy_zzzz, g_xy_0_yyy_xxxx, g_xy_0_yyy_xxxy, g_xy_0_yyy_xxxz, g_xy_0_yyy_xxyy, g_xy_0_yyy_xxyz, g_xy_0_yyy_xxzz, g_xy_0_yyy_xyyy, g_xy_0_yyy_xyyz, g_xy_0_yyy_xyzz, g_xy_0_yyy_xzzz, g_xy_0_yyy_yyyy, g_xy_0_yyy_yyyz, g_xy_0_yyy_yyzz, g_xy_0_yyy_yzzz, g_xy_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyy_xxxx[k] = -g_x_0_yy_xxxx[k] - g_xy_0_yy_xxxx[k] * ab_y + g_xy_0_yy_xxxxy[k];

                g_xy_0_yyy_xxxy[k] = -g_x_0_yy_xxxy[k] - g_xy_0_yy_xxxy[k] * ab_y + g_xy_0_yy_xxxyy[k];

                g_xy_0_yyy_xxxz[k] = -g_x_0_yy_xxxz[k] - g_xy_0_yy_xxxz[k] * ab_y + g_xy_0_yy_xxxyz[k];

                g_xy_0_yyy_xxyy[k] = -g_x_0_yy_xxyy[k] - g_xy_0_yy_xxyy[k] * ab_y + g_xy_0_yy_xxyyy[k];

                g_xy_0_yyy_xxyz[k] = -g_x_0_yy_xxyz[k] - g_xy_0_yy_xxyz[k] * ab_y + g_xy_0_yy_xxyyz[k];

                g_xy_0_yyy_xxzz[k] = -g_x_0_yy_xxzz[k] - g_xy_0_yy_xxzz[k] * ab_y + g_xy_0_yy_xxyzz[k];

                g_xy_0_yyy_xyyy[k] = -g_x_0_yy_xyyy[k] - g_xy_0_yy_xyyy[k] * ab_y + g_xy_0_yy_xyyyy[k];

                g_xy_0_yyy_xyyz[k] = -g_x_0_yy_xyyz[k] - g_xy_0_yy_xyyz[k] * ab_y + g_xy_0_yy_xyyyz[k];

                g_xy_0_yyy_xyzz[k] = -g_x_0_yy_xyzz[k] - g_xy_0_yy_xyzz[k] * ab_y + g_xy_0_yy_xyyzz[k];

                g_xy_0_yyy_xzzz[k] = -g_x_0_yy_xzzz[k] - g_xy_0_yy_xzzz[k] * ab_y + g_xy_0_yy_xyzzz[k];

                g_xy_0_yyy_yyyy[k] = -g_x_0_yy_yyyy[k] - g_xy_0_yy_yyyy[k] * ab_y + g_xy_0_yy_yyyyy[k];

                g_xy_0_yyy_yyyz[k] = -g_x_0_yy_yyyz[k] - g_xy_0_yy_yyyz[k] * ab_y + g_xy_0_yy_yyyyz[k];

                g_xy_0_yyy_yyzz[k] = -g_x_0_yy_yyzz[k] - g_xy_0_yy_yyzz[k] * ab_y + g_xy_0_yy_yyyzz[k];

                g_xy_0_yyy_yzzz[k] = -g_x_0_yy_yzzz[k] - g_xy_0_yy_yzzz[k] * ab_y + g_xy_0_yy_yyzzz[k];

                g_xy_0_yyy_zzzz[k] = -g_x_0_yy_zzzz[k] - g_xy_0_yy_zzzz[k] * ab_y + g_xy_0_yy_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yy_xxxx, g_xy_0_yy_xxxxz, g_xy_0_yy_xxxy, g_xy_0_yy_xxxyz, g_xy_0_yy_xxxz, g_xy_0_yy_xxxzz, g_xy_0_yy_xxyy, g_xy_0_yy_xxyyz, g_xy_0_yy_xxyz, g_xy_0_yy_xxyzz, g_xy_0_yy_xxzz, g_xy_0_yy_xxzzz, g_xy_0_yy_xyyy, g_xy_0_yy_xyyyz, g_xy_0_yy_xyyz, g_xy_0_yy_xyyzz, g_xy_0_yy_xyzz, g_xy_0_yy_xyzzz, g_xy_0_yy_xzzz, g_xy_0_yy_xzzzz, g_xy_0_yy_yyyy, g_xy_0_yy_yyyyz, g_xy_0_yy_yyyz, g_xy_0_yy_yyyzz, g_xy_0_yy_yyzz, g_xy_0_yy_yyzzz, g_xy_0_yy_yzzz, g_xy_0_yy_yzzzz, g_xy_0_yy_zzzz, g_xy_0_yy_zzzzz, g_xy_0_yyz_xxxx, g_xy_0_yyz_xxxy, g_xy_0_yyz_xxxz, g_xy_0_yyz_xxyy, g_xy_0_yyz_xxyz, g_xy_0_yyz_xxzz, g_xy_0_yyz_xyyy, g_xy_0_yyz_xyyz, g_xy_0_yyz_xyzz, g_xy_0_yyz_xzzz, g_xy_0_yyz_yyyy, g_xy_0_yyz_yyyz, g_xy_0_yyz_yyzz, g_xy_0_yyz_yzzz, g_xy_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyz_xxxx[k] = -g_xy_0_yy_xxxx[k] * ab_z + g_xy_0_yy_xxxxz[k];

                g_xy_0_yyz_xxxy[k] = -g_xy_0_yy_xxxy[k] * ab_z + g_xy_0_yy_xxxyz[k];

                g_xy_0_yyz_xxxz[k] = -g_xy_0_yy_xxxz[k] * ab_z + g_xy_0_yy_xxxzz[k];

                g_xy_0_yyz_xxyy[k] = -g_xy_0_yy_xxyy[k] * ab_z + g_xy_0_yy_xxyyz[k];

                g_xy_0_yyz_xxyz[k] = -g_xy_0_yy_xxyz[k] * ab_z + g_xy_0_yy_xxyzz[k];

                g_xy_0_yyz_xxzz[k] = -g_xy_0_yy_xxzz[k] * ab_z + g_xy_0_yy_xxzzz[k];

                g_xy_0_yyz_xyyy[k] = -g_xy_0_yy_xyyy[k] * ab_z + g_xy_0_yy_xyyyz[k];

                g_xy_0_yyz_xyyz[k] = -g_xy_0_yy_xyyz[k] * ab_z + g_xy_0_yy_xyyzz[k];

                g_xy_0_yyz_xyzz[k] = -g_xy_0_yy_xyzz[k] * ab_z + g_xy_0_yy_xyzzz[k];

                g_xy_0_yyz_xzzz[k] = -g_xy_0_yy_xzzz[k] * ab_z + g_xy_0_yy_xzzzz[k];

                g_xy_0_yyz_yyyy[k] = -g_xy_0_yy_yyyy[k] * ab_z + g_xy_0_yy_yyyyz[k];

                g_xy_0_yyz_yyyz[k] = -g_xy_0_yy_yyyz[k] * ab_z + g_xy_0_yy_yyyzz[k];

                g_xy_0_yyz_yyzz[k] = -g_xy_0_yy_yyzz[k] * ab_z + g_xy_0_yy_yyzzz[k];

                g_xy_0_yyz_yzzz[k] = -g_xy_0_yy_yzzz[k] * ab_z + g_xy_0_yy_yzzzz[k];

                g_xy_0_yyz_zzzz[k] = -g_xy_0_yy_zzzz[k] * ab_z + g_xy_0_yy_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yz_xxxx, g_xy_0_yz_xxxxz, g_xy_0_yz_xxxy, g_xy_0_yz_xxxyz, g_xy_0_yz_xxxz, g_xy_0_yz_xxxzz, g_xy_0_yz_xxyy, g_xy_0_yz_xxyyz, g_xy_0_yz_xxyz, g_xy_0_yz_xxyzz, g_xy_0_yz_xxzz, g_xy_0_yz_xxzzz, g_xy_0_yz_xyyy, g_xy_0_yz_xyyyz, g_xy_0_yz_xyyz, g_xy_0_yz_xyyzz, g_xy_0_yz_xyzz, g_xy_0_yz_xyzzz, g_xy_0_yz_xzzz, g_xy_0_yz_xzzzz, g_xy_0_yz_yyyy, g_xy_0_yz_yyyyz, g_xy_0_yz_yyyz, g_xy_0_yz_yyyzz, g_xy_0_yz_yyzz, g_xy_0_yz_yyzzz, g_xy_0_yz_yzzz, g_xy_0_yz_yzzzz, g_xy_0_yz_zzzz, g_xy_0_yz_zzzzz, g_xy_0_yzz_xxxx, g_xy_0_yzz_xxxy, g_xy_0_yzz_xxxz, g_xy_0_yzz_xxyy, g_xy_0_yzz_xxyz, g_xy_0_yzz_xxzz, g_xy_0_yzz_xyyy, g_xy_0_yzz_xyyz, g_xy_0_yzz_xyzz, g_xy_0_yzz_xzzz, g_xy_0_yzz_yyyy, g_xy_0_yzz_yyyz, g_xy_0_yzz_yyzz, g_xy_0_yzz_yzzz, g_xy_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzz_xxxx[k] = -g_xy_0_yz_xxxx[k] * ab_z + g_xy_0_yz_xxxxz[k];

                g_xy_0_yzz_xxxy[k] = -g_xy_0_yz_xxxy[k] * ab_z + g_xy_0_yz_xxxyz[k];

                g_xy_0_yzz_xxxz[k] = -g_xy_0_yz_xxxz[k] * ab_z + g_xy_0_yz_xxxzz[k];

                g_xy_0_yzz_xxyy[k] = -g_xy_0_yz_xxyy[k] * ab_z + g_xy_0_yz_xxyyz[k];

                g_xy_0_yzz_xxyz[k] = -g_xy_0_yz_xxyz[k] * ab_z + g_xy_0_yz_xxyzz[k];

                g_xy_0_yzz_xxzz[k] = -g_xy_0_yz_xxzz[k] * ab_z + g_xy_0_yz_xxzzz[k];

                g_xy_0_yzz_xyyy[k] = -g_xy_0_yz_xyyy[k] * ab_z + g_xy_0_yz_xyyyz[k];

                g_xy_0_yzz_xyyz[k] = -g_xy_0_yz_xyyz[k] * ab_z + g_xy_0_yz_xyyzz[k];

                g_xy_0_yzz_xyzz[k] = -g_xy_0_yz_xyzz[k] * ab_z + g_xy_0_yz_xyzzz[k];

                g_xy_0_yzz_xzzz[k] = -g_xy_0_yz_xzzz[k] * ab_z + g_xy_0_yz_xzzzz[k];

                g_xy_0_yzz_yyyy[k] = -g_xy_0_yz_yyyy[k] * ab_z + g_xy_0_yz_yyyyz[k];

                g_xy_0_yzz_yyyz[k] = -g_xy_0_yz_yyyz[k] * ab_z + g_xy_0_yz_yyyzz[k];

                g_xy_0_yzz_yyzz[k] = -g_xy_0_yz_yyzz[k] * ab_z + g_xy_0_yz_yyzzz[k];

                g_xy_0_yzz_yzzz[k] = -g_xy_0_yz_yzzz[k] * ab_z + g_xy_0_yz_yzzzz[k];

                g_xy_0_yzz_zzzz[k] = -g_xy_0_yz_zzzz[k] * ab_z + g_xy_0_yz_zzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zz_xxxx, g_xy_0_zz_xxxxz, g_xy_0_zz_xxxy, g_xy_0_zz_xxxyz, g_xy_0_zz_xxxz, g_xy_0_zz_xxxzz, g_xy_0_zz_xxyy, g_xy_0_zz_xxyyz, g_xy_0_zz_xxyz, g_xy_0_zz_xxyzz, g_xy_0_zz_xxzz, g_xy_0_zz_xxzzz, g_xy_0_zz_xyyy, g_xy_0_zz_xyyyz, g_xy_0_zz_xyyz, g_xy_0_zz_xyyzz, g_xy_0_zz_xyzz, g_xy_0_zz_xyzzz, g_xy_0_zz_xzzz, g_xy_0_zz_xzzzz, g_xy_0_zz_yyyy, g_xy_0_zz_yyyyz, g_xy_0_zz_yyyz, g_xy_0_zz_yyyzz, g_xy_0_zz_yyzz, g_xy_0_zz_yyzzz, g_xy_0_zz_yzzz, g_xy_0_zz_yzzzz, g_xy_0_zz_zzzz, g_xy_0_zz_zzzzz, g_xy_0_zzz_xxxx, g_xy_0_zzz_xxxy, g_xy_0_zzz_xxxz, g_xy_0_zzz_xxyy, g_xy_0_zzz_xxyz, g_xy_0_zzz_xxzz, g_xy_0_zzz_xyyy, g_xy_0_zzz_xyyz, g_xy_0_zzz_xyzz, g_xy_0_zzz_xzzz, g_xy_0_zzz_yyyy, g_xy_0_zzz_yyyz, g_xy_0_zzz_yyzz, g_xy_0_zzz_yzzz, g_xy_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzz_xxxx[k] = -g_xy_0_zz_xxxx[k] * ab_z + g_xy_0_zz_xxxxz[k];

                g_xy_0_zzz_xxxy[k] = -g_xy_0_zz_xxxy[k] * ab_z + g_xy_0_zz_xxxyz[k];

                g_xy_0_zzz_xxxz[k] = -g_xy_0_zz_xxxz[k] * ab_z + g_xy_0_zz_xxxzz[k];

                g_xy_0_zzz_xxyy[k] = -g_xy_0_zz_xxyy[k] * ab_z + g_xy_0_zz_xxyyz[k];

                g_xy_0_zzz_xxyz[k] = -g_xy_0_zz_xxyz[k] * ab_z + g_xy_0_zz_xxyzz[k];

                g_xy_0_zzz_xxzz[k] = -g_xy_0_zz_xxzz[k] * ab_z + g_xy_0_zz_xxzzz[k];

                g_xy_0_zzz_xyyy[k] = -g_xy_0_zz_xyyy[k] * ab_z + g_xy_0_zz_xyyyz[k];

                g_xy_0_zzz_xyyz[k] = -g_xy_0_zz_xyyz[k] * ab_z + g_xy_0_zz_xyyzz[k];

                g_xy_0_zzz_xyzz[k] = -g_xy_0_zz_xyzz[k] * ab_z + g_xy_0_zz_xyzzz[k];

                g_xy_0_zzz_xzzz[k] = -g_xy_0_zz_xzzz[k] * ab_z + g_xy_0_zz_xzzzz[k];

                g_xy_0_zzz_yyyy[k] = -g_xy_0_zz_yyyy[k] * ab_z + g_xy_0_zz_yyyyz[k];

                g_xy_0_zzz_yyyz[k] = -g_xy_0_zz_yyyz[k] * ab_z + g_xy_0_zz_yyyzz[k];

                g_xy_0_zzz_yyzz[k] = -g_xy_0_zz_yyzz[k] * ab_z + g_xy_0_zz_yyzzz[k];

                g_xy_0_zzz_yzzz[k] = -g_xy_0_zz_yzzz[k] * ab_z + g_xy_0_zz_yzzzz[k];

                g_xy_0_zzz_zzzz[k] = -g_xy_0_zz_zzzz[k] * ab_z + g_xy_0_zz_zzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_xxxx, g_xz_0_xx_xxxxx, g_xz_0_xx_xxxxy, g_xz_0_xx_xxxxz, g_xz_0_xx_xxxy, g_xz_0_xx_xxxyy, g_xz_0_xx_xxxyz, g_xz_0_xx_xxxz, g_xz_0_xx_xxxzz, g_xz_0_xx_xxyy, g_xz_0_xx_xxyyy, g_xz_0_xx_xxyyz, g_xz_0_xx_xxyz, g_xz_0_xx_xxyzz, g_xz_0_xx_xxzz, g_xz_0_xx_xxzzz, g_xz_0_xx_xyyy, g_xz_0_xx_xyyyy, g_xz_0_xx_xyyyz, g_xz_0_xx_xyyz, g_xz_0_xx_xyyzz, g_xz_0_xx_xyzz, g_xz_0_xx_xyzzz, g_xz_0_xx_xzzz, g_xz_0_xx_xzzzz, g_xz_0_xx_yyyy, g_xz_0_xx_yyyz, g_xz_0_xx_yyzz, g_xz_0_xx_yzzz, g_xz_0_xx_zzzz, g_xz_0_xxx_xxxx, g_xz_0_xxx_xxxy, g_xz_0_xxx_xxxz, g_xz_0_xxx_xxyy, g_xz_0_xxx_xxyz, g_xz_0_xxx_xxzz, g_xz_0_xxx_xyyy, g_xz_0_xxx_xyyz, g_xz_0_xxx_xyzz, g_xz_0_xxx_xzzz, g_xz_0_xxx_yyyy, g_xz_0_xxx_yyyz, g_xz_0_xxx_yyzz, g_xz_0_xxx_yzzz, g_xz_0_xxx_zzzz, g_z_0_xx_xxxx, g_z_0_xx_xxxy, g_z_0_xx_xxxz, g_z_0_xx_xxyy, g_z_0_xx_xxyz, g_z_0_xx_xxzz, g_z_0_xx_xyyy, g_z_0_xx_xyyz, g_z_0_xx_xyzz, g_z_0_xx_xzzz, g_z_0_xx_yyyy, g_z_0_xx_yyyz, g_z_0_xx_yyzz, g_z_0_xx_yzzz, g_z_0_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxx_xxxx[k] = -g_z_0_xx_xxxx[k] - g_xz_0_xx_xxxx[k] * ab_x + g_xz_0_xx_xxxxx[k];

                g_xz_0_xxx_xxxy[k] = -g_z_0_xx_xxxy[k] - g_xz_0_xx_xxxy[k] * ab_x + g_xz_0_xx_xxxxy[k];

                g_xz_0_xxx_xxxz[k] = -g_z_0_xx_xxxz[k] - g_xz_0_xx_xxxz[k] * ab_x + g_xz_0_xx_xxxxz[k];

                g_xz_0_xxx_xxyy[k] = -g_z_0_xx_xxyy[k] - g_xz_0_xx_xxyy[k] * ab_x + g_xz_0_xx_xxxyy[k];

                g_xz_0_xxx_xxyz[k] = -g_z_0_xx_xxyz[k] - g_xz_0_xx_xxyz[k] * ab_x + g_xz_0_xx_xxxyz[k];

                g_xz_0_xxx_xxzz[k] = -g_z_0_xx_xxzz[k] - g_xz_0_xx_xxzz[k] * ab_x + g_xz_0_xx_xxxzz[k];

                g_xz_0_xxx_xyyy[k] = -g_z_0_xx_xyyy[k] - g_xz_0_xx_xyyy[k] * ab_x + g_xz_0_xx_xxyyy[k];

                g_xz_0_xxx_xyyz[k] = -g_z_0_xx_xyyz[k] - g_xz_0_xx_xyyz[k] * ab_x + g_xz_0_xx_xxyyz[k];

                g_xz_0_xxx_xyzz[k] = -g_z_0_xx_xyzz[k] - g_xz_0_xx_xyzz[k] * ab_x + g_xz_0_xx_xxyzz[k];

                g_xz_0_xxx_xzzz[k] = -g_z_0_xx_xzzz[k] - g_xz_0_xx_xzzz[k] * ab_x + g_xz_0_xx_xxzzz[k];

                g_xz_0_xxx_yyyy[k] = -g_z_0_xx_yyyy[k] - g_xz_0_xx_yyyy[k] * ab_x + g_xz_0_xx_xyyyy[k];

                g_xz_0_xxx_yyyz[k] = -g_z_0_xx_yyyz[k] - g_xz_0_xx_yyyz[k] * ab_x + g_xz_0_xx_xyyyz[k];

                g_xz_0_xxx_yyzz[k] = -g_z_0_xx_yyzz[k] - g_xz_0_xx_yyzz[k] * ab_x + g_xz_0_xx_xyyzz[k];

                g_xz_0_xxx_yzzz[k] = -g_z_0_xx_yzzz[k] - g_xz_0_xx_yzzz[k] * ab_x + g_xz_0_xx_xyzzz[k];

                g_xz_0_xxx_zzzz[k] = -g_z_0_xx_zzzz[k] - g_xz_0_xx_zzzz[k] * ab_x + g_xz_0_xx_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_xxxx, g_xz_0_xx_xxxxy, g_xz_0_xx_xxxy, g_xz_0_xx_xxxyy, g_xz_0_xx_xxxyz, g_xz_0_xx_xxxz, g_xz_0_xx_xxyy, g_xz_0_xx_xxyyy, g_xz_0_xx_xxyyz, g_xz_0_xx_xxyz, g_xz_0_xx_xxyzz, g_xz_0_xx_xxzz, g_xz_0_xx_xyyy, g_xz_0_xx_xyyyy, g_xz_0_xx_xyyyz, g_xz_0_xx_xyyz, g_xz_0_xx_xyyzz, g_xz_0_xx_xyzz, g_xz_0_xx_xyzzz, g_xz_0_xx_xzzz, g_xz_0_xx_yyyy, g_xz_0_xx_yyyyy, g_xz_0_xx_yyyyz, g_xz_0_xx_yyyz, g_xz_0_xx_yyyzz, g_xz_0_xx_yyzz, g_xz_0_xx_yyzzz, g_xz_0_xx_yzzz, g_xz_0_xx_yzzzz, g_xz_0_xx_zzzz, g_xz_0_xxy_xxxx, g_xz_0_xxy_xxxy, g_xz_0_xxy_xxxz, g_xz_0_xxy_xxyy, g_xz_0_xxy_xxyz, g_xz_0_xxy_xxzz, g_xz_0_xxy_xyyy, g_xz_0_xxy_xyyz, g_xz_0_xxy_xyzz, g_xz_0_xxy_xzzz, g_xz_0_xxy_yyyy, g_xz_0_xxy_yyyz, g_xz_0_xxy_yyzz, g_xz_0_xxy_yzzz, g_xz_0_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxy_xxxx[k] = -g_xz_0_xx_xxxx[k] * ab_y + g_xz_0_xx_xxxxy[k];

                g_xz_0_xxy_xxxy[k] = -g_xz_0_xx_xxxy[k] * ab_y + g_xz_0_xx_xxxyy[k];

                g_xz_0_xxy_xxxz[k] = -g_xz_0_xx_xxxz[k] * ab_y + g_xz_0_xx_xxxyz[k];

                g_xz_0_xxy_xxyy[k] = -g_xz_0_xx_xxyy[k] * ab_y + g_xz_0_xx_xxyyy[k];

                g_xz_0_xxy_xxyz[k] = -g_xz_0_xx_xxyz[k] * ab_y + g_xz_0_xx_xxyyz[k];

                g_xz_0_xxy_xxzz[k] = -g_xz_0_xx_xxzz[k] * ab_y + g_xz_0_xx_xxyzz[k];

                g_xz_0_xxy_xyyy[k] = -g_xz_0_xx_xyyy[k] * ab_y + g_xz_0_xx_xyyyy[k];

                g_xz_0_xxy_xyyz[k] = -g_xz_0_xx_xyyz[k] * ab_y + g_xz_0_xx_xyyyz[k];

                g_xz_0_xxy_xyzz[k] = -g_xz_0_xx_xyzz[k] * ab_y + g_xz_0_xx_xyyzz[k];

                g_xz_0_xxy_xzzz[k] = -g_xz_0_xx_xzzz[k] * ab_y + g_xz_0_xx_xyzzz[k];

                g_xz_0_xxy_yyyy[k] = -g_xz_0_xx_yyyy[k] * ab_y + g_xz_0_xx_yyyyy[k];

                g_xz_0_xxy_yyyz[k] = -g_xz_0_xx_yyyz[k] * ab_y + g_xz_0_xx_yyyyz[k];

                g_xz_0_xxy_yyzz[k] = -g_xz_0_xx_yyzz[k] * ab_y + g_xz_0_xx_yyyzz[k];

                g_xz_0_xxy_yzzz[k] = -g_xz_0_xx_yzzz[k] * ab_y + g_xz_0_xx_yyzzz[k];

                g_xz_0_xxy_zzzz[k] = -g_xz_0_xx_zzzz[k] * ab_y + g_xz_0_xx_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxz_xxxx, g_xz_0_xxz_xxxy, g_xz_0_xxz_xxxz, g_xz_0_xxz_xxyy, g_xz_0_xxz_xxyz, g_xz_0_xxz_xxzz, g_xz_0_xxz_xyyy, g_xz_0_xxz_xyyz, g_xz_0_xxz_xyzz, g_xz_0_xxz_xzzz, g_xz_0_xxz_yyyy, g_xz_0_xxz_yyyz, g_xz_0_xxz_yyzz, g_xz_0_xxz_yzzz, g_xz_0_xxz_zzzz, g_xz_0_xz_xxxx, g_xz_0_xz_xxxxx, g_xz_0_xz_xxxxy, g_xz_0_xz_xxxxz, g_xz_0_xz_xxxy, g_xz_0_xz_xxxyy, g_xz_0_xz_xxxyz, g_xz_0_xz_xxxz, g_xz_0_xz_xxxzz, g_xz_0_xz_xxyy, g_xz_0_xz_xxyyy, g_xz_0_xz_xxyyz, g_xz_0_xz_xxyz, g_xz_0_xz_xxyzz, g_xz_0_xz_xxzz, g_xz_0_xz_xxzzz, g_xz_0_xz_xyyy, g_xz_0_xz_xyyyy, g_xz_0_xz_xyyyz, g_xz_0_xz_xyyz, g_xz_0_xz_xyyzz, g_xz_0_xz_xyzz, g_xz_0_xz_xyzzz, g_xz_0_xz_xzzz, g_xz_0_xz_xzzzz, g_xz_0_xz_yyyy, g_xz_0_xz_yyyz, g_xz_0_xz_yyzz, g_xz_0_xz_yzzz, g_xz_0_xz_zzzz, g_z_0_xz_xxxx, g_z_0_xz_xxxy, g_z_0_xz_xxxz, g_z_0_xz_xxyy, g_z_0_xz_xxyz, g_z_0_xz_xxzz, g_z_0_xz_xyyy, g_z_0_xz_xyyz, g_z_0_xz_xyzz, g_z_0_xz_xzzz, g_z_0_xz_yyyy, g_z_0_xz_yyyz, g_z_0_xz_yyzz, g_z_0_xz_yzzz, g_z_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxz_xxxx[k] = -g_z_0_xz_xxxx[k] - g_xz_0_xz_xxxx[k] * ab_x + g_xz_0_xz_xxxxx[k];

                g_xz_0_xxz_xxxy[k] = -g_z_0_xz_xxxy[k] - g_xz_0_xz_xxxy[k] * ab_x + g_xz_0_xz_xxxxy[k];

                g_xz_0_xxz_xxxz[k] = -g_z_0_xz_xxxz[k] - g_xz_0_xz_xxxz[k] * ab_x + g_xz_0_xz_xxxxz[k];

                g_xz_0_xxz_xxyy[k] = -g_z_0_xz_xxyy[k] - g_xz_0_xz_xxyy[k] * ab_x + g_xz_0_xz_xxxyy[k];

                g_xz_0_xxz_xxyz[k] = -g_z_0_xz_xxyz[k] - g_xz_0_xz_xxyz[k] * ab_x + g_xz_0_xz_xxxyz[k];

                g_xz_0_xxz_xxzz[k] = -g_z_0_xz_xxzz[k] - g_xz_0_xz_xxzz[k] * ab_x + g_xz_0_xz_xxxzz[k];

                g_xz_0_xxz_xyyy[k] = -g_z_0_xz_xyyy[k] - g_xz_0_xz_xyyy[k] * ab_x + g_xz_0_xz_xxyyy[k];

                g_xz_0_xxz_xyyz[k] = -g_z_0_xz_xyyz[k] - g_xz_0_xz_xyyz[k] * ab_x + g_xz_0_xz_xxyyz[k];

                g_xz_0_xxz_xyzz[k] = -g_z_0_xz_xyzz[k] - g_xz_0_xz_xyzz[k] * ab_x + g_xz_0_xz_xxyzz[k];

                g_xz_0_xxz_xzzz[k] = -g_z_0_xz_xzzz[k] - g_xz_0_xz_xzzz[k] * ab_x + g_xz_0_xz_xxzzz[k];

                g_xz_0_xxz_yyyy[k] = -g_z_0_xz_yyyy[k] - g_xz_0_xz_yyyy[k] * ab_x + g_xz_0_xz_xyyyy[k];

                g_xz_0_xxz_yyyz[k] = -g_z_0_xz_yyyz[k] - g_xz_0_xz_yyyz[k] * ab_x + g_xz_0_xz_xyyyz[k];

                g_xz_0_xxz_yyzz[k] = -g_z_0_xz_yyzz[k] - g_xz_0_xz_yyzz[k] * ab_x + g_xz_0_xz_xyyzz[k];

                g_xz_0_xxz_yzzz[k] = -g_z_0_xz_yzzz[k] - g_xz_0_xz_yzzz[k] * ab_x + g_xz_0_xz_xyzzz[k];

                g_xz_0_xxz_zzzz[k] = -g_z_0_xz_zzzz[k] - g_xz_0_xz_zzzz[k] * ab_x + g_xz_0_xz_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xy_xxxx, g_xz_0_xy_xxxxy, g_xz_0_xy_xxxy, g_xz_0_xy_xxxyy, g_xz_0_xy_xxxyz, g_xz_0_xy_xxxz, g_xz_0_xy_xxyy, g_xz_0_xy_xxyyy, g_xz_0_xy_xxyyz, g_xz_0_xy_xxyz, g_xz_0_xy_xxyzz, g_xz_0_xy_xxzz, g_xz_0_xy_xyyy, g_xz_0_xy_xyyyy, g_xz_0_xy_xyyyz, g_xz_0_xy_xyyz, g_xz_0_xy_xyyzz, g_xz_0_xy_xyzz, g_xz_0_xy_xyzzz, g_xz_0_xy_xzzz, g_xz_0_xy_yyyy, g_xz_0_xy_yyyyy, g_xz_0_xy_yyyyz, g_xz_0_xy_yyyz, g_xz_0_xy_yyyzz, g_xz_0_xy_yyzz, g_xz_0_xy_yyzzz, g_xz_0_xy_yzzz, g_xz_0_xy_yzzzz, g_xz_0_xy_zzzz, g_xz_0_xyy_xxxx, g_xz_0_xyy_xxxy, g_xz_0_xyy_xxxz, g_xz_0_xyy_xxyy, g_xz_0_xyy_xxyz, g_xz_0_xyy_xxzz, g_xz_0_xyy_xyyy, g_xz_0_xyy_xyyz, g_xz_0_xyy_xyzz, g_xz_0_xyy_xzzz, g_xz_0_xyy_yyyy, g_xz_0_xyy_yyyz, g_xz_0_xyy_yyzz, g_xz_0_xyy_yzzz, g_xz_0_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyy_xxxx[k] = -g_xz_0_xy_xxxx[k] * ab_y + g_xz_0_xy_xxxxy[k];

                g_xz_0_xyy_xxxy[k] = -g_xz_0_xy_xxxy[k] * ab_y + g_xz_0_xy_xxxyy[k];

                g_xz_0_xyy_xxxz[k] = -g_xz_0_xy_xxxz[k] * ab_y + g_xz_0_xy_xxxyz[k];

                g_xz_0_xyy_xxyy[k] = -g_xz_0_xy_xxyy[k] * ab_y + g_xz_0_xy_xxyyy[k];

                g_xz_0_xyy_xxyz[k] = -g_xz_0_xy_xxyz[k] * ab_y + g_xz_0_xy_xxyyz[k];

                g_xz_0_xyy_xxzz[k] = -g_xz_0_xy_xxzz[k] * ab_y + g_xz_0_xy_xxyzz[k];

                g_xz_0_xyy_xyyy[k] = -g_xz_0_xy_xyyy[k] * ab_y + g_xz_0_xy_xyyyy[k];

                g_xz_0_xyy_xyyz[k] = -g_xz_0_xy_xyyz[k] * ab_y + g_xz_0_xy_xyyyz[k];

                g_xz_0_xyy_xyzz[k] = -g_xz_0_xy_xyzz[k] * ab_y + g_xz_0_xy_xyyzz[k];

                g_xz_0_xyy_xzzz[k] = -g_xz_0_xy_xzzz[k] * ab_y + g_xz_0_xy_xyzzz[k];

                g_xz_0_xyy_yyyy[k] = -g_xz_0_xy_yyyy[k] * ab_y + g_xz_0_xy_yyyyy[k];

                g_xz_0_xyy_yyyz[k] = -g_xz_0_xy_yyyz[k] * ab_y + g_xz_0_xy_yyyyz[k];

                g_xz_0_xyy_yyzz[k] = -g_xz_0_xy_yyzz[k] * ab_y + g_xz_0_xy_yyyzz[k];

                g_xz_0_xyy_yzzz[k] = -g_xz_0_xy_yzzz[k] * ab_y + g_xz_0_xy_yyzzz[k];

                g_xz_0_xyy_zzzz[k] = -g_xz_0_xy_zzzz[k] * ab_y + g_xz_0_xy_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyz_xxxx, g_xz_0_xyz_xxxy, g_xz_0_xyz_xxxz, g_xz_0_xyz_xxyy, g_xz_0_xyz_xxyz, g_xz_0_xyz_xxzz, g_xz_0_xyz_xyyy, g_xz_0_xyz_xyyz, g_xz_0_xyz_xyzz, g_xz_0_xyz_xzzz, g_xz_0_xyz_yyyy, g_xz_0_xyz_yyyz, g_xz_0_xyz_yyzz, g_xz_0_xyz_yzzz, g_xz_0_xyz_zzzz, g_xz_0_xz_xxxx, g_xz_0_xz_xxxxy, g_xz_0_xz_xxxy, g_xz_0_xz_xxxyy, g_xz_0_xz_xxxyz, g_xz_0_xz_xxxz, g_xz_0_xz_xxyy, g_xz_0_xz_xxyyy, g_xz_0_xz_xxyyz, g_xz_0_xz_xxyz, g_xz_0_xz_xxyzz, g_xz_0_xz_xxzz, g_xz_0_xz_xyyy, g_xz_0_xz_xyyyy, g_xz_0_xz_xyyyz, g_xz_0_xz_xyyz, g_xz_0_xz_xyyzz, g_xz_0_xz_xyzz, g_xz_0_xz_xyzzz, g_xz_0_xz_xzzz, g_xz_0_xz_yyyy, g_xz_0_xz_yyyyy, g_xz_0_xz_yyyyz, g_xz_0_xz_yyyz, g_xz_0_xz_yyyzz, g_xz_0_xz_yyzz, g_xz_0_xz_yyzzz, g_xz_0_xz_yzzz, g_xz_0_xz_yzzzz, g_xz_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyz_xxxx[k] = -g_xz_0_xz_xxxx[k] * ab_y + g_xz_0_xz_xxxxy[k];

                g_xz_0_xyz_xxxy[k] = -g_xz_0_xz_xxxy[k] * ab_y + g_xz_0_xz_xxxyy[k];

                g_xz_0_xyz_xxxz[k] = -g_xz_0_xz_xxxz[k] * ab_y + g_xz_0_xz_xxxyz[k];

                g_xz_0_xyz_xxyy[k] = -g_xz_0_xz_xxyy[k] * ab_y + g_xz_0_xz_xxyyy[k];

                g_xz_0_xyz_xxyz[k] = -g_xz_0_xz_xxyz[k] * ab_y + g_xz_0_xz_xxyyz[k];

                g_xz_0_xyz_xxzz[k] = -g_xz_0_xz_xxzz[k] * ab_y + g_xz_0_xz_xxyzz[k];

                g_xz_0_xyz_xyyy[k] = -g_xz_0_xz_xyyy[k] * ab_y + g_xz_0_xz_xyyyy[k];

                g_xz_0_xyz_xyyz[k] = -g_xz_0_xz_xyyz[k] * ab_y + g_xz_0_xz_xyyyz[k];

                g_xz_0_xyz_xyzz[k] = -g_xz_0_xz_xyzz[k] * ab_y + g_xz_0_xz_xyyzz[k];

                g_xz_0_xyz_xzzz[k] = -g_xz_0_xz_xzzz[k] * ab_y + g_xz_0_xz_xyzzz[k];

                g_xz_0_xyz_yyyy[k] = -g_xz_0_xz_yyyy[k] * ab_y + g_xz_0_xz_yyyyy[k];

                g_xz_0_xyz_yyyz[k] = -g_xz_0_xz_yyyz[k] * ab_y + g_xz_0_xz_yyyyz[k];

                g_xz_0_xyz_yyzz[k] = -g_xz_0_xz_yyzz[k] * ab_y + g_xz_0_xz_yyyzz[k];

                g_xz_0_xyz_yzzz[k] = -g_xz_0_xz_yzzz[k] * ab_y + g_xz_0_xz_yyzzz[k];

                g_xz_0_xyz_zzzz[k] = -g_xz_0_xz_zzzz[k] * ab_y + g_xz_0_xz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 377 * ccomps * dcomps);

            auto g_xz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 379 * ccomps * dcomps);

            auto g_xz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 383 * ccomps * dcomps);

            auto g_xz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzz_xxxx, g_xz_0_xzz_xxxy, g_xz_0_xzz_xxxz, g_xz_0_xzz_xxyy, g_xz_0_xzz_xxyz, g_xz_0_xzz_xxzz, g_xz_0_xzz_xyyy, g_xz_0_xzz_xyyz, g_xz_0_xzz_xyzz, g_xz_0_xzz_xzzz, g_xz_0_xzz_yyyy, g_xz_0_xzz_yyyz, g_xz_0_xzz_yyzz, g_xz_0_xzz_yzzz, g_xz_0_xzz_zzzz, g_xz_0_zz_xxxx, g_xz_0_zz_xxxxx, g_xz_0_zz_xxxxy, g_xz_0_zz_xxxxz, g_xz_0_zz_xxxy, g_xz_0_zz_xxxyy, g_xz_0_zz_xxxyz, g_xz_0_zz_xxxz, g_xz_0_zz_xxxzz, g_xz_0_zz_xxyy, g_xz_0_zz_xxyyy, g_xz_0_zz_xxyyz, g_xz_0_zz_xxyz, g_xz_0_zz_xxyzz, g_xz_0_zz_xxzz, g_xz_0_zz_xxzzz, g_xz_0_zz_xyyy, g_xz_0_zz_xyyyy, g_xz_0_zz_xyyyz, g_xz_0_zz_xyyz, g_xz_0_zz_xyyzz, g_xz_0_zz_xyzz, g_xz_0_zz_xyzzz, g_xz_0_zz_xzzz, g_xz_0_zz_xzzzz, g_xz_0_zz_yyyy, g_xz_0_zz_yyyz, g_xz_0_zz_yyzz, g_xz_0_zz_yzzz, g_xz_0_zz_zzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzz_xxxx[k] = -g_z_0_zz_xxxx[k] - g_xz_0_zz_xxxx[k] * ab_x + g_xz_0_zz_xxxxx[k];

                g_xz_0_xzz_xxxy[k] = -g_z_0_zz_xxxy[k] - g_xz_0_zz_xxxy[k] * ab_x + g_xz_0_zz_xxxxy[k];

                g_xz_0_xzz_xxxz[k] = -g_z_0_zz_xxxz[k] - g_xz_0_zz_xxxz[k] * ab_x + g_xz_0_zz_xxxxz[k];

                g_xz_0_xzz_xxyy[k] = -g_z_0_zz_xxyy[k] - g_xz_0_zz_xxyy[k] * ab_x + g_xz_0_zz_xxxyy[k];

                g_xz_0_xzz_xxyz[k] = -g_z_0_zz_xxyz[k] - g_xz_0_zz_xxyz[k] * ab_x + g_xz_0_zz_xxxyz[k];

                g_xz_0_xzz_xxzz[k] = -g_z_0_zz_xxzz[k] - g_xz_0_zz_xxzz[k] * ab_x + g_xz_0_zz_xxxzz[k];

                g_xz_0_xzz_xyyy[k] = -g_z_0_zz_xyyy[k] - g_xz_0_zz_xyyy[k] * ab_x + g_xz_0_zz_xxyyy[k];

                g_xz_0_xzz_xyyz[k] = -g_z_0_zz_xyyz[k] - g_xz_0_zz_xyyz[k] * ab_x + g_xz_0_zz_xxyyz[k];

                g_xz_0_xzz_xyzz[k] = -g_z_0_zz_xyzz[k] - g_xz_0_zz_xyzz[k] * ab_x + g_xz_0_zz_xxyzz[k];

                g_xz_0_xzz_xzzz[k] = -g_z_0_zz_xzzz[k] - g_xz_0_zz_xzzz[k] * ab_x + g_xz_0_zz_xxzzz[k];

                g_xz_0_xzz_yyyy[k] = -g_z_0_zz_yyyy[k] - g_xz_0_zz_yyyy[k] * ab_x + g_xz_0_zz_xyyyy[k];

                g_xz_0_xzz_yyyz[k] = -g_z_0_zz_yyyz[k] - g_xz_0_zz_yyyz[k] * ab_x + g_xz_0_zz_xyyyz[k];

                g_xz_0_xzz_yyzz[k] = -g_z_0_zz_yyzz[k] - g_xz_0_zz_yyzz[k] * ab_x + g_xz_0_zz_xyyzz[k];

                g_xz_0_xzz_yzzz[k] = -g_z_0_zz_yzzz[k] - g_xz_0_zz_yzzz[k] * ab_x + g_xz_0_zz_xyzzz[k];

                g_xz_0_xzz_zzzz[k] = -g_z_0_zz_zzzz[k] - g_xz_0_zz_zzzz[k] * ab_x + g_xz_0_zz_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 391 * ccomps * dcomps);

            auto g_xz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 395 * ccomps * dcomps);

            auto g_xz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 399 * ccomps * dcomps);

            auto g_xz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 401 * ccomps * dcomps);

            auto g_xz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yy_xxxx, g_xz_0_yy_xxxxy, g_xz_0_yy_xxxy, g_xz_0_yy_xxxyy, g_xz_0_yy_xxxyz, g_xz_0_yy_xxxz, g_xz_0_yy_xxyy, g_xz_0_yy_xxyyy, g_xz_0_yy_xxyyz, g_xz_0_yy_xxyz, g_xz_0_yy_xxyzz, g_xz_0_yy_xxzz, g_xz_0_yy_xyyy, g_xz_0_yy_xyyyy, g_xz_0_yy_xyyyz, g_xz_0_yy_xyyz, g_xz_0_yy_xyyzz, g_xz_0_yy_xyzz, g_xz_0_yy_xyzzz, g_xz_0_yy_xzzz, g_xz_0_yy_yyyy, g_xz_0_yy_yyyyy, g_xz_0_yy_yyyyz, g_xz_0_yy_yyyz, g_xz_0_yy_yyyzz, g_xz_0_yy_yyzz, g_xz_0_yy_yyzzz, g_xz_0_yy_yzzz, g_xz_0_yy_yzzzz, g_xz_0_yy_zzzz, g_xz_0_yyy_xxxx, g_xz_0_yyy_xxxy, g_xz_0_yyy_xxxz, g_xz_0_yyy_xxyy, g_xz_0_yyy_xxyz, g_xz_0_yyy_xxzz, g_xz_0_yyy_xyyy, g_xz_0_yyy_xyyz, g_xz_0_yyy_xyzz, g_xz_0_yyy_xzzz, g_xz_0_yyy_yyyy, g_xz_0_yyy_yyyz, g_xz_0_yyy_yyzz, g_xz_0_yyy_yzzz, g_xz_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyy_xxxx[k] = -g_xz_0_yy_xxxx[k] * ab_y + g_xz_0_yy_xxxxy[k];

                g_xz_0_yyy_xxxy[k] = -g_xz_0_yy_xxxy[k] * ab_y + g_xz_0_yy_xxxyy[k];

                g_xz_0_yyy_xxxz[k] = -g_xz_0_yy_xxxz[k] * ab_y + g_xz_0_yy_xxxyz[k];

                g_xz_0_yyy_xxyy[k] = -g_xz_0_yy_xxyy[k] * ab_y + g_xz_0_yy_xxyyy[k];

                g_xz_0_yyy_xxyz[k] = -g_xz_0_yy_xxyz[k] * ab_y + g_xz_0_yy_xxyyz[k];

                g_xz_0_yyy_xxzz[k] = -g_xz_0_yy_xxzz[k] * ab_y + g_xz_0_yy_xxyzz[k];

                g_xz_0_yyy_xyyy[k] = -g_xz_0_yy_xyyy[k] * ab_y + g_xz_0_yy_xyyyy[k];

                g_xz_0_yyy_xyyz[k] = -g_xz_0_yy_xyyz[k] * ab_y + g_xz_0_yy_xyyyz[k];

                g_xz_0_yyy_xyzz[k] = -g_xz_0_yy_xyzz[k] * ab_y + g_xz_0_yy_xyyzz[k];

                g_xz_0_yyy_xzzz[k] = -g_xz_0_yy_xzzz[k] * ab_y + g_xz_0_yy_xyzzz[k];

                g_xz_0_yyy_yyyy[k] = -g_xz_0_yy_yyyy[k] * ab_y + g_xz_0_yy_yyyyy[k];

                g_xz_0_yyy_yyyz[k] = -g_xz_0_yy_yyyz[k] * ab_y + g_xz_0_yy_yyyyz[k];

                g_xz_0_yyy_yyzz[k] = -g_xz_0_yy_yyzz[k] * ab_y + g_xz_0_yy_yyyzz[k];

                g_xz_0_yyy_yzzz[k] = -g_xz_0_yy_yzzz[k] * ab_y + g_xz_0_yy_yyzzz[k];

                g_xz_0_yyy_zzzz[k] = -g_xz_0_yy_zzzz[k] * ab_y + g_xz_0_yy_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 407 * ccomps * dcomps);

            auto g_xz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 409 * ccomps * dcomps);

            auto g_xz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 413 * ccomps * dcomps);

            auto g_xz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyz_xxxx, g_xz_0_yyz_xxxy, g_xz_0_yyz_xxxz, g_xz_0_yyz_xxyy, g_xz_0_yyz_xxyz, g_xz_0_yyz_xxzz, g_xz_0_yyz_xyyy, g_xz_0_yyz_xyyz, g_xz_0_yyz_xyzz, g_xz_0_yyz_xzzz, g_xz_0_yyz_yyyy, g_xz_0_yyz_yyyz, g_xz_0_yyz_yyzz, g_xz_0_yyz_yzzz, g_xz_0_yyz_zzzz, g_xz_0_yz_xxxx, g_xz_0_yz_xxxxy, g_xz_0_yz_xxxy, g_xz_0_yz_xxxyy, g_xz_0_yz_xxxyz, g_xz_0_yz_xxxz, g_xz_0_yz_xxyy, g_xz_0_yz_xxyyy, g_xz_0_yz_xxyyz, g_xz_0_yz_xxyz, g_xz_0_yz_xxyzz, g_xz_0_yz_xxzz, g_xz_0_yz_xyyy, g_xz_0_yz_xyyyy, g_xz_0_yz_xyyyz, g_xz_0_yz_xyyz, g_xz_0_yz_xyyzz, g_xz_0_yz_xyzz, g_xz_0_yz_xyzzz, g_xz_0_yz_xzzz, g_xz_0_yz_yyyy, g_xz_0_yz_yyyyy, g_xz_0_yz_yyyyz, g_xz_0_yz_yyyz, g_xz_0_yz_yyyzz, g_xz_0_yz_yyzz, g_xz_0_yz_yyzzz, g_xz_0_yz_yzzz, g_xz_0_yz_yzzzz, g_xz_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyz_xxxx[k] = -g_xz_0_yz_xxxx[k] * ab_y + g_xz_0_yz_xxxxy[k];

                g_xz_0_yyz_xxxy[k] = -g_xz_0_yz_xxxy[k] * ab_y + g_xz_0_yz_xxxyy[k];

                g_xz_0_yyz_xxxz[k] = -g_xz_0_yz_xxxz[k] * ab_y + g_xz_0_yz_xxxyz[k];

                g_xz_0_yyz_xxyy[k] = -g_xz_0_yz_xxyy[k] * ab_y + g_xz_0_yz_xxyyy[k];

                g_xz_0_yyz_xxyz[k] = -g_xz_0_yz_xxyz[k] * ab_y + g_xz_0_yz_xxyyz[k];

                g_xz_0_yyz_xxzz[k] = -g_xz_0_yz_xxzz[k] * ab_y + g_xz_0_yz_xxyzz[k];

                g_xz_0_yyz_xyyy[k] = -g_xz_0_yz_xyyy[k] * ab_y + g_xz_0_yz_xyyyy[k];

                g_xz_0_yyz_xyyz[k] = -g_xz_0_yz_xyyz[k] * ab_y + g_xz_0_yz_xyyyz[k];

                g_xz_0_yyz_xyzz[k] = -g_xz_0_yz_xyzz[k] * ab_y + g_xz_0_yz_xyyzz[k];

                g_xz_0_yyz_xzzz[k] = -g_xz_0_yz_xzzz[k] * ab_y + g_xz_0_yz_xyzzz[k];

                g_xz_0_yyz_yyyy[k] = -g_xz_0_yz_yyyy[k] * ab_y + g_xz_0_yz_yyyyy[k];

                g_xz_0_yyz_yyyz[k] = -g_xz_0_yz_yyyz[k] * ab_y + g_xz_0_yz_yyyyz[k];

                g_xz_0_yyz_yyzz[k] = -g_xz_0_yz_yyzz[k] * ab_y + g_xz_0_yz_yyyzz[k];

                g_xz_0_yyz_yzzz[k] = -g_xz_0_yz_yzzz[k] * ab_y + g_xz_0_yz_yyzzz[k];

                g_xz_0_yyz_zzzz[k] = -g_xz_0_yz_zzzz[k] * ab_y + g_xz_0_yz_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzz_xxxx, g_xz_0_yzz_xxxy, g_xz_0_yzz_xxxz, g_xz_0_yzz_xxyy, g_xz_0_yzz_xxyz, g_xz_0_yzz_xxzz, g_xz_0_yzz_xyyy, g_xz_0_yzz_xyyz, g_xz_0_yzz_xyzz, g_xz_0_yzz_xzzz, g_xz_0_yzz_yyyy, g_xz_0_yzz_yyyz, g_xz_0_yzz_yyzz, g_xz_0_yzz_yzzz, g_xz_0_yzz_zzzz, g_xz_0_zz_xxxx, g_xz_0_zz_xxxxy, g_xz_0_zz_xxxy, g_xz_0_zz_xxxyy, g_xz_0_zz_xxxyz, g_xz_0_zz_xxxz, g_xz_0_zz_xxyy, g_xz_0_zz_xxyyy, g_xz_0_zz_xxyyz, g_xz_0_zz_xxyz, g_xz_0_zz_xxyzz, g_xz_0_zz_xxzz, g_xz_0_zz_xyyy, g_xz_0_zz_xyyyy, g_xz_0_zz_xyyyz, g_xz_0_zz_xyyz, g_xz_0_zz_xyyzz, g_xz_0_zz_xyzz, g_xz_0_zz_xyzzz, g_xz_0_zz_xzzz, g_xz_0_zz_yyyy, g_xz_0_zz_yyyyy, g_xz_0_zz_yyyyz, g_xz_0_zz_yyyz, g_xz_0_zz_yyyzz, g_xz_0_zz_yyzz, g_xz_0_zz_yyzzz, g_xz_0_zz_yzzz, g_xz_0_zz_yzzzz, g_xz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzz_xxxx[k] = -g_xz_0_zz_xxxx[k] * ab_y + g_xz_0_zz_xxxxy[k];

                g_xz_0_yzz_xxxy[k] = -g_xz_0_zz_xxxy[k] * ab_y + g_xz_0_zz_xxxyy[k];

                g_xz_0_yzz_xxxz[k] = -g_xz_0_zz_xxxz[k] * ab_y + g_xz_0_zz_xxxyz[k];

                g_xz_0_yzz_xxyy[k] = -g_xz_0_zz_xxyy[k] * ab_y + g_xz_0_zz_xxyyy[k];

                g_xz_0_yzz_xxyz[k] = -g_xz_0_zz_xxyz[k] * ab_y + g_xz_0_zz_xxyyz[k];

                g_xz_0_yzz_xxzz[k] = -g_xz_0_zz_xxzz[k] * ab_y + g_xz_0_zz_xxyzz[k];

                g_xz_0_yzz_xyyy[k] = -g_xz_0_zz_xyyy[k] * ab_y + g_xz_0_zz_xyyyy[k];

                g_xz_0_yzz_xyyz[k] = -g_xz_0_zz_xyyz[k] * ab_y + g_xz_0_zz_xyyyz[k];

                g_xz_0_yzz_xyzz[k] = -g_xz_0_zz_xyzz[k] * ab_y + g_xz_0_zz_xyyzz[k];

                g_xz_0_yzz_xzzz[k] = -g_xz_0_zz_xzzz[k] * ab_y + g_xz_0_zz_xyzzz[k];

                g_xz_0_yzz_yyyy[k] = -g_xz_0_zz_yyyy[k] * ab_y + g_xz_0_zz_yyyyy[k];

                g_xz_0_yzz_yyyz[k] = -g_xz_0_zz_yyyz[k] * ab_y + g_xz_0_zz_yyyyz[k];

                g_xz_0_yzz_yyzz[k] = -g_xz_0_zz_yyzz[k] * ab_y + g_xz_0_zz_yyyzz[k];

                g_xz_0_yzz_yzzz[k] = -g_xz_0_zz_yzzz[k] * ab_y + g_xz_0_zz_yyzzz[k];

                g_xz_0_yzz_zzzz[k] = -g_xz_0_zz_zzzz[k] * ab_y + g_xz_0_zz_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xxxx, g_x_0_zz_xxxy, g_x_0_zz_xxxz, g_x_0_zz_xxyy, g_x_0_zz_xxyz, g_x_0_zz_xxzz, g_x_0_zz_xyyy, g_x_0_zz_xyyz, g_x_0_zz_xyzz, g_x_0_zz_xzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyz, g_x_0_zz_yyzz, g_x_0_zz_yzzz, g_x_0_zz_zzzz, g_xz_0_zz_xxxx, g_xz_0_zz_xxxxz, g_xz_0_zz_xxxy, g_xz_0_zz_xxxyz, g_xz_0_zz_xxxz, g_xz_0_zz_xxxzz, g_xz_0_zz_xxyy, g_xz_0_zz_xxyyz, g_xz_0_zz_xxyz, g_xz_0_zz_xxyzz, g_xz_0_zz_xxzz, g_xz_0_zz_xxzzz, g_xz_0_zz_xyyy, g_xz_0_zz_xyyyz, g_xz_0_zz_xyyz, g_xz_0_zz_xyyzz, g_xz_0_zz_xyzz, g_xz_0_zz_xyzzz, g_xz_0_zz_xzzz, g_xz_0_zz_xzzzz, g_xz_0_zz_yyyy, g_xz_0_zz_yyyyz, g_xz_0_zz_yyyz, g_xz_0_zz_yyyzz, g_xz_0_zz_yyzz, g_xz_0_zz_yyzzz, g_xz_0_zz_yzzz, g_xz_0_zz_yzzzz, g_xz_0_zz_zzzz, g_xz_0_zz_zzzzz, g_xz_0_zzz_xxxx, g_xz_0_zzz_xxxy, g_xz_0_zzz_xxxz, g_xz_0_zzz_xxyy, g_xz_0_zzz_xxyz, g_xz_0_zzz_xxzz, g_xz_0_zzz_xyyy, g_xz_0_zzz_xyyz, g_xz_0_zzz_xyzz, g_xz_0_zzz_xzzz, g_xz_0_zzz_yyyy, g_xz_0_zzz_yyyz, g_xz_0_zzz_yyzz, g_xz_0_zzz_yzzz, g_xz_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzz_xxxx[k] = -g_x_0_zz_xxxx[k] - g_xz_0_zz_xxxx[k] * ab_z + g_xz_0_zz_xxxxz[k];

                g_xz_0_zzz_xxxy[k] = -g_x_0_zz_xxxy[k] - g_xz_0_zz_xxxy[k] * ab_z + g_xz_0_zz_xxxyz[k];

                g_xz_0_zzz_xxxz[k] = -g_x_0_zz_xxxz[k] - g_xz_0_zz_xxxz[k] * ab_z + g_xz_0_zz_xxxzz[k];

                g_xz_0_zzz_xxyy[k] = -g_x_0_zz_xxyy[k] - g_xz_0_zz_xxyy[k] * ab_z + g_xz_0_zz_xxyyz[k];

                g_xz_0_zzz_xxyz[k] = -g_x_0_zz_xxyz[k] - g_xz_0_zz_xxyz[k] * ab_z + g_xz_0_zz_xxyzz[k];

                g_xz_0_zzz_xxzz[k] = -g_x_0_zz_xxzz[k] - g_xz_0_zz_xxzz[k] * ab_z + g_xz_0_zz_xxzzz[k];

                g_xz_0_zzz_xyyy[k] = -g_x_0_zz_xyyy[k] - g_xz_0_zz_xyyy[k] * ab_z + g_xz_0_zz_xyyyz[k];

                g_xz_0_zzz_xyyz[k] = -g_x_0_zz_xyyz[k] - g_xz_0_zz_xyyz[k] * ab_z + g_xz_0_zz_xyyzz[k];

                g_xz_0_zzz_xyzz[k] = -g_x_0_zz_xyzz[k] - g_xz_0_zz_xyzz[k] * ab_z + g_xz_0_zz_xyzzz[k];

                g_xz_0_zzz_xzzz[k] = -g_x_0_zz_xzzz[k] - g_xz_0_zz_xzzz[k] * ab_z + g_xz_0_zz_xzzzz[k];

                g_xz_0_zzz_yyyy[k] = -g_x_0_zz_yyyy[k] - g_xz_0_zz_yyyy[k] * ab_z + g_xz_0_zz_yyyyz[k];

                g_xz_0_zzz_yyyz[k] = -g_x_0_zz_yyyz[k] - g_xz_0_zz_yyyz[k] * ab_z + g_xz_0_zz_yyyzz[k];

                g_xz_0_zzz_yyzz[k] = -g_x_0_zz_yyzz[k] - g_xz_0_zz_yyzz[k] * ab_z + g_xz_0_zz_yyzzz[k];

                g_xz_0_zzz_yzzz[k] = -g_x_0_zz_yzzz[k] - g_xz_0_zz_yzzz[k] * ab_z + g_xz_0_zz_yzzzz[k];

                g_xz_0_zzz_zzzz[k] = -g_x_0_zz_zzzz[k] - g_xz_0_zz_zzzz[k] * ab_z + g_xz_0_zz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xx_xxxx, g_yy_0_xx_xxxxx, g_yy_0_xx_xxxxy, g_yy_0_xx_xxxxz, g_yy_0_xx_xxxy, g_yy_0_xx_xxxyy, g_yy_0_xx_xxxyz, g_yy_0_xx_xxxz, g_yy_0_xx_xxxzz, g_yy_0_xx_xxyy, g_yy_0_xx_xxyyy, g_yy_0_xx_xxyyz, g_yy_0_xx_xxyz, g_yy_0_xx_xxyzz, g_yy_0_xx_xxzz, g_yy_0_xx_xxzzz, g_yy_0_xx_xyyy, g_yy_0_xx_xyyyy, g_yy_0_xx_xyyyz, g_yy_0_xx_xyyz, g_yy_0_xx_xyyzz, g_yy_0_xx_xyzz, g_yy_0_xx_xyzzz, g_yy_0_xx_xzzz, g_yy_0_xx_xzzzz, g_yy_0_xx_yyyy, g_yy_0_xx_yyyz, g_yy_0_xx_yyzz, g_yy_0_xx_yzzz, g_yy_0_xx_zzzz, g_yy_0_xxx_xxxx, g_yy_0_xxx_xxxy, g_yy_0_xxx_xxxz, g_yy_0_xxx_xxyy, g_yy_0_xxx_xxyz, g_yy_0_xxx_xxzz, g_yy_0_xxx_xyyy, g_yy_0_xxx_xyyz, g_yy_0_xxx_xyzz, g_yy_0_xxx_xzzz, g_yy_0_xxx_yyyy, g_yy_0_xxx_yyyz, g_yy_0_xxx_yyzz, g_yy_0_xxx_yzzz, g_yy_0_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxx_xxxx[k] = -g_yy_0_xx_xxxx[k] * ab_x + g_yy_0_xx_xxxxx[k];

                g_yy_0_xxx_xxxy[k] = -g_yy_0_xx_xxxy[k] * ab_x + g_yy_0_xx_xxxxy[k];

                g_yy_0_xxx_xxxz[k] = -g_yy_0_xx_xxxz[k] * ab_x + g_yy_0_xx_xxxxz[k];

                g_yy_0_xxx_xxyy[k] = -g_yy_0_xx_xxyy[k] * ab_x + g_yy_0_xx_xxxyy[k];

                g_yy_0_xxx_xxyz[k] = -g_yy_0_xx_xxyz[k] * ab_x + g_yy_0_xx_xxxyz[k];

                g_yy_0_xxx_xxzz[k] = -g_yy_0_xx_xxzz[k] * ab_x + g_yy_0_xx_xxxzz[k];

                g_yy_0_xxx_xyyy[k] = -g_yy_0_xx_xyyy[k] * ab_x + g_yy_0_xx_xxyyy[k];

                g_yy_0_xxx_xyyz[k] = -g_yy_0_xx_xyyz[k] * ab_x + g_yy_0_xx_xxyyz[k];

                g_yy_0_xxx_xyzz[k] = -g_yy_0_xx_xyzz[k] * ab_x + g_yy_0_xx_xxyzz[k];

                g_yy_0_xxx_xzzz[k] = -g_yy_0_xx_xzzz[k] * ab_x + g_yy_0_xx_xxzzz[k];

                g_yy_0_xxx_yyyy[k] = -g_yy_0_xx_yyyy[k] * ab_x + g_yy_0_xx_xyyyy[k];

                g_yy_0_xxx_yyyz[k] = -g_yy_0_xx_yyyz[k] * ab_x + g_yy_0_xx_xyyyz[k];

                g_yy_0_xxx_yyzz[k] = -g_yy_0_xx_yyzz[k] * ab_x + g_yy_0_xx_xyyzz[k];

                g_yy_0_xxx_yzzz[k] = -g_yy_0_xx_yzzz[k] * ab_x + g_yy_0_xx_xyzzz[k];

                g_yy_0_xxx_zzzz[k] = -g_yy_0_xx_zzzz[k] * ab_x + g_yy_0_xx_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxy_xxxx, g_yy_0_xxy_xxxy, g_yy_0_xxy_xxxz, g_yy_0_xxy_xxyy, g_yy_0_xxy_xxyz, g_yy_0_xxy_xxzz, g_yy_0_xxy_xyyy, g_yy_0_xxy_xyyz, g_yy_0_xxy_xyzz, g_yy_0_xxy_xzzz, g_yy_0_xxy_yyyy, g_yy_0_xxy_yyyz, g_yy_0_xxy_yyzz, g_yy_0_xxy_yzzz, g_yy_0_xxy_zzzz, g_yy_0_xy_xxxx, g_yy_0_xy_xxxxx, g_yy_0_xy_xxxxy, g_yy_0_xy_xxxxz, g_yy_0_xy_xxxy, g_yy_0_xy_xxxyy, g_yy_0_xy_xxxyz, g_yy_0_xy_xxxz, g_yy_0_xy_xxxzz, g_yy_0_xy_xxyy, g_yy_0_xy_xxyyy, g_yy_0_xy_xxyyz, g_yy_0_xy_xxyz, g_yy_0_xy_xxyzz, g_yy_0_xy_xxzz, g_yy_0_xy_xxzzz, g_yy_0_xy_xyyy, g_yy_0_xy_xyyyy, g_yy_0_xy_xyyyz, g_yy_0_xy_xyyz, g_yy_0_xy_xyyzz, g_yy_0_xy_xyzz, g_yy_0_xy_xyzzz, g_yy_0_xy_xzzz, g_yy_0_xy_xzzzz, g_yy_0_xy_yyyy, g_yy_0_xy_yyyz, g_yy_0_xy_yyzz, g_yy_0_xy_yzzz, g_yy_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxy_xxxx[k] = -g_yy_0_xy_xxxx[k] * ab_x + g_yy_0_xy_xxxxx[k];

                g_yy_0_xxy_xxxy[k] = -g_yy_0_xy_xxxy[k] * ab_x + g_yy_0_xy_xxxxy[k];

                g_yy_0_xxy_xxxz[k] = -g_yy_0_xy_xxxz[k] * ab_x + g_yy_0_xy_xxxxz[k];

                g_yy_0_xxy_xxyy[k] = -g_yy_0_xy_xxyy[k] * ab_x + g_yy_0_xy_xxxyy[k];

                g_yy_0_xxy_xxyz[k] = -g_yy_0_xy_xxyz[k] * ab_x + g_yy_0_xy_xxxyz[k];

                g_yy_0_xxy_xxzz[k] = -g_yy_0_xy_xxzz[k] * ab_x + g_yy_0_xy_xxxzz[k];

                g_yy_0_xxy_xyyy[k] = -g_yy_0_xy_xyyy[k] * ab_x + g_yy_0_xy_xxyyy[k];

                g_yy_0_xxy_xyyz[k] = -g_yy_0_xy_xyyz[k] * ab_x + g_yy_0_xy_xxyyz[k];

                g_yy_0_xxy_xyzz[k] = -g_yy_0_xy_xyzz[k] * ab_x + g_yy_0_xy_xxyzz[k];

                g_yy_0_xxy_xzzz[k] = -g_yy_0_xy_xzzz[k] * ab_x + g_yy_0_xy_xxzzz[k];

                g_yy_0_xxy_yyyy[k] = -g_yy_0_xy_yyyy[k] * ab_x + g_yy_0_xy_xyyyy[k];

                g_yy_0_xxy_yyyz[k] = -g_yy_0_xy_yyyz[k] * ab_x + g_yy_0_xy_xyyyz[k];

                g_yy_0_xxy_yyzz[k] = -g_yy_0_xy_yyzz[k] * ab_x + g_yy_0_xy_xyyzz[k];

                g_yy_0_xxy_yzzz[k] = -g_yy_0_xy_yzzz[k] * ab_x + g_yy_0_xy_xyzzz[k];

                g_yy_0_xxy_zzzz[k] = -g_yy_0_xy_zzzz[k] * ab_x + g_yy_0_xy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxz_xxxx, g_yy_0_xxz_xxxy, g_yy_0_xxz_xxxz, g_yy_0_xxz_xxyy, g_yy_0_xxz_xxyz, g_yy_0_xxz_xxzz, g_yy_0_xxz_xyyy, g_yy_0_xxz_xyyz, g_yy_0_xxz_xyzz, g_yy_0_xxz_xzzz, g_yy_0_xxz_yyyy, g_yy_0_xxz_yyyz, g_yy_0_xxz_yyzz, g_yy_0_xxz_yzzz, g_yy_0_xxz_zzzz, g_yy_0_xz_xxxx, g_yy_0_xz_xxxxx, g_yy_0_xz_xxxxy, g_yy_0_xz_xxxxz, g_yy_0_xz_xxxy, g_yy_0_xz_xxxyy, g_yy_0_xz_xxxyz, g_yy_0_xz_xxxz, g_yy_0_xz_xxxzz, g_yy_0_xz_xxyy, g_yy_0_xz_xxyyy, g_yy_0_xz_xxyyz, g_yy_0_xz_xxyz, g_yy_0_xz_xxyzz, g_yy_0_xz_xxzz, g_yy_0_xz_xxzzz, g_yy_0_xz_xyyy, g_yy_0_xz_xyyyy, g_yy_0_xz_xyyyz, g_yy_0_xz_xyyz, g_yy_0_xz_xyyzz, g_yy_0_xz_xyzz, g_yy_0_xz_xyzzz, g_yy_0_xz_xzzz, g_yy_0_xz_xzzzz, g_yy_0_xz_yyyy, g_yy_0_xz_yyyz, g_yy_0_xz_yyzz, g_yy_0_xz_yzzz, g_yy_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxz_xxxx[k] = -g_yy_0_xz_xxxx[k] * ab_x + g_yy_0_xz_xxxxx[k];

                g_yy_0_xxz_xxxy[k] = -g_yy_0_xz_xxxy[k] * ab_x + g_yy_0_xz_xxxxy[k];

                g_yy_0_xxz_xxxz[k] = -g_yy_0_xz_xxxz[k] * ab_x + g_yy_0_xz_xxxxz[k];

                g_yy_0_xxz_xxyy[k] = -g_yy_0_xz_xxyy[k] * ab_x + g_yy_0_xz_xxxyy[k];

                g_yy_0_xxz_xxyz[k] = -g_yy_0_xz_xxyz[k] * ab_x + g_yy_0_xz_xxxyz[k];

                g_yy_0_xxz_xxzz[k] = -g_yy_0_xz_xxzz[k] * ab_x + g_yy_0_xz_xxxzz[k];

                g_yy_0_xxz_xyyy[k] = -g_yy_0_xz_xyyy[k] * ab_x + g_yy_0_xz_xxyyy[k];

                g_yy_0_xxz_xyyz[k] = -g_yy_0_xz_xyyz[k] * ab_x + g_yy_0_xz_xxyyz[k];

                g_yy_0_xxz_xyzz[k] = -g_yy_0_xz_xyzz[k] * ab_x + g_yy_0_xz_xxyzz[k];

                g_yy_0_xxz_xzzz[k] = -g_yy_0_xz_xzzz[k] * ab_x + g_yy_0_xz_xxzzz[k];

                g_yy_0_xxz_yyyy[k] = -g_yy_0_xz_yyyy[k] * ab_x + g_yy_0_xz_xyyyy[k];

                g_yy_0_xxz_yyyz[k] = -g_yy_0_xz_yyyz[k] * ab_x + g_yy_0_xz_xyyyz[k];

                g_yy_0_xxz_yyzz[k] = -g_yy_0_xz_yyzz[k] * ab_x + g_yy_0_xz_xyyzz[k];

                g_yy_0_xxz_yzzz[k] = -g_yy_0_xz_yzzz[k] * ab_x + g_yy_0_xz_xyzzz[k];

                g_yy_0_xxz_zzzz[k] = -g_yy_0_xz_zzzz[k] * ab_x + g_yy_0_xz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 503 * ccomps * dcomps);

            auto g_yy_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyy_xxxx, g_yy_0_xyy_xxxy, g_yy_0_xyy_xxxz, g_yy_0_xyy_xxyy, g_yy_0_xyy_xxyz, g_yy_0_xyy_xxzz, g_yy_0_xyy_xyyy, g_yy_0_xyy_xyyz, g_yy_0_xyy_xyzz, g_yy_0_xyy_xzzz, g_yy_0_xyy_yyyy, g_yy_0_xyy_yyyz, g_yy_0_xyy_yyzz, g_yy_0_xyy_yzzz, g_yy_0_xyy_zzzz, g_yy_0_yy_xxxx, g_yy_0_yy_xxxxx, g_yy_0_yy_xxxxy, g_yy_0_yy_xxxxz, g_yy_0_yy_xxxy, g_yy_0_yy_xxxyy, g_yy_0_yy_xxxyz, g_yy_0_yy_xxxz, g_yy_0_yy_xxxzz, g_yy_0_yy_xxyy, g_yy_0_yy_xxyyy, g_yy_0_yy_xxyyz, g_yy_0_yy_xxyz, g_yy_0_yy_xxyzz, g_yy_0_yy_xxzz, g_yy_0_yy_xxzzz, g_yy_0_yy_xyyy, g_yy_0_yy_xyyyy, g_yy_0_yy_xyyyz, g_yy_0_yy_xyyz, g_yy_0_yy_xyyzz, g_yy_0_yy_xyzz, g_yy_0_yy_xyzzz, g_yy_0_yy_xzzz, g_yy_0_yy_xzzzz, g_yy_0_yy_yyyy, g_yy_0_yy_yyyz, g_yy_0_yy_yyzz, g_yy_0_yy_yzzz, g_yy_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyy_xxxx[k] = -g_yy_0_yy_xxxx[k] * ab_x + g_yy_0_yy_xxxxx[k];

                g_yy_0_xyy_xxxy[k] = -g_yy_0_yy_xxxy[k] * ab_x + g_yy_0_yy_xxxxy[k];

                g_yy_0_xyy_xxxz[k] = -g_yy_0_yy_xxxz[k] * ab_x + g_yy_0_yy_xxxxz[k];

                g_yy_0_xyy_xxyy[k] = -g_yy_0_yy_xxyy[k] * ab_x + g_yy_0_yy_xxxyy[k];

                g_yy_0_xyy_xxyz[k] = -g_yy_0_yy_xxyz[k] * ab_x + g_yy_0_yy_xxxyz[k];

                g_yy_0_xyy_xxzz[k] = -g_yy_0_yy_xxzz[k] * ab_x + g_yy_0_yy_xxxzz[k];

                g_yy_0_xyy_xyyy[k] = -g_yy_0_yy_xyyy[k] * ab_x + g_yy_0_yy_xxyyy[k];

                g_yy_0_xyy_xyyz[k] = -g_yy_0_yy_xyyz[k] * ab_x + g_yy_0_yy_xxyyz[k];

                g_yy_0_xyy_xyzz[k] = -g_yy_0_yy_xyzz[k] * ab_x + g_yy_0_yy_xxyzz[k];

                g_yy_0_xyy_xzzz[k] = -g_yy_0_yy_xzzz[k] * ab_x + g_yy_0_yy_xxzzz[k];

                g_yy_0_xyy_yyyy[k] = -g_yy_0_yy_yyyy[k] * ab_x + g_yy_0_yy_xyyyy[k];

                g_yy_0_xyy_yyyz[k] = -g_yy_0_yy_yyyz[k] * ab_x + g_yy_0_yy_xyyyz[k];

                g_yy_0_xyy_yyzz[k] = -g_yy_0_yy_yyzz[k] * ab_x + g_yy_0_yy_xyyzz[k];

                g_yy_0_xyy_yzzz[k] = -g_yy_0_yy_yzzz[k] * ab_x + g_yy_0_yy_xyzzz[k];

                g_yy_0_xyy_zzzz[k] = -g_yy_0_yy_zzzz[k] * ab_x + g_yy_0_yy_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 515 * ccomps * dcomps);

            auto g_yy_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 519 * ccomps * dcomps);

            auto g_yy_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 521 * ccomps * dcomps);

            auto g_yy_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyz_xxxx, g_yy_0_xyz_xxxy, g_yy_0_xyz_xxxz, g_yy_0_xyz_xxyy, g_yy_0_xyz_xxyz, g_yy_0_xyz_xxzz, g_yy_0_xyz_xyyy, g_yy_0_xyz_xyyz, g_yy_0_xyz_xyzz, g_yy_0_xyz_xzzz, g_yy_0_xyz_yyyy, g_yy_0_xyz_yyyz, g_yy_0_xyz_yyzz, g_yy_0_xyz_yzzz, g_yy_0_xyz_zzzz, g_yy_0_yz_xxxx, g_yy_0_yz_xxxxx, g_yy_0_yz_xxxxy, g_yy_0_yz_xxxxz, g_yy_0_yz_xxxy, g_yy_0_yz_xxxyy, g_yy_0_yz_xxxyz, g_yy_0_yz_xxxz, g_yy_0_yz_xxxzz, g_yy_0_yz_xxyy, g_yy_0_yz_xxyyy, g_yy_0_yz_xxyyz, g_yy_0_yz_xxyz, g_yy_0_yz_xxyzz, g_yy_0_yz_xxzz, g_yy_0_yz_xxzzz, g_yy_0_yz_xyyy, g_yy_0_yz_xyyyy, g_yy_0_yz_xyyyz, g_yy_0_yz_xyyz, g_yy_0_yz_xyyzz, g_yy_0_yz_xyzz, g_yy_0_yz_xyzzz, g_yy_0_yz_xzzz, g_yy_0_yz_xzzzz, g_yy_0_yz_yyyy, g_yy_0_yz_yyyz, g_yy_0_yz_yyzz, g_yy_0_yz_yzzz, g_yy_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyz_xxxx[k] = -g_yy_0_yz_xxxx[k] * ab_x + g_yy_0_yz_xxxxx[k];

                g_yy_0_xyz_xxxy[k] = -g_yy_0_yz_xxxy[k] * ab_x + g_yy_0_yz_xxxxy[k];

                g_yy_0_xyz_xxxz[k] = -g_yy_0_yz_xxxz[k] * ab_x + g_yy_0_yz_xxxxz[k];

                g_yy_0_xyz_xxyy[k] = -g_yy_0_yz_xxyy[k] * ab_x + g_yy_0_yz_xxxyy[k];

                g_yy_0_xyz_xxyz[k] = -g_yy_0_yz_xxyz[k] * ab_x + g_yy_0_yz_xxxyz[k];

                g_yy_0_xyz_xxzz[k] = -g_yy_0_yz_xxzz[k] * ab_x + g_yy_0_yz_xxxzz[k];

                g_yy_0_xyz_xyyy[k] = -g_yy_0_yz_xyyy[k] * ab_x + g_yy_0_yz_xxyyy[k];

                g_yy_0_xyz_xyyz[k] = -g_yy_0_yz_xyyz[k] * ab_x + g_yy_0_yz_xxyyz[k];

                g_yy_0_xyz_xyzz[k] = -g_yy_0_yz_xyzz[k] * ab_x + g_yy_0_yz_xxyzz[k];

                g_yy_0_xyz_xzzz[k] = -g_yy_0_yz_xzzz[k] * ab_x + g_yy_0_yz_xxzzz[k];

                g_yy_0_xyz_yyyy[k] = -g_yy_0_yz_yyyy[k] * ab_x + g_yy_0_yz_xyyyy[k];

                g_yy_0_xyz_yyyz[k] = -g_yy_0_yz_yyyz[k] * ab_x + g_yy_0_yz_xyyyz[k];

                g_yy_0_xyz_yyzz[k] = -g_yy_0_yz_yyzz[k] * ab_x + g_yy_0_yz_xyyzz[k];

                g_yy_0_xyz_yzzz[k] = -g_yy_0_yz_yzzz[k] * ab_x + g_yy_0_yz_xyzzz[k];

                g_yy_0_xyz_zzzz[k] = -g_yy_0_yz_zzzz[k] * ab_x + g_yy_0_yz_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 527 * ccomps * dcomps);

            auto g_yy_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 529 * ccomps * dcomps);

            auto g_yy_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 531 * ccomps * dcomps);

            auto g_yy_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 533 * ccomps * dcomps);

            auto g_yy_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzz_xxxx, g_yy_0_xzz_xxxy, g_yy_0_xzz_xxxz, g_yy_0_xzz_xxyy, g_yy_0_xzz_xxyz, g_yy_0_xzz_xxzz, g_yy_0_xzz_xyyy, g_yy_0_xzz_xyyz, g_yy_0_xzz_xyzz, g_yy_0_xzz_xzzz, g_yy_0_xzz_yyyy, g_yy_0_xzz_yyyz, g_yy_0_xzz_yyzz, g_yy_0_xzz_yzzz, g_yy_0_xzz_zzzz, g_yy_0_zz_xxxx, g_yy_0_zz_xxxxx, g_yy_0_zz_xxxxy, g_yy_0_zz_xxxxz, g_yy_0_zz_xxxy, g_yy_0_zz_xxxyy, g_yy_0_zz_xxxyz, g_yy_0_zz_xxxz, g_yy_0_zz_xxxzz, g_yy_0_zz_xxyy, g_yy_0_zz_xxyyy, g_yy_0_zz_xxyyz, g_yy_0_zz_xxyz, g_yy_0_zz_xxyzz, g_yy_0_zz_xxzz, g_yy_0_zz_xxzzz, g_yy_0_zz_xyyy, g_yy_0_zz_xyyyy, g_yy_0_zz_xyyyz, g_yy_0_zz_xyyz, g_yy_0_zz_xyyzz, g_yy_0_zz_xyzz, g_yy_0_zz_xyzzz, g_yy_0_zz_xzzz, g_yy_0_zz_xzzzz, g_yy_0_zz_yyyy, g_yy_0_zz_yyyz, g_yy_0_zz_yyzz, g_yy_0_zz_yzzz, g_yy_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzz_xxxx[k] = -g_yy_0_zz_xxxx[k] * ab_x + g_yy_0_zz_xxxxx[k];

                g_yy_0_xzz_xxxy[k] = -g_yy_0_zz_xxxy[k] * ab_x + g_yy_0_zz_xxxxy[k];

                g_yy_0_xzz_xxxz[k] = -g_yy_0_zz_xxxz[k] * ab_x + g_yy_0_zz_xxxxz[k];

                g_yy_0_xzz_xxyy[k] = -g_yy_0_zz_xxyy[k] * ab_x + g_yy_0_zz_xxxyy[k];

                g_yy_0_xzz_xxyz[k] = -g_yy_0_zz_xxyz[k] * ab_x + g_yy_0_zz_xxxyz[k];

                g_yy_0_xzz_xxzz[k] = -g_yy_0_zz_xxzz[k] * ab_x + g_yy_0_zz_xxxzz[k];

                g_yy_0_xzz_xyyy[k] = -g_yy_0_zz_xyyy[k] * ab_x + g_yy_0_zz_xxyyy[k];

                g_yy_0_xzz_xyyz[k] = -g_yy_0_zz_xyyz[k] * ab_x + g_yy_0_zz_xxyyz[k];

                g_yy_0_xzz_xyzz[k] = -g_yy_0_zz_xyzz[k] * ab_x + g_yy_0_zz_xxyzz[k];

                g_yy_0_xzz_xzzz[k] = -g_yy_0_zz_xzzz[k] * ab_x + g_yy_0_zz_xxzzz[k];

                g_yy_0_xzz_yyyy[k] = -g_yy_0_zz_yyyy[k] * ab_x + g_yy_0_zz_xyyyy[k];

                g_yy_0_xzz_yyyz[k] = -g_yy_0_zz_yyyz[k] * ab_x + g_yy_0_zz_xyyyz[k];

                g_yy_0_xzz_yyzz[k] = -g_yy_0_zz_yyzz[k] * ab_x + g_yy_0_zz_xyyzz[k];

                g_yy_0_xzz_yzzz[k] = -g_yy_0_zz_yzzz[k] * ab_x + g_yy_0_zz_xyzzz[k];

                g_yy_0_xzz_zzzz[k] = -g_yy_0_zz_zzzz[k] * ab_x + g_yy_0_zz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 545 * ccomps * dcomps);

            auto g_yy_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 549 * ccomps * dcomps);

            auto g_yy_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 551 * ccomps * dcomps);

            auto g_yy_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz, g_yy_0_yy_xxxx, g_yy_0_yy_xxxxy, g_yy_0_yy_xxxy, g_yy_0_yy_xxxyy, g_yy_0_yy_xxxyz, g_yy_0_yy_xxxz, g_yy_0_yy_xxyy, g_yy_0_yy_xxyyy, g_yy_0_yy_xxyyz, g_yy_0_yy_xxyz, g_yy_0_yy_xxyzz, g_yy_0_yy_xxzz, g_yy_0_yy_xyyy, g_yy_0_yy_xyyyy, g_yy_0_yy_xyyyz, g_yy_0_yy_xyyz, g_yy_0_yy_xyyzz, g_yy_0_yy_xyzz, g_yy_0_yy_xyzzz, g_yy_0_yy_xzzz, g_yy_0_yy_yyyy, g_yy_0_yy_yyyyy, g_yy_0_yy_yyyyz, g_yy_0_yy_yyyz, g_yy_0_yy_yyyzz, g_yy_0_yy_yyzz, g_yy_0_yy_yyzzz, g_yy_0_yy_yzzz, g_yy_0_yy_yzzzz, g_yy_0_yy_zzzz, g_yy_0_yyy_xxxx, g_yy_0_yyy_xxxy, g_yy_0_yyy_xxxz, g_yy_0_yyy_xxyy, g_yy_0_yyy_xxyz, g_yy_0_yyy_xxzz, g_yy_0_yyy_xyyy, g_yy_0_yyy_xyyz, g_yy_0_yyy_xyzz, g_yy_0_yyy_xzzz, g_yy_0_yyy_yyyy, g_yy_0_yyy_yyyz, g_yy_0_yyy_yyzz, g_yy_0_yyy_yzzz, g_yy_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyy_xxxx[k] = -2.0 * g_y_0_yy_xxxx[k] - g_yy_0_yy_xxxx[k] * ab_y + g_yy_0_yy_xxxxy[k];

                g_yy_0_yyy_xxxy[k] = -2.0 * g_y_0_yy_xxxy[k] - g_yy_0_yy_xxxy[k] * ab_y + g_yy_0_yy_xxxyy[k];

                g_yy_0_yyy_xxxz[k] = -2.0 * g_y_0_yy_xxxz[k] - g_yy_0_yy_xxxz[k] * ab_y + g_yy_0_yy_xxxyz[k];

                g_yy_0_yyy_xxyy[k] = -2.0 * g_y_0_yy_xxyy[k] - g_yy_0_yy_xxyy[k] * ab_y + g_yy_0_yy_xxyyy[k];

                g_yy_0_yyy_xxyz[k] = -2.0 * g_y_0_yy_xxyz[k] - g_yy_0_yy_xxyz[k] * ab_y + g_yy_0_yy_xxyyz[k];

                g_yy_0_yyy_xxzz[k] = -2.0 * g_y_0_yy_xxzz[k] - g_yy_0_yy_xxzz[k] * ab_y + g_yy_0_yy_xxyzz[k];

                g_yy_0_yyy_xyyy[k] = -2.0 * g_y_0_yy_xyyy[k] - g_yy_0_yy_xyyy[k] * ab_y + g_yy_0_yy_xyyyy[k];

                g_yy_0_yyy_xyyz[k] = -2.0 * g_y_0_yy_xyyz[k] - g_yy_0_yy_xyyz[k] * ab_y + g_yy_0_yy_xyyyz[k];

                g_yy_0_yyy_xyzz[k] = -2.0 * g_y_0_yy_xyzz[k] - g_yy_0_yy_xyzz[k] * ab_y + g_yy_0_yy_xyyzz[k];

                g_yy_0_yyy_xzzz[k] = -2.0 * g_y_0_yy_xzzz[k] - g_yy_0_yy_xzzz[k] * ab_y + g_yy_0_yy_xyzzz[k];

                g_yy_0_yyy_yyyy[k] = -2.0 * g_y_0_yy_yyyy[k] - g_yy_0_yy_yyyy[k] * ab_y + g_yy_0_yy_yyyyy[k];

                g_yy_0_yyy_yyyz[k] = -2.0 * g_y_0_yy_yyyz[k] - g_yy_0_yy_yyyz[k] * ab_y + g_yy_0_yy_yyyyz[k];

                g_yy_0_yyy_yyzz[k] = -2.0 * g_y_0_yy_yyzz[k] - g_yy_0_yy_yyzz[k] * ab_y + g_yy_0_yy_yyyzz[k];

                g_yy_0_yyy_yzzz[k] = -2.0 * g_y_0_yy_yzzz[k] - g_yy_0_yy_yzzz[k] * ab_y + g_yy_0_yy_yyzzz[k];

                g_yy_0_yyy_zzzz[k] = -2.0 * g_y_0_yy_zzzz[k] - g_yy_0_yy_zzzz[k] * ab_y + g_yy_0_yy_yzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 557 * ccomps * dcomps);

            auto g_yy_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 559 * ccomps * dcomps);

            auto g_yy_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 563 * ccomps * dcomps);

            auto g_yy_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yy_xxxx, g_yy_0_yy_xxxxz, g_yy_0_yy_xxxy, g_yy_0_yy_xxxyz, g_yy_0_yy_xxxz, g_yy_0_yy_xxxzz, g_yy_0_yy_xxyy, g_yy_0_yy_xxyyz, g_yy_0_yy_xxyz, g_yy_0_yy_xxyzz, g_yy_0_yy_xxzz, g_yy_0_yy_xxzzz, g_yy_0_yy_xyyy, g_yy_0_yy_xyyyz, g_yy_0_yy_xyyz, g_yy_0_yy_xyyzz, g_yy_0_yy_xyzz, g_yy_0_yy_xyzzz, g_yy_0_yy_xzzz, g_yy_0_yy_xzzzz, g_yy_0_yy_yyyy, g_yy_0_yy_yyyyz, g_yy_0_yy_yyyz, g_yy_0_yy_yyyzz, g_yy_0_yy_yyzz, g_yy_0_yy_yyzzz, g_yy_0_yy_yzzz, g_yy_0_yy_yzzzz, g_yy_0_yy_zzzz, g_yy_0_yy_zzzzz, g_yy_0_yyz_xxxx, g_yy_0_yyz_xxxy, g_yy_0_yyz_xxxz, g_yy_0_yyz_xxyy, g_yy_0_yyz_xxyz, g_yy_0_yyz_xxzz, g_yy_0_yyz_xyyy, g_yy_0_yyz_xyyz, g_yy_0_yyz_xyzz, g_yy_0_yyz_xzzz, g_yy_0_yyz_yyyy, g_yy_0_yyz_yyyz, g_yy_0_yyz_yyzz, g_yy_0_yyz_yzzz, g_yy_0_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyz_xxxx[k] = -g_yy_0_yy_xxxx[k] * ab_z + g_yy_0_yy_xxxxz[k];

                g_yy_0_yyz_xxxy[k] = -g_yy_0_yy_xxxy[k] * ab_z + g_yy_0_yy_xxxyz[k];

                g_yy_0_yyz_xxxz[k] = -g_yy_0_yy_xxxz[k] * ab_z + g_yy_0_yy_xxxzz[k];

                g_yy_0_yyz_xxyy[k] = -g_yy_0_yy_xxyy[k] * ab_z + g_yy_0_yy_xxyyz[k];

                g_yy_0_yyz_xxyz[k] = -g_yy_0_yy_xxyz[k] * ab_z + g_yy_0_yy_xxyzz[k];

                g_yy_0_yyz_xxzz[k] = -g_yy_0_yy_xxzz[k] * ab_z + g_yy_0_yy_xxzzz[k];

                g_yy_0_yyz_xyyy[k] = -g_yy_0_yy_xyyy[k] * ab_z + g_yy_0_yy_xyyyz[k];

                g_yy_0_yyz_xyyz[k] = -g_yy_0_yy_xyyz[k] * ab_z + g_yy_0_yy_xyyzz[k];

                g_yy_0_yyz_xyzz[k] = -g_yy_0_yy_xyzz[k] * ab_z + g_yy_0_yy_xyzzz[k];

                g_yy_0_yyz_xzzz[k] = -g_yy_0_yy_xzzz[k] * ab_z + g_yy_0_yy_xzzzz[k];

                g_yy_0_yyz_yyyy[k] = -g_yy_0_yy_yyyy[k] * ab_z + g_yy_0_yy_yyyyz[k];

                g_yy_0_yyz_yyyz[k] = -g_yy_0_yy_yyyz[k] * ab_z + g_yy_0_yy_yyyzz[k];

                g_yy_0_yyz_yyzz[k] = -g_yy_0_yy_yyzz[k] * ab_z + g_yy_0_yy_yyzzz[k];

                g_yy_0_yyz_yzzz[k] = -g_yy_0_yy_yzzz[k] * ab_z + g_yy_0_yy_yzzzz[k];

                g_yy_0_yyz_zzzz[k] = -g_yy_0_yy_zzzz[k] * ab_z + g_yy_0_yy_zzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 575 * ccomps * dcomps);

            auto g_yy_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 579 * ccomps * dcomps);

            auto g_yy_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 581 * ccomps * dcomps);

            auto g_yy_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yz_xxxx, g_yy_0_yz_xxxxz, g_yy_0_yz_xxxy, g_yy_0_yz_xxxyz, g_yy_0_yz_xxxz, g_yy_0_yz_xxxzz, g_yy_0_yz_xxyy, g_yy_0_yz_xxyyz, g_yy_0_yz_xxyz, g_yy_0_yz_xxyzz, g_yy_0_yz_xxzz, g_yy_0_yz_xxzzz, g_yy_0_yz_xyyy, g_yy_0_yz_xyyyz, g_yy_0_yz_xyyz, g_yy_0_yz_xyyzz, g_yy_0_yz_xyzz, g_yy_0_yz_xyzzz, g_yy_0_yz_xzzz, g_yy_0_yz_xzzzz, g_yy_0_yz_yyyy, g_yy_0_yz_yyyyz, g_yy_0_yz_yyyz, g_yy_0_yz_yyyzz, g_yy_0_yz_yyzz, g_yy_0_yz_yyzzz, g_yy_0_yz_yzzz, g_yy_0_yz_yzzzz, g_yy_0_yz_zzzz, g_yy_0_yz_zzzzz, g_yy_0_yzz_xxxx, g_yy_0_yzz_xxxy, g_yy_0_yzz_xxxz, g_yy_0_yzz_xxyy, g_yy_0_yzz_xxyz, g_yy_0_yzz_xxzz, g_yy_0_yzz_xyyy, g_yy_0_yzz_xyyz, g_yy_0_yzz_xyzz, g_yy_0_yzz_xzzz, g_yy_0_yzz_yyyy, g_yy_0_yzz_yyyz, g_yy_0_yzz_yyzz, g_yy_0_yzz_yzzz, g_yy_0_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzz_xxxx[k] = -g_yy_0_yz_xxxx[k] * ab_z + g_yy_0_yz_xxxxz[k];

                g_yy_0_yzz_xxxy[k] = -g_yy_0_yz_xxxy[k] * ab_z + g_yy_0_yz_xxxyz[k];

                g_yy_0_yzz_xxxz[k] = -g_yy_0_yz_xxxz[k] * ab_z + g_yy_0_yz_xxxzz[k];

                g_yy_0_yzz_xxyy[k] = -g_yy_0_yz_xxyy[k] * ab_z + g_yy_0_yz_xxyyz[k];

                g_yy_0_yzz_xxyz[k] = -g_yy_0_yz_xxyz[k] * ab_z + g_yy_0_yz_xxyzz[k];

                g_yy_0_yzz_xxzz[k] = -g_yy_0_yz_xxzz[k] * ab_z + g_yy_0_yz_xxzzz[k];

                g_yy_0_yzz_xyyy[k] = -g_yy_0_yz_xyyy[k] * ab_z + g_yy_0_yz_xyyyz[k];

                g_yy_0_yzz_xyyz[k] = -g_yy_0_yz_xyyz[k] * ab_z + g_yy_0_yz_xyyzz[k];

                g_yy_0_yzz_xyzz[k] = -g_yy_0_yz_xyzz[k] * ab_z + g_yy_0_yz_xyzzz[k];

                g_yy_0_yzz_xzzz[k] = -g_yy_0_yz_xzzz[k] * ab_z + g_yy_0_yz_xzzzz[k];

                g_yy_0_yzz_yyyy[k] = -g_yy_0_yz_yyyy[k] * ab_z + g_yy_0_yz_yyyyz[k];

                g_yy_0_yzz_yyyz[k] = -g_yy_0_yz_yyyz[k] * ab_z + g_yy_0_yz_yyyzz[k];

                g_yy_0_yzz_yyzz[k] = -g_yy_0_yz_yyzz[k] * ab_z + g_yy_0_yz_yyzzz[k];

                g_yy_0_yzz_yzzz[k] = -g_yy_0_yz_yzzz[k] * ab_z + g_yy_0_yz_yzzzz[k];

                g_yy_0_yzz_zzzz[k] = -g_yy_0_yz_zzzz[k] * ab_z + g_yy_0_yz_zzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 587 * ccomps * dcomps);

            auto g_yy_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 589 * ccomps * dcomps);

            auto g_yy_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 593 * ccomps * dcomps);

            auto g_yy_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zz_xxxx, g_yy_0_zz_xxxxz, g_yy_0_zz_xxxy, g_yy_0_zz_xxxyz, g_yy_0_zz_xxxz, g_yy_0_zz_xxxzz, g_yy_0_zz_xxyy, g_yy_0_zz_xxyyz, g_yy_0_zz_xxyz, g_yy_0_zz_xxyzz, g_yy_0_zz_xxzz, g_yy_0_zz_xxzzz, g_yy_0_zz_xyyy, g_yy_0_zz_xyyyz, g_yy_0_zz_xyyz, g_yy_0_zz_xyyzz, g_yy_0_zz_xyzz, g_yy_0_zz_xyzzz, g_yy_0_zz_xzzz, g_yy_0_zz_xzzzz, g_yy_0_zz_yyyy, g_yy_0_zz_yyyyz, g_yy_0_zz_yyyz, g_yy_0_zz_yyyzz, g_yy_0_zz_yyzz, g_yy_0_zz_yyzzz, g_yy_0_zz_yzzz, g_yy_0_zz_yzzzz, g_yy_0_zz_zzzz, g_yy_0_zz_zzzzz, g_yy_0_zzz_xxxx, g_yy_0_zzz_xxxy, g_yy_0_zzz_xxxz, g_yy_0_zzz_xxyy, g_yy_0_zzz_xxyz, g_yy_0_zzz_xxzz, g_yy_0_zzz_xyyy, g_yy_0_zzz_xyyz, g_yy_0_zzz_xyzz, g_yy_0_zzz_xzzz, g_yy_0_zzz_yyyy, g_yy_0_zzz_yyyz, g_yy_0_zzz_yyzz, g_yy_0_zzz_yzzz, g_yy_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzz_xxxx[k] = -g_yy_0_zz_xxxx[k] * ab_z + g_yy_0_zz_xxxxz[k];

                g_yy_0_zzz_xxxy[k] = -g_yy_0_zz_xxxy[k] * ab_z + g_yy_0_zz_xxxyz[k];

                g_yy_0_zzz_xxxz[k] = -g_yy_0_zz_xxxz[k] * ab_z + g_yy_0_zz_xxxzz[k];

                g_yy_0_zzz_xxyy[k] = -g_yy_0_zz_xxyy[k] * ab_z + g_yy_0_zz_xxyyz[k];

                g_yy_0_zzz_xxyz[k] = -g_yy_0_zz_xxyz[k] * ab_z + g_yy_0_zz_xxyzz[k];

                g_yy_0_zzz_xxzz[k] = -g_yy_0_zz_xxzz[k] * ab_z + g_yy_0_zz_xxzzz[k];

                g_yy_0_zzz_xyyy[k] = -g_yy_0_zz_xyyy[k] * ab_z + g_yy_0_zz_xyyyz[k];

                g_yy_0_zzz_xyyz[k] = -g_yy_0_zz_xyyz[k] * ab_z + g_yy_0_zz_xyyzz[k];

                g_yy_0_zzz_xyzz[k] = -g_yy_0_zz_xyzz[k] * ab_z + g_yy_0_zz_xyzzz[k];

                g_yy_0_zzz_xzzz[k] = -g_yy_0_zz_xzzz[k] * ab_z + g_yy_0_zz_xzzzz[k];

                g_yy_0_zzz_yyyy[k] = -g_yy_0_zz_yyyy[k] * ab_z + g_yy_0_zz_yyyyz[k];

                g_yy_0_zzz_yyyz[k] = -g_yy_0_zz_yyyz[k] * ab_z + g_yy_0_zz_yyyzz[k];

                g_yy_0_zzz_yyzz[k] = -g_yy_0_zz_yyzz[k] * ab_z + g_yy_0_zz_yyzzz[k];

                g_yy_0_zzz_yzzz[k] = -g_yy_0_zz_yzzz[k] * ab_z + g_yy_0_zz_yzzzz[k];

                g_yy_0_zzz_zzzz[k] = -g_yy_0_zz_zzzz[k] * ab_z + g_yy_0_zz_zzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xx_xxxx, g_yz_0_xx_xxxxx, g_yz_0_xx_xxxxy, g_yz_0_xx_xxxxz, g_yz_0_xx_xxxy, g_yz_0_xx_xxxyy, g_yz_0_xx_xxxyz, g_yz_0_xx_xxxz, g_yz_0_xx_xxxzz, g_yz_0_xx_xxyy, g_yz_0_xx_xxyyy, g_yz_0_xx_xxyyz, g_yz_0_xx_xxyz, g_yz_0_xx_xxyzz, g_yz_0_xx_xxzz, g_yz_0_xx_xxzzz, g_yz_0_xx_xyyy, g_yz_0_xx_xyyyy, g_yz_0_xx_xyyyz, g_yz_0_xx_xyyz, g_yz_0_xx_xyyzz, g_yz_0_xx_xyzz, g_yz_0_xx_xyzzz, g_yz_0_xx_xzzz, g_yz_0_xx_xzzzz, g_yz_0_xx_yyyy, g_yz_0_xx_yyyz, g_yz_0_xx_yyzz, g_yz_0_xx_yzzz, g_yz_0_xx_zzzz, g_yz_0_xxx_xxxx, g_yz_0_xxx_xxxy, g_yz_0_xxx_xxxz, g_yz_0_xxx_xxyy, g_yz_0_xxx_xxyz, g_yz_0_xxx_xxzz, g_yz_0_xxx_xyyy, g_yz_0_xxx_xyyz, g_yz_0_xxx_xyzz, g_yz_0_xxx_xzzz, g_yz_0_xxx_yyyy, g_yz_0_xxx_yyyz, g_yz_0_xxx_yyzz, g_yz_0_xxx_yzzz, g_yz_0_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxx_xxxx[k] = -g_yz_0_xx_xxxx[k] * ab_x + g_yz_0_xx_xxxxx[k];

                g_yz_0_xxx_xxxy[k] = -g_yz_0_xx_xxxy[k] * ab_x + g_yz_0_xx_xxxxy[k];

                g_yz_0_xxx_xxxz[k] = -g_yz_0_xx_xxxz[k] * ab_x + g_yz_0_xx_xxxxz[k];

                g_yz_0_xxx_xxyy[k] = -g_yz_0_xx_xxyy[k] * ab_x + g_yz_0_xx_xxxyy[k];

                g_yz_0_xxx_xxyz[k] = -g_yz_0_xx_xxyz[k] * ab_x + g_yz_0_xx_xxxyz[k];

                g_yz_0_xxx_xxzz[k] = -g_yz_0_xx_xxzz[k] * ab_x + g_yz_0_xx_xxxzz[k];

                g_yz_0_xxx_xyyy[k] = -g_yz_0_xx_xyyy[k] * ab_x + g_yz_0_xx_xxyyy[k];

                g_yz_0_xxx_xyyz[k] = -g_yz_0_xx_xyyz[k] * ab_x + g_yz_0_xx_xxyyz[k];

                g_yz_0_xxx_xyzz[k] = -g_yz_0_xx_xyzz[k] * ab_x + g_yz_0_xx_xxyzz[k];

                g_yz_0_xxx_xzzz[k] = -g_yz_0_xx_xzzz[k] * ab_x + g_yz_0_xx_xxzzz[k];

                g_yz_0_xxx_yyyy[k] = -g_yz_0_xx_yyyy[k] * ab_x + g_yz_0_xx_xyyyy[k];

                g_yz_0_xxx_yyyz[k] = -g_yz_0_xx_yyyz[k] * ab_x + g_yz_0_xx_xyyyz[k];

                g_yz_0_xxx_yyzz[k] = -g_yz_0_xx_yyzz[k] * ab_x + g_yz_0_xx_xyyzz[k];

                g_yz_0_xxx_yzzz[k] = -g_yz_0_xx_yzzz[k] * ab_x + g_yz_0_xx_xyzzz[k];

                g_yz_0_xxx_zzzz[k] = -g_yz_0_xx_zzzz[k] * ab_x + g_yz_0_xx_xzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxy_xxxx, g_yz_0_xxy_xxxy, g_yz_0_xxy_xxxz, g_yz_0_xxy_xxyy, g_yz_0_xxy_xxyz, g_yz_0_xxy_xxzz, g_yz_0_xxy_xyyy, g_yz_0_xxy_xyyz, g_yz_0_xxy_xyzz, g_yz_0_xxy_xzzz, g_yz_0_xxy_yyyy, g_yz_0_xxy_yyyz, g_yz_0_xxy_yyzz, g_yz_0_xxy_yzzz, g_yz_0_xxy_zzzz, g_yz_0_xy_xxxx, g_yz_0_xy_xxxxx, g_yz_0_xy_xxxxy, g_yz_0_xy_xxxxz, g_yz_0_xy_xxxy, g_yz_0_xy_xxxyy, g_yz_0_xy_xxxyz, g_yz_0_xy_xxxz, g_yz_0_xy_xxxzz, g_yz_0_xy_xxyy, g_yz_0_xy_xxyyy, g_yz_0_xy_xxyyz, g_yz_0_xy_xxyz, g_yz_0_xy_xxyzz, g_yz_0_xy_xxzz, g_yz_0_xy_xxzzz, g_yz_0_xy_xyyy, g_yz_0_xy_xyyyy, g_yz_0_xy_xyyyz, g_yz_0_xy_xyyz, g_yz_0_xy_xyyzz, g_yz_0_xy_xyzz, g_yz_0_xy_xyzzz, g_yz_0_xy_xzzz, g_yz_0_xy_xzzzz, g_yz_0_xy_yyyy, g_yz_0_xy_yyyz, g_yz_0_xy_yyzz, g_yz_0_xy_yzzz, g_yz_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxy_xxxx[k] = -g_yz_0_xy_xxxx[k] * ab_x + g_yz_0_xy_xxxxx[k];

                g_yz_0_xxy_xxxy[k] = -g_yz_0_xy_xxxy[k] * ab_x + g_yz_0_xy_xxxxy[k];

                g_yz_0_xxy_xxxz[k] = -g_yz_0_xy_xxxz[k] * ab_x + g_yz_0_xy_xxxxz[k];

                g_yz_0_xxy_xxyy[k] = -g_yz_0_xy_xxyy[k] * ab_x + g_yz_0_xy_xxxyy[k];

                g_yz_0_xxy_xxyz[k] = -g_yz_0_xy_xxyz[k] * ab_x + g_yz_0_xy_xxxyz[k];

                g_yz_0_xxy_xxzz[k] = -g_yz_0_xy_xxzz[k] * ab_x + g_yz_0_xy_xxxzz[k];

                g_yz_0_xxy_xyyy[k] = -g_yz_0_xy_xyyy[k] * ab_x + g_yz_0_xy_xxyyy[k];

                g_yz_0_xxy_xyyz[k] = -g_yz_0_xy_xyyz[k] * ab_x + g_yz_0_xy_xxyyz[k];

                g_yz_0_xxy_xyzz[k] = -g_yz_0_xy_xyzz[k] * ab_x + g_yz_0_xy_xxyzz[k];

                g_yz_0_xxy_xzzz[k] = -g_yz_0_xy_xzzz[k] * ab_x + g_yz_0_xy_xxzzz[k];

                g_yz_0_xxy_yyyy[k] = -g_yz_0_xy_yyyy[k] * ab_x + g_yz_0_xy_xyyyy[k];

                g_yz_0_xxy_yyyz[k] = -g_yz_0_xy_yyyz[k] * ab_x + g_yz_0_xy_xyyyz[k];

                g_yz_0_xxy_yyzz[k] = -g_yz_0_xy_yyzz[k] * ab_x + g_yz_0_xy_xyyzz[k];

                g_yz_0_xxy_yzzz[k] = -g_yz_0_xy_yzzz[k] * ab_x + g_yz_0_xy_xyzzz[k];

                g_yz_0_xxy_zzzz[k] = -g_yz_0_xy_zzzz[k] * ab_x + g_yz_0_xy_xzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 630 * ccomps * dcomps);

            auto g_yz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 631 * ccomps * dcomps);

            auto g_yz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 632 * ccomps * dcomps);

            auto g_yz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 633 * ccomps * dcomps);

            auto g_yz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 634 * ccomps * dcomps);

            auto g_yz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 635 * ccomps * dcomps);

            auto g_yz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 636 * ccomps * dcomps);

            auto g_yz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 637 * ccomps * dcomps);

            auto g_yz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 638 * ccomps * dcomps);

            auto g_yz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 639 * ccomps * dcomps);

            auto g_yz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 640 * ccomps * dcomps);

            auto g_yz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 641 * ccomps * dcomps);

            auto g_yz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 642 * ccomps * dcomps);

            auto g_yz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 643 * ccomps * dcomps);

            auto g_yz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxz_xxxx, g_yz_0_xxz_xxxy, g_yz_0_xxz_xxxz, g_yz_0_xxz_xxyy, g_yz_0_xxz_xxyz, g_yz_0_xxz_xxzz, g_yz_0_xxz_xyyy, g_yz_0_xxz_xyyz, g_yz_0_xxz_xyzz, g_yz_0_xxz_xzzz, g_yz_0_xxz_yyyy, g_yz_0_xxz_yyyz, g_yz_0_xxz_yyzz, g_yz_0_xxz_yzzz, g_yz_0_xxz_zzzz, g_yz_0_xz_xxxx, g_yz_0_xz_xxxxx, g_yz_0_xz_xxxxy, g_yz_0_xz_xxxxz, g_yz_0_xz_xxxy, g_yz_0_xz_xxxyy, g_yz_0_xz_xxxyz, g_yz_0_xz_xxxz, g_yz_0_xz_xxxzz, g_yz_0_xz_xxyy, g_yz_0_xz_xxyyy, g_yz_0_xz_xxyyz, g_yz_0_xz_xxyz, g_yz_0_xz_xxyzz, g_yz_0_xz_xxzz, g_yz_0_xz_xxzzz, g_yz_0_xz_xyyy, g_yz_0_xz_xyyyy, g_yz_0_xz_xyyyz, g_yz_0_xz_xyyz, g_yz_0_xz_xyyzz, g_yz_0_xz_xyzz, g_yz_0_xz_xyzzz, g_yz_0_xz_xzzz, g_yz_0_xz_xzzzz, g_yz_0_xz_yyyy, g_yz_0_xz_yyyz, g_yz_0_xz_yyzz, g_yz_0_xz_yzzz, g_yz_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxz_xxxx[k] = -g_yz_0_xz_xxxx[k] * ab_x + g_yz_0_xz_xxxxx[k];

                g_yz_0_xxz_xxxy[k] = -g_yz_0_xz_xxxy[k] * ab_x + g_yz_0_xz_xxxxy[k];

                g_yz_0_xxz_xxxz[k] = -g_yz_0_xz_xxxz[k] * ab_x + g_yz_0_xz_xxxxz[k];

                g_yz_0_xxz_xxyy[k] = -g_yz_0_xz_xxyy[k] * ab_x + g_yz_0_xz_xxxyy[k];

                g_yz_0_xxz_xxyz[k] = -g_yz_0_xz_xxyz[k] * ab_x + g_yz_0_xz_xxxyz[k];

                g_yz_0_xxz_xxzz[k] = -g_yz_0_xz_xxzz[k] * ab_x + g_yz_0_xz_xxxzz[k];

                g_yz_0_xxz_xyyy[k] = -g_yz_0_xz_xyyy[k] * ab_x + g_yz_0_xz_xxyyy[k];

                g_yz_0_xxz_xyyz[k] = -g_yz_0_xz_xyyz[k] * ab_x + g_yz_0_xz_xxyyz[k];

                g_yz_0_xxz_xyzz[k] = -g_yz_0_xz_xyzz[k] * ab_x + g_yz_0_xz_xxyzz[k];

                g_yz_0_xxz_xzzz[k] = -g_yz_0_xz_xzzz[k] * ab_x + g_yz_0_xz_xxzzz[k];

                g_yz_0_xxz_yyyy[k] = -g_yz_0_xz_yyyy[k] * ab_x + g_yz_0_xz_xyyyy[k];

                g_yz_0_xxz_yyyz[k] = -g_yz_0_xz_yyyz[k] * ab_x + g_yz_0_xz_xyyyz[k];

                g_yz_0_xxz_yyzz[k] = -g_yz_0_xz_yyzz[k] * ab_x + g_yz_0_xz_xyyzz[k];

                g_yz_0_xxz_yzzz[k] = -g_yz_0_xz_yzzz[k] * ab_x + g_yz_0_xz_xyzzz[k];

                g_yz_0_xxz_zzzz[k] = -g_yz_0_xz_zzzz[k] * ab_x + g_yz_0_xz_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 645 * ccomps * dcomps);

            auto g_yz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 646 * ccomps * dcomps);

            auto g_yz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 647 * ccomps * dcomps);

            auto g_yz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 648 * ccomps * dcomps);

            auto g_yz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 649 * ccomps * dcomps);

            auto g_yz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 650 * ccomps * dcomps);

            auto g_yz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 651 * ccomps * dcomps);

            auto g_yz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 652 * ccomps * dcomps);

            auto g_yz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 653 * ccomps * dcomps);

            auto g_yz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 654 * ccomps * dcomps);

            auto g_yz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 655 * ccomps * dcomps);

            auto g_yz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 656 * ccomps * dcomps);

            auto g_yz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 657 * ccomps * dcomps);

            auto g_yz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 658 * ccomps * dcomps);

            auto g_yz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyy_xxxx, g_yz_0_xyy_xxxy, g_yz_0_xyy_xxxz, g_yz_0_xyy_xxyy, g_yz_0_xyy_xxyz, g_yz_0_xyy_xxzz, g_yz_0_xyy_xyyy, g_yz_0_xyy_xyyz, g_yz_0_xyy_xyzz, g_yz_0_xyy_xzzz, g_yz_0_xyy_yyyy, g_yz_0_xyy_yyyz, g_yz_0_xyy_yyzz, g_yz_0_xyy_yzzz, g_yz_0_xyy_zzzz, g_yz_0_yy_xxxx, g_yz_0_yy_xxxxx, g_yz_0_yy_xxxxy, g_yz_0_yy_xxxxz, g_yz_0_yy_xxxy, g_yz_0_yy_xxxyy, g_yz_0_yy_xxxyz, g_yz_0_yy_xxxz, g_yz_0_yy_xxxzz, g_yz_0_yy_xxyy, g_yz_0_yy_xxyyy, g_yz_0_yy_xxyyz, g_yz_0_yy_xxyz, g_yz_0_yy_xxyzz, g_yz_0_yy_xxzz, g_yz_0_yy_xxzzz, g_yz_0_yy_xyyy, g_yz_0_yy_xyyyy, g_yz_0_yy_xyyyz, g_yz_0_yy_xyyz, g_yz_0_yy_xyyzz, g_yz_0_yy_xyzz, g_yz_0_yy_xyzzz, g_yz_0_yy_xzzz, g_yz_0_yy_xzzzz, g_yz_0_yy_yyyy, g_yz_0_yy_yyyz, g_yz_0_yy_yyzz, g_yz_0_yy_yzzz, g_yz_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyy_xxxx[k] = -g_yz_0_yy_xxxx[k] * ab_x + g_yz_0_yy_xxxxx[k];

                g_yz_0_xyy_xxxy[k] = -g_yz_0_yy_xxxy[k] * ab_x + g_yz_0_yy_xxxxy[k];

                g_yz_0_xyy_xxxz[k] = -g_yz_0_yy_xxxz[k] * ab_x + g_yz_0_yy_xxxxz[k];

                g_yz_0_xyy_xxyy[k] = -g_yz_0_yy_xxyy[k] * ab_x + g_yz_0_yy_xxxyy[k];

                g_yz_0_xyy_xxyz[k] = -g_yz_0_yy_xxyz[k] * ab_x + g_yz_0_yy_xxxyz[k];

                g_yz_0_xyy_xxzz[k] = -g_yz_0_yy_xxzz[k] * ab_x + g_yz_0_yy_xxxzz[k];

                g_yz_0_xyy_xyyy[k] = -g_yz_0_yy_xyyy[k] * ab_x + g_yz_0_yy_xxyyy[k];

                g_yz_0_xyy_xyyz[k] = -g_yz_0_yy_xyyz[k] * ab_x + g_yz_0_yy_xxyyz[k];

                g_yz_0_xyy_xyzz[k] = -g_yz_0_yy_xyzz[k] * ab_x + g_yz_0_yy_xxyzz[k];

                g_yz_0_xyy_xzzz[k] = -g_yz_0_yy_xzzz[k] * ab_x + g_yz_0_yy_xxzzz[k];

                g_yz_0_xyy_yyyy[k] = -g_yz_0_yy_yyyy[k] * ab_x + g_yz_0_yy_xyyyy[k];

                g_yz_0_xyy_yyyz[k] = -g_yz_0_yy_yyyz[k] * ab_x + g_yz_0_yy_xyyyz[k];

                g_yz_0_xyy_yyzz[k] = -g_yz_0_yy_yyzz[k] * ab_x + g_yz_0_yy_xyyzz[k];

                g_yz_0_xyy_yzzz[k] = -g_yz_0_yy_yzzz[k] * ab_x + g_yz_0_yy_xyzzz[k];

                g_yz_0_xyy_zzzz[k] = -g_yz_0_yy_zzzz[k] * ab_x + g_yz_0_yy_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 660 * ccomps * dcomps);

            auto g_yz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 661 * ccomps * dcomps);

            auto g_yz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 662 * ccomps * dcomps);

            auto g_yz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 663 * ccomps * dcomps);

            auto g_yz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 664 * ccomps * dcomps);

            auto g_yz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 665 * ccomps * dcomps);

            auto g_yz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 666 * ccomps * dcomps);

            auto g_yz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 667 * ccomps * dcomps);

            auto g_yz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 668 * ccomps * dcomps);

            auto g_yz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 669 * ccomps * dcomps);

            auto g_yz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 670 * ccomps * dcomps);

            auto g_yz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 671 * ccomps * dcomps);

            auto g_yz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyz_xxxx, g_yz_0_xyz_xxxy, g_yz_0_xyz_xxxz, g_yz_0_xyz_xxyy, g_yz_0_xyz_xxyz, g_yz_0_xyz_xxzz, g_yz_0_xyz_xyyy, g_yz_0_xyz_xyyz, g_yz_0_xyz_xyzz, g_yz_0_xyz_xzzz, g_yz_0_xyz_yyyy, g_yz_0_xyz_yyyz, g_yz_0_xyz_yyzz, g_yz_0_xyz_yzzz, g_yz_0_xyz_zzzz, g_yz_0_yz_xxxx, g_yz_0_yz_xxxxx, g_yz_0_yz_xxxxy, g_yz_0_yz_xxxxz, g_yz_0_yz_xxxy, g_yz_0_yz_xxxyy, g_yz_0_yz_xxxyz, g_yz_0_yz_xxxz, g_yz_0_yz_xxxzz, g_yz_0_yz_xxyy, g_yz_0_yz_xxyyy, g_yz_0_yz_xxyyz, g_yz_0_yz_xxyz, g_yz_0_yz_xxyzz, g_yz_0_yz_xxzz, g_yz_0_yz_xxzzz, g_yz_0_yz_xyyy, g_yz_0_yz_xyyyy, g_yz_0_yz_xyyyz, g_yz_0_yz_xyyz, g_yz_0_yz_xyyzz, g_yz_0_yz_xyzz, g_yz_0_yz_xyzzz, g_yz_0_yz_xzzz, g_yz_0_yz_xzzzz, g_yz_0_yz_yyyy, g_yz_0_yz_yyyz, g_yz_0_yz_yyzz, g_yz_0_yz_yzzz, g_yz_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyz_xxxx[k] = -g_yz_0_yz_xxxx[k] * ab_x + g_yz_0_yz_xxxxx[k];

                g_yz_0_xyz_xxxy[k] = -g_yz_0_yz_xxxy[k] * ab_x + g_yz_0_yz_xxxxy[k];

                g_yz_0_xyz_xxxz[k] = -g_yz_0_yz_xxxz[k] * ab_x + g_yz_0_yz_xxxxz[k];

                g_yz_0_xyz_xxyy[k] = -g_yz_0_yz_xxyy[k] * ab_x + g_yz_0_yz_xxxyy[k];

                g_yz_0_xyz_xxyz[k] = -g_yz_0_yz_xxyz[k] * ab_x + g_yz_0_yz_xxxyz[k];

                g_yz_0_xyz_xxzz[k] = -g_yz_0_yz_xxzz[k] * ab_x + g_yz_0_yz_xxxzz[k];

                g_yz_0_xyz_xyyy[k] = -g_yz_0_yz_xyyy[k] * ab_x + g_yz_0_yz_xxyyy[k];

                g_yz_0_xyz_xyyz[k] = -g_yz_0_yz_xyyz[k] * ab_x + g_yz_0_yz_xxyyz[k];

                g_yz_0_xyz_xyzz[k] = -g_yz_0_yz_xyzz[k] * ab_x + g_yz_0_yz_xxyzz[k];

                g_yz_0_xyz_xzzz[k] = -g_yz_0_yz_xzzz[k] * ab_x + g_yz_0_yz_xxzzz[k];

                g_yz_0_xyz_yyyy[k] = -g_yz_0_yz_yyyy[k] * ab_x + g_yz_0_yz_xyyyy[k];

                g_yz_0_xyz_yyyz[k] = -g_yz_0_yz_yyyz[k] * ab_x + g_yz_0_yz_xyyyz[k];

                g_yz_0_xyz_yyzz[k] = -g_yz_0_yz_yyzz[k] * ab_x + g_yz_0_yz_xyyzz[k];

                g_yz_0_xyz_yzzz[k] = -g_yz_0_yz_yzzz[k] * ab_x + g_yz_0_yz_xyzzz[k];

                g_yz_0_xyz_zzzz[k] = -g_yz_0_yz_zzzz[k] * ab_x + g_yz_0_yz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 677 * ccomps * dcomps);

            auto g_yz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 679 * ccomps * dcomps);

            auto g_yz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 683 * ccomps * dcomps);

            auto g_yz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzz_xxxx, g_yz_0_xzz_xxxy, g_yz_0_xzz_xxxz, g_yz_0_xzz_xxyy, g_yz_0_xzz_xxyz, g_yz_0_xzz_xxzz, g_yz_0_xzz_xyyy, g_yz_0_xzz_xyyz, g_yz_0_xzz_xyzz, g_yz_0_xzz_xzzz, g_yz_0_xzz_yyyy, g_yz_0_xzz_yyyz, g_yz_0_xzz_yyzz, g_yz_0_xzz_yzzz, g_yz_0_xzz_zzzz, g_yz_0_zz_xxxx, g_yz_0_zz_xxxxx, g_yz_0_zz_xxxxy, g_yz_0_zz_xxxxz, g_yz_0_zz_xxxy, g_yz_0_zz_xxxyy, g_yz_0_zz_xxxyz, g_yz_0_zz_xxxz, g_yz_0_zz_xxxzz, g_yz_0_zz_xxyy, g_yz_0_zz_xxyyy, g_yz_0_zz_xxyyz, g_yz_0_zz_xxyz, g_yz_0_zz_xxyzz, g_yz_0_zz_xxzz, g_yz_0_zz_xxzzz, g_yz_0_zz_xyyy, g_yz_0_zz_xyyyy, g_yz_0_zz_xyyyz, g_yz_0_zz_xyyz, g_yz_0_zz_xyyzz, g_yz_0_zz_xyzz, g_yz_0_zz_xyzzz, g_yz_0_zz_xzzz, g_yz_0_zz_xzzzz, g_yz_0_zz_yyyy, g_yz_0_zz_yyyz, g_yz_0_zz_yyzz, g_yz_0_zz_yzzz, g_yz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzz_xxxx[k] = -g_yz_0_zz_xxxx[k] * ab_x + g_yz_0_zz_xxxxx[k];

                g_yz_0_xzz_xxxy[k] = -g_yz_0_zz_xxxy[k] * ab_x + g_yz_0_zz_xxxxy[k];

                g_yz_0_xzz_xxxz[k] = -g_yz_0_zz_xxxz[k] * ab_x + g_yz_0_zz_xxxxz[k];

                g_yz_0_xzz_xxyy[k] = -g_yz_0_zz_xxyy[k] * ab_x + g_yz_0_zz_xxxyy[k];

                g_yz_0_xzz_xxyz[k] = -g_yz_0_zz_xxyz[k] * ab_x + g_yz_0_zz_xxxyz[k];

                g_yz_0_xzz_xxzz[k] = -g_yz_0_zz_xxzz[k] * ab_x + g_yz_0_zz_xxxzz[k];

                g_yz_0_xzz_xyyy[k] = -g_yz_0_zz_xyyy[k] * ab_x + g_yz_0_zz_xxyyy[k];

                g_yz_0_xzz_xyyz[k] = -g_yz_0_zz_xyyz[k] * ab_x + g_yz_0_zz_xxyyz[k];

                g_yz_0_xzz_xyzz[k] = -g_yz_0_zz_xyzz[k] * ab_x + g_yz_0_zz_xxyzz[k];

                g_yz_0_xzz_xzzz[k] = -g_yz_0_zz_xzzz[k] * ab_x + g_yz_0_zz_xxzzz[k];

                g_yz_0_xzz_yyyy[k] = -g_yz_0_zz_yyyy[k] * ab_x + g_yz_0_zz_xyyyy[k];

                g_yz_0_xzz_yyyz[k] = -g_yz_0_zz_yyyz[k] * ab_x + g_yz_0_zz_xyyyz[k];

                g_yz_0_xzz_yyzz[k] = -g_yz_0_zz_yyzz[k] * ab_x + g_yz_0_zz_xyyzz[k];

                g_yz_0_xzz_yzzz[k] = -g_yz_0_zz_yzzz[k] * ab_x + g_yz_0_zz_xyzzz[k];

                g_yz_0_xzz_zzzz[k] = -g_yz_0_zz_zzzz[k] * ab_x + g_yz_0_zz_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 695 * ccomps * dcomps);

            auto g_yz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 699 * ccomps * dcomps);

            auto g_yz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 701 * ccomps * dcomps);

            auto g_yz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yy_xxxx, g_yz_0_yy_xxxxy, g_yz_0_yy_xxxy, g_yz_0_yy_xxxyy, g_yz_0_yy_xxxyz, g_yz_0_yy_xxxz, g_yz_0_yy_xxyy, g_yz_0_yy_xxyyy, g_yz_0_yy_xxyyz, g_yz_0_yy_xxyz, g_yz_0_yy_xxyzz, g_yz_0_yy_xxzz, g_yz_0_yy_xyyy, g_yz_0_yy_xyyyy, g_yz_0_yy_xyyyz, g_yz_0_yy_xyyz, g_yz_0_yy_xyyzz, g_yz_0_yy_xyzz, g_yz_0_yy_xyzzz, g_yz_0_yy_xzzz, g_yz_0_yy_yyyy, g_yz_0_yy_yyyyy, g_yz_0_yy_yyyyz, g_yz_0_yy_yyyz, g_yz_0_yy_yyyzz, g_yz_0_yy_yyzz, g_yz_0_yy_yyzzz, g_yz_0_yy_yzzz, g_yz_0_yy_yzzzz, g_yz_0_yy_zzzz, g_yz_0_yyy_xxxx, g_yz_0_yyy_xxxy, g_yz_0_yyy_xxxz, g_yz_0_yyy_xxyy, g_yz_0_yyy_xxyz, g_yz_0_yyy_xxzz, g_yz_0_yyy_xyyy, g_yz_0_yyy_xyyz, g_yz_0_yyy_xyzz, g_yz_0_yyy_xzzz, g_yz_0_yyy_yyyy, g_yz_0_yyy_yyyz, g_yz_0_yyy_yyzz, g_yz_0_yyy_yzzz, g_yz_0_yyy_zzzz, g_z_0_yy_xxxx, g_z_0_yy_xxxy, g_z_0_yy_xxxz, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxzz, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyzz, g_z_0_yy_xzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyzz, g_z_0_yy_yzzz, g_z_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyy_xxxx[k] = -g_z_0_yy_xxxx[k] - g_yz_0_yy_xxxx[k] * ab_y + g_yz_0_yy_xxxxy[k];

                g_yz_0_yyy_xxxy[k] = -g_z_0_yy_xxxy[k] - g_yz_0_yy_xxxy[k] * ab_y + g_yz_0_yy_xxxyy[k];

                g_yz_0_yyy_xxxz[k] = -g_z_0_yy_xxxz[k] - g_yz_0_yy_xxxz[k] * ab_y + g_yz_0_yy_xxxyz[k];

                g_yz_0_yyy_xxyy[k] = -g_z_0_yy_xxyy[k] - g_yz_0_yy_xxyy[k] * ab_y + g_yz_0_yy_xxyyy[k];

                g_yz_0_yyy_xxyz[k] = -g_z_0_yy_xxyz[k] - g_yz_0_yy_xxyz[k] * ab_y + g_yz_0_yy_xxyyz[k];

                g_yz_0_yyy_xxzz[k] = -g_z_0_yy_xxzz[k] - g_yz_0_yy_xxzz[k] * ab_y + g_yz_0_yy_xxyzz[k];

                g_yz_0_yyy_xyyy[k] = -g_z_0_yy_xyyy[k] - g_yz_0_yy_xyyy[k] * ab_y + g_yz_0_yy_xyyyy[k];

                g_yz_0_yyy_xyyz[k] = -g_z_0_yy_xyyz[k] - g_yz_0_yy_xyyz[k] * ab_y + g_yz_0_yy_xyyyz[k];

                g_yz_0_yyy_xyzz[k] = -g_z_0_yy_xyzz[k] - g_yz_0_yy_xyzz[k] * ab_y + g_yz_0_yy_xyyzz[k];

                g_yz_0_yyy_xzzz[k] = -g_z_0_yy_xzzz[k] - g_yz_0_yy_xzzz[k] * ab_y + g_yz_0_yy_xyzzz[k];

                g_yz_0_yyy_yyyy[k] = -g_z_0_yy_yyyy[k] - g_yz_0_yy_yyyy[k] * ab_y + g_yz_0_yy_yyyyy[k];

                g_yz_0_yyy_yyyz[k] = -g_z_0_yy_yyyz[k] - g_yz_0_yy_yyyz[k] * ab_y + g_yz_0_yy_yyyyz[k];

                g_yz_0_yyy_yyzz[k] = -g_z_0_yy_yyzz[k] - g_yz_0_yy_yyzz[k] * ab_y + g_yz_0_yy_yyyzz[k];

                g_yz_0_yyy_yzzz[k] = -g_z_0_yy_yzzz[k] - g_yz_0_yy_yzzz[k] * ab_y + g_yz_0_yy_yyzzz[k];

                g_yz_0_yyy_zzzz[k] = -g_z_0_yy_zzzz[k] - g_yz_0_yy_zzzz[k] * ab_y + g_yz_0_yy_yzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 707 * ccomps * dcomps);

            auto g_yz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 709 * ccomps * dcomps);

            auto g_yz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 713 * ccomps * dcomps);

            auto g_yz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyz_xxxx, g_yz_0_yyz_xxxy, g_yz_0_yyz_xxxz, g_yz_0_yyz_xxyy, g_yz_0_yyz_xxyz, g_yz_0_yyz_xxzz, g_yz_0_yyz_xyyy, g_yz_0_yyz_xyyz, g_yz_0_yyz_xyzz, g_yz_0_yyz_xzzz, g_yz_0_yyz_yyyy, g_yz_0_yyz_yyyz, g_yz_0_yyz_yyzz, g_yz_0_yyz_yzzz, g_yz_0_yyz_zzzz, g_yz_0_yz_xxxx, g_yz_0_yz_xxxxy, g_yz_0_yz_xxxy, g_yz_0_yz_xxxyy, g_yz_0_yz_xxxyz, g_yz_0_yz_xxxz, g_yz_0_yz_xxyy, g_yz_0_yz_xxyyy, g_yz_0_yz_xxyyz, g_yz_0_yz_xxyz, g_yz_0_yz_xxyzz, g_yz_0_yz_xxzz, g_yz_0_yz_xyyy, g_yz_0_yz_xyyyy, g_yz_0_yz_xyyyz, g_yz_0_yz_xyyz, g_yz_0_yz_xyyzz, g_yz_0_yz_xyzz, g_yz_0_yz_xyzzz, g_yz_0_yz_xzzz, g_yz_0_yz_yyyy, g_yz_0_yz_yyyyy, g_yz_0_yz_yyyyz, g_yz_0_yz_yyyz, g_yz_0_yz_yyyzz, g_yz_0_yz_yyzz, g_yz_0_yz_yyzzz, g_yz_0_yz_yzzz, g_yz_0_yz_yzzzz, g_yz_0_yz_zzzz, g_z_0_yz_xxxx, g_z_0_yz_xxxy, g_z_0_yz_xxxz, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxzz, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyzz, g_z_0_yz_xzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyzz, g_z_0_yz_yzzz, g_z_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyz_xxxx[k] = -g_z_0_yz_xxxx[k] - g_yz_0_yz_xxxx[k] * ab_y + g_yz_0_yz_xxxxy[k];

                g_yz_0_yyz_xxxy[k] = -g_z_0_yz_xxxy[k] - g_yz_0_yz_xxxy[k] * ab_y + g_yz_0_yz_xxxyy[k];

                g_yz_0_yyz_xxxz[k] = -g_z_0_yz_xxxz[k] - g_yz_0_yz_xxxz[k] * ab_y + g_yz_0_yz_xxxyz[k];

                g_yz_0_yyz_xxyy[k] = -g_z_0_yz_xxyy[k] - g_yz_0_yz_xxyy[k] * ab_y + g_yz_0_yz_xxyyy[k];

                g_yz_0_yyz_xxyz[k] = -g_z_0_yz_xxyz[k] - g_yz_0_yz_xxyz[k] * ab_y + g_yz_0_yz_xxyyz[k];

                g_yz_0_yyz_xxzz[k] = -g_z_0_yz_xxzz[k] - g_yz_0_yz_xxzz[k] * ab_y + g_yz_0_yz_xxyzz[k];

                g_yz_0_yyz_xyyy[k] = -g_z_0_yz_xyyy[k] - g_yz_0_yz_xyyy[k] * ab_y + g_yz_0_yz_xyyyy[k];

                g_yz_0_yyz_xyyz[k] = -g_z_0_yz_xyyz[k] - g_yz_0_yz_xyyz[k] * ab_y + g_yz_0_yz_xyyyz[k];

                g_yz_0_yyz_xyzz[k] = -g_z_0_yz_xyzz[k] - g_yz_0_yz_xyzz[k] * ab_y + g_yz_0_yz_xyyzz[k];

                g_yz_0_yyz_xzzz[k] = -g_z_0_yz_xzzz[k] - g_yz_0_yz_xzzz[k] * ab_y + g_yz_0_yz_xyzzz[k];

                g_yz_0_yyz_yyyy[k] = -g_z_0_yz_yyyy[k] - g_yz_0_yz_yyyy[k] * ab_y + g_yz_0_yz_yyyyy[k];

                g_yz_0_yyz_yyyz[k] = -g_z_0_yz_yyyz[k] - g_yz_0_yz_yyyz[k] * ab_y + g_yz_0_yz_yyyyz[k];

                g_yz_0_yyz_yyzz[k] = -g_z_0_yz_yyzz[k] - g_yz_0_yz_yyzz[k] * ab_y + g_yz_0_yz_yyyzz[k];

                g_yz_0_yyz_yzzz[k] = -g_z_0_yz_yzzz[k] - g_yz_0_yz_yzzz[k] * ab_y + g_yz_0_yz_yyzzz[k];

                g_yz_0_yyz_zzzz[k] = -g_z_0_yz_zzzz[k] - g_yz_0_yz_zzzz[k] * ab_y + g_yz_0_yz_yzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 725 * ccomps * dcomps);

            auto g_yz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 727 * ccomps * dcomps);

            auto g_yz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 729 * ccomps * dcomps);

            auto g_yz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 731 * ccomps * dcomps);

            auto g_yz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzz_xxxx, g_yz_0_yzz_xxxy, g_yz_0_yzz_xxxz, g_yz_0_yzz_xxyy, g_yz_0_yzz_xxyz, g_yz_0_yzz_xxzz, g_yz_0_yzz_xyyy, g_yz_0_yzz_xyyz, g_yz_0_yzz_xyzz, g_yz_0_yzz_xzzz, g_yz_0_yzz_yyyy, g_yz_0_yzz_yyyz, g_yz_0_yzz_yyzz, g_yz_0_yzz_yzzz, g_yz_0_yzz_zzzz, g_yz_0_zz_xxxx, g_yz_0_zz_xxxxy, g_yz_0_zz_xxxy, g_yz_0_zz_xxxyy, g_yz_0_zz_xxxyz, g_yz_0_zz_xxxz, g_yz_0_zz_xxyy, g_yz_0_zz_xxyyy, g_yz_0_zz_xxyyz, g_yz_0_zz_xxyz, g_yz_0_zz_xxyzz, g_yz_0_zz_xxzz, g_yz_0_zz_xyyy, g_yz_0_zz_xyyyy, g_yz_0_zz_xyyyz, g_yz_0_zz_xyyz, g_yz_0_zz_xyyzz, g_yz_0_zz_xyzz, g_yz_0_zz_xyzzz, g_yz_0_zz_xzzz, g_yz_0_zz_yyyy, g_yz_0_zz_yyyyy, g_yz_0_zz_yyyyz, g_yz_0_zz_yyyz, g_yz_0_zz_yyyzz, g_yz_0_zz_yyzz, g_yz_0_zz_yyzzz, g_yz_0_zz_yzzz, g_yz_0_zz_yzzzz, g_yz_0_zz_zzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzz_xxxx[k] = -g_z_0_zz_xxxx[k] - g_yz_0_zz_xxxx[k] * ab_y + g_yz_0_zz_xxxxy[k];

                g_yz_0_yzz_xxxy[k] = -g_z_0_zz_xxxy[k] - g_yz_0_zz_xxxy[k] * ab_y + g_yz_0_zz_xxxyy[k];

                g_yz_0_yzz_xxxz[k] = -g_z_0_zz_xxxz[k] - g_yz_0_zz_xxxz[k] * ab_y + g_yz_0_zz_xxxyz[k];

                g_yz_0_yzz_xxyy[k] = -g_z_0_zz_xxyy[k] - g_yz_0_zz_xxyy[k] * ab_y + g_yz_0_zz_xxyyy[k];

                g_yz_0_yzz_xxyz[k] = -g_z_0_zz_xxyz[k] - g_yz_0_zz_xxyz[k] * ab_y + g_yz_0_zz_xxyyz[k];

                g_yz_0_yzz_xxzz[k] = -g_z_0_zz_xxzz[k] - g_yz_0_zz_xxzz[k] * ab_y + g_yz_0_zz_xxyzz[k];

                g_yz_0_yzz_xyyy[k] = -g_z_0_zz_xyyy[k] - g_yz_0_zz_xyyy[k] * ab_y + g_yz_0_zz_xyyyy[k];

                g_yz_0_yzz_xyyz[k] = -g_z_0_zz_xyyz[k] - g_yz_0_zz_xyyz[k] * ab_y + g_yz_0_zz_xyyyz[k];

                g_yz_0_yzz_xyzz[k] = -g_z_0_zz_xyzz[k] - g_yz_0_zz_xyzz[k] * ab_y + g_yz_0_zz_xyyzz[k];

                g_yz_0_yzz_xzzz[k] = -g_z_0_zz_xzzz[k] - g_yz_0_zz_xzzz[k] * ab_y + g_yz_0_zz_xyzzz[k];

                g_yz_0_yzz_yyyy[k] = -g_z_0_zz_yyyy[k] - g_yz_0_zz_yyyy[k] * ab_y + g_yz_0_zz_yyyyy[k];

                g_yz_0_yzz_yyyz[k] = -g_z_0_zz_yyyz[k] - g_yz_0_zz_yyyz[k] * ab_y + g_yz_0_zz_yyyyz[k];

                g_yz_0_yzz_yyzz[k] = -g_z_0_zz_yyzz[k] - g_yz_0_zz_yyzz[k] * ab_y + g_yz_0_zz_yyyzz[k];

                g_yz_0_yzz_yzzz[k] = -g_z_0_zz_yzzz[k] - g_yz_0_zz_yzzz[k] * ab_y + g_yz_0_zz_yyzzz[k];

                g_yz_0_yzz_zzzz[k] = -g_z_0_zz_zzzz[k] - g_yz_0_zz_zzzz[k] * ab_y + g_yz_0_zz_yzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 737 * ccomps * dcomps);

            auto g_yz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 739 * ccomps * dcomps);

            auto g_yz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 743 * ccomps * dcomps);

            auto g_yz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xxxx, g_y_0_zz_xxxy, g_y_0_zz_xxxz, g_y_0_zz_xxyy, g_y_0_zz_xxyz, g_y_0_zz_xxzz, g_y_0_zz_xyyy, g_y_0_zz_xyyz, g_y_0_zz_xyzz, g_y_0_zz_xzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyz, g_y_0_zz_yyzz, g_y_0_zz_yzzz, g_y_0_zz_zzzz, g_yz_0_zz_xxxx, g_yz_0_zz_xxxxz, g_yz_0_zz_xxxy, g_yz_0_zz_xxxyz, g_yz_0_zz_xxxz, g_yz_0_zz_xxxzz, g_yz_0_zz_xxyy, g_yz_0_zz_xxyyz, g_yz_0_zz_xxyz, g_yz_0_zz_xxyzz, g_yz_0_zz_xxzz, g_yz_0_zz_xxzzz, g_yz_0_zz_xyyy, g_yz_0_zz_xyyyz, g_yz_0_zz_xyyz, g_yz_0_zz_xyyzz, g_yz_0_zz_xyzz, g_yz_0_zz_xyzzz, g_yz_0_zz_xzzz, g_yz_0_zz_xzzzz, g_yz_0_zz_yyyy, g_yz_0_zz_yyyyz, g_yz_0_zz_yyyz, g_yz_0_zz_yyyzz, g_yz_0_zz_yyzz, g_yz_0_zz_yyzzz, g_yz_0_zz_yzzz, g_yz_0_zz_yzzzz, g_yz_0_zz_zzzz, g_yz_0_zz_zzzzz, g_yz_0_zzz_xxxx, g_yz_0_zzz_xxxy, g_yz_0_zzz_xxxz, g_yz_0_zzz_xxyy, g_yz_0_zzz_xxyz, g_yz_0_zzz_xxzz, g_yz_0_zzz_xyyy, g_yz_0_zzz_xyyz, g_yz_0_zzz_xyzz, g_yz_0_zzz_xzzz, g_yz_0_zzz_yyyy, g_yz_0_zzz_yyyz, g_yz_0_zzz_yyzz, g_yz_0_zzz_yzzz, g_yz_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzz_xxxx[k] = -g_y_0_zz_xxxx[k] - g_yz_0_zz_xxxx[k] * ab_z + g_yz_0_zz_xxxxz[k];

                g_yz_0_zzz_xxxy[k] = -g_y_0_zz_xxxy[k] - g_yz_0_zz_xxxy[k] * ab_z + g_yz_0_zz_xxxyz[k];

                g_yz_0_zzz_xxxz[k] = -g_y_0_zz_xxxz[k] - g_yz_0_zz_xxxz[k] * ab_z + g_yz_0_zz_xxxzz[k];

                g_yz_0_zzz_xxyy[k] = -g_y_0_zz_xxyy[k] - g_yz_0_zz_xxyy[k] * ab_z + g_yz_0_zz_xxyyz[k];

                g_yz_0_zzz_xxyz[k] = -g_y_0_zz_xxyz[k] - g_yz_0_zz_xxyz[k] * ab_z + g_yz_0_zz_xxyzz[k];

                g_yz_0_zzz_xxzz[k] = -g_y_0_zz_xxzz[k] - g_yz_0_zz_xxzz[k] * ab_z + g_yz_0_zz_xxzzz[k];

                g_yz_0_zzz_xyyy[k] = -g_y_0_zz_xyyy[k] - g_yz_0_zz_xyyy[k] * ab_z + g_yz_0_zz_xyyyz[k];

                g_yz_0_zzz_xyyz[k] = -g_y_0_zz_xyyz[k] - g_yz_0_zz_xyyz[k] * ab_z + g_yz_0_zz_xyyzz[k];

                g_yz_0_zzz_xyzz[k] = -g_y_0_zz_xyzz[k] - g_yz_0_zz_xyzz[k] * ab_z + g_yz_0_zz_xyzzz[k];

                g_yz_0_zzz_xzzz[k] = -g_y_0_zz_xzzz[k] - g_yz_0_zz_xzzz[k] * ab_z + g_yz_0_zz_xzzzz[k];

                g_yz_0_zzz_yyyy[k] = -g_y_0_zz_yyyy[k] - g_yz_0_zz_yyyy[k] * ab_z + g_yz_0_zz_yyyyz[k];

                g_yz_0_zzz_yyyz[k] = -g_y_0_zz_yyyz[k] - g_yz_0_zz_yyyz[k] * ab_z + g_yz_0_zz_yyyzz[k];

                g_yz_0_zzz_yyzz[k] = -g_y_0_zz_yyzz[k] - g_yz_0_zz_yyzz[k] * ab_z + g_yz_0_zz_yyzzz[k];

                g_yz_0_zzz_yzzz[k] = -g_y_0_zz_yzzz[k] - g_yz_0_zz_yzzz[k] * ab_z + g_yz_0_zz_yzzzz[k];

                g_yz_0_zzz_zzzz[k] = -g_y_0_zz_zzzz[k] - g_yz_0_zz_zzzz[k] * ab_z + g_yz_0_zz_zzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxx_xxxx = cbuffer.data(fg_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_xxx_xxxy = cbuffer.data(fg_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_xxx_xxxz = cbuffer.data(fg_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_xxx_xxyy = cbuffer.data(fg_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_xxx_xxyz = cbuffer.data(fg_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_xxx_xxzz = cbuffer.data(fg_geom_20_off + 755 * ccomps * dcomps);

            auto g_zz_0_xxx_xyyy = cbuffer.data(fg_geom_20_off + 756 * ccomps * dcomps);

            auto g_zz_0_xxx_xyyz = cbuffer.data(fg_geom_20_off + 757 * ccomps * dcomps);

            auto g_zz_0_xxx_xyzz = cbuffer.data(fg_geom_20_off + 758 * ccomps * dcomps);

            auto g_zz_0_xxx_xzzz = cbuffer.data(fg_geom_20_off + 759 * ccomps * dcomps);

            auto g_zz_0_xxx_yyyy = cbuffer.data(fg_geom_20_off + 760 * ccomps * dcomps);

            auto g_zz_0_xxx_yyyz = cbuffer.data(fg_geom_20_off + 761 * ccomps * dcomps);

            auto g_zz_0_xxx_yyzz = cbuffer.data(fg_geom_20_off + 762 * ccomps * dcomps);

            auto g_zz_0_xxx_yzzz = cbuffer.data(fg_geom_20_off + 763 * ccomps * dcomps);

            auto g_zz_0_xxx_zzzz = cbuffer.data(fg_geom_20_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xx_xxxx, g_zz_0_xx_xxxxx, g_zz_0_xx_xxxxy, g_zz_0_xx_xxxxz, g_zz_0_xx_xxxy, g_zz_0_xx_xxxyy, g_zz_0_xx_xxxyz, g_zz_0_xx_xxxz, g_zz_0_xx_xxxzz, g_zz_0_xx_xxyy, g_zz_0_xx_xxyyy, g_zz_0_xx_xxyyz, g_zz_0_xx_xxyz, g_zz_0_xx_xxyzz, g_zz_0_xx_xxzz, g_zz_0_xx_xxzzz, g_zz_0_xx_xyyy, g_zz_0_xx_xyyyy, g_zz_0_xx_xyyyz, g_zz_0_xx_xyyz, g_zz_0_xx_xyyzz, g_zz_0_xx_xyzz, g_zz_0_xx_xyzzz, g_zz_0_xx_xzzz, g_zz_0_xx_xzzzz, g_zz_0_xx_yyyy, g_zz_0_xx_yyyz, g_zz_0_xx_yyzz, g_zz_0_xx_yzzz, g_zz_0_xx_zzzz, g_zz_0_xxx_xxxx, g_zz_0_xxx_xxxy, g_zz_0_xxx_xxxz, g_zz_0_xxx_xxyy, g_zz_0_xxx_xxyz, g_zz_0_xxx_xxzz, g_zz_0_xxx_xyyy, g_zz_0_xxx_xyyz, g_zz_0_xxx_xyzz, g_zz_0_xxx_xzzz, g_zz_0_xxx_yyyy, g_zz_0_xxx_yyyz, g_zz_0_xxx_yyzz, g_zz_0_xxx_yzzz, g_zz_0_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxx_xxxx[k] = -g_zz_0_xx_xxxx[k] * ab_x + g_zz_0_xx_xxxxx[k];

                g_zz_0_xxx_xxxy[k] = -g_zz_0_xx_xxxy[k] * ab_x + g_zz_0_xx_xxxxy[k];

                g_zz_0_xxx_xxxz[k] = -g_zz_0_xx_xxxz[k] * ab_x + g_zz_0_xx_xxxxz[k];

                g_zz_0_xxx_xxyy[k] = -g_zz_0_xx_xxyy[k] * ab_x + g_zz_0_xx_xxxyy[k];

                g_zz_0_xxx_xxyz[k] = -g_zz_0_xx_xxyz[k] * ab_x + g_zz_0_xx_xxxyz[k];

                g_zz_0_xxx_xxzz[k] = -g_zz_0_xx_xxzz[k] * ab_x + g_zz_0_xx_xxxzz[k];

                g_zz_0_xxx_xyyy[k] = -g_zz_0_xx_xyyy[k] * ab_x + g_zz_0_xx_xxyyy[k];

                g_zz_0_xxx_xyyz[k] = -g_zz_0_xx_xyyz[k] * ab_x + g_zz_0_xx_xxyyz[k];

                g_zz_0_xxx_xyzz[k] = -g_zz_0_xx_xyzz[k] * ab_x + g_zz_0_xx_xxyzz[k];

                g_zz_0_xxx_xzzz[k] = -g_zz_0_xx_xzzz[k] * ab_x + g_zz_0_xx_xxzzz[k];

                g_zz_0_xxx_yyyy[k] = -g_zz_0_xx_yyyy[k] * ab_x + g_zz_0_xx_xyyyy[k];

                g_zz_0_xxx_yyyz[k] = -g_zz_0_xx_yyyz[k] * ab_x + g_zz_0_xx_xyyyz[k];

                g_zz_0_xxx_yyzz[k] = -g_zz_0_xx_yyzz[k] * ab_x + g_zz_0_xx_xyyzz[k];

                g_zz_0_xxx_yzzz[k] = -g_zz_0_xx_yzzz[k] * ab_x + g_zz_0_xx_xyzzz[k];

                g_zz_0_xxx_zzzz[k] = -g_zz_0_xx_zzzz[k] * ab_x + g_zz_0_xx_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxy_xxxx = cbuffer.data(fg_geom_20_off + 765 * ccomps * dcomps);

            auto g_zz_0_xxy_xxxy = cbuffer.data(fg_geom_20_off + 766 * ccomps * dcomps);

            auto g_zz_0_xxy_xxxz = cbuffer.data(fg_geom_20_off + 767 * ccomps * dcomps);

            auto g_zz_0_xxy_xxyy = cbuffer.data(fg_geom_20_off + 768 * ccomps * dcomps);

            auto g_zz_0_xxy_xxyz = cbuffer.data(fg_geom_20_off + 769 * ccomps * dcomps);

            auto g_zz_0_xxy_xxzz = cbuffer.data(fg_geom_20_off + 770 * ccomps * dcomps);

            auto g_zz_0_xxy_xyyy = cbuffer.data(fg_geom_20_off + 771 * ccomps * dcomps);

            auto g_zz_0_xxy_xyyz = cbuffer.data(fg_geom_20_off + 772 * ccomps * dcomps);

            auto g_zz_0_xxy_xyzz = cbuffer.data(fg_geom_20_off + 773 * ccomps * dcomps);

            auto g_zz_0_xxy_xzzz = cbuffer.data(fg_geom_20_off + 774 * ccomps * dcomps);

            auto g_zz_0_xxy_yyyy = cbuffer.data(fg_geom_20_off + 775 * ccomps * dcomps);

            auto g_zz_0_xxy_yyyz = cbuffer.data(fg_geom_20_off + 776 * ccomps * dcomps);

            auto g_zz_0_xxy_yyzz = cbuffer.data(fg_geom_20_off + 777 * ccomps * dcomps);

            auto g_zz_0_xxy_yzzz = cbuffer.data(fg_geom_20_off + 778 * ccomps * dcomps);

            auto g_zz_0_xxy_zzzz = cbuffer.data(fg_geom_20_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxy_xxxx, g_zz_0_xxy_xxxy, g_zz_0_xxy_xxxz, g_zz_0_xxy_xxyy, g_zz_0_xxy_xxyz, g_zz_0_xxy_xxzz, g_zz_0_xxy_xyyy, g_zz_0_xxy_xyyz, g_zz_0_xxy_xyzz, g_zz_0_xxy_xzzz, g_zz_0_xxy_yyyy, g_zz_0_xxy_yyyz, g_zz_0_xxy_yyzz, g_zz_0_xxy_yzzz, g_zz_0_xxy_zzzz, g_zz_0_xy_xxxx, g_zz_0_xy_xxxxx, g_zz_0_xy_xxxxy, g_zz_0_xy_xxxxz, g_zz_0_xy_xxxy, g_zz_0_xy_xxxyy, g_zz_0_xy_xxxyz, g_zz_0_xy_xxxz, g_zz_0_xy_xxxzz, g_zz_0_xy_xxyy, g_zz_0_xy_xxyyy, g_zz_0_xy_xxyyz, g_zz_0_xy_xxyz, g_zz_0_xy_xxyzz, g_zz_0_xy_xxzz, g_zz_0_xy_xxzzz, g_zz_0_xy_xyyy, g_zz_0_xy_xyyyy, g_zz_0_xy_xyyyz, g_zz_0_xy_xyyz, g_zz_0_xy_xyyzz, g_zz_0_xy_xyzz, g_zz_0_xy_xyzzz, g_zz_0_xy_xzzz, g_zz_0_xy_xzzzz, g_zz_0_xy_yyyy, g_zz_0_xy_yyyz, g_zz_0_xy_yyzz, g_zz_0_xy_yzzz, g_zz_0_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxy_xxxx[k] = -g_zz_0_xy_xxxx[k] * ab_x + g_zz_0_xy_xxxxx[k];

                g_zz_0_xxy_xxxy[k] = -g_zz_0_xy_xxxy[k] * ab_x + g_zz_0_xy_xxxxy[k];

                g_zz_0_xxy_xxxz[k] = -g_zz_0_xy_xxxz[k] * ab_x + g_zz_0_xy_xxxxz[k];

                g_zz_0_xxy_xxyy[k] = -g_zz_0_xy_xxyy[k] * ab_x + g_zz_0_xy_xxxyy[k];

                g_zz_0_xxy_xxyz[k] = -g_zz_0_xy_xxyz[k] * ab_x + g_zz_0_xy_xxxyz[k];

                g_zz_0_xxy_xxzz[k] = -g_zz_0_xy_xxzz[k] * ab_x + g_zz_0_xy_xxxzz[k];

                g_zz_0_xxy_xyyy[k] = -g_zz_0_xy_xyyy[k] * ab_x + g_zz_0_xy_xxyyy[k];

                g_zz_0_xxy_xyyz[k] = -g_zz_0_xy_xyyz[k] * ab_x + g_zz_0_xy_xxyyz[k];

                g_zz_0_xxy_xyzz[k] = -g_zz_0_xy_xyzz[k] * ab_x + g_zz_0_xy_xxyzz[k];

                g_zz_0_xxy_xzzz[k] = -g_zz_0_xy_xzzz[k] * ab_x + g_zz_0_xy_xxzzz[k];

                g_zz_0_xxy_yyyy[k] = -g_zz_0_xy_yyyy[k] * ab_x + g_zz_0_xy_xyyyy[k];

                g_zz_0_xxy_yyyz[k] = -g_zz_0_xy_yyyz[k] * ab_x + g_zz_0_xy_xyyyz[k];

                g_zz_0_xxy_yyzz[k] = -g_zz_0_xy_yyzz[k] * ab_x + g_zz_0_xy_xyyzz[k];

                g_zz_0_xxy_yzzz[k] = -g_zz_0_xy_yzzz[k] * ab_x + g_zz_0_xy_xyzzz[k];

                g_zz_0_xxy_zzzz[k] = -g_zz_0_xy_zzzz[k] * ab_x + g_zz_0_xy_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxz_xxxx = cbuffer.data(fg_geom_20_off + 780 * ccomps * dcomps);

            auto g_zz_0_xxz_xxxy = cbuffer.data(fg_geom_20_off + 781 * ccomps * dcomps);

            auto g_zz_0_xxz_xxxz = cbuffer.data(fg_geom_20_off + 782 * ccomps * dcomps);

            auto g_zz_0_xxz_xxyy = cbuffer.data(fg_geom_20_off + 783 * ccomps * dcomps);

            auto g_zz_0_xxz_xxyz = cbuffer.data(fg_geom_20_off + 784 * ccomps * dcomps);

            auto g_zz_0_xxz_xxzz = cbuffer.data(fg_geom_20_off + 785 * ccomps * dcomps);

            auto g_zz_0_xxz_xyyy = cbuffer.data(fg_geom_20_off + 786 * ccomps * dcomps);

            auto g_zz_0_xxz_xyyz = cbuffer.data(fg_geom_20_off + 787 * ccomps * dcomps);

            auto g_zz_0_xxz_xyzz = cbuffer.data(fg_geom_20_off + 788 * ccomps * dcomps);

            auto g_zz_0_xxz_xzzz = cbuffer.data(fg_geom_20_off + 789 * ccomps * dcomps);

            auto g_zz_0_xxz_yyyy = cbuffer.data(fg_geom_20_off + 790 * ccomps * dcomps);

            auto g_zz_0_xxz_yyyz = cbuffer.data(fg_geom_20_off + 791 * ccomps * dcomps);

            auto g_zz_0_xxz_yyzz = cbuffer.data(fg_geom_20_off + 792 * ccomps * dcomps);

            auto g_zz_0_xxz_yzzz = cbuffer.data(fg_geom_20_off + 793 * ccomps * dcomps);

            auto g_zz_0_xxz_zzzz = cbuffer.data(fg_geom_20_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxz_xxxx, g_zz_0_xxz_xxxy, g_zz_0_xxz_xxxz, g_zz_0_xxz_xxyy, g_zz_0_xxz_xxyz, g_zz_0_xxz_xxzz, g_zz_0_xxz_xyyy, g_zz_0_xxz_xyyz, g_zz_0_xxz_xyzz, g_zz_0_xxz_xzzz, g_zz_0_xxz_yyyy, g_zz_0_xxz_yyyz, g_zz_0_xxz_yyzz, g_zz_0_xxz_yzzz, g_zz_0_xxz_zzzz, g_zz_0_xz_xxxx, g_zz_0_xz_xxxxx, g_zz_0_xz_xxxxy, g_zz_0_xz_xxxxz, g_zz_0_xz_xxxy, g_zz_0_xz_xxxyy, g_zz_0_xz_xxxyz, g_zz_0_xz_xxxz, g_zz_0_xz_xxxzz, g_zz_0_xz_xxyy, g_zz_0_xz_xxyyy, g_zz_0_xz_xxyyz, g_zz_0_xz_xxyz, g_zz_0_xz_xxyzz, g_zz_0_xz_xxzz, g_zz_0_xz_xxzzz, g_zz_0_xz_xyyy, g_zz_0_xz_xyyyy, g_zz_0_xz_xyyyz, g_zz_0_xz_xyyz, g_zz_0_xz_xyyzz, g_zz_0_xz_xyzz, g_zz_0_xz_xyzzz, g_zz_0_xz_xzzz, g_zz_0_xz_xzzzz, g_zz_0_xz_yyyy, g_zz_0_xz_yyyz, g_zz_0_xz_yyzz, g_zz_0_xz_yzzz, g_zz_0_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxz_xxxx[k] = -g_zz_0_xz_xxxx[k] * ab_x + g_zz_0_xz_xxxxx[k];

                g_zz_0_xxz_xxxy[k] = -g_zz_0_xz_xxxy[k] * ab_x + g_zz_0_xz_xxxxy[k];

                g_zz_0_xxz_xxxz[k] = -g_zz_0_xz_xxxz[k] * ab_x + g_zz_0_xz_xxxxz[k];

                g_zz_0_xxz_xxyy[k] = -g_zz_0_xz_xxyy[k] * ab_x + g_zz_0_xz_xxxyy[k];

                g_zz_0_xxz_xxyz[k] = -g_zz_0_xz_xxyz[k] * ab_x + g_zz_0_xz_xxxyz[k];

                g_zz_0_xxz_xxzz[k] = -g_zz_0_xz_xxzz[k] * ab_x + g_zz_0_xz_xxxzz[k];

                g_zz_0_xxz_xyyy[k] = -g_zz_0_xz_xyyy[k] * ab_x + g_zz_0_xz_xxyyy[k];

                g_zz_0_xxz_xyyz[k] = -g_zz_0_xz_xyyz[k] * ab_x + g_zz_0_xz_xxyyz[k];

                g_zz_0_xxz_xyzz[k] = -g_zz_0_xz_xyzz[k] * ab_x + g_zz_0_xz_xxyzz[k];

                g_zz_0_xxz_xzzz[k] = -g_zz_0_xz_xzzz[k] * ab_x + g_zz_0_xz_xxzzz[k];

                g_zz_0_xxz_yyyy[k] = -g_zz_0_xz_yyyy[k] * ab_x + g_zz_0_xz_xyyyy[k];

                g_zz_0_xxz_yyyz[k] = -g_zz_0_xz_yyyz[k] * ab_x + g_zz_0_xz_xyyyz[k];

                g_zz_0_xxz_yyzz[k] = -g_zz_0_xz_yyzz[k] * ab_x + g_zz_0_xz_xyyzz[k];

                g_zz_0_xxz_yzzz[k] = -g_zz_0_xz_yzzz[k] * ab_x + g_zz_0_xz_xyzzz[k];

                g_zz_0_xxz_zzzz[k] = -g_zz_0_xz_zzzz[k] * ab_x + g_zz_0_xz_xzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyy_xxxx = cbuffer.data(fg_geom_20_off + 795 * ccomps * dcomps);

            auto g_zz_0_xyy_xxxy = cbuffer.data(fg_geom_20_off + 796 * ccomps * dcomps);

            auto g_zz_0_xyy_xxxz = cbuffer.data(fg_geom_20_off + 797 * ccomps * dcomps);

            auto g_zz_0_xyy_xxyy = cbuffer.data(fg_geom_20_off + 798 * ccomps * dcomps);

            auto g_zz_0_xyy_xxyz = cbuffer.data(fg_geom_20_off + 799 * ccomps * dcomps);

            auto g_zz_0_xyy_xxzz = cbuffer.data(fg_geom_20_off + 800 * ccomps * dcomps);

            auto g_zz_0_xyy_xyyy = cbuffer.data(fg_geom_20_off + 801 * ccomps * dcomps);

            auto g_zz_0_xyy_xyyz = cbuffer.data(fg_geom_20_off + 802 * ccomps * dcomps);

            auto g_zz_0_xyy_xyzz = cbuffer.data(fg_geom_20_off + 803 * ccomps * dcomps);

            auto g_zz_0_xyy_xzzz = cbuffer.data(fg_geom_20_off + 804 * ccomps * dcomps);

            auto g_zz_0_xyy_yyyy = cbuffer.data(fg_geom_20_off + 805 * ccomps * dcomps);

            auto g_zz_0_xyy_yyyz = cbuffer.data(fg_geom_20_off + 806 * ccomps * dcomps);

            auto g_zz_0_xyy_yyzz = cbuffer.data(fg_geom_20_off + 807 * ccomps * dcomps);

            auto g_zz_0_xyy_yzzz = cbuffer.data(fg_geom_20_off + 808 * ccomps * dcomps);

            auto g_zz_0_xyy_zzzz = cbuffer.data(fg_geom_20_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyy_xxxx, g_zz_0_xyy_xxxy, g_zz_0_xyy_xxxz, g_zz_0_xyy_xxyy, g_zz_0_xyy_xxyz, g_zz_0_xyy_xxzz, g_zz_0_xyy_xyyy, g_zz_0_xyy_xyyz, g_zz_0_xyy_xyzz, g_zz_0_xyy_xzzz, g_zz_0_xyy_yyyy, g_zz_0_xyy_yyyz, g_zz_0_xyy_yyzz, g_zz_0_xyy_yzzz, g_zz_0_xyy_zzzz, g_zz_0_yy_xxxx, g_zz_0_yy_xxxxx, g_zz_0_yy_xxxxy, g_zz_0_yy_xxxxz, g_zz_0_yy_xxxy, g_zz_0_yy_xxxyy, g_zz_0_yy_xxxyz, g_zz_0_yy_xxxz, g_zz_0_yy_xxxzz, g_zz_0_yy_xxyy, g_zz_0_yy_xxyyy, g_zz_0_yy_xxyyz, g_zz_0_yy_xxyz, g_zz_0_yy_xxyzz, g_zz_0_yy_xxzz, g_zz_0_yy_xxzzz, g_zz_0_yy_xyyy, g_zz_0_yy_xyyyy, g_zz_0_yy_xyyyz, g_zz_0_yy_xyyz, g_zz_0_yy_xyyzz, g_zz_0_yy_xyzz, g_zz_0_yy_xyzzz, g_zz_0_yy_xzzz, g_zz_0_yy_xzzzz, g_zz_0_yy_yyyy, g_zz_0_yy_yyyz, g_zz_0_yy_yyzz, g_zz_0_yy_yzzz, g_zz_0_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyy_xxxx[k] = -g_zz_0_yy_xxxx[k] * ab_x + g_zz_0_yy_xxxxx[k];

                g_zz_0_xyy_xxxy[k] = -g_zz_0_yy_xxxy[k] * ab_x + g_zz_0_yy_xxxxy[k];

                g_zz_0_xyy_xxxz[k] = -g_zz_0_yy_xxxz[k] * ab_x + g_zz_0_yy_xxxxz[k];

                g_zz_0_xyy_xxyy[k] = -g_zz_0_yy_xxyy[k] * ab_x + g_zz_0_yy_xxxyy[k];

                g_zz_0_xyy_xxyz[k] = -g_zz_0_yy_xxyz[k] * ab_x + g_zz_0_yy_xxxyz[k];

                g_zz_0_xyy_xxzz[k] = -g_zz_0_yy_xxzz[k] * ab_x + g_zz_0_yy_xxxzz[k];

                g_zz_0_xyy_xyyy[k] = -g_zz_0_yy_xyyy[k] * ab_x + g_zz_0_yy_xxyyy[k];

                g_zz_0_xyy_xyyz[k] = -g_zz_0_yy_xyyz[k] * ab_x + g_zz_0_yy_xxyyz[k];

                g_zz_0_xyy_xyzz[k] = -g_zz_0_yy_xyzz[k] * ab_x + g_zz_0_yy_xxyzz[k];

                g_zz_0_xyy_xzzz[k] = -g_zz_0_yy_xzzz[k] * ab_x + g_zz_0_yy_xxzzz[k];

                g_zz_0_xyy_yyyy[k] = -g_zz_0_yy_yyyy[k] * ab_x + g_zz_0_yy_xyyyy[k];

                g_zz_0_xyy_yyyz[k] = -g_zz_0_yy_yyyz[k] * ab_x + g_zz_0_yy_xyyyz[k];

                g_zz_0_xyy_yyzz[k] = -g_zz_0_yy_yyzz[k] * ab_x + g_zz_0_yy_xyyzz[k];

                g_zz_0_xyy_yzzz[k] = -g_zz_0_yy_yzzz[k] * ab_x + g_zz_0_yy_xyzzz[k];

                g_zz_0_xyy_zzzz[k] = -g_zz_0_yy_zzzz[k] * ab_x + g_zz_0_yy_xzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyz_xxxx = cbuffer.data(fg_geom_20_off + 810 * ccomps * dcomps);

            auto g_zz_0_xyz_xxxy = cbuffer.data(fg_geom_20_off + 811 * ccomps * dcomps);

            auto g_zz_0_xyz_xxxz = cbuffer.data(fg_geom_20_off + 812 * ccomps * dcomps);

            auto g_zz_0_xyz_xxyy = cbuffer.data(fg_geom_20_off + 813 * ccomps * dcomps);

            auto g_zz_0_xyz_xxyz = cbuffer.data(fg_geom_20_off + 814 * ccomps * dcomps);

            auto g_zz_0_xyz_xxzz = cbuffer.data(fg_geom_20_off + 815 * ccomps * dcomps);

            auto g_zz_0_xyz_xyyy = cbuffer.data(fg_geom_20_off + 816 * ccomps * dcomps);

            auto g_zz_0_xyz_xyyz = cbuffer.data(fg_geom_20_off + 817 * ccomps * dcomps);

            auto g_zz_0_xyz_xyzz = cbuffer.data(fg_geom_20_off + 818 * ccomps * dcomps);

            auto g_zz_0_xyz_xzzz = cbuffer.data(fg_geom_20_off + 819 * ccomps * dcomps);

            auto g_zz_0_xyz_yyyy = cbuffer.data(fg_geom_20_off + 820 * ccomps * dcomps);

            auto g_zz_0_xyz_yyyz = cbuffer.data(fg_geom_20_off + 821 * ccomps * dcomps);

            auto g_zz_0_xyz_yyzz = cbuffer.data(fg_geom_20_off + 822 * ccomps * dcomps);

            auto g_zz_0_xyz_yzzz = cbuffer.data(fg_geom_20_off + 823 * ccomps * dcomps);

            auto g_zz_0_xyz_zzzz = cbuffer.data(fg_geom_20_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyz_xxxx, g_zz_0_xyz_xxxy, g_zz_0_xyz_xxxz, g_zz_0_xyz_xxyy, g_zz_0_xyz_xxyz, g_zz_0_xyz_xxzz, g_zz_0_xyz_xyyy, g_zz_0_xyz_xyyz, g_zz_0_xyz_xyzz, g_zz_0_xyz_xzzz, g_zz_0_xyz_yyyy, g_zz_0_xyz_yyyz, g_zz_0_xyz_yyzz, g_zz_0_xyz_yzzz, g_zz_0_xyz_zzzz, g_zz_0_yz_xxxx, g_zz_0_yz_xxxxx, g_zz_0_yz_xxxxy, g_zz_0_yz_xxxxz, g_zz_0_yz_xxxy, g_zz_0_yz_xxxyy, g_zz_0_yz_xxxyz, g_zz_0_yz_xxxz, g_zz_0_yz_xxxzz, g_zz_0_yz_xxyy, g_zz_0_yz_xxyyy, g_zz_0_yz_xxyyz, g_zz_0_yz_xxyz, g_zz_0_yz_xxyzz, g_zz_0_yz_xxzz, g_zz_0_yz_xxzzz, g_zz_0_yz_xyyy, g_zz_0_yz_xyyyy, g_zz_0_yz_xyyyz, g_zz_0_yz_xyyz, g_zz_0_yz_xyyzz, g_zz_0_yz_xyzz, g_zz_0_yz_xyzzz, g_zz_0_yz_xzzz, g_zz_0_yz_xzzzz, g_zz_0_yz_yyyy, g_zz_0_yz_yyyz, g_zz_0_yz_yyzz, g_zz_0_yz_yzzz, g_zz_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyz_xxxx[k] = -g_zz_0_yz_xxxx[k] * ab_x + g_zz_0_yz_xxxxx[k];

                g_zz_0_xyz_xxxy[k] = -g_zz_0_yz_xxxy[k] * ab_x + g_zz_0_yz_xxxxy[k];

                g_zz_0_xyz_xxxz[k] = -g_zz_0_yz_xxxz[k] * ab_x + g_zz_0_yz_xxxxz[k];

                g_zz_0_xyz_xxyy[k] = -g_zz_0_yz_xxyy[k] * ab_x + g_zz_0_yz_xxxyy[k];

                g_zz_0_xyz_xxyz[k] = -g_zz_0_yz_xxyz[k] * ab_x + g_zz_0_yz_xxxyz[k];

                g_zz_0_xyz_xxzz[k] = -g_zz_0_yz_xxzz[k] * ab_x + g_zz_0_yz_xxxzz[k];

                g_zz_0_xyz_xyyy[k] = -g_zz_0_yz_xyyy[k] * ab_x + g_zz_0_yz_xxyyy[k];

                g_zz_0_xyz_xyyz[k] = -g_zz_0_yz_xyyz[k] * ab_x + g_zz_0_yz_xxyyz[k];

                g_zz_0_xyz_xyzz[k] = -g_zz_0_yz_xyzz[k] * ab_x + g_zz_0_yz_xxyzz[k];

                g_zz_0_xyz_xzzz[k] = -g_zz_0_yz_xzzz[k] * ab_x + g_zz_0_yz_xxzzz[k];

                g_zz_0_xyz_yyyy[k] = -g_zz_0_yz_yyyy[k] * ab_x + g_zz_0_yz_xyyyy[k];

                g_zz_0_xyz_yyyz[k] = -g_zz_0_yz_yyyz[k] * ab_x + g_zz_0_yz_xyyyz[k];

                g_zz_0_xyz_yyzz[k] = -g_zz_0_yz_yyzz[k] * ab_x + g_zz_0_yz_xyyzz[k];

                g_zz_0_xyz_yzzz[k] = -g_zz_0_yz_yzzz[k] * ab_x + g_zz_0_yz_xyzzz[k];

                g_zz_0_xyz_zzzz[k] = -g_zz_0_yz_zzzz[k] * ab_x + g_zz_0_yz_xzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzz_xxxx = cbuffer.data(fg_geom_20_off + 825 * ccomps * dcomps);

            auto g_zz_0_xzz_xxxy = cbuffer.data(fg_geom_20_off + 826 * ccomps * dcomps);

            auto g_zz_0_xzz_xxxz = cbuffer.data(fg_geom_20_off + 827 * ccomps * dcomps);

            auto g_zz_0_xzz_xxyy = cbuffer.data(fg_geom_20_off + 828 * ccomps * dcomps);

            auto g_zz_0_xzz_xxyz = cbuffer.data(fg_geom_20_off + 829 * ccomps * dcomps);

            auto g_zz_0_xzz_xxzz = cbuffer.data(fg_geom_20_off + 830 * ccomps * dcomps);

            auto g_zz_0_xzz_xyyy = cbuffer.data(fg_geom_20_off + 831 * ccomps * dcomps);

            auto g_zz_0_xzz_xyyz = cbuffer.data(fg_geom_20_off + 832 * ccomps * dcomps);

            auto g_zz_0_xzz_xyzz = cbuffer.data(fg_geom_20_off + 833 * ccomps * dcomps);

            auto g_zz_0_xzz_xzzz = cbuffer.data(fg_geom_20_off + 834 * ccomps * dcomps);

            auto g_zz_0_xzz_yyyy = cbuffer.data(fg_geom_20_off + 835 * ccomps * dcomps);

            auto g_zz_0_xzz_yyyz = cbuffer.data(fg_geom_20_off + 836 * ccomps * dcomps);

            auto g_zz_0_xzz_yyzz = cbuffer.data(fg_geom_20_off + 837 * ccomps * dcomps);

            auto g_zz_0_xzz_yzzz = cbuffer.data(fg_geom_20_off + 838 * ccomps * dcomps);

            auto g_zz_0_xzz_zzzz = cbuffer.data(fg_geom_20_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzz_xxxx, g_zz_0_xzz_xxxy, g_zz_0_xzz_xxxz, g_zz_0_xzz_xxyy, g_zz_0_xzz_xxyz, g_zz_0_xzz_xxzz, g_zz_0_xzz_xyyy, g_zz_0_xzz_xyyz, g_zz_0_xzz_xyzz, g_zz_0_xzz_xzzz, g_zz_0_xzz_yyyy, g_zz_0_xzz_yyyz, g_zz_0_xzz_yyzz, g_zz_0_xzz_yzzz, g_zz_0_xzz_zzzz, g_zz_0_zz_xxxx, g_zz_0_zz_xxxxx, g_zz_0_zz_xxxxy, g_zz_0_zz_xxxxz, g_zz_0_zz_xxxy, g_zz_0_zz_xxxyy, g_zz_0_zz_xxxyz, g_zz_0_zz_xxxz, g_zz_0_zz_xxxzz, g_zz_0_zz_xxyy, g_zz_0_zz_xxyyy, g_zz_0_zz_xxyyz, g_zz_0_zz_xxyz, g_zz_0_zz_xxyzz, g_zz_0_zz_xxzz, g_zz_0_zz_xxzzz, g_zz_0_zz_xyyy, g_zz_0_zz_xyyyy, g_zz_0_zz_xyyyz, g_zz_0_zz_xyyz, g_zz_0_zz_xyyzz, g_zz_0_zz_xyzz, g_zz_0_zz_xyzzz, g_zz_0_zz_xzzz, g_zz_0_zz_xzzzz, g_zz_0_zz_yyyy, g_zz_0_zz_yyyz, g_zz_0_zz_yyzz, g_zz_0_zz_yzzz, g_zz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzz_xxxx[k] = -g_zz_0_zz_xxxx[k] * ab_x + g_zz_0_zz_xxxxx[k];

                g_zz_0_xzz_xxxy[k] = -g_zz_0_zz_xxxy[k] * ab_x + g_zz_0_zz_xxxxy[k];

                g_zz_0_xzz_xxxz[k] = -g_zz_0_zz_xxxz[k] * ab_x + g_zz_0_zz_xxxxz[k];

                g_zz_0_xzz_xxyy[k] = -g_zz_0_zz_xxyy[k] * ab_x + g_zz_0_zz_xxxyy[k];

                g_zz_0_xzz_xxyz[k] = -g_zz_0_zz_xxyz[k] * ab_x + g_zz_0_zz_xxxyz[k];

                g_zz_0_xzz_xxzz[k] = -g_zz_0_zz_xxzz[k] * ab_x + g_zz_0_zz_xxxzz[k];

                g_zz_0_xzz_xyyy[k] = -g_zz_0_zz_xyyy[k] * ab_x + g_zz_0_zz_xxyyy[k];

                g_zz_0_xzz_xyyz[k] = -g_zz_0_zz_xyyz[k] * ab_x + g_zz_0_zz_xxyyz[k];

                g_zz_0_xzz_xyzz[k] = -g_zz_0_zz_xyzz[k] * ab_x + g_zz_0_zz_xxyzz[k];

                g_zz_0_xzz_xzzz[k] = -g_zz_0_zz_xzzz[k] * ab_x + g_zz_0_zz_xxzzz[k];

                g_zz_0_xzz_yyyy[k] = -g_zz_0_zz_yyyy[k] * ab_x + g_zz_0_zz_xyyyy[k];

                g_zz_0_xzz_yyyz[k] = -g_zz_0_zz_yyyz[k] * ab_x + g_zz_0_zz_xyyyz[k];

                g_zz_0_xzz_yyzz[k] = -g_zz_0_zz_yyzz[k] * ab_x + g_zz_0_zz_xyyzz[k];

                g_zz_0_xzz_yzzz[k] = -g_zz_0_zz_yzzz[k] * ab_x + g_zz_0_zz_xyzzz[k];

                g_zz_0_xzz_zzzz[k] = -g_zz_0_zz_zzzz[k] * ab_x + g_zz_0_zz_xzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyy_xxxx = cbuffer.data(fg_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_yyy_xxxy = cbuffer.data(fg_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_yyy_xxxz = cbuffer.data(fg_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_yyy_xxyy = cbuffer.data(fg_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_yyy_xxyz = cbuffer.data(fg_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_yyy_xxzz = cbuffer.data(fg_geom_20_off + 845 * ccomps * dcomps);

            auto g_zz_0_yyy_xyyy = cbuffer.data(fg_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_yyy_xyyz = cbuffer.data(fg_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_yyy_xyzz = cbuffer.data(fg_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_yyy_xzzz = cbuffer.data(fg_geom_20_off + 849 * ccomps * dcomps);

            auto g_zz_0_yyy_yyyy = cbuffer.data(fg_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_yyy_yyyz = cbuffer.data(fg_geom_20_off + 851 * ccomps * dcomps);

            auto g_zz_0_yyy_yyzz = cbuffer.data(fg_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_yyy_yzzz = cbuffer.data(fg_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_yyy_zzzz = cbuffer.data(fg_geom_20_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yy_xxxx, g_zz_0_yy_xxxxy, g_zz_0_yy_xxxy, g_zz_0_yy_xxxyy, g_zz_0_yy_xxxyz, g_zz_0_yy_xxxz, g_zz_0_yy_xxyy, g_zz_0_yy_xxyyy, g_zz_0_yy_xxyyz, g_zz_0_yy_xxyz, g_zz_0_yy_xxyzz, g_zz_0_yy_xxzz, g_zz_0_yy_xyyy, g_zz_0_yy_xyyyy, g_zz_0_yy_xyyyz, g_zz_0_yy_xyyz, g_zz_0_yy_xyyzz, g_zz_0_yy_xyzz, g_zz_0_yy_xyzzz, g_zz_0_yy_xzzz, g_zz_0_yy_yyyy, g_zz_0_yy_yyyyy, g_zz_0_yy_yyyyz, g_zz_0_yy_yyyz, g_zz_0_yy_yyyzz, g_zz_0_yy_yyzz, g_zz_0_yy_yyzzz, g_zz_0_yy_yzzz, g_zz_0_yy_yzzzz, g_zz_0_yy_zzzz, g_zz_0_yyy_xxxx, g_zz_0_yyy_xxxy, g_zz_0_yyy_xxxz, g_zz_0_yyy_xxyy, g_zz_0_yyy_xxyz, g_zz_0_yyy_xxzz, g_zz_0_yyy_xyyy, g_zz_0_yyy_xyyz, g_zz_0_yyy_xyzz, g_zz_0_yyy_xzzz, g_zz_0_yyy_yyyy, g_zz_0_yyy_yyyz, g_zz_0_yyy_yyzz, g_zz_0_yyy_yzzz, g_zz_0_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyy_xxxx[k] = -g_zz_0_yy_xxxx[k] * ab_y + g_zz_0_yy_xxxxy[k];

                g_zz_0_yyy_xxxy[k] = -g_zz_0_yy_xxxy[k] * ab_y + g_zz_0_yy_xxxyy[k];

                g_zz_0_yyy_xxxz[k] = -g_zz_0_yy_xxxz[k] * ab_y + g_zz_0_yy_xxxyz[k];

                g_zz_0_yyy_xxyy[k] = -g_zz_0_yy_xxyy[k] * ab_y + g_zz_0_yy_xxyyy[k];

                g_zz_0_yyy_xxyz[k] = -g_zz_0_yy_xxyz[k] * ab_y + g_zz_0_yy_xxyyz[k];

                g_zz_0_yyy_xxzz[k] = -g_zz_0_yy_xxzz[k] * ab_y + g_zz_0_yy_xxyzz[k];

                g_zz_0_yyy_xyyy[k] = -g_zz_0_yy_xyyy[k] * ab_y + g_zz_0_yy_xyyyy[k];

                g_zz_0_yyy_xyyz[k] = -g_zz_0_yy_xyyz[k] * ab_y + g_zz_0_yy_xyyyz[k];

                g_zz_0_yyy_xyzz[k] = -g_zz_0_yy_xyzz[k] * ab_y + g_zz_0_yy_xyyzz[k];

                g_zz_0_yyy_xzzz[k] = -g_zz_0_yy_xzzz[k] * ab_y + g_zz_0_yy_xyzzz[k];

                g_zz_0_yyy_yyyy[k] = -g_zz_0_yy_yyyy[k] * ab_y + g_zz_0_yy_yyyyy[k];

                g_zz_0_yyy_yyyz[k] = -g_zz_0_yy_yyyz[k] * ab_y + g_zz_0_yy_yyyyz[k];

                g_zz_0_yyy_yyzz[k] = -g_zz_0_yy_yyzz[k] * ab_y + g_zz_0_yy_yyyzz[k];

                g_zz_0_yyy_yzzz[k] = -g_zz_0_yy_yzzz[k] * ab_y + g_zz_0_yy_yyzzz[k];

                g_zz_0_yyy_zzzz[k] = -g_zz_0_yy_zzzz[k] * ab_y + g_zz_0_yy_yzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyz_xxxx = cbuffer.data(fg_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_yyz_xxxy = cbuffer.data(fg_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_yyz_xxxz = cbuffer.data(fg_geom_20_off + 857 * ccomps * dcomps);

            auto g_zz_0_yyz_xxyy = cbuffer.data(fg_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_yyz_xxyz = cbuffer.data(fg_geom_20_off + 859 * ccomps * dcomps);

            auto g_zz_0_yyz_xxzz = cbuffer.data(fg_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_yyz_xyyy = cbuffer.data(fg_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_yyz_xyyz = cbuffer.data(fg_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_yyz_xyzz = cbuffer.data(fg_geom_20_off + 863 * ccomps * dcomps);

            auto g_zz_0_yyz_xzzz = cbuffer.data(fg_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_yyz_yyyy = cbuffer.data(fg_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_yyz_yyyz = cbuffer.data(fg_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_yyz_yyzz = cbuffer.data(fg_geom_20_off + 867 * ccomps * dcomps);

            auto g_zz_0_yyz_yzzz = cbuffer.data(fg_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_yyz_zzzz = cbuffer.data(fg_geom_20_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyz_xxxx, g_zz_0_yyz_xxxy, g_zz_0_yyz_xxxz, g_zz_0_yyz_xxyy, g_zz_0_yyz_xxyz, g_zz_0_yyz_xxzz, g_zz_0_yyz_xyyy, g_zz_0_yyz_xyyz, g_zz_0_yyz_xyzz, g_zz_0_yyz_xzzz, g_zz_0_yyz_yyyy, g_zz_0_yyz_yyyz, g_zz_0_yyz_yyzz, g_zz_0_yyz_yzzz, g_zz_0_yyz_zzzz, g_zz_0_yz_xxxx, g_zz_0_yz_xxxxy, g_zz_0_yz_xxxy, g_zz_0_yz_xxxyy, g_zz_0_yz_xxxyz, g_zz_0_yz_xxxz, g_zz_0_yz_xxyy, g_zz_0_yz_xxyyy, g_zz_0_yz_xxyyz, g_zz_0_yz_xxyz, g_zz_0_yz_xxyzz, g_zz_0_yz_xxzz, g_zz_0_yz_xyyy, g_zz_0_yz_xyyyy, g_zz_0_yz_xyyyz, g_zz_0_yz_xyyz, g_zz_0_yz_xyyzz, g_zz_0_yz_xyzz, g_zz_0_yz_xyzzz, g_zz_0_yz_xzzz, g_zz_0_yz_yyyy, g_zz_0_yz_yyyyy, g_zz_0_yz_yyyyz, g_zz_0_yz_yyyz, g_zz_0_yz_yyyzz, g_zz_0_yz_yyzz, g_zz_0_yz_yyzzz, g_zz_0_yz_yzzz, g_zz_0_yz_yzzzz, g_zz_0_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyz_xxxx[k] = -g_zz_0_yz_xxxx[k] * ab_y + g_zz_0_yz_xxxxy[k];

                g_zz_0_yyz_xxxy[k] = -g_zz_0_yz_xxxy[k] * ab_y + g_zz_0_yz_xxxyy[k];

                g_zz_0_yyz_xxxz[k] = -g_zz_0_yz_xxxz[k] * ab_y + g_zz_0_yz_xxxyz[k];

                g_zz_0_yyz_xxyy[k] = -g_zz_0_yz_xxyy[k] * ab_y + g_zz_0_yz_xxyyy[k];

                g_zz_0_yyz_xxyz[k] = -g_zz_0_yz_xxyz[k] * ab_y + g_zz_0_yz_xxyyz[k];

                g_zz_0_yyz_xxzz[k] = -g_zz_0_yz_xxzz[k] * ab_y + g_zz_0_yz_xxyzz[k];

                g_zz_0_yyz_xyyy[k] = -g_zz_0_yz_xyyy[k] * ab_y + g_zz_0_yz_xyyyy[k];

                g_zz_0_yyz_xyyz[k] = -g_zz_0_yz_xyyz[k] * ab_y + g_zz_0_yz_xyyyz[k];

                g_zz_0_yyz_xyzz[k] = -g_zz_0_yz_xyzz[k] * ab_y + g_zz_0_yz_xyyzz[k];

                g_zz_0_yyz_xzzz[k] = -g_zz_0_yz_xzzz[k] * ab_y + g_zz_0_yz_xyzzz[k];

                g_zz_0_yyz_yyyy[k] = -g_zz_0_yz_yyyy[k] * ab_y + g_zz_0_yz_yyyyy[k];

                g_zz_0_yyz_yyyz[k] = -g_zz_0_yz_yyyz[k] * ab_y + g_zz_0_yz_yyyyz[k];

                g_zz_0_yyz_yyzz[k] = -g_zz_0_yz_yyzz[k] * ab_y + g_zz_0_yz_yyyzz[k];

                g_zz_0_yyz_yzzz[k] = -g_zz_0_yz_yzzz[k] * ab_y + g_zz_0_yz_yyzzz[k];

                g_zz_0_yyz_zzzz[k] = -g_zz_0_yz_zzzz[k] * ab_y + g_zz_0_yz_yzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzz_xxxx = cbuffer.data(fg_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_yzz_xxxy = cbuffer.data(fg_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_yzz_xxxz = cbuffer.data(fg_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_yzz_xxyy = cbuffer.data(fg_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_yzz_xxyz = cbuffer.data(fg_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_yzz_xxzz = cbuffer.data(fg_geom_20_off + 875 * ccomps * dcomps);

            auto g_zz_0_yzz_xyyy = cbuffer.data(fg_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_yzz_xyyz = cbuffer.data(fg_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_yzz_xyzz = cbuffer.data(fg_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_yzz_xzzz = cbuffer.data(fg_geom_20_off + 879 * ccomps * dcomps);

            auto g_zz_0_yzz_yyyy = cbuffer.data(fg_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_yzz_yyyz = cbuffer.data(fg_geom_20_off + 881 * ccomps * dcomps);

            auto g_zz_0_yzz_yyzz = cbuffer.data(fg_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_yzz_yzzz = cbuffer.data(fg_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_yzz_zzzz = cbuffer.data(fg_geom_20_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzz_xxxx, g_zz_0_yzz_xxxy, g_zz_0_yzz_xxxz, g_zz_0_yzz_xxyy, g_zz_0_yzz_xxyz, g_zz_0_yzz_xxzz, g_zz_0_yzz_xyyy, g_zz_0_yzz_xyyz, g_zz_0_yzz_xyzz, g_zz_0_yzz_xzzz, g_zz_0_yzz_yyyy, g_zz_0_yzz_yyyz, g_zz_0_yzz_yyzz, g_zz_0_yzz_yzzz, g_zz_0_yzz_zzzz, g_zz_0_zz_xxxx, g_zz_0_zz_xxxxy, g_zz_0_zz_xxxy, g_zz_0_zz_xxxyy, g_zz_0_zz_xxxyz, g_zz_0_zz_xxxz, g_zz_0_zz_xxyy, g_zz_0_zz_xxyyy, g_zz_0_zz_xxyyz, g_zz_0_zz_xxyz, g_zz_0_zz_xxyzz, g_zz_0_zz_xxzz, g_zz_0_zz_xyyy, g_zz_0_zz_xyyyy, g_zz_0_zz_xyyyz, g_zz_0_zz_xyyz, g_zz_0_zz_xyyzz, g_zz_0_zz_xyzz, g_zz_0_zz_xyzzz, g_zz_0_zz_xzzz, g_zz_0_zz_yyyy, g_zz_0_zz_yyyyy, g_zz_0_zz_yyyyz, g_zz_0_zz_yyyz, g_zz_0_zz_yyyzz, g_zz_0_zz_yyzz, g_zz_0_zz_yyzzz, g_zz_0_zz_yzzz, g_zz_0_zz_yzzzz, g_zz_0_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzz_xxxx[k] = -g_zz_0_zz_xxxx[k] * ab_y + g_zz_0_zz_xxxxy[k];

                g_zz_0_yzz_xxxy[k] = -g_zz_0_zz_xxxy[k] * ab_y + g_zz_0_zz_xxxyy[k];

                g_zz_0_yzz_xxxz[k] = -g_zz_0_zz_xxxz[k] * ab_y + g_zz_0_zz_xxxyz[k];

                g_zz_0_yzz_xxyy[k] = -g_zz_0_zz_xxyy[k] * ab_y + g_zz_0_zz_xxyyy[k];

                g_zz_0_yzz_xxyz[k] = -g_zz_0_zz_xxyz[k] * ab_y + g_zz_0_zz_xxyyz[k];

                g_zz_0_yzz_xxzz[k] = -g_zz_0_zz_xxzz[k] * ab_y + g_zz_0_zz_xxyzz[k];

                g_zz_0_yzz_xyyy[k] = -g_zz_0_zz_xyyy[k] * ab_y + g_zz_0_zz_xyyyy[k];

                g_zz_0_yzz_xyyz[k] = -g_zz_0_zz_xyyz[k] * ab_y + g_zz_0_zz_xyyyz[k];

                g_zz_0_yzz_xyzz[k] = -g_zz_0_zz_xyzz[k] * ab_y + g_zz_0_zz_xyyzz[k];

                g_zz_0_yzz_xzzz[k] = -g_zz_0_zz_xzzz[k] * ab_y + g_zz_0_zz_xyzzz[k];

                g_zz_0_yzz_yyyy[k] = -g_zz_0_zz_yyyy[k] * ab_y + g_zz_0_zz_yyyyy[k];

                g_zz_0_yzz_yyyz[k] = -g_zz_0_zz_yyyz[k] * ab_y + g_zz_0_zz_yyyyz[k];

                g_zz_0_yzz_yyzz[k] = -g_zz_0_zz_yyzz[k] * ab_y + g_zz_0_zz_yyyzz[k];

                g_zz_0_yzz_yzzz[k] = -g_zz_0_zz_yzzz[k] * ab_y + g_zz_0_zz_yyzzz[k];

                g_zz_0_yzz_zzzz[k] = -g_zz_0_zz_zzzz[k] * ab_y + g_zz_0_zz_yzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzz_xxxx = cbuffer.data(fg_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_zzz_xxxy = cbuffer.data(fg_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_zzz_xxxz = cbuffer.data(fg_geom_20_off + 887 * ccomps * dcomps);

            auto g_zz_0_zzz_xxyy = cbuffer.data(fg_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_zzz_xxyz = cbuffer.data(fg_geom_20_off + 889 * ccomps * dcomps);

            auto g_zz_0_zzz_xxzz = cbuffer.data(fg_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_zzz_xyyy = cbuffer.data(fg_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_zzz_xyyz = cbuffer.data(fg_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_zzz_xyzz = cbuffer.data(fg_geom_20_off + 893 * ccomps * dcomps);

            auto g_zz_0_zzz_xzzz = cbuffer.data(fg_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_zzz_yyyy = cbuffer.data(fg_geom_20_off + 895 * ccomps * dcomps);

            auto g_zz_0_zzz_yyyz = cbuffer.data(fg_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_zzz_yyzz = cbuffer.data(fg_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_zzz_yzzz = cbuffer.data(fg_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_zzz_zzzz = cbuffer.data(fg_geom_20_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_zz_0_zz_xxxx, g_zz_0_zz_xxxxz, g_zz_0_zz_xxxy, g_zz_0_zz_xxxyz, g_zz_0_zz_xxxz, g_zz_0_zz_xxxzz, g_zz_0_zz_xxyy, g_zz_0_zz_xxyyz, g_zz_0_zz_xxyz, g_zz_0_zz_xxyzz, g_zz_0_zz_xxzz, g_zz_0_zz_xxzzz, g_zz_0_zz_xyyy, g_zz_0_zz_xyyyz, g_zz_0_zz_xyyz, g_zz_0_zz_xyyzz, g_zz_0_zz_xyzz, g_zz_0_zz_xyzzz, g_zz_0_zz_xzzz, g_zz_0_zz_xzzzz, g_zz_0_zz_yyyy, g_zz_0_zz_yyyyz, g_zz_0_zz_yyyz, g_zz_0_zz_yyyzz, g_zz_0_zz_yyzz, g_zz_0_zz_yyzzz, g_zz_0_zz_yzzz, g_zz_0_zz_yzzzz, g_zz_0_zz_zzzz, g_zz_0_zz_zzzzz, g_zz_0_zzz_xxxx, g_zz_0_zzz_xxxy, g_zz_0_zzz_xxxz, g_zz_0_zzz_xxyy, g_zz_0_zzz_xxyz, g_zz_0_zzz_xxzz, g_zz_0_zzz_xyyy, g_zz_0_zzz_xyyz, g_zz_0_zzz_xyzz, g_zz_0_zzz_xzzz, g_zz_0_zzz_yyyy, g_zz_0_zzz_yyyz, g_zz_0_zzz_yyzz, g_zz_0_zzz_yzzz, g_zz_0_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzz_xxxx[k] = -2.0 * g_z_0_zz_xxxx[k] - g_zz_0_zz_xxxx[k] * ab_z + g_zz_0_zz_xxxxz[k];

                g_zz_0_zzz_xxxy[k] = -2.0 * g_z_0_zz_xxxy[k] - g_zz_0_zz_xxxy[k] * ab_z + g_zz_0_zz_xxxyz[k];

                g_zz_0_zzz_xxxz[k] = -2.0 * g_z_0_zz_xxxz[k] - g_zz_0_zz_xxxz[k] * ab_z + g_zz_0_zz_xxxzz[k];

                g_zz_0_zzz_xxyy[k] = -2.0 * g_z_0_zz_xxyy[k] - g_zz_0_zz_xxyy[k] * ab_z + g_zz_0_zz_xxyyz[k];

                g_zz_0_zzz_xxyz[k] = -2.0 * g_z_0_zz_xxyz[k] - g_zz_0_zz_xxyz[k] * ab_z + g_zz_0_zz_xxyzz[k];

                g_zz_0_zzz_xxzz[k] = -2.0 * g_z_0_zz_xxzz[k] - g_zz_0_zz_xxzz[k] * ab_z + g_zz_0_zz_xxzzz[k];

                g_zz_0_zzz_xyyy[k] = -2.0 * g_z_0_zz_xyyy[k] - g_zz_0_zz_xyyy[k] * ab_z + g_zz_0_zz_xyyyz[k];

                g_zz_0_zzz_xyyz[k] = -2.0 * g_z_0_zz_xyyz[k] - g_zz_0_zz_xyyz[k] * ab_z + g_zz_0_zz_xyyzz[k];

                g_zz_0_zzz_xyzz[k] = -2.0 * g_z_0_zz_xyzz[k] - g_zz_0_zz_xyzz[k] * ab_z + g_zz_0_zz_xyzzz[k];

                g_zz_0_zzz_xzzz[k] = -2.0 * g_z_0_zz_xzzz[k] - g_zz_0_zz_xzzz[k] * ab_z + g_zz_0_zz_xzzzz[k];

                g_zz_0_zzz_yyyy[k] = -2.0 * g_z_0_zz_yyyy[k] - g_zz_0_zz_yyyy[k] * ab_z + g_zz_0_zz_yyyyz[k];

                g_zz_0_zzz_yyyz[k] = -2.0 * g_z_0_zz_yyyz[k] - g_zz_0_zz_yyyz[k] * ab_z + g_zz_0_zz_yyyzz[k];

                g_zz_0_zzz_yyzz[k] = -2.0 * g_z_0_zz_yyzz[k] - g_zz_0_zz_yyzz[k] * ab_z + g_zz_0_zz_yyzzz[k];

                g_zz_0_zzz_yzzz[k] = -2.0 * g_z_0_zz_yzzz[k] - g_zz_0_zz_yzzz[k] * ab_z + g_zz_0_zz_yzzzz[k];

                g_zz_0_zzz_zzzz[k] = -2.0 * g_z_0_zz_zzzz[k] - g_zz_0_zz_zzzz[k] * ab_z + g_zz_0_zz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

