#include "ElectronRepulsionGeom1100ContrRecFGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_fgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_fgxx,
                                            const size_t idx_geom_01_dgxx,
                                            const size_t idx_geom_10_dgxx,
                                            const size_t idx_geom_11_dgxx,
                                            const size_t idx_geom_11_dhxx,
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

            const auto dg_geom_01_off = idx_geom_01_dgxx + i * dcomps + j;

            auto g_0_x_xx_xxxx = cbuffer.data(dg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxy = cbuffer.data(dg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxz = cbuffer.data(dg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxyy = cbuffer.data(dg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxyz = cbuffer.data(dg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxzz = cbuffer.data(dg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xyyy = cbuffer.data(dg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xyyz = cbuffer.data(dg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xyzz = cbuffer.data(dg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xzzz = cbuffer.data(dg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_yyyy = cbuffer.data(dg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_yyyz = cbuffer.data(dg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_yyzz = cbuffer.data(dg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_yzzz = cbuffer.data(dg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_zzzz = cbuffer.data(dg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xy_xxxx = cbuffer.data(dg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xy_xxxy = cbuffer.data(dg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xy_xxxz = cbuffer.data(dg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xy_xxyy = cbuffer.data(dg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xy_xxyz = cbuffer.data(dg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xy_xxzz = cbuffer.data(dg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xy_xyyy = cbuffer.data(dg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xy_xyyz = cbuffer.data(dg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xy_xyzz = cbuffer.data(dg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xy_xzzz = cbuffer.data(dg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xy_yyyy = cbuffer.data(dg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xy_yyyz = cbuffer.data(dg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xy_yyzz = cbuffer.data(dg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xy_yzzz = cbuffer.data(dg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xy_zzzz = cbuffer.data(dg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xz_xxxx = cbuffer.data(dg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xz_xxxy = cbuffer.data(dg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xz_xxxz = cbuffer.data(dg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xz_xxyy = cbuffer.data(dg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xz_xxyz = cbuffer.data(dg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xz_xxzz = cbuffer.data(dg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xz_xyyy = cbuffer.data(dg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xz_xyyz = cbuffer.data(dg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xz_xyzz = cbuffer.data(dg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xz_xzzz = cbuffer.data(dg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xz_yyyy = cbuffer.data(dg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xz_yyyz = cbuffer.data(dg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xz_yyzz = cbuffer.data(dg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xz_yzzz = cbuffer.data(dg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xz_zzzz = cbuffer.data(dg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yy_xxxx = cbuffer.data(dg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yy_xxxy = cbuffer.data(dg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yy_xxxz = cbuffer.data(dg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yy_xxyy = cbuffer.data(dg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yy_xxyz = cbuffer.data(dg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yy_xxzz = cbuffer.data(dg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_yy_xyyy = cbuffer.data(dg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yy_xyyz = cbuffer.data(dg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yy_xyzz = cbuffer.data(dg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_yy_xzzz = cbuffer.data(dg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_yy_yyyy = cbuffer.data(dg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_yy_yyyz = cbuffer.data(dg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_yy_yyzz = cbuffer.data(dg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_yy_yzzz = cbuffer.data(dg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_yy_zzzz = cbuffer.data(dg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yz_xxxx = cbuffer.data(dg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yz_xxxy = cbuffer.data(dg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yz_xxxz = cbuffer.data(dg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yz_xxyy = cbuffer.data(dg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yz_xxyz = cbuffer.data(dg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yz_xxzz = cbuffer.data(dg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yz_xyyy = cbuffer.data(dg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yz_xyyz = cbuffer.data(dg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yz_xyzz = cbuffer.data(dg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yz_xzzz = cbuffer.data(dg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yz_yyyy = cbuffer.data(dg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yz_yyyz = cbuffer.data(dg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yz_yyzz = cbuffer.data(dg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yz_yzzz = cbuffer.data(dg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yz_zzzz = cbuffer.data(dg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_zz_xxxx = cbuffer.data(dg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_zz_xxxy = cbuffer.data(dg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_zz_xxxz = cbuffer.data(dg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_zz_xxyy = cbuffer.data(dg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_zz_xxyz = cbuffer.data(dg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_zz_xxzz = cbuffer.data(dg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_zz_xyyy = cbuffer.data(dg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_zz_xyyz = cbuffer.data(dg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_zz_xyzz = cbuffer.data(dg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_zz_xzzz = cbuffer.data(dg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zz_yyyy = cbuffer.data(dg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zz_yyyz = cbuffer.data(dg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zz_yyzz = cbuffer.data(dg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zz_yzzz = cbuffer.data(dg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zz_zzzz = cbuffer.data(dg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xx_xxxx = cbuffer.data(dg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xx_xxxy = cbuffer.data(dg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xx_xxxz = cbuffer.data(dg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xx_xxyy = cbuffer.data(dg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xx_xxyz = cbuffer.data(dg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xx_xxzz = cbuffer.data(dg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xx_xyyy = cbuffer.data(dg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xx_xyyz = cbuffer.data(dg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xx_xyzz = cbuffer.data(dg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xx_xzzz = cbuffer.data(dg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xx_yyyy = cbuffer.data(dg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xx_yyyz = cbuffer.data(dg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xx_yyzz = cbuffer.data(dg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xx_yzzz = cbuffer.data(dg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xx_zzzz = cbuffer.data(dg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xy_xxxx = cbuffer.data(dg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xy_xxxy = cbuffer.data(dg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xy_xxxz = cbuffer.data(dg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xy_xxyy = cbuffer.data(dg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xy_xxyz = cbuffer.data(dg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xy_xxzz = cbuffer.data(dg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xy_xyyy = cbuffer.data(dg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xy_xyyz = cbuffer.data(dg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xy_xyzz = cbuffer.data(dg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xy_xzzz = cbuffer.data(dg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xy_yyyy = cbuffer.data(dg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xy_yyyz = cbuffer.data(dg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xy_yyzz = cbuffer.data(dg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xy_yzzz = cbuffer.data(dg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xy_zzzz = cbuffer.data(dg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xz_xxxx = cbuffer.data(dg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xz_xxxy = cbuffer.data(dg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xz_xxxz = cbuffer.data(dg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xz_xxyy = cbuffer.data(dg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xz_xxyz = cbuffer.data(dg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xz_xxzz = cbuffer.data(dg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xz_xyyy = cbuffer.data(dg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xz_xyyz = cbuffer.data(dg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xz_xyzz = cbuffer.data(dg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xz_xzzz = cbuffer.data(dg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xz_yyyy = cbuffer.data(dg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xz_yyyz = cbuffer.data(dg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xz_yyzz = cbuffer.data(dg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xz_yzzz = cbuffer.data(dg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xz_zzzz = cbuffer.data(dg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_yy_xxxx = cbuffer.data(dg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_yy_xxxy = cbuffer.data(dg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_yy_xxxz = cbuffer.data(dg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_yy_xxyy = cbuffer.data(dg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_yy_xxyz = cbuffer.data(dg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_yy_xxzz = cbuffer.data(dg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_yy_xyyy = cbuffer.data(dg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_yy_xyyz = cbuffer.data(dg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_yy_xyzz = cbuffer.data(dg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_yy_xzzz = cbuffer.data(dg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_yy_yyyy = cbuffer.data(dg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_yy_yyyz = cbuffer.data(dg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_yy_yyzz = cbuffer.data(dg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_yy_yzzz = cbuffer.data(dg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_yy_zzzz = cbuffer.data(dg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yz_xxxx = cbuffer.data(dg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yz_xxxy = cbuffer.data(dg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yz_xxxz = cbuffer.data(dg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yz_xxyy = cbuffer.data(dg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yz_xxyz = cbuffer.data(dg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yz_xxzz = cbuffer.data(dg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yz_xyyy = cbuffer.data(dg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yz_xyyz = cbuffer.data(dg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yz_xyzz = cbuffer.data(dg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yz_xzzz = cbuffer.data(dg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yz_yyyy = cbuffer.data(dg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yz_yyyz = cbuffer.data(dg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yz_yyzz = cbuffer.data(dg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yz_yzzz = cbuffer.data(dg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yz_zzzz = cbuffer.data(dg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_zz_xxxx = cbuffer.data(dg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_zz_xxxy = cbuffer.data(dg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_zz_xxxz = cbuffer.data(dg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_zz_xxyy = cbuffer.data(dg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_zz_xxyz = cbuffer.data(dg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_zz_xxzz = cbuffer.data(dg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_zz_xyyy = cbuffer.data(dg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_zz_xyyz = cbuffer.data(dg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_zz_xyzz = cbuffer.data(dg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_zz_xzzz = cbuffer.data(dg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zz_yyyy = cbuffer.data(dg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zz_yyyz = cbuffer.data(dg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zz_yyzz = cbuffer.data(dg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zz_yzzz = cbuffer.data(dg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zz_zzzz = cbuffer.data(dg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xx_xxxx = cbuffer.data(dg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xx_xxxy = cbuffer.data(dg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xx_xxxz = cbuffer.data(dg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xx_xxyy = cbuffer.data(dg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xx_xxyz = cbuffer.data(dg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xx_xxzz = cbuffer.data(dg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xx_xyyy = cbuffer.data(dg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xx_xyyz = cbuffer.data(dg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xx_xyzz = cbuffer.data(dg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xx_xzzz = cbuffer.data(dg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xx_yyyy = cbuffer.data(dg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xx_yyyz = cbuffer.data(dg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xx_yyzz = cbuffer.data(dg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xx_yzzz = cbuffer.data(dg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xx_zzzz = cbuffer.data(dg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xy_xxxx = cbuffer.data(dg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xy_xxxy = cbuffer.data(dg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xy_xxxz = cbuffer.data(dg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xy_xxyy = cbuffer.data(dg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xy_xxyz = cbuffer.data(dg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xy_xxzz = cbuffer.data(dg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xy_xyyy = cbuffer.data(dg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xy_xyyz = cbuffer.data(dg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xy_xyzz = cbuffer.data(dg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xy_xzzz = cbuffer.data(dg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xy_yyyy = cbuffer.data(dg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xy_yyyz = cbuffer.data(dg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xy_yyzz = cbuffer.data(dg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xy_yzzz = cbuffer.data(dg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xy_zzzz = cbuffer.data(dg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xz_xxxx = cbuffer.data(dg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xz_xxxy = cbuffer.data(dg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xz_xxxz = cbuffer.data(dg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xz_xxyy = cbuffer.data(dg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xz_xxyz = cbuffer.data(dg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xz_xxzz = cbuffer.data(dg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xz_xyyy = cbuffer.data(dg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xz_xyyz = cbuffer.data(dg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xz_xyzz = cbuffer.data(dg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xz_xzzz = cbuffer.data(dg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xz_yyyy = cbuffer.data(dg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xz_yyyz = cbuffer.data(dg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xz_yyzz = cbuffer.data(dg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xz_yzzz = cbuffer.data(dg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xz_zzzz = cbuffer.data(dg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_yy_xxxx = cbuffer.data(dg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_yy_xxxy = cbuffer.data(dg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_yy_xxxz = cbuffer.data(dg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_yy_xxyy = cbuffer.data(dg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_yy_xxyz = cbuffer.data(dg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_yy_xxzz = cbuffer.data(dg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_yy_xyyy = cbuffer.data(dg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_yy_xyyz = cbuffer.data(dg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_yy_xyzz = cbuffer.data(dg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_yy_xzzz = cbuffer.data(dg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_yy_yyyy = cbuffer.data(dg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_yy_yyyz = cbuffer.data(dg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_yy_yyzz = cbuffer.data(dg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_yy_yzzz = cbuffer.data(dg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_yy_zzzz = cbuffer.data(dg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yz_xxxx = cbuffer.data(dg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yz_xxxy = cbuffer.data(dg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yz_xxxz = cbuffer.data(dg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yz_xxyy = cbuffer.data(dg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yz_xxyz = cbuffer.data(dg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yz_xxzz = cbuffer.data(dg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yz_xyyy = cbuffer.data(dg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yz_xyyz = cbuffer.data(dg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yz_xyzz = cbuffer.data(dg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yz_xzzz = cbuffer.data(dg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yz_yyyy = cbuffer.data(dg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yz_yyyz = cbuffer.data(dg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_yz_yyzz = cbuffer.data(dg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yz_yzzz = cbuffer.data(dg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yz_zzzz = cbuffer.data(dg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_zz_xxxx = cbuffer.data(dg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_zz_xxxy = cbuffer.data(dg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_zz_xxxz = cbuffer.data(dg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_zz_xxyy = cbuffer.data(dg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_zz_xxyz = cbuffer.data(dg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_zz_xxzz = cbuffer.data(dg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_zz_xyyy = cbuffer.data(dg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_zz_xyyz = cbuffer.data(dg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_zz_xyzz = cbuffer.data(dg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_zz_xzzz = cbuffer.data(dg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zz_yyyy = cbuffer.data(dg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zz_yyyz = cbuffer.data(dg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zz_yyzz = cbuffer.data(dg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zz_yzzz = cbuffer.data(dg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zz_zzzz = cbuffer.data(dg_geom_01_off + 269 * ccomps * dcomps);

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

            const auto dg_geom_11_off = idx_geom_11_dgxx + i * dcomps + j;

            auto g_x_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 419 * ccomps * dcomps);

            auto g_y_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_x_xx_xxxx = cbuffer.data(dg_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_x_xx_xxxy = cbuffer.data(dg_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_x_xx_xxxz = cbuffer.data(dg_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_x_xx_xxyy = cbuffer.data(dg_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_x_xx_xxyz = cbuffer.data(dg_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_x_xx_xxzz = cbuffer.data(dg_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_x_xx_xyyy = cbuffer.data(dg_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_x_xx_xyyz = cbuffer.data(dg_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_x_xx_xyzz = cbuffer.data(dg_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_x_xx_xzzz = cbuffer.data(dg_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_x_xx_yyyy = cbuffer.data(dg_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_x_xx_yyyz = cbuffer.data(dg_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_x_xx_yyzz = cbuffer.data(dg_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_x_xx_yzzz = cbuffer.data(dg_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_x_xx_zzzz = cbuffer.data(dg_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_x_xy_xxxx = cbuffer.data(dg_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_x_xy_xxxy = cbuffer.data(dg_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_x_xy_xxxz = cbuffer.data(dg_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_x_xy_xxyy = cbuffer.data(dg_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_x_xy_xxyz = cbuffer.data(dg_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_x_xy_xxzz = cbuffer.data(dg_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_x_xy_xyyy = cbuffer.data(dg_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_x_xy_xyyz = cbuffer.data(dg_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_x_xy_xyzz = cbuffer.data(dg_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_x_xy_xzzz = cbuffer.data(dg_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_x_xy_yyyy = cbuffer.data(dg_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_x_xy_yyyz = cbuffer.data(dg_geom_11_off + 566 * ccomps * dcomps);

            auto g_z_x_xy_yyzz = cbuffer.data(dg_geom_11_off + 567 * ccomps * dcomps);

            auto g_z_x_xy_yzzz = cbuffer.data(dg_geom_11_off + 568 * ccomps * dcomps);

            auto g_z_x_xy_zzzz = cbuffer.data(dg_geom_11_off + 569 * ccomps * dcomps);

            auto g_z_x_xz_xxxx = cbuffer.data(dg_geom_11_off + 570 * ccomps * dcomps);

            auto g_z_x_xz_xxxy = cbuffer.data(dg_geom_11_off + 571 * ccomps * dcomps);

            auto g_z_x_xz_xxxz = cbuffer.data(dg_geom_11_off + 572 * ccomps * dcomps);

            auto g_z_x_xz_xxyy = cbuffer.data(dg_geom_11_off + 573 * ccomps * dcomps);

            auto g_z_x_xz_xxyz = cbuffer.data(dg_geom_11_off + 574 * ccomps * dcomps);

            auto g_z_x_xz_xxzz = cbuffer.data(dg_geom_11_off + 575 * ccomps * dcomps);

            auto g_z_x_xz_xyyy = cbuffer.data(dg_geom_11_off + 576 * ccomps * dcomps);

            auto g_z_x_xz_xyyz = cbuffer.data(dg_geom_11_off + 577 * ccomps * dcomps);

            auto g_z_x_xz_xyzz = cbuffer.data(dg_geom_11_off + 578 * ccomps * dcomps);

            auto g_z_x_xz_xzzz = cbuffer.data(dg_geom_11_off + 579 * ccomps * dcomps);

            auto g_z_x_xz_yyyy = cbuffer.data(dg_geom_11_off + 580 * ccomps * dcomps);

            auto g_z_x_xz_yyyz = cbuffer.data(dg_geom_11_off + 581 * ccomps * dcomps);

            auto g_z_x_xz_yyzz = cbuffer.data(dg_geom_11_off + 582 * ccomps * dcomps);

            auto g_z_x_xz_yzzz = cbuffer.data(dg_geom_11_off + 583 * ccomps * dcomps);

            auto g_z_x_xz_zzzz = cbuffer.data(dg_geom_11_off + 584 * ccomps * dcomps);

            auto g_z_x_yy_xxxx = cbuffer.data(dg_geom_11_off + 585 * ccomps * dcomps);

            auto g_z_x_yy_xxxy = cbuffer.data(dg_geom_11_off + 586 * ccomps * dcomps);

            auto g_z_x_yy_xxxz = cbuffer.data(dg_geom_11_off + 587 * ccomps * dcomps);

            auto g_z_x_yy_xxyy = cbuffer.data(dg_geom_11_off + 588 * ccomps * dcomps);

            auto g_z_x_yy_xxyz = cbuffer.data(dg_geom_11_off + 589 * ccomps * dcomps);

            auto g_z_x_yy_xxzz = cbuffer.data(dg_geom_11_off + 590 * ccomps * dcomps);

            auto g_z_x_yy_xyyy = cbuffer.data(dg_geom_11_off + 591 * ccomps * dcomps);

            auto g_z_x_yy_xyyz = cbuffer.data(dg_geom_11_off + 592 * ccomps * dcomps);

            auto g_z_x_yy_xyzz = cbuffer.data(dg_geom_11_off + 593 * ccomps * dcomps);

            auto g_z_x_yy_xzzz = cbuffer.data(dg_geom_11_off + 594 * ccomps * dcomps);

            auto g_z_x_yy_yyyy = cbuffer.data(dg_geom_11_off + 595 * ccomps * dcomps);

            auto g_z_x_yy_yyyz = cbuffer.data(dg_geom_11_off + 596 * ccomps * dcomps);

            auto g_z_x_yy_yyzz = cbuffer.data(dg_geom_11_off + 597 * ccomps * dcomps);

            auto g_z_x_yy_yzzz = cbuffer.data(dg_geom_11_off + 598 * ccomps * dcomps);

            auto g_z_x_yy_zzzz = cbuffer.data(dg_geom_11_off + 599 * ccomps * dcomps);

            auto g_z_x_yz_xxxx = cbuffer.data(dg_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_x_yz_xxxy = cbuffer.data(dg_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_x_yz_xxxz = cbuffer.data(dg_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_x_yz_xxyy = cbuffer.data(dg_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_x_yz_xxyz = cbuffer.data(dg_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_x_yz_xxzz = cbuffer.data(dg_geom_11_off + 605 * ccomps * dcomps);

            auto g_z_x_yz_xyyy = cbuffer.data(dg_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_x_yz_xyyz = cbuffer.data(dg_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_x_yz_xyzz = cbuffer.data(dg_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_x_yz_xzzz = cbuffer.data(dg_geom_11_off + 609 * ccomps * dcomps);

            auto g_z_x_yz_yyyy = cbuffer.data(dg_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_x_yz_yyyz = cbuffer.data(dg_geom_11_off + 611 * ccomps * dcomps);

            auto g_z_x_yz_yyzz = cbuffer.data(dg_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_x_yz_yzzz = cbuffer.data(dg_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_x_yz_zzzz = cbuffer.data(dg_geom_11_off + 614 * ccomps * dcomps);

            auto g_z_x_zz_xxxx = cbuffer.data(dg_geom_11_off + 615 * ccomps * dcomps);

            auto g_z_x_zz_xxxy = cbuffer.data(dg_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_x_zz_xxxz = cbuffer.data(dg_geom_11_off + 617 * ccomps * dcomps);

            auto g_z_x_zz_xxyy = cbuffer.data(dg_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_x_zz_xxyz = cbuffer.data(dg_geom_11_off + 619 * ccomps * dcomps);

            auto g_z_x_zz_xxzz = cbuffer.data(dg_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_x_zz_xyyy = cbuffer.data(dg_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_x_zz_xyyz = cbuffer.data(dg_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_x_zz_xyzz = cbuffer.data(dg_geom_11_off + 623 * ccomps * dcomps);

            auto g_z_x_zz_xzzz = cbuffer.data(dg_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_x_zz_yyyy = cbuffer.data(dg_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_x_zz_yyyz = cbuffer.data(dg_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_x_zz_yyzz = cbuffer.data(dg_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_x_zz_yzzz = cbuffer.data(dg_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_x_zz_zzzz = cbuffer.data(dg_geom_11_off + 629 * ccomps * dcomps);

            auto g_z_y_xx_xxxx = cbuffer.data(dg_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_y_xx_xxxy = cbuffer.data(dg_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_y_xx_xxxz = cbuffer.data(dg_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_y_xx_xxyy = cbuffer.data(dg_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_y_xx_xxyz = cbuffer.data(dg_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_y_xx_xxzz = cbuffer.data(dg_geom_11_off + 635 * ccomps * dcomps);

            auto g_z_y_xx_xyyy = cbuffer.data(dg_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_y_xx_xyyz = cbuffer.data(dg_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_y_xx_xyzz = cbuffer.data(dg_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_y_xx_xzzz = cbuffer.data(dg_geom_11_off + 639 * ccomps * dcomps);

            auto g_z_y_xx_yyyy = cbuffer.data(dg_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_y_xx_yyyz = cbuffer.data(dg_geom_11_off + 641 * ccomps * dcomps);

            auto g_z_y_xx_yyzz = cbuffer.data(dg_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_y_xx_yzzz = cbuffer.data(dg_geom_11_off + 643 * ccomps * dcomps);

            auto g_z_y_xx_zzzz = cbuffer.data(dg_geom_11_off + 644 * ccomps * dcomps);

            auto g_z_y_xy_xxxx = cbuffer.data(dg_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_y_xy_xxxy = cbuffer.data(dg_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_y_xy_xxxz = cbuffer.data(dg_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_y_xy_xxyy = cbuffer.data(dg_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_y_xy_xxyz = cbuffer.data(dg_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_y_xy_xxzz = cbuffer.data(dg_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_y_xy_xyyy = cbuffer.data(dg_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_y_xy_xyyz = cbuffer.data(dg_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_y_xy_xyzz = cbuffer.data(dg_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_y_xy_xzzz = cbuffer.data(dg_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_y_xy_yyyy = cbuffer.data(dg_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_y_xy_yyyz = cbuffer.data(dg_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_y_xy_yyzz = cbuffer.data(dg_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_y_xy_yzzz = cbuffer.data(dg_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_y_xy_zzzz = cbuffer.data(dg_geom_11_off + 659 * ccomps * dcomps);

            auto g_z_y_xz_xxxx = cbuffer.data(dg_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_y_xz_xxxy = cbuffer.data(dg_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_y_xz_xxxz = cbuffer.data(dg_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_y_xz_xxyy = cbuffer.data(dg_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_y_xz_xxyz = cbuffer.data(dg_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_y_xz_xxzz = cbuffer.data(dg_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_y_xz_xyyy = cbuffer.data(dg_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_y_xz_xyyz = cbuffer.data(dg_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_y_xz_xyzz = cbuffer.data(dg_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_y_xz_xzzz = cbuffer.data(dg_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_y_xz_yyyy = cbuffer.data(dg_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_y_xz_yyyz = cbuffer.data(dg_geom_11_off + 671 * ccomps * dcomps);

            auto g_z_y_xz_yyzz = cbuffer.data(dg_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_y_xz_yzzz = cbuffer.data(dg_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_y_xz_zzzz = cbuffer.data(dg_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_y_yy_xxxx = cbuffer.data(dg_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_y_yy_xxxy = cbuffer.data(dg_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_y_yy_xxxz = cbuffer.data(dg_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_y_yy_xxyy = cbuffer.data(dg_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_y_yy_xxyz = cbuffer.data(dg_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_y_yy_xxzz = cbuffer.data(dg_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_y_yy_xyyy = cbuffer.data(dg_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_y_yy_xyyz = cbuffer.data(dg_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_y_yy_xyzz = cbuffer.data(dg_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_y_yy_xzzz = cbuffer.data(dg_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_y_yy_yyyy = cbuffer.data(dg_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_y_yy_yyyz = cbuffer.data(dg_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_y_yy_yyzz = cbuffer.data(dg_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_y_yy_yzzz = cbuffer.data(dg_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_y_yy_zzzz = cbuffer.data(dg_geom_11_off + 689 * ccomps * dcomps);

            auto g_z_y_yz_xxxx = cbuffer.data(dg_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_y_yz_xxxy = cbuffer.data(dg_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_y_yz_xxxz = cbuffer.data(dg_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_y_yz_xxyy = cbuffer.data(dg_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_y_yz_xxyz = cbuffer.data(dg_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_y_yz_xxzz = cbuffer.data(dg_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_y_yz_xyyy = cbuffer.data(dg_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_y_yz_xyyz = cbuffer.data(dg_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_y_yz_xyzz = cbuffer.data(dg_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_y_yz_xzzz = cbuffer.data(dg_geom_11_off + 699 * ccomps * dcomps);

            auto g_z_y_yz_yyyy = cbuffer.data(dg_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_y_yz_yyyz = cbuffer.data(dg_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_y_yz_yyzz = cbuffer.data(dg_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_y_yz_yzzz = cbuffer.data(dg_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_y_yz_zzzz = cbuffer.data(dg_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_y_zz_xxxx = cbuffer.data(dg_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_y_zz_xxxy = cbuffer.data(dg_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_y_zz_xxxz = cbuffer.data(dg_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_y_zz_xxyy = cbuffer.data(dg_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_y_zz_xxyz = cbuffer.data(dg_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_y_zz_xxzz = cbuffer.data(dg_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_y_zz_xyyy = cbuffer.data(dg_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_y_zz_xyyz = cbuffer.data(dg_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_y_zz_xyzz = cbuffer.data(dg_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_y_zz_xzzz = cbuffer.data(dg_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_y_zz_yyyy = cbuffer.data(dg_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_y_zz_yyyz = cbuffer.data(dg_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_y_zz_yyzz = cbuffer.data(dg_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_y_zz_yzzz = cbuffer.data(dg_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_y_zz_zzzz = cbuffer.data(dg_geom_11_off + 719 * ccomps * dcomps);

            auto g_z_z_xx_xxxx = cbuffer.data(dg_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_z_xx_xxxy = cbuffer.data(dg_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_z_xx_xxxz = cbuffer.data(dg_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_z_xx_xxyy = cbuffer.data(dg_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_z_xx_xxyz = cbuffer.data(dg_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_z_xx_xxzz = cbuffer.data(dg_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_z_xx_xyyy = cbuffer.data(dg_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_z_xx_xyyz = cbuffer.data(dg_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_z_xx_xyzz = cbuffer.data(dg_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_z_xx_xzzz = cbuffer.data(dg_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_z_xx_yyyy = cbuffer.data(dg_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_z_xx_yyyz = cbuffer.data(dg_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_z_xx_yyzz = cbuffer.data(dg_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_z_xx_yzzz = cbuffer.data(dg_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_z_xx_zzzz = cbuffer.data(dg_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_z_xy_xxxx = cbuffer.data(dg_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_z_xy_xxxy = cbuffer.data(dg_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_z_xy_xxxz = cbuffer.data(dg_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_z_xy_xxyy = cbuffer.data(dg_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_z_xy_xxyz = cbuffer.data(dg_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_z_xy_xxzz = cbuffer.data(dg_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_z_xy_xyyy = cbuffer.data(dg_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_z_xy_xyyz = cbuffer.data(dg_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_z_xy_xyzz = cbuffer.data(dg_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_z_xy_xzzz = cbuffer.data(dg_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_z_xy_yyyy = cbuffer.data(dg_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_z_xy_yyyz = cbuffer.data(dg_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_z_xy_yyzz = cbuffer.data(dg_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_z_xy_yzzz = cbuffer.data(dg_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_z_xy_zzzz = cbuffer.data(dg_geom_11_off + 749 * ccomps * dcomps);

            auto g_z_z_xz_xxxx = cbuffer.data(dg_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_z_xz_xxxy = cbuffer.data(dg_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_z_xz_xxxz = cbuffer.data(dg_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_z_xz_xxyy = cbuffer.data(dg_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_z_xz_xxyz = cbuffer.data(dg_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_z_xz_xxzz = cbuffer.data(dg_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_z_xz_xyyy = cbuffer.data(dg_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_z_xz_xyyz = cbuffer.data(dg_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_z_xz_xyzz = cbuffer.data(dg_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_z_xz_xzzz = cbuffer.data(dg_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_z_xz_yyyy = cbuffer.data(dg_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_z_xz_yyyz = cbuffer.data(dg_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_z_xz_yyzz = cbuffer.data(dg_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_z_xz_yzzz = cbuffer.data(dg_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_z_xz_zzzz = cbuffer.data(dg_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_z_yy_xxxx = cbuffer.data(dg_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_z_yy_xxxy = cbuffer.data(dg_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_z_yy_xxxz = cbuffer.data(dg_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_z_yy_xxyy = cbuffer.data(dg_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_z_yy_xxyz = cbuffer.data(dg_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_z_yy_xxzz = cbuffer.data(dg_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_z_yy_xyyy = cbuffer.data(dg_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_z_yy_xyyz = cbuffer.data(dg_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_z_yy_xyzz = cbuffer.data(dg_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_z_yy_xzzz = cbuffer.data(dg_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_z_yy_yyyy = cbuffer.data(dg_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_z_yy_yyyz = cbuffer.data(dg_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_z_yy_yyzz = cbuffer.data(dg_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_z_yy_yzzz = cbuffer.data(dg_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_z_yy_zzzz = cbuffer.data(dg_geom_11_off + 779 * ccomps * dcomps);

            auto g_z_z_yz_xxxx = cbuffer.data(dg_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_z_yz_xxxy = cbuffer.data(dg_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_z_yz_xxxz = cbuffer.data(dg_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_z_yz_xxyy = cbuffer.data(dg_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_z_yz_xxyz = cbuffer.data(dg_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_z_yz_xxzz = cbuffer.data(dg_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_z_yz_xyyy = cbuffer.data(dg_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_z_yz_xyyz = cbuffer.data(dg_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_z_yz_xyzz = cbuffer.data(dg_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_z_yz_xzzz = cbuffer.data(dg_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_z_yz_yyyy = cbuffer.data(dg_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_z_yz_yyyz = cbuffer.data(dg_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_z_yz_yyzz = cbuffer.data(dg_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_z_yz_yzzz = cbuffer.data(dg_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_z_yz_zzzz = cbuffer.data(dg_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_z_zz_xxxx = cbuffer.data(dg_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_z_zz_xxxy = cbuffer.data(dg_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_z_zz_xxxz = cbuffer.data(dg_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_z_zz_xxyy = cbuffer.data(dg_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_z_zz_xxyz = cbuffer.data(dg_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_z_zz_xxzz = cbuffer.data(dg_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_z_zz_xyyy = cbuffer.data(dg_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_z_zz_xyyz = cbuffer.data(dg_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_z_zz_xyzz = cbuffer.data(dg_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_z_zz_xzzz = cbuffer.data(dg_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_z_zz_yyyy = cbuffer.data(dg_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_z_zz_yyyz = cbuffer.data(dg_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_z_zz_yyzz = cbuffer.data(dg_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_z_zz_yzzz = cbuffer.data(dg_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_z_zz_zzzz = cbuffer.data(dg_geom_11_off + 809 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_11_off = idx_geom_11_dhxx + i * dcomps + j;

            auto g_x_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 269 * ccomps * dcomps);

            auto g_x_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 299 * ccomps * dcomps);

            auto g_x_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 314 * ccomps * dcomps);

            auto g_x_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 323 * ccomps * dcomps);

            auto g_x_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 329 * ccomps * dcomps);

            auto g_x_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 335 * ccomps * dcomps);

            auto g_x_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 341 * ccomps * dcomps);

            auto g_x_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 344 * ccomps * dcomps);

            auto g_x_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 347 * ccomps * dcomps);

            auto g_x_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 353 * ccomps * dcomps);

            auto g_x_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 356 * ccomps * dcomps);

            auto g_x_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 359 * ccomps * dcomps);

            auto g_x_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 363 * ccomps * dcomps);

            auto g_x_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 365 * ccomps * dcomps);

            auto g_x_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 371 * ccomps * dcomps);

            auto g_x_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 374 * ccomps * dcomps);

            auto g_x_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 419 * ccomps * dcomps);

            auto g_y_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 539 * ccomps * dcomps);

            auto g_y_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 569 * ccomps * dcomps);

            auto g_y_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 587 * ccomps * dcomps);

            auto g_y_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 599 * ccomps * dcomps);

            auto g_y_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 629 * ccomps * dcomps);

            auto g_y_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 647 * ccomps * dcomps);

            auto g_y_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 650 * ccomps * dcomps);

            auto g_y_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 653 * ccomps * dcomps);

            auto g_y_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 659 * ccomps * dcomps);

            auto g_y_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 665 * ccomps * dcomps);

            auto g_y_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 671 * ccomps * dcomps);

            auto g_y_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 674 * ccomps * dcomps);

            auto g_y_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 677 * ccomps * dcomps);

            auto g_y_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 683 * ccomps * dcomps);

            auto g_y_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 689 * ccomps * dcomps);

            auto g_y_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 692 * ccomps * dcomps);

            auto g_y_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 695 * ccomps * dcomps);

            auto g_y_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 699 * ccomps * dcomps);

            auto g_y_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 701 * ccomps * dcomps);

            auto g_y_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 704 * ccomps * dcomps);

            auto g_y_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 707 * ccomps * dcomps);

            auto g_y_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 713 * ccomps * dcomps);

            auto g_y_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 719 * ccomps * dcomps);

            auto g_y_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 725 * ccomps * dcomps);

            auto g_y_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 727 * ccomps * dcomps);

            auto g_y_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 731 * ccomps * dcomps);

            auto g_y_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 734 * ccomps * dcomps);

            auto g_y_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 737 * ccomps * dcomps);

            auto g_y_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 743 * ccomps * dcomps);

            auto g_y_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 749 * ccomps * dcomps);

            auto g_y_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 779 * ccomps * dcomps);

            auto g_z_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 809 * ccomps * dcomps);

            auto g_z_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 810 * ccomps * dcomps);

            auto g_z_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 811 * ccomps * dcomps);

            auto g_z_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 812 * ccomps * dcomps);

            auto g_z_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 813 * ccomps * dcomps);

            auto g_z_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 814 * ccomps * dcomps);

            auto g_z_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 815 * ccomps * dcomps);

            auto g_z_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 816 * ccomps * dcomps);

            auto g_z_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 817 * ccomps * dcomps);

            auto g_z_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 818 * ccomps * dcomps);

            auto g_z_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 819 * ccomps * dcomps);

            auto g_z_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 820 * ccomps * dcomps);

            auto g_z_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 821 * ccomps * dcomps);

            auto g_z_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 822 * ccomps * dcomps);

            auto g_z_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 823 * ccomps * dcomps);

            auto g_z_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 824 * ccomps * dcomps);

            auto g_z_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 825 * ccomps * dcomps);

            auto g_z_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 826 * ccomps * dcomps);

            auto g_z_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 827 * ccomps * dcomps);

            auto g_z_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 828 * ccomps * dcomps);

            auto g_z_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 829 * ccomps * dcomps);

            auto g_z_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 830 * ccomps * dcomps);

            auto g_z_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 831 * ccomps * dcomps);

            auto g_z_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 832 * ccomps * dcomps);

            auto g_z_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 833 * ccomps * dcomps);

            auto g_z_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 834 * ccomps * dcomps);

            auto g_z_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 835 * ccomps * dcomps);

            auto g_z_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 836 * ccomps * dcomps);

            auto g_z_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 837 * ccomps * dcomps);

            auto g_z_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 838 * ccomps * dcomps);

            auto g_z_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 839 * ccomps * dcomps);

            auto g_z_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 840 * ccomps * dcomps);

            auto g_z_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 841 * ccomps * dcomps);

            auto g_z_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 842 * ccomps * dcomps);

            auto g_z_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 843 * ccomps * dcomps);

            auto g_z_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 844 * ccomps * dcomps);

            auto g_z_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 845 * ccomps * dcomps);

            auto g_z_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 846 * ccomps * dcomps);

            auto g_z_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 847 * ccomps * dcomps);

            auto g_z_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 848 * ccomps * dcomps);

            auto g_z_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 849 * ccomps * dcomps);

            auto g_z_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 850 * ccomps * dcomps);

            auto g_z_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 851 * ccomps * dcomps);

            auto g_z_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 852 * ccomps * dcomps);

            auto g_z_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 853 * ccomps * dcomps);

            auto g_z_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 854 * ccomps * dcomps);

            auto g_z_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 855 * ccomps * dcomps);

            auto g_z_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 856 * ccomps * dcomps);

            auto g_z_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 857 * ccomps * dcomps);

            auto g_z_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 858 * ccomps * dcomps);

            auto g_z_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 859 * ccomps * dcomps);

            auto g_z_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 860 * ccomps * dcomps);

            auto g_z_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 861 * ccomps * dcomps);

            auto g_z_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 862 * ccomps * dcomps);

            auto g_z_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 863 * ccomps * dcomps);

            auto g_z_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 864 * ccomps * dcomps);

            auto g_z_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 865 * ccomps * dcomps);

            auto g_z_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 866 * ccomps * dcomps);

            auto g_z_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 867 * ccomps * dcomps);

            auto g_z_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 868 * ccomps * dcomps);

            auto g_z_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 869 * ccomps * dcomps);

            auto g_z_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 870 * ccomps * dcomps);

            auto g_z_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 871 * ccomps * dcomps);

            auto g_z_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 872 * ccomps * dcomps);

            auto g_z_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 873 * ccomps * dcomps);

            auto g_z_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 874 * ccomps * dcomps);

            auto g_z_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 875 * ccomps * dcomps);

            auto g_z_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 876 * ccomps * dcomps);

            auto g_z_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 877 * ccomps * dcomps);

            auto g_z_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 878 * ccomps * dcomps);

            auto g_z_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 879 * ccomps * dcomps);

            auto g_z_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 880 * ccomps * dcomps);

            auto g_z_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 881 * ccomps * dcomps);

            auto g_z_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 882 * ccomps * dcomps);

            auto g_z_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 883 * ccomps * dcomps);

            auto g_z_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 884 * ccomps * dcomps);

            auto g_z_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 885 * ccomps * dcomps);

            auto g_z_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 886 * ccomps * dcomps);

            auto g_z_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 887 * ccomps * dcomps);

            auto g_z_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 888 * ccomps * dcomps);

            auto g_z_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 889 * ccomps * dcomps);

            auto g_z_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 890 * ccomps * dcomps);

            auto g_z_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 891 * ccomps * dcomps);

            auto g_z_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 892 * ccomps * dcomps);

            auto g_z_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 893 * ccomps * dcomps);

            auto g_z_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 894 * ccomps * dcomps);

            auto g_z_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 895 * ccomps * dcomps);

            auto g_z_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 896 * ccomps * dcomps);

            auto g_z_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 897 * ccomps * dcomps);

            auto g_z_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 898 * ccomps * dcomps);

            auto g_z_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 899 * ccomps * dcomps);

            auto g_z_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 902 * ccomps * dcomps);

            auto g_z_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 905 * ccomps * dcomps);

            auto g_z_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 911 * ccomps * dcomps);

            auto g_z_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 914 * ccomps * dcomps);

            auto g_z_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 917 * ccomps * dcomps);

            auto g_z_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 923 * ccomps * dcomps);

            auto g_z_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 929 * ccomps * dcomps);

            auto g_z_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 935 * ccomps * dcomps);

            auto g_z_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 941 * ccomps * dcomps);

            auto g_z_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 944 * ccomps * dcomps);

            auto g_z_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 947 * ccomps * dcomps);

            auto g_z_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 953 * ccomps * dcomps);

            auto g_z_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 959 * ccomps * dcomps);

            auto g_z_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 965 * ccomps * dcomps);

            auto g_z_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 971 * ccomps * dcomps);

            auto g_z_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 972 * ccomps * dcomps);

            auto g_z_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 973 * ccomps * dcomps);

            auto g_z_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 974 * ccomps * dcomps);

            auto g_z_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 975 * ccomps * dcomps);

            auto g_z_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 976 * ccomps * dcomps);

            auto g_z_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 977 * ccomps * dcomps);

            auto g_z_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 978 * ccomps * dcomps);

            auto g_z_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 979 * ccomps * dcomps);

            auto g_z_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 980 * ccomps * dcomps);

            auto g_z_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 981 * ccomps * dcomps);

            auto g_z_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 982 * ccomps * dcomps);

            auto g_z_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 983 * ccomps * dcomps);

            auto g_z_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 984 * ccomps * dcomps);

            auto g_z_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 985 * ccomps * dcomps);

            auto g_z_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 986 * ccomps * dcomps);

            auto g_z_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 987 * ccomps * dcomps);

            auto g_z_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 988 * ccomps * dcomps);

            auto g_z_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 989 * ccomps * dcomps);

            auto g_z_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 990 * ccomps * dcomps);

            auto g_z_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 991 * ccomps * dcomps);

            auto g_z_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 992 * ccomps * dcomps);

            auto g_z_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 993 * ccomps * dcomps);

            auto g_z_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 994 * ccomps * dcomps);

            auto g_z_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 995 * ccomps * dcomps);

            auto g_z_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 996 * ccomps * dcomps);

            auto g_z_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 997 * ccomps * dcomps);

            auto g_z_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 998 * ccomps * dcomps);

            auto g_z_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 999 * ccomps * dcomps);

            auto g_z_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 1000 * ccomps * dcomps);

            auto g_z_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 1001 * ccomps * dcomps);

            auto g_z_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 1002 * ccomps * dcomps);

            auto g_z_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 1003 * ccomps * dcomps);

            auto g_z_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 1004 * ccomps * dcomps);

            auto g_z_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 1005 * ccomps * dcomps);

            auto g_z_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 1006 * ccomps * dcomps);

            auto g_z_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 1007 * ccomps * dcomps);

            auto g_z_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 1013 * ccomps * dcomps);

            auto g_z_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 1019 * ccomps * dcomps);

            auto g_z_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 1025 * ccomps * dcomps);

            auto g_z_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 1028 * ccomps * dcomps);

            auto g_z_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 1031 * ccomps * dcomps);

            auto g_z_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 1034 * ccomps * dcomps);

            auto g_z_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 1035 * ccomps * dcomps);

            auto g_z_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 1037 * ccomps * dcomps);

            auto g_z_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 1043 * ccomps * dcomps);

            auto g_z_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 1049 * ccomps * dcomps);

            auto g_z_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 1055 * ccomps * dcomps);

            auto g_z_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 1061 * ccomps * dcomps);

            auto g_z_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 1063 * ccomps * dcomps);

            auto g_z_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 1064 * ccomps * dcomps);

            auto g_z_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 1067 * ccomps * dcomps);

            auto g_z_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 1070 * ccomps * dcomps);

            auto g_z_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 1073 * ccomps * dcomps);

            auto g_z_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 1079 * ccomps * dcomps);

            auto g_z_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 1085 * ccomps * dcomps);

            auto g_z_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 1091 * ccomps * dcomps);

            auto g_z_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 1094 * ccomps * dcomps);

            auto g_z_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 1097 * ccomps * dcomps);

            auto g_z_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 1103 * ccomps * dcomps);

            auto g_z_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 1109 * ccomps * dcomps);

            auto g_z_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 1112 * ccomps * dcomps);

            auto g_z_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 1115 * ccomps * dcomps);

            auto g_z_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 1119 * ccomps * dcomps);

            auto g_z_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 1121 * ccomps * dcomps);

            auto g_z_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 1124 * ccomps * dcomps);

            auto g_z_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 1127 * ccomps * dcomps);

            auto g_z_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 1133 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fgxx

            const auto fg_geom_11_off = idx_geom_11_fgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxx_xxxx = cbuffer.data(fg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxx_xxxy = cbuffer.data(fg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxx_xxxz = cbuffer.data(fg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxx_xxyy = cbuffer.data(fg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxx_xxyz = cbuffer.data(fg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxx_xxzz = cbuffer.data(fg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxx_xyyy = cbuffer.data(fg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxx_xyyz = cbuffer.data(fg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxx_xyzz = cbuffer.data(fg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxx_xzzz = cbuffer.data(fg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxx_yyyy = cbuffer.data(fg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxx_yyyz = cbuffer.data(fg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxx_yyzz = cbuffer.data(fg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxx_yzzz = cbuffer.data(fg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxx_zzzz = cbuffer.data(fg_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_xxxx, g_0_x_xx_xxxy, g_0_x_xx_xxxz, g_0_x_xx_xxyy, g_0_x_xx_xxyz, g_0_x_xx_xxzz, g_0_x_xx_xyyy, g_0_x_xx_xyyz, g_0_x_xx_xyzz, g_0_x_xx_xzzz, g_0_x_xx_yyyy, g_0_x_xx_yyyz, g_0_x_xx_yyzz, g_0_x_xx_yzzz, g_0_x_xx_zzzz, g_x_0_xx_xxxx, g_x_0_xx_xxxy, g_x_0_xx_xxxz, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxzz, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyzz, g_x_0_xx_xzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyzz, g_x_0_xx_yzzz, g_x_0_xx_zzzz, g_x_x_xx_xxxx, g_x_x_xx_xxxxx, g_x_x_xx_xxxxy, g_x_x_xx_xxxxz, g_x_x_xx_xxxy, g_x_x_xx_xxxyy, g_x_x_xx_xxxyz, g_x_x_xx_xxxz, g_x_x_xx_xxxzz, g_x_x_xx_xxyy, g_x_x_xx_xxyyy, g_x_x_xx_xxyyz, g_x_x_xx_xxyz, g_x_x_xx_xxyzz, g_x_x_xx_xxzz, g_x_x_xx_xxzzz, g_x_x_xx_xyyy, g_x_x_xx_xyyyy, g_x_x_xx_xyyyz, g_x_x_xx_xyyz, g_x_x_xx_xyyzz, g_x_x_xx_xyzz, g_x_x_xx_xyzzz, g_x_x_xx_xzzz, g_x_x_xx_xzzzz, g_x_x_xx_yyyy, g_x_x_xx_yyyz, g_x_x_xx_yyzz, g_x_x_xx_yzzz, g_x_x_xx_zzzz, g_x_x_xxx_xxxx, g_x_x_xxx_xxxy, g_x_x_xxx_xxxz, g_x_x_xxx_xxyy, g_x_x_xxx_xxyz, g_x_x_xxx_xxzz, g_x_x_xxx_xyyy, g_x_x_xxx_xyyz, g_x_x_xxx_xyzz, g_x_x_xxx_xzzz, g_x_x_xxx_yyyy, g_x_x_xxx_yyyz, g_x_x_xxx_yyzz, g_x_x_xxx_yzzz, g_x_x_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxx_xxxx[k] = -g_0_x_xx_xxxx[k] + g_x_0_xx_xxxx[k] - g_x_x_xx_xxxx[k] * ab_x + g_x_x_xx_xxxxx[k];

                g_x_x_xxx_xxxy[k] = -g_0_x_xx_xxxy[k] + g_x_0_xx_xxxy[k] - g_x_x_xx_xxxy[k] * ab_x + g_x_x_xx_xxxxy[k];

                g_x_x_xxx_xxxz[k] = -g_0_x_xx_xxxz[k] + g_x_0_xx_xxxz[k] - g_x_x_xx_xxxz[k] * ab_x + g_x_x_xx_xxxxz[k];

                g_x_x_xxx_xxyy[k] = -g_0_x_xx_xxyy[k] + g_x_0_xx_xxyy[k] - g_x_x_xx_xxyy[k] * ab_x + g_x_x_xx_xxxyy[k];

                g_x_x_xxx_xxyz[k] = -g_0_x_xx_xxyz[k] + g_x_0_xx_xxyz[k] - g_x_x_xx_xxyz[k] * ab_x + g_x_x_xx_xxxyz[k];

                g_x_x_xxx_xxzz[k] = -g_0_x_xx_xxzz[k] + g_x_0_xx_xxzz[k] - g_x_x_xx_xxzz[k] * ab_x + g_x_x_xx_xxxzz[k];

                g_x_x_xxx_xyyy[k] = -g_0_x_xx_xyyy[k] + g_x_0_xx_xyyy[k] - g_x_x_xx_xyyy[k] * ab_x + g_x_x_xx_xxyyy[k];

                g_x_x_xxx_xyyz[k] = -g_0_x_xx_xyyz[k] + g_x_0_xx_xyyz[k] - g_x_x_xx_xyyz[k] * ab_x + g_x_x_xx_xxyyz[k];

                g_x_x_xxx_xyzz[k] = -g_0_x_xx_xyzz[k] + g_x_0_xx_xyzz[k] - g_x_x_xx_xyzz[k] * ab_x + g_x_x_xx_xxyzz[k];

                g_x_x_xxx_xzzz[k] = -g_0_x_xx_xzzz[k] + g_x_0_xx_xzzz[k] - g_x_x_xx_xzzz[k] * ab_x + g_x_x_xx_xxzzz[k];

                g_x_x_xxx_yyyy[k] = -g_0_x_xx_yyyy[k] + g_x_0_xx_yyyy[k] - g_x_x_xx_yyyy[k] * ab_x + g_x_x_xx_xyyyy[k];

                g_x_x_xxx_yyyz[k] = -g_0_x_xx_yyyz[k] + g_x_0_xx_yyyz[k] - g_x_x_xx_yyyz[k] * ab_x + g_x_x_xx_xyyyz[k];

                g_x_x_xxx_yyzz[k] = -g_0_x_xx_yyzz[k] + g_x_0_xx_yyzz[k] - g_x_x_xx_yyzz[k] * ab_x + g_x_x_xx_xyyzz[k];

                g_x_x_xxx_yzzz[k] = -g_0_x_xx_yzzz[k] + g_x_0_xx_yzzz[k] - g_x_x_xx_yzzz[k] * ab_x + g_x_x_xx_xyzzz[k];

                g_x_x_xxx_zzzz[k] = -g_0_x_xx_zzzz[k] + g_x_0_xx_zzzz[k] - g_x_x_xx_zzzz[k] * ab_x + g_x_x_xx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxy_xxxx = cbuffer.data(fg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxy_xxxy = cbuffer.data(fg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxy_xxxz = cbuffer.data(fg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xxy_xxyy = cbuffer.data(fg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxy_xxyz = cbuffer.data(fg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xxy_xxzz = cbuffer.data(fg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxy_xyyy = cbuffer.data(fg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxy_xyyz = cbuffer.data(fg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxy_xyzz = cbuffer.data(fg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xxy_xzzz = cbuffer.data(fg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxy_yyyy = cbuffer.data(fg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxy_yyyz = cbuffer.data(fg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxy_yyzz = cbuffer.data(fg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxy_yzzz = cbuffer.data(fg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxy_zzzz = cbuffer.data(fg_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_xxxx, g_x_x_xx_xxxxy, g_x_x_xx_xxxy, g_x_x_xx_xxxyy, g_x_x_xx_xxxyz, g_x_x_xx_xxxz, g_x_x_xx_xxyy, g_x_x_xx_xxyyy, g_x_x_xx_xxyyz, g_x_x_xx_xxyz, g_x_x_xx_xxyzz, g_x_x_xx_xxzz, g_x_x_xx_xyyy, g_x_x_xx_xyyyy, g_x_x_xx_xyyyz, g_x_x_xx_xyyz, g_x_x_xx_xyyzz, g_x_x_xx_xyzz, g_x_x_xx_xyzzz, g_x_x_xx_xzzz, g_x_x_xx_yyyy, g_x_x_xx_yyyyy, g_x_x_xx_yyyyz, g_x_x_xx_yyyz, g_x_x_xx_yyyzz, g_x_x_xx_yyzz, g_x_x_xx_yyzzz, g_x_x_xx_yzzz, g_x_x_xx_yzzzz, g_x_x_xx_zzzz, g_x_x_xxy_xxxx, g_x_x_xxy_xxxy, g_x_x_xxy_xxxz, g_x_x_xxy_xxyy, g_x_x_xxy_xxyz, g_x_x_xxy_xxzz, g_x_x_xxy_xyyy, g_x_x_xxy_xyyz, g_x_x_xxy_xyzz, g_x_x_xxy_xzzz, g_x_x_xxy_yyyy, g_x_x_xxy_yyyz, g_x_x_xxy_yyzz, g_x_x_xxy_yzzz, g_x_x_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxy_xxxx[k] = -g_x_x_xx_xxxx[k] * ab_y + g_x_x_xx_xxxxy[k];

                g_x_x_xxy_xxxy[k] = -g_x_x_xx_xxxy[k] * ab_y + g_x_x_xx_xxxyy[k];

                g_x_x_xxy_xxxz[k] = -g_x_x_xx_xxxz[k] * ab_y + g_x_x_xx_xxxyz[k];

                g_x_x_xxy_xxyy[k] = -g_x_x_xx_xxyy[k] * ab_y + g_x_x_xx_xxyyy[k];

                g_x_x_xxy_xxyz[k] = -g_x_x_xx_xxyz[k] * ab_y + g_x_x_xx_xxyyz[k];

                g_x_x_xxy_xxzz[k] = -g_x_x_xx_xxzz[k] * ab_y + g_x_x_xx_xxyzz[k];

                g_x_x_xxy_xyyy[k] = -g_x_x_xx_xyyy[k] * ab_y + g_x_x_xx_xyyyy[k];

                g_x_x_xxy_xyyz[k] = -g_x_x_xx_xyyz[k] * ab_y + g_x_x_xx_xyyyz[k];

                g_x_x_xxy_xyzz[k] = -g_x_x_xx_xyzz[k] * ab_y + g_x_x_xx_xyyzz[k];

                g_x_x_xxy_xzzz[k] = -g_x_x_xx_xzzz[k] * ab_y + g_x_x_xx_xyzzz[k];

                g_x_x_xxy_yyyy[k] = -g_x_x_xx_yyyy[k] * ab_y + g_x_x_xx_yyyyy[k];

                g_x_x_xxy_yyyz[k] = -g_x_x_xx_yyyz[k] * ab_y + g_x_x_xx_yyyyz[k];

                g_x_x_xxy_yyzz[k] = -g_x_x_xx_yyzz[k] * ab_y + g_x_x_xx_yyyzz[k];

                g_x_x_xxy_yzzz[k] = -g_x_x_xx_yzzz[k] * ab_y + g_x_x_xx_yyzzz[k];

                g_x_x_xxy_zzzz[k] = -g_x_x_xx_zzzz[k] * ab_y + g_x_x_xx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxz_xxxx = cbuffer.data(fg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xxz_xxxy = cbuffer.data(fg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xxz_xxxz = cbuffer.data(fg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xxz_xxyy = cbuffer.data(fg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xxz_xxyz = cbuffer.data(fg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xxz_xxzz = cbuffer.data(fg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xxz_xyyy = cbuffer.data(fg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xxz_xyyz = cbuffer.data(fg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xxz_xyzz = cbuffer.data(fg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xxz_xzzz = cbuffer.data(fg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xxz_yyyy = cbuffer.data(fg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xxz_yyyz = cbuffer.data(fg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xxz_yyzz = cbuffer.data(fg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xxz_yzzz = cbuffer.data(fg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xxz_zzzz = cbuffer.data(fg_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_xxxx, g_x_x_xx_xxxxz, g_x_x_xx_xxxy, g_x_x_xx_xxxyz, g_x_x_xx_xxxz, g_x_x_xx_xxxzz, g_x_x_xx_xxyy, g_x_x_xx_xxyyz, g_x_x_xx_xxyz, g_x_x_xx_xxyzz, g_x_x_xx_xxzz, g_x_x_xx_xxzzz, g_x_x_xx_xyyy, g_x_x_xx_xyyyz, g_x_x_xx_xyyz, g_x_x_xx_xyyzz, g_x_x_xx_xyzz, g_x_x_xx_xyzzz, g_x_x_xx_xzzz, g_x_x_xx_xzzzz, g_x_x_xx_yyyy, g_x_x_xx_yyyyz, g_x_x_xx_yyyz, g_x_x_xx_yyyzz, g_x_x_xx_yyzz, g_x_x_xx_yyzzz, g_x_x_xx_yzzz, g_x_x_xx_yzzzz, g_x_x_xx_zzzz, g_x_x_xx_zzzzz, g_x_x_xxz_xxxx, g_x_x_xxz_xxxy, g_x_x_xxz_xxxz, g_x_x_xxz_xxyy, g_x_x_xxz_xxyz, g_x_x_xxz_xxzz, g_x_x_xxz_xyyy, g_x_x_xxz_xyyz, g_x_x_xxz_xyzz, g_x_x_xxz_xzzz, g_x_x_xxz_yyyy, g_x_x_xxz_yyyz, g_x_x_xxz_yyzz, g_x_x_xxz_yzzz, g_x_x_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxz_xxxx[k] = -g_x_x_xx_xxxx[k] * ab_z + g_x_x_xx_xxxxz[k];

                g_x_x_xxz_xxxy[k] = -g_x_x_xx_xxxy[k] * ab_z + g_x_x_xx_xxxyz[k];

                g_x_x_xxz_xxxz[k] = -g_x_x_xx_xxxz[k] * ab_z + g_x_x_xx_xxxzz[k];

                g_x_x_xxz_xxyy[k] = -g_x_x_xx_xxyy[k] * ab_z + g_x_x_xx_xxyyz[k];

                g_x_x_xxz_xxyz[k] = -g_x_x_xx_xxyz[k] * ab_z + g_x_x_xx_xxyzz[k];

                g_x_x_xxz_xxzz[k] = -g_x_x_xx_xxzz[k] * ab_z + g_x_x_xx_xxzzz[k];

                g_x_x_xxz_xyyy[k] = -g_x_x_xx_xyyy[k] * ab_z + g_x_x_xx_xyyyz[k];

                g_x_x_xxz_xyyz[k] = -g_x_x_xx_xyyz[k] * ab_z + g_x_x_xx_xyyzz[k];

                g_x_x_xxz_xyzz[k] = -g_x_x_xx_xyzz[k] * ab_z + g_x_x_xx_xyzzz[k];

                g_x_x_xxz_xzzz[k] = -g_x_x_xx_xzzz[k] * ab_z + g_x_x_xx_xzzzz[k];

                g_x_x_xxz_yyyy[k] = -g_x_x_xx_yyyy[k] * ab_z + g_x_x_xx_yyyyz[k];

                g_x_x_xxz_yyyz[k] = -g_x_x_xx_yyyz[k] * ab_z + g_x_x_xx_yyyzz[k];

                g_x_x_xxz_yyzz[k] = -g_x_x_xx_yyzz[k] * ab_z + g_x_x_xx_yyzzz[k];

                g_x_x_xxz_yzzz[k] = -g_x_x_xx_yzzz[k] * ab_z + g_x_x_xx_yzzzz[k];

                g_x_x_xxz_zzzz[k] = -g_x_x_xx_zzzz[k] * ab_z + g_x_x_xx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyy_xxxx = cbuffer.data(fg_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xyy_xxxy = cbuffer.data(fg_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xyy_xxxz = cbuffer.data(fg_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xyy_xxyy = cbuffer.data(fg_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xyy_xxyz = cbuffer.data(fg_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xyy_xxzz = cbuffer.data(fg_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xyy_xyyy = cbuffer.data(fg_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xyy_xyyz = cbuffer.data(fg_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xyy_xyzz = cbuffer.data(fg_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xyy_xzzz = cbuffer.data(fg_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xyy_yyyy = cbuffer.data(fg_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xyy_yyyz = cbuffer.data(fg_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xyy_yyzz = cbuffer.data(fg_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xyy_yzzz = cbuffer.data(fg_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xyy_zzzz = cbuffer.data(fg_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xy_xxxx, g_x_x_xy_xxxxy, g_x_x_xy_xxxy, g_x_x_xy_xxxyy, g_x_x_xy_xxxyz, g_x_x_xy_xxxz, g_x_x_xy_xxyy, g_x_x_xy_xxyyy, g_x_x_xy_xxyyz, g_x_x_xy_xxyz, g_x_x_xy_xxyzz, g_x_x_xy_xxzz, g_x_x_xy_xyyy, g_x_x_xy_xyyyy, g_x_x_xy_xyyyz, g_x_x_xy_xyyz, g_x_x_xy_xyyzz, g_x_x_xy_xyzz, g_x_x_xy_xyzzz, g_x_x_xy_xzzz, g_x_x_xy_yyyy, g_x_x_xy_yyyyy, g_x_x_xy_yyyyz, g_x_x_xy_yyyz, g_x_x_xy_yyyzz, g_x_x_xy_yyzz, g_x_x_xy_yyzzz, g_x_x_xy_yzzz, g_x_x_xy_yzzzz, g_x_x_xy_zzzz, g_x_x_xyy_xxxx, g_x_x_xyy_xxxy, g_x_x_xyy_xxxz, g_x_x_xyy_xxyy, g_x_x_xyy_xxyz, g_x_x_xyy_xxzz, g_x_x_xyy_xyyy, g_x_x_xyy_xyyz, g_x_x_xyy_xyzz, g_x_x_xyy_xzzz, g_x_x_xyy_yyyy, g_x_x_xyy_yyyz, g_x_x_xyy_yyzz, g_x_x_xyy_yzzz, g_x_x_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyy_xxxx[k] = -g_x_x_xy_xxxx[k] * ab_y + g_x_x_xy_xxxxy[k];

                g_x_x_xyy_xxxy[k] = -g_x_x_xy_xxxy[k] * ab_y + g_x_x_xy_xxxyy[k];

                g_x_x_xyy_xxxz[k] = -g_x_x_xy_xxxz[k] * ab_y + g_x_x_xy_xxxyz[k];

                g_x_x_xyy_xxyy[k] = -g_x_x_xy_xxyy[k] * ab_y + g_x_x_xy_xxyyy[k];

                g_x_x_xyy_xxyz[k] = -g_x_x_xy_xxyz[k] * ab_y + g_x_x_xy_xxyyz[k];

                g_x_x_xyy_xxzz[k] = -g_x_x_xy_xxzz[k] * ab_y + g_x_x_xy_xxyzz[k];

                g_x_x_xyy_xyyy[k] = -g_x_x_xy_xyyy[k] * ab_y + g_x_x_xy_xyyyy[k];

                g_x_x_xyy_xyyz[k] = -g_x_x_xy_xyyz[k] * ab_y + g_x_x_xy_xyyyz[k];

                g_x_x_xyy_xyzz[k] = -g_x_x_xy_xyzz[k] * ab_y + g_x_x_xy_xyyzz[k];

                g_x_x_xyy_xzzz[k] = -g_x_x_xy_xzzz[k] * ab_y + g_x_x_xy_xyzzz[k];

                g_x_x_xyy_yyyy[k] = -g_x_x_xy_yyyy[k] * ab_y + g_x_x_xy_yyyyy[k];

                g_x_x_xyy_yyyz[k] = -g_x_x_xy_yyyz[k] * ab_y + g_x_x_xy_yyyyz[k];

                g_x_x_xyy_yyzz[k] = -g_x_x_xy_yyzz[k] * ab_y + g_x_x_xy_yyyzz[k];

                g_x_x_xyy_yzzz[k] = -g_x_x_xy_yzzz[k] * ab_y + g_x_x_xy_yyzzz[k];

                g_x_x_xyy_zzzz[k] = -g_x_x_xy_zzzz[k] * ab_y + g_x_x_xy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyz_xxxx = cbuffer.data(fg_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xyz_xxxy = cbuffer.data(fg_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xyz_xxxz = cbuffer.data(fg_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_xyz_xxyy = cbuffer.data(fg_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_xyz_xxyz = cbuffer.data(fg_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_xyz_xxzz = cbuffer.data(fg_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_xyz_xyyy = cbuffer.data(fg_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_xyz_xyyz = cbuffer.data(fg_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_xyz_xyzz = cbuffer.data(fg_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_xyz_xzzz = cbuffer.data(fg_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_xyz_yyyy = cbuffer.data(fg_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_xyz_yyyz = cbuffer.data(fg_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_xyz_yyzz = cbuffer.data(fg_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_xyz_yzzz = cbuffer.data(fg_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_xyz_zzzz = cbuffer.data(fg_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyz_xxxx, g_x_x_xyz_xxxy, g_x_x_xyz_xxxz, g_x_x_xyz_xxyy, g_x_x_xyz_xxyz, g_x_x_xyz_xxzz, g_x_x_xyz_xyyy, g_x_x_xyz_xyyz, g_x_x_xyz_xyzz, g_x_x_xyz_xzzz, g_x_x_xyz_yyyy, g_x_x_xyz_yyyz, g_x_x_xyz_yyzz, g_x_x_xyz_yzzz, g_x_x_xyz_zzzz, g_x_x_xz_xxxx, g_x_x_xz_xxxxy, g_x_x_xz_xxxy, g_x_x_xz_xxxyy, g_x_x_xz_xxxyz, g_x_x_xz_xxxz, g_x_x_xz_xxyy, g_x_x_xz_xxyyy, g_x_x_xz_xxyyz, g_x_x_xz_xxyz, g_x_x_xz_xxyzz, g_x_x_xz_xxzz, g_x_x_xz_xyyy, g_x_x_xz_xyyyy, g_x_x_xz_xyyyz, g_x_x_xz_xyyz, g_x_x_xz_xyyzz, g_x_x_xz_xyzz, g_x_x_xz_xyzzz, g_x_x_xz_xzzz, g_x_x_xz_yyyy, g_x_x_xz_yyyyy, g_x_x_xz_yyyyz, g_x_x_xz_yyyz, g_x_x_xz_yyyzz, g_x_x_xz_yyzz, g_x_x_xz_yyzzz, g_x_x_xz_yzzz, g_x_x_xz_yzzzz, g_x_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyz_xxxx[k] = -g_x_x_xz_xxxx[k] * ab_y + g_x_x_xz_xxxxy[k];

                g_x_x_xyz_xxxy[k] = -g_x_x_xz_xxxy[k] * ab_y + g_x_x_xz_xxxyy[k];

                g_x_x_xyz_xxxz[k] = -g_x_x_xz_xxxz[k] * ab_y + g_x_x_xz_xxxyz[k];

                g_x_x_xyz_xxyy[k] = -g_x_x_xz_xxyy[k] * ab_y + g_x_x_xz_xxyyy[k];

                g_x_x_xyz_xxyz[k] = -g_x_x_xz_xxyz[k] * ab_y + g_x_x_xz_xxyyz[k];

                g_x_x_xyz_xxzz[k] = -g_x_x_xz_xxzz[k] * ab_y + g_x_x_xz_xxyzz[k];

                g_x_x_xyz_xyyy[k] = -g_x_x_xz_xyyy[k] * ab_y + g_x_x_xz_xyyyy[k];

                g_x_x_xyz_xyyz[k] = -g_x_x_xz_xyyz[k] * ab_y + g_x_x_xz_xyyyz[k];

                g_x_x_xyz_xyzz[k] = -g_x_x_xz_xyzz[k] * ab_y + g_x_x_xz_xyyzz[k];

                g_x_x_xyz_xzzz[k] = -g_x_x_xz_xzzz[k] * ab_y + g_x_x_xz_xyzzz[k];

                g_x_x_xyz_yyyy[k] = -g_x_x_xz_yyyy[k] * ab_y + g_x_x_xz_yyyyy[k];

                g_x_x_xyz_yyyz[k] = -g_x_x_xz_yyyz[k] * ab_y + g_x_x_xz_yyyyz[k];

                g_x_x_xyz_yyzz[k] = -g_x_x_xz_yyzz[k] * ab_y + g_x_x_xz_yyyzz[k];

                g_x_x_xyz_yzzz[k] = -g_x_x_xz_yzzz[k] * ab_y + g_x_x_xz_yyzzz[k];

                g_x_x_xyz_zzzz[k] = -g_x_x_xz_zzzz[k] * ab_y + g_x_x_xz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzz_xxxx = cbuffer.data(fg_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_xzz_xxxy = cbuffer.data(fg_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_xzz_xxxz = cbuffer.data(fg_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_xzz_xxyy = cbuffer.data(fg_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_xzz_xxyz = cbuffer.data(fg_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_xzz_xxzz = cbuffer.data(fg_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_xzz_xyyy = cbuffer.data(fg_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_xzz_xyyz = cbuffer.data(fg_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_xzz_xyzz = cbuffer.data(fg_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_xzz_xzzz = cbuffer.data(fg_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_xzz_yyyy = cbuffer.data(fg_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_xzz_yyyz = cbuffer.data(fg_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_xzz_yyzz = cbuffer.data(fg_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_xzz_yzzz = cbuffer.data(fg_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_xzz_zzzz = cbuffer.data(fg_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xz_xxxx, g_x_x_xz_xxxxz, g_x_x_xz_xxxy, g_x_x_xz_xxxyz, g_x_x_xz_xxxz, g_x_x_xz_xxxzz, g_x_x_xz_xxyy, g_x_x_xz_xxyyz, g_x_x_xz_xxyz, g_x_x_xz_xxyzz, g_x_x_xz_xxzz, g_x_x_xz_xxzzz, g_x_x_xz_xyyy, g_x_x_xz_xyyyz, g_x_x_xz_xyyz, g_x_x_xz_xyyzz, g_x_x_xz_xyzz, g_x_x_xz_xyzzz, g_x_x_xz_xzzz, g_x_x_xz_xzzzz, g_x_x_xz_yyyy, g_x_x_xz_yyyyz, g_x_x_xz_yyyz, g_x_x_xz_yyyzz, g_x_x_xz_yyzz, g_x_x_xz_yyzzz, g_x_x_xz_yzzz, g_x_x_xz_yzzzz, g_x_x_xz_zzzz, g_x_x_xz_zzzzz, g_x_x_xzz_xxxx, g_x_x_xzz_xxxy, g_x_x_xzz_xxxz, g_x_x_xzz_xxyy, g_x_x_xzz_xxyz, g_x_x_xzz_xxzz, g_x_x_xzz_xyyy, g_x_x_xzz_xyyz, g_x_x_xzz_xyzz, g_x_x_xzz_xzzz, g_x_x_xzz_yyyy, g_x_x_xzz_yyyz, g_x_x_xzz_yyzz, g_x_x_xzz_yzzz, g_x_x_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzz_xxxx[k] = -g_x_x_xz_xxxx[k] * ab_z + g_x_x_xz_xxxxz[k];

                g_x_x_xzz_xxxy[k] = -g_x_x_xz_xxxy[k] * ab_z + g_x_x_xz_xxxyz[k];

                g_x_x_xzz_xxxz[k] = -g_x_x_xz_xxxz[k] * ab_z + g_x_x_xz_xxxzz[k];

                g_x_x_xzz_xxyy[k] = -g_x_x_xz_xxyy[k] * ab_z + g_x_x_xz_xxyyz[k];

                g_x_x_xzz_xxyz[k] = -g_x_x_xz_xxyz[k] * ab_z + g_x_x_xz_xxyzz[k];

                g_x_x_xzz_xxzz[k] = -g_x_x_xz_xxzz[k] * ab_z + g_x_x_xz_xxzzz[k];

                g_x_x_xzz_xyyy[k] = -g_x_x_xz_xyyy[k] * ab_z + g_x_x_xz_xyyyz[k];

                g_x_x_xzz_xyyz[k] = -g_x_x_xz_xyyz[k] * ab_z + g_x_x_xz_xyyzz[k];

                g_x_x_xzz_xyzz[k] = -g_x_x_xz_xyzz[k] * ab_z + g_x_x_xz_xyzzz[k];

                g_x_x_xzz_xzzz[k] = -g_x_x_xz_xzzz[k] * ab_z + g_x_x_xz_xzzzz[k];

                g_x_x_xzz_yyyy[k] = -g_x_x_xz_yyyy[k] * ab_z + g_x_x_xz_yyyyz[k];

                g_x_x_xzz_yyyz[k] = -g_x_x_xz_yyyz[k] * ab_z + g_x_x_xz_yyyzz[k];

                g_x_x_xzz_yyzz[k] = -g_x_x_xz_yyzz[k] * ab_z + g_x_x_xz_yyzzz[k];

                g_x_x_xzz_yzzz[k] = -g_x_x_xz_yzzz[k] * ab_z + g_x_x_xz_yzzzz[k];

                g_x_x_xzz_zzzz[k] = -g_x_x_xz_zzzz[k] * ab_z + g_x_x_xz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyy_xxxx = cbuffer.data(fg_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_yyy_xxxy = cbuffer.data(fg_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_yyy_xxxz = cbuffer.data(fg_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_yyy_xxyy = cbuffer.data(fg_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_yyy_xxyz = cbuffer.data(fg_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_yyy_xxzz = cbuffer.data(fg_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_yyy_xyyy = cbuffer.data(fg_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_yyy_xyyz = cbuffer.data(fg_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_yyy_xyzz = cbuffer.data(fg_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_yyy_xzzz = cbuffer.data(fg_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yyy_yyyy = cbuffer.data(fg_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yyy_yyyz = cbuffer.data(fg_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_yyy_yyzz = cbuffer.data(fg_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yyy_yzzz = cbuffer.data(fg_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yyy_zzzz = cbuffer.data(fg_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yy_xxxx, g_x_x_yy_xxxxy, g_x_x_yy_xxxy, g_x_x_yy_xxxyy, g_x_x_yy_xxxyz, g_x_x_yy_xxxz, g_x_x_yy_xxyy, g_x_x_yy_xxyyy, g_x_x_yy_xxyyz, g_x_x_yy_xxyz, g_x_x_yy_xxyzz, g_x_x_yy_xxzz, g_x_x_yy_xyyy, g_x_x_yy_xyyyy, g_x_x_yy_xyyyz, g_x_x_yy_xyyz, g_x_x_yy_xyyzz, g_x_x_yy_xyzz, g_x_x_yy_xyzzz, g_x_x_yy_xzzz, g_x_x_yy_yyyy, g_x_x_yy_yyyyy, g_x_x_yy_yyyyz, g_x_x_yy_yyyz, g_x_x_yy_yyyzz, g_x_x_yy_yyzz, g_x_x_yy_yyzzz, g_x_x_yy_yzzz, g_x_x_yy_yzzzz, g_x_x_yy_zzzz, g_x_x_yyy_xxxx, g_x_x_yyy_xxxy, g_x_x_yyy_xxxz, g_x_x_yyy_xxyy, g_x_x_yyy_xxyz, g_x_x_yyy_xxzz, g_x_x_yyy_xyyy, g_x_x_yyy_xyyz, g_x_x_yyy_xyzz, g_x_x_yyy_xzzz, g_x_x_yyy_yyyy, g_x_x_yyy_yyyz, g_x_x_yyy_yyzz, g_x_x_yyy_yzzz, g_x_x_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyy_xxxx[k] = -g_x_x_yy_xxxx[k] * ab_y + g_x_x_yy_xxxxy[k];

                g_x_x_yyy_xxxy[k] = -g_x_x_yy_xxxy[k] * ab_y + g_x_x_yy_xxxyy[k];

                g_x_x_yyy_xxxz[k] = -g_x_x_yy_xxxz[k] * ab_y + g_x_x_yy_xxxyz[k];

                g_x_x_yyy_xxyy[k] = -g_x_x_yy_xxyy[k] * ab_y + g_x_x_yy_xxyyy[k];

                g_x_x_yyy_xxyz[k] = -g_x_x_yy_xxyz[k] * ab_y + g_x_x_yy_xxyyz[k];

                g_x_x_yyy_xxzz[k] = -g_x_x_yy_xxzz[k] * ab_y + g_x_x_yy_xxyzz[k];

                g_x_x_yyy_xyyy[k] = -g_x_x_yy_xyyy[k] * ab_y + g_x_x_yy_xyyyy[k];

                g_x_x_yyy_xyyz[k] = -g_x_x_yy_xyyz[k] * ab_y + g_x_x_yy_xyyyz[k];

                g_x_x_yyy_xyzz[k] = -g_x_x_yy_xyzz[k] * ab_y + g_x_x_yy_xyyzz[k];

                g_x_x_yyy_xzzz[k] = -g_x_x_yy_xzzz[k] * ab_y + g_x_x_yy_xyzzz[k];

                g_x_x_yyy_yyyy[k] = -g_x_x_yy_yyyy[k] * ab_y + g_x_x_yy_yyyyy[k];

                g_x_x_yyy_yyyz[k] = -g_x_x_yy_yyyz[k] * ab_y + g_x_x_yy_yyyyz[k];

                g_x_x_yyy_yyzz[k] = -g_x_x_yy_yyzz[k] * ab_y + g_x_x_yy_yyyzz[k];

                g_x_x_yyy_yzzz[k] = -g_x_x_yy_yzzz[k] * ab_y + g_x_x_yy_yyzzz[k];

                g_x_x_yyy_zzzz[k] = -g_x_x_yy_zzzz[k] * ab_y + g_x_x_yy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyz_xxxx = cbuffer.data(fg_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_yyz_xxxy = cbuffer.data(fg_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_yyz_xxxz = cbuffer.data(fg_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_x_yyz_xxyy = cbuffer.data(fg_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_yyz_xxyz = cbuffer.data(fg_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_yyz_xxzz = cbuffer.data(fg_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_yyz_xyyy = cbuffer.data(fg_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_x_yyz_xyyz = cbuffer.data(fg_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_yyz_xyzz = cbuffer.data(fg_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_x_yyz_xzzz = cbuffer.data(fg_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_yyz_yyyy = cbuffer.data(fg_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_yyz_yyyz = cbuffer.data(fg_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_yyz_yyzz = cbuffer.data(fg_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_yyz_yzzz = cbuffer.data(fg_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_yyz_zzzz = cbuffer.data(fg_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyz_xxxx, g_x_x_yyz_xxxy, g_x_x_yyz_xxxz, g_x_x_yyz_xxyy, g_x_x_yyz_xxyz, g_x_x_yyz_xxzz, g_x_x_yyz_xyyy, g_x_x_yyz_xyyz, g_x_x_yyz_xyzz, g_x_x_yyz_xzzz, g_x_x_yyz_yyyy, g_x_x_yyz_yyyz, g_x_x_yyz_yyzz, g_x_x_yyz_yzzz, g_x_x_yyz_zzzz, g_x_x_yz_xxxx, g_x_x_yz_xxxxy, g_x_x_yz_xxxy, g_x_x_yz_xxxyy, g_x_x_yz_xxxyz, g_x_x_yz_xxxz, g_x_x_yz_xxyy, g_x_x_yz_xxyyy, g_x_x_yz_xxyyz, g_x_x_yz_xxyz, g_x_x_yz_xxyzz, g_x_x_yz_xxzz, g_x_x_yz_xyyy, g_x_x_yz_xyyyy, g_x_x_yz_xyyyz, g_x_x_yz_xyyz, g_x_x_yz_xyyzz, g_x_x_yz_xyzz, g_x_x_yz_xyzzz, g_x_x_yz_xzzz, g_x_x_yz_yyyy, g_x_x_yz_yyyyy, g_x_x_yz_yyyyz, g_x_x_yz_yyyz, g_x_x_yz_yyyzz, g_x_x_yz_yyzz, g_x_x_yz_yyzzz, g_x_x_yz_yzzz, g_x_x_yz_yzzzz, g_x_x_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyz_xxxx[k] = -g_x_x_yz_xxxx[k] * ab_y + g_x_x_yz_xxxxy[k];

                g_x_x_yyz_xxxy[k] = -g_x_x_yz_xxxy[k] * ab_y + g_x_x_yz_xxxyy[k];

                g_x_x_yyz_xxxz[k] = -g_x_x_yz_xxxz[k] * ab_y + g_x_x_yz_xxxyz[k];

                g_x_x_yyz_xxyy[k] = -g_x_x_yz_xxyy[k] * ab_y + g_x_x_yz_xxyyy[k];

                g_x_x_yyz_xxyz[k] = -g_x_x_yz_xxyz[k] * ab_y + g_x_x_yz_xxyyz[k];

                g_x_x_yyz_xxzz[k] = -g_x_x_yz_xxzz[k] * ab_y + g_x_x_yz_xxyzz[k];

                g_x_x_yyz_xyyy[k] = -g_x_x_yz_xyyy[k] * ab_y + g_x_x_yz_xyyyy[k];

                g_x_x_yyz_xyyz[k] = -g_x_x_yz_xyyz[k] * ab_y + g_x_x_yz_xyyyz[k];

                g_x_x_yyz_xyzz[k] = -g_x_x_yz_xyzz[k] * ab_y + g_x_x_yz_xyyzz[k];

                g_x_x_yyz_xzzz[k] = -g_x_x_yz_xzzz[k] * ab_y + g_x_x_yz_xyzzz[k];

                g_x_x_yyz_yyyy[k] = -g_x_x_yz_yyyy[k] * ab_y + g_x_x_yz_yyyyy[k];

                g_x_x_yyz_yyyz[k] = -g_x_x_yz_yyyz[k] * ab_y + g_x_x_yz_yyyyz[k];

                g_x_x_yyz_yyzz[k] = -g_x_x_yz_yyzz[k] * ab_y + g_x_x_yz_yyyzz[k];

                g_x_x_yyz_yzzz[k] = -g_x_x_yz_yzzz[k] * ab_y + g_x_x_yz_yyzzz[k];

                g_x_x_yyz_zzzz[k] = -g_x_x_yz_zzzz[k] * ab_y + g_x_x_yz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzz_xxxx = cbuffer.data(fg_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_yzz_xxxy = cbuffer.data(fg_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_yzz_xxxz = cbuffer.data(fg_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_yzz_xxyy = cbuffer.data(fg_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_yzz_xxyz = cbuffer.data(fg_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_yzz_xxzz = cbuffer.data(fg_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_x_yzz_xyyy = cbuffer.data(fg_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_x_yzz_xyyz = cbuffer.data(fg_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_x_yzz_xyzz = cbuffer.data(fg_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_x_yzz_xzzz = cbuffer.data(fg_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_x_yzz_yyyy = cbuffer.data(fg_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_x_yzz_yyyz = cbuffer.data(fg_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_x_yzz_yyzz = cbuffer.data(fg_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_x_yzz_yzzz = cbuffer.data(fg_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_x_yzz_zzzz = cbuffer.data(fg_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzz_xxxx, g_x_x_yzz_xxxy, g_x_x_yzz_xxxz, g_x_x_yzz_xxyy, g_x_x_yzz_xxyz, g_x_x_yzz_xxzz, g_x_x_yzz_xyyy, g_x_x_yzz_xyyz, g_x_x_yzz_xyzz, g_x_x_yzz_xzzz, g_x_x_yzz_yyyy, g_x_x_yzz_yyyz, g_x_x_yzz_yyzz, g_x_x_yzz_yzzz, g_x_x_yzz_zzzz, g_x_x_zz_xxxx, g_x_x_zz_xxxxy, g_x_x_zz_xxxy, g_x_x_zz_xxxyy, g_x_x_zz_xxxyz, g_x_x_zz_xxxz, g_x_x_zz_xxyy, g_x_x_zz_xxyyy, g_x_x_zz_xxyyz, g_x_x_zz_xxyz, g_x_x_zz_xxyzz, g_x_x_zz_xxzz, g_x_x_zz_xyyy, g_x_x_zz_xyyyy, g_x_x_zz_xyyyz, g_x_x_zz_xyyz, g_x_x_zz_xyyzz, g_x_x_zz_xyzz, g_x_x_zz_xyzzz, g_x_x_zz_xzzz, g_x_x_zz_yyyy, g_x_x_zz_yyyyy, g_x_x_zz_yyyyz, g_x_x_zz_yyyz, g_x_x_zz_yyyzz, g_x_x_zz_yyzz, g_x_x_zz_yyzzz, g_x_x_zz_yzzz, g_x_x_zz_yzzzz, g_x_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzz_xxxx[k] = -g_x_x_zz_xxxx[k] * ab_y + g_x_x_zz_xxxxy[k];

                g_x_x_yzz_xxxy[k] = -g_x_x_zz_xxxy[k] * ab_y + g_x_x_zz_xxxyy[k];

                g_x_x_yzz_xxxz[k] = -g_x_x_zz_xxxz[k] * ab_y + g_x_x_zz_xxxyz[k];

                g_x_x_yzz_xxyy[k] = -g_x_x_zz_xxyy[k] * ab_y + g_x_x_zz_xxyyy[k];

                g_x_x_yzz_xxyz[k] = -g_x_x_zz_xxyz[k] * ab_y + g_x_x_zz_xxyyz[k];

                g_x_x_yzz_xxzz[k] = -g_x_x_zz_xxzz[k] * ab_y + g_x_x_zz_xxyzz[k];

                g_x_x_yzz_xyyy[k] = -g_x_x_zz_xyyy[k] * ab_y + g_x_x_zz_xyyyy[k];

                g_x_x_yzz_xyyz[k] = -g_x_x_zz_xyyz[k] * ab_y + g_x_x_zz_xyyyz[k];

                g_x_x_yzz_xyzz[k] = -g_x_x_zz_xyzz[k] * ab_y + g_x_x_zz_xyyzz[k];

                g_x_x_yzz_xzzz[k] = -g_x_x_zz_xzzz[k] * ab_y + g_x_x_zz_xyzzz[k];

                g_x_x_yzz_yyyy[k] = -g_x_x_zz_yyyy[k] * ab_y + g_x_x_zz_yyyyy[k];

                g_x_x_yzz_yyyz[k] = -g_x_x_zz_yyyz[k] * ab_y + g_x_x_zz_yyyyz[k];

                g_x_x_yzz_yyzz[k] = -g_x_x_zz_yyzz[k] * ab_y + g_x_x_zz_yyyzz[k];

                g_x_x_yzz_yzzz[k] = -g_x_x_zz_yzzz[k] * ab_y + g_x_x_zz_yyzzz[k];

                g_x_x_yzz_zzzz[k] = -g_x_x_zz_zzzz[k] * ab_y + g_x_x_zz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzz_xxxx = cbuffer.data(fg_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_x_zzz_xxxy = cbuffer.data(fg_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_x_zzz_xxxz = cbuffer.data(fg_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_x_zzz_xxyy = cbuffer.data(fg_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_x_zzz_xxyz = cbuffer.data(fg_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_x_zzz_xxzz = cbuffer.data(fg_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_x_zzz_xyyy = cbuffer.data(fg_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_x_zzz_xyyz = cbuffer.data(fg_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_x_zzz_xyzz = cbuffer.data(fg_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_x_zzz_xzzz = cbuffer.data(fg_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_x_zzz_yyyy = cbuffer.data(fg_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_x_zzz_yyyz = cbuffer.data(fg_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_x_zzz_yyzz = cbuffer.data(fg_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_x_zzz_yzzz = cbuffer.data(fg_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_x_zzz_zzzz = cbuffer.data(fg_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zz_xxxx, g_x_x_zz_xxxxz, g_x_x_zz_xxxy, g_x_x_zz_xxxyz, g_x_x_zz_xxxz, g_x_x_zz_xxxzz, g_x_x_zz_xxyy, g_x_x_zz_xxyyz, g_x_x_zz_xxyz, g_x_x_zz_xxyzz, g_x_x_zz_xxzz, g_x_x_zz_xxzzz, g_x_x_zz_xyyy, g_x_x_zz_xyyyz, g_x_x_zz_xyyz, g_x_x_zz_xyyzz, g_x_x_zz_xyzz, g_x_x_zz_xyzzz, g_x_x_zz_xzzz, g_x_x_zz_xzzzz, g_x_x_zz_yyyy, g_x_x_zz_yyyyz, g_x_x_zz_yyyz, g_x_x_zz_yyyzz, g_x_x_zz_yyzz, g_x_x_zz_yyzzz, g_x_x_zz_yzzz, g_x_x_zz_yzzzz, g_x_x_zz_zzzz, g_x_x_zz_zzzzz, g_x_x_zzz_xxxx, g_x_x_zzz_xxxy, g_x_x_zzz_xxxz, g_x_x_zzz_xxyy, g_x_x_zzz_xxyz, g_x_x_zzz_xxzz, g_x_x_zzz_xyyy, g_x_x_zzz_xyyz, g_x_x_zzz_xyzz, g_x_x_zzz_xzzz, g_x_x_zzz_yyyy, g_x_x_zzz_yyyz, g_x_x_zzz_yyzz, g_x_x_zzz_yzzz, g_x_x_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzz_xxxx[k] = -g_x_x_zz_xxxx[k] * ab_z + g_x_x_zz_xxxxz[k];

                g_x_x_zzz_xxxy[k] = -g_x_x_zz_xxxy[k] * ab_z + g_x_x_zz_xxxyz[k];

                g_x_x_zzz_xxxz[k] = -g_x_x_zz_xxxz[k] * ab_z + g_x_x_zz_xxxzz[k];

                g_x_x_zzz_xxyy[k] = -g_x_x_zz_xxyy[k] * ab_z + g_x_x_zz_xxyyz[k];

                g_x_x_zzz_xxyz[k] = -g_x_x_zz_xxyz[k] * ab_z + g_x_x_zz_xxyzz[k];

                g_x_x_zzz_xxzz[k] = -g_x_x_zz_xxzz[k] * ab_z + g_x_x_zz_xxzzz[k];

                g_x_x_zzz_xyyy[k] = -g_x_x_zz_xyyy[k] * ab_z + g_x_x_zz_xyyyz[k];

                g_x_x_zzz_xyyz[k] = -g_x_x_zz_xyyz[k] * ab_z + g_x_x_zz_xyyzz[k];

                g_x_x_zzz_xyzz[k] = -g_x_x_zz_xyzz[k] * ab_z + g_x_x_zz_xyzzz[k];

                g_x_x_zzz_xzzz[k] = -g_x_x_zz_xzzz[k] * ab_z + g_x_x_zz_xzzzz[k];

                g_x_x_zzz_yyyy[k] = -g_x_x_zz_yyyy[k] * ab_z + g_x_x_zz_yyyyz[k];

                g_x_x_zzz_yyyz[k] = -g_x_x_zz_yyyz[k] * ab_z + g_x_x_zz_yyyzz[k];

                g_x_x_zzz_yyzz[k] = -g_x_x_zz_yyzz[k] * ab_z + g_x_x_zz_yyzzz[k];

                g_x_x_zzz_yzzz[k] = -g_x_x_zz_yzzz[k] * ab_z + g_x_x_zz_yzzzz[k];

                g_x_x_zzz_zzzz[k] = -g_x_x_zz_zzzz[k] * ab_z + g_x_x_zz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxx_xxxx = cbuffer.data(fg_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xxx_xxxy = cbuffer.data(fg_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xxx_xxxz = cbuffer.data(fg_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xxx_xxyy = cbuffer.data(fg_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xxx_xxyz = cbuffer.data(fg_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xxx_xxzz = cbuffer.data(fg_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_xxx_xyyy = cbuffer.data(fg_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xxx_xyyz = cbuffer.data(fg_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xxx_xyzz = cbuffer.data(fg_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xxx_xzzz = cbuffer.data(fg_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_xxx_yyyy = cbuffer.data(fg_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_xxx_yyyz = cbuffer.data(fg_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_xxx_yyzz = cbuffer.data(fg_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_xxx_yzzz = cbuffer.data(fg_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_xxx_zzzz = cbuffer.data(fg_geom_11_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_xxxx, g_0_y_xx_xxxy, g_0_y_xx_xxxz, g_0_y_xx_xxyy, g_0_y_xx_xxyz, g_0_y_xx_xxzz, g_0_y_xx_xyyy, g_0_y_xx_xyyz, g_0_y_xx_xyzz, g_0_y_xx_xzzz, g_0_y_xx_yyyy, g_0_y_xx_yyyz, g_0_y_xx_yyzz, g_0_y_xx_yzzz, g_0_y_xx_zzzz, g_x_y_xx_xxxx, g_x_y_xx_xxxxx, g_x_y_xx_xxxxy, g_x_y_xx_xxxxz, g_x_y_xx_xxxy, g_x_y_xx_xxxyy, g_x_y_xx_xxxyz, g_x_y_xx_xxxz, g_x_y_xx_xxxzz, g_x_y_xx_xxyy, g_x_y_xx_xxyyy, g_x_y_xx_xxyyz, g_x_y_xx_xxyz, g_x_y_xx_xxyzz, g_x_y_xx_xxzz, g_x_y_xx_xxzzz, g_x_y_xx_xyyy, g_x_y_xx_xyyyy, g_x_y_xx_xyyyz, g_x_y_xx_xyyz, g_x_y_xx_xyyzz, g_x_y_xx_xyzz, g_x_y_xx_xyzzz, g_x_y_xx_xzzz, g_x_y_xx_xzzzz, g_x_y_xx_yyyy, g_x_y_xx_yyyz, g_x_y_xx_yyzz, g_x_y_xx_yzzz, g_x_y_xx_zzzz, g_x_y_xxx_xxxx, g_x_y_xxx_xxxy, g_x_y_xxx_xxxz, g_x_y_xxx_xxyy, g_x_y_xxx_xxyz, g_x_y_xxx_xxzz, g_x_y_xxx_xyyy, g_x_y_xxx_xyyz, g_x_y_xxx_xyzz, g_x_y_xxx_xzzz, g_x_y_xxx_yyyy, g_x_y_xxx_yyyz, g_x_y_xxx_yyzz, g_x_y_xxx_yzzz, g_x_y_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxx_xxxx[k] = -g_0_y_xx_xxxx[k] - g_x_y_xx_xxxx[k] * ab_x + g_x_y_xx_xxxxx[k];

                g_x_y_xxx_xxxy[k] = -g_0_y_xx_xxxy[k] - g_x_y_xx_xxxy[k] * ab_x + g_x_y_xx_xxxxy[k];

                g_x_y_xxx_xxxz[k] = -g_0_y_xx_xxxz[k] - g_x_y_xx_xxxz[k] * ab_x + g_x_y_xx_xxxxz[k];

                g_x_y_xxx_xxyy[k] = -g_0_y_xx_xxyy[k] - g_x_y_xx_xxyy[k] * ab_x + g_x_y_xx_xxxyy[k];

                g_x_y_xxx_xxyz[k] = -g_0_y_xx_xxyz[k] - g_x_y_xx_xxyz[k] * ab_x + g_x_y_xx_xxxyz[k];

                g_x_y_xxx_xxzz[k] = -g_0_y_xx_xxzz[k] - g_x_y_xx_xxzz[k] * ab_x + g_x_y_xx_xxxzz[k];

                g_x_y_xxx_xyyy[k] = -g_0_y_xx_xyyy[k] - g_x_y_xx_xyyy[k] * ab_x + g_x_y_xx_xxyyy[k];

                g_x_y_xxx_xyyz[k] = -g_0_y_xx_xyyz[k] - g_x_y_xx_xyyz[k] * ab_x + g_x_y_xx_xxyyz[k];

                g_x_y_xxx_xyzz[k] = -g_0_y_xx_xyzz[k] - g_x_y_xx_xyzz[k] * ab_x + g_x_y_xx_xxyzz[k];

                g_x_y_xxx_xzzz[k] = -g_0_y_xx_xzzz[k] - g_x_y_xx_xzzz[k] * ab_x + g_x_y_xx_xxzzz[k];

                g_x_y_xxx_yyyy[k] = -g_0_y_xx_yyyy[k] - g_x_y_xx_yyyy[k] * ab_x + g_x_y_xx_xyyyy[k];

                g_x_y_xxx_yyyz[k] = -g_0_y_xx_yyyz[k] - g_x_y_xx_yyyz[k] * ab_x + g_x_y_xx_xyyyz[k];

                g_x_y_xxx_yyzz[k] = -g_0_y_xx_yyzz[k] - g_x_y_xx_yyzz[k] * ab_x + g_x_y_xx_xyyzz[k];

                g_x_y_xxx_yzzz[k] = -g_0_y_xx_yzzz[k] - g_x_y_xx_yzzz[k] * ab_x + g_x_y_xx_xyzzz[k];

                g_x_y_xxx_zzzz[k] = -g_0_y_xx_zzzz[k] - g_x_y_xx_zzzz[k] * ab_x + g_x_y_xx_xzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxy_xxxx = cbuffer.data(fg_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_xxy_xxxy = cbuffer.data(fg_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_xxy_xxxz = cbuffer.data(fg_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_xxy_xxyy = cbuffer.data(fg_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xxy_xxyz = cbuffer.data(fg_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xxy_xxzz = cbuffer.data(fg_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xxy_xyyy = cbuffer.data(fg_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xxy_xyyz = cbuffer.data(fg_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xxy_xyzz = cbuffer.data(fg_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_xxy_xzzz = cbuffer.data(fg_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xxy_yyyy = cbuffer.data(fg_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xxy_yyyz = cbuffer.data(fg_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xxy_yyzz = cbuffer.data(fg_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xxy_yzzz = cbuffer.data(fg_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xxy_zzzz = cbuffer.data(fg_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_xxxx, g_0_y_xy_xxxy, g_0_y_xy_xxxz, g_0_y_xy_xxyy, g_0_y_xy_xxyz, g_0_y_xy_xxzz, g_0_y_xy_xyyy, g_0_y_xy_xyyz, g_0_y_xy_xyzz, g_0_y_xy_xzzz, g_0_y_xy_yyyy, g_0_y_xy_yyyz, g_0_y_xy_yyzz, g_0_y_xy_yzzz, g_0_y_xy_zzzz, g_x_y_xxy_xxxx, g_x_y_xxy_xxxy, g_x_y_xxy_xxxz, g_x_y_xxy_xxyy, g_x_y_xxy_xxyz, g_x_y_xxy_xxzz, g_x_y_xxy_xyyy, g_x_y_xxy_xyyz, g_x_y_xxy_xyzz, g_x_y_xxy_xzzz, g_x_y_xxy_yyyy, g_x_y_xxy_yyyz, g_x_y_xxy_yyzz, g_x_y_xxy_yzzz, g_x_y_xxy_zzzz, g_x_y_xy_xxxx, g_x_y_xy_xxxxx, g_x_y_xy_xxxxy, g_x_y_xy_xxxxz, g_x_y_xy_xxxy, g_x_y_xy_xxxyy, g_x_y_xy_xxxyz, g_x_y_xy_xxxz, g_x_y_xy_xxxzz, g_x_y_xy_xxyy, g_x_y_xy_xxyyy, g_x_y_xy_xxyyz, g_x_y_xy_xxyz, g_x_y_xy_xxyzz, g_x_y_xy_xxzz, g_x_y_xy_xxzzz, g_x_y_xy_xyyy, g_x_y_xy_xyyyy, g_x_y_xy_xyyyz, g_x_y_xy_xyyz, g_x_y_xy_xyyzz, g_x_y_xy_xyzz, g_x_y_xy_xyzzz, g_x_y_xy_xzzz, g_x_y_xy_xzzzz, g_x_y_xy_yyyy, g_x_y_xy_yyyz, g_x_y_xy_yyzz, g_x_y_xy_yzzz, g_x_y_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxy_xxxx[k] = -g_0_y_xy_xxxx[k] - g_x_y_xy_xxxx[k] * ab_x + g_x_y_xy_xxxxx[k];

                g_x_y_xxy_xxxy[k] = -g_0_y_xy_xxxy[k] - g_x_y_xy_xxxy[k] * ab_x + g_x_y_xy_xxxxy[k];

                g_x_y_xxy_xxxz[k] = -g_0_y_xy_xxxz[k] - g_x_y_xy_xxxz[k] * ab_x + g_x_y_xy_xxxxz[k];

                g_x_y_xxy_xxyy[k] = -g_0_y_xy_xxyy[k] - g_x_y_xy_xxyy[k] * ab_x + g_x_y_xy_xxxyy[k];

                g_x_y_xxy_xxyz[k] = -g_0_y_xy_xxyz[k] - g_x_y_xy_xxyz[k] * ab_x + g_x_y_xy_xxxyz[k];

                g_x_y_xxy_xxzz[k] = -g_0_y_xy_xxzz[k] - g_x_y_xy_xxzz[k] * ab_x + g_x_y_xy_xxxzz[k];

                g_x_y_xxy_xyyy[k] = -g_0_y_xy_xyyy[k] - g_x_y_xy_xyyy[k] * ab_x + g_x_y_xy_xxyyy[k];

                g_x_y_xxy_xyyz[k] = -g_0_y_xy_xyyz[k] - g_x_y_xy_xyyz[k] * ab_x + g_x_y_xy_xxyyz[k];

                g_x_y_xxy_xyzz[k] = -g_0_y_xy_xyzz[k] - g_x_y_xy_xyzz[k] * ab_x + g_x_y_xy_xxyzz[k];

                g_x_y_xxy_xzzz[k] = -g_0_y_xy_xzzz[k] - g_x_y_xy_xzzz[k] * ab_x + g_x_y_xy_xxzzz[k];

                g_x_y_xxy_yyyy[k] = -g_0_y_xy_yyyy[k] - g_x_y_xy_yyyy[k] * ab_x + g_x_y_xy_xyyyy[k];

                g_x_y_xxy_yyyz[k] = -g_0_y_xy_yyyz[k] - g_x_y_xy_yyyz[k] * ab_x + g_x_y_xy_xyyyz[k];

                g_x_y_xxy_yyzz[k] = -g_0_y_xy_yyzz[k] - g_x_y_xy_yyzz[k] * ab_x + g_x_y_xy_xyyzz[k];

                g_x_y_xxy_yzzz[k] = -g_0_y_xy_yzzz[k] - g_x_y_xy_yzzz[k] * ab_x + g_x_y_xy_xyzzz[k];

                g_x_y_xxy_zzzz[k] = -g_0_y_xy_zzzz[k] - g_x_y_xy_zzzz[k] * ab_x + g_x_y_xy_xzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxz_xxxx = cbuffer.data(fg_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xxz_xxxy = cbuffer.data(fg_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xxz_xxxz = cbuffer.data(fg_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xxz_xxyy = cbuffer.data(fg_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xxz_xxyz = cbuffer.data(fg_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xxz_xxzz = cbuffer.data(fg_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_xxz_xyyy = cbuffer.data(fg_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xxz_xyyz = cbuffer.data(fg_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xxz_xyzz = cbuffer.data(fg_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_xxz_xzzz = cbuffer.data(fg_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_xxz_yyyy = cbuffer.data(fg_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_xxz_yyyz = cbuffer.data(fg_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_xxz_yyzz = cbuffer.data(fg_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_xxz_yzzz = cbuffer.data(fg_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_xxz_zzzz = cbuffer.data(fg_geom_11_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xx_xxxx, g_x_y_xx_xxxxz, g_x_y_xx_xxxy, g_x_y_xx_xxxyz, g_x_y_xx_xxxz, g_x_y_xx_xxxzz, g_x_y_xx_xxyy, g_x_y_xx_xxyyz, g_x_y_xx_xxyz, g_x_y_xx_xxyzz, g_x_y_xx_xxzz, g_x_y_xx_xxzzz, g_x_y_xx_xyyy, g_x_y_xx_xyyyz, g_x_y_xx_xyyz, g_x_y_xx_xyyzz, g_x_y_xx_xyzz, g_x_y_xx_xyzzz, g_x_y_xx_xzzz, g_x_y_xx_xzzzz, g_x_y_xx_yyyy, g_x_y_xx_yyyyz, g_x_y_xx_yyyz, g_x_y_xx_yyyzz, g_x_y_xx_yyzz, g_x_y_xx_yyzzz, g_x_y_xx_yzzz, g_x_y_xx_yzzzz, g_x_y_xx_zzzz, g_x_y_xx_zzzzz, g_x_y_xxz_xxxx, g_x_y_xxz_xxxy, g_x_y_xxz_xxxz, g_x_y_xxz_xxyy, g_x_y_xxz_xxyz, g_x_y_xxz_xxzz, g_x_y_xxz_xyyy, g_x_y_xxz_xyyz, g_x_y_xxz_xyzz, g_x_y_xxz_xzzz, g_x_y_xxz_yyyy, g_x_y_xxz_yyyz, g_x_y_xxz_yyzz, g_x_y_xxz_yzzz, g_x_y_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxz_xxxx[k] = -g_x_y_xx_xxxx[k] * ab_z + g_x_y_xx_xxxxz[k];

                g_x_y_xxz_xxxy[k] = -g_x_y_xx_xxxy[k] * ab_z + g_x_y_xx_xxxyz[k];

                g_x_y_xxz_xxxz[k] = -g_x_y_xx_xxxz[k] * ab_z + g_x_y_xx_xxxzz[k];

                g_x_y_xxz_xxyy[k] = -g_x_y_xx_xxyy[k] * ab_z + g_x_y_xx_xxyyz[k];

                g_x_y_xxz_xxyz[k] = -g_x_y_xx_xxyz[k] * ab_z + g_x_y_xx_xxyzz[k];

                g_x_y_xxz_xxzz[k] = -g_x_y_xx_xxzz[k] * ab_z + g_x_y_xx_xxzzz[k];

                g_x_y_xxz_xyyy[k] = -g_x_y_xx_xyyy[k] * ab_z + g_x_y_xx_xyyyz[k];

                g_x_y_xxz_xyyz[k] = -g_x_y_xx_xyyz[k] * ab_z + g_x_y_xx_xyyzz[k];

                g_x_y_xxz_xyzz[k] = -g_x_y_xx_xyzz[k] * ab_z + g_x_y_xx_xyzzz[k];

                g_x_y_xxz_xzzz[k] = -g_x_y_xx_xzzz[k] * ab_z + g_x_y_xx_xzzzz[k];

                g_x_y_xxz_yyyy[k] = -g_x_y_xx_yyyy[k] * ab_z + g_x_y_xx_yyyyz[k];

                g_x_y_xxz_yyyz[k] = -g_x_y_xx_yyyz[k] * ab_z + g_x_y_xx_yyyzz[k];

                g_x_y_xxz_yyzz[k] = -g_x_y_xx_yyzz[k] * ab_z + g_x_y_xx_yyzzz[k];

                g_x_y_xxz_yzzz[k] = -g_x_y_xx_yzzz[k] * ab_z + g_x_y_xx_yzzzz[k];

                g_x_y_xxz_zzzz[k] = -g_x_y_xx_zzzz[k] * ab_z + g_x_y_xx_zzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyy_xxxx = cbuffer.data(fg_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_xyy_xxxy = cbuffer.data(fg_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_xyy_xxxz = cbuffer.data(fg_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_xyy_xxyy = cbuffer.data(fg_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_xyy_xxyz = cbuffer.data(fg_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_xyy_xxzz = cbuffer.data(fg_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_xyy_xyyy = cbuffer.data(fg_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_xyy_xyyz = cbuffer.data(fg_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_xyy_xyzz = cbuffer.data(fg_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_xyy_xzzz = cbuffer.data(fg_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_xyy_yyyy = cbuffer.data(fg_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_xyy_yyyz = cbuffer.data(fg_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_xyy_yyzz = cbuffer.data(fg_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_xyy_yzzz = cbuffer.data(fg_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_xyy_zzzz = cbuffer.data(fg_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxxx, g_0_y_yy_xxxy, g_0_y_yy_xxxz, g_0_y_yy_xxyy, g_0_y_yy_xxyz, g_0_y_yy_xxzz, g_0_y_yy_xyyy, g_0_y_yy_xyyz, g_0_y_yy_xyzz, g_0_y_yy_xzzz, g_0_y_yy_yyyy, g_0_y_yy_yyyz, g_0_y_yy_yyzz, g_0_y_yy_yzzz, g_0_y_yy_zzzz, g_x_y_xyy_xxxx, g_x_y_xyy_xxxy, g_x_y_xyy_xxxz, g_x_y_xyy_xxyy, g_x_y_xyy_xxyz, g_x_y_xyy_xxzz, g_x_y_xyy_xyyy, g_x_y_xyy_xyyz, g_x_y_xyy_xyzz, g_x_y_xyy_xzzz, g_x_y_xyy_yyyy, g_x_y_xyy_yyyz, g_x_y_xyy_yyzz, g_x_y_xyy_yzzz, g_x_y_xyy_zzzz, g_x_y_yy_xxxx, g_x_y_yy_xxxxx, g_x_y_yy_xxxxy, g_x_y_yy_xxxxz, g_x_y_yy_xxxy, g_x_y_yy_xxxyy, g_x_y_yy_xxxyz, g_x_y_yy_xxxz, g_x_y_yy_xxxzz, g_x_y_yy_xxyy, g_x_y_yy_xxyyy, g_x_y_yy_xxyyz, g_x_y_yy_xxyz, g_x_y_yy_xxyzz, g_x_y_yy_xxzz, g_x_y_yy_xxzzz, g_x_y_yy_xyyy, g_x_y_yy_xyyyy, g_x_y_yy_xyyyz, g_x_y_yy_xyyz, g_x_y_yy_xyyzz, g_x_y_yy_xyzz, g_x_y_yy_xyzzz, g_x_y_yy_xzzz, g_x_y_yy_xzzzz, g_x_y_yy_yyyy, g_x_y_yy_yyyz, g_x_y_yy_yyzz, g_x_y_yy_yzzz, g_x_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyy_xxxx[k] = -g_0_y_yy_xxxx[k] - g_x_y_yy_xxxx[k] * ab_x + g_x_y_yy_xxxxx[k];

                g_x_y_xyy_xxxy[k] = -g_0_y_yy_xxxy[k] - g_x_y_yy_xxxy[k] * ab_x + g_x_y_yy_xxxxy[k];

                g_x_y_xyy_xxxz[k] = -g_0_y_yy_xxxz[k] - g_x_y_yy_xxxz[k] * ab_x + g_x_y_yy_xxxxz[k];

                g_x_y_xyy_xxyy[k] = -g_0_y_yy_xxyy[k] - g_x_y_yy_xxyy[k] * ab_x + g_x_y_yy_xxxyy[k];

                g_x_y_xyy_xxyz[k] = -g_0_y_yy_xxyz[k] - g_x_y_yy_xxyz[k] * ab_x + g_x_y_yy_xxxyz[k];

                g_x_y_xyy_xxzz[k] = -g_0_y_yy_xxzz[k] - g_x_y_yy_xxzz[k] * ab_x + g_x_y_yy_xxxzz[k];

                g_x_y_xyy_xyyy[k] = -g_0_y_yy_xyyy[k] - g_x_y_yy_xyyy[k] * ab_x + g_x_y_yy_xxyyy[k];

                g_x_y_xyy_xyyz[k] = -g_0_y_yy_xyyz[k] - g_x_y_yy_xyyz[k] * ab_x + g_x_y_yy_xxyyz[k];

                g_x_y_xyy_xyzz[k] = -g_0_y_yy_xyzz[k] - g_x_y_yy_xyzz[k] * ab_x + g_x_y_yy_xxyzz[k];

                g_x_y_xyy_xzzz[k] = -g_0_y_yy_xzzz[k] - g_x_y_yy_xzzz[k] * ab_x + g_x_y_yy_xxzzz[k];

                g_x_y_xyy_yyyy[k] = -g_0_y_yy_yyyy[k] - g_x_y_yy_yyyy[k] * ab_x + g_x_y_yy_xyyyy[k];

                g_x_y_xyy_yyyz[k] = -g_0_y_yy_yyyz[k] - g_x_y_yy_yyyz[k] * ab_x + g_x_y_yy_xyyyz[k];

                g_x_y_xyy_yyzz[k] = -g_0_y_yy_yyzz[k] - g_x_y_yy_yyzz[k] * ab_x + g_x_y_yy_xyyzz[k];

                g_x_y_xyy_yzzz[k] = -g_0_y_yy_yzzz[k] - g_x_y_yy_yzzz[k] * ab_x + g_x_y_yy_xyzzz[k];

                g_x_y_xyy_zzzz[k] = -g_0_y_yy_zzzz[k] - g_x_y_yy_zzzz[k] * ab_x + g_x_y_yy_xzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyz_xxxx = cbuffer.data(fg_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_xyz_xxxy = cbuffer.data(fg_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_xyz_xxxz = cbuffer.data(fg_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_xyz_xxyy = cbuffer.data(fg_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_xyz_xxyz = cbuffer.data(fg_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_xyz_xxzz = cbuffer.data(fg_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_y_xyz_xyyy = cbuffer.data(fg_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_xyz_xyyz = cbuffer.data(fg_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_xyz_xyzz = cbuffer.data(fg_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_xyz_xzzz = cbuffer.data(fg_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_xyz_yyyy = cbuffer.data(fg_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_xyz_yyyz = cbuffer.data(fg_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_y_xyz_yyzz = cbuffer.data(fg_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_xyz_yzzz = cbuffer.data(fg_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_y_xyz_zzzz = cbuffer.data(fg_geom_11_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xy_xxxx, g_x_y_xy_xxxxz, g_x_y_xy_xxxy, g_x_y_xy_xxxyz, g_x_y_xy_xxxz, g_x_y_xy_xxxzz, g_x_y_xy_xxyy, g_x_y_xy_xxyyz, g_x_y_xy_xxyz, g_x_y_xy_xxyzz, g_x_y_xy_xxzz, g_x_y_xy_xxzzz, g_x_y_xy_xyyy, g_x_y_xy_xyyyz, g_x_y_xy_xyyz, g_x_y_xy_xyyzz, g_x_y_xy_xyzz, g_x_y_xy_xyzzz, g_x_y_xy_xzzz, g_x_y_xy_xzzzz, g_x_y_xy_yyyy, g_x_y_xy_yyyyz, g_x_y_xy_yyyz, g_x_y_xy_yyyzz, g_x_y_xy_yyzz, g_x_y_xy_yyzzz, g_x_y_xy_yzzz, g_x_y_xy_yzzzz, g_x_y_xy_zzzz, g_x_y_xy_zzzzz, g_x_y_xyz_xxxx, g_x_y_xyz_xxxy, g_x_y_xyz_xxxz, g_x_y_xyz_xxyy, g_x_y_xyz_xxyz, g_x_y_xyz_xxzz, g_x_y_xyz_xyyy, g_x_y_xyz_xyyz, g_x_y_xyz_xyzz, g_x_y_xyz_xzzz, g_x_y_xyz_yyyy, g_x_y_xyz_yyyz, g_x_y_xyz_yyzz, g_x_y_xyz_yzzz, g_x_y_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyz_xxxx[k] = -g_x_y_xy_xxxx[k] * ab_z + g_x_y_xy_xxxxz[k];

                g_x_y_xyz_xxxy[k] = -g_x_y_xy_xxxy[k] * ab_z + g_x_y_xy_xxxyz[k];

                g_x_y_xyz_xxxz[k] = -g_x_y_xy_xxxz[k] * ab_z + g_x_y_xy_xxxzz[k];

                g_x_y_xyz_xxyy[k] = -g_x_y_xy_xxyy[k] * ab_z + g_x_y_xy_xxyyz[k];

                g_x_y_xyz_xxyz[k] = -g_x_y_xy_xxyz[k] * ab_z + g_x_y_xy_xxyzz[k];

                g_x_y_xyz_xxzz[k] = -g_x_y_xy_xxzz[k] * ab_z + g_x_y_xy_xxzzz[k];

                g_x_y_xyz_xyyy[k] = -g_x_y_xy_xyyy[k] * ab_z + g_x_y_xy_xyyyz[k];

                g_x_y_xyz_xyyz[k] = -g_x_y_xy_xyyz[k] * ab_z + g_x_y_xy_xyyzz[k];

                g_x_y_xyz_xyzz[k] = -g_x_y_xy_xyzz[k] * ab_z + g_x_y_xy_xyzzz[k];

                g_x_y_xyz_xzzz[k] = -g_x_y_xy_xzzz[k] * ab_z + g_x_y_xy_xzzzz[k];

                g_x_y_xyz_yyyy[k] = -g_x_y_xy_yyyy[k] * ab_z + g_x_y_xy_yyyyz[k];

                g_x_y_xyz_yyyz[k] = -g_x_y_xy_yyyz[k] * ab_z + g_x_y_xy_yyyzz[k];

                g_x_y_xyz_yyzz[k] = -g_x_y_xy_yyzz[k] * ab_z + g_x_y_xy_yyzzz[k];

                g_x_y_xyz_yzzz[k] = -g_x_y_xy_yzzz[k] * ab_z + g_x_y_xy_yzzzz[k];

                g_x_y_xyz_zzzz[k] = -g_x_y_xy_zzzz[k] * ab_z + g_x_y_xy_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzz_xxxx = cbuffer.data(fg_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_xzz_xxxy = cbuffer.data(fg_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_xzz_xxxz = cbuffer.data(fg_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_y_xzz_xxyy = cbuffer.data(fg_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_xzz_xxyz = cbuffer.data(fg_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_xzz_xxzz = cbuffer.data(fg_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_y_xzz_xyyy = cbuffer.data(fg_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_xzz_xyyz = cbuffer.data(fg_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_xzz_xyzz = cbuffer.data(fg_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_y_xzz_xzzz = cbuffer.data(fg_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_xzz_yyyy = cbuffer.data(fg_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_xzz_yyyz = cbuffer.data(fg_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_xzz_yyzz = cbuffer.data(fg_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_xzz_yzzz = cbuffer.data(fg_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_xzz_zzzz = cbuffer.data(fg_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xz_xxxx, g_x_y_xz_xxxxz, g_x_y_xz_xxxy, g_x_y_xz_xxxyz, g_x_y_xz_xxxz, g_x_y_xz_xxxzz, g_x_y_xz_xxyy, g_x_y_xz_xxyyz, g_x_y_xz_xxyz, g_x_y_xz_xxyzz, g_x_y_xz_xxzz, g_x_y_xz_xxzzz, g_x_y_xz_xyyy, g_x_y_xz_xyyyz, g_x_y_xz_xyyz, g_x_y_xz_xyyzz, g_x_y_xz_xyzz, g_x_y_xz_xyzzz, g_x_y_xz_xzzz, g_x_y_xz_xzzzz, g_x_y_xz_yyyy, g_x_y_xz_yyyyz, g_x_y_xz_yyyz, g_x_y_xz_yyyzz, g_x_y_xz_yyzz, g_x_y_xz_yyzzz, g_x_y_xz_yzzz, g_x_y_xz_yzzzz, g_x_y_xz_zzzz, g_x_y_xz_zzzzz, g_x_y_xzz_xxxx, g_x_y_xzz_xxxy, g_x_y_xzz_xxxz, g_x_y_xzz_xxyy, g_x_y_xzz_xxyz, g_x_y_xzz_xxzz, g_x_y_xzz_xyyy, g_x_y_xzz_xyyz, g_x_y_xzz_xyzz, g_x_y_xzz_xzzz, g_x_y_xzz_yyyy, g_x_y_xzz_yyyz, g_x_y_xzz_yyzz, g_x_y_xzz_yzzz, g_x_y_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzz_xxxx[k] = -g_x_y_xz_xxxx[k] * ab_z + g_x_y_xz_xxxxz[k];

                g_x_y_xzz_xxxy[k] = -g_x_y_xz_xxxy[k] * ab_z + g_x_y_xz_xxxyz[k];

                g_x_y_xzz_xxxz[k] = -g_x_y_xz_xxxz[k] * ab_z + g_x_y_xz_xxxzz[k];

                g_x_y_xzz_xxyy[k] = -g_x_y_xz_xxyy[k] * ab_z + g_x_y_xz_xxyyz[k];

                g_x_y_xzz_xxyz[k] = -g_x_y_xz_xxyz[k] * ab_z + g_x_y_xz_xxyzz[k];

                g_x_y_xzz_xxzz[k] = -g_x_y_xz_xxzz[k] * ab_z + g_x_y_xz_xxzzz[k];

                g_x_y_xzz_xyyy[k] = -g_x_y_xz_xyyy[k] * ab_z + g_x_y_xz_xyyyz[k];

                g_x_y_xzz_xyyz[k] = -g_x_y_xz_xyyz[k] * ab_z + g_x_y_xz_xyyzz[k];

                g_x_y_xzz_xyzz[k] = -g_x_y_xz_xyzz[k] * ab_z + g_x_y_xz_xyzzz[k];

                g_x_y_xzz_xzzz[k] = -g_x_y_xz_xzzz[k] * ab_z + g_x_y_xz_xzzzz[k];

                g_x_y_xzz_yyyy[k] = -g_x_y_xz_yyyy[k] * ab_z + g_x_y_xz_yyyyz[k];

                g_x_y_xzz_yyyz[k] = -g_x_y_xz_yyyz[k] * ab_z + g_x_y_xz_yyyzz[k];

                g_x_y_xzz_yyzz[k] = -g_x_y_xz_yyzz[k] * ab_z + g_x_y_xz_yyzzz[k];

                g_x_y_xzz_yzzz[k] = -g_x_y_xz_yzzz[k] * ab_z + g_x_y_xz_yzzzz[k];

                g_x_y_xzz_zzzz[k] = -g_x_y_xz_zzzz[k] * ab_z + g_x_y_xz_zzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyy_xxxx = cbuffer.data(fg_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_yyy_xxxy = cbuffer.data(fg_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_yyy_xxxz = cbuffer.data(fg_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_yyy_xxyy = cbuffer.data(fg_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_yyy_xxyz = cbuffer.data(fg_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_yyy_xxzz = cbuffer.data(fg_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_y_yyy_xyyy = cbuffer.data(fg_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_yyy_xyyz = cbuffer.data(fg_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_yyy_xyzz = cbuffer.data(fg_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_yyy_xzzz = cbuffer.data(fg_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_yyy_yyyy = cbuffer.data(fg_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_yyy_yyyz = cbuffer.data(fg_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_y_yyy_yyzz = cbuffer.data(fg_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_y_yyy_yzzz = cbuffer.data(fg_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_y_yyy_zzzz = cbuffer.data(fg_geom_11_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xxxx, g_x_0_yy_xxxy, g_x_0_yy_xxxz, g_x_0_yy_xxyy, g_x_0_yy_xxyz, g_x_0_yy_xxzz, g_x_0_yy_xyyy, g_x_0_yy_xyyz, g_x_0_yy_xyzz, g_x_0_yy_xzzz, g_x_0_yy_yyyy, g_x_0_yy_yyyz, g_x_0_yy_yyzz, g_x_0_yy_yzzz, g_x_0_yy_zzzz, g_x_y_yy_xxxx, g_x_y_yy_xxxxy, g_x_y_yy_xxxy, g_x_y_yy_xxxyy, g_x_y_yy_xxxyz, g_x_y_yy_xxxz, g_x_y_yy_xxyy, g_x_y_yy_xxyyy, g_x_y_yy_xxyyz, g_x_y_yy_xxyz, g_x_y_yy_xxyzz, g_x_y_yy_xxzz, g_x_y_yy_xyyy, g_x_y_yy_xyyyy, g_x_y_yy_xyyyz, g_x_y_yy_xyyz, g_x_y_yy_xyyzz, g_x_y_yy_xyzz, g_x_y_yy_xyzzz, g_x_y_yy_xzzz, g_x_y_yy_yyyy, g_x_y_yy_yyyyy, g_x_y_yy_yyyyz, g_x_y_yy_yyyz, g_x_y_yy_yyyzz, g_x_y_yy_yyzz, g_x_y_yy_yyzzz, g_x_y_yy_yzzz, g_x_y_yy_yzzzz, g_x_y_yy_zzzz, g_x_y_yyy_xxxx, g_x_y_yyy_xxxy, g_x_y_yyy_xxxz, g_x_y_yyy_xxyy, g_x_y_yyy_xxyz, g_x_y_yyy_xxzz, g_x_y_yyy_xyyy, g_x_y_yyy_xyyz, g_x_y_yyy_xyzz, g_x_y_yyy_xzzz, g_x_y_yyy_yyyy, g_x_y_yyy_yyyz, g_x_y_yyy_yyzz, g_x_y_yyy_yzzz, g_x_y_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyy_xxxx[k] = g_x_0_yy_xxxx[k] - g_x_y_yy_xxxx[k] * ab_y + g_x_y_yy_xxxxy[k];

                g_x_y_yyy_xxxy[k] = g_x_0_yy_xxxy[k] - g_x_y_yy_xxxy[k] * ab_y + g_x_y_yy_xxxyy[k];

                g_x_y_yyy_xxxz[k] = g_x_0_yy_xxxz[k] - g_x_y_yy_xxxz[k] * ab_y + g_x_y_yy_xxxyz[k];

                g_x_y_yyy_xxyy[k] = g_x_0_yy_xxyy[k] - g_x_y_yy_xxyy[k] * ab_y + g_x_y_yy_xxyyy[k];

                g_x_y_yyy_xxyz[k] = g_x_0_yy_xxyz[k] - g_x_y_yy_xxyz[k] * ab_y + g_x_y_yy_xxyyz[k];

                g_x_y_yyy_xxzz[k] = g_x_0_yy_xxzz[k] - g_x_y_yy_xxzz[k] * ab_y + g_x_y_yy_xxyzz[k];

                g_x_y_yyy_xyyy[k] = g_x_0_yy_xyyy[k] - g_x_y_yy_xyyy[k] * ab_y + g_x_y_yy_xyyyy[k];

                g_x_y_yyy_xyyz[k] = g_x_0_yy_xyyz[k] - g_x_y_yy_xyyz[k] * ab_y + g_x_y_yy_xyyyz[k];

                g_x_y_yyy_xyzz[k] = g_x_0_yy_xyzz[k] - g_x_y_yy_xyzz[k] * ab_y + g_x_y_yy_xyyzz[k];

                g_x_y_yyy_xzzz[k] = g_x_0_yy_xzzz[k] - g_x_y_yy_xzzz[k] * ab_y + g_x_y_yy_xyzzz[k];

                g_x_y_yyy_yyyy[k] = g_x_0_yy_yyyy[k] - g_x_y_yy_yyyy[k] * ab_y + g_x_y_yy_yyyyy[k];

                g_x_y_yyy_yyyz[k] = g_x_0_yy_yyyz[k] - g_x_y_yy_yyyz[k] * ab_y + g_x_y_yy_yyyyz[k];

                g_x_y_yyy_yyzz[k] = g_x_0_yy_yyzz[k] - g_x_y_yy_yyzz[k] * ab_y + g_x_y_yy_yyyzz[k];

                g_x_y_yyy_yzzz[k] = g_x_0_yy_yzzz[k] - g_x_y_yy_yzzz[k] * ab_y + g_x_y_yy_yyzzz[k];

                g_x_y_yyy_zzzz[k] = g_x_0_yy_zzzz[k] - g_x_y_yy_zzzz[k] * ab_y + g_x_y_yy_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyz_xxxx = cbuffer.data(fg_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_y_yyz_xxxy = cbuffer.data(fg_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_y_yyz_xxxz = cbuffer.data(fg_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_y_yyz_xxyy = cbuffer.data(fg_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_y_yyz_xxyz = cbuffer.data(fg_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_y_yyz_xxzz = cbuffer.data(fg_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_y_yyz_xyyy = cbuffer.data(fg_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_y_yyz_xyyz = cbuffer.data(fg_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_y_yyz_xyzz = cbuffer.data(fg_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_y_yyz_xzzz = cbuffer.data(fg_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_y_yyz_yyyy = cbuffer.data(fg_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_y_yyz_yyyz = cbuffer.data(fg_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_y_yyz_yyzz = cbuffer.data(fg_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_y_yyz_yzzz = cbuffer.data(fg_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_y_yyz_zzzz = cbuffer.data(fg_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yy_xxxx, g_x_y_yy_xxxxz, g_x_y_yy_xxxy, g_x_y_yy_xxxyz, g_x_y_yy_xxxz, g_x_y_yy_xxxzz, g_x_y_yy_xxyy, g_x_y_yy_xxyyz, g_x_y_yy_xxyz, g_x_y_yy_xxyzz, g_x_y_yy_xxzz, g_x_y_yy_xxzzz, g_x_y_yy_xyyy, g_x_y_yy_xyyyz, g_x_y_yy_xyyz, g_x_y_yy_xyyzz, g_x_y_yy_xyzz, g_x_y_yy_xyzzz, g_x_y_yy_xzzz, g_x_y_yy_xzzzz, g_x_y_yy_yyyy, g_x_y_yy_yyyyz, g_x_y_yy_yyyz, g_x_y_yy_yyyzz, g_x_y_yy_yyzz, g_x_y_yy_yyzzz, g_x_y_yy_yzzz, g_x_y_yy_yzzzz, g_x_y_yy_zzzz, g_x_y_yy_zzzzz, g_x_y_yyz_xxxx, g_x_y_yyz_xxxy, g_x_y_yyz_xxxz, g_x_y_yyz_xxyy, g_x_y_yyz_xxyz, g_x_y_yyz_xxzz, g_x_y_yyz_xyyy, g_x_y_yyz_xyyz, g_x_y_yyz_xyzz, g_x_y_yyz_xzzz, g_x_y_yyz_yyyy, g_x_y_yyz_yyyz, g_x_y_yyz_yyzz, g_x_y_yyz_yzzz, g_x_y_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyz_xxxx[k] = -g_x_y_yy_xxxx[k] * ab_z + g_x_y_yy_xxxxz[k];

                g_x_y_yyz_xxxy[k] = -g_x_y_yy_xxxy[k] * ab_z + g_x_y_yy_xxxyz[k];

                g_x_y_yyz_xxxz[k] = -g_x_y_yy_xxxz[k] * ab_z + g_x_y_yy_xxxzz[k];

                g_x_y_yyz_xxyy[k] = -g_x_y_yy_xxyy[k] * ab_z + g_x_y_yy_xxyyz[k];

                g_x_y_yyz_xxyz[k] = -g_x_y_yy_xxyz[k] * ab_z + g_x_y_yy_xxyzz[k];

                g_x_y_yyz_xxzz[k] = -g_x_y_yy_xxzz[k] * ab_z + g_x_y_yy_xxzzz[k];

                g_x_y_yyz_xyyy[k] = -g_x_y_yy_xyyy[k] * ab_z + g_x_y_yy_xyyyz[k];

                g_x_y_yyz_xyyz[k] = -g_x_y_yy_xyyz[k] * ab_z + g_x_y_yy_xyyzz[k];

                g_x_y_yyz_xyzz[k] = -g_x_y_yy_xyzz[k] * ab_z + g_x_y_yy_xyzzz[k];

                g_x_y_yyz_xzzz[k] = -g_x_y_yy_xzzz[k] * ab_z + g_x_y_yy_xzzzz[k];

                g_x_y_yyz_yyyy[k] = -g_x_y_yy_yyyy[k] * ab_z + g_x_y_yy_yyyyz[k];

                g_x_y_yyz_yyyz[k] = -g_x_y_yy_yyyz[k] * ab_z + g_x_y_yy_yyyzz[k];

                g_x_y_yyz_yyzz[k] = -g_x_y_yy_yyzz[k] * ab_z + g_x_y_yy_yyzzz[k];

                g_x_y_yyz_yzzz[k] = -g_x_y_yy_yzzz[k] * ab_z + g_x_y_yy_yzzzz[k];

                g_x_y_yyz_zzzz[k] = -g_x_y_yy_zzzz[k] * ab_z + g_x_y_yy_zzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzz_xxxx = cbuffer.data(fg_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_y_yzz_xxxy = cbuffer.data(fg_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_y_yzz_xxxz = cbuffer.data(fg_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_y_yzz_xxyy = cbuffer.data(fg_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_y_yzz_xxyz = cbuffer.data(fg_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_y_yzz_xxzz = cbuffer.data(fg_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_y_yzz_xyyy = cbuffer.data(fg_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_y_yzz_xyyz = cbuffer.data(fg_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_y_yzz_xyzz = cbuffer.data(fg_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_y_yzz_xzzz = cbuffer.data(fg_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_y_yzz_yyyy = cbuffer.data(fg_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_y_yzz_yyyz = cbuffer.data(fg_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_y_yzz_yyzz = cbuffer.data(fg_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_y_yzz_yzzz = cbuffer.data(fg_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_y_yzz_zzzz = cbuffer.data(fg_geom_11_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yz_xxxx, g_x_y_yz_xxxxz, g_x_y_yz_xxxy, g_x_y_yz_xxxyz, g_x_y_yz_xxxz, g_x_y_yz_xxxzz, g_x_y_yz_xxyy, g_x_y_yz_xxyyz, g_x_y_yz_xxyz, g_x_y_yz_xxyzz, g_x_y_yz_xxzz, g_x_y_yz_xxzzz, g_x_y_yz_xyyy, g_x_y_yz_xyyyz, g_x_y_yz_xyyz, g_x_y_yz_xyyzz, g_x_y_yz_xyzz, g_x_y_yz_xyzzz, g_x_y_yz_xzzz, g_x_y_yz_xzzzz, g_x_y_yz_yyyy, g_x_y_yz_yyyyz, g_x_y_yz_yyyz, g_x_y_yz_yyyzz, g_x_y_yz_yyzz, g_x_y_yz_yyzzz, g_x_y_yz_yzzz, g_x_y_yz_yzzzz, g_x_y_yz_zzzz, g_x_y_yz_zzzzz, g_x_y_yzz_xxxx, g_x_y_yzz_xxxy, g_x_y_yzz_xxxz, g_x_y_yzz_xxyy, g_x_y_yzz_xxyz, g_x_y_yzz_xxzz, g_x_y_yzz_xyyy, g_x_y_yzz_xyyz, g_x_y_yzz_xyzz, g_x_y_yzz_xzzz, g_x_y_yzz_yyyy, g_x_y_yzz_yyyz, g_x_y_yzz_yyzz, g_x_y_yzz_yzzz, g_x_y_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzz_xxxx[k] = -g_x_y_yz_xxxx[k] * ab_z + g_x_y_yz_xxxxz[k];

                g_x_y_yzz_xxxy[k] = -g_x_y_yz_xxxy[k] * ab_z + g_x_y_yz_xxxyz[k];

                g_x_y_yzz_xxxz[k] = -g_x_y_yz_xxxz[k] * ab_z + g_x_y_yz_xxxzz[k];

                g_x_y_yzz_xxyy[k] = -g_x_y_yz_xxyy[k] * ab_z + g_x_y_yz_xxyyz[k];

                g_x_y_yzz_xxyz[k] = -g_x_y_yz_xxyz[k] * ab_z + g_x_y_yz_xxyzz[k];

                g_x_y_yzz_xxzz[k] = -g_x_y_yz_xxzz[k] * ab_z + g_x_y_yz_xxzzz[k];

                g_x_y_yzz_xyyy[k] = -g_x_y_yz_xyyy[k] * ab_z + g_x_y_yz_xyyyz[k];

                g_x_y_yzz_xyyz[k] = -g_x_y_yz_xyyz[k] * ab_z + g_x_y_yz_xyyzz[k];

                g_x_y_yzz_xyzz[k] = -g_x_y_yz_xyzz[k] * ab_z + g_x_y_yz_xyzzz[k];

                g_x_y_yzz_xzzz[k] = -g_x_y_yz_xzzz[k] * ab_z + g_x_y_yz_xzzzz[k];

                g_x_y_yzz_yyyy[k] = -g_x_y_yz_yyyy[k] * ab_z + g_x_y_yz_yyyyz[k];

                g_x_y_yzz_yyyz[k] = -g_x_y_yz_yyyz[k] * ab_z + g_x_y_yz_yyyzz[k];

                g_x_y_yzz_yyzz[k] = -g_x_y_yz_yyzz[k] * ab_z + g_x_y_yz_yyzzz[k];

                g_x_y_yzz_yzzz[k] = -g_x_y_yz_yzzz[k] * ab_z + g_x_y_yz_yzzzz[k];

                g_x_y_yzz_zzzz[k] = -g_x_y_yz_zzzz[k] * ab_z + g_x_y_yz_zzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzz_xxxx = cbuffer.data(fg_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_y_zzz_xxxy = cbuffer.data(fg_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_y_zzz_xxxz = cbuffer.data(fg_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_y_zzz_xxyy = cbuffer.data(fg_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_y_zzz_xxyz = cbuffer.data(fg_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_y_zzz_xxzz = cbuffer.data(fg_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_y_zzz_xyyy = cbuffer.data(fg_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_y_zzz_xyyz = cbuffer.data(fg_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_y_zzz_xyzz = cbuffer.data(fg_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_y_zzz_xzzz = cbuffer.data(fg_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_y_zzz_yyyy = cbuffer.data(fg_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_y_zzz_yyyz = cbuffer.data(fg_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_y_zzz_yyzz = cbuffer.data(fg_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_y_zzz_yzzz = cbuffer.data(fg_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_y_zzz_zzzz = cbuffer.data(fg_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zz_xxxx, g_x_y_zz_xxxxz, g_x_y_zz_xxxy, g_x_y_zz_xxxyz, g_x_y_zz_xxxz, g_x_y_zz_xxxzz, g_x_y_zz_xxyy, g_x_y_zz_xxyyz, g_x_y_zz_xxyz, g_x_y_zz_xxyzz, g_x_y_zz_xxzz, g_x_y_zz_xxzzz, g_x_y_zz_xyyy, g_x_y_zz_xyyyz, g_x_y_zz_xyyz, g_x_y_zz_xyyzz, g_x_y_zz_xyzz, g_x_y_zz_xyzzz, g_x_y_zz_xzzz, g_x_y_zz_xzzzz, g_x_y_zz_yyyy, g_x_y_zz_yyyyz, g_x_y_zz_yyyz, g_x_y_zz_yyyzz, g_x_y_zz_yyzz, g_x_y_zz_yyzzz, g_x_y_zz_yzzz, g_x_y_zz_yzzzz, g_x_y_zz_zzzz, g_x_y_zz_zzzzz, g_x_y_zzz_xxxx, g_x_y_zzz_xxxy, g_x_y_zzz_xxxz, g_x_y_zzz_xxyy, g_x_y_zzz_xxyz, g_x_y_zzz_xxzz, g_x_y_zzz_xyyy, g_x_y_zzz_xyyz, g_x_y_zzz_xyzz, g_x_y_zzz_xzzz, g_x_y_zzz_yyyy, g_x_y_zzz_yyyz, g_x_y_zzz_yyzz, g_x_y_zzz_yzzz, g_x_y_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzz_xxxx[k] = -g_x_y_zz_xxxx[k] * ab_z + g_x_y_zz_xxxxz[k];

                g_x_y_zzz_xxxy[k] = -g_x_y_zz_xxxy[k] * ab_z + g_x_y_zz_xxxyz[k];

                g_x_y_zzz_xxxz[k] = -g_x_y_zz_xxxz[k] * ab_z + g_x_y_zz_xxxzz[k];

                g_x_y_zzz_xxyy[k] = -g_x_y_zz_xxyy[k] * ab_z + g_x_y_zz_xxyyz[k];

                g_x_y_zzz_xxyz[k] = -g_x_y_zz_xxyz[k] * ab_z + g_x_y_zz_xxyzz[k];

                g_x_y_zzz_xxzz[k] = -g_x_y_zz_xxzz[k] * ab_z + g_x_y_zz_xxzzz[k];

                g_x_y_zzz_xyyy[k] = -g_x_y_zz_xyyy[k] * ab_z + g_x_y_zz_xyyyz[k];

                g_x_y_zzz_xyyz[k] = -g_x_y_zz_xyyz[k] * ab_z + g_x_y_zz_xyyzz[k];

                g_x_y_zzz_xyzz[k] = -g_x_y_zz_xyzz[k] * ab_z + g_x_y_zz_xyzzz[k];

                g_x_y_zzz_xzzz[k] = -g_x_y_zz_xzzz[k] * ab_z + g_x_y_zz_xzzzz[k];

                g_x_y_zzz_yyyy[k] = -g_x_y_zz_yyyy[k] * ab_z + g_x_y_zz_yyyyz[k];

                g_x_y_zzz_yyyz[k] = -g_x_y_zz_yyyz[k] * ab_z + g_x_y_zz_yyyzz[k];

                g_x_y_zzz_yyzz[k] = -g_x_y_zz_yyzz[k] * ab_z + g_x_y_zz_yyzzz[k];

                g_x_y_zzz_yzzz[k] = -g_x_y_zz_yzzz[k] * ab_z + g_x_y_zz_yzzzz[k];

                g_x_y_zzz_zzzz[k] = -g_x_y_zz_zzzz[k] * ab_z + g_x_y_zz_zzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxx_xxxx = cbuffer.data(fg_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_xxx_xxxy = cbuffer.data(fg_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_xxx_xxxz = cbuffer.data(fg_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_xxx_xxyy = cbuffer.data(fg_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_xxx_xxyz = cbuffer.data(fg_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_xxx_xxzz = cbuffer.data(fg_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_z_xxx_xyyy = cbuffer.data(fg_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_xxx_xyyz = cbuffer.data(fg_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_xxx_xyzz = cbuffer.data(fg_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_xxx_xzzz = cbuffer.data(fg_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_xxx_yyyy = cbuffer.data(fg_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_xxx_yyyz = cbuffer.data(fg_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_z_xxx_yyzz = cbuffer.data(fg_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_xxx_yzzz = cbuffer.data(fg_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_xxx_zzzz = cbuffer.data(fg_geom_11_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_xxxx, g_0_z_xx_xxxy, g_0_z_xx_xxxz, g_0_z_xx_xxyy, g_0_z_xx_xxyz, g_0_z_xx_xxzz, g_0_z_xx_xyyy, g_0_z_xx_xyyz, g_0_z_xx_xyzz, g_0_z_xx_xzzz, g_0_z_xx_yyyy, g_0_z_xx_yyyz, g_0_z_xx_yyzz, g_0_z_xx_yzzz, g_0_z_xx_zzzz, g_x_z_xx_xxxx, g_x_z_xx_xxxxx, g_x_z_xx_xxxxy, g_x_z_xx_xxxxz, g_x_z_xx_xxxy, g_x_z_xx_xxxyy, g_x_z_xx_xxxyz, g_x_z_xx_xxxz, g_x_z_xx_xxxzz, g_x_z_xx_xxyy, g_x_z_xx_xxyyy, g_x_z_xx_xxyyz, g_x_z_xx_xxyz, g_x_z_xx_xxyzz, g_x_z_xx_xxzz, g_x_z_xx_xxzzz, g_x_z_xx_xyyy, g_x_z_xx_xyyyy, g_x_z_xx_xyyyz, g_x_z_xx_xyyz, g_x_z_xx_xyyzz, g_x_z_xx_xyzz, g_x_z_xx_xyzzz, g_x_z_xx_xzzz, g_x_z_xx_xzzzz, g_x_z_xx_yyyy, g_x_z_xx_yyyz, g_x_z_xx_yyzz, g_x_z_xx_yzzz, g_x_z_xx_zzzz, g_x_z_xxx_xxxx, g_x_z_xxx_xxxy, g_x_z_xxx_xxxz, g_x_z_xxx_xxyy, g_x_z_xxx_xxyz, g_x_z_xxx_xxzz, g_x_z_xxx_xyyy, g_x_z_xxx_xyyz, g_x_z_xxx_xyzz, g_x_z_xxx_xzzz, g_x_z_xxx_yyyy, g_x_z_xxx_yyyz, g_x_z_xxx_yyzz, g_x_z_xxx_yzzz, g_x_z_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxx_xxxx[k] = -g_0_z_xx_xxxx[k] - g_x_z_xx_xxxx[k] * ab_x + g_x_z_xx_xxxxx[k];

                g_x_z_xxx_xxxy[k] = -g_0_z_xx_xxxy[k] - g_x_z_xx_xxxy[k] * ab_x + g_x_z_xx_xxxxy[k];

                g_x_z_xxx_xxxz[k] = -g_0_z_xx_xxxz[k] - g_x_z_xx_xxxz[k] * ab_x + g_x_z_xx_xxxxz[k];

                g_x_z_xxx_xxyy[k] = -g_0_z_xx_xxyy[k] - g_x_z_xx_xxyy[k] * ab_x + g_x_z_xx_xxxyy[k];

                g_x_z_xxx_xxyz[k] = -g_0_z_xx_xxyz[k] - g_x_z_xx_xxyz[k] * ab_x + g_x_z_xx_xxxyz[k];

                g_x_z_xxx_xxzz[k] = -g_0_z_xx_xxzz[k] - g_x_z_xx_xxzz[k] * ab_x + g_x_z_xx_xxxzz[k];

                g_x_z_xxx_xyyy[k] = -g_0_z_xx_xyyy[k] - g_x_z_xx_xyyy[k] * ab_x + g_x_z_xx_xxyyy[k];

                g_x_z_xxx_xyyz[k] = -g_0_z_xx_xyyz[k] - g_x_z_xx_xyyz[k] * ab_x + g_x_z_xx_xxyyz[k];

                g_x_z_xxx_xyzz[k] = -g_0_z_xx_xyzz[k] - g_x_z_xx_xyzz[k] * ab_x + g_x_z_xx_xxyzz[k];

                g_x_z_xxx_xzzz[k] = -g_0_z_xx_xzzz[k] - g_x_z_xx_xzzz[k] * ab_x + g_x_z_xx_xxzzz[k];

                g_x_z_xxx_yyyy[k] = -g_0_z_xx_yyyy[k] - g_x_z_xx_yyyy[k] * ab_x + g_x_z_xx_xyyyy[k];

                g_x_z_xxx_yyyz[k] = -g_0_z_xx_yyyz[k] - g_x_z_xx_yyyz[k] * ab_x + g_x_z_xx_xyyyz[k];

                g_x_z_xxx_yyzz[k] = -g_0_z_xx_yyzz[k] - g_x_z_xx_yyzz[k] * ab_x + g_x_z_xx_xyyzz[k];

                g_x_z_xxx_yzzz[k] = -g_0_z_xx_yzzz[k] - g_x_z_xx_yzzz[k] * ab_x + g_x_z_xx_xyzzz[k];

                g_x_z_xxx_zzzz[k] = -g_0_z_xx_zzzz[k] - g_x_z_xx_zzzz[k] * ab_x + g_x_z_xx_xzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxy_xxxx = cbuffer.data(fg_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_xxy_xxxy = cbuffer.data(fg_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_xxy_xxxz = cbuffer.data(fg_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_z_xxy_xxyy = cbuffer.data(fg_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_xxy_xxyz = cbuffer.data(fg_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_xxy_xxzz = cbuffer.data(fg_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_xxy_xyyy = cbuffer.data(fg_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_xxy_xyyz = cbuffer.data(fg_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_xxy_xyzz = cbuffer.data(fg_geom_11_off + 323 * ccomps * dcomps);

            auto g_x_z_xxy_xzzz = cbuffer.data(fg_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_z_xxy_yyyy = cbuffer.data(fg_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_z_xxy_yyyz = cbuffer.data(fg_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_z_xxy_yyzz = cbuffer.data(fg_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_z_xxy_yzzz = cbuffer.data(fg_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_z_xxy_zzzz = cbuffer.data(fg_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xx_xxxx, g_x_z_xx_xxxxy, g_x_z_xx_xxxy, g_x_z_xx_xxxyy, g_x_z_xx_xxxyz, g_x_z_xx_xxxz, g_x_z_xx_xxyy, g_x_z_xx_xxyyy, g_x_z_xx_xxyyz, g_x_z_xx_xxyz, g_x_z_xx_xxyzz, g_x_z_xx_xxzz, g_x_z_xx_xyyy, g_x_z_xx_xyyyy, g_x_z_xx_xyyyz, g_x_z_xx_xyyz, g_x_z_xx_xyyzz, g_x_z_xx_xyzz, g_x_z_xx_xyzzz, g_x_z_xx_xzzz, g_x_z_xx_yyyy, g_x_z_xx_yyyyy, g_x_z_xx_yyyyz, g_x_z_xx_yyyz, g_x_z_xx_yyyzz, g_x_z_xx_yyzz, g_x_z_xx_yyzzz, g_x_z_xx_yzzz, g_x_z_xx_yzzzz, g_x_z_xx_zzzz, g_x_z_xxy_xxxx, g_x_z_xxy_xxxy, g_x_z_xxy_xxxz, g_x_z_xxy_xxyy, g_x_z_xxy_xxyz, g_x_z_xxy_xxzz, g_x_z_xxy_xyyy, g_x_z_xxy_xyyz, g_x_z_xxy_xyzz, g_x_z_xxy_xzzz, g_x_z_xxy_yyyy, g_x_z_xxy_yyyz, g_x_z_xxy_yyzz, g_x_z_xxy_yzzz, g_x_z_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxy_xxxx[k] = -g_x_z_xx_xxxx[k] * ab_y + g_x_z_xx_xxxxy[k];

                g_x_z_xxy_xxxy[k] = -g_x_z_xx_xxxy[k] * ab_y + g_x_z_xx_xxxyy[k];

                g_x_z_xxy_xxxz[k] = -g_x_z_xx_xxxz[k] * ab_y + g_x_z_xx_xxxyz[k];

                g_x_z_xxy_xxyy[k] = -g_x_z_xx_xxyy[k] * ab_y + g_x_z_xx_xxyyy[k];

                g_x_z_xxy_xxyz[k] = -g_x_z_xx_xxyz[k] * ab_y + g_x_z_xx_xxyyz[k];

                g_x_z_xxy_xxzz[k] = -g_x_z_xx_xxzz[k] * ab_y + g_x_z_xx_xxyzz[k];

                g_x_z_xxy_xyyy[k] = -g_x_z_xx_xyyy[k] * ab_y + g_x_z_xx_xyyyy[k];

                g_x_z_xxy_xyyz[k] = -g_x_z_xx_xyyz[k] * ab_y + g_x_z_xx_xyyyz[k];

                g_x_z_xxy_xyzz[k] = -g_x_z_xx_xyzz[k] * ab_y + g_x_z_xx_xyyzz[k];

                g_x_z_xxy_xzzz[k] = -g_x_z_xx_xzzz[k] * ab_y + g_x_z_xx_xyzzz[k];

                g_x_z_xxy_yyyy[k] = -g_x_z_xx_yyyy[k] * ab_y + g_x_z_xx_yyyyy[k];

                g_x_z_xxy_yyyz[k] = -g_x_z_xx_yyyz[k] * ab_y + g_x_z_xx_yyyyz[k];

                g_x_z_xxy_yyzz[k] = -g_x_z_xx_yyzz[k] * ab_y + g_x_z_xx_yyyzz[k];

                g_x_z_xxy_yzzz[k] = -g_x_z_xx_yzzz[k] * ab_y + g_x_z_xx_yyzzz[k];

                g_x_z_xxy_zzzz[k] = -g_x_z_xx_zzzz[k] * ab_y + g_x_z_xx_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxz_xxxx = cbuffer.data(fg_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_z_xxz_xxxy = cbuffer.data(fg_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_z_xxz_xxxz = cbuffer.data(fg_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_z_xxz_xxyy = cbuffer.data(fg_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_z_xxz_xxyz = cbuffer.data(fg_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_z_xxz_xxzz = cbuffer.data(fg_geom_11_off + 335 * ccomps * dcomps);

            auto g_x_z_xxz_xyyy = cbuffer.data(fg_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_xxz_xyyz = cbuffer.data(fg_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_xxz_xyzz = cbuffer.data(fg_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_xxz_xzzz = cbuffer.data(fg_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_xxz_yyyy = cbuffer.data(fg_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_xxz_yyyz = cbuffer.data(fg_geom_11_off + 341 * ccomps * dcomps);

            auto g_x_z_xxz_yyzz = cbuffer.data(fg_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_xxz_yzzz = cbuffer.data(fg_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_xxz_zzzz = cbuffer.data(fg_geom_11_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_xxxx, g_0_z_xz_xxxy, g_0_z_xz_xxxz, g_0_z_xz_xxyy, g_0_z_xz_xxyz, g_0_z_xz_xxzz, g_0_z_xz_xyyy, g_0_z_xz_xyyz, g_0_z_xz_xyzz, g_0_z_xz_xzzz, g_0_z_xz_yyyy, g_0_z_xz_yyyz, g_0_z_xz_yyzz, g_0_z_xz_yzzz, g_0_z_xz_zzzz, g_x_z_xxz_xxxx, g_x_z_xxz_xxxy, g_x_z_xxz_xxxz, g_x_z_xxz_xxyy, g_x_z_xxz_xxyz, g_x_z_xxz_xxzz, g_x_z_xxz_xyyy, g_x_z_xxz_xyyz, g_x_z_xxz_xyzz, g_x_z_xxz_xzzz, g_x_z_xxz_yyyy, g_x_z_xxz_yyyz, g_x_z_xxz_yyzz, g_x_z_xxz_yzzz, g_x_z_xxz_zzzz, g_x_z_xz_xxxx, g_x_z_xz_xxxxx, g_x_z_xz_xxxxy, g_x_z_xz_xxxxz, g_x_z_xz_xxxy, g_x_z_xz_xxxyy, g_x_z_xz_xxxyz, g_x_z_xz_xxxz, g_x_z_xz_xxxzz, g_x_z_xz_xxyy, g_x_z_xz_xxyyy, g_x_z_xz_xxyyz, g_x_z_xz_xxyz, g_x_z_xz_xxyzz, g_x_z_xz_xxzz, g_x_z_xz_xxzzz, g_x_z_xz_xyyy, g_x_z_xz_xyyyy, g_x_z_xz_xyyyz, g_x_z_xz_xyyz, g_x_z_xz_xyyzz, g_x_z_xz_xyzz, g_x_z_xz_xyzzz, g_x_z_xz_xzzz, g_x_z_xz_xzzzz, g_x_z_xz_yyyy, g_x_z_xz_yyyz, g_x_z_xz_yyzz, g_x_z_xz_yzzz, g_x_z_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxz_xxxx[k] = -g_0_z_xz_xxxx[k] - g_x_z_xz_xxxx[k] * ab_x + g_x_z_xz_xxxxx[k];

                g_x_z_xxz_xxxy[k] = -g_0_z_xz_xxxy[k] - g_x_z_xz_xxxy[k] * ab_x + g_x_z_xz_xxxxy[k];

                g_x_z_xxz_xxxz[k] = -g_0_z_xz_xxxz[k] - g_x_z_xz_xxxz[k] * ab_x + g_x_z_xz_xxxxz[k];

                g_x_z_xxz_xxyy[k] = -g_0_z_xz_xxyy[k] - g_x_z_xz_xxyy[k] * ab_x + g_x_z_xz_xxxyy[k];

                g_x_z_xxz_xxyz[k] = -g_0_z_xz_xxyz[k] - g_x_z_xz_xxyz[k] * ab_x + g_x_z_xz_xxxyz[k];

                g_x_z_xxz_xxzz[k] = -g_0_z_xz_xxzz[k] - g_x_z_xz_xxzz[k] * ab_x + g_x_z_xz_xxxzz[k];

                g_x_z_xxz_xyyy[k] = -g_0_z_xz_xyyy[k] - g_x_z_xz_xyyy[k] * ab_x + g_x_z_xz_xxyyy[k];

                g_x_z_xxz_xyyz[k] = -g_0_z_xz_xyyz[k] - g_x_z_xz_xyyz[k] * ab_x + g_x_z_xz_xxyyz[k];

                g_x_z_xxz_xyzz[k] = -g_0_z_xz_xyzz[k] - g_x_z_xz_xyzz[k] * ab_x + g_x_z_xz_xxyzz[k];

                g_x_z_xxz_xzzz[k] = -g_0_z_xz_xzzz[k] - g_x_z_xz_xzzz[k] * ab_x + g_x_z_xz_xxzzz[k];

                g_x_z_xxz_yyyy[k] = -g_0_z_xz_yyyy[k] - g_x_z_xz_yyyy[k] * ab_x + g_x_z_xz_xyyyy[k];

                g_x_z_xxz_yyyz[k] = -g_0_z_xz_yyyz[k] - g_x_z_xz_yyyz[k] * ab_x + g_x_z_xz_xyyyz[k];

                g_x_z_xxz_yyzz[k] = -g_0_z_xz_yyzz[k] - g_x_z_xz_yyzz[k] * ab_x + g_x_z_xz_xyyzz[k];

                g_x_z_xxz_yzzz[k] = -g_0_z_xz_yzzz[k] - g_x_z_xz_yzzz[k] * ab_x + g_x_z_xz_xyzzz[k];

                g_x_z_xxz_zzzz[k] = -g_0_z_xz_zzzz[k] - g_x_z_xz_zzzz[k] * ab_x + g_x_z_xz_xzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyy_xxxx = cbuffer.data(fg_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_xyy_xxxy = cbuffer.data(fg_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_xyy_xxxz = cbuffer.data(fg_geom_11_off + 347 * ccomps * dcomps);

            auto g_x_z_xyy_xxyy = cbuffer.data(fg_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_xyy_xxyz = cbuffer.data(fg_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_xyy_xxzz = cbuffer.data(fg_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_xyy_xyyy = cbuffer.data(fg_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_xyy_xyyz = cbuffer.data(fg_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_xyy_xyzz = cbuffer.data(fg_geom_11_off + 353 * ccomps * dcomps);

            auto g_x_z_xyy_xzzz = cbuffer.data(fg_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_xyy_yyyy = cbuffer.data(fg_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_xyy_yyyz = cbuffer.data(fg_geom_11_off + 356 * ccomps * dcomps);

            auto g_x_z_xyy_yyzz = cbuffer.data(fg_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_xyy_yzzz = cbuffer.data(fg_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_xyy_zzzz = cbuffer.data(fg_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xy_xxxx, g_x_z_xy_xxxxy, g_x_z_xy_xxxy, g_x_z_xy_xxxyy, g_x_z_xy_xxxyz, g_x_z_xy_xxxz, g_x_z_xy_xxyy, g_x_z_xy_xxyyy, g_x_z_xy_xxyyz, g_x_z_xy_xxyz, g_x_z_xy_xxyzz, g_x_z_xy_xxzz, g_x_z_xy_xyyy, g_x_z_xy_xyyyy, g_x_z_xy_xyyyz, g_x_z_xy_xyyz, g_x_z_xy_xyyzz, g_x_z_xy_xyzz, g_x_z_xy_xyzzz, g_x_z_xy_xzzz, g_x_z_xy_yyyy, g_x_z_xy_yyyyy, g_x_z_xy_yyyyz, g_x_z_xy_yyyz, g_x_z_xy_yyyzz, g_x_z_xy_yyzz, g_x_z_xy_yyzzz, g_x_z_xy_yzzz, g_x_z_xy_yzzzz, g_x_z_xy_zzzz, g_x_z_xyy_xxxx, g_x_z_xyy_xxxy, g_x_z_xyy_xxxz, g_x_z_xyy_xxyy, g_x_z_xyy_xxyz, g_x_z_xyy_xxzz, g_x_z_xyy_xyyy, g_x_z_xyy_xyyz, g_x_z_xyy_xyzz, g_x_z_xyy_xzzz, g_x_z_xyy_yyyy, g_x_z_xyy_yyyz, g_x_z_xyy_yyzz, g_x_z_xyy_yzzz, g_x_z_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyy_xxxx[k] = -g_x_z_xy_xxxx[k] * ab_y + g_x_z_xy_xxxxy[k];

                g_x_z_xyy_xxxy[k] = -g_x_z_xy_xxxy[k] * ab_y + g_x_z_xy_xxxyy[k];

                g_x_z_xyy_xxxz[k] = -g_x_z_xy_xxxz[k] * ab_y + g_x_z_xy_xxxyz[k];

                g_x_z_xyy_xxyy[k] = -g_x_z_xy_xxyy[k] * ab_y + g_x_z_xy_xxyyy[k];

                g_x_z_xyy_xxyz[k] = -g_x_z_xy_xxyz[k] * ab_y + g_x_z_xy_xxyyz[k];

                g_x_z_xyy_xxzz[k] = -g_x_z_xy_xxzz[k] * ab_y + g_x_z_xy_xxyzz[k];

                g_x_z_xyy_xyyy[k] = -g_x_z_xy_xyyy[k] * ab_y + g_x_z_xy_xyyyy[k];

                g_x_z_xyy_xyyz[k] = -g_x_z_xy_xyyz[k] * ab_y + g_x_z_xy_xyyyz[k];

                g_x_z_xyy_xyzz[k] = -g_x_z_xy_xyzz[k] * ab_y + g_x_z_xy_xyyzz[k];

                g_x_z_xyy_xzzz[k] = -g_x_z_xy_xzzz[k] * ab_y + g_x_z_xy_xyzzz[k];

                g_x_z_xyy_yyyy[k] = -g_x_z_xy_yyyy[k] * ab_y + g_x_z_xy_yyyyy[k];

                g_x_z_xyy_yyyz[k] = -g_x_z_xy_yyyz[k] * ab_y + g_x_z_xy_yyyyz[k];

                g_x_z_xyy_yyzz[k] = -g_x_z_xy_yyzz[k] * ab_y + g_x_z_xy_yyyzz[k];

                g_x_z_xyy_yzzz[k] = -g_x_z_xy_yzzz[k] * ab_y + g_x_z_xy_yyzzz[k];

                g_x_z_xyy_zzzz[k] = -g_x_z_xy_zzzz[k] * ab_y + g_x_z_xy_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyz_xxxx = cbuffer.data(fg_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_xyz_xxxy = cbuffer.data(fg_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_xyz_xxxz = cbuffer.data(fg_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_xyz_xxyy = cbuffer.data(fg_geom_11_off + 363 * ccomps * dcomps);

            auto g_x_z_xyz_xxyz = cbuffer.data(fg_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_xyz_xxzz = cbuffer.data(fg_geom_11_off + 365 * ccomps * dcomps);

            auto g_x_z_xyz_xyyy = cbuffer.data(fg_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_xyz_xyyz = cbuffer.data(fg_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_xyz_xyzz = cbuffer.data(fg_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_xyz_xzzz = cbuffer.data(fg_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_xyz_yyyy = cbuffer.data(fg_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_xyz_yyyz = cbuffer.data(fg_geom_11_off + 371 * ccomps * dcomps);

            auto g_x_z_xyz_yyzz = cbuffer.data(fg_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_xyz_yzzz = cbuffer.data(fg_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_xyz_zzzz = cbuffer.data(fg_geom_11_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyz_xxxx, g_x_z_xyz_xxxy, g_x_z_xyz_xxxz, g_x_z_xyz_xxyy, g_x_z_xyz_xxyz, g_x_z_xyz_xxzz, g_x_z_xyz_xyyy, g_x_z_xyz_xyyz, g_x_z_xyz_xyzz, g_x_z_xyz_xzzz, g_x_z_xyz_yyyy, g_x_z_xyz_yyyz, g_x_z_xyz_yyzz, g_x_z_xyz_yzzz, g_x_z_xyz_zzzz, g_x_z_xz_xxxx, g_x_z_xz_xxxxy, g_x_z_xz_xxxy, g_x_z_xz_xxxyy, g_x_z_xz_xxxyz, g_x_z_xz_xxxz, g_x_z_xz_xxyy, g_x_z_xz_xxyyy, g_x_z_xz_xxyyz, g_x_z_xz_xxyz, g_x_z_xz_xxyzz, g_x_z_xz_xxzz, g_x_z_xz_xyyy, g_x_z_xz_xyyyy, g_x_z_xz_xyyyz, g_x_z_xz_xyyz, g_x_z_xz_xyyzz, g_x_z_xz_xyzz, g_x_z_xz_xyzzz, g_x_z_xz_xzzz, g_x_z_xz_yyyy, g_x_z_xz_yyyyy, g_x_z_xz_yyyyz, g_x_z_xz_yyyz, g_x_z_xz_yyyzz, g_x_z_xz_yyzz, g_x_z_xz_yyzzz, g_x_z_xz_yzzz, g_x_z_xz_yzzzz, g_x_z_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyz_xxxx[k] = -g_x_z_xz_xxxx[k] * ab_y + g_x_z_xz_xxxxy[k];

                g_x_z_xyz_xxxy[k] = -g_x_z_xz_xxxy[k] * ab_y + g_x_z_xz_xxxyy[k];

                g_x_z_xyz_xxxz[k] = -g_x_z_xz_xxxz[k] * ab_y + g_x_z_xz_xxxyz[k];

                g_x_z_xyz_xxyy[k] = -g_x_z_xz_xxyy[k] * ab_y + g_x_z_xz_xxyyy[k];

                g_x_z_xyz_xxyz[k] = -g_x_z_xz_xxyz[k] * ab_y + g_x_z_xz_xxyyz[k];

                g_x_z_xyz_xxzz[k] = -g_x_z_xz_xxzz[k] * ab_y + g_x_z_xz_xxyzz[k];

                g_x_z_xyz_xyyy[k] = -g_x_z_xz_xyyy[k] * ab_y + g_x_z_xz_xyyyy[k];

                g_x_z_xyz_xyyz[k] = -g_x_z_xz_xyyz[k] * ab_y + g_x_z_xz_xyyyz[k];

                g_x_z_xyz_xyzz[k] = -g_x_z_xz_xyzz[k] * ab_y + g_x_z_xz_xyyzz[k];

                g_x_z_xyz_xzzz[k] = -g_x_z_xz_xzzz[k] * ab_y + g_x_z_xz_xyzzz[k];

                g_x_z_xyz_yyyy[k] = -g_x_z_xz_yyyy[k] * ab_y + g_x_z_xz_yyyyy[k];

                g_x_z_xyz_yyyz[k] = -g_x_z_xz_yyyz[k] * ab_y + g_x_z_xz_yyyyz[k];

                g_x_z_xyz_yyzz[k] = -g_x_z_xz_yyzz[k] * ab_y + g_x_z_xz_yyyzz[k];

                g_x_z_xyz_yzzz[k] = -g_x_z_xz_yzzz[k] * ab_y + g_x_z_xz_yyzzz[k];

                g_x_z_xyz_zzzz[k] = -g_x_z_xz_zzzz[k] * ab_y + g_x_z_xz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzz_xxxx = cbuffer.data(fg_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_xzz_xxxy = cbuffer.data(fg_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_xzz_xxxz = cbuffer.data(fg_geom_11_off + 377 * ccomps * dcomps);

            auto g_x_z_xzz_xxyy = cbuffer.data(fg_geom_11_off + 378 * ccomps * dcomps);

            auto g_x_z_xzz_xxyz = cbuffer.data(fg_geom_11_off + 379 * ccomps * dcomps);

            auto g_x_z_xzz_xxzz = cbuffer.data(fg_geom_11_off + 380 * ccomps * dcomps);

            auto g_x_z_xzz_xyyy = cbuffer.data(fg_geom_11_off + 381 * ccomps * dcomps);

            auto g_x_z_xzz_xyyz = cbuffer.data(fg_geom_11_off + 382 * ccomps * dcomps);

            auto g_x_z_xzz_xyzz = cbuffer.data(fg_geom_11_off + 383 * ccomps * dcomps);

            auto g_x_z_xzz_xzzz = cbuffer.data(fg_geom_11_off + 384 * ccomps * dcomps);

            auto g_x_z_xzz_yyyy = cbuffer.data(fg_geom_11_off + 385 * ccomps * dcomps);

            auto g_x_z_xzz_yyyz = cbuffer.data(fg_geom_11_off + 386 * ccomps * dcomps);

            auto g_x_z_xzz_yyzz = cbuffer.data(fg_geom_11_off + 387 * ccomps * dcomps);

            auto g_x_z_xzz_yzzz = cbuffer.data(fg_geom_11_off + 388 * ccomps * dcomps);

            auto g_x_z_xzz_zzzz = cbuffer.data(fg_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxxx, g_0_z_zz_xxxy, g_0_z_zz_xxxz, g_0_z_zz_xxyy, g_0_z_zz_xxyz, g_0_z_zz_xxzz, g_0_z_zz_xyyy, g_0_z_zz_xyyz, g_0_z_zz_xyzz, g_0_z_zz_xzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyz, g_0_z_zz_yyzz, g_0_z_zz_yzzz, g_0_z_zz_zzzz, g_x_z_xzz_xxxx, g_x_z_xzz_xxxy, g_x_z_xzz_xxxz, g_x_z_xzz_xxyy, g_x_z_xzz_xxyz, g_x_z_xzz_xxzz, g_x_z_xzz_xyyy, g_x_z_xzz_xyyz, g_x_z_xzz_xyzz, g_x_z_xzz_xzzz, g_x_z_xzz_yyyy, g_x_z_xzz_yyyz, g_x_z_xzz_yyzz, g_x_z_xzz_yzzz, g_x_z_xzz_zzzz, g_x_z_zz_xxxx, g_x_z_zz_xxxxx, g_x_z_zz_xxxxy, g_x_z_zz_xxxxz, g_x_z_zz_xxxy, g_x_z_zz_xxxyy, g_x_z_zz_xxxyz, g_x_z_zz_xxxz, g_x_z_zz_xxxzz, g_x_z_zz_xxyy, g_x_z_zz_xxyyy, g_x_z_zz_xxyyz, g_x_z_zz_xxyz, g_x_z_zz_xxyzz, g_x_z_zz_xxzz, g_x_z_zz_xxzzz, g_x_z_zz_xyyy, g_x_z_zz_xyyyy, g_x_z_zz_xyyyz, g_x_z_zz_xyyz, g_x_z_zz_xyyzz, g_x_z_zz_xyzz, g_x_z_zz_xyzzz, g_x_z_zz_xzzz, g_x_z_zz_xzzzz, g_x_z_zz_yyyy, g_x_z_zz_yyyz, g_x_z_zz_yyzz, g_x_z_zz_yzzz, g_x_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzz_xxxx[k] = -g_0_z_zz_xxxx[k] - g_x_z_zz_xxxx[k] * ab_x + g_x_z_zz_xxxxx[k];

                g_x_z_xzz_xxxy[k] = -g_0_z_zz_xxxy[k] - g_x_z_zz_xxxy[k] * ab_x + g_x_z_zz_xxxxy[k];

                g_x_z_xzz_xxxz[k] = -g_0_z_zz_xxxz[k] - g_x_z_zz_xxxz[k] * ab_x + g_x_z_zz_xxxxz[k];

                g_x_z_xzz_xxyy[k] = -g_0_z_zz_xxyy[k] - g_x_z_zz_xxyy[k] * ab_x + g_x_z_zz_xxxyy[k];

                g_x_z_xzz_xxyz[k] = -g_0_z_zz_xxyz[k] - g_x_z_zz_xxyz[k] * ab_x + g_x_z_zz_xxxyz[k];

                g_x_z_xzz_xxzz[k] = -g_0_z_zz_xxzz[k] - g_x_z_zz_xxzz[k] * ab_x + g_x_z_zz_xxxzz[k];

                g_x_z_xzz_xyyy[k] = -g_0_z_zz_xyyy[k] - g_x_z_zz_xyyy[k] * ab_x + g_x_z_zz_xxyyy[k];

                g_x_z_xzz_xyyz[k] = -g_0_z_zz_xyyz[k] - g_x_z_zz_xyyz[k] * ab_x + g_x_z_zz_xxyyz[k];

                g_x_z_xzz_xyzz[k] = -g_0_z_zz_xyzz[k] - g_x_z_zz_xyzz[k] * ab_x + g_x_z_zz_xxyzz[k];

                g_x_z_xzz_xzzz[k] = -g_0_z_zz_xzzz[k] - g_x_z_zz_xzzz[k] * ab_x + g_x_z_zz_xxzzz[k];

                g_x_z_xzz_yyyy[k] = -g_0_z_zz_yyyy[k] - g_x_z_zz_yyyy[k] * ab_x + g_x_z_zz_xyyyy[k];

                g_x_z_xzz_yyyz[k] = -g_0_z_zz_yyyz[k] - g_x_z_zz_yyyz[k] * ab_x + g_x_z_zz_xyyyz[k];

                g_x_z_xzz_yyzz[k] = -g_0_z_zz_yyzz[k] - g_x_z_zz_yyzz[k] * ab_x + g_x_z_zz_xyyzz[k];

                g_x_z_xzz_yzzz[k] = -g_0_z_zz_yzzz[k] - g_x_z_zz_yzzz[k] * ab_x + g_x_z_zz_xyzzz[k];

                g_x_z_xzz_zzzz[k] = -g_0_z_zz_zzzz[k] - g_x_z_zz_zzzz[k] * ab_x + g_x_z_zz_xzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyy_xxxx = cbuffer.data(fg_geom_11_off + 390 * ccomps * dcomps);

            auto g_x_z_yyy_xxxy = cbuffer.data(fg_geom_11_off + 391 * ccomps * dcomps);

            auto g_x_z_yyy_xxxz = cbuffer.data(fg_geom_11_off + 392 * ccomps * dcomps);

            auto g_x_z_yyy_xxyy = cbuffer.data(fg_geom_11_off + 393 * ccomps * dcomps);

            auto g_x_z_yyy_xxyz = cbuffer.data(fg_geom_11_off + 394 * ccomps * dcomps);

            auto g_x_z_yyy_xxzz = cbuffer.data(fg_geom_11_off + 395 * ccomps * dcomps);

            auto g_x_z_yyy_xyyy = cbuffer.data(fg_geom_11_off + 396 * ccomps * dcomps);

            auto g_x_z_yyy_xyyz = cbuffer.data(fg_geom_11_off + 397 * ccomps * dcomps);

            auto g_x_z_yyy_xyzz = cbuffer.data(fg_geom_11_off + 398 * ccomps * dcomps);

            auto g_x_z_yyy_xzzz = cbuffer.data(fg_geom_11_off + 399 * ccomps * dcomps);

            auto g_x_z_yyy_yyyy = cbuffer.data(fg_geom_11_off + 400 * ccomps * dcomps);

            auto g_x_z_yyy_yyyz = cbuffer.data(fg_geom_11_off + 401 * ccomps * dcomps);

            auto g_x_z_yyy_yyzz = cbuffer.data(fg_geom_11_off + 402 * ccomps * dcomps);

            auto g_x_z_yyy_yzzz = cbuffer.data(fg_geom_11_off + 403 * ccomps * dcomps);

            auto g_x_z_yyy_zzzz = cbuffer.data(fg_geom_11_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yy_xxxx, g_x_z_yy_xxxxy, g_x_z_yy_xxxy, g_x_z_yy_xxxyy, g_x_z_yy_xxxyz, g_x_z_yy_xxxz, g_x_z_yy_xxyy, g_x_z_yy_xxyyy, g_x_z_yy_xxyyz, g_x_z_yy_xxyz, g_x_z_yy_xxyzz, g_x_z_yy_xxzz, g_x_z_yy_xyyy, g_x_z_yy_xyyyy, g_x_z_yy_xyyyz, g_x_z_yy_xyyz, g_x_z_yy_xyyzz, g_x_z_yy_xyzz, g_x_z_yy_xyzzz, g_x_z_yy_xzzz, g_x_z_yy_yyyy, g_x_z_yy_yyyyy, g_x_z_yy_yyyyz, g_x_z_yy_yyyz, g_x_z_yy_yyyzz, g_x_z_yy_yyzz, g_x_z_yy_yyzzz, g_x_z_yy_yzzz, g_x_z_yy_yzzzz, g_x_z_yy_zzzz, g_x_z_yyy_xxxx, g_x_z_yyy_xxxy, g_x_z_yyy_xxxz, g_x_z_yyy_xxyy, g_x_z_yyy_xxyz, g_x_z_yyy_xxzz, g_x_z_yyy_xyyy, g_x_z_yyy_xyyz, g_x_z_yyy_xyzz, g_x_z_yyy_xzzz, g_x_z_yyy_yyyy, g_x_z_yyy_yyyz, g_x_z_yyy_yyzz, g_x_z_yyy_yzzz, g_x_z_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyy_xxxx[k] = -g_x_z_yy_xxxx[k] * ab_y + g_x_z_yy_xxxxy[k];

                g_x_z_yyy_xxxy[k] = -g_x_z_yy_xxxy[k] * ab_y + g_x_z_yy_xxxyy[k];

                g_x_z_yyy_xxxz[k] = -g_x_z_yy_xxxz[k] * ab_y + g_x_z_yy_xxxyz[k];

                g_x_z_yyy_xxyy[k] = -g_x_z_yy_xxyy[k] * ab_y + g_x_z_yy_xxyyy[k];

                g_x_z_yyy_xxyz[k] = -g_x_z_yy_xxyz[k] * ab_y + g_x_z_yy_xxyyz[k];

                g_x_z_yyy_xxzz[k] = -g_x_z_yy_xxzz[k] * ab_y + g_x_z_yy_xxyzz[k];

                g_x_z_yyy_xyyy[k] = -g_x_z_yy_xyyy[k] * ab_y + g_x_z_yy_xyyyy[k];

                g_x_z_yyy_xyyz[k] = -g_x_z_yy_xyyz[k] * ab_y + g_x_z_yy_xyyyz[k];

                g_x_z_yyy_xyzz[k] = -g_x_z_yy_xyzz[k] * ab_y + g_x_z_yy_xyyzz[k];

                g_x_z_yyy_xzzz[k] = -g_x_z_yy_xzzz[k] * ab_y + g_x_z_yy_xyzzz[k];

                g_x_z_yyy_yyyy[k] = -g_x_z_yy_yyyy[k] * ab_y + g_x_z_yy_yyyyy[k];

                g_x_z_yyy_yyyz[k] = -g_x_z_yy_yyyz[k] * ab_y + g_x_z_yy_yyyyz[k];

                g_x_z_yyy_yyzz[k] = -g_x_z_yy_yyzz[k] * ab_y + g_x_z_yy_yyyzz[k];

                g_x_z_yyy_yzzz[k] = -g_x_z_yy_yzzz[k] * ab_y + g_x_z_yy_yyzzz[k];

                g_x_z_yyy_zzzz[k] = -g_x_z_yy_zzzz[k] * ab_y + g_x_z_yy_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyz_xxxx = cbuffer.data(fg_geom_11_off + 405 * ccomps * dcomps);

            auto g_x_z_yyz_xxxy = cbuffer.data(fg_geom_11_off + 406 * ccomps * dcomps);

            auto g_x_z_yyz_xxxz = cbuffer.data(fg_geom_11_off + 407 * ccomps * dcomps);

            auto g_x_z_yyz_xxyy = cbuffer.data(fg_geom_11_off + 408 * ccomps * dcomps);

            auto g_x_z_yyz_xxyz = cbuffer.data(fg_geom_11_off + 409 * ccomps * dcomps);

            auto g_x_z_yyz_xxzz = cbuffer.data(fg_geom_11_off + 410 * ccomps * dcomps);

            auto g_x_z_yyz_xyyy = cbuffer.data(fg_geom_11_off + 411 * ccomps * dcomps);

            auto g_x_z_yyz_xyyz = cbuffer.data(fg_geom_11_off + 412 * ccomps * dcomps);

            auto g_x_z_yyz_xyzz = cbuffer.data(fg_geom_11_off + 413 * ccomps * dcomps);

            auto g_x_z_yyz_xzzz = cbuffer.data(fg_geom_11_off + 414 * ccomps * dcomps);

            auto g_x_z_yyz_yyyy = cbuffer.data(fg_geom_11_off + 415 * ccomps * dcomps);

            auto g_x_z_yyz_yyyz = cbuffer.data(fg_geom_11_off + 416 * ccomps * dcomps);

            auto g_x_z_yyz_yyzz = cbuffer.data(fg_geom_11_off + 417 * ccomps * dcomps);

            auto g_x_z_yyz_yzzz = cbuffer.data(fg_geom_11_off + 418 * ccomps * dcomps);

            auto g_x_z_yyz_zzzz = cbuffer.data(fg_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyz_xxxx, g_x_z_yyz_xxxy, g_x_z_yyz_xxxz, g_x_z_yyz_xxyy, g_x_z_yyz_xxyz, g_x_z_yyz_xxzz, g_x_z_yyz_xyyy, g_x_z_yyz_xyyz, g_x_z_yyz_xyzz, g_x_z_yyz_xzzz, g_x_z_yyz_yyyy, g_x_z_yyz_yyyz, g_x_z_yyz_yyzz, g_x_z_yyz_yzzz, g_x_z_yyz_zzzz, g_x_z_yz_xxxx, g_x_z_yz_xxxxy, g_x_z_yz_xxxy, g_x_z_yz_xxxyy, g_x_z_yz_xxxyz, g_x_z_yz_xxxz, g_x_z_yz_xxyy, g_x_z_yz_xxyyy, g_x_z_yz_xxyyz, g_x_z_yz_xxyz, g_x_z_yz_xxyzz, g_x_z_yz_xxzz, g_x_z_yz_xyyy, g_x_z_yz_xyyyy, g_x_z_yz_xyyyz, g_x_z_yz_xyyz, g_x_z_yz_xyyzz, g_x_z_yz_xyzz, g_x_z_yz_xyzzz, g_x_z_yz_xzzz, g_x_z_yz_yyyy, g_x_z_yz_yyyyy, g_x_z_yz_yyyyz, g_x_z_yz_yyyz, g_x_z_yz_yyyzz, g_x_z_yz_yyzz, g_x_z_yz_yyzzz, g_x_z_yz_yzzz, g_x_z_yz_yzzzz, g_x_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyz_xxxx[k] = -g_x_z_yz_xxxx[k] * ab_y + g_x_z_yz_xxxxy[k];

                g_x_z_yyz_xxxy[k] = -g_x_z_yz_xxxy[k] * ab_y + g_x_z_yz_xxxyy[k];

                g_x_z_yyz_xxxz[k] = -g_x_z_yz_xxxz[k] * ab_y + g_x_z_yz_xxxyz[k];

                g_x_z_yyz_xxyy[k] = -g_x_z_yz_xxyy[k] * ab_y + g_x_z_yz_xxyyy[k];

                g_x_z_yyz_xxyz[k] = -g_x_z_yz_xxyz[k] * ab_y + g_x_z_yz_xxyyz[k];

                g_x_z_yyz_xxzz[k] = -g_x_z_yz_xxzz[k] * ab_y + g_x_z_yz_xxyzz[k];

                g_x_z_yyz_xyyy[k] = -g_x_z_yz_xyyy[k] * ab_y + g_x_z_yz_xyyyy[k];

                g_x_z_yyz_xyyz[k] = -g_x_z_yz_xyyz[k] * ab_y + g_x_z_yz_xyyyz[k];

                g_x_z_yyz_xyzz[k] = -g_x_z_yz_xyzz[k] * ab_y + g_x_z_yz_xyyzz[k];

                g_x_z_yyz_xzzz[k] = -g_x_z_yz_xzzz[k] * ab_y + g_x_z_yz_xyzzz[k];

                g_x_z_yyz_yyyy[k] = -g_x_z_yz_yyyy[k] * ab_y + g_x_z_yz_yyyyy[k];

                g_x_z_yyz_yyyz[k] = -g_x_z_yz_yyyz[k] * ab_y + g_x_z_yz_yyyyz[k];

                g_x_z_yyz_yyzz[k] = -g_x_z_yz_yyzz[k] * ab_y + g_x_z_yz_yyyzz[k];

                g_x_z_yyz_yzzz[k] = -g_x_z_yz_yzzz[k] * ab_y + g_x_z_yz_yyzzz[k];

                g_x_z_yyz_zzzz[k] = -g_x_z_yz_zzzz[k] * ab_y + g_x_z_yz_yzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzz_xxxx = cbuffer.data(fg_geom_11_off + 420 * ccomps * dcomps);

            auto g_x_z_yzz_xxxy = cbuffer.data(fg_geom_11_off + 421 * ccomps * dcomps);

            auto g_x_z_yzz_xxxz = cbuffer.data(fg_geom_11_off + 422 * ccomps * dcomps);

            auto g_x_z_yzz_xxyy = cbuffer.data(fg_geom_11_off + 423 * ccomps * dcomps);

            auto g_x_z_yzz_xxyz = cbuffer.data(fg_geom_11_off + 424 * ccomps * dcomps);

            auto g_x_z_yzz_xxzz = cbuffer.data(fg_geom_11_off + 425 * ccomps * dcomps);

            auto g_x_z_yzz_xyyy = cbuffer.data(fg_geom_11_off + 426 * ccomps * dcomps);

            auto g_x_z_yzz_xyyz = cbuffer.data(fg_geom_11_off + 427 * ccomps * dcomps);

            auto g_x_z_yzz_xyzz = cbuffer.data(fg_geom_11_off + 428 * ccomps * dcomps);

            auto g_x_z_yzz_xzzz = cbuffer.data(fg_geom_11_off + 429 * ccomps * dcomps);

            auto g_x_z_yzz_yyyy = cbuffer.data(fg_geom_11_off + 430 * ccomps * dcomps);

            auto g_x_z_yzz_yyyz = cbuffer.data(fg_geom_11_off + 431 * ccomps * dcomps);

            auto g_x_z_yzz_yyzz = cbuffer.data(fg_geom_11_off + 432 * ccomps * dcomps);

            auto g_x_z_yzz_yzzz = cbuffer.data(fg_geom_11_off + 433 * ccomps * dcomps);

            auto g_x_z_yzz_zzzz = cbuffer.data(fg_geom_11_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzz_xxxx, g_x_z_yzz_xxxy, g_x_z_yzz_xxxz, g_x_z_yzz_xxyy, g_x_z_yzz_xxyz, g_x_z_yzz_xxzz, g_x_z_yzz_xyyy, g_x_z_yzz_xyyz, g_x_z_yzz_xyzz, g_x_z_yzz_xzzz, g_x_z_yzz_yyyy, g_x_z_yzz_yyyz, g_x_z_yzz_yyzz, g_x_z_yzz_yzzz, g_x_z_yzz_zzzz, g_x_z_zz_xxxx, g_x_z_zz_xxxxy, g_x_z_zz_xxxy, g_x_z_zz_xxxyy, g_x_z_zz_xxxyz, g_x_z_zz_xxxz, g_x_z_zz_xxyy, g_x_z_zz_xxyyy, g_x_z_zz_xxyyz, g_x_z_zz_xxyz, g_x_z_zz_xxyzz, g_x_z_zz_xxzz, g_x_z_zz_xyyy, g_x_z_zz_xyyyy, g_x_z_zz_xyyyz, g_x_z_zz_xyyz, g_x_z_zz_xyyzz, g_x_z_zz_xyzz, g_x_z_zz_xyzzz, g_x_z_zz_xzzz, g_x_z_zz_yyyy, g_x_z_zz_yyyyy, g_x_z_zz_yyyyz, g_x_z_zz_yyyz, g_x_z_zz_yyyzz, g_x_z_zz_yyzz, g_x_z_zz_yyzzz, g_x_z_zz_yzzz, g_x_z_zz_yzzzz, g_x_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzz_xxxx[k] = -g_x_z_zz_xxxx[k] * ab_y + g_x_z_zz_xxxxy[k];

                g_x_z_yzz_xxxy[k] = -g_x_z_zz_xxxy[k] * ab_y + g_x_z_zz_xxxyy[k];

                g_x_z_yzz_xxxz[k] = -g_x_z_zz_xxxz[k] * ab_y + g_x_z_zz_xxxyz[k];

                g_x_z_yzz_xxyy[k] = -g_x_z_zz_xxyy[k] * ab_y + g_x_z_zz_xxyyy[k];

                g_x_z_yzz_xxyz[k] = -g_x_z_zz_xxyz[k] * ab_y + g_x_z_zz_xxyyz[k];

                g_x_z_yzz_xxzz[k] = -g_x_z_zz_xxzz[k] * ab_y + g_x_z_zz_xxyzz[k];

                g_x_z_yzz_xyyy[k] = -g_x_z_zz_xyyy[k] * ab_y + g_x_z_zz_xyyyy[k];

                g_x_z_yzz_xyyz[k] = -g_x_z_zz_xyyz[k] * ab_y + g_x_z_zz_xyyyz[k];

                g_x_z_yzz_xyzz[k] = -g_x_z_zz_xyzz[k] * ab_y + g_x_z_zz_xyyzz[k];

                g_x_z_yzz_xzzz[k] = -g_x_z_zz_xzzz[k] * ab_y + g_x_z_zz_xyzzz[k];

                g_x_z_yzz_yyyy[k] = -g_x_z_zz_yyyy[k] * ab_y + g_x_z_zz_yyyyy[k];

                g_x_z_yzz_yyyz[k] = -g_x_z_zz_yyyz[k] * ab_y + g_x_z_zz_yyyyz[k];

                g_x_z_yzz_yyzz[k] = -g_x_z_zz_yyzz[k] * ab_y + g_x_z_zz_yyyzz[k];

                g_x_z_yzz_yzzz[k] = -g_x_z_zz_yzzz[k] * ab_y + g_x_z_zz_yyzzz[k];

                g_x_z_yzz_zzzz[k] = -g_x_z_zz_zzzz[k] * ab_y + g_x_z_zz_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzz_xxxx = cbuffer.data(fg_geom_11_off + 435 * ccomps * dcomps);

            auto g_x_z_zzz_xxxy = cbuffer.data(fg_geom_11_off + 436 * ccomps * dcomps);

            auto g_x_z_zzz_xxxz = cbuffer.data(fg_geom_11_off + 437 * ccomps * dcomps);

            auto g_x_z_zzz_xxyy = cbuffer.data(fg_geom_11_off + 438 * ccomps * dcomps);

            auto g_x_z_zzz_xxyz = cbuffer.data(fg_geom_11_off + 439 * ccomps * dcomps);

            auto g_x_z_zzz_xxzz = cbuffer.data(fg_geom_11_off + 440 * ccomps * dcomps);

            auto g_x_z_zzz_xyyy = cbuffer.data(fg_geom_11_off + 441 * ccomps * dcomps);

            auto g_x_z_zzz_xyyz = cbuffer.data(fg_geom_11_off + 442 * ccomps * dcomps);

            auto g_x_z_zzz_xyzz = cbuffer.data(fg_geom_11_off + 443 * ccomps * dcomps);

            auto g_x_z_zzz_xzzz = cbuffer.data(fg_geom_11_off + 444 * ccomps * dcomps);

            auto g_x_z_zzz_yyyy = cbuffer.data(fg_geom_11_off + 445 * ccomps * dcomps);

            auto g_x_z_zzz_yyyz = cbuffer.data(fg_geom_11_off + 446 * ccomps * dcomps);

            auto g_x_z_zzz_yyzz = cbuffer.data(fg_geom_11_off + 447 * ccomps * dcomps);

            auto g_x_z_zzz_yzzz = cbuffer.data(fg_geom_11_off + 448 * ccomps * dcomps);

            auto g_x_z_zzz_zzzz = cbuffer.data(fg_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xxxx, g_x_0_zz_xxxy, g_x_0_zz_xxxz, g_x_0_zz_xxyy, g_x_0_zz_xxyz, g_x_0_zz_xxzz, g_x_0_zz_xyyy, g_x_0_zz_xyyz, g_x_0_zz_xyzz, g_x_0_zz_xzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyz, g_x_0_zz_yyzz, g_x_0_zz_yzzz, g_x_0_zz_zzzz, g_x_z_zz_xxxx, g_x_z_zz_xxxxz, g_x_z_zz_xxxy, g_x_z_zz_xxxyz, g_x_z_zz_xxxz, g_x_z_zz_xxxzz, g_x_z_zz_xxyy, g_x_z_zz_xxyyz, g_x_z_zz_xxyz, g_x_z_zz_xxyzz, g_x_z_zz_xxzz, g_x_z_zz_xxzzz, g_x_z_zz_xyyy, g_x_z_zz_xyyyz, g_x_z_zz_xyyz, g_x_z_zz_xyyzz, g_x_z_zz_xyzz, g_x_z_zz_xyzzz, g_x_z_zz_xzzz, g_x_z_zz_xzzzz, g_x_z_zz_yyyy, g_x_z_zz_yyyyz, g_x_z_zz_yyyz, g_x_z_zz_yyyzz, g_x_z_zz_yyzz, g_x_z_zz_yyzzz, g_x_z_zz_yzzz, g_x_z_zz_yzzzz, g_x_z_zz_zzzz, g_x_z_zz_zzzzz, g_x_z_zzz_xxxx, g_x_z_zzz_xxxy, g_x_z_zzz_xxxz, g_x_z_zzz_xxyy, g_x_z_zzz_xxyz, g_x_z_zzz_xxzz, g_x_z_zzz_xyyy, g_x_z_zzz_xyyz, g_x_z_zzz_xyzz, g_x_z_zzz_xzzz, g_x_z_zzz_yyyy, g_x_z_zzz_yyyz, g_x_z_zzz_yyzz, g_x_z_zzz_yzzz, g_x_z_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzz_xxxx[k] = g_x_0_zz_xxxx[k] - g_x_z_zz_xxxx[k] * ab_z + g_x_z_zz_xxxxz[k];

                g_x_z_zzz_xxxy[k] = g_x_0_zz_xxxy[k] - g_x_z_zz_xxxy[k] * ab_z + g_x_z_zz_xxxyz[k];

                g_x_z_zzz_xxxz[k] = g_x_0_zz_xxxz[k] - g_x_z_zz_xxxz[k] * ab_z + g_x_z_zz_xxxzz[k];

                g_x_z_zzz_xxyy[k] = g_x_0_zz_xxyy[k] - g_x_z_zz_xxyy[k] * ab_z + g_x_z_zz_xxyyz[k];

                g_x_z_zzz_xxyz[k] = g_x_0_zz_xxyz[k] - g_x_z_zz_xxyz[k] * ab_z + g_x_z_zz_xxyzz[k];

                g_x_z_zzz_xxzz[k] = g_x_0_zz_xxzz[k] - g_x_z_zz_xxzz[k] * ab_z + g_x_z_zz_xxzzz[k];

                g_x_z_zzz_xyyy[k] = g_x_0_zz_xyyy[k] - g_x_z_zz_xyyy[k] * ab_z + g_x_z_zz_xyyyz[k];

                g_x_z_zzz_xyyz[k] = g_x_0_zz_xyyz[k] - g_x_z_zz_xyyz[k] * ab_z + g_x_z_zz_xyyzz[k];

                g_x_z_zzz_xyzz[k] = g_x_0_zz_xyzz[k] - g_x_z_zz_xyzz[k] * ab_z + g_x_z_zz_xyzzz[k];

                g_x_z_zzz_xzzz[k] = g_x_0_zz_xzzz[k] - g_x_z_zz_xzzz[k] * ab_z + g_x_z_zz_xzzzz[k];

                g_x_z_zzz_yyyy[k] = g_x_0_zz_yyyy[k] - g_x_z_zz_yyyy[k] * ab_z + g_x_z_zz_yyyyz[k];

                g_x_z_zzz_yyyz[k] = g_x_0_zz_yyyz[k] - g_x_z_zz_yyyz[k] * ab_z + g_x_z_zz_yyyzz[k];

                g_x_z_zzz_yyzz[k] = g_x_0_zz_yyzz[k] - g_x_z_zz_yyzz[k] * ab_z + g_x_z_zz_yyzzz[k];

                g_x_z_zzz_yzzz[k] = g_x_0_zz_yzzz[k] - g_x_z_zz_yzzz[k] * ab_z + g_x_z_zz_yzzzz[k];

                g_x_z_zzz_zzzz[k] = g_x_0_zz_zzzz[k] - g_x_z_zz_zzzz[k] * ab_z + g_x_z_zz_zzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxx_xxxx = cbuffer.data(fg_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_x_xxx_xxxy = cbuffer.data(fg_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_x_xxx_xxxz = cbuffer.data(fg_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_x_xxx_xxyy = cbuffer.data(fg_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_x_xxx_xxyz = cbuffer.data(fg_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_x_xxx_xxzz = cbuffer.data(fg_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_x_xxx_xyyy = cbuffer.data(fg_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_x_xxx_xyyz = cbuffer.data(fg_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_x_xxx_xyzz = cbuffer.data(fg_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_x_xxx_xzzz = cbuffer.data(fg_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_x_xxx_yyyy = cbuffer.data(fg_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_x_xxx_yyyz = cbuffer.data(fg_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_x_xxx_yyzz = cbuffer.data(fg_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_x_xxx_yzzz = cbuffer.data(fg_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_x_xxx_zzzz = cbuffer.data(fg_geom_11_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_xxxx, g_y_0_xx_xxxy, g_y_0_xx_xxxz, g_y_0_xx_xxyy, g_y_0_xx_xxyz, g_y_0_xx_xxzz, g_y_0_xx_xyyy, g_y_0_xx_xyyz, g_y_0_xx_xyzz, g_y_0_xx_xzzz, g_y_0_xx_yyyy, g_y_0_xx_yyyz, g_y_0_xx_yyzz, g_y_0_xx_yzzz, g_y_0_xx_zzzz, g_y_x_xx_xxxx, g_y_x_xx_xxxxx, g_y_x_xx_xxxxy, g_y_x_xx_xxxxz, g_y_x_xx_xxxy, g_y_x_xx_xxxyy, g_y_x_xx_xxxyz, g_y_x_xx_xxxz, g_y_x_xx_xxxzz, g_y_x_xx_xxyy, g_y_x_xx_xxyyy, g_y_x_xx_xxyyz, g_y_x_xx_xxyz, g_y_x_xx_xxyzz, g_y_x_xx_xxzz, g_y_x_xx_xxzzz, g_y_x_xx_xyyy, g_y_x_xx_xyyyy, g_y_x_xx_xyyyz, g_y_x_xx_xyyz, g_y_x_xx_xyyzz, g_y_x_xx_xyzz, g_y_x_xx_xyzzz, g_y_x_xx_xzzz, g_y_x_xx_xzzzz, g_y_x_xx_yyyy, g_y_x_xx_yyyz, g_y_x_xx_yyzz, g_y_x_xx_yzzz, g_y_x_xx_zzzz, g_y_x_xxx_xxxx, g_y_x_xxx_xxxy, g_y_x_xxx_xxxz, g_y_x_xxx_xxyy, g_y_x_xxx_xxyz, g_y_x_xxx_xxzz, g_y_x_xxx_xyyy, g_y_x_xxx_xyyz, g_y_x_xxx_xyzz, g_y_x_xxx_xzzz, g_y_x_xxx_yyyy, g_y_x_xxx_yyyz, g_y_x_xxx_yyzz, g_y_x_xxx_yzzz, g_y_x_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxx_xxxx[k] = g_y_0_xx_xxxx[k] - g_y_x_xx_xxxx[k] * ab_x + g_y_x_xx_xxxxx[k];

                g_y_x_xxx_xxxy[k] = g_y_0_xx_xxxy[k] - g_y_x_xx_xxxy[k] * ab_x + g_y_x_xx_xxxxy[k];

                g_y_x_xxx_xxxz[k] = g_y_0_xx_xxxz[k] - g_y_x_xx_xxxz[k] * ab_x + g_y_x_xx_xxxxz[k];

                g_y_x_xxx_xxyy[k] = g_y_0_xx_xxyy[k] - g_y_x_xx_xxyy[k] * ab_x + g_y_x_xx_xxxyy[k];

                g_y_x_xxx_xxyz[k] = g_y_0_xx_xxyz[k] - g_y_x_xx_xxyz[k] * ab_x + g_y_x_xx_xxxyz[k];

                g_y_x_xxx_xxzz[k] = g_y_0_xx_xxzz[k] - g_y_x_xx_xxzz[k] * ab_x + g_y_x_xx_xxxzz[k];

                g_y_x_xxx_xyyy[k] = g_y_0_xx_xyyy[k] - g_y_x_xx_xyyy[k] * ab_x + g_y_x_xx_xxyyy[k];

                g_y_x_xxx_xyyz[k] = g_y_0_xx_xyyz[k] - g_y_x_xx_xyyz[k] * ab_x + g_y_x_xx_xxyyz[k];

                g_y_x_xxx_xyzz[k] = g_y_0_xx_xyzz[k] - g_y_x_xx_xyzz[k] * ab_x + g_y_x_xx_xxyzz[k];

                g_y_x_xxx_xzzz[k] = g_y_0_xx_xzzz[k] - g_y_x_xx_xzzz[k] * ab_x + g_y_x_xx_xxzzz[k];

                g_y_x_xxx_yyyy[k] = g_y_0_xx_yyyy[k] - g_y_x_xx_yyyy[k] * ab_x + g_y_x_xx_xyyyy[k];

                g_y_x_xxx_yyyz[k] = g_y_0_xx_yyyz[k] - g_y_x_xx_yyyz[k] * ab_x + g_y_x_xx_xyyyz[k];

                g_y_x_xxx_yyzz[k] = g_y_0_xx_yyzz[k] - g_y_x_xx_yyzz[k] * ab_x + g_y_x_xx_xyyzz[k];

                g_y_x_xxx_yzzz[k] = g_y_0_xx_yzzz[k] - g_y_x_xx_yzzz[k] * ab_x + g_y_x_xx_xyzzz[k];

                g_y_x_xxx_zzzz[k] = g_y_0_xx_zzzz[k] - g_y_x_xx_zzzz[k] * ab_x + g_y_x_xx_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxy_xxxx = cbuffer.data(fg_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_x_xxy_xxxy = cbuffer.data(fg_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_x_xxy_xxxz = cbuffer.data(fg_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_x_xxy_xxyy = cbuffer.data(fg_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_x_xxy_xxyz = cbuffer.data(fg_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_x_xxy_xxzz = cbuffer.data(fg_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_x_xxy_xyyy = cbuffer.data(fg_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_x_xxy_xyyz = cbuffer.data(fg_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_x_xxy_xyzz = cbuffer.data(fg_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_x_xxy_xzzz = cbuffer.data(fg_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_x_xxy_yyyy = cbuffer.data(fg_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_x_xxy_yyyz = cbuffer.data(fg_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_x_xxy_yyzz = cbuffer.data(fg_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_x_xxy_yzzz = cbuffer.data(fg_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_x_xxy_zzzz = cbuffer.data(fg_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xxxx, g_y_0_xy_xxxy, g_y_0_xy_xxxz, g_y_0_xy_xxyy, g_y_0_xy_xxyz, g_y_0_xy_xxzz, g_y_0_xy_xyyy, g_y_0_xy_xyyz, g_y_0_xy_xyzz, g_y_0_xy_xzzz, g_y_0_xy_yyyy, g_y_0_xy_yyyz, g_y_0_xy_yyzz, g_y_0_xy_yzzz, g_y_0_xy_zzzz, g_y_x_xxy_xxxx, g_y_x_xxy_xxxy, g_y_x_xxy_xxxz, g_y_x_xxy_xxyy, g_y_x_xxy_xxyz, g_y_x_xxy_xxzz, g_y_x_xxy_xyyy, g_y_x_xxy_xyyz, g_y_x_xxy_xyzz, g_y_x_xxy_xzzz, g_y_x_xxy_yyyy, g_y_x_xxy_yyyz, g_y_x_xxy_yyzz, g_y_x_xxy_yzzz, g_y_x_xxy_zzzz, g_y_x_xy_xxxx, g_y_x_xy_xxxxx, g_y_x_xy_xxxxy, g_y_x_xy_xxxxz, g_y_x_xy_xxxy, g_y_x_xy_xxxyy, g_y_x_xy_xxxyz, g_y_x_xy_xxxz, g_y_x_xy_xxxzz, g_y_x_xy_xxyy, g_y_x_xy_xxyyy, g_y_x_xy_xxyyz, g_y_x_xy_xxyz, g_y_x_xy_xxyzz, g_y_x_xy_xxzz, g_y_x_xy_xxzzz, g_y_x_xy_xyyy, g_y_x_xy_xyyyy, g_y_x_xy_xyyyz, g_y_x_xy_xyyz, g_y_x_xy_xyyzz, g_y_x_xy_xyzz, g_y_x_xy_xyzzz, g_y_x_xy_xzzz, g_y_x_xy_xzzzz, g_y_x_xy_yyyy, g_y_x_xy_yyyz, g_y_x_xy_yyzz, g_y_x_xy_yzzz, g_y_x_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxy_xxxx[k] = g_y_0_xy_xxxx[k] - g_y_x_xy_xxxx[k] * ab_x + g_y_x_xy_xxxxx[k];

                g_y_x_xxy_xxxy[k] = g_y_0_xy_xxxy[k] - g_y_x_xy_xxxy[k] * ab_x + g_y_x_xy_xxxxy[k];

                g_y_x_xxy_xxxz[k] = g_y_0_xy_xxxz[k] - g_y_x_xy_xxxz[k] * ab_x + g_y_x_xy_xxxxz[k];

                g_y_x_xxy_xxyy[k] = g_y_0_xy_xxyy[k] - g_y_x_xy_xxyy[k] * ab_x + g_y_x_xy_xxxyy[k];

                g_y_x_xxy_xxyz[k] = g_y_0_xy_xxyz[k] - g_y_x_xy_xxyz[k] * ab_x + g_y_x_xy_xxxyz[k];

                g_y_x_xxy_xxzz[k] = g_y_0_xy_xxzz[k] - g_y_x_xy_xxzz[k] * ab_x + g_y_x_xy_xxxzz[k];

                g_y_x_xxy_xyyy[k] = g_y_0_xy_xyyy[k] - g_y_x_xy_xyyy[k] * ab_x + g_y_x_xy_xxyyy[k];

                g_y_x_xxy_xyyz[k] = g_y_0_xy_xyyz[k] - g_y_x_xy_xyyz[k] * ab_x + g_y_x_xy_xxyyz[k];

                g_y_x_xxy_xyzz[k] = g_y_0_xy_xyzz[k] - g_y_x_xy_xyzz[k] * ab_x + g_y_x_xy_xxyzz[k];

                g_y_x_xxy_xzzz[k] = g_y_0_xy_xzzz[k] - g_y_x_xy_xzzz[k] * ab_x + g_y_x_xy_xxzzz[k];

                g_y_x_xxy_yyyy[k] = g_y_0_xy_yyyy[k] - g_y_x_xy_yyyy[k] * ab_x + g_y_x_xy_xyyyy[k];

                g_y_x_xxy_yyyz[k] = g_y_0_xy_yyyz[k] - g_y_x_xy_yyyz[k] * ab_x + g_y_x_xy_xyyyz[k];

                g_y_x_xxy_yyzz[k] = g_y_0_xy_yyzz[k] - g_y_x_xy_yyzz[k] * ab_x + g_y_x_xy_xyyzz[k];

                g_y_x_xxy_yzzz[k] = g_y_0_xy_yzzz[k] - g_y_x_xy_yzzz[k] * ab_x + g_y_x_xy_xyzzz[k];

                g_y_x_xxy_zzzz[k] = g_y_0_xy_zzzz[k] - g_y_x_xy_zzzz[k] * ab_x + g_y_x_xy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxz_xxxx = cbuffer.data(fg_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_x_xxz_xxxy = cbuffer.data(fg_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_x_xxz_xxxz = cbuffer.data(fg_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_x_xxz_xxyy = cbuffer.data(fg_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_x_xxz_xxyz = cbuffer.data(fg_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_x_xxz_xxzz = cbuffer.data(fg_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_x_xxz_xyyy = cbuffer.data(fg_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_x_xxz_xyyz = cbuffer.data(fg_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_x_xxz_xyzz = cbuffer.data(fg_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_x_xxz_xzzz = cbuffer.data(fg_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_x_xxz_yyyy = cbuffer.data(fg_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_x_xxz_yyyz = cbuffer.data(fg_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_x_xxz_yyzz = cbuffer.data(fg_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_x_xxz_yzzz = cbuffer.data(fg_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_x_xxz_zzzz = cbuffer.data(fg_geom_11_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xx_xxxx, g_y_x_xx_xxxxz, g_y_x_xx_xxxy, g_y_x_xx_xxxyz, g_y_x_xx_xxxz, g_y_x_xx_xxxzz, g_y_x_xx_xxyy, g_y_x_xx_xxyyz, g_y_x_xx_xxyz, g_y_x_xx_xxyzz, g_y_x_xx_xxzz, g_y_x_xx_xxzzz, g_y_x_xx_xyyy, g_y_x_xx_xyyyz, g_y_x_xx_xyyz, g_y_x_xx_xyyzz, g_y_x_xx_xyzz, g_y_x_xx_xyzzz, g_y_x_xx_xzzz, g_y_x_xx_xzzzz, g_y_x_xx_yyyy, g_y_x_xx_yyyyz, g_y_x_xx_yyyz, g_y_x_xx_yyyzz, g_y_x_xx_yyzz, g_y_x_xx_yyzzz, g_y_x_xx_yzzz, g_y_x_xx_yzzzz, g_y_x_xx_zzzz, g_y_x_xx_zzzzz, g_y_x_xxz_xxxx, g_y_x_xxz_xxxy, g_y_x_xxz_xxxz, g_y_x_xxz_xxyy, g_y_x_xxz_xxyz, g_y_x_xxz_xxzz, g_y_x_xxz_xyyy, g_y_x_xxz_xyyz, g_y_x_xxz_xyzz, g_y_x_xxz_xzzz, g_y_x_xxz_yyyy, g_y_x_xxz_yyyz, g_y_x_xxz_yyzz, g_y_x_xxz_yzzz, g_y_x_xxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxz_xxxx[k] = -g_y_x_xx_xxxx[k] * ab_z + g_y_x_xx_xxxxz[k];

                g_y_x_xxz_xxxy[k] = -g_y_x_xx_xxxy[k] * ab_z + g_y_x_xx_xxxyz[k];

                g_y_x_xxz_xxxz[k] = -g_y_x_xx_xxxz[k] * ab_z + g_y_x_xx_xxxzz[k];

                g_y_x_xxz_xxyy[k] = -g_y_x_xx_xxyy[k] * ab_z + g_y_x_xx_xxyyz[k];

                g_y_x_xxz_xxyz[k] = -g_y_x_xx_xxyz[k] * ab_z + g_y_x_xx_xxyzz[k];

                g_y_x_xxz_xxzz[k] = -g_y_x_xx_xxzz[k] * ab_z + g_y_x_xx_xxzzz[k];

                g_y_x_xxz_xyyy[k] = -g_y_x_xx_xyyy[k] * ab_z + g_y_x_xx_xyyyz[k];

                g_y_x_xxz_xyyz[k] = -g_y_x_xx_xyyz[k] * ab_z + g_y_x_xx_xyyzz[k];

                g_y_x_xxz_xyzz[k] = -g_y_x_xx_xyzz[k] * ab_z + g_y_x_xx_xyzzz[k];

                g_y_x_xxz_xzzz[k] = -g_y_x_xx_xzzz[k] * ab_z + g_y_x_xx_xzzzz[k];

                g_y_x_xxz_yyyy[k] = -g_y_x_xx_yyyy[k] * ab_z + g_y_x_xx_yyyyz[k];

                g_y_x_xxz_yyyz[k] = -g_y_x_xx_yyyz[k] * ab_z + g_y_x_xx_yyyzz[k];

                g_y_x_xxz_yyzz[k] = -g_y_x_xx_yyzz[k] * ab_z + g_y_x_xx_yyzzz[k];

                g_y_x_xxz_yzzz[k] = -g_y_x_xx_yzzz[k] * ab_z + g_y_x_xx_yzzzz[k];

                g_y_x_xxz_zzzz[k] = -g_y_x_xx_zzzz[k] * ab_z + g_y_x_xx_zzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyy_xxxx = cbuffer.data(fg_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_x_xyy_xxxy = cbuffer.data(fg_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_x_xyy_xxxz = cbuffer.data(fg_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_x_xyy_xxyy = cbuffer.data(fg_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_x_xyy_xxyz = cbuffer.data(fg_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_x_xyy_xxzz = cbuffer.data(fg_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_x_xyy_xyyy = cbuffer.data(fg_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_x_xyy_xyyz = cbuffer.data(fg_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_x_xyy_xyzz = cbuffer.data(fg_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_x_xyy_xzzz = cbuffer.data(fg_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_x_xyy_yyyy = cbuffer.data(fg_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_x_xyy_yyyz = cbuffer.data(fg_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_x_xyy_yyzz = cbuffer.data(fg_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_x_xyy_yzzz = cbuffer.data(fg_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_x_xyy_zzzz = cbuffer.data(fg_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz, g_y_x_xyy_xxxx, g_y_x_xyy_xxxy, g_y_x_xyy_xxxz, g_y_x_xyy_xxyy, g_y_x_xyy_xxyz, g_y_x_xyy_xxzz, g_y_x_xyy_xyyy, g_y_x_xyy_xyyz, g_y_x_xyy_xyzz, g_y_x_xyy_xzzz, g_y_x_xyy_yyyy, g_y_x_xyy_yyyz, g_y_x_xyy_yyzz, g_y_x_xyy_yzzz, g_y_x_xyy_zzzz, g_y_x_yy_xxxx, g_y_x_yy_xxxxx, g_y_x_yy_xxxxy, g_y_x_yy_xxxxz, g_y_x_yy_xxxy, g_y_x_yy_xxxyy, g_y_x_yy_xxxyz, g_y_x_yy_xxxz, g_y_x_yy_xxxzz, g_y_x_yy_xxyy, g_y_x_yy_xxyyy, g_y_x_yy_xxyyz, g_y_x_yy_xxyz, g_y_x_yy_xxyzz, g_y_x_yy_xxzz, g_y_x_yy_xxzzz, g_y_x_yy_xyyy, g_y_x_yy_xyyyy, g_y_x_yy_xyyyz, g_y_x_yy_xyyz, g_y_x_yy_xyyzz, g_y_x_yy_xyzz, g_y_x_yy_xyzzz, g_y_x_yy_xzzz, g_y_x_yy_xzzzz, g_y_x_yy_yyyy, g_y_x_yy_yyyz, g_y_x_yy_yyzz, g_y_x_yy_yzzz, g_y_x_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyy_xxxx[k] = g_y_0_yy_xxxx[k] - g_y_x_yy_xxxx[k] * ab_x + g_y_x_yy_xxxxx[k];

                g_y_x_xyy_xxxy[k] = g_y_0_yy_xxxy[k] - g_y_x_yy_xxxy[k] * ab_x + g_y_x_yy_xxxxy[k];

                g_y_x_xyy_xxxz[k] = g_y_0_yy_xxxz[k] - g_y_x_yy_xxxz[k] * ab_x + g_y_x_yy_xxxxz[k];

                g_y_x_xyy_xxyy[k] = g_y_0_yy_xxyy[k] - g_y_x_yy_xxyy[k] * ab_x + g_y_x_yy_xxxyy[k];

                g_y_x_xyy_xxyz[k] = g_y_0_yy_xxyz[k] - g_y_x_yy_xxyz[k] * ab_x + g_y_x_yy_xxxyz[k];

                g_y_x_xyy_xxzz[k] = g_y_0_yy_xxzz[k] - g_y_x_yy_xxzz[k] * ab_x + g_y_x_yy_xxxzz[k];

                g_y_x_xyy_xyyy[k] = g_y_0_yy_xyyy[k] - g_y_x_yy_xyyy[k] * ab_x + g_y_x_yy_xxyyy[k];

                g_y_x_xyy_xyyz[k] = g_y_0_yy_xyyz[k] - g_y_x_yy_xyyz[k] * ab_x + g_y_x_yy_xxyyz[k];

                g_y_x_xyy_xyzz[k] = g_y_0_yy_xyzz[k] - g_y_x_yy_xyzz[k] * ab_x + g_y_x_yy_xxyzz[k];

                g_y_x_xyy_xzzz[k] = g_y_0_yy_xzzz[k] - g_y_x_yy_xzzz[k] * ab_x + g_y_x_yy_xxzzz[k];

                g_y_x_xyy_yyyy[k] = g_y_0_yy_yyyy[k] - g_y_x_yy_yyyy[k] * ab_x + g_y_x_yy_xyyyy[k];

                g_y_x_xyy_yyyz[k] = g_y_0_yy_yyyz[k] - g_y_x_yy_yyyz[k] * ab_x + g_y_x_yy_xyyyz[k];

                g_y_x_xyy_yyzz[k] = g_y_0_yy_yyzz[k] - g_y_x_yy_yyzz[k] * ab_x + g_y_x_yy_xyyzz[k];

                g_y_x_xyy_yzzz[k] = g_y_0_yy_yzzz[k] - g_y_x_yy_yzzz[k] * ab_x + g_y_x_yy_xyzzz[k];

                g_y_x_xyy_zzzz[k] = g_y_0_yy_zzzz[k] - g_y_x_yy_zzzz[k] * ab_x + g_y_x_yy_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyz_xxxx = cbuffer.data(fg_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_x_xyz_xxxy = cbuffer.data(fg_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_x_xyz_xxxz = cbuffer.data(fg_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_x_xyz_xxyy = cbuffer.data(fg_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_x_xyz_xxyz = cbuffer.data(fg_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_x_xyz_xxzz = cbuffer.data(fg_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_x_xyz_xyyy = cbuffer.data(fg_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_x_xyz_xyyz = cbuffer.data(fg_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_x_xyz_xyzz = cbuffer.data(fg_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_x_xyz_xzzz = cbuffer.data(fg_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_x_xyz_yyyy = cbuffer.data(fg_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_x_xyz_yyyz = cbuffer.data(fg_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_x_xyz_yyzz = cbuffer.data(fg_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_x_xyz_yzzz = cbuffer.data(fg_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_x_xyz_zzzz = cbuffer.data(fg_geom_11_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xy_xxxx, g_y_x_xy_xxxxz, g_y_x_xy_xxxy, g_y_x_xy_xxxyz, g_y_x_xy_xxxz, g_y_x_xy_xxxzz, g_y_x_xy_xxyy, g_y_x_xy_xxyyz, g_y_x_xy_xxyz, g_y_x_xy_xxyzz, g_y_x_xy_xxzz, g_y_x_xy_xxzzz, g_y_x_xy_xyyy, g_y_x_xy_xyyyz, g_y_x_xy_xyyz, g_y_x_xy_xyyzz, g_y_x_xy_xyzz, g_y_x_xy_xyzzz, g_y_x_xy_xzzz, g_y_x_xy_xzzzz, g_y_x_xy_yyyy, g_y_x_xy_yyyyz, g_y_x_xy_yyyz, g_y_x_xy_yyyzz, g_y_x_xy_yyzz, g_y_x_xy_yyzzz, g_y_x_xy_yzzz, g_y_x_xy_yzzzz, g_y_x_xy_zzzz, g_y_x_xy_zzzzz, g_y_x_xyz_xxxx, g_y_x_xyz_xxxy, g_y_x_xyz_xxxz, g_y_x_xyz_xxyy, g_y_x_xyz_xxyz, g_y_x_xyz_xxzz, g_y_x_xyz_xyyy, g_y_x_xyz_xyyz, g_y_x_xyz_xyzz, g_y_x_xyz_xzzz, g_y_x_xyz_yyyy, g_y_x_xyz_yyyz, g_y_x_xyz_yyzz, g_y_x_xyz_yzzz, g_y_x_xyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyz_xxxx[k] = -g_y_x_xy_xxxx[k] * ab_z + g_y_x_xy_xxxxz[k];

                g_y_x_xyz_xxxy[k] = -g_y_x_xy_xxxy[k] * ab_z + g_y_x_xy_xxxyz[k];

                g_y_x_xyz_xxxz[k] = -g_y_x_xy_xxxz[k] * ab_z + g_y_x_xy_xxxzz[k];

                g_y_x_xyz_xxyy[k] = -g_y_x_xy_xxyy[k] * ab_z + g_y_x_xy_xxyyz[k];

                g_y_x_xyz_xxyz[k] = -g_y_x_xy_xxyz[k] * ab_z + g_y_x_xy_xxyzz[k];

                g_y_x_xyz_xxzz[k] = -g_y_x_xy_xxzz[k] * ab_z + g_y_x_xy_xxzzz[k];

                g_y_x_xyz_xyyy[k] = -g_y_x_xy_xyyy[k] * ab_z + g_y_x_xy_xyyyz[k];

                g_y_x_xyz_xyyz[k] = -g_y_x_xy_xyyz[k] * ab_z + g_y_x_xy_xyyzz[k];

                g_y_x_xyz_xyzz[k] = -g_y_x_xy_xyzz[k] * ab_z + g_y_x_xy_xyzzz[k];

                g_y_x_xyz_xzzz[k] = -g_y_x_xy_xzzz[k] * ab_z + g_y_x_xy_xzzzz[k];

                g_y_x_xyz_yyyy[k] = -g_y_x_xy_yyyy[k] * ab_z + g_y_x_xy_yyyyz[k];

                g_y_x_xyz_yyyz[k] = -g_y_x_xy_yyyz[k] * ab_z + g_y_x_xy_yyyzz[k];

                g_y_x_xyz_yyzz[k] = -g_y_x_xy_yyzz[k] * ab_z + g_y_x_xy_yyzzz[k];

                g_y_x_xyz_yzzz[k] = -g_y_x_xy_yzzz[k] * ab_z + g_y_x_xy_yzzzz[k];

                g_y_x_xyz_zzzz[k] = -g_y_x_xy_zzzz[k] * ab_z + g_y_x_xy_zzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzz_xxxx = cbuffer.data(fg_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_x_xzz_xxxy = cbuffer.data(fg_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_x_xzz_xxxz = cbuffer.data(fg_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_x_xzz_xxyy = cbuffer.data(fg_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_x_xzz_xxyz = cbuffer.data(fg_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_x_xzz_xxzz = cbuffer.data(fg_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_x_xzz_xyyy = cbuffer.data(fg_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_x_xzz_xyyz = cbuffer.data(fg_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_x_xzz_xyzz = cbuffer.data(fg_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_x_xzz_xzzz = cbuffer.data(fg_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_x_xzz_yyyy = cbuffer.data(fg_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_x_xzz_yyyz = cbuffer.data(fg_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_x_xzz_yyzz = cbuffer.data(fg_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_x_xzz_yzzz = cbuffer.data(fg_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_x_xzz_zzzz = cbuffer.data(fg_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xz_xxxx, g_y_x_xz_xxxxz, g_y_x_xz_xxxy, g_y_x_xz_xxxyz, g_y_x_xz_xxxz, g_y_x_xz_xxxzz, g_y_x_xz_xxyy, g_y_x_xz_xxyyz, g_y_x_xz_xxyz, g_y_x_xz_xxyzz, g_y_x_xz_xxzz, g_y_x_xz_xxzzz, g_y_x_xz_xyyy, g_y_x_xz_xyyyz, g_y_x_xz_xyyz, g_y_x_xz_xyyzz, g_y_x_xz_xyzz, g_y_x_xz_xyzzz, g_y_x_xz_xzzz, g_y_x_xz_xzzzz, g_y_x_xz_yyyy, g_y_x_xz_yyyyz, g_y_x_xz_yyyz, g_y_x_xz_yyyzz, g_y_x_xz_yyzz, g_y_x_xz_yyzzz, g_y_x_xz_yzzz, g_y_x_xz_yzzzz, g_y_x_xz_zzzz, g_y_x_xz_zzzzz, g_y_x_xzz_xxxx, g_y_x_xzz_xxxy, g_y_x_xzz_xxxz, g_y_x_xzz_xxyy, g_y_x_xzz_xxyz, g_y_x_xzz_xxzz, g_y_x_xzz_xyyy, g_y_x_xzz_xyyz, g_y_x_xzz_xyzz, g_y_x_xzz_xzzz, g_y_x_xzz_yyyy, g_y_x_xzz_yyyz, g_y_x_xzz_yyzz, g_y_x_xzz_yzzz, g_y_x_xzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzz_xxxx[k] = -g_y_x_xz_xxxx[k] * ab_z + g_y_x_xz_xxxxz[k];

                g_y_x_xzz_xxxy[k] = -g_y_x_xz_xxxy[k] * ab_z + g_y_x_xz_xxxyz[k];

                g_y_x_xzz_xxxz[k] = -g_y_x_xz_xxxz[k] * ab_z + g_y_x_xz_xxxzz[k];

                g_y_x_xzz_xxyy[k] = -g_y_x_xz_xxyy[k] * ab_z + g_y_x_xz_xxyyz[k];

                g_y_x_xzz_xxyz[k] = -g_y_x_xz_xxyz[k] * ab_z + g_y_x_xz_xxyzz[k];

                g_y_x_xzz_xxzz[k] = -g_y_x_xz_xxzz[k] * ab_z + g_y_x_xz_xxzzz[k];

                g_y_x_xzz_xyyy[k] = -g_y_x_xz_xyyy[k] * ab_z + g_y_x_xz_xyyyz[k];

                g_y_x_xzz_xyyz[k] = -g_y_x_xz_xyyz[k] * ab_z + g_y_x_xz_xyyzz[k];

                g_y_x_xzz_xyzz[k] = -g_y_x_xz_xyzz[k] * ab_z + g_y_x_xz_xyzzz[k];

                g_y_x_xzz_xzzz[k] = -g_y_x_xz_xzzz[k] * ab_z + g_y_x_xz_xzzzz[k];

                g_y_x_xzz_yyyy[k] = -g_y_x_xz_yyyy[k] * ab_z + g_y_x_xz_yyyyz[k];

                g_y_x_xzz_yyyz[k] = -g_y_x_xz_yyyz[k] * ab_z + g_y_x_xz_yyyzz[k];

                g_y_x_xzz_yyzz[k] = -g_y_x_xz_yyzz[k] * ab_z + g_y_x_xz_yyzzz[k];

                g_y_x_xzz_yzzz[k] = -g_y_x_xz_yzzz[k] * ab_z + g_y_x_xz_yzzzz[k];

                g_y_x_xzz_zzzz[k] = -g_y_x_xz_zzzz[k] * ab_z + g_y_x_xz_zzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyy_xxxx = cbuffer.data(fg_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_x_yyy_xxxy = cbuffer.data(fg_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_x_yyy_xxxz = cbuffer.data(fg_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_x_yyy_xxyy = cbuffer.data(fg_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_x_yyy_xxyz = cbuffer.data(fg_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_x_yyy_xxzz = cbuffer.data(fg_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_x_yyy_xyyy = cbuffer.data(fg_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_x_yyy_xyyz = cbuffer.data(fg_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_x_yyy_xyzz = cbuffer.data(fg_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_x_yyy_xzzz = cbuffer.data(fg_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_x_yyy_yyyy = cbuffer.data(fg_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_x_yyy_yyyz = cbuffer.data(fg_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_x_yyy_yyzz = cbuffer.data(fg_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_x_yyy_yzzz = cbuffer.data(fg_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_x_yyy_zzzz = cbuffer.data(fg_geom_11_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_xxxx, g_0_x_yy_xxxy, g_0_x_yy_xxxz, g_0_x_yy_xxyy, g_0_x_yy_xxyz, g_0_x_yy_xxzz, g_0_x_yy_xyyy, g_0_x_yy_xyyz, g_0_x_yy_xyzz, g_0_x_yy_xzzz, g_0_x_yy_yyyy, g_0_x_yy_yyyz, g_0_x_yy_yyzz, g_0_x_yy_yzzz, g_0_x_yy_zzzz, g_y_x_yy_xxxx, g_y_x_yy_xxxxy, g_y_x_yy_xxxy, g_y_x_yy_xxxyy, g_y_x_yy_xxxyz, g_y_x_yy_xxxz, g_y_x_yy_xxyy, g_y_x_yy_xxyyy, g_y_x_yy_xxyyz, g_y_x_yy_xxyz, g_y_x_yy_xxyzz, g_y_x_yy_xxzz, g_y_x_yy_xyyy, g_y_x_yy_xyyyy, g_y_x_yy_xyyyz, g_y_x_yy_xyyz, g_y_x_yy_xyyzz, g_y_x_yy_xyzz, g_y_x_yy_xyzzz, g_y_x_yy_xzzz, g_y_x_yy_yyyy, g_y_x_yy_yyyyy, g_y_x_yy_yyyyz, g_y_x_yy_yyyz, g_y_x_yy_yyyzz, g_y_x_yy_yyzz, g_y_x_yy_yyzzz, g_y_x_yy_yzzz, g_y_x_yy_yzzzz, g_y_x_yy_zzzz, g_y_x_yyy_xxxx, g_y_x_yyy_xxxy, g_y_x_yyy_xxxz, g_y_x_yyy_xxyy, g_y_x_yyy_xxyz, g_y_x_yyy_xxzz, g_y_x_yyy_xyyy, g_y_x_yyy_xyyz, g_y_x_yyy_xyzz, g_y_x_yyy_xzzz, g_y_x_yyy_yyyy, g_y_x_yyy_yyyz, g_y_x_yyy_yyzz, g_y_x_yyy_yzzz, g_y_x_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyy_xxxx[k] = -g_0_x_yy_xxxx[k] - g_y_x_yy_xxxx[k] * ab_y + g_y_x_yy_xxxxy[k];

                g_y_x_yyy_xxxy[k] = -g_0_x_yy_xxxy[k] - g_y_x_yy_xxxy[k] * ab_y + g_y_x_yy_xxxyy[k];

                g_y_x_yyy_xxxz[k] = -g_0_x_yy_xxxz[k] - g_y_x_yy_xxxz[k] * ab_y + g_y_x_yy_xxxyz[k];

                g_y_x_yyy_xxyy[k] = -g_0_x_yy_xxyy[k] - g_y_x_yy_xxyy[k] * ab_y + g_y_x_yy_xxyyy[k];

                g_y_x_yyy_xxyz[k] = -g_0_x_yy_xxyz[k] - g_y_x_yy_xxyz[k] * ab_y + g_y_x_yy_xxyyz[k];

                g_y_x_yyy_xxzz[k] = -g_0_x_yy_xxzz[k] - g_y_x_yy_xxzz[k] * ab_y + g_y_x_yy_xxyzz[k];

                g_y_x_yyy_xyyy[k] = -g_0_x_yy_xyyy[k] - g_y_x_yy_xyyy[k] * ab_y + g_y_x_yy_xyyyy[k];

                g_y_x_yyy_xyyz[k] = -g_0_x_yy_xyyz[k] - g_y_x_yy_xyyz[k] * ab_y + g_y_x_yy_xyyyz[k];

                g_y_x_yyy_xyzz[k] = -g_0_x_yy_xyzz[k] - g_y_x_yy_xyzz[k] * ab_y + g_y_x_yy_xyyzz[k];

                g_y_x_yyy_xzzz[k] = -g_0_x_yy_xzzz[k] - g_y_x_yy_xzzz[k] * ab_y + g_y_x_yy_xyzzz[k];

                g_y_x_yyy_yyyy[k] = -g_0_x_yy_yyyy[k] - g_y_x_yy_yyyy[k] * ab_y + g_y_x_yy_yyyyy[k];

                g_y_x_yyy_yyyz[k] = -g_0_x_yy_yyyz[k] - g_y_x_yy_yyyz[k] * ab_y + g_y_x_yy_yyyyz[k];

                g_y_x_yyy_yyzz[k] = -g_0_x_yy_yyzz[k] - g_y_x_yy_yyzz[k] * ab_y + g_y_x_yy_yyyzz[k];

                g_y_x_yyy_yzzz[k] = -g_0_x_yy_yzzz[k] - g_y_x_yy_yzzz[k] * ab_y + g_y_x_yy_yyzzz[k];

                g_y_x_yyy_zzzz[k] = -g_0_x_yy_zzzz[k] - g_y_x_yy_zzzz[k] * ab_y + g_y_x_yy_yzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyz_xxxx = cbuffer.data(fg_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_x_yyz_xxxy = cbuffer.data(fg_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_x_yyz_xxxz = cbuffer.data(fg_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_x_yyz_xxyy = cbuffer.data(fg_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_x_yyz_xxyz = cbuffer.data(fg_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_x_yyz_xxzz = cbuffer.data(fg_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_x_yyz_xyyy = cbuffer.data(fg_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_x_yyz_xyyz = cbuffer.data(fg_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_x_yyz_xyzz = cbuffer.data(fg_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_x_yyz_xzzz = cbuffer.data(fg_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_x_yyz_yyyy = cbuffer.data(fg_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_x_yyz_yyyz = cbuffer.data(fg_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_x_yyz_yyzz = cbuffer.data(fg_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_x_yyz_yzzz = cbuffer.data(fg_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_x_yyz_zzzz = cbuffer.data(fg_geom_11_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yy_xxxx, g_y_x_yy_xxxxz, g_y_x_yy_xxxy, g_y_x_yy_xxxyz, g_y_x_yy_xxxz, g_y_x_yy_xxxzz, g_y_x_yy_xxyy, g_y_x_yy_xxyyz, g_y_x_yy_xxyz, g_y_x_yy_xxyzz, g_y_x_yy_xxzz, g_y_x_yy_xxzzz, g_y_x_yy_xyyy, g_y_x_yy_xyyyz, g_y_x_yy_xyyz, g_y_x_yy_xyyzz, g_y_x_yy_xyzz, g_y_x_yy_xyzzz, g_y_x_yy_xzzz, g_y_x_yy_xzzzz, g_y_x_yy_yyyy, g_y_x_yy_yyyyz, g_y_x_yy_yyyz, g_y_x_yy_yyyzz, g_y_x_yy_yyzz, g_y_x_yy_yyzzz, g_y_x_yy_yzzz, g_y_x_yy_yzzzz, g_y_x_yy_zzzz, g_y_x_yy_zzzzz, g_y_x_yyz_xxxx, g_y_x_yyz_xxxy, g_y_x_yyz_xxxz, g_y_x_yyz_xxyy, g_y_x_yyz_xxyz, g_y_x_yyz_xxzz, g_y_x_yyz_xyyy, g_y_x_yyz_xyyz, g_y_x_yyz_xyzz, g_y_x_yyz_xzzz, g_y_x_yyz_yyyy, g_y_x_yyz_yyyz, g_y_x_yyz_yyzz, g_y_x_yyz_yzzz, g_y_x_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyz_xxxx[k] = -g_y_x_yy_xxxx[k] * ab_z + g_y_x_yy_xxxxz[k];

                g_y_x_yyz_xxxy[k] = -g_y_x_yy_xxxy[k] * ab_z + g_y_x_yy_xxxyz[k];

                g_y_x_yyz_xxxz[k] = -g_y_x_yy_xxxz[k] * ab_z + g_y_x_yy_xxxzz[k];

                g_y_x_yyz_xxyy[k] = -g_y_x_yy_xxyy[k] * ab_z + g_y_x_yy_xxyyz[k];

                g_y_x_yyz_xxyz[k] = -g_y_x_yy_xxyz[k] * ab_z + g_y_x_yy_xxyzz[k];

                g_y_x_yyz_xxzz[k] = -g_y_x_yy_xxzz[k] * ab_z + g_y_x_yy_xxzzz[k];

                g_y_x_yyz_xyyy[k] = -g_y_x_yy_xyyy[k] * ab_z + g_y_x_yy_xyyyz[k];

                g_y_x_yyz_xyyz[k] = -g_y_x_yy_xyyz[k] * ab_z + g_y_x_yy_xyyzz[k];

                g_y_x_yyz_xyzz[k] = -g_y_x_yy_xyzz[k] * ab_z + g_y_x_yy_xyzzz[k];

                g_y_x_yyz_xzzz[k] = -g_y_x_yy_xzzz[k] * ab_z + g_y_x_yy_xzzzz[k];

                g_y_x_yyz_yyyy[k] = -g_y_x_yy_yyyy[k] * ab_z + g_y_x_yy_yyyyz[k];

                g_y_x_yyz_yyyz[k] = -g_y_x_yy_yyyz[k] * ab_z + g_y_x_yy_yyyzz[k];

                g_y_x_yyz_yyzz[k] = -g_y_x_yy_yyzz[k] * ab_z + g_y_x_yy_yyzzz[k];

                g_y_x_yyz_yzzz[k] = -g_y_x_yy_yzzz[k] * ab_z + g_y_x_yy_yzzzz[k];

                g_y_x_yyz_zzzz[k] = -g_y_x_yy_zzzz[k] * ab_z + g_y_x_yy_zzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzz_xxxx = cbuffer.data(fg_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_x_yzz_xxxy = cbuffer.data(fg_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_x_yzz_xxxz = cbuffer.data(fg_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_x_yzz_xxyy = cbuffer.data(fg_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_x_yzz_xxyz = cbuffer.data(fg_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_x_yzz_xxzz = cbuffer.data(fg_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_x_yzz_xyyy = cbuffer.data(fg_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_x_yzz_xyyz = cbuffer.data(fg_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_x_yzz_xyzz = cbuffer.data(fg_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_x_yzz_xzzz = cbuffer.data(fg_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_x_yzz_yyyy = cbuffer.data(fg_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_x_yzz_yyyz = cbuffer.data(fg_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_x_yzz_yyzz = cbuffer.data(fg_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_x_yzz_yzzz = cbuffer.data(fg_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_x_yzz_zzzz = cbuffer.data(fg_geom_11_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yz_xxxx, g_y_x_yz_xxxxz, g_y_x_yz_xxxy, g_y_x_yz_xxxyz, g_y_x_yz_xxxz, g_y_x_yz_xxxzz, g_y_x_yz_xxyy, g_y_x_yz_xxyyz, g_y_x_yz_xxyz, g_y_x_yz_xxyzz, g_y_x_yz_xxzz, g_y_x_yz_xxzzz, g_y_x_yz_xyyy, g_y_x_yz_xyyyz, g_y_x_yz_xyyz, g_y_x_yz_xyyzz, g_y_x_yz_xyzz, g_y_x_yz_xyzzz, g_y_x_yz_xzzz, g_y_x_yz_xzzzz, g_y_x_yz_yyyy, g_y_x_yz_yyyyz, g_y_x_yz_yyyz, g_y_x_yz_yyyzz, g_y_x_yz_yyzz, g_y_x_yz_yyzzz, g_y_x_yz_yzzz, g_y_x_yz_yzzzz, g_y_x_yz_zzzz, g_y_x_yz_zzzzz, g_y_x_yzz_xxxx, g_y_x_yzz_xxxy, g_y_x_yzz_xxxz, g_y_x_yzz_xxyy, g_y_x_yzz_xxyz, g_y_x_yzz_xxzz, g_y_x_yzz_xyyy, g_y_x_yzz_xyyz, g_y_x_yzz_xyzz, g_y_x_yzz_xzzz, g_y_x_yzz_yyyy, g_y_x_yzz_yyyz, g_y_x_yzz_yyzz, g_y_x_yzz_yzzz, g_y_x_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzz_xxxx[k] = -g_y_x_yz_xxxx[k] * ab_z + g_y_x_yz_xxxxz[k];

                g_y_x_yzz_xxxy[k] = -g_y_x_yz_xxxy[k] * ab_z + g_y_x_yz_xxxyz[k];

                g_y_x_yzz_xxxz[k] = -g_y_x_yz_xxxz[k] * ab_z + g_y_x_yz_xxxzz[k];

                g_y_x_yzz_xxyy[k] = -g_y_x_yz_xxyy[k] * ab_z + g_y_x_yz_xxyyz[k];

                g_y_x_yzz_xxyz[k] = -g_y_x_yz_xxyz[k] * ab_z + g_y_x_yz_xxyzz[k];

                g_y_x_yzz_xxzz[k] = -g_y_x_yz_xxzz[k] * ab_z + g_y_x_yz_xxzzz[k];

                g_y_x_yzz_xyyy[k] = -g_y_x_yz_xyyy[k] * ab_z + g_y_x_yz_xyyyz[k];

                g_y_x_yzz_xyyz[k] = -g_y_x_yz_xyyz[k] * ab_z + g_y_x_yz_xyyzz[k];

                g_y_x_yzz_xyzz[k] = -g_y_x_yz_xyzz[k] * ab_z + g_y_x_yz_xyzzz[k];

                g_y_x_yzz_xzzz[k] = -g_y_x_yz_xzzz[k] * ab_z + g_y_x_yz_xzzzz[k];

                g_y_x_yzz_yyyy[k] = -g_y_x_yz_yyyy[k] * ab_z + g_y_x_yz_yyyyz[k];

                g_y_x_yzz_yyyz[k] = -g_y_x_yz_yyyz[k] * ab_z + g_y_x_yz_yyyzz[k];

                g_y_x_yzz_yyzz[k] = -g_y_x_yz_yyzz[k] * ab_z + g_y_x_yz_yyzzz[k];

                g_y_x_yzz_yzzz[k] = -g_y_x_yz_yzzz[k] * ab_z + g_y_x_yz_yzzzz[k];

                g_y_x_yzz_zzzz[k] = -g_y_x_yz_zzzz[k] * ab_z + g_y_x_yz_zzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzz_xxxx = cbuffer.data(fg_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_x_zzz_xxxy = cbuffer.data(fg_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_x_zzz_xxxz = cbuffer.data(fg_geom_11_off + 587 * ccomps * dcomps);

            auto g_y_x_zzz_xxyy = cbuffer.data(fg_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_x_zzz_xxyz = cbuffer.data(fg_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_x_zzz_xxzz = cbuffer.data(fg_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_x_zzz_xyyy = cbuffer.data(fg_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_x_zzz_xyyz = cbuffer.data(fg_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_x_zzz_xyzz = cbuffer.data(fg_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_x_zzz_xzzz = cbuffer.data(fg_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_x_zzz_yyyy = cbuffer.data(fg_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_x_zzz_yyyz = cbuffer.data(fg_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_x_zzz_yyzz = cbuffer.data(fg_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_x_zzz_yzzz = cbuffer.data(fg_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_x_zzz_zzzz = cbuffer.data(fg_geom_11_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zz_xxxx, g_y_x_zz_xxxxz, g_y_x_zz_xxxy, g_y_x_zz_xxxyz, g_y_x_zz_xxxz, g_y_x_zz_xxxzz, g_y_x_zz_xxyy, g_y_x_zz_xxyyz, g_y_x_zz_xxyz, g_y_x_zz_xxyzz, g_y_x_zz_xxzz, g_y_x_zz_xxzzz, g_y_x_zz_xyyy, g_y_x_zz_xyyyz, g_y_x_zz_xyyz, g_y_x_zz_xyyzz, g_y_x_zz_xyzz, g_y_x_zz_xyzzz, g_y_x_zz_xzzz, g_y_x_zz_xzzzz, g_y_x_zz_yyyy, g_y_x_zz_yyyyz, g_y_x_zz_yyyz, g_y_x_zz_yyyzz, g_y_x_zz_yyzz, g_y_x_zz_yyzzz, g_y_x_zz_yzzz, g_y_x_zz_yzzzz, g_y_x_zz_zzzz, g_y_x_zz_zzzzz, g_y_x_zzz_xxxx, g_y_x_zzz_xxxy, g_y_x_zzz_xxxz, g_y_x_zzz_xxyy, g_y_x_zzz_xxyz, g_y_x_zzz_xxzz, g_y_x_zzz_xyyy, g_y_x_zzz_xyyz, g_y_x_zzz_xyzz, g_y_x_zzz_xzzz, g_y_x_zzz_yyyy, g_y_x_zzz_yyyz, g_y_x_zzz_yyzz, g_y_x_zzz_yzzz, g_y_x_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzz_xxxx[k] = -g_y_x_zz_xxxx[k] * ab_z + g_y_x_zz_xxxxz[k];

                g_y_x_zzz_xxxy[k] = -g_y_x_zz_xxxy[k] * ab_z + g_y_x_zz_xxxyz[k];

                g_y_x_zzz_xxxz[k] = -g_y_x_zz_xxxz[k] * ab_z + g_y_x_zz_xxxzz[k];

                g_y_x_zzz_xxyy[k] = -g_y_x_zz_xxyy[k] * ab_z + g_y_x_zz_xxyyz[k];

                g_y_x_zzz_xxyz[k] = -g_y_x_zz_xxyz[k] * ab_z + g_y_x_zz_xxyzz[k];

                g_y_x_zzz_xxzz[k] = -g_y_x_zz_xxzz[k] * ab_z + g_y_x_zz_xxzzz[k];

                g_y_x_zzz_xyyy[k] = -g_y_x_zz_xyyy[k] * ab_z + g_y_x_zz_xyyyz[k];

                g_y_x_zzz_xyyz[k] = -g_y_x_zz_xyyz[k] * ab_z + g_y_x_zz_xyyzz[k];

                g_y_x_zzz_xyzz[k] = -g_y_x_zz_xyzz[k] * ab_z + g_y_x_zz_xyzzz[k];

                g_y_x_zzz_xzzz[k] = -g_y_x_zz_xzzz[k] * ab_z + g_y_x_zz_xzzzz[k];

                g_y_x_zzz_yyyy[k] = -g_y_x_zz_yyyy[k] * ab_z + g_y_x_zz_yyyyz[k];

                g_y_x_zzz_yyyz[k] = -g_y_x_zz_yyyz[k] * ab_z + g_y_x_zz_yyyzz[k];

                g_y_x_zzz_yyzz[k] = -g_y_x_zz_yyzz[k] * ab_z + g_y_x_zz_yyzzz[k];

                g_y_x_zzz_yzzz[k] = -g_y_x_zz_yzzz[k] * ab_z + g_y_x_zz_yzzzz[k];

                g_y_x_zzz_zzzz[k] = -g_y_x_zz_zzzz[k] * ab_z + g_y_x_zz_zzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxx_xxxx = cbuffer.data(fg_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_y_xxx_xxxy = cbuffer.data(fg_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_y_xxx_xxxz = cbuffer.data(fg_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_y_xxx_xxyy = cbuffer.data(fg_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_y_xxx_xxyz = cbuffer.data(fg_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_y_xxx_xxzz = cbuffer.data(fg_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_y_xxx_xyyy = cbuffer.data(fg_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_y_xxx_xyyz = cbuffer.data(fg_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_y_xxx_xyzz = cbuffer.data(fg_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_y_xxx_xzzz = cbuffer.data(fg_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_y_xxx_yyyy = cbuffer.data(fg_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_y_xxx_yyyz = cbuffer.data(fg_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_y_xxx_yyzz = cbuffer.data(fg_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_y_xxx_yzzz = cbuffer.data(fg_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_y_xxx_zzzz = cbuffer.data(fg_geom_11_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xx_xxxx, g_y_y_xx_xxxxx, g_y_y_xx_xxxxy, g_y_y_xx_xxxxz, g_y_y_xx_xxxy, g_y_y_xx_xxxyy, g_y_y_xx_xxxyz, g_y_y_xx_xxxz, g_y_y_xx_xxxzz, g_y_y_xx_xxyy, g_y_y_xx_xxyyy, g_y_y_xx_xxyyz, g_y_y_xx_xxyz, g_y_y_xx_xxyzz, g_y_y_xx_xxzz, g_y_y_xx_xxzzz, g_y_y_xx_xyyy, g_y_y_xx_xyyyy, g_y_y_xx_xyyyz, g_y_y_xx_xyyz, g_y_y_xx_xyyzz, g_y_y_xx_xyzz, g_y_y_xx_xyzzz, g_y_y_xx_xzzz, g_y_y_xx_xzzzz, g_y_y_xx_yyyy, g_y_y_xx_yyyz, g_y_y_xx_yyzz, g_y_y_xx_yzzz, g_y_y_xx_zzzz, g_y_y_xxx_xxxx, g_y_y_xxx_xxxy, g_y_y_xxx_xxxz, g_y_y_xxx_xxyy, g_y_y_xxx_xxyz, g_y_y_xxx_xxzz, g_y_y_xxx_xyyy, g_y_y_xxx_xyyz, g_y_y_xxx_xyzz, g_y_y_xxx_xzzz, g_y_y_xxx_yyyy, g_y_y_xxx_yyyz, g_y_y_xxx_yyzz, g_y_y_xxx_yzzz, g_y_y_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxx_xxxx[k] = -g_y_y_xx_xxxx[k] * ab_x + g_y_y_xx_xxxxx[k];

                g_y_y_xxx_xxxy[k] = -g_y_y_xx_xxxy[k] * ab_x + g_y_y_xx_xxxxy[k];

                g_y_y_xxx_xxxz[k] = -g_y_y_xx_xxxz[k] * ab_x + g_y_y_xx_xxxxz[k];

                g_y_y_xxx_xxyy[k] = -g_y_y_xx_xxyy[k] * ab_x + g_y_y_xx_xxxyy[k];

                g_y_y_xxx_xxyz[k] = -g_y_y_xx_xxyz[k] * ab_x + g_y_y_xx_xxxyz[k];

                g_y_y_xxx_xxzz[k] = -g_y_y_xx_xxzz[k] * ab_x + g_y_y_xx_xxxzz[k];

                g_y_y_xxx_xyyy[k] = -g_y_y_xx_xyyy[k] * ab_x + g_y_y_xx_xxyyy[k];

                g_y_y_xxx_xyyz[k] = -g_y_y_xx_xyyz[k] * ab_x + g_y_y_xx_xxyyz[k];

                g_y_y_xxx_xyzz[k] = -g_y_y_xx_xyzz[k] * ab_x + g_y_y_xx_xxyzz[k];

                g_y_y_xxx_xzzz[k] = -g_y_y_xx_xzzz[k] * ab_x + g_y_y_xx_xxzzz[k];

                g_y_y_xxx_yyyy[k] = -g_y_y_xx_yyyy[k] * ab_x + g_y_y_xx_xyyyy[k];

                g_y_y_xxx_yyyz[k] = -g_y_y_xx_yyyz[k] * ab_x + g_y_y_xx_xyyyz[k];

                g_y_y_xxx_yyzz[k] = -g_y_y_xx_yyzz[k] * ab_x + g_y_y_xx_xyyzz[k];

                g_y_y_xxx_yzzz[k] = -g_y_y_xx_yzzz[k] * ab_x + g_y_y_xx_xyzzz[k];

                g_y_y_xxx_zzzz[k] = -g_y_y_xx_zzzz[k] * ab_x + g_y_y_xx_xzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxy_xxxx = cbuffer.data(fg_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_y_xxy_xxxy = cbuffer.data(fg_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_y_xxy_xxxz = cbuffer.data(fg_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_y_xxy_xxyy = cbuffer.data(fg_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_y_xxy_xxyz = cbuffer.data(fg_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_y_xxy_xxzz = cbuffer.data(fg_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_y_xxy_xyyy = cbuffer.data(fg_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_y_xxy_xyyz = cbuffer.data(fg_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_y_xxy_xyzz = cbuffer.data(fg_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_y_xxy_xzzz = cbuffer.data(fg_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_y_xxy_yyyy = cbuffer.data(fg_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_y_xxy_yyyz = cbuffer.data(fg_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_y_xxy_yyzz = cbuffer.data(fg_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_y_xxy_yzzz = cbuffer.data(fg_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_y_xxy_zzzz = cbuffer.data(fg_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxy_xxxx, g_y_y_xxy_xxxy, g_y_y_xxy_xxxz, g_y_y_xxy_xxyy, g_y_y_xxy_xxyz, g_y_y_xxy_xxzz, g_y_y_xxy_xyyy, g_y_y_xxy_xyyz, g_y_y_xxy_xyzz, g_y_y_xxy_xzzz, g_y_y_xxy_yyyy, g_y_y_xxy_yyyz, g_y_y_xxy_yyzz, g_y_y_xxy_yzzz, g_y_y_xxy_zzzz, g_y_y_xy_xxxx, g_y_y_xy_xxxxx, g_y_y_xy_xxxxy, g_y_y_xy_xxxxz, g_y_y_xy_xxxy, g_y_y_xy_xxxyy, g_y_y_xy_xxxyz, g_y_y_xy_xxxz, g_y_y_xy_xxxzz, g_y_y_xy_xxyy, g_y_y_xy_xxyyy, g_y_y_xy_xxyyz, g_y_y_xy_xxyz, g_y_y_xy_xxyzz, g_y_y_xy_xxzz, g_y_y_xy_xxzzz, g_y_y_xy_xyyy, g_y_y_xy_xyyyy, g_y_y_xy_xyyyz, g_y_y_xy_xyyz, g_y_y_xy_xyyzz, g_y_y_xy_xyzz, g_y_y_xy_xyzzz, g_y_y_xy_xzzz, g_y_y_xy_xzzzz, g_y_y_xy_yyyy, g_y_y_xy_yyyz, g_y_y_xy_yyzz, g_y_y_xy_yzzz, g_y_y_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxy_xxxx[k] = -g_y_y_xy_xxxx[k] * ab_x + g_y_y_xy_xxxxx[k];

                g_y_y_xxy_xxxy[k] = -g_y_y_xy_xxxy[k] * ab_x + g_y_y_xy_xxxxy[k];

                g_y_y_xxy_xxxz[k] = -g_y_y_xy_xxxz[k] * ab_x + g_y_y_xy_xxxxz[k];

                g_y_y_xxy_xxyy[k] = -g_y_y_xy_xxyy[k] * ab_x + g_y_y_xy_xxxyy[k];

                g_y_y_xxy_xxyz[k] = -g_y_y_xy_xxyz[k] * ab_x + g_y_y_xy_xxxyz[k];

                g_y_y_xxy_xxzz[k] = -g_y_y_xy_xxzz[k] * ab_x + g_y_y_xy_xxxzz[k];

                g_y_y_xxy_xyyy[k] = -g_y_y_xy_xyyy[k] * ab_x + g_y_y_xy_xxyyy[k];

                g_y_y_xxy_xyyz[k] = -g_y_y_xy_xyyz[k] * ab_x + g_y_y_xy_xxyyz[k];

                g_y_y_xxy_xyzz[k] = -g_y_y_xy_xyzz[k] * ab_x + g_y_y_xy_xxyzz[k];

                g_y_y_xxy_xzzz[k] = -g_y_y_xy_xzzz[k] * ab_x + g_y_y_xy_xxzzz[k];

                g_y_y_xxy_yyyy[k] = -g_y_y_xy_yyyy[k] * ab_x + g_y_y_xy_xyyyy[k];

                g_y_y_xxy_yyyz[k] = -g_y_y_xy_yyyz[k] * ab_x + g_y_y_xy_xyyyz[k];

                g_y_y_xxy_yyzz[k] = -g_y_y_xy_yyzz[k] * ab_x + g_y_y_xy_xyyzz[k];

                g_y_y_xxy_yzzz[k] = -g_y_y_xy_yzzz[k] * ab_x + g_y_y_xy_xyzzz[k];

                g_y_y_xxy_zzzz[k] = -g_y_y_xy_zzzz[k] * ab_x + g_y_y_xy_xzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxz_xxxx = cbuffer.data(fg_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_y_xxz_xxxy = cbuffer.data(fg_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_y_xxz_xxxz = cbuffer.data(fg_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_y_xxz_xxyy = cbuffer.data(fg_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_y_xxz_xxyz = cbuffer.data(fg_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_y_xxz_xxzz = cbuffer.data(fg_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_y_xxz_xyyy = cbuffer.data(fg_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_y_xxz_xyyz = cbuffer.data(fg_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_y_xxz_xyzz = cbuffer.data(fg_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_y_xxz_xzzz = cbuffer.data(fg_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_y_xxz_yyyy = cbuffer.data(fg_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_y_xxz_yyyz = cbuffer.data(fg_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_y_xxz_yyzz = cbuffer.data(fg_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_y_xxz_yzzz = cbuffer.data(fg_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_y_xxz_zzzz = cbuffer.data(fg_geom_11_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxz_xxxx, g_y_y_xxz_xxxy, g_y_y_xxz_xxxz, g_y_y_xxz_xxyy, g_y_y_xxz_xxyz, g_y_y_xxz_xxzz, g_y_y_xxz_xyyy, g_y_y_xxz_xyyz, g_y_y_xxz_xyzz, g_y_y_xxz_xzzz, g_y_y_xxz_yyyy, g_y_y_xxz_yyyz, g_y_y_xxz_yyzz, g_y_y_xxz_yzzz, g_y_y_xxz_zzzz, g_y_y_xz_xxxx, g_y_y_xz_xxxxx, g_y_y_xz_xxxxy, g_y_y_xz_xxxxz, g_y_y_xz_xxxy, g_y_y_xz_xxxyy, g_y_y_xz_xxxyz, g_y_y_xz_xxxz, g_y_y_xz_xxxzz, g_y_y_xz_xxyy, g_y_y_xz_xxyyy, g_y_y_xz_xxyyz, g_y_y_xz_xxyz, g_y_y_xz_xxyzz, g_y_y_xz_xxzz, g_y_y_xz_xxzzz, g_y_y_xz_xyyy, g_y_y_xz_xyyyy, g_y_y_xz_xyyyz, g_y_y_xz_xyyz, g_y_y_xz_xyyzz, g_y_y_xz_xyzz, g_y_y_xz_xyzzz, g_y_y_xz_xzzz, g_y_y_xz_xzzzz, g_y_y_xz_yyyy, g_y_y_xz_yyyz, g_y_y_xz_yyzz, g_y_y_xz_yzzz, g_y_y_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxz_xxxx[k] = -g_y_y_xz_xxxx[k] * ab_x + g_y_y_xz_xxxxx[k];

                g_y_y_xxz_xxxy[k] = -g_y_y_xz_xxxy[k] * ab_x + g_y_y_xz_xxxxy[k];

                g_y_y_xxz_xxxz[k] = -g_y_y_xz_xxxz[k] * ab_x + g_y_y_xz_xxxxz[k];

                g_y_y_xxz_xxyy[k] = -g_y_y_xz_xxyy[k] * ab_x + g_y_y_xz_xxxyy[k];

                g_y_y_xxz_xxyz[k] = -g_y_y_xz_xxyz[k] * ab_x + g_y_y_xz_xxxyz[k];

                g_y_y_xxz_xxzz[k] = -g_y_y_xz_xxzz[k] * ab_x + g_y_y_xz_xxxzz[k];

                g_y_y_xxz_xyyy[k] = -g_y_y_xz_xyyy[k] * ab_x + g_y_y_xz_xxyyy[k];

                g_y_y_xxz_xyyz[k] = -g_y_y_xz_xyyz[k] * ab_x + g_y_y_xz_xxyyz[k];

                g_y_y_xxz_xyzz[k] = -g_y_y_xz_xyzz[k] * ab_x + g_y_y_xz_xxyzz[k];

                g_y_y_xxz_xzzz[k] = -g_y_y_xz_xzzz[k] * ab_x + g_y_y_xz_xxzzz[k];

                g_y_y_xxz_yyyy[k] = -g_y_y_xz_yyyy[k] * ab_x + g_y_y_xz_xyyyy[k];

                g_y_y_xxz_yyyz[k] = -g_y_y_xz_yyyz[k] * ab_x + g_y_y_xz_xyyyz[k];

                g_y_y_xxz_yyzz[k] = -g_y_y_xz_yyzz[k] * ab_x + g_y_y_xz_xyyzz[k];

                g_y_y_xxz_yzzz[k] = -g_y_y_xz_yzzz[k] * ab_x + g_y_y_xz_xyzzz[k];

                g_y_y_xxz_zzzz[k] = -g_y_y_xz_zzzz[k] * ab_x + g_y_y_xz_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyy_xxxx = cbuffer.data(fg_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_y_xyy_xxxy = cbuffer.data(fg_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_y_xyy_xxxz = cbuffer.data(fg_geom_11_off + 647 * ccomps * dcomps);

            auto g_y_y_xyy_xxyy = cbuffer.data(fg_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_y_xyy_xxyz = cbuffer.data(fg_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_y_xyy_xxzz = cbuffer.data(fg_geom_11_off + 650 * ccomps * dcomps);

            auto g_y_y_xyy_xyyy = cbuffer.data(fg_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_y_xyy_xyyz = cbuffer.data(fg_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_y_xyy_xyzz = cbuffer.data(fg_geom_11_off + 653 * ccomps * dcomps);

            auto g_y_y_xyy_xzzz = cbuffer.data(fg_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_y_xyy_yyyy = cbuffer.data(fg_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_y_xyy_yyyz = cbuffer.data(fg_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_y_xyy_yyzz = cbuffer.data(fg_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_y_xyy_yzzz = cbuffer.data(fg_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_y_xyy_zzzz = cbuffer.data(fg_geom_11_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyy_xxxx, g_y_y_xyy_xxxy, g_y_y_xyy_xxxz, g_y_y_xyy_xxyy, g_y_y_xyy_xxyz, g_y_y_xyy_xxzz, g_y_y_xyy_xyyy, g_y_y_xyy_xyyz, g_y_y_xyy_xyzz, g_y_y_xyy_xzzz, g_y_y_xyy_yyyy, g_y_y_xyy_yyyz, g_y_y_xyy_yyzz, g_y_y_xyy_yzzz, g_y_y_xyy_zzzz, g_y_y_yy_xxxx, g_y_y_yy_xxxxx, g_y_y_yy_xxxxy, g_y_y_yy_xxxxz, g_y_y_yy_xxxy, g_y_y_yy_xxxyy, g_y_y_yy_xxxyz, g_y_y_yy_xxxz, g_y_y_yy_xxxzz, g_y_y_yy_xxyy, g_y_y_yy_xxyyy, g_y_y_yy_xxyyz, g_y_y_yy_xxyz, g_y_y_yy_xxyzz, g_y_y_yy_xxzz, g_y_y_yy_xxzzz, g_y_y_yy_xyyy, g_y_y_yy_xyyyy, g_y_y_yy_xyyyz, g_y_y_yy_xyyz, g_y_y_yy_xyyzz, g_y_y_yy_xyzz, g_y_y_yy_xyzzz, g_y_y_yy_xzzz, g_y_y_yy_xzzzz, g_y_y_yy_yyyy, g_y_y_yy_yyyz, g_y_y_yy_yyzz, g_y_y_yy_yzzz, g_y_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyy_xxxx[k] = -g_y_y_yy_xxxx[k] * ab_x + g_y_y_yy_xxxxx[k];

                g_y_y_xyy_xxxy[k] = -g_y_y_yy_xxxy[k] * ab_x + g_y_y_yy_xxxxy[k];

                g_y_y_xyy_xxxz[k] = -g_y_y_yy_xxxz[k] * ab_x + g_y_y_yy_xxxxz[k];

                g_y_y_xyy_xxyy[k] = -g_y_y_yy_xxyy[k] * ab_x + g_y_y_yy_xxxyy[k];

                g_y_y_xyy_xxyz[k] = -g_y_y_yy_xxyz[k] * ab_x + g_y_y_yy_xxxyz[k];

                g_y_y_xyy_xxzz[k] = -g_y_y_yy_xxzz[k] * ab_x + g_y_y_yy_xxxzz[k];

                g_y_y_xyy_xyyy[k] = -g_y_y_yy_xyyy[k] * ab_x + g_y_y_yy_xxyyy[k];

                g_y_y_xyy_xyyz[k] = -g_y_y_yy_xyyz[k] * ab_x + g_y_y_yy_xxyyz[k];

                g_y_y_xyy_xyzz[k] = -g_y_y_yy_xyzz[k] * ab_x + g_y_y_yy_xxyzz[k];

                g_y_y_xyy_xzzz[k] = -g_y_y_yy_xzzz[k] * ab_x + g_y_y_yy_xxzzz[k];

                g_y_y_xyy_yyyy[k] = -g_y_y_yy_yyyy[k] * ab_x + g_y_y_yy_xyyyy[k];

                g_y_y_xyy_yyyz[k] = -g_y_y_yy_yyyz[k] * ab_x + g_y_y_yy_xyyyz[k];

                g_y_y_xyy_yyzz[k] = -g_y_y_yy_yyzz[k] * ab_x + g_y_y_yy_xyyzz[k];

                g_y_y_xyy_yzzz[k] = -g_y_y_yy_yzzz[k] * ab_x + g_y_y_yy_xyzzz[k];

                g_y_y_xyy_zzzz[k] = -g_y_y_yy_zzzz[k] * ab_x + g_y_y_yy_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyz_xxxx = cbuffer.data(fg_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_y_xyz_xxxy = cbuffer.data(fg_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_y_xyz_xxxz = cbuffer.data(fg_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_y_xyz_xxyy = cbuffer.data(fg_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_y_xyz_xxyz = cbuffer.data(fg_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_y_xyz_xxzz = cbuffer.data(fg_geom_11_off + 665 * ccomps * dcomps);

            auto g_y_y_xyz_xyyy = cbuffer.data(fg_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_y_xyz_xyyz = cbuffer.data(fg_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_y_xyz_xyzz = cbuffer.data(fg_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_y_xyz_xzzz = cbuffer.data(fg_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_y_xyz_yyyy = cbuffer.data(fg_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_y_xyz_yyyz = cbuffer.data(fg_geom_11_off + 671 * ccomps * dcomps);

            auto g_y_y_xyz_yyzz = cbuffer.data(fg_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_y_xyz_yzzz = cbuffer.data(fg_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_y_xyz_zzzz = cbuffer.data(fg_geom_11_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyz_xxxx, g_y_y_xyz_xxxy, g_y_y_xyz_xxxz, g_y_y_xyz_xxyy, g_y_y_xyz_xxyz, g_y_y_xyz_xxzz, g_y_y_xyz_xyyy, g_y_y_xyz_xyyz, g_y_y_xyz_xyzz, g_y_y_xyz_xzzz, g_y_y_xyz_yyyy, g_y_y_xyz_yyyz, g_y_y_xyz_yyzz, g_y_y_xyz_yzzz, g_y_y_xyz_zzzz, g_y_y_yz_xxxx, g_y_y_yz_xxxxx, g_y_y_yz_xxxxy, g_y_y_yz_xxxxz, g_y_y_yz_xxxy, g_y_y_yz_xxxyy, g_y_y_yz_xxxyz, g_y_y_yz_xxxz, g_y_y_yz_xxxzz, g_y_y_yz_xxyy, g_y_y_yz_xxyyy, g_y_y_yz_xxyyz, g_y_y_yz_xxyz, g_y_y_yz_xxyzz, g_y_y_yz_xxzz, g_y_y_yz_xxzzz, g_y_y_yz_xyyy, g_y_y_yz_xyyyy, g_y_y_yz_xyyyz, g_y_y_yz_xyyz, g_y_y_yz_xyyzz, g_y_y_yz_xyzz, g_y_y_yz_xyzzz, g_y_y_yz_xzzz, g_y_y_yz_xzzzz, g_y_y_yz_yyyy, g_y_y_yz_yyyz, g_y_y_yz_yyzz, g_y_y_yz_yzzz, g_y_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyz_xxxx[k] = -g_y_y_yz_xxxx[k] * ab_x + g_y_y_yz_xxxxx[k];

                g_y_y_xyz_xxxy[k] = -g_y_y_yz_xxxy[k] * ab_x + g_y_y_yz_xxxxy[k];

                g_y_y_xyz_xxxz[k] = -g_y_y_yz_xxxz[k] * ab_x + g_y_y_yz_xxxxz[k];

                g_y_y_xyz_xxyy[k] = -g_y_y_yz_xxyy[k] * ab_x + g_y_y_yz_xxxyy[k];

                g_y_y_xyz_xxyz[k] = -g_y_y_yz_xxyz[k] * ab_x + g_y_y_yz_xxxyz[k];

                g_y_y_xyz_xxzz[k] = -g_y_y_yz_xxzz[k] * ab_x + g_y_y_yz_xxxzz[k];

                g_y_y_xyz_xyyy[k] = -g_y_y_yz_xyyy[k] * ab_x + g_y_y_yz_xxyyy[k];

                g_y_y_xyz_xyyz[k] = -g_y_y_yz_xyyz[k] * ab_x + g_y_y_yz_xxyyz[k];

                g_y_y_xyz_xyzz[k] = -g_y_y_yz_xyzz[k] * ab_x + g_y_y_yz_xxyzz[k];

                g_y_y_xyz_xzzz[k] = -g_y_y_yz_xzzz[k] * ab_x + g_y_y_yz_xxzzz[k];

                g_y_y_xyz_yyyy[k] = -g_y_y_yz_yyyy[k] * ab_x + g_y_y_yz_xyyyy[k];

                g_y_y_xyz_yyyz[k] = -g_y_y_yz_yyyz[k] * ab_x + g_y_y_yz_xyyyz[k];

                g_y_y_xyz_yyzz[k] = -g_y_y_yz_yyzz[k] * ab_x + g_y_y_yz_xyyzz[k];

                g_y_y_xyz_yzzz[k] = -g_y_y_yz_yzzz[k] * ab_x + g_y_y_yz_xyzzz[k];

                g_y_y_xyz_zzzz[k] = -g_y_y_yz_zzzz[k] * ab_x + g_y_y_yz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzz_xxxx = cbuffer.data(fg_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_y_xzz_xxxy = cbuffer.data(fg_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_y_xzz_xxxz = cbuffer.data(fg_geom_11_off + 677 * ccomps * dcomps);

            auto g_y_y_xzz_xxyy = cbuffer.data(fg_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_y_xzz_xxyz = cbuffer.data(fg_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_y_xzz_xxzz = cbuffer.data(fg_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_y_xzz_xyyy = cbuffer.data(fg_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_y_xzz_xyyz = cbuffer.data(fg_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_y_xzz_xyzz = cbuffer.data(fg_geom_11_off + 683 * ccomps * dcomps);

            auto g_y_y_xzz_xzzz = cbuffer.data(fg_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_y_xzz_yyyy = cbuffer.data(fg_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_y_xzz_yyyz = cbuffer.data(fg_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_y_xzz_yyzz = cbuffer.data(fg_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_y_xzz_yzzz = cbuffer.data(fg_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_y_xzz_zzzz = cbuffer.data(fg_geom_11_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzz_xxxx, g_y_y_xzz_xxxy, g_y_y_xzz_xxxz, g_y_y_xzz_xxyy, g_y_y_xzz_xxyz, g_y_y_xzz_xxzz, g_y_y_xzz_xyyy, g_y_y_xzz_xyyz, g_y_y_xzz_xyzz, g_y_y_xzz_xzzz, g_y_y_xzz_yyyy, g_y_y_xzz_yyyz, g_y_y_xzz_yyzz, g_y_y_xzz_yzzz, g_y_y_xzz_zzzz, g_y_y_zz_xxxx, g_y_y_zz_xxxxx, g_y_y_zz_xxxxy, g_y_y_zz_xxxxz, g_y_y_zz_xxxy, g_y_y_zz_xxxyy, g_y_y_zz_xxxyz, g_y_y_zz_xxxz, g_y_y_zz_xxxzz, g_y_y_zz_xxyy, g_y_y_zz_xxyyy, g_y_y_zz_xxyyz, g_y_y_zz_xxyz, g_y_y_zz_xxyzz, g_y_y_zz_xxzz, g_y_y_zz_xxzzz, g_y_y_zz_xyyy, g_y_y_zz_xyyyy, g_y_y_zz_xyyyz, g_y_y_zz_xyyz, g_y_y_zz_xyyzz, g_y_y_zz_xyzz, g_y_y_zz_xyzzz, g_y_y_zz_xzzz, g_y_y_zz_xzzzz, g_y_y_zz_yyyy, g_y_y_zz_yyyz, g_y_y_zz_yyzz, g_y_y_zz_yzzz, g_y_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzz_xxxx[k] = -g_y_y_zz_xxxx[k] * ab_x + g_y_y_zz_xxxxx[k];

                g_y_y_xzz_xxxy[k] = -g_y_y_zz_xxxy[k] * ab_x + g_y_y_zz_xxxxy[k];

                g_y_y_xzz_xxxz[k] = -g_y_y_zz_xxxz[k] * ab_x + g_y_y_zz_xxxxz[k];

                g_y_y_xzz_xxyy[k] = -g_y_y_zz_xxyy[k] * ab_x + g_y_y_zz_xxxyy[k];

                g_y_y_xzz_xxyz[k] = -g_y_y_zz_xxyz[k] * ab_x + g_y_y_zz_xxxyz[k];

                g_y_y_xzz_xxzz[k] = -g_y_y_zz_xxzz[k] * ab_x + g_y_y_zz_xxxzz[k];

                g_y_y_xzz_xyyy[k] = -g_y_y_zz_xyyy[k] * ab_x + g_y_y_zz_xxyyy[k];

                g_y_y_xzz_xyyz[k] = -g_y_y_zz_xyyz[k] * ab_x + g_y_y_zz_xxyyz[k];

                g_y_y_xzz_xyzz[k] = -g_y_y_zz_xyzz[k] * ab_x + g_y_y_zz_xxyzz[k];

                g_y_y_xzz_xzzz[k] = -g_y_y_zz_xzzz[k] * ab_x + g_y_y_zz_xxzzz[k];

                g_y_y_xzz_yyyy[k] = -g_y_y_zz_yyyy[k] * ab_x + g_y_y_zz_xyyyy[k];

                g_y_y_xzz_yyyz[k] = -g_y_y_zz_yyyz[k] * ab_x + g_y_y_zz_xyyyz[k];

                g_y_y_xzz_yyzz[k] = -g_y_y_zz_yyzz[k] * ab_x + g_y_y_zz_xyyzz[k];

                g_y_y_xzz_yzzz[k] = -g_y_y_zz_yzzz[k] * ab_x + g_y_y_zz_xyzzz[k];

                g_y_y_xzz_zzzz[k] = -g_y_y_zz_zzzz[k] * ab_x + g_y_y_zz_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyy_xxxx = cbuffer.data(fg_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_y_yyy_xxxy = cbuffer.data(fg_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_y_yyy_xxxz = cbuffer.data(fg_geom_11_off + 692 * ccomps * dcomps);

            auto g_y_y_yyy_xxyy = cbuffer.data(fg_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_y_yyy_xxyz = cbuffer.data(fg_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_y_yyy_xxzz = cbuffer.data(fg_geom_11_off + 695 * ccomps * dcomps);

            auto g_y_y_yyy_xyyy = cbuffer.data(fg_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_y_yyy_xyyz = cbuffer.data(fg_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_y_yyy_xyzz = cbuffer.data(fg_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_y_yyy_xzzz = cbuffer.data(fg_geom_11_off + 699 * ccomps * dcomps);

            auto g_y_y_yyy_yyyy = cbuffer.data(fg_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_y_yyy_yyyz = cbuffer.data(fg_geom_11_off + 701 * ccomps * dcomps);

            auto g_y_y_yyy_yyzz = cbuffer.data(fg_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_y_yyy_yzzz = cbuffer.data(fg_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_y_yyy_zzzz = cbuffer.data(fg_geom_11_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxxx, g_0_y_yy_xxxy, g_0_y_yy_xxxz, g_0_y_yy_xxyy, g_0_y_yy_xxyz, g_0_y_yy_xxzz, g_0_y_yy_xyyy, g_0_y_yy_xyyz, g_0_y_yy_xyzz, g_0_y_yy_xzzz, g_0_y_yy_yyyy, g_0_y_yy_yyyz, g_0_y_yy_yyzz, g_0_y_yy_yzzz, g_0_y_yy_zzzz, g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz, g_y_y_yy_xxxx, g_y_y_yy_xxxxy, g_y_y_yy_xxxy, g_y_y_yy_xxxyy, g_y_y_yy_xxxyz, g_y_y_yy_xxxz, g_y_y_yy_xxyy, g_y_y_yy_xxyyy, g_y_y_yy_xxyyz, g_y_y_yy_xxyz, g_y_y_yy_xxyzz, g_y_y_yy_xxzz, g_y_y_yy_xyyy, g_y_y_yy_xyyyy, g_y_y_yy_xyyyz, g_y_y_yy_xyyz, g_y_y_yy_xyyzz, g_y_y_yy_xyzz, g_y_y_yy_xyzzz, g_y_y_yy_xzzz, g_y_y_yy_yyyy, g_y_y_yy_yyyyy, g_y_y_yy_yyyyz, g_y_y_yy_yyyz, g_y_y_yy_yyyzz, g_y_y_yy_yyzz, g_y_y_yy_yyzzz, g_y_y_yy_yzzz, g_y_y_yy_yzzzz, g_y_y_yy_zzzz, g_y_y_yyy_xxxx, g_y_y_yyy_xxxy, g_y_y_yyy_xxxz, g_y_y_yyy_xxyy, g_y_y_yyy_xxyz, g_y_y_yyy_xxzz, g_y_y_yyy_xyyy, g_y_y_yyy_xyyz, g_y_y_yyy_xyzz, g_y_y_yyy_xzzz, g_y_y_yyy_yyyy, g_y_y_yyy_yyyz, g_y_y_yyy_yyzz, g_y_y_yyy_yzzz, g_y_y_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyy_xxxx[k] = -g_0_y_yy_xxxx[k] + g_y_0_yy_xxxx[k] - g_y_y_yy_xxxx[k] * ab_y + g_y_y_yy_xxxxy[k];

                g_y_y_yyy_xxxy[k] = -g_0_y_yy_xxxy[k] + g_y_0_yy_xxxy[k] - g_y_y_yy_xxxy[k] * ab_y + g_y_y_yy_xxxyy[k];

                g_y_y_yyy_xxxz[k] = -g_0_y_yy_xxxz[k] + g_y_0_yy_xxxz[k] - g_y_y_yy_xxxz[k] * ab_y + g_y_y_yy_xxxyz[k];

                g_y_y_yyy_xxyy[k] = -g_0_y_yy_xxyy[k] + g_y_0_yy_xxyy[k] - g_y_y_yy_xxyy[k] * ab_y + g_y_y_yy_xxyyy[k];

                g_y_y_yyy_xxyz[k] = -g_0_y_yy_xxyz[k] + g_y_0_yy_xxyz[k] - g_y_y_yy_xxyz[k] * ab_y + g_y_y_yy_xxyyz[k];

                g_y_y_yyy_xxzz[k] = -g_0_y_yy_xxzz[k] + g_y_0_yy_xxzz[k] - g_y_y_yy_xxzz[k] * ab_y + g_y_y_yy_xxyzz[k];

                g_y_y_yyy_xyyy[k] = -g_0_y_yy_xyyy[k] + g_y_0_yy_xyyy[k] - g_y_y_yy_xyyy[k] * ab_y + g_y_y_yy_xyyyy[k];

                g_y_y_yyy_xyyz[k] = -g_0_y_yy_xyyz[k] + g_y_0_yy_xyyz[k] - g_y_y_yy_xyyz[k] * ab_y + g_y_y_yy_xyyyz[k];

                g_y_y_yyy_xyzz[k] = -g_0_y_yy_xyzz[k] + g_y_0_yy_xyzz[k] - g_y_y_yy_xyzz[k] * ab_y + g_y_y_yy_xyyzz[k];

                g_y_y_yyy_xzzz[k] = -g_0_y_yy_xzzz[k] + g_y_0_yy_xzzz[k] - g_y_y_yy_xzzz[k] * ab_y + g_y_y_yy_xyzzz[k];

                g_y_y_yyy_yyyy[k] = -g_0_y_yy_yyyy[k] + g_y_0_yy_yyyy[k] - g_y_y_yy_yyyy[k] * ab_y + g_y_y_yy_yyyyy[k];

                g_y_y_yyy_yyyz[k] = -g_0_y_yy_yyyz[k] + g_y_0_yy_yyyz[k] - g_y_y_yy_yyyz[k] * ab_y + g_y_y_yy_yyyyz[k];

                g_y_y_yyy_yyzz[k] = -g_0_y_yy_yyzz[k] + g_y_0_yy_yyzz[k] - g_y_y_yy_yyzz[k] * ab_y + g_y_y_yy_yyyzz[k];

                g_y_y_yyy_yzzz[k] = -g_0_y_yy_yzzz[k] + g_y_0_yy_yzzz[k] - g_y_y_yy_yzzz[k] * ab_y + g_y_y_yy_yyzzz[k];

                g_y_y_yyy_zzzz[k] = -g_0_y_yy_zzzz[k] + g_y_0_yy_zzzz[k] - g_y_y_yy_zzzz[k] * ab_y + g_y_y_yy_yzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyz_xxxx = cbuffer.data(fg_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_y_yyz_xxxy = cbuffer.data(fg_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_y_yyz_xxxz = cbuffer.data(fg_geom_11_off + 707 * ccomps * dcomps);

            auto g_y_y_yyz_xxyy = cbuffer.data(fg_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_y_yyz_xxyz = cbuffer.data(fg_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_y_yyz_xxzz = cbuffer.data(fg_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_y_yyz_xyyy = cbuffer.data(fg_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_y_yyz_xyyz = cbuffer.data(fg_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_y_yyz_xyzz = cbuffer.data(fg_geom_11_off + 713 * ccomps * dcomps);

            auto g_y_y_yyz_xzzz = cbuffer.data(fg_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_y_yyz_yyyy = cbuffer.data(fg_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_y_yyz_yyyz = cbuffer.data(fg_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_y_yyz_yyzz = cbuffer.data(fg_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_y_yyz_yzzz = cbuffer.data(fg_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_y_yyz_zzzz = cbuffer.data(fg_geom_11_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yy_xxxx, g_y_y_yy_xxxxz, g_y_y_yy_xxxy, g_y_y_yy_xxxyz, g_y_y_yy_xxxz, g_y_y_yy_xxxzz, g_y_y_yy_xxyy, g_y_y_yy_xxyyz, g_y_y_yy_xxyz, g_y_y_yy_xxyzz, g_y_y_yy_xxzz, g_y_y_yy_xxzzz, g_y_y_yy_xyyy, g_y_y_yy_xyyyz, g_y_y_yy_xyyz, g_y_y_yy_xyyzz, g_y_y_yy_xyzz, g_y_y_yy_xyzzz, g_y_y_yy_xzzz, g_y_y_yy_xzzzz, g_y_y_yy_yyyy, g_y_y_yy_yyyyz, g_y_y_yy_yyyz, g_y_y_yy_yyyzz, g_y_y_yy_yyzz, g_y_y_yy_yyzzz, g_y_y_yy_yzzz, g_y_y_yy_yzzzz, g_y_y_yy_zzzz, g_y_y_yy_zzzzz, g_y_y_yyz_xxxx, g_y_y_yyz_xxxy, g_y_y_yyz_xxxz, g_y_y_yyz_xxyy, g_y_y_yyz_xxyz, g_y_y_yyz_xxzz, g_y_y_yyz_xyyy, g_y_y_yyz_xyyz, g_y_y_yyz_xyzz, g_y_y_yyz_xzzz, g_y_y_yyz_yyyy, g_y_y_yyz_yyyz, g_y_y_yyz_yyzz, g_y_y_yyz_yzzz, g_y_y_yyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyz_xxxx[k] = -g_y_y_yy_xxxx[k] * ab_z + g_y_y_yy_xxxxz[k];

                g_y_y_yyz_xxxy[k] = -g_y_y_yy_xxxy[k] * ab_z + g_y_y_yy_xxxyz[k];

                g_y_y_yyz_xxxz[k] = -g_y_y_yy_xxxz[k] * ab_z + g_y_y_yy_xxxzz[k];

                g_y_y_yyz_xxyy[k] = -g_y_y_yy_xxyy[k] * ab_z + g_y_y_yy_xxyyz[k];

                g_y_y_yyz_xxyz[k] = -g_y_y_yy_xxyz[k] * ab_z + g_y_y_yy_xxyzz[k];

                g_y_y_yyz_xxzz[k] = -g_y_y_yy_xxzz[k] * ab_z + g_y_y_yy_xxzzz[k];

                g_y_y_yyz_xyyy[k] = -g_y_y_yy_xyyy[k] * ab_z + g_y_y_yy_xyyyz[k];

                g_y_y_yyz_xyyz[k] = -g_y_y_yy_xyyz[k] * ab_z + g_y_y_yy_xyyzz[k];

                g_y_y_yyz_xyzz[k] = -g_y_y_yy_xyzz[k] * ab_z + g_y_y_yy_xyzzz[k];

                g_y_y_yyz_xzzz[k] = -g_y_y_yy_xzzz[k] * ab_z + g_y_y_yy_xzzzz[k];

                g_y_y_yyz_yyyy[k] = -g_y_y_yy_yyyy[k] * ab_z + g_y_y_yy_yyyyz[k];

                g_y_y_yyz_yyyz[k] = -g_y_y_yy_yyyz[k] * ab_z + g_y_y_yy_yyyzz[k];

                g_y_y_yyz_yyzz[k] = -g_y_y_yy_yyzz[k] * ab_z + g_y_y_yy_yyzzz[k];

                g_y_y_yyz_yzzz[k] = -g_y_y_yy_yzzz[k] * ab_z + g_y_y_yy_yzzzz[k];

                g_y_y_yyz_zzzz[k] = -g_y_y_yy_zzzz[k] * ab_z + g_y_y_yy_zzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzz_xxxx = cbuffer.data(fg_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_y_yzz_xxxy = cbuffer.data(fg_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_y_yzz_xxxz = cbuffer.data(fg_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_y_yzz_xxyy = cbuffer.data(fg_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_y_yzz_xxyz = cbuffer.data(fg_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_y_yzz_xxzz = cbuffer.data(fg_geom_11_off + 725 * ccomps * dcomps);

            auto g_y_y_yzz_xyyy = cbuffer.data(fg_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_y_yzz_xyyz = cbuffer.data(fg_geom_11_off + 727 * ccomps * dcomps);

            auto g_y_y_yzz_xyzz = cbuffer.data(fg_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_y_yzz_xzzz = cbuffer.data(fg_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_y_yzz_yyyy = cbuffer.data(fg_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_y_yzz_yyyz = cbuffer.data(fg_geom_11_off + 731 * ccomps * dcomps);

            auto g_y_y_yzz_yyzz = cbuffer.data(fg_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_y_yzz_yzzz = cbuffer.data(fg_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_y_yzz_zzzz = cbuffer.data(fg_geom_11_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yz_xxxx, g_y_y_yz_xxxxz, g_y_y_yz_xxxy, g_y_y_yz_xxxyz, g_y_y_yz_xxxz, g_y_y_yz_xxxzz, g_y_y_yz_xxyy, g_y_y_yz_xxyyz, g_y_y_yz_xxyz, g_y_y_yz_xxyzz, g_y_y_yz_xxzz, g_y_y_yz_xxzzz, g_y_y_yz_xyyy, g_y_y_yz_xyyyz, g_y_y_yz_xyyz, g_y_y_yz_xyyzz, g_y_y_yz_xyzz, g_y_y_yz_xyzzz, g_y_y_yz_xzzz, g_y_y_yz_xzzzz, g_y_y_yz_yyyy, g_y_y_yz_yyyyz, g_y_y_yz_yyyz, g_y_y_yz_yyyzz, g_y_y_yz_yyzz, g_y_y_yz_yyzzz, g_y_y_yz_yzzz, g_y_y_yz_yzzzz, g_y_y_yz_zzzz, g_y_y_yz_zzzzz, g_y_y_yzz_xxxx, g_y_y_yzz_xxxy, g_y_y_yzz_xxxz, g_y_y_yzz_xxyy, g_y_y_yzz_xxyz, g_y_y_yzz_xxzz, g_y_y_yzz_xyyy, g_y_y_yzz_xyyz, g_y_y_yzz_xyzz, g_y_y_yzz_xzzz, g_y_y_yzz_yyyy, g_y_y_yzz_yyyz, g_y_y_yzz_yyzz, g_y_y_yzz_yzzz, g_y_y_yzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzz_xxxx[k] = -g_y_y_yz_xxxx[k] * ab_z + g_y_y_yz_xxxxz[k];

                g_y_y_yzz_xxxy[k] = -g_y_y_yz_xxxy[k] * ab_z + g_y_y_yz_xxxyz[k];

                g_y_y_yzz_xxxz[k] = -g_y_y_yz_xxxz[k] * ab_z + g_y_y_yz_xxxzz[k];

                g_y_y_yzz_xxyy[k] = -g_y_y_yz_xxyy[k] * ab_z + g_y_y_yz_xxyyz[k];

                g_y_y_yzz_xxyz[k] = -g_y_y_yz_xxyz[k] * ab_z + g_y_y_yz_xxyzz[k];

                g_y_y_yzz_xxzz[k] = -g_y_y_yz_xxzz[k] * ab_z + g_y_y_yz_xxzzz[k];

                g_y_y_yzz_xyyy[k] = -g_y_y_yz_xyyy[k] * ab_z + g_y_y_yz_xyyyz[k];

                g_y_y_yzz_xyyz[k] = -g_y_y_yz_xyyz[k] * ab_z + g_y_y_yz_xyyzz[k];

                g_y_y_yzz_xyzz[k] = -g_y_y_yz_xyzz[k] * ab_z + g_y_y_yz_xyzzz[k];

                g_y_y_yzz_xzzz[k] = -g_y_y_yz_xzzz[k] * ab_z + g_y_y_yz_xzzzz[k];

                g_y_y_yzz_yyyy[k] = -g_y_y_yz_yyyy[k] * ab_z + g_y_y_yz_yyyyz[k];

                g_y_y_yzz_yyyz[k] = -g_y_y_yz_yyyz[k] * ab_z + g_y_y_yz_yyyzz[k];

                g_y_y_yzz_yyzz[k] = -g_y_y_yz_yyzz[k] * ab_z + g_y_y_yz_yyzzz[k];

                g_y_y_yzz_yzzz[k] = -g_y_y_yz_yzzz[k] * ab_z + g_y_y_yz_yzzzz[k];

                g_y_y_yzz_zzzz[k] = -g_y_y_yz_zzzz[k] * ab_z + g_y_y_yz_zzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzz_xxxx = cbuffer.data(fg_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_y_zzz_xxxy = cbuffer.data(fg_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_y_zzz_xxxz = cbuffer.data(fg_geom_11_off + 737 * ccomps * dcomps);

            auto g_y_y_zzz_xxyy = cbuffer.data(fg_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_y_zzz_xxyz = cbuffer.data(fg_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_y_zzz_xxzz = cbuffer.data(fg_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_y_zzz_xyyy = cbuffer.data(fg_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_y_zzz_xyyz = cbuffer.data(fg_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_y_zzz_xyzz = cbuffer.data(fg_geom_11_off + 743 * ccomps * dcomps);

            auto g_y_y_zzz_xzzz = cbuffer.data(fg_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_y_zzz_yyyy = cbuffer.data(fg_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_y_zzz_yyyz = cbuffer.data(fg_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_y_zzz_yyzz = cbuffer.data(fg_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_y_zzz_yzzz = cbuffer.data(fg_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_y_zzz_zzzz = cbuffer.data(fg_geom_11_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zz_xxxx, g_y_y_zz_xxxxz, g_y_y_zz_xxxy, g_y_y_zz_xxxyz, g_y_y_zz_xxxz, g_y_y_zz_xxxzz, g_y_y_zz_xxyy, g_y_y_zz_xxyyz, g_y_y_zz_xxyz, g_y_y_zz_xxyzz, g_y_y_zz_xxzz, g_y_y_zz_xxzzz, g_y_y_zz_xyyy, g_y_y_zz_xyyyz, g_y_y_zz_xyyz, g_y_y_zz_xyyzz, g_y_y_zz_xyzz, g_y_y_zz_xyzzz, g_y_y_zz_xzzz, g_y_y_zz_xzzzz, g_y_y_zz_yyyy, g_y_y_zz_yyyyz, g_y_y_zz_yyyz, g_y_y_zz_yyyzz, g_y_y_zz_yyzz, g_y_y_zz_yyzzz, g_y_y_zz_yzzz, g_y_y_zz_yzzzz, g_y_y_zz_zzzz, g_y_y_zz_zzzzz, g_y_y_zzz_xxxx, g_y_y_zzz_xxxy, g_y_y_zzz_xxxz, g_y_y_zzz_xxyy, g_y_y_zzz_xxyz, g_y_y_zzz_xxzz, g_y_y_zzz_xyyy, g_y_y_zzz_xyyz, g_y_y_zzz_xyzz, g_y_y_zzz_xzzz, g_y_y_zzz_yyyy, g_y_y_zzz_yyyz, g_y_y_zzz_yyzz, g_y_y_zzz_yzzz, g_y_y_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzz_xxxx[k] = -g_y_y_zz_xxxx[k] * ab_z + g_y_y_zz_xxxxz[k];

                g_y_y_zzz_xxxy[k] = -g_y_y_zz_xxxy[k] * ab_z + g_y_y_zz_xxxyz[k];

                g_y_y_zzz_xxxz[k] = -g_y_y_zz_xxxz[k] * ab_z + g_y_y_zz_xxxzz[k];

                g_y_y_zzz_xxyy[k] = -g_y_y_zz_xxyy[k] * ab_z + g_y_y_zz_xxyyz[k];

                g_y_y_zzz_xxyz[k] = -g_y_y_zz_xxyz[k] * ab_z + g_y_y_zz_xxyzz[k];

                g_y_y_zzz_xxzz[k] = -g_y_y_zz_xxzz[k] * ab_z + g_y_y_zz_xxzzz[k];

                g_y_y_zzz_xyyy[k] = -g_y_y_zz_xyyy[k] * ab_z + g_y_y_zz_xyyyz[k];

                g_y_y_zzz_xyyz[k] = -g_y_y_zz_xyyz[k] * ab_z + g_y_y_zz_xyyzz[k];

                g_y_y_zzz_xyzz[k] = -g_y_y_zz_xyzz[k] * ab_z + g_y_y_zz_xyzzz[k];

                g_y_y_zzz_xzzz[k] = -g_y_y_zz_xzzz[k] * ab_z + g_y_y_zz_xzzzz[k];

                g_y_y_zzz_yyyy[k] = -g_y_y_zz_yyyy[k] * ab_z + g_y_y_zz_yyyyz[k];

                g_y_y_zzz_yyyz[k] = -g_y_y_zz_yyyz[k] * ab_z + g_y_y_zz_yyyzz[k];

                g_y_y_zzz_yyzz[k] = -g_y_y_zz_yyzz[k] * ab_z + g_y_y_zz_yyzzz[k];

                g_y_y_zzz_yzzz[k] = -g_y_y_zz_yzzz[k] * ab_z + g_y_y_zz_yzzzz[k];

                g_y_y_zzz_zzzz[k] = -g_y_y_zz_zzzz[k] * ab_z + g_y_y_zz_zzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxx_xxxx = cbuffer.data(fg_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_z_xxx_xxxy = cbuffer.data(fg_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_z_xxx_xxxz = cbuffer.data(fg_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_z_xxx_xxyy = cbuffer.data(fg_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_z_xxx_xxyz = cbuffer.data(fg_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_z_xxx_xxzz = cbuffer.data(fg_geom_11_off + 755 * ccomps * dcomps);

            auto g_y_z_xxx_xyyy = cbuffer.data(fg_geom_11_off + 756 * ccomps * dcomps);

            auto g_y_z_xxx_xyyz = cbuffer.data(fg_geom_11_off + 757 * ccomps * dcomps);

            auto g_y_z_xxx_xyzz = cbuffer.data(fg_geom_11_off + 758 * ccomps * dcomps);

            auto g_y_z_xxx_xzzz = cbuffer.data(fg_geom_11_off + 759 * ccomps * dcomps);

            auto g_y_z_xxx_yyyy = cbuffer.data(fg_geom_11_off + 760 * ccomps * dcomps);

            auto g_y_z_xxx_yyyz = cbuffer.data(fg_geom_11_off + 761 * ccomps * dcomps);

            auto g_y_z_xxx_yyzz = cbuffer.data(fg_geom_11_off + 762 * ccomps * dcomps);

            auto g_y_z_xxx_yzzz = cbuffer.data(fg_geom_11_off + 763 * ccomps * dcomps);

            auto g_y_z_xxx_zzzz = cbuffer.data(fg_geom_11_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xx_xxxx, g_y_z_xx_xxxxx, g_y_z_xx_xxxxy, g_y_z_xx_xxxxz, g_y_z_xx_xxxy, g_y_z_xx_xxxyy, g_y_z_xx_xxxyz, g_y_z_xx_xxxz, g_y_z_xx_xxxzz, g_y_z_xx_xxyy, g_y_z_xx_xxyyy, g_y_z_xx_xxyyz, g_y_z_xx_xxyz, g_y_z_xx_xxyzz, g_y_z_xx_xxzz, g_y_z_xx_xxzzz, g_y_z_xx_xyyy, g_y_z_xx_xyyyy, g_y_z_xx_xyyyz, g_y_z_xx_xyyz, g_y_z_xx_xyyzz, g_y_z_xx_xyzz, g_y_z_xx_xyzzz, g_y_z_xx_xzzz, g_y_z_xx_xzzzz, g_y_z_xx_yyyy, g_y_z_xx_yyyz, g_y_z_xx_yyzz, g_y_z_xx_yzzz, g_y_z_xx_zzzz, g_y_z_xxx_xxxx, g_y_z_xxx_xxxy, g_y_z_xxx_xxxz, g_y_z_xxx_xxyy, g_y_z_xxx_xxyz, g_y_z_xxx_xxzz, g_y_z_xxx_xyyy, g_y_z_xxx_xyyz, g_y_z_xxx_xyzz, g_y_z_xxx_xzzz, g_y_z_xxx_yyyy, g_y_z_xxx_yyyz, g_y_z_xxx_yyzz, g_y_z_xxx_yzzz, g_y_z_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxx_xxxx[k] = -g_y_z_xx_xxxx[k] * ab_x + g_y_z_xx_xxxxx[k];

                g_y_z_xxx_xxxy[k] = -g_y_z_xx_xxxy[k] * ab_x + g_y_z_xx_xxxxy[k];

                g_y_z_xxx_xxxz[k] = -g_y_z_xx_xxxz[k] * ab_x + g_y_z_xx_xxxxz[k];

                g_y_z_xxx_xxyy[k] = -g_y_z_xx_xxyy[k] * ab_x + g_y_z_xx_xxxyy[k];

                g_y_z_xxx_xxyz[k] = -g_y_z_xx_xxyz[k] * ab_x + g_y_z_xx_xxxyz[k];

                g_y_z_xxx_xxzz[k] = -g_y_z_xx_xxzz[k] * ab_x + g_y_z_xx_xxxzz[k];

                g_y_z_xxx_xyyy[k] = -g_y_z_xx_xyyy[k] * ab_x + g_y_z_xx_xxyyy[k];

                g_y_z_xxx_xyyz[k] = -g_y_z_xx_xyyz[k] * ab_x + g_y_z_xx_xxyyz[k];

                g_y_z_xxx_xyzz[k] = -g_y_z_xx_xyzz[k] * ab_x + g_y_z_xx_xxyzz[k];

                g_y_z_xxx_xzzz[k] = -g_y_z_xx_xzzz[k] * ab_x + g_y_z_xx_xxzzz[k];

                g_y_z_xxx_yyyy[k] = -g_y_z_xx_yyyy[k] * ab_x + g_y_z_xx_xyyyy[k];

                g_y_z_xxx_yyyz[k] = -g_y_z_xx_yyyz[k] * ab_x + g_y_z_xx_xyyyz[k];

                g_y_z_xxx_yyzz[k] = -g_y_z_xx_yyzz[k] * ab_x + g_y_z_xx_xyyzz[k];

                g_y_z_xxx_yzzz[k] = -g_y_z_xx_yzzz[k] * ab_x + g_y_z_xx_xyzzz[k];

                g_y_z_xxx_zzzz[k] = -g_y_z_xx_zzzz[k] * ab_x + g_y_z_xx_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxy_xxxx = cbuffer.data(fg_geom_11_off + 765 * ccomps * dcomps);

            auto g_y_z_xxy_xxxy = cbuffer.data(fg_geom_11_off + 766 * ccomps * dcomps);

            auto g_y_z_xxy_xxxz = cbuffer.data(fg_geom_11_off + 767 * ccomps * dcomps);

            auto g_y_z_xxy_xxyy = cbuffer.data(fg_geom_11_off + 768 * ccomps * dcomps);

            auto g_y_z_xxy_xxyz = cbuffer.data(fg_geom_11_off + 769 * ccomps * dcomps);

            auto g_y_z_xxy_xxzz = cbuffer.data(fg_geom_11_off + 770 * ccomps * dcomps);

            auto g_y_z_xxy_xyyy = cbuffer.data(fg_geom_11_off + 771 * ccomps * dcomps);

            auto g_y_z_xxy_xyyz = cbuffer.data(fg_geom_11_off + 772 * ccomps * dcomps);

            auto g_y_z_xxy_xyzz = cbuffer.data(fg_geom_11_off + 773 * ccomps * dcomps);

            auto g_y_z_xxy_xzzz = cbuffer.data(fg_geom_11_off + 774 * ccomps * dcomps);

            auto g_y_z_xxy_yyyy = cbuffer.data(fg_geom_11_off + 775 * ccomps * dcomps);

            auto g_y_z_xxy_yyyz = cbuffer.data(fg_geom_11_off + 776 * ccomps * dcomps);

            auto g_y_z_xxy_yyzz = cbuffer.data(fg_geom_11_off + 777 * ccomps * dcomps);

            auto g_y_z_xxy_yzzz = cbuffer.data(fg_geom_11_off + 778 * ccomps * dcomps);

            auto g_y_z_xxy_zzzz = cbuffer.data(fg_geom_11_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxy_xxxx, g_y_z_xxy_xxxy, g_y_z_xxy_xxxz, g_y_z_xxy_xxyy, g_y_z_xxy_xxyz, g_y_z_xxy_xxzz, g_y_z_xxy_xyyy, g_y_z_xxy_xyyz, g_y_z_xxy_xyzz, g_y_z_xxy_xzzz, g_y_z_xxy_yyyy, g_y_z_xxy_yyyz, g_y_z_xxy_yyzz, g_y_z_xxy_yzzz, g_y_z_xxy_zzzz, g_y_z_xy_xxxx, g_y_z_xy_xxxxx, g_y_z_xy_xxxxy, g_y_z_xy_xxxxz, g_y_z_xy_xxxy, g_y_z_xy_xxxyy, g_y_z_xy_xxxyz, g_y_z_xy_xxxz, g_y_z_xy_xxxzz, g_y_z_xy_xxyy, g_y_z_xy_xxyyy, g_y_z_xy_xxyyz, g_y_z_xy_xxyz, g_y_z_xy_xxyzz, g_y_z_xy_xxzz, g_y_z_xy_xxzzz, g_y_z_xy_xyyy, g_y_z_xy_xyyyy, g_y_z_xy_xyyyz, g_y_z_xy_xyyz, g_y_z_xy_xyyzz, g_y_z_xy_xyzz, g_y_z_xy_xyzzz, g_y_z_xy_xzzz, g_y_z_xy_xzzzz, g_y_z_xy_yyyy, g_y_z_xy_yyyz, g_y_z_xy_yyzz, g_y_z_xy_yzzz, g_y_z_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxy_xxxx[k] = -g_y_z_xy_xxxx[k] * ab_x + g_y_z_xy_xxxxx[k];

                g_y_z_xxy_xxxy[k] = -g_y_z_xy_xxxy[k] * ab_x + g_y_z_xy_xxxxy[k];

                g_y_z_xxy_xxxz[k] = -g_y_z_xy_xxxz[k] * ab_x + g_y_z_xy_xxxxz[k];

                g_y_z_xxy_xxyy[k] = -g_y_z_xy_xxyy[k] * ab_x + g_y_z_xy_xxxyy[k];

                g_y_z_xxy_xxyz[k] = -g_y_z_xy_xxyz[k] * ab_x + g_y_z_xy_xxxyz[k];

                g_y_z_xxy_xxzz[k] = -g_y_z_xy_xxzz[k] * ab_x + g_y_z_xy_xxxzz[k];

                g_y_z_xxy_xyyy[k] = -g_y_z_xy_xyyy[k] * ab_x + g_y_z_xy_xxyyy[k];

                g_y_z_xxy_xyyz[k] = -g_y_z_xy_xyyz[k] * ab_x + g_y_z_xy_xxyyz[k];

                g_y_z_xxy_xyzz[k] = -g_y_z_xy_xyzz[k] * ab_x + g_y_z_xy_xxyzz[k];

                g_y_z_xxy_xzzz[k] = -g_y_z_xy_xzzz[k] * ab_x + g_y_z_xy_xxzzz[k];

                g_y_z_xxy_yyyy[k] = -g_y_z_xy_yyyy[k] * ab_x + g_y_z_xy_xyyyy[k];

                g_y_z_xxy_yyyz[k] = -g_y_z_xy_yyyz[k] * ab_x + g_y_z_xy_xyyyz[k];

                g_y_z_xxy_yyzz[k] = -g_y_z_xy_yyzz[k] * ab_x + g_y_z_xy_xyyzz[k];

                g_y_z_xxy_yzzz[k] = -g_y_z_xy_yzzz[k] * ab_x + g_y_z_xy_xyzzz[k];

                g_y_z_xxy_zzzz[k] = -g_y_z_xy_zzzz[k] * ab_x + g_y_z_xy_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxz_xxxx = cbuffer.data(fg_geom_11_off + 780 * ccomps * dcomps);

            auto g_y_z_xxz_xxxy = cbuffer.data(fg_geom_11_off + 781 * ccomps * dcomps);

            auto g_y_z_xxz_xxxz = cbuffer.data(fg_geom_11_off + 782 * ccomps * dcomps);

            auto g_y_z_xxz_xxyy = cbuffer.data(fg_geom_11_off + 783 * ccomps * dcomps);

            auto g_y_z_xxz_xxyz = cbuffer.data(fg_geom_11_off + 784 * ccomps * dcomps);

            auto g_y_z_xxz_xxzz = cbuffer.data(fg_geom_11_off + 785 * ccomps * dcomps);

            auto g_y_z_xxz_xyyy = cbuffer.data(fg_geom_11_off + 786 * ccomps * dcomps);

            auto g_y_z_xxz_xyyz = cbuffer.data(fg_geom_11_off + 787 * ccomps * dcomps);

            auto g_y_z_xxz_xyzz = cbuffer.data(fg_geom_11_off + 788 * ccomps * dcomps);

            auto g_y_z_xxz_xzzz = cbuffer.data(fg_geom_11_off + 789 * ccomps * dcomps);

            auto g_y_z_xxz_yyyy = cbuffer.data(fg_geom_11_off + 790 * ccomps * dcomps);

            auto g_y_z_xxz_yyyz = cbuffer.data(fg_geom_11_off + 791 * ccomps * dcomps);

            auto g_y_z_xxz_yyzz = cbuffer.data(fg_geom_11_off + 792 * ccomps * dcomps);

            auto g_y_z_xxz_yzzz = cbuffer.data(fg_geom_11_off + 793 * ccomps * dcomps);

            auto g_y_z_xxz_zzzz = cbuffer.data(fg_geom_11_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxz_xxxx, g_y_z_xxz_xxxy, g_y_z_xxz_xxxz, g_y_z_xxz_xxyy, g_y_z_xxz_xxyz, g_y_z_xxz_xxzz, g_y_z_xxz_xyyy, g_y_z_xxz_xyyz, g_y_z_xxz_xyzz, g_y_z_xxz_xzzz, g_y_z_xxz_yyyy, g_y_z_xxz_yyyz, g_y_z_xxz_yyzz, g_y_z_xxz_yzzz, g_y_z_xxz_zzzz, g_y_z_xz_xxxx, g_y_z_xz_xxxxx, g_y_z_xz_xxxxy, g_y_z_xz_xxxxz, g_y_z_xz_xxxy, g_y_z_xz_xxxyy, g_y_z_xz_xxxyz, g_y_z_xz_xxxz, g_y_z_xz_xxxzz, g_y_z_xz_xxyy, g_y_z_xz_xxyyy, g_y_z_xz_xxyyz, g_y_z_xz_xxyz, g_y_z_xz_xxyzz, g_y_z_xz_xxzz, g_y_z_xz_xxzzz, g_y_z_xz_xyyy, g_y_z_xz_xyyyy, g_y_z_xz_xyyyz, g_y_z_xz_xyyz, g_y_z_xz_xyyzz, g_y_z_xz_xyzz, g_y_z_xz_xyzzz, g_y_z_xz_xzzz, g_y_z_xz_xzzzz, g_y_z_xz_yyyy, g_y_z_xz_yyyz, g_y_z_xz_yyzz, g_y_z_xz_yzzz, g_y_z_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxz_xxxx[k] = -g_y_z_xz_xxxx[k] * ab_x + g_y_z_xz_xxxxx[k];

                g_y_z_xxz_xxxy[k] = -g_y_z_xz_xxxy[k] * ab_x + g_y_z_xz_xxxxy[k];

                g_y_z_xxz_xxxz[k] = -g_y_z_xz_xxxz[k] * ab_x + g_y_z_xz_xxxxz[k];

                g_y_z_xxz_xxyy[k] = -g_y_z_xz_xxyy[k] * ab_x + g_y_z_xz_xxxyy[k];

                g_y_z_xxz_xxyz[k] = -g_y_z_xz_xxyz[k] * ab_x + g_y_z_xz_xxxyz[k];

                g_y_z_xxz_xxzz[k] = -g_y_z_xz_xxzz[k] * ab_x + g_y_z_xz_xxxzz[k];

                g_y_z_xxz_xyyy[k] = -g_y_z_xz_xyyy[k] * ab_x + g_y_z_xz_xxyyy[k];

                g_y_z_xxz_xyyz[k] = -g_y_z_xz_xyyz[k] * ab_x + g_y_z_xz_xxyyz[k];

                g_y_z_xxz_xyzz[k] = -g_y_z_xz_xyzz[k] * ab_x + g_y_z_xz_xxyzz[k];

                g_y_z_xxz_xzzz[k] = -g_y_z_xz_xzzz[k] * ab_x + g_y_z_xz_xxzzz[k];

                g_y_z_xxz_yyyy[k] = -g_y_z_xz_yyyy[k] * ab_x + g_y_z_xz_xyyyy[k];

                g_y_z_xxz_yyyz[k] = -g_y_z_xz_yyyz[k] * ab_x + g_y_z_xz_xyyyz[k];

                g_y_z_xxz_yyzz[k] = -g_y_z_xz_yyzz[k] * ab_x + g_y_z_xz_xyyzz[k];

                g_y_z_xxz_yzzz[k] = -g_y_z_xz_yzzz[k] * ab_x + g_y_z_xz_xyzzz[k];

                g_y_z_xxz_zzzz[k] = -g_y_z_xz_zzzz[k] * ab_x + g_y_z_xz_xzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyy_xxxx = cbuffer.data(fg_geom_11_off + 795 * ccomps * dcomps);

            auto g_y_z_xyy_xxxy = cbuffer.data(fg_geom_11_off + 796 * ccomps * dcomps);

            auto g_y_z_xyy_xxxz = cbuffer.data(fg_geom_11_off + 797 * ccomps * dcomps);

            auto g_y_z_xyy_xxyy = cbuffer.data(fg_geom_11_off + 798 * ccomps * dcomps);

            auto g_y_z_xyy_xxyz = cbuffer.data(fg_geom_11_off + 799 * ccomps * dcomps);

            auto g_y_z_xyy_xxzz = cbuffer.data(fg_geom_11_off + 800 * ccomps * dcomps);

            auto g_y_z_xyy_xyyy = cbuffer.data(fg_geom_11_off + 801 * ccomps * dcomps);

            auto g_y_z_xyy_xyyz = cbuffer.data(fg_geom_11_off + 802 * ccomps * dcomps);

            auto g_y_z_xyy_xyzz = cbuffer.data(fg_geom_11_off + 803 * ccomps * dcomps);

            auto g_y_z_xyy_xzzz = cbuffer.data(fg_geom_11_off + 804 * ccomps * dcomps);

            auto g_y_z_xyy_yyyy = cbuffer.data(fg_geom_11_off + 805 * ccomps * dcomps);

            auto g_y_z_xyy_yyyz = cbuffer.data(fg_geom_11_off + 806 * ccomps * dcomps);

            auto g_y_z_xyy_yyzz = cbuffer.data(fg_geom_11_off + 807 * ccomps * dcomps);

            auto g_y_z_xyy_yzzz = cbuffer.data(fg_geom_11_off + 808 * ccomps * dcomps);

            auto g_y_z_xyy_zzzz = cbuffer.data(fg_geom_11_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyy_xxxx, g_y_z_xyy_xxxy, g_y_z_xyy_xxxz, g_y_z_xyy_xxyy, g_y_z_xyy_xxyz, g_y_z_xyy_xxzz, g_y_z_xyy_xyyy, g_y_z_xyy_xyyz, g_y_z_xyy_xyzz, g_y_z_xyy_xzzz, g_y_z_xyy_yyyy, g_y_z_xyy_yyyz, g_y_z_xyy_yyzz, g_y_z_xyy_yzzz, g_y_z_xyy_zzzz, g_y_z_yy_xxxx, g_y_z_yy_xxxxx, g_y_z_yy_xxxxy, g_y_z_yy_xxxxz, g_y_z_yy_xxxy, g_y_z_yy_xxxyy, g_y_z_yy_xxxyz, g_y_z_yy_xxxz, g_y_z_yy_xxxzz, g_y_z_yy_xxyy, g_y_z_yy_xxyyy, g_y_z_yy_xxyyz, g_y_z_yy_xxyz, g_y_z_yy_xxyzz, g_y_z_yy_xxzz, g_y_z_yy_xxzzz, g_y_z_yy_xyyy, g_y_z_yy_xyyyy, g_y_z_yy_xyyyz, g_y_z_yy_xyyz, g_y_z_yy_xyyzz, g_y_z_yy_xyzz, g_y_z_yy_xyzzz, g_y_z_yy_xzzz, g_y_z_yy_xzzzz, g_y_z_yy_yyyy, g_y_z_yy_yyyz, g_y_z_yy_yyzz, g_y_z_yy_yzzz, g_y_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyy_xxxx[k] = -g_y_z_yy_xxxx[k] * ab_x + g_y_z_yy_xxxxx[k];

                g_y_z_xyy_xxxy[k] = -g_y_z_yy_xxxy[k] * ab_x + g_y_z_yy_xxxxy[k];

                g_y_z_xyy_xxxz[k] = -g_y_z_yy_xxxz[k] * ab_x + g_y_z_yy_xxxxz[k];

                g_y_z_xyy_xxyy[k] = -g_y_z_yy_xxyy[k] * ab_x + g_y_z_yy_xxxyy[k];

                g_y_z_xyy_xxyz[k] = -g_y_z_yy_xxyz[k] * ab_x + g_y_z_yy_xxxyz[k];

                g_y_z_xyy_xxzz[k] = -g_y_z_yy_xxzz[k] * ab_x + g_y_z_yy_xxxzz[k];

                g_y_z_xyy_xyyy[k] = -g_y_z_yy_xyyy[k] * ab_x + g_y_z_yy_xxyyy[k];

                g_y_z_xyy_xyyz[k] = -g_y_z_yy_xyyz[k] * ab_x + g_y_z_yy_xxyyz[k];

                g_y_z_xyy_xyzz[k] = -g_y_z_yy_xyzz[k] * ab_x + g_y_z_yy_xxyzz[k];

                g_y_z_xyy_xzzz[k] = -g_y_z_yy_xzzz[k] * ab_x + g_y_z_yy_xxzzz[k];

                g_y_z_xyy_yyyy[k] = -g_y_z_yy_yyyy[k] * ab_x + g_y_z_yy_xyyyy[k];

                g_y_z_xyy_yyyz[k] = -g_y_z_yy_yyyz[k] * ab_x + g_y_z_yy_xyyyz[k];

                g_y_z_xyy_yyzz[k] = -g_y_z_yy_yyzz[k] * ab_x + g_y_z_yy_xyyzz[k];

                g_y_z_xyy_yzzz[k] = -g_y_z_yy_yzzz[k] * ab_x + g_y_z_yy_xyzzz[k];

                g_y_z_xyy_zzzz[k] = -g_y_z_yy_zzzz[k] * ab_x + g_y_z_yy_xzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyz_xxxx = cbuffer.data(fg_geom_11_off + 810 * ccomps * dcomps);

            auto g_y_z_xyz_xxxy = cbuffer.data(fg_geom_11_off + 811 * ccomps * dcomps);

            auto g_y_z_xyz_xxxz = cbuffer.data(fg_geom_11_off + 812 * ccomps * dcomps);

            auto g_y_z_xyz_xxyy = cbuffer.data(fg_geom_11_off + 813 * ccomps * dcomps);

            auto g_y_z_xyz_xxyz = cbuffer.data(fg_geom_11_off + 814 * ccomps * dcomps);

            auto g_y_z_xyz_xxzz = cbuffer.data(fg_geom_11_off + 815 * ccomps * dcomps);

            auto g_y_z_xyz_xyyy = cbuffer.data(fg_geom_11_off + 816 * ccomps * dcomps);

            auto g_y_z_xyz_xyyz = cbuffer.data(fg_geom_11_off + 817 * ccomps * dcomps);

            auto g_y_z_xyz_xyzz = cbuffer.data(fg_geom_11_off + 818 * ccomps * dcomps);

            auto g_y_z_xyz_xzzz = cbuffer.data(fg_geom_11_off + 819 * ccomps * dcomps);

            auto g_y_z_xyz_yyyy = cbuffer.data(fg_geom_11_off + 820 * ccomps * dcomps);

            auto g_y_z_xyz_yyyz = cbuffer.data(fg_geom_11_off + 821 * ccomps * dcomps);

            auto g_y_z_xyz_yyzz = cbuffer.data(fg_geom_11_off + 822 * ccomps * dcomps);

            auto g_y_z_xyz_yzzz = cbuffer.data(fg_geom_11_off + 823 * ccomps * dcomps);

            auto g_y_z_xyz_zzzz = cbuffer.data(fg_geom_11_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyz_xxxx, g_y_z_xyz_xxxy, g_y_z_xyz_xxxz, g_y_z_xyz_xxyy, g_y_z_xyz_xxyz, g_y_z_xyz_xxzz, g_y_z_xyz_xyyy, g_y_z_xyz_xyyz, g_y_z_xyz_xyzz, g_y_z_xyz_xzzz, g_y_z_xyz_yyyy, g_y_z_xyz_yyyz, g_y_z_xyz_yyzz, g_y_z_xyz_yzzz, g_y_z_xyz_zzzz, g_y_z_yz_xxxx, g_y_z_yz_xxxxx, g_y_z_yz_xxxxy, g_y_z_yz_xxxxz, g_y_z_yz_xxxy, g_y_z_yz_xxxyy, g_y_z_yz_xxxyz, g_y_z_yz_xxxz, g_y_z_yz_xxxzz, g_y_z_yz_xxyy, g_y_z_yz_xxyyy, g_y_z_yz_xxyyz, g_y_z_yz_xxyz, g_y_z_yz_xxyzz, g_y_z_yz_xxzz, g_y_z_yz_xxzzz, g_y_z_yz_xyyy, g_y_z_yz_xyyyy, g_y_z_yz_xyyyz, g_y_z_yz_xyyz, g_y_z_yz_xyyzz, g_y_z_yz_xyzz, g_y_z_yz_xyzzz, g_y_z_yz_xzzz, g_y_z_yz_xzzzz, g_y_z_yz_yyyy, g_y_z_yz_yyyz, g_y_z_yz_yyzz, g_y_z_yz_yzzz, g_y_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyz_xxxx[k] = -g_y_z_yz_xxxx[k] * ab_x + g_y_z_yz_xxxxx[k];

                g_y_z_xyz_xxxy[k] = -g_y_z_yz_xxxy[k] * ab_x + g_y_z_yz_xxxxy[k];

                g_y_z_xyz_xxxz[k] = -g_y_z_yz_xxxz[k] * ab_x + g_y_z_yz_xxxxz[k];

                g_y_z_xyz_xxyy[k] = -g_y_z_yz_xxyy[k] * ab_x + g_y_z_yz_xxxyy[k];

                g_y_z_xyz_xxyz[k] = -g_y_z_yz_xxyz[k] * ab_x + g_y_z_yz_xxxyz[k];

                g_y_z_xyz_xxzz[k] = -g_y_z_yz_xxzz[k] * ab_x + g_y_z_yz_xxxzz[k];

                g_y_z_xyz_xyyy[k] = -g_y_z_yz_xyyy[k] * ab_x + g_y_z_yz_xxyyy[k];

                g_y_z_xyz_xyyz[k] = -g_y_z_yz_xyyz[k] * ab_x + g_y_z_yz_xxyyz[k];

                g_y_z_xyz_xyzz[k] = -g_y_z_yz_xyzz[k] * ab_x + g_y_z_yz_xxyzz[k];

                g_y_z_xyz_xzzz[k] = -g_y_z_yz_xzzz[k] * ab_x + g_y_z_yz_xxzzz[k];

                g_y_z_xyz_yyyy[k] = -g_y_z_yz_yyyy[k] * ab_x + g_y_z_yz_xyyyy[k];

                g_y_z_xyz_yyyz[k] = -g_y_z_yz_yyyz[k] * ab_x + g_y_z_yz_xyyyz[k];

                g_y_z_xyz_yyzz[k] = -g_y_z_yz_yyzz[k] * ab_x + g_y_z_yz_xyyzz[k];

                g_y_z_xyz_yzzz[k] = -g_y_z_yz_yzzz[k] * ab_x + g_y_z_yz_xyzzz[k];

                g_y_z_xyz_zzzz[k] = -g_y_z_yz_zzzz[k] * ab_x + g_y_z_yz_xzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzz_xxxx = cbuffer.data(fg_geom_11_off + 825 * ccomps * dcomps);

            auto g_y_z_xzz_xxxy = cbuffer.data(fg_geom_11_off + 826 * ccomps * dcomps);

            auto g_y_z_xzz_xxxz = cbuffer.data(fg_geom_11_off + 827 * ccomps * dcomps);

            auto g_y_z_xzz_xxyy = cbuffer.data(fg_geom_11_off + 828 * ccomps * dcomps);

            auto g_y_z_xzz_xxyz = cbuffer.data(fg_geom_11_off + 829 * ccomps * dcomps);

            auto g_y_z_xzz_xxzz = cbuffer.data(fg_geom_11_off + 830 * ccomps * dcomps);

            auto g_y_z_xzz_xyyy = cbuffer.data(fg_geom_11_off + 831 * ccomps * dcomps);

            auto g_y_z_xzz_xyyz = cbuffer.data(fg_geom_11_off + 832 * ccomps * dcomps);

            auto g_y_z_xzz_xyzz = cbuffer.data(fg_geom_11_off + 833 * ccomps * dcomps);

            auto g_y_z_xzz_xzzz = cbuffer.data(fg_geom_11_off + 834 * ccomps * dcomps);

            auto g_y_z_xzz_yyyy = cbuffer.data(fg_geom_11_off + 835 * ccomps * dcomps);

            auto g_y_z_xzz_yyyz = cbuffer.data(fg_geom_11_off + 836 * ccomps * dcomps);

            auto g_y_z_xzz_yyzz = cbuffer.data(fg_geom_11_off + 837 * ccomps * dcomps);

            auto g_y_z_xzz_yzzz = cbuffer.data(fg_geom_11_off + 838 * ccomps * dcomps);

            auto g_y_z_xzz_zzzz = cbuffer.data(fg_geom_11_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzz_xxxx, g_y_z_xzz_xxxy, g_y_z_xzz_xxxz, g_y_z_xzz_xxyy, g_y_z_xzz_xxyz, g_y_z_xzz_xxzz, g_y_z_xzz_xyyy, g_y_z_xzz_xyyz, g_y_z_xzz_xyzz, g_y_z_xzz_xzzz, g_y_z_xzz_yyyy, g_y_z_xzz_yyyz, g_y_z_xzz_yyzz, g_y_z_xzz_yzzz, g_y_z_xzz_zzzz, g_y_z_zz_xxxx, g_y_z_zz_xxxxx, g_y_z_zz_xxxxy, g_y_z_zz_xxxxz, g_y_z_zz_xxxy, g_y_z_zz_xxxyy, g_y_z_zz_xxxyz, g_y_z_zz_xxxz, g_y_z_zz_xxxzz, g_y_z_zz_xxyy, g_y_z_zz_xxyyy, g_y_z_zz_xxyyz, g_y_z_zz_xxyz, g_y_z_zz_xxyzz, g_y_z_zz_xxzz, g_y_z_zz_xxzzz, g_y_z_zz_xyyy, g_y_z_zz_xyyyy, g_y_z_zz_xyyyz, g_y_z_zz_xyyz, g_y_z_zz_xyyzz, g_y_z_zz_xyzz, g_y_z_zz_xyzzz, g_y_z_zz_xzzz, g_y_z_zz_xzzzz, g_y_z_zz_yyyy, g_y_z_zz_yyyz, g_y_z_zz_yyzz, g_y_z_zz_yzzz, g_y_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzz_xxxx[k] = -g_y_z_zz_xxxx[k] * ab_x + g_y_z_zz_xxxxx[k];

                g_y_z_xzz_xxxy[k] = -g_y_z_zz_xxxy[k] * ab_x + g_y_z_zz_xxxxy[k];

                g_y_z_xzz_xxxz[k] = -g_y_z_zz_xxxz[k] * ab_x + g_y_z_zz_xxxxz[k];

                g_y_z_xzz_xxyy[k] = -g_y_z_zz_xxyy[k] * ab_x + g_y_z_zz_xxxyy[k];

                g_y_z_xzz_xxyz[k] = -g_y_z_zz_xxyz[k] * ab_x + g_y_z_zz_xxxyz[k];

                g_y_z_xzz_xxzz[k] = -g_y_z_zz_xxzz[k] * ab_x + g_y_z_zz_xxxzz[k];

                g_y_z_xzz_xyyy[k] = -g_y_z_zz_xyyy[k] * ab_x + g_y_z_zz_xxyyy[k];

                g_y_z_xzz_xyyz[k] = -g_y_z_zz_xyyz[k] * ab_x + g_y_z_zz_xxyyz[k];

                g_y_z_xzz_xyzz[k] = -g_y_z_zz_xyzz[k] * ab_x + g_y_z_zz_xxyzz[k];

                g_y_z_xzz_xzzz[k] = -g_y_z_zz_xzzz[k] * ab_x + g_y_z_zz_xxzzz[k];

                g_y_z_xzz_yyyy[k] = -g_y_z_zz_yyyy[k] * ab_x + g_y_z_zz_xyyyy[k];

                g_y_z_xzz_yyyz[k] = -g_y_z_zz_yyyz[k] * ab_x + g_y_z_zz_xyyyz[k];

                g_y_z_xzz_yyzz[k] = -g_y_z_zz_yyzz[k] * ab_x + g_y_z_zz_xyyzz[k];

                g_y_z_xzz_yzzz[k] = -g_y_z_zz_yzzz[k] * ab_x + g_y_z_zz_xyzzz[k];

                g_y_z_xzz_zzzz[k] = -g_y_z_zz_zzzz[k] * ab_x + g_y_z_zz_xzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyy_xxxx = cbuffer.data(fg_geom_11_off + 840 * ccomps * dcomps);

            auto g_y_z_yyy_xxxy = cbuffer.data(fg_geom_11_off + 841 * ccomps * dcomps);

            auto g_y_z_yyy_xxxz = cbuffer.data(fg_geom_11_off + 842 * ccomps * dcomps);

            auto g_y_z_yyy_xxyy = cbuffer.data(fg_geom_11_off + 843 * ccomps * dcomps);

            auto g_y_z_yyy_xxyz = cbuffer.data(fg_geom_11_off + 844 * ccomps * dcomps);

            auto g_y_z_yyy_xxzz = cbuffer.data(fg_geom_11_off + 845 * ccomps * dcomps);

            auto g_y_z_yyy_xyyy = cbuffer.data(fg_geom_11_off + 846 * ccomps * dcomps);

            auto g_y_z_yyy_xyyz = cbuffer.data(fg_geom_11_off + 847 * ccomps * dcomps);

            auto g_y_z_yyy_xyzz = cbuffer.data(fg_geom_11_off + 848 * ccomps * dcomps);

            auto g_y_z_yyy_xzzz = cbuffer.data(fg_geom_11_off + 849 * ccomps * dcomps);

            auto g_y_z_yyy_yyyy = cbuffer.data(fg_geom_11_off + 850 * ccomps * dcomps);

            auto g_y_z_yyy_yyyz = cbuffer.data(fg_geom_11_off + 851 * ccomps * dcomps);

            auto g_y_z_yyy_yyzz = cbuffer.data(fg_geom_11_off + 852 * ccomps * dcomps);

            auto g_y_z_yyy_yzzz = cbuffer.data(fg_geom_11_off + 853 * ccomps * dcomps);

            auto g_y_z_yyy_zzzz = cbuffer.data(fg_geom_11_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_xxxx, g_0_z_yy_xxxy, g_0_z_yy_xxxz, g_0_z_yy_xxyy, g_0_z_yy_xxyz, g_0_z_yy_xxzz, g_0_z_yy_xyyy, g_0_z_yy_xyyz, g_0_z_yy_xyzz, g_0_z_yy_xzzz, g_0_z_yy_yyyy, g_0_z_yy_yyyz, g_0_z_yy_yyzz, g_0_z_yy_yzzz, g_0_z_yy_zzzz, g_y_z_yy_xxxx, g_y_z_yy_xxxxy, g_y_z_yy_xxxy, g_y_z_yy_xxxyy, g_y_z_yy_xxxyz, g_y_z_yy_xxxz, g_y_z_yy_xxyy, g_y_z_yy_xxyyy, g_y_z_yy_xxyyz, g_y_z_yy_xxyz, g_y_z_yy_xxyzz, g_y_z_yy_xxzz, g_y_z_yy_xyyy, g_y_z_yy_xyyyy, g_y_z_yy_xyyyz, g_y_z_yy_xyyz, g_y_z_yy_xyyzz, g_y_z_yy_xyzz, g_y_z_yy_xyzzz, g_y_z_yy_xzzz, g_y_z_yy_yyyy, g_y_z_yy_yyyyy, g_y_z_yy_yyyyz, g_y_z_yy_yyyz, g_y_z_yy_yyyzz, g_y_z_yy_yyzz, g_y_z_yy_yyzzz, g_y_z_yy_yzzz, g_y_z_yy_yzzzz, g_y_z_yy_zzzz, g_y_z_yyy_xxxx, g_y_z_yyy_xxxy, g_y_z_yyy_xxxz, g_y_z_yyy_xxyy, g_y_z_yyy_xxyz, g_y_z_yyy_xxzz, g_y_z_yyy_xyyy, g_y_z_yyy_xyyz, g_y_z_yyy_xyzz, g_y_z_yyy_xzzz, g_y_z_yyy_yyyy, g_y_z_yyy_yyyz, g_y_z_yyy_yyzz, g_y_z_yyy_yzzz, g_y_z_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyy_xxxx[k] = -g_0_z_yy_xxxx[k] - g_y_z_yy_xxxx[k] * ab_y + g_y_z_yy_xxxxy[k];

                g_y_z_yyy_xxxy[k] = -g_0_z_yy_xxxy[k] - g_y_z_yy_xxxy[k] * ab_y + g_y_z_yy_xxxyy[k];

                g_y_z_yyy_xxxz[k] = -g_0_z_yy_xxxz[k] - g_y_z_yy_xxxz[k] * ab_y + g_y_z_yy_xxxyz[k];

                g_y_z_yyy_xxyy[k] = -g_0_z_yy_xxyy[k] - g_y_z_yy_xxyy[k] * ab_y + g_y_z_yy_xxyyy[k];

                g_y_z_yyy_xxyz[k] = -g_0_z_yy_xxyz[k] - g_y_z_yy_xxyz[k] * ab_y + g_y_z_yy_xxyyz[k];

                g_y_z_yyy_xxzz[k] = -g_0_z_yy_xxzz[k] - g_y_z_yy_xxzz[k] * ab_y + g_y_z_yy_xxyzz[k];

                g_y_z_yyy_xyyy[k] = -g_0_z_yy_xyyy[k] - g_y_z_yy_xyyy[k] * ab_y + g_y_z_yy_xyyyy[k];

                g_y_z_yyy_xyyz[k] = -g_0_z_yy_xyyz[k] - g_y_z_yy_xyyz[k] * ab_y + g_y_z_yy_xyyyz[k];

                g_y_z_yyy_xyzz[k] = -g_0_z_yy_xyzz[k] - g_y_z_yy_xyzz[k] * ab_y + g_y_z_yy_xyyzz[k];

                g_y_z_yyy_xzzz[k] = -g_0_z_yy_xzzz[k] - g_y_z_yy_xzzz[k] * ab_y + g_y_z_yy_xyzzz[k];

                g_y_z_yyy_yyyy[k] = -g_0_z_yy_yyyy[k] - g_y_z_yy_yyyy[k] * ab_y + g_y_z_yy_yyyyy[k];

                g_y_z_yyy_yyyz[k] = -g_0_z_yy_yyyz[k] - g_y_z_yy_yyyz[k] * ab_y + g_y_z_yy_yyyyz[k];

                g_y_z_yyy_yyzz[k] = -g_0_z_yy_yyzz[k] - g_y_z_yy_yyzz[k] * ab_y + g_y_z_yy_yyyzz[k];

                g_y_z_yyy_yzzz[k] = -g_0_z_yy_yzzz[k] - g_y_z_yy_yzzz[k] * ab_y + g_y_z_yy_yyzzz[k];

                g_y_z_yyy_zzzz[k] = -g_0_z_yy_zzzz[k] - g_y_z_yy_zzzz[k] * ab_y + g_y_z_yy_yzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyz_xxxx = cbuffer.data(fg_geom_11_off + 855 * ccomps * dcomps);

            auto g_y_z_yyz_xxxy = cbuffer.data(fg_geom_11_off + 856 * ccomps * dcomps);

            auto g_y_z_yyz_xxxz = cbuffer.data(fg_geom_11_off + 857 * ccomps * dcomps);

            auto g_y_z_yyz_xxyy = cbuffer.data(fg_geom_11_off + 858 * ccomps * dcomps);

            auto g_y_z_yyz_xxyz = cbuffer.data(fg_geom_11_off + 859 * ccomps * dcomps);

            auto g_y_z_yyz_xxzz = cbuffer.data(fg_geom_11_off + 860 * ccomps * dcomps);

            auto g_y_z_yyz_xyyy = cbuffer.data(fg_geom_11_off + 861 * ccomps * dcomps);

            auto g_y_z_yyz_xyyz = cbuffer.data(fg_geom_11_off + 862 * ccomps * dcomps);

            auto g_y_z_yyz_xyzz = cbuffer.data(fg_geom_11_off + 863 * ccomps * dcomps);

            auto g_y_z_yyz_xzzz = cbuffer.data(fg_geom_11_off + 864 * ccomps * dcomps);

            auto g_y_z_yyz_yyyy = cbuffer.data(fg_geom_11_off + 865 * ccomps * dcomps);

            auto g_y_z_yyz_yyyz = cbuffer.data(fg_geom_11_off + 866 * ccomps * dcomps);

            auto g_y_z_yyz_yyzz = cbuffer.data(fg_geom_11_off + 867 * ccomps * dcomps);

            auto g_y_z_yyz_yzzz = cbuffer.data(fg_geom_11_off + 868 * ccomps * dcomps);

            auto g_y_z_yyz_zzzz = cbuffer.data(fg_geom_11_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_xxxx, g_0_z_yz_xxxy, g_0_z_yz_xxxz, g_0_z_yz_xxyy, g_0_z_yz_xxyz, g_0_z_yz_xxzz, g_0_z_yz_xyyy, g_0_z_yz_xyyz, g_0_z_yz_xyzz, g_0_z_yz_xzzz, g_0_z_yz_yyyy, g_0_z_yz_yyyz, g_0_z_yz_yyzz, g_0_z_yz_yzzz, g_0_z_yz_zzzz, g_y_z_yyz_xxxx, g_y_z_yyz_xxxy, g_y_z_yyz_xxxz, g_y_z_yyz_xxyy, g_y_z_yyz_xxyz, g_y_z_yyz_xxzz, g_y_z_yyz_xyyy, g_y_z_yyz_xyyz, g_y_z_yyz_xyzz, g_y_z_yyz_xzzz, g_y_z_yyz_yyyy, g_y_z_yyz_yyyz, g_y_z_yyz_yyzz, g_y_z_yyz_yzzz, g_y_z_yyz_zzzz, g_y_z_yz_xxxx, g_y_z_yz_xxxxy, g_y_z_yz_xxxy, g_y_z_yz_xxxyy, g_y_z_yz_xxxyz, g_y_z_yz_xxxz, g_y_z_yz_xxyy, g_y_z_yz_xxyyy, g_y_z_yz_xxyyz, g_y_z_yz_xxyz, g_y_z_yz_xxyzz, g_y_z_yz_xxzz, g_y_z_yz_xyyy, g_y_z_yz_xyyyy, g_y_z_yz_xyyyz, g_y_z_yz_xyyz, g_y_z_yz_xyyzz, g_y_z_yz_xyzz, g_y_z_yz_xyzzz, g_y_z_yz_xzzz, g_y_z_yz_yyyy, g_y_z_yz_yyyyy, g_y_z_yz_yyyyz, g_y_z_yz_yyyz, g_y_z_yz_yyyzz, g_y_z_yz_yyzz, g_y_z_yz_yyzzz, g_y_z_yz_yzzz, g_y_z_yz_yzzzz, g_y_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyz_xxxx[k] = -g_0_z_yz_xxxx[k] - g_y_z_yz_xxxx[k] * ab_y + g_y_z_yz_xxxxy[k];

                g_y_z_yyz_xxxy[k] = -g_0_z_yz_xxxy[k] - g_y_z_yz_xxxy[k] * ab_y + g_y_z_yz_xxxyy[k];

                g_y_z_yyz_xxxz[k] = -g_0_z_yz_xxxz[k] - g_y_z_yz_xxxz[k] * ab_y + g_y_z_yz_xxxyz[k];

                g_y_z_yyz_xxyy[k] = -g_0_z_yz_xxyy[k] - g_y_z_yz_xxyy[k] * ab_y + g_y_z_yz_xxyyy[k];

                g_y_z_yyz_xxyz[k] = -g_0_z_yz_xxyz[k] - g_y_z_yz_xxyz[k] * ab_y + g_y_z_yz_xxyyz[k];

                g_y_z_yyz_xxzz[k] = -g_0_z_yz_xxzz[k] - g_y_z_yz_xxzz[k] * ab_y + g_y_z_yz_xxyzz[k];

                g_y_z_yyz_xyyy[k] = -g_0_z_yz_xyyy[k] - g_y_z_yz_xyyy[k] * ab_y + g_y_z_yz_xyyyy[k];

                g_y_z_yyz_xyyz[k] = -g_0_z_yz_xyyz[k] - g_y_z_yz_xyyz[k] * ab_y + g_y_z_yz_xyyyz[k];

                g_y_z_yyz_xyzz[k] = -g_0_z_yz_xyzz[k] - g_y_z_yz_xyzz[k] * ab_y + g_y_z_yz_xyyzz[k];

                g_y_z_yyz_xzzz[k] = -g_0_z_yz_xzzz[k] - g_y_z_yz_xzzz[k] * ab_y + g_y_z_yz_xyzzz[k];

                g_y_z_yyz_yyyy[k] = -g_0_z_yz_yyyy[k] - g_y_z_yz_yyyy[k] * ab_y + g_y_z_yz_yyyyy[k];

                g_y_z_yyz_yyyz[k] = -g_0_z_yz_yyyz[k] - g_y_z_yz_yyyz[k] * ab_y + g_y_z_yz_yyyyz[k];

                g_y_z_yyz_yyzz[k] = -g_0_z_yz_yyzz[k] - g_y_z_yz_yyzz[k] * ab_y + g_y_z_yz_yyyzz[k];

                g_y_z_yyz_yzzz[k] = -g_0_z_yz_yzzz[k] - g_y_z_yz_yzzz[k] * ab_y + g_y_z_yz_yyzzz[k];

                g_y_z_yyz_zzzz[k] = -g_0_z_yz_zzzz[k] - g_y_z_yz_zzzz[k] * ab_y + g_y_z_yz_yzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzz_xxxx = cbuffer.data(fg_geom_11_off + 870 * ccomps * dcomps);

            auto g_y_z_yzz_xxxy = cbuffer.data(fg_geom_11_off + 871 * ccomps * dcomps);

            auto g_y_z_yzz_xxxz = cbuffer.data(fg_geom_11_off + 872 * ccomps * dcomps);

            auto g_y_z_yzz_xxyy = cbuffer.data(fg_geom_11_off + 873 * ccomps * dcomps);

            auto g_y_z_yzz_xxyz = cbuffer.data(fg_geom_11_off + 874 * ccomps * dcomps);

            auto g_y_z_yzz_xxzz = cbuffer.data(fg_geom_11_off + 875 * ccomps * dcomps);

            auto g_y_z_yzz_xyyy = cbuffer.data(fg_geom_11_off + 876 * ccomps * dcomps);

            auto g_y_z_yzz_xyyz = cbuffer.data(fg_geom_11_off + 877 * ccomps * dcomps);

            auto g_y_z_yzz_xyzz = cbuffer.data(fg_geom_11_off + 878 * ccomps * dcomps);

            auto g_y_z_yzz_xzzz = cbuffer.data(fg_geom_11_off + 879 * ccomps * dcomps);

            auto g_y_z_yzz_yyyy = cbuffer.data(fg_geom_11_off + 880 * ccomps * dcomps);

            auto g_y_z_yzz_yyyz = cbuffer.data(fg_geom_11_off + 881 * ccomps * dcomps);

            auto g_y_z_yzz_yyzz = cbuffer.data(fg_geom_11_off + 882 * ccomps * dcomps);

            auto g_y_z_yzz_yzzz = cbuffer.data(fg_geom_11_off + 883 * ccomps * dcomps);

            auto g_y_z_yzz_zzzz = cbuffer.data(fg_geom_11_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxxx, g_0_z_zz_xxxy, g_0_z_zz_xxxz, g_0_z_zz_xxyy, g_0_z_zz_xxyz, g_0_z_zz_xxzz, g_0_z_zz_xyyy, g_0_z_zz_xyyz, g_0_z_zz_xyzz, g_0_z_zz_xzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyz, g_0_z_zz_yyzz, g_0_z_zz_yzzz, g_0_z_zz_zzzz, g_y_z_yzz_xxxx, g_y_z_yzz_xxxy, g_y_z_yzz_xxxz, g_y_z_yzz_xxyy, g_y_z_yzz_xxyz, g_y_z_yzz_xxzz, g_y_z_yzz_xyyy, g_y_z_yzz_xyyz, g_y_z_yzz_xyzz, g_y_z_yzz_xzzz, g_y_z_yzz_yyyy, g_y_z_yzz_yyyz, g_y_z_yzz_yyzz, g_y_z_yzz_yzzz, g_y_z_yzz_zzzz, g_y_z_zz_xxxx, g_y_z_zz_xxxxy, g_y_z_zz_xxxy, g_y_z_zz_xxxyy, g_y_z_zz_xxxyz, g_y_z_zz_xxxz, g_y_z_zz_xxyy, g_y_z_zz_xxyyy, g_y_z_zz_xxyyz, g_y_z_zz_xxyz, g_y_z_zz_xxyzz, g_y_z_zz_xxzz, g_y_z_zz_xyyy, g_y_z_zz_xyyyy, g_y_z_zz_xyyyz, g_y_z_zz_xyyz, g_y_z_zz_xyyzz, g_y_z_zz_xyzz, g_y_z_zz_xyzzz, g_y_z_zz_xzzz, g_y_z_zz_yyyy, g_y_z_zz_yyyyy, g_y_z_zz_yyyyz, g_y_z_zz_yyyz, g_y_z_zz_yyyzz, g_y_z_zz_yyzz, g_y_z_zz_yyzzz, g_y_z_zz_yzzz, g_y_z_zz_yzzzz, g_y_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzz_xxxx[k] = -g_0_z_zz_xxxx[k] - g_y_z_zz_xxxx[k] * ab_y + g_y_z_zz_xxxxy[k];

                g_y_z_yzz_xxxy[k] = -g_0_z_zz_xxxy[k] - g_y_z_zz_xxxy[k] * ab_y + g_y_z_zz_xxxyy[k];

                g_y_z_yzz_xxxz[k] = -g_0_z_zz_xxxz[k] - g_y_z_zz_xxxz[k] * ab_y + g_y_z_zz_xxxyz[k];

                g_y_z_yzz_xxyy[k] = -g_0_z_zz_xxyy[k] - g_y_z_zz_xxyy[k] * ab_y + g_y_z_zz_xxyyy[k];

                g_y_z_yzz_xxyz[k] = -g_0_z_zz_xxyz[k] - g_y_z_zz_xxyz[k] * ab_y + g_y_z_zz_xxyyz[k];

                g_y_z_yzz_xxzz[k] = -g_0_z_zz_xxzz[k] - g_y_z_zz_xxzz[k] * ab_y + g_y_z_zz_xxyzz[k];

                g_y_z_yzz_xyyy[k] = -g_0_z_zz_xyyy[k] - g_y_z_zz_xyyy[k] * ab_y + g_y_z_zz_xyyyy[k];

                g_y_z_yzz_xyyz[k] = -g_0_z_zz_xyyz[k] - g_y_z_zz_xyyz[k] * ab_y + g_y_z_zz_xyyyz[k];

                g_y_z_yzz_xyzz[k] = -g_0_z_zz_xyzz[k] - g_y_z_zz_xyzz[k] * ab_y + g_y_z_zz_xyyzz[k];

                g_y_z_yzz_xzzz[k] = -g_0_z_zz_xzzz[k] - g_y_z_zz_xzzz[k] * ab_y + g_y_z_zz_xyzzz[k];

                g_y_z_yzz_yyyy[k] = -g_0_z_zz_yyyy[k] - g_y_z_zz_yyyy[k] * ab_y + g_y_z_zz_yyyyy[k];

                g_y_z_yzz_yyyz[k] = -g_0_z_zz_yyyz[k] - g_y_z_zz_yyyz[k] * ab_y + g_y_z_zz_yyyyz[k];

                g_y_z_yzz_yyzz[k] = -g_0_z_zz_yyzz[k] - g_y_z_zz_yyzz[k] * ab_y + g_y_z_zz_yyyzz[k];

                g_y_z_yzz_yzzz[k] = -g_0_z_zz_yzzz[k] - g_y_z_zz_yzzz[k] * ab_y + g_y_z_zz_yyzzz[k];

                g_y_z_yzz_zzzz[k] = -g_0_z_zz_zzzz[k] - g_y_z_zz_zzzz[k] * ab_y + g_y_z_zz_yzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzz_xxxx = cbuffer.data(fg_geom_11_off + 885 * ccomps * dcomps);

            auto g_y_z_zzz_xxxy = cbuffer.data(fg_geom_11_off + 886 * ccomps * dcomps);

            auto g_y_z_zzz_xxxz = cbuffer.data(fg_geom_11_off + 887 * ccomps * dcomps);

            auto g_y_z_zzz_xxyy = cbuffer.data(fg_geom_11_off + 888 * ccomps * dcomps);

            auto g_y_z_zzz_xxyz = cbuffer.data(fg_geom_11_off + 889 * ccomps * dcomps);

            auto g_y_z_zzz_xxzz = cbuffer.data(fg_geom_11_off + 890 * ccomps * dcomps);

            auto g_y_z_zzz_xyyy = cbuffer.data(fg_geom_11_off + 891 * ccomps * dcomps);

            auto g_y_z_zzz_xyyz = cbuffer.data(fg_geom_11_off + 892 * ccomps * dcomps);

            auto g_y_z_zzz_xyzz = cbuffer.data(fg_geom_11_off + 893 * ccomps * dcomps);

            auto g_y_z_zzz_xzzz = cbuffer.data(fg_geom_11_off + 894 * ccomps * dcomps);

            auto g_y_z_zzz_yyyy = cbuffer.data(fg_geom_11_off + 895 * ccomps * dcomps);

            auto g_y_z_zzz_yyyz = cbuffer.data(fg_geom_11_off + 896 * ccomps * dcomps);

            auto g_y_z_zzz_yyzz = cbuffer.data(fg_geom_11_off + 897 * ccomps * dcomps);

            auto g_y_z_zzz_yzzz = cbuffer.data(fg_geom_11_off + 898 * ccomps * dcomps);

            auto g_y_z_zzz_zzzz = cbuffer.data(fg_geom_11_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xxxx, g_y_0_zz_xxxy, g_y_0_zz_xxxz, g_y_0_zz_xxyy, g_y_0_zz_xxyz, g_y_0_zz_xxzz, g_y_0_zz_xyyy, g_y_0_zz_xyyz, g_y_0_zz_xyzz, g_y_0_zz_xzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyz, g_y_0_zz_yyzz, g_y_0_zz_yzzz, g_y_0_zz_zzzz, g_y_z_zz_xxxx, g_y_z_zz_xxxxz, g_y_z_zz_xxxy, g_y_z_zz_xxxyz, g_y_z_zz_xxxz, g_y_z_zz_xxxzz, g_y_z_zz_xxyy, g_y_z_zz_xxyyz, g_y_z_zz_xxyz, g_y_z_zz_xxyzz, g_y_z_zz_xxzz, g_y_z_zz_xxzzz, g_y_z_zz_xyyy, g_y_z_zz_xyyyz, g_y_z_zz_xyyz, g_y_z_zz_xyyzz, g_y_z_zz_xyzz, g_y_z_zz_xyzzz, g_y_z_zz_xzzz, g_y_z_zz_xzzzz, g_y_z_zz_yyyy, g_y_z_zz_yyyyz, g_y_z_zz_yyyz, g_y_z_zz_yyyzz, g_y_z_zz_yyzz, g_y_z_zz_yyzzz, g_y_z_zz_yzzz, g_y_z_zz_yzzzz, g_y_z_zz_zzzz, g_y_z_zz_zzzzz, g_y_z_zzz_xxxx, g_y_z_zzz_xxxy, g_y_z_zzz_xxxz, g_y_z_zzz_xxyy, g_y_z_zzz_xxyz, g_y_z_zzz_xxzz, g_y_z_zzz_xyyy, g_y_z_zzz_xyyz, g_y_z_zzz_xyzz, g_y_z_zzz_xzzz, g_y_z_zzz_yyyy, g_y_z_zzz_yyyz, g_y_z_zzz_yyzz, g_y_z_zzz_yzzz, g_y_z_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzz_xxxx[k] = g_y_0_zz_xxxx[k] - g_y_z_zz_xxxx[k] * ab_z + g_y_z_zz_xxxxz[k];

                g_y_z_zzz_xxxy[k] = g_y_0_zz_xxxy[k] - g_y_z_zz_xxxy[k] * ab_z + g_y_z_zz_xxxyz[k];

                g_y_z_zzz_xxxz[k] = g_y_0_zz_xxxz[k] - g_y_z_zz_xxxz[k] * ab_z + g_y_z_zz_xxxzz[k];

                g_y_z_zzz_xxyy[k] = g_y_0_zz_xxyy[k] - g_y_z_zz_xxyy[k] * ab_z + g_y_z_zz_xxyyz[k];

                g_y_z_zzz_xxyz[k] = g_y_0_zz_xxyz[k] - g_y_z_zz_xxyz[k] * ab_z + g_y_z_zz_xxyzz[k];

                g_y_z_zzz_xxzz[k] = g_y_0_zz_xxzz[k] - g_y_z_zz_xxzz[k] * ab_z + g_y_z_zz_xxzzz[k];

                g_y_z_zzz_xyyy[k] = g_y_0_zz_xyyy[k] - g_y_z_zz_xyyy[k] * ab_z + g_y_z_zz_xyyyz[k];

                g_y_z_zzz_xyyz[k] = g_y_0_zz_xyyz[k] - g_y_z_zz_xyyz[k] * ab_z + g_y_z_zz_xyyzz[k];

                g_y_z_zzz_xyzz[k] = g_y_0_zz_xyzz[k] - g_y_z_zz_xyzz[k] * ab_z + g_y_z_zz_xyzzz[k];

                g_y_z_zzz_xzzz[k] = g_y_0_zz_xzzz[k] - g_y_z_zz_xzzz[k] * ab_z + g_y_z_zz_xzzzz[k];

                g_y_z_zzz_yyyy[k] = g_y_0_zz_yyyy[k] - g_y_z_zz_yyyy[k] * ab_z + g_y_z_zz_yyyyz[k];

                g_y_z_zzz_yyyz[k] = g_y_0_zz_yyyz[k] - g_y_z_zz_yyyz[k] * ab_z + g_y_z_zz_yyyzz[k];

                g_y_z_zzz_yyzz[k] = g_y_0_zz_yyzz[k] - g_y_z_zz_yyzz[k] * ab_z + g_y_z_zz_yyzzz[k];

                g_y_z_zzz_yzzz[k] = g_y_0_zz_yzzz[k] - g_y_z_zz_yzzz[k] * ab_z + g_y_z_zz_yzzzz[k];

                g_y_z_zzz_zzzz[k] = g_y_0_zz_zzzz[k] - g_y_z_zz_zzzz[k] * ab_z + g_y_z_zz_zzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxx_xxxx = cbuffer.data(fg_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_x_xxx_xxxy = cbuffer.data(fg_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_x_xxx_xxxz = cbuffer.data(fg_geom_11_off + 902 * ccomps * dcomps);

            auto g_z_x_xxx_xxyy = cbuffer.data(fg_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_x_xxx_xxyz = cbuffer.data(fg_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_x_xxx_xxzz = cbuffer.data(fg_geom_11_off + 905 * ccomps * dcomps);

            auto g_z_x_xxx_xyyy = cbuffer.data(fg_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_x_xxx_xyyz = cbuffer.data(fg_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_x_xxx_xyzz = cbuffer.data(fg_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_x_xxx_xzzz = cbuffer.data(fg_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_x_xxx_yyyy = cbuffer.data(fg_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_x_xxx_yyyz = cbuffer.data(fg_geom_11_off + 911 * ccomps * dcomps);

            auto g_z_x_xxx_yyzz = cbuffer.data(fg_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_x_xxx_yzzz = cbuffer.data(fg_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_x_xxx_zzzz = cbuffer.data(fg_geom_11_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_xxxx, g_z_0_xx_xxxy, g_z_0_xx_xxxz, g_z_0_xx_xxyy, g_z_0_xx_xxyz, g_z_0_xx_xxzz, g_z_0_xx_xyyy, g_z_0_xx_xyyz, g_z_0_xx_xyzz, g_z_0_xx_xzzz, g_z_0_xx_yyyy, g_z_0_xx_yyyz, g_z_0_xx_yyzz, g_z_0_xx_yzzz, g_z_0_xx_zzzz, g_z_x_xx_xxxx, g_z_x_xx_xxxxx, g_z_x_xx_xxxxy, g_z_x_xx_xxxxz, g_z_x_xx_xxxy, g_z_x_xx_xxxyy, g_z_x_xx_xxxyz, g_z_x_xx_xxxz, g_z_x_xx_xxxzz, g_z_x_xx_xxyy, g_z_x_xx_xxyyy, g_z_x_xx_xxyyz, g_z_x_xx_xxyz, g_z_x_xx_xxyzz, g_z_x_xx_xxzz, g_z_x_xx_xxzzz, g_z_x_xx_xyyy, g_z_x_xx_xyyyy, g_z_x_xx_xyyyz, g_z_x_xx_xyyz, g_z_x_xx_xyyzz, g_z_x_xx_xyzz, g_z_x_xx_xyzzz, g_z_x_xx_xzzz, g_z_x_xx_xzzzz, g_z_x_xx_yyyy, g_z_x_xx_yyyz, g_z_x_xx_yyzz, g_z_x_xx_yzzz, g_z_x_xx_zzzz, g_z_x_xxx_xxxx, g_z_x_xxx_xxxy, g_z_x_xxx_xxxz, g_z_x_xxx_xxyy, g_z_x_xxx_xxyz, g_z_x_xxx_xxzz, g_z_x_xxx_xyyy, g_z_x_xxx_xyyz, g_z_x_xxx_xyzz, g_z_x_xxx_xzzz, g_z_x_xxx_yyyy, g_z_x_xxx_yyyz, g_z_x_xxx_yyzz, g_z_x_xxx_yzzz, g_z_x_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxx_xxxx[k] = g_z_0_xx_xxxx[k] - g_z_x_xx_xxxx[k] * ab_x + g_z_x_xx_xxxxx[k];

                g_z_x_xxx_xxxy[k] = g_z_0_xx_xxxy[k] - g_z_x_xx_xxxy[k] * ab_x + g_z_x_xx_xxxxy[k];

                g_z_x_xxx_xxxz[k] = g_z_0_xx_xxxz[k] - g_z_x_xx_xxxz[k] * ab_x + g_z_x_xx_xxxxz[k];

                g_z_x_xxx_xxyy[k] = g_z_0_xx_xxyy[k] - g_z_x_xx_xxyy[k] * ab_x + g_z_x_xx_xxxyy[k];

                g_z_x_xxx_xxyz[k] = g_z_0_xx_xxyz[k] - g_z_x_xx_xxyz[k] * ab_x + g_z_x_xx_xxxyz[k];

                g_z_x_xxx_xxzz[k] = g_z_0_xx_xxzz[k] - g_z_x_xx_xxzz[k] * ab_x + g_z_x_xx_xxxzz[k];

                g_z_x_xxx_xyyy[k] = g_z_0_xx_xyyy[k] - g_z_x_xx_xyyy[k] * ab_x + g_z_x_xx_xxyyy[k];

                g_z_x_xxx_xyyz[k] = g_z_0_xx_xyyz[k] - g_z_x_xx_xyyz[k] * ab_x + g_z_x_xx_xxyyz[k];

                g_z_x_xxx_xyzz[k] = g_z_0_xx_xyzz[k] - g_z_x_xx_xyzz[k] * ab_x + g_z_x_xx_xxyzz[k];

                g_z_x_xxx_xzzz[k] = g_z_0_xx_xzzz[k] - g_z_x_xx_xzzz[k] * ab_x + g_z_x_xx_xxzzz[k];

                g_z_x_xxx_yyyy[k] = g_z_0_xx_yyyy[k] - g_z_x_xx_yyyy[k] * ab_x + g_z_x_xx_xyyyy[k];

                g_z_x_xxx_yyyz[k] = g_z_0_xx_yyyz[k] - g_z_x_xx_yyyz[k] * ab_x + g_z_x_xx_xyyyz[k];

                g_z_x_xxx_yyzz[k] = g_z_0_xx_yyzz[k] - g_z_x_xx_yyzz[k] * ab_x + g_z_x_xx_xyyzz[k];

                g_z_x_xxx_yzzz[k] = g_z_0_xx_yzzz[k] - g_z_x_xx_yzzz[k] * ab_x + g_z_x_xx_xyzzz[k];

                g_z_x_xxx_zzzz[k] = g_z_0_xx_zzzz[k] - g_z_x_xx_zzzz[k] * ab_x + g_z_x_xx_xzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxy_xxxx = cbuffer.data(fg_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_x_xxy_xxxy = cbuffer.data(fg_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_x_xxy_xxxz = cbuffer.data(fg_geom_11_off + 917 * ccomps * dcomps);

            auto g_z_x_xxy_xxyy = cbuffer.data(fg_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_x_xxy_xxyz = cbuffer.data(fg_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_x_xxy_xxzz = cbuffer.data(fg_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_x_xxy_xyyy = cbuffer.data(fg_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_x_xxy_xyyz = cbuffer.data(fg_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_x_xxy_xyzz = cbuffer.data(fg_geom_11_off + 923 * ccomps * dcomps);

            auto g_z_x_xxy_xzzz = cbuffer.data(fg_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_x_xxy_yyyy = cbuffer.data(fg_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_x_xxy_yyyz = cbuffer.data(fg_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_x_xxy_yyzz = cbuffer.data(fg_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_x_xxy_yzzz = cbuffer.data(fg_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_x_xxy_zzzz = cbuffer.data(fg_geom_11_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xx_xxxx, g_z_x_xx_xxxxy, g_z_x_xx_xxxy, g_z_x_xx_xxxyy, g_z_x_xx_xxxyz, g_z_x_xx_xxxz, g_z_x_xx_xxyy, g_z_x_xx_xxyyy, g_z_x_xx_xxyyz, g_z_x_xx_xxyz, g_z_x_xx_xxyzz, g_z_x_xx_xxzz, g_z_x_xx_xyyy, g_z_x_xx_xyyyy, g_z_x_xx_xyyyz, g_z_x_xx_xyyz, g_z_x_xx_xyyzz, g_z_x_xx_xyzz, g_z_x_xx_xyzzz, g_z_x_xx_xzzz, g_z_x_xx_yyyy, g_z_x_xx_yyyyy, g_z_x_xx_yyyyz, g_z_x_xx_yyyz, g_z_x_xx_yyyzz, g_z_x_xx_yyzz, g_z_x_xx_yyzzz, g_z_x_xx_yzzz, g_z_x_xx_yzzzz, g_z_x_xx_zzzz, g_z_x_xxy_xxxx, g_z_x_xxy_xxxy, g_z_x_xxy_xxxz, g_z_x_xxy_xxyy, g_z_x_xxy_xxyz, g_z_x_xxy_xxzz, g_z_x_xxy_xyyy, g_z_x_xxy_xyyz, g_z_x_xxy_xyzz, g_z_x_xxy_xzzz, g_z_x_xxy_yyyy, g_z_x_xxy_yyyz, g_z_x_xxy_yyzz, g_z_x_xxy_yzzz, g_z_x_xxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxy_xxxx[k] = -g_z_x_xx_xxxx[k] * ab_y + g_z_x_xx_xxxxy[k];

                g_z_x_xxy_xxxy[k] = -g_z_x_xx_xxxy[k] * ab_y + g_z_x_xx_xxxyy[k];

                g_z_x_xxy_xxxz[k] = -g_z_x_xx_xxxz[k] * ab_y + g_z_x_xx_xxxyz[k];

                g_z_x_xxy_xxyy[k] = -g_z_x_xx_xxyy[k] * ab_y + g_z_x_xx_xxyyy[k];

                g_z_x_xxy_xxyz[k] = -g_z_x_xx_xxyz[k] * ab_y + g_z_x_xx_xxyyz[k];

                g_z_x_xxy_xxzz[k] = -g_z_x_xx_xxzz[k] * ab_y + g_z_x_xx_xxyzz[k];

                g_z_x_xxy_xyyy[k] = -g_z_x_xx_xyyy[k] * ab_y + g_z_x_xx_xyyyy[k];

                g_z_x_xxy_xyyz[k] = -g_z_x_xx_xyyz[k] * ab_y + g_z_x_xx_xyyyz[k];

                g_z_x_xxy_xyzz[k] = -g_z_x_xx_xyzz[k] * ab_y + g_z_x_xx_xyyzz[k];

                g_z_x_xxy_xzzz[k] = -g_z_x_xx_xzzz[k] * ab_y + g_z_x_xx_xyzzz[k];

                g_z_x_xxy_yyyy[k] = -g_z_x_xx_yyyy[k] * ab_y + g_z_x_xx_yyyyy[k];

                g_z_x_xxy_yyyz[k] = -g_z_x_xx_yyyz[k] * ab_y + g_z_x_xx_yyyyz[k];

                g_z_x_xxy_yyzz[k] = -g_z_x_xx_yyzz[k] * ab_y + g_z_x_xx_yyyzz[k];

                g_z_x_xxy_yzzz[k] = -g_z_x_xx_yzzz[k] * ab_y + g_z_x_xx_yyzzz[k];

                g_z_x_xxy_zzzz[k] = -g_z_x_xx_zzzz[k] * ab_y + g_z_x_xx_yzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxz_xxxx = cbuffer.data(fg_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_x_xxz_xxxy = cbuffer.data(fg_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_x_xxz_xxxz = cbuffer.data(fg_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_x_xxz_xxyy = cbuffer.data(fg_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_x_xxz_xxyz = cbuffer.data(fg_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_x_xxz_xxzz = cbuffer.data(fg_geom_11_off + 935 * ccomps * dcomps);

            auto g_z_x_xxz_xyyy = cbuffer.data(fg_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_x_xxz_xyyz = cbuffer.data(fg_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_x_xxz_xyzz = cbuffer.data(fg_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_x_xxz_xzzz = cbuffer.data(fg_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_x_xxz_yyyy = cbuffer.data(fg_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_x_xxz_yyyz = cbuffer.data(fg_geom_11_off + 941 * ccomps * dcomps);

            auto g_z_x_xxz_yyzz = cbuffer.data(fg_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_x_xxz_yzzz = cbuffer.data(fg_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_x_xxz_zzzz = cbuffer.data(fg_geom_11_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xxxx, g_z_0_xz_xxxy, g_z_0_xz_xxxz, g_z_0_xz_xxyy, g_z_0_xz_xxyz, g_z_0_xz_xxzz, g_z_0_xz_xyyy, g_z_0_xz_xyyz, g_z_0_xz_xyzz, g_z_0_xz_xzzz, g_z_0_xz_yyyy, g_z_0_xz_yyyz, g_z_0_xz_yyzz, g_z_0_xz_yzzz, g_z_0_xz_zzzz, g_z_x_xxz_xxxx, g_z_x_xxz_xxxy, g_z_x_xxz_xxxz, g_z_x_xxz_xxyy, g_z_x_xxz_xxyz, g_z_x_xxz_xxzz, g_z_x_xxz_xyyy, g_z_x_xxz_xyyz, g_z_x_xxz_xyzz, g_z_x_xxz_xzzz, g_z_x_xxz_yyyy, g_z_x_xxz_yyyz, g_z_x_xxz_yyzz, g_z_x_xxz_yzzz, g_z_x_xxz_zzzz, g_z_x_xz_xxxx, g_z_x_xz_xxxxx, g_z_x_xz_xxxxy, g_z_x_xz_xxxxz, g_z_x_xz_xxxy, g_z_x_xz_xxxyy, g_z_x_xz_xxxyz, g_z_x_xz_xxxz, g_z_x_xz_xxxzz, g_z_x_xz_xxyy, g_z_x_xz_xxyyy, g_z_x_xz_xxyyz, g_z_x_xz_xxyz, g_z_x_xz_xxyzz, g_z_x_xz_xxzz, g_z_x_xz_xxzzz, g_z_x_xz_xyyy, g_z_x_xz_xyyyy, g_z_x_xz_xyyyz, g_z_x_xz_xyyz, g_z_x_xz_xyyzz, g_z_x_xz_xyzz, g_z_x_xz_xyzzz, g_z_x_xz_xzzz, g_z_x_xz_xzzzz, g_z_x_xz_yyyy, g_z_x_xz_yyyz, g_z_x_xz_yyzz, g_z_x_xz_yzzz, g_z_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxz_xxxx[k] = g_z_0_xz_xxxx[k] - g_z_x_xz_xxxx[k] * ab_x + g_z_x_xz_xxxxx[k];

                g_z_x_xxz_xxxy[k] = g_z_0_xz_xxxy[k] - g_z_x_xz_xxxy[k] * ab_x + g_z_x_xz_xxxxy[k];

                g_z_x_xxz_xxxz[k] = g_z_0_xz_xxxz[k] - g_z_x_xz_xxxz[k] * ab_x + g_z_x_xz_xxxxz[k];

                g_z_x_xxz_xxyy[k] = g_z_0_xz_xxyy[k] - g_z_x_xz_xxyy[k] * ab_x + g_z_x_xz_xxxyy[k];

                g_z_x_xxz_xxyz[k] = g_z_0_xz_xxyz[k] - g_z_x_xz_xxyz[k] * ab_x + g_z_x_xz_xxxyz[k];

                g_z_x_xxz_xxzz[k] = g_z_0_xz_xxzz[k] - g_z_x_xz_xxzz[k] * ab_x + g_z_x_xz_xxxzz[k];

                g_z_x_xxz_xyyy[k] = g_z_0_xz_xyyy[k] - g_z_x_xz_xyyy[k] * ab_x + g_z_x_xz_xxyyy[k];

                g_z_x_xxz_xyyz[k] = g_z_0_xz_xyyz[k] - g_z_x_xz_xyyz[k] * ab_x + g_z_x_xz_xxyyz[k];

                g_z_x_xxz_xyzz[k] = g_z_0_xz_xyzz[k] - g_z_x_xz_xyzz[k] * ab_x + g_z_x_xz_xxyzz[k];

                g_z_x_xxz_xzzz[k] = g_z_0_xz_xzzz[k] - g_z_x_xz_xzzz[k] * ab_x + g_z_x_xz_xxzzz[k];

                g_z_x_xxz_yyyy[k] = g_z_0_xz_yyyy[k] - g_z_x_xz_yyyy[k] * ab_x + g_z_x_xz_xyyyy[k];

                g_z_x_xxz_yyyz[k] = g_z_0_xz_yyyz[k] - g_z_x_xz_yyyz[k] * ab_x + g_z_x_xz_xyyyz[k];

                g_z_x_xxz_yyzz[k] = g_z_0_xz_yyzz[k] - g_z_x_xz_yyzz[k] * ab_x + g_z_x_xz_xyyzz[k];

                g_z_x_xxz_yzzz[k] = g_z_0_xz_yzzz[k] - g_z_x_xz_yzzz[k] * ab_x + g_z_x_xz_xyzzz[k];

                g_z_x_xxz_zzzz[k] = g_z_0_xz_zzzz[k] - g_z_x_xz_zzzz[k] * ab_x + g_z_x_xz_xzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyy_xxxx = cbuffer.data(fg_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_x_xyy_xxxy = cbuffer.data(fg_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_x_xyy_xxxz = cbuffer.data(fg_geom_11_off + 947 * ccomps * dcomps);

            auto g_z_x_xyy_xxyy = cbuffer.data(fg_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_x_xyy_xxyz = cbuffer.data(fg_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_x_xyy_xxzz = cbuffer.data(fg_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_x_xyy_xyyy = cbuffer.data(fg_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_x_xyy_xyyz = cbuffer.data(fg_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_x_xyy_xyzz = cbuffer.data(fg_geom_11_off + 953 * ccomps * dcomps);

            auto g_z_x_xyy_xzzz = cbuffer.data(fg_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_x_xyy_yyyy = cbuffer.data(fg_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_x_xyy_yyyz = cbuffer.data(fg_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_x_xyy_yyzz = cbuffer.data(fg_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_x_xyy_yzzz = cbuffer.data(fg_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_x_xyy_zzzz = cbuffer.data(fg_geom_11_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xy_xxxx, g_z_x_xy_xxxxy, g_z_x_xy_xxxy, g_z_x_xy_xxxyy, g_z_x_xy_xxxyz, g_z_x_xy_xxxz, g_z_x_xy_xxyy, g_z_x_xy_xxyyy, g_z_x_xy_xxyyz, g_z_x_xy_xxyz, g_z_x_xy_xxyzz, g_z_x_xy_xxzz, g_z_x_xy_xyyy, g_z_x_xy_xyyyy, g_z_x_xy_xyyyz, g_z_x_xy_xyyz, g_z_x_xy_xyyzz, g_z_x_xy_xyzz, g_z_x_xy_xyzzz, g_z_x_xy_xzzz, g_z_x_xy_yyyy, g_z_x_xy_yyyyy, g_z_x_xy_yyyyz, g_z_x_xy_yyyz, g_z_x_xy_yyyzz, g_z_x_xy_yyzz, g_z_x_xy_yyzzz, g_z_x_xy_yzzz, g_z_x_xy_yzzzz, g_z_x_xy_zzzz, g_z_x_xyy_xxxx, g_z_x_xyy_xxxy, g_z_x_xyy_xxxz, g_z_x_xyy_xxyy, g_z_x_xyy_xxyz, g_z_x_xyy_xxzz, g_z_x_xyy_xyyy, g_z_x_xyy_xyyz, g_z_x_xyy_xyzz, g_z_x_xyy_xzzz, g_z_x_xyy_yyyy, g_z_x_xyy_yyyz, g_z_x_xyy_yyzz, g_z_x_xyy_yzzz, g_z_x_xyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyy_xxxx[k] = -g_z_x_xy_xxxx[k] * ab_y + g_z_x_xy_xxxxy[k];

                g_z_x_xyy_xxxy[k] = -g_z_x_xy_xxxy[k] * ab_y + g_z_x_xy_xxxyy[k];

                g_z_x_xyy_xxxz[k] = -g_z_x_xy_xxxz[k] * ab_y + g_z_x_xy_xxxyz[k];

                g_z_x_xyy_xxyy[k] = -g_z_x_xy_xxyy[k] * ab_y + g_z_x_xy_xxyyy[k];

                g_z_x_xyy_xxyz[k] = -g_z_x_xy_xxyz[k] * ab_y + g_z_x_xy_xxyyz[k];

                g_z_x_xyy_xxzz[k] = -g_z_x_xy_xxzz[k] * ab_y + g_z_x_xy_xxyzz[k];

                g_z_x_xyy_xyyy[k] = -g_z_x_xy_xyyy[k] * ab_y + g_z_x_xy_xyyyy[k];

                g_z_x_xyy_xyyz[k] = -g_z_x_xy_xyyz[k] * ab_y + g_z_x_xy_xyyyz[k];

                g_z_x_xyy_xyzz[k] = -g_z_x_xy_xyzz[k] * ab_y + g_z_x_xy_xyyzz[k];

                g_z_x_xyy_xzzz[k] = -g_z_x_xy_xzzz[k] * ab_y + g_z_x_xy_xyzzz[k];

                g_z_x_xyy_yyyy[k] = -g_z_x_xy_yyyy[k] * ab_y + g_z_x_xy_yyyyy[k];

                g_z_x_xyy_yyyz[k] = -g_z_x_xy_yyyz[k] * ab_y + g_z_x_xy_yyyyz[k];

                g_z_x_xyy_yyzz[k] = -g_z_x_xy_yyzz[k] * ab_y + g_z_x_xy_yyyzz[k];

                g_z_x_xyy_yzzz[k] = -g_z_x_xy_yzzz[k] * ab_y + g_z_x_xy_yyzzz[k];

                g_z_x_xyy_zzzz[k] = -g_z_x_xy_zzzz[k] * ab_y + g_z_x_xy_yzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyz_xxxx = cbuffer.data(fg_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_x_xyz_xxxy = cbuffer.data(fg_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_x_xyz_xxxz = cbuffer.data(fg_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_x_xyz_xxyy = cbuffer.data(fg_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_x_xyz_xxyz = cbuffer.data(fg_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_x_xyz_xxzz = cbuffer.data(fg_geom_11_off + 965 * ccomps * dcomps);

            auto g_z_x_xyz_xyyy = cbuffer.data(fg_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_x_xyz_xyyz = cbuffer.data(fg_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_x_xyz_xyzz = cbuffer.data(fg_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_x_xyz_xzzz = cbuffer.data(fg_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_x_xyz_yyyy = cbuffer.data(fg_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_x_xyz_yyyz = cbuffer.data(fg_geom_11_off + 971 * ccomps * dcomps);

            auto g_z_x_xyz_yyzz = cbuffer.data(fg_geom_11_off + 972 * ccomps * dcomps);

            auto g_z_x_xyz_yzzz = cbuffer.data(fg_geom_11_off + 973 * ccomps * dcomps);

            auto g_z_x_xyz_zzzz = cbuffer.data(fg_geom_11_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyz_xxxx, g_z_x_xyz_xxxy, g_z_x_xyz_xxxz, g_z_x_xyz_xxyy, g_z_x_xyz_xxyz, g_z_x_xyz_xxzz, g_z_x_xyz_xyyy, g_z_x_xyz_xyyz, g_z_x_xyz_xyzz, g_z_x_xyz_xzzz, g_z_x_xyz_yyyy, g_z_x_xyz_yyyz, g_z_x_xyz_yyzz, g_z_x_xyz_yzzz, g_z_x_xyz_zzzz, g_z_x_xz_xxxx, g_z_x_xz_xxxxy, g_z_x_xz_xxxy, g_z_x_xz_xxxyy, g_z_x_xz_xxxyz, g_z_x_xz_xxxz, g_z_x_xz_xxyy, g_z_x_xz_xxyyy, g_z_x_xz_xxyyz, g_z_x_xz_xxyz, g_z_x_xz_xxyzz, g_z_x_xz_xxzz, g_z_x_xz_xyyy, g_z_x_xz_xyyyy, g_z_x_xz_xyyyz, g_z_x_xz_xyyz, g_z_x_xz_xyyzz, g_z_x_xz_xyzz, g_z_x_xz_xyzzz, g_z_x_xz_xzzz, g_z_x_xz_yyyy, g_z_x_xz_yyyyy, g_z_x_xz_yyyyz, g_z_x_xz_yyyz, g_z_x_xz_yyyzz, g_z_x_xz_yyzz, g_z_x_xz_yyzzz, g_z_x_xz_yzzz, g_z_x_xz_yzzzz, g_z_x_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyz_xxxx[k] = -g_z_x_xz_xxxx[k] * ab_y + g_z_x_xz_xxxxy[k];

                g_z_x_xyz_xxxy[k] = -g_z_x_xz_xxxy[k] * ab_y + g_z_x_xz_xxxyy[k];

                g_z_x_xyz_xxxz[k] = -g_z_x_xz_xxxz[k] * ab_y + g_z_x_xz_xxxyz[k];

                g_z_x_xyz_xxyy[k] = -g_z_x_xz_xxyy[k] * ab_y + g_z_x_xz_xxyyy[k];

                g_z_x_xyz_xxyz[k] = -g_z_x_xz_xxyz[k] * ab_y + g_z_x_xz_xxyyz[k];

                g_z_x_xyz_xxzz[k] = -g_z_x_xz_xxzz[k] * ab_y + g_z_x_xz_xxyzz[k];

                g_z_x_xyz_xyyy[k] = -g_z_x_xz_xyyy[k] * ab_y + g_z_x_xz_xyyyy[k];

                g_z_x_xyz_xyyz[k] = -g_z_x_xz_xyyz[k] * ab_y + g_z_x_xz_xyyyz[k];

                g_z_x_xyz_xyzz[k] = -g_z_x_xz_xyzz[k] * ab_y + g_z_x_xz_xyyzz[k];

                g_z_x_xyz_xzzz[k] = -g_z_x_xz_xzzz[k] * ab_y + g_z_x_xz_xyzzz[k];

                g_z_x_xyz_yyyy[k] = -g_z_x_xz_yyyy[k] * ab_y + g_z_x_xz_yyyyy[k];

                g_z_x_xyz_yyyz[k] = -g_z_x_xz_yyyz[k] * ab_y + g_z_x_xz_yyyyz[k];

                g_z_x_xyz_yyzz[k] = -g_z_x_xz_yyzz[k] * ab_y + g_z_x_xz_yyyzz[k];

                g_z_x_xyz_yzzz[k] = -g_z_x_xz_yzzz[k] * ab_y + g_z_x_xz_yyzzz[k];

                g_z_x_xyz_zzzz[k] = -g_z_x_xz_zzzz[k] * ab_y + g_z_x_xz_yzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzz_xxxx = cbuffer.data(fg_geom_11_off + 975 * ccomps * dcomps);

            auto g_z_x_xzz_xxxy = cbuffer.data(fg_geom_11_off + 976 * ccomps * dcomps);

            auto g_z_x_xzz_xxxz = cbuffer.data(fg_geom_11_off + 977 * ccomps * dcomps);

            auto g_z_x_xzz_xxyy = cbuffer.data(fg_geom_11_off + 978 * ccomps * dcomps);

            auto g_z_x_xzz_xxyz = cbuffer.data(fg_geom_11_off + 979 * ccomps * dcomps);

            auto g_z_x_xzz_xxzz = cbuffer.data(fg_geom_11_off + 980 * ccomps * dcomps);

            auto g_z_x_xzz_xyyy = cbuffer.data(fg_geom_11_off + 981 * ccomps * dcomps);

            auto g_z_x_xzz_xyyz = cbuffer.data(fg_geom_11_off + 982 * ccomps * dcomps);

            auto g_z_x_xzz_xyzz = cbuffer.data(fg_geom_11_off + 983 * ccomps * dcomps);

            auto g_z_x_xzz_xzzz = cbuffer.data(fg_geom_11_off + 984 * ccomps * dcomps);

            auto g_z_x_xzz_yyyy = cbuffer.data(fg_geom_11_off + 985 * ccomps * dcomps);

            auto g_z_x_xzz_yyyz = cbuffer.data(fg_geom_11_off + 986 * ccomps * dcomps);

            auto g_z_x_xzz_yyzz = cbuffer.data(fg_geom_11_off + 987 * ccomps * dcomps);

            auto g_z_x_xzz_yzzz = cbuffer.data(fg_geom_11_off + 988 * ccomps * dcomps);

            auto g_z_x_xzz_zzzz = cbuffer.data(fg_geom_11_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_z_x_xzz_xxxx, g_z_x_xzz_xxxy, g_z_x_xzz_xxxz, g_z_x_xzz_xxyy, g_z_x_xzz_xxyz, g_z_x_xzz_xxzz, g_z_x_xzz_xyyy, g_z_x_xzz_xyyz, g_z_x_xzz_xyzz, g_z_x_xzz_xzzz, g_z_x_xzz_yyyy, g_z_x_xzz_yyyz, g_z_x_xzz_yyzz, g_z_x_xzz_yzzz, g_z_x_xzz_zzzz, g_z_x_zz_xxxx, g_z_x_zz_xxxxx, g_z_x_zz_xxxxy, g_z_x_zz_xxxxz, g_z_x_zz_xxxy, g_z_x_zz_xxxyy, g_z_x_zz_xxxyz, g_z_x_zz_xxxz, g_z_x_zz_xxxzz, g_z_x_zz_xxyy, g_z_x_zz_xxyyy, g_z_x_zz_xxyyz, g_z_x_zz_xxyz, g_z_x_zz_xxyzz, g_z_x_zz_xxzz, g_z_x_zz_xxzzz, g_z_x_zz_xyyy, g_z_x_zz_xyyyy, g_z_x_zz_xyyyz, g_z_x_zz_xyyz, g_z_x_zz_xyyzz, g_z_x_zz_xyzz, g_z_x_zz_xyzzz, g_z_x_zz_xzzz, g_z_x_zz_xzzzz, g_z_x_zz_yyyy, g_z_x_zz_yyyz, g_z_x_zz_yyzz, g_z_x_zz_yzzz, g_z_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzz_xxxx[k] = g_z_0_zz_xxxx[k] - g_z_x_zz_xxxx[k] * ab_x + g_z_x_zz_xxxxx[k];

                g_z_x_xzz_xxxy[k] = g_z_0_zz_xxxy[k] - g_z_x_zz_xxxy[k] * ab_x + g_z_x_zz_xxxxy[k];

                g_z_x_xzz_xxxz[k] = g_z_0_zz_xxxz[k] - g_z_x_zz_xxxz[k] * ab_x + g_z_x_zz_xxxxz[k];

                g_z_x_xzz_xxyy[k] = g_z_0_zz_xxyy[k] - g_z_x_zz_xxyy[k] * ab_x + g_z_x_zz_xxxyy[k];

                g_z_x_xzz_xxyz[k] = g_z_0_zz_xxyz[k] - g_z_x_zz_xxyz[k] * ab_x + g_z_x_zz_xxxyz[k];

                g_z_x_xzz_xxzz[k] = g_z_0_zz_xxzz[k] - g_z_x_zz_xxzz[k] * ab_x + g_z_x_zz_xxxzz[k];

                g_z_x_xzz_xyyy[k] = g_z_0_zz_xyyy[k] - g_z_x_zz_xyyy[k] * ab_x + g_z_x_zz_xxyyy[k];

                g_z_x_xzz_xyyz[k] = g_z_0_zz_xyyz[k] - g_z_x_zz_xyyz[k] * ab_x + g_z_x_zz_xxyyz[k];

                g_z_x_xzz_xyzz[k] = g_z_0_zz_xyzz[k] - g_z_x_zz_xyzz[k] * ab_x + g_z_x_zz_xxyzz[k];

                g_z_x_xzz_xzzz[k] = g_z_0_zz_xzzz[k] - g_z_x_zz_xzzz[k] * ab_x + g_z_x_zz_xxzzz[k];

                g_z_x_xzz_yyyy[k] = g_z_0_zz_yyyy[k] - g_z_x_zz_yyyy[k] * ab_x + g_z_x_zz_xyyyy[k];

                g_z_x_xzz_yyyz[k] = g_z_0_zz_yyyz[k] - g_z_x_zz_yyyz[k] * ab_x + g_z_x_zz_xyyyz[k];

                g_z_x_xzz_yyzz[k] = g_z_0_zz_yyzz[k] - g_z_x_zz_yyzz[k] * ab_x + g_z_x_zz_xyyzz[k];

                g_z_x_xzz_yzzz[k] = g_z_0_zz_yzzz[k] - g_z_x_zz_yzzz[k] * ab_x + g_z_x_zz_xyzzz[k];

                g_z_x_xzz_zzzz[k] = g_z_0_zz_zzzz[k] - g_z_x_zz_zzzz[k] * ab_x + g_z_x_zz_xzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyy_xxxx = cbuffer.data(fg_geom_11_off + 990 * ccomps * dcomps);

            auto g_z_x_yyy_xxxy = cbuffer.data(fg_geom_11_off + 991 * ccomps * dcomps);

            auto g_z_x_yyy_xxxz = cbuffer.data(fg_geom_11_off + 992 * ccomps * dcomps);

            auto g_z_x_yyy_xxyy = cbuffer.data(fg_geom_11_off + 993 * ccomps * dcomps);

            auto g_z_x_yyy_xxyz = cbuffer.data(fg_geom_11_off + 994 * ccomps * dcomps);

            auto g_z_x_yyy_xxzz = cbuffer.data(fg_geom_11_off + 995 * ccomps * dcomps);

            auto g_z_x_yyy_xyyy = cbuffer.data(fg_geom_11_off + 996 * ccomps * dcomps);

            auto g_z_x_yyy_xyyz = cbuffer.data(fg_geom_11_off + 997 * ccomps * dcomps);

            auto g_z_x_yyy_xyzz = cbuffer.data(fg_geom_11_off + 998 * ccomps * dcomps);

            auto g_z_x_yyy_xzzz = cbuffer.data(fg_geom_11_off + 999 * ccomps * dcomps);

            auto g_z_x_yyy_yyyy = cbuffer.data(fg_geom_11_off + 1000 * ccomps * dcomps);

            auto g_z_x_yyy_yyyz = cbuffer.data(fg_geom_11_off + 1001 * ccomps * dcomps);

            auto g_z_x_yyy_yyzz = cbuffer.data(fg_geom_11_off + 1002 * ccomps * dcomps);

            auto g_z_x_yyy_yzzz = cbuffer.data(fg_geom_11_off + 1003 * ccomps * dcomps);

            auto g_z_x_yyy_zzzz = cbuffer.data(fg_geom_11_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yy_xxxx, g_z_x_yy_xxxxy, g_z_x_yy_xxxy, g_z_x_yy_xxxyy, g_z_x_yy_xxxyz, g_z_x_yy_xxxz, g_z_x_yy_xxyy, g_z_x_yy_xxyyy, g_z_x_yy_xxyyz, g_z_x_yy_xxyz, g_z_x_yy_xxyzz, g_z_x_yy_xxzz, g_z_x_yy_xyyy, g_z_x_yy_xyyyy, g_z_x_yy_xyyyz, g_z_x_yy_xyyz, g_z_x_yy_xyyzz, g_z_x_yy_xyzz, g_z_x_yy_xyzzz, g_z_x_yy_xzzz, g_z_x_yy_yyyy, g_z_x_yy_yyyyy, g_z_x_yy_yyyyz, g_z_x_yy_yyyz, g_z_x_yy_yyyzz, g_z_x_yy_yyzz, g_z_x_yy_yyzzz, g_z_x_yy_yzzz, g_z_x_yy_yzzzz, g_z_x_yy_zzzz, g_z_x_yyy_xxxx, g_z_x_yyy_xxxy, g_z_x_yyy_xxxz, g_z_x_yyy_xxyy, g_z_x_yyy_xxyz, g_z_x_yyy_xxzz, g_z_x_yyy_xyyy, g_z_x_yyy_xyyz, g_z_x_yyy_xyzz, g_z_x_yyy_xzzz, g_z_x_yyy_yyyy, g_z_x_yyy_yyyz, g_z_x_yyy_yyzz, g_z_x_yyy_yzzz, g_z_x_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyy_xxxx[k] = -g_z_x_yy_xxxx[k] * ab_y + g_z_x_yy_xxxxy[k];

                g_z_x_yyy_xxxy[k] = -g_z_x_yy_xxxy[k] * ab_y + g_z_x_yy_xxxyy[k];

                g_z_x_yyy_xxxz[k] = -g_z_x_yy_xxxz[k] * ab_y + g_z_x_yy_xxxyz[k];

                g_z_x_yyy_xxyy[k] = -g_z_x_yy_xxyy[k] * ab_y + g_z_x_yy_xxyyy[k];

                g_z_x_yyy_xxyz[k] = -g_z_x_yy_xxyz[k] * ab_y + g_z_x_yy_xxyyz[k];

                g_z_x_yyy_xxzz[k] = -g_z_x_yy_xxzz[k] * ab_y + g_z_x_yy_xxyzz[k];

                g_z_x_yyy_xyyy[k] = -g_z_x_yy_xyyy[k] * ab_y + g_z_x_yy_xyyyy[k];

                g_z_x_yyy_xyyz[k] = -g_z_x_yy_xyyz[k] * ab_y + g_z_x_yy_xyyyz[k];

                g_z_x_yyy_xyzz[k] = -g_z_x_yy_xyzz[k] * ab_y + g_z_x_yy_xyyzz[k];

                g_z_x_yyy_xzzz[k] = -g_z_x_yy_xzzz[k] * ab_y + g_z_x_yy_xyzzz[k];

                g_z_x_yyy_yyyy[k] = -g_z_x_yy_yyyy[k] * ab_y + g_z_x_yy_yyyyy[k];

                g_z_x_yyy_yyyz[k] = -g_z_x_yy_yyyz[k] * ab_y + g_z_x_yy_yyyyz[k];

                g_z_x_yyy_yyzz[k] = -g_z_x_yy_yyzz[k] * ab_y + g_z_x_yy_yyyzz[k];

                g_z_x_yyy_yzzz[k] = -g_z_x_yy_yzzz[k] * ab_y + g_z_x_yy_yyzzz[k];

                g_z_x_yyy_zzzz[k] = -g_z_x_yy_zzzz[k] * ab_y + g_z_x_yy_yzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyz_xxxx = cbuffer.data(fg_geom_11_off + 1005 * ccomps * dcomps);

            auto g_z_x_yyz_xxxy = cbuffer.data(fg_geom_11_off + 1006 * ccomps * dcomps);

            auto g_z_x_yyz_xxxz = cbuffer.data(fg_geom_11_off + 1007 * ccomps * dcomps);

            auto g_z_x_yyz_xxyy = cbuffer.data(fg_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_x_yyz_xxyz = cbuffer.data(fg_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_x_yyz_xxzz = cbuffer.data(fg_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_x_yyz_xyyy = cbuffer.data(fg_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_x_yyz_xyyz = cbuffer.data(fg_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_x_yyz_xyzz = cbuffer.data(fg_geom_11_off + 1013 * ccomps * dcomps);

            auto g_z_x_yyz_xzzz = cbuffer.data(fg_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_x_yyz_yyyy = cbuffer.data(fg_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_x_yyz_yyyz = cbuffer.data(fg_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_x_yyz_yyzz = cbuffer.data(fg_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_x_yyz_yzzz = cbuffer.data(fg_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_x_yyz_zzzz = cbuffer.data(fg_geom_11_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyz_xxxx, g_z_x_yyz_xxxy, g_z_x_yyz_xxxz, g_z_x_yyz_xxyy, g_z_x_yyz_xxyz, g_z_x_yyz_xxzz, g_z_x_yyz_xyyy, g_z_x_yyz_xyyz, g_z_x_yyz_xyzz, g_z_x_yyz_xzzz, g_z_x_yyz_yyyy, g_z_x_yyz_yyyz, g_z_x_yyz_yyzz, g_z_x_yyz_yzzz, g_z_x_yyz_zzzz, g_z_x_yz_xxxx, g_z_x_yz_xxxxy, g_z_x_yz_xxxy, g_z_x_yz_xxxyy, g_z_x_yz_xxxyz, g_z_x_yz_xxxz, g_z_x_yz_xxyy, g_z_x_yz_xxyyy, g_z_x_yz_xxyyz, g_z_x_yz_xxyz, g_z_x_yz_xxyzz, g_z_x_yz_xxzz, g_z_x_yz_xyyy, g_z_x_yz_xyyyy, g_z_x_yz_xyyyz, g_z_x_yz_xyyz, g_z_x_yz_xyyzz, g_z_x_yz_xyzz, g_z_x_yz_xyzzz, g_z_x_yz_xzzz, g_z_x_yz_yyyy, g_z_x_yz_yyyyy, g_z_x_yz_yyyyz, g_z_x_yz_yyyz, g_z_x_yz_yyyzz, g_z_x_yz_yyzz, g_z_x_yz_yyzzz, g_z_x_yz_yzzz, g_z_x_yz_yzzzz, g_z_x_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyz_xxxx[k] = -g_z_x_yz_xxxx[k] * ab_y + g_z_x_yz_xxxxy[k];

                g_z_x_yyz_xxxy[k] = -g_z_x_yz_xxxy[k] * ab_y + g_z_x_yz_xxxyy[k];

                g_z_x_yyz_xxxz[k] = -g_z_x_yz_xxxz[k] * ab_y + g_z_x_yz_xxxyz[k];

                g_z_x_yyz_xxyy[k] = -g_z_x_yz_xxyy[k] * ab_y + g_z_x_yz_xxyyy[k];

                g_z_x_yyz_xxyz[k] = -g_z_x_yz_xxyz[k] * ab_y + g_z_x_yz_xxyyz[k];

                g_z_x_yyz_xxzz[k] = -g_z_x_yz_xxzz[k] * ab_y + g_z_x_yz_xxyzz[k];

                g_z_x_yyz_xyyy[k] = -g_z_x_yz_xyyy[k] * ab_y + g_z_x_yz_xyyyy[k];

                g_z_x_yyz_xyyz[k] = -g_z_x_yz_xyyz[k] * ab_y + g_z_x_yz_xyyyz[k];

                g_z_x_yyz_xyzz[k] = -g_z_x_yz_xyzz[k] * ab_y + g_z_x_yz_xyyzz[k];

                g_z_x_yyz_xzzz[k] = -g_z_x_yz_xzzz[k] * ab_y + g_z_x_yz_xyzzz[k];

                g_z_x_yyz_yyyy[k] = -g_z_x_yz_yyyy[k] * ab_y + g_z_x_yz_yyyyy[k];

                g_z_x_yyz_yyyz[k] = -g_z_x_yz_yyyz[k] * ab_y + g_z_x_yz_yyyyz[k];

                g_z_x_yyz_yyzz[k] = -g_z_x_yz_yyzz[k] * ab_y + g_z_x_yz_yyyzz[k];

                g_z_x_yyz_yzzz[k] = -g_z_x_yz_yzzz[k] * ab_y + g_z_x_yz_yyzzz[k];

                g_z_x_yyz_zzzz[k] = -g_z_x_yz_zzzz[k] * ab_y + g_z_x_yz_yzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzz_xxxx = cbuffer.data(fg_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_x_yzz_xxxy = cbuffer.data(fg_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_x_yzz_xxxz = cbuffer.data(fg_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_x_yzz_xxyy = cbuffer.data(fg_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_x_yzz_xxyz = cbuffer.data(fg_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_x_yzz_xxzz = cbuffer.data(fg_geom_11_off + 1025 * ccomps * dcomps);

            auto g_z_x_yzz_xyyy = cbuffer.data(fg_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_x_yzz_xyyz = cbuffer.data(fg_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_x_yzz_xyzz = cbuffer.data(fg_geom_11_off + 1028 * ccomps * dcomps);

            auto g_z_x_yzz_xzzz = cbuffer.data(fg_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_x_yzz_yyyy = cbuffer.data(fg_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_x_yzz_yyyz = cbuffer.data(fg_geom_11_off + 1031 * ccomps * dcomps);

            auto g_z_x_yzz_yyzz = cbuffer.data(fg_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_x_yzz_yzzz = cbuffer.data(fg_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_x_yzz_zzzz = cbuffer.data(fg_geom_11_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzz_xxxx, g_z_x_yzz_xxxy, g_z_x_yzz_xxxz, g_z_x_yzz_xxyy, g_z_x_yzz_xxyz, g_z_x_yzz_xxzz, g_z_x_yzz_xyyy, g_z_x_yzz_xyyz, g_z_x_yzz_xyzz, g_z_x_yzz_xzzz, g_z_x_yzz_yyyy, g_z_x_yzz_yyyz, g_z_x_yzz_yyzz, g_z_x_yzz_yzzz, g_z_x_yzz_zzzz, g_z_x_zz_xxxx, g_z_x_zz_xxxxy, g_z_x_zz_xxxy, g_z_x_zz_xxxyy, g_z_x_zz_xxxyz, g_z_x_zz_xxxz, g_z_x_zz_xxyy, g_z_x_zz_xxyyy, g_z_x_zz_xxyyz, g_z_x_zz_xxyz, g_z_x_zz_xxyzz, g_z_x_zz_xxzz, g_z_x_zz_xyyy, g_z_x_zz_xyyyy, g_z_x_zz_xyyyz, g_z_x_zz_xyyz, g_z_x_zz_xyyzz, g_z_x_zz_xyzz, g_z_x_zz_xyzzz, g_z_x_zz_xzzz, g_z_x_zz_yyyy, g_z_x_zz_yyyyy, g_z_x_zz_yyyyz, g_z_x_zz_yyyz, g_z_x_zz_yyyzz, g_z_x_zz_yyzz, g_z_x_zz_yyzzz, g_z_x_zz_yzzz, g_z_x_zz_yzzzz, g_z_x_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzz_xxxx[k] = -g_z_x_zz_xxxx[k] * ab_y + g_z_x_zz_xxxxy[k];

                g_z_x_yzz_xxxy[k] = -g_z_x_zz_xxxy[k] * ab_y + g_z_x_zz_xxxyy[k];

                g_z_x_yzz_xxxz[k] = -g_z_x_zz_xxxz[k] * ab_y + g_z_x_zz_xxxyz[k];

                g_z_x_yzz_xxyy[k] = -g_z_x_zz_xxyy[k] * ab_y + g_z_x_zz_xxyyy[k];

                g_z_x_yzz_xxyz[k] = -g_z_x_zz_xxyz[k] * ab_y + g_z_x_zz_xxyyz[k];

                g_z_x_yzz_xxzz[k] = -g_z_x_zz_xxzz[k] * ab_y + g_z_x_zz_xxyzz[k];

                g_z_x_yzz_xyyy[k] = -g_z_x_zz_xyyy[k] * ab_y + g_z_x_zz_xyyyy[k];

                g_z_x_yzz_xyyz[k] = -g_z_x_zz_xyyz[k] * ab_y + g_z_x_zz_xyyyz[k];

                g_z_x_yzz_xyzz[k] = -g_z_x_zz_xyzz[k] * ab_y + g_z_x_zz_xyyzz[k];

                g_z_x_yzz_xzzz[k] = -g_z_x_zz_xzzz[k] * ab_y + g_z_x_zz_xyzzz[k];

                g_z_x_yzz_yyyy[k] = -g_z_x_zz_yyyy[k] * ab_y + g_z_x_zz_yyyyy[k];

                g_z_x_yzz_yyyz[k] = -g_z_x_zz_yyyz[k] * ab_y + g_z_x_zz_yyyyz[k];

                g_z_x_yzz_yyzz[k] = -g_z_x_zz_yyzz[k] * ab_y + g_z_x_zz_yyyzz[k];

                g_z_x_yzz_yzzz[k] = -g_z_x_zz_yzzz[k] * ab_y + g_z_x_zz_yyzzz[k];

                g_z_x_yzz_zzzz[k] = -g_z_x_zz_zzzz[k] * ab_y + g_z_x_zz_yzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzz_xxxx = cbuffer.data(fg_geom_11_off + 1035 * ccomps * dcomps);

            auto g_z_x_zzz_xxxy = cbuffer.data(fg_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_x_zzz_xxxz = cbuffer.data(fg_geom_11_off + 1037 * ccomps * dcomps);

            auto g_z_x_zzz_xxyy = cbuffer.data(fg_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_x_zzz_xxyz = cbuffer.data(fg_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_x_zzz_xxzz = cbuffer.data(fg_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_x_zzz_xyyy = cbuffer.data(fg_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_x_zzz_xyyz = cbuffer.data(fg_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_x_zzz_xyzz = cbuffer.data(fg_geom_11_off + 1043 * ccomps * dcomps);

            auto g_z_x_zzz_xzzz = cbuffer.data(fg_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_x_zzz_yyyy = cbuffer.data(fg_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_x_zzz_yyyz = cbuffer.data(fg_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_x_zzz_yyzz = cbuffer.data(fg_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_x_zzz_yzzz = cbuffer.data(fg_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_x_zzz_zzzz = cbuffer.data(fg_geom_11_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_xxxx, g_0_x_zz_xxxy, g_0_x_zz_xxxz, g_0_x_zz_xxyy, g_0_x_zz_xxyz, g_0_x_zz_xxzz, g_0_x_zz_xyyy, g_0_x_zz_xyyz, g_0_x_zz_xyzz, g_0_x_zz_xzzz, g_0_x_zz_yyyy, g_0_x_zz_yyyz, g_0_x_zz_yyzz, g_0_x_zz_yzzz, g_0_x_zz_zzzz, g_z_x_zz_xxxx, g_z_x_zz_xxxxz, g_z_x_zz_xxxy, g_z_x_zz_xxxyz, g_z_x_zz_xxxz, g_z_x_zz_xxxzz, g_z_x_zz_xxyy, g_z_x_zz_xxyyz, g_z_x_zz_xxyz, g_z_x_zz_xxyzz, g_z_x_zz_xxzz, g_z_x_zz_xxzzz, g_z_x_zz_xyyy, g_z_x_zz_xyyyz, g_z_x_zz_xyyz, g_z_x_zz_xyyzz, g_z_x_zz_xyzz, g_z_x_zz_xyzzz, g_z_x_zz_xzzz, g_z_x_zz_xzzzz, g_z_x_zz_yyyy, g_z_x_zz_yyyyz, g_z_x_zz_yyyz, g_z_x_zz_yyyzz, g_z_x_zz_yyzz, g_z_x_zz_yyzzz, g_z_x_zz_yzzz, g_z_x_zz_yzzzz, g_z_x_zz_zzzz, g_z_x_zz_zzzzz, g_z_x_zzz_xxxx, g_z_x_zzz_xxxy, g_z_x_zzz_xxxz, g_z_x_zzz_xxyy, g_z_x_zzz_xxyz, g_z_x_zzz_xxzz, g_z_x_zzz_xyyy, g_z_x_zzz_xyyz, g_z_x_zzz_xyzz, g_z_x_zzz_xzzz, g_z_x_zzz_yyyy, g_z_x_zzz_yyyz, g_z_x_zzz_yyzz, g_z_x_zzz_yzzz, g_z_x_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzz_xxxx[k] = -g_0_x_zz_xxxx[k] - g_z_x_zz_xxxx[k] * ab_z + g_z_x_zz_xxxxz[k];

                g_z_x_zzz_xxxy[k] = -g_0_x_zz_xxxy[k] - g_z_x_zz_xxxy[k] * ab_z + g_z_x_zz_xxxyz[k];

                g_z_x_zzz_xxxz[k] = -g_0_x_zz_xxxz[k] - g_z_x_zz_xxxz[k] * ab_z + g_z_x_zz_xxxzz[k];

                g_z_x_zzz_xxyy[k] = -g_0_x_zz_xxyy[k] - g_z_x_zz_xxyy[k] * ab_z + g_z_x_zz_xxyyz[k];

                g_z_x_zzz_xxyz[k] = -g_0_x_zz_xxyz[k] - g_z_x_zz_xxyz[k] * ab_z + g_z_x_zz_xxyzz[k];

                g_z_x_zzz_xxzz[k] = -g_0_x_zz_xxzz[k] - g_z_x_zz_xxzz[k] * ab_z + g_z_x_zz_xxzzz[k];

                g_z_x_zzz_xyyy[k] = -g_0_x_zz_xyyy[k] - g_z_x_zz_xyyy[k] * ab_z + g_z_x_zz_xyyyz[k];

                g_z_x_zzz_xyyz[k] = -g_0_x_zz_xyyz[k] - g_z_x_zz_xyyz[k] * ab_z + g_z_x_zz_xyyzz[k];

                g_z_x_zzz_xyzz[k] = -g_0_x_zz_xyzz[k] - g_z_x_zz_xyzz[k] * ab_z + g_z_x_zz_xyzzz[k];

                g_z_x_zzz_xzzz[k] = -g_0_x_zz_xzzz[k] - g_z_x_zz_xzzz[k] * ab_z + g_z_x_zz_xzzzz[k];

                g_z_x_zzz_yyyy[k] = -g_0_x_zz_yyyy[k] - g_z_x_zz_yyyy[k] * ab_z + g_z_x_zz_yyyyz[k];

                g_z_x_zzz_yyyz[k] = -g_0_x_zz_yyyz[k] - g_z_x_zz_yyyz[k] * ab_z + g_z_x_zz_yyyzz[k];

                g_z_x_zzz_yyzz[k] = -g_0_x_zz_yyzz[k] - g_z_x_zz_yyzz[k] * ab_z + g_z_x_zz_yyzzz[k];

                g_z_x_zzz_yzzz[k] = -g_0_x_zz_yzzz[k] - g_z_x_zz_yzzz[k] * ab_z + g_z_x_zz_yzzzz[k];

                g_z_x_zzz_zzzz[k] = -g_0_x_zz_zzzz[k] - g_z_x_zz_zzzz[k] * ab_z + g_z_x_zz_zzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxx_xxxx = cbuffer.data(fg_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_y_xxx_xxxy = cbuffer.data(fg_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_y_xxx_xxxz = cbuffer.data(fg_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_y_xxx_xxyy = cbuffer.data(fg_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_y_xxx_xxyz = cbuffer.data(fg_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_y_xxx_xxzz = cbuffer.data(fg_geom_11_off + 1055 * ccomps * dcomps);

            auto g_z_y_xxx_xyyy = cbuffer.data(fg_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_y_xxx_xyyz = cbuffer.data(fg_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_y_xxx_xyzz = cbuffer.data(fg_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_y_xxx_xzzz = cbuffer.data(fg_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_y_xxx_yyyy = cbuffer.data(fg_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_y_xxx_yyyz = cbuffer.data(fg_geom_11_off + 1061 * ccomps * dcomps);

            auto g_z_y_xxx_yyzz = cbuffer.data(fg_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_y_xxx_yzzz = cbuffer.data(fg_geom_11_off + 1063 * ccomps * dcomps);

            auto g_z_y_xxx_zzzz = cbuffer.data(fg_geom_11_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xx_xxxx, g_z_y_xx_xxxxx, g_z_y_xx_xxxxy, g_z_y_xx_xxxxz, g_z_y_xx_xxxy, g_z_y_xx_xxxyy, g_z_y_xx_xxxyz, g_z_y_xx_xxxz, g_z_y_xx_xxxzz, g_z_y_xx_xxyy, g_z_y_xx_xxyyy, g_z_y_xx_xxyyz, g_z_y_xx_xxyz, g_z_y_xx_xxyzz, g_z_y_xx_xxzz, g_z_y_xx_xxzzz, g_z_y_xx_xyyy, g_z_y_xx_xyyyy, g_z_y_xx_xyyyz, g_z_y_xx_xyyz, g_z_y_xx_xyyzz, g_z_y_xx_xyzz, g_z_y_xx_xyzzz, g_z_y_xx_xzzz, g_z_y_xx_xzzzz, g_z_y_xx_yyyy, g_z_y_xx_yyyz, g_z_y_xx_yyzz, g_z_y_xx_yzzz, g_z_y_xx_zzzz, g_z_y_xxx_xxxx, g_z_y_xxx_xxxy, g_z_y_xxx_xxxz, g_z_y_xxx_xxyy, g_z_y_xxx_xxyz, g_z_y_xxx_xxzz, g_z_y_xxx_xyyy, g_z_y_xxx_xyyz, g_z_y_xxx_xyzz, g_z_y_xxx_xzzz, g_z_y_xxx_yyyy, g_z_y_xxx_yyyz, g_z_y_xxx_yyzz, g_z_y_xxx_yzzz, g_z_y_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxx_xxxx[k] = -g_z_y_xx_xxxx[k] * ab_x + g_z_y_xx_xxxxx[k];

                g_z_y_xxx_xxxy[k] = -g_z_y_xx_xxxy[k] * ab_x + g_z_y_xx_xxxxy[k];

                g_z_y_xxx_xxxz[k] = -g_z_y_xx_xxxz[k] * ab_x + g_z_y_xx_xxxxz[k];

                g_z_y_xxx_xxyy[k] = -g_z_y_xx_xxyy[k] * ab_x + g_z_y_xx_xxxyy[k];

                g_z_y_xxx_xxyz[k] = -g_z_y_xx_xxyz[k] * ab_x + g_z_y_xx_xxxyz[k];

                g_z_y_xxx_xxzz[k] = -g_z_y_xx_xxzz[k] * ab_x + g_z_y_xx_xxxzz[k];

                g_z_y_xxx_xyyy[k] = -g_z_y_xx_xyyy[k] * ab_x + g_z_y_xx_xxyyy[k];

                g_z_y_xxx_xyyz[k] = -g_z_y_xx_xyyz[k] * ab_x + g_z_y_xx_xxyyz[k];

                g_z_y_xxx_xyzz[k] = -g_z_y_xx_xyzz[k] * ab_x + g_z_y_xx_xxyzz[k];

                g_z_y_xxx_xzzz[k] = -g_z_y_xx_xzzz[k] * ab_x + g_z_y_xx_xxzzz[k];

                g_z_y_xxx_yyyy[k] = -g_z_y_xx_yyyy[k] * ab_x + g_z_y_xx_xyyyy[k];

                g_z_y_xxx_yyyz[k] = -g_z_y_xx_yyyz[k] * ab_x + g_z_y_xx_xyyyz[k];

                g_z_y_xxx_yyzz[k] = -g_z_y_xx_yyzz[k] * ab_x + g_z_y_xx_xyyzz[k];

                g_z_y_xxx_yzzz[k] = -g_z_y_xx_yzzz[k] * ab_x + g_z_y_xx_xyzzz[k];

                g_z_y_xxx_zzzz[k] = -g_z_y_xx_zzzz[k] * ab_x + g_z_y_xx_xzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxy_xxxx = cbuffer.data(fg_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_y_xxy_xxxy = cbuffer.data(fg_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_y_xxy_xxxz = cbuffer.data(fg_geom_11_off + 1067 * ccomps * dcomps);

            auto g_z_y_xxy_xxyy = cbuffer.data(fg_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_y_xxy_xxyz = cbuffer.data(fg_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_y_xxy_xxzz = cbuffer.data(fg_geom_11_off + 1070 * ccomps * dcomps);

            auto g_z_y_xxy_xyyy = cbuffer.data(fg_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_y_xxy_xyyz = cbuffer.data(fg_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_y_xxy_xyzz = cbuffer.data(fg_geom_11_off + 1073 * ccomps * dcomps);

            auto g_z_y_xxy_xzzz = cbuffer.data(fg_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_y_xxy_yyyy = cbuffer.data(fg_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_y_xxy_yyyz = cbuffer.data(fg_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_y_xxy_yyzz = cbuffer.data(fg_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_y_xxy_yzzz = cbuffer.data(fg_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_y_xxy_zzzz = cbuffer.data(fg_geom_11_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxy_xxxx, g_z_y_xxy_xxxy, g_z_y_xxy_xxxz, g_z_y_xxy_xxyy, g_z_y_xxy_xxyz, g_z_y_xxy_xxzz, g_z_y_xxy_xyyy, g_z_y_xxy_xyyz, g_z_y_xxy_xyzz, g_z_y_xxy_xzzz, g_z_y_xxy_yyyy, g_z_y_xxy_yyyz, g_z_y_xxy_yyzz, g_z_y_xxy_yzzz, g_z_y_xxy_zzzz, g_z_y_xy_xxxx, g_z_y_xy_xxxxx, g_z_y_xy_xxxxy, g_z_y_xy_xxxxz, g_z_y_xy_xxxy, g_z_y_xy_xxxyy, g_z_y_xy_xxxyz, g_z_y_xy_xxxz, g_z_y_xy_xxxzz, g_z_y_xy_xxyy, g_z_y_xy_xxyyy, g_z_y_xy_xxyyz, g_z_y_xy_xxyz, g_z_y_xy_xxyzz, g_z_y_xy_xxzz, g_z_y_xy_xxzzz, g_z_y_xy_xyyy, g_z_y_xy_xyyyy, g_z_y_xy_xyyyz, g_z_y_xy_xyyz, g_z_y_xy_xyyzz, g_z_y_xy_xyzz, g_z_y_xy_xyzzz, g_z_y_xy_xzzz, g_z_y_xy_xzzzz, g_z_y_xy_yyyy, g_z_y_xy_yyyz, g_z_y_xy_yyzz, g_z_y_xy_yzzz, g_z_y_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxy_xxxx[k] = -g_z_y_xy_xxxx[k] * ab_x + g_z_y_xy_xxxxx[k];

                g_z_y_xxy_xxxy[k] = -g_z_y_xy_xxxy[k] * ab_x + g_z_y_xy_xxxxy[k];

                g_z_y_xxy_xxxz[k] = -g_z_y_xy_xxxz[k] * ab_x + g_z_y_xy_xxxxz[k];

                g_z_y_xxy_xxyy[k] = -g_z_y_xy_xxyy[k] * ab_x + g_z_y_xy_xxxyy[k];

                g_z_y_xxy_xxyz[k] = -g_z_y_xy_xxyz[k] * ab_x + g_z_y_xy_xxxyz[k];

                g_z_y_xxy_xxzz[k] = -g_z_y_xy_xxzz[k] * ab_x + g_z_y_xy_xxxzz[k];

                g_z_y_xxy_xyyy[k] = -g_z_y_xy_xyyy[k] * ab_x + g_z_y_xy_xxyyy[k];

                g_z_y_xxy_xyyz[k] = -g_z_y_xy_xyyz[k] * ab_x + g_z_y_xy_xxyyz[k];

                g_z_y_xxy_xyzz[k] = -g_z_y_xy_xyzz[k] * ab_x + g_z_y_xy_xxyzz[k];

                g_z_y_xxy_xzzz[k] = -g_z_y_xy_xzzz[k] * ab_x + g_z_y_xy_xxzzz[k];

                g_z_y_xxy_yyyy[k] = -g_z_y_xy_yyyy[k] * ab_x + g_z_y_xy_xyyyy[k];

                g_z_y_xxy_yyyz[k] = -g_z_y_xy_yyyz[k] * ab_x + g_z_y_xy_xyyyz[k];

                g_z_y_xxy_yyzz[k] = -g_z_y_xy_yyzz[k] * ab_x + g_z_y_xy_xyyzz[k];

                g_z_y_xxy_yzzz[k] = -g_z_y_xy_yzzz[k] * ab_x + g_z_y_xy_xyzzz[k];

                g_z_y_xxy_zzzz[k] = -g_z_y_xy_zzzz[k] * ab_x + g_z_y_xy_xzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxz_xxxx = cbuffer.data(fg_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_y_xxz_xxxy = cbuffer.data(fg_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_y_xxz_xxxz = cbuffer.data(fg_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_y_xxz_xxyy = cbuffer.data(fg_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_y_xxz_xxyz = cbuffer.data(fg_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_y_xxz_xxzz = cbuffer.data(fg_geom_11_off + 1085 * ccomps * dcomps);

            auto g_z_y_xxz_xyyy = cbuffer.data(fg_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_y_xxz_xyyz = cbuffer.data(fg_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_y_xxz_xyzz = cbuffer.data(fg_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_y_xxz_xzzz = cbuffer.data(fg_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_y_xxz_yyyy = cbuffer.data(fg_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_y_xxz_yyyz = cbuffer.data(fg_geom_11_off + 1091 * ccomps * dcomps);

            auto g_z_y_xxz_yyzz = cbuffer.data(fg_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_y_xxz_yzzz = cbuffer.data(fg_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_y_xxz_zzzz = cbuffer.data(fg_geom_11_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxz_xxxx, g_z_y_xxz_xxxy, g_z_y_xxz_xxxz, g_z_y_xxz_xxyy, g_z_y_xxz_xxyz, g_z_y_xxz_xxzz, g_z_y_xxz_xyyy, g_z_y_xxz_xyyz, g_z_y_xxz_xyzz, g_z_y_xxz_xzzz, g_z_y_xxz_yyyy, g_z_y_xxz_yyyz, g_z_y_xxz_yyzz, g_z_y_xxz_yzzz, g_z_y_xxz_zzzz, g_z_y_xz_xxxx, g_z_y_xz_xxxxx, g_z_y_xz_xxxxy, g_z_y_xz_xxxxz, g_z_y_xz_xxxy, g_z_y_xz_xxxyy, g_z_y_xz_xxxyz, g_z_y_xz_xxxz, g_z_y_xz_xxxzz, g_z_y_xz_xxyy, g_z_y_xz_xxyyy, g_z_y_xz_xxyyz, g_z_y_xz_xxyz, g_z_y_xz_xxyzz, g_z_y_xz_xxzz, g_z_y_xz_xxzzz, g_z_y_xz_xyyy, g_z_y_xz_xyyyy, g_z_y_xz_xyyyz, g_z_y_xz_xyyz, g_z_y_xz_xyyzz, g_z_y_xz_xyzz, g_z_y_xz_xyzzz, g_z_y_xz_xzzz, g_z_y_xz_xzzzz, g_z_y_xz_yyyy, g_z_y_xz_yyyz, g_z_y_xz_yyzz, g_z_y_xz_yzzz, g_z_y_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxz_xxxx[k] = -g_z_y_xz_xxxx[k] * ab_x + g_z_y_xz_xxxxx[k];

                g_z_y_xxz_xxxy[k] = -g_z_y_xz_xxxy[k] * ab_x + g_z_y_xz_xxxxy[k];

                g_z_y_xxz_xxxz[k] = -g_z_y_xz_xxxz[k] * ab_x + g_z_y_xz_xxxxz[k];

                g_z_y_xxz_xxyy[k] = -g_z_y_xz_xxyy[k] * ab_x + g_z_y_xz_xxxyy[k];

                g_z_y_xxz_xxyz[k] = -g_z_y_xz_xxyz[k] * ab_x + g_z_y_xz_xxxyz[k];

                g_z_y_xxz_xxzz[k] = -g_z_y_xz_xxzz[k] * ab_x + g_z_y_xz_xxxzz[k];

                g_z_y_xxz_xyyy[k] = -g_z_y_xz_xyyy[k] * ab_x + g_z_y_xz_xxyyy[k];

                g_z_y_xxz_xyyz[k] = -g_z_y_xz_xyyz[k] * ab_x + g_z_y_xz_xxyyz[k];

                g_z_y_xxz_xyzz[k] = -g_z_y_xz_xyzz[k] * ab_x + g_z_y_xz_xxyzz[k];

                g_z_y_xxz_xzzz[k] = -g_z_y_xz_xzzz[k] * ab_x + g_z_y_xz_xxzzz[k];

                g_z_y_xxz_yyyy[k] = -g_z_y_xz_yyyy[k] * ab_x + g_z_y_xz_xyyyy[k];

                g_z_y_xxz_yyyz[k] = -g_z_y_xz_yyyz[k] * ab_x + g_z_y_xz_xyyyz[k];

                g_z_y_xxz_yyzz[k] = -g_z_y_xz_yyzz[k] * ab_x + g_z_y_xz_xyyzz[k];

                g_z_y_xxz_yzzz[k] = -g_z_y_xz_yzzz[k] * ab_x + g_z_y_xz_xyzzz[k];

                g_z_y_xxz_zzzz[k] = -g_z_y_xz_zzzz[k] * ab_x + g_z_y_xz_xzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyy_xxxx = cbuffer.data(fg_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_y_xyy_xxxy = cbuffer.data(fg_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_y_xyy_xxxz = cbuffer.data(fg_geom_11_off + 1097 * ccomps * dcomps);

            auto g_z_y_xyy_xxyy = cbuffer.data(fg_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_y_xyy_xxyz = cbuffer.data(fg_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_y_xyy_xxzz = cbuffer.data(fg_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_y_xyy_xyyy = cbuffer.data(fg_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_y_xyy_xyyz = cbuffer.data(fg_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_y_xyy_xyzz = cbuffer.data(fg_geom_11_off + 1103 * ccomps * dcomps);

            auto g_z_y_xyy_xzzz = cbuffer.data(fg_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_y_xyy_yyyy = cbuffer.data(fg_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_y_xyy_yyyz = cbuffer.data(fg_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_y_xyy_yyzz = cbuffer.data(fg_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_y_xyy_yzzz = cbuffer.data(fg_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_y_xyy_zzzz = cbuffer.data(fg_geom_11_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyy_xxxx, g_z_y_xyy_xxxy, g_z_y_xyy_xxxz, g_z_y_xyy_xxyy, g_z_y_xyy_xxyz, g_z_y_xyy_xxzz, g_z_y_xyy_xyyy, g_z_y_xyy_xyyz, g_z_y_xyy_xyzz, g_z_y_xyy_xzzz, g_z_y_xyy_yyyy, g_z_y_xyy_yyyz, g_z_y_xyy_yyzz, g_z_y_xyy_yzzz, g_z_y_xyy_zzzz, g_z_y_yy_xxxx, g_z_y_yy_xxxxx, g_z_y_yy_xxxxy, g_z_y_yy_xxxxz, g_z_y_yy_xxxy, g_z_y_yy_xxxyy, g_z_y_yy_xxxyz, g_z_y_yy_xxxz, g_z_y_yy_xxxzz, g_z_y_yy_xxyy, g_z_y_yy_xxyyy, g_z_y_yy_xxyyz, g_z_y_yy_xxyz, g_z_y_yy_xxyzz, g_z_y_yy_xxzz, g_z_y_yy_xxzzz, g_z_y_yy_xyyy, g_z_y_yy_xyyyy, g_z_y_yy_xyyyz, g_z_y_yy_xyyz, g_z_y_yy_xyyzz, g_z_y_yy_xyzz, g_z_y_yy_xyzzz, g_z_y_yy_xzzz, g_z_y_yy_xzzzz, g_z_y_yy_yyyy, g_z_y_yy_yyyz, g_z_y_yy_yyzz, g_z_y_yy_yzzz, g_z_y_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyy_xxxx[k] = -g_z_y_yy_xxxx[k] * ab_x + g_z_y_yy_xxxxx[k];

                g_z_y_xyy_xxxy[k] = -g_z_y_yy_xxxy[k] * ab_x + g_z_y_yy_xxxxy[k];

                g_z_y_xyy_xxxz[k] = -g_z_y_yy_xxxz[k] * ab_x + g_z_y_yy_xxxxz[k];

                g_z_y_xyy_xxyy[k] = -g_z_y_yy_xxyy[k] * ab_x + g_z_y_yy_xxxyy[k];

                g_z_y_xyy_xxyz[k] = -g_z_y_yy_xxyz[k] * ab_x + g_z_y_yy_xxxyz[k];

                g_z_y_xyy_xxzz[k] = -g_z_y_yy_xxzz[k] * ab_x + g_z_y_yy_xxxzz[k];

                g_z_y_xyy_xyyy[k] = -g_z_y_yy_xyyy[k] * ab_x + g_z_y_yy_xxyyy[k];

                g_z_y_xyy_xyyz[k] = -g_z_y_yy_xyyz[k] * ab_x + g_z_y_yy_xxyyz[k];

                g_z_y_xyy_xyzz[k] = -g_z_y_yy_xyzz[k] * ab_x + g_z_y_yy_xxyzz[k];

                g_z_y_xyy_xzzz[k] = -g_z_y_yy_xzzz[k] * ab_x + g_z_y_yy_xxzzz[k];

                g_z_y_xyy_yyyy[k] = -g_z_y_yy_yyyy[k] * ab_x + g_z_y_yy_xyyyy[k];

                g_z_y_xyy_yyyz[k] = -g_z_y_yy_yyyz[k] * ab_x + g_z_y_yy_xyyyz[k];

                g_z_y_xyy_yyzz[k] = -g_z_y_yy_yyzz[k] * ab_x + g_z_y_yy_xyyzz[k];

                g_z_y_xyy_yzzz[k] = -g_z_y_yy_yzzz[k] * ab_x + g_z_y_yy_xyzzz[k];

                g_z_y_xyy_zzzz[k] = -g_z_y_yy_zzzz[k] * ab_x + g_z_y_yy_xzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyz_xxxx = cbuffer.data(fg_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_y_xyz_xxxy = cbuffer.data(fg_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_y_xyz_xxxz = cbuffer.data(fg_geom_11_off + 1112 * ccomps * dcomps);

            auto g_z_y_xyz_xxyy = cbuffer.data(fg_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_y_xyz_xxyz = cbuffer.data(fg_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_y_xyz_xxzz = cbuffer.data(fg_geom_11_off + 1115 * ccomps * dcomps);

            auto g_z_y_xyz_xyyy = cbuffer.data(fg_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_y_xyz_xyyz = cbuffer.data(fg_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_y_xyz_xyzz = cbuffer.data(fg_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_y_xyz_xzzz = cbuffer.data(fg_geom_11_off + 1119 * ccomps * dcomps);

            auto g_z_y_xyz_yyyy = cbuffer.data(fg_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_y_xyz_yyyz = cbuffer.data(fg_geom_11_off + 1121 * ccomps * dcomps);

            auto g_z_y_xyz_yyzz = cbuffer.data(fg_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_y_xyz_yzzz = cbuffer.data(fg_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_y_xyz_zzzz = cbuffer.data(fg_geom_11_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyz_xxxx, g_z_y_xyz_xxxy, g_z_y_xyz_xxxz, g_z_y_xyz_xxyy, g_z_y_xyz_xxyz, g_z_y_xyz_xxzz, g_z_y_xyz_xyyy, g_z_y_xyz_xyyz, g_z_y_xyz_xyzz, g_z_y_xyz_xzzz, g_z_y_xyz_yyyy, g_z_y_xyz_yyyz, g_z_y_xyz_yyzz, g_z_y_xyz_yzzz, g_z_y_xyz_zzzz, g_z_y_yz_xxxx, g_z_y_yz_xxxxx, g_z_y_yz_xxxxy, g_z_y_yz_xxxxz, g_z_y_yz_xxxy, g_z_y_yz_xxxyy, g_z_y_yz_xxxyz, g_z_y_yz_xxxz, g_z_y_yz_xxxzz, g_z_y_yz_xxyy, g_z_y_yz_xxyyy, g_z_y_yz_xxyyz, g_z_y_yz_xxyz, g_z_y_yz_xxyzz, g_z_y_yz_xxzz, g_z_y_yz_xxzzz, g_z_y_yz_xyyy, g_z_y_yz_xyyyy, g_z_y_yz_xyyyz, g_z_y_yz_xyyz, g_z_y_yz_xyyzz, g_z_y_yz_xyzz, g_z_y_yz_xyzzz, g_z_y_yz_xzzz, g_z_y_yz_xzzzz, g_z_y_yz_yyyy, g_z_y_yz_yyyz, g_z_y_yz_yyzz, g_z_y_yz_yzzz, g_z_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyz_xxxx[k] = -g_z_y_yz_xxxx[k] * ab_x + g_z_y_yz_xxxxx[k];

                g_z_y_xyz_xxxy[k] = -g_z_y_yz_xxxy[k] * ab_x + g_z_y_yz_xxxxy[k];

                g_z_y_xyz_xxxz[k] = -g_z_y_yz_xxxz[k] * ab_x + g_z_y_yz_xxxxz[k];

                g_z_y_xyz_xxyy[k] = -g_z_y_yz_xxyy[k] * ab_x + g_z_y_yz_xxxyy[k];

                g_z_y_xyz_xxyz[k] = -g_z_y_yz_xxyz[k] * ab_x + g_z_y_yz_xxxyz[k];

                g_z_y_xyz_xxzz[k] = -g_z_y_yz_xxzz[k] * ab_x + g_z_y_yz_xxxzz[k];

                g_z_y_xyz_xyyy[k] = -g_z_y_yz_xyyy[k] * ab_x + g_z_y_yz_xxyyy[k];

                g_z_y_xyz_xyyz[k] = -g_z_y_yz_xyyz[k] * ab_x + g_z_y_yz_xxyyz[k];

                g_z_y_xyz_xyzz[k] = -g_z_y_yz_xyzz[k] * ab_x + g_z_y_yz_xxyzz[k];

                g_z_y_xyz_xzzz[k] = -g_z_y_yz_xzzz[k] * ab_x + g_z_y_yz_xxzzz[k];

                g_z_y_xyz_yyyy[k] = -g_z_y_yz_yyyy[k] * ab_x + g_z_y_yz_xyyyy[k];

                g_z_y_xyz_yyyz[k] = -g_z_y_yz_yyyz[k] * ab_x + g_z_y_yz_xyyyz[k];

                g_z_y_xyz_yyzz[k] = -g_z_y_yz_yyzz[k] * ab_x + g_z_y_yz_xyyzz[k];

                g_z_y_xyz_yzzz[k] = -g_z_y_yz_yzzz[k] * ab_x + g_z_y_yz_xyzzz[k];

                g_z_y_xyz_zzzz[k] = -g_z_y_yz_zzzz[k] * ab_x + g_z_y_yz_xzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzz_xxxx = cbuffer.data(fg_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_y_xzz_xxxy = cbuffer.data(fg_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_y_xzz_xxxz = cbuffer.data(fg_geom_11_off + 1127 * ccomps * dcomps);

            auto g_z_y_xzz_xxyy = cbuffer.data(fg_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_y_xzz_xxyz = cbuffer.data(fg_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_y_xzz_xxzz = cbuffer.data(fg_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_y_xzz_xyyy = cbuffer.data(fg_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_y_xzz_xyyz = cbuffer.data(fg_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_y_xzz_xyzz = cbuffer.data(fg_geom_11_off + 1133 * ccomps * dcomps);

            auto g_z_y_xzz_xzzz = cbuffer.data(fg_geom_11_off + 1134 * ccomps * dcomps);

            auto g_z_y_xzz_yyyy = cbuffer.data(fg_geom_11_off + 1135 * ccomps * dcomps);

            auto g_z_y_xzz_yyyz = cbuffer.data(fg_geom_11_off + 1136 * ccomps * dcomps);

            auto g_z_y_xzz_yyzz = cbuffer.data(fg_geom_11_off + 1137 * ccomps * dcomps);

            auto g_z_y_xzz_yzzz = cbuffer.data(fg_geom_11_off + 1138 * ccomps * dcomps);

            auto g_z_y_xzz_zzzz = cbuffer.data(fg_geom_11_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzz_xxxx, g_z_y_xzz_xxxy, g_z_y_xzz_xxxz, g_z_y_xzz_xxyy, g_z_y_xzz_xxyz, g_z_y_xzz_xxzz, g_z_y_xzz_xyyy, g_z_y_xzz_xyyz, g_z_y_xzz_xyzz, g_z_y_xzz_xzzz, g_z_y_xzz_yyyy, g_z_y_xzz_yyyz, g_z_y_xzz_yyzz, g_z_y_xzz_yzzz, g_z_y_xzz_zzzz, g_z_y_zz_xxxx, g_z_y_zz_xxxxx, g_z_y_zz_xxxxy, g_z_y_zz_xxxxz, g_z_y_zz_xxxy, g_z_y_zz_xxxyy, g_z_y_zz_xxxyz, g_z_y_zz_xxxz, g_z_y_zz_xxxzz, g_z_y_zz_xxyy, g_z_y_zz_xxyyy, g_z_y_zz_xxyyz, g_z_y_zz_xxyz, g_z_y_zz_xxyzz, g_z_y_zz_xxzz, g_z_y_zz_xxzzz, g_z_y_zz_xyyy, g_z_y_zz_xyyyy, g_z_y_zz_xyyyz, g_z_y_zz_xyyz, g_z_y_zz_xyyzz, g_z_y_zz_xyzz, g_z_y_zz_xyzzz, g_z_y_zz_xzzz, g_z_y_zz_xzzzz, g_z_y_zz_yyyy, g_z_y_zz_yyyz, g_z_y_zz_yyzz, g_z_y_zz_yzzz, g_z_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzz_xxxx[k] = -g_z_y_zz_xxxx[k] * ab_x + g_z_y_zz_xxxxx[k];

                g_z_y_xzz_xxxy[k] = -g_z_y_zz_xxxy[k] * ab_x + g_z_y_zz_xxxxy[k];

                g_z_y_xzz_xxxz[k] = -g_z_y_zz_xxxz[k] * ab_x + g_z_y_zz_xxxxz[k];

                g_z_y_xzz_xxyy[k] = -g_z_y_zz_xxyy[k] * ab_x + g_z_y_zz_xxxyy[k];

                g_z_y_xzz_xxyz[k] = -g_z_y_zz_xxyz[k] * ab_x + g_z_y_zz_xxxyz[k];

                g_z_y_xzz_xxzz[k] = -g_z_y_zz_xxzz[k] * ab_x + g_z_y_zz_xxxzz[k];

                g_z_y_xzz_xyyy[k] = -g_z_y_zz_xyyy[k] * ab_x + g_z_y_zz_xxyyy[k];

                g_z_y_xzz_xyyz[k] = -g_z_y_zz_xyyz[k] * ab_x + g_z_y_zz_xxyyz[k];

                g_z_y_xzz_xyzz[k] = -g_z_y_zz_xyzz[k] * ab_x + g_z_y_zz_xxyzz[k];

                g_z_y_xzz_xzzz[k] = -g_z_y_zz_xzzz[k] * ab_x + g_z_y_zz_xxzzz[k];

                g_z_y_xzz_yyyy[k] = -g_z_y_zz_yyyy[k] * ab_x + g_z_y_zz_xyyyy[k];

                g_z_y_xzz_yyyz[k] = -g_z_y_zz_yyyz[k] * ab_x + g_z_y_zz_xyyyz[k];

                g_z_y_xzz_yyzz[k] = -g_z_y_zz_yyzz[k] * ab_x + g_z_y_zz_xyyzz[k];

                g_z_y_xzz_yzzz[k] = -g_z_y_zz_yzzz[k] * ab_x + g_z_y_zz_xyzzz[k];

                g_z_y_xzz_zzzz[k] = -g_z_y_zz_zzzz[k] * ab_x + g_z_y_zz_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyy_xxxx = cbuffer.data(fg_geom_11_off + 1140 * ccomps * dcomps);

            auto g_z_y_yyy_xxxy = cbuffer.data(fg_geom_11_off + 1141 * ccomps * dcomps);

            auto g_z_y_yyy_xxxz = cbuffer.data(fg_geom_11_off + 1142 * ccomps * dcomps);

            auto g_z_y_yyy_xxyy = cbuffer.data(fg_geom_11_off + 1143 * ccomps * dcomps);

            auto g_z_y_yyy_xxyz = cbuffer.data(fg_geom_11_off + 1144 * ccomps * dcomps);

            auto g_z_y_yyy_xxzz = cbuffer.data(fg_geom_11_off + 1145 * ccomps * dcomps);

            auto g_z_y_yyy_xyyy = cbuffer.data(fg_geom_11_off + 1146 * ccomps * dcomps);

            auto g_z_y_yyy_xyyz = cbuffer.data(fg_geom_11_off + 1147 * ccomps * dcomps);

            auto g_z_y_yyy_xyzz = cbuffer.data(fg_geom_11_off + 1148 * ccomps * dcomps);

            auto g_z_y_yyy_xzzz = cbuffer.data(fg_geom_11_off + 1149 * ccomps * dcomps);

            auto g_z_y_yyy_yyyy = cbuffer.data(fg_geom_11_off + 1150 * ccomps * dcomps);

            auto g_z_y_yyy_yyyz = cbuffer.data(fg_geom_11_off + 1151 * ccomps * dcomps);

            auto g_z_y_yyy_yyzz = cbuffer.data(fg_geom_11_off + 1152 * ccomps * dcomps);

            auto g_z_y_yyy_yzzz = cbuffer.data(fg_geom_11_off + 1153 * ccomps * dcomps);

            auto g_z_y_yyy_zzzz = cbuffer.data(fg_geom_11_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_xxxx, g_z_0_yy_xxxy, g_z_0_yy_xxxz, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxzz, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyzz, g_z_0_yy_xzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyzz, g_z_0_yy_yzzz, g_z_0_yy_zzzz, g_z_y_yy_xxxx, g_z_y_yy_xxxxy, g_z_y_yy_xxxy, g_z_y_yy_xxxyy, g_z_y_yy_xxxyz, g_z_y_yy_xxxz, g_z_y_yy_xxyy, g_z_y_yy_xxyyy, g_z_y_yy_xxyyz, g_z_y_yy_xxyz, g_z_y_yy_xxyzz, g_z_y_yy_xxzz, g_z_y_yy_xyyy, g_z_y_yy_xyyyy, g_z_y_yy_xyyyz, g_z_y_yy_xyyz, g_z_y_yy_xyyzz, g_z_y_yy_xyzz, g_z_y_yy_xyzzz, g_z_y_yy_xzzz, g_z_y_yy_yyyy, g_z_y_yy_yyyyy, g_z_y_yy_yyyyz, g_z_y_yy_yyyz, g_z_y_yy_yyyzz, g_z_y_yy_yyzz, g_z_y_yy_yyzzz, g_z_y_yy_yzzz, g_z_y_yy_yzzzz, g_z_y_yy_zzzz, g_z_y_yyy_xxxx, g_z_y_yyy_xxxy, g_z_y_yyy_xxxz, g_z_y_yyy_xxyy, g_z_y_yyy_xxyz, g_z_y_yyy_xxzz, g_z_y_yyy_xyyy, g_z_y_yyy_xyyz, g_z_y_yyy_xyzz, g_z_y_yyy_xzzz, g_z_y_yyy_yyyy, g_z_y_yyy_yyyz, g_z_y_yyy_yyzz, g_z_y_yyy_yzzz, g_z_y_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyy_xxxx[k] = g_z_0_yy_xxxx[k] - g_z_y_yy_xxxx[k] * ab_y + g_z_y_yy_xxxxy[k];

                g_z_y_yyy_xxxy[k] = g_z_0_yy_xxxy[k] - g_z_y_yy_xxxy[k] * ab_y + g_z_y_yy_xxxyy[k];

                g_z_y_yyy_xxxz[k] = g_z_0_yy_xxxz[k] - g_z_y_yy_xxxz[k] * ab_y + g_z_y_yy_xxxyz[k];

                g_z_y_yyy_xxyy[k] = g_z_0_yy_xxyy[k] - g_z_y_yy_xxyy[k] * ab_y + g_z_y_yy_xxyyy[k];

                g_z_y_yyy_xxyz[k] = g_z_0_yy_xxyz[k] - g_z_y_yy_xxyz[k] * ab_y + g_z_y_yy_xxyyz[k];

                g_z_y_yyy_xxzz[k] = g_z_0_yy_xxzz[k] - g_z_y_yy_xxzz[k] * ab_y + g_z_y_yy_xxyzz[k];

                g_z_y_yyy_xyyy[k] = g_z_0_yy_xyyy[k] - g_z_y_yy_xyyy[k] * ab_y + g_z_y_yy_xyyyy[k];

                g_z_y_yyy_xyyz[k] = g_z_0_yy_xyyz[k] - g_z_y_yy_xyyz[k] * ab_y + g_z_y_yy_xyyyz[k];

                g_z_y_yyy_xyzz[k] = g_z_0_yy_xyzz[k] - g_z_y_yy_xyzz[k] * ab_y + g_z_y_yy_xyyzz[k];

                g_z_y_yyy_xzzz[k] = g_z_0_yy_xzzz[k] - g_z_y_yy_xzzz[k] * ab_y + g_z_y_yy_xyzzz[k];

                g_z_y_yyy_yyyy[k] = g_z_0_yy_yyyy[k] - g_z_y_yy_yyyy[k] * ab_y + g_z_y_yy_yyyyy[k];

                g_z_y_yyy_yyyz[k] = g_z_0_yy_yyyz[k] - g_z_y_yy_yyyz[k] * ab_y + g_z_y_yy_yyyyz[k];

                g_z_y_yyy_yyzz[k] = g_z_0_yy_yyzz[k] - g_z_y_yy_yyzz[k] * ab_y + g_z_y_yy_yyyzz[k];

                g_z_y_yyy_yzzz[k] = g_z_0_yy_yzzz[k] - g_z_y_yy_yzzz[k] * ab_y + g_z_y_yy_yyzzz[k];

                g_z_y_yyy_zzzz[k] = g_z_0_yy_zzzz[k] - g_z_y_yy_zzzz[k] * ab_y + g_z_y_yy_yzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyz_xxxx = cbuffer.data(fg_geom_11_off + 1155 * ccomps * dcomps);

            auto g_z_y_yyz_xxxy = cbuffer.data(fg_geom_11_off + 1156 * ccomps * dcomps);

            auto g_z_y_yyz_xxxz = cbuffer.data(fg_geom_11_off + 1157 * ccomps * dcomps);

            auto g_z_y_yyz_xxyy = cbuffer.data(fg_geom_11_off + 1158 * ccomps * dcomps);

            auto g_z_y_yyz_xxyz = cbuffer.data(fg_geom_11_off + 1159 * ccomps * dcomps);

            auto g_z_y_yyz_xxzz = cbuffer.data(fg_geom_11_off + 1160 * ccomps * dcomps);

            auto g_z_y_yyz_xyyy = cbuffer.data(fg_geom_11_off + 1161 * ccomps * dcomps);

            auto g_z_y_yyz_xyyz = cbuffer.data(fg_geom_11_off + 1162 * ccomps * dcomps);

            auto g_z_y_yyz_xyzz = cbuffer.data(fg_geom_11_off + 1163 * ccomps * dcomps);

            auto g_z_y_yyz_xzzz = cbuffer.data(fg_geom_11_off + 1164 * ccomps * dcomps);

            auto g_z_y_yyz_yyyy = cbuffer.data(fg_geom_11_off + 1165 * ccomps * dcomps);

            auto g_z_y_yyz_yyyz = cbuffer.data(fg_geom_11_off + 1166 * ccomps * dcomps);

            auto g_z_y_yyz_yyzz = cbuffer.data(fg_geom_11_off + 1167 * ccomps * dcomps);

            auto g_z_y_yyz_yzzz = cbuffer.data(fg_geom_11_off + 1168 * ccomps * dcomps);

            auto g_z_y_yyz_zzzz = cbuffer.data(fg_geom_11_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xxxx, g_z_0_yz_xxxy, g_z_0_yz_xxxz, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxzz, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyzz, g_z_0_yz_xzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyzz, g_z_0_yz_yzzz, g_z_0_yz_zzzz, g_z_y_yyz_xxxx, g_z_y_yyz_xxxy, g_z_y_yyz_xxxz, g_z_y_yyz_xxyy, g_z_y_yyz_xxyz, g_z_y_yyz_xxzz, g_z_y_yyz_xyyy, g_z_y_yyz_xyyz, g_z_y_yyz_xyzz, g_z_y_yyz_xzzz, g_z_y_yyz_yyyy, g_z_y_yyz_yyyz, g_z_y_yyz_yyzz, g_z_y_yyz_yzzz, g_z_y_yyz_zzzz, g_z_y_yz_xxxx, g_z_y_yz_xxxxy, g_z_y_yz_xxxy, g_z_y_yz_xxxyy, g_z_y_yz_xxxyz, g_z_y_yz_xxxz, g_z_y_yz_xxyy, g_z_y_yz_xxyyy, g_z_y_yz_xxyyz, g_z_y_yz_xxyz, g_z_y_yz_xxyzz, g_z_y_yz_xxzz, g_z_y_yz_xyyy, g_z_y_yz_xyyyy, g_z_y_yz_xyyyz, g_z_y_yz_xyyz, g_z_y_yz_xyyzz, g_z_y_yz_xyzz, g_z_y_yz_xyzzz, g_z_y_yz_xzzz, g_z_y_yz_yyyy, g_z_y_yz_yyyyy, g_z_y_yz_yyyyz, g_z_y_yz_yyyz, g_z_y_yz_yyyzz, g_z_y_yz_yyzz, g_z_y_yz_yyzzz, g_z_y_yz_yzzz, g_z_y_yz_yzzzz, g_z_y_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyz_xxxx[k] = g_z_0_yz_xxxx[k] - g_z_y_yz_xxxx[k] * ab_y + g_z_y_yz_xxxxy[k];

                g_z_y_yyz_xxxy[k] = g_z_0_yz_xxxy[k] - g_z_y_yz_xxxy[k] * ab_y + g_z_y_yz_xxxyy[k];

                g_z_y_yyz_xxxz[k] = g_z_0_yz_xxxz[k] - g_z_y_yz_xxxz[k] * ab_y + g_z_y_yz_xxxyz[k];

                g_z_y_yyz_xxyy[k] = g_z_0_yz_xxyy[k] - g_z_y_yz_xxyy[k] * ab_y + g_z_y_yz_xxyyy[k];

                g_z_y_yyz_xxyz[k] = g_z_0_yz_xxyz[k] - g_z_y_yz_xxyz[k] * ab_y + g_z_y_yz_xxyyz[k];

                g_z_y_yyz_xxzz[k] = g_z_0_yz_xxzz[k] - g_z_y_yz_xxzz[k] * ab_y + g_z_y_yz_xxyzz[k];

                g_z_y_yyz_xyyy[k] = g_z_0_yz_xyyy[k] - g_z_y_yz_xyyy[k] * ab_y + g_z_y_yz_xyyyy[k];

                g_z_y_yyz_xyyz[k] = g_z_0_yz_xyyz[k] - g_z_y_yz_xyyz[k] * ab_y + g_z_y_yz_xyyyz[k];

                g_z_y_yyz_xyzz[k] = g_z_0_yz_xyzz[k] - g_z_y_yz_xyzz[k] * ab_y + g_z_y_yz_xyyzz[k];

                g_z_y_yyz_xzzz[k] = g_z_0_yz_xzzz[k] - g_z_y_yz_xzzz[k] * ab_y + g_z_y_yz_xyzzz[k];

                g_z_y_yyz_yyyy[k] = g_z_0_yz_yyyy[k] - g_z_y_yz_yyyy[k] * ab_y + g_z_y_yz_yyyyy[k];

                g_z_y_yyz_yyyz[k] = g_z_0_yz_yyyz[k] - g_z_y_yz_yyyz[k] * ab_y + g_z_y_yz_yyyyz[k];

                g_z_y_yyz_yyzz[k] = g_z_0_yz_yyzz[k] - g_z_y_yz_yyzz[k] * ab_y + g_z_y_yz_yyyzz[k];

                g_z_y_yyz_yzzz[k] = g_z_0_yz_yzzz[k] - g_z_y_yz_yzzz[k] * ab_y + g_z_y_yz_yyzzz[k];

                g_z_y_yyz_zzzz[k] = g_z_0_yz_zzzz[k] - g_z_y_yz_zzzz[k] * ab_y + g_z_y_yz_yzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzz_xxxx = cbuffer.data(fg_geom_11_off + 1170 * ccomps * dcomps);

            auto g_z_y_yzz_xxxy = cbuffer.data(fg_geom_11_off + 1171 * ccomps * dcomps);

            auto g_z_y_yzz_xxxz = cbuffer.data(fg_geom_11_off + 1172 * ccomps * dcomps);

            auto g_z_y_yzz_xxyy = cbuffer.data(fg_geom_11_off + 1173 * ccomps * dcomps);

            auto g_z_y_yzz_xxyz = cbuffer.data(fg_geom_11_off + 1174 * ccomps * dcomps);

            auto g_z_y_yzz_xxzz = cbuffer.data(fg_geom_11_off + 1175 * ccomps * dcomps);

            auto g_z_y_yzz_xyyy = cbuffer.data(fg_geom_11_off + 1176 * ccomps * dcomps);

            auto g_z_y_yzz_xyyz = cbuffer.data(fg_geom_11_off + 1177 * ccomps * dcomps);

            auto g_z_y_yzz_xyzz = cbuffer.data(fg_geom_11_off + 1178 * ccomps * dcomps);

            auto g_z_y_yzz_xzzz = cbuffer.data(fg_geom_11_off + 1179 * ccomps * dcomps);

            auto g_z_y_yzz_yyyy = cbuffer.data(fg_geom_11_off + 1180 * ccomps * dcomps);

            auto g_z_y_yzz_yyyz = cbuffer.data(fg_geom_11_off + 1181 * ccomps * dcomps);

            auto g_z_y_yzz_yyzz = cbuffer.data(fg_geom_11_off + 1182 * ccomps * dcomps);

            auto g_z_y_yzz_yzzz = cbuffer.data(fg_geom_11_off + 1183 * ccomps * dcomps);

            auto g_z_y_yzz_zzzz = cbuffer.data(fg_geom_11_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_z_y_yzz_xxxx, g_z_y_yzz_xxxy, g_z_y_yzz_xxxz, g_z_y_yzz_xxyy, g_z_y_yzz_xxyz, g_z_y_yzz_xxzz, g_z_y_yzz_xyyy, g_z_y_yzz_xyyz, g_z_y_yzz_xyzz, g_z_y_yzz_xzzz, g_z_y_yzz_yyyy, g_z_y_yzz_yyyz, g_z_y_yzz_yyzz, g_z_y_yzz_yzzz, g_z_y_yzz_zzzz, g_z_y_zz_xxxx, g_z_y_zz_xxxxy, g_z_y_zz_xxxy, g_z_y_zz_xxxyy, g_z_y_zz_xxxyz, g_z_y_zz_xxxz, g_z_y_zz_xxyy, g_z_y_zz_xxyyy, g_z_y_zz_xxyyz, g_z_y_zz_xxyz, g_z_y_zz_xxyzz, g_z_y_zz_xxzz, g_z_y_zz_xyyy, g_z_y_zz_xyyyy, g_z_y_zz_xyyyz, g_z_y_zz_xyyz, g_z_y_zz_xyyzz, g_z_y_zz_xyzz, g_z_y_zz_xyzzz, g_z_y_zz_xzzz, g_z_y_zz_yyyy, g_z_y_zz_yyyyy, g_z_y_zz_yyyyz, g_z_y_zz_yyyz, g_z_y_zz_yyyzz, g_z_y_zz_yyzz, g_z_y_zz_yyzzz, g_z_y_zz_yzzz, g_z_y_zz_yzzzz, g_z_y_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzz_xxxx[k] = g_z_0_zz_xxxx[k] - g_z_y_zz_xxxx[k] * ab_y + g_z_y_zz_xxxxy[k];

                g_z_y_yzz_xxxy[k] = g_z_0_zz_xxxy[k] - g_z_y_zz_xxxy[k] * ab_y + g_z_y_zz_xxxyy[k];

                g_z_y_yzz_xxxz[k] = g_z_0_zz_xxxz[k] - g_z_y_zz_xxxz[k] * ab_y + g_z_y_zz_xxxyz[k];

                g_z_y_yzz_xxyy[k] = g_z_0_zz_xxyy[k] - g_z_y_zz_xxyy[k] * ab_y + g_z_y_zz_xxyyy[k];

                g_z_y_yzz_xxyz[k] = g_z_0_zz_xxyz[k] - g_z_y_zz_xxyz[k] * ab_y + g_z_y_zz_xxyyz[k];

                g_z_y_yzz_xxzz[k] = g_z_0_zz_xxzz[k] - g_z_y_zz_xxzz[k] * ab_y + g_z_y_zz_xxyzz[k];

                g_z_y_yzz_xyyy[k] = g_z_0_zz_xyyy[k] - g_z_y_zz_xyyy[k] * ab_y + g_z_y_zz_xyyyy[k];

                g_z_y_yzz_xyyz[k] = g_z_0_zz_xyyz[k] - g_z_y_zz_xyyz[k] * ab_y + g_z_y_zz_xyyyz[k];

                g_z_y_yzz_xyzz[k] = g_z_0_zz_xyzz[k] - g_z_y_zz_xyzz[k] * ab_y + g_z_y_zz_xyyzz[k];

                g_z_y_yzz_xzzz[k] = g_z_0_zz_xzzz[k] - g_z_y_zz_xzzz[k] * ab_y + g_z_y_zz_xyzzz[k];

                g_z_y_yzz_yyyy[k] = g_z_0_zz_yyyy[k] - g_z_y_zz_yyyy[k] * ab_y + g_z_y_zz_yyyyy[k];

                g_z_y_yzz_yyyz[k] = g_z_0_zz_yyyz[k] - g_z_y_zz_yyyz[k] * ab_y + g_z_y_zz_yyyyz[k];

                g_z_y_yzz_yyzz[k] = g_z_0_zz_yyzz[k] - g_z_y_zz_yyzz[k] * ab_y + g_z_y_zz_yyyzz[k];

                g_z_y_yzz_yzzz[k] = g_z_0_zz_yzzz[k] - g_z_y_zz_yzzz[k] * ab_y + g_z_y_zz_yyzzz[k];

                g_z_y_yzz_zzzz[k] = g_z_0_zz_zzzz[k] - g_z_y_zz_zzzz[k] * ab_y + g_z_y_zz_yzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzz_xxxx = cbuffer.data(fg_geom_11_off + 1185 * ccomps * dcomps);

            auto g_z_y_zzz_xxxy = cbuffer.data(fg_geom_11_off + 1186 * ccomps * dcomps);

            auto g_z_y_zzz_xxxz = cbuffer.data(fg_geom_11_off + 1187 * ccomps * dcomps);

            auto g_z_y_zzz_xxyy = cbuffer.data(fg_geom_11_off + 1188 * ccomps * dcomps);

            auto g_z_y_zzz_xxyz = cbuffer.data(fg_geom_11_off + 1189 * ccomps * dcomps);

            auto g_z_y_zzz_xxzz = cbuffer.data(fg_geom_11_off + 1190 * ccomps * dcomps);

            auto g_z_y_zzz_xyyy = cbuffer.data(fg_geom_11_off + 1191 * ccomps * dcomps);

            auto g_z_y_zzz_xyyz = cbuffer.data(fg_geom_11_off + 1192 * ccomps * dcomps);

            auto g_z_y_zzz_xyzz = cbuffer.data(fg_geom_11_off + 1193 * ccomps * dcomps);

            auto g_z_y_zzz_xzzz = cbuffer.data(fg_geom_11_off + 1194 * ccomps * dcomps);

            auto g_z_y_zzz_yyyy = cbuffer.data(fg_geom_11_off + 1195 * ccomps * dcomps);

            auto g_z_y_zzz_yyyz = cbuffer.data(fg_geom_11_off + 1196 * ccomps * dcomps);

            auto g_z_y_zzz_yyzz = cbuffer.data(fg_geom_11_off + 1197 * ccomps * dcomps);

            auto g_z_y_zzz_yzzz = cbuffer.data(fg_geom_11_off + 1198 * ccomps * dcomps);

            auto g_z_y_zzz_zzzz = cbuffer.data(fg_geom_11_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_xxxx, g_0_y_zz_xxxy, g_0_y_zz_xxxz, g_0_y_zz_xxyy, g_0_y_zz_xxyz, g_0_y_zz_xxzz, g_0_y_zz_xyyy, g_0_y_zz_xyyz, g_0_y_zz_xyzz, g_0_y_zz_xzzz, g_0_y_zz_yyyy, g_0_y_zz_yyyz, g_0_y_zz_yyzz, g_0_y_zz_yzzz, g_0_y_zz_zzzz, g_z_y_zz_xxxx, g_z_y_zz_xxxxz, g_z_y_zz_xxxy, g_z_y_zz_xxxyz, g_z_y_zz_xxxz, g_z_y_zz_xxxzz, g_z_y_zz_xxyy, g_z_y_zz_xxyyz, g_z_y_zz_xxyz, g_z_y_zz_xxyzz, g_z_y_zz_xxzz, g_z_y_zz_xxzzz, g_z_y_zz_xyyy, g_z_y_zz_xyyyz, g_z_y_zz_xyyz, g_z_y_zz_xyyzz, g_z_y_zz_xyzz, g_z_y_zz_xyzzz, g_z_y_zz_xzzz, g_z_y_zz_xzzzz, g_z_y_zz_yyyy, g_z_y_zz_yyyyz, g_z_y_zz_yyyz, g_z_y_zz_yyyzz, g_z_y_zz_yyzz, g_z_y_zz_yyzzz, g_z_y_zz_yzzz, g_z_y_zz_yzzzz, g_z_y_zz_zzzz, g_z_y_zz_zzzzz, g_z_y_zzz_xxxx, g_z_y_zzz_xxxy, g_z_y_zzz_xxxz, g_z_y_zzz_xxyy, g_z_y_zzz_xxyz, g_z_y_zzz_xxzz, g_z_y_zzz_xyyy, g_z_y_zzz_xyyz, g_z_y_zzz_xyzz, g_z_y_zzz_xzzz, g_z_y_zzz_yyyy, g_z_y_zzz_yyyz, g_z_y_zzz_yyzz, g_z_y_zzz_yzzz, g_z_y_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzz_xxxx[k] = -g_0_y_zz_xxxx[k] - g_z_y_zz_xxxx[k] * ab_z + g_z_y_zz_xxxxz[k];

                g_z_y_zzz_xxxy[k] = -g_0_y_zz_xxxy[k] - g_z_y_zz_xxxy[k] * ab_z + g_z_y_zz_xxxyz[k];

                g_z_y_zzz_xxxz[k] = -g_0_y_zz_xxxz[k] - g_z_y_zz_xxxz[k] * ab_z + g_z_y_zz_xxxzz[k];

                g_z_y_zzz_xxyy[k] = -g_0_y_zz_xxyy[k] - g_z_y_zz_xxyy[k] * ab_z + g_z_y_zz_xxyyz[k];

                g_z_y_zzz_xxyz[k] = -g_0_y_zz_xxyz[k] - g_z_y_zz_xxyz[k] * ab_z + g_z_y_zz_xxyzz[k];

                g_z_y_zzz_xxzz[k] = -g_0_y_zz_xxzz[k] - g_z_y_zz_xxzz[k] * ab_z + g_z_y_zz_xxzzz[k];

                g_z_y_zzz_xyyy[k] = -g_0_y_zz_xyyy[k] - g_z_y_zz_xyyy[k] * ab_z + g_z_y_zz_xyyyz[k];

                g_z_y_zzz_xyyz[k] = -g_0_y_zz_xyyz[k] - g_z_y_zz_xyyz[k] * ab_z + g_z_y_zz_xyyzz[k];

                g_z_y_zzz_xyzz[k] = -g_0_y_zz_xyzz[k] - g_z_y_zz_xyzz[k] * ab_z + g_z_y_zz_xyzzz[k];

                g_z_y_zzz_xzzz[k] = -g_0_y_zz_xzzz[k] - g_z_y_zz_xzzz[k] * ab_z + g_z_y_zz_xzzzz[k];

                g_z_y_zzz_yyyy[k] = -g_0_y_zz_yyyy[k] - g_z_y_zz_yyyy[k] * ab_z + g_z_y_zz_yyyyz[k];

                g_z_y_zzz_yyyz[k] = -g_0_y_zz_yyyz[k] - g_z_y_zz_yyyz[k] * ab_z + g_z_y_zz_yyyzz[k];

                g_z_y_zzz_yyzz[k] = -g_0_y_zz_yyzz[k] - g_z_y_zz_yyzz[k] * ab_z + g_z_y_zz_yyzzz[k];

                g_z_y_zzz_yzzz[k] = -g_0_y_zz_yzzz[k] - g_z_y_zz_yzzz[k] * ab_z + g_z_y_zz_yzzzz[k];

                g_z_y_zzz_zzzz[k] = -g_0_y_zz_zzzz[k] - g_z_y_zz_zzzz[k] * ab_z + g_z_y_zz_zzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxx_xxxx = cbuffer.data(fg_geom_11_off + 1200 * ccomps * dcomps);

            auto g_z_z_xxx_xxxy = cbuffer.data(fg_geom_11_off + 1201 * ccomps * dcomps);

            auto g_z_z_xxx_xxxz = cbuffer.data(fg_geom_11_off + 1202 * ccomps * dcomps);

            auto g_z_z_xxx_xxyy = cbuffer.data(fg_geom_11_off + 1203 * ccomps * dcomps);

            auto g_z_z_xxx_xxyz = cbuffer.data(fg_geom_11_off + 1204 * ccomps * dcomps);

            auto g_z_z_xxx_xxzz = cbuffer.data(fg_geom_11_off + 1205 * ccomps * dcomps);

            auto g_z_z_xxx_xyyy = cbuffer.data(fg_geom_11_off + 1206 * ccomps * dcomps);

            auto g_z_z_xxx_xyyz = cbuffer.data(fg_geom_11_off + 1207 * ccomps * dcomps);

            auto g_z_z_xxx_xyzz = cbuffer.data(fg_geom_11_off + 1208 * ccomps * dcomps);

            auto g_z_z_xxx_xzzz = cbuffer.data(fg_geom_11_off + 1209 * ccomps * dcomps);

            auto g_z_z_xxx_yyyy = cbuffer.data(fg_geom_11_off + 1210 * ccomps * dcomps);

            auto g_z_z_xxx_yyyz = cbuffer.data(fg_geom_11_off + 1211 * ccomps * dcomps);

            auto g_z_z_xxx_yyzz = cbuffer.data(fg_geom_11_off + 1212 * ccomps * dcomps);

            auto g_z_z_xxx_yzzz = cbuffer.data(fg_geom_11_off + 1213 * ccomps * dcomps);

            auto g_z_z_xxx_zzzz = cbuffer.data(fg_geom_11_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xx_xxxx, g_z_z_xx_xxxxx, g_z_z_xx_xxxxy, g_z_z_xx_xxxxz, g_z_z_xx_xxxy, g_z_z_xx_xxxyy, g_z_z_xx_xxxyz, g_z_z_xx_xxxz, g_z_z_xx_xxxzz, g_z_z_xx_xxyy, g_z_z_xx_xxyyy, g_z_z_xx_xxyyz, g_z_z_xx_xxyz, g_z_z_xx_xxyzz, g_z_z_xx_xxzz, g_z_z_xx_xxzzz, g_z_z_xx_xyyy, g_z_z_xx_xyyyy, g_z_z_xx_xyyyz, g_z_z_xx_xyyz, g_z_z_xx_xyyzz, g_z_z_xx_xyzz, g_z_z_xx_xyzzz, g_z_z_xx_xzzz, g_z_z_xx_xzzzz, g_z_z_xx_yyyy, g_z_z_xx_yyyz, g_z_z_xx_yyzz, g_z_z_xx_yzzz, g_z_z_xx_zzzz, g_z_z_xxx_xxxx, g_z_z_xxx_xxxy, g_z_z_xxx_xxxz, g_z_z_xxx_xxyy, g_z_z_xxx_xxyz, g_z_z_xxx_xxzz, g_z_z_xxx_xyyy, g_z_z_xxx_xyyz, g_z_z_xxx_xyzz, g_z_z_xxx_xzzz, g_z_z_xxx_yyyy, g_z_z_xxx_yyyz, g_z_z_xxx_yyzz, g_z_z_xxx_yzzz, g_z_z_xxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxx_xxxx[k] = -g_z_z_xx_xxxx[k] * ab_x + g_z_z_xx_xxxxx[k];

                g_z_z_xxx_xxxy[k] = -g_z_z_xx_xxxy[k] * ab_x + g_z_z_xx_xxxxy[k];

                g_z_z_xxx_xxxz[k] = -g_z_z_xx_xxxz[k] * ab_x + g_z_z_xx_xxxxz[k];

                g_z_z_xxx_xxyy[k] = -g_z_z_xx_xxyy[k] * ab_x + g_z_z_xx_xxxyy[k];

                g_z_z_xxx_xxyz[k] = -g_z_z_xx_xxyz[k] * ab_x + g_z_z_xx_xxxyz[k];

                g_z_z_xxx_xxzz[k] = -g_z_z_xx_xxzz[k] * ab_x + g_z_z_xx_xxxzz[k];

                g_z_z_xxx_xyyy[k] = -g_z_z_xx_xyyy[k] * ab_x + g_z_z_xx_xxyyy[k];

                g_z_z_xxx_xyyz[k] = -g_z_z_xx_xyyz[k] * ab_x + g_z_z_xx_xxyyz[k];

                g_z_z_xxx_xyzz[k] = -g_z_z_xx_xyzz[k] * ab_x + g_z_z_xx_xxyzz[k];

                g_z_z_xxx_xzzz[k] = -g_z_z_xx_xzzz[k] * ab_x + g_z_z_xx_xxzzz[k];

                g_z_z_xxx_yyyy[k] = -g_z_z_xx_yyyy[k] * ab_x + g_z_z_xx_xyyyy[k];

                g_z_z_xxx_yyyz[k] = -g_z_z_xx_yyyz[k] * ab_x + g_z_z_xx_xyyyz[k];

                g_z_z_xxx_yyzz[k] = -g_z_z_xx_yyzz[k] * ab_x + g_z_z_xx_xyyzz[k];

                g_z_z_xxx_yzzz[k] = -g_z_z_xx_yzzz[k] * ab_x + g_z_z_xx_xyzzz[k];

                g_z_z_xxx_zzzz[k] = -g_z_z_xx_zzzz[k] * ab_x + g_z_z_xx_xzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxy_xxxx = cbuffer.data(fg_geom_11_off + 1215 * ccomps * dcomps);

            auto g_z_z_xxy_xxxy = cbuffer.data(fg_geom_11_off + 1216 * ccomps * dcomps);

            auto g_z_z_xxy_xxxz = cbuffer.data(fg_geom_11_off + 1217 * ccomps * dcomps);

            auto g_z_z_xxy_xxyy = cbuffer.data(fg_geom_11_off + 1218 * ccomps * dcomps);

            auto g_z_z_xxy_xxyz = cbuffer.data(fg_geom_11_off + 1219 * ccomps * dcomps);

            auto g_z_z_xxy_xxzz = cbuffer.data(fg_geom_11_off + 1220 * ccomps * dcomps);

            auto g_z_z_xxy_xyyy = cbuffer.data(fg_geom_11_off + 1221 * ccomps * dcomps);

            auto g_z_z_xxy_xyyz = cbuffer.data(fg_geom_11_off + 1222 * ccomps * dcomps);

            auto g_z_z_xxy_xyzz = cbuffer.data(fg_geom_11_off + 1223 * ccomps * dcomps);

            auto g_z_z_xxy_xzzz = cbuffer.data(fg_geom_11_off + 1224 * ccomps * dcomps);

            auto g_z_z_xxy_yyyy = cbuffer.data(fg_geom_11_off + 1225 * ccomps * dcomps);

            auto g_z_z_xxy_yyyz = cbuffer.data(fg_geom_11_off + 1226 * ccomps * dcomps);

            auto g_z_z_xxy_yyzz = cbuffer.data(fg_geom_11_off + 1227 * ccomps * dcomps);

            auto g_z_z_xxy_yzzz = cbuffer.data(fg_geom_11_off + 1228 * ccomps * dcomps);

            auto g_z_z_xxy_zzzz = cbuffer.data(fg_geom_11_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxy_xxxx, g_z_z_xxy_xxxy, g_z_z_xxy_xxxz, g_z_z_xxy_xxyy, g_z_z_xxy_xxyz, g_z_z_xxy_xxzz, g_z_z_xxy_xyyy, g_z_z_xxy_xyyz, g_z_z_xxy_xyzz, g_z_z_xxy_xzzz, g_z_z_xxy_yyyy, g_z_z_xxy_yyyz, g_z_z_xxy_yyzz, g_z_z_xxy_yzzz, g_z_z_xxy_zzzz, g_z_z_xy_xxxx, g_z_z_xy_xxxxx, g_z_z_xy_xxxxy, g_z_z_xy_xxxxz, g_z_z_xy_xxxy, g_z_z_xy_xxxyy, g_z_z_xy_xxxyz, g_z_z_xy_xxxz, g_z_z_xy_xxxzz, g_z_z_xy_xxyy, g_z_z_xy_xxyyy, g_z_z_xy_xxyyz, g_z_z_xy_xxyz, g_z_z_xy_xxyzz, g_z_z_xy_xxzz, g_z_z_xy_xxzzz, g_z_z_xy_xyyy, g_z_z_xy_xyyyy, g_z_z_xy_xyyyz, g_z_z_xy_xyyz, g_z_z_xy_xyyzz, g_z_z_xy_xyzz, g_z_z_xy_xyzzz, g_z_z_xy_xzzz, g_z_z_xy_xzzzz, g_z_z_xy_yyyy, g_z_z_xy_yyyz, g_z_z_xy_yyzz, g_z_z_xy_yzzz, g_z_z_xy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxy_xxxx[k] = -g_z_z_xy_xxxx[k] * ab_x + g_z_z_xy_xxxxx[k];

                g_z_z_xxy_xxxy[k] = -g_z_z_xy_xxxy[k] * ab_x + g_z_z_xy_xxxxy[k];

                g_z_z_xxy_xxxz[k] = -g_z_z_xy_xxxz[k] * ab_x + g_z_z_xy_xxxxz[k];

                g_z_z_xxy_xxyy[k] = -g_z_z_xy_xxyy[k] * ab_x + g_z_z_xy_xxxyy[k];

                g_z_z_xxy_xxyz[k] = -g_z_z_xy_xxyz[k] * ab_x + g_z_z_xy_xxxyz[k];

                g_z_z_xxy_xxzz[k] = -g_z_z_xy_xxzz[k] * ab_x + g_z_z_xy_xxxzz[k];

                g_z_z_xxy_xyyy[k] = -g_z_z_xy_xyyy[k] * ab_x + g_z_z_xy_xxyyy[k];

                g_z_z_xxy_xyyz[k] = -g_z_z_xy_xyyz[k] * ab_x + g_z_z_xy_xxyyz[k];

                g_z_z_xxy_xyzz[k] = -g_z_z_xy_xyzz[k] * ab_x + g_z_z_xy_xxyzz[k];

                g_z_z_xxy_xzzz[k] = -g_z_z_xy_xzzz[k] * ab_x + g_z_z_xy_xxzzz[k];

                g_z_z_xxy_yyyy[k] = -g_z_z_xy_yyyy[k] * ab_x + g_z_z_xy_xyyyy[k];

                g_z_z_xxy_yyyz[k] = -g_z_z_xy_yyyz[k] * ab_x + g_z_z_xy_xyyyz[k];

                g_z_z_xxy_yyzz[k] = -g_z_z_xy_yyzz[k] * ab_x + g_z_z_xy_xyyzz[k];

                g_z_z_xxy_yzzz[k] = -g_z_z_xy_yzzz[k] * ab_x + g_z_z_xy_xyzzz[k];

                g_z_z_xxy_zzzz[k] = -g_z_z_xy_zzzz[k] * ab_x + g_z_z_xy_xzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxz_xxxx = cbuffer.data(fg_geom_11_off + 1230 * ccomps * dcomps);

            auto g_z_z_xxz_xxxy = cbuffer.data(fg_geom_11_off + 1231 * ccomps * dcomps);

            auto g_z_z_xxz_xxxz = cbuffer.data(fg_geom_11_off + 1232 * ccomps * dcomps);

            auto g_z_z_xxz_xxyy = cbuffer.data(fg_geom_11_off + 1233 * ccomps * dcomps);

            auto g_z_z_xxz_xxyz = cbuffer.data(fg_geom_11_off + 1234 * ccomps * dcomps);

            auto g_z_z_xxz_xxzz = cbuffer.data(fg_geom_11_off + 1235 * ccomps * dcomps);

            auto g_z_z_xxz_xyyy = cbuffer.data(fg_geom_11_off + 1236 * ccomps * dcomps);

            auto g_z_z_xxz_xyyz = cbuffer.data(fg_geom_11_off + 1237 * ccomps * dcomps);

            auto g_z_z_xxz_xyzz = cbuffer.data(fg_geom_11_off + 1238 * ccomps * dcomps);

            auto g_z_z_xxz_xzzz = cbuffer.data(fg_geom_11_off + 1239 * ccomps * dcomps);

            auto g_z_z_xxz_yyyy = cbuffer.data(fg_geom_11_off + 1240 * ccomps * dcomps);

            auto g_z_z_xxz_yyyz = cbuffer.data(fg_geom_11_off + 1241 * ccomps * dcomps);

            auto g_z_z_xxz_yyzz = cbuffer.data(fg_geom_11_off + 1242 * ccomps * dcomps);

            auto g_z_z_xxz_yzzz = cbuffer.data(fg_geom_11_off + 1243 * ccomps * dcomps);

            auto g_z_z_xxz_zzzz = cbuffer.data(fg_geom_11_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxz_xxxx, g_z_z_xxz_xxxy, g_z_z_xxz_xxxz, g_z_z_xxz_xxyy, g_z_z_xxz_xxyz, g_z_z_xxz_xxzz, g_z_z_xxz_xyyy, g_z_z_xxz_xyyz, g_z_z_xxz_xyzz, g_z_z_xxz_xzzz, g_z_z_xxz_yyyy, g_z_z_xxz_yyyz, g_z_z_xxz_yyzz, g_z_z_xxz_yzzz, g_z_z_xxz_zzzz, g_z_z_xz_xxxx, g_z_z_xz_xxxxx, g_z_z_xz_xxxxy, g_z_z_xz_xxxxz, g_z_z_xz_xxxy, g_z_z_xz_xxxyy, g_z_z_xz_xxxyz, g_z_z_xz_xxxz, g_z_z_xz_xxxzz, g_z_z_xz_xxyy, g_z_z_xz_xxyyy, g_z_z_xz_xxyyz, g_z_z_xz_xxyz, g_z_z_xz_xxyzz, g_z_z_xz_xxzz, g_z_z_xz_xxzzz, g_z_z_xz_xyyy, g_z_z_xz_xyyyy, g_z_z_xz_xyyyz, g_z_z_xz_xyyz, g_z_z_xz_xyyzz, g_z_z_xz_xyzz, g_z_z_xz_xyzzz, g_z_z_xz_xzzz, g_z_z_xz_xzzzz, g_z_z_xz_yyyy, g_z_z_xz_yyyz, g_z_z_xz_yyzz, g_z_z_xz_yzzz, g_z_z_xz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxz_xxxx[k] = -g_z_z_xz_xxxx[k] * ab_x + g_z_z_xz_xxxxx[k];

                g_z_z_xxz_xxxy[k] = -g_z_z_xz_xxxy[k] * ab_x + g_z_z_xz_xxxxy[k];

                g_z_z_xxz_xxxz[k] = -g_z_z_xz_xxxz[k] * ab_x + g_z_z_xz_xxxxz[k];

                g_z_z_xxz_xxyy[k] = -g_z_z_xz_xxyy[k] * ab_x + g_z_z_xz_xxxyy[k];

                g_z_z_xxz_xxyz[k] = -g_z_z_xz_xxyz[k] * ab_x + g_z_z_xz_xxxyz[k];

                g_z_z_xxz_xxzz[k] = -g_z_z_xz_xxzz[k] * ab_x + g_z_z_xz_xxxzz[k];

                g_z_z_xxz_xyyy[k] = -g_z_z_xz_xyyy[k] * ab_x + g_z_z_xz_xxyyy[k];

                g_z_z_xxz_xyyz[k] = -g_z_z_xz_xyyz[k] * ab_x + g_z_z_xz_xxyyz[k];

                g_z_z_xxz_xyzz[k] = -g_z_z_xz_xyzz[k] * ab_x + g_z_z_xz_xxyzz[k];

                g_z_z_xxz_xzzz[k] = -g_z_z_xz_xzzz[k] * ab_x + g_z_z_xz_xxzzz[k];

                g_z_z_xxz_yyyy[k] = -g_z_z_xz_yyyy[k] * ab_x + g_z_z_xz_xyyyy[k];

                g_z_z_xxz_yyyz[k] = -g_z_z_xz_yyyz[k] * ab_x + g_z_z_xz_xyyyz[k];

                g_z_z_xxz_yyzz[k] = -g_z_z_xz_yyzz[k] * ab_x + g_z_z_xz_xyyzz[k];

                g_z_z_xxz_yzzz[k] = -g_z_z_xz_yzzz[k] * ab_x + g_z_z_xz_xyzzz[k];

                g_z_z_xxz_zzzz[k] = -g_z_z_xz_zzzz[k] * ab_x + g_z_z_xz_xzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyy_xxxx = cbuffer.data(fg_geom_11_off + 1245 * ccomps * dcomps);

            auto g_z_z_xyy_xxxy = cbuffer.data(fg_geom_11_off + 1246 * ccomps * dcomps);

            auto g_z_z_xyy_xxxz = cbuffer.data(fg_geom_11_off + 1247 * ccomps * dcomps);

            auto g_z_z_xyy_xxyy = cbuffer.data(fg_geom_11_off + 1248 * ccomps * dcomps);

            auto g_z_z_xyy_xxyz = cbuffer.data(fg_geom_11_off + 1249 * ccomps * dcomps);

            auto g_z_z_xyy_xxzz = cbuffer.data(fg_geom_11_off + 1250 * ccomps * dcomps);

            auto g_z_z_xyy_xyyy = cbuffer.data(fg_geom_11_off + 1251 * ccomps * dcomps);

            auto g_z_z_xyy_xyyz = cbuffer.data(fg_geom_11_off + 1252 * ccomps * dcomps);

            auto g_z_z_xyy_xyzz = cbuffer.data(fg_geom_11_off + 1253 * ccomps * dcomps);

            auto g_z_z_xyy_xzzz = cbuffer.data(fg_geom_11_off + 1254 * ccomps * dcomps);

            auto g_z_z_xyy_yyyy = cbuffer.data(fg_geom_11_off + 1255 * ccomps * dcomps);

            auto g_z_z_xyy_yyyz = cbuffer.data(fg_geom_11_off + 1256 * ccomps * dcomps);

            auto g_z_z_xyy_yyzz = cbuffer.data(fg_geom_11_off + 1257 * ccomps * dcomps);

            auto g_z_z_xyy_yzzz = cbuffer.data(fg_geom_11_off + 1258 * ccomps * dcomps);

            auto g_z_z_xyy_zzzz = cbuffer.data(fg_geom_11_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyy_xxxx, g_z_z_xyy_xxxy, g_z_z_xyy_xxxz, g_z_z_xyy_xxyy, g_z_z_xyy_xxyz, g_z_z_xyy_xxzz, g_z_z_xyy_xyyy, g_z_z_xyy_xyyz, g_z_z_xyy_xyzz, g_z_z_xyy_xzzz, g_z_z_xyy_yyyy, g_z_z_xyy_yyyz, g_z_z_xyy_yyzz, g_z_z_xyy_yzzz, g_z_z_xyy_zzzz, g_z_z_yy_xxxx, g_z_z_yy_xxxxx, g_z_z_yy_xxxxy, g_z_z_yy_xxxxz, g_z_z_yy_xxxy, g_z_z_yy_xxxyy, g_z_z_yy_xxxyz, g_z_z_yy_xxxz, g_z_z_yy_xxxzz, g_z_z_yy_xxyy, g_z_z_yy_xxyyy, g_z_z_yy_xxyyz, g_z_z_yy_xxyz, g_z_z_yy_xxyzz, g_z_z_yy_xxzz, g_z_z_yy_xxzzz, g_z_z_yy_xyyy, g_z_z_yy_xyyyy, g_z_z_yy_xyyyz, g_z_z_yy_xyyz, g_z_z_yy_xyyzz, g_z_z_yy_xyzz, g_z_z_yy_xyzzz, g_z_z_yy_xzzz, g_z_z_yy_xzzzz, g_z_z_yy_yyyy, g_z_z_yy_yyyz, g_z_z_yy_yyzz, g_z_z_yy_yzzz, g_z_z_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyy_xxxx[k] = -g_z_z_yy_xxxx[k] * ab_x + g_z_z_yy_xxxxx[k];

                g_z_z_xyy_xxxy[k] = -g_z_z_yy_xxxy[k] * ab_x + g_z_z_yy_xxxxy[k];

                g_z_z_xyy_xxxz[k] = -g_z_z_yy_xxxz[k] * ab_x + g_z_z_yy_xxxxz[k];

                g_z_z_xyy_xxyy[k] = -g_z_z_yy_xxyy[k] * ab_x + g_z_z_yy_xxxyy[k];

                g_z_z_xyy_xxyz[k] = -g_z_z_yy_xxyz[k] * ab_x + g_z_z_yy_xxxyz[k];

                g_z_z_xyy_xxzz[k] = -g_z_z_yy_xxzz[k] * ab_x + g_z_z_yy_xxxzz[k];

                g_z_z_xyy_xyyy[k] = -g_z_z_yy_xyyy[k] * ab_x + g_z_z_yy_xxyyy[k];

                g_z_z_xyy_xyyz[k] = -g_z_z_yy_xyyz[k] * ab_x + g_z_z_yy_xxyyz[k];

                g_z_z_xyy_xyzz[k] = -g_z_z_yy_xyzz[k] * ab_x + g_z_z_yy_xxyzz[k];

                g_z_z_xyy_xzzz[k] = -g_z_z_yy_xzzz[k] * ab_x + g_z_z_yy_xxzzz[k];

                g_z_z_xyy_yyyy[k] = -g_z_z_yy_yyyy[k] * ab_x + g_z_z_yy_xyyyy[k];

                g_z_z_xyy_yyyz[k] = -g_z_z_yy_yyyz[k] * ab_x + g_z_z_yy_xyyyz[k];

                g_z_z_xyy_yyzz[k] = -g_z_z_yy_yyzz[k] * ab_x + g_z_z_yy_xyyzz[k];

                g_z_z_xyy_yzzz[k] = -g_z_z_yy_yzzz[k] * ab_x + g_z_z_yy_xyzzz[k];

                g_z_z_xyy_zzzz[k] = -g_z_z_yy_zzzz[k] * ab_x + g_z_z_yy_xzzzz[k];
            }

            /// Set up 1260-1275 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyz_xxxx = cbuffer.data(fg_geom_11_off + 1260 * ccomps * dcomps);

            auto g_z_z_xyz_xxxy = cbuffer.data(fg_geom_11_off + 1261 * ccomps * dcomps);

            auto g_z_z_xyz_xxxz = cbuffer.data(fg_geom_11_off + 1262 * ccomps * dcomps);

            auto g_z_z_xyz_xxyy = cbuffer.data(fg_geom_11_off + 1263 * ccomps * dcomps);

            auto g_z_z_xyz_xxyz = cbuffer.data(fg_geom_11_off + 1264 * ccomps * dcomps);

            auto g_z_z_xyz_xxzz = cbuffer.data(fg_geom_11_off + 1265 * ccomps * dcomps);

            auto g_z_z_xyz_xyyy = cbuffer.data(fg_geom_11_off + 1266 * ccomps * dcomps);

            auto g_z_z_xyz_xyyz = cbuffer.data(fg_geom_11_off + 1267 * ccomps * dcomps);

            auto g_z_z_xyz_xyzz = cbuffer.data(fg_geom_11_off + 1268 * ccomps * dcomps);

            auto g_z_z_xyz_xzzz = cbuffer.data(fg_geom_11_off + 1269 * ccomps * dcomps);

            auto g_z_z_xyz_yyyy = cbuffer.data(fg_geom_11_off + 1270 * ccomps * dcomps);

            auto g_z_z_xyz_yyyz = cbuffer.data(fg_geom_11_off + 1271 * ccomps * dcomps);

            auto g_z_z_xyz_yyzz = cbuffer.data(fg_geom_11_off + 1272 * ccomps * dcomps);

            auto g_z_z_xyz_yzzz = cbuffer.data(fg_geom_11_off + 1273 * ccomps * dcomps);

            auto g_z_z_xyz_zzzz = cbuffer.data(fg_geom_11_off + 1274 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyz_xxxx, g_z_z_xyz_xxxy, g_z_z_xyz_xxxz, g_z_z_xyz_xxyy, g_z_z_xyz_xxyz, g_z_z_xyz_xxzz, g_z_z_xyz_xyyy, g_z_z_xyz_xyyz, g_z_z_xyz_xyzz, g_z_z_xyz_xzzz, g_z_z_xyz_yyyy, g_z_z_xyz_yyyz, g_z_z_xyz_yyzz, g_z_z_xyz_yzzz, g_z_z_xyz_zzzz, g_z_z_yz_xxxx, g_z_z_yz_xxxxx, g_z_z_yz_xxxxy, g_z_z_yz_xxxxz, g_z_z_yz_xxxy, g_z_z_yz_xxxyy, g_z_z_yz_xxxyz, g_z_z_yz_xxxz, g_z_z_yz_xxxzz, g_z_z_yz_xxyy, g_z_z_yz_xxyyy, g_z_z_yz_xxyyz, g_z_z_yz_xxyz, g_z_z_yz_xxyzz, g_z_z_yz_xxzz, g_z_z_yz_xxzzz, g_z_z_yz_xyyy, g_z_z_yz_xyyyy, g_z_z_yz_xyyyz, g_z_z_yz_xyyz, g_z_z_yz_xyyzz, g_z_z_yz_xyzz, g_z_z_yz_xyzzz, g_z_z_yz_xzzz, g_z_z_yz_xzzzz, g_z_z_yz_yyyy, g_z_z_yz_yyyz, g_z_z_yz_yyzz, g_z_z_yz_yzzz, g_z_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyz_xxxx[k] = -g_z_z_yz_xxxx[k] * ab_x + g_z_z_yz_xxxxx[k];

                g_z_z_xyz_xxxy[k] = -g_z_z_yz_xxxy[k] * ab_x + g_z_z_yz_xxxxy[k];

                g_z_z_xyz_xxxz[k] = -g_z_z_yz_xxxz[k] * ab_x + g_z_z_yz_xxxxz[k];

                g_z_z_xyz_xxyy[k] = -g_z_z_yz_xxyy[k] * ab_x + g_z_z_yz_xxxyy[k];

                g_z_z_xyz_xxyz[k] = -g_z_z_yz_xxyz[k] * ab_x + g_z_z_yz_xxxyz[k];

                g_z_z_xyz_xxzz[k] = -g_z_z_yz_xxzz[k] * ab_x + g_z_z_yz_xxxzz[k];

                g_z_z_xyz_xyyy[k] = -g_z_z_yz_xyyy[k] * ab_x + g_z_z_yz_xxyyy[k];

                g_z_z_xyz_xyyz[k] = -g_z_z_yz_xyyz[k] * ab_x + g_z_z_yz_xxyyz[k];

                g_z_z_xyz_xyzz[k] = -g_z_z_yz_xyzz[k] * ab_x + g_z_z_yz_xxyzz[k];

                g_z_z_xyz_xzzz[k] = -g_z_z_yz_xzzz[k] * ab_x + g_z_z_yz_xxzzz[k];

                g_z_z_xyz_yyyy[k] = -g_z_z_yz_yyyy[k] * ab_x + g_z_z_yz_xyyyy[k];

                g_z_z_xyz_yyyz[k] = -g_z_z_yz_yyyz[k] * ab_x + g_z_z_yz_xyyyz[k];

                g_z_z_xyz_yyzz[k] = -g_z_z_yz_yyzz[k] * ab_x + g_z_z_yz_xyyzz[k];

                g_z_z_xyz_yzzz[k] = -g_z_z_yz_yzzz[k] * ab_x + g_z_z_yz_xyzzz[k];

                g_z_z_xyz_zzzz[k] = -g_z_z_yz_zzzz[k] * ab_x + g_z_z_yz_xzzzz[k];
            }

            /// Set up 1275-1290 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzz_xxxx = cbuffer.data(fg_geom_11_off + 1275 * ccomps * dcomps);

            auto g_z_z_xzz_xxxy = cbuffer.data(fg_geom_11_off + 1276 * ccomps * dcomps);

            auto g_z_z_xzz_xxxz = cbuffer.data(fg_geom_11_off + 1277 * ccomps * dcomps);

            auto g_z_z_xzz_xxyy = cbuffer.data(fg_geom_11_off + 1278 * ccomps * dcomps);

            auto g_z_z_xzz_xxyz = cbuffer.data(fg_geom_11_off + 1279 * ccomps * dcomps);

            auto g_z_z_xzz_xxzz = cbuffer.data(fg_geom_11_off + 1280 * ccomps * dcomps);

            auto g_z_z_xzz_xyyy = cbuffer.data(fg_geom_11_off + 1281 * ccomps * dcomps);

            auto g_z_z_xzz_xyyz = cbuffer.data(fg_geom_11_off + 1282 * ccomps * dcomps);

            auto g_z_z_xzz_xyzz = cbuffer.data(fg_geom_11_off + 1283 * ccomps * dcomps);

            auto g_z_z_xzz_xzzz = cbuffer.data(fg_geom_11_off + 1284 * ccomps * dcomps);

            auto g_z_z_xzz_yyyy = cbuffer.data(fg_geom_11_off + 1285 * ccomps * dcomps);

            auto g_z_z_xzz_yyyz = cbuffer.data(fg_geom_11_off + 1286 * ccomps * dcomps);

            auto g_z_z_xzz_yyzz = cbuffer.data(fg_geom_11_off + 1287 * ccomps * dcomps);

            auto g_z_z_xzz_yzzz = cbuffer.data(fg_geom_11_off + 1288 * ccomps * dcomps);

            auto g_z_z_xzz_zzzz = cbuffer.data(fg_geom_11_off + 1289 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzz_xxxx, g_z_z_xzz_xxxy, g_z_z_xzz_xxxz, g_z_z_xzz_xxyy, g_z_z_xzz_xxyz, g_z_z_xzz_xxzz, g_z_z_xzz_xyyy, g_z_z_xzz_xyyz, g_z_z_xzz_xyzz, g_z_z_xzz_xzzz, g_z_z_xzz_yyyy, g_z_z_xzz_yyyz, g_z_z_xzz_yyzz, g_z_z_xzz_yzzz, g_z_z_xzz_zzzz, g_z_z_zz_xxxx, g_z_z_zz_xxxxx, g_z_z_zz_xxxxy, g_z_z_zz_xxxxz, g_z_z_zz_xxxy, g_z_z_zz_xxxyy, g_z_z_zz_xxxyz, g_z_z_zz_xxxz, g_z_z_zz_xxxzz, g_z_z_zz_xxyy, g_z_z_zz_xxyyy, g_z_z_zz_xxyyz, g_z_z_zz_xxyz, g_z_z_zz_xxyzz, g_z_z_zz_xxzz, g_z_z_zz_xxzzz, g_z_z_zz_xyyy, g_z_z_zz_xyyyy, g_z_z_zz_xyyyz, g_z_z_zz_xyyz, g_z_z_zz_xyyzz, g_z_z_zz_xyzz, g_z_z_zz_xyzzz, g_z_z_zz_xzzz, g_z_z_zz_xzzzz, g_z_z_zz_yyyy, g_z_z_zz_yyyz, g_z_z_zz_yyzz, g_z_z_zz_yzzz, g_z_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzz_xxxx[k] = -g_z_z_zz_xxxx[k] * ab_x + g_z_z_zz_xxxxx[k];

                g_z_z_xzz_xxxy[k] = -g_z_z_zz_xxxy[k] * ab_x + g_z_z_zz_xxxxy[k];

                g_z_z_xzz_xxxz[k] = -g_z_z_zz_xxxz[k] * ab_x + g_z_z_zz_xxxxz[k];

                g_z_z_xzz_xxyy[k] = -g_z_z_zz_xxyy[k] * ab_x + g_z_z_zz_xxxyy[k];

                g_z_z_xzz_xxyz[k] = -g_z_z_zz_xxyz[k] * ab_x + g_z_z_zz_xxxyz[k];

                g_z_z_xzz_xxzz[k] = -g_z_z_zz_xxzz[k] * ab_x + g_z_z_zz_xxxzz[k];

                g_z_z_xzz_xyyy[k] = -g_z_z_zz_xyyy[k] * ab_x + g_z_z_zz_xxyyy[k];

                g_z_z_xzz_xyyz[k] = -g_z_z_zz_xyyz[k] * ab_x + g_z_z_zz_xxyyz[k];

                g_z_z_xzz_xyzz[k] = -g_z_z_zz_xyzz[k] * ab_x + g_z_z_zz_xxyzz[k];

                g_z_z_xzz_xzzz[k] = -g_z_z_zz_xzzz[k] * ab_x + g_z_z_zz_xxzzz[k];

                g_z_z_xzz_yyyy[k] = -g_z_z_zz_yyyy[k] * ab_x + g_z_z_zz_xyyyy[k];

                g_z_z_xzz_yyyz[k] = -g_z_z_zz_yyyz[k] * ab_x + g_z_z_zz_xyyyz[k];

                g_z_z_xzz_yyzz[k] = -g_z_z_zz_yyzz[k] * ab_x + g_z_z_zz_xyyzz[k];

                g_z_z_xzz_yzzz[k] = -g_z_z_zz_yzzz[k] * ab_x + g_z_z_zz_xyzzz[k];

                g_z_z_xzz_zzzz[k] = -g_z_z_zz_zzzz[k] * ab_x + g_z_z_zz_xzzzz[k];
            }

            /// Set up 1290-1305 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyy_xxxx = cbuffer.data(fg_geom_11_off + 1290 * ccomps * dcomps);

            auto g_z_z_yyy_xxxy = cbuffer.data(fg_geom_11_off + 1291 * ccomps * dcomps);

            auto g_z_z_yyy_xxxz = cbuffer.data(fg_geom_11_off + 1292 * ccomps * dcomps);

            auto g_z_z_yyy_xxyy = cbuffer.data(fg_geom_11_off + 1293 * ccomps * dcomps);

            auto g_z_z_yyy_xxyz = cbuffer.data(fg_geom_11_off + 1294 * ccomps * dcomps);

            auto g_z_z_yyy_xxzz = cbuffer.data(fg_geom_11_off + 1295 * ccomps * dcomps);

            auto g_z_z_yyy_xyyy = cbuffer.data(fg_geom_11_off + 1296 * ccomps * dcomps);

            auto g_z_z_yyy_xyyz = cbuffer.data(fg_geom_11_off + 1297 * ccomps * dcomps);

            auto g_z_z_yyy_xyzz = cbuffer.data(fg_geom_11_off + 1298 * ccomps * dcomps);

            auto g_z_z_yyy_xzzz = cbuffer.data(fg_geom_11_off + 1299 * ccomps * dcomps);

            auto g_z_z_yyy_yyyy = cbuffer.data(fg_geom_11_off + 1300 * ccomps * dcomps);

            auto g_z_z_yyy_yyyz = cbuffer.data(fg_geom_11_off + 1301 * ccomps * dcomps);

            auto g_z_z_yyy_yyzz = cbuffer.data(fg_geom_11_off + 1302 * ccomps * dcomps);

            auto g_z_z_yyy_yzzz = cbuffer.data(fg_geom_11_off + 1303 * ccomps * dcomps);

            auto g_z_z_yyy_zzzz = cbuffer.data(fg_geom_11_off + 1304 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yy_xxxx, g_z_z_yy_xxxxy, g_z_z_yy_xxxy, g_z_z_yy_xxxyy, g_z_z_yy_xxxyz, g_z_z_yy_xxxz, g_z_z_yy_xxyy, g_z_z_yy_xxyyy, g_z_z_yy_xxyyz, g_z_z_yy_xxyz, g_z_z_yy_xxyzz, g_z_z_yy_xxzz, g_z_z_yy_xyyy, g_z_z_yy_xyyyy, g_z_z_yy_xyyyz, g_z_z_yy_xyyz, g_z_z_yy_xyyzz, g_z_z_yy_xyzz, g_z_z_yy_xyzzz, g_z_z_yy_xzzz, g_z_z_yy_yyyy, g_z_z_yy_yyyyy, g_z_z_yy_yyyyz, g_z_z_yy_yyyz, g_z_z_yy_yyyzz, g_z_z_yy_yyzz, g_z_z_yy_yyzzz, g_z_z_yy_yzzz, g_z_z_yy_yzzzz, g_z_z_yy_zzzz, g_z_z_yyy_xxxx, g_z_z_yyy_xxxy, g_z_z_yyy_xxxz, g_z_z_yyy_xxyy, g_z_z_yyy_xxyz, g_z_z_yyy_xxzz, g_z_z_yyy_xyyy, g_z_z_yyy_xyyz, g_z_z_yyy_xyzz, g_z_z_yyy_xzzz, g_z_z_yyy_yyyy, g_z_z_yyy_yyyz, g_z_z_yyy_yyzz, g_z_z_yyy_yzzz, g_z_z_yyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyy_xxxx[k] = -g_z_z_yy_xxxx[k] * ab_y + g_z_z_yy_xxxxy[k];

                g_z_z_yyy_xxxy[k] = -g_z_z_yy_xxxy[k] * ab_y + g_z_z_yy_xxxyy[k];

                g_z_z_yyy_xxxz[k] = -g_z_z_yy_xxxz[k] * ab_y + g_z_z_yy_xxxyz[k];

                g_z_z_yyy_xxyy[k] = -g_z_z_yy_xxyy[k] * ab_y + g_z_z_yy_xxyyy[k];

                g_z_z_yyy_xxyz[k] = -g_z_z_yy_xxyz[k] * ab_y + g_z_z_yy_xxyyz[k];

                g_z_z_yyy_xxzz[k] = -g_z_z_yy_xxzz[k] * ab_y + g_z_z_yy_xxyzz[k];

                g_z_z_yyy_xyyy[k] = -g_z_z_yy_xyyy[k] * ab_y + g_z_z_yy_xyyyy[k];

                g_z_z_yyy_xyyz[k] = -g_z_z_yy_xyyz[k] * ab_y + g_z_z_yy_xyyyz[k];

                g_z_z_yyy_xyzz[k] = -g_z_z_yy_xyzz[k] * ab_y + g_z_z_yy_xyyzz[k];

                g_z_z_yyy_xzzz[k] = -g_z_z_yy_xzzz[k] * ab_y + g_z_z_yy_xyzzz[k];

                g_z_z_yyy_yyyy[k] = -g_z_z_yy_yyyy[k] * ab_y + g_z_z_yy_yyyyy[k];

                g_z_z_yyy_yyyz[k] = -g_z_z_yy_yyyz[k] * ab_y + g_z_z_yy_yyyyz[k];

                g_z_z_yyy_yyzz[k] = -g_z_z_yy_yyzz[k] * ab_y + g_z_z_yy_yyyzz[k];

                g_z_z_yyy_yzzz[k] = -g_z_z_yy_yzzz[k] * ab_y + g_z_z_yy_yyzzz[k];

                g_z_z_yyy_zzzz[k] = -g_z_z_yy_zzzz[k] * ab_y + g_z_z_yy_yzzzz[k];
            }

            /// Set up 1305-1320 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyz_xxxx = cbuffer.data(fg_geom_11_off + 1305 * ccomps * dcomps);

            auto g_z_z_yyz_xxxy = cbuffer.data(fg_geom_11_off + 1306 * ccomps * dcomps);

            auto g_z_z_yyz_xxxz = cbuffer.data(fg_geom_11_off + 1307 * ccomps * dcomps);

            auto g_z_z_yyz_xxyy = cbuffer.data(fg_geom_11_off + 1308 * ccomps * dcomps);

            auto g_z_z_yyz_xxyz = cbuffer.data(fg_geom_11_off + 1309 * ccomps * dcomps);

            auto g_z_z_yyz_xxzz = cbuffer.data(fg_geom_11_off + 1310 * ccomps * dcomps);

            auto g_z_z_yyz_xyyy = cbuffer.data(fg_geom_11_off + 1311 * ccomps * dcomps);

            auto g_z_z_yyz_xyyz = cbuffer.data(fg_geom_11_off + 1312 * ccomps * dcomps);

            auto g_z_z_yyz_xyzz = cbuffer.data(fg_geom_11_off + 1313 * ccomps * dcomps);

            auto g_z_z_yyz_xzzz = cbuffer.data(fg_geom_11_off + 1314 * ccomps * dcomps);

            auto g_z_z_yyz_yyyy = cbuffer.data(fg_geom_11_off + 1315 * ccomps * dcomps);

            auto g_z_z_yyz_yyyz = cbuffer.data(fg_geom_11_off + 1316 * ccomps * dcomps);

            auto g_z_z_yyz_yyzz = cbuffer.data(fg_geom_11_off + 1317 * ccomps * dcomps);

            auto g_z_z_yyz_yzzz = cbuffer.data(fg_geom_11_off + 1318 * ccomps * dcomps);

            auto g_z_z_yyz_zzzz = cbuffer.data(fg_geom_11_off + 1319 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyz_xxxx, g_z_z_yyz_xxxy, g_z_z_yyz_xxxz, g_z_z_yyz_xxyy, g_z_z_yyz_xxyz, g_z_z_yyz_xxzz, g_z_z_yyz_xyyy, g_z_z_yyz_xyyz, g_z_z_yyz_xyzz, g_z_z_yyz_xzzz, g_z_z_yyz_yyyy, g_z_z_yyz_yyyz, g_z_z_yyz_yyzz, g_z_z_yyz_yzzz, g_z_z_yyz_zzzz, g_z_z_yz_xxxx, g_z_z_yz_xxxxy, g_z_z_yz_xxxy, g_z_z_yz_xxxyy, g_z_z_yz_xxxyz, g_z_z_yz_xxxz, g_z_z_yz_xxyy, g_z_z_yz_xxyyy, g_z_z_yz_xxyyz, g_z_z_yz_xxyz, g_z_z_yz_xxyzz, g_z_z_yz_xxzz, g_z_z_yz_xyyy, g_z_z_yz_xyyyy, g_z_z_yz_xyyyz, g_z_z_yz_xyyz, g_z_z_yz_xyyzz, g_z_z_yz_xyzz, g_z_z_yz_xyzzz, g_z_z_yz_xzzz, g_z_z_yz_yyyy, g_z_z_yz_yyyyy, g_z_z_yz_yyyyz, g_z_z_yz_yyyz, g_z_z_yz_yyyzz, g_z_z_yz_yyzz, g_z_z_yz_yyzzz, g_z_z_yz_yzzz, g_z_z_yz_yzzzz, g_z_z_yz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyz_xxxx[k] = -g_z_z_yz_xxxx[k] * ab_y + g_z_z_yz_xxxxy[k];

                g_z_z_yyz_xxxy[k] = -g_z_z_yz_xxxy[k] * ab_y + g_z_z_yz_xxxyy[k];

                g_z_z_yyz_xxxz[k] = -g_z_z_yz_xxxz[k] * ab_y + g_z_z_yz_xxxyz[k];

                g_z_z_yyz_xxyy[k] = -g_z_z_yz_xxyy[k] * ab_y + g_z_z_yz_xxyyy[k];

                g_z_z_yyz_xxyz[k] = -g_z_z_yz_xxyz[k] * ab_y + g_z_z_yz_xxyyz[k];

                g_z_z_yyz_xxzz[k] = -g_z_z_yz_xxzz[k] * ab_y + g_z_z_yz_xxyzz[k];

                g_z_z_yyz_xyyy[k] = -g_z_z_yz_xyyy[k] * ab_y + g_z_z_yz_xyyyy[k];

                g_z_z_yyz_xyyz[k] = -g_z_z_yz_xyyz[k] * ab_y + g_z_z_yz_xyyyz[k];

                g_z_z_yyz_xyzz[k] = -g_z_z_yz_xyzz[k] * ab_y + g_z_z_yz_xyyzz[k];

                g_z_z_yyz_xzzz[k] = -g_z_z_yz_xzzz[k] * ab_y + g_z_z_yz_xyzzz[k];

                g_z_z_yyz_yyyy[k] = -g_z_z_yz_yyyy[k] * ab_y + g_z_z_yz_yyyyy[k];

                g_z_z_yyz_yyyz[k] = -g_z_z_yz_yyyz[k] * ab_y + g_z_z_yz_yyyyz[k];

                g_z_z_yyz_yyzz[k] = -g_z_z_yz_yyzz[k] * ab_y + g_z_z_yz_yyyzz[k];

                g_z_z_yyz_yzzz[k] = -g_z_z_yz_yzzz[k] * ab_y + g_z_z_yz_yyzzz[k];

                g_z_z_yyz_zzzz[k] = -g_z_z_yz_zzzz[k] * ab_y + g_z_z_yz_yzzzz[k];
            }

            /// Set up 1320-1335 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzz_xxxx = cbuffer.data(fg_geom_11_off + 1320 * ccomps * dcomps);

            auto g_z_z_yzz_xxxy = cbuffer.data(fg_geom_11_off + 1321 * ccomps * dcomps);

            auto g_z_z_yzz_xxxz = cbuffer.data(fg_geom_11_off + 1322 * ccomps * dcomps);

            auto g_z_z_yzz_xxyy = cbuffer.data(fg_geom_11_off + 1323 * ccomps * dcomps);

            auto g_z_z_yzz_xxyz = cbuffer.data(fg_geom_11_off + 1324 * ccomps * dcomps);

            auto g_z_z_yzz_xxzz = cbuffer.data(fg_geom_11_off + 1325 * ccomps * dcomps);

            auto g_z_z_yzz_xyyy = cbuffer.data(fg_geom_11_off + 1326 * ccomps * dcomps);

            auto g_z_z_yzz_xyyz = cbuffer.data(fg_geom_11_off + 1327 * ccomps * dcomps);

            auto g_z_z_yzz_xyzz = cbuffer.data(fg_geom_11_off + 1328 * ccomps * dcomps);

            auto g_z_z_yzz_xzzz = cbuffer.data(fg_geom_11_off + 1329 * ccomps * dcomps);

            auto g_z_z_yzz_yyyy = cbuffer.data(fg_geom_11_off + 1330 * ccomps * dcomps);

            auto g_z_z_yzz_yyyz = cbuffer.data(fg_geom_11_off + 1331 * ccomps * dcomps);

            auto g_z_z_yzz_yyzz = cbuffer.data(fg_geom_11_off + 1332 * ccomps * dcomps);

            auto g_z_z_yzz_yzzz = cbuffer.data(fg_geom_11_off + 1333 * ccomps * dcomps);

            auto g_z_z_yzz_zzzz = cbuffer.data(fg_geom_11_off + 1334 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzz_xxxx, g_z_z_yzz_xxxy, g_z_z_yzz_xxxz, g_z_z_yzz_xxyy, g_z_z_yzz_xxyz, g_z_z_yzz_xxzz, g_z_z_yzz_xyyy, g_z_z_yzz_xyyz, g_z_z_yzz_xyzz, g_z_z_yzz_xzzz, g_z_z_yzz_yyyy, g_z_z_yzz_yyyz, g_z_z_yzz_yyzz, g_z_z_yzz_yzzz, g_z_z_yzz_zzzz, g_z_z_zz_xxxx, g_z_z_zz_xxxxy, g_z_z_zz_xxxy, g_z_z_zz_xxxyy, g_z_z_zz_xxxyz, g_z_z_zz_xxxz, g_z_z_zz_xxyy, g_z_z_zz_xxyyy, g_z_z_zz_xxyyz, g_z_z_zz_xxyz, g_z_z_zz_xxyzz, g_z_z_zz_xxzz, g_z_z_zz_xyyy, g_z_z_zz_xyyyy, g_z_z_zz_xyyyz, g_z_z_zz_xyyz, g_z_z_zz_xyyzz, g_z_z_zz_xyzz, g_z_z_zz_xyzzz, g_z_z_zz_xzzz, g_z_z_zz_yyyy, g_z_z_zz_yyyyy, g_z_z_zz_yyyyz, g_z_z_zz_yyyz, g_z_z_zz_yyyzz, g_z_z_zz_yyzz, g_z_z_zz_yyzzz, g_z_z_zz_yzzz, g_z_z_zz_yzzzz, g_z_z_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzz_xxxx[k] = -g_z_z_zz_xxxx[k] * ab_y + g_z_z_zz_xxxxy[k];

                g_z_z_yzz_xxxy[k] = -g_z_z_zz_xxxy[k] * ab_y + g_z_z_zz_xxxyy[k];

                g_z_z_yzz_xxxz[k] = -g_z_z_zz_xxxz[k] * ab_y + g_z_z_zz_xxxyz[k];

                g_z_z_yzz_xxyy[k] = -g_z_z_zz_xxyy[k] * ab_y + g_z_z_zz_xxyyy[k];

                g_z_z_yzz_xxyz[k] = -g_z_z_zz_xxyz[k] * ab_y + g_z_z_zz_xxyyz[k];

                g_z_z_yzz_xxzz[k] = -g_z_z_zz_xxzz[k] * ab_y + g_z_z_zz_xxyzz[k];

                g_z_z_yzz_xyyy[k] = -g_z_z_zz_xyyy[k] * ab_y + g_z_z_zz_xyyyy[k];

                g_z_z_yzz_xyyz[k] = -g_z_z_zz_xyyz[k] * ab_y + g_z_z_zz_xyyyz[k];

                g_z_z_yzz_xyzz[k] = -g_z_z_zz_xyzz[k] * ab_y + g_z_z_zz_xyyzz[k];

                g_z_z_yzz_xzzz[k] = -g_z_z_zz_xzzz[k] * ab_y + g_z_z_zz_xyzzz[k];

                g_z_z_yzz_yyyy[k] = -g_z_z_zz_yyyy[k] * ab_y + g_z_z_zz_yyyyy[k];

                g_z_z_yzz_yyyz[k] = -g_z_z_zz_yyyz[k] * ab_y + g_z_z_zz_yyyyz[k];

                g_z_z_yzz_yyzz[k] = -g_z_z_zz_yyzz[k] * ab_y + g_z_z_zz_yyyzz[k];

                g_z_z_yzz_yzzz[k] = -g_z_z_zz_yzzz[k] * ab_y + g_z_z_zz_yyzzz[k];

                g_z_z_yzz_zzzz[k] = -g_z_z_zz_zzzz[k] * ab_y + g_z_z_zz_yzzzz[k];
            }

            /// Set up 1335-1350 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzz_xxxx = cbuffer.data(fg_geom_11_off + 1335 * ccomps * dcomps);

            auto g_z_z_zzz_xxxy = cbuffer.data(fg_geom_11_off + 1336 * ccomps * dcomps);

            auto g_z_z_zzz_xxxz = cbuffer.data(fg_geom_11_off + 1337 * ccomps * dcomps);

            auto g_z_z_zzz_xxyy = cbuffer.data(fg_geom_11_off + 1338 * ccomps * dcomps);

            auto g_z_z_zzz_xxyz = cbuffer.data(fg_geom_11_off + 1339 * ccomps * dcomps);

            auto g_z_z_zzz_xxzz = cbuffer.data(fg_geom_11_off + 1340 * ccomps * dcomps);

            auto g_z_z_zzz_xyyy = cbuffer.data(fg_geom_11_off + 1341 * ccomps * dcomps);

            auto g_z_z_zzz_xyyz = cbuffer.data(fg_geom_11_off + 1342 * ccomps * dcomps);

            auto g_z_z_zzz_xyzz = cbuffer.data(fg_geom_11_off + 1343 * ccomps * dcomps);

            auto g_z_z_zzz_xzzz = cbuffer.data(fg_geom_11_off + 1344 * ccomps * dcomps);

            auto g_z_z_zzz_yyyy = cbuffer.data(fg_geom_11_off + 1345 * ccomps * dcomps);

            auto g_z_z_zzz_yyyz = cbuffer.data(fg_geom_11_off + 1346 * ccomps * dcomps);

            auto g_z_z_zzz_yyzz = cbuffer.data(fg_geom_11_off + 1347 * ccomps * dcomps);

            auto g_z_z_zzz_yzzz = cbuffer.data(fg_geom_11_off + 1348 * ccomps * dcomps);

            auto g_z_z_zzz_zzzz = cbuffer.data(fg_geom_11_off + 1349 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxxx, g_0_z_zz_xxxy, g_0_z_zz_xxxz, g_0_z_zz_xxyy, g_0_z_zz_xxyz, g_0_z_zz_xxzz, g_0_z_zz_xyyy, g_0_z_zz_xyyz, g_0_z_zz_xyzz, g_0_z_zz_xzzz, g_0_z_zz_yyyy, g_0_z_zz_yyyz, g_0_z_zz_yyzz, g_0_z_zz_yzzz, g_0_z_zz_zzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_z_z_zz_xxxx, g_z_z_zz_xxxxz, g_z_z_zz_xxxy, g_z_z_zz_xxxyz, g_z_z_zz_xxxz, g_z_z_zz_xxxzz, g_z_z_zz_xxyy, g_z_z_zz_xxyyz, g_z_z_zz_xxyz, g_z_z_zz_xxyzz, g_z_z_zz_xxzz, g_z_z_zz_xxzzz, g_z_z_zz_xyyy, g_z_z_zz_xyyyz, g_z_z_zz_xyyz, g_z_z_zz_xyyzz, g_z_z_zz_xyzz, g_z_z_zz_xyzzz, g_z_z_zz_xzzz, g_z_z_zz_xzzzz, g_z_z_zz_yyyy, g_z_z_zz_yyyyz, g_z_z_zz_yyyz, g_z_z_zz_yyyzz, g_z_z_zz_yyzz, g_z_z_zz_yyzzz, g_z_z_zz_yzzz, g_z_z_zz_yzzzz, g_z_z_zz_zzzz, g_z_z_zz_zzzzz, g_z_z_zzz_xxxx, g_z_z_zzz_xxxy, g_z_z_zzz_xxxz, g_z_z_zzz_xxyy, g_z_z_zzz_xxyz, g_z_z_zzz_xxzz, g_z_z_zzz_xyyy, g_z_z_zzz_xyyz, g_z_z_zzz_xyzz, g_z_z_zzz_xzzz, g_z_z_zzz_yyyy, g_z_z_zzz_yyyz, g_z_z_zzz_yyzz, g_z_z_zzz_yzzz, g_z_z_zzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzz_xxxx[k] = -g_0_z_zz_xxxx[k] + g_z_0_zz_xxxx[k] - g_z_z_zz_xxxx[k] * ab_z + g_z_z_zz_xxxxz[k];

                g_z_z_zzz_xxxy[k] = -g_0_z_zz_xxxy[k] + g_z_0_zz_xxxy[k] - g_z_z_zz_xxxy[k] * ab_z + g_z_z_zz_xxxyz[k];

                g_z_z_zzz_xxxz[k] = -g_0_z_zz_xxxz[k] + g_z_0_zz_xxxz[k] - g_z_z_zz_xxxz[k] * ab_z + g_z_z_zz_xxxzz[k];

                g_z_z_zzz_xxyy[k] = -g_0_z_zz_xxyy[k] + g_z_0_zz_xxyy[k] - g_z_z_zz_xxyy[k] * ab_z + g_z_z_zz_xxyyz[k];

                g_z_z_zzz_xxyz[k] = -g_0_z_zz_xxyz[k] + g_z_0_zz_xxyz[k] - g_z_z_zz_xxyz[k] * ab_z + g_z_z_zz_xxyzz[k];

                g_z_z_zzz_xxzz[k] = -g_0_z_zz_xxzz[k] + g_z_0_zz_xxzz[k] - g_z_z_zz_xxzz[k] * ab_z + g_z_z_zz_xxzzz[k];

                g_z_z_zzz_xyyy[k] = -g_0_z_zz_xyyy[k] + g_z_0_zz_xyyy[k] - g_z_z_zz_xyyy[k] * ab_z + g_z_z_zz_xyyyz[k];

                g_z_z_zzz_xyyz[k] = -g_0_z_zz_xyyz[k] + g_z_0_zz_xyyz[k] - g_z_z_zz_xyyz[k] * ab_z + g_z_z_zz_xyyzz[k];

                g_z_z_zzz_xyzz[k] = -g_0_z_zz_xyzz[k] + g_z_0_zz_xyzz[k] - g_z_z_zz_xyzz[k] * ab_z + g_z_z_zz_xyzzz[k];

                g_z_z_zzz_xzzz[k] = -g_0_z_zz_xzzz[k] + g_z_0_zz_xzzz[k] - g_z_z_zz_xzzz[k] * ab_z + g_z_z_zz_xzzzz[k];

                g_z_z_zzz_yyyy[k] = -g_0_z_zz_yyyy[k] + g_z_0_zz_yyyy[k] - g_z_z_zz_yyyy[k] * ab_z + g_z_z_zz_yyyyz[k];

                g_z_z_zzz_yyyz[k] = -g_0_z_zz_yyyz[k] + g_z_0_zz_yyyz[k] - g_z_z_zz_yyyz[k] * ab_z + g_z_z_zz_yyyzz[k];

                g_z_z_zzz_yyzz[k] = -g_0_z_zz_yyzz[k] + g_z_0_zz_yyzz[k] - g_z_z_zz_yyzz[k] * ab_z + g_z_z_zz_yyzzz[k];

                g_z_z_zzz_yzzz[k] = -g_0_z_zz_yzzz[k] + g_z_0_zz_yzzz[k] - g_z_z_zz_yzzz[k] * ab_z + g_z_z_zz_yzzzz[k];

                g_z_z_zzz_zzzz[k] = -g_0_z_zz_zzzz[k] + g_z_0_zz_zzzz[k] - g_z_z_zz_zzzz[k] * ab_z + g_z_z_zz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

