#include "ElectronRepulsionGeom2000ContrRecHDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_hdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_hdxx,
                                            const size_t idx_geom_10_gdxx,
                                            const size_t idx_geom_20_gdxx,
                                            const size_t idx_geom_20_gfxx,
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
            /// Set up components of auxilary buffer : GDSS

            const auto gd_geom_10_off = idx_geom_10_gdxx + i * dcomps + j;

            auto g_x_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_xxxx_xx = cbuffer.data(gd_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_xxxx_xy = cbuffer.data(gd_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_xxxx_xz = cbuffer.data(gd_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_xxxx_yy = cbuffer.data(gd_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_xxxx_yz = cbuffer.data(gd_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_xxxx_zz = cbuffer.data(gd_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_xxxy_xx = cbuffer.data(gd_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_xxxy_xy = cbuffer.data(gd_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_xxxy_xz = cbuffer.data(gd_geom_10_off + 188 * ccomps * dcomps);

            auto g_z_0_xxxy_yy = cbuffer.data(gd_geom_10_off + 189 * ccomps * dcomps);

            auto g_z_0_xxxy_yz = cbuffer.data(gd_geom_10_off + 190 * ccomps * dcomps);

            auto g_z_0_xxxy_zz = cbuffer.data(gd_geom_10_off + 191 * ccomps * dcomps);

            auto g_z_0_xxxz_xx = cbuffer.data(gd_geom_10_off + 192 * ccomps * dcomps);

            auto g_z_0_xxxz_xy = cbuffer.data(gd_geom_10_off + 193 * ccomps * dcomps);

            auto g_z_0_xxxz_xz = cbuffer.data(gd_geom_10_off + 194 * ccomps * dcomps);

            auto g_z_0_xxxz_yy = cbuffer.data(gd_geom_10_off + 195 * ccomps * dcomps);

            auto g_z_0_xxxz_yz = cbuffer.data(gd_geom_10_off + 196 * ccomps * dcomps);

            auto g_z_0_xxxz_zz = cbuffer.data(gd_geom_10_off + 197 * ccomps * dcomps);

            auto g_z_0_xxyy_xx = cbuffer.data(gd_geom_10_off + 198 * ccomps * dcomps);

            auto g_z_0_xxyy_xy = cbuffer.data(gd_geom_10_off + 199 * ccomps * dcomps);

            auto g_z_0_xxyy_xz = cbuffer.data(gd_geom_10_off + 200 * ccomps * dcomps);

            auto g_z_0_xxyy_yy = cbuffer.data(gd_geom_10_off + 201 * ccomps * dcomps);

            auto g_z_0_xxyy_yz = cbuffer.data(gd_geom_10_off + 202 * ccomps * dcomps);

            auto g_z_0_xxyy_zz = cbuffer.data(gd_geom_10_off + 203 * ccomps * dcomps);

            auto g_z_0_xxyz_xx = cbuffer.data(gd_geom_10_off + 204 * ccomps * dcomps);

            auto g_z_0_xxyz_xy = cbuffer.data(gd_geom_10_off + 205 * ccomps * dcomps);

            auto g_z_0_xxyz_xz = cbuffer.data(gd_geom_10_off + 206 * ccomps * dcomps);

            auto g_z_0_xxyz_yy = cbuffer.data(gd_geom_10_off + 207 * ccomps * dcomps);

            auto g_z_0_xxyz_yz = cbuffer.data(gd_geom_10_off + 208 * ccomps * dcomps);

            auto g_z_0_xxyz_zz = cbuffer.data(gd_geom_10_off + 209 * ccomps * dcomps);

            auto g_z_0_xxzz_xx = cbuffer.data(gd_geom_10_off + 210 * ccomps * dcomps);

            auto g_z_0_xxzz_xy = cbuffer.data(gd_geom_10_off + 211 * ccomps * dcomps);

            auto g_z_0_xxzz_xz = cbuffer.data(gd_geom_10_off + 212 * ccomps * dcomps);

            auto g_z_0_xxzz_yy = cbuffer.data(gd_geom_10_off + 213 * ccomps * dcomps);

            auto g_z_0_xxzz_yz = cbuffer.data(gd_geom_10_off + 214 * ccomps * dcomps);

            auto g_z_0_xxzz_zz = cbuffer.data(gd_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_xyyy_xx = cbuffer.data(gd_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_xyyy_xy = cbuffer.data(gd_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_xyyy_xz = cbuffer.data(gd_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_xyyy_yy = cbuffer.data(gd_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_xyyy_yz = cbuffer.data(gd_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_xyyy_zz = cbuffer.data(gd_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_xyyz_xx = cbuffer.data(gd_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_xyyz_xy = cbuffer.data(gd_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_xyyz_xz = cbuffer.data(gd_geom_10_off + 224 * ccomps * dcomps);

            auto g_z_0_xyyz_yy = cbuffer.data(gd_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_xyyz_yz = cbuffer.data(gd_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_xyyz_zz = cbuffer.data(gd_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_xyzz_xx = cbuffer.data(gd_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_xyzz_xy = cbuffer.data(gd_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_xyzz_xz = cbuffer.data(gd_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_xyzz_yy = cbuffer.data(gd_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_xyzz_yz = cbuffer.data(gd_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_xyzz_zz = cbuffer.data(gd_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_xzzz_xx = cbuffer.data(gd_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_xzzz_xy = cbuffer.data(gd_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_xzzz_xz = cbuffer.data(gd_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_xzzz_yy = cbuffer.data(gd_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_xzzz_yz = cbuffer.data(gd_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_xzzz_zz = cbuffer.data(gd_geom_10_off + 239 * ccomps * dcomps);

            auto g_z_0_yyyy_xx = cbuffer.data(gd_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_yyyy_xy = cbuffer.data(gd_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_yyyy_xz = cbuffer.data(gd_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_yyyy_yy = cbuffer.data(gd_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_yyyy_yz = cbuffer.data(gd_geom_10_off + 244 * ccomps * dcomps);

            auto g_z_0_yyyy_zz = cbuffer.data(gd_geom_10_off + 245 * ccomps * dcomps);

            auto g_z_0_yyyz_xx = cbuffer.data(gd_geom_10_off + 246 * ccomps * dcomps);

            auto g_z_0_yyyz_xy = cbuffer.data(gd_geom_10_off + 247 * ccomps * dcomps);

            auto g_z_0_yyyz_xz = cbuffer.data(gd_geom_10_off + 248 * ccomps * dcomps);

            auto g_z_0_yyyz_yy = cbuffer.data(gd_geom_10_off + 249 * ccomps * dcomps);

            auto g_z_0_yyyz_yz = cbuffer.data(gd_geom_10_off + 250 * ccomps * dcomps);

            auto g_z_0_yyyz_zz = cbuffer.data(gd_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_yyzz_xx = cbuffer.data(gd_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_yyzz_xy = cbuffer.data(gd_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_yyzz_xz = cbuffer.data(gd_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_yyzz_yy = cbuffer.data(gd_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_yyzz_yz = cbuffer.data(gd_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_yyzz_zz = cbuffer.data(gd_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_yzzz_xx = cbuffer.data(gd_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_yzzz_xy = cbuffer.data(gd_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_yzzz_xz = cbuffer.data(gd_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_yzzz_yy = cbuffer.data(gd_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_yzzz_yz = cbuffer.data(gd_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_yzzz_zz = cbuffer.data(gd_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_zzzz_xx = cbuffer.data(gd_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_zzzz_xy = cbuffer.data(gd_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_zzzz_xz = cbuffer.data(gd_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_zzzz_yy = cbuffer.data(gd_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_zzzz_yz = cbuffer.data(gd_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_zzzz_zz = cbuffer.data(gd_geom_10_off + 269 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GDSS

            const auto gd_geom_20_off = idx_geom_20_gdxx + i * dcomps + j;

            auto g_xx_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 191 * ccomps * dcomps);

            auto g_xz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 197 * ccomps * dcomps);

            auto g_xz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 269 * ccomps * dcomps);

            auto g_yy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 275 * ccomps * dcomps);

            auto g_yy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 281 * ccomps * dcomps);

            auto g_yy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 284 * ccomps * dcomps);

            auto g_yy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 287 * ccomps * dcomps);

            auto g_yy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 293 * ccomps * dcomps);

            auto g_yy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 359 * ccomps * dcomps);

            auto g_yz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 365 * ccomps * dcomps);

            auto g_yz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 371 * ccomps * dcomps);

            auto g_yz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 374 * ccomps * dcomps);

            auto g_yz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 377 * ccomps * dcomps);

            auto g_yz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 383 * ccomps * dcomps);

            auto g_yz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 389 * ccomps * dcomps);

            auto g_yz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 395 * ccomps * dcomps);

            auto g_yz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 419 * ccomps * dcomps);

            auto g_yz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 449 * ccomps * dcomps);

            auto g_zz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 455 * ccomps * dcomps);

            auto g_zz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 461 * ccomps * dcomps);

            auto g_zz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 464 * ccomps * dcomps);

            auto g_zz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 467 * ccomps * dcomps);

            auto g_zz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 473 * ccomps * dcomps);

            auto g_zz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 479 * ccomps * dcomps);

            auto g_zz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 485 * ccomps * dcomps);

            auto g_zz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 491 * ccomps * dcomps);

            auto g_zz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 494 * ccomps * dcomps);

            auto g_zz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 497 * ccomps * dcomps);

            auto g_zz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 509 * ccomps * dcomps);

            auto g_zz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 539 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GFSS

            const auto gf_geom_20_off = idx_geom_20_gfxx + i * dcomps + j;

            auto g_xx_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 139 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 251 * ccomps * dcomps);

            auto g_xy_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 269 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 279 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 377 * ccomps * dcomps);

            auto g_xz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 379 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 383 * ccomps * dcomps);

            auto g_xz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 389 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 391 * ccomps * dcomps);

            auto g_xz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 395 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 399 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 401 * ccomps * dcomps);

            auto g_xz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 404 * ccomps * dcomps);

            auto g_xz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 407 * ccomps * dcomps);

            auto g_xz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 409 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 413 * ccomps * dcomps);

            auto g_xz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 419 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 447 * ccomps * dcomps);

            auto g_xz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 503 * ccomps * dcomps);

            auto g_yy_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 509 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 515 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 519 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 521 * ccomps * dcomps);

            auto g_yy_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 524 * ccomps * dcomps);

            auto g_yy_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 527 * ccomps * dcomps);

            auto g_yy_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 529 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 531 * ccomps * dcomps);

            auto g_yy_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 533 * ccomps * dcomps);

            auto g_yy_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 539 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 545 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 549 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 551 * ccomps * dcomps);

            auto g_yy_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 554 * ccomps * dcomps);

            auto g_yy_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 557 * ccomps * dcomps);

            auto g_yy_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 559 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 563 * ccomps * dcomps);

            auto g_yy_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 569 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 575 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 579 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 581 * ccomps * dcomps);

            auto g_yy_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 584 * ccomps * dcomps);

            auto g_yy_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 587 * ccomps * dcomps);

            auto g_yy_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 589 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 593 * ccomps * dcomps);

            auto g_yy_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 629 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 630 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 631 * ccomps * dcomps);

            auto g_yz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 632 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 633 * ccomps * dcomps);

            auto g_yz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 634 * ccomps * dcomps);

            auto g_yz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 635 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 636 * ccomps * dcomps);

            auto g_yz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 637 * ccomps * dcomps);

            auto g_yz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 638 * ccomps * dcomps);

            auto g_yz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 639 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 640 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 641 * ccomps * dcomps);

            auto g_yz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 642 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 643 * ccomps * dcomps);

            auto g_yz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 644 * ccomps * dcomps);

            auto g_yz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 645 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 646 * ccomps * dcomps);

            auto g_yz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 647 * ccomps * dcomps);

            auto g_yz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 648 * ccomps * dcomps);

            auto g_yz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 649 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 650 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 651 * ccomps * dcomps);

            auto g_yz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 652 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 653 * ccomps * dcomps);

            auto g_yz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 654 * ccomps * dcomps);

            auto g_yz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 655 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 656 * ccomps * dcomps);

            auto g_yz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 657 * ccomps * dcomps);

            auto g_yz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 658 * ccomps * dcomps);

            auto g_yz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 659 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 660 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 661 * ccomps * dcomps);

            auto g_yz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 662 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 663 * ccomps * dcomps);

            auto g_yz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 664 * ccomps * dcomps);

            auto g_yz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 665 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 666 * ccomps * dcomps);

            auto g_yz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 667 * ccomps * dcomps);

            auto g_yz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 668 * ccomps * dcomps);

            auto g_yz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 669 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 670 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 671 * ccomps * dcomps);

            auto g_yz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 674 * ccomps * dcomps);

            auto g_yz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 677 * ccomps * dcomps);

            auto g_yz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 679 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 683 * ccomps * dcomps);

            auto g_yz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 689 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 695 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 699 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 701 * ccomps * dcomps);

            auto g_yz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 704 * ccomps * dcomps);

            auto g_yz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 707 * ccomps * dcomps);

            auto g_yz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 709 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 713 * ccomps * dcomps);

            auto g_yz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 719 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 725 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 727 * ccomps * dcomps);

            auto g_yz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 729 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 731 * ccomps * dcomps);

            auto g_yz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 734 * ccomps * dcomps);

            auto g_yz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 737 * ccomps * dcomps);

            auto g_yz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 739 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 743 * ccomps * dcomps);

            auto g_yz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxx = cbuffer.data(gf_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxy = cbuffer.data(gf_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_xxxx_xxz = cbuffer.data(gf_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyy = cbuffer.data(gf_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_xxxx_xyz = cbuffer.data(gf_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_xxxx_xzz = cbuffer.data(gf_geom_20_off + 755 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyy = cbuffer.data(gf_geom_20_off + 756 * ccomps * dcomps);

            auto g_zz_0_xxxx_yyz = cbuffer.data(gf_geom_20_off + 757 * ccomps * dcomps);

            auto g_zz_0_xxxx_yzz = cbuffer.data(gf_geom_20_off + 758 * ccomps * dcomps);

            auto g_zz_0_xxxx_zzz = cbuffer.data(gf_geom_20_off + 759 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxx = cbuffer.data(gf_geom_20_off + 760 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxy = cbuffer.data(gf_geom_20_off + 761 * ccomps * dcomps);

            auto g_zz_0_xxxy_xxz = cbuffer.data(gf_geom_20_off + 762 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyy = cbuffer.data(gf_geom_20_off + 763 * ccomps * dcomps);

            auto g_zz_0_xxxy_xyz = cbuffer.data(gf_geom_20_off + 764 * ccomps * dcomps);

            auto g_zz_0_xxxy_xzz = cbuffer.data(gf_geom_20_off + 765 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyy = cbuffer.data(gf_geom_20_off + 766 * ccomps * dcomps);

            auto g_zz_0_xxxy_yyz = cbuffer.data(gf_geom_20_off + 767 * ccomps * dcomps);

            auto g_zz_0_xxxy_yzz = cbuffer.data(gf_geom_20_off + 768 * ccomps * dcomps);

            auto g_zz_0_xxxy_zzz = cbuffer.data(gf_geom_20_off + 769 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxx = cbuffer.data(gf_geom_20_off + 770 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxy = cbuffer.data(gf_geom_20_off + 771 * ccomps * dcomps);

            auto g_zz_0_xxxz_xxz = cbuffer.data(gf_geom_20_off + 772 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyy = cbuffer.data(gf_geom_20_off + 773 * ccomps * dcomps);

            auto g_zz_0_xxxz_xyz = cbuffer.data(gf_geom_20_off + 774 * ccomps * dcomps);

            auto g_zz_0_xxxz_xzz = cbuffer.data(gf_geom_20_off + 775 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyy = cbuffer.data(gf_geom_20_off + 776 * ccomps * dcomps);

            auto g_zz_0_xxxz_yyz = cbuffer.data(gf_geom_20_off + 777 * ccomps * dcomps);

            auto g_zz_0_xxxz_yzz = cbuffer.data(gf_geom_20_off + 778 * ccomps * dcomps);

            auto g_zz_0_xxxz_zzz = cbuffer.data(gf_geom_20_off + 779 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxx = cbuffer.data(gf_geom_20_off + 780 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxy = cbuffer.data(gf_geom_20_off + 781 * ccomps * dcomps);

            auto g_zz_0_xxyy_xxz = cbuffer.data(gf_geom_20_off + 782 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyy = cbuffer.data(gf_geom_20_off + 783 * ccomps * dcomps);

            auto g_zz_0_xxyy_xyz = cbuffer.data(gf_geom_20_off + 784 * ccomps * dcomps);

            auto g_zz_0_xxyy_xzz = cbuffer.data(gf_geom_20_off + 785 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyy = cbuffer.data(gf_geom_20_off + 786 * ccomps * dcomps);

            auto g_zz_0_xxyy_yyz = cbuffer.data(gf_geom_20_off + 787 * ccomps * dcomps);

            auto g_zz_0_xxyy_yzz = cbuffer.data(gf_geom_20_off + 788 * ccomps * dcomps);

            auto g_zz_0_xxyy_zzz = cbuffer.data(gf_geom_20_off + 789 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxx = cbuffer.data(gf_geom_20_off + 790 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxy = cbuffer.data(gf_geom_20_off + 791 * ccomps * dcomps);

            auto g_zz_0_xxyz_xxz = cbuffer.data(gf_geom_20_off + 792 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyy = cbuffer.data(gf_geom_20_off + 793 * ccomps * dcomps);

            auto g_zz_0_xxyz_xyz = cbuffer.data(gf_geom_20_off + 794 * ccomps * dcomps);

            auto g_zz_0_xxyz_xzz = cbuffer.data(gf_geom_20_off + 795 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyy = cbuffer.data(gf_geom_20_off + 796 * ccomps * dcomps);

            auto g_zz_0_xxyz_yyz = cbuffer.data(gf_geom_20_off + 797 * ccomps * dcomps);

            auto g_zz_0_xxyz_yzz = cbuffer.data(gf_geom_20_off + 798 * ccomps * dcomps);

            auto g_zz_0_xxyz_zzz = cbuffer.data(gf_geom_20_off + 799 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxx = cbuffer.data(gf_geom_20_off + 800 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxy = cbuffer.data(gf_geom_20_off + 801 * ccomps * dcomps);

            auto g_zz_0_xxzz_xxz = cbuffer.data(gf_geom_20_off + 802 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyy = cbuffer.data(gf_geom_20_off + 803 * ccomps * dcomps);

            auto g_zz_0_xxzz_xyz = cbuffer.data(gf_geom_20_off + 804 * ccomps * dcomps);

            auto g_zz_0_xxzz_xzz = cbuffer.data(gf_geom_20_off + 805 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyy = cbuffer.data(gf_geom_20_off + 806 * ccomps * dcomps);

            auto g_zz_0_xxzz_yyz = cbuffer.data(gf_geom_20_off + 807 * ccomps * dcomps);

            auto g_zz_0_xxzz_yzz = cbuffer.data(gf_geom_20_off + 808 * ccomps * dcomps);

            auto g_zz_0_xxzz_zzz = cbuffer.data(gf_geom_20_off + 809 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxx = cbuffer.data(gf_geom_20_off + 810 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxy = cbuffer.data(gf_geom_20_off + 811 * ccomps * dcomps);

            auto g_zz_0_xyyy_xxz = cbuffer.data(gf_geom_20_off + 812 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyy = cbuffer.data(gf_geom_20_off + 813 * ccomps * dcomps);

            auto g_zz_0_xyyy_xyz = cbuffer.data(gf_geom_20_off + 814 * ccomps * dcomps);

            auto g_zz_0_xyyy_xzz = cbuffer.data(gf_geom_20_off + 815 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyy = cbuffer.data(gf_geom_20_off + 816 * ccomps * dcomps);

            auto g_zz_0_xyyy_yyz = cbuffer.data(gf_geom_20_off + 817 * ccomps * dcomps);

            auto g_zz_0_xyyy_yzz = cbuffer.data(gf_geom_20_off + 818 * ccomps * dcomps);

            auto g_zz_0_xyyy_zzz = cbuffer.data(gf_geom_20_off + 819 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxx = cbuffer.data(gf_geom_20_off + 820 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxy = cbuffer.data(gf_geom_20_off + 821 * ccomps * dcomps);

            auto g_zz_0_xyyz_xxz = cbuffer.data(gf_geom_20_off + 822 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyy = cbuffer.data(gf_geom_20_off + 823 * ccomps * dcomps);

            auto g_zz_0_xyyz_xyz = cbuffer.data(gf_geom_20_off + 824 * ccomps * dcomps);

            auto g_zz_0_xyyz_xzz = cbuffer.data(gf_geom_20_off + 825 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyy = cbuffer.data(gf_geom_20_off + 826 * ccomps * dcomps);

            auto g_zz_0_xyyz_yyz = cbuffer.data(gf_geom_20_off + 827 * ccomps * dcomps);

            auto g_zz_0_xyyz_yzz = cbuffer.data(gf_geom_20_off + 828 * ccomps * dcomps);

            auto g_zz_0_xyyz_zzz = cbuffer.data(gf_geom_20_off + 829 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxx = cbuffer.data(gf_geom_20_off + 830 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxy = cbuffer.data(gf_geom_20_off + 831 * ccomps * dcomps);

            auto g_zz_0_xyzz_xxz = cbuffer.data(gf_geom_20_off + 832 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyy = cbuffer.data(gf_geom_20_off + 833 * ccomps * dcomps);

            auto g_zz_0_xyzz_xyz = cbuffer.data(gf_geom_20_off + 834 * ccomps * dcomps);

            auto g_zz_0_xyzz_xzz = cbuffer.data(gf_geom_20_off + 835 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyy = cbuffer.data(gf_geom_20_off + 836 * ccomps * dcomps);

            auto g_zz_0_xyzz_yyz = cbuffer.data(gf_geom_20_off + 837 * ccomps * dcomps);

            auto g_zz_0_xyzz_yzz = cbuffer.data(gf_geom_20_off + 838 * ccomps * dcomps);

            auto g_zz_0_xyzz_zzz = cbuffer.data(gf_geom_20_off + 839 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxx = cbuffer.data(gf_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxy = cbuffer.data(gf_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_xzzz_xxz = cbuffer.data(gf_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyy = cbuffer.data(gf_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_xzzz_xyz = cbuffer.data(gf_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_xzzz_xzz = cbuffer.data(gf_geom_20_off + 845 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyy = cbuffer.data(gf_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_xzzz_yyz = cbuffer.data(gf_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_xzzz_yzz = cbuffer.data(gf_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_xzzz_zzz = cbuffer.data(gf_geom_20_off + 849 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxx = cbuffer.data(gf_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxy = cbuffer.data(gf_geom_20_off + 851 * ccomps * dcomps);

            auto g_zz_0_yyyy_xxz = cbuffer.data(gf_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyy = cbuffer.data(gf_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_yyyy_xyz = cbuffer.data(gf_geom_20_off + 854 * ccomps * dcomps);

            auto g_zz_0_yyyy_xzz = cbuffer.data(gf_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyy = cbuffer.data(gf_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_yyyy_yyz = cbuffer.data(gf_geom_20_off + 857 * ccomps * dcomps);

            auto g_zz_0_yyyy_yzz = cbuffer.data(gf_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_yyyy_zzz = cbuffer.data(gf_geom_20_off + 859 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxx = cbuffer.data(gf_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxy = cbuffer.data(gf_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_yyyz_xxz = cbuffer.data(gf_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyy = cbuffer.data(gf_geom_20_off + 863 * ccomps * dcomps);

            auto g_zz_0_yyyz_xyz = cbuffer.data(gf_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_yyyz_xzz = cbuffer.data(gf_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyy = cbuffer.data(gf_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_yyyz_yyz = cbuffer.data(gf_geom_20_off + 867 * ccomps * dcomps);

            auto g_zz_0_yyyz_yzz = cbuffer.data(gf_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_yyyz_zzz = cbuffer.data(gf_geom_20_off + 869 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxx = cbuffer.data(gf_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxy = cbuffer.data(gf_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_yyzz_xxz = cbuffer.data(gf_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyy = cbuffer.data(gf_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_yyzz_xyz = cbuffer.data(gf_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_yyzz_xzz = cbuffer.data(gf_geom_20_off + 875 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyy = cbuffer.data(gf_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_yyzz_yyz = cbuffer.data(gf_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_yyzz_yzz = cbuffer.data(gf_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_yyzz_zzz = cbuffer.data(gf_geom_20_off + 879 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxx = cbuffer.data(gf_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxy = cbuffer.data(gf_geom_20_off + 881 * ccomps * dcomps);

            auto g_zz_0_yzzz_xxz = cbuffer.data(gf_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyy = cbuffer.data(gf_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_yzzz_xyz = cbuffer.data(gf_geom_20_off + 884 * ccomps * dcomps);

            auto g_zz_0_yzzz_xzz = cbuffer.data(gf_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyy = cbuffer.data(gf_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_yzzz_yyz = cbuffer.data(gf_geom_20_off + 887 * ccomps * dcomps);

            auto g_zz_0_yzzz_yzz = cbuffer.data(gf_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_yzzz_zzz = cbuffer.data(gf_geom_20_off + 889 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxx = cbuffer.data(gf_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxy = cbuffer.data(gf_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_zzzz_xxz = cbuffer.data(gf_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyy = cbuffer.data(gf_geom_20_off + 893 * ccomps * dcomps);

            auto g_zz_0_zzzz_xyz = cbuffer.data(gf_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_zzzz_xzz = cbuffer.data(gf_geom_20_off + 895 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyy = cbuffer.data(gf_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_zzzz_yyz = cbuffer.data(gf_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_zzzz_yzz = cbuffer.data(gf_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_zzzz_zzz = cbuffer.data(gf_geom_20_off + 899 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hdxx

            const auto hd_geom_20_off = idx_geom_20_hdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xx, g_x_0_xxxx_xy, g_x_0_xxxx_xz, g_x_0_xxxx_yy, g_x_0_xxxx_yz, g_x_0_xxxx_zz, g_xx_0_xxxx_xx, g_xx_0_xxxx_xxx, g_xx_0_xxxx_xxy, g_xx_0_xxxx_xxz, g_xx_0_xxxx_xy, g_xx_0_xxxx_xyy, g_xx_0_xxxx_xyz, g_xx_0_xxxx_xz, g_xx_0_xxxx_xzz, g_xx_0_xxxx_yy, g_xx_0_xxxx_yz, g_xx_0_xxxx_zz, g_xx_0_xxxxx_xx, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_yy, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxx_xx[k] = -2.0 * g_x_0_xxxx_xx[k] - g_xx_0_xxxx_xx[k] * ab_x + g_xx_0_xxxx_xxx[k];

                g_xx_0_xxxxx_xy[k] = -2.0 * g_x_0_xxxx_xy[k] - g_xx_0_xxxx_xy[k] * ab_x + g_xx_0_xxxx_xxy[k];

                g_xx_0_xxxxx_xz[k] = -2.0 * g_x_0_xxxx_xz[k] - g_xx_0_xxxx_xz[k] * ab_x + g_xx_0_xxxx_xxz[k];

                g_xx_0_xxxxx_yy[k] = -2.0 * g_x_0_xxxx_yy[k] - g_xx_0_xxxx_yy[k] * ab_x + g_xx_0_xxxx_xyy[k];

                g_xx_0_xxxxx_yz[k] = -2.0 * g_x_0_xxxx_yz[k] - g_xx_0_xxxx_yz[k] * ab_x + g_xx_0_xxxx_xyz[k];

                g_xx_0_xxxxx_zz[k] = -2.0 * g_x_0_xxxx_zz[k] - g_xx_0_xxxx_zz[k] * ab_x + g_xx_0_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_xx, g_xx_0_xxxx_xxy, g_xx_0_xxxx_xy, g_xx_0_xxxx_xyy, g_xx_0_xxxx_xyz, g_xx_0_xxxx_xz, g_xx_0_xxxx_yy, g_xx_0_xxxx_yyy, g_xx_0_xxxx_yyz, g_xx_0_xxxx_yz, g_xx_0_xxxx_yzz, g_xx_0_xxxx_zz, g_xx_0_xxxxy_xx, g_xx_0_xxxxy_xy, g_xx_0_xxxxy_xz, g_xx_0_xxxxy_yy, g_xx_0_xxxxy_yz, g_xx_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxy_xx[k] = -g_xx_0_xxxx_xx[k] * ab_y + g_xx_0_xxxx_xxy[k];

                g_xx_0_xxxxy_xy[k] = -g_xx_0_xxxx_xy[k] * ab_y + g_xx_0_xxxx_xyy[k];

                g_xx_0_xxxxy_xz[k] = -g_xx_0_xxxx_xz[k] * ab_y + g_xx_0_xxxx_xyz[k];

                g_xx_0_xxxxy_yy[k] = -g_xx_0_xxxx_yy[k] * ab_y + g_xx_0_xxxx_yyy[k];

                g_xx_0_xxxxy_yz[k] = -g_xx_0_xxxx_yz[k] * ab_y + g_xx_0_xxxx_yyz[k];

                g_xx_0_xxxxy_zz[k] = -g_xx_0_xxxx_zz[k] * ab_y + g_xx_0_xxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxx_xx, g_xx_0_xxxx_xxz, g_xx_0_xxxx_xy, g_xx_0_xxxx_xyz, g_xx_0_xxxx_xz, g_xx_0_xxxx_xzz, g_xx_0_xxxx_yy, g_xx_0_xxxx_yyz, g_xx_0_xxxx_yz, g_xx_0_xxxx_yzz, g_xx_0_xxxx_zz, g_xx_0_xxxx_zzz, g_xx_0_xxxxz_xx, g_xx_0_xxxxz_xy, g_xx_0_xxxxz_xz, g_xx_0_xxxxz_yy, g_xx_0_xxxxz_yz, g_xx_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxz_xx[k] = -g_xx_0_xxxx_xx[k] * ab_z + g_xx_0_xxxx_xxz[k];

                g_xx_0_xxxxz_xy[k] = -g_xx_0_xxxx_xy[k] * ab_z + g_xx_0_xxxx_xyz[k];

                g_xx_0_xxxxz_xz[k] = -g_xx_0_xxxx_xz[k] * ab_z + g_xx_0_xxxx_xzz[k];

                g_xx_0_xxxxz_yy[k] = -g_xx_0_xxxx_yy[k] * ab_z + g_xx_0_xxxx_yyz[k];

                g_xx_0_xxxxz_yz[k] = -g_xx_0_xxxx_yz[k] * ab_z + g_xx_0_xxxx_yzz[k];

                g_xx_0_xxxxz_zz[k] = -g_xx_0_xxxx_zz[k] * ab_z + g_xx_0_xxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxy_xx, g_xx_0_xxxy_xxy, g_xx_0_xxxy_xy, g_xx_0_xxxy_xyy, g_xx_0_xxxy_xyz, g_xx_0_xxxy_xz, g_xx_0_xxxy_yy, g_xx_0_xxxy_yyy, g_xx_0_xxxy_yyz, g_xx_0_xxxy_yz, g_xx_0_xxxy_yzz, g_xx_0_xxxy_zz, g_xx_0_xxxyy_xx, g_xx_0_xxxyy_xy, g_xx_0_xxxyy_xz, g_xx_0_xxxyy_yy, g_xx_0_xxxyy_yz, g_xx_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyy_xx[k] = -g_xx_0_xxxy_xx[k] * ab_y + g_xx_0_xxxy_xxy[k];

                g_xx_0_xxxyy_xy[k] = -g_xx_0_xxxy_xy[k] * ab_y + g_xx_0_xxxy_xyy[k];

                g_xx_0_xxxyy_xz[k] = -g_xx_0_xxxy_xz[k] * ab_y + g_xx_0_xxxy_xyz[k];

                g_xx_0_xxxyy_yy[k] = -g_xx_0_xxxy_yy[k] * ab_y + g_xx_0_xxxy_yyy[k];

                g_xx_0_xxxyy_yz[k] = -g_xx_0_xxxy_yz[k] * ab_y + g_xx_0_xxxy_yyz[k];

                g_xx_0_xxxyy_zz[k] = -g_xx_0_xxxy_zz[k] * ab_y + g_xx_0_xxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyz_xx, g_xx_0_xxxyz_xy, g_xx_0_xxxyz_xz, g_xx_0_xxxyz_yy, g_xx_0_xxxyz_yz, g_xx_0_xxxyz_zz, g_xx_0_xxxz_xx, g_xx_0_xxxz_xxy, g_xx_0_xxxz_xy, g_xx_0_xxxz_xyy, g_xx_0_xxxz_xyz, g_xx_0_xxxz_xz, g_xx_0_xxxz_yy, g_xx_0_xxxz_yyy, g_xx_0_xxxz_yyz, g_xx_0_xxxz_yz, g_xx_0_xxxz_yzz, g_xx_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyz_xx[k] = -g_xx_0_xxxz_xx[k] * ab_y + g_xx_0_xxxz_xxy[k];

                g_xx_0_xxxyz_xy[k] = -g_xx_0_xxxz_xy[k] * ab_y + g_xx_0_xxxz_xyy[k];

                g_xx_0_xxxyz_xz[k] = -g_xx_0_xxxz_xz[k] * ab_y + g_xx_0_xxxz_xyz[k];

                g_xx_0_xxxyz_yy[k] = -g_xx_0_xxxz_yy[k] * ab_y + g_xx_0_xxxz_yyy[k];

                g_xx_0_xxxyz_yz[k] = -g_xx_0_xxxz_yz[k] * ab_y + g_xx_0_xxxz_yyz[k];

                g_xx_0_xxxyz_zz[k] = -g_xx_0_xxxz_zz[k] * ab_y + g_xx_0_xxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxz_xx, g_xx_0_xxxz_xxz, g_xx_0_xxxz_xy, g_xx_0_xxxz_xyz, g_xx_0_xxxz_xz, g_xx_0_xxxz_xzz, g_xx_0_xxxz_yy, g_xx_0_xxxz_yyz, g_xx_0_xxxz_yz, g_xx_0_xxxz_yzz, g_xx_0_xxxz_zz, g_xx_0_xxxz_zzz, g_xx_0_xxxzz_xx, g_xx_0_xxxzz_xy, g_xx_0_xxxzz_xz, g_xx_0_xxxzz_yy, g_xx_0_xxxzz_yz, g_xx_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzz_xx[k] = -g_xx_0_xxxz_xx[k] * ab_z + g_xx_0_xxxz_xxz[k];

                g_xx_0_xxxzz_xy[k] = -g_xx_0_xxxz_xy[k] * ab_z + g_xx_0_xxxz_xyz[k];

                g_xx_0_xxxzz_xz[k] = -g_xx_0_xxxz_xz[k] * ab_z + g_xx_0_xxxz_xzz[k];

                g_xx_0_xxxzz_yy[k] = -g_xx_0_xxxz_yy[k] * ab_z + g_xx_0_xxxz_yyz[k];

                g_xx_0_xxxzz_yz[k] = -g_xx_0_xxxz_yz[k] * ab_z + g_xx_0_xxxz_yzz[k];

                g_xx_0_xxxzz_zz[k] = -g_xx_0_xxxz_zz[k] * ab_z + g_xx_0_xxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyy_xx, g_xx_0_xxyy_xxy, g_xx_0_xxyy_xy, g_xx_0_xxyy_xyy, g_xx_0_xxyy_xyz, g_xx_0_xxyy_xz, g_xx_0_xxyy_yy, g_xx_0_xxyy_yyy, g_xx_0_xxyy_yyz, g_xx_0_xxyy_yz, g_xx_0_xxyy_yzz, g_xx_0_xxyy_zz, g_xx_0_xxyyy_xx, g_xx_0_xxyyy_xy, g_xx_0_xxyyy_xz, g_xx_0_xxyyy_yy, g_xx_0_xxyyy_yz, g_xx_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyy_xx[k] = -g_xx_0_xxyy_xx[k] * ab_y + g_xx_0_xxyy_xxy[k];

                g_xx_0_xxyyy_xy[k] = -g_xx_0_xxyy_xy[k] * ab_y + g_xx_0_xxyy_xyy[k];

                g_xx_0_xxyyy_xz[k] = -g_xx_0_xxyy_xz[k] * ab_y + g_xx_0_xxyy_xyz[k];

                g_xx_0_xxyyy_yy[k] = -g_xx_0_xxyy_yy[k] * ab_y + g_xx_0_xxyy_yyy[k];

                g_xx_0_xxyyy_yz[k] = -g_xx_0_xxyy_yz[k] * ab_y + g_xx_0_xxyy_yyz[k];

                g_xx_0_xxyyy_zz[k] = -g_xx_0_xxyy_zz[k] * ab_y + g_xx_0_xxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyz_xx, g_xx_0_xxyyz_xy, g_xx_0_xxyyz_xz, g_xx_0_xxyyz_yy, g_xx_0_xxyyz_yz, g_xx_0_xxyyz_zz, g_xx_0_xxyz_xx, g_xx_0_xxyz_xxy, g_xx_0_xxyz_xy, g_xx_0_xxyz_xyy, g_xx_0_xxyz_xyz, g_xx_0_xxyz_xz, g_xx_0_xxyz_yy, g_xx_0_xxyz_yyy, g_xx_0_xxyz_yyz, g_xx_0_xxyz_yz, g_xx_0_xxyz_yzz, g_xx_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyz_xx[k] = -g_xx_0_xxyz_xx[k] * ab_y + g_xx_0_xxyz_xxy[k];

                g_xx_0_xxyyz_xy[k] = -g_xx_0_xxyz_xy[k] * ab_y + g_xx_0_xxyz_xyy[k];

                g_xx_0_xxyyz_xz[k] = -g_xx_0_xxyz_xz[k] * ab_y + g_xx_0_xxyz_xyz[k];

                g_xx_0_xxyyz_yy[k] = -g_xx_0_xxyz_yy[k] * ab_y + g_xx_0_xxyz_yyy[k];

                g_xx_0_xxyyz_yz[k] = -g_xx_0_xxyz_yz[k] * ab_y + g_xx_0_xxyz_yyz[k];

                g_xx_0_xxyyz_zz[k] = -g_xx_0_xxyz_zz[k] * ab_y + g_xx_0_xxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzz_xx, g_xx_0_xxyzz_xy, g_xx_0_xxyzz_xz, g_xx_0_xxyzz_yy, g_xx_0_xxyzz_yz, g_xx_0_xxyzz_zz, g_xx_0_xxzz_xx, g_xx_0_xxzz_xxy, g_xx_0_xxzz_xy, g_xx_0_xxzz_xyy, g_xx_0_xxzz_xyz, g_xx_0_xxzz_xz, g_xx_0_xxzz_yy, g_xx_0_xxzz_yyy, g_xx_0_xxzz_yyz, g_xx_0_xxzz_yz, g_xx_0_xxzz_yzz, g_xx_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzz_xx[k] = -g_xx_0_xxzz_xx[k] * ab_y + g_xx_0_xxzz_xxy[k];

                g_xx_0_xxyzz_xy[k] = -g_xx_0_xxzz_xy[k] * ab_y + g_xx_0_xxzz_xyy[k];

                g_xx_0_xxyzz_xz[k] = -g_xx_0_xxzz_xz[k] * ab_y + g_xx_0_xxzz_xyz[k];

                g_xx_0_xxyzz_yy[k] = -g_xx_0_xxzz_yy[k] * ab_y + g_xx_0_xxzz_yyy[k];

                g_xx_0_xxyzz_yz[k] = -g_xx_0_xxzz_yz[k] * ab_y + g_xx_0_xxzz_yyz[k];

                g_xx_0_xxyzz_zz[k] = -g_xx_0_xxzz_zz[k] * ab_y + g_xx_0_xxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzz_xx, g_xx_0_xxzz_xxz, g_xx_0_xxzz_xy, g_xx_0_xxzz_xyz, g_xx_0_xxzz_xz, g_xx_0_xxzz_xzz, g_xx_0_xxzz_yy, g_xx_0_xxzz_yyz, g_xx_0_xxzz_yz, g_xx_0_xxzz_yzz, g_xx_0_xxzz_zz, g_xx_0_xxzz_zzz, g_xx_0_xxzzz_xx, g_xx_0_xxzzz_xy, g_xx_0_xxzzz_xz, g_xx_0_xxzzz_yy, g_xx_0_xxzzz_yz, g_xx_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzz_xx[k] = -g_xx_0_xxzz_xx[k] * ab_z + g_xx_0_xxzz_xxz[k];

                g_xx_0_xxzzz_xy[k] = -g_xx_0_xxzz_xy[k] * ab_z + g_xx_0_xxzz_xyz[k];

                g_xx_0_xxzzz_xz[k] = -g_xx_0_xxzz_xz[k] * ab_z + g_xx_0_xxzz_xzz[k];

                g_xx_0_xxzzz_yy[k] = -g_xx_0_xxzz_yy[k] * ab_z + g_xx_0_xxzz_yyz[k];

                g_xx_0_xxzzz_yz[k] = -g_xx_0_xxzz_yz[k] * ab_z + g_xx_0_xxzz_yzz[k];

                g_xx_0_xxzzz_zz[k] = -g_xx_0_xxzz_zz[k] * ab_z + g_xx_0_xxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyy_xx, g_xx_0_xyyy_xxy, g_xx_0_xyyy_xy, g_xx_0_xyyy_xyy, g_xx_0_xyyy_xyz, g_xx_0_xyyy_xz, g_xx_0_xyyy_yy, g_xx_0_xyyy_yyy, g_xx_0_xyyy_yyz, g_xx_0_xyyy_yz, g_xx_0_xyyy_yzz, g_xx_0_xyyy_zz, g_xx_0_xyyyy_xx, g_xx_0_xyyyy_xy, g_xx_0_xyyyy_xz, g_xx_0_xyyyy_yy, g_xx_0_xyyyy_yz, g_xx_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyy_xx[k] = -g_xx_0_xyyy_xx[k] * ab_y + g_xx_0_xyyy_xxy[k];

                g_xx_0_xyyyy_xy[k] = -g_xx_0_xyyy_xy[k] * ab_y + g_xx_0_xyyy_xyy[k];

                g_xx_0_xyyyy_xz[k] = -g_xx_0_xyyy_xz[k] * ab_y + g_xx_0_xyyy_xyz[k];

                g_xx_0_xyyyy_yy[k] = -g_xx_0_xyyy_yy[k] * ab_y + g_xx_0_xyyy_yyy[k];

                g_xx_0_xyyyy_yz[k] = -g_xx_0_xyyy_yz[k] * ab_y + g_xx_0_xyyy_yyz[k];

                g_xx_0_xyyyy_zz[k] = -g_xx_0_xyyy_zz[k] * ab_y + g_xx_0_xyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyz_xx, g_xx_0_xyyyz_xy, g_xx_0_xyyyz_xz, g_xx_0_xyyyz_yy, g_xx_0_xyyyz_yz, g_xx_0_xyyyz_zz, g_xx_0_xyyz_xx, g_xx_0_xyyz_xxy, g_xx_0_xyyz_xy, g_xx_0_xyyz_xyy, g_xx_0_xyyz_xyz, g_xx_0_xyyz_xz, g_xx_0_xyyz_yy, g_xx_0_xyyz_yyy, g_xx_0_xyyz_yyz, g_xx_0_xyyz_yz, g_xx_0_xyyz_yzz, g_xx_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyz_xx[k] = -g_xx_0_xyyz_xx[k] * ab_y + g_xx_0_xyyz_xxy[k];

                g_xx_0_xyyyz_xy[k] = -g_xx_0_xyyz_xy[k] * ab_y + g_xx_0_xyyz_xyy[k];

                g_xx_0_xyyyz_xz[k] = -g_xx_0_xyyz_xz[k] * ab_y + g_xx_0_xyyz_xyz[k];

                g_xx_0_xyyyz_yy[k] = -g_xx_0_xyyz_yy[k] * ab_y + g_xx_0_xyyz_yyy[k];

                g_xx_0_xyyyz_yz[k] = -g_xx_0_xyyz_yz[k] * ab_y + g_xx_0_xyyz_yyz[k];

                g_xx_0_xyyyz_zz[k] = -g_xx_0_xyyz_zz[k] * ab_y + g_xx_0_xyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzz_xx, g_xx_0_xyyzz_xy, g_xx_0_xyyzz_xz, g_xx_0_xyyzz_yy, g_xx_0_xyyzz_yz, g_xx_0_xyyzz_zz, g_xx_0_xyzz_xx, g_xx_0_xyzz_xxy, g_xx_0_xyzz_xy, g_xx_0_xyzz_xyy, g_xx_0_xyzz_xyz, g_xx_0_xyzz_xz, g_xx_0_xyzz_yy, g_xx_0_xyzz_yyy, g_xx_0_xyzz_yyz, g_xx_0_xyzz_yz, g_xx_0_xyzz_yzz, g_xx_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzz_xx[k] = -g_xx_0_xyzz_xx[k] * ab_y + g_xx_0_xyzz_xxy[k];

                g_xx_0_xyyzz_xy[k] = -g_xx_0_xyzz_xy[k] * ab_y + g_xx_0_xyzz_xyy[k];

                g_xx_0_xyyzz_xz[k] = -g_xx_0_xyzz_xz[k] * ab_y + g_xx_0_xyzz_xyz[k];

                g_xx_0_xyyzz_yy[k] = -g_xx_0_xyzz_yy[k] * ab_y + g_xx_0_xyzz_yyy[k];

                g_xx_0_xyyzz_yz[k] = -g_xx_0_xyzz_yz[k] * ab_y + g_xx_0_xyzz_yyz[k];

                g_xx_0_xyyzz_zz[k] = -g_xx_0_xyzz_zz[k] * ab_y + g_xx_0_xyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzz_xx, g_xx_0_xyzzz_xy, g_xx_0_xyzzz_xz, g_xx_0_xyzzz_yy, g_xx_0_xyzzz_yz, g_xx_0_xyzzz_zz, g_xx_0_xzzz_xx, g_xx_0_xzzz_xxy, g_xx_0_xzzz_xy, g_xx_0_xzzz_xyy, g_xx_0_xzzz_xyz, g_xx_0_xzzz_xz, g_xx_0_xzzz_yy, g_xx_0_xzzz_yyy, g_xx_0_xzzz_yyz, g_xx_0_xzzz_yz, g_xx_0_xzzz_yzz, g_xx_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzz_xx[k] = -g_xx_0_xzzz_xx[k] * ab_y + g_xx_0_xzzz_xxy[k];

                g_xx_0_xyzzz_xy[k] = -g_xx_0_xzzz_xy[k] * ab_y + g_xx_0_xzzz_xyy[k];

                g_xx_0_xyzzz_xz[k] = -g_xx_0_xzzz_xz[k] * ab_y + g_xx_0_xzzz_xyz[k];

                g_xx_0_xyzzz_yy[k] = -g_xx_0_xzzz_yy[k] * ab_y + g_xx_0_xzzz_yyy[k];

                g_xx_0_xyzzz_yz[k] = -g_xx_0_xzzz_yz[k] * ab_y + g_xx_0_xzzz_yyz[k];

                g_xx_0_xyzzz_zz[k] = -g_xx_0_xzzz_zz[k] * ab_y + g_xx_0_xzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzz_xx, g_xx_0_xzzz_xxz, g_xx_0_xzzz_xy, g_xx_0_xzzz_xyz, g_xx_0_xzzz_xz, g_xx_0_xzzz_xzz, g_xx_0_xzzz_yy, g_xx_0_xzzz_yyz, g_xx_0_xzzz_yz, g_xx_0_xzzz_yzz, g_xx_0_xzzz_zz, g_xx_0_xzzz_zzz, g_xx_0_xzzzz_xx, g_xx_0_xzzzz_xy, g_xx_0_xzzzz_xz, g_xx_0_xzzzz_yy, g_xx_0_xzzzz_yz, g_xx_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzz_xx[k] = -g_xx_0_xzzz_xx[k] * ab_z + g_xx_0_xzzz_xxz[k];

                g_xx_0_xzzzz_xy[k] = -g_xx_0_xzzz_xy[k] * ab_z + g_xx_0_xzzz_xyz[k];

                g_xx_0_xzzzz_xz[k] = -g_xx_0_xzzz_xz[k] * ab_z + g_xx_0_xzzz_xzz[k];

                g_xx_0_xzzzz_yy[k] = -g_xx_0_xzzz_yy[k] * ab_z + g_xx_0_xzzz_yyz[k];

                g_xx_0_xzzzz_yz[k] = -g_xx_0_xzzz_yz[k] * ab_z + g_xx_0_xzzz_yzz[k];

                g_xx_0_xzzzz_zz[k] = -g_xx_0_xzzz_zz[k] * ab_z + g_xx_0_xzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyy_xx, g_xx_0_yyyy_xxy, g_xx_0_yyyy_xy, g_xx_0_yyyy_xyy, g_xx_0_yyyy_xyz, g_xx_0_yyyy_xz, g_xx_0_yyyy_yy, g_xx_0_yyyy_yyy, g_xx_0_yyyy_yyz, g_xx_0_yyyy_yz, g_xx_0_yyyy_yzz, g_xx_0_yyyy_zz, g_xx_0_yyyyy_xx, g_xx_0_yyyyy_xy, g_xx_0_yyyyy_xz, g_xx_0_yyyyy_yy, g_xx_0_yyyyy_yz, g_xx_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyy_xx[k] = -g_xx_0_yyyy_xx[k] * ab_y + g_xx_0_yyyy_xxy[k];

                g_xx_0_yyyyy_xy[k] = -g_xx_0_yyyy_xy[k] * ab_y + g_xx_0_yyyy_xyy[k];

                g_xx_0_yyyyy_xz[k] = -g_xx_0_yyyy_xz[k] * ab_y + g_xx_0_yyyy_xyz[k];

                g_xx_0_yyyyy_yy[k] = -g_xx_0_yyyy_yy[k] * ab_y + g_xx_0_yyyy_yyy[k];

                g_xx_0_yyyyy_yz[k] = -g_xx_0_yyyy_yz[k] * ab_y + g_xx_0_yyyy_yyz[k];

                g_xx_0_yyyyy_zz[k] = -g_xx_0_yyyy_zz[k] * ab_y + g_xx_0_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyz_xx, g_xx_0_yyyyz_xy, g_xx_0_yyyyz_xz, g_xx_0_yyyyz_yy, g_xx_0_yyyyz_yz, g_xx_0_yyyyz_zz, g_xx_0_yyyz_xx, g_xx_0_yyyz_xxy, g_xx_0_yyyz_xy, g_xx_0_yyyz_xyy, g_xx_0_yyyz_xyz, g_xx_0_yyyz_xz, g_xx_0_yyyz_yy, g_xx_0_yyyz_yyy, g_xx_0_yyyz_yyz, g_xx_0_yyyz_yz, g_xx_0_yyyz_yzz, g_xx_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyz_xx[k] = -g_xx_0_yyyz_xx[k] * ab_y + g_xx_0_yyyz_xxy[k];

                g_xx_0_yyyyz_xy[k] = -g_xx_0_yyyz_xy[k] * ab_y + g_xx_0_yyyz_xyy[k];

                g_xx_0_yyyyz_xz[k] = -g_xx_0_yyyz_xz[k] * ab_y + g_xx_0_yyyz_xyz[k];

                g_xx_0_yyyyz_yy[k] = -g_xx_0_yyyz_yy[k] * ab_y + g_xx_0_yyyz_yyy[k];

                g_xx_0_yyyyz_yz[k] = -g_xx_0_yyyz_yz[k] * ab_y + g_xx_0_yyyz_yyz[k];

                g_xx_0_yyyyz_zz[k] = -g_xx_0_yyyz_zz[k] * ab_y + g_xx_0_yyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzz_xx, g_xx_0_yyyzz_xy, g_xx_0_yyyzz_xz, g_xx_0_yyyzz_yy, g_xx_0_yyyzz_yz, g_xx_0_yyyzz_zz, g_xx_0_yyzz_xx, g_xx_0_yyzz_xxy, g_xx_0_yyzz_xy, g_xx_0_yyzz_xyy, g_xx_0_yyzz_xyz, g_xx_0_yyzz_xz, g_xx_0_yyzz_yy, g_xx_0_yyzz_yyy, g_xx_0_yyzz_yyz, g_xx_0_yyzz_yz, g_xx_0_yyzz_yzz, g_xx_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzz_xx[k] = -g_xx_0_yyzz_xx[k] * ab_y + g_xx_0_yyzz_xxy[k];

                g_xx_0_yyyzz_xy[k] = -g_xx_0_yyzz_xy[k] * ab_y + g_xx_0_yyzz_xyy[k];

                g_xx_0_yyyzz_xz[k] = -g_xx_0_yyzz_xz[k] * ab_y + g_xx_0_yyzz_xyz[k];

                g_xx_0_yyyzz_yy[k] = -g_xx_0_yyzz_yy[k] * ab_y + g_xx_0_yyzz_yyy[k];

                g_xx_0_yyyzz_yz[k] = -g_xx_0_yyzz_yz[k] * ab_y + g_xx_0_yyzz_yyz[k];

                g_xx_0_yyyzz_zz[k] = -g_xx_0_yyzz_zz[k] * ab_y + g_xx_0_yyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzz_xx, g_xx_0_yyzzz_xy, g_xx_0_yyzzz_xz, g_xx_0_yyzzz_yy, g_xx_0_yyzzz_yz, g_xx_0_yyzzz_zz, g_xx_0_yzzz_xx, g_xx_0_yzzz_xxy, g_xx_0_yzzz_xy, g_xx_0_yzzz_xyy, g_xx_0_yzzz_xyz, g_xx_0_yzzz_xz, g_xx_0_yzzz_yy, g_xx_0_yzzz_yyy, g_xx_0_yzzz_yyz, g_xx_0_yzzz_yz, g_xx_0_yzzz_yzz, g_xx_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzz_xx[k] = -g_xx_0_yzzz_xx[k] * ab_y + g_xx_0_yzzz_xxy[k];

                g_xx_0_yyzzz_xy[k] = -g_xx_0_yzzz_xy[k] * ab_y + g_xx_0_yzzz_xyy[k];

                g_xx_0_yyzzz_xz[k] = -g_xx_0_yzzz_xz[k] * ab_y + g_xx_0_yzzz_xyz[k];

                g_xx_0_yyzzz_yy[k] = -g_xx_0_yzzz_yy[k] * ab_y + g_xx_0_yzzz_yyy[k];

                g_xx_0_yyzzz_yz[k] = -g_xx_0_yzzz_yz[k] * ab_y + g_xx_0_yzzz_yyz[k];

                g_xx_0_yyzzz_zz[k] = -g_xx_0_yzzz_zz[k] * ab_y + g_xx_0_yzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzz_xx, g_xx_0_yzzzz_xy, g_xx_0_yzzzz_xz, g_xx_0_yzzzz_yy, g_xx_0_yzzzz_yz, g_xx_0_yzzzz_zz, g_xx_0_zzzz_xx, g_xx_0_zzzz_xxy, g_xx_0_zzzz_xy, g_xx_0_zzzz_xyy, g_xx_0_zzzz_xyz, g_xx_0_zzzz_xz, g_xx_0_zzzz_yy, g_xx_0_zzzz_yyy, g_xx_0_zzzz_yyz, g_xx_0_zzzz_yz, g_xx_0_zzzz_yzz, g_xx_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzz_xx[k] = -g_xx_0_zzzz_xx[k] * ab_y + g_xx_0_zzzz_xxy[k];

                g_xx_0_yzzzz_xy[k] = -g_xx_0_zzzz_xy[k] * ab_y + g_xx_0_zzzz_xyy[k];

                g_xx_0_yzzzz_xz[k] = -g_xx_0_zzzz_xz[k] * ab_y + g_xx_0_zzzz_xyz[k];

                g_xx_0_yzzzz_yy[k] = -g_xx_0_zzzz_yy[k] * ab_y + g_xx_0_zzzz_yyy[k];

                g_xx_0_yzzzz_yz[k] = -g_xx_0_zzzz_yz[k] * ab_y + g_xx_0_zzzz_yyz[k];

                g_xx_0_yzzzz_zz[k] = -g_xx_0_zzzz_zz[k] * ab_y + g_xx_0_zzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzz_xx, g_xx_0_zzzz_xxz, g_xx_0_zzzz_xy, g_xx_0_zzzz_xyz, g_xx_0_zzzz_xz, g_xx_0_zzzz_xzz, g_xx_0_zzzz_yy, g_xx_0_zzzz_yyz, g_xx_0_zzzz_yz, g_xx_0_zzzz_yzz, g_xx_0_zzzz_zz, g_xx_0_zzzz_zzz, g_xx_0_zzzzz_xx, g_xx_0_zzzzz_xy, g_xx_0_zzzzz_xz, g_xx_0_zzzzz_yy, g_xx_0_zzzzz_yz, g_xx_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzz_xx[k] = -g_xx_0_zzzz_xx[k] * ab_z + g_xx_0_zzzz_xxz[k];

                g_xx_0_zzzzz_xy[k] = -g_xx_0_zzzz_xy[k] * ab_z + g_xx_0_zzzz_xyz[k];

                g_xx_0_zzzzz_xz[k] = -g_xx_0_zzzz_xz[k] * ab_z + g_xx_0_zzzz_xzz[k];

                g_xx_0_zzzzz_yy[k] = -g_xx_0_zzzz_yy[k] * ab_z + g_xx_0_zzzz_yyz[k];

                g_xx_0_zzzzz_yz[k] = -g_xx_0_zzzz_yz[k] * ab_z + g_xx_0_zzzz_yzz[k];

                g_xx_0_zzzzz_zz[k] = -g_xx_0_zzzz_zz[k] * ab_z + g_xx_0_zzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_xx, g_xy_0_xxxx_xxx, g_xy_0_xxxx_xxy, g_xy_0_xxxx_xxz, g_xy_0_xxxx_xy, g_xy_0_xxxx_xyy, g_xy_0_xxxx_xyz, g_xy_0_xxxx_xz, g_xy_0_xxxx_xzz, g_xy_0_xxxx_yy, g_xy_0_xxxx_yz, g_xy_0_xxxx_zz, g_xy_0_xxxxx_xx, g_xy_0_xxxxx_xy, g_xy_0_xxxxx_xz, g_xy_0_xxxxx_yy, g_xy_0_xxxxx_yz, g_xy_0_xxxxx_zz, g_y_0_xxxx_xx, g_y_0_xxxx_xy, g_y_0_xxxx_xz, g_y_0_xxxx_yy, g_y_0_xxxx_yz, g_y_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxx_xx[k] = -g_y_0_xxxx_xx[k] - g_xy_0_xxxx_xx[k] * ab_x + g_xy_0_xxxx_xxx[k];

                g_xy_0_xxxxx_xy[k] = -g_y_0_xxxx_xy[k] - g_xy_0_xxxx_xy[k] * ab_x + g_xy_0_xxxx_xxy[k];

                g_xy_0_xxxxx_xz[k] = -g_y_0_xxxx_xz[k] - g_xy_0_xxxx_xz[k] * ab_x + g_xy_0_xxxx_xxz[k];

                g_xy_0_xxxxx_yy[k] = -g_y_0_xxxx_yy[k] - g_xy_0_xxxx_yy[k] * ab_x + g_xy_0_xxxx_xyy[k];

                g_xy_0_xxxxx_yz[k] = -g_y_0_xxxx_yz[k] - g_xy_0_xxxx_yz[k] * ab_x + g_xy_0_xxxx_xyz[k];

                g_xy_0_xxxxx_zz[k] = -g_y_0_xxxx_zz[k] - g_xy_0_xxxx_zz[k] * ab_x + g_xy_0_xxxx_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_xx, g_xy_0_xxxxy_xy, g_xy_0_xxxxy_xz, g_xy_0_xxxxy_yy, g_xy_0_xxxxy_yz, g_xy_0_xxxxy_zz, g_xy_0_xxxy_xx, g_xy_0_xxxy_xxx, g_xy_0_xxxy_xxy, g_xy_0_xxxy_xxz, g_xy_0_xxxy_xy, g_xy_0_xxxy_xyy, g_xy_0_xxxy_xyz, g_xy_0_xxxy_xz, g_xy_0_xxxy_xzz, g_xy_0_xxxy_yy, g_xy_0_xxxy_yz, g_xy_0_xxxy_zz, g_y_0_xxxy_xx, g_y_0_xxxy_xy, g_y_0_xxxy_xz, g_y_0_xxxy_yy, g_y_0_xxxy_yz, g_y_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxy_xx[k] = -g_y_0_xxxy_xx[k] - g_xy_0_xxxy_xx[k] * ab_x + g_xy_0_xxxy_xxx[k];

                g_xy_0_xxxxy_xy[k] = -g_y_0_xxxy_xy[k] - g_xy_0_xxxy_xy[k] * ab_x + g_xy_0_xxxy_xxy[k];

                g_xy_0_xxxxy_xz[k] = -g_y_0_xxxy_xz[k] - g_xy_0_xxxy_xz[k] * ab_x + g_xy_0_xxxy_xxz[k];

                g_xy_0_xxxxy_yy[k] = -g_y_0_xxxy_yy[k] - g_xy_0_xxxy_yy[k] * ab_x + g_xy_0_xxxy_xyy[k];

                g_xy_0_xxxxy_yz[k] = -g_y_0_xxxy_yz[k] - g_xy_0_xxxy_yz[k] * ab_x + g_xy_0_xxxy_xyz[k];

                g_xy_0_xxxxy_zz[k] = -g_y_0_xxxy_zz[k] - g_xy_0_xxxy_zz[k] * ab_x + g_xy_0_xxxy_xzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxx_xx, g_xy_0_xxxx_xxz, g_xy_0_xxxx_xy, g_xy_0_xxxx_xyz, g_xy_0_xxxx_xz, g_xy_0_xxxx_xzz, g_xy_0_xxxx_yy, g_xy_0_xxxx_yyz, g_xy_0_xxxx_yz, g_xy_0_xxxx_yzz, g_xy_0_xxxx_zz, g_xy_0_xxxx_zzz, g_xy_0_xxxxz_xx, g_xy_0_xxxxz_xy, g_xy_0_xxxxz_xz, g_xy_0_xxxxz_yy, g_xy_0_xxxxz_yz, g_xy_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxz_xx[k] = -g_xy_0_xxxx_xx[k] * ab_z + g_xy_0_xxxx_xxz[k];

                g_xy_0_xxxxz_xy[k] = -g_xy_0_xxxx_xy[k] * ab_z + g_xy_0_xxxx_xyz[k];

                g_xy_0_xxxxz_xz[k] = -g_xy_0_xxxx_xz[k] * ab_z + g_xy_0_xxxx_xzz[k];

                g_xy_0_xxxxz_yy[k] = -g_xy_0_xxxx_yy[k] * ab_z + g_xy_0_xxxx_yyz[k];

                g_xy_0_xxxxz_yz[k] = -g_xy_0_xxxx_yz[k] * ab_z + g_xy_0_xxxx_yzz[k];

                g_xy_0_xxxxz_zz[k] = -g_xy_0_xxxx_zz[k] * ab_z + g_xy_0_xxxx_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_xx, g_xy_0_xxxyy_xy, g_xy_0_xxxyy_xz, g_xy_0_xxxyy_yy, g_xy_0_xxxyy_yz, g_xy_0_xxxyy_zz, g_xy_0_xxyy_xx, g_xy_0_xxyy_xxx, g_xy_0_xxyy_xxy, g_xy_0_xxyy_xxz, g_xy_0_xxyy_xy, g_xy_0_xxyy_xyy, g_xy_0_xxyy_xyz, g_xy_0_xxyy_xz, g_xy_0_xxyy_xzz, g_xy_0_xxyy_yy, g_xy_0_xxyy_yz, g_xy_0_xxyy_zz, g_y_0_xxyy_xx, g_y_0_xxyy_xy, g_y_0_xxyy_xz, g_y_0_xxyy_yy, g_y_0_xxyy_yz, g_y_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyy_xx[k] = -g_y_0_xxyy_xx[k] - g_xy_0_xxyy_xx[k] * ab_x + g_xy_0_xxyy_xxx[k];

                g_xy_0_xxxyy_xy[k] = -g_y_0_xxyy_xy[k] - g_xy_0_xxyy_xy[k] * ab_x + g_xy_0_xxyy_xxy[k];

                g_xy_0_xxxyy_xz[k] = -g_y_0_xxyy_xz[k] - g_xy_0_xxyy_xz[k] * ab_x + g_xy_0_xxyy_xxz[k];

                g_xy_0_xxxyy_yy[k] = -g_y_0_xxyy_yy[k] - g_xy_0_xxyy_yy[k] * ab_x + g_xy_0_xxyy_xyy[k];

                g_xy_0_xxxyy_yz[k] = -g_y_0_xxyy_yz[k] - g_xy_0_xxyy_yz[k] * ab_x + g_xy_0_xxyy_xyz[k];

                g_xy_0_xxxyy_zz[k] = -g_y_0_xxyy_zz[k] - g_xy_0_xxyy_zz[k] * ab_x + g_xy_0_xxyy_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_xx, g_xy_0_xxxy_xxz, g_xy_0_xxxy_xy, g_xy_0_xxxy_xyz, g_xy_0_xxxy_xz, g_xy_0_xxxy_xzz, g_xy_0_xxxy_yy, g_xy_0_xxxy_yyz, g_xy_0_xxxy_yz, g_xy_0_xxxy_yzz, g_xy_0_xxxy_zz, g_xy_0_xxxy_zzz, g_xy_0_xxxyz_xx, g_xy_0_xxxyz_xy, g_xy_0_xxxyz_xz, g_xy_0_xxxyz_yy, g_xy_0_xxxyz_yz, g_xy_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyz_xx[k] = -g_xy_0_xxxy_xx[k] * ab_z + g_xy_0_xxxy_xxz[k];

                g_xy_0_xxxyz_xy[k] = -g_xy_0_xxxy_xy[k] * ab_z + g_xy_0_xxxy_xyz[k];

                g_xy_0_xxxyz_xz[k] = -g_xy_0_xxxy_xz[k] * ab_z + g_xy_0_xxxy_xzz[k];

                g_xy_0_xxxyz_yy[k] = -g_xy_0_xxxy_yy[k] * ab_z + g_xy_0_xxxy_yyz[k];

                g_xy_0_xxxyz_yz[k] = -g_xy_0_xxxy_yz[k] * ab_z + g_xy_0_xxxy_yzz[k];

                g_xy_0_xxxyz_zz[k] = -g_xy_0_xxxy_zz[k] * ab_z + g_xy_0_xxxy_zzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxz_xx, g_xy_0_xxxz_xxz, g_xy_0_xxxz_xy, g_xy_0_xxxz_xyz, g_xy_0_xxxz_xz, g_xy_0_xxxz_xzz, g_xy_0_xxxz_yy, g_xy_0_xxxz_yyz, g_xy_0_xxxz_yz, g_xy_0_xxxz_yzz, g_xy_0_xxxz_zz, g_xy_0_xxxz_zzz, g_xy_0_xxxzz_xx, g_xy_0_xxxzz_xy, g_xy_0_xxxzz_xz, g_xy_0_xxxzz_yy, g_xy_0_xxxzz_yz, g_xy_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzz_xx[k] = -g_xy_0_xxxz_xx[k] * ab_z + g_xy_0_xxxz_xxz[k];

                g_xy_0_xxxzz_xy[k] = -g_xy_0_xxxz_xy[k] * ab_z + g_xy_0_xxxz_xyz[k];

                g_xy_0_xxxzz_xz[k] = -g_xy_0_xxxz_xz[k] * ab_z + g_xy_0_xxxz_xzz[k];

                g_xy_0_xxxzz_yy[k] = -g_xy_0_xxxz_yy[k] * ab_z + g_xy_0_xxxz_yyz[k];

                g_xy_0_xxxzz_yz[k] = -g_xy_0_xxxz_yz[k] * ab_z + g_xy_0_xxxz_yzz[k];

                g_xy_0_xxxzz_zz[k] = -g_xy_0_xxxz_zz[k] * ab_z + g_xy_0_xxxz_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_xx, g_xy_0_xxyyy_xy, g_xy_0_xxyyy_xz, g_xy_0_xxyyy_yy, g_xy_0_xxyyy_yz, g_xy_0_xxyyy_zz, g_xy_0_xyyy_xx, g_xy_0_xyyy_xxx, g_xy_0_xyyy_xxy, g_xy_0_xyyy_xxz, g_xy_0_xyyy_xy, g_xy_0_xyyy_xyy, g_xy_0_xyyy_xyz, g_xy_0_xyyy_xz, g_xy_0_xyyy_xzz, g_xy_0_xyyy_yy, g_xy_0_xyyy_yz, g_xy_0_xyyy_zz, g_y_0_xyyy_xx, g_y_0_xyyy_xy, g_y_0_xyyy_xz, g_y_0_xyyy_yy, g_y_0_xyyy_yz, g_y_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyy_xx[k] = -g_y_0_xyyy_xx[k] - g_xy_0_xyyy_xx[k] * ab_x + g_xy_0_xyyy_xxx[k];

                g_xy_0_xxyyy_xy[k] = -g_y_0_xyyy_xy[k] - g_xy_0_xyyy_xy[k] * ab_x + g_xy_0_xyyy_xxy[k];

                g_xy_0_xxyyy_xz[k] = -g_y_0_xyyy_xz[k] - g_xy_0_xyyy_xz[k] * ab_x + g_xy_0_xyyy_xxz[k];

                g_xy_0_xxyyy_yy[k] = -g_y_0_xyyy_yy[k] - g_xy_0_xyyy_yy[k] * ab_x + g_xy_0_xyyy_xyy[k];

                g_xy_0_xxyyy_yz[k] = -g_y_0_xyyy_yz[k] - g_xy_0_xyyy_yz[k] * ab_x + g_xy_0_xyyy_xyz[k];

                g_xy_0_xxyyy_zz[k] = -g_y_0_xyyy_zz[k] - g_xy_0_xyyy_zz[k] * ab_x + g_xy_0_xyyy_xzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_xx, g_xy_0_xxyy_xxz, g_xy_0_xxyy_xy, g_xy_0_xxyy_xyz, g_xy_0_xxyy_xz, g_xy_0_xxyy_xzz, g_xy_0_xxyy_yy, g_xy_0_xxyy_yyz, g_xy_0_xxyy_yz, g_xy_0_xxyy_yzz, g_xy_0_xxyy_zz, g_xy_0_xxyy_zzz, g_xy_0_xxyyz_xx, g_xy_0_xxyyz_xy, g_xy_0_xxyyz_xz, g_xy_0_xxyyz_yy, g_xy_0_xxyyz_yz, g_xy_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyz_xx[k] = -g_xy_0_xxyy_xx[k] * ab_z + g_xy_0_xxyy_xxz[k];

                g_xy_0_xxyyz_xy[k] = -g_xy_0_xxyy_xy[k] * ab_z + g_xy_0_xxyy_xyz[k];

                g_xy_0_xxyyz_xz[k] = -g_xy_0_xxyy_xz[k] * ab_z + g_xy_0_xxyy_xzz[k];

                g_xy_0_xxyyz_yy[k] = -g_xy_0_xxyy_yy[k] * ab_z + g_xy_0_xxyy_yyz[k];

                g_xy_0_xxyyz_yz[k] = -g_xy_0_xxyy_yz[k] * ab_z + g_xy_0_xxyy_yzz[k];

                g_xy_0_xxyyz_zz[k] = -g_xy_0_xxyy_zz[k] * ab_z + g_xy_0_xxyy_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyz_xx, g_xy_0_xxyz_xxz, g_xy_0_xxyz_xy, g_xy_0_xxyz_xyz, g_xy_0_xxyz_xz, g_xy_0_xxyz_xzz, g_xy_0_xxyz_yy, g_xy_0_xxyz_yyz, g_xy_0_xxyz_yz, g_xy_0_xxyz_yzz, g_xy_0_xxyz_zz, g_xy_0_xxyz_zzz, g_xy_0_xxyzz_xx, g_xy_0_xxyzz_xy, g_xy_0_xxyzz_xz, g_xy_0_xxyzz_yy, g_xy_0_xxyzz_yz, g_xy_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzz_xx[k] = -g_xy_0_xxyz_xx[k] * ab_z + g_xy_0_xxyz_xxz[k];

                g_xy_0_xxyzz_xy[k] = -g_xy_0_xxyz_xy[k] * ab_z + g_xy_0_xxyz_xyz[k];

                g_xy_0_xxyzz_xz[k] = -g_xy_0_xxyz_xz[k] * ab_z + g_xy_0_xxyz_xzz[k];

                g_xy_0_xxyzz_yy[k] = -g_xy_0_xxyz_yy[k] * ab_z + g_xy_0_xxyz_yyz[k];

                g_xy_0_xxyzz_yz[k] = -g_xy_0_xxyz_yz[k] * ab_z + g_xy_0_xxyz_yzz[k];

                g_xy_0_xxyzz_zz[k] = -g_xy_0_xxyz_zz[k] * ab_z + g_xy_0_xxyz_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzz_xx, g_xy_0_xxzz_xxz, g_xy_0_xxzz_xy, g_xy_0_xxzz_xyz, g_xy_0_xxzz_xz, g_xy_0_xxzz_xzz, g_xy_0_xxzz_yy, g_xy_0_xxzz_yyz, g_xy_0_xxzz_yz, g_xy_0_xxzz_yzz, g_xy_0_xxzz_zz, g_xy_0_xxzz_zzz, g_xy_0_xxzzz_xx, g_xy_0_xxzzz_xy, g_xy_0_xxzzz_xz, g_xy_0_xxzzz_yy, g_xy_0_xxzzz_yz, g_xy_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzz_xx[k] = -g_xy_0_xxzz_xx[k] * ab_z + g_xy_0_xxzz_xxz[k];

                g_xy_0_xxzzz_xy[k] = -g_xy_0_xxzz_xy[k] * ab_z + g_xy_0_xxzz_xyz[k];

                g_xy_0_xxzzz_xz[k] = -g_xy_0_xxzz_xz[k] * ab_z + g_xy_0_xxzz_xzz[k];

                g_xy_0_xxzzz_yy[k] = -g_xy_0_xxzz_yy[k] * ab_z + g_xy_0_xxzz_yyz[k];

                g_xy_0_xxzzz_yz[k] = -g_xy_0_xxzz_yz[k] * ab_z + g_xy_0_xxzz_yzz[k];

                g_xy_0_xxzzz_zz[k] = -g_xy_0_xxzz_zz[k] * ab_z + g_xy_0_xxzz_zzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_xx, g_xy_0_xyyyy_xy, g_xy_0_xyyyy_xz, g_xy_0_xyyyy_yy, g_xy_0_xyyyy_yz, g_xy_0_xyyyy_zz, g_xy_0_yyyy_xx, g_xy_0_yyyy_xxx, g_xy_0_yyyy_xxy, g_xy_0_yyyy_xxz, g_xy_0_yyyy_xy, g_xy_0_yyyy_xyy, g_xy_0_yyyy_xyz, g_xy_0_yyyy_xz, g_xy_0_yyyy_xzz, g_xy_0_yyyy_yy, g_xy_0_yyyy_yz, g_xy_0_yyyy_zz, g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyy_xx[k] = -g_y_0_yyyy_xx[k] - g_xy_0_yyyy_xx[k] * ab_x + g_xy_0_yyyy_xxx[k];

                g_xy_0_xyyyy_xy[k] = -g_y_0_yyyy_xy[k] - g_xy_0_yyyy_xy[k] * ab_x + g_xy_0_yyyy_xxy[k];

                g_xy_0_xyyyy_xz[k] = -g_y_0_yyyy_xz[k] - g_xy_0_yyyy_xz[k] * ab_x + g_xy_0_yyyy_xxz[k];

                g_xy_0_xyyyy_yy[k] = -g_y_0_yyyy_yy[k] - g_xy_0_yyyy_yy[k] * ab_x + g_xy_0_yyyy_xyy[k];

                g_xy_0_xyyyy_yz[k] = -g_y_0_yyyy_yz[k] - g_xy_0_yyyy_yz[k] * ab_x + g_xy_0_yyyy_xyz[k];

                g_xy_0_xyyyy_zz[k] = -g_y_0_yyyy_zz[k] - g_xy_0_yyyy_zz[k] * ab_x + g_xy_0_yyyy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_xx, g_xy_0_xyyy_xxz, g_xy_0_xyyy_xy, g_xy_0_xyyy_xyz, g_xy_0_xyyy_xz, g_xy_0_xyyy_xzz, g_xy_0_xyyy_yy, g_xy_0_xyyy_yyz, g_xy_0_xyyy_yz, g_xy_0_xyyy_yzz, g_xy_0_xyyy_zz, g_xy_0_xyyy_zzz, g_xy_0_xyyyz_xx, g_xy_0_xyyyz_xy, g_xy_0_xyyyz_xz, g_xy_0_xyyyz_yy, g_xy_0_xyyyz_yz, g_xy_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyz_xx[k] = -g_xy_0_xyyy_xx[k] * ab_z + g_xy_0_xyyy_xxz[k];

                g_xy_0_xyyyz_xy[k] = -g_xy_0_xyyy_xy[k] * ab_z + g_xy_0_xyyy_xyz[k];

                g_xy_0_xyyyz_xz[k] = -g_xy_0_xyyy_xz[k] * ab_z + g_xy_0_xyyy_xzz[k];

                g_xy_0_xyyyz_yy[k] = -g_xy_0_xyyy_yy[k] * ab_z + g_xy_0_xyyy_yyz[k];

                g_xy_0_xyyyz_yz[k] = -g_xy_0_xyyy_yz[k] * ab_z + g_xy_0_xyyy_yzz[k];

                g_xy_0_xyyyz_zz[k] = -g_xy_0_xyyy_zz[k] * ab_z + g_xy_0_xyyy_zzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyz_xx, g_xy_0_xyyz_xxz, g_xy_0_xyyz_xy, g_xy_0_xyyz_xyz, g_xy_0_xyyz_xz, g_xy_0_xyyz_xzz, g_xy_0_xyyz_yy, g_xy_0_xyyz_yyz, g_xy_0_xyyz_yz, g_xy_0_xyyz_yzz, g_xy_0_xyyz_zz, g_xy_0_xyyz_zzz, g_xy_0_xyyzz_xx, g_xy_0_xyyzz_xy, g_xy_0_xyyzz_xz, g_xy_0_xyyzz_yy, g_xy_0_xyyzz_yz, g_xy_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzz_xx[k] = -g_xy_0_xyyz_xx[k] * ab_z + g_xy_0_xyyz_xxz[k];

                g_xy_0_xyyzz_xy[k] = -g_xy_0_xyyz_xy[k] * ab_z + g_xy_0_xyyz_xyz[k];

                g_xy_0_xyyzz_xz[k] = -g_xy_0_xyyz_xz[k] * ab_z + g_xy_0_xyyz_xzz[k];

                g_xy_0_xyyzz_yy[k] = -g_xy_0_xyyz_yy[k] * ab_z + g_xy_0_xyyz_yyz[k];

                g_xy_0_xyyzz_yz[k] = -g_xy_0_xyyz_yz[k] * ab_z + g_xy_0_xyyz_yzz[k];

                g_xy_0_xyyzz_zz[k] = -g_xy_0_xyyz_zz[k] * ab_z + g_xy_0_xyyz_zzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzz_xx, g_xy_0_xyzz_xxz, g_xy_0_xyzz_xy, g_xy_0_xyzz_xyz, g_xy_0_xyzz_xz, g_xy_0_xyzz_xzz, g_xy_0_xyzz_yy, g_xy_0_xyzz_yyz, g_xy_0_xyzz_yz, g_xy_0_xyzz_yzz, g_xy_0_xyzz_zz, g_xy_0_xyzz_zzz, g_xy_0_xyzzz_xx, g_xy_0_xyzzz_xy, g_xy_0_xyzzz_xz, g_xy_0_xyzzz_yy, g_xy_0_xyzzz_yz, g_xy_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzz_xx[k] = -g_xy_0_xyzz_xx[k] * ab_z + g_xy_0_xyzz_xxz[k];

                g_xy_0_xyzzz_xy[k] = -g_xy_0_xyzz_xy[k] * ab_z + g_xy_0_xyzz_xyz[k];

                g_xy_0_xyzzz_xz[k] = -g_xy_0_xyzz_xz[k] * ab_z + g_xy_0_xyzz_xzz[k];

                g_xy_0_xyzzz_yy[k] = -g_xy_0_xyzz_yy[k] * ab_z + g_xy_0_xyzz_yyz[k];

                g_xy_0_xyzzz_yz[k] = -g_xy_0_xyzz_yz[k] * ab_z + g_xy_0_xyzz_yzz[k];

                g_xy_0_xyzzz_zz[k] = -g_xy_0_xyzz_zz[k] * ab_z + g_xy_0_xyzz_zzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzz_xx, g_xy_0_xzzz_xxz, g_xy_0_xzzz_xy, g_xy_0_xzzz_xyz, g_xy_0_xzzz_xz, g_xy_0_xzzz_xzz, g_xy_0_xzzz_yy, g_xy_0_xzzz_yyz, g_xy_0_xzzz_yz, g_xy_0_xzzz_yzz, g_xy_0_xzzz_zz, g_xy_0_xzzz_zzz, g_xy_0_xzzzz_xx, g_xy_0_xzzzz_xy, g_xy_0_xzzzz_xz, g_xy_0_xzzzz_yy, g_xy_0_xzzzz_yz, g_xy_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzz_xx[k] = -g_xy_0_xzzz_xx[k] * ab_z + g_xy_0_xzzz_xxz[k];

                g_xy_0_xzzzz_xy[k] = -g_xy_0_xzzz_xy[k] * ab_z + g_xy_0_xzzz_xyz[k];

                g_xy_0_xzzzz_xz[k] = -g_xy_0_xzzz_xz[k] * ab_z + g_xy_0_xzzz_xzz[k];

                g_xy_0_xzzzz_yy[k] = -g_xy_0_xzzz_yy[k] * ab_z + g_xy_0_xzzz_yyz[k];

                g_xy_0_xzzzz_yz[k] = -g_xy_0_xzzz_yz[k] * ab_z + g_xy_0_xzzz_yzz[k];

                g_xy_0_xzzzz_zz[k] = -g_xy_0_xzzz_zz[k] * ab_z + g_xy_0_xzzz_zzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_xx, g_x_0_yyyy_xy, g_x_0_yyyy_xz, g_x_0_yyyy_yy, g_x_0_yyyy_yz, g_x_0_yyyy_zz, g_xy_0_yyyy_xx, g_xy_0_yyyy_xxy, g_xy_0_yyyy_xy, g_xy_0_yyyy_xyy, g_xy_0_yyyy_xyz, g_xy_0_yyyy_xz, g_xy_0_yyyy_yy, g_xy_0_yyyy_yyy, g_xy_0_yyyy_yyz, g_xy_0_yyyy_yz, g_xy_0_yyyy_yzz, g_xy_0_yyyy_zz, g_xy_0_yyyyy_xx, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_yy, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyy_xx[k] = -g_x_0_yyyy_xx[k] - g_xy_0_yyyy_xx[k] * ab_y + g_xy_0_yyyy_xxy[k];

                g_xy_0_yyyyy_xy[k] = -g_x_0_yyyy_xy[k] - g_xy_0_yyyy_xy[k] * ab_y + g_xy_0_yyyy_xyy[k];

                g_xy_0_yyyyy_xz[k] = -g_x_0_yyyy_xz[k] - g_xy_0_yyyy_xz[k] * ab_y + g_xy_0_yyyy_xyz[k];

                g_xy_0_yyyyy_yy[k] = -g_x_0_yyyy_yy[k] - g_xy_0_yyyy_yy[k] * ab_y + g_xy_0_yyyy_yyy[k];

                g_xy_0_yyyyy_yz[k] = -g_x_0_yyyy_yz[k] - g_xy_0_yyyy_yz[k] * ab_y + g_xy_0_yyyy_yyz[k];

                g_xy_0_yyyyy_zz[k] = -g_x_0_yyyy_zz[k] - g_xy_0_yyyy_zz[k] * ab_y + g_xy_0_yyyy_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyy_xx, g_xy_0_yyyy_xxz, g_xy_0_yyyy_xy, g_xy_0_yyyy_xyz, g_xy_0_yyyy_xz, g_xy_0_yyyy_xzz, g_xy_0_yyyy_yy, g_xy_0_yyyy_yyz, g_xy_0_yyyy_yz, g_xy_0_yyyy_yzz, g_xy_0_yyyy_zz, g_xy_0_yyyy_zzz, g_xy_0_yyyyz_xx, g_xy_0_yyyyz_xy, g_xy_0_yyyyz_xz, g_xy_0_yyyyz_yy, g_xy_0_yyyyz_yz, g_xy_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyz_xx[k] = -g_xy_0_yyyy_xx[k] * ab_z + g_xy_0_yyyy_xxz[k];

                g_xy_0_yyyyz_xy[k] = -g_xy_0_yyyy_xy[k] * ab_z + g_xy_0_yyyy_xyz[k];

                g_xy_0_yyyyz_xz[k] = -g_xy_0_yyyy_xz[k] * ab_z + g_xy_0_yyyy_xzz[k];

                g_xy_0_yyyyz_yy[k] = -g_xy_0_yyyy_yy[k] * ab_z + g_xy_0_yyyy_yyz[k];

                g_xy_0_yyyyz_yz[k] = -g_xy_0_yyyy_yz[k] * ab_z + g_xy_0_yyyy_yzz[k];

                g_xy_0_yyyyz_zz[k] = -g_xy_0_yyyy_zz[k] * ab_z + g_xy_0_yyyy_zzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyz_xx, g_xy_0_yyyz_xxz, g_xy_0_yyyz_xy, g_xy_0_yyyz_xyz, g_xy_0_yyyz_xz, g_xy_0_yyyz_xzz, g_xy_0_yyyz_yy, g_xy_0_yyyz_yyz, g_xy_0_yyyz_yz, g_xy_0_yyyz_yzz, g_xy_0_yyyz_zz, g_xy_0_yyyz_zzz, g_xy_0_yyyzz_xx, g_xy_0_yyyzz_xy, g_xy_0_yyyzz_xz, g_xy_0_yyyzz_yy, g_xy_0_yyyzz_yz, g_xy_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzz_xx[k] = -g_xy_0_yyyz_xx[k] * ab_z + g_xy_0_yyyz_xxz[k];

                g_xy_0_yyyzz_xy[k] = -g_xy_0_yyyz_xy[k] * ab_z + g_xy_0_yyyz_xyz[k];

                g_xy_0_yyyzz_xz[k] = -g_xy_0_yyyz_xz[k] * ab_z + g_xy_0_yyyz_xzz[k];

                g_xy_0_yyyzz_yy[k] = -g_xy_0_yyyz_yy[k] * ab_z + g_xy_0_yyyz_yyz[k];

                g_xy_0_yyyzz_yz[k] = -g_xy_0_yyyz_yz[k] * ab_z + g_xy_0_yyyz_yzz[k];

                g_xy_0_yyyzz_zz[k] = -g_xy_0_yyyz_zz[k] * ab_z + g_xy_0_yyyz_zzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzz_xx, g_xy_0_yyzz_xxz, g_xy_0_yyzz_xy, g_xy_0_yyzz_xyz, g_xy_0_yyzz_xz, g_xy_0_yyzz_xzz, g_xy_0_yyzz_yy, g_xy_0_yyzz_yyz, g_xy_0_yyzz_yz, g_xy_0_yyzz_yzz, g_xy_0_yyzz_zz, g_xy_0_yyzz_zzz, g_xy_0_yyzzz_xx, g_xy_0_yyzzz_xy, g_xy_0_yyzzz_xz, g_xy_0_yyzzz_yy, g_xy_0_yyzzz_yz, g_xy_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzz_xx[k] = -g_xy_0_yyzz_xx[k] * ab_z + g_xy_0_yyzz_xxz[k];

                g_xy_0_yyzzz_xy[k] = -g_xy_0_yyzz_xy[k] * ab_z + g_xy_0_yyzz_xyz[k];

                g_xy_0_yyzzz_xz[k] = -g_xy_0_yyzz_xz[k] * ab_z + g_xy_0_yyzz_xzz[k];

                g_xy_0_yyzzz_yy[k] = -g_xy_0_yyzz_yy[k] * ab_z + g_xy_0_yyzz_yyz[k];

                g_xy_0_yyzzz_yz[k] = -g_xy_0_yyzz_yz[k] * ab_z + g_xy_0_yyzz_yzz[k];

                g_xy_0_yyzzz_zz[k] = -g_xy_0_yyzz_zz[k] * ab_z + g_xy_0_yyzz_zzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzz_xx, g_xy_0_yzzz_xxz, g_xy_0_yzzz_xy, g_xy_0_yzzz_xyz, g_xy_0_yzzz_xz, g_xy_0_yzzz_xzz, g_xy_0_yzzz_yy, g_xy_0_yzzz_yyz, g_xy_0_yzzz_yz, g_xy_0_yzzz_yzz, g_xy_0_yzzz_zz, g_xy_0_yzzz_zzz, g_xy_0_yzzzz_xx, g_xy_0_yzzzz_xy, g_xy_0_yzzzz_xz, g_xy_0_yzzzz_yy, g_xy_0_yzzzz_yz, g_xy_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzz_xx[k] = -g_xy_0_yzzz_xx[k] * ab_z + g_xy_0_yzzz_xxz[k];

                g_xy_0_yzzzz_xy[k] = -g_xy_0_yzzz_xy[k] * ab_z + g_xy_0_yzzz_xyz[k];

                g_xy_0_yzzzz_xz[k] = -g_xy_0_yzzz_xz[k] * ab_z + g_xy_0_yzzz_xzz[k];

                g_xy_0_yzzzz_yy[k] = -g_xy_0_yzzz_yy[k] * ab_z + g_xy_0_yzzz_yyz[k];

                g_xy_0_yzzzz_yz[k] = -g_xy_0_yzzz_yz[k] * ab_z + g_xy_0_yzzz_yzz[k];

                g_xy_0_yzzzz_zz[k] = -g_xy_0_yzzz_zz[k] * ab_z + g_xy_0_yzzz_zzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzz_xx, g_xy_0_zzzz_xxz, g_xy_0_zzzz_xy, g_xy_0_zzzz_xyz, g_xy_0_zzzz_xz, g_xy_0_zzzz_xzz, g_xy_0_zzzz_yy, g_xy_0_zzzz_yyz, g_xy_0_zzzz_yz, g_xy_0_zzzz_yzz, g_xy_0_zzzz_zz, g_xy_0_zzzz_zzz, g_xy_0_zzzzz_xx, g_xy_0_zzzzz_xy, g_xy_0_zzzzz_xz, g_xy_0_zzzzz_yy, g_xy_0_zzzzz_yz, g_xy_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzz_xx[k] = -g_xy_0_zzzz_xx[k] * ab_z + g_xy_0_zzzz_xxz[k];

                g_xy_0_zzzzz_xy[k] = -g_xy_0_zzzz_xy[k] * ab_z + g_xy_0_zzzz_xyz[k];

                g_xy_0_zzzzz_xz[k] = -g_xy_0_zzzz_xz[k] * ab_z + g_xy_0_zzzz_xzz[k];

                g_xy_0_zzzzz_yy[k] = -g_xy_0_zzzz_yy[k] * ab_z + g_xy_0_zzzz_yyz[k];

                g_xy_0_zzzzz_yz[k] = -g_xy_0_zzzz_yz[k] * ab_z + g_xy_0_zzzz_yzz[k];

                g_xy_0_zzzzz_zz[k] = -g_xy_0_zzzz_zz[k] * ab_z + g_xy_0_zzzz_zzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_xx, g_xz_0_xxxx_xxx, g_xz_0_xxxx_xxy, g_xz_0_xxxx_xxz, g_xz_0_xxxx_xy, g_xz_0_xxxx_xyy, g_xz_0_xxxx_xyz, g_xz_0_xxxx_xz, g_xz_0_xxxx_xzz, g_xz_0_xxxx_yy, g_xz_0_xxxx_yz, g_xz_0_xxxx_zz, g_xz_0_xxxxx_xx, g_xz_0_xxxxx_xy, g_xz_0_xxxxx_xz, g_xz_0_xxxxx_yy, g_xz_0_xxxxx_yz, g_xz_0_xxxxx_zz, g_z_0_xxxx_xx, g_z_0_xxxx_xy, g_z_0_xxxx_xz, g_z_0_xxxx_yy, g_z_0_xxxx_yz, g_z_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxx_xx[k] = -g_z_0_xxxx_xx[k] - g_xz_0_xxxx_xx[k] * ab_x + g_xz_0_xxxx_xxx[k];

                g_xz_0_xxxxx_xy[k] = -g_z_0_xxxx_xy[k] - g_xz_0_xxxx_xy[k] * ab_x + g_xz_0_xxxx_xxy[k];

                g_xz_0_xxxxx_xz[k] = -g_z_0_xxxx_xz[k] - g_xz_0_xxxx_xz[k] * ab_x + g_xz_0_xxxx_xxz[k];

                g_xz_0_xxxxx_yy[k] = -g_z_0_xxxx_yy[k] - g_xz_0_xxxx_yy[k] * ab_x + g_xz_0_xxxx_xyy[k];

                g_xz_0_xxxxx_yz[k] = -g_z_0_xxxx_yz[k] - g_xz_0_xxxx_yz[k] * ab_x + g_xz_0_xxxx_xyz[k];

                g_xz_0_xxxxx_zz[k] = -g_z_0_xxxx_zz[k] - g_xz_0_xxxx_zz[k] * ab_x + g_xz_0_xxxx_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxx_xx, g_xz_0_xxxx_xxy, g_xz_0_xxxx_xy, g_xz_0_xxxx_xyy, g_xz_0_xxxx_xyz, g_xz_0_xxxx_xz, g_xz_0_xxxx_yy, g_xz_0_xxxx_yyy, g_xz_0_xxxx_yyz, g_xz_0_xxxx_yz, g_xz_0_xxxx_yzz, g_xz_0_xxxx_zz, g_xz_0_xxxxy_xx, g_xz_0_xxxxy_xy, g_xz_0_xxxxy_xz, g_xz_0_xxxxy_yy, g_xz_0_xxxxy_yz, g_xz_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxy_xx[k] = -g_xz_0_xxxx_xx[k] * ab_y + g_xz_0_xxxx_xxy[k];

                g_xz_0_xxxxy_xy[k] = -g_xz_0_xxxx_xy[k] * ab_y + g_xz_0_xxxx_xyy[k];

                g_xz_0_xxxxy_xz[k] = -g_xz_0_xxxx_xz[k] * ab_y + g_xz_0_xxxx_xyz[k];

                g_xz_0_xxxxy_yy[k] = -g_xz_0_xxxx_yy[k] * ab_y + g_xz_0_xxxx_yyy[k];

                g_xz_0_xxxxy_yz[k] = -g_xz_0_xxxx_yz[k] * ab_y + g_xz_0_xxxx_yyz[k];

                g_xz_0_xxxxy_zz[k] = -g_xz_0_xxxx_zz[k] * ab_y + g_xz_0_xxxx_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxz_xx, g_xz_0_xxxxz_xy, g_xz_0_xxxxz_xz, g_xz_0_xxxxz_yy, g_xz_0_xxxxz_yz, g_xz_0_xxxxz_zz, g_xz_0_xxxz_xx, g_xz_0_xxxz_xxx, g_xz_0_xxxz_xxy, g_xz_0_xxxz_xxz, g_xz_0_xxxz_xy, g_xz_0_xxxz_xyy, g_xz_0_xxxz_xyz, g_xz_0_xxxz_xz, g_xz_0_xxxz_xzz, g_xz_0_xxxz_yy, g_xz_0_xxxz_yz, g_xz_0_xxxz_zz, g_z_0_xxxz_xx, g_z_0_xxxz_xy, g_z_0_xxxz_xz, g_z_0_xxxz_yy, g_z_0_xxxz_yz, g_z_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxz_xx[k] = -g_z_0_xxxz_xx[k] - g_xz_0_xxxz_xx[k] * ab_x + g_xz_0_xxxz_xxx[k];

                g_xz_0_xxxxz_xy[k] = -g_z_0_xxxz_xy[k] - g_xz_0_xxxz_xy[k] * ab_x + g_xz_0_xxxz_xxy[k];

                g_xz_0_xxxxz_xz[k] = -g_z_0_xxxz_xz[k] - g_xz_0_xxxz_xz[k] * ab_x + g_xz_0_xxxz_xxz[k];

                g_xz_0_xxxxz_yy[k] = -g_z_0_xxxz_yy[k] - g_xz_0_xxxz_yy[k] * ab_x + g_xz_0_xxxz_xyy[k];

                g_xz_0_xxxxz_yz[k] = -g_z_0_xxxz_yz[k] - g_xz_0_xxxz_yz[k] * ab_x + g_xz_0_xxxz_xyz[k];

                g_xz_0_xxxxz_zz[k] = -g_z_0_xxxz_zz[k] - g_xz_0_xxxz_zz[k] * ab_x + g_xz_0_xxxz_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxy_xx, g_xz_0_xxxy_xxy, g_xz_0_xxxy_xy, g_xz_0_xxxy_xyy, g_xz_0_xxxy_xyz, g_xz_0_xxxy_xz, g_xz_0_xxxy_yy, g_xz_0_xxxy_yyy, g_xz_0_xxxy_yyz, g_xz_0_xxxy_yz, g_xz_0_xxxy_yzz, g_xz_0_xxxy_zz, g_xz_0_xxxyy_xx, g_xz_0_xxxyy_xy, g_xz_0_xxxyy_xz, g_xz_0_xxxyy_yy, g_xz_0_xxxyy_yz, g_xz_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyy_xx[k] = -g_xz_0_xxxy_xx[k] * ab_y + g_xz_0_xxxy_xxy[k];

                g_xz_0_xxxyy_xy[k] = -g_xz_0_xxxy_xy[k] * ab_y + g_xz_0_xxxy_xyy[k];

                g_xz_0_xxxyy_xz[k] = -g_xz_0_xxxy_xz[k] * ab_y + g_xz_0_xxxy_xyz[k];

                g_xz_0_xxxyy_yy[k] = -g_xz_0_xxxy_yy[k] * ab_y + g_xz_0_xxxy_yyy[k];

                g_xz_0_xxxyy_yz[k] = -g_xz_0_xxxy_yz[k] * ab_y + g_xz_0_xxxy_yyz[k];

                g_xz_0_xxxyy_zz[k] = -g_xz_0_xxxy_zz[k] * ab_y + g_xz_0_xxxy_yzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyz_xx, g_xz_0_xxxyz_xy, g_xz_0_xxxyz_xz, g_xz_0_xxxyz_yy, g_xz_0_xxxyz_yz, g_xz_0_xxxyz_zz, g_xz_0_xxxz_xx, g_xz_0_xxxz_xxy, g_xz_0_xxxz_xy, g_xz_0_xxxz_xyy, g_xz_0_xxxz_xyz, g_xz_0_xxxz_xz, g_xz_0_xxxz_yy, g_xz_0_xxxz_yyy, g_xz_0_xxxz_yyz, g_xz_0_xxxz_yz, g_xz_0_xxxz_yzz, g_xz_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyz_xx[k] = -g_xz_0_xxxz_xx[k] * ab_y + g_xz_0_xxxz_xxy[k];

                g_xz_0_xxxyz_xy[k] = -g_xz_0_xxxz_xy[k] * ab_y + g_xz_0_xxxz_xyy[k];

                g_xz_0_xxxyz_xz[k] = -g_xz_0_xxxz_xz[k] * ab_y + g_xz_0_xxxz_xyz[k];

                g_xz_0_xxxyz_yy[k] = -g_xz_0_xxxz_yy[k] * ab_y + g_xz_0_xxxz_yyy[k];

                g_xz_0_xxxyz_yz[k] = -g_xz_0_xxxz_yz[k] * ab_y + g_xz_0_xxxz_yyz[k];

                g_xz_0_xxxyz_zz[k] = -g_xz_0_xxxz_zz[k] * ab_y + g_xz_0_xxxz_yzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzz_xx, g_xz_0_xxxzz_xy, g_xz_0_xxxzz_xz, g_xz_0_xxxzz_yy, g_xz_0_xxxzz_yz, g_xz_0_xxxzz_zz, g_xz_0_xxzz_xx, g_xz_0_xxzz_xxx, g_xz_0_xxzz_xxy, g_xz_0_xxzz_xxz, g_xz_0_xxzz_xy, g_xz_0_xxzz_xyy, g_xz_0_xxzz_xyz, g_xz_0_xxzz_xz, g_xz_0_xxzz_xzz, g_xz_0_xxzz_yy, g_xz_0_xxzz_yz, g_xz_0_xxzz_zz, g_z_0_xxzz_xx, g_z_0_xxzz_xy, g_z_0_xxzz_xz, g_z_0_xxzz_yy, g_z_0_xxzz_yz, g_z_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzz_xx[k] = -g_z_0_xxzz_xx[k] - g_xz_0_xxzz_xx[k] * ab_x + g_xz_0_xxzz_xxx[k];

                g_xz_0_xxxzz_xy[k] = -g_z_0_xxzz_xy[k] - g_xz_0_xxzz_xy[k] * ab_x + g_xz_0_xxzz_xxy[k];

                g_xz_0_xxxzz_xz[k] = -g_z_0_xxzz_xz[k] - g_xz_0_xxzz_xz[k] * ab_x + g_xz_0_xxzz_xxz[k];

                g_xz_0_xxxzz_yy[k] = -g_z_0_xxzz_yy[k] - g_xz_0_xxzz_yy[k] * ab_x + g_xz_0_xxzz_xyy[k];

                g_xz_0_xxxzz_yz[k] = -g_z_0_xxzz_yz[k] - g_xz_0_xxzz_yz[k] * ab_x + g_xz_0_xxzz_xyz[k];

                g_xz_0_xxxzz_zz[k] = -g_z_0_xxzz_zz[k] - g_xz_0_xxzz_zz[k] * ab_x + g_xz_0_xxzz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyy_xx, g_xz_0_xxyy_xxy, g_xz_0_xxyy_xy, g_xz_0_xxyy_xyy, g_xz_0_xxyy_xyz, g_xz_0_xxyy_xz, g_xz_0_xxyy_yy, g_xz_0_xxyy_yyy, g_xz_0_xxyy_yyz, g_xz_0_xxyy_yz, g_xz_0_xxyy_yzz, g_xz_0_xxyy_zz, g_xz_0_xxyyy_xx, g_xz_0_xxyyy_xy, g_xz_0_xxyyy_xz, g_xz_0_xxyyy_yy, g_xz_0_xxyyy_yz, g_xz_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyy_xx[k] = -g_xz_0_xxyy_xx[k] * ab_y + g_xz_0_xxyy_xxy[k];

                g_xz_0_xxyyy_xy[k] = -g_xz_0_xxyy_xy[k] * ab_y + g_xz_0_xxyy_xyy[k];

                g_xz_0_xxyyy_xz[k] = -g_xz_0_xxyy_xz[k] * ab_y + g_xz_0_xxyy_xyz[k];

                g_xz_0_xxyyy_yy[k] = -g_xz_0_xxyy_yy[k] * ab_y + g_xz_0_xxyy_yyy[k];

                g_xz_0_xxyyy_yz[k] = -g_xz_0_xxyy_yz[k] * ab_y + g_xz_0_xxyy_yyz[k];

                g_xz_0_xxyyy_zz[k] = -g_xz_0_xxyy_zz[k] * ab_y + g_xz_0_xxyy_yzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyz_xx, g_xz_0_xxyyz_xy, g_xz_0_xxyyz_xz, g_xz_0_xxyyz_yy, g_xz_0_xxyyz_yz, g_xz_0_xxyyz_zz, g_xz_0_xxyz_xx, g_xz_0_xxyz_xxy, g_xz_0_xxyz_xy, g_xz_0_xxyz_xyy, g_xz_0_xxyz_xyz, g_xz_0_xxyz_xz, g_xz_0_xxyz_yy, g_xz_0_xxyz_yyy, g_xz_0_xxyz_yyz, g_xz_0_xxyz_yz, g_xz_0_xxyz_yzz, g_xz_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyz_xx[k] = -g_xz_0_xxyz_xx[k] * ab_y + g_xz_0_xxyz_xxy[k];

                g_xz_0_xxyyz_xy[k] = -g_xz_0_xxyz_xy[k] * ab_y + g_xz_0_xxyz_xyy[k];

                g_xz_0_xxyyz_xz[k] = -g_xz_0_xxyz_xz[k] * ab_y + g_xz_0_xxyz_xyz[k];

                g_xz_0_xxyyz_yy[k] = -g_xz_0_xxyz_yy[k] * ab_y + g_xz_0_xxyz_yyy[k];

                g_xz_0_xxyyz_yz[k] = -g_xz_0_xxyz_yz[k] * ab_y + g_xz_0_xxyz_yyz[k];

                g_xz_0_xxyyz_zz[k] = -g_xz_0_xxyz_zz[k] * ab_y + g_xz_0_xxyz_yzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzz_xx, g_xz_0_xxyzz_xy, g_xz_0_xxyzz_xz, g_xz_0_xxyzz_yy, g_xz_0_xxyzz_yz, g_xz_0_xxyzz_zz, g_xz_0_xxzz_xx, g_xz_0_xxzz_xxy, g_xz_0_xxzz_xy, g_xz_0_xxzz_xyy, g_xz_0_xxzz_xyz, g_xz_0_xxzz_xz, g_xz_0_xxzz_yy, g_xz_0_xxzz_yyy, g_xz_0_xxzz_yyz, g_xz_0_xxzz_yz, g_xz_0_xxzz_yzz, g_xz_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzz_xx[k] = -g_xz_0_xxzz_xx[k] * ab_y + g_xz_0_xxzz_xxy[k];

                g_xz_0_xxyzz_xy[k] = -g_xz_0_xxzz_xy[k] * ab_y + g_xz_0_xxzz_xyy[k];

                g_xz_0_xxyzz_xz[k] = -g_xz_0_xxzz_xz[k] * ab_y + g_xz_0_xxzz_xyz[k];

                g_xz_0_xxyzz_yy[k] = -g_xz_0_xxzz_yy[k] * ab_y + g_xz_0_xxzz_yyy[k];

                g_xz_0_xxyzz_yz[k] = -g_xz_0_xxzz_yz[k] * ab_y + g_xz_0_xxzz_yyz[k];

                g_xz_0_xxyzz_zz[k] = -g_xz_0_xxzz_zz[k] * ab_y + g_xz_0_xxzz_yzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzz_xx, g_xz_0_xxzzz_xy, g_xz_0_xxzzz_xz, g_xz_0_xxzzz_yy, g_xz_0_xxzzz_yz, g_xz_0_xxzzz_zz, g_xz_0_xzzz_xx, g_xz_0_xzzz_xxx, g_xz_0_xzzz_xxy, g_xz_0_xzzz_xxz, g_xz_0_xzzz_xy, g_xz_0_xzzz_xyy, g_xz_0_xzzz_xyz, g_xz_0_xzzz_xz, g_xz_0_xzzz_xzz, g_xz_0_xzzz_yy, g_xz_0_xzzz_yz, g_xz_0_xzzz_zz, g_z_0_xzzz_xx, g_z_0_xzzz_xy, g_z_0_xzzz_xz, g_z_0_xzzz_yy, g_z_0_xzzz_yz, g_z_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzz_xx[k] = -g_z_0_xzzz_xx[k] - g_xz_0_xzzz_xx[k] * ab_x + g_xz_0_xzzz_xxx[k];

                g_xz_0_xxzzz_xy[k] = -g_z_0_xzzz_xy[k] - g_xz_0_xzzz_xy[k] * ab_x + g_xz_0_xzzz_xxy[k];

                g_xz_0_xxzzz_xz[k] = -g_z_0_xzzz_xz[k] - g_xz_0_xzzz_xz[k] * ab_x + g_xz_0_xzzz_xxz[k];

                g_xz_0_xxzzz_yy[k] = -g_z_0_xzzz_yy[k] - g_xz_0_xzzz_yy[k] * ab_x + g_xz_0_xzzz_xyy[k];

                g_xz_0_xxzzz_yz[k] = -g_z_0_xzzz_yz[k] - g_xz_0_xzzz_yz[k] * ab_x + g_xz_0_xzzz_xyz[k];

                g_xz_0_xxzzz_zz[k] = -g_z_0_xzzz_zz[k] - g_xz_0_xzzz_zz[k] * ab_x + g_xz_0_xzzz_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyy_xx, g_xz_0_xyyy_xxy, g_xz_0_xyyy_xy, g_xz_0_xyyy_xyy, g_xz_0_xyyy_xyz, g_xz_0_xyyy_xz, g_xz_0_xyyy_yy, g_xz_0_xyyy_yyy, g_xz_0_xyyy_yyz, g_xz_0_xyyy_yz, g_xz_0_xyyy_yzz, g_xz_0_xyyy_zz, g_xz_0_xyyyy_xx, g_xz_0_xyyyy_xy, g_xz_0_xyyyy_xz, g_xz_0_xyyyy_yy, g_xz_0_xyyyy_yz, g_xz_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyy_xx[k] = -g_xz_0_xyyy_xx[k] * ab_y + g_xz_0_xyyy_xxy[k];

                g_xz_0_xyyyy_xy[k] = -g_xz_0_xyyy_xy[k] * ab_y + g_xz_0_xyyy_xyy[k];

                g_xz_0_xyyyy_xz[k] = -g_xz_0_xyyy_xz[k] * ab_y + g_xz_0_xyyy_xyz[k];

                g_xz_0_xyyyy_yy[k] = -g_xz_0_xyyy_yy[k] * ab_y + g_xz_0_xyyy_yyy[k];

                g_xz_0_xyyyy_yz[k] = -g_xz_0_xyyy_yz[k] * ab_y + g_xz_0_xyyy_yyz[k];

                g_xz_0_xyyyy_zz[k] = -g_xz_0_xyyy_zz[k] * ab_y + g_xz_0_xyyy_yzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyz_xx, g_xz_0_xyyyz_xy, g_xz_0_xyyyz_xz, g_xz_0_xyyyz_yy, g_xz_0_xyyyz_yz, g_xz_0_xyyyz_zz, g_xz_0_xyyz_xx, g_xz_0_xyyz_xxy, g_xz_0_xyyz_xy, g_xz_0_xyyz_xyy, g_xz_0_xyyz_xyz, g_xz_0_xyyz_xz, g_xz_0_xyyz_yy, g_xz_0_xyyz_yyy, g_xz_0_xyyz_yyz, g_xz_0_xyyz_yz, g_xz_0_xyyz_yzz, g_xz_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyz_xx[k] = -g_xz_0_xyyz_xx[k] * ab_y + g_xz_0_xyyz_xxy[k];

                g_xz_0_xyyyz_xy[k] = -g_xz_0_xyyz_xy[k] * ab_y + g_xz_0_xyyz_xyy[k];

                g_xz_0_xyyyz_xz[k] = -g_xz_0_xyyz_xz[k] * ab_y + g_xz_0_xyyz_xyz[k];

                g_xz_0_xyyyz_yy[k] = -g_xz_0_xyyz_yy[k] * ab_y + g_xz_0_xyyz_yyy[k];

                g_xz_0_xyyyz_yz[k] = -g_xz_0_xyyz_yz[k] * ab_y + g_xz_0_xyyz_yyz[k];

                g_xz_0_xyyyz_zz[k] = -g_xz_0_xyyz_zz[k] * ab_y + g_xz_0_xyyz_yzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzz_xx, g_xz_0_xyyzz_xy, g_xz_0_xyyzz_xz, g_xz_0_xyyzz_yy, g_xz_0_xyyzz_yz, g_xz_0_xyyzz_zz, g_xz_0_xyzz_xx, g_xz_0_xyzz_xxy, g_xz_0_xyzz_xy, g_xz_0_xyzz_xyy, g_xz_0_xyzz_xyz, g_xz_0_xyzz_xz, g_xz_0_xyzz_yy, g_xz_0_xyzz_yyy, g_xz_0_xyzz_yyz, g_xz_0_xyzz_yz, g_xz_0_xyzz_yzz, g_xz_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzz_xx[k] = -g_xz_0_xyzz_xx[k] * ab_y + g_xz_0_xyzz_xxy[k];

                g_xz_0_xyyzz_xy[k] = -g_xz_0_xyzz_xy[k] * ab_y + g_xz_0_xyzz_xyy[k];

                g_xz_0_xyyzz_xz[k] = -g_xz_0_xyzz_xz[k] * ab_y + g_xz_0_xyzz_xyz[k];

                g_xz_0_xyyzz_yy[k] = -g_xz_0_xyzz_yy[k] * ab_y + g_xz_0_xyzz_yyy[k];

                g_xz_0_xyyzz_yz[k] = -g_xz_0_xyzz_yz[k] * ab_y + g_xz_0_xyzz_yyz[k];

                g_xz_0_xyyzz_zz[k] = -g_xz_0_xyzz_zz[k] * ab_y + g_xz_0_xyzz_yzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzz_xx, g_xz_0_xyzzz_xy, g_xz_0_xyzzz_xz, g_xz_0_xyzzz_yy, g_xz_0_xyzzz_yz, g_xz_0_xyzzz_zz, g_xz_0_xzzz_xx, g_xz_0_xzzz_xxy, g_xz_0_xzzz_xy, g_xz_0_xzzz_xyy, g_xz_0_xzzz_xyz, g_xz_0_xzzz_xz, g_xz_0_xzzz_yy, g_xz_0_xzzz_yyy, g_xz_0_xzzz_yyz, g_xz_0_xzzz_yz, g_xz_0_xzzz_yzz, g_xz_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzz_xx[k] = -g_xz_0_xzzz_xx[k] * ab_y + g_xz_0_xzzz_xxy[k];

                g_xz_0_xyzzz_xy[k] = -g_xz_0_xzzz_xy[k] * ab_y + g_xz_0_xzzz_xyy[k];

                g_xz_0_xyzzz_xz[k] = -g_xz_0_xzzz_xz[k] * ab_y + g_xz_0_xzzz_xyz[k];

                g_xz_0_xyzzz_yy[k] = -g_xz_0_xzzz_yy[k] * ab_y + g_xz_0_xzzz_yyy[k];

                g_xz_0_xyzzz_yz[k] = -g_xz_0_xzzz_yz[k] * ab_y + g_xz_0_xzzz_yyz[k];

                g_xz_0_xyzzz_zz[k] = -g_xz_0_xzzz_zz[k] * ab_y + g_xz_0_xzzz_yzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzz_xx, g_xz_0_xzzzz_xy, g_xz_0_xzzzz_xz, g_xz_0_xzzzz_yy, g_xz_0_xzzzz_yz, g_xz_0_xzzzz_zz, g_xz_0_zzzz_xx, g_xz_0_zzzz_xxx, g_xz_0_zzzz_xxy, g_xz_0_zzzz_xxz, g_xz_0_zzzz_xy, g_xz_0_zzzz_xyy, g_xz_0_zzzz_xyz, g_xz_0_zzzz_xz, g_xz_0_zzzz_xzz, g_xz_0_zzzz_yy, g_xz_0_zzzz_yz, g_xz_0_zzzz_zz, g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzz_xx[k] = -g_z_0_zzzz_xx[k] - g_xz_0_zzzz_xx[k] * ab_x + g_xz_0_zzzz_xxx[k];

                g_xz_0_xzzzz_xy[k] = -g_z_0_zzzz_xy[k] - g_xz_0_zzzz_xy[k] * ab_x + g_xz_0_zzzz_xxy[k];

                g_xz_0_xzzzz_xz[k] = -g_z_0_zzzz_xz[k] - g_xz_0_zzzz_xz[k] * ab_x + g_xz_0_zzzz_xxz[k];

                g_xz_0_xzzzz_yy[k] = -g_z_0_zzzz_yy[k] - g_xz_0_zzzz_yy[k] * ab_x + g_xz_0_zzzz_xyy[k];

                g_xz_0_xzzzz_yz[k] = -g_z_0_zzzz_yz[k] - g_xz_0_zzzz_yz[k] * ab_x + g_xz_0_zzzz_xyz[k];

                g_xz_0_xzzzz_zz[k] = -g_z_0_zzzz_zz[k] - g_xz_0_zzzz_zz[k] * ab_x + g_xz_0_zzzz_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyy_xx, g_xz_0_yyyy_xxy, g_xz_0_yyyy_xy, g_xz_0_yyyy_xyy, g_xz_0_yyyy_xyz, g_xz_0_yyyy_xz, g_xz_0_yyyy_yy, g_xz_0_yyyy_yyy, g_xz_0_yyyy_yyz, g_xz_0_yyyy_yz, g_xz_0_yyyy_yzz, g_xz_0_yyyy_zz, g_xz_0_yyyyy_xx, g_xz_0_yyyyy_xy, g_xz_0_yyyyy_xz, g_xz_0_yyyyy_yy, g_xz_0_yyyyy_yz, g_xz_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyy_xx[k] = -g_xz_0_yyyy_xx[k] * ab_y + g_xz_0_yyyy_xxy[k];

                g_xz_0_yyyyy_xy[k] = -g_xz_0_yyyy_xy[k] * ab_y + g_xz_0_yyyy_xyy[k];

                g_xz_0_yyyyy_xz[k] = -g_xz_0_yyyy_xz[k] * ab_y + g_xz_0_yyyy_xyz[k];

                g_xz_0_yyyyy_yy[k] = -g_xz_0_yyyy_yy[k] * ab_y + g_xz_0_yyyy_yyy[k];

                g_xz_0_yyyyy_yz[k] = -g_xz_0_yyyy_yz[k] * ab_y + g_xz_0_yyyy_yyz[k];

                g_xz_0_yyyyy_zz[k] = -g_xz_0_yyyy_zz[k] * ab_y + g_xz_0_yyyy_yzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyz_xx, g_xz_0_yyyyz_xy, g_xz_0_yyyyz_xz, g_xz_0_yyyyz_yy, g_xz_0_yyyyz_yz, g_xz_0_yyyyz_zz, g_xz_0_yyyz_xx, g_xz_0_yyyz_xxy, g_xz_0_yyyz_xy, g_xz_0_yyyz_xyy, g_xz_0_yyyz_xyz, g_xz_0_yyyz_xz, g_xz_0_yyyz_yy, g_xz_0_yyyz_yyy, g_xz_0_yyyz_yyz, g_xz_0_yyyz_yz, g_xz_0_yyyz_yzz, g_xz_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyz_xx[k] = -g_xz_0_yyyz_xx[k] * ab_y + g_xz_0_yyyz_xxy[k];

                g_xz_0_yyyyz_xy[k] = -g_xz_0_yyyz_xy[k] * ab_y + g_xz_0_yyyz_xyy[k];

                g_xz_0_yyyyz_xz[k] = -g_xz_0_yyyz_xz[k] * ab_y + g_xz_0_yyyz_xyz[k];

                g_xz_0_yyyyz_yy[k] = -g_xz_0_yyyz_yy[k] * ab_y + g_xz_0_yyyz_yyy[k];

                g_xz_0_yyyyz_yz[k] = -g_xz_0_yyyz_yz[k] * ab_y + g_xz_0_yyyz_yyz[k];

                g_xz_0_yyyyz_zz[k] = -g_xz_0_yyyz_zz[k] * ab_y + g_xz_0_yyyz_yzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzz_xx, g_xz_0_yyyzz_xy, g_xz_0_yyyzz_xz, g_xz_0_yyyzz_yy, g_xz_0_yyyzz_yz, g_xz_0_yyyzz_zz, g_xz_0_yyzz_xx, g_xz_0_yyzz_xxy, g_xz_0_yyzz_xy, g_xz_0_yyzz_xyy, g_xz_0_yyzz_xyz, g_xz_0_yyzz_xz, g_xz_0_yyzz_yy, g_xz_0_yyzz_yyy, g_xz_0_yyzz_yyz, g_xz_0_yyzz_yz, g_xz_0_yyzz_yzz, g_xz_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzz_xx[k] = -g_xz_0_yyzz_xx[k] * ab_y + g_xz_0_yyzz_xxy[k];

                g_xz_0_yyyzz_xy[k] = -g_xz_0_yyzz_xy[k] * ab_y + g_xz_0_yyzz_xyy[k];

                g_xz_0_yyyzz_xz[k] = -g_xz_0_yyzz_xz[k] * ab_y + g_xz_0_yyzz_xyz[k];

                g_xz_0_yyyzz_yy[k] = -g_xz_0_yyzz_yy[k] * ab_y + g_xz_0_yyzz_yyy[k];

                g_xz_0_yyyzz_yz[k] = -g_xz_0_yyzz_yz[k] * ab_y + g_xz_0_yyzz_yyz[k];

                g_xz_0_yyyzz_zz[k] = -g_xz_0_yyzz_zz[k] * ab_y + g_xz_0_yyzz_yzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzz_xx, g_xz_0_yyzzz_xy, g_xz_0_yyzzz_xz, g_xz_0_yyzzz_yy, g_xz_0_yyzzz_yz, g_xz_0_yyzzz_zz, g_xz_0_yzzz_xx, g_xz_0_yzzz_xxy, g_xz_0_yzzz_xy, g_xz_0_yzzz_xyy, g_xz_0_yzzz_xyz, g_xz_0_yzzz_xz, g_xz_0_yzzz_yy, g_xz_0_yzzz_yyy, g_xz_0_yzzz_yyz, g_xz_0_yzzz_yz, g_xz_0_yzzz_yzz, g_xz_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzz_xx[k] = -g_xz_0_yzzz_xx[k] * ab_y + g_xz_0_yzzz_xxy[k];

                g_xz_0_yyzzz_xy[k] = -g_xz_0_yzzz_xy[k] * ab_y + g_xz_0_yzzz_xyy[k];

                g_xz_0_yyzzz_xz[k] = -g_xz_0_yzzz_xz[k] * ab_y + g_xz_0_yzzz_xyz[k];

                g_xz_0_yyzzz_yy[k] = -g_xz_0_yzzz_yy[k] * ab_y + g_xz_0_yzzz_yyy[k];

                g_xz_0_yyzzz_yz[k] = -g_xz_0_yzzz_yz[k] * ab_y + g_xz_0_yzzz_yyz[k];

                g_xz_0_yyzzz_zz[k] = -g_xz_0_yzzz_zz[k] * ab_y + g_xz_0_yzzz_yzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzz_xx, g_xz_0_yzzzz_xy, g_xz_0_yzzzz_xz, g_xz_0_yzzzz_yy, g_xz_0_yzzzz_yz, g_xz_0_yzzzz_zz, g_xz_0_zzzz_xx, g_xz_0_zzzz_xxy, g_xz_0_zzzz_xy, g_xz_0_zzzz_xyy, g_xz_0_zzzz_xyz, g_xz_0_zzzz_xz, g_xz_0_zzzz_yy, g_xz_0_zzzz_yyy, g_xz_0_zzzz_yyz, g_xz_0_zzzz_yz, g_xz_0_zzzz_yzz, g_xz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzz_xx[k] = -g_xz_0_zzzz_xx[k] * ab_y + g_xz_0_zzzz_xxy[k];

                g_xz_0_yzzzz_xy[k] = -g_xz_0_zzzz_xy[k] * ab_y + g_xz_0_zzzz_xyy[k];

                g_xz_0_yzzzz_xz[k] = -g_xz_0_zzzz_xz[k] * ab_y + g_xz_0_zzzz_xyz[k];

                g_xz_0_yzzzz_yy[k] = -g_xz_0_zzzz_yy[k] * ab_y + g_xz_0_zzzz_yyy[k];

                g_xz_0_yzzzz_yz[k] = -g_xz_0_zzzz_yz[k] * ab_y + g_xz_0_zzzz_yyz[k];

                g_xz_0_yzzzz_zz[k] = -g_xz_0_zzzz_zz[k] * ab_y + g_xz_0_zzzz_yzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_xx, g_x_0_zzzz_xy, g_x_0_zzzz_xz, g_x_0_zzzz_yy, g_x_0_zzzz_yz, g_x_0_zzzz_zz, g_xz_0_zzzz_xx, g_xz_0_zzzz_xxz, g_xz_0_zzzz_xy, g_xz_0_zzzz_xyz, g_xz_0_zzzz_xz, g_xz_0_zzzz_xzz, g_xz_0_zzzz_yy, g_xz_0_zzzz_yyz, g_xz_0_zzzz_yz, g_xz_0_zzzz_yzz, g_xz_0_zzzz_zz, g_xz_0_zzzz_zzz, g_xz_0_zzzzz_xx, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_yy, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzz_xx[k] = -g_x_0_zzzz_xx[k] - g_xz_0_zzzz_xx[k] * ab_z + g_xz_0_zzzz_xxz[k];

                g_xz_0_zzzzz_xy[k] = -g_x_0_zzzz_xy[k] - g_xz_0_zzzz_xy[k] * ab_z + g_xz_0_zzzz_xyz[k];

                g_xz_0_zzzzz_xz[k] = -g_x_0_zzzz_xz[k] - g_xz_0_zzzz_xz[k] * ab_z + g_xz_0_zzzz_xzz[k];

                g_xz_0_zzzzz_yy[k] = -g_x_0_zzzz_yy[k] - g_xz_0_zzzz_yy[k] * ab_z + g_xz_0_zzzz_yyz[k];

                g_xz_0_zzzzz_yz[k] = -g_x_0_zzzz_yz[k] - g_xz_0_zzzz_yz[k] * ab_z + g_xz_0_zzzz_yzz[k];

                g_xz_0_zzzzz_zz[k] = -g_x_0_zzzz_zz[k] - g_xz_0_zzzz_zz[k] * ab_z + g_xz_0_zzzz_zzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxx_xx, g_yy_0_xxxx_xxx, g_yy_0_xxxx_xxy, g_yy_0_xxxx_xxz, g_yy_0_xxxx_xy, g_yy_0_xxxx_xyy, g_yy_0_xxxx_xyz, g_yy_0_xxxx_xz, g_yy_0_xxxx_xzz, g_yy_0_xxxx_yy, g_yy_0_xxxx_yz, g_yy_0_xxxx_zz, g_yy_0_xxxxx_xx, g_yy_0_xxxxx_xy, g_yy_0_xxxxx_xz, g_yy_0_xxxxx_yy, g_yy_0_xxxxx_yz, g_yy_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxx_xx[k] = -g_yy_0_xxxx_xx[k] * ab_x + g_yy_0_xxxx_xxx[k];

                g_yy_0_xxxxx_xy[k] = -g_yy_0_xxxx_xy[k] * ab_x + g_yy_0_xxxx_xxy[k];

                g_yy_0_xxxxx_xz[k] = -g_yy_0_xxxx_xz[k] * ab_x + g_yy_0_xxxx_xxz[k];

                g_yy_0_xxxxx_yy[k] = -g_yy_0_xxxx_yy[k] * ab_x + g_yy_0_xxxx_xyy[k];

                g_yy_0_xxxxx_yz[k] = -g_yy_0_xxxx_yz[k] * ab_x + g_yy_0_xxxx_xyz[k];

                g_yy_0_xxxxx_zz[k] = -g_yy_0_xxxx_zz[k] * ab_x + g_yy_0_xxxx_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxy_xx, g_yy_0_xxxxy_xy, g_yy_0_xxxxy_xz, g_yy_0_xxxxy_yy, g_yy_0_xxxxy_yz, g_yy_0_xxxxy_zz, g_yy_0_xxxy_xx, g_yy_0_xxxy_xxx, g_yy_0_xxxy_xxy, g_yy_0_xxxy_xxz, g_yy_0_xxxy_xy, g_yy_0_xxxy_xyy, g_yy_0_xxxy_xyz, g_yy_0_xxxy_xz, g_yy_0_xxxy_xzz, g_yy_0_xxxy_yy, g_yy_0_xxxy_yz, g_yy_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxy_xx[k] = -g_yy_0_xxxy_xx[k] * ab_x + g_yy_0_xxxy_xxx[k];

                g_yy_0_xxxxy_xy[k] = -g_yy_0_xxxy_xy[k] * ab_x + g_yy_0_xxxy_xxy[k];

                g_yy_0_xxxxy_xz[k] = -g_yy_0_xxxy_xz[k] * ab_x + g_yy_0_xxxy_xxz[k];

                g_yy_0_xxxxy_yy[k] = -g_yy_0_xxxy_yy[k] * ab_x + g_yy_0_xxxy_xyy[k];

                g_yy_0_xxxxy_yz[k] = -g_yy_0_xxxy_yz[k] * ab_x + g_yy_0_xxxy_xyz[k];

                g_yy_0_xxxxy_zz[k] = -g_yy_0_xxxy_zz[k] * ab_x + g_yy_0_xxxy_xzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxz_xx, g_yy_0_xxxxz_xy, g_yy_0_xxxxz_xz, g_yy_0_xxxxz_yy, g_yy_0_xxxxz_yz, g_yy_0_xxxxz_zz, g_yy_0_xxxz_xx, g_yy_0_xxxz_xxx, g_yy_0_xxxz_xxy, g_yy_0_xxxz_xxz, g_yy_0_xxxz_xy, g_yy_0_xxxz_xyy, g_yy_0_xxxz_xyz, g_yy_0_xxxz_xz, g_yy_0_xxxz_xzz, g_yy_0_xxxz_yy, g_yy_0_xxxz_yz, g_yy_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxz_xx[k] = -g_yy_0_xxxz_xx[k] * ab_x + g_yy_0_xxxz_xxx[k];

                g_yy_0_xxxxz_xy[k] = -g_yy_0_xxxz_xy[k] * ab_x + g_yy_0_xxxz_xxy[k];

                g_yy_0_xxxxz_xz[k] = -g_yy_0_xxxz_xz[k] * ab_x + g_yy_0_xxxz_xxz[k];

                g_yy_0_xxxxz_yy[k] = -g_yy_0_xxxz_yy[k] * ab_x + g_yy_0_xxxz_xyy[k];

                g_yy_0_xxxxz_yz[k] = -g_yy_0_xxxz_yz[k] * ab_x + g_yy_0_xxxz_xyz[k];

                g_yy_0_xxxxz_zz[k] = -g_yy_0_xxxz_zz[k] * ab_x + g_yy_0_xxxz_xzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyy_xx, g_yy_0_xxxyy_xy, g_yy_0_xxxyy_xz, g_yy_0_xxxyy_yy, g_yy_0_xxxyy_yz, g_yy_0_xxxyy_zz, g_yy_0_xxyy_xx, g_yy_0_xxyy_xxx, g_yy_0_xxyy_xxy, g_yy_0_xxyy_xxz, g_yy_0_xxyy_xy, g_yy_0_xxyy_xyy, g_yy_0_xxyy_xyz, g_yy_0_xxyy_xz, g_yy_0_xxyy_xzz, g_yy_0_xxyy_yy, g_yy_0_xxyy_yz, g_yy_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyy_xx[k] = -g_yy_0_xxyy_xx[k] * ab_x + g_yy_0_xxyy_xxx[k];

                g_yy_0_xxxyy_xy[k] = -g_yy_0_xxyy_xy[k] * ab_x + g_yy_0_xxyy_xxy[k];

                g_yy_0_xxxyy_xz[k] = -g_yy_0_xxyy_xz[k] * ab_x + g_yy_0_xxyy_xxz[k];

                g_yy_0_xxxyy_yy[k] = -g_yy_0_xxyy_yy[k] * ab_x + g_yy_0_xxyy_xyy[k];

                g_yy_0_xxxyy_yz[k] = -g_yy_0_xxyy_yz[k] * ab_x + g_yy_0_xxyy_xyz[k];

                g_yy_0_xxxyy_zz[k] = -g_yy_0_xxyy_zz[k] * ab_x + g_yy_0_xxyy_xzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyz_xx, g_yy_0_xxxyz_xy, g_yy_0_xxxyz_xz, g_yy_0_xxxyz_yy, g_yy_0_xxxyz_yz, g_yy_0_xxxyz_zz, g_yy_0_xxyz_xx, g_yy_0_xxyz_xxx, g_yy_0_xxyz_xxy, g_yy_0_xxyz_xxz, g_yy_0_xxyz_xy, g_yy_0_xxyz_xyy, g_yy_0_xxyz_xyz, g_yy_0_xxyz_xz, g_yy_0_xxyz_xzz, g_yy_0_xxyz_yy, g_yy_0_xxyz_yz, g_yy_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyz_xx[k] = -g_yy_0_xxyz_xx[k] * ab_x + g_yy_0_xxyz_xxx[k];

                g_yy_0_xxxyz_xy[k] = -g_yy_0_xxyz_xy[k] * ab_x + g_yy_0_xxyz_xxy[k];

                g_yy_0_xxxyz_xz[k] = -g_yy_0_xxyz_xz[k] * ab_x + g_yy_0_xxyz_xxz[k];

                g_yy_0_xxxyz_yy[k] = -g_yy_0_xxyz_yy[k] * ab_x + g_yy_0_xxyz_xyy[k];

                g_yy_0_xxxyz_yz[k] = -g_yy_0_xxyz_yz[k] * ab_x + g_yy_0_xxyz_xyz[k];

                g_yy_0_xxxyz_zz[k] = -g_yy_0_xxyz_zz[k] * ab_x + g_yy_0_xxyz_xzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzz_xx, g_yy_0_xxxzz_xy, g_yy_0_xxxzz_xz, g_yy_0_xxxzz_yy, g_yy_0_xxxzz_yz, g_yy_0_xxxzz_zz, g_yy_0_xxzz_xx, g_yy_0_xxzz_xxx, g_yy_0_xxzz_xxy, g_yy_0_xxzz_xxz, g_yy_0_xxzz_xy, g_yy_0_xxzz_xyy, g_yy_0_xxzz_xyz, g_yy_0_xxzz_xz, g_yy_0_xxzz_xzz, g_yy_0_xxzz_yy, g_yy_0_xxzz_yz, g_yy_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzz_xx[k] = -g_yy_0_xxzz_xx[k] * ab_x + g_yy_0_xxzz_xxx[k];

                g_yy_0_xxxzz_xy[k] = -g_yy_0_xxzz_xy[k] * ab_x + g_yy_0_xxzz_xxy[k];

                g_yy_0_xxxzz_xz[k] = -g_yy_0_xxzz_xz[k] * ab_x + g_yy_0_xxzz_xxz[k];

                g_yy_0_xxxzz_yy[k] = -g_yy_0_xxzz_yy[k] * ab_x + g_yy_0_xxzz_xyy[k];

                g_yy_0_xxxzz_yz[k] = -g_yy_0_xxzz_yz[k] * ab_x + g_yy_0_xxzz_xyz[k];

                g_yy_0_xxxzz_zz[k] = -g_yy_0_xxzz_zz[k] * ab_x + g_yy_0_xxzz_xzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyy_xx, g_yy_0_xxyyy_xy, g_yy_0_xxyyy_xz, g_yy_0_xxyyy_yy, g_yy_0_xxyyy_yz, g_yy_0_xxyyy_zz, g_yy_0_xyyy_xx, g_yy_0_xyyy_xxx, g_yy_0_xyyy_xxy, g_yy_0_xyyy_xxz, g_yy_0_xyyy_xy, g_yy_0_xyyy_xyy, g_yy_0_xyyy_xyz, g_yy_0_xyyy_xz, g_yy_0_xyyy_xzz, g_yy_0_xyyy_yy, g_yy_0_xyyy_yz, g_yy_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyy_xx[k] = -g_yy_0_xyyy_xx[k] * ab_x + g_yy_0_xyyy_xxx[k];

                g_yy_0_xxyyy_xy[k] = -g_yy_0_xyyy_xy[k] * ab_x + g_yy_0_xyyy_xxy[k];

                g_yy_0_xxyyy_xz[k] = -g_yy_0_xyyy_xz[k] * ab_x + g_yy_0_xyyy_xxz[k];

                g_yy_0_xxyyy_yy[k] = -g_yy_0_xyyy_yy[k] * ab_x + g_yy_0_xyyy_xyy[k];

                g_yy_0_xxyyy_yz[k] = -g_yy_0_xyyy_yz[k] * ab_x + g_yy_0_xyyy_xyz[k];

                g_yy_0_xxyyy_zz[k] = -g_yy_0_xyyy_zz[k] * ab_x + g_yy_0_xyyy_xzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyz_xx, g_yy_0_xxyyz_xy, g_yy_0_xxyyz_xz, g_yy_0_xxyyz_yy, g_yy_0_xxyyz_yz, g_yy_0_xxyyz_zz, g_yy_0_xyyz_xx, g_yy_0_xyyz_xxx, g_yy_0_xyyz_xxy, g_yy_0_xyyz_xxz, g_yy_0_xyyz_xy, g_yy_0_xyyz_xyy, g_yy_0_xyyz_xyz, g_yy_0_xyyz_xz, g_yy_0_xyyz_xzz, g_yy_0_xyyz_yy, g_yy_0_xyyz_yz, g_yy_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyz_xx[k] = -g_yy_0_xyyz_xx[k] * ab_x + g_yy_0_xyyz_xxx[k];

                g_yy_0_xxyyz_xy[k] = -g_yy_0_xyyz_xy[k] * ab_x + g_yy_0_xyyz_xxy[k];

                g_yy_0_xxyyz_xz[k] = -g_yy_0_xyyz_xz[k] * ab_x + g_yy_0_xyyz_xxz[k];

                g_yy_0_xxyyz_yy[k] = -g_yy_0_xyyz_yy[k] * ab_x + g_yy_0_xyyz_xyy[k];

                g_yy_0_xxyyz_yz[k] = -g_yy_0_xyyz_yz[k] * ab_x + g_yy_0_xyyz_xyz[k];

                g_yy_0_xxyyz_zz[k] = -g_yy_0_xyyz_zz[k] * ab_x + g_yy_0_xyyz_xzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzz_xx, g_yy_0_xxyzz_xy, g_yy_0_xxyzz_xz, g_yy_0_xxyzz_yy, g_yy_0_xxyzz_yz, g_yy_0_xxyzz_zz, g_yy_0_xyzz_xx, g_yy_0_xyzz_xxx, g_yy_0_xyzz_xxy, g_yy_0_xyzz_xxz, g_yy_0_xyzz_xy, g_yy_0_xyzz_xyy, g_yy_0_xyzz_xyz, g_yy_0_xyzz_xz, g_yy_0_xyzz_xzz, g_yy_0_xyzz_yy, g_yy_0_xyzz_yz, g_yy_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzz_xx[k] = -g_yy_0_xyzz_xx[k] * ab_x + g_yy_0_xyzz_xxx[k];

                g_yy_0_xxyzz_xy[k] = -g_yy_0_xyzz_xy[k] * ab_x + g_yy_0_xyzz_xxy[k];

                g_yy_0_xxyzz_xz[k] = -g_yy_0_xyzz_xz[k] * ab_x + g_yy_0_xyzz_xxz[k];

                g_yy_0_xxyzz_yy[k] = -g_yy_0_xyzz_yy[k] * ab_x + g_yy_0_xyzz_xyy[k];

                g_yy_0_xxyzz_yz[k] = -g_yy_0_xyzz_yz[k] * ab_x + g_yy_0_xyzz_xyz[k];

                g_yy_0_xxyzz_zz[k] = -g_yy_0_xyzz_zz[k] * ab_x + g_yy_0_xyzz_xzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzz_xx, g_yy_0_xxzzz_xy, g_yy_0_xxzzz_xz, g_yy_0_xxzzz_yy, g_yy_0_xxzzz_yz, g_yy_0_xxzzz_zz, g_yy_0_xzzz_xx, g_yy_0_xzzz_xxx, g_yy_0_xzzz_xxy, g_yy_0_xzzz_xxz, g_yy_0_xzzz_xy, g_yy_0_xzzz_xyy, g_yy_0_xzzz_xyz, g_yy_0_xzzz_xz, g_yy_0_xzzz_xzz, g_yy_0_xzzz_yy, g_yy_0_xzzz_yz, g_yy_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzz_xx[k] = -g_yy_0_xzzz_xx[k] * ab_x + g_yy_0_xzzz_xxx[k];

                g_yy_0_xxzzz_xy[k] = -g_yy_0_xzzz_xy[k] * ab_x + g_yy_0_xzzz_xxy[k];

                g_yy_0_xxzzz_xz[k] = -g_yy_0_xzzz_xz[k] * ab_x + g_yy_0_xzzz_xxz[k];

                g_yy_0_xxzzz_yy[k] = -g_yy_0_xzzz_yy[k] * ab_x + g_yy_0_xzzz_xyy[k];

                g_yy_0_xxzzz_yz[k] = -g_yy_0_xzzz_yz[k] * ab_x + g_yy_0_xzzz_xyz[k];

                g_yy_0_xxzzz_zz[k] = -g_yy_0_xzzz_zz[k] * ab_x + g_yy_0_xzzz_xzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyy_xx, g_yy_0_xyyyy_xy, g_yy_0_xyyyy_xz, g_yy_0_xyyyy_yy, g_yy_0_xyyyy_yz, g_yy_0_xyyyy_zz, g_yy_0_yyyy_xx, g_yy_0_yyyy_xxx, g_yy_0_yyyy_xxy, g_yy_0_yyyy_xxz, g_yy_0_yyyy_xy, g_yy_0_yyyy_xyy, g_yy_0_yyyy_xyz, g_yy_0_yyyy_xz, g_yy_0_yyyy_xzz, g_yy_0_yyyy_yy, g_yy_0_yyyy_yz, g_yy_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyy_xx[k] = -g_yy_0_yyyy_xx[k] * ab_x + g_yy_0_yyyy_xxx[k];

                g_yy_0_xyyyy_xy[k] = -g_yy_0_yyyy_xy[k] * ab_x + g_yy_0_yyyy_xxy[k];

                g_yy_0_xyyyy_xz[k] = -g_yy_0_yyyy_xz[k] * ab_x + g_yy_0_yyyy_xxz[k];

                g_yy_0_xyyyy_yy[k] = -g_yy_0_yyyy_yy[k] * ab_x + g_yy_0_yyyy_xyy[k];

                g_yy_0_xyyyy_yz[k] = -g_yy_0_yyyy_yz[k] * ab_x + g_yy_0_yyyy_xyz[k];

                g_yy_0_xyyyy_zz[k] = -g_yy_0_yyyy_zz[k] * ab_x + g_yy_0_yyyy_xzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyz_xx, g_yy_0_xyyyz_xy, g_yy_0_xyyyz_xz, g_yy_0_xyyyz_yy, g_yy_0_xyyyz_yz, g_yy_0_xyyyz_zz, g_yy_0_yyyz_xx, g_yy_0_yyyz_xxx, g_yy_0_yyyz_xxy, g_yy_0_yyyz_xxz, g_yy_0_yyyz_xy, g_yy_0_yyyz_xyy, g_yy_0_yyyz_xyz, g_yy_0_yyyz_xz, g_yy_0_yyyz_xzz, g_yy_0_yyyz_yy, g_yy_0_yyyz_yz, g_yy_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyz_xx[k] = -g_yy_0_yyyz_xx[k] * ab_x + g_yy_0_yyyz_xxx[k];

                g_yy_0_xyyyz_xy[k] = -g_yy_0_yyyz_xy[k] * ab_x + g_yy_0_yyyz_xxy[k];

                g_yy_0_xyyyz_xz[k] = -g_yy_0_yyyz_xz[k] * ab_x + g_yy_0_yyyz_xxz[k];

                g_yy_0_xyyyz_yy[k] = -g_yy_0_yyyz_yy[k] * ab_x + g_yy_0_yyyz_xyy[k];

                g_yy_0_xyyyz_yz[k] = -g_yy_0_yyyz_yz[k] * ab_x + g_yy_0_yyyz_xyz[k];

                g_yy_0_xyyyz_zz[k] = -g_yy_0_yyyz_zz[k] * ab_x + g_yy_0_yyyz_xzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzz_xx, g_yy_0_xyyzz_xy, g_yy_0_xyyzz_xz, g_yy_0_xyyzz_yy, g_yy_0_xyyzz_yz, g_yy_0_xyyzz_zz, g_yy_0_yyzz_xx, g_yy_0_yyzz_xxx, g_yy_0_yyzz_xxy, g_yy_0_yyzz_xxz, g_yy_0_yyzz_xy, g_yy_0_yyzz_xyy, g_yy_0_yyzz_xyz, g_yy_0_yyzz_xz, g_yy_0_yyzz_xzz, g_yy_0_yyzz_yy, g_yy_0_yyzz_yz, g_yy_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzz_xx[k] = -g_yy_0_yyzz_xx[k] * ab_x + g_yy_0_yyzz_xxx[k];

                g_yy_0_xyyzz_xy[k] = -g_yy_0_yyzz_xy[k] * ab_x + g_yy_0_yyzz_xxy[k];

                g_yy_0_xyyzz_xz[k] = -g_yy_0_yyzz_xz[k] * ab_x + g_yy_0_yyzz_xxz[k];

                g_yy_0_xyyzz_yy[k] = -g_yy_0_yyzz_yy[k] * ab_x + g_yy_0_yyzz_xyy[k];

                g_yy_0_xyyzz_yz[k] = -g_yy_0_yyzz_yz[k] * ab_x + g_yy_0_yyzz_xyz[k];

                g_yy_0_xyyzz_zz[k] = -g_yy_0_yyzz_zz[k] * ab_x + g_yy_0_yyzz_xzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzz_xx, g_yy_0_xyzzz_xy, g_yy_0_xyzzz_xz, g_yy_0_xyzzz_yy, g_yy_0_xyzzz_yz, g_yy_0_xyzzz_zz, g_yy_0_yzzz_xx, g_yy_0_yzzz_xxx, g_yy_0_yzzz_xxy, g_yy_0_yzzz_xxz, g_yy_0_yzzz_xy, g_yy_0_yzzz_xyy, g_yy_0_yzzz_xyz, g_yy_0_yzzz_xz, g_yy_0_yzzz_xzz, g_yy_0_yzzz_yy, g_yy_0_yzzz_yz, g_yy_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzz_xx[k] = -g_yy_0_yzzz_xx[k] * ab_x + g_yy_0_yzzz_xxx[k];

                g_yy_0_xyzzz_xy[k] = -g_yy_0_yzzz_xy[k] * ab_x + g_yy_0_yzzz_xxy[k];

                g_yy_0_xyzzz_xz[k] = -g_yy_0_yzzz_xz[k] * ab_x + g_yy_0_yzzz_xxz[k];

                g_yy_0_xyzzz_yy[k] = -g_yy_0_yzzz_yy[k] * ab_x + g_yy_0_yzzz_xyy[k];

                g_yy_0_xyzzz_yz[k] = -g_yy_0_yzzz_yz[k] * ab_x + g_yy_0_yzzz_xyz[k];

                g_yy_0_xyzzz_zz[k] = -g_yy_0_yzzz_zz[k] * ab_x + g_yy_0_yzzz_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzz_xx, g_yy_0_xzzzz_xy, g_yy_0_xzzzz_xz, g_yy_0_xzzzz_yy, g_yy_0_xzzzz_yz, g_yy_0_xzzzz_zz, g_yy_0_zzzz_xx, g_yy_0_zzzz_xxx, g_yy_0_zzzz_xxy, g_yy_0_zzzz_xxz, g_yy_0_zzzz_xy, g_yy_0_zzzz_xyy, g_yy_0_zzzz_xyz, g_yy_0_zzzz_xz, g_yy_0_zzzz_xzz, g_yy_0_zzzz_yy, g_yy_0_zzzz_yz, g_yy_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzz_xx[k] = -g_yy_0_zzzz_xx[k] * ab_x + g_yy_0_zzzz_xxx[k];

                g_yy_0_xzzzz_xy[k] = -g_yy_0_zzzz_xy[k] * ab_x + g_yy_0_zzzz_xxy[k];

                g_yy_0_xzzzz_xz[k] = -g_yy_0_zzzz_xz[k] * ab_x + g_yy_0_zzzz_xxz[k];

                g_yy_0_xzzzz_yy[k] = -g_yy_0_zzzz_yy[k] * ab_x + g_yy_0_zzzz_xyy[k];

                g_yy_0_xzzzz_yz[k] = -g_yy_0_zzzz_yz[k] * ab_x + g_yy_0_zzzz_xyz[k];

                g_yy_0_xzzzz_zz[k] = -g_yy_0_zzzz_zz[k] * ab_x + g_yy_0_zzzz_xzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz, g_yy_0_yyyy_xx, g_yy_0_yyyy_xxy, g_yy_0_yyyy_xy, g_yy_0_yyyy_xyy, g_yy_0_yyyy_xyz, g_yy_0_yyyy_xz, g_yy_0_yyyy_yy, g_yy_0_yyyy_yyy, g_yy_0_yyyy_yyz, g_yy_0_yyyy_yz, g_yy_0_yyyy_yzz, g_yy_0_yyyy_zz, g_yy_0_yyyyy_xx, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_yy, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyy_xx[k] = -2.0 * g_y_0_yyyy_xx[k] - g_yy_0_yyyy_xx[k] * ab_y + g_yy_0_yyyy_xxy[k];

                g_yy_0_yyyyy_xy[k] = -2.0 * g_y_0_yyyy_xy[k] - g_yy_0_yyyy_xy[k] * ab_y + g_yy_0_yyyy_xyy[k];

                g_yy_0_yyyyy_xz[k] = -2.0 * g_y_0_yyyy_xz[k] - g_yy_0_yyyy_xz[k] * ab_y + g_yy_0_yyyy_xyz[k];

                g_yy_0_yyyyy_yy[k] = -2.0 * g_y_0_yyyy_yy[k] - g_yy_0_yyyy_yy[k] * ab_y + g_yy_0_yyyy_yyy[k];

                g_yy_0_yyyyy_yz[k] = -2.0 * g_y_0_yyyy_yz[k] - g_yy_0_yyyy_yz[k] * ab_y + g_yy_0_yyyy_yyz[k];

                g_yy_0_yyyyy_zz[k] = -2.0 * g_y_0_yyyy_zz[k] - g_yy_0_yyyy_zz[k] * ab_y + g_yy_0_yyyy_yzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyy_xx, g_yy_0_yyyy_xxz, g_yy_0_yyyy_xy, g_yy_0_yyyy_xyz, g_yy_0_yyyy_xz, g_yy_0_yyyy_xzz, g_yy_0_yyyy_yy, g_yy_0_yyyy_yyz, g_yy_0_yyyy_yz, g_yy_0_yyyy_yzz, g_yy_0_yyyy_zz, g_yy_0_yyyy_zzz, g_yy_0_yyyyz_xx, g_yy_0_yyyyz_xy, g_yy_0_yyyyz_xz, g_yy_0_yyyyz_yy, g_yy_0_yyyyz_yz, g_yy_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyz_xx[k] = -g_yy_0_yyyy_xx[k] * ab_z + g_yy_0_yyyy_xxz[k];

                g_yy_0_yyyyz_xy[k] = -g_yy_0_yyyy_xy[k] * ab_z + g_yy_0_yyyy_xyz[k];

                g_yy_0_yyyyz_xz[k] = -g_yy_0_yyyy_xz[k] * ab_z + g_yy_0_yyyy_xzz[k];

                g_yy_0_yyyyz_yy[k] = -g_yy_0_yyyy_yy[k] * ab_z + g_yy_0_yyyy_yyz[k];

                g_yy_0_yyyyz_yz[k] = -g_yy_0_yyyy_yz[k] * ab_z + g_yy_0_yyyy_yzz[k];

                g_yy_0_yyyyz_zz[k] = -g_yy_0_yyyy_zz[k] * ab_z + g_yy_0_yyyy_zzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyz_xx, g_yy_0_yyyz_xxz, g_yy_0_yyyz_xy, g_yy_0_yyyz_xyz, g_yy_0_yyyz_xz, g_yy_0_yyyz_xzz, g_yy_0_yyyz_yy, g_yy_0_yyyz_yyz, g_yy_0_yyyz_yz, g_yy_0_yyyz_yzz, g_yy_0_yyyz_zz, g_yy_0_yyyz_zzz, g_yy_0_yyyzz_xx, g_yy_0_yyyzz_xy, g_yy_0_yyyzz_xz, g_yy_0_yyyzz_yy, g_yy_0_yyyzz_yz, g_yy_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzz_xx[k] = -g_yy_0_yyyz_xx[k] * ab_z + g_yy_0_yyyz_xxz[k];

                g_yy_0_yyyzz_xy[k] = -g_yy_0_yyyz_xy[k] * ab_z + g_yy_0_yyyz_xyz[k];

                g_yy_0_yyyzz_xz[k] = -g_yy_0_yyyz_xz[k] * ab_z + g_yy_0_yyyz_xzz[k];

                g_yy_0_yyyzz_yy[k] = -g_yy_0_yyyz_yy[k] * ab_z + g_yy_0_yyyz_yyz[k];

                g_yy_0_yyyzz_yz[k] = -g_yy_0_yyyz_yz[k] * ab_z + g_yy_0_yyyz_yzz[k];

                g_yy_0_yyyzz_zz[k] = -g_yy_0_yyyz_zz[k] * ab_z + g_yy_0_yyyz_zzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzz_xx, g_yy_0_yyzz_xxz, g_yy_0_yyzz_xy, g_yy_0_yyzz_xyz, g_yy_0_yyzz_xz, g_yy_0_yyzz_xzz, g_yy_0_yyzz_yy, g_yy_0_yyzz_yyz, g_yy_0_yyzz_yz, g_yy_0_yyzz_yzz, g_yy_0_yyzz_zz, g_yy_0_yyzz_zzz, g_yy_0_yyzzz_xx, g_yy_0_yyzzz_xy, g_yy_0_yyzzz_xz, g_yy_0_yyzzz_yy, g_yy_0_yyzzz_yz, g_yy_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzz_xx[k] = -g_yy_0_yyzz_xx[k] * ab_z + g_yy_0_yyzz_xxz[k];

                g_yy_0_yyzzz_xy[k] = -g_yy_0_yyzz_xy[k] * ab_z + g_yy_0_yyzz_xyz[k];

                g_yy_0_yyzzz_xz[k] = -g_yy_0_yyzz_xz[k] * ab_z + g_yy_0_yyzz_xzz[k];

                g_yy_0_yyzzz_yy[k] = -g_yy_0_yyzz_yy[k] * ab_z + g_yy_0_yyzz_yyz[k];

                g_yy_0_yyzzz_yz[k] = -g_yy_0_yyzz_yz[k] * ab_z + g_yy_0_yyzz_yzz[k];

                g_yy_0_yyzzz_zz[k] = -g_yy_0_yyzz_zz[k] * ab_z + g_yy_0_yyzz_zzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzz_xx, g_yy_0_yzzz_xxz, g_yy_0_yzzz_xy, g_yy_0_yzzz_xyz, g_yy_0_yzzz_xz, g_yy_0_yzzz_xzz, g_yy_0_yzzz_yy, g_yy_0_yzzz_yyz, g_yy_0_yzzz_yz, g_yy_0_yzzz_yzz, g_yy_0_yzzz_zz, g_yy_0_yzzz_zzz, g_yy_0_yzzzz_xx, g_yy_0_yzzzz_xy, g_yy_0_yzzzz_xz, g_yy_0_yzzzz_yy, g_yy_0_yzzzz_yz, g_yy_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzz_xx[k] = -g_yy_0_yzzz_xx[k] * ab_z + g_yy_0_yzzz_xxz[k];

                g_yy_0_yzzzz_xy[k] = -g_yy_0_yzzz_xy[k] * ab_z + g_yy_0_yzzz_xyz[k];

                g_yy_0_yzzzz_xz[k] = -g_yy_0_yzzz_xz[k] * ab_z + g_yy_0_yzzz_xzz[k];

                g_yy_0_yzzzz_yy[k] = -g_yy_0_yzzz_yy[k] * ab_z + g_yy_0_yzzz_yyz[k];

                g_yy_0_yzzzz_yz[k] = -g_yy_0_yzzz_yz[k] * ab_z + g_yy_0_yzzz_yzz[k];

                g_yy_0_yzzzz_zz[k] = -g_yy_0_yzzz_zz[k] * ab_z + g_yy_0_yzzz_zzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzz_xx, g_yy_0_zzzz_xxz, g_yy_0_zzzz_xy, g_yy_0_zzzz_xyz, g_yy_0_zzzz_xz, g_yy_0_zzzz_xzz, g_yy_0_zzzz_yy, g_yy_0_zzzz_yyz, g_yy_0_zzzz_yz, g_yy_0_zzzz_yzz, g_yy_0_zzzz_zz, g_yy_0_zzzz_zzz, g_yy_0_zzzzz_xx, g_yy_0_zzzzz_xy, g_yy_0_zzzzz_xz, g_yy_0_zzzzz_yy, g_yy_0_zzzzz_yz, g_yy_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzz_xx[k] = -g_yy_0_zzzz_xx[k] * ab_z + g_yy_0_zzzz_xxz[k];

                g_yy_0_zzzzz_xy[k] = -g_yy_0_zzzz_xy[k] * ab_z + g_yy_0_zzzz_xyz[k];

                g_yy_0_zzzzz_xz[k] = -g_yy_0_zzzz_xz[k] * ab_z + g_yy_0_zzzz_xzz[k];

                g_yy_0_zzzzz_yy[k] = -g_yy_0_zzzz_yy[k] * ab_z + g_yy_0_zzzz_yyz[k];

                g_yy_0_zzzzz_yz[k] = -g_yy_0_zzzz_yz[k] * ab_z + g_yy_0_zzzz_yzz[k];

                g_yy_0_zzzzz_zz[k] = -g_yy_0_zzzz_zz[k] * ab_z + g_yy_0_zzzz_zzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxx_xx, g_yz_0_xxxx_xxx, g_yz_0_xxxx_xxy, g_yz_0_xxxx_xxz, g_yz_0_xxxx_xy, g_yz_0_xxxx_xyy, g_yz_0_xxxx_xyz, g_yz_0_xxxx_xz, g_yz_0_xxxx_xzz, g_yz_0_xxxx_yy, g_yz_0_xxxx_yz, g_yz_0_xxxx_zz, g_yz_0_xxxxx_xx, g_yz_0_xxxxx_xy, g_yz_0_xxxxx_xz, g_yz_0_xxxxx_yy, g_yz_0_xxxxx_yz, g_yz_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxx_xx[k] = -g_yz_0_xxxx_xx[k] * ab_x + g_yz_0_xxxx_xxx[k];

                g_yz_0_xxxxx_xy[k] = -g_yz_0_xxxx_xy[k] * ab_x + g_yz_0_xxxx_xxy[k];

                g_yz_0_xxxxx_xz[k] = -g_yz_0_xxxx_xz[k] * ab_x + g_yz_0_xxxx_xxz[k];

                g_yz_0_xxxxx_yy[k] = -g_yz_0_xxxx_yy[k] * ab_x + g_yz_0_xxxx_xyy[k];

                g_yz_0_xxxxx_yz[k] = -g_yz_0_xxxx_yz[k] * ab_x + g_yz_0_xxxx_xyz[k];

                g_yz_0_xxxxx_zz[k] = -g_yz_0_xxxx_zz[k] * ab_x + g_yz_0_xxxx_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxy_xx, g_yz_0_xxxxy_xy, g_yz_0_xxxxy_xz, g_yz_0_xxxxy_yy, g_yz_0_xxxxy_yz, g_yz_0_xxxxy_zz, g_yz_0_xxxy_xx, g_yz_0_xxxy_xxx, g_yz_0_xxxy_xxy, g_yz_0_xxxy_xxz, g_yz_0_xxxy_xy, g_yz_0_xxxy_xyy, g_yz_0_xxxy_xyz, g_yz_0_xxxy_xz, g_yz_0_xxxy_xzz, g_yz_0_xxxy_yy, g_yz_0_xxxy_yz, g_yz_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxy_xx[k] = -g_yz_0_xxxy_xx[k] * ab_x + g_yz_0_xxxy_xxx[k];

                g_yz_0_xxxxy_xy[k] = -g_yz_0_xxxy_xy[k] * ab_x + g_yz_0_xxxy_xxy[k];

                g_yz_0_xxxxy_xz[k] = -g_yz_0_xxxy_xz[k] * ab_x + g_yz_0_xxxy_xxz[k];

                g_yz_0_xxxxy_yy[k] = -g_yz_0_xxxy_yy[k] * ab_x + g_yz_0_xxxy_xyy[k];

                g_yz_0_xxxxy_yz[k] = -g_yz_0_xxxy_yz[k] * ab_x + g_yz_0_xxxy_xyz[k];

                g_yz_0_xxxxy_zz[k] = -g_yz_0_xxxy_zz[k] * ab_x + g_yz_0_xxxy_xzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxz_xx, g_yz_0_xxxxz_xy, g_yz_0_xxxxz_xz, g_yz_0_xxxxz_yy, g_yz_0_xxxxz_yz, g_yz_0_xxxxz_zz, g_yz_0_xxxz_xx, g_yz_0_xxxz_xxx, g_yz_0_xxxz_xxy, g_yz_0_xxxz_xxz, g_yz_0_xxxz_xy, g_yz_0_xxxz_xyy, g_yz_0_xxxz_xyz, g_yz_0_xxxz_xz, g_yz_0_xxxz_xzz, g_yz_0_xxxz_yy, g_yz_0_xxxz_yz, g_yz_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxz_xx[k] = -g_yz_0_xxxz_xx[k] * ab_x + g_yz_0_xxxz_xxx[k];

                g_yz_0_xxxxz_xy[k] = -g_yz_0_xxxz_xy[k] * ab_x + g_yz_0_xxxz_xxy[k];

                g_yz_0_xxxxz_xz[k] = -g_yz_0_xxxz_xz[k] * ab_x + g_yz_0_xxxz_xxz[k];

                g_yz_0_xxxxz_yy[k] = -g_yz_0_xxxz_yy[k] * ab_x + g_yz_0_xxxz_xyy[k];

                g_yz_0_xxxxz_yz[k] = -g_yz_0_xxxz_yz[k] * ab_x + g_yz_0_xxxz_xyz[k];

                g_yz_0_xxxxz_zz[k] = -g_yz_0_xxxz_zz[k] * ab_x + g_yz_0_xxxz_xzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 524 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyy_xx, g_yz_0_xxxyy_xy, g_yz_0_xxxyy_xz, g_yz_0_xxxyy_yy, g_yz_0_xxxyy_yz, g_yz_0_xxxyy_zz, g_yz_0_xxyy_xx, g_yz_0_xxyy_xxx, g_yz_0_xxyy_xxy, g_yz_0_xxyy_xxz, g_yz_0_xxyy_xy, g_yz_0_xxyy_xyy, g_yz_0_xxyy_xyz, g_yz_0_xxyy_xz, g_yz_0_xxyy_xzz, g_yz_0_xxyy_yy, g_yz_0_xxyy_yz, g_yz_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyy_xx[k] = -g_yz_0_xxyy_xx[k] * ab_x + g_yz_0_xxyy_xxx[k];

                g_yz_0_xxxyy_xy[k] = -g_yz_0_xxyy_xy[k] * ab_x + g_yz_0_xxyy_xxy[k];

                g_yz_0_xxxyy_xz[k] = -g_yz_0_xxyy_xz[k] * ab_x + g_yz_0_xxyy_xxz[k];

                g_yz_0_xxxyy_yy[k] = -g_yz_0_xxyy_yy[k] * ab_x + g_yz_0_xxyy_xyy[k];

                g_yz_0_xxxyy_yz[k] = -g_yz_0_xxyy_yz[k] * ab_x + g_yz_0_xxyy_xyz[k];

                g_yz_0_xxxyy_zz[k] = -g_yz_0_xxyy_zz[k] * ab_x + g_yz_0_xxyy_xzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyz_xx, g_yz_0_xxxyz_xy, g_yz_0_xxxyz_xz, g_yz_0_xxxyz_yy, g_yz_0_xxxyz_yz, g_yz_0_xxxyz_zz, g_yz_0_xxyz_xx, g_yz_0_xxyz_xxx, g_yz_0_xxyz_xxy, g_yz_0_xxyz_xxz, g_yz_0_xxyz_xy, g_yz_0_xxyz_xyy, g_yz_0_xxyz_xyz, g_yz_0_xxyz_xz, g_yz_0_xxyz_xzz, g_yz_0_xxyz_yy, g_yz_0_xxyz_yz, g_yz_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyz_xx[k] = -g_yz_0_xxyz_xx[k] * ab_x + g_yz_0_xxyz_xxx[k];

                g_yz_0_xxxyz_xy[k] = -g_yz_0_xxyz_xy[k] * ab_x + g_yz_0_xxyz_xxy[k];

                g_yz_0_xxxyz_xz[k] = -g_yz_0_xxyz_xz[k] * ab_x + g_yz_0_xxyz_xxz[k];

                g_yz_0_xxxyz_yy[k] = -g_yz_0_xxyz_yy[k] * ab_x + g_yz_0_xxyz_xyy[k];

                g_yz_0_xxxyz_yz[k] = -g_yz_0_xxyz_yz[k] * ab_x + g_yz_0_xxyz_xyz[k];

                g_yz_0_xxxyz_zz[k] = -g_yz_0_xxyz_zz[k] * ab_x + g_yz_0_xxyz_xzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzz_xx, g_yz_0_xxxzz_xy, g_yz_0_xxxzz_xz, g_yz_0_xxxzz_yy, g_yz_0_xxxzz_yz, g_yz_0_xxxzz_zz, g_yz_0_xxzz_xx, g_yz_0_xxzz_xxx, g_yz_0_xxzz_xxy, g_yz_0_xxzz_xxz, g_yz_0_xxzz_xy, g_yz_0_xxzz_xyy, g_yz_0_xxzz_xyz, g_yz_0_xxzz_xz, g_yz_0_xxzz_xzz, g_yz_0_xxzz_yy, g_yz_0_xxzz_yz, g_yz_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzz_xx[k] = -g_yz_0_xxzz_xx[k] * ab_x + g_yz_0_xxzz_xxx[k];

                g_yz_0_xxxzz_xy[k] = -g_yz_0_xxzz_xy[k] * ab_x + g_yz_0_xxzz_xxy[k];

                g_yz_0_xxxzz_xz[k] = -g_yz_0_xxzz_xz[k] * ab_x + g_yz_0_xxzz_xxz[k];

                g_yz_0_xxxzz_yy[k] = -g_yz_0_xxzz_yy[k] * ab_x + g_yz_0_xxzz_xyy[k];

                g_yz_0_xxxzz_yz[k] = -g_yz_0_xxzz_yz[k] * ab_x + g_yz_0_xxzz_xyz[k];

                g_yz_0_xxxzz_zz[k] = -g_yz_0_xxzz_zz[k] * ab_x + g_yz_0_xxzz_xzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 540 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 541 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 542 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 543 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 544 * ccomps * dcomps);

            auto g_yz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyy_xx, g_yz_0_xxyyy_xy, g_yz_0_xxyyy_xz, g_yz_0_xxyyy_yy, g_yz_0_xxyyy_yz, g_yz_0_xxyyy_zz, g_yz_0_xyyy_xx, g_yz_0_xyyy_xxx, g_yz_0_xyyy_xxy, g_yz_0_xyyy_xxz, g_yz_0_xyyy_xy, g_yz_0_xyyy_xyy, g_yz_0_xyyy_xyz, g_yz_0_xyyy_xz, g_yz_0_xyyy_xzz, g_yz_0_xyyy_yy, g_yz_0_xyyy_yz, g_yz_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyy_xx[k] = -g_yz_0_xyyy_xx[k] * ab_x + g_yz_0_xyyy_xxx[k];

                g_yz_0_xxyyy_xy[k] = -g_yz_0_xyyy_xy[k] * ab_x + g_yz_0_xyyy_xxy[k];

                g_yz_0_xxyyy_xz[k] = -g_yz_0_xyyy_xz[k] * ab_x + g_yz_0_xyyy_xxz[k];

                g_yz_0_xxyyy_yy[k] = -g_yz_0_xyyy_yy[k] * ab_x + g_yz_0_xyyy_xyy[k];

                g_yz_0_xxyyy_yz[k] = -g_yz_0_xyyy_yz[k] * ab_x + g_yz_0_xyyy_xyz[k];

                g_yz_0_xxyyy_zz[k] = -g_yz_0_xyyy_zz[k] * ab_x + g_yz_0_xyyy_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 546 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 547 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 548 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 549 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 550 * ccomps * dcomps);

            auto g_yz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyz_xx, g_yz_0_xxyyz_xy, g_yz_0_xxyyz_xz, g_yz_0_xxyyz_yy, g_yz_0_xxyyz_yz, g_yz_0_xxyyz_zz, g_yz_0_xyyz_xx, g_yz_0_xyyz_xxx, g_yz_0_xyyz_xxy, g_yz_0_xyyz_xxz, g_yz_0_xyyz_xy, g_yz_0_xyyz_xyy, g_yz_0_xyyz_xyz, g_yz_0_xyyz_xz, g_yz_0_xyyz_xzz, g_yz_0_xyyz_yy, g_yz_0_xyyz_yz, g_yz_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyz_xx[k] = -g_yz_0_xyyz_xx[k] * ab_x + g_yz_0_xyyz_xxx[k];

                g_yz_0_xxyyz_xy[k] = -g_yz_0_xyyz_xy[k] * ab_x + g_yz_0_xyyz_xxy[k];

                g_yz_0_xxyyz_xz[k] = -g_yz_0_xyyz_xz[k] * ab_x + g_yz_0_xyyz_xxz[k];

                g_yz_0_xxyyz_yy[k] = -g_yz_0_xyyz_yy[k] * ab_x + g_yz_0_xyyz_xyy[k];

                g_yz_0_xxyyz_yz[k] = -g_yz_0_xyyz_yz[k] * ab_x + g_yz_0_xyyz_xyz[k];

                g_yz_0_xxyyz_zz[k] = -g_yz_0_xyyz_zz[k] * ab_x + g_yz_0_xyyz_xzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 552 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 553 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 554 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 555 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 556 * ccomps * dcomps);

            auto g_yz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzz_xx, g_yz_0_xxyzz_xy, g_yz_0_xxyzz_xz, g_yz_0_xxyzz_yy, g_yz_0_xxyzz_yz, g_yz_0_xxyzz_zz, g_yz_0_xyzz_xx, g_yz_0_xyzz_xxx, g_yz_0_xyzz_xxy, g_yz_0_xyzz_xxz, g_yz_0_xyzz_xy, g_yz_0_xyzz_xyy, g_yz_0_xyzz_xyz, g_yz_0_xyzz_xz, g_yz_0_xyzz_xzz, g_yz_0_xyzz_yy, g_yz_0_xyzz_yz, g_yz_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzz_xx[k] = -g_yz_0_xyzz_xx[k] * ab_x + g_yz_0_xyzz_xxx[k];

                g_yz_0_xxyzz_xy[k] = -g_yz_0_xyzz_xy[k] * ab_x + g_yz_0_xyzz_xxy[k];

                g_yz_0_xxyzz_xz[k] = -g_yz_0_xyzz_xz[k] * ab_x + g_yz_0_xyzz_xxz[k];

                g_yz_0_xxyzz_yy[k] = -g_yz_0_xyzz_yy[k] * ab_x + g_yz_0_xyzz_xyy[k];

                g_yz_0_xxyzz_yz[k] = -g_yz_0_xyzz_yz[k] * ab_x + g_yz_0_xyzz_xyz[k];

                g_yz_0_xxyzz_zz[k] = -g_yz_0_xyzz_zz[k] * ab_x + g_yz_0_xyzz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 558 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 559 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 560 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 561 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 562 * ccomps * dcomps);

            auto g_yz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzz_xx, g_yz_0_xxzzz_xy, g_yz_0_xxzzz_xz, g_yz_0_xxzzz_yy, g_yz_0_xxzzz_yz, g_yz_0_xxzzz_zz, g_yz_0_xzzz_xx, g_yz_0_xzzz_xxx, g_yz_0_xzzz_xxy, g_yz_0_xzzz_xxz, g_yz_0_xzzz_xy, g_yz_0_xzzz_xyy, g_yz_0_xzzz_xyz, g_yz_0_xzzz_xz, g_yz_0_xzzz_xzz, g_yz_0_xzzz_yy, g_yz_0_xzzz_yz, g_yz_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzz_xx[k] = -g_yz_0_xzzz_xx[k] * ab_x + g_yz_0_xzzz_xxx[k];

                g_yz_0_xxzzz_xy[k] = -g_yz_0_xzzz_xy[k] * ab_x + g_yz_0_xzzz_xxy[k];

                g_yz_0_xxzzz_xz[k] = -g_yz_0_xzzz_xz[k] * ab_x + g_yz_0_xzzz_xxz[k];

                g_yz_0_xxzzz_yy[k] = -g_yz_0_xzzz_yy[k] * ab_x + g_yz_0_xzzz_xyy[k];

                g_yz_0_xxzzz_yz[k] = -g_yz_0_xzzz_yz[k] * ab_x + g_yz_0_xzzz_xyz[k];

                g_yz_0_xxzzz_zz[k] = -g_yz_0_xzzz_zz[k] * ab_x + g_yz_0_xzzz_xzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 564 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 565 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 566 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 567 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 568 * ccomps * dcomps);

            auto g_yz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyy_xx, g_yz_0_xyyyy_xy, g_yz_0_xyyyy_xz, g_yz_0_xyyyy_yy, g_yz_0_xyyyy_yz, g_yz_0_xyyyy_zz, g_yz_0_yyyy_xx, g_yz_0_yyyy_xxx, g_yz_0_yyyy_xxy, g_yz_0_yyyy_xxz, g_yz_0_yyyy_xy, g_yz_0_yyyy_xyy, g_yz_0_yyyy_xyz, g_yz_0_yyyy_xz, g_yz_0_yyyy_xzz, g_yz_0_yyyy_yy, g_yz_0_yyyy_yz, g_yz_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyy_xx[k] = -g_yz_0_yyyy_xx[k] * ab_x + g_yz_0_yyyy_xxx[k];

                g_yz_0_xyyyy_xy[k] = -g_yz_0_yyyy_xy[k] * ab_x + g_yz_0_yyyy_xxy[k];

                g_yz_0_xyyyy_xz[k] = -g_yz_0_yyyy_xz[k] * ab_x + g_yz_0_yyyy_xxz[k];

                g_yz_0_xyyyy_yy[k] = -g_yz_0_yyyy_yy[k] * ab_x + g_yz_0_yyyy_xyy[k];

                g_yz_0_xyyyy_yz[k] = -g_yz_0_yyyy_yz[k] * ab_x + g_yz_0_yyyy_xyz[k];

                g_yz_0_xyyyy_zz[k] = -g_yz_0_yyyy_zz[k] * ab_x + g_yz_0_yyyy_xzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 570 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 571 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 572 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 573 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 574 * ccomps * dcomps);

            auto g_yz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyz_xx, g_yz_0_xyyyz_xy, g_yz_0_xyyyz_xz, g_yz_0_xyyyz_yy, g_yz_0_xyyyz_yz, g_yz_0_xyyyz_zz, g_yz_0_yyyz_xx, g_yz_0_yyyz_xxx, g_yz_0_yyyz_xxy, g_yz_0_yyyz_xxz, g_yz_0_yyyz_xy, g_yz_0_yyyz_xyy, g_yz_0_yyyz_xyz, g_yz_0_yyyz_xz, g_yz_0_yyyz_xzz, g_yz_0_yyyz_yy, g_yz_0_yyyz_yz, g_yz_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyz_xx[k] = -g_yz_0_yyyz_xx[k] * ab_x + g_yz_0_yyyz_xxx[k];

                g_yz_0_xyyyz_xy[k] = -g_yz_0_yyyz_xy[k] * ab_x + g_yz_0_yyyz_xxy[k];

                g_yz_0_xyyyz_xz[k] = -g_yz_0_yyyz_xz[k] * ab_x + g_yz_0_yyyz_xxz[k];

                g_yz_0_xyyyz_yy[k] = -g_yz_0_yyyz_yy[k] * ab_x + g_yz_0_yyyz_xyy[k];

                g_yz_0_xyyyz_yz[k] = -g_yz_0_yyyz_yz[k] * ab_x + g_yz_0_yyyz_xyz[k];

                g_yz_0_xyyyz_zz[k] = -g_yz_0_yyyz_zz[k] * ab_x + g_yz_0_yyyz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 576 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 577 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 578 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 579 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 580 * ccomps * dcomps);

            auto g_yz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzz_xx, g_yz_0_xyyzz_xy, g_yz_0_xyyzz_xz, g_yz_0_xyyzz_yy, g_yz_0_xyyzz_yz, g_yz_0_xyyzz_zz, g_yz_0_yyzz_xx, g_yz_0_yyzz_xxx, g_yz_0_yyzz_xxy, g_yz_0_yyzz_xxz, g_yz_0_yyzz_xy, g_yz_0_yyzz_xyy, g_yz_0_yyzz_xyz, g_yz_0_yyzz_xz, g_yz_0_yyzz_xzz, g_yz_0_yyzz_yy, g_yz_0_yyzz_yz, g_yz_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzz_xx[k] = -g_yz_0_yyzz_xx[k] * ab_x + g_yz_0_yyzz_xxx[k];

                g_yz_0_xyyzz_xy[k] = -g_yz_0_yyzz_xy[k] * ab_x + g_yz_0_yyzz_xxy[k];

                g_yz_0_xyyzz_xz[k] = -g_yz_0_yyzz_xz[k] * ab_x + g_yz_0_yyzz_xxz[k];

                g_yz_0_xyyzz_yy[k] = -g_yz_0_yyzz_yy[k] * ab_x + g_yz_0_yyzz_xyy[k];

                g_yz_0_xyyzz_yz[k] = -g_yz_0_yyzz_yz[k] * ab_x + g_yz_0_yyzz_xyz[k];

                g_yz_0_xyyzz_zz[k] = -g_yz_0_yyzz_zz[k] * ab_x + g_yz_0_yyzz_xzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 582 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 583 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 584 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 585 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 586 * ccomps * dcomps);

            auto g_yz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzz_xx, g_yz_0_xyzzz_xy, g_yz_0_xyzzz_xz, g_yz_0_xyzzz_yy, g_yz_0_xyzzz_yz, g_yz_0_xyzzz_zz, g_yz_0_yzzz_xx, g_yz_0_yzzz_xxx, g_yz_0_yzzz_xxy, g_yz_0_yzzz_xxz, g_yz_0_yzzz_xy, g_yz_0_yzzz_xyy, g_yz_0_yzzz_xyz, g_yz_0_yzzz_xz, g_yz_0_yzzz_xzz, g_yz_0_yzzz_yy, g_yz_0_yzzz_yz, g_yz_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzz_xx[k] = -g_yz_0_yzzz_xx[k] * ab_x + g_yz_0_yzzz_xxx[k];

                g_yz_0_xyzzz_xy[k] = -g_yz_0_yzzz_xy[k] * ab_x + g_yz_0_yzzz_xxy[k];

                g_yz_0_xyzzz_xz[k] = -g_yz_0_yzzz_xz[k] * ab_x + g_yz_0_yzzz_xxz[k];

                g_yz_0_xyzzz_yy[k] = -g_yz_0_yzzz_yy[k] * ab_x + g_yz_0_yzzz_xyy[k];

                g_yz_0_xyzzz_yz[k] = -g_yz_0_yzzz_yz[k] * ab_x + g_yz_0_yzzz_xyz[k];

                g_yz_0_xyzzz_zz[k] = -g_yz_0_yzzz_zz[k] * ab_x + g_yz_0_yzzz_xzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 588 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 589 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 590 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 591 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 592 * ccomps * dcomps);

            auto g_yz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzz_xx, g_yz_0_xzzzz_xy, g_yz_0_xzzzz_xz, g_yz_0_xzzzz_yy, g_yz_0_xzzzz_yz, g_yz_0_xzzzz_zz, g_yz_0_zzzz_xx, g_yz_0_zzzz_xxx, g_yz_0_zzzz_xxy, g_yz_0_zzzz_xxz, g_yz_0_zzzz_xy, g_yz_0_zzzz_xyy, g_yz_0_zzzz_xyz, g_yz_0_zzzz_xz, g_yz_0_zzzz_xzz, g_yz_0_zzzz_yy, g_yz_0_zzzz_yz, g_yz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzz_xx[k] = -g_yz_0_zzzz_xx[k] * ab_x + g_yz_0_zzzz_xxx[k];

                g_yz_0_xzzzz_xy[k] = -g_yz_0_zzzz_xy[k] * ab_x + g_yz_0_zzzz_xxy[k];

                g_yz_0_xzzzz_xz[k] = -g_yz_0_zzzz_xz[k] * ab_x + g_yz_0_zzzz_xxz[k];

                g_yz_0_xzzzz_yy[k] = -g_yz_0_zzzz_yy[k] * ab_x + g_yz_0_zzzz_xyy[k];

                g_yz_0_xzzzz_yz[k] = -g_yz_0_zzzz_yz[k] * ab_x + g_yz_0_zzzz_xyz[k];

                g_yz_0_xzzzz_zz[k] = -g_yz_0_zzzz_zz[k] * ab_x + g_yz_0_zzzz_xzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 594 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 595 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 596 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 597 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 598 * ccomps * dcomps);

            auto g_yz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyy_xx, g_yz_0_yyyy_xxy, g_yz_0_yyyy_xy, g_yz_0_yyyy_xyy, g_yz_0_yyyy_xyz, g_yz_0_yyyy_xz, g_yz_0_yyyy_yy, g_yz_0_yyyy_yyy, g_yz_0_yyyy_yyz, g_yz_0_yyyy_yz, g_yz_0_yyyy_yzz, g_yz_0_yyyy_zz, g_yz_0_yyyyy_xx, g_yz_0_yyyyy_xy, g_yz_0_yyyyy_xz, g_yz_0_yyyyy_yy, g_yz_0_yyyyy_yz, g_yz_0_yyyyy_zz, g_z_0_yyyy_xx, g_z_0_yyyy_xy, g_z_0_yyyy_xz, g_z_0_yyyy_yy, g_z_0_yyyy_yz, g_z_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyy_xx[k] = -g_z_0_yyyy_xx[k] - g_yz_0_yyyy_xx[k] * ab_y + g_yz_0_yyyy_xxy[k];

                g_yz_0_yyyyy_xy[k] = -g_z_0_yyyy_xy[k] - g_yz_0_yyyy_xy[k] * ab_y + g_yz_0_yyyy_xyy[k];

                g_yz_0_yyyyy_xz[k] = -g_z_0_yyyy_xz[k] - g_yz_0_yyyy_xz[k] * ab_y + g_yz_0_yyyy_xyz[k];

                g_yz_0_yyyyy_yy[k] = -g_z_0_yyyy_yy[k] - g_yz_0_yyyy_yy[k] * ab_y + g_yz_0_yyyy_yyy[k];

                g_yz_0_yyyyy_yz[k] = -g_z_0_yyyy_yz[k] - g_yz_0_yyyy_yz[k] * ab_y + g_yz_0_yyyy_yyz[k];

                g_yz_0_yyyyy_zz[k] = -g_z_0_yyyy_zz[k] - g_yz_0_yyyy_zz[k] * ab_y + g_yz_0_yyyy_yzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyz_xx, g_yz_0_yyyyz_xy, g_yz_0_yyyyz_xz, g_yz_0_yyyyz_yy, g_yz_0_yyyyz_yz, g_yz_0_yyyyz_zz, g_yz_0_yyyz_xx, g_yz_0_yyyz_xxy, g_yz_0_yyyz_xy, g_yz_0_yyyz_xyy, g_yz_0_yyyz_xyz, g_yz_0_yyyz_xz, g_yz_0_yyyz_yy, g_yz_0_yyyz_yyy, g_yz_0_yyyz_yyz, g_yz_0_yyyz_yz, g_yz_0_yyyz_yzz, g_yz_0_yyyz_zz, g_z_0_yyyz_xx, g_z_0_yyyz_xy, g_z_0_yyyz_xz, g_z_0_yyyz_yy, g_z_0_yyyz_yz, g_z_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyz_xx[k] = -g_z_0_yyyz_xx[k] - g_yz_0_yyyz_xx[k] * ab_y + g_yz_0_yyyz_xxy[k];

                g_yz_0_yyyyz_xy[k] = -g_z_0_yyyz_xy[k] - g_yz_0_yyyz_xy[k] * ab_y + g_yz_0_yyyz_xyy[k];

                g_yz_0_yyyyz_xz[k] = -g_z_0_yyyz_xz[k] - g_yz_0_yyyz_xz[k] * ab_y + g_yz_0_yyyz_xyz[k];

                g_yz_0_yyyyz_yy[k] = -g_z_0_yyyz_yy[k] - g_yz_0_yyyz_yy[k] * ab_y + g_yz_0_yyyz_yyy[k];

                g_yz_0_yyyyz_yz[k] = -g_z_0_yyyz_yz[k] - g_yz_0_yyyz_yz[k] * ab_y + g_yz_0_yyyz_yyz[k];

                g_yz_0_yyyyz_zz[k] = -g_z_0_yyyz_zz[k] - g_yz_0_yyyz_zz[k] * ab_y + g_yz_0_yyyz_yzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzz_xx, g_yz_0_yyyzz_xy, g_yz_0_yyyzz_xz, g_yz_0_yyyzz_yy, g_yz_0_yyyzz_yz, g_yz_0_yyyzz_zz, g_yz_0_yyzz_xx, g_yz_0_yyzz_xxy, g_yz_0_yyzz_xy, g_yz_0_yyzz_xyy, g_yz_0_yyzz_xyz, g_yz_0_yyzz_xz, g_yz_0_yyzz_yy, g_yz_0_yyzz_yyy, g_yz_0_yyzz_yyz, g_yz_0_yyzz_yz, g_yz_0_yyzz_yzz, g_yz_0_yyzz_zz, g_z_0_yyzz_xx, g_z_0_yyzz_xy, g_z_0_yyzz_xz, g_z_0_yyzz_yy, g_z_0_yyzz_yz, g_z_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzz_xx[k] = -g_z_0_yyzz_xx[k] - g_yz_0_yyzz_xx[k] * ab_y + g_yz_0_yyzz_xxy[k];

                g_yz_0_yyyzz_xy[k] = -g_z_0_yyzz_xy[k] - g_yz_0_yyzz_xy[k] * ab_y + g_yz_0_yyzz_xyy[k];

                g_yz_0_yyyzz_xz[k] = -g_z_0_yyzz_xz[k] - g_yz_0_yyzz_xz[k] * ab_y + g_yz_0_yyzz_xyz[k];

                g_yz_0_yyyzz_yy[k] = -g_z_0_yyzz_yy[k] - g_yz_0_yyzz_yy[k] * ab_y + g_yz_0_yyzz_yyy[k];

                g_yz_0_yyyzz_yz[k] = -g_z_0_yyzz_yz[k] - g_yz_0_yyzz_yz[k] * ab_y + g_yz_0_yyzz_yyz[k];

                g_yz_0_yyyzz_zz[k] = -g_z_0_yyzz_zz[k] - g_yz_0_yyzz_zz[k] * ab_y + g_yz_0_yyzz_yzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzz_xx, g_yz_0_yyzzz_xy, g_yz_0_yyzzz_xz, g_yz_0_yyzzz_yy, g_yz_0_yyzzz_yz, g_yz_0_yyzzz_zz, g_yz_0_yzzz_xx, g_yz_0_yzzz_xxy, g_yz_0_yzzz_xy, g_yz_0_yzzz_xyy, g_yz_0_yzzz_xyz, g_yz_0_yzzz_xz, g_yz_0_yzzz_yy, g_yz_0_yzzz_yyy, g_yz_0_yzzz_yyz, g_yz_0_yzzz_yz, g_yz_0_yzzz_yzz, g_yz_0_yzzz_zz, g_z_0_yzzz_xx, g_z_0_yzzz_xy, g_z_0_yzzz_xz, g_z_0_yzzz_yy, g_z_0_yzzz_yz, g_z_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzz_xx[k] = -g_z_0_yzzz_xx[k] - g_yz_0_yzzz_xx[k] * ab_y + g_yz_0_yzzz_xxy[k];

                g_yz_0_yyzzz_xy[k] = -g_z_0_yzzz_xy[k] - g_yz_0_yzzz_xy[k] * ab_y + g_yz_0_yzzz_xyy[k];

                g_yz_0_yyzzz_xz[k] = -g_z_0_yzzz_xz[k] - g_yz_0_yzzz_xz[k] * ab_y + g_yz_0_yzzz_xyz[k];

                g_yz_0_yyzzz_yy[k] = -g_z_0_yzzz_yy[k] - g_yz_0_yzzz_yy[k] * ab_y + g_yz_0_yzzz_yyy[k];

                g_yz_0_yyzzz_yz[k] = -g_z_0_yzzz_yz[k] - g_yz_0_yzzz_yz[k] * ab_y + g_yz_0_yzzz_yyz[k];

                g_yz_0_yyzzz_zz[k] = -g_z_0_yzzz_zz[k] - g_yz_0_yzzz_zz[k] * ab_y + g_yz_0_yzzz_yzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzz_xx, g_yz_0_yzzzz_xy, g_yz_0_yzzzz_xz, g_yz_0_yzzzz_yy, g_yz_0_yzzzz_yz, g_yz_0_yzzzz_zz, g_yz_0_zzzz_xx, g_yz_0_zzzz_xxy, g_yz_0_zzzz_xy, g_yz_0_zzzz_xyy, g_yz_0_zzzz_xyz, g_yz_0_zzzz_xz, g_yz_0_zzzz_yy, g_yz_0_zzzz_yyy, g_yz_0_zzzz_yyz, g_yz_0_zzzz_yz, g_yz_0_zzzz_yzz, g_yz_0_zzzz_zz, g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzz_xx[k] = -g_z_0_zzzz_xx[k] - g_yz_0_zzzz_xx[k] * ab_y + g_yz_0_zzzz_xxy[k];

                g_yz_0_yzzzz_xy[k] = -g_z_0_zzzz_xy[k] - g_yz_0_zzzz_xy[k] * ab_y + g_yz_0_zzzz_xyy[k];

                g_yz_0_yzzzz_xz[k] = -g_z_0_zzzz_xz[k] - g_yz_0_zzzz_xz[k] * ab_y + g_yz_0_zzzz_xyz[k];

                g_yz_0_yzzzz_yy[k] = -g_z_0_zzzz_yy[k] - g_yz_0_zzzz_yy[k] * ab_y + g_yz_0_zzzz_yyy[k];

                g_yz_0_yzzzz_yz[k] = -g_z_0_zzzz_yz[k] - g_yz_0_zzzz_yz[k] * ab_y + g_yz_0_zzzz_yyz[k];

                g_yz_0_yzzzz_zz[k] = -g_z_0_zzzz_zz[k] - g_yz_0_zzzz_zz[k] * ab_y + g_yz_0_zzzz_yzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_xx, g_y_0_zzzz_xy, g_y_0_zzzz_xz, g_y_0_zzzz_yy, g_y_0_zzzz_yz, g_y_0_zzzz_zz, g_yz_0_zzzz_xx, g_yz_0_zzzz_xxz, g_yz_0_zzzz_xy, g_yz_0_zzzz_xyz, g_yz_0_zzzz_xz, g_yz_0_zzzz_xzz, g_yz_0_zzzz_yy, g_yz_0_zzzz_yyz, g_yz_0_zzzz_yz, g_yz_0_zzzz_yzz, g_yz_0_zzzz_zz, g_yz_0_zzzz_zzz, g_yz_0_zzzzz_xx, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_yy, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzz_xx[k] = -g_y_0_zzzz_xx[k] - g_yz_0_zzzz_xx[k] * ab_z + g_yz_0_zzzz_xxz[k];

                g_yz_0_zzzzz_xy[k] = -g_y_0_zzzz_xy[k] - g_yz_0_zzzz_xy[k] * ab_z + g_yz_0_zzzz_xyz[k];

                g_yz_0_zzzzz_xz[k] = -g_y_0_zzzz_xz[k] - g_yz_0_zzzz_xz[k] * ab_z + g_yz_0_zzzz_xzz[k];

                g_yz_0_zzzzz_yy[k] = -g_y_0_zzzz_yy[k] - g_yz_0_zzzz_yy[k] * ab_z + g_yz_0_zzzz_yyz[k];

                g_yz_0_zzzzz_yz[k] = -g_y_0_zzzz_yz[k] - g_yz_0_zzzz_yz[k] * ab_z + g_yz_0_zzzz_yzz[k];

                g_yz_0_zzzzz_zz[k] = -g_y_0_zzzz_zz[k] - g_yz_0_zzzz_zz[k] * ab_z + g_yz_0_zzzz_zzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxx_xx, g_zz_0_xxxx_xxx, g_zz_0_xxxx_xxy, g_zz_0_xxxx_xxz, g_zz_0_xxxx_xy, g_zz_0_xxxx_xyy, g_zz_0_xxxx_xyz, g_zz_0_xxxx_xz, g_zz_0_xxxx_xzz, g_zz_0_xxxx_yy, g_zz_0_xxxx_yz, g_zz_0_xxxx_zz, g_zz_0_xxxxx_xx, g_zz_0_xxxxx_xy, g_zz_0_xxxxx_xz, g_zz_0_xxxxx_yy, g_zz_0_xxxxx_yz, g_zz_0_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxx_xx[k] = -g_zz_0_xxxx_xx[k] * ab_x + g_zz_0_xxxx_xxx[k];

                g_zz_0_xxxxx_xy[k] = -g_zz_0_xxxx_xy[k] * ab_x + g_zz_0_xxxx_xxy[k];

                g_zz_0_xxxxx_xz[k] = -g_zz_0_xxxx_xz[k] * ab_x + g_zz_0_xxxx_xxz[k];

                g_zz_0_xxxxx_yy[k] = -g_zz_0_xxxx_yy[k] * ab_x + g_zz_0_xxxx_xyy[k];

                g_zz_0_xxxxx_yz[k] = -g_zz_0_xxxx_yz[k] * ab_x + g_zz_0_xxxx_xyz[k];

                g_zz_0_xxxxx_zz[k] = -g_zz_0_xxxx_zz[k] * ab_x + g_zz_0_xxxx_xzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxy_xx, g_zz_0_xxxxy_xy, g_zz_0_xxxxy_xz, g_zz_0_xxxxy_yy, g_zz_0_xxxxy_yz, g_zz_0_xxxxy_zz, g_zz_0_xxxy_xx, g_zz_0_xxxy_xxx, g_zz_0_xxxy_xxy, g_zz_0_xxxy_xxz, g_zz_0_xxxy_xy, g_zz_0_xxxy_xyy, g_zz_0_xxxy_xyz, g_zz_0_xxxy_xz, g_zz_0_xxxy_xzz, g_zz_0_xxxy_yy, g_zz_0_xxxy_yz, g_zz_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxy_xx[k] = -g_zz_0_xxxy_xx[k] * ab_x + g_zz_0_xxxy_xxx[k];

                g_zz_0_xxxxy_xy[k] = -g_zz_0_xxxy_xy[k] * ab_x + g_zz_0_xxxy_xxy[k];

                g_zz_0_xxxxy_xz[k] = -g_zz_0_xxxy_xz[k] * ab_x + g_zz_0_xxxy_xxz[k];

                g_zz_0_xxxxy_yy[k] = -g_zz_0_xxxy_yy[k] * ab_x + g_zz_0_xxxy_xyy[k];

                g_zz_0_xxxxy_yz[k] = -g_zz_0_xxxy_yz[k] * ab_x + g_zz_0_xxxy_xyz[k];

                g_zz_0_xxxxy_zz[k] = -g_zz_0_xxxy_zz[k] * ab_x + g_zz_0_xxxy_xzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxz_xx, g_zz_0_xxxxz_xy, g_zz_0_xxxxz_xz, g_zz_0_xxxxz_yy, g_zz_0_xxxxz_yz, g_zz_0_xxxxz_zz, g_zz_0_xxxz_xx, g_zz_0_xxxz_xxx, g_zz_0_xxxz_xxy, g_zz_0_xxxz_xxz, g_zz_0_xxxz_xy, g_zz_0_xxxz_xyy, g_zz_0_xxxz_xyz, g_zz_0_xxxz_xz, g_zz_0_xxxz_xzz, g_zz_0_xxxz_yy, g_zz_0_xxxz_yz, g_zz_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxz_xx[k] = -g_zz_0_xxxz_xx[k] * ab_x + g_zz_0_xxxz_xxx[k];

                g_zz_0_xxxxz_xy[k] = -g_zz_0_xxxz_xy[k] * ab_x + g_zz_0_xxxz_xxy[k];

                g_zz_0_xxxxz_xz[k] = -g_zz_0_xxxz_xz[k] * ab_x + g_zz_0_xxxz_xxz[k];

                g_zz_0_xxxxz_yy[k] = -g_zz_0_xxxz_yy[k] * ab_x + g_zz_0_xxxz_xyy[k];

                g_zz_0_xxxxz_yz[k] = -g_zz_0_xxxz_yz[k] * ab_x + g_zz_0_xxxz_xyz[k];

                g_zz_0_xxxxz_zz[k] = -g_zz_0_xxxz_zz[k] * ab_x + g_zz_0_xxxz_xzz[k];
            }

            /// Set up 648-654 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 648 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 649 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 650 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 651 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 652 * ccomps * dcomps);

            auto g_zz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyy_xx, g_zz_0_xxxyy_xy, g_zz_0_xxxyy_xz, g_zz_0_xxxyy_yy, g_zz_0_xxxyy_yz, g_zz_0_xxxyy_zz, g_zz_0_xxyy_xx, g_zz_0_xxyy_xxx, g_zz_0_xxyy_xxy, g_zz_0_xxyy_xxz, g_zz_0_xxyy_xy, g_zz_0_xxyy_xyy, g_zz_0_xxyy_xyz, g_zz_0_xxyy_xz, g_zz_0_xxyy_xzz, g_zz_0_xxyy_yy, g_zz_0_xxyy_yz, g_zz_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyy_xx[k] = -g_zz_0_xxyy_xx[k] * ab_x + g_zz_0_xxyy_xxx[k];

                g_zz_0_xxxyy_xy[k] = -g_zz_0_xxyy_xy[k] * ab_x + g_zz_0_xxyy_xxy[k];

                g_zz_0_xxxyy_xz[k] = -g_zz_0_xxyy_xz[k] * ab_x + g_zz_0_xxyy_xxz[k];

                g_zz_0_xxxyy_yy[k] = -g_zz_0_xxyy_yy[k] * ab_x + g_zz_0_xxyy_xyy[k];

                g_zz_0_xxxyy_yz[k] = -g_zz_0_xxyy_yz[k] * ab_x + g_zz_0_xxyy_xyz[k];

                g_zz_0_xxxyy_zz[k] = -g_zz_0_xxyy_zz[k] * ab_x + g_zz_0_xxyy_xzz[k];
            }

            /// Set up 654-660 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 654 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 655 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 656 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 657 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 658 * ccomps * dcomps);

            auto g_zz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyz_xx, g_zz_0_xxxyz_xy, g_zz_0_xxxyz_xz, g_zz_0_xxxyz_yy, g_zz_0_xxxyz_yz, g_zz_0_xxxyz_zz, g_zz_0_xxyz_xx, g_zz_0_xxyz_xxx, g_zz_0_xxyz_xxy, g_zz_0_xxyz_xxz, g_zz_0_xxyz_xy, g_zz_0_xxyz_xyy, g_zz_0_xxyz_xyz, g_zz_0_xxyz_xz, g_zz_0_xxyz_xzz, g_zz_0_xxyz_yy, g_zz_0_xxyz_yz, g_zz_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyz_xx[k] = -g_zz_0_xxyz_xx[k] * ab_x + g_zz_0_xxyz_xxx[k];

                g_zz_0_xxxyz_xy[k] = -g_zz_0_xxyz_xy[k] * ab_x + g_zz_0_xxyz_xxy[k];

                g_zz_0_xxxyz_xz[k] = -g_zz_0_xxyz_xz[k] * ab_x + g_zz_0_xxyz_xxz[k];

                g_zz_0_xxxyz_yy[k] = -g_zz_0_xxyz_yy[k] * ab_x + g_zz_0_xxyz_xyy[k];

                g_zz_0_xxxyz_yz[k] = -g_zz_0_xxyz_yz[k] * ab_x + g_zz_0_xxyz_xyz[k];

                g_zz_0_xxxyz_zz[k] = -g_zz_0_xxyz_zz[k] * ab_x + g_zz_0_xxyz_xzz[k];
            }

            /// Set up 660-666 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 660 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 661 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 662 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 663 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 664 * ccomps * dcomps);

            auto g_zz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzz_xx, g_zz_0_xxxzz_xy, g_zz_0_xxxzz_xz, g_zz_0_xxxzz_yy, g_zz_0_xxxzz_yz, g_zz_0_xxxzz_zz, g_zz_0_xxzz_xx, g_zz_0_xxzz_xxx, g_zz_0_xxzz_xxy, g_zz_0_xxzz_xxz, g_zz_0_xxzz_xy, g_zz_0_xxzz_xyy, g_zz_0_xxzz_xyz, g_zz_0_xxzz_xz, g_zz_0_xxzz_xzz, g_zz_0_xxzz_yy, g_zz_0_xxzz_yz, g_zz_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzz_xx[k] = -g_zz_0_xxzz_xx[k] * ab_x + g_zz_0_xxzz_xxx[k];

                g_zz_0_xxxzz_xy[k] = -g_zz_0_xxzz_xy[k] * ab_x + g_zz_0_xxzz_xxy[k];

                g_zz_0_xxxzz_xz[k] = -g_zz_0_xxzz_xz[k] * ab_x + g_zz_0_xxzz_xxz[k];

                g_zz_0_xxxzz_yy[k] = -g_zz_0_xxzz_yy[k] * ab_x + g_zz_0_xxzz_xyy[k];

                g_zz_0_xxxzz_yz[k] = -g_zz_0_xxzz_yz[k] * ab_x + g_zz_0_xxzz_xyz[k];

                g_zz_0_xxxzz_zz[k] = -g_zz_0_xxzz_zz[k] * ab_x + g_zz_0_xxzz_xzz[k];
            }

            /// Set up 666-672 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 666 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 667 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 668 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 669 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 670 * ccomps * dcomps);

            auto g_zz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyy_xx, g_zz_0_xxyyy_xy, g_zz_0_xxyyy_xz, g_zz_0_xxyyy_yy, g_zz_0_xxyyy_yz, g_zz_0_xxyyy_zz, g_zz_0_xyyy_xx, g_zz_0_xyyy_xxx, g_zz_0_xyyy_xxy, g_zz_0_xyyy_xxz, g_zz_0_xyyy_xy, g_zz_0_xyyy_xyy, g_zz_0_xyyy_xyz, g_zz_0_xyyy_xz, g_zz_0_xyyy_xzz, g_zz_0_xyyy_yy, g_zz_0_xyyy_yz, g_zz_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyy_xx[k] = -g_zz_0_xyyy_xx[k] * ab_x + g_zz_0_xyyy_xxx[k];

                g_zz_0_xxyyy_xy[k] = -g_zz_0_xyyy_xy[k] * ab_x + g_zz_0_xyyy_xxy[k];

                g_zz_0_xxyyy_xz[k] = -g_zz_0_xyyy_xz[k] * ab_x + g_zz_0_xyyy_xxz[k];

                g_zz_0_xxyyy_yy[k] = -g_zz_0_xyyy_yy[k] * ab_x + g_zz_0_xyyy_xyy[k];

                g_zz_0_xxyyy_yz[k] = -g_zz_0_xyyy_yz[k] * ab_x + g_zz_0_xyyy_xyz[k];

                g_zz_0_xxyyy_zz[k] = -g_zz_0_xyyy_zz[k] * ab_x + g_zz_0_xyyy_xzz[k];
            }

            /// Set up 672-678 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 672 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 673 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 674 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 675 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 676 * ccomps * dcomps);

            auto g_zz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyz_xx, g_zz_0_xxyyz_xy, g_zz_0_xxyyz_xz, g_zz_0_xxyyz_yy, g_zz_0_xxyyz_yz, g_zz_0_xxyyz_zz, g_zz_0_xyyz_xx, g_zz_0_xyyz_xxx, g_zz_0_xyyz_xxy, g_zz_0_xyyz_xxz, g_zz_0_xyyz_xy, g_zz_0_xyyz_xyy, g_zz_0_xyyz_xyz, g_zz_0_xyyz_xz, g_zz_0_xyyz_xzz, g_zz_0_xyyz_yy, g_zz_0_xyyz_yz, g_zz_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyz_xx[k] = -g_zz_0_xyyz_xx[k] * ab_x + g_zz_0_xyyz_xxx[k];

                g_zz_0_xxyyz_xy[k] = -g_zz_0_xyyz_xy[k] * ab_x + g_zz_0_xyyz_xxy[k];

                g_zz_0_xxyyz_xz[k] = -g_zz_0_xyyz_xz[k] * ab_x + g_zz_0_xyyz_xxz[k];

                g_zz_0_xxyyz_yy[k] = -g_zz_0_xyyz_yy[k] * ab_x + g_zz_0_xyyz_xyy[k];

                g_zz_0_xxyyz_yz[k] = -g_zz_0_xyyz_yz[k] * ab_x + g_zz_0_xyyz_xyz[k];

                g_zz_0_xxyyz_zz[k] = -g_zz_0_xyyz_zz[k] * ab_x + g_zz_0_xyyz_xzz[k];
            }

            /// Set up 678-684 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 678 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 679 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 680 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 681 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 682 * ccomps * dcomps);

            auto g_zz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzz_xx, g_zz_0_xxyzz_xy, g_zz_0_xxyzz_xz, g_zz_0_xxyzz_yy, g_zz_0_xxyzz_yz, g_zz_0_xxyzz_zz, g_zz_0_xyzz_xx, g_zz_0_xyzz_xxx, g_zz_0_xyzz_xxy, g_zz_0_xyzz_xxz, g_zz_0_xyzz_xy, g_zz_0_xyzz_xyy, g_zz_0_xyzz_xyz, g_zz_0_xyzz_xz, g_zz_0_xyzz_xzz, g_zz_0_xyzz_yy, g_zz_0_xyzz_yz, g_zz_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzz_xx[k] = -g_zz_0_xyzz_xx[k] * ab_x + g_zz_0_xyzz_xxx[k];

                g_zz_0_xxyzz_xy[k] = -g_zz_0_xyzz_xy[k] * ab_x + g_zz_0_xyzz_xxy[k];

                g_zz_0_xxyzz_xz[k] = -g_zz_0_xyzz_xz[k] * ab_x + g_zz_0_xyzz_xxz[k];

                g_zz_0_xxyzz_yy[k] = -g_zz_0_xyzz_yy[k] * ab_x + g_zz_0_xyzz_xyy[k];

                g_zz_0_xxyzz_yz[k] = -g_zz_0_xyzz_yz[k] * ab_x + g_zz_0_xyzz_xyz[k];

                g_zz_0_xxyzz_zz[k] = -g_zz_0_xyzz_zz[k] * ab_x + g_zz_0_xyzz_xzz[k];
            }

            /// Set up 684-690 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 684 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 685 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 686 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 687 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 688 * ccomps * dcomps);

            auto g_zz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzz_xx, g_zz_0_xxzzz_xy, g_zz_0_xxzzz_xz, g_zz_0_xxzzz_yy, g_zz_0_xxzzz_yz, g_zz_0_xxzzz_zz, g_zz_0_xzzz_xx, g_zz_0_xzzz_xxx, g_zz_0_xzzz_xxy, g_zz_0_xzzz_xxz, g_zz_0_xzzz_xy, g_zz_0_xzzz_xyy, g_zz_0_xzzz_xyz, g_zz_0_xzzz_xz, g_zz_0_xzzz_xzz, g_zz_0_xzzz_yy, g_zz_0_xzzz_yz, g_zz_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzz_xx[k] = -g_zz_0_xzzz_xx[k] * ab_x + g_zz_0_xzzz_xxx[k];

                g_zz_0_xxzzz_xy[k] = -g_zz_0_xzzz_xy[k] * ab_x + g_zz_0_xzzz_xxy[k];

                g_zz_0_xxzzz_xz[k] = -g_zz_0_xzzz_xz[k] * ab_x + g_zz_0_xzzz_xxz[k];

                g_zz_0_xxzzz_yy[k] = -g_zz_0_xzzz_yy[k] * ab_x + g_zz_0_xzzz_xyy[k];

                g_zz_0_xxzzz_yz[k] = -g_zz_0_xzzz_yz[k] * ab_x + g_zz_0_xzzz_xyz[k];

                g_zz_0_xxzzz_zz[k] = -g_zz_0_xzzz_zz[k] * ab_x + g_zz_0_xzzz_xzz[k];
            }

            /// Set up 690-696 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 690 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 691 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 692 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 693 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 694 * ccomps * dcomps);

            auto g_zz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyy_xx, g_zz_0_xyyyy_xy, g_zz_0_xyyyy_xz, g_zz_0_xyyyy_yy, g_zz_0_xyyyy_yz, g_zz_0_xyyyy_zz, g_zz_0_yyyy_xx, g_zz_0_yyyy_xxx, g_zz_0_yyyy_xxy, g_zz_0_yyyy_xxz, g_zz_0_yyyy_xy, g_zz_0_yyyy_xyy, g_zz_0_yyyy_xyz, g_zz_0_yyyy_xz, g_zz_0_yyyy_xzz, g_zz_0_yyyy_yy, g_zz_0_yyyy_yz, g_zz_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyy_xx[k] = -g_zz_0_yyyy_xx[k] * ab_x + g_zz_0_yyyy_xxx[k];

                g_zz_0_xyyyy_xy[k] = -g_zz_0_yyyy_xy[k] * ab_x + g_zz_0_yyyy_xxy[k];

                g_zz_0_xyyyy_xz[k] = -g_zz_0_yyyy_xz[k] * ab_x + g_zz_0_yyyy_xxz[k];

                g_zz_0_xyyyy_yy[k] = -g_zz_0_yyyy_yy[k] * ab_x + g_zz_0_yyyy_xyy[k];

                g_zz_0_xyyyy_yz[k] = -g_zz_0_yyyy_yz[k] * ab_x + g_zz_0_yyyy_xyz[k];

                g_zz_0_xyyyy_zz[k] = -g_zz_0_yyyy_zz[k] * ab_x + g_zz_0_yyyy_xzz[k];
            }

            /// Set up 696-702 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 696 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 697 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 698 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 699 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 700 * ccomps * dcomps);

            auto g_zz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyz_xx, g_zz_0_xyyyz_xy, g_zz_0_xyyyz_xz, g_zz_0_xyyyz_yy, g_zz_0_xyyyz_yz, g_zz_0_xyyyz_zz, g_zz_0_yyyz_xx, g_zz_0_yyyz_xxx, g_zz_0_yyyz_xxy, g_zz_0_yyyz_xxz, g_zz_0_yyyz_xy, g_zz_0_yyyz_xyy, g_zz_0_yyyz_xyz, g_zz_0_yyyz_xz, g_zz_0_yyyz_xzz, g_zz_0_yyyz_yy, g_zz_0_yyyz_yz, g_zz_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyz_xx[k] = -g_zz_0_yyyz_xx[k] * ab_x + g_zz_0_yyyz_xxx[k];

                g_zz_0_xyyyz_xy[k] = -g_zz_0_yyyz_xy[k] * ab_x + g_zz_0_yyyz_xxy[k];

                g_zz_0_xyyyz_xz[k] = -g_zz_0_yyyz_xz[k] * ab_x + g_zz_0_yyyz_xxz[k];

                g_zz_0_xyyyz_yy[k] = -g_zz_0_yyyz_yy[k] * ab_x + g_zz_0_yyyz_xyy[k];

                g_zz_0_xyyyz_yz[k] = -g_zz_0_yyyz_yz[k] * ab_x + g_zz_0_yyyz_xyz[k];

                g_zz_0_xyyyz_zz[k] = -g_zz_0_yyyz_zz[k] * ab_x + g_zz_0_yyyz_xzz[k];
            }

            /// Set up 702-708 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 702 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 703 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 704 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 705 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 706 * ccomps * dcomps);

            auto g_zz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzz_xx, g_zz_0_xyyzz_xy, g_zz_0_xyyzz_xz, g_zz_0_xyyzz_yy, g_zz_0_xyyzz_yz, g_zz_0_xyyzz_zz, g_zz_0_yyzz_xx, g_zz_0_yyzz_xxx, g_zz_0_yyzz_xxy, g_zz_0_yyzz_xxz, g_zz_0_yyzz_xy, g_zz_0_yyzz_xyy, g_zz_0_yyzz_xyz, g_zz_0_yyzz_xz, g_zz_0_yyzz_xzz, g_zz_0_yyzz_yy, g_zz_0_yyzz_yz, g_zz_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzz_xx[k] = -g_zz_0_yyzz_xx[k] * ab_x + g_zz_0_yyzz_xxx[k];

                g_zz_0_xyyzz_xy[k] = -g_zz_0_yyzz_xy[k] * ab_x + g_zz_0_yyzz_xxy[k];

                g_zz_0_xyyzz_xz[k] = -g_zz_0_yyzz_xz[k] * ab_x + g_zz_0_yyzz_xxz[k];

                g_zz_0_xyyzz_yy[k] = -g_zz_0_yyzz_yy[k] * ab_x + g_zz_0_yyzz_xyy[k];

                g_zz_0_xyyzz_yz[k] = -g_zz_0_yyzz_yz[k] * ab_x + g_zz_0_yyzz_xyz[k];

                g_zz_0_xyyzz_zz[k] = -g_zz_0_yyzz_zz[k] * ab_x + g_zz_0_yyzz_xzz[k];
            }

            /// Set up 708-714 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 708 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 709 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 710 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 711 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 712 * ccomps * dcomps);

            auto g_zz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzz_xx, g_zz_0_xyzzz_xy, g_zz_0_xyzzz_xz, g_zz_0_xyzzz_yy, g_zz_0_xyzzz_yz, g_zz_0_xyzzz_zz, g_zz_0_yzzz_xx, g_zz_0_yzzz_xxx, g_zz_0_yzzz_xxy, g_zz_0_yzzz_xxz, g_zz_0_yzzz_xy, g_zz_0_yzzz_xyy, g_zz_0_yzzz_xyz, g_zz_0_yzzz_xz, g_zz_0_yzzz_xzz, g_zz_0_yzzz_yy, g_zz_0_yzzz_yz, g_zz_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzz_xx[k] = -g_zz_0_yzzz_xx[k] * ab_x + g_zz_0_yzzz_xxx[k];

                g_zz_0_xyzzz_xy[k] = -g_zz_0_yzzz_xy[k] * ab_x + g_zz_0_yzzz_xxy[k];

                g_zz_0_xyzzz_xz[k] = -g_zz_0_yzzz_xz[k] * ab_x + g_zz_0_yzzz_xxz[k];

                g_zz_0_xyzzz_yy[k] = -g_zz_0_yzzz_yy[k] * ab_x + g_zz_0_yzzz_xyy[k];

                g_zz_0_xyzzz_yz[k] = -g_zz_0_yzzz_yz[k] * ab_x + g_zz_0_yzzz_xyz[k];

                g_zz_0_xyzzz_zz[k] = -g_zz_0_yzzz_zz[k] * ab_x + g_zz_0_yzzz_xzz[k];
            }

            /// Set up 714-720 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 714 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 715 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 716 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 717 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 718 * ccomps * dcomps);

            auto g_zz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzz_xx, g_zz_0_xzzzz_xy, g_zz_0_xzzzz_xz, g_zz_0_xzzzz_yy, g_zz_0_xzzzz_yz, g_zz_0_xzzzz_zz, g_zz_0_zzzz_xx, g_zz_0_zzzz_xxx, g_zz_0_zzzz_xxy, g_zz_0_zzzz_xxz, g_zz_0_zzzz_xy, g_zz_0_zzzz_xyy, g_zz_0_zzzz_xyz, g_zz_0_zzzz_xz, g_zz_0_zzzz_xzz, g_zz_0_zzzz_yy, g_zz_0_zzzz_yz, g_zz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzz_xx[k] = -g_zz_0_zzzz_xx[k] * ab_x + g_zz_0_zzzz_xxx[k];

                g_zz_0_xzzzz_xy[k] = -g_zz_0_zzzz_xy[k] * ab_x + g_zz_0_zzzz_xxy[k];

                g_zz_0_xzzzz_xz[k] = -g_zz_0_zzzz_xz[k] * ab_x + g_zz_0_zzzz_xxz[k];

                g_zz_0_xzzzz_yy[k] = -g_zz_0_zzzz_yy[k] * ab_x + g_zz_0_zzzz_xyy[k];

                g_zz_0_xzzzz_yz[k] = -g_zz_0_zzzz_yz[k] * ab_x + g_zz_0_zzzz_xyz[k];

                g_zz_0_xzzzz_zz[k] = -g_zz_0_zzzz_zz[k] * ab_x + g_zz_0_zzzz_xzz[k];
            }

            /// Set up 720-726 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 720 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 721 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 722 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 723 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 724 * ccomps * dcomps);

            auto g_zz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyy_xx, g_zz_0_yyyy_xxy, g_zz_0_yyyy_xy, g_zz_0_yyyy_xyy, g_zz_0_yyyy_xyz, g_zz_0_yyyy_xz, g_zz_0_yyyy_yy, g_zz_0_yyyy_yyy, g_zz_0_yyyy_yyz, g_zz_0_yyyy_yz, g_zz_0_yyyy_yzz, g_zz_0_yyyy_zz, g_zz_0_yyyyy_xx, g_zz_0_yyyyy_xy, g_zz_0_yyyyy_xz, g_zz_0_yyyyy_yy, g_zz_0_yyyyy_yz, g_zz_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyy_xx[k] = -g_zz_0_yyyy_xx[k] * ab_y + g_zz_0_yyyy_xxy[k];

                g_zz_0_yyyyy_xy[k] = -g_zz_0_yyyy_xy[k] * ab_y + g_zz_0_yyyy_xyy[k];

                g_zz_0_yyyyy_xz[k] = -g_zz_0_yyyy_xz[k] * ab_y + g_zz_0_yyyy_xyz[k];

                g_zz_0_yyyyy_yy[k] = -g_zz_0_yyyy_yy[k] * ab_y + g_zz_0_yyyy_yyy[k];

                g_zz_0_yyyyy_yz[k] = -g_zz_0_yyyy_yz[k] * ab_y + g_zz_0_yyyy_yyz[k];

                g_zz_0_yyyyy_zz[k] = -g_zz_0_yyyy_zz[k] * ab_y + g_zz_0_yyyy_yzz[k];
            }

            /// Set up 726-732 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 726 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 727 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 728 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 729 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 730 * ccomps * dcomps);

            auto g_zz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyz_xx, g_zz_0_yyyyz_xy, g_zz_0_yyyyz_xz, g_zz_0_yyyyz_yy, g_zz_0_yyyyz_yz, g_zz_0_yyyyz_zz, g_zz_0_yyyz_xx, g_zz_0_yyyz_xxy, g_zz_0_yyyz_xy, g_zz_0_yyyz_xyy, g_zz_0_yyyz_xyz, g_zz_0_yyyz_xz, g_zz_0_yyyz_yy, g_zz_0_yyyz_yyy, g_zz_0_yyyz_yyz, g_zz_0_yyyz_yz, g_zz_0_yyyz_yzz, g_zz_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyz_xx[k] = -g_zz_0_yyyz_xx[k] * ab_y + g_zz_0_yyyz_xxy[k];

                g_zz_0_yyyyz_xy[k] = -g_zz_0_yyyz_xy[k] * ab_y + g_zz_0_yyyz_xyy[k];

                g_zz_0_yyyyz_xz[k] = -g_zz_0_yyyz_xz[k] * ab_y + g_zz_0_yyyz_xyz[k];

                g_zz_0_yyyyz_yy[k] = -g_zz_0_yyyz_yy[k] * ab_y + g_zz_0_yyyz_yyy[k];

                g_zz_0_yyyyz_yz[k] = -g_zz_0_yyyz_yz[k] * ab_y + g_zz_0_yyyz_yyz[k];

                g_zz_0_yyyyz_zz[k] = -g_zz_0_yyyz_zz[k] * ab_y + g_zz_0_yyyz_yzz[k];
            }

            /// Set up 732-738 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 732 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 733 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 734 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 735 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 736 * ccomps * dcomps);

            auto g_zz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzz_xx, g_zz_0_yyyzz_xy, g_zz_0_yyyzz_xz, g_zz_0_yyyzz_yy, g_zz_0_yyyzz_yz, g_zz_0_yyyzz_zz, g_zz_0_yyzz_xx, g_zz_0_yyzz_xxy, g_zz_0_yyzz_xy, g_zz_0_yyzz_xyy, g_zz_0_yyzz_xyz, g_zz_0_yyzz_xz, g_zz_0_yyzz_yy, g_zz_0_yyzz_yyy, g_zz_0_yyzz_yyz, g_zz_0_yyzz_yz, g_zz_0_yyzz_yzz, g_zz_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzz_xx[k] = -g_zz_0_yyzz_xx[k] * ab_y + g_zz_0_yyzz_xxy[k];

                g_zz_0_yyyzz_xy[k] = -g_zz_0_yyzz_xy[k] * ab_y + g_zz_0_yyzz_xyy[k];

                g_zz_0_yyyzz_xz[k] = -g_zz_0_yyzz_xz[k] * ab_y + g_zz_0_yyzz_xyz[k];

                g_zz_0_yyyzz_yy[k] = -g_zz_0_yyzz_yy[k] * ab_y + g_zz_0_yyzz_yyy[k];

                g_zz_0_yyyzz_yz[k] = -g_zz_0_yyzz_yz[k] * ab_y + g_zz_0_yyzz_yyz[k];

                g_zz_0_yyyzz_zz[k] = -g_zz_0_yyzz_zz[k] * ab_y + g_zz_0_yyzz_yzz[k];
            }

            /// Set up 738-744 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 738 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 739 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 740 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 741 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 742 * ccomps * dcomps);

            auto g_zz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzz_xx, g_zz_0_yyzzz_xy, g_zz_0_yyzzz_xz, g_zz_0_yyzzz_yy, g_zz_0_yyzzz_yz, g_zz_0_yyzzz_zz, g_zz_0_yzzz_xx, g_zz_0_yzzz_xxy, g_zz_0_yzzz_xy, g_zz_0_yzzz_xyy, g_zz_0_yzzz_xyz, g_zz_0_yzzz_xz, g_zz_0_yzzz_yy, g_zz_0_yzzz_yyy, g_zz_0_yzzz_yyz, g_zz_0_yzzz_yz, g_zz_0_yzzz_yzz, g_zz_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzz_xx[k] = -g_zz_0_yzzz_xx[k] * ab_y + g_zz_0_yzzz_xxy[k];

                g_zz_0_yyzzz_xy[k] = -g_zz_0_yzzz_xy[k] * ab_y + g_zz_0_yzzz_xyy[k];

                g_zz_0_yyzzz_xz[k] = -g_zz_0_yzzz_xz[k] * ab_y + g_zz_0_yzzz_xyz[k];

                g_zz_0_yyzzz_yy[k] = -g_zz_0_yzzz_yy[k] * ab_y + g_zz_0_yzzz_yyy[k];

                g_zz_0_yyzzz_yz[k] = -g_zz_0_yzzz_yz[k] * ab_y + g_zz_0_yzzz_yyz[k];

                g_zz_0_yyzzz_zz[k] = -g_zz_0_yzzz_zz[k] * ab_y + g_zz_0_yzzz_yzz[k];
            }

            /// Set up 744-750 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 744 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 745 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 746 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 747 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 748 * ccomps * dcomps);

            auto g_zz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzz_xx, g_zz_0_yzzzz_xy, g_zz_0_yzzzz_xz, g_zz_0_yzzzz_yy, g_zz_0_yzzzz_yz, g_zz_0_yzzzz_zz, g_zz_0_zzzz_xx, g_zz_0_zzzz_xxy, g_zz_0_zzzz_xy, g_zz_0_zzzz_xyy, g_zz_0_zzzz_xyz, g_zz_0_zzzz_xz, g_zz_0_zzzz_yy, g_zz_0_zzzz_yyy, g_zz_0_zzzz_yyz, g_zz_0_zzzz_yz, g_zz_0_zzzz_yzz, g_zz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzz_xx[k] = -g_zz_0_zzzz_xx[k] * ab_y + g_zz_0_zzzz_xxy[k];

                g_zz_0_yzzzz_xy[k] = -g_zz_0_zzzz_xy[k] * ab_y + g_zz_0_zzzz_xyy[k];

                g_zz_0_yzzzz_xz[k] = -g_zz_0_zzzz_xz[k] * ab_y + g_zz_0_zzzz_xyz[k];

                g_zz_0_yzzzz_yy[k] = -g_zz_0_zzzz_yy[k] * ab_y + g_zz_0_zzzz_yyy[k];

                g_zz_0_yzzzz_yz[k] = -g_zz_0_zzzz_yz[k] * ab_y + g_zz_0_zzzz_yyz[k];

                g_zz_0_yzzzz_zz[k] = -g_zz_0_zzzz_zz[k] * ab_y + g_zz_0_zzzz_yzz[k];
            }

            /// Set up 750-756 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz, g_zz_0_zzzz_xx, g_zz_0_zzzz_xxz, g_zz_0_zzzz_xy, g_zz_0_zzzz_xyz, g_zz_0_zzzz_xz, g_zz_0_zzzz_xzz, g_zz_0_zzzz_yy, g_zz_0_zzzz_yyz, g_zz_0_zzzz_yz, g_zz_0_zzzz_yzz, g_zz_0_zzzz_zz, g_zz_0_zzzz_zzz, g_zz_0_zzzzz_xx, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_yy, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzz_xx[k] = -2.0 * g_z_0_zzzz_xx[k] - g_zz_0_zzzz_xx[k] * ab_z + g_zz_0_zzzz_xxz[k];

                g_zz_0_zzzzz_xy[k] = -2.0 * g_z_0_zzzz_xy[k] - g_zz_0_zzzz_xy[k] * ab_z + g_zz_0_zzzz_xyz[k];

                g_zz_0_zzzzz_xz[k] = -2.0 * g_z_0_zzzz_xz[k] - g_zz_0_zzzz_xz[k] * ab_z + g_zz_0_zzzz_xzz[k];

                g_zz_0_zzzzz_yy[k] = -2.0 * g_z_0_zzzz_yy[k] - g_zz_0_zzzz_yy[k] * ab_z + g_zz_0_zzzz_yyz[k];

                g_zz_0_zzzzz_yz[k] = -2.0 * g_z_0_zzzz_yz[k] - g_zz_0_zzzz_yz[k] * ab_z + g_zz_0_zzzz_yzz[k];

                g_zz_0_zzzzz_zz[k] = -2.0 * g_z_0_zzzz_zz[k] - g_zz_0_zzzz_zz[k] * ab_z + g_zz_0_zzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

