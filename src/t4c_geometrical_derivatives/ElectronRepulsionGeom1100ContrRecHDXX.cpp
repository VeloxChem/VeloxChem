#include "ElectronRepulsionGeom1100ContrRecHDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_hdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_hdxx,
                                            const size_t idx_geom_01_gdxx,
                                            const size_t idx_geom_10_gdxx,
                                            const size_t idx_geom_11_gdxx,
                                            const size_t idx_geom_11_gfxx,
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

            const auto gd_geom_01_off = idx_geom_01_gdxx + i * dcomps + j;

            auto g_0_x_xxxx_xx = cbuffer.data(gd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xy = cbuffer.data(gd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xz = cbuffer.data(gd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_yy = cbuffer.data(gd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_yz = cbuffer.data(gd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_zz = cbuffer.data(gd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxy_xx = cbuffer.data(gd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxy_xy = cbuffer.data(gd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxy_xz = cbuffer.data(gd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxy_yy = cbuffer.data(gd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxy_yz = cbuffer.data(gd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxy_zz = cbuffer.data(gd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxz_xx = cbuffer.data(gd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxz_xy = cbuffer.data(gd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxz_xz = cbuffer.data(gd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxz_yy = cbuffer.data(gd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxz_yz = cbuffer.data(gd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxz_zz = cbuffer.data(gd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxyy_xx = cbuffer.data(gd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxyy_xy = cbuffer.data(gd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxyy_xz = cbuffer.data(gd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxyy_yy = cbuffer.data(gd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxyy_yz = cbuffer.data(gd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxyy_zz = cbuffer.data(gd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxyz_xx = cbuffer.data(gd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxyz_xy = cbuffer.data(gd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxyz_xz = cbuffer.data(gd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxyz_yy = cbuffer.data(gd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxyz_yz = cbuffer.data(gd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxyz_zz = cbuffer.data(gd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxzz_xx = cbuffer.data(gd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxzz_xy = cbuffer.data(gd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxzz_xz = cbuffer.data(gd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxzz_yy = cbuffer.data(gd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxzz_yz = cbuffer.data(gd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxzz_zz = cbuffer.data(gd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xyyy_xx = cbuffer.data(gd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyyy_xy = cbuffer.data(gd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyyy_xz = cbuffer.data(gd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyyy_yy = cbuffer.data(gd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyyy_yz = cbuffer.data(gd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyyy_zz = cbuffer.data(gd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xyyz_xx = cbuffer.data(gd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xyyz_xy = cbuffer.data(gd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xyyz_xz = cbuffer.data(gd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyyz_yy = cbuffer.data(gd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyyz_yz = cbuffer.data(gd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyyz_zz = cbuffer.data(gd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyzz_xx = cbuffer.data(gd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyzz_xy = cbuffer.data(gd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyzz_xz = cbuffer.data(gd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyzz_yy = cbuffer.data(gd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyzz_yz = cbuffer.data(gd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyzz_zz = cbuffer.data(gd_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xzzz_xx = cbuffer.data(gd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xzzz_xy = cbuffer.data(gd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xzzz_xz = cbuffer.data(gd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xzzz_yy = cbuffer.data(gd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xzzz_yz = cbuffer.data(gd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xzzz_zz = cbuffer.data(gd_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yyyy_xx = cbuffer.data(gd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yyyy_xy = cbuffer.data(gd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yyyy_xz = cbuffer.data(gd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyyy_yy = cbuffer.data(gd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyyy_yz = cbuffer.data(gd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyyy_zz = cbuffer.data(gd_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yyyz_xx = cbuffer.data(gd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyyz_xy = cbuffer.data(gd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyyz_xz = cbuffer.data(gd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyyz_yy = cbuffer.data(gd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyyz_yz = cbuffer.data(gd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyyz_zz = cbuffer.data(gd_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yyzz_xx = cbuffer.data(gd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyzz_xy = cbuffer.data(gd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyzz_xz = cbuffer.data(gd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyzz_yy = cbuffer.data(gd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyzz_yz = cbuffer.data(gd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyzz_zz = cbuffer.data(gd_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yzzz_xx = cbuffer.data(gd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yzzz_xy = cbuffer.data(gd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzzz_xz = cbuffer.data(gd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yzzz_yy = cbuffer.data(gd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yzzz_yz = cbuffer.data(gd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yzzz_zz = cbuffer.data(gd_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_zzzz_xx = cbuffer.data(gd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zzzz_xy = cbuffer.data(gd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zzzz_xz = cbuffer.data(gd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zzzz_yy = cbuffer.data(gd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zzzz_yz = cbuffer.data(gd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zzzz_zz = cbuffer.data(gd_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xxxx_xx = cbuffer.data(gd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxxx_xy = cbuffer.data(gd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxxx_xz = cbuffer.data(gd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xxxx_yy = cbuffer.data(gd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xxxx_yz = cbuffer.data(gd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xxxx_zz = cbuffer.data(gd_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xxxy_xx = cbuffer.data(gd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xxxy_xy = cbuffer.data(gd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xxxy_xz = cbuffer.data(gd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xxxy_yy = cbuffer.data(gd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxxy_yz = cbuffer.data(gd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxxy_zz = cbuffer.data(gd_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xxxz_xx = cbuffer.data(gd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxxz_xy = cbuffer.data(gd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxxz_xz = cbuffer.data(gd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxxz_yy = cbuffer.data(gd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxxz_yz = cbuffer.data(gd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxxz_zz = cbuffer.data(gd_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xxyy_xx = cbuffer.data(gd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxyy_xy = cbuffer.data(gd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxyy_xz = cbuffer.data(gd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxyy_yy = cbuffer.data(gd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxyy_yz = cbuffer.data(gd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxyy_zz = cbuffer.data(gd_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xxyz_xx = cbuffer.data(gd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxyz_xy = cbuffer.data(gd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxyz_xz = cbuffer.data(gd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxyz_yy = cbuffer.data(gd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxyz_yz = cbuffer.data(gd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxyz_zz = cbuffer.data(gd_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xxzz_xx = cbuffer.data(gd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxzz_xy = cbuffer.data(gd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxzz_xz = cbuffer.data(gd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxzz_yy = cbuffer.data(gd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxzz_yz = cbuffer.data(gd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxzz_zz = cbuffer.data(gd_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xyyy_xx = cbuffer.data(gd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xyyy_xy = cbuffer.data(gd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xyyy_xz = cbuffer.data(gd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xyyy_yy = cbuffer.data(gd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyyy_yz = cbuffer.data(gd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyyy_zz = cbuffer.data(gd_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xyyz_xx = cbuffer.data(gd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyyz_xy = cbuffer.data(gd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyyz_xz = cbuffer.data(gd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyyz_yy = cbuffer.data(gd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyyz_yz = cbuffer.data(gd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyyz_zz = cbuffer.data(gd_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xyzz_xx = cbuffer.data(gd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyzz_xy = cbuffer.data(gd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyzz_xz = cbuffer.data(gd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyzz_yy = cbuffer.data(gd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyzz_yz = cbuffer.data(gd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyzz_zz = cbuffer.data(gd_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xzzz_xx = cbuffer.data(gd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xzzz_xy = cbuffer.data(gd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xzzz_xz = cbuffer.data(gd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xzzz_yy = cbuffer.data(gd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xzzz_yz = cbuffer.data(gd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xzzz_zz = cbuffer.data(gd_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yyyy_xx = cbuffer.data(gd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yyyy_xy = cbuffer.data(gd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yyyy_xz = cbuffer.data(gd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yyyy_yy = cbuffer.data(gd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yyyy_yz = cbuffer.data(gd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yyyy_zz = cbuffer.data(gd_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yyyz_xx = cbuffer.data(gd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yyyz_xy = cbuffer.data(gd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yyyz_xz = cbuffer.data(gd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yyyz_yy = cbuffer.data(gd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyyz_yz = cbuffer.data(gd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyyz_zz = cbuffer.data(gd_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yyzz_xx = cbuffer.data(gd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yyzz_xy = cbuffer.data(gd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yyzz_xz = cbuffer.data(gd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_yyzz_yy = cbuffer.data(gd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_yyzz_yz = cbuffer.data(gd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_yyzz_zz = cbuffer.data(gd_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_yzzz_xx = cbuffer.data(gd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_yzzz_xy = cbuffer.data(gd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_yzzz_xz = cbuffer.data(gd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_yzzz_yy = cbuffer.data(gd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_yzzz_yz = cbuffer.data(gd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_yzzz_zz = cbuffer.data(gd_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_zzzz_xx = cbuffer.data(gd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zzzz_xy = cbuffer.data(gd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zzzz_xz = cbuffer.data(gd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zzzz_yy = cbuffer.data(gd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zzzz_yz = cbuffer.data(gd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zzzz_zz = cbuffer.data(gd_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xxxx_xx = cbuffer.data(gd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xxxx_xy = cbuffer.data(gd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xxxx_xz = cbuffer.data(gd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xxxx_yy = cbuffer.data(gd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xxxx_yz = cbuffer.data(gd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xxxx_zz = cbuffer.data(gd_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xxxy_xx = cbuffer.data(gd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xxxy_xy = cbuffer.data(gd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xxxy_xz = cbuffer.data(gd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xxxy_yy = cbuffer.data(gd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xxxy_yz = cbuffer.data(gd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xxxy_zz = cbuffer.data(gd_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xxxz_xx = cbuffer.data(gd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xxxz_xy = cbuffer.data(gd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xxxz_xz = cbuffer.data(gd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xxxz_yy = cbuffer.data(gd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xxxz_yz = cbuffer.data(gd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xxxz_zz = cbuffer.data(gd_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xxyy_xx = cbuffer.data(gd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xxyy_xy = cbuffer.data(gd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxyy_xz = cbuffer.data(gd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxyy_yy = cbuffer.data(gd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxyy_yz = cbuffer.data(gd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxyy_zz = cbuffer.data(gd_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xxyz_xx = cbuffer.data(gd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxyz_xy = cbuffer.data(gd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxyz_xz = cbuffer.data(gd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxyz_yy = cbuffer.data(gd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxyz_yz = cbuffer.data(gd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxyz_zz = cbuffer.data(gd_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xxzz_xx = cbuffer.data(gd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxzz_xy = cbuffer.data(gd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxzz_xz = cbuffer.data(gd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xxzz_yy = cbuffer.data(gd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xxzz_yz = cbuffer.data(gd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xxzz_zz = cbuffer.data(gd_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xyyy_xx = cbuffer.data(gd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xyyy_xy = cbuffer.data(gd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xyyy_xz = cbuffer.data(gd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xyyy_yy = cbuffer.data(gd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xyyy_yz = cbuffer.data(gd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xyyy_zz = cbuffer.data(gd_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xyyz_xx = cbuffer.data(gd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xyyz_xy = cbuffer.data(gd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xyyz_xz = cbuffer.data(gd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xyyz_yy = cbuffer.data(gd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xyyz_yz = cbuffer.data(gd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xyyz_zz = cbuffer.data(gd_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_xyzz_xx = cbuffer.data(gd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xyzz_xy = cbuffer.data(gd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xyzz_xz = cbuffer.data(gd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_xyzz_yy = cbuffer.data(gd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xyzz_yz = cbuffer.data(gd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xyzz_zz = cbuffer.data(gd_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_xzzz_xx = cbuffer.data(gd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xzzz_xy = cbuffer.data(gd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xzzz_xz = cbuffer.data(gd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_xzzz_yy = cbuffer.data(gd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xzzz_yz = cbuffer.data(gd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xzzz_zz = cbuffer.data(gd_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yyyy_xx = cbuffer.data(gd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yyyy_xy = cbuffer.data(gd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yyyy_xz = cbuffer.data(gd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yyyy_yy = cbuffer.data(gd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yyyy_yz = cbuffer.data(gd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yyyy_zz = cbuffer.data(gd_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yyyz_xx = cbuffer.data(gd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yyyz_xy = cbuffer.data(gd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yyyz_xz = cbuffer.data(gd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yyyz_yy = cbuffer.data(gd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yyyz_yz = cbuffer.data(gd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yyyz_zz = cbuffer.data(gd_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_yyzz_xx = cbuffer.data(gd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yyzz_xy = cbuffer.data(gd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yyzz_xz = cbuffer.data(gd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_yyzz_yy = cbuffer.data(gd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_yyzz_yz = cbuffer.data(gd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_yyzz_zz = cbuffer.data(gd_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_yzzz_xx = cbuffer.data(gd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_yzzz_xy = cbuffer.data(gd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_yzzz_xz = cbuffer.data(gd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_yzzz_yy = cbuffer.data(gd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_yzzz_yz = cbuffer.data(gd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_yzzz_zz = cbuffer.data(gd_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_zzzz_xx = cbuffer.data(gd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zzzz_xy = cbuffer.data(gd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zzzz_xz = cbuffer.data(gd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zzzz_yy = cbuffer.data(gd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zzzz_yz = cbuffer.data(gd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zzzz_zz = cbuffer.data(gd_geom_01_off + 269 * ccomps * dcomps);

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

            const auto gd_geom_11_off = idx_geom_11_gdxx + i * dcomps + j;

            auto g_x_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 419 * ccomps * dcomps);

            auto g_y_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 566 * ccomps * dcomps);

            auto g_z_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 567 * ccomps * dcomps);

            auto g_z_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 568 * ccomps * dcomps);

            auto g_z_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 569 * ccomps * dcomps);

            auto g_z_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 570 * ccomps * dcomps);

            auto g_z_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 571 * ccomps * dcomps);

            auto g_z_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 572 * ccomps * dcomps);

            auto g_z_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 573 * ccomps * dcomps);

            auto g_z_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 574 * ccomps * dcomps);

            auto g_z_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 575 * ccomps * dcomps);

            auto g_z_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 576 * ccomps * dcomps);

            auto g_z_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 577 * ccomps * dcomps);

            auto g_z_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 578 * ccomps * dcomps);

            auto g_z_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 579 * ccomps * dcomps);

            auto g_z_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 580 * ccomps * dcomps);

            auto g_z_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 581 * ccomps * dcomps);

            auto g_z_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 582 * ccomps * dcomps);

            auto g_z_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 583 * ccomps * dcomps);

            auto g_z_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 584 * ccomps * dcomps);

            auto g_z_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 585 * ccomps * dcomps);

            auto g_z_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 586 * ccomps * dcomps);

            auto g_z_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 587 * ccomps * dcomps);

            auto g_z_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 588 * ccomps * dcomps);

            auto g_z_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 589 * ccomps * dcomps);

            auto g_z_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 590 * ccomps * dcomps);

            auto g_z_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 591 * ccomps * dcomps);

            auto g_z_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 592 * ccomps * dcomps);

            auto g_z_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 593 * ccomps * dcomps);

            auto g_z_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 594 * ccomps * dcomps);

            auto g_z_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 595 * ccomps * dcomps);

            auto g_z_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 596 * ccomps * dcomps);

            auto g_z_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 597 * ccomps * dcomps);

            auto g_z_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 598 * ccomps * dcomps);

            auto g_z_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 599 * ccomps * dcomps);

            auto g_z_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 605 * ccomps * dcomps);

            auto g_z_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 609 * ccomps * dcomps);

            auto g_z_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 611 * ccomps * dcomps);

            auto g_z_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 614 * ccomps * dcomps);

            auto g_z_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 615 * ccomps * dcomps);

            auto g_z_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 617 * ccomps * dcomps);

            auto g_z_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 619 * ccomps * dcomps);

            auto g_z_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 623 * ccomps * dcomps);

            auto g_z_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 629 * ccomps * dcomps);

            auto g_z_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 635 * ccomps * dcomps);

            auto g_z_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 639 * ccomps * dcomps);

            auto g_z_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 641 * ccomps * dcomps);

            auto g_z_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 643 * ccomps * dcomps);

            auto g_z_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 644 * ccomps * dcomps);

            auto g_z_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 659 * ccomps * dcomps);

            auto g_z_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 671 * ccomps * dcomps);

            auto g_z_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 689 * ccomps * dcomps);

            auto g_z_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 699 * ccomps * dcomps);

            auto g_z_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 719 * ccomps * dcomps);

            auto g_z_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 749 * ccomps * dcomps);

            auto g_z_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 779 * ccomps * dcomps);

            auto g_z_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 809 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GFSS

            const auto gf_geom_11_off = idx_geom_11_gfxx + i * dcomps + j;

            auto g_x_x_xxxx_xxx = cbuffer.data(gf_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxx_xxy = cbuffer.data(gf_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxx_xxz = cbuffer.data(gf_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxxx_xyy = cbuffer.data(gf_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxx_xyz = cbuffer.data(gf_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxx_xzz = cbuffer.data(gf_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxxx_yyy = cbuffer.data(gf_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxx_yyz = cbuffer.data(gf_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxx_yzz = cbuffer.data(gf_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxxx_zzz = cbuffer.data(gf_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxxy_xxx = cbuffer.data(gf_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxxy_xxy = cbuffer.data(gf_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxxy_xxz = cbuffer.data(gf_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxxy_xyy = cbuffer.data(gf_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxxy_xyz = cbuffer.data(gf_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxxy_xzz = cbuffer.data(gf_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxxy_yyy = cbuffer.data(gf_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxxy_yyz = cbuffer.data(gf_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xxxy_yzz = cbuffer.data(gf_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxxy_zzz = cbuffer.data(gf_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xxxz_xxx = cbuffer.data(gf_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxxz_xxy = cbuffer.data(gf_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxxz_xxz = cbuffer.data(gf_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxxz_xyy = cbuffer.data(gf_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xxxz_xyz = cbuffer.data(gf_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxxz_xzz = cbuffer.data(gf_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxxz_yyy = cbuffer.data(gf_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxxz_yyz = cbuffer.data(gf_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxxz_yzz = cbuffer.data(gf_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxxz_zzz = cbuffer.data(gf_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xxyy_xxx = cbuffer.data(gf_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xxyy_xxy = cbuffer.data(gf_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xxyy_xxz = cbuffer.data(gf_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xxyy_xyy = cbuffer.data(gf_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xxyy_xyz = cbuffer.data(gf_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xxyy_xzz = cbuffer.data(gf_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xxyy_yyy = cbuffer.data(gf_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xxyy_yyz = cbuffer.data(gf_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xxyy_yzz = cbuffer.data(gf_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xxyy_zzz = cbuffer.data(gf_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xxyz_xxx = cbuffer.data(gf_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xxyz_xxy = cbuffer.data(gf_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xxyz_xxz = cbuffer.data(gf_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xxyz_xyy = cbuffer.data(gf_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xxyz_xyz = cbuffer.data(gf_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xxyz_xzz = cbuffer.data(gf_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xxyz_yyy = cbuffer.data(gf_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xxyz_yyz = cbuffer.data(gf_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xxyz_yzz = cbuffer.data(gf_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xxyz_zzz = cbuffer.data(gf_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xxzz_xxx = cbuffer.data(gf_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xxzz_xxy = cbuffer.data(gf_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xxzz_xxz = cbuffer.data(gf_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xxzz_xyy = cbuffer.data(gf_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xxzz_xyz = cbuffer.data(gf_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xxzz_xzz = cbuffer.data(gf_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xxzz_yyy = cbuffer.data(gf_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xxzz_yyz = cbuffer.data(gf_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xxzz_yzz = cbuffer.data(gf_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xxzz_zzz = cbuffer.data(gf_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_xyyy_xxx = cbuffer.data(gf_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xyyy_xxy = cbuffer.data(gf_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xyyy_xxz = cbuffer.data(gf_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_xyyy_xyy = cbuffer.data(gf_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_xyyy_xyz = cbuffer.data(gf_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_xyyy_xzz = cbuffer.data(gf_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_xyyy_yyy = cbuffer.data(gf_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_xyyy_yyz = cbuffer.data(gf_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_xyyy_yzz = cbuffer.data(gf_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_xyyy_zzz = cbuffer.data(gf_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_xyyz_xxx = cbuffer.data(gf_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_xyyz_xxy = cbuffer.data(gf_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_xyyz_xxz = cbuffer.data(gf_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_xyyz_xyy = cbuffer.data(gf_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_xyyz_xyz = cbuffer.data(gf_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_xyyz_xzz = cbuffer.data(gf_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_xyyz_yyy = cbuffer.data(gf_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_xyyz_yyz = cbuffer.data(gf_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_xyyz_yzz = cbuffer.data(gf_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_xyyz_zzz = cbuffer.data(gf_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_xyzz_xxx = cbuffer.data(gf_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_xyzz_xxy = cbuffer.data(gf_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_xyzz_xxz = cbuffer.data(gf_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_xyzz_xyy = cbuffer.data(gf_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_xyzz_xyz = cbuffer.data(gf_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_xyzz_xzz = cbuffer.data(gf_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_xyzz_yyy = cbuffer.data(gf_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_xyzz_yyz = cbuffer.data(gf_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_xyzz_yzz = cbuffer.data(gf_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_xyzz_zzz = cbuffer.data(gf_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_x_xzzz_xxx = cbuffer.data(gf_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_xzzz_xxy = cbuffer.data(gf_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_xzzz_xxz = cbuffer.data(gf_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_xzzz_xyy = cbuffer.data(gf_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_xzzz_xyz = cbuffer.data(gf_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_xzzz_xzz = cbuffer.data(gf_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_xzzz_yyy = cbuffer.data(gf_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_xzzz_yyz = cbuffer.data(gf_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_xzzz_yzz = cbuffer.data(gf_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_xzzz_zzz = cbuffer.data(gf_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yyyy_xxx = cbuffer.data(gf_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yyyy_xxy = cbuffer.data(gf_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_yyyy_xxz = cbuffer.data(gf_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yyyy_xyy = cbuffer.data(gf_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yyyy_xyz = cbuffer.data(gf_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_x_yyyy_xzz = cbuffer.data(gf_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_yyyy_yyy = cbuffer.data(gf_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_yyyy_yyz = cbuffer.data(gf_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_x_yyyy_yzz = cbuffer.data(gf_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_yyyy_zzz = cbuffer.data(gf_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_yyyz_xxx = cbuffer.data(gf_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_yyyz_xxy = cbuffer.data(gf_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_x_yyyz_xxz = cbuffer.data(gf_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_yyyz_xyy = cbuffer.data(gf_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_x_yyyz_xyz = cbuffer.data(gf_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_yyyz_xzz = cbuffer.data(gf_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_yyyz_yyy = cbuffer.data(gf_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_yyyz_yyz = cbuffer.data(gf_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_yyyz_yzz = cbuffer.data(gf_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_yyyz_zzz = cbuffer.data(gf_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_x_yyzz_xxx = cbuffer.data(gf_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_yyzz_xxy = cbuffer.data(gf_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_yyzz_xxz = cbuffer.data(gf_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_yyzz_xyy = cbuffer.data(gf_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_yyzz_xyz = cbuffer.data(gf_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_yyzz_xzz = cbuffer.data(gf_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_x_yyzz_yyy = cbuffer.data(gf_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_x_yyzz_yyz = cbuffer.data(gf_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_x_yyzz_yzz = cbuffer.data(gf_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_x_yyzz_zzz = cbuffer.data(gf_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_x_yzzz_xxx = cbuffer.data(gf_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_x_yzzz_xxy = cbuffer.data(gf_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_x_yzzz_xxz = cbuffer.data(gf_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_x_yzzz_xyy = cbuffer.data(gf_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_x_yzzz_xyz = cbuffer.data(gf_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_x_yzzz_xzz = cbuffer.data(gf_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_x_yzzz_yyy = cbuffer.data(gf_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_x_yzzz_yyz = cbuffer.data(gf_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_x_yzzz_yzz = cbuffer.data(gf_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_x_yzzz_zzz = cbuffer.data(gf_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_x_zzzz_xxx = cbuffer.data(gf_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_x_zzzz_xxy = cbuffer.data(gf_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_x_zzzz_xxz = cbuffer.data(gf_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_x_zzzz_xyy = cbuffer.data(gf_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_x_zzzz_xyz = cbuffer.data(gf_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_x_zzzz_xzz = cbuffer.data(gf_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_x_zzzz_yyy = cbuffer.data(gf_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_x_zzzz_yyz = cbuffer.data(gf_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_x_zzzz_yzz = cbuffer.data(gf_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_x_zzzz_zzz = cbuffer.data(gf_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_xxxx_xxx = cbuffer.data(gf_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xxxx_xxy = cbuffer.data(gf_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xxxx_xxz = cbuffer.data(gf_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xxxx_xyy = cbuffer.data(gf_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xxxx_xyz = cbuffer.data(gf_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xxxx_xzz = cbuffer.data(gf_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_xxxx_yyy = cbuffer.data(gf_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xxxx_yyz = cbuffer.data(gf_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xxxx_yzz = cbuffer.data(gf_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xxxx_zzz = cbuffer.data(gf_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_xxxy_xxx = cbuffer.data(gf_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_xxxy_xxy = cbuffer.data(gf_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_xxxy_xxz = cbuffer.data(gf_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_xxxy_xyy = cbuffer.data(gf_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_xxxy_xyz = cbuffer.data(gf_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_xxxy_xzz = cbuffer.data(gf_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_xxxy_yyy = cbuffer.data(gf_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_xxxy_yyz = cbuffer.data(gf_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_xxxy_yzz = cbuffer.data(gf_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xxxy_zzz = cbuffer.data(gf_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xxxz_xxx = cbuffer.data(gf_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xxxz_xxy = cbuffer.data(gf_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xxxz_xxz = cbuffer.data(gf_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xxxz_xyy = cbuffer.data(gf_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_xxxz_xyz = cbuffer.data(gf_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xxxz_xzz = cbuffer.data(gf_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xxxz_yyy = cbuffer.data(gf_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xxxz_yyz = cbuffer.data(gf_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xxxz_yzz = cbuffer.data(gf_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xxxz_zzz = cbuffer.data(gf_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_y_xxyy_xxx = cbuffer.data(gf_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xxyy_xxy = cbuffer.data(gf_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xxyy_xxz = cbuffer.data(gf_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xxyy_xyy = cbuffer.data(gf_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xxyy_xyz = cbuffer.data(gf_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xxyy_xzz = cbuffer.data(gf_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_xxyy_yyy = cbuffer.data(gf_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xxyy_yyz = cbuffer.data(gf_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xxyy_yzz = cbuffer.data(gf_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_xxyy_zzz = cbuffer.data(gf_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_xxyz_xxx = cbuffer.data(gf_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_xxyz_xxy = cbuffer.data(gf_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_xxyz_xxz = cbuffer.data(gf_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_xxyz_xyy = cbuffer.data(gf_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_xxyz_xyz = cbuffer.data(gf_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_xxyz_xzz = cbuffer.data(gf_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_xxyz_yyy = cbuffer.data(gf_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_xxyz_yyz = cbuffer.data(gf_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_xxyz_yzz = cbuffer.data(gf_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_xxyz_zzz = cbuffer.data(gf_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_xxzz_xxx = cbuffer.data(gf_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_xxzz_xxy = cbuffer.data(gf_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_xxzz_xxz = cbuffer.data(gf_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_xxzz_xyy = cbuffer.data(gf_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_xxzz_xyz = cbuffer.data(gf_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_xxzz_xzz = cbuffer.data(gf_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_xxzz_yyy = cbuffer.data(gf_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_xxzz_yyz = cbuffer.data(gf_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_xxzz_yzz = cbuffer.data(gf_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_xxzz_zzz = cbuffer.data(gf_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_y_xyyy_xxx = cbuffer.data(gf_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_xyyy_xxy = cbuffer.data(gf_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_xyyy_xxz = cbuffer.data(gf_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_xyyy_xyy = cbuffer.data(gf_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_xyyy_xyz = cbuffer.data(gf_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_xyyy_xzz = cbuffer.data(gf_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_y_xyyy_yyy = cbuffer.data(gf_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_xyyy_yyz = cbuffer.data(gf_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_xyyy_yzz = cbuffer.data(gf_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_xyyy_zzz = cbuffer.data(gf_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_xyyz_xxx = cbuffer.data(gf_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_xyyz_xxy = cbuffer.data(gf_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_y_xyyz_xxz = cbuffer.data(gf_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_xyyz_xyy = cbuffer.data(gf_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_y_xyyz_xyz = cbuffer.data(gf_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_y_xyyz_xzz = cbuffer.data(gf_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_xyyz_yyy = cbuffer.data(gf_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_xyyz_yyz = cbuffer.data(gf_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_y_xyyz_yzz = cbuffer.data(gf_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_xyyz_zzz = cbuffer.data(gf_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_xyzz_xxx = cbuffer.data(gf_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_y_xyzz_xxy = cbuffer.data(gf_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_xyzz_xxz = cbuffer.data(gf_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_xyzz_xyy = cbuffer.data(gf_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_y_xyzz_xyz = cbuffer.data(gf_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_xyzz_xzz = cbuffer.data(gf_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_xyzz_yyy = cbuffer.data(gf_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_xyzz_yyz = cbuffer.data(gf_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_xyzz_yzz = cbuffer.data(gf_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_xyzz_zzz = cbuffer.data(gf_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_y_xzzz_xxx = cbuffer.data(gf_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_xzzz_xxy = cbuffer.data(gf_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_xzzz_xxz = cbuffer.data(gf_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_xzzz_xyy = cbuffer.data(gf_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_xzzz_xyz = cbuffer.data(gf_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_xzzz_xzz = cbuffer.data(gf_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_y_xzzz_yyy = cbuffer.data(gf_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_xzzz_yyz = cbuffer.data(gf_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_xzzz_yzz = cbuffer.data(gf_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_xzzz_zzz = cbuffer.data(gf_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_yyyy_xxx = cbuffer.data(gf_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_yyyy_xxy = cbuffer.data(gf_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_y_yyyy_xxz = cbuffer.data(gf_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_y_yyyy_xyy = cbuffer.data(gf_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_y_yyyy_xyz = cbuffer.data(gf_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_y_yyyy_xzz = cbuffer.data(gf_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_y_yyyy_yyy = cbuffer.data(gf_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_y_yyyy_yyz = cbuffer.data(gf_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_y_yyyy_yzz = cbuffer.data(gf_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_y_yyyy_zzz = cbuffer.data(gf_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_y_yyyz_xxx = cbuffer.data(gf_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_y_yyyz_xxy = cbuffer.data(gf_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_y_yyyz_xxz = cbuffer.data(gf_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_y_yyyz_xyy = cbuffer.data(gf_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_y_yyyz_xyz = cbuffer.data(gf_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_y_yyyz_xzz = cbuffer.data(gf_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_y_yyyz_yyy = cbuffer.data(gf_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_y_yyyz_yyz = cbuffer.data(gf_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_y_yyyz_yzz = cbuffer.data(gf_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_y_yyyz_zzz = cbuffer.data(gf_geom_11_off + 269 * ccomps * dcomps);

            auto g_x_y_yyzz_xxx = cbuffer.data(gf_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_y_yyzz_xxy = cbuffer.data(gf_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_y_yyzz_xxz = cbuffer.data(gf_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_y_yyzz_xyy = cbuffer.data(gf_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_y_yyzz_xyz = cbuffer.data(gf_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_y_yyzz_xzz = cbuffer.data(gf_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_y_yyzz_yyy = cbuffer.data(gf_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_y_yyzz_yyz = cbuffer.data(gf_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_y_yyzz_yzz = cbuffer.data(gf_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_y_yyzz_zzz = cbuffer.data(gf_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_y_yzzz_xxx = cbuffer.data(gf_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_y_yzzz_xxy = cbuffer.data(gf_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_y_yzzz_xxz = cbuffer.data(gf_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_y_yzzz_xyy = cbuffer.data(gf_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_y_yzzz_xyz = cbuffer.data(gf_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_y_yzzz_xzz = cbuffer.data(gf_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_y_yzzz_yyy = cbuffer.data(gf_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_y_yzzz_yyz = cbuffer.data(gf_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_y_yzzz_yzz = cbuffer.data(gf_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_y_yzzz_zzz = cbuffer.data(gf_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_y_zzzz_xxx = cbuffer.data(gf_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_y_zzzz_xxy = cbuffer.data(gf_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_y_zzzz_xxz = cbuffer.data(gf_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_y_zzzz_xyy = cbuffer.data(gf_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_y_zzzz_xyz = cbuffer.data(gf_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_y_zzzz_xzz = cbuffer.data(gf_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_y_zzzz_yyy = cbuffer.data(gf_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_y_zzzz_yyz = cbuffer.data(gf_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_y_zzzz_yzz = cbuffer.data(gf_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_y_zzzz_zzz = cbuffer.data(gf_geom_11_off + 299 * ccomps * dcomps);

            auto g_x_z_xxxx_xxx = cbuffer.data(gf_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_xxxx_xxy = cbuffer.data(gf_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_xxxx_xxz = cbuffer.data(gf_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_xxxx_xyy = cbuffer.data(gf_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_xxxx_xyz = cbuffer.data(gf_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_xxxx_xzz = cbuffer.data(gf_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_z_xxxx_yyy = cbuffer.data(gf_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_xxxx_yyz = cbuffer.data(gf_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_xxxx_yzz = cbuffer.data(gf_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_xxxx_zzz = cbuffer.data(gf_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_xxxy_xxx = cbuffer.data(gf_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_xxxy_xxy = cbuffer.data(gf_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_z_xxxy_xxz = cbuffer.data(gf_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_xxxy_xyy = cbuffer.data(gf_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_xxxy_xyz = cbuffer.data(gf_geom_11_off + 314 * ccomps * dcomps);

            auto g_x_z_xxxy_xzz = cbuffer.data(gf_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_xxxy_yyy = cbuffer.data(gf_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_xxxy_yyz = cbuffer.data(gf_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_z_xxxy_yzz = cbuffer.data(gf_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_xxxy_zzz = cbuffer.data(gf_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_xxxz_xxx = cbuffer.data(gf_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_xxxz_xxy = cbuffer.data(gf_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_xxxz_xxz = cbuffer.data(gf_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_xxxz_xyy = cbuffer.data(gf_geom_11_off + 323 * ccomps * dcomps);

            auto g_x_z_xxxz_xyz = cbuffer.data(gf_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_z_xxxz_xzz = cbuffer.data(gf_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_z_xxxz_yyy = cbuffer.data(gf_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_z_xxxz_yyz = cbuffer.data(gf_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_z_xxxz_yzz = cbuffer.data(gf_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_z_xxxz_zzz = cbuffer.data(gf_geom_11_off + 329 * ccomps * dcomps);

            auto g_x_z_xxyy_xxx = cbuffer.data(gf_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_z_xxyy_xxy = cbuffer.data(gf_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_z_xxyy_xxz = cbuffer.data(gf_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_z_xxyy_xyy = cbuffer.data(gf_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_z_xxyy_xyz = cbuffer.data(gf_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_z_xxyy_xzz = cbuffer.data(gf_geom_11_off + 335 * ccomps * dcomps);

            auto g_x_z_xxyy_yyy = cbuffer.data(gf_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_xxyy_yyz = cbuffer.data(gf_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_xxyy_yzz = cbuffer.data(gf_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_xxyy_zzz = cbuffer.data(gf_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_xxyz_xxx = cbuffer.data(gf_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_xxyz_xxy = cbuffer.data(gf_geom_11_off + 341 * ccomps * dcomps);

            auto g_x_z_xxyz_xxz = cbuffer.data(gf_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_xxyz_xyy = cbuffer.data(gf_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_xxyz_xyz = cbuffer.data(gf_geom_11_off + 344 * ccomps * dcomps);

            auto g_x_z_xxyz_xzz = cbuffer.data(gf_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_xxyz_yyy = cbuffer.data(gf_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_xxyz_yyz = cbuffer.data(gf_geom_11_off + 347 * ccomps * dcomps);

            auto g_x_z_xxyz_yzz = cbuffer.data(gf_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_xxyz_zzz = cbuffer.data(gf_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_xxzz_xxx = cbuffer.data(gf_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_xxzz_xxy = cbuffer.data(gf_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_xxzz_xxz = cbuffer.data(gf_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_xxzz_xyy = cbuffer.data(gf_geom_11_off + 353 * ccomps * dcomps);

            auto g_x_z_xxzz_xyz = cbuffer.data(gf_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_xxzz_xzz = cbuffer.data(gf_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_xxzz_yyy = cbuffer.data(gf_geom_11_off + 356 * ccomps * dcomps);

            auto g_x_z_xxzz_yyz = cbuffer.data(gf_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_xxzz_yzz = cbuffer.data(gf_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_xxzz_zzz = cbuffer.data(gf_geom_11_off + 359 * ccomps * dcomps);

            auto g_x_z_xyyy_xxx = cbuffer.data(gf_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_xyyy_xxy = cbuffer.data(gf_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_xyyy_xxz = cbuffer.data(gf_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_xyyy_xyy = cbuffer.data(gf_geom_11_off + 363 * ccomps * dcomps);

            auto g_x_z_xyyy_xyz = cbuffer.data(gf_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_xyyy_xzz = cbuffer.data(gf_geom_11_off + 365 * ccomps * dcomps);

            auto g_x_z_xyyy_yyy = cbuffer.data(gf_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_xyyy_yyz = cbuffer.data(gf_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_xyyy_yzz = cbuffer.data(gf_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_xyyy_zzz = cbuffer.data(gf_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_xyyz_xxx = cbuffer.data(gf_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_xyyz_xxy = cbuffer.data(gf_geom_11_off + 371 * ccomps * dcomps);

            auto g_x_z_xyyz_xxz = cbuffer.data(gf_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_xyyz_xyy = cbuffer.data(gf_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_xyyz_xyz = cbuffer.data(gf_geom_11_off + 374 * ccomps * dcomps);

            auto g_x_z_xyyz_xzz = cbuffer.data(gf_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_xyyz_yyy = cbuffer.data(gf_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_xyyz_yyz = cbuffer.data(gf_geom_11_off + 377 * ccomps * dcomps);

            auto g_x_z_xyyz_yzz = cbuffer.data(gf_geom_11_off + 378 * ccomps * dcomps);

            auto g_x_z_xyyz_zzz = cbuffer.data(gf_geom_11_off + 379 * ccomps * dcomps);

            auto g_x_z_xyzz_xxx = cbuffer.data(gf_geom_11_off + 380 * ccomps * dcomps);

            auto g_x_z_xyzz_xxy = cbuffer.data(gf_geom_11_off + 381 * ccomps * dcomps);

            auto g_x_z_xyzz_xxz = cbuffer.data(gf_geom_11_off + 382 * ccomps * dcomps);

            auto g_x_z_xyzz_xyy = cbuffer.data(gf_geom_11_off + 383 * ccomps * dcomps);

            auto g_x_z_xyzz_xyz = cbuffer.data(gf_geom_11_off + 384 * ccomps * dcomps);

            auto g_x_z_xyzz_xzz = cbuffer.data(gf_geom_11_off + 385 * ccomps * dcomps);

            auto g_x_z_xyzz_yyy = cbuffer.data(gf_geom_11_off + 386 * ccomps * dcomps);

            auto g_x_z_xyzz_yyz = cbuffer.data(gf_geom_11_off + 387 * ccomps * dcomps);

            auto g_x_z_xyzz_yzz = cbuffer.data(gf_geom_11_off + 388 * ccomps * dcomps);

            auto g_x_z_xyzz_zzz = cbuffer.data(gf_geom_11_off + 389 * ccomps * dcomps);

            auto g_x_z_xzzz_xxx = cbuffer.data(gf_geom_11_off + 390 * ccomps * dcomps);

            auto g_x_z_xzzz_xxy = cbuffer.data(gf_geom_11_off + 391 * ccomps * dcomps);

            auto g_x_z_xzzz_xxz = cbuffer.data(gf_geom_11_off + 392 * ccomps * dcomps);

            auto g_x_z_xzzz_xyy = cbuffer.data(gf_geom_11_off + 393 * ccomps * dcomps);

            auto g_x_z_xzzz_xyz = cbuffer.data(gf_geom_11_off + 394 * ccomps * dcomps);

            auto g_x_z_xzzz_xzz = cbuffer.data(gf_geom_11_off + 395 * ccomps * dcomps);

            auto g_x_z_xzzz_yyy = cbuffer.data(gf_geom_11_off + 396 * ccomps * dcomps);

            auto g_x_z_xzzz_yyz = cbuffer.data(gf_geom_11_off + 397 * ccomps * dcomps);

            auto g_x_z_xzzz_yzz = cbuffer.data(gf_geom_11_off + 398 * ccomps * dcomps);

            auto g_x_z_xzzz_zzz = cbuffer.data(gf_geom_11_off + 399 * ccomps * dcomps);

            auto g_x_z_yyyy_xxx = cbuffer.data(gf_geom_11_off + 400 * ccomps * dcomps);

            auto g_x_z_yyyy_xxy = cbuffer.data(gf_geom_11_off + 401 * ccomps * dcomps);

            auto g_x_z_yyyy_xxz = cbuffer.data(gf_geom_11_off + 402 * ccomps * dcomps);

            auto g_x_z_yyyy_xyy = cbuffer.data(gf_geom_11_off + 403 * ccomps * dcomps);

            auto g_x_z_yyyy_xyz = cbuffer.data(gf_geom_11_off + 404 * ccomps * dcomps);

            auto g_x_z_yyyy_xzz = cbuffer.data(gf_geom_11_off + 405 * ccomps * dcomps);

            auto g_x_z_yyyy_yyy = cbuffer.data(gf_geom_11_off + 406 * ccomps * dcomps);

            auto g_x_z_yyyy_yyz = cbuffer.data(gf_geom_11_off + 407 * ccomps * dcomps);

            auto g_x_z_yyyy_yzz = cbuffer.data(gf_geom_11_off + 408 * ccomps * dcomps);

            auto g_x_z_yyyy_zzz = cbuffer.data(gf_geom_11_off + 409 * ccomps * dcomps);

            auto g_x_z_yyyz_xxx = cbuffer.data(gf_geom_11_off + 410 * ccomps * dcomps);

            auto g_x_z_yyyz_xxy = cbuffer.data(gf_geom_11_off + 411 * ccomps * dcomps);

            auto g_x_z_yyyz_xxz = cbuffer.data(gf_geom_11_off + 412 * ccomps * dcomps);

            auto g_x_z_yyyz_xyy = cbuffer.data(gf_geom_11_off + 413 * ccomps * dcomps);

            auto g_x_z_yyyz_xyz = cbuffer.data(gf_geom_11_off + 414 * ccomps * dcomps);

            auto g_x_z_yyyz_xzz = cbuffer.data(gf_geom_11_off + 415 * ccomps * dcomps);

            auto g_x_z_yyyz_yyy = cbuffer.data(gf_geom_11_off + 416 * ccomps * dcomps);

            auto g_x_z_yyyz_yyz = cbuffer.data(gf_geom_11_off + 417 * ccomps * dcomps);

            auto g_x_z_yyyz_yzz = cbuffer.data(gf_geom_11_off + 418 * ccomps * dcomps);

            auto g_x_z_yyyz_zzz = cbuffer.data(gf_geom_11_off + 419 * ccomps * dcomps);

            auto g_x_z_yyzz_xxx = cbuffer.data(gf_geom_11_off + 420 * ccomps * dcomps);

            auto g_x_z_yyzz_xxy = cbuffer.data(gf_geom_11_off + 421 * ccomps * dcomps);

            auto g_x_z_yyzz_xxz = cbuffer.data(gf_geom_11_off + 422 * ccomps * dcomps);

            auto g_x_z_yyzz_xyy = cbuffer.data(gf_geom_11_off + 423 * ccomps * dcomps);

            auto g_x_z_yyzz_xyz = cbuffer.data(gf_geom_11_off + 424 * ccomps * dcomps);

            auto g_x_z_yyzz_xzz = cbuffer.data(gf_geom_11_off + 425 * ccomps * dcomps);

            auto g_x_z_yyzz_yyy = cbuffer.data(gf_geom_11_off + 426 * ccomps * dcomps);

            auto g_x_z_yyzz_yyz = cbuffer.data(gf_geom_11_off + 427 * ccomps * dcomps);

            auto g_x_z_yyzz_yzz = cbuffer.data(gf_geom_11_off + 428 * ccomps * dcomps);

            auto g_x_z_yyzz_zzz = cbuffer.data(gf_geom_11_off + 429 * ccomps * dcomps);

            auto g_x_z_yzzz_xxx = cbuffer.data(gf_geom_11_off + 430 * ccomps * dcomps);

            auto g_x_z_yzzz_xxy = cbuffer.data(gf_geom_11_off + 431 * ccomps * dcomps);

            auto g_x_z_yzzz_xxz = cbuffer.data(gf_geom_11_off + 432 * ccomps * dcomps);

            auto g_x_z_yzzz_xyy = cbuffer.data(gf_geom_11_off + 433 * ccomps * dcomps);

            auto g_x_z_yzzz_xyz = cbuffer.data(gf_geom_11_off + 434 * ccomps * dcomps);

            auto g_x_z_yzzz_xzz = cbuffer.data(gf_geom_11_off + 435 * ccomps * dcomps);

            auto g_x_z_yzzz_yyy = cbuffer.data(gf_geom_11_off + 436 * ccomps * dcomps);

            auto g_x_z_yzzz_yyz = cbuffer.data(gf_geom_11_off + 437 * ccomps * dcomps);

            auto g_x_z_yzzz_yzz = cbuffer.data(gf_geom_11_off + 438 * ccomps * dcomps);

            auto g_x_z_yzzz_zzz = cbuffer.data(gf_geom_11_off + 439 * ccomps * dcomps);

            auto g_x_z_zzzz_xxx = cbuffer.data(gf_geom_11_off + 440 * ccomps * dcomps);

            auto g_x_z_zzzz_xxy = cbuffer.data(gf_geom_11_off + 441 * ccomps * dcomps);

            auto g_x_z_zzzz_xxz = cbuffer.data(gf_geom_11_off + 442 * ccomps * dcomps);

            auto g_x_z_zzzz_xyy = cbuffer.data(gf_geom_11_off + 443 * ccomps * dcomps);

            auto g_x_z_zzzz_xyz = cbuffer.data(gf_geom_11_off + 444 * ccomps * dcomps);

            auto g_x_z_zzzz_xzz = cbuffer.data(gf_geom_11_off + 445 * ccomps * dcomps);

            auto g_x_z_zzzz_yyy = cbuffer.data(gf_geom_11_off + 446 * ccomps * dcomps);

            auto g_x_z_zzzz_yyz = cbuffer.data(gf_geom_11_off + 447 * ccomps * dcomps);

            auto g_x_z_zzzz_yzz = cbuffer.data(gf_geom_11_off + 448 * ccomps * dcomps);

            auto g_x_z_zzzz_zzz = cbuffer.data(gf_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_x_xxxx_xxx = cbuffer.data(gf_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_x_xxxx_xxy = cbuffer.data(gf_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_x_xxxx_xxz = cbuffer.data(gf_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_x_xxxx_xyy = cbuffer.data(gf_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_x_xxxx_xyz = cbuffer.data(gf_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_x_xxxx_xzz = cbuffer.data(gf_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_x_xxxx_yyy = cbuffer.data(gf_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_x_xxxx_yyz = cbuffer.data(gf_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_x_xxxx_yzz = cbuffer.data(gf_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_x_xxxx_zzz = cbuffer.data(gf_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_x_xxxy_xxx = cbuffer.data(gf_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_x_xxxy_xxy = cbuffer.data(gf_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_x_xxxy_xxz = cbuffer.data(gf_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_x_xxxy_xyy = cbuffer.data(gf_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_x_xxxy_xyz = cbuffer.data(gf_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_x_xxxy_xzz = cbuffer.data(gf_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_x_xxxy_yyy = cbuffer.data(gf_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_x_xxxy_yyz = cbuffer.data(gf_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_x_xxxy_yzz = cbuffer.data(gf_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_x_xxxy_zzz = cbuffer.data(gf_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_x_xxxz_xxx = cbuffer.data(gf_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_x_xxxz_xxy = cbuffer.data(gf_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_x_xxxz_xxz = cbuffer.data(gf_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_x_xxxz_xyy = cbuffer.data(gf_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_x_xxxz_xyz = cbuffer.data(gf_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_x_xxxz_xzz = cbuffer.data(gf_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_x_xxxz_yyy = cbuffer.data(gf_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_x_xxxz_yyz = cbuffer.data(gf_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_x_xxxz_yzz = cbuffer.data(gf_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_x_xxxz_zzz = cbuffer.data(gf_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_x_xxyy_xxx = cbuffer.data(gf_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_x_xxyy_xxy = cbuffer.data(gf_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_x_xxyy_xxz = cbuffer.data(gf_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_x_xxyy_xyy = cbuffer.data(gf_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_x_xxyy_xyz = cbuffer.data(gf_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_x_xxyy_xzz = cbuffer.data(gf_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_x_xxyy_yyy = cbuffer.data(gf_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_x_xxyy_yyz = cbuffer.data(gf_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_x_xxyy_yzz = cbuffer.data(gf_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_x_xxyy_zzz = cbuffer.data(gf_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_x_xxyz_xxx = cbuffer.data(gf_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_x_xxyz_xxy = cbuffer.data(gf_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_x_xxyz_xxz = cbuffer.data(gf_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_x_xxyz_xyy = cbuffer.data(gf_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_x_xxyz_xyz = cbuffer.data(gf_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_x_xxyz_xzz = cbuffer.data(gf_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_x_xxyz_yyy = cbuffer.data(gf_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_x_xxyz_yyz = cbuffer.data(gf_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_x_xxyz_yzz = cbuffer.data(gf_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_x_xxyz_zzz = cbuffer.data(gf_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_x_xxzz_xxx = cbuffer.data(gf_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_x_xxzz_xxy = cbuffer.data(gf_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_x_xxzz_xxz = cbuffer.data(gf_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_x_xxzz_xyy = cbuffer.data(gf_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_x_xxzz_xyz = cbuffer.data(gf_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_x_xxzz_xzz = cbuffer.data(gf_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_x_xxzz_yyy = cbuffer.data(gf_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_x_xxzz_yyz = cbuffer.data(gf_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_x_xxzz_yzz = cbuffer.data(gf_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_x_xxzz_zzz = cbuffer.data(gf_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_x_xyyy_xxx = cbuffer.data(gf_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_x_xyyy_xxy = cbuffer.data(gf_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_x_xyyy_xxz = cbuffer.data(gf_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_x_xyyy_xyy = cbuffer.data(gf_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_x_xyyy_xyz = cbuffer.data(gf_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_x_xyyy_xzz = cbuffer.data(gf_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_x_xyyy_yyy = cbuffer.data(gf_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_x_xyyy_yyz = cbuffer.data(gf_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_x_xyyy_yzz = cbuffer.data(gf_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_x_xyyy_zzz = cbuffer.data(gf_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_x_xyyz_xxx = cbuffer.data(gf_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_x_xyyz_xxy = cbuffer.data(gf_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_x_xyyz_xxz = cbuffer.data(gf_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_x_xyyz_xyy = cbuffer.data(gf_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_x_xyyz_xyz = cbuffer.data(gf_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_x_xyyz_xzz = cbuffer.data(gf_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_x_xyyz_yyy = cbuffer.data(gf_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_x_xyyz_yyz = cbuffer.data(gf_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_x_xyyz_yzz = cbuffer.data(gf_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_x_xyyz_zzz = cbuffer.data(gf_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_x_xyzz_xxx = cbuffer.data(gf_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_x_xyzz_xxy = cbuffer.data(gf_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_x_xyzz_xxz = cbuffer.data(gf_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_x_xyzz_xyy = cbuffer.data(gf_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_x_xyzz_xyz = cbuffer.data(gf_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_x_xyzz_xzz = cbuffer.data(gf_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_x_xyzz_yyy = cbuffer.data(gf_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_x_xyzz_yyz = cbuffer.data(gf_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_x_xyzz_yzz = cbuffer.data(gf_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_x_xyzz_zzz = cbuffer.data(gf_geom_11_off + 539 * ccomps * dcomps);

            auto g_y_x_xzzz_xxx = cbuffer.data(gf_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_x_xzzz_xxy = cbuffer.data(gf_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_x_xzzz_xxz = cbuffer.data(gf_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_x_xzzz_xyy = cbuffer.data(gf_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_x_xzzz_xyz = cbuffer.data(gf_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_x_xzzz_xzz = cbuffer.data(gf_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_x_xzzz_yyy = cbuffer.data(gf_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_x_xzzz_yyz = cbuffer.data(gf_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_x_xzzz_yzz = cbuffer.data(gf_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_x_xzzz_zzz = cbuffer.data(gf_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_x_yyyy_xxx = cbuffer.data(gf_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_x_yyyy_xxy = cbuffer.data(gf_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_x_yyyy_xxz = cbuffer.data(gf_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_x_yyyy_xyy = cbuffer.data(gf_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_x_yyyy_xyz = cbuffer.data(gf_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_x_yyyy_xzz = cbuffer.data(gf_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_x_yyyy_yyy = cbuffer.data(gf_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_x_yyyy_yyz = cbuffer.data(gf_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_x_yyyy_yzz = cbuffer.data(gf_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_x_yyyy_zzz = cbuffer.data(gf_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_x_yyyz_xxx = cbuffer.data(gf_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_x_yyyz_xxy = cbuffer.data(gf_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_x_yyyz_xxz = cbuffer.data(gf_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_x_yyyz_xyy = cbuffer.data(gf_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_x_yyyz_xyz = cbuffer.data(gf_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_x_yyyz_xzz = cbuffer.data(gf_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_x_yyyz_yyy = cbuffer.data(gf_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_x_yyyz_yyz = cbuffer.data(gf_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_x_yyyz_yzz = cbuffer.data(gf_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_x_yyyz_zzz = cbuffer.data(gf_geom_11_off + 569 * ccomps * dcomps);

            auto g_y_x_yyzz_xxx = cbuffer.data(gf_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_x_yyzz_xxy = cbuffer.data(gf_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_x_yyzz_xxz = cbuffer.data(gf_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_x_yyzz_xyy = cbuffer.data(gf_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_x_yyzz_xyz = cbuffer.data(gf_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_x_yyzz_xzz = cbuffer.data(gf_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_x_yyzz_yyy = cbuffer.data(gf_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_x_yyzz_yyz = cbuffer.data(gf_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_x_yyzz_yzz = cbuffer.data(gf_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_x_yyzz_zzz = cbuffer.data(gf_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_x_yzzz_xxx = cbuffer.data(gf_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_x_yzzz_xxy = cbuffer.data(gf_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_x_yzzz_xxz = cbuffer.data(gf_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_x_yzzz_xyy = cbuffer.data(gf_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_x_yzzz_xyz = cbuffer.data(gf_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_x_yzzz_xzz = cbuffer.data(gf_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_x_yzzz_yyy = cbuffer.data(gf_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_x_yzzz_yyz = cbuffer.data(gf_geom_11_off + 587 * ccomps * dcomps);

            auto g_y_x_yzzz_yzz = cbuffer.data(gf_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_x_yzzz_zzz = cbuffer.data(gf_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_x_zzzz_xxx = cbuffer.data(gf_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_x_zzzz_xxy = cbuffer.data(gf_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_x_zzzz_xxz = cbuffer.data(gf_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_x_zzzz_xyy = cbuffer.data(gf_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_x_zzzz_xyz = cbuffer.data(gf_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_x_zzzz_xzz = cbuffer.data(gf_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_x_zzzz_yyy = cbuffer.data(gf_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_x_zzzz_yyz = cbuffer.data(gf_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_x_zzzz_yzz = cbuffer.data(gf_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_x_zzzz_zzz = cbuffer.data(gf_geom_11_off + 599 * ccomps * dcomps);

            auto g_y_y_xxxx_xxx = cbuffer.data(gf_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_y_xxxx_xxy = cbuffer.data(gf_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_y_xxxx_xxz = cbuffer.data(gf_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_y_xxxx_xyy = cbuffer.data(gf_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_y_xxxx_xyz = cbuffer.data(gf_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_y_xxxx_xzz = cbuffer.data(gf_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_y_xxxx_yyy = cbuffer.data(gf_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_y_xxxx_yyz = cbuffer.data(gf_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_y_xxxx_yzz = cbuffer.data(gf_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_y_xxxx_zzz = cbuffer.data(gf_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_y_xxxy_xxx = cbuffer.data(gf_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_y_xxxy_xxy = cbuffer.data(gf_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_y_xxxy_xxz = cbuffer.data(gf_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_y_xxxy_xyy = cbuffer.data(gf_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_y_xxxy_xyz = cbuffer.data(gf_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_y_xxxy_xzz = cbuffer.data(gf_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_y_xxxy_yyy = cbuffer.data(gf_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_y_xxxy_yyz = cbuffer.data(gf_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_y_xxxy_yzz = cbuffer.data(gf_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_y_xxxy_zzz = cbuffer.data(gf_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_y_xxxz_xxx = cbuffer.data(gf_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_y_xxxz_xxy = cbuffer.data(gf_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_y_xxxz_xxz = cbuffer.data(gf_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_y_xxxz_xyy = cbuffer.data(gf_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_y_xxxz_xyz = cbuffer.data(gf_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_y_xxxz_xzz = cbuffer.data(gf_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_y_xxxz_yyy = cbuffer.data(gf_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_y_xxxz_yyz = cbuffer.data(gf_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_y_xxxz_yzz = cbuffer.data(gf_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_y_xxxz_zzz = cbuffer.data(gf_geom_11_off + 629 * ccomps * dcomps);

            auto g_y_y_xxyy_xxx = cbuffer.data(gf_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_y_xxyy_xxy = cbuffer.data(gf_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_y_xxyy_xxz = cbuffer.data(gf_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_y_xxyy_xyy = cbuffer.data(gf_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_y_xxyy_xyz = cbuffer.data(gf_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_y_xxyy_xzz = cbuffer.data(gf_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_y_xxyy_yyy = cbuffer.data(gf_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_y_xxyy_yyz = cbuffer.data(gf_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_y_xxyy_yzz = cbuffer.data(gf_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_y_xxyy_zzz = cbuffer.data(gf_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_y_xxyz_xxx = cbuffer.data(gf_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_y_xxyz_xxy = cbuffer.data(gf_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_y_xxyz_xxz = cbuffer.data(gf_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_y_xxyz_xyy = cbuffer.data(gf_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_y_xxyz_xyz = cbuffer.data(gf_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_y_xxyz_xzz = cbuffer.data(gf_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_y_xxyz_yyy = cbuffer.data(gf_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_y_xxyz_yyz = cbuffer.data(gf_geom_11_off + 647 * ccomps * dcomps);

            auto g_y_y_xxyz_yzz = cbuffer.data(gf_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_y_xxyz_zzz = cbuffer.data(gf_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_y_xxzz_xxx = cbuffer.data(gf_geom_11_off + 650 * ccomps * dcomps);

            auto g_y_y_xxzz_xxy = cbuffer.data(gf_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_y_xxzz_xxz = cbuffer.data(gf_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_y_xxzz_xyy = cbuffer.data(gf_geom_11_off + 653 * ccomps * dcomps);

            auto g_y_y_xxzz_xyz = cbuffer.data(gf_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_y_xxzz_xzz = cbuffer.data(gf_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_y_xxzz_yyy = cbuffer.data(gf_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_y_xxzz_yyz = cbuffer.data(gf_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_y_xxzz_yzz = cbuffer.data(gf_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_y_xxzz_zzz = cbuffer.data(gf_geom_11_off + 659 * ccomps * dcomps);

            auto g_y_y_xyyy_xxx = cbuffer.data(gf_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_y_xyyy_xxy = cbuffer.data(gf_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_y_xyyy_xxz = cbuffer.data(gf_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_y_xyyy_xyy = cbuffer.data(gf_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_y_xyyy_xyz = cbuffer.data(gf_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_y_xyyy_xzz = cbuffer.data(gf_geom_11_off + 665 * ccomps * dcomps);

            auto g_y_y_xyyy_yyy = cbuffer.data(gf_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_y_xyyy_yyz = cbuffer.data(gf_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_y_xyyy_yzz = cbuffer.data(gf_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_y_xyyy_zzz = cbuffer.data(gf_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_y_xyyz_xxx = cbuffer.data(gf_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_y_xyyz_xxy = cbuffer.data(gf_geom_11_off + 671 * ccomps * dcomps);

            auto g_y_y_xyyz_xxz = cbuffer.data(gf_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_y_xyyz_xyy = cbuffer.data(gf_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_y_xyyz_xyz = cbuffer.data(gf_geom_11_off + 674 * ccomps * dcomps);

            auto g_y_y_xyyz_xzz = cbuffer.data(gf_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_y_xyyz_yyy = cbuffer.data(gf_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_y_xyyz_yyz = cbuffer.data(gf_geom_11_off + 677 * ccomps * dcomps);

            auto g_y_y_xyyz_yzz = cbuffer.data(gf_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_y_xyyz_zzz = cbuffer.data(gf_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_y_xyzz_xxx = cbuffer.data(gf_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_y_xyzz_xxy = cbuffer.data(gf_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_y_xyzz_xxz = cbuffer.data(gf_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_y_xyzz_xyy = cbuffer.data(gf_geom_11_off + 683 * ccomps * dcomps);

            auto g_y_y_xyzz_xyz = cbuffer.data(gf_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_y_xyzz_xzz = cbuffer.data(gf_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_y_xyzz_yyy = cbuffer.data(gf_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_y_xyzz_yyz = cbuffer.data(gf_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_y_xyzz_yzz = cbuffer.data(gf_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_y_xyzz_zzz = cbuffer.data(gf_geom_11_off + 689 * ccomps * dcomps);

            auto g_y_y_xzzz_xxx = cbuffer.data(gf_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_y_xzzz_xxy = cbuffer.data(gf_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_y_xzzz_xxz = cbuffer.data(gf_geom_11_off + 692 * ccomps * dcomps);

            auto g_y_y_xzzz_xyy = cbuffer.data(gf_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_y_xzzz_xyz = cbuffer.data(gf_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_y_xzzz_xzz = cbuffer.data(gf_geom_11_off + 695 * ccomps * dcomps);

            auto g_y_y_xzzz_yyy = cbuffer.data(gf_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_y_xzzz_yyz = cbuffer.data(gf_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_y_xzzz_yzz = cbuffer.data(gf_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_y_xzzz_zzz = cbuffer.data(gf_geom_11_off + 699 * ccomps * dcomps);

            auto g_y_y_yyyy_xxx = cbuffer.data(gf_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_y_yyyy_xxy = cbuffer.data(gf_geom_11_off + 701 * ccomps * dcomps);

            auto g_y_y_yyyy_xxz = cbuffer.data(gf_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_y_yyyy_xyy = cbuffer.data(gf_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_y_yyyy_xyz = cbuffer.data(gf_geom_11_off + 704 * ccomps * dcomps);

            auto g_y_y_yyyy_xzz = cbuffer.data(gf_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_y_yyyy_yyy = cbuffer.data(gf_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_y_yyyy_yyz = cbuffer.data(gf_geom_11_off + 707 * ccomps * dcomps);

            auto g_y_y_yyyy_yzz = cbuffer.data(gf_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_y_yyyy_zzz = cbuffer.data(gf_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_y_yyyz_xxx = cbuffer.data(gf_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_y_yyyz_xxy = cbuffer.data(gf_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_y_yyyz_xxz = cbuffer.data(gf_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_y_yyyz_xyy = cbuffer.data(gf_geom_11_off + 713 * ccomps * dcomps);

            auto g_y_y_yyyz_xyz = cbuffer.data(gf_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_y_yyyz_xzz = cbuffer.data(gf_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_y_yyyz_yyy = cbuffer.data(gf_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_y_yyyz_yyz = cbuffer.data(gf_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_y_yyyz_yzz = cbuffer.data(gf_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_y_yyyz_zzz = cbuffer.data(gf_geom_11_off + 719 * ccomps * dcomps);

            auto g_y_y_yyzz_xxx = cbuffer.data(gf_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_y_yyzz_xxy = cbuffer.data(gf_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_y_yyzz_xxz = cbuffer.data(gf_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_y_yyzz_xyy = cbuffer.data(gf_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_y_yyzz_xyz = cbuffer.data(gf_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_y_yyzz_xzz = cbuffer.data(gf_geom_11_off + 725 * ccomps * dcomps);

            auto g_y_y_yyzz_yyy = cbuffer.data(gf_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_y_yyzz_yyz = cbuffer.data(gf_geom_11_off + 727 * ccomps * dcomps);

            auto g_y_y_yyzz_yzz = cbuffer.data(gf_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_y_yyzz_zzz = cbuffer.data(gf_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_y_yzzz_xxx = cbuffer.data(gf_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_y_yzzz_xxy = cbuffer.data(gf_geom_11_off + 731 * ccomps * dcomps);

            auto g_y_y_yzzz_xxz = cbuffer.data(gf_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_y_yzzz_xyy = cbuffer.data(gf_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_y_yzzz_xyz = cbuffer.data(gf_geom_11_off + 734 * ccomps * dcomps);

            auto g_y_y_yzzz_xzz = cbuffer.data(gf_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_y_yzzz_yyy = cbuffer.data(gf_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_y_yzzz_yyz = cbuffer.data(gf_geom_11_off + 737 * ccomps * dcomps);

            auto g_y_y_yzzz_yzz = cbuffer.data(gf_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_y_yzzz_zzz = cbuffer.data(gf_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_y_zzzz_xxx = cbuffer.data(gf_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_y_zzzz_xxy = cbuffer.data(gf_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_y_zzzz_xxz = cbuffer.data(gf_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_y_zzzz_xyy = cbuffer.data(gf_geom_11_off + 743 * ccomps * dcomps);

            auto g_y_y_zzzz_xyz = cbuffer.data(gf_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_y_zzzz_xzz = cbuffer.data(gf_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_y_zzzz_yyy = cbuffer.data(gf_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_y_zzzz_yyz = cbuffer.data(gf_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_y_zzzz_yzz = cbuffer.data(gf_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_y_zzzz_zzz = cbuffer.data(gf_geom_11_off + 749 * ccomps * dcomps);

            auto g_y_z_xxxx_xxx = cbuffer.data(gf_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_z_xxxx_xxy = cbuffer.data(gf_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_z_xxxx_xxz = cbuffer.data(gf_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_z_xxxx_xyy = cbuffer.data(gf_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_z_xxxx_xyz = cbuffer.data(gf_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_z_xxxx_xzz = cbuffer.data(gf_geom_11_off + 755 * ccomps * dcomps);

            auto g_y_z_xxxx_yyy = cbuffer.data(gf_geom_11_off + 756 * ccomps * dcomps);

            auto g_y_z_xxxx_yyz = cbuffer.data(gf_geom_11_off + 757 * ccomps * dcomps);

            auto g_y_z_xxxx_yzz = cbuffer.data(gf_geom_11_off + 758 * ccomps * dcomps);

            auto g_y_z_xxxx_zzz = cbuffer.data(gf_geom_11_off + 759 * ccomps * dcomps);

            auto g_y_z_xxxy_xxx = cbuffer.data(gf_geom_11_off + 760 * ccomps * dcomps);

            auto g_y_z_xxxy_xxy = cbuffer.data(gf_geom_11_off + 761 * ccomps * dcomps);

            auto g_y_z_xxxy_xxz = cbuffer.data(gf_geom_11_off + 762 * ccomps * dcomps);

            auto g_y_z_xxxy_xyy = cbuffer.data(gf_geom_11_off + 763 * ccomps * dcomps);

            auto g_y_z_xxxy_xyz = cbuffer.data(gf_geom_11_off + 764 * ccomps * dcomps);

            auto g_y_z_xxxy_xzz = cbuffer.data(gf_geom_11_off + 765 * ccomps * dcomps);

            auto g_y_z_xxxy_yyy = cbuffer.data(gf_geom_11_off + 766 * ccomps * dcomps);

            auto g_y_z_xxxy_yyz = cbuffer.data(gf_geom_11_off + 767 * ccomps * dcomps);

            auto g_y_z_xxxy_yzz = cbuffer.data(gf_geom_11_off + 768 * ccomps * dcomps);

            auto g_y_z_xxxy_zzz = cbuffer.data(gf_geom_11_off + 769 * ccomps * dcomps);

            auto g_y_z_xxxz_xxx = cbuffer.data(gf_geom_11_off + 770 * ccomps * dcomps);

            auto g_y_z_xxxz_xxy = cbuffer.data(gf_geom_11_off + 771 * ccomps * dcomps);

            auto g_y_z_xxxz_xxz = cbuffer.data(gf_geom_11_off + 772 * ccomps * dcomps);

            auto g_y_z_xxxz_xyy = cbuffer.data(gf_geom_11_off + 773 * ccomps * dcomps);

            auto g_y_z_xxxz_xyz = cbuffer.data(gf_geom_11_off + 774 * ccomps * dcomps);

            auto g_y_z_xxxz_xzz = cbuffer.data(gf_geom_11_off + 775 * ccomps * dcomps);

            auto g_y_z_xxxz_yyy = cbuffer.data(gf_geom_11_off + 776 * ccomps * dcomps);

            auto g_y_z_xxxz_yyz = cbuffer.data(gf_geom_11_off + 777 * ccomps * dcomps);

            auto g_y_z_xxxz_yzz = cbuffer.data(gf_geom_11_off + 778 * ccomps * dcomps);

            auto g_y_z_xxxz_zzz = cbuffer.data(gf_geom_11_off + 779 * ccomps * dcomps);

            auto g_y_z_xxyy_xxx = cbuffer.data(gf_geom_11_off + 780 * ccomps * dcomps);

            auto g_y_z_xxyy_xxy = cbuffer.data(gf_geom_11_off + 781 * ccomps * dcomps);

            auto g_y_z_xxyy_xxz = cbuffer.data(gf_geom_11_off + 782 * ccomps * dcomps);

            auto g_y_z_xxyy_xyy = cbuffer.data(gf_geom_11_off + 783 * ccomps * dcomps);

            auto g_y_z_xxyy_xyz = cbuffer.data(gf_geom_11_off + 784 * ccomps * dcomps);

            auto g_y_z_xxyy_xzz = cbuffer.data(gf_geom_11_off + 785 * ccomps * dcomps);

            auto g_y_z_xxyy_yyy = cbuffer.data(gf_geom_11_off + 786 * ccomps * dcomps);

            auto g_y_z_xxyy_yyz = cbuffer.data(gf_geom_11_off + 787 * ccomps * dcomps);

            auto g_y_z_xxyy_yzz = cbuffer.data(gf_geom_11_off + 788 * ccomps * dcomps);

            auto g_y_z_xxyy_zzz = cbuffer.data(gf_geom_11_off + 789 * ccomps * dcomps);

            auto g_y_z_xxyz_xxx = cbuffer.data(gf_geom_11_off + 790 * ccomps * dcomps);

            auto g_y_z_xxyz_xxy = cbuffer.data(gf_geom_11_off + 791 * ccomps * dcomps);

            auto g_y_z_xxyz_xxz = cbuffer.data(gf_geom_11_off + 792 * ccomps * dcomps);

            auto g_y_z_xxyz_xyy = cbuffer.data(gf_geom_11_off + 793 * ccomps * dcomps);

            auto g_y_z_xxyz_xyz = cbuffer.data(gf_geom_11_off + 794 * ccomps * dcomps);

            auto g_y_z_xxyz_xzz = cbuffer.data(gf_geom_11_off + 795 * ccomps * dcomps);

            auto g_y_z_xxyz_yyy = cbuffer.data(gf_geom_11_off + 796 * ccomps * dcomps);

            auto g_y_z_xxyz_yyz = cbuffer.data(gf_geom_11_off + 797 * ccomps * dcomps);

            auto g_y_z_xxyz_yzz = cbuffer.data(gf_geom_11_off + 798 * ccomps * dcomps);

            auto g_y_z_xxyz_zzz = cbuffer.data(gf_geom_11_off + 799 * ccomps * dcomps);

            auto g_y_z_xxzz_xxx = cbuffer.data(gf_geom_11_off + 800 * ccomps * dcomps);

            auto g_y_z_xxzz_xxy = cbuffer.data(gf_geom_11_off + 801 * ccomps * dcomps);

            auto g_y_z_xxzz_xxz = cbuffer.data(gf_geom_11_off + 802 * ccomps * dcomps);

            auto g_y_z_xxzz_xyy = cbuffer.data(gf_geom_11_off + 803 * ccomps * dcomps);

            auto g_y_z_xxzz_xyz = cbuffer.data(gf_geom_11_off + 804 * ccomps * dcomps);

            auto g_y_z_xxzz_xzz = cbuffer.data(gf_geom_11_off + 805 * ccomps * dcomps);

            auto g_y_z_xxzz_yyy = cbuffer.data(gf_geom_11_off + 806 * ccomps * dcomps);

            auto g_y_z_xxzz_yyz = cbuffer.data(gf_geom_11_off + 807 * ccomps * dcomps);

            auto g_y_z_xxzz_yzz = cbuffer.data(gf_geom_11_off + 808 * ccomps * dcomps);

            auto g_y_z_xxzz_zzz = cbuffer.data(gf_geom_11_off + 809 * ccomps * dcomps);

            auto g_y_z_xyyy_xxx = cbuffer.data(gf_geom_11_off + 810 * ccomps * dcomps);

            auto g_y_z_xyyy_xxy = cbuffer.data(gf_geom_11_off + 811 * ccomps * dcomps);

            auto g_y_z_xyyy_xxz = cbuffer.data(gf_geom_11_off + 812 * ccomps * dcomps);

            auto g_y_z_xyyy_xyy = cbuffer.data(gf_geom_11_off + 813 * ccomps * dcomps);

            auto g_y_z_xyyy_xyz = cbuffer.data(gf_geom_11_off + 814 * ccomps * dcomps);

            auto g_y_z_xyyy_xzz = cbuffer.data(gf_geom_11_off + 815 * ccomps * dcomps);

            auto g_y_z_xyyy_yyy = cbuffer.data(gf_geom_11_off + 816 * ccomps * dcomps);

            auto g_y_z_xyyy_yyz = cbuffer.data(gf_geom_11_off + 817 * ccomps * dcomps);

            auto g_y_z_xyyy_yzz = cbuffer.data(gf_geom_11_off + 818 * ccomps * dcomps);

            auto g_y_z_xyyy_zzz = cbuffer.data(gf_geom_11_off + 819 * ccomps * dcomps);

            auto g_y_z_xyyz_xxx = cbuffer.data(gf_geom_11_off + 820 * ccomps * dcomps);

            auto g_y_z_xyyz_xxy = cbuffer.data(gf_geom_11_off + 821 * ccomps * dcomps);

            auto g_y_z_xyyz_xxz = cbuffer.data(gf_geom_11_off + 822 * ccomps * dcomps);

            auto g_y_z_xyyz_xyy = cbuffer.data(gf_geom_11_off + 823 * ccomps * dcomps);

            auto g_y_z_xyyz_xyz = cbuffer.data(gf_geom_11_off + 824 * ccomps * dcomps);

            auto g_y_z_xyyz_xzz = cbuffer.data(gf_geom_11_off + 825 * ccomps * dcomps);

            auto g_y_z_xyyz_yyy = cbuffer.data(gf_geom_11_off + 826 * ccomps * dcomps);

            auto g_y_z_xyyz_yyz = cbuffer.data(gf_geom_11_off + 827 * ccomps * dcomps);

            auto g_y_z_xyyz_yzz = cbuffer.data(gf_geom_11_off + 828 * ccomps * dcomps);

            auto g_y_z_xyyz_zzz = cbuffer.data(gf_geom_11_off + 829 * ccomps * dcomps);

            auto g_y_z_xyzz_xxx = cbuffer.data(gf_geom_11_off + 830 * ccomps * dcomps);

            auto g_y_z_xyzz_xxy = cbuffer.data(gf_geom_11_off + 831 * ccomps * dcomps);

            auto g_y_z_xyzz_xxz = cbuffer.data(gf_geom_11_off + 832 * ccomps * dcomps);

            auto g_y_z_xyzz_xyy = cbuffer.data(gf_geom_11_off + 833 * ccomps * dcomps);

            auto g_y_z_xyzz_xyz = cbuffer.data(gf_geom_11_off + 834 * ccomps * dcomps);

            auto g_y_z_xyzz_xzz = cbuffer.data(gf_geom_11_off + 835 * ccomps * dcomps);

            auto g_y_z_xyzz_yyy = cbuffer.data(gf_geom_11_off + 836 * ccomps * dcomps);

            auto g_y_z_xyzz_yyz = cbuffer.data(gf_geom_11_off + 837 * ccomps * dcomps);

            auto g_y_z_xyzz_yzz = cbuffer.data(gf_geom_11_off + 838 * ccomps * dcomps);

            auto g_y_z_xyzz_zzz = cbuffer.data(gf_geom_11_off + 839 * ccomps * dcomps);

            auto g_y_z_xzzz_xxx = cbuffer.data(gf_geom_11_off + 840 * ccomps * dcomps);

            auto g_y_z_xzzz_xxy = cbuffer.data(gf_geom_11_off + 841 * ccomps * dcomps);

            auto g_y_z_xzzz_xxz = cbuffer.data(gf_geom_11_off + 842 * ccomps * dcomps);

            auto g_y_z_xzzz_xyy = cbuffer.data(gf_geom_11_off + 843 * ccomps * dcomps);

            auto g_y_z_xzzz_xyz = cbuffer.data(gf_geom_11_off + 844 * ccomps * dcomps);

            auto g_y_z_xzzz_xzz = cbuffer.data(gf_geom_11_off + 845 * ccomps * dcomps);

            auto g_y_z_xzzz_yyy = cbuffer.data(gf_geom_11_off + 846 * ccomps * dcomps);

            auto g_y_z_xzzz_yyz = cbuffer.data(gf_geom_11_off + 847 * ccomps * dcomps);

            auto g_y_z_xzzz_yzz = cbuffer.data(gf_geom_11_off + 848 * ccomps * dcomps);

            auto g_y_z_xzzz_zzz = cbuffer.data(gf_geom_11_off + 849 * ccomps * dcomps);

            auto g_y_z_yyyy_xxx = cbuffer.data(gf_geom_11_off + 850 * ccomps * dcomps);

            auto g_y_z_yyyy_xxy = cbuffer.data(gf_geom_11_off + 851 * ccomps * dcomps);

            auto g_y_z_yyyy_xxz = cbuffer.data(gf_geom_11_off + 852 * ccomps * dcomps);

            auto g_y_z_yyyy_xyy = cbuffer.data(gf_geom_11_off + 853 * ccomps * dcomps);

            auto g_y_z_yyyy_xyz = cbuffer.data(gf_geom_11_off + 854 * ccomps * dcomps);

            auto g_y_z_yyyy_xzz = cbuffer.data(gf_geom_11_off + 855 * ccomps * dcomps);

            auto g_y_z_yyyy_yyy = cbuffer.data(gf_geom_11_off + 856 * ccomps * dcomps);

            auto g_y_z_yyyy_yyz = cbuffer.data(gf_geom_11_off + 857 * ccomps * dcomps);

            auto g_y_z_yyyy_yzz = cbuffer.data(gf_geom_11_off + 858 * ccomps * dcomps);

            auto g_y_z_yyyy_zzz = cbuffer.data(gf_geom_11_off + 859 * ccomps * dcomps);

            auto g_y_z_yyyz_xxx = cbuffer.data(gf_geom_11_off + 860 * ccomps * dcomps);

            auto g_y_z_yyyz_xxy = cbuffer.data(gf_geom_11_off + 861 * ccomps * dcomps);

            auto g_y_z_yyyz_xxz = cbuffer.data(gf_geom_11_off + 862 * ccomps * dcomps);

            auto g_y_z_yyyz_xyy = cbuffer.data(gf_geom_11_off + 863 * ccomps * dcomps);

            auto g_y_z_yyyz_xyz = cbuffer.data(gf_geom_11_off + 864 * ccomps * dcomps);

            auto g_y_z_yyyz_xzz = cbuffer.data(gf_geom_11_off + 865 * ccomps * dcomps);

            auto g_y_z_yyyz_yyy = cbuffer.data(gf_geom_11_off + 866 * ccomps * dcomps);

            auto g_y_z_yyyz_yyz = cbuffer.data(gf_geom_11_off + 867 * ccomps * dcomps);

            auto g_y_z_yyyz_yzz = cbuffer.data(gf_geom_11_off + 868 * ccomps * dcomps);

            auto g_y_z_yyyz_zzz = cbuffer.data(gf_geom_11_off + 869 * ccomps * dcomps);

            auto g_y_z_yyzz_xxx = cbuffer.data(gf_geom_11_off + 870 * ccomps * dcomps);

            auto g_y_z_yyzz_xxy = cbuffer.data(gf_geom_11_off + 871 * ccomps * dcomps);

            auto g_y_z_yyzz_xxz = cbuffer.data(gf_geom_11_off + 872 * ccomps * dcomps);

            auto g_y_z_yyzz_xyy = cbuffer.data(gf_geom_11_off + 873 * ccomps * dcomps);

            auto g_y_z_yyzz_xyz = cbuffer.data(gf_geom_11_off + 874 * ccomps * dcomps);

            auto g_y_z_yyzz_xzz = cbuffer.data(gf_geom_11_off + 875 * ccomps * dcomps);

            auto g_y_z_yyzz_yyy = cbuffer.data(gf_geom_11_off + 876 * ccomps * dcomps);

            auto g_y_z_yyzz_yyz = cbuffer.data(gf_geom_11_off + 877 * ccomps * dcomps);

            auto g_y_z_yyzz_yzz = cbuffer.data(gf_geom_11_off + 878 * ccomps * dcomps);

            auto g_y_z_yyzz_zzz = cbuffer.data(gf_geom_11_off + 879 * ccomps * dcomps);

            auto g_y_z_yzzz_xxx = cbuffer.data(gf_geom_11_off + 880 * ccomps * dcomps);

            auto g_y_z_yzzz_xxy = cbuffer.data(gf_geom_11_off + 881 * ccomps * dcomps);

            auto g_y_z_yzzz_xxz = cbuffer.data(gf_geom_11_off + 882 * ccomps * dcomps);

            auto g_y_z_yzzz_xyy = cbuffer.data(gf_geom_11_off + 883 * ccomps * dcomps);

            auto g_y_z_yzzz_xyz = cbuffer.data(gf_geom_11_off + 884 * ccomps * dcomps);

            auto g_y_z_yzzz_xzz = cbuffer.data(gf_geom_11_off + 885 * ccomps * dcomps);

            auto g_y_z_yzzz_yyy = cbuffer.data(gf_geom_11_off + 886 * ccomps * dcomps);

            auto g_y_z_yzzz_yyz = cbuffer.data(gf_geom_11_off + 887 * ccomps * dcomps);

            auto g_y_z_yzzz_yzz = cbuffer.data(gf_geom_11_off + 888 * ccomps * dcomps);

            auto g_y_z_yzzz_zzz = cbuffer.data(gf_geom_11_off + 889 * ccomps * dcomps);

            auto g_y_z_zzzz_xxx = cbuffer.data(gf_geom_11_off + 890 * ccomps * dcomps);

            auto g_y_z_zzzz_xxy = cbuffer.data(gf_geom_11_off + 891 * ccomps * dcomps);

            auto g_y_z_zzzz_xxz = cbuffer.data(gf_geom_11_off + 892 * ccomps * dcomps);

            auto g_y_z_zzzz_xyy = cbuffer.data(gf_geom_11_off + 893 * ccomps * dcomps);

            auto g_y_z_zzzz_xyz = cbuffer.data(gf_geom_11_off + 894 * ccomps * dcomps);

            auto g_y_z_zzzz_xzz = cbuffer.data(gf_geom_11_off + 895 * ccomps * dcomps);

            auto g_y_z_zzzz_yyy = cbuffer.data(gf_geom_11_off + 896 * ccomps * dcomps);

            auto g_y_z_zzzz_yyz = cbuffer.data(gf_geom_11_off + 897 * ccomps * dcomps);

            auto g_y_z_zzzz_yzz = cbuffer.data(gf_geom_11_off + 898 * ccomps * dcomps);

            auto g_y_z_zzzz_zzz = cbuffer.data(gf_geom_11_off + 899 * ccomps * dcomps);

            auto g_z_x_xxxx_xxx = cbuffer.data(gf_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_x_xxxx_xxy = cbuffer.data(gf_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_x_xxxx_xxz = cbuffer.data(gf_geom_11_off + 902 * ccomps * dcomps);

            auto g_z_x_xxxx_xyy = cbuffer.data(gf_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_x_xxxx_xyz = cbuffer.data(gf_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_x_xxxx_xzz = cbuffer.data(gf_geom_11_off + 905 * ccomps * dcomps);

            auto g_z_x_xxxx_yyy = cbuffer.data(gf_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_x_xxxx_yyz = cbuffer.data(gf_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_x_xxxx_yzz = cbuffer.data(gf_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_x_xxxx_zzz = cbuffer.data(gf_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_x_xxxy_xxx = cbuffer.data(gf_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_x_xxxy_xxy = cbuffer.data(gf_geom_11_off + 911 * ccomps * dcomps);

            auto g_z_x_xxxy_xxz = cbuffer.data(gf_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_x_xxxy_xyy = cbuffer.data(gf_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_x_xxxy_xyz = cbuffer.data(gf_geom_11_off + 914 * ccomps * dcomps);

            auto g_z_x_xxxy_xzz = cbuffer.data(gf_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_x_xxxy_yyy = cbuffer.data(gf_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_x_xxxy_yyz = cbuffer.data(gf_geom_11_off + 917 * ccomps * dcomps);

            auto g_z_x_xxxy_yzz = cbuffer.data(gf_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_x_xxxy_zzz = cbuffer.data(gf_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_x_xxxz_xxx = cbuffer.data(gf_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_x_xxxz_xxy = cbuffer.data(gf_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_x_xxxz_xxz = cbuffer.data(gf_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_x_xxxz_xyy = cbuffer.data(gf_geom_11_off + 923 * ccomps * dcomps);

            auto g_z_x_xxxz_xyz = cbuffer.data(gf_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_x_xxxz_xzz = cbuffer.data(gf_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_x_xxxz_yyy = cbuffer.data(gf_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_x_xxxz_yyz = cbuffer.data(gf_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_x_xxxz_yzz = cbuffer.data(gf_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_x_xxxz_zzz = cbuffer.data(gf_geom_11_off + 929 * ccomps * dcomps);

            auto g_z_x_xxyy_xxx = cbuffer.data(gf_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_x_xxyy_xxy = cbuffer.data(gf_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_x_xxyy_xxz = cbuffer.data(gf_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_x_xxyy_xyy = cbuffer.data(gf_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_x_xxyy_xyz = cbuffer.data(gf_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_x_xxyy_xzz = cbuffer.data(gf_geom_11_off + 935 * ccomps * dcomps);

            auto g_z_x_xxyy_yyy = cbuffer.data(gf_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_x_xxyy_yyz = cbuffer.data(gf_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_x_xxyy_yzz = cbuffer.data(gf_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_x_xxyy_zzz = cbuffer.data(gf_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_x_xxyz_xxx = cbuffer.data(gf_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_x_xxyz_xxy = cbuffer.data(gf_geom_11_off + 941 * ccomps * dcomps);

            auto g_z_x_xxyz_xxz = cbuffer.data(gf_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_x_xxyz_xyy = cbuffer.data(gf_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_x_xxyz_xyz = cbuffer.data(gf_geom_11_off + 944 * ccomps * dcomps);

            auto g_z_x_xxyz_xzz = cbuffer.data(gf_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_x_xxyz_yyy = cbuffer.data(gf_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_x_xxyz_yyz = cbuffer.data(gf_geom_11_off + 947 * ccomps * dcomps);

            auto g_z_x_xxyz_yzz = cbuffer.data(gf_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_x_xxyz_zzz = cbuffer.data(gf_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_x_xxzz_xxx = cbuffer.data(gf_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_x_xxzz_xxy = cbuffer.data(gf_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_x_xxzz_xxz = cbuffer.data(gf_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_x_xxzz_xyy = cbuffer.data(gf_geom_11_off + 953 * ccomps * dcomps);

            auto g_z_x_xxzz_xyz = cbuffer.data(gf_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_x_xxzz_xzz = cbuffer.data(gf_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_x_xxzz_yyy = cbuffer.data(gf_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_x_xxzz_yyz = cbuffer.data(gf_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_x_xxzz_yzz = cbuffer.data(gf_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_x_xxzz_zzz = cbuffer.data(gf_geom_11_off + 959 * ccomps * dcomps);

            auto g_z_x_xyyy_xxx = cbuffer.data(gf_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_x_xyyy_xxy = cbuffer.data(gf_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_x_xyyy_xxz = cbuffer.data(gf_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_x_xyyy_xyy = cbuffer.data(gf_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_x_xyyy_xyz = cbuffer.data(gf_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_x_xyyy_xzz = cbuffer.data(gf_geom_11_off + 965 * ccomps * dcomps);

            auto g_z_x_xyyy_yyy = cbuffer.data(gf_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_x_xyyy_yyz = cbuffer.data(gf_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_x_xyyy_yzz = cbuffer.data(gf_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_x_xyyy_zzz = cbuffer.data(gf_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_x_xyyz_xxx = cbuffer.data(gf_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_x_xyyz_xxy = cbuffer.data(gf_geom_11_off + 971 * ccomps * dcomps);

            auto g_z_x_xyyz_xxz = cbuffer.data(gf_geom_11_off + 972 * ccomps * dcomps);

            auto g_z_x_xyyz_xyy = cbuffer.data(gf_geom_11_off + 973 * ccomps * dcomps);

            auto g_z_x_xyyz_xyz = cbuffer.data(gf_geom_11_off + 974 * ccomps * dcomps);

            auto g_z_x_xyyz_xzz = cbuffer.data(gf_geom_11_off + 975 * ccomps * dcomps);

            auto g_z_x_xyyz_yyy = cbuffer.data(gf_geom_11_off + 976 * ccomps * dcomps);

            auto g_z_x_xyyz_yyz = cbuffer.data(gf_geom_11_off + 977 * ccomps * dcomps);

            auto g_z_x_xyyz_yzz = cbuffer.data(gf_geom_11_off + 978 * ccomps * dcomps);

            auto g_z_x_xyyz_zzz = cbuffer.data(gf_geom_11_off + 979 * ccomps * dcomps);

            auto g_z_x_xyzz_xxx = cbuffer.data(gf_geom_11_off + 980 * ccomps * dcomps);

            auto g_z_x_xyzz_xxy = cbuffer.data(gf_geom_11_off + 981 * ccomps * dcomps);

            auto g_z_x_xyzz_xxz = cbuffer.data(gf_geom_11_off + 982 * ccomps * dcomps);

            auto g_z_x_xyzz_xyy = cbuffer.data(gf_geom_11_off + 983 * ccomps * dcomps);

            auto g_z_x_xyzz_xyz = cbuffer.data(gf_geom_11_off + 984 * ccomps * dcomps);

            auto g_z_x_xyzz_xzz = cbuffer.data(gf_geom_11_off + 985 * ccomps * dcomps);

            auto g_z_x_xyzz_yyy = cbuffer.data(gf_geom_11_off + 986 * ccomps * dcomps);

            auto g_z_x_xyzz_yyz = cbuffer.data(gf_geom_11_off + 987 * ccomps * dcomps);

            auto g_z_x_xyzz_yzz = cbuffer.data(gf_geom_11_off + 988 * ccomps * dcomps);

            auto g_z_x_xyzz_zzz = cbuffer.data(gf_geom_11_off + 989 * ccomps * dcomps);

            auto g_z_x_xzzz_xxx = cbuffer.data(gf_geom_11_off + 990 * ccomps * dcomps);

            auto g_z_x_xzzz_xxy = cbuffer.data(gf_geom_11_off + 991 * ccomps * dcomps);

            auto g_z_x_xzzz_xxz = cbuffer.data(gf_geom_11_off + 992 * ccomps * dcomps);

            auto g_z_x_xzzz_xyy = cbuffer.data(gf_geom_11_off + 993 * ccomps * dcomps);

            auto g_z_x_xzzz_xyz = cbuffer.data(gf_geom_11_off + 994 * ccomps * dcomps);

            auto g_z_x_xzzz_xzz = cbuffer.data(gf_geom_11_off + 995 * ccomps * dcomps);

            auto g_z_x_xzzz_yyy = cbuffer.data(gf_geom_11_off + 996 * ccomps * dcomps);

            auto g_z_x_xzzz_yyz = cbuffer.data(gf_geom_11_off + 997 * ccomps * dcomps);

            auto g_z_x_xzzz_yzz = cbuffer.data(gf_geom_11_off + 998 * ccomps * dcomps);

            auto g_z_x_xzzz_zzz = cbuffer.data(gf_geom_11_off + 999 * ccomps * dcomps);

            auto g_z_x_yyyy_xxx = cbuffer.data(gf_geom_11_off + 1000 * ccomps * dcomps);

            auto g_z_x_yyyy_xxy = cbuffer.data(gf_geom_11_off + 1001 * ccomps * dcomps);

            auto g_z_x_yyyy_xxz = cbuffer.data(gf_geom_11_off + 1002 * ccomps * dcomps);

            auto g_z_x_yyyy_xyy = cbuffer.data(gf_geom_11_off + 1003 * ccomps * dcomps);

            auto g_z_x_yyyy_xyz = cbuffer.data(gf_geom_11_off + 1004 * ccomps * dcomps);

            auto g_z_x_yyyy_xzz = cbuffer.data(gf_geom_11_off + 1005 * ccomps * dcomps);

            auto g_z_x_yyyy_yyy = cbuffer.data(gf_geom_11_off + 1006 * ccomps * dcomps);

            auto g_z_x_yyyy_yyz = cbuffer.data(gf_geom_11_off + 1007 * ccomps * dcomps);

            auto g_z_x_yyyy_yzz = cbuffer.data(gf_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_x_yyyy_zzz = cbuffer.data(gf_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_x_yyyz_xxx = cbuffer.data(gf_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_x_yyyz_xxy = cbuffer.data(gf_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_x_yyyz_xxz = cbuffer.data(gf_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_x_yyyz_xyy = cbuffer.data(gf_geom_11_off + 1013 * ccomps * dcomps);

            auto g_z_x_yyyz_xyz = cbuffer.data(gf_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_x_yyyz_xzz = cbuffer.data(gf_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_x_yyyz_yyy = cbuffer.data(gf_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_x_yyyz_yyz = cbuffer.data(gf_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_x_yyyz_yzz = cbuffer.data(gf_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_x_yyyz_zzz = cbuffer.data(gf_geom_11_off + 1019 * ccomps * dcomps);

            auto g_z_x_yyzz_xxx = cbuffer.data(gf_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_x_yyzz_xxy = cbuffer.data(gf_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_x_yyzz_xxz = cbuffer.data(gf_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_x_yyzz_xyy = cbuffer.data(gf_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_x_yyzz_xyz = cbuffer.data(gf_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_x_yyzz_xzz = cbuffer.data(gf_geom_11_off + 1025 * ccomps * dcomps);

            auto g_z_x_yyzz_yyy = cbuffer.data(gf_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_x_yyzz_yyz = cbuffer.data(gf_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_x_yyzz_yzz = cbuffer.data(gf_geom_11_off + 1028 * ccomps * dcomps);

            auto g_z_x_yyzz_zzz = cbuffer.data(gf_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_x_yzzz_xxx = cbuffer.data(gf_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_x_yzzz_xxy = cbuffer.data(gf_geom_11_off + 1031 * ccomps * dcomps);

            auto g_z_x_yzzz_xxz = cbuffer.data(gf_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_x_yzzz_xyy = cbuffer.data(gf_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_x_yzzz_xyz = cbuffer.data(gf_geom_11_off + 1034 * ccomps * dcomps);

            auto g_z_x_yzzz_xzz = cbuffer.data(gf_geom_11_off + 1035 * ccomps * dcomps);

            auto g_z_x_yzzz_yyy = cbuffer.data(gf_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_x_yzzz_yyz = cbuffer.data(gf_geom_11_off + 1037 * ccomps * dcomps);

            auto g_z_x_yzzz_yzz = cbuffer.data(gf_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_x_yzzz_zzz = cbuffer.data(gf_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_x_zzzz_xxx = cbuffer.data(gf_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_x_zzzz_xxy = cbuffer.data(gf_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_x_zzzz_xxz = cbuffer.data(gf_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_x_zzzz_xyy = cbuffer.data(gf_geom_11_off + 1043 * ccomps * dcomps);

            auto g_z_x_zzzz_xyz = cbuffer.data(gf_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_x_zzzz_xzz = cbuffer.data(gf_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_x_zzzz_yyy = cbuffer.data(gf_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_x_zzzz_yyz = cbuffer.data(gf_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_x_zzzz_yzz = cbuffer.data(gf_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_x_zzzz_zzz = cbuffer.data(gf_geom_11_off + 1049 * ccomps * dcomps);

            auto g_z_y_xxxx_xxx = cbuffer.data(gf_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_y_xxxx_xxy = cbuffer.data(gf_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_y_xxxx_xxz = cbuffer.data(gf_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_y_xxxx_xyy = cbuffer.data(gf_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_y_xxxx_xyz = cbuffer.data(gf_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_y_xxxx_xzz = cbuffer.data(gf_geom_11_off + 1055 * ccomps * dcomps);

            auto g_z_y_xxxx_yyy = cbuffer.data(gf_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_y_xxxx_yyz = cbuffer.data(gf_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_y_xxxx_yzz = cbuffer.data(gf_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_y_xxxx_zzz = cbuffer.data(gf_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_y_xxxy_xxx = cbuffer.data(gf_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_y_xxxy_xxy = cbuffer.data(gf_geom_11_off + 1061 * ccomps * dcomps);

            auto g_z_y_xxxy_xxz = cbuffer.data(gf_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_y_xxxy_xyy = cbuffer.data(gf_geom_11_off + 1063 * ccomps * dcomps);

            auto g_z_y_xxxy_xyz = cbuffer.data(gf_geom_11_off + 1064 * ccomps * dcomps);

            auto g_z_y_xxxy_xzz = cbuffer.data(gf_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_y_xxxy_yyy = cbuffer.data(gf_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_y_xxxy_yyz = cbuffer.data(gf_geom_11_off + 1067 * ccomps * dcomps);

            auto g_z_y_xxxy_yzz = cbuffer.data(gf_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_y_xxxy_zzz = cbuffer.data(gf_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_y_xxxz_xxx = cbuffer.data(gf_geom_11_off + 1070 * ccomps * dcomps);

            auto g_z_y_xxxz_xxy = cbuffer.data(gf_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_y_xxxz_xxz = cbuffer.data(gf_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_y_xxxz_xyy = cbuffer.data(gf_geom_11_off + 1073 * ccomps * dcomps);

            auto g_z_y_xxxz_xyz = cbuffer.data(gf_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_y_xxxz_xzz = cbuffer.data(gf_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_y_xxxz_yyy = cbuffer.data(gf_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_y_xxxz_yyz = cbuffer.data(gf_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_y_xxxz_yzz = cbuffer.data(gf_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_y_xxxz_zzz = cbuffer.data(gf_geom_11_off + 1079 * ccomps * dcomps);

            auto g_z_y_xxyy_xxx = cbuffer.data(gf_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_y_xxyy_xxy = cbuffer.data(gf_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_y_xxyy_xxz = cbuffer.data(gf_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_y_xxyy_xyy = cbuffer.data(gf_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_y_xxyy_xyz = cbuffer.data(gf_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_y_xxyy_xzz = cbuffer.data(gf_geom_11_off + 1085 * ccomps * dcomps);

            auto g_z_y_xxyy_yyy = cbuffer.data(gf_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_y_xxyy_yyz = cbuffer.data(gf_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_y_xxyy_yzz = cbuffer.data(gf_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_y_xxyy_zzz = cbuffer.data(gf_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_y_xxyz_xxx = cbuffer.data(gf_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_y_xxyz_xxy = cbuffer.data(gf_geom_11_off + 1091 * ccomps * dcomps);

            auto g_z_y_xxyz_xxz = cbuffer.data(gf_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_y_xxyz_xyy = cbuffer.data(gf_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_y_xxyz_xyz = cbuffer.data(gf_geom_11_off + 1094 * ccomps * dcomps);

            auto g_z_y_xxyz_xzz = cbuffer.data(gf_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_y_xxyz_yyy = cbuffer.data(gf_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_y_xxyz_yyz = cbuffer.data(gf_geom_11_off + 1097 * ccomps * dcomps);

            auto g_z_y_xxyz_yzz = cbuffer.data(gf_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_y_xxyz_zzz = cbuffer.data(gf_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_y_xxzz_xxx = cbuffer.data(gf_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_y_xxzz_xxy = cbuffer.data(gf_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_y_xxzz_xxz = cbuffer.data(gf_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_y_xxzz_xyy = cbuffer.data(gf_geom_11_off + 1103 * ccomps * dcomps);

            auto g_z_y_xxzz_xyz = cbuffer.data(gf_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_y_xxzz_xzz = cbuffer.data(gf_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_y_xxzz_yyy = cbuffer.data(gf_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_y_xxzz_yyz = cbuffer.data(gf_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_y_xxzz_yzz = cbuffer.data(gf_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_y_xxzz_zzz = cbuffer.data(gf_geom_11_off + 1109 * ccomps * dcomps);

            auto g_z_y_xyyy_xxx = cbuffer.data(gf_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_y_xyyy_xxy = cbuffer.data(gf_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_y_xyyy_xxz = cbuffer.data(gf_geom_11_off + 1112 * ccomps * dcomps);

            auto g_z_y_xyyy_xyy = cbuffer.data(gf_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_y_xyyy_xyz = cbuffer.data(gf_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_y_xyyy_xzz = cbuffer.data(gf_geom_11_off + 1115 * ccomps * dcomps);

            auto g_z_y_xyyy_yyy = cbuffer.data(gf_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_y_xyyy_yyz = cbuffer.data(gf_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_y_xyyy_yzz = cbuffer.data(gf_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_y_xyyy_zzz = cbuffer.data(gf_geom_11_off + 1119 * ccomps * dcomps);

            auto g_z_y_xyyz_xxx = cbuffer.data(gf_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_y_xyyz_xxy = cbuffer.data(gf_geom_11_off + 1121 * ccomps * dcomps);

            auto g_z_y_xyyz_xxz = cbuffer.data(gf_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_y_xyyz_xyy = cbuffer.data(gf_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_y_xyyz_xyz = cbuffer.data(gf_geom_11_off + 1124 * ccomps * dcomps);

            auto g_z_y_xyyz_xzz = cbuffer.data(gf_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_y_xyyz_yyy = cbuffer.data(gf_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_y_xyyz_yyz = cbuffer.data(gf_geom_11_off + 1127 * ccomps * dcomps);

            auto g_z_y_xyyz_yzz = cbuffer.data(gf_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_y_xyyz_zzz = cbuffer.data(gf_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_y_xyzz_xxx = cbuffer.data(gf_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_y_xyzz_xxy = cbuffer.data(gf_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_y_xyzz_xxz = cbuffer.data(gf_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_y_xyzz_xyy = cbuffer.data(gf_geom_11_off + 1133 * ccomps * dcomps);

            auto g_z_y_xyzz_xyz = cbuffer.data(gf_geom_11_off + 1134 * ccomps * dcomps);

            auto g_z_y_xyzz_xzz = cbuffer.data(gf_geom_11_off + 1135 * ccomps * dcomps);

            auto g_z_y_xyzz_yyy = cbuffer.data(gf_geom_11_off + 1136 * ccomps * dcomps);

            auto g_z_y_xyzz_yyz = cbuffer.data(gf_geom_11_off + 1137 * ccomps * dcomps);

            auto g_z_y_xyzz_yzz = cbuffer.data(gf_geom_11_off + 1138 * ccomps * dcomps);

            auto g_z_y_xyzz_zzz = cbuffer.data(gf_geom_11_off + 1139 * ccomps * dcomps);

            auto g_z_y_xzzz_xxx = cbuffer.data(gf_geom_11_off + 1140 * ccomps * dcomps);

            auto g_z_y_xzzz_xxy = cbuffer.data(gf_geom_11_off + 1141 * ccomps * dcomps);

            auto g_z_y_xzzz_xxz = cbuffer.data(gf_geom_11_off + 1142 * ccomps * dcomps);

            auto g_z_y_xzzz_xyy = cbuffer.data(gf_geom_11_off + 1143 * ccomps * dcomps);

            auto g_z_y_xzzz_xyz = cbuffer.data(gf_geom_11_off + 1144 * ccomps * dcomps);

            auto g_z_y_xzzz_xzz = cbuffer.data(gf_geom_11_off + 1145 * ccomps * dcomps);

            auto g_z_y_xzzz_yyy = cbuffer.data(gf_geom_11_off + 1146 * ccomps * dcomps);

            auto g_z_y_xzzz_yyz = cbuffer.data(gf_geom_11_off + 1147 * ccomps * dcomps);

            auto g_z_y_xzzz_yzz = cbuffer.data(gf_geom_11_off + 1148 * ccomps * dcomps);

            auto g_z_y_xzzz_zzz = cbuffer.data(gf_geom_11_off + 1149 * ccomps * dcomps);

            auto g_z_y_yyyy_xxx = cbuffer.data(gf_geom_11_off + 1150 * ccomps * dcomps);

            auto g_z_y_yyyy_xxy = cbuffer.data(gf_geom_11_off + 1151 * ccomps * dcomps);

            auto g_z_y_yyyy_xxz = cbuffer.data(gf_geom_11_off + 1152 * ccomps * dcomps);

            auto g_z_y_yyyy_xyy = cbuffer.data(gf_geom_11_off + 1153 * ccomps * dcomps);

            auto g_z_y_yyyy_xyz = cbuffer.data(gf_geom_11_off + 1154 * ccomps * dcomps);

            auto g_z_y_yyyy_xzz = cbuffer.data(gf_geom_11_off + 1155 * ccomps * dcomps);

            auto g_z_y_yyyy_yyy = cbuffer.data(gf_geom_11_off + 1156 * ccomps * dcomps);

            auto g_z_y_yyyy_yyz = cbuffer.data(gf_geom_11_off + 1157 * ccomps * dcomps);

            auto g_z_y_yyyy_yzz = cbuffer.data(gf_geom_11_off + 1158 * ccomps * dcomps);

            auto g_z_y_yyyy_zzz = cbuffer.data(gf_geom_11_off + 1159 * ccomps * dcomps);

            auto g_z_y_yyyz_xxx = cbuffer.data(gf_geom_11_off + 1160 * ccomps * dcomps);

            auto g_z_y_yyyz_xxy = cbuffer.data(gf_geom_11_off + 1161 * ccomps * dcomps);

            auto g_z_y_yyyz_xxz = cbuffer.data(gf_geom_11_off + 1162 * ccomps * dcomps);

            auto g_z_y_yyyz_xyy = cbuffer.data(gf_geom_11_off + 1163 * ccomps * dcomps);

            auto g_z_y_yyyz_xyz = cbuffer.data(gf_geom_11_off + 1164 * ccomps * dcomps);

            auto g_z_y_yyyz_xzz = cbuffer.data(gf_geom_11_off + 1165 * ccomps * dcomps);

            auto g_z_y_yyyz_yyy = cbuffer.data(gf_geom_11_off + 1166 * ccomps * dcomps);

            auto g_z_y_yyyz_yyz = cbuffer.data(gf_geom_11_off + 1167 * ccomps * dcomps);

            auto g_z_y_yyyz_yzz = cbuffer.data(gf_geom_11_off + 1168 * ccomps * dcomps);

            auto g_z_y_yyyz_zzz = cbuffer.data(gf_geom_11_off + 1169 * ccomps * dcomps);

            auto g_z_y_yyzz_xxx = cbuffer.data(gf_geom_11_off + 1170 * ccomps * dcomps);

            auto g_z_y_yyzz_xxy = cbuffer.data(gf_geom_11_off + 1171 * ccomps * dcomps);

            auto g_z_y_yyzz_xxz = cbuffer.data(gf_geom_11_off + 1172 * ccomps * dcomps);

            auto g_z_y_yyzz_xyy = cbuffer.data(gf_geom_11_off + 1173 * ccomps * dcomps);

            auto g_z_y_yyzz_xyz = cbuffer.data(gf_geom_11_off + 1174 * ccomps * dcomps);

            auto g_z_y_yyzz_xzz = cbuffer.data(gf_geom_11_off + 1175 * ccomps * dcomps);

            auto g_z_y_yyzz_yyy = cbuffer.data(gf_geom_11_off + 1176 * ccomps * dcomps);

            auto g_z_y_yyzz_yyz = cbuffer.data(gf_geom_11_off + 1177 * ccomps * dcomps);

            auto g_z_y_yyzz_yzz = cbuffer.data(gf_geom_11_off + 1178 * ccomps * dcomps);

            auto g_z_y_yyzz_zzz = cbuffer.data(gf_geom_11_off + 1179 * ccomps * dcomps);

            auto g_z_y_yzzz_xxx = cbuffer.data(gf_geom_11_off + 1180 * ccomps * dcomps);

            auto g_z_y_yzzz_xxy = cbuffer.data(gf_geom_11_off + 1181 * ccomps * dcomps);

            auto g_z_y_yzzz_xxz = cbuffer.data(gf_geom_11_off + 1182 * ccomps * dcomps);

            auto g_z_y_yzzz_xyy = cbuffer.data(gf_geom_11_off + 1183 * ccomps * dcomps);

            auto g_z_y_yzzz_xyz = cbuffer.data(gf_geom_11_off + 1184 * ccomps * dcomps);

            auto g_z_y_yzzz_xzz = cbuffer.data(gf_geom_11_off + 1185 * ccomps * dcomps);

            auto g_z_y_yzzz_yyy = cbuffer.data(gf_geom_11_off + 1186 * ccomps * dcomps);

            auto g_z_y_yzzz_yyz = cbuffer.data(gf_geom_11_off + 1187 * ccomps * dcomps);

            auto g_z_y_yzzz_yzz = cbuffer.data(gf_geom_11_off + 1188 * ccomps * dcomps);

            auto g_z_y_yzzz_zzz = cbuffer.data(gf_geom_11_off + 1189 * ccomps * dcomps);

            auto g_z_y_zzzz_xxx = cbuffer.data(gf_geom_11_off + 1190 * ccomps * dcomps);

            auto g_z_y_zzzz_xxy = cbuffer.data(gf_geom_11_off + 1191 * ccomps * dcomps);

            auto g_z_y_zzzz_xxz = cbuffer.data(gf_geom_11_off + 1192 * ccomps * dcomps);

            auto g_z_y_zzzz_xyy = cbuffer.data(gf_geom_11_off + 1193 * ccomps * dcomps);

            auto g_z_y_zzzz_xyz = cbuffer.data(gf_geom_11_off + 1194 * ccomps * dcomps);

            auto g_z_y_zzzz_xzz = cbuffer.data(gf_geom_11_off + 1195 * ccomps * dcomps);

            auto g_z_y_zzzz_yyy = cbuffer.data(gf_geom_11_off + 1196 * ccomps * dcomps);

            auto g_z_y_zzzz_yyz = cbuffer.data(gf_geom_11_off + 1197 * ccomps * dcomps);

            auto g_z_y_zzzz_yzz = cbuffer.data(gf_geom_11_off + 1198 * ccomps * dcomps);

            auto g_z_y_zzzz_zzz = cbuffer.data(gf_geom_11_off + 1199 * ccomps * dcomps);

            auto g_z_z_xxxx_xxx = cbuffer.data(gf_geom_11_off + 1200 * ccomps * dcomps);

            auto g_z_z_xxxx_xxy = cbuffer.data(gf_geom_11_off + 1201 * ccomps * dcomps);

            auto g_z_z_xxxx_xxz = cbuffer.data(gf_geom_11_off + 1202 * ccomps * dcomps);

            auto g_z_z_xxxx_xyy = cbuffer.data(gf_geom_11_off + 1203 * ccomps * dcomps);

            auto g_z_z_xxxx_xyz = cbuffer.data(gf_geom_11_off + 1204 * ccomps * dcomps);

            auto g_z_z_xxxx_xzz = cbuffer.data(gf_geom_11_off + 1205 * ccomps * dcomps);

            auto g_z_z_xxxx_yyy = cbuffer.data(gf_geom_11_off + 1206 * ccomps * dcomps);

            auto g_z_z_xxxx_yyz = cbuffer.data(gf_geom_11_off + 1207 * ccomps * dcomps);

            auto g_z_z_xxxx_yzz = cbuffer.data(gf_geom_11_off + 1208 * ccomps * dcomps);

            auto g_z_z_xxxx_zzz = cbuffer.data(gf_geom_11_off + 1209 * ccomps * dcomps);

            auto g_z_z_xxxy_xxx = cbuffer.data(gf_geom_11_off + 1210 * ccomps * dcomps);

            auto g_z_z_xxxy_xxy = cbuffer.data(gf_geom_11_off + 1211 * ccomps * dcomps);

            auto g_z_z_xxxy_xxz = cbuffer.data(gf_geom_11_off + 1212 * ccomps * dcomps);

            auto g_z_z_xxxy_xyy = cbuffer.data(gf_geom_11_off + 1213 * ccomps * dcomps);

            auto g_z_z_xxxy_xyz = cbuffer.data(gf_geom_11_off + 1214 * ccomps * dcomps);

            auto g_z_z_xxxy_xzz = cbuffer.data(gf_geom_11_off + 1215 * ccomps * dcomps);

            auto g_z_z_xxxy_yyy = cbuffer.data(gf_geom_11_off + 1216 * ccomps * dcomps);

            auto g_z_z_xxxy_yyz = cbuffer.data(gf_geom_11_off + 1217 * ccomps * dcomps);

            auto g_z_z_xxxy_yzz = cbuffer.data(gf_geom_11_off + 1218 * ccomps * dcomps);

            auto g_z_z_xxxy_zzz = cbuffer.data(gf_geom_11_off + 1219 * ccomps * dcomps);

            auto g_z_z_xxxz_xxx = cbuffer.data(gf_geom_11_off + 1220 * ccomps * dcomps);

            auto g_z_z_xxxz_xxy = cbuffer.data(gf_geom_11_off + 1221 * ccomps * dcomps);

            auto g_z_z_xxxz_xxz = cbuffer.data(gf_geom_11_off + 1222 * ccomps * dcomps);

            auto g_z_z_xxxz_xyy = cbuffer.data(gf_geom_11_off + 1223 * ccomps * dcomps);

            auto g_z_z_xxxz_xyz = cbuffer.data(gf_geom_11_off + 1224 * ccomps * dcomps);

            auto g_z_z_xxxz_xzz = cbuffer.data(gf_geom_11_off + 1225 * ccomps * dcomps);

            auto g_z_z_xxxz_yyy = cbuffer.data(gf_geom_11_off + 1226 * ccomps * dcomps);

            auto g_z_z_xxxz_yyz = cbuffer.data(gf_geom_11_off + 1227 * ccomps * dcomps);

            auto g_z_z_xxxz_yzz = cbuffer.data(gf_geom_11_off + 1228 * ccomps * dcomps);

            auto g_z_z_xxxz_zzz = cbuffer.data(gf_geom_11_off + 1229 * ccomps * dcomps);

            auto g_z_z_xxyy_xxx = cbuffer.data(gf_geom_11_off + 1230 * ccomps * dcomps);

            auto g_z_z_xxyy_xxy = cbuffer.data(gf_geom_11_off + 1231 * ccomps * dcomps);

            auto g_z_z_xxyy_xxz = cbuffer.data(gf_geom_11_off + 1232 * ccomps * dcomps);

            auto g_z_z_xxyy_xyy = cbuffer.data(gf_geom_11_off + 1233 * ccomps * dcomps);

            auto g_z_z_xxyy_xyz = cbuffer.data(gf_geom_11_off + 1234 * ccomps * dcomps);

            auto g_z_z_xxyy_xzz = cbuffer.data(gf_geom_11_off + 1235 * ccomps * dcomps);

            auto g_z_z_xxyy_yyy = cbuffer.data(gf_geom_11_off + 1236 * ccomps * dcomps);

            auto g_z_z_xxyy_yyz = cbuffer.data(gf_geom_11_off + 1237 * ccomps * dcomps);

            auto g_z_z_xxyy_yzz = cbuffer.data(gf_geom_11_off + 1238 * ccomps * dcomps);

            auto g_z_z_xxyy_zzz = cbuffer.data(gf_geom_11_off + 1239 * ccomps * dcomps);

            auto g_z_z_xxyz_xxx = cbuffer.data(gf_geom_11_off + 1240 * ccomps * dcomps);

            auto g_z_z_xxyz_xxy = cbuffer.data(gf_geom_11_off + 1241 * ccomps * dcomps);

            auto g_z_z_xxyz_xxz = cbuffer.data(gf_geom_11_off + 1242 * ccomps * dcomps);

            auto g_z_z_xxyz_xyy = cbuffer.data(gf_geom_11_off + 1243 * ccomps * dcomps);

            auto g_z_z_xxyz_xyz = cbuffer.data(gf_geom_11_off + 1244 * ccomps * dcomps);

            auto g_z_z_xxyz_xzz = cbuffer.data(gf_geom_11_off + 1245 * ccomps * dcomps);

            auto g_z_z_xxyz_yyy = cbuffer.data(gf_geom_11_off + 1246 * ccomps * dcomps);

            auto g_z_z_xxyz_yyz = cbuffer.data(gf_geom_11_off + 1247 * ccomps * dcomps);

            auto g_z_z_xxyz_yzz = cbuffer.data(gf_geom_11_off + 1248 * ccomps * dcomps);

            auto g_z_z_xxyz_zzz = cbuffer.data(gf_geom_11_off + 1249 * ccomps * dcomps);

            auto g_z_z_xxzz_xxx = cbuffer.data(gf_geom_11_off + 1250 * ccomps * dcomps);

            auto g_z_z_xxzz_xxy = cbuffer.data(gf_geom_11_off + 1251 * ccomps * dcomps);

            auto g_z_z_xxzz_xxz = cbuffer.data(gf_geom_11_off + 1252 * ccomps * dcomps);

            auto g_z_z_xxzz_xyy = cbuffer.data(gf_geom_11_off + 1253 * ccomps * dcomps);

            auto g_z_z_xxzz_xyz = cbuffer.data(gf_geom_11_off + 1254 * ccomps * dcomps);

            auto g_z_z_xxzz_xzz = cbuffer.data(gf_geom_11_off + 1255 * ccomps * dcomps);

            auto g_z_z_xxzz_yyy = cbuffer.data(gf_geom_11_off + 1256 * ccomps * dcomps);

            auto g_z_z_xxzz_yyz = cbuffer.data(gf_geom_11_off + 1257 * ccomps * dcomps);

            auto g_z_z_xxzz_yzz = cbuffer.data(gf_geom_11_off + 1258 * ccomps * dcomps);

            auto g_z_z_xxzz_zzz = cbuffer.data(gf_geom_11_off + 1259 * ccomps * dcomps);

            auto g_z_z_xyyy_xxx = cbuffer.data(gf_geom_11_off + 1260 * ccomps * dcomps);

            auto g_z_z_xyyy_xxy = cbuffer.data(gf_geom_11_off + 1261 * ccomps * dcomps);

            auto g_z_z_xyyy_xxz = cbuffer.data(gf_geom_11_off + 1262 * ccomps * dcomps);

            auto g_z_z_xyyy_xyy = cbuffer.data(gf_geom_11_off + 1263 * ccomps * dcomps);

            auto g_z_z_xyyy_xyz = cbuffer.data(gf_geom_11_off + 1264 * ccomps * dcomps);

            auto g_z_z_xyyy_xzz = cbuffer.data(gf_geom_11_off + 1265 * ccomps * dcomps);

            auto g_z_z_xyyy_yyy = cbuffer.data(gf_geom_11_off + 1266 * ccomps * dcomps);

            auto g_z_z_xyyy_yyz = cbuffer.data(gf_geom_11_off + 1267 * ccomps * dcomps);

            auto g_z_z_xyyy_yzz = cbuffer.data(gf_geom_11_off + 1268 * ccomps * dcomps);

            auto g_z_z_xyyy_zzz = cbuffer.data(gf_geom_11_off + 1269 * ccomps * dcomps);

            auto g_z_z_xyyz_xxx = cbuffer.data(gf_geom_11_off + 1270 * ccomps * dcomps);

            auto g_z_z_xyyz_xxy = cbuffer.data(gf_geom_11_off + 1271 * ccomps * dcomps);

            auto g_z_z_xyyz_xxz = cbuffer.data(gf_geom_11_off + 1272 * ccomps * dcomps);

            auto g_z_z_xyyz_xyy = cbuffer.data(gf_geom_11_off + 1273 * ccomps * dcomps);

            auto g_z_z_xyyz_xyz = cbuffer.data(gf_geom_11_off + 1274 * ccomps * dcomps);

            auto g_z_z_xyyz_xzz = cbuffer.data(gf_geom_11_off + 1275 * ccomps * dcomps);

            auto g_z_z_xyyz_yyy = cbuffer.data(gf_geom_11_off + 1276 * ccomps * dcomps);

            auto g_z_z_xyyz_yyz = cbuffer.data(gf_geom_11_off + 1277 * ccomps * dcomps);

            auto g_z_z_xyyz_yzz = cbuffer.data(gf_geom_11_off + 1278 * ccomps * dcomps);

            auto g_z_z_xyyz_zzz = cbuffer.data(gf_geom_11_off + 1279 * ccomps * dcomps);

            auto g_z_z_xyzz_xxx = cbuffer.data(gf_geom_11_off + 1280 * ccomps * dcomps);

            auto g_z_z_xyzz_xxy = cbuffer.data(gf_geom_11_off + 1281 * ccomps * dcomps);

            auto g_z_z_xyzz_xxz = cbuffer.data(gf_geom_11_off + 1282 * ccomps * dcomps);

            auto g_z_z_xyzz_xyy = cbuffer.data(gf_geom_11_off + 1283 * ccomps * dcomps);

            auto g_z_z_xyzz_xyz = cbuffer.data(gf_geom_11_off + 1284 * ccomps * dcomps);

            auto g_z_z_xyzz_xzz = cbuffer.data(gf_geom_11_off + 1285 * ccomps * dcomps);

            auto g_z_z_xyzz_yyy = cbuffer.data(gf_geom_11_off + 1286 * ccomps * dcomps);

            auto g_z_z_xyzz_yyz = cbuffer.data(gf_geom_11_off + 1287 * ccomps * dcomps);

            auto g_z_z_xyzz_yzz = cbuffer.data(gf_geom_11_off + 1288 * ccomps * dcomps);

            auto g_z_z_xyzz_zzz = cbuffer.data(gf_geom_11_off + 1289 * ccomps * dcomps);

            auto g_z_z_xzzz_xxx = cbuffer.data(gf_geom_11_off + 1290 * ccomps * dcomps);

            auto g_z_z_xzzz_xxy = cbuffer.data(gf_geom_11_off + 1291 * ccomps * dcomps);

            auto g_z_z_xzzz_xxz = cbuffer.data(gf_geom_11_off + 1292 * ccomps * dcomps);

            auto g_z_z_xzzz_xyy = cbuffer.data(gf_geom_11_off + 1293 * ccomps * dcomps);

            auto g_z_z_xzzz_xyz = cbuffer.data(gf_geom_11_off + 1294 * ccomps * dcomps);

            auto g_z_z_xzzz_xzz = cbuffer.data(gf_geom_11_off + 1295 * ccomps * dcomps);

            auto g_z_z_xzzz_yyy = cbuffer.data(gf_geom_11_off + 1296 * ccomps * dcomps);

            auto g_z_z_xzzz_yyz = cbuffer.data(gf_geom_11_off + 1297 * ccomps * dcomps);

            auto g_z_z_xzzz_yzz = cbuffer.data(gf_geom_11_off + 1298 * ccomps * dcomps);

            auto g_z_z_xzzz_zzz = cbuffer.data(gf_geom_11_off + 1299 * ccomps * dcomps);

            auto g_z_z_yyyy_xxx = cbuffer.data(gf_geom_11_off + 1300 * ccomps * dcomps);

            auto g_z_z_yyyy_xxy = cbuffer.data(gf_geom_11_off + 1301 * ccomps * dcomps);

            auto g_z_z_yyyy_xxz = cbuffer.data(gf_geom_11_off + 1302 * ccomps * dcomps);

            auto g_z_z_yyyy_xyy = cbuffer.data(gf_geom_11_off + 1303 * ccomps * dcomps);

            auto g_z_z_yyyy_xyz = cbuffer.data(gf_geom_11_off + 1304 * ccomps * dcomps);

            auto g_z_z_yyyy_xzz = cbuffer.data(gf_geom_11_off + 1305 * ccomps * dcomps);

            auto g_z_z_yyyy_yyy = cbuffer.data(gf_geom_11_off + 1306 * ccomps * dcomps);

            auto g_z_z_yyyy_yyz = cbuffer.data(gf_geom_11_off + 1307 * ccomps * dcomps);

            auto g_z_z_yyyy_yzz = cbuffer.data(gf_geom_11_off + 1308 * ccomps * dcomps);

            auto g_z_z_yyyy_zzz = cbuffer.data(gf_geom_11_off + 1309 * ccomps * dcomps);

            auto g_z_z_yyyz_xxx = cbuffer.data(gf_geom_11_off + 1310 * ccomps * dcomps);

            auto g_z_z_yyyz_xxy = cbuffer.data(gf_geom_11_off + 1311 * ccomps * dcomps);

            auto g_z_z_yyyz_xxz = cbuffer.data(gf_geom_11_off + 1312 * ccomps * dcomps);

            auto g_z_z_yyyz_xyy = cbuffer.data(gf_geom_11_off + 1313 * ccomps * dcomps);

            auto g_z_z_yyyz_xyz = cbuffer.data(gf_geom_11_off + 1314 * ccomps * dcomps);

            auto g_z_z_yyyz_xzz = cbuffer.data(gf_geom_11_off + 1315 * ccomps * dcomps);

            auto g_z_z_yyyz_yyy = cbuffer.data(gf_geom_11_off + 1316 * ccomps * dcomps);

            auto g_z_z_yyyz_yyz = cbuffer.data(gf_geom_11_off + 1317 * ccomps * dcomps);

            auto g_z_z_yyyz_yzz = cbuffer.data(gf_geom_11_off + 1318 * ccomps * dcomps);

            auto g_z_z_yyyz_zzz = cbuffer.data(gf_geom_11_off + 1319 * ccomps * dcomps);

            auto g_z_z_yyzz_xxx = cbuffer.data(gf_geom_11_off + 1320 * ccomps * dcomps);

            auto g_z_z_yyzz_xxy = cbuffer.data(gf_geom_11_off + 1321 * ccomps * dcomps);

            auto g_z_z_yyzz_xxz = cbuffer.data(gf_geom_11_off + 1322 * ccomps * dcomps);

            auto g_z_z_yyzz_xyy = cbuffer.data(gf_geom_11_off + 1323 * ccomps * dcomps);

            auto g_z_z_yyzz_xyz = cbuffer.data(gf_geom_11_off + 1324 * ccomps * dcomps);

            auto g_z_z_yyzz_xzz = cbuffer.data(gf_geom_11_off + 1325 * ccomps * dcomps);

            auto g_z_z_yyzz_yyy = cbuffer.data(gf_geom_11_off + 1326 * ccomps * dcomps);

            auto g_z_z_yyzz_yyz = cbuffer.data(gf_geom_11_off + 1327 * ccomps * dcomps);

            auto g_z_z_yyzz_yzz = cbuffer.data(gf_geom_11_off + 1328 * ccomps * dcomps);

            auto g_z_z_yyzz_zzz = cbuffer.data(gf_geom_11_off + 1329 * ccomps * dcomps);

            auto g_z_z_yzzz_xxx = cbuffer.data(gf_geom_11_off + 1330 * ccomps * dcomps);

            auto g_z_z_yzzz_xxy = cbuffer.data(gf_geom_11_off + 1331 * ccomps * dcomps);

            auto g_z_z_yzzz_xxz = cbuffer.data(gf_geom_11_off + 1332 * ccomps * dcomps);

            auto g_z_z_yzzz_xyy = cbuffer.data(gf_geom_11_off + 1333 * ccomps * dcomps);

            auto g_z_z_yzzz_xyz = cbuffer.data(gf_geom_11_off + 1334 * ccomps * dcomps);

            auto g_z_z_yzzz_xzz = cbuffer.data(gf_geom_11_off + 1335 * ccomps * dcomps);

            auto g_z_z_yzzz_yyy = cbuffer.data(gf_geom_11_off + 1336 * ccomps * dcomps);

            auto g_z_z_yzzz_yyz = cbuffer.data(gf_geom_11_off + 1337 * ccomps * dcomps);

            auto g_z_z_yzzz_yzz = cbuffer.data(gf_geom_11_off + 1338 * ccomps * dcomps);

            auto g_z_z_yzzz_zzz = cbuffer.data(gf_geom_11_off + 1339 * ccomps * dcomps);

            auto g_z_z_zzzz_xxx = cbuffer.data(gf_geom_11_off + 1340 * ccomps * dcomps);

            auto g_z_z_zzzz_xxy = cbuffer.data(gf_geom_11_off + 1341 * ccomps * dcomps);

            auto g_z_z_zzzz_xxz = cbuffer.data(gf_geom_11_off + 1342 * ccomps * dcomps);

            auto g_z_z_zzzz_xyy = cbuffer.data(gf_geom_11_off + 1343 * ccomps * dcomps);

            auto g_z_z_zzzz_xyz = cbuffer.data(gf_geom_11_off + 1344 * ccomps * dcomps);

            auto g_z_z_zzzz_xzz = cbuffer.data(gf_geom_11_off + 1345 * ccomps * dcomps);

            auto g_z_z_zzzz_yyy = cbuffer.data(gf_geom_11_off + 1346 * ccomps * dcomps);

            auto g_z_z_zzzz_yyz = cbuffer.data(gf_geom_11_off + 1347 * ccomps * dcomps);

            auto g_z_z_zzzz_yzz = cbuffer.data(gf_geom_11_off + 1348 * ccomps * dcomps);

            auto g_z_z_zzzz_zzz = cbuffer.data(gf_geom_11_off + 1349 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hdxx

            const auto hd_geom_11_off = idx_geom_11_hdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxx_xx = cbuffer.data(hd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxxx_xy = cbuffer.data(hd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxxx_xz = cbuffer.data(hd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxxxx_yy = cbuffer.data(hd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxxx_yz = cbuffer.data(hd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxxx_zz = cbuffer.data(hd_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xx, g_0_x_xxxx_xy, g_0_x_xxxx_xz, g_0_x_xxxx_yy, g_0_x_xxxx_yz, g_0_x_xxxx_zz, g_x_0_xxxx_xx, g_x_0_xxxx_xy, g_x_0_xxxx_xz, g_x_0_xxxx_yy, g_x_0_xxxx_yz, g_x_0_xxxx_zz, g_x_x_xxxx_xx, g_x_x_xxxx_xxx, g_x_x_xxxx_xxy, g_x_x_xxxx_xxz, g_x_x_xxxx_xy, g_x_x_xxxx_xyy, g_x_x_xxxx_xyz, g_x_x_xxxx_xz, g_x_x_xxxx_xzz, g_x_x_xxxx_yy, g_x_x_xxxx_yz, g_x_x_xxxx_zz, g_x_x_xxxxx_xx, g_x_x_xxxxx_xy, g_x_x_xxxxx_xz, g_x_x_xxxxx_yy, g_x_x_xxxxx_yz, g_x_x_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxx_xx[k] = -g_0_x_xxxx_xx[k] + g_x_0_xxxx_xx[k] - g_x_x_xxxx_xx[k] * ab_x + g_x_x_xxxx_xxx[k];

                g_x_x_xxxxx_xy[k] = -g_0_x_xxxx_xy[k] + g_x_0_xxxx_xy[k] - g_x_x_xxxx_xy[k] * ab_x + g_x_x_xxxx_xxy[k];

                g_x_x_xxxxx_xz[k] = -g_0_x_xxxx_xz[k] + g_x_0_xxxx_xz[k] - g_x_x_xxxx_xz[k] * ab_x + g_x_x_xxxx_xxz[k];

                g_x_x_xxxxx_yy[k] = -g_0_x_xxxx_yy[k] + g_x_0_xxxx_yy[k] - g_x_x_xxxx_yy[k] * ab_x + g_x_x_xxxx_xyy[k];

                g_x_x_xxxxx_yz[k] = -g_0_x_xxxx_yz[k] + g_x_0_xxxx_yz[k] - g_x_x_xxxx_yz[k] * ab_x + g_x_x_xxxx_xyz[k];

                g_x_x_xxxxx_zz[k] = -g_0_x_xxxx_zz[k] + g_x_0_xxxx_zz[k] - g_x_x_xxxx_zz[k] * ab_x + g_x_x_xxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxy_xx = cbuffer.data(hd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxxy_xy = cbuffer.data(hd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxxy_xz = cbuffer.data(hd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxxxy_yy = cbuffer.data(hd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxxxy_yz = cbuffer.data(hd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxxxy_zz = cbuffer.data(hd_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxx_xx, g_x_x_xxxx_xxy, g_x_x_xxxx_xy, g_x_x_xxxx_xyy, g_x_x_xxxx_xyz, g_x_x_xxxx_xz, g_x_x_xxxx_yy, g_x_x_xxxx_yyy, g_x_x_xxxx_yyz, g_x_x_xxxx_yz, g_x_x_xxxx_yzz, g_x_x_xxxx_zz, g_x_x_xxxxy_xx, g_x_x_xxxxy_xy, g_x_x_xxxxy_xz, g_x_x_xxxxy_yy, g_x_x_xxxxy_yz, g_x_x_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxy_xx[k] = -g_x_x_xxxx_xx[k] * ab_y + g_x_x_xxxx_xxy[k];

                g_x_x_xxxxy_xy[k] = -g_x_x_xxxx_xy[k] * ab_y + g_x_x_xxxx_xyy[k];

                g_x_x_xxxxy_xz[k] = -g_x_x_xxxx_xz[k] * ab_y + g_x_x_xxxx_xyz[k];

                g_x_x_xxxxy_yy[k] = -g_x_x_xxxx_yy[k] * ab_y + g_x_x_xxxx_yyy[k];

                g_x_x_xxxxy_yz[k] = -g_x_x_xxxx_yz[k] * ab_y + g_x_x_xxxx_yyz[k];

                g_x_x_xxxxy_zz[k] = -g_x_x_xxxx_zz[k] * ab_y + g_x_x_xxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxxz_xx = cbuffer.data(hd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxxxz_xy = cbuffer.data(hd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxxxz_xz = cbuffer.data(hd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxxxz_yy = cbuffer.data(hd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxxxz_yz = cbuffer.data(hd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxxxz_zz = cbuffer.data(hd_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxx_xx, g_x_x_xxxx_xxz, g_x_x_xxxx_xy, g_x_x_xxxx_xyz, g_x_x_xxxx_xz, g_x_x_xxxx_xzz, g_x_x_xxxx_yy, g_x_x_xxxx_yyz, g_x_x_xxxx_yz, g_x_x_xxxx_yzz, g_x_x_xxxx_zz, g_x_x_xxxx_zzz, g_x_x_xxxxz_xx, g_x_x_xxxxz_xy, g_x_x_xxxxz_xz, g_x_x_xxxxz_yy, g_x_x_xxxxz_yz, g_x_x_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxxz_xx[k] = -g_x_x_xxxx_xx[k] * ab_z + g_x_x_xxxx_xxz[k];

                g_x_x_xxxxz_xy[k] = -g_x_x_xxxx_xy[k] * ab_z + g_x_x_xxxx_xyz[k];

                g_x_x_xxxxz_xz[k] = -g_x_x_xxxx_xz[k] * ab_z + g_x_x_xxxx_xzz[k];

                g_x_x_xxxxz_yy[k] = -g_x_x_xxxx_yy[k] * ab_z + g_x_x_xxxx_yyz[k];

                g_x_x_xxxxz_yz[k] = -g_x_x_xxxx_yz[k] * ab_z + g_x_x_xxxx_yzz[k];

                g_x_x_xxxxz_zz[k] = -g_x_x_xxxx_zz[k] * ab_z + g_x_x_xxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxyy_xx = cbuffer.data(hd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxxyy_xy = cbuffer.data(hd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xxxyy_xz = cbuffer.data(hd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxxyy_yy = cbuffer.data(hd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxxyy_yz = cbuffer.data(hd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxxyy_zz = cbuffer.data(hd_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxy_xx, g_x_x_xxxy_xxy, g_x_x_xxxy_xy, g_x_x_xxxy_xyy, g_x_x_xxxy_xyz, g_x_x_xxxy_xz, g_x_x_xxxy_yy, g_x_x_xxxy_yyy, g_x_x_xxxy_yyz, g_x_x_xxxy_yz, g_x_x_xxxy_yzz, g_x_x_xxxy_zz, g_x_x_xxxyy_xx, g_x_x_xxxyy_xy, g_x_x_xxxyy_xz, g_x_x_xxxyy_yy, g_x_x_xxxyy_yz, g_x_x_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxyy_xx[k] = -g_x_x_xxxy_xx[k] * ab_y + g_x_x_xxxy_xxy[k];

                g_x_x_xxxyy_xy[k] = -g_x_x_xxxy_xy[k] * ab_y + g_x_x_xxxy_xyy[k];

                g_x_x_xxxyy_xz[k] = -g_x_x_xxxy_xz[k] * ab_y + g_x_x_xxxy_xyz[k];

                g_x_x_xxxyy_yy[k] = -g_x_x_xxxy_yy[k] * ab_y + g_x_x_xxxy_yyy[k];

                g_x_x_xxxyy_yz[k] = -g_x_x_xxxy_yz[k] * ab_y + g_x_x_xxxy_yyz[k];

                g_x_x_xxxyy_zz[k] = -g_x_x_xxxy_zz[k] * ab_y + g_x_x_xxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxyz_xx = cbuffer.data(hd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxxyz_xy = cbuffer.data(hd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxxyz_xz = cbuffer.data(hd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxxyz_yy = cbuffer.data(hd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxxyz_yz = cbuffer.data(hd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxxyz_zz = cbuffer.data(hd_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxyz_xx, g_x_x_xxxyz_xy, g_x_x_xxxyz_xz, g_x_x_xxxyz_yy, g_x_x_xxxyz_yz, g_x_x_xxxyz_zz, g_x_x_xxxz_xx, g_x_x_xxxz_xxy, g_x_x_xxxz_xy, g_x_x_xxxz_xyy, g_x_x_xxxz_xyz, g_x_x_xxxz_xz, g_x_x_xxxz_yy, g_x_x_xxxz_yyy, g_x_x_xxxz_yyz, g_x_x_xxxz_yz, g_x_x_xxxz_yzz, g_x_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxyz_xx[k] = -g_x_x_xxxz_xx[k] * ab_y + g_x_x_xxxz_xxy[k];

                g_x_x_xxxyz_xy[k] = -g_x_x_xxxz_xy[k] * ab_y + g_x_x_xxxz_xyy[k];

                g_x_x_xxxyz_xz[k] = -g_x_x_xxxz_xz[k] * ab_y + g_x_x_xxxz_xyz[k];

                g_x_x_xxxyz_yy[k] = -g_x_x_xxxz_yy[k] * ab_y + g_x_x_xxxz_yyy[k];

                g_x_x_xxxyz_yz[k] = -g_x_x_xxxz_yz[k] * ab_y + g_x_x_xxxz_yyz[k];

                g_x_x_xxxyz_zz[k] = -g_x_x_xxxz_zz[k] * ab_y + g_x_x_xxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxzz_xx = cbuffer.data(hd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xxxzz_xy = cbuffer.data(hd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xxxzz_xz = cbuffer.data(hd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xxxzz_yy = cbuffer.data(hd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xxxzz_yz = cbuffer.data(hd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xxxzz_zz = cbuffer.data(hd_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxxz_xx, g_x_x_xxxz_xxz, g_x_x_xxxz_xy, g_x_x_xxxz_xyz, g_x_x_xxxz_xz, g_x_x_xxxz_xzz, g_x_x_xxxz_yy, g_x_x_xxxz_yyz, g_x_x_xxxz_yz, g_x_x_xxxz_yzz, g_x_x_xxxz_zz, g_x_x_xxxz_zzz, g_x_x_xxxzz_xx, g_x_x_xxxzz_xy, g_x_x_xxxzz_xz, g_x_x_xxxzz_yy, g_x_x_xxxzz_yz, g_x_x_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxzz_xx[k] = -g_x_x_xxxz_xx[k] * ab_z + g_x_x_xxxz_xxz[k];

                g_x_x_xxxzz_xy[k] = -g_x_x_xxxz_xy[k] * ab_z + g_x_x_xxxz_xyz[k];

                g_x_x_xxxzz_xz[k] = -g_x_x_xxxz_xz[k] * ab_z + g_x_x_xxxz_xzz[k];

                g_x_x_xxxzz_yy[k] = -g_x_x_xxxz_yy[k] * ab_z + g_x_x_xxxz_yyz[k];

                g_x_x_xxxzz_yz[k] = -g_x_x_xxxz_yz[k] * ab_z + g_x_x_xxxz_yzz[k];

                g_x_x_xxxzz_zz[k] = -g_x_x_xxxz_zz[k] * ab_z + g_x_x_xxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyyy_xx = cbuffer.data(hd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xxyyy_xy = cbuffer.data(hd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xxyyy_xz = cbuffer.data(hd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xxyyy_yy = cbuffer.data(hd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xxyyy_yz = cbuffer.data(hd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xxyyy_zz = cbuffer.data(hd_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyy_xx, g_x_x_xxyy_xxy, g_x_x_xxyy_xy, g_x_x_xxyy_xyy, g_x_x_xxyy_xyz, g_x_x_xxyy_xz, g_x_x_xxyy_yy, g_x_x_xxyy_yyy, g_x_x_xxyy_yyz, g_x_x_xxyy_yz, g_x_x_xxyy_yzz, g_x_x_xxyy_zz, g_x_x_xxyyy_xx, g_x_x_xxyyy_xy, g_x_x_xxyyy_xz, g_x_x_xxyyy_yy, g_x_x_xxyyy_yz, g_x_x_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyyy_xx[k] = -g_x_x_xxyy_xx[k] * ab_y + g_x_x_xxyy_xxy[k];

                g_x_x_xxyyy_xy[k] = -g_x_x_xxyy_xy[k] * ab_y + g_x_x_xxyy_xyy[k];

                g_x_x_xxyyy_xz[k] = -g_x_x_xxyy_xz[k] * ab_y + g_x_x_xxyy_xyz[k];

                g_x_x_xxyyy_yy[k] = -g_x_x_xxyy_yy[k] * ab_y + g_x_x_xxyy_yyy[k];

                g_x_x_xxyyy_yz[k] = -g_x_x_xxyy_yz[k] * ab_y + g_x_x_xxyy_yyz[k];

                g_x_x_xxyyy_zz[k] = -g_x_x_xxyy_zz[k] * ab_y + g_x_x_xxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyyz_xx = cbuffer.data(hd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xxyyz_xy = cbuffer.data(hd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xxyyz_xz = cbuffer.data(hd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xxyyz_yy = cbuffer.data(hd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xxyyz_yz = cbuffer.data(hd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xxyyz_zz = cbuffer.data(hd_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyyz_xx, g_x_x_xxyyz_xy, g_x_x_xxyyz_xz, g_x_x_xxyyz_yy, g_x_x_xxyyz_yz, g_x_x_xxyyz_zz, g_x_x_xxyz_xx, g_x_x_xxyz_xxy, g_x_x_xxyz_xy, g_x_x_xxyz_xyy, g_x_x_xxyz_xyz, g_x_x_xxyz_xz, g_x_x_xxyz_yy, g_x_x_xxyz_yyy, g_x_x_xxyz_yyz, g_x_x_xxyz_yz, g_x_x_xxyz_yzz, g_x_x_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyyz_xx[k] = -g_x_x_xxyz_xx[k] * ab_y + g_x_x_xxyz_xxy[k];

                g_x_x_xxyyz_xy[k] = -g_x_x_xxyz_xy[k] * ab_y + g_x_x_xxyz_xyy[k];

                g_x_x_xxyyz_xz[k] = -g_x_x_xxyz_xz[k] * ab_y + g_x_x_xxyz_xyz[k];

                g_x_x_xxyyz_yy[k] = -g_x_x_xxyz_yy[k] * ab_y + g_x_x_xxyz_yyy[k];

                g_x_x_xxyyz_yz[k] = -g_x_x_xxyz_yz[k] * ab_y + g_x_x_xxyz_yyz[k];

                g_x_x_xxyyz_zz[k] = -g_x_x_xxyz_zz[k] * ab_y + g_x_x_xxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyzz_xx = cbuffer.data(hd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xxyzz_xy = cbuffer.data(hd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xxyzz_xz = cbuffer.data(hd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xxyzz_yy = cbuffer.data(hd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xxyzz_yz = cbuffer.data(hd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xxyzz_zz = cbuffer.data(hd_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyzz_xx, g_x_x_xxyzz_xy, g_x_x_xxyzz_xz, g_x_x_xxyzz_yy, g_x_x_xxyzz_yz, g_x_x_xxyzz_zz, g_x_x_xxzz_xx, g_x_x_xxzz_xxy, g_x_x_xxzz_xy, g_x_x_xxzz_xyy, g_x_x_xxzz_xyz, g_x_x_xxzz_xz, g_x_x_xxzz_yy, g_x_x_xxzz_yyy, g_x_x_xxzz_yyz, g_x_x_xxzz_yz, g_x_x_xxzz_yzz, g_x_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyzz_xx[k] = -g_x_x_xxzz_xx[k] * ab_y + g_x_x_xxzz_xxy[k];

                g_x_x_xxyzz_xy[k] = -g_x_x_xxzz_xy[k] * ab_y + g_x_x_xxzz_xyy[k];

                g_x_x_xxyzz_xz[k] = -g_x_x_xxzz_xz[k] * ab_y + g_x_x_xxzz_xyz[k];

                g_x_x_xxyzz_yy[k] = -g_x_x_xxzz_yy[k] * ab_y + g_x_x_xxzz_yyy[k];

                g_x_x_xxyzz_yz[k] = -g_x_x_xxzz_yz[k] * ab_y + g_x_x_xxzz_yyz[k];

                g_x_x_xxyzz_zz[k] = -g_x_x_xxzz_zz[k] * ab_y + g_x_x_xxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxzzz_xx = cbuffer.data(hd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xxzzz_xy = cbuffer.data(hd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xxzzz_xz = cbuffer.data(hd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xxzzz_yy = cbuffer.data(hd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xxzzz_yz = cbuffer.data(hd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xxzzz_zz = cbuffer.data(hd_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxzz_xx, g_x_x_xxzz_xxz, g_x_x_xxzz_xy, g_x_x_xxzz_xyz, g_x_x_xxzz_xz, g_x_x_xxzz_xzz, g_x_x_xxzz_yy, g_x_x_xxzz_yyz, g_x_x_xxzz_yz, g_x_x_xxzz_yzz, g_x_x_xxzz_zz, g_x_x_xxzz_zzz, g_x_x_xxzzz_xx, g_x_x_xxzzz_xy, g_x_x_xxzzz_xz, g_x_x_xxzzz_yy, g_x_x_xxzzz_yz, g_x_x_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxzzz_xx[k] = -g_x_x_xxzz_xx[k] * ab_z + g_x_x_xxzz_xxz[k];

                g_x_x_xxzzz_xy[k] = -g_x_x_xxzz_xy[k] * ab_z + g_x_x_xxzz_xyz[k];

                g_x_x_xxzzz_xz[k] = -g_x_x_xxzz_xz[k] * ab_z + g_x_x_xxzz_xzz[k];

                g_x_x_xxzzz_yy[k] = -g_x_x_xxzz_yy[k] * ab_z + g_x_x_xxzz_yyz[k];

                g_x_x_xxzzz_yz[k] = -g_x_x_xxzz_yz[k] * ab_z + g_x_x_xxzz_yzz[k];

                g_x_x_xxzzz_zz[k] = -g_x_x_xxzz_zz[k] * ab_z + g_x_x_xxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyyy_xx = cbuffer.data(hd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xyyyy_xy = cbuffer.data(hd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xyyyy_xz = cbuffer.data(hd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_xyyyy_yy = cbuffer.data(hd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_xyyyy_yz = cbuffer.data(hd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_xyyyy_zz = cbuffer.data(hd_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyy_xx, g_x_x_xyyy_xxy, g_x_x_xyyy_xy, g_x_x_xyyy_xyy, g_x_x_xyyy_xyz, g_x_x_xyyy_xz, g_x_x_xyyy_yy, g_x_x_xyyy_yyy, g_x_x_xyyy_yyz, g_x_x_xyyy_yz, g_x_x_xyyy_yzz, g_x_x_xyyy_zz, g_x_x_xyyyy_xx, g_x_x_xyyyy_xy, g_x_x_xyyyy_xz, g_x_x_xyyyy_yy, g_x_x_xyyyy_yz, g_x_x_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyyy_xx[k] = -g_x_x_xyyy_xx[k] * ab_y + g_x_x_xyyy_xxy[k];

                g_x_x_xyyyy_xy[k] = -g_x_x_xyyy_xy[k] * ab_y + g_x_x_xyyy_xyy[k];

                g_x_x_xyyyy_xz[k] = -g_x_x_xyyy_xz[k] * ab_y + g_x_x_xyyy_xyz[k];

                g_x_x_xyyyy_yy[k] = -g_x_x_xyyy_yy[k] * ab_y + g_x_x_xyyy_yyy[k];

                g_x_x_xyyyy_yz[k] = -g_x_x_xyyy_yz[k] * ab_y + g_x_x_xyyy_yyz[k];

                g_x_x_xyyyy_zz[k] = -g_x_x_xyyy_zz[k] * ab_y + g_x_x_xyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyyz_xx = cbuffer.data(hd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_xyyyz_xy = cbuffer.data(hd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_xyyyz_xz = cbuffer.data(hd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_xyyyz_yy = cbuffer.data(hd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_xyyyz_yz = cbuffer.data(hd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_xyyyz_zz = cbuffer.data(hd_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyyz_xx, g_x_x_xyyyz_xy, g_x_x_xyyyz_xz, g_x_x_xyyyz_yy, g_x_x_xyyyz_yz, g_x_x_xyyyz_zz, g_x_x_xyyz_xx, g_x_x_xyyz_xxy, g_x_x_xyyz_xy, g_x_x_xyyz_xyy, g_x_x_xyyz_xyz, g_x_x_xyyz_xz, g_x_x_xyyz_yy, g_x_x_xyyz_yyy, g_x_x_xyyz_yyz, g_x_x_xyyz_yz, g_x_x_xyyz_yzz, g_x_x_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyyz_xx[k] = -g_x_x_xyyz_xx[k] * ab_y + g_x_x_xyyz_xxy[k];

                g_x_x_xyyyz_xy[k] = -g_x_x_xyyz_xy[k] * ab_y + g_x_x_xyyz_xyy[k];

                g_x_x_xyyyz_xz[k] = -g_x_x_xyyz_xz[k] * ab_y + g_x_x_xyyz_xyz[k];

                g_x_x_xyyyz_yy[k] = -g_x_x_xyyz_yy[k] * ab_y + g_x_x_xyyz_yyy[k];

                g_x_x_xyyyz_yz[k] = -g_x_x_xyyz_yz[k] * ab_y + g_x_x_xyyz_yyz[k];

                g_x_x_xyyyz_zz[k] = -g_x_x_xyyz_zz[k] * ab_y + g_x_x_xyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyzz_xx = cbuffer.data(hd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_xyyzz_xy = cbuffer.data(hd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_xyyzz_xz = cbuffer.data(hd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_xyyzz_yy = cbuffer.data(hd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_xyyzz_yz = cbuffer.data(hd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_xyyzz_zz = cbuffer.data(hd_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyzz_xx, g_x_x_xyyzz_xy, g_x_x_xyyzz_xz, g_x_x_xyyzz_yy, g_x_x_xyyzz_yz, g_x_x_xyyzz_zz, g_x_x_xyzz_xx, g_x_x_xyzz_xxy, g_x_x_xyzz_xy, g_x_x_xyzz_xyy, g_x_x_xyzz_xyz, g_x_x_xyzz_xz, g_x_x_xyzz_yy, g_x_x_xyzz_yyy, g_x_x_xyzz_yyz, g_x_x_xyzz_yz, g_x_x_xyzz_yzz, g_x_x_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyzz_xx[k] = -g_x_x_xyzz_xx[k] * ab_y + g_x_x_xyzz_xxy[k];

                g_x_x_xyyzz_xy[k] = -g_x_x_xyzz_xy[k] * ab_y + g_x_x_xyzz_xyy[k];

                g_x_x_xyyzz_xz[k] = -g_x_x_xyzz_xz[k] * ab_y + g_x_x_xyzz_xyz[k];

                g_x_x_xyyzz_yy[k] = -g_x_x_xyzz_yy[k] * ab_y + g_x_x_xyzz_yyy[k];

                g_x_x_xyyzz_yz[k] = -g_x_x_xyzz_yz[k] * ab_y + g_x_x_xyzz_yyz[k];

                g_x_x_xyyzz_zz[k] = -g_x_x_xyzz_zz[k] * ab_y + g_x_x_xyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyzzz_xx = cbuffer.data(hd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_xyzzz_xy = cbuffer.data(hd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_xyzzz_xz = cbuffer.data(hd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_xyzzz_yy = cbuffer.data(hd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_xyzzz_yz = cbuffer.data(hd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_xyzzz_zz = cbuffer.data(hd_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyzzz_xx, g_x_x_xyzzz_xy, g_x_x_xyzzz_xz, g_x_x_xyzzz_yy, g_x_x_xyzzz_yz, g_x_x_xyzzz_zz, g_x_x_xzzz_xx, g_x_x_xzzz_xxy, g_x_x_xzzz_xy, g_x_x_xzzz_xyy, g_x_x_xzzz_xyz, g_x_x_xzzz_xz, g_x_x_xzzz_yy, g_x_x_xzzz_yyy, g_x_x_xzzz_yyz, g_x_x_xzzz_yz, g_x_x_xzzz_yzz, g_x_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyzzz_xx[k] = -g_x_x_xzzz_xx[k] * ab_y + g_x_x_xzzz_xxy[k];

                g_x_x_xyzzz_xy[k] = -g_x_x_xzzz_xy[k] * ab_y + g_x_x_xzzz_xyy[k];

                g_x_x_xyzzz_xz[k] = -g_x_x_xzzz_xz[k] * ab_y + g_x_x_xzzz_xyz[k];

                g_x_x_xyzzz_yy[k] = -g_x_x_xzzz_yy[k] * ab_y + g_x_x_xzzz_yyy[k];

                g_x_x_xyzzz_yz[k] = -g_x_x_xzzz_yz[k] * ab_y + g_x_x_xzzz_yyz[k];

                g_x_x_xyzzz_zz[k] = -g_x_x_xzzz_zz[k] * ab_y + g_x_x_xzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzzzz_xx = cbuffer.data(hd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_xzzzz_xy = cbuffer.data(hd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_xzzzz_xz = cbuffer.data(hd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_xzzzz_yy = cbuffer.data(hd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_xzzzz_yz = cbuffer.data(hd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_xzzzz_zz = cbuffer.data(hd_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xzzz_xx, g_x_x_xzzz_xxz, g_x_x_xzzz_xy, g_x_x_xzzz_xyz, g_x_x_xzzz_xz, g_x_x_xzzz_xzz, g_x_x_xzzz_yy, g_x_x_xzzz_yyz, g_x_x_xzzz_yz, g_x_x_xzzz_yzz, g_x_x_xzzz_zz, g_x_x_xzzz_zzz, g_x_x_xzzzz_xx, g_x_x_xzzzz_xy, g_x_x_xzzzz_xz, g_x_x_xzzzz_yy, g_x_x_xzzzz_yz, g_x_x_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzzzz_xx[k] = -g_x_x_xzzz_xx[k] * ab_z + g_x_x_xzzz_xxz[k];

                g_x_x_xzzzz_xy[k] = -g_x_x_xzzz_xy[k] * ab_z + g_x_x_xzzz_xyz[k];

                g_x_x_xzzzz_xz[k] = -g_x_x_xzzz_xz[k] * ab_z + g_x_x_xzzz_xzz[k];

                g_x_x_xzzzz_yy[k] = -g_x_x_xzzz_yy[k] * ab_z + g_x_x_xzzz_yyz[k];

                g_x_x_xzzzz_yz[k] = -g_x_x_xzzz_yz[k] * ab_z + g_x_x_xzzz_yzz[k];

                g_x_x_xzzzz_zz[k] = -g_x_x_xzzz_zz[k] * ab_z + g_x_x_xzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyyy_xx = cbuffer.data(hd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_yyyyy_xy = cbuffer.data(hd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_yyyyy_xz = cbuffer.data(hd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_yyyyy_yy = cbuffer.data(hd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_yyyyy_yz = cbuffer.data(hd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_yyyyy_zz = cbuffer.data(hd_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyy_xx, g_x_x_yyyy_xxy, g_x_x_yyyy_xy, g_x_x_yyyy_xyy, g_x_x_yyyy_xyz, g_x_x_yyyy_xz, g_x_x_yyyy_yy, g_x_x_yyyy_yyy, g_x_x_yyyy_yyz, g_x_x_yyyy_yz, g_x_x_yyyy_yzz, g_x_x_yyyy_zz, g_x_x_yyyyy_xx, g_x_x_yyyyy_xy, g_x_x_yyyyy_xz, g_x_x_yyyyy_yy, g_x_x_yyyyy_yz, g_x_x_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyyy_xx[k] = -g_x_x_yyyy_xx[k] * ab_y + g_x_x_yyyy_xxy[k];

                g_x_x_yyyyy_xy[k] = -g_x_x_yyyy_xy[k] * ab_y + g_x_x_yyyy_xyy[k];

                g_x_x_yyyyy_xz[k] = -g_x_x_yyyy_xz[k] * ab_y + g_x_x_yyyy_xyz[k];

                g_x_x_yyyyy_yy[k] = -g_x_x_yyyy_yy[k] * ab_y + g_x_x_yyyy_yyy[k];

                g_x_x_yyyyy_yz[k] = -g_x_x_yyyy_yz[k] * ab_y + g_x_x_yyyy_yyz[k];

                g_x_x_yyyyy_zz[k] = -g_x_x_yyyy_zz[k] * ab_y + g_x_x_yyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyyz_xx = cbuffer.data(hd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_yyyyz_xy = cbuffer.data(hd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_yyyyz_xz = cbuffer.data(hd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_yyyyz_yy = cbuffer.data(hd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yyyyz_yz = cbuffer.data(hd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yyyyz_zz = cbuffer.data(hd_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyyz_xx, g_x_x_yyyyz_xy, g_x_x_yyyyz_xz, g_x_x_yyyyz_yy, g_x_x_yyyyz_yz, g_x_x_yyyyz_zz, g_x_x_yyyz_xx, g_x_x_yyyz_xxy, g_x_x_yyyz_xy, g_x_x_yyyz_xyy, g_x_x_yyyz_xyz, g_x_x_yyyz_xz, g_x_x_yyyz_yy, g_x_x_yyyz_yyy, g_x_x_yyyz_yyz, g_x_x_yyyz_yz, g_x_x_yyyz_yzz, g_x_x_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyyz_xx[k] = -g_x_x_yyyz_xx[k] * ab_y + g_x_x_yyyz_xxy[k];

                g_x_x_yyyyz_xy[k] = -g_x_x_yyyz_xy[k] * ab_y + g_x_x_yyyz_xyy[k];

                g_x_x_yyyyz_xz[k] = -g_x_x_yyyz_xz[k] * ab_y + g_x_x_yyyz_xyz[k];

                g_x_x_yyyyz_yy[k] = -g_x_x_yyyz_yy[k] * ab_y + g_x_x_yyyz_yyy[k];

                g_x_x_yyyyz_yz[k] = -g_x_x_yyyz_yz[k] * ab_y + g_x_x_yyyz_yyz[k];

                g_x_x_yyyyz_zz[k] = -g_x_x_yyyz_zz[k] * ab_y + g_x_x_yyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyzz_xx = cbuffer.data(hd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yyyzz_xy = cbuffer.data(hd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yyyzz_xz = cbuffer.data(hd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_x_yyyzz_yy = cbuffer.data(hd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_yyyzz_yz = cbuffer.data(hd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_yyyzz_zz = cbuffer.data(hd_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyzz_xx, g_x_x_yyyzz_xy, g_x_x_yyyzz_xz, g_x_x_yyyzz_yy, g_x_x_yyyzz_yz, g_x_x_yyyzz_zz, g_x_x_yyzz_xx, g_x_x_yyzz_xxy, g_x_x_yyzz_xy, g_x_x_yyzz_xyy, g_x_x_yyzz_xyz, g_x_x_yyzz_xz, g_x_x_yyzz_yy, g_x_x_yyzz_yyy, g_x_x_yyzz_yyz, g_x_x_yyzz_yz, g_x_x_yyzz_yzz, g_x_x_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyzz_xx[k] = -g_x_x_yyzz_xx[k] * ab_y + g_x_x_yyzz_xxy[k];

                g_x_x_yyyzz_xy[k] = -g_x_x_yyzz_xy[k] * ab_y + g_x_x_yyzz_xyy[k];

                g_x_x_yyyzz_xz[k] = -g_x_x_yyzz_xz[k] * ab_y + g_x_x_yyzz_xyz[k];

                g_x_x_yyyzz_yy[k] = -g_x_x_yyzz_yy[k] * ab_y + g_x_x_yyzz_yyy[k];

                g_x_x_yyyzz_yz[k] = -g_x_x_yyzz_yz[k] * ab_y + g_x_x_yyzz_yyz[k];

                g_x_x_yyyzz_zz[k] = -g_x_x_yyzz_zz[k] * ab_y + g_x_x_yyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyzzz_xx = cbuffer.data(hd_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_yyzzz_xy = cbuffer.data(hd_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_yyzzz_xz = cbuffer.data(hd_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_yyzzz_yy = cbuffer.data(hd_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_x_yyzzz_yz = cbuffer.data(hd_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_yyzzz_zz = cbuffer.data(hd_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyzzz_xx, g_x_x_yyzzz_xy, g_x_x_yyzzz_xz, g_x_x_yyzzz_yy, g_x_x_yyzzz_yz, g_x_x_yyzzz_zz, g_x_x_yzzz_xx, g_x_x_yzzz_xxy, g_x_x_yzzz_xy, g_x_x_yzzz_xyy, g_x_x_yzzz_xyz, g_x_x_yzzz_xz, g_x_x_yzzz_yy, g_x_x_yzzz_yyy, g_x_x_yzzz_yyz, g_x_x_yzzz_yz, g_x_x_yzzz_yzz, g_x_x_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyzzz_xx[k] = -g_x_x_yzzz_xx[k] * ab_y + g_x_x_yzzz_xxy[k];

                g_x_x_yyzzz_xy[k] = -g_x_x_yzzz_xy[k] * ab_y + g_x_x_yzzz_xyy[k];

                g_x_x_yyzzz_xz[k] = -g_x_x_yzzz_xz[k] * ab_y + g_x_x_yzzz_xyz[k];

                g_x_x_yyzzz_yy[k] = -g_x_x_yzzz_yy[k] * ab_y + g_x_x_yzzz_yyy[k];

                g_x_x_yyzzz_yz[k] = -g_x_x_yzzz_yz[k] * ab_y + g_x_x_yzzz_yyz[k];

                g_x_x_yyzzz_zz[k] = -g_x_x_yzzz_zz[k] * ab_y + g_x_x_yzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzzzz_xx = cbuffer.data(hd_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_yzzzz_xy = cbuffer.data(hd_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_yzzzz_xz = cbuffer.data(hd_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_yzzzz_yy = cbuffer.data(hd_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_yzzzz_yz = cbuffer.data(hd_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_yzzzz_zz = cbuffer.data(hd_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzzzz_xx, g_x_x_yzzzz_xy, g_x_x_yzzzz_xz, g_x_x_yzzzz_yy, g_x_x_yzzzz_yz, g_x_x_yzzzz_zz, g_x_x_zzzz_xx, g_x_x_zzzz_xxy, g_x_x_zzzz_xy, g_x_x_zzzz_xyy, g_x_x_zzzz_xyz, g_x_x_zzzz_xz, g_x_x_zzzz_yy, g_x_x_zzzz_yyy, g_x_x_zzzz_yyz, g_x_x_zzzz_yz, g_x_x_zzzz_yzz, g_x_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzzzz_xx[k] = -g_x_x_zzzz_xx[k] * ab_y + g_x_x_zzzz_xxy[k];

                g_x_x_yzzzz_xy[k] = -g_x_x_zzzz_xy[k] * ab_y + g_x_x_zzzz_xyy[k];

                g_x_x_yzzzz_xz[k] = -g_x_x_zzzz_xz[k] * ab_y + g_x_x_zzzz_xyz[k];

                g_x_x_yzzzz_yy[k] = -g_x_x_zzzz_yy[k] * ab_y + g_x_x_zzzz_yyy[k];

                g_x_x_yzzzz_yz[k] = -g_x_x_zzzz_yz[k] * ab_y + g_x_x_zzzz_yyz[k];

                g_x_x_yzzzz_zz[k] = -g_x_x_zzzz_zz[k] * ab_y + g_x_x_zzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzzzz_xx = cbuffer.data(hd_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_zzzzz_xy = cbuffer.data(hd_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_zzzzz_xz = cbuffer.data(hd_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_zzzzz_yy = cbuffer.data(hd_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_zzzzz_yz = cbuffer.data(hd_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_zzzzz_zz = cbuffer.data(hd_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zzzz_xx, g_x_x_zzzz_xxz, g_x_x_zzzz_xy, g_x_x_zzzz_xyz, g_x_x_zzzz_xz, g_x_x_zzzz_xzz, g_x_x_zzzz_yy, g_x_x_zzzz_yyz, g_x_x_zzzz_yz, g_x_x_zzzz_yzz, g_x_x_zzzz_zz, g_x_x_zzzz_zzz, g_x_x_zzzzz_xx, g_x_x_zzzzz_xy, g_x_x_zzzzz_xz, g_x_x_zzzzz_yy, g_x_x_zzzzz_yz, g_x_x_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzzzz_xx[k] = -g_x_x_zzzz_xx[k] * ab_z + g_x_x_zzzz_xxz[k];

                g_x_x_zzzzz_xy[k] = -g_x_x_zzzz_xy[k] * ab_z + g_x_x_zzzz_xyz[k];

                g_x_x_zzzzz_xz[k] = -g_x_x_zzzz_xz[k] * ab_z + g_x_x_zzzz_xzz[k];

                g_x_x_zzzzz_yy[k] = -g_x_x_zzzz_yy[k] * ab_z + g_x_x_zzzz_yyz[k];

                g_x_x_zzzzz_yz[k] = -g_x_x_zzzz_yz[k] * ab_z + g_x_x_zzzz_yzz[k];

                g_x_x_zzzzz_zz[k] = -g_x_x_zzzz_zz[k] * ab_z + g_x_x_zzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxx_xx = cbuffer.data(hd_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xxxxx_xy = cbuffer.data(hd_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xxxxx_xz = cbuffer.data(hd_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xxxxx_yy = cbuffer.data(hd_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xxxxx_yz = cbuffer.data(hd_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xxxxx_zz = cbuffer.data(hd_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_xx, g_0_y_xxxx_xy, g_0_y_xxxx_xz, g_0_y_xxxx_yy, g_0_y_xxxx_yz, g_0_y_xxxx_zz, g_x_y_xxxx_xx, g_x_y_xxxx_xxx, g_x_y_xxxx_xxy, g_x_y_xxxx_xxz, g_x_y_xxxx_xy, g_x_y_xxxx_xyy, g_x_y_xxxx_xyz, g_x_y_xxxx_xz, g_x_y_xxxx_xzz, g_x_y_xxxx_yy, g_x_y_xxxx_yz, g_x_y_xxxx_zz, g_x_y_xxxxx_xx, g_x_y_xxxxx_xy, g_x_y_xxxxx_xz, g_x_y_xxxxx_yy, g_x_y_xxxxx_yz, g_x_y_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxx_xx[k] = -g_0_y_xxxx_xx[k] - g_x_y_xxxx_xx[k] * ab_x + g_x_y_xxxx_xxx[k];

                g_x_y_xxxxx_xy[k] = -g_0_y_xxxx_xy[k] - g_x_y_xxxx_xy[k] * ab_x + g_x_y_xxxx_xxy[k];

                g_x_y_xxxxx_xz[k] = -g_0_y_xxxx_xz[k] - g_x_y_xxxx_xz[k] * ab_x + g_x_y_xxxx_xxz[k];

                g_x_y_xxxxx_yy[k] = -g_0_y_xxxx_yy[k] - g_x_y_xxxx_yy[k] * ab_x + g_x_y_xxxx_xyy[k];

                g_x_y_xxxxx_yz[k] = -g_0_y_xxxx_yz[k] - g_x_y_xxxx_yz[k] * ab_x + g_x_y_xxxx_xyz[k];

                g_x_y_xxxxx_zz[k] = -g_0_y_xxxx_zz[k] - g_x_y_xxxx_zz[k] * ab_x + g_x_y_xxxx_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxy_xx = cbuffer.data(hd_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xxxxy_xy = cbuffer.data(hd_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xxxxy_xz = cbuffer.data(hd_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xxxxy_yy = cbuffer.data(hd_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xxxxy_yz = cbuffer.data(hd_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xxxxy_zz = cbuffer.data(hd_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_xx, g_0_y_xxxy_xy, g_0_y_xxxy_xz, g_0_y_xxxy_yy, g_0_y_xxxy_yz, g_0_y_xxxy_zz, g_x_y_xxxxy_xx, g_x_y_xxxxy_xy, g_x_y_xxxxy_xz, g_x_y_xxxxy_yy, g_x_y_xxxxy_yz, g_x_y_xxxxy_zz, g_x_y_xxxy_xx, g_x_y_xxxy_xxx, g_x_y_xxxy_xxy, g_x_y_xxxy_xxz, g_x_y_xxxy_xy, g_x_y_xxxy_xyy, g_x_y_xxxy_xyz, g_x_y_xxxy_xz, g_x_y_xxxy_xzz, g_x_y_xxxy_yy, g_x_y_xxxy_yz, g_x_y_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxy_xx[k] = -g_0_y_xxxy_xx[k] - g_x_y_xxxy_xx[k] * ab_x + g_x_y_xxxy_xxx[k];

                g_x_y_xxxxy_xy[k] = -g_0_y_xxxy_xy[k] - g_x_y_xxxy_xy[k] * ab_x + g_x_y_xxxy_xxy[k];

                g_x_y_xxxxy_xz[k] = -g_0_y_xxxy_xz[k] - g_x_y_xxxy_xz[k] * ab_x + g_x_y_xxxy_xxz[k];

                g_x_y_xxxxy_yy[k] = -g_0_y_xxxy_yy[k] - g_x_y_xxxy_yy[k] * ab_x + g_x_y_xxxy_xyy[k];

                g_x_y_xxxxy_yz[k] = -g_0_y_xxxy_yz[k] - g_x_y_xxxy_yz[k] * ab_x + g_x_y_xxxy_xyz[k];

                g_x_y_xxxxy_zz[k] = -g_0_y_xxxy_zz[k] - g_x_y_xxxy_zz[k] * ab_x + g_x_y_xxxy_xzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxxz_xx = cbuffer.data(hd_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xxxxz_xy = cbuffer.data(hd_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_xxxxz_xz = cbuffer.data(hd_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xxxxz_yy = cbuffer.data(hd_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xxxxz_yz = cbuffer.data(hd_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xxxxz_zz = cbuffer.data(hd_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxx_xx, g_x_y_xxxx_xxz, g_x_y_xxxx_xy, g_x_y_xxxx_xyz, g_x_y_xxxx_xz, g_x_y_xxxx_xzz, g_x_y_xxxx_yy, g_x_y_xxxx_yyz, g_x_y_xxxx_yz, g_x_y_xxxx_yzz, g_x_y_xxxx_zz, g_x_y_xxxx_zzz, g_x_y_xxxxz_xx, g_x_y_xxxxz_xy, g_x_y_xxxxz_xz, g_x_y_xxxxz_yy, g_x_y_xxxxz_yz, g_x_y_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxxz_xx[k] = -g_x_y_xxxx_xx[k] * ab_z + g_x_y_xxxx_xxz[k];

                g_x_y_xxxxz_xy[k] = -g_x_y_xxxx_xy[k] * ab_z + g_x_y_xxxx_xyz[k];

                g_x_y_xxxxz_xz[k] = -g_x_y_xxxx_xz[k] * ab_z + g_x_y_xxxx_xzz[k];

                g_x_y_xxxxz_yy[k] = -g_x_y_xxxx_yy[k] * ab_z + g_x_y_xxxx_yyz[k];

                g_x_y_xxxxz_yz[k] = -g_x_y_xxxx_yz[k] * ab_z + g_x_y_xxxx_yzz[k];

                g_x_y_xxxxz_zz[k] = -g_x_y_xxxx_zz[k] * ab_z + g_x_y_xxxx_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxyy_xx = cbuffer.data(hd_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xxxyy_xy = cbuffer.data(hd_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xxxyy_xz = cbuffer.data(hd_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_xxxyy_yy = cbuffer.data(hd_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xxxyy_yz = cbuffer.data(hd_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xxxyy_zz = cbuffer.data(hd_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_xx, g_0_y_xxyy_xy, g_0_y_xxyy_xz, g_0_y_xxyy_yy, g_0_y_xxyy_yz, g_0_y_xxyy_zz, g_x_y_xxxyy_xx, g_x_y_xxxyy_xy, g_x_y_xxxyy_xz, g_x_y_xxxyy_yy, g_x_y_xxxyy_yz, g_x_y_xxxyy_zz, g_x_y_xxyy_xx, g_x_y_xxyy_xxx, g_x_y_xxyy_xxy, g_x_y_xxyy_xxz, g_x_y_xxyy_xy, g_x_y_xxyy_xyy, g_x_y_xxyy_xyz, g_x_y_xxyy_xz, g_x_y_xxyy_xzz, g_x_y_xxyy_yy, g_x_y_xxyy_yz, g_x_y_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxyy_xx[k] = -g_0_y_xxyy_xx[k] - g_x_y_xxyy_xx[k] * ab_x + g_x_y_xxyy_xxx[k];

                g_x_y_xxxyy_xy[k] = -g_0_y_xxyy_xy[k] - g_x_y_xxyy_xy[k] * ab_x + g_x_y_xxyy_xxy[k];

                g_x_y_xxxyy_xz[k] = -g_0_y_xxyy_xz[k] - g_x_y_xxyy_xz[k] * ab_x + g_x_y_xxyy_xxz[k];

                g_x_y_xxxyy_yy[k] = -g_0_y_xxyy_yy[k] - g_x_y_xxyy_yy[k] * ab_x + g_x_y_xxyy_xyy[k];

                g_x_y_xxxyy_yz[k] = -g_0_y_xxyy_yz[k] - g_x_y_xxyy_yz[k] * ab_x + g_x_y_xxyy_xyz[k];

                g_x_y_xxxyy_zz[k] = -g_0_y_xxyy_zz[k] - g_x_y_xxyy_zz[k] * ab_x + g_x_y_xxyy_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxyz_xx = cbuffer.data(hd_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xxxyz_xy = cbuffer.data(hd_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xxxyz_xz = cbuffer.data(hd_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xxxyz_yy = cbuffer.data(hd_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xxxyz_yz = cbuffer.data(hd_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xxxyz_zz = cbuffer.data(hd_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxy_xx, g_x_y_xxxy_xxz, g_x_y_xxxy_xy, g_x_y_xxxy_xyz, g_x_y_xxxy_xz, g_x_y_xxxy_xzz, g_x_y_xxxy_yy, g_x_y_xxxy_yyz, g_x_y_xxxy_yz, g_x_y_xxxy_yzz, g_x_y_xxxy_zz, g_x_y_xxxy_zzz, g_x_y_xxxyz_xx, g_x_y_xxxyz_xy, g_x_y_xxxyz_xz, g_x_y_xxxyz_yy, g_x_y_xxxyz_yz, g_x_y_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxyz_xx[k] = -g_x_y_xxxy_xx[k] * ab_z + g_x_y_xxxy_xxz[k];

                g_x_y_xxxyz_xy[k] = -g_x_y_xxxy_xy[k] * ab_z + g_x_y_xxxy_xyz[k];

                g_x_y_xxxyz_xz[k] = -g_x_y_xxxy_xz[k] * ab_z + g_x_y_xxxy_xzz[k];

                g_x_y_xxxyz_yy[k] = -g_x_y_xxxy_yy[k] * ab_z + g_x_y_xxxy_yyz[k];

                g_x_y_xxxyz_yz[k] = -g_x_y_xxxy_yz[k] * ab_z + g_x_y_xxxy_yzz[k];

                g_x_y_xxxyz_zz[k] = -g_x_y_xxxy_zz[k] * ab_z + g_x_y_xxxy_zzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxzz_xx = cbuffer.data(hd_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xxxzz_xy = cbuffer.data(hd_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xxxzz_xz = cbuffer.data(hd_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xxxzz_yy = cbuffer.data(hd_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_xxxzz_yz = cbuffer.data(hd_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_xxxzz_zz = cbuffer.data(hd_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxxz_xx, g_x_y_xxxz_xxz, g_x_y_xxxz_xy, g_x_y_xxxz_xyz, g_x_y_xxxz_xz, g_x_y_xxxz_xzz, g_x_y_xxxz_yy, g_x_y_xxxz_yyz, g_x_y_xxxz_yz, g_x_y_xxxz_yzz, g_x_y_xxxz_zz, g_x_y_xxxz_zzz, g_x_y_xxxzz_xx, g_x_y_xxxzz_xy, g_x_y_xxxzz_xz, g_x_y_xxxzz_yy, g_x_y_xxxzz_yz, g_x_y_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxzz_xx[k] = -g_x_y_xxxz_xx[k] * ab_z + g_x_y_xxxz_xxz[k];

                g_x_y_xxxzz_xy[k] = -g_x_y_xxxz_xy[k] * ab_z + g_x_y_xxxz_xyz[k];

                g_x_y_xxxzz_xz[k] = -g_x_y_xxxz_xz[k] * ab_z + g_x_y_xxxz_xzz[k];

                g_x_y_xxxzz_yy[k] = -g_x_y_xxxz_yy[k] * ab_z + g_x_y_xxxz_yyz[k];

                g_x_y_xxxzz_yz[k] = -g_x_y_xxxz_yz[k] * ab_z + g_x_y_xxxz_yzz[k];

                g_x_y_xxxzz_zz[k] = -g_x_y_xxxz_zz[k] * ab_z + g_x_y_xxxz_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyyy_xx = cbuffer.data(hd_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_xxyyy_xy = cbuffer.data(hd_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_xxyyy_xz = cbuffer.data(hd_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_xxyyy_yy = cbuffer.data(hd_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_xxyyy_yz = cbuffer.data(hd_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_xxyyy_zz = cbuffer.data(hd_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_xx, g_0_y_xyyy_xy, g_0_y_xyyy_xz, g_0_y_xyyy_yy, g_0_y_xyyy_yz, g_0_y_xyyy_zz, g_x_y_xxyyy_xx, g_x_y_xxyyy_xy, g_x_y_xxyyy_xz, g_x_y_xxyyy_yy, g_x_y_xxyyy_yz, g_x_y_xxyyy_zz, g_x_y_xyyy_xx, g_x_y_xyyy_xxx, g_x_y_xyyy_xxy, g_x_y_xyyy_xxz, g_x_y_xyyy_xy, g_x_y_xyyy_xyy, g_x_y_xyyy_xyz, g_x_y_xyyy_xz, g_x_y_xyyy_xzz, g_x_y_xyyy_yy, g_x_y_xyyy_yz, g_x_y_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyyy_xx[k] = -g_0_y_xyyy_xx[k] - g_x_y_xyyy_xx[k] * ab_x + g_x_y_xyyy_xxx[k];

                g_x_y_xxyyy_xy[k] = -g_0_y_xyyy_xy[k] - g_x_y_xyyy_xy[k] * ab_x + g_x_y_xyyy_xxy[k];

                g_x_y_xxyyy_xz[k] = -g_0_y_xyyy_xz[k] - g_x_y_xyyy_xz[k] * ab_x + g_x_y_xyyy_xxz[k];

                g_x_y_xxyyy_yy[k] = -g_0_y_xyyy_yy[k] - g_x_y_xyyy_yy[k] * ab_x + g_x_y_xyyy_xyy[k];

                g_x_y_xxyyy_yz[k] = -g_0_y_xyyy_yz[k] - g_x_y_xyyy_yz[k] * ab_x + g_x_y_xyyy_xyz[k];

                g_x_y_xxyyy_zz[k] = -g_0_y_xyyy_zz[k] - g_x_y_xyyy_zz[k] * ab_x + g_x_y_xyyy_xzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyyz_xx = cbuffer.data(hd_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xxyyz_xy = cbuffer.data(hd_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xxyyz_xz = cbuffer.data(hd_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xxyyz_yy = cbuffer.data(hd_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xxyyz_yz = cbuffer.data(hd_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xxyyz_zz = cbuffer.data(hd_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxyy_xx, g_x_y_xxyy_xxz, g_x_y_xxyy_xy, g_x_y_xxyy_xyz, g_x_y_xxyy_xz, g_x_y_xxyy_xzz, g_x_y_xxyy_yy, g_x_y_xxyy_yyz, g_x_y_xxyy_yz, g_x_y_xxyy_yzz, g_x_y_xxyy_zz, g_x_y_xxyy_zzz, g_x_y_xxyyz_xx, g_x_y_xxyyz_xy, g_x_y_xxyyz_xz, g_x_y_xxyyz_yy, g_x_y_xxyyz_yz, g_x_y_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyyz_xx[k] = -g_x_y_xxyy_xx[k] * ab_z + g_x_y_xxyy_xxz[k];

                g_x_y_xxyyz_xy[k] = -g_x_y_xxyy_xy[k] * ab_z + g_x_y_xxyy_xyz[k];

                g_x_y_xxyyz_xz[k] = -g_x_y_xxyy_xz[k] * ab_z + g_x_y_xxyy_xzz[k];

                g_x_y_xxyyz_yy[k] = -g_x_y_xxyy_yy[k] * ab_z + g_x_y_xxyy_yyz[k];

                g_x_y_xxyyz_yz[k] = -g_x_y_xxyy_yz[k] * ab_z + g_x_y_xxyy_yzz[k];

                g_x_y_xxyyz_zz[k] = -g_x_y_xxyy_zz[k] * ab_z + g_x_y_xxyy_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyzz_xx = cbuffer.data(hd_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xxyzz_xy = cbuffer.data(hd_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xxyzz_xz = cbuffer.data(hd_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xxyzz_yy = cbuffer.data(hd_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xxyzz_yz = cbuffer.data(hd_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xxyzz_zz = cbuffer.data(hd_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxyz_xx, g_x_y_xxyz_xxz, g_x_y_xxyz_xy, g_x_y_xxyz_xyz, g_x_y_xxyz_xz, g_x_y_xxyz_xzz, g_x_y_xxyz_yy, g_x_y_xxyz_yyz, g_x_y_xxyz_yz, g_x_y_xxyz_yzz, g_x_y_xxyz_zz, g_x_y_xxyz_zzz, g_x_y_xxyzz_xx, g_x_y_xxyzz_xy, g_x_y_xxyzz_xz, g_x_y_xxyzz_yy, g_x_y_xxyzz_yz, g_x_y_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyzz_xx[k] = -g_x_y_xxyz_xx[k] * ab_z + g_x_y_xxyz_xxz[k];

                g_x_y_xxyzz_xy[k] = -g_x_y_xxyz_xy[k] * ab_z + g_x_y_xxyz_xyz[k];

                g_x_y_xxyzz_xz[k] = -g_x_y_xxyz_xz[k] * ab_z + g_x_y_xxyz_xzz[k];

                g_x_y_xxyzz_yy[k] = -g_x_y_xxyz_yy[k] * ab_z + g_x_y_xxyz_yyz[k];

                g_x_y_xxyzz_yz[k] = -g_x_y_xxyz_yz[k] * ab_z + g_x_y_xxyz_yzz[k];

                g_x_y_xxyzz_zz[k] = -g_x_y_xxyz_zz[k] * ab_z + g_x_y_xxyz_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxzzz_xx = cbuffer.data(hd_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xxzzz_xy = cbuffer.data(hd_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xxzzz_xz = cbuffer.data(hd_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xxzzz_yy = cbuffer.data(hd_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xxzzz_yz = cbuffer.data(hd_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xxzzz_zz = cbuffer.data(hd_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxzz_xx, g_x_y_xxzz_xxz, g_x_y_xxzz_xy, g_x_y_xxzz_xyz, g_x_y_xxzz_xz, g_x_y_xxzz_xzz, g_x_y_xxzz_yy, g_x_y_xxzz_yyz, g_x_y_xxzz_yz, g_x_y_xxzz_yzz, g_x_y_xxzz_zz, g_x_y_xxzz_zzz, g_x_y_xxzzz_xx, g_x_y_xxzzz_xy, g_x_y_xxzzz_xz, g_x_y_xxzzz_yy, g_x_y_xxzzz_yz, g_x_y_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxzzz_xx[k] = -g_x_y_xxzz_xx[k] * ab_z + g_x_y_xxzz_xxz[k];

                g_x_y_xxzzz_xy[k] = -g_x_y_xxzz_xy[k] * ab_z + g_x_y_xxzz_xyz[k];

                g_x_y_xxzzz_xz[k] = -g_x_y_xxzz_xz[k] * ab_z + g_x_y_xxzz_xzz[k];

                g_x_y_xxzzz_yy[k] = -g_x_y_xxzz_yy[k] * ab_z + g_x_y_xxzz_yyz[k];

                g_x_y_xxzzz_yz[k] = -g_x_y_xxzz_yz[k] * ab_z + g_x_y_xxzz_yzz[k];

                g_x_y_xxzzz_zz[k] = -g_x_y_xxzz_zz[k] * ab_z + g_x_y_xxzz_zzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyyy_xx = cbuffer.data(hd_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xyyyy_xy = cbuffer.data(hd_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xyyyy_xz = cbuffer.data(hd_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_xyyyy_yy = cbuffer.data(hd_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_xyyyy_yz = cbuffer.data(hd_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_xyyyy_zz = cbuffer.data(hd_geom_11_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xx, g_0_y_yyyy_xy, g_0_y_yyyy_xz, g_0_y_yyyy_yy, g_0_y_yyyy_yz, g_0_y_yyyy_zz, g_x_y_xyyyy_xx, g_x_y_xyyyy_xy, g_x_y_xyyyy_xz, g_x_y_xyyyy_yy, g_x_y_xyyyy_yz, g_x_y_xyyyy_zz, g_x_y_yyyy_xx, g_x_y_yyyy_xxx, g_x_y_yyyy_xxy, g_x_y_yyyy_xxz, g_x_y_yyyy_xy, g_x_y_yyyy_xyy, g_x_y_yyyy_xyz, g_x_y_yyyy_xz, g_x_y_yyyy_xzz, g_x_y_yyyy_yy, g_x_y_yyyy_yz, g_x_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyyy_xx[k] = -g_0_y_yyyy_xx[k] - g_x_y_yyyy_xx[k] * ab_x + g_x_y_yyyy_xxx[k];

                g_x_y_xyyyy_xy[k] = -g_0_y_yyyy_xy[k] - g_x_y_yyyy_xy[k] * ab_x + g_x_y_yyyy_xxy[k];

                g_x_y_xyyyy_xz[k] = -g_0_y_yyyy_xz[k] - g_x_y_yyyy_xz[k] * ab_x + g_x_y_yyyy_xxz[k];

                g_x_y_xyyyy_yy[k] = -g_0_y_yyyy_yy[k] - g_x_y_yyyy_yy[k] * ab_x + g_x_y_yyyy_xyy[k];

                g_x_y_xyyyy_yz[k] = -g_0_y_yyyy_yz[k] - g_x_y_yyyy_yz[k] * ab_x + g_x_y_yyyy_xyz[k];

                g_x_y_xyyyy_zz[k] = -g_0_y_yyyy_zz[k] - g_x_y_yyyy_zz[k] * ab_x + g_x_y_yyyy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyyz_xx = cbuffer.data(hd_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_xyyyz_xy = cbuffer.data(hd_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_xyyyz_xz = cbuffer.data(hd_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_xyyyz_yy = cbuffer.data(hd_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_xyyyz_yz = cbuffer.data(hd_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_xyyyz_zz = cbuffer.data(hd_geom_11_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyyy_xx, g_x_y_xyyy_xxz, g_x_y_xyyy_xy, g_x_y_xyyy_xyz, g_x_y_xyyy_xz, g_x_y_xyyy_xzz, g_x_y_xyyy_yy, g_x_y_xyyy_yyz, g_x_y_xyyy_yz, g_x_y_xyyy_yzz, g_x_y_xyyy_zz, g_x_y_xyyy_zzz, g_x_y_xyyyz_xx, g_x_y_xyyyz_xy, g_x_y_xyyyz_xz, g_x_y_xyyyz_yy, g_x_y_xyyyz_yz, g_x_y_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyyz_xx[k] = -g_x_y_xyyy_xx[k] * ab_z + g_x_y_xyyy_xxz[k];

                g_x_y_xyyyz_xy[k] = -g_x_y_xyyy_xy[k] * ab_z + g_x_y_xyyy_xyz[k];

                g_x_y_xyyyz_xz[k] = -g_x_y_xyyy_xz[k] * ab_z + g_x_y_xyyy_xzz[k];

                g_x_y_xyyyz_yy[k] = -g_x_y_xyyy_yy[k] * ab_z + g_x_y_xyyy_yyz[k];

                g_x_y_xyyyz_yz[k] = -g_x_y_xyyy_yz[k] * ab_z + g_x_y_xyyy_yzz[k];

                g_x_y_xyyyz_zz[k] = -g_x_y_xyyy_zz[k] * ab_z + g_x_y_xyyy_zzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyzz_xx = cbuffer.data(hd_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_xyyzz_xy = cbuffer.data(hd_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_xyyzz_xz = cbuffer.data(hd_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_xyyzz_yy = cbuffer.data(hd_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_xyyzz_yz = cbuffer.data(hd_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_xyyzz_zz = cbuffer.data(hd_geom_11_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyyz_xx, g_x_y_xyyz_xxz, g_x_y_xyyz_xy, g_x_y_xyyz_xyz, g_x_y_xyyz_xz, g_x_y_xyyz_xzz, g_x_y_xyyz_yy, g_x_y_xyyz_yyz, g_x_y_xyyz_yz, g_x_y_xyyz_yzz, g_x_y_xyyz_zz, g_x_y_xyyz_zzz, g_x_y_xyyzz_xx, g_x_y_xyyzz_xy, g_x_y_xyyzz_xz, g_x_y_xyyzz_yy, g_x_y_xyyzz_yz, g_x_y_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyzz_xx[k] = -g_x_y_xyyz_xx[k] * ab_z + g_x_y_xyyz_xxz[k];

                g_x_y_xyyzz_xy[k] = -g_x_y_xyyz_xy[k] * ab_z + g_x_y_xyyz_xyz[k];

                g_x_y_xyyzz_xz[k] = -g_x_y_xyyz_xz[k] * ab_z + g_x_y_xyyz_xzz[k];

                g_x_y_xyyzz_yy[k] = -g_x_y_xyyz_yy[k] * ab_z + g_x_y_xyyz_yyz[k];

                g_x_y_xyyzz_yz[k] = -g_x_y_xyyz_yz[k] * ab_z + g_x_y_xyyz_yzz[k];

                g_x_y_xyyzz_zz[k] = -g_x_y_xyyz_zz[k] * ab_z + g_x_y_xyyz_zzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyzzz_xx = cbuffer.data(hd_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_xyzzz_xy = cbuffer.data(hd_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_xyzzz_xz = cbuffer.data(hd_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_xyzzz_yy = cbuffer.data(hd_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_xyzzz_yz = cbuffer.data(hd_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_xyzzz_zz = cbuffer.data(hd_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyzz_xx, g_x_y_xyzz_xxz, g_x_y_xyzz_xy, g_x_y_xyzz_xyz, g_x_y_xyzz_xz, g_x_y_xyzz_xzz, g_x_y_xyzz_yy, g_x_y_xyzz_yyz, g_x_y_xyzz_yz, g_x_y_xyzz_yzz, g_x_y_xyzz_zz, g_x_y_xyzz_zzz, g_x_y_xyzzz_xx, g_x_y_xyzzz_xy, g_x_y_xyzzz_xz, g_x_y_xyzzz_yy, g_x_y_xyzzz_yz, g_x_y_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyzzz_xx[k] = -g_x_y_xyzz_xx[k] * ab_z + g_x_y_xyzz_xxz[k];

                g_x_y_xyzzz_xy[k] = -g_x_y_xyzz_xy[k] * ab_z + g_x_y_xyzz_xyz[k];

                g_x_y_xyzzz_xz[k] = -g_x_y_xyzz_xz[k] * ab_z + g_x_y_xyzz_xzz[k];

                g_x_y_xyzzz_yy[k] = -g_x_y_xyzz_yy[k] * ab_z + g_x_y_xyzz_yyz[k];

                g_x_y_xyzzz_yz[k] = -g_x_y_xyzz_yz[k] * ab_z + g_x_y_xyzz_yzz[k];

                g_x_y_xyzzz_zz[k] = -g_x_y_xyzz_zz[k] * ab_z + g_x_y_xyzz_zzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzzzz_xx = cbuffer.data(hd_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_xzzzz_xy = cbuffer.data(hd_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_xzzzz_xz = cbuffer.data(hd_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_xzzzz_yy = cbuffer.data(hd_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_xzzzz_yz = cbuffer.data(hd_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_xzzzz_zz = cbuffer.data(hd_geom_11_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xzzz_xx, g_x_y_xzzz_xxz, g_x_y_xzzz_xy, g_x_y_xzzz_xyz, g_x_y_xzzz_xz, g_x_y_xzzz_xzz, g_x_y_xzzz_yy, g_x_y_xzzz_yyz, g_x_y_xzzz_yz, g_x_y_xzzz_yzz, g_x_y_xzzz_zz, g_x_y_xzzz_zzz, g_x_y_xzzzz_xx, g_x_y_xzzzz_xy, g_x_y_xzzzz_xz, g_x_y_xzzzz_yy, g_x_y_xzzzz_yz, g_x_y_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzzzz_xx[k] = -g_x_y_xzzz_xx[k] * ab_z + g_x_y_xzzz_xxz[k];

                g_x_y_xzzzz_xy[k] = -g_x_y_xzzz_xy[k] * ab_z + g_x_y_xzzz_xyz[k];

                g_x_y_xzzzz_xz[k] = -g_x_y_xzzz_xz[k] * ab_z + g_x_y_xzzz_xzz[k];

                g_x_y_xzzzz_yy[k] = -g_x_y_xzzz_yy[k] * ab_z + g_x_y_xzzz_yyz[k];

                g_x_y_xzzzz_yz[k] = -g_x_y_xzzz_yz[k] * ab_z + g_x_y_xzzz_yzz[k];

                g_x_y_xzzzz_zz[k] = -g_x_y_xzzz_zz[k] * ab_z + g_x_y_xzzz_zzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyyy_xx = cbuffer.data(hd_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_yyyyy_xy = cbuffer.data(hd_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_yyyyy_xz = cbuffer.data(hd_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_yyyyy_yy = cbuffer.data(hd_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_yyyyy_yz = cbuffer.data(hd_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_yyyyy_zz = cbuffer.data(hd_geom_11_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_xx, g_x_0_yyyy_xy, g_x_0_yyyy_xz, g_x_0_yyyy_yy, g_x_0_yyyy_yz, g_x_0_yyyy_zz, g_x_y_yyyy_xx, g_x_y_yyyy_xxy, g_x_y_yyyy_xy, g_x_y_yyyy_xyy, g_x_y_yyyy_xyz, g_x_y_yyyy_xz, g_x_y_yyyy_yy, g_x_y_yyyy_yyy, g_x_y_yyyy_yyz, g_x_y_yyyy_yz, g_x_y_yyyy_yzz, g_x_y_yyyy_zz, g_x_y_yyyyy_xx, g_x_y_yyyyy_xy, g_x_y_yyyyy_xz, g_x_y_yyyyy_yy, g_x_y_yyyyy_yz, g_x_y_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyyy_xx[k] = g_x_0_yyyy_xx[k] - g_x_y_yyyy_xx[k] * ab_y + g_x_y_yyyy_xxy[k];

                g_x_y_yyyyy_xy[k] = g_x_0_yyyy_xy[k] - g_x_y_yyyy_xy[k] * ab_y + g_x_y_yyyy_xyy[k];

                g_x_y_yyyyy_xz[k] = g_x_0_yyyy_xz[k] - g_x_y_yyyy_xz[k] * ab_y + g_x_y_yyyy_xyz[k];

                g_x_y_yyyyy_yy[k] = g_x_0_yyyy_yy[k] - g_x_y_yyyy_yy[k] * ab_y + g_x_y_yyyy_yyy[k];

                g_x_y_yyyyy_yz[k] = g_x_0_yyyy_yz[k] - g_x_y_yyyy_yz[k] * ab_y + g_x_y_yyyy_yyz[k];

                g_x_y_yyyyy_zz[k] = g_x_0_yyyy_zz[k] - g_x_y_yyyy_zz[k] * ab_y + g_x_y_yyyy_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyyz_xx = cbuffer.data(hd_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_yyyyz_xy = cbuffer.data(hd_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_y_yyyyz_xz = cbuffer.data(hd_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_y_yyyyz_yy = cbuffer.data(hd_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_yyyyz_yz = cbuffer.data(hd_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_yyyyz_zz = cbuffer.data(hd_geom_11_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyyy_xx, g_x_y_yyyy_xxz, g_x_y_yyyy_xy, g_x_y_yyyy_xyz, g_x_y_yyyy_xz, g_x_y_yyyy_xzz, g_x_y_yyyy_yy, g_x_y_yyyy_yyz, g_x_y_yyyy_yz, g_x_y_yyyy_yzz, g_x_y_yyyy_zz, g_x_y_yyyy_zzz, g_x_y_yyyyz_xx, g_x_y_yyyyz_xy, g_x_y_yyyyz_xz, g_x_y_yyyyz_yy, g_x_y_yyyyz_yz, g_x_y_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyyz_xx[k] = -g_x_y_yyyy_xx[k] * ab_z + g_x_y_yyyy_xxz[k];

                g_x_y_yyyyz_xy[k] = -g_x_y_yyyy_xy[k] * ab_z + g_x_y_yyyy_xyz[k];

                g_x_y_yyyyz_xz[k] = -g_x_y_yyyy_xz[k] * ab_z + g_x_y_yyyy_xzz[k];

                g_x_y_yyyyz_yy[k] = -g_x_y_yyyy_yy[k] * ab_z + g_x_y_yyyy_yyz[k];

                g_x_y_yyyyz_yz[k] = -g_x_y_yyyy_yz[k] * ab_z + g_x_y_yyyy_yzz[k];

                g_x_y_yyyyz_zz[k] = -g_x_y_yyyy_zz[k] * ab_z + g_x_y_yyyy_zzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyzz_xx = cbuffer.data(hd_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_yyyzz_xy = cbuffer.data(hd_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_yyyzz_xz = cbuffer.data(hd_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_y_yyyzz_yy = cbuffer.data(hd_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_yyyzz_yz = cbuffer.data(hd_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_yyyzz_zz = cbuffer.data(hd_geom_11_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyyz_xx, g_x_y_yyyz_xxz, g_x_y_yyyz_xy, g_x_y_yyyz_xyz, g_x_y_yyyz_xz, g_x_y_yyyz_xzz, g_x_y_yyyz_yy, g_x_y_yyyz_yyz, g_x_y_yyyz_yz, g_x_y_yyyz_yzz, g_x_y_yyyz_zz, g_x_y_yyyz_zzz, g_x_y_yyyzz_xx, g_x_y_yyyzz_xy, g_x_y_yyyzz_xz, g_x_y_yyyzz_yy, g_x_y_yyyzz_yz, g_x_y_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyzz_xx[k] = -g_x_y_yyyz_xx[k] * ab_z + g_x_y_yyyz_xxz[k];

                g_x_y_yyyzz_xy[k] = -g_x_y_yyyz_xy[k] * ab_z + g_x_y_yyyz_xyz[k];

                g_x_y_yyyzz_xz[k] = -g_x_y_yyyz_xz[k] * ab_z + g_x_y_yyyz_xzz[k];

                g_x_y_yyyzz_yy[k] = -g_x_y_yyyz_yy[k] * ab_z + g_x_y_yyyz_yyz[k];

                g_x_y_yyyzz_yz[k] = -g_x_y_yyyz_yz[k] * ab_z + g_x_y_yyyz_yzz[k];

                g_x_y_yyyzz_zz[k] = -g_x_y_yyyz_zz[k] * ab_z + g_x_y_yyyz_zzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyzzz_xx = cbuffer.data(hd_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_yyzzz_xy = cbuffer.data(hd_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_yyzzz_xz = cbuffer.data(hd_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_yyzzz_yy = cbuffer.data(hd_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_yyzzz_yz = cbuffer.data(hd_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_yyzzz_zz = cbuffer.data(hd_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyzz_xx, g_x_y_yyzz_xxz, g_x_y_yyzz_xy, g_x_y_yyzz_xyz, g_x_y_yyzz_xz, g_x_y_yyzz_xzz, g_x_y_yyzz_yy, g_x_y_yyzz_yyz, g_x_y_yyzz_yz, g_x_y_yyzz_yzz, g_x_y_yyzz_zz, g_x_y_yyzz_zzz, g_x_y_yyzzz_xx, g_x_y_yyzzz_xy, g_x_y_yyzzz_xz, g_x_y_yyzzz_yy, g_x_y_yyzzz_yz, g_x_y_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyzzz_xx[k] = -g_x_y_yyzz_xx[k] * ab_z + g_x_y_yyzz_xxz[k];

                g_x_y_yyzzz_xy[k] = -g_x_y_yyzz_xy[k] * ab_z + g_x_y_yyzz_xyz[k];

                g_x_y_yyzzz_xz[k] = -g_x_y_yyzz_xz[k] * ab_z + g_x_y_yyzz_xzz[k];

                g_x_y_yyzzz_yy[k] = -g_x_y_yyzz_yy[k] * ab_z + g_x_y_yyzz_yyz[k];

                g_x_y_yyzzz_yz[k] = -g_x_y_yyzz_yz[k] * ab_z + g_x_y_yyzz_yzz[k];

                g_x_y_yyzzz_zz[k] = -g_x_y_yyzz_zz[k] * ab_z + g_x_y_yyzz_zzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzzzz_xx = cbuffer.data(hd_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_yzzzz_xy = cbuffer.data(hd_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_yzzzz_xz = cbuffer.data(hd_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_yzzzz_yy = cbuffer.data(hd_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_yzzzz_yz = cbuffer.data(hd_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_yzzzz_zz = cbuffer.data(hd_geom_11_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yzzz_xx, g_x_y_yzzz_xxz, g_x_y_yzzz_xy, g_x_y_yzzz_xyz, g_x_y_yzzz_xz, g_x_y_yzzz_xzz, g_x_y_yzzz_yy, g_x_y_yzzz_yyz, g_x_y_yzzz_yz, g_x_y_yzzz_yzz, g_x_y_yzzz_zz, g_x_y_yzzz_zzz, g_x_y_yzzzz_xx, g_x_y_yzzzz_xy, g_x_y_yzzzz_xz, g_x_y_yzzzz_yy, g_x_y_yzzzz_yz, g_x_y_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzzzz_xx[k] = -g_x_y_yzzz_xx[k] * ab_z + g_x_y_yzzz_xxz[k];

                g_x_y_yzzzz_xy[k] = -g_x_y_yzzz_xy[k] * ab_z + g_x_y_yzzz_xyz[k];

                g_x_y_yzzzz_xz[k] = -g_x_y_yzzz_xz[k] * ab_z + g_x_y_yzzz_xzz[k];

                g_x_y_yzzzz_yy[k] = -g_x_y_yzzz_yy[k] * ab_z + g_x_y_yzzz_yyz[k];

                g_x_y_yzzzz_yz[k] = -g_x_y_yzzz_yz[k] * ab_z + g_x_y_yzzz_yzz[k];

                g_x_y_yzzzz_zz[k] = -g_x_y_yzzz_zz[k] * ab_z + g_x_y_yzzz_zzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzzzz_xx = cbuffer.data(hd_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_zzzzz_xy = cbuffer.data(hd_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_zzzzz_xz = cbuffer.data(hd_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_zzzzz_yy = cbuffer.data(hd_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_zzzzz_yz = cbuffer.data(hd_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_zzzzz_zz = cbuffer.data(hd_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zzzz_xx, g_x_y_zzzz_xxz, g_x_y_zzzz_xy, g_x_y_zzzz_xyz, g_x_y_zzzz_xz, g_x_y_zzzz_xzz, g_x_y_zzzz_yy, g_x_y_zzzz_yyz, g_x_y_zzzz_yz, g_x_y_zzzz_yzz, g_x_y_zzzz_zz, g_x_y_zzzz_zzz, g_x_y_zzzzz_xx, g_x_y_zzzzz_xy, g_x_y_zzzzz_xz, g_x_y_zzzzz_yy, g_x_y_zzzzz_yz, g_x_y_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzzzz_xx[k] = -g_x_y_zzzz_xx[k] * ab_z + g_x_y_zzzz_xxz[k];

                g_x_y_zzzzz_xy[k] = -g_x_y_zzzz_xy[k] * ab_z + g_x_y_zzzz_xyz[k];

                g_x_y_zzzzz_xz[k] = -g_x_y_zzzz_xz[k] * ab_z + g_x_y_zzzz_xzz[k];

                g_x_y_zzzzz_yy[k] = -g_x_y_zzzz_yy[k] * ab_z + g_x_y_zzzz_yyz[k];

                g_x_y_zzzzz_yz[k] = -g_x_y_zzzz_yz[k] * ab_z + g_x_y_zzzz_yzz[k];

                g_x_y_zzzzz_zz[k] = -g_x_y_zzzz_zz[k] * ab_z + g_x_y_zzzz_zzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxx_xx = cbuffer.data(hd_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_xxxxx_xy = cbuffer.data(hd_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_xxxxx_xz = cbuffer.data(hd_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_xxxxx_yy = cbuffer.data(hd_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_xxxxx_yz = cbuffer.data(hd_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_xxxxx_zz = cbuffer.data(hd_geom_11_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_xx, g_0_z_xxxx_xy, g_0_z_xxxx_xz, g_0_z_xxxx_yy, g_0_z_xxxx_yz, g_0_z_xxxx_zz, g_x_z_xxxx_xx, g_x_z_xxxx_xxx, g_x_z_xxxx_xxy, g_x_z_xxxx_xxz, g_x_z_xxxx_xy, g_x_z_xxxx_xyy, g_x_z_xxxx_xyz, g_x_z_xxxx_xz, g_x_z_xxxx_xzz, g_x_z_xxxx_yy, g_x_z_xxxx_yz, g_x_z_xxxx_zz, g_x_z_xxxxx_xx, g_x_z_xxxxx_xy, g_x_z_xxxxx_xz, g_x_z_xxxxx_yy, g_x_z_xxxxx_yz, g_x_z_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxx_xx[k] = -g_0_z_xxxx_xx[k] - g_x_z_xxxx_xx[k] * ab_x + g_x_z_xxxx_xxx[k];

                g_x_z_xxxxx_xy[k] = -g_0_z_xxxx_xy[k] - g_x_z_xxxx_xy[k] * ab_x + g_x_z_xxxx_xxy[k];

                g_x_z_xxxxx_xz[k] = -g_0_z_xxxx_xz[k] - g_x_z_xxxx_xz[k] * ab_x + g_x_z_xxxx_xxz[k];

                g_x_z_xxxxx_yy[k] = -g_0_z_xxxx_yy[k] - g_x_z_xxxx_yy[k] * ab_x + g_x_z_xxxx_xyy[k];

                g_x_z_xxxxx_yz[k] = -g_0_z_xxxx_yz[k] - g_x_z_xxxx_yz[k] * ab_x + g_x_z_xxxx_xyz[k];

                g_x_z_xxxxx_zz[k] = -g_0_z_xxxx_zz[k] - g_x_z_xxxx_zz[k] * ab_x + g_x_z_xxxx_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxy_xx = cbuffer.data(hd_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_xxxxy_xy = cbuffer.data(hd_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_xxxxy_xz = cbuffer.data(hd_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_xxxxy_yy = cbuffer.data(hd_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_xxxxy_yz = cbuffer.data(hd_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_xxxxy_zz = cbuffer.data(hd_geom_11_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxx_xx, g_x_z_xxxx_xxy, g_x_z_xxxx_xy, g_x_z_xxxx_xyy, g_x_z_xxxx_xyz, g_x_z_xxxx_xz, g_x_z_xxxx_yy, g_x_z_xxxx_yyy, g_x_z_xxxx_yyz, g_x_z_xxxx_yz, g_x_z_xxxx_yzz, g_x_z_xxxx_zz, g_x_z_xxxxy_xx, g_x_z_xxxxy_xy, g_x_z_xxxxy_xz, g_x_z_xxxxy_yy, g_x_z_xxxxy_yz, g_x_z_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxy_xx[k] = -g_x_z_xxxx_xx[k] * ab_y + g_x_z_xxxx_xxy[k];

                g_x_z_xxxxy_xy[k] = -g_x_z_xxxx_xy[k] * ab_y + g_x_z_xxxx_xyy[k];

                g_x_z_xxxxy_xz[k] = -g_x_z_xxxx_xz[k] * ab_y + g_x_z_xxxx_xyz[k];

                g_x_z_xxxxy_yy[k] = -g_x_z_xxxx_yy[k] * ab_y + g_x_z_xxxx_yyy[k];

                g_x_z_xxxxy_yz[k] = -g_x_z_xxxx_yz[k] * ab_y + g_x_z_xxxx_yyz[k];

                g_x_z_xxxxy_zz[k] = -g_x_z_xxxx_zz[k] * ab_y + g_x_z_xxxx_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxxz_xx = cbuffer.data(hd_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_xxxxz_xy = cbuffer.data(hd_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_xxxxz_xz = cbuffer.data(hd_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_xxxxz_yy = cbuffer.data(hd_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_xxxxz_yz = cbuffer.data(hd_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_xxxxz_zz = cbuffer.data(hd_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_xx, g_0_z_xxxz_xy, g_0_z_xxxz_xz, g_0_z_xxxz_yy, g_0_z_xxxz_yz, g_0_z_xxxz_zz, g_x_z_xxxxz_xx, g_x_z_xxxxz_xy, g_x_z_xxxxz_xz, g_x_z_xxxxz_yy, g_x_z_xxxxz_yz, g_x_z_xxxxz_zz, g_x_z_xxxz_xx, g_x_z_xxxz_xxx, g_x_z_xxxz_xxy, g_x_z_xxxz_xxz, g_x_z_xxxz_xy, g_x_z_xxxz_xyy, g_x_z_xxxz_xyz, g_x_z_xxxz_xz, g_x_z_xxxz_xzz, g_x_z_xxxz_yy, g_x_z_xxxz_yz, g_x_z_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxxz_xx[k] = -g_0_z_xxxz_xx[k] - g_x_z_xxxz_xx[k] * ab_x + g_x_z_xxxz_xxx[k];

                g_x_z_xxxxz_xy[k] = -g_0_z_xxxz_xy[k] - g_x_z_xxxz_xy[k] * ab_x + g_x_z_xxxz_xxy[k];

                g_x_z_xxxxz_xz[k] = -g_0_z_xxxz_xz[k] - g_x_z_xxxz_xz[k] * ab_x + g_x_z_xxxz_xxz[k];

                g_x_z_xxxxz_yy[k] = -g_0_z_xxxz_yy[k] - g_x_z_xxxz_yy[k] * ab_x + g_x_z_xxxz_xyy[k];

                g_x_z_xxxxz_yz[k] = -g_0_z_xxxz_yz[k] - g_x_z_xxxz_yz[k] * ab_x + g_x_z_xxxz_xyz[k];

                g_x_z_xxxxz_zz[k] = -g_0_z_xxxz_zz[k] - g_x_z_xxxz_zz[k] * ab_x + g_x_z_xxxz_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxyy_xx = cbuffer.data(hd_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_z_xxxyy_xy = cbuffer.data(hd_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_z_xxxyy_xz = cbuffer.data(hd_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_z_xxxyy_yy = cbuffer.data(hd_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_z_xxxyy_yz = cbuffer.data(hd_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_z_xxxyy_zz = cbuffer.data(hd_geom_11_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxy_xx, g_x_z_xxxy_xxy, g_x_z_xxxy_xy, g_x_z_xxxy_xyy, g_x_z_xxxy_xyz, g_x_z_xxxy_xz, g_x_z_xxxy_yy, g_x_z_xxxy_yyy, g_x_z_xxxy_yyz, g_x_z_xxxy_yz, g_x_z_xxxy_yzz, g_x_z_xxxy_zz, g_x_z_xxxyy_xx, g_x_z_xxxyy_xy, g_x_z_xxxyy_xz, g_x_z_xxxyy_yy, g_x_z_xxxyy_yz, g_x_z_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxyy_xx[k] = -g_x_z_xxxy_xx[k] * ab_y + g_x_z_xxxy_xxy[k];

                g_x_z_xxxyy_xy[k] = -g_x_z_xxxy_xy[k] * ab_y + g_x_z_xxxy_xyy[k];

                g_x_z_xxxyy_xz[k] = -g_x_z_xxxy_xz[k] * ab_y + g_x_z_xxxy_xyz[k];

                g_x_z_xxxyy_yy[k] = -g_x_z_xxxy_yy[k] * ab_y + g_x_z_xxxy_yyy[k];

                g_x_z_xxxyy_yz[k] = -g_x_z_xxxy_yz[k] * ab_y + g_x_z_xxxy_yyz[k];

                g_x_z_xxxyy_zz[k] = -g_x_z_xxxy_zz[k] * ab_y + g_x_z_xxxy_yzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxyz_xx = cbuffer.data(hd_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_z_xxxyz_xy = cbuffer.data(hd_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_z_xxxyz_xz = cbuffer.data(hd_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_z_xxxyz_yy = cbuffer.data(hd_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_z_xxxyz_yz = cbuffer.data(hd_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_z_xxxyz_zz = cbuffer.data(hd_geom_11_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxxyz_xx, g_x_z_xxxyz_xy, g_x_z_xxxyz_xz, g_x_z_xxxyz_yy, g_x_z_xxxyz_yz, g_x_z_xxxyz_zz, g_x_z_xxxz_xx, g_x_z_xxxz_xxy, g_x_z_xxxz_xy, g_x_z_xxxz_xyy, g_x_z_xxxz_xyz, g_x_z_xxxz_xz, g_x_z_xxxz_yy, g_x_z_xxxz_yyy, g_x_z_xxxz_yyz, g_x_z_xxxz_yz, g_x_z_xxxz_yzz, g_x_z_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxyz_xx[k] = -g_x_z_xxxz_xx[k] * ab_y + g_x_z_xxxz_xxy[k];

                g_x_z_xxxyz_xy[k] = -g_x_z_xxxz_xy[k] * ab_y + g_x_z_xxxz_xyy[k];

                g_x_z_xxxyz_xz[k] = -g_x_z_xxxz_xz[k] * ab_y + g_x_z_xxxz_xyz[k];

                g_x_z_xxxyz_yy[k] = -g_x_z_xxxz_yy[k] * ab_y + g_x_z_xxxz_yyy[k];

                g_x_z_xxxyz_yz[k] = -g_x_z_xxxz_yz[k] * ab_y + g_x_z_xxxz_yyz[k];

                g_x_z_xxxyz_zz[k] = -g_x_z_xxxz_zz[k] * ab_y + g_x_z_xxxz_yzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxzz_xx = cbuffer.data(hd_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_z_xxxzz_xy = cbuffer.data(hd_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_z_xxxzz_xz = cbuffer.data(hd_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_z_xxxzz_yy = cbuffer.data(hd_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_z_xxxzz_yz = cbuffer.data(hd_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_z_xxxzz_zz = cbuffer.data(hd_geom_11_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_xx, g_0_z_xxzz_xy, g_0_z_xxzz_xz, g_0_z_xxzz_yy, g_0_z_xxzz_yz, g_0_z_xxzz_zz, g_x_z_xxxzz_xx, g_x_z_xxxzz_xy, g_x_z_xxxzz_xz, g_x_z_xxxzz_yy, g_x_z_xxxzz_yz, g_x_z_xxxzz_zz, g_x_z_xxzz_xx, g_x_z_xxzz_xxx, g_x_z_xxzz_xxy, g_x_z_xxzz_xxz, g_x_z_xxzz_xy, g_x_z_xxzz_xyy, g_x_z_xxzz_xyz, g_x_z_xxzz_xz, g_x_z_xxzz_xzz, g_x_z_xxzz_yy, g_x_z_xxzz_yz, g_x_z_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxzz_xx[k] = -g_0_z_xxzz_xx[k] - g_x_z_xxzz_xx[k] * ab_x + g_x_z_xxzz_xxx[k];

                g_x_z_xxxzz_xy[k] = -g_0_z_xxzz_xy[k] - g_x_z_xxzz_xy[k] * ab_x + g_x_z_xxzz_xxy[k];

                g_x_z_xxxzz_xz[k] = -g_0_z_xxzz_xz[k] - g_x_z_xxzz_xz[k] * ab_x + g_x_z_xxzz_xxz[k];

                g_x_z_xxxzz_yy[k] = -g_0_z_xxzz_yy[k] - g_x_z_xxzz_yy[k] * ab_x + g_x_z_xxzz_xyy[k];

                g_x_z_xxxzz_yz[k] = -g_0_z_xxzz_yz[k] - g_x_z_xxzz_yz[k] * ab_x + g_x_z_xxzz_xyz[k];

                g_x_z_xxxzz_zz[k] = -g_0_z_xxzz_zz[k] - g_x_z_xxzz_zz[k] * ab_x + g_x_z_xxzz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyyy_xx = cbuffer.data(hd_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_z_xxyyy_xy = cbuffer.data(hd_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_z_xxyyy_xz = cbuffer.data(hd_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_z_xxyyy_yy = cbuffer.data(hd_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_z_xxyyy_yz = cbuffer.data(hd_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_z_xxyyy_zz = cbuffer.data(hd_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyy_xx, g_x_z_xxyy_xxy, g_x_z_xxyy_xy, g_x_z_xxyy_xyy, g_x_z_xxyy_xyz, g_x_z_xxyy_xz, g_x_z_xxyy_yy, g_x_z_xxyy_yyy, g_x_z_xxyy_yyz, g_x_z_xxyy_yz, g_x_z_xxyy_yzz, g_x_z_xxyy_zz, g_x_z_xxyyy_xx, g_x_z_xxyyy_xy, g_x_z_xxyyy_xz, g_x_z_xxyyy_yy, g_x_z_xxyyy_yz, g_x_z_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyyy_xx[k] = -g_x_z_xxyy_xx[k] * ab_y + g_x_z_xxyy_xxy[k];

                g_x_z_xxyyy_xy[k] = -g_x_z_xxyy_xy[k] * ab_y + g_x_z_xxyy_xyy[k];

                g_x_z_xxyyy_xz[k] = -g_x_z_xxyy_xz[k] * ab_y + g_x_z_xxyy_xyz[k];

                g_x_z_xxyyy_yy[k] = -g_x_z_xxyy_yy[k] * ab_y + g_x_z_xxyy_yyy[k];

                g_x_z_xxyyy_yz[k] = -g_x_z_xxyy_yz[k] * ab_y + g_x_z_xxyy_yyz[k];

                g_x_z_xxyyy_zz[k] = -g_x_z_xxyy_zz[k] * ab_y + g_x_z_xxyy_yzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyyz_xx = cbuffer.data(hd_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_z_xxyyz_xy = cbuffer.data(hd_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_z_xxyyz_xz = cbuffer.data(hd_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_z_xxyyz_yy = cbuffer.data(hd_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_z_xxyyz_yz = cbuffer.data(hd_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_z_xxyyz_zz = cbuffer.data(hd_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyyz_xx, g_x_z_xxyyz_xy, g_x_z_xxyyz_xz, g_x_z_xxyyz_yy, g_x_z_xxyyz_yz, g_x_z_xxyyz_zz, g_x_z_xxyz_xx, g_x_z_xxyz_xxy, g_x_z_xxyz_xy, g_x_z_xxyz_xyy, g_x_z_xxyz_xyz, g_x_z_xxyz_xz, g_x_z_xxyz_yy, g_x_z_xxyz_yyy, g_x_z_xxyz_yyz, g_x_z_xxyz_yz, g_x_z_xxyz_yzz, g_x_z_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyyz_xx[k] = -g_x_z_xxyz_xx[k] * ab_y + g_x_z_xxyz_xxy[k];

                g_x_z_xxyyz_xy[k] = -g_x_z_xxyz_xy[k] * ab_y + g_x_z_xxyz_xyy[k];

                g_x_z_xxyyz_xz[k] = -g_x_z_xxyz_xz[k] * ab_y + g_x_z_xxyz_xyz[k];

                g_x_z_xxyyz_yy[k] = -g_x_z_xxyz_yy[k] * ab_y + g_x_z_xxyz_yyy[k];

                g_x_z_xxyyz_yz[k] = -g_x_z_xxyz_yz[k] * ab_y + g_x_z_xxyz_yyz[k];

                g_x_z_xxyyz_zz[k] = -g_x_z_xxyz_zz[k] * ab_y + g_x_z_xxyz_yzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyzz_xx = cbuffer.data(hd_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_xxyzz_xy = cbuffer.data(hd_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_xxyzz_xz = cbuffer.data(hd_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_xxyzz_yy = cbuffer.data(hd_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_xxyzz_yz = cbuffer.data(hd_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_xxyzz_zz = cbuffer.data(hd_geom_11_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyzz_xx, g_x_z_xxyzz_xy, g_x_z_xxyzz_xz, g_x_z_xxyzz_yy, g_x_z_xxyzz_yz, g_x_z_xxyzz_zz, g_x_z_xxzz_xx, g_x_z_xxzz_xxy, g_x_z_xxzz_xy, g_x_z_xxzz_xyy, g_x_z_xxzz_xyz, g_x_z_xxzz_xz, g_x_z_xxzz_yy, g_x_z_xxzz_yyy, g_x_z_xxzz_yyz, g_x_z_xxzz_yz, g_x_z_xxzz_yzz, g_x_z_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyzz_xx[k] = -g_x_z_xxzz_xx[k] * ab_y + g_x_z_xxzz_xxy[k];

                g_x_z_xxyzz_xy[k] = -g_x_z_xxzz_xy[k] * ab_y + g_x_z_xxzz_xyy[k];

                g_x_z_xxyzz_xz[k] = -g_x_z_xxzz_xz[k] * ab_y + g_x_z_xxzz_xyz[k];

                g_x_z_xxyzz_yy[k] = -g_x_z_xxzz_yy[k] * ab_y + g_x_z_xxzz_yyy[k];

                g_x_z_xxyzz_yz[k] = -g_x_z_xxzz_yz[k] * ab_y + g_x_z_xxzz_yyz[k];

                g_x_z_xxyzz_zz[k] = -g_x_z_xxzz_zz[k] * ab_y + g_x_z_xxzz_yzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxzzz_xx = cbuffer.data(hd_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_xxzzz_xy = cbuffer.data(hd_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_xxzzz_xz = cbuffer.data(hd_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_xxzzz_yy = cbuffer.data(hd_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_xxzzz_yz = cbuffer.data(hd_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_xxzzz_zz = cbuffer.data(hd_geom_11_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_xx, g_0_z_xzzz_xy, g_0_z_xzzz_xz, g_0_z_xzzz_yy, g_0_z_xzzz_yz, g_0_z_xzzz_zz, g_x_z_xxzzz_xx, g_x_z_xxzzz_xy, g_x_z_xxzzz_xz, g_x_z_xxzzz_yy, g_x_z_xxzzz_yz, g_x_z_xxzzz_zz, g_x_z_xzzz_xx, g_x_z_xzzz_xxx, g_x_z_xzzz_xxy, g_x_z_xzzz_xxz, g_x_z_xzzz_xy, g_x_z_xzzz_xyy, g_x_z_xzzz_xyz, g_x_z_xzzz_xz, g_x_z_xzzz_xzz, g_x_z_xzzz_yy, g_x_z_xzzz_yz, g_x_z_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxzzz_xx[k] = -g_0_z_xzzz_xx[k] - g_x_z_xzzz_xx[k] * ab_x + g_x_z_xzzz_xxx[k];

                g_x_z_xxzzz_xy[k] = -g_0_z_xzzz_xy[k] - g_x_z_xzzz_xy[k] * ab_x + g_x_z_xzzz_xxy[k];

                g_x_z_xxzzz_xz[k] = -g_0_z_xzzz_xz[k] - g_x_z_xzzz_xz[k] * ab_x + g_x_z_xzzz_xxz[k];

                g_x_z_xxzzz_yy[k] = -g_0_z_xzzz_yy[k] - g_x_z_xzzz_yy[k] * ab_x + g_x_z_xzzz_xyy[k];

                g_x_z_xxzzz_yz[k] = -g_0_z_xzzz_yz[k] - g_x_z_xzzz_yz[k] * ab_x + g_x_z_xzzz_xyz[k];

                g_x_z_xxzzz_zz[k] = -g_0_z_xzzz_zz[k] - g_x_z_xzzz_zz[k] * ab_x + g_x_z_xzzz_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyyy_xx = cbuffer.data(hd_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_xyyyy_xy = cbuffer.data(hd_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_xyyyy_xz = cbuffer.data(hd_geom_11_off + 314 * ccomps * dcomps);

            auto g_x_z_xyyyy_yy = cbuffer.data(hd_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_xyyyy_yz = cbuffer.data(hd_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_xyyyy_zz = cbuffer.data(hd_geom_11_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyy_xx, g_x_z_xyyy_xxy, g_x_z_xyyy_xy, g_x_z_xyyy_xyy, g_x_z_xyyy_xyz, g_x_z_xyyy_xz, g_x_z_xyyy_yy, g_x_z_xyyy_yyy, g_x_z_xyyy_yyz, g_x_z_xyyy_yz, g_x_z_xyyy_yzz, g_x_z_xyyy_zz, g_x_z_xyyyy_xx, g_x_z_xyyyy_xy, g_x_z_xyyyy_xz, g_x_z_xyyyy_yy, g_x_z_xyyyy_yz, g_x_z_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyyy_xx[k] = -g_x_z_xyyy_xx[k] * ab_y + g_x_z_xyyy_xxy[k];

                g_x_z_xyyyy_xy[k] = -g_x_z_xyyy_xy[k] * ab_y + g_x_z_xyyy_xyy[k];

                g_x_z_xyyyy_xz[k] = -g_x_z_xyyy_xz[k] * ab_y + g_x_z_xyyy_xyz[k];

                g_x_z_xyyyy_yy[k] = -g_x_z_xyyy_yy[k] * ab_y + g_x_z_xyyy_yyy[k];

                g_x_z_xyyyy_yz[k] = -g_x_z_xyyy_yz[k] * ab_y + g_x_z_xyyy_yyz[k];

                g_x_z_xyyyy_zz[k] = -g_x_z_xyyy_zz[k] * ab_y + g_x_z_xyyy_yzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyyz_xx = cbuffer.data(hd_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_xyyyz_xy = cbuffer.data(hd_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_xyyyz_xz = cbuffer.data(hd_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_xyyyz_yy = cbuffer.data(hd_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_xyyyz_yz = cbuffer.data(hd_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_xyyyz_zz = cbuffer.data(hd_geom_11_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyyz_xx, g_x_z_xyyyz_xy, g_x_z_xyyyz_xz, g_x_z_xyyyz_yy, g_x_z_xyyyz_yz, g_x_z_xyyyz_zz, g_x_z_xyyz_xx, g_x_z_xyyz_xxy, g_x_z_xyyz_xy, g_x_z_xyyz_xyy, g_x_z_xyyz_xyz, g_x_z_xyyz_xz, g_x_z_xyyz_yy, g_x_z_xyyz_yyy, g_x_z_xyyz_yyz, g_x_z_xyyz_yz, g_x_z_xyyz_yzz, g_x_z_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyyz_xx[k] = -g_x_z_xyyz_xx[k] * ab_y + g_x_z_xyyz_xxy[k];

                g_x_z_xyyyz_xy[k] = -g_x_z_xyyz_xy[k] * ab_y + g_x_z_xyyz_xyy[k];

                g_x_z_xyyyz_xz[k] = -g_x_z_xyyz_xz[k] * ab_y + g_x_z_xyyz_xyz[k];

                g_x_z_xyyyz_yy[k] = -g_x_z_xyyz_yy[k] * ab_y + g_x_z_xyyz_yyy[k];

                g_x_z_xyyyz_yz[k] = -g_x_z_xyyz_yz[k] * ab_y + g_x_z_xyyz_yyz[k];

                g_x_z_xyyyz_zz[k] = -g_x_z_xyyz_zz[k] * ab_y + g_x_z_xyyz_yzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyzz_xx = cbuffer.data(hd_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_z_xyyzz_xy = cbuffer.data(hd_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_z_xyyzz_xz = cbuffer.data(hd_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_z_xyyzz_yy = cbuffer.data(hd_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_z_xyyzz_yz = cbuffer.data(hd_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_z_xyyzz_zz = cbuffer.data(hd_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyzz_xx, g_x_z_xyyzz_xy, g_x_z_xyyzz_xz, g_x_z_xyyzz_yy, g_x_z_xyyzz_yz, g_x_z_xyyzz_zz, g_x_z_xyzz_xx, g_x_z_xyzz_xxy, g_x_z_xyzz_xy, g_x_z_xyzz_xyy, g_x_z_xyzz_xyz, g_x_z_xyzz_xz, g_x_z_xyzz_yy, g_x_z_xyzz_yyy, g_x_z_xyzz_yyz, g_x_z_xyzz_yz, g_x_z_xyzz_yzz, g_x_z_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyzz_xx[k] = -g_x_z_xyzz_xx[k] * ab_y + g_x_z_xyzz_xxy[k];

                g_x_z_xyyzz_xy[k] = -g_x_z_xyzz_xy[k] * ab_y + g_x_z_xyzz_xyy[k];

                g_x_z_xyyzz_xz[k] = -g_x_z_xyzz_xz[k] * ab_y + g_x_z_xyzz_xyz[k];

                g_x_z_xyyzz_yy[k] = -g_x_z_xyzz_yy[k] * ab_y + g_x_z_xyzz_yyy[k];

                g_x_z_xyyzz_yz[k] = -g_x_z_xyzz_yz[k] * ab_y + g_x_z_xyzz_yyz[k];

                g_x_z_xyyzz_zz[k] = -g_x_z_xyzz_zz[k] * ab_y + g_x_z_xyzz_yzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyzzz_xx = cbuffer.data(hd_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_z_xyzzz_xy = cbuffer.data(hd_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_z_xyzzz_xz = cbuffer.data(hd_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_z_xyzzz_yy = cbuffer.data(hd_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_z_xyzzz_yz = cbuffer.data(hd_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_z_xyzzz_zz = cbuffer.data(hd_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyzzz_xx, g_x_z_xyzzz_xy, g_x_z_xyzzz_xz, g_x_z_xyzzz_yy, g_x_z_xyzzz_yz, g_x_z_xyzzz_zz, g_x_z_xzzz_xx, g_x_z_xzzz_xxy, g_x_z_xzzz_xy, g_x_z_xzzz_xyy, g_x_z_xzzz_xyz, g_x_z_xzzz_xz, g_x_z_xzzz_yy, g_x_z_xzzz_yyy, g_x_z_xzzz_yyz, g_x_z_xzzz_yz, g_x_z_xzzz_yzz, g_x_z_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyzzz_xx[k] = -g_x_z_xzzz_xx[k] * ab_y + g_x_z_xzzz_xxy[k];

                g_x_z_xyzzz_xy[k] = -g_x_z_xzzz_xy[k] * ab_y + g_x_z_xzzz_xyy[k];

                g_x_z_xyzzz_xz[k] = -g_x_z_xzzz_xz[k] * ab_y + g_x_z_xzzz_xyz[k];

                g_x_z_xyzzz_yy[k] = -g_x_z_xzzz_yy[k] * ab_y + g_x_z_xzzz_yyy[k];

                g_x_z_xyzzz_yz[k] = -g_x_z_xzzz_yz[k] * ab_y + g_x_z_xzzz_yyz[k];

                g_x_z_xyzzz_zz[k] = -g_x_z_xzzz_zz[k] * ab_y + g_x_z_xzzz_yzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzzzz_xx = cbuffer.data(hd_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_xzzzz_xy = cbuffer.data(hd_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_xzzzz_xz = cbuffer.data(hd_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_xzzzz_yy = cbuffer.data(hd_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_xzzzz_yz = cbuffer.data(hd_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_xzzzz_zz = cbuffer.data(hd_geom_11_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xx, g_0_z_zzzz_xy, g_0_z_zzzz_xz, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_zz, g_x_z_xzzzz_xx, g_x_z_xzzzz_xy, g_x_z_xzzzz_xz, g_x_z_xzzzz_yy, g_x_z_xzzzz_yz, g_x_z_xzzzz_zz, g_x_z_zzzz_xx, g_x_z_zzzz_xxx, g_x_z_zzzz_xxy, g_x_z_zzzz_xxz, g_x_z_zzzz_xy, g_x_z_zzzz_xyy, g_x_z_zzzz_xyz, g_x_z_zzzz_xz, g_x_z_zzzz_xzz, g_x_z_zzzz_yy, g_x_z_zzzz_yz, g_x_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzzzz_xx[k] = -g_0_z_zzzz_xx[k] - g_x_z_zzzz_xx[k] * ab_x + g_x_z_zzzz_xxx[k];

                g_x_z_xzzzz_xy[k] = -g_0_z_zzzz_xy[k] - g_x_z_zzzz_xy[k] * ab_x + g_x_z_zzzz_xxy[k];

                g_x_z_xzzzz_xz[k] = -g_0_z_zzzz_xz[k] - g_x_z_zzzz_xz[k] * ab_x + g_x_z_zzzz_xxz[k];

                g_x_z_xzzzz_yy[k] = -g_0_z_zzzz_yy[k] - g_x_z_zzzz_yy[k] * ab_x + g_x_z_zzzz_xyy[k];

                g_x_z_xzzzz_yz[k] = -g_0_z_zzzz_yz[k] - g_x_z_zzzz_yz[k] * ab_x + g_x_z_zzzz_xyz[k];

                g_x_z_xzzzz_zz[k] = -g_0_z_zzzz_zz[k] - g_x_z_zzzz_zz[k] * ab_x + g_x_z_zzzz_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyyy_xx = cbuffer.data(hd_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_yyyyy_xy = cbuffer.data(hd_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_yyyyy_xz = cbuffer.data(hd_geom_11_off + 344 * ccomps * dcomps);

            auto g_x_z_yyyyy_yy = cbuffer.data(hd_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_yyyyy_yz = cbuffer.data(hd_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_yyyyy_zz = cbuffer.data(hd_geom_11_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyy_xx, g_x_z_yyyy_xxy, g_x_z_yyyy_xy, g_x_z_yyyy_xyy, g_x_z_yyyy_xyz, g_x_z_yyyy_xz, g_x_z_yyyy_yy, g_x_z_yyyy_yyy, g_x_z_yyyy_yyz, g_x_z_yyyy_yz, g_x_z_yyyy_yzz, g_x_z_yyyy_zz, g_x_z_yyyyy_xx, g_x_z_yyyyy_xy, g_x_z_yyyyy_xz, g_x_z_yyyyy_yy, g_x_z_yyyyy_yz, g_x_z_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyyy_xx[k] = -g_x_z_yyyy_xx[k] * ab_y + g_x_z_yyyy_xxy[k];

                g_x_z_yyyyy_xy[k] = -g_x_z_yyyy_xy[k] * ab_y + g_x_z_yyyy_xyy[k];

                g_x_z_yyyyy_xz[k] = -g_x_z_yyyy_xz[k] * ab_y + g_x_z_yyyy_xyz[k];

                g_x_z_yyyyy_yy[k] = -g_x_z_yyyy_yy[k] * ab_y + g_x_z_yyyy_yyy[k];

                g_x_z_yyyyy_yz[k] = -g_x_z_yyyy_yz[k] * ab_y + g_x_z_yyyy_yyz[k];

                g_x_z_yyyyy_zz[k] = -g_x_z_yyyy_zz[k] * ab_y + g_x_z_yyyy_yzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyyz_xx = cbuffer.data(hd_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_yyyyz_xy = cbuffer.data(hd_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_yyyyz_xz = cbuffer.data(hd_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_yyyyz_yy = cbuffer.data(hd_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_yyyyz_yz = cbuffer.data(hd_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_yyyyz_zz = cbuffer.data(hd_geom_11_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyyz_xx, g_x_z_yyyyz_xy, g_x_z_yyyyz_xz, g_x_z_yyyyz_yy, g_x_z_yyyyz_yz, g_x_z_yyyyz_zz, g_x_z_yyyz_xx, g_x_z_yyyz_xxy, g_x_z_yyyz_xy, g_x_z_yyyz_xyy, g_x_z_yyyz_xyz, g_x_z_yyyz_xz, g_x_z_yyyz_yy, g_x_z_yyyz_yyy, g_x_z_yyyz_yyz, g_x_z_yyyz_yz, g_x_z_yyyz_yzz, g_x_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyyz_xx[k] = -g_x_z_yyyz_xx[k] * ab_y + g_x_z_yyyz_xxy[k];

                g_x_z_yyyyz_xy[k] = -g_x_z_yyyz_xy[k] * ab_y + g_x_z_yyyz_xyy[k];

                g_x_z_yyyyz_xz[k] = -g_x_z_yyyz_xz[k] * ab_y + g_x_z_yyyz_xyz[k];

                g_x_z_yyyyz_yy[k] = -g_x_z_yyyz_yy[k] * ab_y + g_x_z_yyyz_yyy[k];

                g_x_z_yyyyz_yz[k] = -g_x_z_yyyz_yz[k] * ab_y + g_x_z_yyyz_yyz[k];

                g_x_z_yyyyz_zz[k] = -g_x_z_yyyz_zz[k] * ab_y + g_x_z_yyyz_yzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyzz_xx = cbuffer.data(hd_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_yyyzz_xy = cbuffer.data(hd_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_yyyzz_xz = cbuffer.data(hd_geom_11_off + 356 * ccomps * dcomps);

            auto g_x_z_yyyzz_yy = cbuffer.data(hd_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_yyyzz_yz = cbuffer.data(hd_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_yyyzz_zz = cbuffer.data(hd_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyzz_xx, g_x_z_yyyzz_xy, g_x_z_yyyzz_xz, g_x_z_yyyzz_yy, g_x_z_yyyzz_yz, g_x_z_yyyzz_zz, g_x_z_yyzz_xx, g_x_z_yyzz_xxy, g_x_z_yyzz_xy, g_x_z_yyzz_xyy, g_x_z_yyzz_xyz, g_x_z_yyzz_xz, g_x_z_yyzz_yy, g_x_z_yyzz_yyy, g_x_z_yyzz_yyz, g_x_z_yyzz_yz, g_x_z_yyzz_yzz, g_x_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyzz_xx[k] = -g_x_z_yyzz_xx[k] * ab_y + g_x_z_yyzz_xxy[k];

                g_x_z_yyyzz_xy[k] = -g_x_z_yyzz_xy[k] * ab_y + g_x_z_yyzz_xyy[k];

                g_x_z_yyyzz_xz[k] = -g_x_z_yyzz_xz[k] * ab_y + g_x_z_yyzz_xyz[k];

                g_x_z_yyyzz_yy[k] = -g_x_z_yyzz_yy[k] * ab_y + g_x_z_yyzz_yyy[k];

                g_x_z_yyyzz_yz[k] = -g_x_z_yyzz_yz[k] * ab_y + g_x_z_yyzz_yyz[k];

                g_x_z_yyyzz_zz[k] = -g_x_z_yyzz_zz[k] * ab_y + g_x_z_yyzz_yzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyzzz_xx = cbuffer.data(hd_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_yyzzz_xy = cbuffer.data(hd_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_yyzzz_xz = cbuffer.data(hd_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_yyzzz_yy = cbuffer.data(hd_geom_11_off + 363 * ccomps * dcomps);

            auto g_x_z_yyzzz_yz = cbuffer.data(hd_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_yyzzz_zz = cbuffer.data(hd_geom_11_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyzzz_xx, g_x_z_yyzzz_xy, g_x_z_yyzzz_xz, g_x_z_yyzzz_yy, g_x_z_yyzzz_yz, g_x_z_yyzzz_zz, g_x_z_yzzz_xx, g_x_z_yzzz_xxy, g_x_z_yzzz_xy, g_x_z_yzzz_xyy, g_x_z_yzzz_xyz, g_x_z_yzzz_xz, g_x_z_yzzz_yy, g_x_z_yzzz_yyy, g_x_z_yzzz_yyz, g_x_z_yzzz_yz, g_x_z_yzzz_yzz, g_x_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyzzz_xx[k] = -g_x_z_yzzz_xx[k] * ab_y + g_x_z_yzzz_xxy[k];

                g_x_z_yyzzz_xy[k] = -g_x_z_yzzz_xy[k] * ab_y + g_x_z_yzzz_xyy[k];

                g_x_z_yyzzz_xz[k] = -g_x_z_yzzz_xz[k] * ab_y + g_x_z_yzzz_xyz[k];

                g_x_z_yyzzz_yy[k] = -g_x_z_yzzz_yy[k] * ab_y + g_x_z_yzzz_yyy[k];

                g_x_z_yyzzz_yz[k] = -g_x_z_yzzz_yz[k] * ab_y + g_x_z_yzzz_yyz[k];

                g_x_z_yyzzz_zz[k] = -g_x_z_yzzz_zz[k] * ab_y + g_x_z_yzzz_yzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzzzz_xx = cbuffer.data(hd_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_yzzzz_xy = cbuffer.data(hd_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_yzzzz_xz = cbuffer.data(hd_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_yzzzz_yy = cbuffer.data(hd_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_yzzzz_yz = cbuffer.data(hd_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_yzzzz_zz = cbuffer.data(hd_geom_11_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzzzz_xx, g_x_z_yzzzz_xy, g_x_z_yzzzz_xz, g_x_z_yzzzz_yy, g_x_z_yzzzz_yz, g_x_z_yzzzz_zz, g_x_z_zzzz_xx, g_x_z_zzzz_xxy, g_x_z_zzzz_xy, g_x_z_zzzz_xyy, g_x_z_zzzz_xyz, g_x_z_zzzz_xz, g_x_z_zzzz_yy, g_x_z_zzzz_yyy, g_x_z_zzzz_yyz, g_x_z_zzzz_yz, g_x_z_zzzz_yzz, g_x_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzzzz_xx[k] = -g_x_z_zzzz_xx[k] * ab_y + g_x_z_zzzz_xxy[k];

                g_x_z_yzzzz_xy[k] = -g_x_z_zzzz_xy[k] * ab_y + g_x_z_zzzz_xyy[k];

                g_x_z_yzzzz_xz[k] = -g_x_z_zzzz_xz[k] * ab_y + g_x_z_zzzz_xyz[k];

                g_x_z_yzzzz_yy[k] = -g_x_z_zzzz_yy[k] * ab_y + g_x_z_zzzz_yyy[k];

                g_x_z_yzzzz_yz[k] = -g_x_z_zzzz_yz[k] * ab_y + g_x_z_zzzz_yyz[k];

                g_x_z_yzzzz_zz[k] = -g_x_z_zzzz_zz[k] * ab_y + g_x_z_zzzz_yzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzzzz_xx = cbuffer.data(hd_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_zzzzz_xy = cbuffer.data(hd_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_zzzzz_xz = cbuffer.data(hd_geom_11_off + 374 * ccomps * dcomps);

            auto g_x_z_zzzzz_yy = cbuffer.data(hd_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_zzzzz_yz = cbuffer.data(hd_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_zzzzz_zz = cbuffer.data(hd_geom_11_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_xx, g_x_0_zzzz_xy, g_x_0_zzzz_xz, g_x_0_zzzz_yy, g_x_0_zzzz_yz, g_x_0_zzzz_zz, g_x_z_zzzz_xx, g_x_z_zzzz_xxz, g_x_z_zzzz_xy, g_x_z_zzzz_xyz, g_x_z_zzzz_xz, g_x_z_zzzz_xzz, g_x_z_zzzz_yy, g_x_z_zzzz_yyz, g_x_z_zzzz_yz, g_x_z_zzzz_yzz, g_x_z_zzzz_zz, g_x_z_zzzz_zzz, g_x_z_zzzzz_xx, g_x_z_zzzzz_xy, g_x_z_zzzzz_xz, g_x_z_zzzzz_yy, g_x_z_zzzzz_yz, g_x_z_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzzzz_xx[k] = g_x_0_zzzz_xx[k] - g_x_z_zzzz_xx[k] * ab_z + g_x_z_zzzz_xxz[k];

                g_x_z_zzzzz_xy[k] = g_x_0_zzzz_xy[k] - g_x_z_zzzz_xy[k] * ab_z + g_x_z_zzzz_xyz[k];

                g_x_z_zzzzz_xz[k] = g_x_0_zzzz_xz[k] - g_x_z_zzzz_xz[k] * ab_z + g_x_z_zzzz_xzz[k];

                g_x_z_zzzzz_yy[k] = g_x_0_zzzz_yy[k] - g_x_z_zzzz_yy[k] * ab_z + g_x_z_zzzz_yyz[k];

                g_x_z_zzzzz_yz[k] = g_x_0_zzzz_yz[k] - g_x_z_zzzz_yz[k] * ab_z + g_x_z_zzzz_yzz[k];

                g_x_z_zzzzz_zz[k] = g_x_0_zzzz_zz[k] - g_x_z_zzzz_zz[k] * ab_z + g_x_z_zzzz_zzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxx_xx = cbuffer.data(hd_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_x_xxxxx_xy = cbuffer.data(hd_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_x_xxxxx_xz = cbuffer.data(hd_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_x_xxxxx_yy = cbuffer.data(hd_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_x_xxxxx_yz = cbuffer.data(hd_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_x_xxxxx_zz = cbuffer.data(hd_geom_11_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxx_xx, g_y_0_xxxx_xy, g_y_0_xxxx_xz, g_y_0_xxxx_yy, g_y_0_xxxx_yz, g_y_0_xxxx_zz, g_y_x_xxxx_xx, g_y_x_xxxx_xxx, g_y_x_xxxx_xxy, g_y_x_xxxx_xxz, g_y_x_xxxx_xy, g_y_x_xxxx_xyy, g_y_x_xxxx_xyz, g_y_x_xxxx_xz, g_y_x_xxxx_xzz, g_y_x_xxxx_yy, g_y_x_xxxx_yz, g_y_x_xxxx_zz, g_y_x_xxxxx_xx, g_y_x_xxxxx_xy, g_y_x_xxxxx_xz, g_y_x_xxxxx_yy, g_y_x_xxxxx_yz, g_y_x_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxx_xx[k] = g_y_0_xxxx_xx[k] - g_y_x_xxxx_xx[k] * ab_x + g_y_x_xxxx_xxx[k];

                g_y_x_xxxxx_xy[k] = g_y_0_xxxx_xy[k] - g_y_x_xxxx_xy[k] * ab_x + g_y_x_xxxx_xxy[k];

                g_y_x_xxxxx_xz[k] = g_y_0_xxxx_xz[k] - g_y_x_xxxx_xz[k] * ab_x + g_y_x_xxxx_xxz[k];

                g_y_x_xxxxx_yy[k] = g_y_0_xxxx_yy[k] - g_y_x_xxxx_yy[k] * ab_x + g_y_x_xxxx_xyy[k];

                g_y_x_xxxxx_yz[k] = g_y_0_xxxx_yz[k] - g_y_x_xxxx_yz[k] * ab_x + g_y_x_xxxx_xyz[k];

                g_y_x_xxxxx_zz[k] = g_y_0_xxxx_zz[k] - g_y_x_xxxx_zz[k] * ab_x + g_y_x_xxxx_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxy_xx = cbuffer.data(hd_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_x_xxxxy_xy = cbuffer.data(hd_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_x_xxxxy_xz = cbuffer.data(hd_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_x_xxxxy_yy = cbuffer.data(hd_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_x_xxxxy_yz = cbuffer.data(hd_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_x_xxxxy_zz = cbuffer.data(hd_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_xx, g_y_0_xxxy_xy, g_y_0_xxxy_xz, g_y_0_xxxy_yy, g_y_0_xxxy_yz, g_y_0_xxxy_zz, g_y_x_xxxxy_xx, g_y_x_xxxxy_xy, g_y_x_xxxxy_xz, g_y_x_xxxxy_yy, g_y_x_xxxxy_yz, g_y_x_xxxxy_zz, g_y_x_xxxy_xx, g_y_x_xxxy_xxx, g_y_x_xxxy_xxy, g_y_x_xxxy_xxz, g_y_x_xxxy_xy, g_y_x_xxxy_xyy, g_y_x_xxxy_xyz, g_y_x_xxxy_xz, g_y_x_xxxy_xzz, g_y_x_xxxy_yy, g_y_x_xxxy_yz, g_y_x_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxy_xx[k] = g_y_0_xxxy_xx[k] - g_y_x_xxxy_xx[k] * ab_x + g_y_x_xxxy_xxx[k];

                g_y_x_xxxxy_xy[k] = g_y_0_xxxy_xy[k] - g_y_x_xxxy_xy[k] * ab_x + g_y_x_xxxy_xxy[k];

                g_y_x_xxxxy_xz[k] = g_y_0_xxxy_xz[k] - g_y_x_xxxy_xz[k] * ab_x + g_y_x_xxxy_xxz[k];

                g_y_x_xxxxy_yy[k] = g_y_0_xxxy_yy[k] - g_y_x_xxxy_yy[k] * ab_x + g_y_x_xxxy_xyy[k];

                g_y_x_xxxxy_yz[k] = g_y_0_xxxy_yz[k] - g_y_x_xxxy_yz[k] * ab_x + g_y_x_xxxy_xyz[k];

                g_y_x_xxxxy_zz[k] = g_y_0_xxxy_zz[k] - g_y_x_xxxy_zz[k] * ab_x + g_y_x_xxxy_xzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxxz_xx = cbuffer.data(hd_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_x_xxxxz_xy = cbuffer.data(hd_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_x_xxxxz_xz = cbuffer.data(hd_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_x_xxxxz_yy = cbuffer.data(hd_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_x_xxxxz_yz = cbuffer.data(hd_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_x_xxxxz_zz = cbuffer.data(hd_geom_11_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxx_xx, g_y_x_xxxx_xxz, g_y_x_xxxx_xy, g_y_x_xxxx_xyz, g_y_x_xxxx_xz, g_y_x_xxxx_xzz, g_y_x_xxxx_yy, g_y_x_xxxx_yyz, g_y_x_xxxx_yz, g_y_x_xxxx_yzz, g_y_x_xxxx_zz, g_y_x_xxxx_zzz, g_y_x_xxxxz_xx, g_y_x_xxxxz_xy, g_y_x_xxxxz_xz, g_y_x_xxxxz_yy, g_y_x_xxxxz_yz, g_y_x_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxxz_xx[k] = -g_y_x_xxxx_xx[k] * ab_z + g_y_x_xxxx_xxz[k];

                g_y_x_xxxxz_xy[k] = -g_y_x_xxxx_xy[k] * ab_z + g_y_x_xxxx_xyz[k];

                g_y_x_xxxxz_xz[k] = -g_y_x_xxxx_xz[k] * ab_z + g_y_x_xxxx_xzz[k];

                g_y_x_xxxxz_yy[k] = -g_y_x_xxxx_yy[k] * ab_z + g_y_x_xxxx_yyz[k];

                g_y_x_xxxxz_yz[k] = -g_y_x_xxxx_yz[k] * ab_z + g_y_x_xxxx_yzz[k];

                g_y_x_xxxxz_zz[k] = -g_y_x_xxxx_zz[k] * ab_z + g_y_x_xxxx_zzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxyy_xx = cbuffer.data(hd_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_x_xxxyy_xy = cbuffer.data(hd_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_x_xxxyy_xz = cbuffer.data(hd_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_x_xxxyy_yy = cbuffer.data(hd_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_x_xxxyy_yz = cbuffer.data(hd_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_x_xxxyy_zz = cbuffer.data(hd_geom_11_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_xx, g_y_0_xxyy_xy, g_y_0_xxyy_xz, g_y_0_xxyy_yy, g_y_0_xxyy_yz, g_y_0_xxyy_zz, g_y_x_xxxyy_xx, g_y_x_xxxyy_xy, g_y_x_xxxyy_xz, g_y_x_xxxyy_yy, g_y_x_xxxyy_yz, g_y_x_xxxyy_zz, g_y_x_xxyy_xx, g_y_x_xxyy_xxx, g_y_x_xxyy_xxy, g_y_x_xxyy_xxz, g_y_x_xxyy_xy, g_y_x_xxyy_xyy, g_y_x_xxyy_xyz, g_y_x_xxyy_xz, g_y_x_xxyy_xzz, g_y_x_xxyy_yy, g_y_x_xxyy_yz, g_y_x_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxyy_xx[k] = g_y_0_xxyy_xx[k] - g_y_x_xxyy_xx[k] * ab_x + g_y_x_xxyy_xxx[k];

                g_y_x_xxxyy_xy[k] = g_y_0_xxyy_xy[k] - g_y_x_xxyy_xy[k] * ab_x + g_y_x_xxyy_xxy[k];

                g_y_x_xxxyy_xz[k] = g_y_0_xxyy_xz[k] - g_y_x_xxyy_xz[k] * ab_x + g_y_x_xxyy_xxz[k];

                g_y_x_xxxyy_yy[k] = g_y_0_xxyy_yy[k] - g_y_x_xxyy_yy[k] * ab_x + g_y_x_xxyy_xyy[k];

                g_y_x_xxxyy_yz[k] = g_y_0_xxyy_yz[k] - g_y_x_xxyy_yz[k] * ab_x + g_y_x_xxyy_xyz[k];

                g_y_x_xxxyy_zz[k] = g_y_0_xxyy_zz[k] - g_y_x_xxyy_zz[k] * ab_x + g_y_x_xxyy_xzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxyz_xx = cbuffer.data(hd_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_x_xxxyz_xy = cbuffer.data(hd_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_x_xxxyz_xz = cbuffer.data(hd_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_x_xxxyz_yy = cbuffer.data(hd_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_x_xxxyz_yz = cbuffer.data(hd_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_x_xxxyz_zz = cbuffer.data(hd_geom_11_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxy_xx, g_y_x_xxxy_xxz, g_y_x_xxxy_xy, g_y_x_xxxy_xyz, g_y_x_xxxy_xz, g_y_x_xxxy_xzz, g_y_x_xxxy_yy, g_y_x_xxxy_yyz, g_y_x_xxxy_yz, g_y_x_xxxy_yzz, g_y_x_xxxy_zz, g_y_x_xxxy_zzz, g_y_x_xxxyz_xx, g_y_x_xxxyz_xy, g_y_x_xxxyz_xz, g_y_x_xxxyz_yy, g_y_x_xxxyz_yz, g_y_x_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxyz_xx[k] = -g_y_x_xxxy_xx[k] * ab_z + g_y_x_xxxy_xxz[k];

                g_y_x_xxxyz_xy[k] = -g_y_x_xxxy_xy[k] * ab_z + g_y_x_xxxy_xyz[k];

                g_y_x_xxxyz_xz[k] = -g_y_x_xxxy_xz[k] * ab_z + g_y_x_xxxy_xzz[k];

                g_y_x_xxxyz_yy[k] = -g_y_x_xxxy_yy[k] * ab_z + g_y_x_xxxy_yyz[k];

                g_y_x_xxxyz_yz[k] = -g_y_x_xxxy_yz[k] * ab_z + g_y_x_xxxy_yzz[k];

                g_y_x_xxxyz_zz[k] = -g_y_x_xxxy_zz[k] * ab_z + g_y_x_xxxy_zzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxzz_xx = cbuffer.data(hd_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_x_xxxzz_xy = cbuffer.data(hd_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_x_xxxzz_xz = cbuffer.data(hd_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_x_xxxzz_yy = cbuffer.data(hd_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_x_xxxzz_yz = cbuffer.data(hd_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_x_xxxzz_zz = cbuffer.data(hd_geom_11_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxxz_xx, g_y_x_xxxz_xxz, g_y_x_xxxz_xy, g_y_x_xxxz_xyz, g_y_x_xxxz_xz, g_y_x_xxxz_xzz, g_y_x_xxxz_yy, g_y_x_xxxz_yyz, g_y_x_xxxz_yz, g_y_x_xxxz_yzz, g_y_x_xxxz_zz, g_y_x_xxxz_zzz, g_y_x_xxxzz_xx, g_y_x_xxxzz_xy, g_y_x_xxxzz_xz, g_y_x_xxxzz_yy, g_y_x_xxxzz_yz, g_y_x_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxzz_xx[k] = -g_y_x_xxxz_xx[k] * ab_z + g_y_x_xxxz_xxz[k];

                g_y_x_xxxzz_xy[k] = -g_y_x_xxxz_xy[k] * ab_z + g_y_x_xxxz_xyz[k];

                g_y_x_xxxzz_xz[k] = -g_y_x_xxxz_xz[k] * ab_z + g_y_x_xxxz_xzz[k];

                g_y_x_xxxzz_yy[k] = -g_y_x_xxxz_yy[k] * ab_z + g_y_x_xxxz_yyz[k];

                g_y_x_xxxzz_yz[k] = -g_y_x_xxxz_yz[k] * ab_z + g_y_x_xxxz_yzz[k];

                g_y_x_xxxzz_zz[k] = -g_y_x_xxxz_zz[k] * ab_z + g_y_x_xxxz_zzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyyy_xx = cbuffer.data(hd_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_x_xxyyy_xy = cbuffer.data(hd_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_x_xxyyy_xz = cbuffer.data(hd_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_x_xxyyy_yy = cbuffer.data(hd_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_x_xxyyy_yz = cbuffer.data(hd_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_x_xxyyy_zz = cbuffer.data(hd_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_xx, g_y_0_xyyy_xy, g_y_0_xyyy_xz, g_y_0_xyyy_yy, g_y_0_xyyy_yz, g_y_0_xyyy_zz, g_y_x_xxyyy_xx, g_y_x_xxyyy_xy, g_y_x_xxyyy_xz, g_y_x_xxyyy_yy, g_y_x_xxyyy_yz, g_y_x_xxyyy_zz, g_y_x_xyyy_xx, g_y_x_xyyy_xxx, g_y_x_xyyy_xxy, g_y_x_xyyy_xxz, g_y_x_xyyy_xy, g_y_x_xyyy_xyy, g_y_x_xyyy_xyz, g_y_x_xyyy_xz, g_y_x_xyyy_xzz, g_y_x_xyyy_yy, g_y_x_xyyy_yz, g_y_x_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyyy_xx[k] = g_y_0_xyyy_xx[k] - g_y_x_xyyy_xx[k] * ab_x + g_y_x_xyyy_xxx[k];

                g_y_x_xxyyy_xy[k] = g_y_0_xyyy_xy[k] - g_y_x_xyyy_xy[k] * ab_x + g_y_x_xyyy_xxy[k];

                g_y_x_xxyyy_xz[k] = g_y_0_xyyy_xz[k] - g_y_x_xyyy_xz[k] * ab_x + g_y_x_xyyy_xxz[k];

                g_y_x_xxyyy_yy[k] = g_y_0_xyyy_yy[k] - g_y_x_xyyy_yy[k] * ab_x + g_y_x_xyyy_xyy[k];

                g_y_x_xxyyy_yz[k] = g_y_0_xyyy_yz[k] - g_y_x_xyyy_yz[k] * ab_x + g_y_x_xyyy_xyz[k];

                g_y_x_xxyyy_zz[k] = g_y_0_xyyy_zz[k] - g_y_x_xyyy_zz[k] * ab_x + g_y_x_xyyy_xzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyyz_xx = cbuffer.data(hd_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_x_xxyyz_xy = cbuffer.data(hd_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_x_xxyyz_xz = cbuffer.data(hd_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_x_xxyyz_yy = cbuffer.data(hd_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_x_xxyyz_yz = cbuffer.data(hd_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_x_xxyyz_zz = cbuffer.data(hd_geom_11_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxyy_xx, g_y_x_xxyy_xxz, g_y_x_xxyy_xy, g_y_x_xxyy_xyz, g_y_x_xxyy_xz, g_y_x_xxyy_xzz, g_y_x_xxyy_yy, g_y_x_xxyy_yyz, g_y_x_xxyy_yz, g_y_x_xxyy_yzz, g_y_x_xxyy_zz, g_y_x_xxyy_zzz, g_y_x_xxyyz_xx, g_y_x_xxyyz_xy, g_y_x_xxyyz_xz, g_y_x_xxyyz_yy, g_y_x_xxyyz_yz, g_y_x_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyyz_xx[k] = -g_y_x_xxyy_xx[k] * ab_z + g_y_x_xxyy_xxz[k];

                g_y_x_xxyyz_xy[k] = -g_y_x_xxyy_xy[k] * ab_z + g_y_x_xxyy_xyz[k];

                g_y_x_xxyyz_xz[k] = -g_y_x_xxyy_xz[k] * ab_z + g_y_x_xxyy_xzz[k];

                g_y_x_xxyyz_yy[k] = -g_y_x_xxyy_yy[k] * ab_z + g_y_x_xxyy_yyz[k];

                g_y_x_xxyyz_yz[k] = -g_y_x_xxyy_yz[k] * ab_z + g_y_x_xxyy_yzz[k];

                g_y_x_xxyyz_zz[k] = -g_y_x_xxyy_zz[k] * ab_z + g_y_x_xxyy_zzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyzz_xx = cbuffer.data(hd_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_x_xxyzz_xy = cbuffer.data(hd_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_x_xxyzz_xz = cbuffer.data(hd_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_x_xxyzz_yy = cbuffer.data(hd_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_x_xxyzz_yz = cbuffer.data(hd_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_x_xxyzz_zz = cbuffer.data(hd_geom_11_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxyz_xx, g_y_x_xxyz_xxz, g_y_x_xxyz_xy, g_y_x_xxyz_xyz, g_y_x_xxyz_xz, g_y_x_xxyz_xzz, g_y_x_xxyz_yy, g_y_x_xxyz_yyz, g_y_x_xxyz_yz, g_y_x_xxyz_yzz, g_y_x_xxyz_zz, g_y_x_xxyz_zzz, g_y_x_xxyzz_xx, g_y_x_xxyzz_xy, g_y_x_xxyzz_xz, g_y_x_xxyzz_yy, g_y_x_xxyzz_yz, g_y_x_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyzz_xx[k] = -g_y_x_xxyz_xx[k] * ab_z + g_y_x_xxyz_xxz[k];

                g_y_x_xxyzz_xy[k] = -g_y_x_xxyz_xy[k] * ab_z + g_y_x_xxyz_xyz[k];

                g_y_x_xxyzz_xz[k] = -g_y_x_xxyz_xz[k] * ab_z + g_y_x_xxyz_xzz[k];

                g_y_x_xxyzz_yy[k] = -g_y_x_xxyz_yy[k] * ab_z + g_y_x_xxyz_yyz[k];

                g_y_x_xxyzz_yz[k] = -g_y_x_xxyz_yz[k] * ab_z + g_y_x_xxyz_yzz[k];

                g_y_x_xxyzz_zz[k] = -g_y_x_xxyz_zz[k] * ab_z + g_y_x_xxyz_zzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxzzz_xx = cbuffer.data(hd_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_x_xxzzz_xy = cbuffer.data(hd_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_x_xxzzz_xz = cbuffer.data(hd_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_x_xxzzz_yy = cbuffer.data(hd_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_x_xxzzz_yz = cbuffer.data(hd_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_x_xxzzz_zz = cbuffer.data(hd_geom_11_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxzz_xx, g_y_x_xxzz_xxz, g_y_x_xxzz_xy, g_y_x_xxzz_xyz, g_y_x_xxzz_xz, g_y_x_xxzz_xzz, g_y_x_xxzz_yy, g_y_x_xxzz_yyz, g_y_x_xxzz_yz, g_y_x_xxzz_yzz, g_y_x_xxzz_zz, g_y_x_xxzz_zzz, g_y_x_xxzzz_xx, g_y_x_xxzzz_xy, g_y_x_xxzzz_xz, g_y_x_xxzzz_yy, g_y_x_xxzzz_yz, g_y_x_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxzzz_xx[k] = -g_y_x_xxzz_xx[k] * ab_z + g_y_x_xxzz_xxz[k];

                g_y_x_xxzzz_xy[k] = -g_y_x_xxzz_xy[k] * ab_z + g_y_x_xxzz_xyz[k];

                g_y_x_xxzzz_xz[k] = -g_y_x_xxzz_xz[k] * ab_z + g_y_x_xxzz_xzz[k];

                g_y_x_xxzzz_yy[k] = -g_y_x_xxzz_yy[k] * ab_z + g_y_x_xxzz_yyz[k];

                g_y_x_xxzzz_yz[k] = -g_y_x_xxzz_yz[k] * ab_z + g_y_x_xxzz_yzz[k];

                g_y_x_xxzzz_zz[k] = -g_y_x_xxzz_zz[k] * ab_z + g_y_x_xxzz_zzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyyy_xx = cbuffer.data(hd_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_x_xyyyy_xy = cbuffer.data(hd_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_x_xyyyy_xz = cbuffer.data(hd_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_x_xyyyy_yy = cbuffer.data(hd_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_x_xyyyy_yz = cbuffer.data(hd_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_x_xyyyy_zz = cbuffer.data(hd_geom_11_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz, g_y_x_xyyyy_xx, g_y_x_xyyyy_xy, g_y_x_xyyyy_xz, g_y_x_xyyyy_yy, g_y_x_xyyyy_yz, g_y_x_xyyyy_zz, g_y_x_yyyy_xx, g_y_x_yyyy_xxx, g_y_x_yyyy_xxy, g_y_x_yyyy_xxz, g_y_x_yyyy_xy, g_y_x_yyyy_xyy, g_y_x_yyyy_xyz, g_y_x_yyyy_xz, g_y_x_yyyy_xzz, g_y_x_yyyy_yy, g_y_x_yyyy_yz, g_y_x_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyyy_xx[k] = g_y_0_yyyy_xx[k] - g_y_x_yyyy_xx[k] * ab_x + g_y_x_yyyy_xxx[k];

                g_y_x_xyyyy_xy[k] = g_y_0_yyyy_xy[k] - g_y_x_yyyy_xy[k] * ab_x + g_y_x_yyyy_xxy[k];

                g_y_x_xyyyy_xz[k] = g_y_0_yyyy_xz[k] - g_y_x_yyyy_xz[k] * ab_x + g_y_x_yyyy_xxz[k];

                g_y_x_xyyyy_yy[k] = g_y_0_yyyy_yy[k] - g_y_x_yyyy_yy[k] * ab_x + g_y_x_yyyy_xyy[k];

                g_y_x_xyyyy_yz[k] = g_y_0_yyyy_yz[k] - g_y_x_yyyy_yz[k] * ab_x + g_y_x_yyyy_xyz[k];

                g_y_x_xyyyy_zz[k] = g_y_0_yyyy_zz[k] - g_y_x_yyyy_zz[k] * ab_x + g_y_x_yyyy_xzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyyz_xx = cbuffer.data(hd_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_x_xyyyz_xy = cbuffer.data(hd_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_x_xyyyz_xz = cbuffer.data(hd_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_x_xyyyz_yy = cbuffer.data(hd_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_x_xyyyz_yz = cbuffer.data(hd_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_x_xyyyz_zz = cbuffer.data(hd_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyyy_xx, g_y_x_xyyy_xxz, g_y_x_xyyy_xy, g_y_x_xyyy_xyz, g_y_x_xyyy_xz, g_y_x_xyyy_xzz, g_y_x_xyyy_yy, g_y_x_xyyy_yyz, g_y_x_xyyy_yz, g_y_x_xyyy_yzz, g_y_x_xyyy_zz, g_y_x_xyyy_zzz, g_y_x_xyyyz_xx, g_y_x_xyyyz_xy, g_y_x_xyyyz_xz, g_y_x_xyyyz_yy, g_y_x_xyyyz_yz, g_y_x_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyyz_xx[k] = -g_y_x_xyyy_xx[k] * ab_z + g_y_x_xyyy_xxz[k];

                g_y_x_xyyyz_xy[k] = -g_y_x_xyyy_xy[k] * ab_z + g_y_x_xyyy_xyz[k];

                g_y_x_xyyyz_xz[k] = -g_y_x_xyyy_xz[k] * ab_z + g_y_x_xyyy_xzz[k];

                g_y_x_xyyyz_yy[k] = -g_y_x_xyyy_yy[k] * ab_z + g_y_x_xyyy_yyz[k];

                g_y_x_xyyyz_yz[k] = -g_y_x_xyyy_yz[k] * ab_z + g_y_x_xyyy_yzz[k];

                g_y_x_xyyyz_zz[k] = -g_y_x_xyyy_zz[k] * ab_z + g_y_x_xyyy_zzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyzz_xx = cbuffer.data(hd_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_x_xyyzz_xy = cbuffer.data(hd_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_x_xyyzz_xz = cbuffer.data(hd_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_x_xyyzz_yy = cbuffer.data(hd_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_x_xyyzz_yz = cbuffer.data(hd_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_x_xyyzz_zz = cbuffer.data(hd_geom_11_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyyz_xx, g_y_x_xyyz_xxz, g_y_x_xyyz_xy, g_y_x_xyyz_xyz, g_y_x_xyyz_xz, g_y_x_xyyz_xzz, g_y_x_xyyz_yy, g_y_x_xyyz_yyz, g_y_x_xyyz_yz, g_y_x_xyyz_yzz, g_y_x_xyyz_zz, g_y_x_xyyz_zzz, g_y_x_xyyzz_xx, g_y_x_xyyzz_xy, g_y_x_xyyzz_xz, g_y_x_xyyzz_yy, g_y_x_xyyzz_yz, g_y_x_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyzz_xx[k] = -g_y_x_xyyz_xx[k] * ab_z + g_y_x_xyyz_xxz[k];

                g_y_x_xyyzz_xy[k] = -g_y_x_xyyz_xy[k] * ab_z + g_y_x_xyyz_xyz[k];

                g_y_x_xyyzz_xz[k] = -g_y_x_xyyz_xz[k] * ab_z + g_y_x_xyyz_xzz[k];

                g_y_x_xyyzz_yy[k] = -g_y_x_xyyz_yy[k] * ab_z + g_y_x_xyyz_yyz[k];

                g_y_x_xyyzz_yz[k] = -g_y_x_xyyz_yz[k] * ab_z + g_y_x_xyyz_yzz[k];

                g_y_x_xyyzz_zz[k] = -g_y_x_xyyz_zz[k] * ab_z + g_y_x_xyyz_zzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyzzz_xx = cbuffer.data(hd_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_x_xyzzz_xy = cbuffer.data(hd_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_x_xyzzz_xz = cbuffer.data(hd_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_x_xyzzz_yy = cbuffer.data(hd_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_x_xyzzz_yz = cbuffer.data(hd_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_x_xyzzz_zz = cbuffer.data(hd_geom_11_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyzz_xx, g_y_x_xyzz_xxz, g_y_x_xyzz_xy, g_y_x_xyzz_xyz, g_y_x_xyzz_xz, g_y_x_xyzz_xzz, g_y_x_xyzz_yy, g_y_x_xyzz_yyz, g_y_x_xyzz_yz, g_y_x_xyzz_yzz, g_y_x_xyzz_zz, g_y_x_xyzz_zzz, g_y_x_xyzzz_xx, g_y_x_xyzzz_xy, g_y_x_xyzzz_xz, g_y_x_xyzzz_yy, g_y_x_xyzzz_yz, g_y_x_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyzzz_xx[k] = -g_y_x_xyzz_xx[k] * ab_z + g_y_x_xyzz_xxz[k];

                g_y_x_xyzzz_xy[k] = -g_y_x_xyzz_xy[k] * ab_z + g_y_x_xyzz_xyz[k];

                g_y_x_xyzzz_xz[k] = -g_y_x_xyzz_xz[k] * ab_z + g_y_x_xyzz_xzz[k];

                g_y_x_xyzzz_yy[k] = -g_y_x_xyzz_yy[k] * ab_z + g_y_x_xyzz_yyz[k];

                g_y_x_xyzzz_yz[k] = -g_y_x_xyzz_yz[k] * ab_z + g_y_x_xyzz_yzz[k];

                g_y_x_xyzzz_zz[k] = -g_y_x_xyzz_zz[k] * ab_z + g_y_x_xyzz_zzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzzzz_xx = cbuffer.data(hd_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_x_xzzzz_xy = cbuffer.data(hd_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_x_xzzzz_xz = cbuffer.data(hd_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_x_xzzzz_yy = cbuffer.data(hd_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_x_xzzzz_yz = cbuffer.data(hd_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_x_xzzzz_zz = cbuffer.data(hd_geom_11_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xzzz_xx, g_y_x_xzzz_xxz, g_y_x_xzzz_xy, g_y_x_xzzz_xyz, g_y_x_xzzz_xz, g_y_x_xzzz_xzz, g_y_x_xzzz_yy, g_y_x_xzzz_yyz, g_y_x_xzzz_yz, g_y_x_xzzz_yzz, g_y_x_xzzz_zz, g_y_x_xzzz_zzz, g_y_x_xzzzz_xx, g_y_x_xzzzz_xy, g_y_x_xzzzz_xz, g_y_x_xzzzz_yy, g_y_x_xzzzz_yz, g_y_x_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzzzz_xx[k] = -g_y_x_xzzz_xx[k] * ab_z + g_y_x_xzzz_xxz[k];

                g_y_x_xzzzz_xy[k] = -g_y_x_xzzz_xy[k] * ab_z + g_y_x_xzzz_xyz[k];

                g_y_x_xzzzz_xz[k] = -g_y_x_xzzz_xz[k] * ab_z + g_y_x_xzzz_xzz[k];

                g_y_x_xzzzz_yy[k] = -g_y_x_xzzz_yy[k] * ab_z + g_y_x_xzzz_yyz[k];

                g_y_x_xzzzz_yz[k] = -g_y_x_xzzz_yz[k] * ab_z + g_y_x_xzzz_yzz[k];

                g_y_x_xzzzz_zz[k] = -g_y_x_xzzz_zz[k] * ab_z + g_y_x_xzzz_zzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyyy_xx = cbuffer.data(hd_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_x_yyyyy_xy = cbuffer.data(hd_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_x_yyyyy_xz = cbuffer.data(hd_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_x_yyyyy_yy = cbuffer.data(hd_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_x_yyyyy_yz = cbuffer.data(hd_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_x_yyyyy_zz = cbuffer.data(hd_geom_11_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_xx, g_0_x_yyyy_xy, g_0_x_yyyy_xz, g_0_x_yyyy_yy, g_0_x_yyyy_yz, g_0_x_yyyy_zz, g_y_x_yyyy_xx, g_y_x_yyyy_xxy, g_y_x_yyyy_xy, g_y_x_yyyy_xyy, g_y_x_yyyy_xyz, g_y_x_yyyy_xz, g_y_x_yyyy_yy, g_y_x_yyyy_yyy, g_y_x_yyyy_yyz, g_y_x_yyyy_yz, g_y_x_yyyy_yzz, g_y_x_yyyy_zz, g_y_x_yyyyy_xx, g_y_x_yyyyy_xy, g_y_x_yyyyy_xz, g_y_x_yyyyy_yy, g_y_x_yyyyy_yz, g_y_x_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyyy_xx[k] = -g_0_x_yyyy_xx[k] - g_y_x_yyyy_xx[k] * ab_y + g_y_x_yyyy_xxy[k];

                g_y_x_yyyyy_xy[k] = -g_0_x_yyyy_xy[k] - g_y_x_yyyy_xy[k] * ab_y + g_y_x_yyyy_xyy[k];

                g_y_x_yyyyy_xz[k] = -g_0_x_yyyy_xz[k] - g_y_x_yyyy_xz[k] * ab_y + g_y_x_yyyy_xyz[k];

                g_y_x_yyyyy_yy[k] = -g_0_x_yyyy_yy[k] - g_y_x_yyyy_yy[k] * ab_y + g_y_x_yyyy_yyy[k];

                g_y_x_yyyyy_yz[k] = -g_0_x_yyyy_yz[k] - g_y_x_yyyy_yz[k] * ab_y + g_y_x_yyyy_yyz[k];

                g_y_x_yyyyy_zz[k] = -g_0_x_yyyy_zz[k] - g_y_x_yyyy_zz[k] * ab_y + g_y_x_yyyy_yzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyyz_xx = cbuffer.data(hd_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_x_yyyyz_xy = cbuffer.data(hd_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_x_yyyyz_xz = cbuffer.data(hd_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_x_yyyyz_yy = cbuffer.data(hd_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_x_yyyyz_yz = cbuffer.data(hd_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_x_yyyyz_zz = cbuffer.data(hd_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyyy_xx, g_y_x_yyyy_xxz, g_y_x_yyyy_xy, g_y_x_yyyy_xyz, g_y_x_yyyy_xz, g_y_x_yyyy_xzz, g_y_x_yyyy_yy, g_y_x_yyyy_yyz, g_y_x_yyyy_yz, g_y_x_yyyy_yzz, g_y_x_yyyy_zz, g_y_x_yyyy_zzz, g_y_x_yyyyz_xx, g_y_x_yyyyz_xy, g_y_x_yyyyz_xz, g_y_x_yyyyz_yy, g_y_x_yyyyz_yz, g_y_x_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyyz_xx[k] = -g_y_x_yyyy_xx[k] * ab_z + g_y_x_yyyy_xxz[k];

                g_y_x_yyyyz_xy[k] = -g_y_x_yyyy_xy[k] * ab_z + g_y_x_yyyy_xyz[k];

                g_y_x_yyyyz_xz[k] = -g_y_x_yyyy_xz[k] * ab_z + g_y_x_yyyy_xzz[k];

                g_y_x_yyyyz_yy[k] = -g_y_x_yyyy_yy[k] * ab_z + g_y_x_yyyy_yyz[k];

                g_y_x_yyyyz_yz[k] = -g_y_x_yyyy_yz[k] * ab_z + g_y_x_yyyy_yzz[k];

                g_y_x_yyyyz_zz[k] = -g_y_x_yyyy_zz[k] * ab_z + g_y_x_yyyy_zzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyzz_xx = cbuffer.data(hd_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_x_yyyzz_xy = cbuffer.data(hd_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_x_yyyzz_xz = cbuffer.data(hd_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_x_yyyzz_yy = cbuffer.data(hd_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_x_yyyzz_yz = cbuffer.data(hd_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_x_yyyzz_zz = cbuffer.data(hd_geom_11_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyyz_xx, g_y_x_yyyz_xxz, g_y_x_yyyz_xy, g_y_x_yyyz_xyz, g_y_x_yyyz_xz, g_y_x_yyyz_xzz, g_y_x_yyyz_yy, g_y_x_yyyz_yyz, g_y_x_yyyz_yz, g_y_x_yyyz_yzz, g_y_x_yyyz_zz, g_y_x_yyyz_zzz, g_y_x_yyyzz_xx, g_y_x_yyyzz_xy, g_y_x_yyyzz_xz, g_y_x_yyyzz_yy, g_y_x_yyyzz_yz, g_y_x_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyzz_xx[k] = -g_y_x_yyyz_xx[k] * ab_z + g_y_x_yyyz_xxz[k];

                g_y_x_yyyzz_xy[k] = -g_y_x_yyyz_xy[k] * ab_z + g_y_x_yyyz_xyz[k];

                g_y_x_yyyzz_xz[k] = -g_y_x_yyyz_xz[k] * ab_z + g_y_x_yyyz_xzz[k];

                g_y_x_yyyzz_yy[k] = -g_y_x_yyyz_yy[k] * ab_z + g_y_x_yyyz_yyz[k];

                g_y_x_yyyzz_yz[k] = -g_y_x_yyyz_yz[k] * ab_z + g_y_x_yyyz_yzz[k];

                g_y_x_yyyzz_zz[k] = -g_y_x_yyyz_zz[k] * ab_z + g_y_x_yyyz_zzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyzzz_xx = cbuffer.data(hd_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_x_yyzzz_xy = cbuffer.data(hd_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_x_yyzzz_xz = cbuffer.data(hd_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_x_yyzzz_yy = cbuffer.data(hd_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_x_yyzzz_yz = cbuffer.data(hd_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_x_yyzzz_zz = cbuffer.data(hd_geom_11_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyzz_xx, g_y_x_yyzz_xxz, g_y_x_yyzz_xy, g_y_x_yyzz_xyz, g_y_x_yyzz_xz, g_y_x_yyzz_xzz, g_y_x_yyzz_yy, g_y_x_yyzz_yyz, g_y_x_yyzz_yz, g_y_x_yyzz_yzz, g_y_x_yyzz_zz, g_y_x_yyzz_zzz, g_y_x_yyzzz_xx, g_y_x_yyzzz_xy, g_y_x_yyzzz_xz, g_y_x_yyzzz_yy, g_y_x_yyzzz_yz, g_y_x_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyzzz_xx[k] = -g_y_x_yyzz_xx[k] * ab_z + g_y_x_yyzz_xxz[k];

                g_y_x_yyzzz_xy[k] = -g_y_x_yyzz_xy[k] * ab_z + g_y_x_yyzz_xyz[k];

                g_y_x_yyzzz_xz[k] = -g_y_x_yyzz_xz[k] * ab_z + g_y_x_yyzz_xzz[k];

                g_y_x_yyzzz_yy[k] = -g_y_x_yyzz_yy[k] * ab_z + g_y_x_yyzz_yyz[k];

                g_y_x_yyzzz_yz[k] = -g_y_x_yyzz_yz[k] * ab_z + g_y_x_yyzz_yzz[k];

                g_y_x_yyzzz_zz[k] = -g_y_x_yyzz_zz[k] * ab_z + g_y_x_yyzz_zzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzzzz_xx = cbuffer.data(hd_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_x_yzzzz_xy = cbuffer.data(hd_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_x_yzzzz_xz = cbuffer.data(hd_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_x_yzzzz_yy = cbuffer.data(hd_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_x_yzzzz_yz = cbuffer.data(hd_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_x_yzzzz_zz = cbuffer.data(hd_geom_11_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yzzz_xx, g_y_x_yzzz_xxz, g_y_x_yzzz_xy, g_y_x_yzzz_xyz, g_y_x_yzzz_xz, g_y_x_yzzz_xzz, g_y_x_yzzz_yy, g_y_x_yzzz_yyz, g_y_x_yzzz_yz, g_y_x_yzzz_yzz, g_y_x_yzzz_zz, g_y_x_yzzz_zzz, g_y_x_yzzzz_xx, g_y_x_yzzzz_xy, g_y_x_yzzzz_xz, g_y_x_yzzzz_yy, g_y_x_yzzzz_yz, g_y_x_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzzzz_xx[k] = -g_y_x_yzzz_xx[k] * ab_z + g_y_x_yzzz_xxz[k];

                g_y_x_yzzzz_xy[k] = -g_y_x_yzzz_xy[k] * ab_z + g_y_x_yzzz_xyz[k];

                g_y_x_yzzzz_xz[k] = -g_y_x_yzzz_xz[k] * ab_z + g_y_x_yzzz_xzz[k];

                g_y_x_yzzzz_yy[k] = -g_y_x_yzzz_yy[k] * ab_z + g_y_x_yzzz_yyz[k];

                g_y_x_yzzzz_yz[k] = -g_y_x_yzzz_yz[k] * ab_z + g_y_x_yzzz_yzz[k];

                g_y_x_yzzzz_zz[k] = -g_y_x_yzzz_zz[k] * ab_z + g_y_x_yzzz_zzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzzzz_xx = cbuffer.data(hd_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_x_zzzzz_xy = cbuffer.data(hd_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_x_zzzzz_xz = cbuffer.data(hd_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_x_zzzzz_yy = cbuffer.data(hd_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_x_zzzzz_yz = cbuffer.data(hd_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_x_zzzzz_zz = cbuffer.data(hd_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zzzz_xx, g_y_x_zzzz_xxz, g_y_x_zzzz_xy, g_y_x_zzzz_xyz, g_y_x_zzzz_xz, g_y_x_zzzz_xzz, g_y_x_zzzz_yy, g_y_x_zzzz_yyz, g_y_x_zzzz_yz, g_y_x_zzzz_yzz, g_y_x_zzzz_zz, g_y_x_zzzz_zzz, g_y_x_zzzzz_xx, g_y_x_zzzzz_xy, g_y_x_zzzzz_xz, g_y_x_zzzzz_yy, g_y_x_zzzzz_yz, g_y_x_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzzzz_xx[k] = -g_y_x_zzzz_xx[k] * ab_z + g_y_x_zzzz_xxz[k];

                g_y_x_zzzzz_xy[k] = -g_y_x_zzzz_xy[k] * ab_z + g_y_x_zzzz_xyz[k];

                g_y_x_zzzzz_xz[k] = -g_y_x_zzzz_xz[k] * ab_z + g_y_x_zzzz_xzz[k];

                g_y_x_zzzzz_yy[k] = -g_y_x_zzzz_yy[k] * ab_z + g_y_x_zzzz_yyz[k];

                g_y_x_zzzzz_yz[k] = -g_y_x_zzzz_yz[k] * ab_z + g_y_x_zzzz_yzz[k];

                g_y_x_zzzzz_zz[k] = -g_y_x_zzzz_zz[k] * ab_z + g_y_x_zzzz_zzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxx_xx = cbuffer.data(hd_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_y_xxxxx_xy = cbuffer.data(hd_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_y_xxxxx_xz = cbuffer.data(hd_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_y_xxxxx_yy = cbuffer.data(hd_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_y_xxxxx_yz = cbuffer.data(hd_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_y_xxxxx_zz = cbuffer.data(hd_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxx_xx, g_y_y_xxxx_xxx, g_y_y_xxxx_xxy, g_y_y_xxxx_xxz, g_y_y_xxxx_xy, g_y_y_xxxx_xyy, g_y_y_xxxx_xyz, g_y_y_xxxx_xz, g_y_y_xxxx_xzz, g_y_y_xxxx_yy, g_y_y_xxxx_yz, g_y_y_xxxx_zz, g_y_y_xxxxx_xx, g_y_y_xxxxx_xy, g_y_y_xxxxx_xz, g_y_y_xxxxx_yy, g_y_y_xxxxx_yz, g_y_y_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxx_xx[k] = -g_y_y_xxxx_xx[k] * ab_x + g_y_y_xxxx_xxx[k];

                g_y_y_xxxxx_xy[k] = -g_y_y_xxxx_xy[k] * ab_x + g_y_y_xxxx_xxy[k];

                g_y_y_xxxxx_xz[k] = -g_y_y_xxxx_xz[k] * ab_x + g_y_y_xxxx_xxz[k];

                g_y_y_xxxxx_yy[k] = -g_y_y_xxxx_yy[k] * ab_x + g_y_y_xxxx_xyy[k];

                g_y_y_xxxxx_yz[k] = -g_y_y_xxxx_yz[k] * ab_x + g_y_y_xxxx_xyz[k];

                g_y_y_xxxxx_zz[k] = -g_y_y_xxxx_zz[k] * ab_x + g_y_y_xxxx_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxy_xx = cbuffer.data(hd_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_y_xxxxy_xy = cbuffer.data(hd_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_y_xxxxy_xz = cbuffer.data(hd_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_y_xxxxy_yy = cbuffer.data(hd_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_y_xxxxy_yz = cbuffer.data(hd_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_y_xxxxy_zz = cbuffer.data(hd_geom_11_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxxy_xx, g_y_y_xxxxy_xy, g_y_y_xxxxy_xz, g_y_y_xxxxy_yy, g_y_y_xxxxy_yz, g_y_y_xxxxy_zz, g_y_y_xxxy_xx, g_y_y_xxxy_xxx, g_y_y_xxxy_xxy, g_y_y_xxxy_xxz, g_y_y_xxxy_xy, g_y_y_xxxy_xyy, g_y_y_xxxy_xyz, g_y_y_xxxy_xz, g_y_y_xxxy_xzz, g_y_y_xxxy_yy, g_y_y_xxxy_yz, g_y_y_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxy_xx[k] = -g_y_y_xxxy_xx[k] * ab_x + g_y_y_xxxy_xxx[k];

                g_y_y_xxxxy_xy[k] = -g_y_y_xxxy_xy[k] * ab_x + g_y_y_xxxy_xxy[k];

                g_y_y_xxxxy_xz[k] = -g_y_y_xxxy_xz[k] * ab_x + g_y_y_xxxy_xxz[k];

                g_y_y_xxxxy_yy[k] = -g_y_y_xxxy_yy[k] * ab_x + g_y_y_xxxy_xyy[k];

                g_y_y_xxxxy_yz[k] = -g_y_y_xxxy_yz[k] * ab_x + g_y_y_xxxy_xyz[k];

                g_y_y_xxxxy_zz[k] = -g_y_y_xxxy_zz[k] * ab_x + g_y_y_xxxy_xzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxxz_xx = cbuffer.data(hd_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_y_xxxxz_xy = cbuffer.data(hd_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_y_xxxxz_xz = cbuffer.data(hd_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_y_xxxxz_yy = cbuffer.data(hd_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_y_xxxxz_yz = cbuffer.data(hd_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_y_xxxxz_zz = cbuffer.data(hd_geom_11_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxxz_xx, g_y_y_xxxxz_xy, g_y_y_xxxxz_xz, g_y_y_xxxxz_yy, g_y_y_xxxxz_yz, g_y_y_xxxxz_zz, g_y_y_xxxz_xx, g_y_y_xxxz_xxx, g_y_y_xxxz_xxy, g_y_y_xxxz_xxz, g_y_y_xxxz_xy, g_y_y_xxxz_xyy, g_y_y_xxxz_xyz, g_y_y_xxxz_xz, g_y_y_xxxz_xzz, g_y_y_xxxz_yy, g_y_y_xxxz_yz, g_y_y_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxxz_xx[k] = -g_y_y_xxxz_xx[k] * ab_x + g_y_y_xxxz_xxx[k];

                g_y_y_xxxxz_xy[k] = -g_y_y_xxxz_xy[k] * ab_x + g_y_y_xxxz_xxy[k];

                g_y_y_xxxxz_xz[k] = -g_y_y_xxxz_xz[k] * ab_x + g_y_y_xxxz_xxz[k];

                g_y_y_xxxxz_yy[k] = -g_y_y_xxxz_yy[k] * ab_x + g_y_y_xxxz_xyy[k];

                g_y_y_xxxxz_yz[k] = -g_y_y_xxxz_yz[k] * ab_x + g_y_y_xxxz_xyz[k];

                g_y_y_xxxxz_zz[k] = -g_y_y_xxxz_zz[k] * ab_x + g_y_y_xxxz_xzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxyy_xx = cbuffer.data(hd_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_y_xxxyy_xy = cbuffer.data(hd_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_y_xxxyy_xz = cbuffer.data(hd_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_y_xxxyy_yy = cbuffer.data(hd_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_y_xxxyy_yz = cbuffer.data(hd_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_y_xxxyy_zz = cbuffer.data(hd_geom_11_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxyy_xx, g_y_y_xxxyy_xy, g_y_y_xxxyy_xz, g_y_y_xxxyy_yy, g_y_y_xxxyy_yz, g_y_y_xxxyy_zz, g_y_y_xxyy_xx, g_y_y_xxyy_xxx, g_y_y_xxyy_xxy, g_y_y_xxyy_xxz, g_y_y_xxyy_xy, g_y_y_xxyy_xyy, g_y_y_xxyy_xyz, g_y_y_xxyy_xz, g_y_y_xxyy_xzz, g_y_y_xxyy_yy, g_y_y_xxyy_yz, g_y_y_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxyy_xx[k] = -g_y_y_xxyy_xx[k] * ab_x + g_y_y_xxyy_xxx[k];

                g_y_y_xxxyy_xy[k] = -g_y_y_xxyy_xy[k] * ab_x + g_y_y_xxyy_xxy[k];

                g_y_y_xxxyy_xz[k] = -g_y_y_xxyy_xz[k] * ab_x + g_y_y_xxyy_xxz[k];

                g_y_y_xxxyy_yy[k] = -g_y_y_xxyy_yy[k] * ab_x + g_y_y_xxyy_xyy[k];

                g_y_y_xxxyy_yz[k] = -g_y_y_xxyy_yz[k] * ab_x + g_y_y_xxyy_xyz[k];

                g_y_y_xxxyy_zz[k] = -g_y_y_xxyy_zz[k] * ab_x + g_y_y_xxyy_xzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxyz_xx = cbuffer.data(hd_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_y_xxxyz_xy = cbuffer.data(hd_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_y_xxxyz_xz = cbuffer.data(hd_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_y_xxxyz_yy = cbuffer.data(hd_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_y_xxxyz_yz = cbuffer.data(hd_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_y_xxxyz_zz = cbuffer.data(hd_geom_11_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxyz_xx, g_y_y_xxxyz_xy, g_y_y_xxxyz_xz, g_y_y_xxxyz_yy, g_y_y_xxxyz_yz, g_y_y_xxxyz_zz, g_y_y_xxyz_xx, g_y_y_xxyz_xxx, g_y_y_xxyz_xxy, g_y_y_xxyz_xxz, g_y_y_xxyz_xy, g_y_y_xxyz_xyy, g_y_y_xxyz_xyz, g_y_y_xxyz_xz, g_y_y_xxyz_xzz, g_y_y_xxyz_yy, g_y_y_xxyz_yz, g_y_y_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxyz_xx[k] = -g_y_y_xxyz_xx[k] * ab_x + g_y_y_xxyz_xxx[k];

                g_y_y_xxxyz_xy[k] = -g_y_y_xxyz_xy[k] * ab_x + g_y_y_xxyz_xxy[k];

                g_y_y_xxxyz_xz[k] = -g_y_y_xxyz_xz[k] * ab_x + g_y_y_xxyz_xxz[k];

                g_y_y_xxxyz_yy[k] = -g_y_y_xxyz_yy[k] * ab_x + g_y_y_xxyz_xyy[k];

                g_y_y_xxxyz_yz[k] = -g_y_y_xxyz_yz[k] * ab_x + g_y_y_xxyz_xyz[k];

                g_y_y_xxxyz_zz[k] = -g_y_y_xxyz_zz[k] * ab_x + g_y_y_xxyz_xzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxzz_xx = cbuffer.data(hd_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_y_xxxzz_xy = cbuffer.data(hd_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_y_xxxzz_xz = cbuffer.data(hd_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_y_xxxzz_yy = cbuffer.data(hd_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_y_xxxzz_yz = cbuffer.data(hd_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_y_xxxzz_zz = cbuffer.data(hd_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxzz_xx, g_y_y_xxxzz_xy, g_y_y_xxxzz_xz, g_y_y_xxxzz_yy, g_y_y_xxxzz_yz, g_y_y_xxxzz_zz, g_y_y_xxzz_xx, g_y_y_xxzz_xxx, g_y_y_xxzz_xxy, g_y_y_xxzz_xxz, g_y_y_xxzz_xy, g_y_y_xxzz_xyy, g_y_y_xxzz_xyz, g_y_y_xxzz_xz, g_y_y_xxzz_xzz, g_y_y_xxzz_yy, g_y_y_xxzz_yz, g_y_y_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxzz_xx[k] = -g_y_y_xxzz_xx[k] * ab_x + g_y_y_xxzz_xxx[k];

                g_y_y_xxxzz_xy[k] = -g_y_y_xxzz_xy[k] * ab_x + g_y_y_xxzz_xxy[k];

                g_y_y_xxxzz_xz[k] = -g_y_y_xxzz_xz[k] * ab_x + g_y_y_xxzz_xxz[k];

                g_y_y_xxxzz_yy[k] = -g_y_y_xxzz_yy[k] * ab_x + g_y_y_xxzz_xyy[k];

                g_y_y_xxxzz_yz[k] = -g_y_y_xxzz_yz[k] * ab_x + g_y_y_xxzz_xyz[k];

                g_y_y_xxxzz_zz[k] = -g_y_y_xxzz_zz[k] * ab_x + g_y_y_xxzz_xzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyyy_xx = cbuffer.data(hd_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_y_xxyyy_xy = cbuffer.data(hd_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_y_xxyyy_xz = cbuffer.data(hd_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_y_xxyyy_yy = cbuffer.data(hd_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_y_xxyyy_yz = cbuffer.data(hd_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_y_xxyyy_zz = cbuffer.data(hd_geom_11_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyyy_xx, g_y_y_xxyyy_xy, g_y_y_xxyyy_xz, g_y_y_xxyyy_yy, g_y_y_xxyyy_yz, g_y_y_xxyyy_zz, g_y_y_xyyy_xx, g_y_y_xyyy_xxx, g_y_y_xyyy_xxy, g_y_y_xyyy_xxz, g_y_y_xyyy_xy, g_y_y_xyyy_xyy, g_y_y_xyyy_xyz, g_y_y_xyyy_xz, g_y_y_xyyy_xzz, g_y_y_xyyy_yy, g_y_y_xyyy_yz, g_y_y_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyyy_xx[k] = -g_y_y_xyyy_xx[k] * ab_x + g_y_y_xyyy_xxx[k];

                g_y_y_xxyyy_xy[k] = -g_y_y_xyyy_xy[k] * ab_x + g_y_y_xyyy_xxy[k];

                g_y_y_xxyyy_xz[k] = -g_y_y_xyyy_xz[k] * ab_x + g_y_y_xyyy_xxz[k];

                g_y_y_xxyyy_yy[k] = -g_y_y_xyyy_yy[k] * ab_x + g_y_y_xyyy_xyy[k];

                g_y_y_xxyyy_yz[k] = -g_y_y_xyyy_yz[k] * ab_x + g_y_y_xyyy_xyz[k];

                g_y_y_xxyyy_zz[k] = -g_y_y_xyyy_zz[k] * ab_x + g_y_y_xyyy_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyyz_xx = cbuffer.data(hd_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_y_xxyyz_xy = cbuffer.data(hd_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_y_xxyyz_xz = cbuffer.data(hd_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_y_xxyyz_yy = cbuffer.data(hd_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_y_xxyyz_yz = cbuffer.data(hd_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_y_xxyyz_zz = cbuffer.data(hd_geom_11_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyyz_xx, g_y_y_xxyyz_xy, g_y_y_xxyyz_xz, g_y_y_xxyyz_yy, g_y_y_xxyyz_yz, g_y_y_xxyyz_zz, g_y_y_xyyz_xx, g_y_y_xyyz_xxx, g_y_y_xyyz_xxy, g_y_y_xyyz_xxz, g_y_y_xyyz_xy, g_y_y_xyyz_xyy, g_y_y_xyyz_xyz, g_y_y_xyyz_xz, g_y_y_xyyz_xzz, g_y_y_xyyz_yy, g_y_y_xyyz_yz, g_y_y_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyyz_xx[k] = -g_y_y_xyyz_xx[k] * ab_x + g_y_y_xyyz_xxx[k];

                g_y_y_xxyyz_xy[k] = -g_y_y_xyyz_xy[k] * ab_x + g_y_y_xyyz_xxy[k];

                g_y_y_xxyyz_xz[k] = -g_y_y_xyyz_xz[k] * ab_x + g_y_y_xyyz_xxz[k];

                g_y_y_xxyyz_yy[k] = -g_y_y_xyyz_yy[k] * ab_x + g_y_y_xyyz_xyy[k];

                g_y_y_xxyyz_yz[k] = -g_y_y_xyyz_yz[k] * ab_x + g_y_y_xyyz_xyz[k];

                g_y_y_xxyyz_zz[k] = -g_y_y_xyyz_zz[k] * ab_x + g_y_y_xyyz_xzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyzz_xx = cbuffer.data(hd_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_y_xxyzz_xy = cbuffer.data(hd_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_y_xxyzz_xz = cbuffer.data(hd_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_y_xxyzz_yy = cbuffer.data(hd_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_y_xxyzz_yz = cbuffer.data(hd_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_y_xxyzz_zz = cbuffer.data(hd_geom_11_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyzz_xx, g_y_y_xxyzz_xy, g_y_y_xxyzz_xz, g_y_y_xxyzz_yy, g_y_y_xxyzz_yz, g_y_y_xxyzz_zz, g_y_y_xyzz_xx, g_y_y_xyzz_xxx, g_y_y_xyzz_xxy, g_y_y_xyzz_xxz, g_y_y_xyzz_xy, g_y_y_xyzz_xyy, g_y_y_xyzz_xyz, g_y_y_xyzz_xz, g_y_y_xyzz_xzz, g_y_y_xyzz_yy, g_y_y_xyzz_yz, g_y_y_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyzz_xx[k] = -g_y_y_xyzz_xx[k] * ab_x + g_y_y_xyzz_xxx[k];

                g_y_y_xxyzz_xy[k] = -g_y_y_xyzz_xy[k] * ab_x + g_y_y_xyzz_xxy[k];

                g_y_y_xxyzz_xz[k] = -g_y_y_xyzz_xz[k] * ab_x + g_y_y_xyzz_xxz[k];

                g_y_y_xxyzz_yy[k] = -g_y_y_xyzz_yy[k] * ab_x + g_y_y_xyzz_xyy[k];

                g_y_y_xxyzz_yz[k] = -g_y_y_xyzz_yz[k] * ab_x + g_y_y_xyzz_xyz[k];

                g_y_y_xxyzz_zz[k] = -g_y_y_xyzz_zz[k] * ab_x + g_y_y_xyzz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxzzz_xx = cbuffer.data(hd_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_y_xxzzz_xy = cbuffer.data(hd_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_y_xxzzz_xz = cbuffer.data(hd_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_y_xxzzz_yy = cbuffer.data(hd_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_y_xxzzz_yz = cbuffer.data(hd_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_y_xxzzz_zz = cbuffer.data(hd_geom_11_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxzzz_xx, g_y_y_xxzzz_xy, g_y_y_xxzzz_xz, g_y_y_xxzzz_yy, g_y_y_xxzzz_yz, g_y_y_xxzzz_zz, g_y_y_xzzz_xx, g_y_y_xzzz_xxx, g_y_y_xzzz_xxy, g_y_y_xzzz_xxz, g_y_y_xzzz_xy, g_y_y_xzzz_xyy, g_y_y_xzzz_xyz, g_y_y_xzzz_xz, g_y_y_xzzz_xzz, g_y_y_xzzz_yy, g_y_y_xzzz_yz, g_y_y_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxzzz_xx[k] = -g_y_y_xzzz_xx[k] * ab_x + g_y_y_xzzz_xxx[k];

                g_y_y_xxzzz_xy[k] = -g_y_y_xzzz_xy[k] * ab_x + g_y_y_xzzz_xxy[k];

                g_y_y_xxzzz_xz[k] = -g_y_y_xzzz_xz[k] * ab_x + g_y_y_xzzz_xxz[k];

                g_y_y_xxzzz_yy[k] = -g_y_y_xzzz_yy[k] * ab_x + g_y_y_xzzz_xyy[k];

                g_y_y_xxzzz_yz[k] = -g_y_y_xzzz_yz[k] * ab_x + g_y_y_xzzz_xyz[k];

                g_y_y_xxzzz_zz[k] = -g_y_y_xzzz_zz[k] * ab_x + g_y_y_xzzz_xzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyyy_xx = cbuffer.data(hd_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_y_xyyyy_xy = cbuffer.data(hd_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_y_xyyyy_xz = cbuffer.data(hd_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_y_xyyyy_yy = cbuffer.data(hd_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_y_xyyyy_yz = cbuffer.data(hd_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_y_xyyyy_zz = cbuffer.data(hd_geom_11_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyyy_xx, g_y_y_xyyyy_xy, g_y_y_xyyyy_xz, g_y_y_xyyyy_yy, g_y_y_xyyyy_yz, g_y_y_xyyyy_zz, g_y_y_yyyy_xx, g_y_y_yyyy_xxx, g_y_y_yyyy_xxy, g_y_y_yyyy_xxz, g_y_y_yyyy_xy, g_y_y_yyyy_xyy, g_y_y_yyyy_xyz, g_y_y_yyyy_xz, g_y_y_yyyy_xzz, g_y_y_yyyy_yy, g_y_y_yyyy_yz, g_y_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyyy_xx[k] = -g_y_y_yyyy_xx[k] * ab_x + g_y_y_yyyy_xxx[k];

                g_y_y_xyyyy_xy[k] = -g_y_y_yyyy_xy[k] * ab_x + g_y_y_yyyy_xxy[k];

                g_y_y_xyyyy_xz[k] = -g_y_y_yyyy_xz[k] * ab_x + g_y_y_yyyy_xxz[k];

                g_y_y_xyyyy_yy[k] = -g_y_y_yyyy_yy[k] * ab_x + g_y_y_yyyy_xyy[k];

                g_y_y_xyyyy_yz[k] = -g_y_y_yyyy_yz[k] * ab_x + g_y_y_yyyy_xyz[k];

                g_y_y_xyyyy_zz[k] = -g_y_y_yyyy_zz[k] * ab_x + g_y_y_yyyy_xzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyyz_xx = cbuffer.data(hd_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_y_xyyyz_xy = cbuffer.data(hd_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_y_xyyyz_xz = cbuffer.data(hd_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_y_xyyyz_yy = cbuffer.data(hd_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_y_xyyyz_yz = cbuffer.data(hd_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_y_xyyyz_zz = cbuffer.data(hd_geom_11_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyyz_xx, g_y_y_xyyyz_xy, g_y_y_xyyyz_xz, g_y_y_xyyyz_yy, g_y_y_xyyyz_yz, g_y_y_xyyyz_zz, g_y_y_yyyz_xx, g_y_y_yyyz_xxx, g_y_y_yyyz_xxy, g_y_y_yyyz_xxz, g_y_y_yyyz_xy, g_y_y_yyyz_xyy, g_y_y_yyyz_xyz, g_y_y_yyyz_xz, g_y_y_yyyz_xzz, g_y_y_yyyz_yy, g_y_y_yyyz_yz, g_y_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyyz_xx[k] = -g_y_y_yyyz_xx[k] * ab_x + g_y_y_yyyz_xxx[k];

                g_y_y_xyyyz_xy[k] = -g_y_y_yyyz_xy[k] * ab_x + g_y_y_yyyz_xxy[k];

                g_y_y_xyyyz_xz[k] = -g_y_y_yyyz_xz[k] * ab_x + g_y_y_yyyz_xxz[k];

                g_y_y_xyyyz_yy[k] = -g_y_y_yyyz_yy[k] * ab_x + g_y_y_yyyz_xyy[k];

                g_y_y_xyyyz_yz[k] = -g_y_y_yyyz_yz[k] * ab_x + g_y_y_yyyz_xyz[k];

                g_y_y_xyyyz_zz[k] = -g_y_y_yyyz_zz[k] * ab_x + g_y_y_yyyz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyzz_xx = cbuffer.data(hd_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_y_xyyzz_xy = cbuffer.data(hd_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_y_xyyzz_xz = cbuffer.data(hd_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_y_xyyzz_yy = cbuffer.data(hd_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_y_xyyzz_yz = cbuffer.data(hd_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_y_xyyzz_zz = cbuffer.data(hd_geom_11_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyzz_xx, g_y_y_xyyzz_xy, g_y_y_xyyzz_xz, g_y_y_xyyzz_yy, g_y_y_xyyzz_yz, g_y_y_xyyzz_zz, g_y_y_yyzz_xx, g_y_y_yyzz_xxx, g_y_y_yyzz_xxy, g_y_y_yyzz_xxz, g_y_y_yyzz_xy, g_y_y_yyzz_xyy, g_y_y_yyzz_xyz, g_y_y_yyzz_xz, g_y_y_yyzz_xzz, g_y_y_yyzz_yy, g_y_y_yyzz_yz, g_y_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyzz_xx[k] = -g_y_y_yyzz_xx[k] * ab_x + g_y_y_yyzz_xxx[k];

                g_y_y_xyyzz_xy[k] = -g_y_y_yyzz_xy[k] * ab_x + g_y_y_yyzz_xxy[k];

                g_y_y_xyyzz_xz[k] = -g_y_y_yyzz_xz[k] * ab_x + g_y_y_yyzz_xxz[k];

                g_y_y_xyyzz_yy[k] = -g_y_y_yyzz_yy[k] * ab_x + g_y_y_yyzz_xyy[k];

                g_y_y_xyyzz_yz[k] = -g_y_y_yyzz_yz[k] * ab_x + g_y_y_yyzz_xyz[k];

                g_y_y_xyyzz_zz[k] = -g_y_y_yyzz_zz[k] * ab_x + g_y_y_yyzz_xzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyzzz_xx = cbuffer.data(hd_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_y_xyzzz_xy = cbuffer.data(hd_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_y_xyzzz_xz = cbuffer.data(hd_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_y_xyzzz_yy = cbuffer.data(hd_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_y_xyzzz_yz = cbuffer.data(hd_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_y_xyzzz_zz = cbuffer.data(hd_geom_11_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyzzz_xx, g_y_y_xyzzz_xy, g_y_y_xyzzz_xz, g_y_y_xyzzz_yy, g_y_y_xyzzz_yz, g_y_y_xyzzz_zz, g_y_y_yzzz_xx, g_y_y_yzzz_xxx, g_y_y_yzzz_xxy, g_y_y_yzzz_xxz, g_y_y_yzzz_xy, g_y_y_yzzz_xyy, g_y_y_yzzz_xyz, g_y_y_yzzz_xz, g_y_y_yzzz_xzz, g_y_y_yzzz_yy, g_y_y_yzzz_yz, g_y_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyzzz_xx[k] = -g_y_y_yzzz_xx[k] * ab_x + g_y_y_yzzz_xxx[k];

                g_y_y_xyzzz_xy[k] = -g_y_y_yzzz_xy[k] * ab_x + g_y_y_yzzz_xxy[k];

                g_y_y_xyzzz_xz[k] = -g_y_y_yzzz_xz[k] * ab_x + g_y_y_yzzz_xxz[k];

                g_y_y_xyzzz_yy[k] = -g_y_y_yzzz_yy[k] * ab_x + g_y_y_yzzz_xyy[k];

                g_y_y_xyzzz_yz[k] = -g_y_y_yzzz_yz[k] * ab_x + g_y_y_yzzz_xyz[k];

                g_y_y_xyzzz_zz[k] = -g_y_y_yzzz_zz[k] * ab_x + g_y_y_yzzz_xzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzzzz_xx = cbuffer.data(hd_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_y_xzzzz_xy = cbuffer.data(hd_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_y_xzzzz_xz = cbuffer.data(hd_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_y_xzzzz_yy = cbuffer.data(hd_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_y_xzzzz_yz = cbuffer.data(hd_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_y_xzzzz_zz = cbuffer.data(hd_geom_11_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzzzz_xx, g_y_y_xzzzz_xy, g_y_y_xzzzz_xz, g_y_y_xzzzz_yy, g_y_y_xzzzz_yz, g_y_y_xzzzz_zz, g_y_y_zzzz_xx, g_y_y_zzzz_xxx, g_y_y_zzzz_xxy, g_y_y_zzzz_xxz, g_y_y_zzzz_xy, g_y_y_zzzz_xyy, g_y_y_zzzz_xyz, g_y_y_zzzz_xz, g_y_y_zzzz_xzz, g_y_y_zzzz_yy, g_y_y_zzzz_yz, g_y_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzzzz_xx[k] = -g_y_y_zzzz_xx[k] * ab_x + g_y_y_zzzz_xxx[k];

                g_y_y_xzzzz_xy[k] = -g_y_y_zzzz_xy[k] * ab_x + g_y_y_zzzz_xxy[k];

                g_y_y_xzzzz_xz[k] = -g_y_y_zzzz_xz[k] * ab_x + g_y_y_zzzz_xxz[k];

                g_y_y_xzzzz_yy[k] = -g_y_y_zzzz_yy[k] * ab_x + g_y_y_zzzz_xyy[k];

                g_y_y_xzzzz_yz[k] = -g_y_y_zzzz_yz[k] * ab_x + g_y_y_zzzz_xyz[k];

                g_y_y_xzzzz_zz[k] = -g_y_y_zzzz_zz[k] * ab_x + g_y_y_zzzz_xzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyyy_xx = cbuffer.data(hd_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_y_yyyyy_xy = cbuffer.data(hd_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_y_yyyyy_xz = cbuffer.data(hd_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_y_yyyyy_yy = cbuffer.data(hd_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_y_yyyyy_yz = cbuffer.data(hd_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_y_yyyyy_zz = cbuffer.data(hd_geom_11_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xx, g_0_y_yyyy_xy, g_0_y_yyyy_xz, g_0_y_yyyy_yy, g_0_y_yyyy_yz, g_0_y_yyyy_zz, g_y_0_yyyy_xx, g_y_0_yyyy_xy, g_y_0_yyyy_xz, g_y_0_yyyy_yy, g_y_0_yyyy_yz, g_y_0_yyyy_zz, g_y_y_yyyy_xx, g_y_y_yyyy_xxy, g_y_y_yyyy_xy, g_y_y_yyyy_xyy, g_y_y_yyyy_xyz, g_y_y_yyyy_xz, g_y_y_yyyy_yy, g_y_y_yyyy_yyy, g_y_y_yyyy_yyz, g_y_y_yyyy_yz, g_y_y_yyyy_yzz, g_y_y_yyyy_zz, g_y_y_yyyyy_xx, g_y_y_yyyyy_xy, g_y_y_yyyyy_xz, g_y_y_yyyyy_yy, g_y_y_yyyyy_yz, g_y_y_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyyy_xx[k] = -g_0_y_yyyy_xx[k] + g_y_0_yyyy_xx[k] - g_y_y_yyyy_xx[k] * ab_y + g_y_y_yyyy_xxy[k];

                g_y_y_yyyyy_xy[k] = -g_0_y_yyyy_xy[k] + g_y_0_yyyy_xy[k] - g_y_y_yyyy_xy[k] * ab_y + g_y_y_yyyy_xyy[k];

                g_y_y_yyyyy_xz[k] = -g_0_y_yyyy_xz[k] + g_y_0_yyyy_xz[k] - g_y_y_yyyy_xz[k] * ab_y + g_y_y_yyyy_xyz[k];

                g_y_y_yyyyy_yy[k] = -g_0_y_yyyy_yy[k] + g_y_0_yyyy_yy[k] - g_y_y_yyyy_yy[k] * ab_y + g_y_y_yyyy_yyy[k];

                g_y_y_yyyyy_yz[k] = -g_0_y_yyyy_yz[k] + g_y_0_yyyy_yz[k] - g_y_y_yyyy_yz[k] * ab_y + g_y_y_yyyy_yyz[k];

                g_y_y_yyyyy_zz[k] = -g_0_y_yyyy_zz[k] + g_y_0_yyyy_zz[k] - g_y_y_yyyy_zz[k] * ab_y + g_y_y_yyyy_yzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyyz_xx = cbuffer.data(hd_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_y_yyyyz_xy = cbuffer.data(hd_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_y_yyyyz_xz = cbuffer.data(hd_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_y_yyyyz_yy = cbuffer.data(hd_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_y_yyyyz_yz = cbuffer.data(hd_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_y_yyyyz_zz = cbuffer.data(hd_geom_11_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyyy_xx, g_y_y_yyyy_xxz, g_y_y_yyyy_xy, g_y_y_yyyy_xyz, g_y_y_yyyy_xz, g_y_y_yyyy_xzz, g_y_y_yyyy_yy, g_y_y_yyyy_yyz, g_y_y_yyyy_yz, g_y_y_yyyy_yzz, g_y_y_yyyy_zz, g_y_y_yyyy_zzz, g_y_y_yyyyz_xx, g_y_y_yyyyz_xy, g_y_y_yyyyz_xz, g_y_y_yyyyz_yy, g_y_y_yyyyz_yz, g_y_y_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyyz_xx[k] = -g_y_y_yyyy_xx[k] * ab_z + g_y_y_yyyy_xxz[k];

                g_y_y_yyyyz_xy[k] = -g_y_y_yyyy_xy[k] * ab_z + g_y_y_yyyy_xyz[k];

                g_y_y_yyyyz_xz[k] = -g_y_y_yyyy_xz[k] * ab_z + g_y_y_yyyy_xzz[k];

                g_y_y_yyyyz_yy[k] = -g_y_y_yyyy_yy[k] * ab_z + g_y_y_yyyy_yyz[k];

                g_y_y_yyyyz_yz[k] = -g_y_y_yyyy_yz[k] * ab_z + g_y_y_yyyy_yzz[k];

                g_y_y_yyyyz_zz[k] = -g_y_y_yyyy_zz[k] * ab_z + g_y_y_yyyy_zzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyzz_xx = cbuffer.data(hd_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_y_yyyzz_xy = cbuffer.data(hd_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_y_yyyzz_xz = cbuffer.data(hd_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_y_yyyzz_yy = cbuffer.data(hd_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_y_yyyzz_yz = cbuffer.data(hd_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_y_yyyzz_zz = cbuffer.data(hd_geom_11_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyyz_xx, g_y_y_yyyz_xxz, g_y_y_yyyz_xy, g_y_y_yyyz_xyz, g_y_y_yyyz_xz, g_y_y_yyyz_xzz, g_y_y_yyyz_yy, g_y_y_yyyz_yyz, g_y_y_yyyz_yz, g_y_y_yyyz_yzz, g_y_y_yyyz_zz, g_y_y_yyyz_zzz, g_y_y_yyyzz_xx, g_y_y_yyyzz_xy, g_y_y_yyyzz_xz, g_y_y_yyyzz_yy, g_y_y_yyyzz_yz, g_y_y_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyzz_xx[k] = -g_y_y_yyyz_xx[k] * ab_z + g_y_y_yyyz_xxz[k];

                g_y_y_yyyzz_xy[k] = -g_y_y_yyyz_xy[k] * ab_z + g_y_y_yyyz_xyz[k];

                g_y_y_yyyzz_xz[k] = -g_y_y_yyyz_xz[k] * ab_z + g_y_y_yyyz_xzz[k];

                g_y_y_yyyzz_yy[k] = -g_y_y_yyyz_yy[k] * ab_z + g_y_y_yyyz_yyz[k];

                g_y_y_yyyzz_yz[k] = -g_y_y_yyyz_yz[k] * ab_z + g_y_y_yyyz_yzz[k];

                g_y_y_yyyzz_zz[k] = -g_y_y_yyyz_zz[k] * ab_z + g_y_y_yyyz_zzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyzzz_xx = cbuffer.data(hd_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_y_yyzzz_xy = cbuffer.data(hd_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_y_yyzzz_xz = cbuffer.data(hd_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_y_yyzzz_yy = cbuffer.data(hd_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_y_yyzzz_yz = cbuffer.data(hd_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_y_yyzzz_zz = cbuffer.data(hd_geom_11_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyzz_xx, g_y_y_yyzz_xxz, g_y_y_yyzz_xy, g_y_y_yyzz_xyz, g_y_y_yyzz_xz, g_y_y_yyzz_xzz, g_y_y_yyzz_yy, g_y_y_yyzz_yyz, g_y_y_yyzz_yz, g_y_y_yyzz_yzz, g_y_y_yyzz_zz, g_y_y_yyzz_zzz, g_y_y_yyzzz_xx, g_y_y_yyzzz_xy, g_y_y_yyzzz_xz, g_y_y_yyzzz_yy, g_y_y_yyzzz_yz, g_y_y_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyzzz_xx[k] = -g_y_y_yyzz_xx[k] * ab_z + g_y_y_yyzz_xxz[k];

                g_y_y_yyzzz_xy[k] = -g_y_y_yyzz_xy[k] * ab_z + g_y_y_yyzz_xyz[k];

                g_y_y_yyzzz_xz[k] = -g_y_y_yyzz_xz[k] * ab_z + g_y_y_yyzz_xzz[k];

                g_y_y_yyzzz_yy[k] = -g_y_y_yyzz_yy[k] * ab_z + g_y_y_yyzz_yyz[k];

                g_y_y_yyzzz_yz[k] = -g_y_y_yyzz_yz[k] * ab_z + g_y_y_yyzz_yzz[k];

                g_y_y_yyzzz_zz[k] = -g_y_y_yyzz_zz[k] * ab_z + g_y_y_yyzz_zzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzzzz_xx = cbuffer.data(hd_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_y_yzzzz_xy = cbuffer.data(hd_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_y_yzzzz_xz = cbuffer.data(hd_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_y_yzzzz_yy = cbuffer.data(hd_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_y_yzzzz_yz = cbuffer.data(hd_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_y_yzzzz_zz = cbuffer.data(hd_geom_11_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yzzz_xx, g_y_y_yzzz_xxz, g_y_y_yzzz_xy, g_y_y_yzzz_xyz, g_y_y_yzzz_xz, g_y_y_yzzz_xzz, g_y_y_yzzz_yy, g_y_y_yzzz_yyz, g_y_y_yzzz_yz, g_y_y_yzzz_yzz, g_y_y_yzzz_zz, g_y_y_yzzz_zzz, g_y_y_yzzzz_xx, g_y_y_yzzzz_xy, g_y_y_yzzzz_xz, g_y_y_yzzzz_yy, g_y_y_yzzzz_yz, g_y_y_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzzzz_xx[k] = -g_y_y_yzzz_xx[k] * ab_z + g_y_y_yzzz_xxz[k];

                g_y_y_yzzzz_xy[k] = -g_y_y_yzzz_xy[k] * ab_z + g_y_y_yzzz_xyz[k];

                g_y_y_yzzzz_xz[k] = -g_y_y_yzzz_xz[k] * ab_z + g_y_y_yzzz_xzz[k];

                g_y_y_yzzzz_yy[k] = -g_y_y_yzzz_yy[k] * ab_z + g_y_y_yzzz_yyz[k];

                g_y_y_yzzzz_yz[k] = -g_y_y_yzzz_yz[k] * ab_z + g_y_y_yzzz_yzz[k];

                g_y_y_yzzzz_zz[k] = -g_y_y_yzzz_zz[k] * ab_z + g_y_y_yzzz_zzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzzzz_xx = cbuffer.data(hd_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_y_zzzzz_xy = cbuffer.data(hd_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_y_zzzzz_xz = cbuffer.data(hd_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_y_zzzzz_yy = cbuffer.data(hd_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_y_zzzzz_yz = cbuffer.data(hd_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_y_zzzzz_zz = cbuffer.data(hd_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zzzz_xx, g_y_y_zzzz_xxz, g_y_y_zzzz_xy, g_y_y_zzzz_xyz, g_y_y_zzzz_xz, g_y_y_zzzz_xzz, g_y_y_zzzz_yy, g_y_y_zzzz_yyz, g_y_y_zzzz_yz, g_y_y_zzzz_yzz, g_y_y_zzzz_zz, g_y_y_zzzz_zzz, g_y_y_zzzzz_xx, g_y_y_zzzzz_xy, g_y_y_zzzzz_xz, g_y_y_zzzzz_yy, g_y_y_zzzzz_yz, g_y_y_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzzzz_xx[k] = -g_y_y_zzzz_xx[k] * ab_z + g_y_y_zzzz_xxz[k];

                g_y_y_zzzzz_xy[k] = -g_y_y_zzzz_xy[k] * ab_z + g_y_y_zzzz_xyz[k];

                g_y_y_zzzzz_xz[k] = -g_y_y_zzzz_xz[k] * ab_z + g_y_y_zzzz_xzz[k];

                g_y_y_zzzzz_yy[k] = -g_y_y_zzzz_yy[k] * ab_z + g_y_y_zzzz_yyz[k];

                g_y_y_zzzzz_yz[k] = -g_y_y_zzzz_yz[k] * ab_z + g_y_y_zzzz_yzz[k];

                g_y_y_zzzzz_zz[k] = -g_y_y_zzzz_zz[k] * ab_z + g_y_y_zzzz_zzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxx_xx = cbuffer.data(hd_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_z_xxxxx_xy = cbuffer.data(hd_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_z_xxxxx_xz = cbuffer.data(hd_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_z_xxxxx_yy = cbuffer.data(hd_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_z_xxxxx_yz = cbuffer.data(hd_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_z_xxxxx_zz = cbuffer.data(hd_geom_11_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxx_xx, g_y_z_xxxx_xxx, g_y_z_xxxx_xxy, g_y_z_xxxx_xxz, g_y_z_xxxx_xy, g_y_z_xxxx_xyy, g_y_z_xxxx_xyz, g_y_z_xxxx_xz, g_y_z_xxxx_xzz, g_y_z_xxxx_yy, g_y_z_xxxx_yz, g_y_z_xxxx_zz, g_y_z_xxxxx_xx, g_y_z_xxxxx_xy, g_y_z_xxxxx_xz, g_y_z_xxxxx_yy, g_y_z_xxxxx_yz, g_y_z_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxx_xx[k] = -g_y_z_xxxx_xx[k] * ab_x + g_y_z_xxxx_xxx[k];

                g_y_z_xxxxx_xy[k] = -g_y_z_xxxx_xy[k] * ab_x + g_y_z_xxxx_xxy[k];

                g_y_z_xxxxx_xz[k] = -g_y_z_xxxx_xz[k] * ab_x + g_y_z_xxxx_xxz[k];

                g_y_z_xxxxx_yy[k] = -g_y_z_xxxx_yy[k] * ab_x + g_y_z_xxxx_xyy[k];

                g_y_z_xxxxx_yz[k] = -g_y_z_xxxx_yz[k] * ab_x + g_y_z_xxxx_xyz[k];

                g_y_z_xxxxx_zz[k] = -g_y_z_xxxx_zz[k] * ab_x + g_y_z_xxxx_xzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxy_xx = cbuffer.data(hd_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_z_xxxxy_xy = cbuffer.data(hd_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_z_xxxxy_xz = cbuffer.data(hd_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_z_xxxxy_yy = cbuffer.data(hd_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_z_xxxxy_yz = cbuffer.data(hd_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_z_xxxxy_zz = cbuffer.data(hd_geom_11_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxxy_xx, g_y_z_xxxxy_xy, g_y_z_xxxxy_xz, g_y_z_xxxxy_yy, g_y_z_xxxxy_yz, g_y_z_xxxxy_zz, g_y_z_xxxy_xx, g_y_z_xxxy_xxx, g_y_z_xxxy_xxy, g_y_z_xxxy_xxz, g_y_z_xxxy_xy, g_y_z_xxxy_xyy, g_y_z_xxxy_xyz, g_y_z_xxxy_xz, g_y_z_xxxy_xzz, g_y_z_xxxy_yy, g_y_z_xxxy_yz, g_y_z_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxy_xx[k] = -g_y_z_xxxy_xx[k] * ab_x + g_y_z_xxxy_xxx[k];

                g_y_z_xxxxy_xy[k] = -g_y_z_xxxy_xy[k] * ab_x + g_y_z_xxxy_xxy[k];

                g_y_z_xxxxy_xz[k] = -g_y_z_xxxy_xz[k] * ab_x + g_y_z_xxxy_xxz[k];

                g_y_z_xxxxy_yy[k] = -g_y_z_xxxy_yy[k] * ab_x + g_y_z_xxxy_xyy[k];

                g_y_z_xxxxy_yz[k] = -g_y_z_xxxy_yz[k] * ab_x + g_y_z_xxxy_xyz[k];

                g_y_z_xxxxy_zz[k] = -g_y_z_xxxy_zz[k] * ab_x + g_y_z_xxxy_xzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxxz_xx = cbuffer.data(hd_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_z_xxxxz_xy = cbuffer.data(hd_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_z_xxxxz_xz = cbuffer.data(hd_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_z_xxxxz_yy = cbuffer.data(hd_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_z_xxxxz_yz = cbuffer.data(hd_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_z_xxxxz_zz = cbuffer.data(hd_geom_11_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxxz_xx, g_y_z_xxxxz_xy, g_y_z_xxxxz_xz, g_y_z_xxxxz_yy, g_y_z_xxxxz_yz, g_y_z_xxxxz_zz, g_y_z_xxxz_xx, g_y_z_xxxz_xxx, g_y_z_xxxz_xxy, g_y_z_xxxz_xxz, g_y_z_xxxz_xy, g_y_z_xxxz_xyy, g_y_z_xxxz_xyz, g_y_z_xxxz_xz, g_y_z_xxxz_xzz, g_y_z_xxxz_yy, g_y_z_xxxz_yz, g_y_z_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxxz_xx[k] = -g_y_z_xxxz_xx[k] * ab_x + g_y_z_xxxz_xxx[k];

                g_y_z_xxxxz_xy[k] = -g_y_z_xxxz_xy[k] * ab_x + g_y_z_xxxz_xxy[k];

                g_y_z_xxxxz_xz[k] = -g_y_z_xxxz_xz[k] * ab_x + g_y_z_xxxz_xxz[k];

                g_y_z_xxxxz_yy[k] = -g_y_z_xxxz_yy[k] * ab_x + g_y_z_xxxz_xyy[k];

                g_y_z_xxxxz_yz[k] = -g_y_z_xxxz_yz[k] * ab_x + g_y_z_xxxz_xyz[k];

                g_y_z_xxxxz_zz[k] = -g_y_z_xxxz_zz[k] * ab_x + g_y_z_xxxz_xzz[k];
            }

            /// Set up 648-654 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxyy_xx = cbuffer.data(hd_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_z_xxxyy_xy = cbuffer.data(hd_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_z_xxxyy_xz = cbuffer.data(hd_geom_11_off + 650 * ccomps * dcomps);

            auto g_y_z_xxxyy_yy = cbuffer.data(hd_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_z_xxxyy_yz = cbuffer.data(hd_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_z_xxxyy_zz = cbuffer.data(hd_geom_11_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxyy_xx, g_y_z_xxxyy_xy, g_y_z_xxxyy_xz, g_y_z_xxxyy_yy, g_y_z_xxxyy_yz, g_y_z_xxxyy_zz, g_y_z_xxyy_xx, g_y_z_xxyy_xxx, g_y_z_xxyy_xxy, g_y_z_xxyy_xxz, g_y_z_xxyy_xy, g_y_z_xxyy_xyy, g_y_z_xxyy_xyz, g_y_z_xxyy_xz, g_y_z_xxyy_xzz, g_y_z_xxyy_yy, g_y_z_xxyy_yz, g_y_z_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxyy_xx[k] = -g_y_z_xxyy_xx[k] * ab_x + g_y_z_xxyy_xxx[k];

                g_y_z_xxxyy_xy[k] = -g_y_z_xxyy_xy[k] * ab_x + g_y_z_xxyy_xxy[k];

                g_y_z_xxxyy_xz[k] = -g_y_z_xxyy_xz[k] * ab_x + g_y_z_xxyy_xxz[k];

                g_y_z_xxxyy_yy[k] = -g_y_z_xxyy_yy[k] * ab_x + g_y_z_xxyy_xyy[k];

                g_y_z_xxxyy_yz[k] = -g_y_z_xxyy_yz[k] * ab_x + g_y_z_xxyy_xyz[k];

                g_y_z_xxxyy_zz[k] = -g_y_z_xxyy_zz[k] * ab_x + g_y_z_xxyy_xzz[k];
            }

            /// Set up 654-660 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxyz_xx = cbuffer.data(hd_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_z_xxxyz_xy = cbuffer.data(hd_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_z_xxxyz_xz = cbuffer.data(hd_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_z_xxxyz_yy = cbuffer.data(hd_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_z_xxxyz_yz = cbuffer.data(hd_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_z_xxxyz_zz = cbuffer.data(hd_geom_11_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxyz_xx, g_y_z_xxxyz_xy, g_y_z_xxxyz_xz, g_y_z_xxxyz_yy, g_y_z_xxxyz_yz, g_y_z_xxxyz_zz, g_y_z_xxyz_xx, g_y_z_xxyz_xxx, g_y_z_xxyz_xxy, g_y_z_xxyz_xxz, g_y_z_xxyz_xy, g_y_z_xxyz_xyy, g_y_z_xxyz_xyz, g_y_z_xxyz_xz, g_y_z_xxyz_xzz, g_y_z_xxyz_yy, g_y_z_xxyz_yz, g_y_z_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxyz_xx[k] = -g_y_z_xxyz_xx[k] * ab_x + g_y_z_xxyz_xxx[k];

                g_y_z_xxxyz_xy[k] = -g_y_z_xxyz_xy[k] * ab_x + g_y_z_xxyz_xxy[k];

                g_y_z_xxxyz_xz[k] = -g_y_z_xxyz_xz[k] * ab_x + g_y_z_xxyz_xxz[k];

                g_y_z_xxxyz_yy[k] = -g_y_z_xxyz_yy[k] * ab_x + g_y_z_xxyz_xyy[k];

                g_y_z_xxxyz_yz[k] = -g_y_z_xxyz_yz[k] * ab_x + g_y_z_xxyz_xyz[k];

                g_y_z_xxxyz_zz[k] = -g_y_z_xxyz_zz[k] * ab_x + g_y_z_xxyz_xzz[k];
            }

            /// Set up 660-666 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxzz_xx = cbuffer.data(hd_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_z_xxxzz_xy = cbuffer.data(hd_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_z_xxxzz_xz = cbuffer.data(hd_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_z_xxxzz_yy = cbuffer.data(hd_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_z_xxxzz_yz = cbuffer.data(hd_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_z_xxxzz_zz = cbuffer.data(hd_geom_11_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxzz_xx, g_y_z_xxxzz_xy, g_y_z_xxxzz_xz, g_y_z_xxxzz_yy, g_y_z_xxxzz_yz, g_y_z_xxxzz_zz, g_y_z_xxzz_xx, g_y_z_xxzz_xxx, g_y_z_xxzz_xxy, g_y_z_xxzz_xxz, g_y_z_xxzz_xy, g_y_z_xxzz_xyy, g_y_z_xxzz_xyz, g_y_z_xxzz_xz, g_y_z_xxzz_xzz, g_y_z_xxzz_yy, g_y_z_xxzz_yz, g_y_z_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxzz_xx[k] = -g_y_z_xxzz_xx[k] * ab_x + g_y_z_xxzz_xxx[k];

                g_y_z_xxxzz_xy[k] = -g_y_z_xxzz_xy[k] * ab_x + g_y_z_xxzz_xxy[k];

                g_y_z_xxxzz_xz[k] = -g_y_z_xxzz_xz[k] * ab_x + g_y_z_xxzz_xxz[k];

                g_y_z_xxxzz_yy[k] = -g_y_z_xxzz_yy[k] * ab_x + g_y_z_xxzz_xyy[k];

                g_y_z_xxxzz_yz[k] = -g_y_z_xxzz_yz[k] * ab_x + g_y_z_xxzz_xyz[k];

                g_y_z_xxxzz_zz[k] = -g_y_z_xxzz_zz[k] * ab_x + g_y_z_xxzz_xzz[k];
            }

            /// Set up 666-672 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyyy_xx = cbuffer.data(hd_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_z_xxyyy_xy = cbuffer.data(hd_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_z_xxyyy_xz = cbuffer.data(hd_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_z_xxyyy_yy = cbuffer.data(hd_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_z_xxyyy_yz = cbuffer.data(hd_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_z_xxyyy_zz = cbuffer.data(hd_geom_11_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyyy_xx, g_y_z_xxyyy_xy, g_y_z_xxyyy_xz, g_y_z_xxyyy_yy, g_y_z_xxyyy_yz, g_y_z_xxyyy_zz, g_y_z_xyyy_xx, g_y_z_xyyy_xxx, g_y_z_xyyy_xxy, g_y_z_xyyy_xxz, g_y_z_xyyy_xy, g_y_z_xyyy_xyy, g_y_z_xyyy_xyz, g_y_z_xyyy_xz, g_y_z_xyyy_xzz, g_y_z_xyyy_yy, g_y_z_xyyy_yz, g_y_z_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyyy_xx[k] = -g_y_z_xyyy_xx[k] * ab_x + g_y_z_xyyy_xxx[k];

                g_y_z_xxyyy_xy[k] = -g_y_z_xyyy_xy[k] * ab_x + g_y_z_xyyy_xxy[k];

                g_y_z_xxyyy_xz[k] = -g_y_z_xyyy_xz[k] * ab_x + g_y_z_xyyy_xxz[k];

                g_y_z_xxyyy_yy[k] = -g_y_z_xyyy_yy[k] * ab_x + g_y_z_xyyy_xyy[k];

                g_y_z_xxyyy_yz[k] = -g_y_z_xyyy_yz[k] * ab_x + g_y_z_xyyy_xyz[k];

                g_y_z_xxyyy_zz[k] = -g_y_z_xyyy_zz[k] * ab_x + g_y_z_xyyy_xzz[k];
            }

            /// Set up 672-678 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyyz_xx = cbuffer.data(hd_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_z_xxyyz_xy = cbuffer.data(hd_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_z_xxyyz_xz = cbuffer.data(hd_geom_11_off + 674 * ccomps * dcomps);

            auto g_y_z_xxyyz_yy = cbuffer.data(hd_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_z_xxyyz_yz = cbuffer.data(hd_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_z_xxyyz_zz = cbuffer.data(hd_geom_11_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyyz_xx, g_y_z_xxyyz_xy, g_y_z_xxyyz_xz, g_y_z_xxyyz_yy, g_y_z_xxyyz_yz, g_y_z_xxyyz_zz, g_y_z_xyyz_xx, g_y_z_xyyz_xxx, g_y_z_xyyz_xxy, g_y_z_xyyz_xxz, g_y_z_xyyz_xy, g_y_z_xyyz_xyy, g_y_z_xyyz_xyz, g_y_z_xyyz_xz, g_y_z_xyyz_xzz, g_y_z_xyyz_yy, g_y_z_xyyz_yz, g_y_z_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyyz_xx[k] = -g_y_z_xyyz_xx[k] * ab_x + g_y_z_xyyz_xxx[k];

                g_y_z_xxyyz_xy[k] = -g_y_z_xyyz_xy[k] * ab_x + g_y_z_xyyz_xxy[k];

                g_y_z_xxyyz_xz[k] = -g_y_z_xyyz_xz[k] * ab_x + g_y_z_xyyz_xxz[k];

                g_y_z_xxyyz_yy[k] = -g_y_z_xyyz_yy[k] * ab_x + g_y_z_xyyz_xyy[k];

                g_y_z_xxyyz_yz[k] = -g_y_z_xyyz_yz[k] * ab_x + g_y_z_xyyz_xyz[k];

                g_y_z_xxyyz_zz[k] = -g_y_z_xyyz_zz[k] * ab_x + g_y_z_xyyz_xzz[k];
            }

            /// Set up 678-684 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyzz_xx = cbuffer.data(hd_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_z_xxyzz_xy = cbuffer.data(hd_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_z_xxyzz_xz = cbuffer.data(hd_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_z_xxyzz_yy = cbuffer.data(hd_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_z_xxyzz_yz = cbuffer.data(hd_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_z_xxyzz_zz = cbuffer.data(hd_geom_11_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyzz_xx, g_y_z_xxyzz_xy, g_y_z_xxyzz_xz, g_y_z_xxyzz_yy, g_y_z_xxyzz_yz, g_y_z_xxyzz_zz, g_y_z_xyzz_xx, g_y_z_xyzz_xxx, g_y_z_xyzz_xxy, g_y_z_xyzz_xxz, g_y_z_xyzz_xy, g_y_z_xyzz_xyy, g_y_z_xyzz_xyz, g_y_z_xyzz_xz, g_y_z_xyzz_xzz, g_y_z_xyzz_yy, g_y_z_xyzz_yz, g_y_z_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyzz_xx[k] = -g_y_z_xyzz_xx[k] * ab_x + g_y_z_xyzz_xxx[k];

                g_y_z_xxyzz_xy[k] = -g_y_z_xyzz_xy[k] * ab_x + g_y_z_xyzz_xxy[k];

                g_y_z_xxyzz_xz[k] = -g_y_z_xyzz_xz[k] * ab_x + g_y_z_xyzz_xxz[k];

                g_y_z_xxyzz_yy[k] = -g_y_z_xyzz_yy[k] * ab_x + g_y_z_xyzz_xyy[k];

                g_y_z_xxyzz_yz[k] = -g_y_z_xyzz_yz[k] * ab_x + g_y_z_xyzz_xyz[k];

                g_y_z_xxyzz_zz[k] = -g_y_z_xyzz_zz[k] * ab_x + g_y_z_xyzz_xzz[k];
            }

            /// Set up 684-690 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxzzz_xx = cbuffer.data(hd_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_z_xxzzz_xy = cbuffer.data(hd_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_z_xxzzz_xz = cbuffer.data(hd_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_z_xxzzz_yy = cbuffer.data(hd_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_z_xxzzz_yz = cbuffer.data(hd_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_z_xxzzz_zz = cbuffer.data(hd_geom_11_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxzzz_xx, g_y_z_xxzzz_xy, g_y_z_xxzzz_xz, g_y_z_xxzzz_yy, g_y_z_xxzzz_yz, g_y_z_xxzzz_zz, g_y_z_xzzz_xx, g_y_z_xzzz_xxx, g_y_z_xzzz_xxy, g_y_z_xzzz_xxz, g_y_z_xzzz_xy, g_y_z_xzzz_xyy, g_y_z_xzzz_xyz, g_y_z_xzzz_xz, g_y_z_xzzz_xzz, g_y_z_xzzz_yy, g_y_z_xzzz_yz, g_y_z_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxzzz_xx[k] = -g_y_z_xzzz_xx[k] * ab_x + g_y_z_xzzz_xxx[k];

                g_y_z_xxzzz_xy[k] = -g_y_z_xzzz_xy[k] * ab_x + g_y_z_xzzz_xxy[k];

                g_y_z_xxzzz_xz[k] = -g_y_z_xzzz_xz[k] * ab_x + g_y_z_xzzz_xxz[k];

                g_y_z_xxzzz_yy[k] = -g_y_z_xzzz_yy[k] * ab_x + g_y_z_xzzz_xyy[k];

                g_y_z_xxzzz_yz[k] = -g_y_z_xzzz_yz[k] * ab_x + g_y_z_xzzz_xyz[k];

                g_y_z_xxzzz_zz[k] = -g_y_z_xzzz_zz[k] * ab_x + g_y_z_xzzz_xzz[k];
            }

            /// Set up 690-696 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyyy_xx = cbuffer.data(hd_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_z_xyyyy_xy = cbuffer.data(hd_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_z_xyyyy_xz = cbuffer.data(hd_geom_11_off + 692 * ccomps * dcomps);

            auto g_y_z_xyyyy_yy = cbuffer.data(hd_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_z_xyyyy_yz = cbuffer.data(hd_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_z_xyyyy_zz = cbuffer.data(hd_geom_11_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyyy_xx, g_y_z_xyyyy_xy, g_y_z_xyyyy_xz, g_y_z_xyyyy_yy, g_y_z_xyyyy_yz, g_y_z_xyyyy_zz, g_y_z_yyyy_xx, g_y_z_yyyy_xxx, g_y_z_yyyy_xxy, g_y_z_yyyy_xxz, g_y_z_yyyy_xy, g_y_z_yyyy_xyy, g_y_z_yyyy_xyz, g_y_z_yyyy_xz, g_y_z_yyyy_xzz, g_y_z_yyyy_yy, g_y_z_yyyy_yz, g_y_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyyy_xx[k] = -g_y_z_yyyy_xx[k] * ab_x + g_y_z_yyyy_xxx[k];

                g_y_z_xyyyy_xy[k] = -g_y_z_yyyy_xy[k] * ab_x + g_y_z_yyyy_xxy[k];

                g_y_z_xyyyy_xz[k] = -g_y_z_yyyy_xz[k] * ab_x + g_y_z_yyyy_xxz[k];

                g_y_z_xyyyy_yy[k] = -g_y_z_yyyy_yy[k] * ab_x + g_y_z_yyyy_xyy[k];

                g_y_z_xyyyy_yz[k] = -g_y_z_yyyy_yz[k] * ab_x + g_y_z_yyyy_xyz[k];

                g_y_z_xyyyy_zz[k] = -g_y_z_yyyy_zz[k] * ab_x + g_y_z_yyyy_xzz[k];
            }

            /// Set up 696-702 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyyz_xx = cbuffer.data(hd_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_z_xyyyz_xy = cbuffer.data(hd_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_z_xyyyz_xz = cbuffer.data(hd_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_z_xyyyz_yy = cbuffer.data(hd_geom_11_off + 699 * ccomps * dcomps);

            auto g_y_z_xyyyz_yz = cbuffer.data(hd_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_z_xyyyz_zz = cbuffer.data(hd_geom_11_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyyz_xx, g_y_z_xyyyz_xy, g_y_z_xyyyz_xz, g_y_z_xyyyz_yy, g_y_z_xyyyz_yz, g_y_z_xyyyz_zz, g_y_z_yyyz_xx, g_y_z_yyyz_xxx, g_y_z_yyyz_xxy, g_y_z_yyyz_xxz, g_y_z_yyyz_xy, g_y_z_yyyz_xyy, g_y_z_yyyz_xyz, g_y_z_yyyz_xz, g_y_z_yyyz_xzz, g_y_z_yyyz_yy, g_y_z_yyyz_yz, g_y_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyyz_xx[k] = -g_y_z_yyyz_xx[k] * ab_x + g_y_z_yyyz_xxx[k];

                g_y_z_xyyyz_xy[k] = -g_y_z_yyyz_xy[k] * ab_x + g_y_z_yyyz_xxy[k];

                g_y_z_xyyyz_xz[k] = -g_y_z_yyyz_xz[k] * ab_x + g_y_z_yyyz_xxz[k];

                g_y_z_xyyyz_yy[k] = -g_y_z_yyyz_yy[k] * ab_x + g_y_z_yyyz_xyy[k];

                g_y_z_xyyyz_yz[k] = -g_y_z_yyyz_yz[k] * ab_x + g_y_z_yyyz_xyz[k];

                g_y_z_xyyyz_zz[k] = -g_y_z_yyyz_zz[k] * ab_x + g_y_z_yyyz_xzz[k];
            }

            /// Set up 702-708 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyzz_xx = cbuffer.data(hd_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_z_xyyzz_xy = cbuffer.data(hd_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_z_xyyzz_xz = cbuffer.data(hd_geom_11_off + 704 * ccomps * dcomps);

            auto g_y_z_xyyzz_yy = cbuffer.data(hd_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_z_xyyzz_yz = cbuffer.data(hd_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_z_xyyzz_zz = cbuffer.data(hd_geom_11_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyzz_xx, g_y_z_xyyzz_xy, g_y_z_xyyzz_xz, g_y_z_xyyzz_yy, g_y_z_xyyzz_yz, g_y_z_xyyzz_zz, g_y_z_yyzz_xx, g_y_z_yyzz_xxx, g_y_z_yyzz_xxy, g_y_z_yyzz_xxz, g_y_z_yyzz_xy, g_y_z_yyzz_xyy, g_y_z_yyzz_xyz, g_y_z_yyzz_xz, g_y_z_yyzz_xzz, g_y_z_yyzz_yy, g_y_z_yyzz_yz, g_y_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyzz_xx[k] = -g_y_z_yyzz_xx[k] * ab_x + g_y_z_yyzz_xxx[k];

                g_y_z_xyyzz_xy[k] = -g_y_z_yyzz_xy[k] * ab_x + g_y_z_yyzz_xxy[k];

                g_y_z_xyyzz_xz[k] = -g_y_z_yyzz_xz[k] * ab_x + g_y_z_yyzz_xxz[k];

                g_y_z_xyyzz_yy[k] = -g_y_z_yyzz_yy[k] * ab_x + g_y_z_yyzz_xyy[k];

                g_y_z_xyyzz_yz[k] = -g_y_z_yyzz_yz[k] * ab_x + g_y_z_yyzz_xyz[k];

                g_y_z_xyyzz_zz[k] = -g_y_z_yyzz_zz[k] * ab_x + g_y_z_yyzz_xzz[k];
            }

            /// Set up 708-714 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyzzz_xx = cbuffer.data(hd_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_z_xyzzz_xy = cbuffer.data(hd_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_z_xyzzz_xz = cbuffer.data(hd_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_z_xyzzz_yy = cbuffer.data(hd_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_z_xyzzz_yz = cbuffer.data(hd_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_z_xyzzz_zz = cbuffer.data(hd_geom_11_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyzzz_xx, g_y_z_xyzzz_xy, g_y_z_xyzzz_xz, g_y_z_xyzzz_yy, g_y_z_xyzzz_yz, g_y_z_xyzzz_zz, g_y_z_yzzz_xx, g_y_z_yzzz_xxx, g_y_z_yzzz_xxy, g_y_z_yzzz_xxz, g_y_z_yzzz_xy, g_y_z_yzzz_xyy, g_y_z_yzzz_xyz, g_y_z_yzzz_xz, g_y_z_yzzz_xzz, g_y_z_yzzz_yy, g_y_z_yzzz_yz, g_y_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyzzz_xx[k] = -g_y_z_yzzz_xx[k] * ab_x + g_y_z_yzzz_xxx[k];

                g_y_z_xyzzz_xy[k] = -g_y_z_yzzz_xy[k] * ab_x + g_y_z_yzzz_xxy[k];

                g_y_z_xyzzz_xz[k] = -g_y_z_yzzz_xz[k] * ab_x + g_y_z_yzzz_xxz[k];

                g_y_z_xyzzz_yy[k] = -g_y_z_yzzz_yy[k] * ab_x + g_y_z_yzzz_xyy[k];

                g_y_z_xyzzz_yz[k] = -g_y_z_yzzz_yz[k] * ab_x + g_y_z_yzzz_xyz[k];

                g_y_z_xyzzz_zz[k] = -g_y_z_yzzz_zz[k] * ab_x + g_y_z_yzzz_xzz[k];
            }

            /// Set up 714-720 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzzzz_xx = cbuffer.data(hd_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_z_xzzzz_xy = cbuffer.data(hd_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_z_xzzzz_xz = cbuffer.data(hd_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_z_xzzzz_yy = cbuffer.data(hd_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_z_xzzzz_yz = cbuffer.data(hd_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_z_xzzzz_zz = cbuffer.data(hd_geom_11_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzzzz_xx, g_y_z_xzzzz_xy, g_y_z_xzzzz_xz, g_y_z_xzzzz_yy, g_y_z_xzzzz_yz, g_y_z_xzzzz_zz, g_y_z_zzzz_xx, g_y_z_zzzz_xxx, g_y_z_zzzz_xxy, g_y_z_zzzz_xxz, g_y_z_zzzz_xy, g_y_z_zzzz_xyy, g_y_z_zzzz_xyz, g_y_z_zzzz_xz, g_y_z_zzzz_xzz, g_y_z_zzzz_yy, g_y_z_zzzz_yz, g_y_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzzzz_xx[k] = -g_y_z_zzzz_xx[k] * ab_x + g_y_z_zzzz_xxx[k];

                g_y_z_xzzzz_xy[k] = -g_y_z_zzzz_xy[k] * ab_x + g_y_z_zzzz_xxy[k];

                g_y_z_xzzzz_xz[k] = -g_y_z_zzzz_xz[k] * ab_x + g_y_z_zzzz_xxz[k];

                g_y_z_xzzzz_yy[k] = -g_y_z_zzzz_yy[k] * ab_x + g_y_z_zzzz_xyy[k];

                g_y_z_xzzzz_yz[k] = -g_y_z_zzzz_yz[k] * ab_x + g_y_z_zzzz_xyz[k];

                g_y_z_xzzzz_zz[k] = -g_y_z_zzzz_zz[k] * ab_x + g_y_z_zzzz_xzz[k];
            }

            /// Set up 720-726 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyyy_xx = cbuffer.data(hd_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_z_yyyyy_xy = cbuffer.data(hd_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_z_yyyyy_xz = cbuffer.data(hd_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_z_yyyyy_yy = cbuffer.data(hd_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_z_yyyyy_yz = cbuffer.data(hd_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_z_yyyyy_zz = cbuffer.data(hd_geom_11_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_xx, g_0_z_yyyy_xy, g_0_z_yyyy_xz, g_0_z_yyyy_yy, g_0_z_yyyy_yz, g_0_z_yyyy_zz, g_y_z_yyyy_xx, g_y_z_yyyy_xxy, g_y_z_yyyy_xy, g_y_z_yyyy_xyy, g_y_z_yyyy_xyz, g_y_z_yyyy_xz, g_y_z_yyyy_yy, g_y_z_yyyy_yyy, g_y_z_yyyy_yyz, g_y_z_yyyy_yz, g_y_z_yyyy_yzz, g_y_z_yyyy_zz, g_y_z_yyyyy_xx, g_y_z_yyyyy_xy, g_y_z_yyyyy_xz, g_y_z_yyyyy_yy, g_y_z_yyyyy_yz, g_y_z_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyyy_xx[k] = -g_0_z_yyyy_xx[k] - g_y_z_yyyy_xx[k] * ab_y + g_y_z_yyyy_xxy[k];

                g_y_z_yyyyy_xy[k] = -g_0_z_yyyy_xy[k] - g_y_z_yyyy_xy[k] * ab_y + g_y_z_yyyy_xyy[k];

                g_y_z_yyyyy_xz[k] = -g_0_z_yyyy_xz[k] - g_y_z_yyyy_xz[k] * ab_y + g_y_z_yyyy_xyz[k];

                g_y_z_yyyyy_yy[k] = -g_0_z_yyyy_yy[k] - g_y_z_yyyy_yy[k] * ab_y + g_y_z_yyyy_yyy[k];

                g_y_z_yyyyy_yz[k] = -g_0_z_yyyy_yz[k] - g_y_z_yyyy_yz[k] * ab_y + g_y_z_yyyy_yyz[k];

                g_y_z_yyyyy_zz[k] = -g_0_z_yyyy_zz[k] - g_y_z_yyyy_zz[k] * ab_y + g_y_z_yyyy_yzz[k];
            }

            /// Set up 726-732 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyyz_xx = cbuffer.data(hd_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_z_yyyyz_xy = cbuffer.data(hd_geom_11_off + 727 * ccomps * dcomps);

            auto g_y_z_yyyyz_xz = cbuffer.data(hd_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_z_yyyyz_yy = cbuffer.data(hd_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_z_yyyyz_yz = cbuffer.data(hd_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_z_yyyyz_zz = cbuffer.data(hd_geom_11_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_xx, g_0_z_yyyz_xy, g_0_z_yyyz_xz, g_0_z_yyyz_yy, g_0_z_yyyz_yz, g_0_z_yyyz_zz, g_y_z_yyyyz_xx, g_y_z_yyyyz_xy, g_y_z_yyyyz_xz, g_y_z_yyyyz_yy, g_y_z_yyyyz_yz, g_y_z_yyyyz_zz, g_y_z_yyyz_xx, g_y_z_yyyz_xxy, g_y_z_yyyz_xy, g_y_z_yyyz_xyy, g_y_z_yyyz_xyz, g_y_z_yyyz_xz, g_y_z_yyyz_yy, g_y_z_yyyz_yyy, g_y_z_yyyz_yyz, g_y_z_yyyz_yz, g_y_z_yyyz_yzz, g_y_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyyz_xx[k] = -g_0_z_yyyz_xx[k] - g_y_z_yyyz_xx[k] * ab_y + g_y_z_yyyz_xxy[k];

                g_y_z_yyyyz_xy[k] = -g_0_z_yyyz_xy[k] - g_y_z_yyyz_xy[k] * ab_y + g_y_z_yyyz_xyy[k];

                g_y_z_yyyyz_xz[k] = -g_0_z_yyyz_xz[k] - g_y_z_yyyz_xz[k] * ab_y + g_y_z_yyyz_xyz[k];

                g_y_z_yyyyz_yy[k] = -g_0_z_yyyz_yy[k] - g_y_z_yyyz_yy[k] * ab_y + g_y_z_yyyz_yyy[k];

                g_y_z_yyyyz_yz[k] = -g_0_z_yyyz_yz[k] - g_y_z_yyyz_yz[k] * ab_y + g_y_z_yyyz_yyz[k];

                g_y_z_yyyyz_zz[k] = -g_0_z_yyyz_zz[k] - g_y_z_yyyz_zz[k] * ab_y + g_y_z_yyyz_yzz[k];
            }

            /// Set up 732-738 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyzz_xx = cbuffer.data(hd_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_z_yyyzz_xy = cbuffer.data(hd_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_z_yyyzz_xz = cbuffer.data(hd_geom_11_off + 734 * ccomps * dcomps);

            auto g_y_z_yyyzz_yy = cbuffer.data(hd_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_z_yyyzz_yz = cbuffer.data(hd_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_z_yyyzz_zz = cbuffer.data(hd_geom_11_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_xx, g_0_z_yyzz_xy, g_0_z_yyzz_xz, g_0_z_yyzz_yy, g_0_z_yyzz_yz, g_0_z_yyzz_zz, g_y_z_yyyzz_xx, g_y_z_yyyzz_xy, g_y_z_yyyzz_xz, g_y_z_yyyzz_yy, g_y_z_yyyzz_yz, g_y_z_yyyzz_zz, g_y_z_yyzz_xx, g_y_z_yyzz_xxy, g_y_z_yyzz_xy, g_y_z_yyzz_xyy, g_y_z_yyzz_xyz, g_y_z_yyzz_xz, g_y_z_yyzz_yy, g_y_z_yyzz_yyy, g_y_z_yyzz_yyz, g_y_z_yyzz_yz, g_y_z_yyzz_yzz, g_y_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyzz_xx[k] = -g_0_z_yyzz_xx[k] - g_y_z_yyzz_xx[k] * ab_y + g_y_z_yyzz_xxy[k];

                g_y_z_yyyzz_xy[k] = -g_0_z_yyzz_xy[k] - g_y_z_yyzz_xy[k] * ab_y + g_y_z_yyzz_xyy[k];

                g_y_z_yyyzz_xz[k] = -g_0_z_yyzz_xz[k] - g_y_z_yyzz_xz[k] * ab_y + g_y_z_yyzz_xyz[k];

                g_y_z_yyyzz_yy[k] = -g_0_z_yyzz_yy[k] - g_y_z_yyzz_yy[k] * ab_y + g_y_z_yyzz_yyy[k];

                g_y_z_yyyzz_yz[k] = -g_0_z_yyzz_yz[k] - g_y_z_yyzz_yz[k] * ab_y + g_y_z_yyzz_yyz[k];

                g_y_z_yyyzz_zz[k] = -g_0_z_yyzz_zz[k] - g_y_z_yyzz_zz[k] * ab_y + g_y_z_yyzz_yzz[k];
            }

            /// Set up 738-744 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyzzz_xx = cbuffer.data(hd_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_z_yyzzz_xy = cbuffer.data(hd_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_z_yyzzz_xz = cbuffer.data(hd_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_z_yyzzz_yy = cbuffer.data(hd_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_z_yyzzz_yz = cbuffer.data(hd_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_z_yyzzz_zz = cbuffer.data(hd_geom_11_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_xx, g_0_z_yzzz_xy, g_0_z_yzzz_xz, g_0_z_yzzz_yy, g_0_z_yzzz_yz, g_0_z_yzzz_zz, g_y_z_yyzzz_xx, g_y_z_yyzzz_xy, g_y_z_yyzzz_xz, g_y_z_yyzzz_yy, g_y_z_yyzzz_yz, g_y_z_yyzzz_zz, g_y_z_yzzz_xx, g_y_z_yzzz_xxy, g_y_z_yzzz_xy, g_y_z_yzzz_xyy, g_y_z_yzzz_xyz, g_y_z_yzzz_xz, g_y_z_yzzz_yy, g_y_z_yzzz_yyy, g_y_z_yzzz_yyz, g_y_z_yzzz_yz, g_y_z_yzzz_yzz, g_y_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyzzz_xx[k] = -g_0_z_yzzz_xx[k] - g_y_z_yzzz_xx[k] * ab_y + g_y_z_yzzz_xxy[k];

                g_y_z_yyzzz_xy[k] = -g_0_z_yzzz_xy[k] - g_y_z_yzzz_xy[k] * ab_y + g_y_z_yzzz_xyy[k];

                g_y_z_yyzzz_xz[k] = -g_0_z_yzzz_xz[k] - g_y_z_yzzz_xz[k] * ab_y + g_y_z_yzzz_xyz[k];

                g_y_z_yyzzz_yy[k] = -g_0_z_yzzz_yy[k] - g_y_z_yzzz_yy[k] * ab_y + g_y_z_yzzz_yyy[k];

                g_y_z_yyzzz_yz[k] = -g_0_z_yzzz_yz[k] - g_y_z_yzzz_yz[k] * ab_y + g_y_z_yzzz_yyz[k];

                g_y_z_yyzzz_zz[k] = -g_0_z_yzzz_zz[k] - g_y_z_yzzz_zz[k] * ab_y + g_y_z_yzzz_yzz[k];
            }

            /// Set up 744-750 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzzzz_xx = cbuffer.data(hd_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_z_yzzzz_xy = cbuffer.data(hd_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_z_yzzzz_xz = cbuffer.data(hd_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_z_yzzzz_yy = cbuffer.data(hd_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_z_yzzzz_yz = cbuffer.data(hd_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_z_yzzzz_zz = cbuffer.data(hd_geom_11_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xx, g_0_z_zzzz_xy, g_0_z_zzzz_xz, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_zz, g_y_z_yzzzz_xx, g_y_z_yzzzz_xy, g_y_z_yzzzz_xz, g_y_z_yzzzz_yy, g_y_z_yzzzz_yz, g_y_z_yzzzz_zz, g_y_z_zzzz_xx, g_y_z_zzzz_xxy, g_y_z_zzzz_xy, g_y_z_zzzz_xyy, g_y_z_zzzz_xyz, g_y_z_zzzz_xz, g_y_z_zzzz_yy, g_y_z_zzzz_yyy, g_y_z_zzzz_yyz, g_y_z_zzzz_yz, g_y_z_zzzz_yzz, g_y_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzzzz_xx[k] = -g_0_z_zzzz_xx[k] - g_y_z_zzzz_xx[k] * ab_y + g_y_z_zzzz_xxy[k];

                g_y_z_yzzzz_xy[k] = -g_0_z_zzzz_xy[k] - g_y_z_zzzz_xy[k] * ab_y + g_y_z_zzzz_xyy[k];

                g_y_z_yzzzz_xz[k] = -g_0_z_zzzz_xz[k] - g_y_z_zzzz_xz[k] * ab_y + g_y_z_zzzz_xyz[k];

                g_y_z_yzzzz_yy[k] = -g_0_z_zzzz_yy[k] - g_y_z_zzzz_yy[k] * ab_y + g_y_z_zzzz_yyy[k];

                g_y_z_yzzzz_yz[k] = -g_0_z_zzzz_yz[k] - g_y_z_zzzz_yz[k] * ab_y + g_y_z_zzzz_yyz[k];

                g_y_z_yzzzz_zz[k] = -g_0_z_zzzz_zz[k] - g_y_z_zzzz_zz[k] * ab_y + g_y_z_zzzz_yzz[k];
            }

            /// Set up 750-756 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzzzz_xx = cbuffer.data(hd_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_z_zzzzz_xy = cbuffer.data(hd_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_z_zzzzz_xz = cbuffer.data(hd_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_z_zzzzz_yy = cbuffer.data(hd_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_z_zzzzz_yz = cbuffer.data(hd_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_z_zzzzz_zz = cbuffer.data(hd_geom_11_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_xx, g_y_0_zzzz_xy, g_y_0_zzzz_xz, g_y_0_zzzz_yy, g_y_0_zzzz_yz, g_y_0_zzzz_zz, g_y_z_zzzz_xx, g_y_z_zzzz_xxz, g_y_z_zzzz_xy, g_y_z_zzzz_xyz, g_y_z_zzzz_xz, g_y_z_zzzz_xzz, g_y_z_zzzz_yy, g_y_z_zzzz_yyz, g_y_z_zzzz_yz, g_y_z_zzzz_yzz, g_y_z_zzzz_zz, g_y_z_zzzz_zzz, g_y_z_zzzzz_xx, g_y_z_zzzzz_xy, g_y_z_zzzzz_xz, g_y_z_zzzzz_yy, g_y_z_zzzzz_yz, g_y_z_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzzzz_xx[k] = g_y_0_zzzz_xx[k] - g_y_z_zzzz_xx[k] * ab_z + g_y_z_zzzz_xxz[k];

                g_y_z_zzzzz_xy[k] = g_y_0_zzzz_xy[k] - g_y_z_zzzz_xy[k] * ab_z + g_y_z_zzzz_xyz[k];

                g_y_z_zzzzz_xz[k] = g_y_0_zzzz_xz[k] - g_y_z_zzzz_xz[k] * ab_z + g_y_z_zzzz_xzz[k];

                g_y_z_zzzzz_yy[k] = g_y_0_zzzz_yy[k] - g_y_z_zzzz_yy[k] * ab_z + g_y_z_zzzz_yyz[k];

                g_y_z_zzzzz_yz[k] = g_y_0_zzzz_yz[k] - g_y_z_zzzz_yz[k] * ab_z + g_y_z_zzzz_yzz[k];

                g_y_z_zzzzz_zz[k] = g_y_0_zzzz_zz[k] - g_y_z_zzzz_zz[k] * ab_z + g_y_z_zzzz_zzz[k];
            }

            /// Set up 756-762 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxx_xx = cbuffer.data(hd_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_x_xxxxx_xy = cbuffer.data(hd_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_x_xxxxx_xz = cbuffer.data(hd_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_x_xxxxx_yy = cbuffer.data(hd_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_x_xxxxx_yz = cbuffer.data(hd_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_x_xxxxx_zz = cbuffer.data(hd_geom_11_off + 761 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxx_xx, g_z_0_xxxx_xy, g_z_0_xxxx_xz, g_z_0_xxxx_yy, g_z_0_xxxx_yz, g_z_0_xxxx_zz, g_z_x_xxxx_xx, g_z_x_xxxx_xxx, g_z_x_xxxx_xxy, g_z_x_xxxx_xxz, g_z_x_xxxx_xy, g_z_x_xxxx_xyy, g_z_x_xxxx_xyz, g_z_x_xxxx_xz, g_z_x_xxxx_xzz, g_z_x_xxxx_yy, g_z_x_xxxx_yz, g_z_x_xxxx_zz, g_z_x_xxxxx_xx, g_z_x_xxxxx_xy, g_z_x_xxxxx_xz, g_z_x_xxxxx_yy, g_z_x_xxxxx_yz, g_z_x_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxx_xx[k] = g_z_0_xxxx_xx[k] - g_z_x_xxxx_xx[k] * ab_x + g_z_x_xxxx_xxx[k];

                g_z_x_xxxxx_xy[k] = g_z_0_xxxx_xy[k] - g_z_x_xxxx_xy[k] * ab_x + g_z_x_xxxx_xxy[k];

                g_z_x_xxxxx_xz[k] = g_z_0_xxxx_xz[k] - g_z_x_xxxx_xz[k] * ab_x + g_z_x_xxxx_xxz[k];

                g_z_x_xxxxx_yy[k] = g_z_0_xxxx_yy[k] - g_z_x_xxxx_yy[k] * ab_x + g_z_x_xxxx_xyy[k];

                g_z_x_xxxxx_yz[k] = g_z_0_xxxx_yz[k] - g_z_x_xxxx_yz[k] * ab_x + g_z_x_xxxx_xyz[k];

                g_z_x_xxxxx_zz[k] = g_z_0_xxxx_zz[k] - g_z_x_xxxx_zz[k] * ab_x + g_z_x_xxxx_xzz[k];
            }

            /// Set up 762-768 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxy_xx = cbuffer.data(hd_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_x_xxxxy_xy = cbuffer.data(hd_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_x_xxxxy_xz = cbuffer.data(hd_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_x_xxxxy_yy = cbuffer.data(hd_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_x_xxxxy_yz = cbuffer.data(hd_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_x_xxxxy_zz = cbuffer.data(hd_geom_11_off + 767 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxx_xx, g_z_x_xxxx_xxy, g_z_x_xxxx_xy, g_z_x_xxxx_xyy, g_z_x_xxxx_xyz, g_z_x_xxxx_xz, g_z_x_xxxx_yy, g_z_x_xxxx_yyy, g_z_x_xxxx_yyz, g_z_x_xxxx_yz, g_z_x_xxxx_yzz, g_z_x_xxxx_zz, g_z_x_xxxxy_xx, g_z_x_xxxxy_xy, g_z_x_xxxxy_xz, g_z_x_xxxxy_yy, g_z_x_xxxxy_yz, g_z_x_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxy_xx[k] = -g_z_x_xxxx_xx[k] * ab_y + g_z_x_xxxx_xxy[k];

                g_z_x_xxxxy_xy[k] = -g_z_x_xxxx_xy[k] * ab_y + g_z_x_xxxx_xyy[k];

                g_z_x_xxxxy_xz[k] = -g_z_x_xxxx_xz[k] * ab_y + g_z_x_xxxx_xyz[k];

                g_z_x_xxxxy_yy[k] = -g_z_x_xxxx_yy[k] * ab_y + g_z_x_xxxx_yyy[k];

                g_z_x_xxxxy_yz[k] = -g_z_x_xxxx_yz[k] * ab_y + g_z_x_xxxx_yyz[k];

                g_z_x_xxxxy_zz[k] = -g_z_x_xxxx_zz[k] * ab_y + g_z_x_xxxx_yzz[k];
            }

            /// Set up 768-774 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxxz_xx = cbuffer.data(hd_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_x_xxxxz_xy = cbuffer.data(hd_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_x_xxxxz_xz = cbuffer.data(hd_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_x_xxxxz_yy = cbuffer.data(hd_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_x_xxxxz_yz = cbuffer.data(hd_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_x_xxxxz_zz = cbuffer.data(hd_geom_11_off + 773 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_xx, g_z_0_xxxz_xy, g_z_0_xxxz_xz, g_z_0_xxxz_yy, g_z_0_xxxz_yz, g_z_0_xxxz_zz, g_z_x_xxxxz_xx, g_z_x_xxxxz_xy, g_z_x_xxxxz_xz, g_z_x_xxxxz_yy, g_z_x_xxxxz_yz, g_z_x_xxxxz_zz, g_z_x_xxxz_xx, g_z_x_xxxz_xxx, g_z_x_xxxz_xxy, g_z_x_xxxz_xxz, g_z_x_xxxz_xy, g_z_x_xxxz_xyy, g_z_x_xxxz_xyz, g_z_x_xxxz_xz, g_z_x_xxxz_xzz, g_z_x_xxxz_yy, g_z_x_xxxz_yz, g_z_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxxz_xx[k] = g_z_0_xxxz_xx[k] - g_z_x_xxxz_xx[k] * ab_x + g_z_x_xxxz_xxx[k];

                g_z_x_xxxxz_xy[k] = g_z_0_xxxz_xy[k] - g_z_x_xxxz_xy[k] * ab_x + g_z_x_xxxz_xxy[k];

                g_z_x_xxxxz_xz[k] = g_z_0_xxxz_xz[k] - g_z_x_xxxz_xz[k] * ab_x + g_z_x_xxxz_xxz[k];

                g_z_x_xxxxz_yy[k] = g_z_0_xxxz_yy[k] - g_z_x_xxxz_yy[k] * ab_x + g_z_x_xxxz_xyy[k];

                g_z_x_xxxxz_yz[k] = g_z_0_xxxz_yz[k] - g_z_x_xxxz_yz[k] * ab_x + g_z_x_xxxz_xyz[k];

                g_z_x_xxxxz_zz[k] = g_z_0_xxxz_zz[k] - g_z_x_xxxz_zz[k] * ab_x + g_z_x_xxxz_xzz[k];
            }

            /// Set up 774-780 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxyy_xx = cbuffer.data(hd_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_x_xxxyy_xy = cbuffer.data(hd_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_x_xxxyy_xz = cbuffer.data(hd_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_x_xxxyy_yy = cbuffer.data(hd_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_x_xxxyy_yz = cbuffer.data(hd_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_x_xxxyy_zz = cbuffer.data(hd_geom_11_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxy_xx, g_z_x_xxxy_xxy, g_z_x_xxxy_xy, g_z_x_xxxy_xyy, g_z_x_xxxy_xyz, g_z_x_xxxy_xz, g_z_x_xxxy_yy, g_z_x_xxxy_yyy, g_z_x_xxxy_yyz, g_z_x_xxxy_yz, g_z_x_xxxy_yzz, g_z_x_xxxy_zz, g_z_x_xxxyy_xx, g_z_x_xxxyy_xy, g_z_x_xxxyy_xz, g_z_x_xxxyy_yy, g_z_x_xxxyy_yz, g_z_x_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxyy_xx[k] = -g_z_x_xxxy_xx[k] * ab_y + g_z_x_xxxy_xxy[k];

                g_z_x_xxxyy_xy[k] = -g_z_x_xxxy_xy[k] * ab_y + g_z_x_xxxy_xyy[k];

                g_z_x_xxxyy_xz[k] = -g_z_x_xxxy_xz[k] * ab_y + g_z_x_xxxy_xyz[k];

                g_z_x_xxxyy_yy[k] = -g_z_x_xxxy_yy[k] * ab_y + g_z_x_xxxy_yyy[k];

                g_z_x_xxxyy_yz[k] = -g_z_x_xxxy_yz[k] * ab_y + g_z_x_xxxy_yyz[k];

                g_z_x_xxxyy_zz[k] = -g_z_x_xxxy_zz[k] * ab_y + g_z_x_xxxy_yzz[k];
            }

            /// Set up 780-786 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxyz_xx = cbuffer.data(hd_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_x_xxxyz_xy = cbuffer.data(hd_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_x_xxxyz_xz = cbuffer.data(hd_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_x_xxxyz_yy = cbuffer.data(hd_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_x_xxxyz_yz = cbuffer.data(hd_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_x_xxxyz_zz = cbuffer.data(hd_geom_11_off + 785 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxxyz_xx, g_z_x_xxxyz_xy, g_z_x_xxxyz_xz, g_z_x_xxxyz_yy, g_z_x_xxxyz_yz, g_z_x_xxxyz_zz, g_z_x_xxxz_xx, g_z_x_xxxz_xxy, g_z_x_xxxz_xy, g_z_x_xxxz_xyy, g_z_x_xxxz_xyz, g_z_x_xxxz_xz, g_z_x_xxxz_yy, g_z_x_xxxz_yyy, g_z_x_xxxz_yyz, g_z_x_xxxz_yz, g_z_x_xxxz_yzz, g_z_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxyz_xx[k] = -g_z_x_xxxz_xx[k] * ab_y + g_z_x_xxxz_xxy[k];

                g_z_x_xxxyz_xy[k] = -g_z_x_xxxz_xy[k] * ab_y + g_z_x_xxxz_xyy[k];

                g_z_x_xxxyz_xz[k] = -g_z_x_xxxz_xz[k] * ab_y + g_z_x_xxxz_xyz[k];

                g_z_x_xxxyz_yy[k] = -g_z_x_xxxz_yy[k] * ab_y + g_z_x_xxxz_yyy[k];

                g_z_x_xxxyz_yz[k] = -g_z_x_xxxz_yz[k] * ab_y + g_z_x_xxxz_yyz[k];

                g_z_x_xxxyz_zz[k] = -g_z_x_xxxz_zz[k] * ab_y + g_z_x_xxxz_yzz[k];
            }

            /// Set up 786-792 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxzz_xx = cbuffer.data(hd_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_x_xxxzz_xy = cbuffer.data(hd_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_x_xxxzz_xz = cbuffer.data(hd_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_x_xxxzz_yy = cbuffer.data(hd_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_x_xxxzz_yz = cbuffer.data(hd_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_x_xxxzz_zz = cbuffer.data(hd_geom_11_off + 791 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_xx, g_z_0_xxzz_xy, g_z_0_xxzz_xz, g_z_0_xxzz_yy, g_z_0_xxzz_yz, g_z_0_xxzz_zz, g_z_x_xxxzz_xx, g_z_x_xxxzz_xy, g_z_x_xxxzz_xz, g_z_x_xxxzz_yy, g_z_x_xxxzz_yz, g_z_x_xxxzz_zz, g_z_x_xxzz_xx, g_z_x_xxzz_xxx, g_z_x_xxzz_xxy, g_z_x_xxzz_xxz, g_z_x_xxzz_xy, g_z_x_xxzz_xyy, g_z_x_xxzz_xyz, g_z_x_xxzz_xz, g_z_x_xxzz_xzz, g_z_x_xxzz_yy, g_z_x_xxzz_yz, g_z_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxzz_xx[k] = g_z_0_xxzz_xx[k] - g_z_x_xxzz_xx[k] * ab_x + g_z_x_xxzz_xxx[k];

                g_z_x_xxxzz_xy[k] = g_z_0_xxzz_xy[k] - g_z_x_xxzz_xy[k] * ab_x + g_z_x_xxzz_xxy[k];

                g_z_x_xxxzz_xz[k] = g_z_0_xxzz_xz[k] - g_z_x_xxzz_xz[k] * ab_x + g_z_x_xxzz_xxz[k];

                g_z_x_xxxzz_yy[k] = g_z_0_xxzz_yy[k] - g_z_x_xxzz_yy[k] * ab_x + g_z_x_xxzz_xyy[k];

                g_z_x_xxxzz_yz[k] = g_z_0_xxzz_yz[k] - g_z_x_xxzz_yz[k] * ab_x + g_z_x_xxzz_xyz[k];

                g_z_x_xxxzz_zz[k] = g_z_0_xxzz_zz[k] - g_z_x_xxzz_zz[k] * ab_x + g_z_x_xxzz_xzz[k];
            }

            /// Set up 792-798 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyyy_xx = cbuffer.data(hd_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_x_xxyyy_xy = cbuffer.data(hd_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_x_xxyyy_xz = cbuffer.data(hd_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_x_xxyyy_yy = cbuffer.data(hd_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_x_xxyyy_yz = cbuffer.data(hd_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_x_xxyyy_zz = cbuffer.data(hd_geom_11_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyy_xx, g_z_x_xxyy_xxy, g_z_x_xxyy_xy, g_z_x_xxyy_xyy, g_z_x_xxyy_xyz, g_z_x_xxyy_xz, g_z_x_xxyy_yy, g_z_x_xxyy_yyy, g_z_x_xxyy_yyz, g_z_x_xxyy_yz, g_z_x_xxyy_yzz, g_z_x_xxyy_zz, g_z_x_xxyyy_xx, g_z_x_xxyyy_xy, g_z_x_xxyyy_xz, g_z_x_xxyyy_yy, g_z_x_xxyyy_yz, g_z_x_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyyy_xx[k] = -g_z_x_xxyy_xx[k] * ab_y + g_z_x_xxyy_xxy[k];

                g_z_x_xxyyy_xy[k] = -g_z_x_xxyy_xy[k] * ab_y + g_z_x_xxyy_xyy[k];

                g_z_x_xxyyy_xz[k] = -g_z_x_xxyy_xz[k] * ab_y + g_z_x_xxyy_xyz[k];

                g_z_x_xxyyy_yy[k] = -g_z_x_xxyy_yy[k] * ab_y + g_z_x_xxyy_yyy[k];

                g_z_x_xxyyy_yz[k] = -g_z_x_xxyy_yz[k] * ab_y + g_z_x_xxyy_yyz[k];

                g_z_x_xxyyy_zz[k] = -g_z_x_xxyy_zz[k] * ab_y + g_z_x_xxyy_yzz[k];
            }

            /// Set up 798-804 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyyz_xx = cbuffer.data(hd_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_x_xxyyz_xy = cbuffer.data(hd_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_x_xxyyz_xz = cbuffer.data(hd_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_x_xxyyz_yy = cbuffer.data(hd_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_x_xxyyz_yz = cbuffer.data(hd_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_x_xxyyz_zz = cbuffer.data(hd_geom_11_off + 803 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyyz_xx, g_z_x_xxyyz_xy, g_z_x_xxyyz_xz, g_z_x_xxyyz_yy, g_z_x_xxyyz_yz, g_z_x_xxyyz_zz, g_z_x_xxyz_xx, g_z_x_xxyz_xxy, g_z_x_xxyz_xy, g_z_x_xxyz_xyy, g_z_x_xxyz_xyz, g_z_x_xxyz_xz, g_z_x_xxyz_yy, g_z_x_xxyz_yyy, g_z_x_xxyz_yyz, g_z_x_xxyz_yz, g_z_x_xxyz_yzz, g_z_x_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyyz_xx[k] = -g_z_x_xxyz_xx[k] * ab_y + g_z_x_xxyz_xxy[k];

                g_z_x_xxyyz_xy[k] = -g_z_x_xxyz_xy[k] * ab_y + g_z_x_xxyz_xyy[k];

                g_z_x_xxyyz_xz[k] = -g_z_x_xxyz_xz[k] * ab_y + g_z_x_xxyz_xyz[k];

                g_z_x_xxyyz_yy[k] = -g_z_x_xxyz_yy[k] * ab_y + g_z_x_xxyz_yyy[k];

                g_z_x_xxyyz_yz[k] = -g_z_x_xxyz_yz[k] * ab_y + g_z_x_xxyz_yyz[k];

                g_z_x_xxyyz_zz[k] = -g_z_x_xxyz_zz[k] * ab_y + g_z_x_xxyz_yzz[k];
            }

            /// Set up 804-810 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyzz_xx = cbuffer.data(hd_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_x_xxyzz_xy = cbuffer.data(hd_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_x_xxyzz_xz = cbuffer.data(hd_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_x_xxyzz_yy = cbuffer.data(hd_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_x_xxyzz_yz = cbuffer.data(hd_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_x_xxyzz_zz = cbuffer.data(hd_geom_11_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyzz_xx, g_z_x_xxyzz_xy, g_z_x_xxyzz_xz, g_z_x_xxyzz_yy, g_z_x_xxyzz_yz, g_z_x_xxyzz_zz, g_z_x_xxzz_xx, g_z_x_xxzz_xxy, g_z_x_xxzz_xy, g_z_x_xxzz_xyy, g_z_x_xxzz_xyz, g_z_x_xxzz_xz, g_z_x_xxzz_yy, g_z_x_xxzz_yyy, g_z_x_xxzz_yyz, g_z_x_xxzz_yz, g_z_x_xxzz_yzz, g_z_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyzz_xx[k] = -g_z_x_xxzz_xx[k] * ab_y + g_z_x_xxzz_xxy[k];

                g_z_x_xxyzz_xy[k] = -g_z_x_xxzz_xy[k] * ab_y + g_z_x_xxzz_xyy[k];

                g_z_x_xxyzz_xz[k] = -g_z_x_xxzz_xz[k] * ab_y + g_z_x_xxzz_xyz[k];

                g_z_x_xxyzz_yy[k] = -g_z_x_xxzz_yy[k] * ab_y + g_z_x_xxzz_yyy[k];

                g_z_x_xxyzz_yz[k] = -g_z_x_xxzz_yz[k] * ab_y + g_z_x_xxzz_yyz[k];

                g_z_x_xxyzz_zz[k] = -g_z_x_xxzz_zz[k] * ab_y + g_z_x_xxzz_yzz[k];
            }

            /// Set up 810-816 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxzzz_xx = cbuffer.data(hd_geom_11_off + 810 * ccomps * dcomps);

            auto g_z_x_xxzzz_xy = cbuffer.data(hd_geom_11_off + 811 * ccomps * dcomps);

            auto g_z_x_xxzzz_xz = cbuffer.data(hd_geom_11_off + 812 * ccomps * dcomps);

            auto g_z_x_xxzzz_yy = cbuffer.data(hd_geom_11_off + 813 * ccomps * dcomps);

            auto g_z_x_xxzzz_yz = cbuffer.data(hd_geom_11_off + 814 * ccomps * dcomps);

            auto g_z_x_xxzzz_zz = cbuffer.data(hd_geom_11_off + 815 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_xx, g_z_0_xzzz_xy, g_z_0_xzzz_xz, g_z_0_xzzz_yy, g_z_0_xzzz_yz, g_z_0_xzzz_zz, g_z_x_xxzzz_xx, g_z_x_xxzzz_xy, g_z_x_xxzzz_xz, g_z_x_xxzzz_yy, g_z_x_xxzzz_yz, g_z_x_xxzzz_zz, g_z_x_xzzz_xx, g_z_x_xzzz_xxx, g_z_x_xzzz_xxy, g_z_x_xzzz_xxz, g_z_x_xzzz_xy, g_z_x_xzzz_xyy, g_z_x_xzzz_xyz, g_z_x_xzzz_xz, g_z_x_xzzz_xzz, g_z_x_xzzz_yy, g_z_x_xzzz_yz, g_z_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxzzz_xx[k] = g_z_0_xzzz_xx[k] - g_z_x_xzzz_xx[k] * ab_x + g_z_x_xzzz_xxx[k];

                g_z_x_xxzzz_xy[k] = g_z_0_xzzz_xy[k] - g_z_x_xzzz_xy[k] * ab_x + g_z_x_xzzz_xxy[k];

                g_z_x_xxzzz_xz[k] = g_z_0_xzzz_xz[k] - g_z_x_xzzz_xz[k] * ab_x + g_z_x_xzzz_xxz[k];

                g_z_x_xxzzz_yy[k] = g_z_0_xzzz_yy[k] - g_z_x_xzzz_yy[k] * ab_x + g_z_x_xzzz_xyy[k];

                g_z_x_xxzzz_yz[k] = g_z_0_xzzz_yz[k] - g_z_x_xzzz_yz[k] * ab_x + g_z_x_xzzz_xyz[k];

                g_z_x_xxzzz_zz[k] = g_z_0_xzzz_zz[k] - g_z_x_xzzz_zz[k] * ab_x + g_z_x_xzzz_xzz[k];
            }

            /// Set up 816-822 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyyy_xx = cbuffer.data(hd_geom_11_off + 816 * ccomps * dcomps);

            auto g_z_x_xyyyy_xy = cbuffer.data(hd_geom_11_off + 817 * ccomps * dcomps);

            auto g_z_x_xyyyy_xz = cbuffer.data(hd_geom_11_off + 818 * ccomps * dcomps);

            auto g_z_x_xyyyy_yy = cbuffer.data(hd_geom_11_off + 819 * ccomps * dcomps);

            auto g_z_x_xyyyy_yz = cbuffer.data(hd_geom_11_off + 820 * ccomps * dcomps);

            auto g_z_x_xyyyy_zz = cbuffer.data(hd_geom_11_off + 821 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyy_xx, g_z_x_xyyy_xxy, g_z_x_xyyy_xy, g_z_x_xyyy_xyy, g_z_x_xyyy_xyz, g_z_x_xyyy_xz, g_z_x_xyyy_yy, g_z_x_xyyy_yyy, g_z_x_xyyy_yyz, g_z_x_xyyy_yz, g_z_x_xyyy_yzz, g_z_x_xyyy_zz, g_z_x_xyyyy_xx, g_z_x_xyyyy_xy, g_z_x_xyyyy_xz, g_z_x_xyyyy_yy, g_z_x_xyyyy_yz, g_z_x_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyyy_xx[k] = -g_z_x_xyyy_xx[k] * ab_y + g_z_x_xyyy_xxy[k];

                g_z_x_xyyyy_xy[k] = -g_z_x_xyyy_xy[k] * ab_y + g_z_x_xyyy_xyy[k];

                g_z_x_xyyyy_xz[k] = -g_z_x_xyyy_xz[k] * ab_y + g_z_x_xyyy_xyz[k];

                g_z_x_xyyyy_yy[k] = -g_z_x_xyyy_yy[k] * ab_y + g_z_x_xyyy_yyy[k];

                g_z_x_xyyyy_yz[k] = -g_z_x_xyyy_yz[k] * ab_y + g_z_x_xyyy_yyz[k];

                g_z_x_xyyyy_zz[k] = -g_z_x_xyyy_zz[k] * ab_y + g_z_x_xyyy_yzz[k];
            }

            /// Set up 822-828 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyyz_xx = cbuffer.data(hd_geom_11_off + 822 * ccomps * dcomps);

            auto g_z_x_xyyyz_xy = cbuffer.data(hd_geom_11_off + 823 * ccomps * dcomps);

            auto g_z_x_xyyyz_xz = cbuffer.data(hd_geom_11_off + 824 * ccomps * dcomps);

            auto g_z_x_xyyyz_yy = cbuffer.data(hd_geom_11_off + 825 * ccomps * dcomps);

            auto g_z_x_xyyyz_yz = cbuffer.data(hd_geom_11_off + 826 * ccomps * dcomps);

            auto g_z_x_xyyyz_zz = cbuffer.data(hd_geom_11_off + 827 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyyz_xx, g_z_x_xyyyz_xy, g_z_x_xyyyz_xz, g_z_x_xyyyz_yy, g_z_x_xyyyz_yz, g_z_x_xyyyz_zz, g_z_x_xyyz_xx, g_z_x_xyyz_xxy, g_z_x_xyyz_xy, g_z_x_xyyz_xyy, g_z_x_xyyz_xyz, g_z_x_xyyz_xz, g_z_x_xyyz_yy, g_z_x_xyyz_yyy, g_z_x_xyyz_yyz, g_z_x_xyyz_yz, g_z_x_xyyz_yzz, g_z_x_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyyz_xx[k] = -g_z_x_xyyz_xx[k] * ab_y + g_z_x_xyyz_xxy[k];

                g_z_x_xyyyz_xy[k] = -g_z_x_xyyz_xy[k] * ab_y + g_z_x_xyyz_xyy[k];

                g_z_x_xyyyz_xz[k] = -g_z_x_xyyz_xz[k] * ab_y + g_z_x_xyyz_xyz[k];

                g_z_x_xyyyz_yy[k] = -g_z_x_xyyz_yy[k] * ab_y + g_z_x_xyyz_yyy[k];

                g_z_x_xyyyz_yz[k] = -g_z_x_xyyz_yz[k] * ab_y + g_z_x_xyyz_yyz[k];

                g_z_x_xyyyz_zz[k] = -g_z_x_xyyz_zz[k] * ab_y + g_z_x_xyyz_yzz[k];
            }

            /// Set up 828-834 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyzz_xx = cbuffer.data(hd_geom_11_off + 828 * ccomps * dcomps);

            auto g_z_x_xyyzz_xy = cbuffer.data(hd_geom_11_off + 829 * ccomps * dcomps);

            auto g_z_x_xyyzz_xz = cbuffer.data(hd_geom_11_off + 830 * ccomps * dcomps);

            auto g_z_x_xyyzz_yy = cbuffer.data(hd_geom_11_off + 831 * ccomps * dcomps);

            auto g_z_x_xyyzz_yz = cbuffer.data(hd_geom_11_off + 832 * ccomps * dcomps);

            auto g_z_x_xyyzz_zz = cbuffer.data(hd_geom_11_off + 833 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyzz_xx, g_z_x_xyyzz_xy, g_z_x_xyyzz_xz, g_z_x_xyyzz_yy, g_z_x_xyyzz_yz, g_z_x_xyyzz_zz, g_z_x_xyzz_xx, g_z_x_xyzz_xxy, g_z_x_xyzz_xy, g_z_x_xyzz_xyy, g_z_x_xyzz_xyz, g_z_x_xyzz_xz, g_z_x_xyzz_yy, g_z_x_xyzz_yyy, g_z_x_xyzz_yyz, g_z_x_xyzz_yz, g_z_x_xyzz_yzz, g_z_x_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyzz_xx[k] = -g_z_x_xyzz_xx[k] * ab_y + g_z_x_xyzz_xxy[k];

                g_z_x_xyyzz_xy[k] = -g_z_x_xyzz_xy[k] * ab_y + g_z_x_xyzz_xyy[k];

                g_z_x_xyyzz_xz[k] = -g_z_x_xyzz_xz[k] * ab_y + g_z_x_xyzz_xyz[k];

                g_z_x_xyyzz_yy[k] = -g_z_x_xyzz_yy[k] * ab_y + g_z_x_xyzz_yyy[k];

                g_z_x_xyyzz_yz[k] = -g_z_x_xyzz_yz[k] * ab_y + g_z_x_xyzz_yyz[k];

                g_z_x_xyyzz_zz[k] = -g_z_x_xyzz_zz[k] * ab_y + g_z_x_xyzz_yzz[k];
            }

            /// Set up 834-840 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyzzz_xx = cbuffer.data(hd_geom_11_off + 834 * ccomps * dcomps);

            auto g_z_x_xyzzz_xy = cbuffer.data(hd_geom_11_off + 835 * ccomps * dcomps);

            auto g_z_x_xyzzz_xz = cbuffer.data(hd_geom_11_off + 836 * ccomps * dcomps);

            auto g_z_x_xyzzz_yy = cbuffer.data(hd_geom_11_off + 837 * ccomps * dcomps);

            auto g_z_x_xyzzz_yz = cbuffer.data(hd_geom_11_off + 838 * ccomps * dcomps);

            auto g_z_x_xyzzz_zz = cbuffer.data(hd_geom_11_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyzzz_xx, g_z_x_xyzzz_xy, g_z_x_xyzzz_xz, g_z_x_xyzzz_yy, g_z_x_xyzzz_yz, g_z_x_xyzzz_zz, g_z_x_xzzz_xx, g_z_x_xzzz_xxy, g_z_x_xzzz_xy, g_z_x_xzzz_xyy, g_z_x_xzzz_xyz, g_z_x_xzzz_xz, g_z_x_xzzz_yy, g_z_x_xzzz_yyy, g_z_x_xzzz_yyz, g_z_x_xzzz_yz, g_z_x_xzzz_yzz, g_z_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyzzz_xx[k] = -g_z_x_xzzz_xx[k] * ab_y + g_z_x_xzzz_xxy[k];

                g_z_x_xyzzz_xy[k] = -g_z_x_xzzz_xy[k] * ab_y + g_z_x_xzzz_xyy[k];

                g_z_x_xyzzz_xz[k] = -g_z_x_xzzz_xz[k] * ab_y + g_z_x_xzzz_xyz[k];

                g_z_x_xyzzz_yy[k] = -g_z_x_xzzz_yy[k] * ab_y + g_z_x_xzzz_yyy[k];

                g_z_x_xyzzz_yz[k] = -g_z_x_xzzz_yz[k] * ab_y + g_z_x_xzzz_yyz[k];

                g_z_x_xyzzz_zz[k] = -g_z_x_xzzz_zz[k] * ab_y + g_z_x_xzzz_yzz[k];
            }

            /// Set up 840-846 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzzzz_xx = cbuffer.data(hd_geom_11_off + 840 * ccomps * dcomps);

            auto g_z_x_xzzzz_xy = cbuffer.data(hd_geom_11_off + 841 * ccomps * dcomps);

            auto g_z_x_xzzzz_xz = cbuffer.data(hd_geom_11_off + 842 * ccomps * dcomps);

            auto g_z_x_xzzzz_yy = cbuffer.data(hd_geom_11_off + 843 * ccomps * dcomps);

            auto g_z_x_xzzzz_yz = cbuffer.data(hd_geom_11_off + 844 * ccomps * dcomps);

            auto g_z_x_xzzzz_zz = cbuffer.data(hd_geom_11_off + 845 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz, g_z_x_xzzzz_xx, g_z_x_xzzzz_xy, g_z_x_xzzzz_xz, g_z_x_xzzzz_yy, g_z_x_xzzzz_yz, g_z_x_xzzzz_zz, g_z_x_zzzz_xx, g_z_x_zzzz_xxx, g_z_x_zzzz_xxy, g_z_x_zzzz_xxz, g_z_x_zzzz_xy, g_z_x_zzzz_xyy, g_z_x_zzzz_xyz, g_z_x_zzzz_xz, g_z_x_zzzz_xzz, g_z_x_zzzz_yy, g_z_x_zzzz_yz, g_z_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzzzz_xx[k] = g_z_0_zzzz_xx[k] - g_z_x_zzzz_xx[k] * ab_x + g_z_x_zzzz_xxx[k];

                g_z_x_xzzzz_xy[k] = g_z_0_zzzz_xy[k] - g_z_x_zzzz_xy[k] * ab_x + g_z_x_zzzz_xxy[k];

                g_z_x_xzzzz_xz[k] = g_z_0_zzzz_xz[k] - g_z_x_zzzz_xz[k] * ab_x + g_z_x_zzzz_xxz[k];

                g_z_x_xzzzz_yy[k] = g_z_0_zzzz_yy[k] - g_z_x_zzzz_yy[k] * ab_x + g_z_x_zzzz_xyy[k];

                g_z_x_xzzzz_yz[k] = g_z_0_zzzz_yz[k] - g_z_x_zzzz_yz[k] * ab_x + g_z_x_zzzz_xyz[k];

                g_z_x_xzzzz_zz[k] = g_z_0_zzzz_zz[k] - g_z_x_zzzz_zz[k] * ab_x + g_z_x_zzzz_xzz[k];
            }

            /// Set up 846-852 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyyy_xx = cbuffer.data(hd_geom_11_off + 846 * ccomps * dcomps);

            auto g_z_x_yyyyy_xy = cbuffer.data(hd_geom_11_off + 847 * ccomps * dcomps);

            auto g_z_x_yyyyy_xz = cbuffer.data(hd_geom_11_off + 848 * ccomps * dcomps);

            auto g_z_x_yyyyy_yy = cbuffer.data(hd_geom_11_off + 849 * ccomps * dcomps);

            auto g_z_x_yyyyy_yz = cbuffer.data(hd_geom_11_off + 850 * ccomps * dcomps);

            auto g_z_x_yyyyy_zz = cbuffer.data(hd_geom_11_off + 851 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyy_xx, g_z_x_yyyy_xxy, g_z_x_yyyy_xy, g_z_x_yyyy_xyy, g_z_x_yyyy_xyz, g_z_x_yyyy_xz, g_z_x_yyyy_yy, g_z_x_yyyy_yyy, g_z_x_yyyy_yyz, g_z_x_yyyy_yz, g_z_x_yyyy_yzz, g_z_x_yyyy_zz, g_z_x_yyyyy_xx, g_z_x_yyyyy_xy, g_z_x_yyyyy_xz, g_z_x_yyyyy_yy, g_z_x_yyyyy_yz, g_z_x_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyyy_xx[k] = -g_z_x_yyyy_xx[k] * ab_y + g_z_x_yyyy_xxy[k];

                g_z_x_yyyyy_xy[k] = -g_z_x_yyyy_xy[k] * ab_y + g_z_x_yyyy_xyy[k];

                g_z_x_yyyyy_xz[k] = -g_z_x_yyyy_xz[k] * ab_y + g_z_x_yyyy_xyz[k];

                g_z_x_yyyyy_yy[k] = -g_z_x_yyyy_yy[k] * ab_y + g_z_x_yyyy_yyy[k];

                g_z_x_yyyyy_yz[k] = -g_z_x_yyyy_yz[k] * ab_y + g_z_x_yyyy_yyz[k];

                g_z_x_yyyyy_zz[k] = -g_z_x_yyyy_zz[k] * ab_y + g_z_x_yyyy_yzz[k];
            }

            /// Set up 852-858 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyyz_xx = cbuffer.data(hd_geom_11_off + 852 * ccomps * dcomps);

            auto g_z_x_yyyyz_xy = cbuffer.data(hd_geom_11_off + 853 * ccomps * dcomps);

            auto g_z_x_yyyyz_xz = cbuffer.data(hd_geom_11_off + 854 * ccomps * dcomps);

            auto g_z_x_yyyyz_yy = cbuffer.data(hd_geom_11_off + 855 * ccomps * dcomps);

            auto g_z_x_yyyyz_yz = cbuffer.data(hd_geom_11_off + 856 * ccomps * dcomps);

            auto g_z_x_yyyyz_zz = cbuffer.data(hd_geom_11_off + 857 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyyz_xx, g_z_x_yyyyz_xy, g_z_x_yyyyz_xz, g_z_x_yyyyz_yy, g_z_x_yyyyz_yz, g_z_x_yyyyz_zz, g_z_x_yyyz_xx, g_z_x_yyyz_xxy, g_z_x_yyyz_xy, g_z_x_yyyz_xyy, g_z_x_yyyz_xyz, g_z_x_yyyz_xz, g_z_x_yyyz_yy, g_z_x_yyyz_yyy, g_z_x_yyyz_yyz, g_z_x_yyyz_yz, g_z_x_yyyz_yzz, g_z_x_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyyz_xx[k] = -g_z_x_yyyz_xx[k] * ab_y + g_z_x_yyyz_xxy[k];

                g_z_x_yyyyz_xy[k] = -g_z_x_yyyz_xy[k] * ab_y + g_z_x_yyyz_xyy[k];

                g_z_x_yyyyz_xz[k] = -g_z_x_yyyz_xz[k] * ab_y + g_z_x_yyyz_xyz[k];

                g_z_x_yyyyz_yy[k] = -g_z_x_yyyz_yy[k] * ab_y + g_z_x_yyyz_yyy[k];

                g_z_x_yyyyz_yz[k] = -g_z_x_yyyz_yz[k] * ab_y + g_z_x_yyyz_yyz[k];

                g_z_x_yyyyz_zz[k] = -g_z_x_yyyz_zz[k] * ab_y + g_z_x_yyyz_yzz[k];
            }

            /// Set up 858-864 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyzz_xx = cbuffer.data(hd_geom_11_off + 858 * ccomps * dcomps);

            auto g_z_x_yyyzz_xy = cbuffer.data(hd_geom_11_off + 859 * ccomps * dcomps);

            auto g_z_x_yyyzz_xz = cbuffer.data(hd_geom_11_off + 860 * ccomps * dcomps);

            auto g_z_x_yyyzz_yy = cbuffer.data(hd_geom_11_off + 861 * ccomps * dcomps);

            auto g_z_x_yyyzz_yz = cbuffer.data(hd_geom_11_off + 862 * ccomps * dcomps);

            auto g_z_x_yyyzz_zz = cbuffer.data(hd_geom_11_off + 863 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyzz_xx, g_z_x_yyyzz_xy, g_z_x_yyyzz_xz, g_z_x_yyyzz_yy, g_z_x_yyyzz_yz, g_z_x_yyyzz_zz, g_z_x_yyzz_xx, g_z_x_yyzz_xxy, g_z_x_yyzz_xy, g_z_x_yyzz_xyy, g_z_x_yyzz_xyz, g_z_x_yyzz_xz, g_z_x_yyzz_yy, g_z_x_yyzz_yyy, g_z_x_yyzz_yyz, g_z_x_yyzz_yz, g_z_x_yyzz_yzz, g_z_x_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyzz_xx[k] = -g_z_x_yyzz_xx[k] * ab_y + g_z_x_yyzz_xxy[k];

                g_z_x_yyyzz_xy[k] = -g_z_x_yyzz_xy[k] * ab_y + g_z_x_yyzz_xyy[k];

                g_z_x_yyyzz_xz[k] = -g_z_x_yyzz_xz[k] * ab_y + g_z_x_yyzz_xyz[k];

                g_z_x_yyyzz_yy[k] = -g_z_x_yyzz_yy[k] * ab_y + g_z_x_yyzz_yyy[k];

                g_z_x_yyyzz_yz[k] = -g_z_x_yyzz_yz[k] * ab_y + g_z_x_yyzz_yyz[k];

                g_z_x_yyyzz_zz[k] = -g_z_x_yyzz_zz[k] * ab_y + g_z_x_yyzz_yzz[k];
            }

            /// Set up 864-870 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyzzz_xx = cbuffer.data(hd_geom_11_off + 864 * ccomps * dcomps);

            auto g_z_x_yyzzz_xy = cbuffer.data(hd_geom_11_off + 865 * ccomps * dcomps);

            auto g_z_x_yyzzz_xz = cbuffer.data(hd_geom_11_off + 866 * ccomps * dcomps);

            auto g_z_x_yyzzz_yy = cbuffer.data(hd_geom_11_off + 867 * ccomps * dcomps);

            auto g_z_x_yyzzz_yz = cbuffer.data(hd_geom_11_off + 868 * ccomps * dcomps);

            auto g_z_x_yyzzz_zz = cbuffer.data(hd_geom_11_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyzzz_xx, g_z_x_yyzzz_xy, g_z_x_yyzzz_xz, g_z_x_yyzzz_yy, g_z_x_yyzzz_yz, g_z_x_yyzzz_zz, g_z_x_yzzz_xx, g_z_x_yzzz_xxy, g_z_x_yzzz_xy, g_z_x_yzzz_xyy, g_z_x_yzzz_xyz, g_z_x_yzzz_xz, g_z_x_yzzz_yy, g_z_x_yzzz_yyy, g_z_x_yzzz_yyz, g_z_x_yzzz_yz, g_z_x_yzzz_yzz, g_z_x_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyzzz_xx[k] = -g_z_x_yzzz_xx[k] * ab_y + g_z_x_yzzz_xxy[k];

                g_z_x_yyzzz_xy[k] = -g_z_x_yzzz_xy[k] * ab_y + g_z_x_yzzz_xyy[k];

                g_z_x_yyzzz_xz[k] = -g_z_x_yzzz_xz[k] * ab_y + g_z_x_yzzz_xyz[k];

                g_z_x_yyzzz_yy[k] = -g_z_x_yzzz_yy[k] * ab_y + g_z_x_yzzz_yyy[k];

                g_z_x_yyzzz_yz[k] = -g_z_x_yzzz_yz[k] * ab_y + g_z_x_yzzz_yyz[k];

                g_z_x_yyzzz_zz[k] = -g_z_x_yzzz_zz[k] * ab_y + g_z_x_yzzz_yzz[k];
            }

            /// Set up 870-876 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzzzz_xx = cbuffer.data(hd_geom_11_off + 870 * ccomps * dcomps);

            auto g_z_x_yzzzz_xy = cbuffer.data(hd_geom_11_off + 871 * ccomps * dcomps);

            auto g_z_x_yzzzz_xz = cbuffer.data(hd_geom_11_off + 872 * ccomps * dcomps);

            auto g_z_x_yzzzz_yy = cbuffer.data(hd_geom_11_off + 873 * ccomps * dcomps);

            auto g_z_x_yzzzz_yz = cbuffer.data(hd_geom_11_off + 874 * ccomps * dcomps);

            auto g_z_x_yzzzz_zz = cbuffer.data(hd_geom_11_off + 875 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzzzz_xx, g_z_x_yzzzz_xy, g_z_x_yzzzz_xz, g_z_x_yzzzz_yy, g_z_x_yzzzz_yz, g_z_x_yzzzz_zz, g_z_x_zzzz_xx, g_z_x_zzzz_xxy, g_z_x_zzzz_xy, g_z_x_zzzz_xyy, g_z_x_zzzz_xyz, g_z_x_zzzz_xz, g_z_x_zzzz_yy, g_z_x_zzzz_yyy, g_z_x_zzzz_yyz, g_z_x_zzzz_yz, g_z_x_zzzz_yzz, g_z_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzzzz_xx[k] = -g_z_x_zzzz_xx[k] * ab_y + g_z_x_zzzz_xxy[k];

                g_z_x_yzzzz_xy[k] = -g_z_x_zzzz_xy[k] * ab_y + g_z_x_zzzz_xyy[k];

                g_z_x_yzzzz_xz[k] = -g_z_x_zzzz_xz[k] * ab_y + g_z_x_zzzz_xyz[k];

                g_z_x_yzzzz_yy[k] = -g_z_x_zzzz_yy[k] * ab_y + g_z_x_zzzz_yyy[k];

                g_z_x_yzzzz_yz[k] = -g_z_x_zzzz_yz[k] * ab_y + g_z_x_zzzz_yyz[k];

                g_z_x_yzzzz_zz[k] = -g_z_x_zzzz_zz[k] * ab_y + g_z_x_zzzz_yzz[k];
            }

            /// Set up 876-882 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzzzz_xx = cbuffer.data(hd_geom_11_off + 876 * ccomps * dcomps);

            auto g_z_x_zzzzz_xy = cbuffer.data(hd_geom_11_off + 877 * ccomps * dcomps);

            auto g_z_x_zzzzz_xz = cbuffer.data(hd_geom_11_off + 878 * ccomps * dcomps);

            auto g_z_x_zzzzz_yy = cbuffer.data(hd_geom_11_off + 879 * ccomps * dcomps);

            auto g_z_x_zzzzz_yz = cbuffer.data(hd_geom_11_off + 880 * ccomps * dcomps);

            auto g_z_x_zzzzz_zz = cbuffer.data(hd_geom_11_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_xx, g_0_x_zzzz_xy, g_0_x_zzzz_xz, g_0_x_zzzz_yy, g_0_x_zzzz_yz, g_0_x_zzzz_zz, g_z_x_zzzz_xx, g_z_x_zzzz_xxz, g_z_x_zzzz_xy, g_z_x_zzzz_xyz, g_z_x_zzzz_xz, g_z_x_zzzz_xzz, g_z_x_zzzz_yy, g_z_x_zzzz_yyz, g_z_x_zzzz_yz, g_z_x_zzzz_yzz, g_z_x_zzzz_zz, g_z_x_zzzz_zzz, g_z_x_zzzzz_xx, g_z_x_zzzzz_xy, g_z_x_zzzzz_xz, g_z_x_zzzzz_yy, g_z_x_zzzzz_yz, g_z_x_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzzzz_xx[k] = -g_0_x_zzzz_xx[k] - g_z_x_zzzz_xx[k] * ab_z + g_z_x_zzzz_xxz[k];

                g_z_x_zzzzz_xy[k] = -g_0_x_zzzz_xy[k] - g_z_x_zzzz_xy[k] * ab_z + g_z_x_zzzz_xyz[k];

                g_z_x_zzzzz_xz[k] = -g_0_x_zzzz_xz[k] - g_z_x_zzzz_xz[k] * ab_z + g_z_x_zzzz_xzz[k];

                g_z_x_zzzzz_yy[k] = -g_0_x_zzzz_yy[k] - g_z_x_zzzz_yy[k] * ab_z + g_z_x_zzzz_yyz[k];

                g_z_x_zzzzz_yz[k] = -g_0_x_zzzz_yz[k] - g_z_x_zzzz_yz[k] * ab_z + g_z_x_zzzz_yzz[k];

                g_z_x_zzzzz_zz[k] = -g_0_x_zzzz_zz[k] - g_z_x_zzzz_zz[k] * ab_z + g_z_x_zzzz_zzz[k];
            }

            /// Set up 882-888 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxx_xx = cbuffer.data(hd_geom_11_off + 882 * ccomps * dcomps);

            auto g_z_y_xxxxx_xy = cbuffer.data(hd_geom_11_off + 883 * ccomps * dcomps);

            auto g_z_y_xxxxx_xz = cbuffer.data(hd_geom_11_off + 884 * ccomps * dcomps);

            auto g_z_y_xxxxx_yy = cbuffer.data(hd_geom_11_off + 885 * ccomps * dcomps);

            auto g_z_y_xxxxx_yz = cbuffer.data(hd_geom_11_off + 886 * ccomps * dcomps);

            auto g_z_y_xxxxx_zz = cbuffer.data(hd_geom_11_off + 887 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxx_xx, g_z_y_xxxx_xxx, g_z_y_xxxx_xxy, g_z_y_xxxx_xxz, g_z_y_xxxx_xy, g_z_y_xxxx_xyy, g_z_y_xxxx_xyz, g_z_y_xxxx_xz, g_z_y_xxxx_xzz, g_z_y_xxxx_yy, g_z_y_xxxx_yz, g_z_y_xxxx_zz, g_z_y_xxxxx_xx, g_z_y_xxxxx_xy, g_z_y_xxxxx_xz, g_z_y_xxxxx_yy, g_z_y_xxxxx_yz, g_z_y_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxx_xx[k] = -g_z_y_xxxx_xx[k] * ab_x + g_z_y_xxxx_xxx[k];

                g_z_y_xxxxx_xy[k] = -g_z_y_xxxx_xy[k] * ab_x + g_z_y_xxxx_xxy[k];

                g_z_y_xxxxx_xz[k] = -g_z_y_xxxx_xz[k] * ab_x + g_z_y_xxxx_xxz[k];

                g_z_y_xxxxx_yy[k] = -g_z_y_xxxx_yy[k] * ab_x + g_z_y_xxxx_xyy[k];

                g_z_y_xxxxx_yz[k] = -g_z_y_xxxx_yz[k] * ab_x + g_z_y_xxxx_xyz[k];

                g_z_y_xxxxx_zz[k] = -g_z_y_xxxx_zz[k] * ab_x + g_z_y_xxxx_xzz[k];
            }

            /// Set up 888-894 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxy_xx = cbuffer.data(hd_geom_11_off + 888 * ccomps * dcomps);

            auto g_z_y_xxxxy_xy = cbuffer.data(hd_geom_11_off + 889 * ccomps * dcomps);

            auto g_z_y_xxxxy_xz = cbuffer.data(hd_geom_11_off + 890 * ccomps * dcomps);

            auto g_z_y_xxxxy_yy = cbuffer.data(hd_geom_11_off + 891 * ccomps * dcomps);

            auto g_z_y_xxxxy_yz = cbuffer.data(hd_geom_11_off + 892 * ccomps * dcomps);

            auto g_z_y_xxxxy_zz = cbuffer.data(hd_geom_11_off + 893 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxxy_xx, g_z_y_xxxxy_xy, g_z_y_xxxxy_xz, g_z_y_xxxxy_yy, g_z_y_xxxxy_yz, g_z_y_xxxxy_zz, g_z_y_xxxy_xx, g_z_y_xxxy_xxx, g_z_y_xxxy_xxy, g_z_y_xxxy_xxz, g_z_y_xxxy_xy, g_z_y_xxxy_xyy, g_z_y_xxxy_xyz, g_z_y_xxxy_xz, g_z_y_xxxy_xzz, g_z_y_xxxy_yy, g_z_y_xxxy_yz, g_z_y_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxy_xx[k] = -g_z_y_xxxy_xx[k] * ab_x + g_z_y_xxxy_xxx[k];

                g_z_y_xxxxy_xy[k] = -g_z_y_xxxy_xy[k] * ab_x + g_z_y_xxxy_xxy[k];

                g_z_y_xxxxy_xz[k] = -g_z_y_xxxy_xz[k] * ab_x + g_z_y_xxxy_xxz[k];

                g_z_y_xxxxy_yy[k] = -g_z_y_xxxy_yy[k] * ab_x + g_z_y_xxxy_xyy[k];

                g_z_y_xxxxy_yz[k] = -g_z_y_xxxy_yz[k] * ab_x + g_z_y_xxxy_xyz[k];

                g_z_y_xxxxy_zz[k] = -g_z_y_xxxy_zz[k] * ab_x + g_z_y_xxxy_xzz[k];
            }

            /// Set up 894-900 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxxz_xx = cbuffer.data(hd_geom_11_off + 894 * ccomps * dcomps);

            auto g_z_y_xxxxz_xy = cbuffer.data(hd_geom_11_off + 895 * ccomps * dcomps);

            auto g_z_y_xxxxz_xz = cbuffer.data(hd_geom_11_off + 896 * ccomps * dcomps);

            auto g_z_y_xxxxz_yy = cbuffer.data(hd_geom_11_off + 897 * ccomps * dcomps);

            auto g_z_y_xxxxz_yz = cbuffer.data(hd_geom_11_off + 898 * ccomps * dcomps);

            auto g_z_y_xxxxz_zz = cbuffer.data(hd_geom_11_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxxz_xx, g_z_y_xxxxz_xy, g_z_y_xxxxz_xz, g_z_y_xxxxz_yy, g_z_y_xxxxz_yz, g_z_y_xxxxz_zz, g_z_y_xxxz_xx, g_z_y_xxxz_xxx, g_z_y_xxxz_xxy, g_z_y_xxxz_xxz, g_z_y_xxxz_xy, g_z_y_xxxz_xyy, g_z_y_xxxz_xyz, g_z_y_xxxz_xz, g_z_y_xxxz_xzz, g_z_y_xxxz_yy, g_z_y_xxxz_yz, g_z_y_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxxz_xx[k] = -g_z_y_xxxz_xx[k] * ab_x + g_z_y_xxxz_xxx[k];

                g_z_y_xxxxz_xy[k] = -g_z_y_xxxz_xy[k] * ab_x + g_z_y_xxxz_xxy[k];

                g_z_y_xxxxz_xz[k] = -g_z_y_xxxz_xz[k] * ab_x + g_z_y_xxxz_xxz[k];

                g_z_y_xxxxz_yy[k] = -g_z_y_xxxz_yy[k] * ab_x + g_z_y_xxxz_xyy[k];

                g_z_y_xxxxz_yz[k] = -g_z_y_xxxz_yz[k] * ab_x + g_z_y_xxxz_xyz[k];

                g_z_y_xxxxz_zz[k] = -g_z_y_xxxz_zz[k] * ab_x + g_z_y_xxxz_xzz[k];
            }

            /// Set up 900-906 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxyy_xx = cbuffer.data(hd_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_y_xxxyy_xy = cbuffer.data(hd_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_y_xxxyy_xz = cbuffer.data(hd_geom_11_off + 902 * ccomps * dcomps);

            auto g_z_y_xxxyy_yy = cbuffer.data(hd_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_y_xxxyy_yz = cbuffer.data(hd_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_y_xxxyy_zz = cbuffer.data(hd_geom_11_off + 905 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxyy_xx, g_z_y_xxxyy_xy, g_z_y_xxxyy_xz, g_z_y_xxxyy_yy, g_z_y_xxxyy_yz, g_z_y_xxxyy_zz, g_z_y_xxyy_xx, g_z_y_xxyy_xxx, g_z_y_xxyy_xxy, g_z_y_xxyy_xxz, g_z_y_xxyy_xy, g_z_y_xxyy_xyy, g_z_y_xxyy_xyz, g_z_y_xxyy_xz, g_z_y_xxyy_xzz, g_z_y_xxyy_yy, g_z_y_xxyy_yz, g_z_y_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxyy_xx[k] = -g_z_y_xxyy_xx[k] * ab_x + g_z_y_xxyy_xxx[k];

                g_z_y_xxxyy_xy[k] = -g_z_y_xxyy_xy[k] * ab_x + g_z_y_xxyy_xxy[k];

                g_z_y_xxxyy_xz[k] = -g_z_y_xxyy_xz[k] * ab_x + g_z_y_xxyy_xxz[k];

                g_z_y_xxxyy_yy[k] = -g_z_y_xxyy_yy[k] * ab_x + g_z_y_xxyy_xyy[k];

                g_z_y_xxxyy_yz[k] = -g_z_y_xxyy_yz[k] * ab_x + g_z_y_xxyy_xyz[k];

                g_z_y_xxxyy_zz[k] = -g_z_y_xxyy_zz[k] * ab_x + g_z_y_xxyy_xzz[k];
            }

            /// Set up 906-912 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxyz_xx = cbuffer.data(hd_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_y_xxxyz_xy = cbuffer.data(hd_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_y_xxxyz_xz = cbuffer.data(hd_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_y_xxxyz_yy = cbuffer.data(hd_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_y_xxxyz_yz = cbuffer.data(hd_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_y_xxxyz_zz = cbuffer.data(hd_geom_11_off + 911 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxyz_xx, g_z_y_xxxyz_xy, g_z_y_xxxyz_xz, g_z_y_xxxyz_yy, g_z_y_xxxyz_yz, g_z_y_xxxyz_zz, g_z_y_xxyz_xx, g_z_y_xxyz_xxx, g_z_y_xxyz_xxy, g_z_y_xxyz_xxz, g_z_y_xxyz_xy, g_z_y_xxyz_xyy, g_z_y_xxyz_xyz, g_z_y_xxyz_xz, g_z_y_xxyz_xzz, g_z_y_xxyz_yy, g_z_y_xxyz_yz, g_z_y_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxyz_xx[k] = -g_z_y_xxyz_xx[k] * ab_x + g_z_y_xxyz_xxx[k];

                g_z_y_xxxyz_xy[k] = -g_z_y_xxyz_xy[k] * ab_x + g_z_y_xxyz_xxy[k];

                g_z_y_xxxyz_xz[k] = -g_z_y_xxyz_xz[k] * ab_x + g_z_y_xxyz_xxz[k];

                g_z_y_xxxyz_yy[k] = -g_z_y_xxyz_yy[k] * ab_x + g_z_y_xxyz_xyy[k];

                g_z_y_xxxyz_yz[k] = -g_z_y_xxyz_yz[k] * ab_x + g_z_y_xxyz_xyz[k];

                g_z_y_xxxyz_zz[k] = -g_z_y_xxyz_zz[k] * ab_x + g_z_y_xxyz_xzz[k];
            }

            /// Set up 912-918 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxzz_xx = cbuffer.data(hd_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_y_xxxzz_xy = cbuffer.data(hd_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_y_xxxzz_xz = cbuffer.data(hd_geom_11_off + 914 * ccomps * dcomps);

            auto g_z_y_xxxzz_yy = cbuffer.data(hd_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_y_xxxzz_yz = cbuffer.data(hd_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_y_xxxzz_zz = cbuffer.data(hd_geom_11_off + 917 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxzz_xx, g_z_y_xxxzz_xy, g_z_y_xxxzz_xz, g_z_y_xxxzz_yy, g_z_y_xxxzz_yz, g_z_y_xxxzz_zz, g_z_y_xxzz_xx, g_z_y_xxzz_xxx, g_z_y_xxzz_xxy, g_z_y_xxzz_xxz, g_z_y_xxzz_xy, g_z_y_xxzz_xyy, g_z_y_xxzz_xyz, g_z_y_xxzz_xz, g_z_y_xxzz_xzz, g_z_y_xxzz_yy, g_z_y_xxzz_yz, g_z_y_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxzz_xx[k] = -g_z_y_xxzz_xx[k] * ab_x + g_z_y_xxzz_xxx[k];

                g_z_y_xxxzz_xy[k] = -g_z_y_xxzz_xy[k] * ab_x + g_z_y_xxzz_xxy[k];

                g_z_y_xxxzz_xz[k] = -g_z_y_xxzz_xz[k] * ab_x + g_z_y_xxzz_xxz[k];

                g_z_y_xxxzz_yy[k] = -g_z_y_xxzz_yy[k] * ab_x + g_z_y_xxzz_xyy[k];

                g_z_y_xxxzz_yz[k] = -g_z_y_xxzz_yz[k] * ab_x + g_z_y_xxzz_xyz[k];

                g_z_y_xxxzz_zz[k] = -g_z_y_xxzz_zz[k] * ab_x + g_z_y_xxzz_xzz[k];
            }

            /// Set up 918-924 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyyy_xx = cbuffer.data(hd_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_y_xxyyy_xy = cbuffer.data(hd_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_y_xxyyy_xz = cbuffer.data(hd_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_y_xxyyy_yy = cbuffer.data(hd_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_y_xxyyy_yz = cbuffer.data(hd_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_y_xxyyy_zz = cbuffer.data(hd_geom_11_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyyy_xx, g_z_y_xxyyy_xy, g_z_y_xxyyy_xz, g_z_y_xxyyy_yy, g_z_y_xxyyy_yz, g_z_y_xxyyy_zz, g_z_y_xyyy_xx, g_z_y_xyyy_xxx, g_z_y_xyyy_xxy, g_z_y_xyyy_xxz, g_z_y_xyyy_xy, g_z_y_xyyy_xyy, g_z_y_xyyy_xyz, g_z_y_xyyy_xz, g_z_y_xyyy_xzz, g_z_y_xyyy_yy, g_z_y_xyyy_yz, g_z_y_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyyy_xx[k] = -g_z_y_xyyy_xx[k] * ab_x + g_z_y_xyyy_xxx[k];

                g_z_y_xxyyy_xy[k] = -g_z_y_xyyy_xy[k] * ab_x + g_z_y_xyyy_xxy[k];

                g_z_y_xxyyy_xz[k] = -g_z_y_xyyy_xz[k] * ab_x + g_z_y_xyyy_xxz[k];

                g_z_y_xxyyy_yy[k] = -g_z_y_xyyy_yy[k] * ab_x + g_z_y_xyyy_xyy[k];

                g_z_y_xxyyy_yz[k] = -g_z_y_xyyy_yz[k] * ab_x + g_z_y_xyyy_xyz[k];

                g_z_y_xxyyy_zz[k] = -g_z_y_xyyy_zz[k] * ab_x + g_z_y_xyyy_xzz[k];
            }

            /// Set up 924-930 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyyz_xx = cbuffer.data(hd_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_y_xxyyz_xy = cbuffer.data(hd_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_y_xxyyz_xz = cbuffer.data(hd_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_y_xxyyz_yy = cbuffer.data(hd_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_y_xxyyz_yz = cbuffer.data(hd_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_y_xxyyz_zz = cbuffer.data(hd_geom_11_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyyz_xx, g_z_y_xxyyz_xy, g_z_y_xxyyz_xz, g_z_y_xxyyz_yy, g_z_y_xxyyz_yz, g_z_y_xxyyz_zz, g_z_y_xyyz_xx, g_z_y_xyyz_xxx, g_z_y_xyyz_xxy, g_z_y_xyyz_xxz, g_z_y_xyyz_xy, g_z_y_xyyz_xyy, g_z_y_xyyz_xyz, g_z_y_xyyz_xz, g_z_y_xyyz_xzz, g_z_y_xyyz_yy, g_z_y_xyyz_yz, g_z_y_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyyz_xx[k] = -g_z_y_xyyz_xx[k] * ab_x + g_z_y_xyyz_xxx[k];

                g_z_y_xxyyz_xy[k] = -g_z_y_xyyz_xy[k] * ab_x + g_z_y_xyyz_xxy[k];

                g_z_y_xxyyz_xz[k] = -g_z_y_xyyz_xz[k] * ab_x + g_z_y_xyyz_xxz[k];

                g_z_y_xxyyz_yy[k] = -g_z_y_xyyz_yy[k] * ab_x + g_z_y_xyyz_xyy[k];

                g_z_y_xxyyz_yz[k] = -g_z_y_xyyz_yz[k] * ab_x + g_z_y_xyyz_xyz[k];

                g_z_y_xxyyz_zz[k] = -g_z_y_xyyz_zz[k] * ab_x + g_z_y_xyyz_xzz[k];
            }

            /// Set up 930-936 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyzz_xx = cbuffer.data(hd_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_y_xxyzz_xy = cbuffer.data(hd_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_y_xxyzz_xz = cbuffer.data(hd_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_y_xxyzz_yy = cbuffer.data(hd_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_y_xxyzz_yz = cbuffer.data(hd_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_y_xxyzz_zz = cbuffer.data(hd_geom_11_off + 935 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyzz_xx, g_z_y_xxyzz_xy, g_z_y_xxyzz_xz, g_z_y_xxyzz_yy, g_z_y_xxyzz_yz, g_z_y_xxyzz_zz, g_z_y_xyzz_xx, g_z_y_xyzz_xxx, g_z_y_xyzz_xxy, g_z_y_xyzz_xxz, g_z_y_xyzz_xy, g_z_y_xyzz_xyy, g_z_y_xyzz_xyz, g_z_y_xyzz_xz, g_z_y_xyzz_xzz, g_z_y_xyzz_yy, g_z_y_xyzz_yz, g_z_y_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyzz_xx[k] = -g_z_y_xyzz_xx[k] * ab_x + g_z_y_xyzz_xxx[k];

                g_z_y_xxyzz_xy[k] = -g_z_y_xyzz_xy[k] * ab_x + g_z_y_xyzz_xxy[k];

                g_z_y_xxyzz_xz[k] = -g_z_y_xyzz_xz[k] * ab_x + g_z_y_xyzz_xxz[k];

                g_z_y_xxyzz_yy[k] = -g_z_y_xyzz_yy[k] * ab_x + g_z_y_xyzz_xyy[k];

                g_z_y_xxyzz_yz[k] = -g_z_y_xyzz_yz[k] * ab_x + g_z_y_xyzz_xyz[k];

                g_z_y_xxyzz_zz[k] = -g_z_y_xyzz_zz[k] * ab_x + g_z_y_xyzz_xzz[k];
            }

            /// Set up 936-942 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxzzz_xx = cbuffer.data(hd_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_y_xxzzz_xy = cbuffer.data(hd_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_y_xxzzz_xz = cbuffer.data(hd_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_y_xxzzz_yy = cbuffer.data(hd_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_y_xxzzz_yz = cbuffer.data(hd_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_y_xxzzz_zz = cbuffer.data(hd_geom_11_off + 941 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxzzz_xx, g_z_y_xxzzz_xy, g_z_y_xxzzz_xz, g_z_y_xxzzz_yy, g_z_y_xxzzz_yz, g_z_y_xxzzz_zz, g_z_y_xzzz_xx, g_z_y_xzzz_xxx, g_z_y_xzzz_xxy, g_z_y_xzzz_xxz, g_z_y_xzzz_xy, g_z_y_xzzz_xyy, g_z_y_xzzz_xyz, g_z_y_xzzz_xz, g_z_y_xzzz_xzz, g_z_y_xzzz_yy, g_z_y_xzzz_yz, g_z_y_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxzzz_xx[k] = -g_z_y_xzzz_xx[k] * ab_x + g_z_y_xzzz_xxx[k];

                g_z_y_xxzzz_xy[k] = -g_z_y_xzzz_xy[k] * ab_x + g_z_y_xzzz_xxy[k];

                g_z_y_xxzzz_xz[k] = -g_z_y_xzzz_xz[k] * ab_x + g_z_y_xzzz_xxz[k];

                g_z_y_xxzzz_yy[k] = -g_z_y_xzzz_yy[k] * ab_x + g_z_y_xzzz_xyy[k];

                g_z_y_xxzzz_yz[k] = -g_z_y_xzzz_yz[k] * ab_x + g_z_y_xzzz_xyz[k];

                g_z_y_xxzzz_zz[k] = -g_z_y_xzzz_zz[k] * ab_x + g_z_y_xzzz_xzz[k];
            }

            /// Set up 942-948 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyyy_xx = cbuffer.data(hd_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_y_xyyyy_xy = cbuffer.data(hd_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_y_xyyyy_xz = cbuffer.data(hd_geom_11_off + 944 * ccomps * dcomps);

            auto g_z_y_xyyyy_yy = cbuffer.data(hd_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_y_xyyyy_yz = cbuffer.data(hd_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_y_xyyyy_zz = cbuffer.data(hd_geom_11_off + 947 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyyy_xx, g_z_y_xyyyy_xy, g_z_y_xyyyy_xz, g_z_y_xyyyy_yy, g_z_y_xyyyy_yz, g_z_y_xyyyy_zz, g_z_y_yyyy_xx, g_z_y_yyyy_xxx, g_z_y_yyyy_xxy, g_z_y_yyyy_xxz, g_z_y_yyyy_xy, g_z_y_yyyy_xyy, g_z_y_yyyy_xyz, g_z_y_yyyy_xz, g_z_y_yyyy_xzz, g_z_y_yyyy_yy, g_z_y_yyyy_yz, g_z_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyyy_xx[k] = -g_z_y_yyyy_xx[k] * ab_x + g_z_y_yyyy_xxx[k];

                g_z_y_xyyyy_xy[k] = -g_z_y_yyyy_xy[k] * ab_x + g_z_y_yyyy_xxy[k];

                g_z_y_xyyyy_xz[k] = -g_z_y_yyyy_xz[k] * ab_x + g_z_y_yyyy_xxz[k];

                g_z_y_xyyyy_yy[k] = -g_z_y_yyyy_yy[k] * ab_x + g_z_y_yyyy_xyy[k];

                g_z_y_xyyyy_yz[k] = -g_z_y_yyyy_yz[k] * ab_x + g_z_y_yyyy_xyz[k];

                g_z_y_xyyyy_zz[k] = -g_z_y_yyyy_zz[k] * ab_x + g_z_y_yyyy_xzz[k];
            }

            /// Set up 948-954 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyyz_xx = cbuffer.data(hd_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_y_xyyyz_xy = cbuffer.data(hd_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_y_xyyyz_xz = cbuffer.data(hd_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_y_xyyyz_yy = cbuffer.data(hd_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_y_xyyyz_yz = cbuffer.data(hd_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_y_xyyyz_zz = cbuffer.data(hd_geom_11_off + 953 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyyz_xx, g_z_y_xyyyz_xy, g_z_y_xyyyz_xz, g_z_y_xyyyz_yy, g_z_y_xyyyz_yz, g_z_y_xyyyz_zz, g_z_y_yyyz_xx, g_z_y_yyyz_xxx, g_z_y_yyyz_xxy, g_z_y_yyyz_xxz, g_z_y_yyyz_xy, g_z_y_yyyz_xyy, g_z_y_yyyz_xyz, g_z_y_yyyz_xz, g_z_y_yyyz_xzz, g_z_y_yyyz_yy, g_z_y_yyyz_yz, g_z_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyyz_xx[k] = -g_z_y_yyyz_xx[k] * ab_x + g_z_y_yyyz_xxx[k];

                g_z_y_xyyyz_xy[k] = -g_z_y_yyyz_xy[k] * ab_x + g_z_y_yyyz_xxy[k];

                g_z_y_xyyyz_xz[k] = -g_z_y_yyyz_xz[k] * ab_x + g_z_y_yyyz_xxz[k];

                g_z_y_xyyyz_yy[k] = -g_z_y_yyyz_yy[k] * ab_x + g_z_y_yyyz_xyy[k];

                g_z_y_xyyyz_yz[k] = -g_z_y_yyyz_yz[k] * ab_x + g_z_y_yyyz_xyz[k];

                g_z_y_xyyyz_zz[k] = -g_z_y_yyyz_zz[k] * ab_x + g_z_y_yyyz_xzz[k];
            }

            /// Set up 954-960 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyzz_xx = cbuffer.data(hd_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_y_xyyzz_xy = cbuffer.data(hd_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_y_xyyzz_xz = cbuffer.data(hd_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_y_xyyzz_yy = cbuffer.data(hd_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_y_xyyzz_yz = cbuffer.data(hd_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_y_xyyzz_zz = cbuffer.data(hd_geom_11_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyzz_xx, g_z_y_xyyzz_xy, g_z_y_xyyzz_xz, g_z_y_xyyzz_yy, g_z_y_xyyzz_yz, g_z_y_xyyzz_zz, g_z_y_yyzz_xx, g_z_y_yyzz_xxx, g_z_y_yyzz_xxy, g_z_y_yyzz_xxz, g_z_y_yyzz_xy, g_z_y_yyzz_xyy, g_z_y_yyzz_xyz, g_z_y_yyzz_xz, g_z_y_yyzz_xzz, g_z_y_yyzz_yy, g_z_y_yyzz_yz, g_z_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyzz_xx[k] = -g_z_y_yyzz_xx[k] * ab_x + g_z_y_yyzz_xxx[k];

                g_z_y_xyyzz_xy[k] = -g_z_y_yyzz_xy[k] * ab_x + g_z_y_yyzz_xxy[k];

                g_z_y_xyyzz_xz[k] = -g_z_y_yyzz_xz[k] * ab_x + g_z_y_yyzz_xxz[k];

                g_z_y_xyyzz_yy[k] = -g_z_y_yyzz_yy[k] * ab_x + g_z_y_yyzz_xyy[k];

                g_z_y_xyyzz_yz[k] = -g_z_y_yyzz_yz[k] * ab_x + g_z_y_yyzz_xyz[k];

                g_z_y_xyyzz_zz[k] = -g_z_y_yyzz_zz[k] * ab_x + g_z_y_yyzz_xzz[k];
            }

            /// Set up 960-966 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyzzz_xx = cbuffer.data(hd_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_y_xyzzz_xy = cbuffer.data(hd_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_y_xyzzz_xz = cbuffer.data(hd_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_y_xyzzz_yy = cbuffer.data(hd_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_y_xyzzz_yz = cbuffer.data(hd_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_y_xyzzz_zz = cbuffer.data(hd_geom_11_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyzzz_xx, g_z_y_xyzzz_xy, g_z_y_xyzzz_xz, g_z_y_xyzzz_yy, g_z_y_xyzzz_yz, g_z_y_xyzzz_zz, g_z_y_yzzz_xx, g_z_y_yzzz_xxx, g_z_y_yzzz_xxy, g_z_y_yzzz_xxz, g_z_y_yzzz_xy, g_z_y_yzzz_xyy, g_z_y_yzzz_xyz, g_z_y_yzzz_xz, g_z_y_yzzz_xzz, g_z_y_yzzz_yy, g_z_y_yzzz_yz, g_z_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyzzz_xx[k] = -g_z_y_yzzz_xx[k] * ab_x + g_z_y_yzzz_xxx[k];

                g_z_y_xyzzz_xy[k] = -g_z_y_yzzz_xy[k] * ab_x + g_z_y_yzzz_xxy[k];

                g_z_y_xyzzz_xz[k] = -g_z_y_yzzz_xz[k] * ab_x + g_z_y_yzzz_xxz[k];

                g_z_y_xyzzz_yy[k] = -g_z_y_yzzz_yy[k] * ab_x + g_z_y_yzzz_xyy[k];

                g_z_y_xyzzz_yz[k] = -g_z_y_yzzz_yz[k] * ab_x + g_z_y_yzzz_xyz[k];

                g_z_y_xyzzz_zz[k] = -g_z_y_yzzz_zz[k] * ab_x + g_z_y_yzzz_xzz[k];
            }

            /// Set up 966-972 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzzzz_xx = cbuffer.data(hd_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_y_xzzzz_xy = cbuffer.data(hd_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_y_xzzzz_xz = cbuffer.data(hd_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_y_xzzzz_yy = cbuffer.data(hd_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_y_xzzzz_yz = cbuffer.data(hd_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_y_xzzzz_zz = cbuffer.data(hd_geom_11_off + 971 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzzzz_xx, g_z_y_xzzzz_xy, g_z_y_xzzzz_xz, g_z_y_xzzzz_yy, g_z_y_xzzzz_yz, g_z_y_xzzzz_zz, g_z_y_zzzz_xx, g_z_y_zzzz_xxx, g_z_y_zzzz_xxy, g_z_y_zzzz_xxz, g_z_y_zzzz_xy, g_z_y_zzzz_xyy, g_z_y_zzzz_xyz, g_z_y_zzzz_xz, g_z_y_zzzz_xzz, g_z_y_zzzz_yy, g_z_y_zzzz_yz, g_z_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzzzz_xx[k] = -g_z_y_zzzz_xx[k] * ab_x + g_z_y_zzzz_xxx[k];

                g_z_y_xzzzz_xy[k] = -g_z_y_zzzz_xy[k] * ab_x + g_z_y_zzzz_xxy[k];

                g_z_y_xzzzz_xz[k] = -g_z_y_zzzz_xz[k] * ab_x + g_z_y_zzzz_xxz[k];

                g_z_y_xzzzz_yy[k] = -g_z_y_zzzz_yy[k] * ab_x + g_z_y_zzzz_xyy[k];

                g_z_y_xzzzz_yz[k] = -g_z_y_zzzz_yz[k] * ab_x + g_z_y_zzzz_xyz[k];

                g_z_y_xzzzz_zz[k] = -g_z_y_zzzz_zz[k] * ab_x + g_z_y_zzzz_xzz[k];
            }

            /// Set up 972-978 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyyy_xx = cbuffer.data(hd_geom_11_off + 972 * ccomps * dcomps);

            auto g_z_y_yyyyy_xy = cbuffer.data(hd_geom_11_off + 973 * ccomps * dcomps);

            auto g_z_y_yyyyy_xz = cbuffer.data(hd_geom_11_off + 974 * ccomps * dcomps);

            auto g_z_y_yyyyy_yy = cbuffer.data(hd_geom_11_off + 975 * ccomps * dcomps);

            auto g_z_y_yyyyy_yz = cbuffer.data(hd_geom_11_off + 976 * ccomps * dcomps);

            auto g_z_y_yyyyy_zz = cbuffer.data(hd_geom_11_off + 977 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyy_xx, g_z_0_yyyy_xy, g_z_0_yyyy_xz, g_z_0_yyyy_yy, g_z_0_yyyy_yz, g_z_0_yyyy_zz, g_z_y_yyyy_xx, g_z_y_yyyy_xxy, g_z_y_yyyy_xy, g_z_y_yyyy_xyy, g_z_y_yyyy_xyz, g_z_y_yyyy_xz, g_z_y_yyyy_yy, g_z_y_yyyy_yyy, g_z_y_yyyy_yyz, g_z_y_yyyy_yz, g_z_y_yyyy_yzz, g_z_y_yyyy_zz, g_z_y_yyyyy_xx, g_z_y_yyyyy_xy, g_z_y_yyyyy_xz, g_z_y_yyyyy_yy, g_z_y_yyyyy_yz, g_z_y_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyyy_xx[k] = g_z_0_yyyy_xx[k] - g_z_y_yyyy_xx[k] * ab_y + g_z_y_yyyy_xxy[k];

                g_z_y_yyyyy_xy[k] = g_z_0_yyyy_xy[k] - g_z_y_yyyy_xy[k] * ab_y + g_z_y_yyyy_xyy[k];

                g_z_y_yyyyy_xz[k] = g_z_0_yyyy_xz[k] - g_z_y_yyyy_xz[k] * ab_y + g_z_y_yyyy_xyz[k];

                g_z_y_yyyyy_yy[k] = g_z_0_yyyy_yy[k] - g_z_y_yyyy_yy[k] * ab_y + g_z_y_yyyy_yyy[k];

                g_z_y_yyyyy_yz[k] = g_z_0_yyyy_yz[k] - g_z_y_yyyy_yz[k] * ab_y + g_z_y_yyyy_yyz[k];

                g_z_y_yyyyy_zz[k] = g_z_0_yyyy_zz[k] - g_z_y_yyyy_zz[k] * ab_y + g_z_y_yyyy_yzz[k];
            }

            /// Set up 978-984 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyyz_xx = cbuffer.data(hd_geom_11_off + 978 * ccomps * dcomps);

            auto g_z_y_yyyyz_xy = cbuffer.data(hd_geom_11_off + 979 * ccomps * dcomps);

            auto g_z_y_yyyyz_xz = cbuffer.data(hd_geom_11_off + 980 * ccomps * dcomps);

            auto g_z_y_yyyyz_yy = cbuffer.data(hd_geom_11_off + 981 * ccomps * dcomps);

            auto g_z_y_yyyyz_yz = cbuffer.data(hd_geom_11_off + 982 * ccomps * dcomps);

            auto g_z_y_yyyyz_zz = cbuffer.data(hd_geom_11_off + 983 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_xx, g_z_0_yyyz_xy, g_z_0_yyyz_xz, g_z_0_yyyz_yy, g_z_0_yyyz_yz, g_z_0_yyyz_zz, g_z_y_yyyyz_xx, g_z_y_yyyyz_xy, g_z_y_yyyyz_xz, g_z_y_yyyyz_yy, g_z_y_yyyyz_yz, g_z_y_yyyyz_zz, g_z_y_yyyz_xx, g_z_y_yyyz_xxy, g_z_y_yyyz_xy, g_z_y_yyyz_xyy, g_z_y_yyyz_xyz, g_z_y_yyyz_xz, g_z_y_yyyz_yy, g_z_y_yyyz_yyy, g_z_y_yyyz_yyz, g_z_y_yyyz_yz, g_z_y_yyyz_yzz, g_z_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyyz_xx[k] = g_z_0_yyyz_xx[k] - g_z_y_yyyz_xx[k] * ab_y + g_z_y_yyyz_xxy[k];

                g_z_y_yyyyz_xy[k] = g_z_0_yyyz_xy[k] - g_z_y_yyyz_xy[k] * ab_y + g_z_y_yyyz_xyy[k];

                g_z_y_yyyyz_xz[k] = g_z_0_yyyz_xz[k] - g_z_y_yyyz_xz[k] * ab_y + g_z_y_yyyz_xyz[k];

                g_z_y_yyyyz_yy[k] = g_z_0_yyyz_yy[k] - g_z_y_yyyz_yy[k] * ab_y + g_z_y_yyyz_yyy[k];

                g_z_y_yyyyz_yz[k] = g_z_0_yyyz_yz[k] - g_z_y_yyyz_yz[k] * ab_y + g_z_y_yyyz_yyz[k];

                g_z_y_yyyyz_zz[k] = g_z_0_yyyz_zz[k] - g_z_y_yyyz_zz[k] * ab_y + g_z_y_yyyz_yzz[k];
            }

            /// Set up 984-990 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyzz_xx = cbuffer.data(hd_geom_11_off + 984 * ccomps * dcomps);

            auto g_z_y_yyyzz_xy = cbuffer.data(hd_geom_11_off + 985 * ccomps * dcomps);

            auto g_z_y_yyyzz_xz = cbuffer.data(hd_geom_11_off + 986 * ccomps * dcomps);

            auto g_z_y_yyyzz_yy = cbuffer.data(hd_geom_11_off + 987 * ccomps * dcomps);

            auto g_z_y_yyyzz_yz = cbuffer.data(hd_geom_11_off + 988 * ccomps * dcomps);

            auto g_z_y_yyyzz_zz = cbuffer.data(hd_geom_11_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_xx, g_z_0_yyzz_xy, g_z_0_yyzz_xz, g_z_0_yyzz_yy, g_z_0_yyzz_yz, g_z_0_yyzz_zz, g_z_y_yyyzz_xx, g_z_y_yyyzz_xy, g_z_y_yyyzz_xz, g_z_y_yyyzz_yy, g_z_y_yyyzz_yz, g_z_y_yyyzz_zz, g_z_y_yyzz_xx, g_z_y_yyzz_xxy, g_z_y_yyzz_xy, g_z_y_yyzz_xyy, g_z_y_yyzz_xyz, g_z_y_yyzz_xz, g_z_y_yyzz_yy, g_z_y_yyzz_yyy, g_z_y_yyzz_yyz, g_z_y_yyzz_yz, g_z_y_yyzz_yzz, g_z_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyzz_xx[k] = g_z_0_yyzz_xx[k] - g_z_y_yyzz_xx[k] * ab_y + g_z_y_yyzz_xxy[k];

                g_z_y_yyyzz_xy[k] = g_z_0_yyzz_xy[k] - g_z_y_yyzz_xy[k] * ab_y + g_z_y_yyzz_xyy[k];

                g_z_y_yyyzz_xz[k] = g_z_0_yyzz_xz[k] - g_z_y_yyzz_xz[k] * ab_y + g_z_y_yyzz_xyz[k];

                g_z_y_yyyzz_yy[k] = g_z_0_yyzz_yy[k] - g_z_y_yyzz_yy[k] * ab_y + g_z_y_yyzz_yyy[k];

                g_z_y_yyyzz_yz[k] = g_z_0_yyzz_yz[k] - g_z_y_yyzz_yz[k] * ab_y + g_z_y_yyzz_yyz[k];

                g_z_y_yyyzz_zz[k] = g_z_0_yyzz_zz[k] - g_z_y_yyzz_zz[k] * ab_y + g_z_y_yyzz_yzz[k];
            }

            /// Set up 990-996 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyzzz_xx = cbuffer.data(hd_geom_11_off + 990 * ccomps * dcomps);

            auto g_z_y_yyzzz_xy = cbuffer.data(hd_geom_11_off + 991 * ccomps * dcomps);

            auto g_z_y_yyzzz_xz = cbuffer.data(hd_geom_11_off + 992 * ccomps * dcomps);

            auto g_z_y_yyzzz_yy = cbuffer.data(hd_geom_11_off + 993 * ccomps * dcomps);

            auto g_z_y_yyzzz_yz = cbuffer.data(hd_geom_11_off + 994 * ccomps * dcomps);

            auto g_z_y_yyzzz_zz = cbuffer.data(hd_geom_11_off + 995 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_xx, g_z_0_yzzz_xy, g_z_0_yzzz_xz, g_z_0_yzzz_yy, g_z_0_yzzz_yz, g_z_0_yzzz_zz, g_z_y_yyzzz_xx, g_z_y_yyzzz_xy, g_z_y_yyzzz_xz, g_z_y_yyzzz_yy, g_z_y_yyzzz_yz, g_z_y_yyzzz_zz, g_z_y_yzzz_xx, g_z_y_yzzz_xxy, g_z_y_yzzz_xy, g_z_y_yzzz_xyy, g_z_y_yzzz_xyz, g_z_y_yzzz_xz, g_z_y_yzzz_yy, g_z_y_yzzz_yyy, g_z_y_yzzz_yyz, g_z_y_yzzz_yz, g_z_y_yzzz_yzz, g_z_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyzzz_xx[k] = g_z_0_yzzz_xx[k] - g_z_y_yzzz_xx[k] * ab_y + g_z_y_yzzz_xxy[k];

                g_z_y_yyzzz_xy[k] = g_z_0_yzzz_xy[k] - g_z_y_yzzz_xy[k] * ab_y + g_z_y_yzzz_xyy[k];

                g_z_y_yyzzz_xz[k] = g_z_0_yzzz_xz[k] - g_z_y_yzzz_xz[k] * ab_y + g_z_y_yzzz_xyz[k];

                g_z_y_yyzzz_yy[k] = g_z_0_yzzz_yy[k] - g_z_y_yzzz_yy[k] * ab_y + g_z_y_yzzz_yyy[k];

                g_z_y_yyzzz_yz[k] = g_z_0_yzzz_yz[k] - g_z_y_yzzz_yz[k] * ab_y + g_z_y_yzzz_yyz[k];

                g_z_y_yyzzz_zz[k] = g_z_0_yzzz_zz[k] - g_z_y_yzzz_zz[k] * ab_y + g_z_y_yzzz_yzz[k];
            }

            /// Set up 996-1002 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzzzz_xx = cbuffer.data(hd_geom_11_off + 996 * ccomps * dcomps);

            auto g_z_y_yzzzz_xy = cbuffer.data(hd_geom_11_off + 997 * ccomps * dcomps);

            auto g_z_y_yzzzz_xz = cbuffer.data(hd_geom_11_off + 998 * ccomps * dcomps);

            auto g_z_y_yzzzz_yy = cbuffer.data(hd_geom_11_off + 999 * ccomps * dcomps);

            auto g_z_y_yzzzz_yz = cbuffer.data(hd_geom_11_off + 1000 * ccomps * dcomps);

            auto g_z_y_yzzzz_zz = cbuffer.data(hd_geom_11_off + 1001 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz, g_z_y_yzzzz_xx, g_z_y_yzzzz_xy, g_z_y_yzzzz_xz, g_z_y_yzzzz_yy, g_z_y_yzzzz_yz, g_z_y_yzzzz_zz, g_z_y_zzzz_xx, g_z_y_zzzz_xxy, g_z_y_zzzz_xy, g_z_y_zzzz_xyy, g_z_y_zzzz_xyz, g_z_y_zzzz_xz, g_z_y_zzzz_yy, g_z_y_zzzz_yyy, g_z_y_zzzz_yyz, g_z_y_zzzz_yz, g_z_y_zzzz_yzz, g_z_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzzzz_xx[k] = g_z_0_zzzz_xx[k] - g_z_y_zzzz_xx[k] * ab_y + g_z_y_zzzz_xxy[k];

                g_z_y_yzzzz_xy[k] = g_z_0_zzzz_xy[k] - g_z_y_zzzz_xy[k] * ab_y + g_z_y_zzzz_xyy[k];

                g_z_y_yzzzz_xz[k] = g_z_0_zzzz_xz[k] - g_z_y_zzzz_xz[k] * ab_y + g_z_y_zzzz_xyz[k];

                g_z_y_yzzzz_yy[k] = g_z_0_zzzz_yy[k] - g_z_y_zzzz_yy[k] * ab_y + g_z_y_zzzz_yyy[k];

                g_z_y_yzzzz_yz[k] = g_z_0_zzzz_yz[k] - g_z_y_zzzz_yz[k] * ab_y + g_z_y_zzzz_yyz[k];

                g_z_y_yzzzz_zz[k] = g_z_0_zzzz_zz[k] - g_z_y_zzzz_zz[k] * ab_y + g_z_y_zzzz_yzz[k];
            }

            /// Set up 1002-1008 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzzzz_xx = cbuffer.data(hd_geom_11_off + 1002 * ccomps * dcomps);

            auto g_z_y_zzzzz_xy = cbuffer.data(hd_geom_11_off + 1003 * ccomps * dcomps);

            auto g_z_y_zzzzz_xz = cbuffer.data(hd_geom_11_off + 1004 * ccomps * dcomps);

            auto g_z_y_zzzzz_yy = cbuffer.data(hd_geom_11_off + 1005 * ccomps * dcomps);

            auto g_z_y_zzzzz_yz = cbuffer.data(hd_geom_11_off + 1006 * ccomps * dcomps);

            auto g_z_y_zzzzz_zz = cbuffer.data(hd_geom_11_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_xx, g_0_y_zzzz_xy, g_0_y_zzzz_xz, g_0_y_zzzz_yy, g_0_y_zzzz_yz, g_0_y_zzzz_zz, g_z_y_zzzz_xx, g_z_y_zzzz_xxz, g_z_y_zzzz_xy, g_z_y_zzzz_xyz, g_z_y_zzzz_xz, g_z_y_zzzz_xzz, g_z_y_zzzz_yy, g_z_y_zzzz_yyz, g_z_y_zzzz_yz, g_z_y_zzzz_yzz, g_z_y_zzzz_zz, g_z_y_zzzz_zzz, g_z_y_zzzzz_xx, g_z_y_zzzzz_xy, g_z_y_zzzzz_xz, g_z_y_zzzzz_yy, g_z_y_zzzzz_yz, g_z_y_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzzzz_xx[k] = -g_0_y_zzzz_xx[k] - g_z_y_zzzz_xx[k] * ab_z + g_z_y_zzzz_xxz[k];

                g_z_y_zzzzz_xy[k] = -g_0_y_zzzz_xy[k] - g_z_y_zzzz_xy[k] * ab_z + g_z_y_zzzz_xyz[k];

                g_z_y_zzzzz_xz[k] = -g_0_y_zzzz_xz[k] - g_z_y_zzzz_xz[k] * ab_z + g_z_y_zzzz_xzz[k];

                g_z_y_zzzzz_yy[k] = -g_0_y_zzzz_yy[k] - g_z_y_zzzz_yy[k] * ab_z + g_z_y_zzzz_yyz[k];

                g_z_y_zzzzz_yz[k] = -g_0_y_zzzz_yz[k] - g_z_y_zzzz_yz[k] * ab_z + g_z_y_zzzz_yzz[k];

                g_z_y_zzzzz_zz[k] = -g_0_y_zzzz_zz[k] - g_z_y_zzzz_zz[k] * ab_z + g_z_y_zzzz_zzz[k];
            }

            /// Set up 1008-1014 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxx_xx = cbuffer.data(hd_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_z_xxxxx_xy = cbuffer.data(hd_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_z_xxxxx_xz = cbuffer.data(hd_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_z_xxxxx_yy = cbuffer.data(hd_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_z_xxxxx_yz = cbuffer.data(hd_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_z_xxxxx_zz = cbuffer.data(hd_geom_11_off + 1013 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxx_xx, g_z_z_xxxx_xxx, g_z_z_xxxx_xxy, g_z_z_xxxx_xxz, g_z_z_xxxx_xy, g_z_z_xxxx_xyy, g_z_z_xxxx_xyz, g_z_z_xxxx_xz, g_z_z_xxxx_xzz, g_z_z_xxxx_yy, g_z_z_xxxx_yz, g_z_z_xxxx_zz, g_z_z_xxxxx_xx, g_z_z_xxxxx_xy, g_z_z_xxxxx_xz, g_z_z_xxxxx_yy, g_z_z_xxxxx_yz, g_z_z_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxx_xx[k] = -g_z_z_xxxx_xx[k] * ab_x + g_z_z_xxxx_xxx[k];

                g_z_z_xxxxx_xy[k] = -g_z_z_xxxx_xy[k] * ab_x + g_z_z_xxxx_xxy[k];

                g_z_z_xxxxx_xz[k] = -g_z_z_xxxx_xz[k] * ab_x + g_z_z_xxxx_xxz[k];

                g_z_z_xxxxx_yy[k] = -g_z_z_xxxx_yy[k] * ab_x + g_z_z_xxxx_xyy[k];

                g_z_z_xxxxx_yz[k] = -g_z_z_xxxx_yz[k] * ab_x + g_z_z_xxxx_xyz[k];

                g_z_z_xxxxx_zz[k] = -g_z_z_xxxx_zz[k] * ab_x + g_z_z_xxxx_xzz[k];
            }

            /// Set up 1014-1020 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxy_xx = cbuffer.data(hd_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_z_xxxxy_xy = cbuffer.data(hd_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_z_xxxxy_xz = cbuffer.data(hd_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_z_xxxxy_yy = cbuffer.data(hd_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_z_xxxxy_yz = cbuffer.data(hd_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_z_xxxxy_zz = cbuffer.data(hd_geom_11_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxxy_xx, g_z_z_xxxxy_xy, g_z_z_xxxxy_xz, g_z_z_xxxxy_yy, g_z_z_xxxxy_yz, g_z_z_xxxxy_zz, g_z_z_xxxy_xx, g_z_z_xxxy_xxx, g_z_z_xxxy_xxy, g_z_z_xxxy_xxz, g_z_z_xxxy_xy, g_z_z_xxxy_xyy, g_z_z_xxxy_xyz, g_z_z_xxxy_xz, g_z_z_xxxy_xzz, g_z_z_xxxy_yy, g_z_z_xxxy_yz, g_z_z_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxy_xx[k] = -g_z_z_xxxy_xx[k] * ab_x + g_z_z_xxxy_xxx[k];

                g_z_z_xxxxy_xy[k] = -g_z_z_xxxy_xy[k] * ab_x + g_z_z_xxxy_xxy[k];

                g_z_z_xxxxy_xz[k] = -g_z_z_xxxy_xz[k] * ab_x + g_z_z_xxxy_xxz[k];

                g_z_z_xxxxy_yy[k] = -g_z_z_xxxy_yy[k] * ab_x + g_z_z_xxxy_xyy[k];

                g_z_z_xxxxy_yz[k] = -g_z_z_xxxy_yz[k] * ab_x + g_z_z_xxxy_xyz[k];

                g_z_z_xxxxy_zz[k] = -g_z_z_xxxy_zz[k] * ab_x + g_z_z_xxxy_xzz[k];
            }

            /// Set up 1020-1026 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxxz_xx = cbuffer.data(hd_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_z_xxxxz_xy = cbuffer.data(hd_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_z_xxxxz_xz = cbuffer.data(hd_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_z_xxxxz_yy = cbuffer.data(hd_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_z_xxxxz_yz = cbuffer.data(hd_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_z_xxxxz_zz = cbuffer.data(hd_geom_11_off + 1025 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxxz_xx, g_z_z_xxxxz_xy, g_z_z_xxxxz_xz, g_z_z_xxxxz_yy, g_z_z_xxxxz_yz, g_z_z_xxxxz_zz, g_z_z_xxxz_xx, g_z_z_xxxz_xxx, g_z_z_xxxz_xxy, g_z_z_xxxz_xxz, g_z_z_xxxz_xy, g_z_z_xxxz_xyy, g_z_z_xxxz_xyz, g_z_z_xxxz_xz, g_z_z_xxxz_xzz, g_z_z_xxxz_yy, g_z_z_xxxz_yz, g_z_z_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxxz_xx[k] = -g_z_z_xxxz_xx[k] * ab_x + g_z_z_xxxz_xxx[k];

                g_z_z_xxxxz_xy[k] = -g_z_z_xxxz_xy[k] * ab_x + g_z_z_xxxz_xxy[k];

                g_z_z_xxxxz_xz[k] = -g_z_z_xxxz_xz[k] * ab_x + g_z_z_xxxz_xxz[k];

                g_z_z_xxxxz_yy[k] = -g_z_z_xxxz_yy[k] * ab_x + g_z_z_xxxz_xyy[k];

                g_z_z_xxxxz_yz[k] = -g_z_z_xxxz_yz[k] * ab_x + g_z_z_xxxz_xyz[k];

                g_z_z_xxxxz_zz[k] = -g_z_z_xxxz_zz[k] * ab_x + g_z_z_xxxz_xzz[k];
            }

            /// Set up 1026-1032 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxyy_xx = cbuffer.data(hd_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_z_xxxyy_xy = cbuffer.data(hd_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_z_xxxyy_xz = cbuffer.data(hd_geom_11_off + 1028 * ccomps * dcomps);

            auto g_z_z_xxxyy_yy = cbuffer.data(hd_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_z_xxxyy_yz = cbuffer.data(hd_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_z_xxxyy_zz = cbuffer.data(hd_geom_11_off + 1031 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxyy_xx, g_z_z_xxxyy_xy, g_z_z_xxxyy_xz, g_z_z_xxxyy_yy, g_z_z_xxxyy_yz, g_z_z_xxxyy_zz, g_z_z_xxyy_xx, g_z_z_xxyy_xxx, g_z_z_xxyy_xxy, g_z_z_xxyy_xxz, g_z_z_xxyy_xy, g_z_z_xxyy_xyy, g_z_z_xxyy_xyz, g_z_z_xxyy_xz, g_z_z_xxyy_xzz, g_z_z_xxyy_yy, g_z_z_xxyy_yz, g_z_z_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxyy_xx[k] = -g_z_z_xxyy_xx[k] * ab_x + g_z_z_xxyy_xxx[k];

                g_z_z_xxxyy_xy[k] = -g_z_z_xxyy_xy[k] * ab_x + g_z_z_xxyy_xxy[k];

                g_z_z_xxxyy_xz[k] = -g_z_z_xxyy_xz[k] * ab_x + g_z_z_xxyy_xxz[k];

                g_z_z_xxxyy_yy[k] = -g_z_z_xxyy_yy[k] * ab_x + g_z_z_xxyy_xyy[k];

                g_z_z_xxxyy_yz[k] = -g_z_z_xxyy_yz[k] * ab_x + g_z_z_xxyy_xyz[k];

                g_z_z_xxxyy_zz[k] = -g_z_z_xxyy_zz[k] * ab_x + g_z_z_xxyy_xzz[k];
            }

            /// Set up 1032-1038 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxyz_xx = cbuffer.data(hd_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_z_xxxyz_xy = cbuffer.data(hd_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_z_xxxyz_xz = cbuffer.data(hd_geom_11_off + 1034 * ccomps * dcomps);

            auto g_z_z_xxxyz_yy = cbuffer.data(hd_geom_11_off + 1035 * ccomps * dcomps);

            auto g_z_z_xxxyz_yz = cbuffer.data(hd_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_z_xxxyz_zz = cbuffer.data(hd_geom_11_off + 1037 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxyz_xx, g_z_z_xxxyz_xy, g_z_z_xxxyz_xz, g_z_z_xxxyz_yy, g_z_z_xxxyz_yz, g_z_z_xxxyz_zz, g_z_z_xxyz_xx, g_z_z_xxyz_xxx, g_z_z_xxyz_xxy, g_z_z_xxyz_xxz, g_z_z_xxyz_xy, g_z_z_xxyz_xyy, g_z_z_xxyz_xyz, g_z_z_xxyz_xz, g_z_z_xxyz_xzz, g_z_z_xxyz_yy, g_z_z_xxyz_yz, g_z_z_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxyz_xx[k] = -g_z_z_xxyz_xx[k] * ab_x + g_z_z_xxyz_xxx[k];

                g_z_z_xxxyz_xy[k] = -g_z_z_xxyz_xy[k] * ab_x + g_z_z_xxyz_xxy[k];

                g_z_z_xxxyz_xz[k] = -g_z_z_xxyz_xz[k] * ab_x + g_z_z_xxyz_xxz[k];

                g_z_z_xxxyz_yy[k] = -g_z_z_xxyz_yy[k] * ab_x + g_z_z_xxyz_xyy[k];

                g_z_z_xxxyz_yz[k] = -g_z_z_xxyz_yz[k] * ab_x + g_z_z_xxyz_xyz[k];

                g_z_z_xxxyz_zz[k] = -g_z_z_xxyz_zz[k] * ab_x + g_z_z_xxyz_xzz[k];
            }

            /// Set up 1038-1044 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxzz_xx = cbuffer.data(hd_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_z_xxxzz_xy = cbuffer.data(hd_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_z_xxxzz_xz = cbuffer.data(hd_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_z_xxxzz_yy = cbuffer.data(hd_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_z_xxxzz_yz = cbuffer.data(hd_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_z_xxxzz_zz = cbuffer.data(hd_geom_11_off + 1043 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxzz_xx, g_z_z_xxxzz_xy, g_z_z_xxxzz_xz, g_z_z_xxxzz_yy, g_z_z_xxxzz_yz, g_z_z_xxxzz_zz, g_z_z_xxzz_xx, g_z_z_xxzz_xxx, g_z_z_xxzz_xxy, g_z_z_xxzz_xxz, g_z_z_xxzz_xy, g_z_z_xxzz_xyy, g_z_z_xxzz_xyz, g_z_z_xxzz_xz, g_z_z_xxzz_xzz, g_z_z_xxzz_yy, g_z_z_xxzz_yz, g_z_z_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxzz_xx[k] = -g_z_z_xxzz_xx[k] * ab_x + g_z_z_xxzz_xxx[k];

                g_z_z_xxxzz_xy[k] = -g_z_z_xxzz_xy[k] * ab_x + g_z_z_xxzz_xxy[k];

                g_z_z_xxxzz_xz[k] = -g_z_z_xxzz_xz[k] * ab_x + g_z_z_xxzz_xxz[k];

                g_z_z_xxxzz_yy[k] = -g_z_z_xxzz_yy[k] * ab_x + g_z_z_xxzz_xyy[k];

                g_z_z_xxxzz_yz[k] = -g_z_z_xxzz_yz[k] * ab_x + g_z_z_xxzz_xyz[k];

                g_z_z_xxxzz_zz[k] = -g_z_z_xxzz_zz[k] * ab_x + g_z_z_xxzz_xzz[k];
            }

            /// Set up 1044-1050 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyyy_xx = cbuffer.data(hd_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_z_xxyyy_xy = cbuffer.data(hd_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_z_xxyyy_xz = cbuffer.data(hd_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_z_xxyyy_yy = cbuffer.data(hd_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_z_xxyyy_yz = cbuffer.data(hd_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_z_xxyyy_zz = cbuffer.data(hd_geom_11_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyyy_xx, g_z_z_xxyyy_xy, g_z_z_xxyyy_xz, g_z_z_xxyyy_yy, g_z_z_xxyyy_yz, g_z_z_xxyyy_zz, g_z_z_xyyy_xx, g_z_z_xyyy_xxx, g_z_z_xyyy_xxy, g_z_z_xyyy_xxz, g_z_z_xyyy_xy, g_z_z_xyyy_xyy, g_z_z_xyyy_xyz, g_z_z_xyyy_xz, g_z_z_xyyy_xzz, g_z_z_xyyy_yy, g_z_z_xyyy_yz, g_z_z_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyyy_xx[k] = -g_z_z_xyyy_xx[k] * ab_x + g_z_z_xyyy_xxx[k];

                g_z_z_xxyyy_xy[k] = -g_z_z_xyyy_xy[k] * ab_x + g_z_z_xyyy_xxy[k];

                g_z_z_xxyyy_xz[k] = -g_z_z_xyyy_xz[k] * ab_x + g_z_z_xyyy_xxz[k];

                g_z_z_xxyyy_yy[k] = -g_z_z_xyyy_yy[k] * ab_x + g_z_z_xyyy_xyy[k];

                g_z_z_xxyyy_yz[k] = -g_z_z_xyyy_yz[k] * ab_x + g_z_z_xyyy_xyz[k];

                g_z_z_xxyyy_zz[k] = -g_z_z_xyyy_zz[k] * ab_x + g_z_z_xyyy_xzz[k];
            }

            /// Set up 1050-1056 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyyz_xx = cbuffer.data(hd_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_z_xxyyz_xy = cbuffer.data(hd_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_z_xxyyz_xz = cbuffer.data(hd_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_z_xxyyz_yy = cbuffer.data(hd_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_z_xxyyz_yz = cbuffer.data(hd_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_z_xxyyz_zz = cbuffer.data(hd_geom_11_off + 1055 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyyz_xx, g_z_z_xxyyz_xy, g_z_z_xxyyz_xz, g_z_z_xxyyz_yy, g_z_z_xxyyz_yz, g_z_z_xxyyz_zz, g_z_z_xyyz_xx, g_z_z_xyyz_xxx, g_z_z_xyyz_xxy, g_z_z_xyyz_xxz, g_z_z_xyyz_xy, g_z_z_xyyz_xyy, g_z_z_xyyz_xyz, g_z_z_xyyz_xz, g_z_z_xyyz_xzz, g_z_z_xyyz_yy, g_z_z_xyyz_yz, g_z_z_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyyz_xx[k] = -g_z_z_xyyz_xx[k] * ab_x + g_z_z_xyyz_xxx[k];

                g_z_z_xxyyz_xy[k] = -g_z_z_xyyz_xy[k] * ab_x + g_z_z_xyyz_xxy[k];

                g_z_z_xxyyz_xz[k] = -g_z_z_xyyz_xz[k] * ab_x + g_z_z_xyyz_xxz[k];

                g_z_z_xxyyz_yy[k] = -g_z_z_xyyz_yy[k] * ab_x + g_z_z_xyyz_xyy[k];

                g_z_z_xxyyz_yz[k] = -g_z_z_xyyz_yz[k] * ab_x + g_z_z_xyyz_xyz[k];

                g_z_z_xxyyz_zz[k] = -g_z_z_xyyz_zz[k] * ab_x + g_z_z_xyyz_xzz[k];
            }

            /// Set up 1056-1062 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyzz_xx = cbuffer.data(hd_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_z_xxyzz_xy = cbuffer.data(hd_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_z_xxyzz_xz = cbuffer.data(hd_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_z_xxyzz_yy = cbuffer.data(hd_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_z_xxyzz_yz = cbuffer.data(hd_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_z_xxyzz_zz = cbuffer.data(hd_geom_11_off + 1061 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyzz_xx, g_z_z_xxyzz_xy, g_z_z_xxyzz_xz, g_z_z_xxyzz_yy, g_z_z_xxyzz_yz, g_z_z_xxyzz_zz, g_z_z_xyzz_xx, g_z_z_xyzz_xxx, g_z_z_xyzz_xxy, g_z_z_xyzz_xxz, g_z_z_xyzz_xy, g_z_z_xyzz_xyy, g_z_z_xyzz_xyz, g_z_z_xyzz_xz, g_z_z_xyzz_xzz, g_z_z_xyzz_yy, g_z_z_xyzz_yz, g_z_z_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyzz_xx[k] = -g_z_z_xyzz_xx[k] * ab_x + g_z_z_xyzz_xxx[k];

                g_z_z_xxyzz_xy[k] = -g_z_z_xyzz_xy[k] * ab_x + g_z_z_xyzz_xxy[k];

                g_z_z_xxyzz_xz[k] = -g_z_z_xyzz_xz[k] * ab_x + g_z_z_xyzz_xxz[k];

                g_z_z_xxyzz_yy[k] = -g_z_z_xyzz_yy[k] * ab_x + g_z_z_xyzz_xyy[k];

                g_z_z_xxyzz_yz[k] = -g_z_z_xyzz_yz[k] * ab_x + g_z_z_xyzz_xyz[k];

                g_z_z_xxyzz_zz[k] = -g_z_z_xyzz_zz[k] * ab_x + g_z_z_xyzz_xzz[k];
            }

            /// Set up 1062-1068 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxzzz_xx = cbuffer.data(hd_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_z_xxzzz_xy = cbuffer.data(hd_geom_11_off + 1063 * ccomps * dcomps);

            auto g_z_z_xxzzz_xz = cbuffer.data(hd_geom_11_off + 1064 * ccomps * dcomps);

            auto g_z_z_xxzzz_yy = cbuffer.data(hd_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_z_xxzzz_yz = cbuffer.data(hd_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_z_xxzzz_zz = cbuffer.data(hd_geom_11_off + 1067 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxzzz_xx, g_z_z_xxzzz_xy, g_z_z_xxzzz_xz, g_z_z_xxzzz_yy, g_z_z_xxzzz_yz, g_z_z_xxzzz_zz, g_z_z_xzzz_xx, g_z_z_xzzz_xxx, g_z_z_xzzz_xxy, g_z_z_xzzz_xxz, g_z_z_xzzz_xy, g_z_z_xzzz_xyy, g_z_z_xzzz_xyz, g_z_z_xzzz_xz, g_z_z_xzzz_xzz, g_z_z_xzzz_yy, g_z_z_xzzz_yz, g_z_z_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxzzz_xx[k] = -g_z_z_xzzz_xx[k] * ab_x + g_z_z_xzzz_xxx[k];

                g_z_z_xxzzz_xy[k] = -g_z_z_xzzz_xy[k] * ab_x + g_z_z_xzzz_xxy[k];

                g_z_z_xxzzz_xz[k] = -g_z_z_xzzz_xz[k] * ab_x + g_z_z_xzzz_xxz[k];

                g_z_z_xxzzz_yy[k] = -g_z_z_xzzz_yy[k] * ab_x + g_z_z_xzzz_xyy[k];

                g_z_z_xxzzz_yz[k] = -g_z_z_xzzz_yz[k] * ab_x + g_z_z_xzzz_xyz[k];

                g_z_z_xxzzz_zz[k] = -g_z_z_xzzz_zz[k] * ab_x + g_z_z_xzzz_xzz[k];
            }

            /// Set up 1068-1074 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyyy_xx = cbuffer.data(hd_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_z_xyyyy_xy = cbuffer.data(hd_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_z_xyyyy_xz = cbuffer.data(hd_geom_11_off + 1070 * ccomps * dcomps);

            auto g_z_z_xyyyy_yy = cbuffer.data(hd_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_z_xyyyy_yz = cbuffer.data(hd_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_z_xyyyy_zz = cbuffer.data(hd_geom_11_off + 1073 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyyy_xx, g_z_z_xyyyy_xy, g_z_z_xyyyy_xz, g_z_z_xyyyy_yy, g_z_z_xyyyy_yz, g_z_z_xyyyy_zz, g_z_z_yyyy_xx, g_z_z_yyyy_xxx, g_z_z_yyyy_xxy, g_z_z_yyyy_xxz, g_z_z_yyyy_xy, g_z_z_yyyy_xyy, g_z_z_yyyy_xyz, g_z_z_yyyy_xz, g_z_z_yyyy_xzz, g_z_z_yyyy_yy, g_z_z_yyyy_yz, g_z_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyyy_xx[k] = -g_z_z_yyyy_xx[k] * ab_x + g_z_z_yyyy_xxx[k];

                g_z_z_xyyyy_xy[k] = -g_z_z_yyyy_xy[k] * ab_x + g_z_z_yyyy_xxy[k];

                g_z_z_xyyyy_xz[k] = -g_z_z_yyyy_xz[k] * ab_x + g_z_z_yyyy_xxz[k];

                g_z_z_xyyyy_yy[k] = -g_z_z_yyyy_yy[k] * ab_x + g_z_z_yyyy_xyy[k];

                g_z_z_xyyyy_yz[k] = -g_z_z_yyyy_yz[k] * ab_x + g_z_z_yyyy_xyz[k];

                g_z_z_xyyyy_zz[k] = -g_z_z_yyyy_zz[k] * ab_x + g_z_z_yyyy_xzz[k];
            }

            /// Set up 1074-1080 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyyz_xx = cbuffer.data(hd_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_z_xyyyz_xy = cbuffer.data(hd_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_z_xyyyz_xz = cbuffer.data(hd_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_z_xyyyz_yy = cbuffer.data(hd_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_z_xyyyz_yz = cbuffer.data(hd_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_z_xyyyz_zz = cbuffer.data(hd_geom_11_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyyz_xx, g_z_z_xyyyz_xy, g_z_z_xyyyz_xz, g_z_z_xyyyz_yy, g_z_z_xyyyz_yz, g_z_z_xyyyz_zz, g_z_z_yyyz_xx, g_z_z_yyyz_xxx, g_z_z_yyyz_xxy, g_z_z_yyyz_xxz, g_z_z_yyyz_xy, g_z_z_yyyz_xyy, g_z_z_yyyz_xyz, g_z_z_yyyz_xz, g_z_z_yyyz_xzz, g_z_z_yyyz_yy, g_z_z_yyyz_yz, g_z_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyyz_xx[k] = -g_z_z_yyyz_xx[k] * ab_x + g_z_z_yyyz_xxx[k];

                g_z_z_xyyyz_xy[k] = -g_z_z_yyyz_xy[k] * ab_x + g_z_z_yyyz_xxy[k];

                g_z_z_xyyyz_xz[k] = -g_z_z_yyyz_xz[k] * ab_x + g_z_z_yyyz_xxz[k];

                g_z_z_xyyyz_yy[k] = -g_z_z_yyyz_yy[k] * ab_x + g_z_z_yyyz_xyy[k];

                g_z_z_xyyyz_yz[k] = -g_z_z_yyyz_yz[k] * ab_x + g_z_z_yyyz_xyz[k];

                g_z_z_xyyyz_zz[k] = -g_z_z_yyyz_zz[k] * ab_x + g_z_z_yyyz_xzz[k];
            }

            /// Set up 1080-1086 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyzz_xx = cbuffer.data(hd_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_z_xyyzz_xy = cbuffer.data(hd_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_z_xyyzz_xz = cbuffer.data(hd_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_z_xyyzz_yy = cbuffer.data(hd_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_z_xyyzz_yz = cbuffer.data(hd_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_z_xyyzz_zz = cbuffer.data(hd_geom_11_off + 1085 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyzz_xx, g_z_z_xyyzz_xy, g_z_z_xyyzz_xz, g_z_z_xyyzz_yy, g_z_z_xyyzz_yz, g_z_z_xyyzz_zz, g_z_z_yyzz_xx, g_z_z_yyzz_xxx, g_z_z_yyzz_xxy, g_z_z_yyzz_xxz, g_z_z_yyzz_xy, g_z_z_yyzz_xyy, g_z_z_yyzz_xyz, g_z_z_yyzz_xz, g_z_z_yyzz_xzz, g_z_z_yyzz_yy, g_z_z_yyzz_yz, g_z_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyzz_xx[k] = -g_z_z_yyzz_xx[k] * ab_x + g_z_z_yyzz_xxx[k];

                g_z_z_xyyzz_xy[k] = -g_z_z_yyzz_xy[k] * ab_x + g_z_z_yyzz_xxy[k];

                g_z_z_xyyzz_xz[k] = -g_z_z_yyzz_xz[k] * ab_x + g_z_z_yyzz_xxz[k];

                g_z_z_xyyzz_yy[k] = -g_z_z_yyzz_yy[k] * ab_x + g_z_z_yyzz_xyy[k];

                g_z_z_xyyzz_yz[k] = -g_z_z_yyzz_yz[k] * ab_x + g_z_z_yyzz_xyz[k];

                g_z_z_xyyzz_zz[k] = -g_z_z_yyzz_zz[k] * ab_x + g_z_z_yyzz_xzz[k];
            }

            /// Set up 1086-1092 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyzzz_xx = cbuffer.data(hd_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_z_xyzzz_xy = cbuffer.data(hd_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_z_xyzzz_xz = cbuffer.data(hd_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_z_xyzzz_yy = cbuffer.data(hd_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_z_xyzzz_yz = cbuffer.data(hd_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_z_xyzzz_zz = cbuffer.data(hd_geom_11_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyzzz_xx, g_z_z_xyzzz_xy, g_z_z_xyzzz_xz, g_z_z_xyzzz_yy, g_z_z_xyzzz_yz, g_z_z_xyzzz_zz, g_z_z_yzzz_xx, g_z_z_yzzz_xxx, g_z_z_yzzz_xxy, g_z_z_yzzz_xxz, g_z_z_yzzz_xy, g_z_z_yzzz_xyy, g_z_z_yzzz_xyz, g_z_z_yzzz_xz, g_z_z_yzzz_xzz, g_z_z_yzzz_yy, g_z_z_yzzz_yz, g_z_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyzzz_xx[k] = -g_z_z_yzzz_xx[k] * ab_x + g_z_z_yzzz_xxx[k];

                g_z_z_xyzzz_xy[k] = -g_z_z_yzzz_xy[k] * ab_x + g_z_z_yzzz_xxy[k];

                g_z_z_xyzzz_xz[k] = -g_z_z_yzzz_xz[k] * ab_x + g_z_z_yzzz_xxz[k];

                g_z_z_xyzzz_yy[k] = -g_z_z_yzzz_yy[k] * ab_x + g_z_z_yzzz_xyy[k];

                g_z_z_xyzzz_yz[k] = -g_z_z_yzzz_yz[k] * ab_x + g_z_z_yzzz_xyz[k];

                g_z_z_xyzzz_zz[k] = -g_z_z_yzzz_zz[k] * ab_x + g_z_z_yzzz_xzz[k];
            }

            /// Set up 1092-1098 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzzzz_xx = cbuffer.data(hd_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_z_xzzzz_xy = cbuffer.data(hd_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_z_xzzzz_xz = cbuffer.data(hd_geom_11_off + 1094 * ccomps * dcomps);

            auto g_z_z_xzzzz_yy = cbuffer.data(hd_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_z_xzzzz_yz = cbuffer.data(hd_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_z_xzzzz_zz = cbuffer.data(hd_geom_11_off + 1097 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzzzz_xx, g_z_z_xzzzz_xy, g_z_z_xzzzz_xz, g_z_z_xzzzz_yy, g_z_z_xzzzz_yz, g_z_z_xzzzz_zz, g_z_z_zzzz_xx, g_z_z_zzzz_xxx, g_z_z_zzzz_xxy, g_z_z_zzzz_xxz, g_z_z_zzzz_xy, g_z_z_zzzz_xyy, g_z_z_zzzz_xyz, g_z_z_zzzz_xz, g_z_z_zzzz_xzz, g_z_z_zzzz_yy, g_z_z_zzzz_yz, g_z_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzzzz_xx[k] = -g_z_z_zzzz_xx[k] * ab_x + g_z_z_zzzz_xxx[k];

                g_z_z_xzzzz_xy[k] = -g_z_z_zzzz_xy[k] * ab_x + g_z_z_zzzz_xxy[k];

                g_z_z_xzzzz_xz[k] = -g_z_z_zzzz_xz[k] * ab_x + g_z_z_zzzz_xxz[k];

                g_z_z_xzzzz_yy[k] = -g_z_z_zzzz_yy[k] * ab_x + g_z_z_zzzz_xyy[k];

                g_z_z_xzzzz_yz[k] = -g_z_z_zzzz_yz[k] * ab_x + g_z_z_zzzz_xyz[k];

                g_z_z_xzzzz_zz[k] = -g_z_z_zzzz_zz[k] * ab_x + g_z_z_zzzz_xzz[k];
            }

            /// Set up 1098-1104 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyyy_xx = cbuffer.data(hd_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_z_yyyyy_xy = cbuffer.data(hd_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_z_yyyyy_xz = cbuffer.data(hd_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_z_yyyyy_yy = cbuffer.data(hd_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_z_yyyyy_yz = cbuffer.data(hd_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_z_yyyyy_zz = cbuffer.data(hd_geom_11_off + 1103 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyy_xx, g_z_z_yyyy_xxy, g_z_z_yyyy_xy, g_z_z_yyyy_xyy, g_z_z_yyyy_xyz, g_z_z_yyyy_xz, g_z_z_yyyy_yy, g_z_z_yyyy_yyy, g_z_z_yyyy_yyz, g_z_z_yyyy_yz, g_z_z_yyyy_yzz, g_z_z_yyyy_zz, g_z_z_yyyyy_xx, g_z_z_yyyyy_xy, g_z_z_yyyyy_xz, g_z_z_yyyyy_yy, g_z_z_yyyyy_yz, g_z_z_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyyy_xx[k] = -g_z_z_yyyy_xx[k] * ab_y + g_z_z_yyyy_xxy[k];

                g_z_z_yyyyy_xy[k] = -g_z_z_yyyy_xy[k] * ab_y + g_z_z_yyyy_xyy[k];

                g_z_z_yyyyy_xz[k] = -g_z_z_yyyy_xz[k] * ab_y + g_z_z_yyyy_xyz[k];

                g_z_z_yyyyy_yy[k] = -g_z_z_yyyy_yy[k] * ab_y + g_z_z_yyyy_yyy[k];

                g_z_z_yyyyy_yz[k] = -g_z_z_yyyy_yz[k] * ab_y + g_z_z_yyyy_yyz[k];

                g_z_z_yyyyy_zz[k] = -g_z_z_yyyy_zz[k] * ab_y + g_z_z_yyyy_yzz[k];
            }

            /// Set up 1104-1110 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyyz_xx = cbuffer.data(hd_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_z_yyyyz_xy = cbuffer.data(hd_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_z_yyyyz_xz = cbuffer.data(hd_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_z_yyyyz_yy = cbuffer.data(hd_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_z_yyyyz_yz = cbuffer.data(hd_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_z_yyyyz_zz = cbuffer.data(hd_geom_11_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyyz_xx, g_z_z_yyyyz_xy, g_z_z_yyyyz_xz, g_z_z_yyyyz_yy, g_z_z_yyyyz_yz, g_z_z_yyyyz_zz, g_z_z_yyyz_xx, g_z_z_yyyz_xxy, g_z_z_yyyz_xy, g_z_z_yyyz_xyy, g_z_z_yyyz_xyz, g_z_z_yyyz_xz, g_z_z_yyyz_yy, g_z_z_yyyz_yyy, g_z_z_yyyz_yyz, g_z_z_yyyz_yz, g_z_z_yyyz_yzz, g_z_z_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyyz_xx[k] = -g_z_z_yyyz_xx[k] * ab_y + g_z_z_yyyz_xxy[k];

                g_z_z_yyyyz_xy[k] = -g_z_z_yyyz_xy[k] * ab_y + g_z_z_yyyz_xyy[k];

                g_z_z_yyyyz_xz[k] = -g_z_z_yyyz_xz[k] * ab_y + g_z_z_yyyz_xyz[k];

                g_z_z_yyyyz_yy[k] = -g_z_z_yyyz_yy[k] * ab_y + g_z_z_yyyz_yyy[k];

                g_z_z_yyyyz_yz[k] = -g_z_z_yyyz_yz[k] * ab_y + g_z_z_yyyz_yyz[k];

                g_z_z_yyyyz_zz[k] = -g_z_z_yyyz_zz[k] * ab_y + g_z_z_yyyz_yzz[k];
            }

            /// Set up 1110-1116 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyzz_xx = cbuffer.data(hd_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_z_yyyzz_xy = cbuffer.data(hd_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_z_yyyzz_xz = cbuffer.data(hd_geom_11_off + 1112 * ccomps * dcomps);

            auto g_z_z_yyyzz_yy = cbuffer.data(hd_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_z_yyyzz_yz = cbuffer.data(hd_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_z_yyyzz_zz = cbuffer.data(hd_geom_11_off + 1115 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyzz_xx, g_z_z_yyyzz_xy, g_z_z_yyyzz_xz, g_z_z_yyyzz_yy, g_z_z_yyyzz_yz, g_z_z_yyyzz_zz, g_z_z_yyzz_xx, g_z_z_yyzz_xxy, g_z_z_yyzz_xy, g_z_z_yyzz_xyy, g_z_z_yyzz_xyz, g_z_z_yyzz_xz, g_z_z_yyzz_yy, g_z_z_yyzz_yyy, g_z_z_yyzz_yyz, g_z_z_yyzz_yz, g_z_z_yyzz_yzz, g_z_z_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyzz_xx[k] = -g_z_z_yyzz_xx[k] * ab_y + g_z_z_yyzz_xxy[k];

                g_z_z_yyyzz_xy[k] = -g_z_z_yyzz_xy[k] * ab_y + g_z_z_yyzz_xyy[k];

                g_z_z_yyyzz_xz[k] = -g_z_z_yyzz_xz[k] * ab_y + g_z_z_yyzz_xyz[k];

                g_z_z_yyyzz_yy[k] = -g_z_z_yyzz_yy[k] * ab_y + g_z_z_yyzz_yyy[k];

                g_z_z_yyyzz_yz[k] = -g_z_z_yyzz_yz[k] * ab_y + g_z_z_yyzz_yyz[k];

                g_z_z_yyyzz_zz[k] = -g_z_z_yyzz_zz[k] * ab_y + g_z_z_yyzz_yzz[k];
            }

            /// Set up 1116-1122 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyzzz_xx = cbuffer.data(hd_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_z_yyzzz_xy = cbuffer.data(hd_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_z_yyzzz_xz = cbuffer.data(hd_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_z_yyzzz_yy = cbuffer.data(hd_geom_11_off + 1119 * ccomps * dcomps);

            auto g_z_z_yyzzz_yz = cbuffer.data(hd_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_z_yyzzz_zz = cbuffer.data(hd_geom_11_off + 1121 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyzzz_xx, g_z_z_yyzzz_xy, g_z_z_yyzzz_xz, g_z_z_yyzzz_yy, g_z_z_yyzzz_yz, g_z_z_yyzzz_zz, g_z_z_yzzz_xx, g_z_z_yzzz_xxy, g_z_z_yzzz_xy, g_z_z_yzzz_xyy, g_z_z_yzzz_xyz, g_z_z_yzzz_xz, g_z_z_yzzz_yy, g_z_z_yzzz_yyy, g_z_z_yzzz_yyz, g_z_z_yzzz_yz, g_z_z_yzzz_yzz, g_z_z_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyzzz_xx[k] = -g_z_z_yzzz_xx[k] * ab_y + g_z_z_yzzz_xxy[k];

                g_z_z_yyzzz_xy[k] = -g_z_z_yzzz_xy[k] * ab_y + g_z_z_yzzz_xyy[k];

                g_z_z_yyzzz_xz[k] = -g_z_z_yzzz_xz[k] * ab_y + g_z_z_yzzz_xyz[k];

                g_z_z_yyzzz_yy[k] = -g_z_z_yzzz_yy[k] * ab_y + g_z_z_yzzz_yyy[k];

                g_z_z_yyzzz_yz[k] = -g_z_z_yzzz_yz[k] * ab_y + g_z_z_yzzz_yyz[k];

                g_z_z_yyzzz_zz[k] = -g_z_z_yzzz_zz[k] * ab_y + g_z_z_yzzz_yzz[k];
            }

            /// Set up 1122-1128 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzzzz_xx = cbuffer.data(hd_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_z_yzzzz_xy = cbuffer.data(hd_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_z_yzzzz_xz = cbuffer.data(hd_geom_11_off + 1124 * ccomps * dcomps);

            auto g_z_z_yzzzz_yy = cbuffer.data(hd_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_z_yzzzz_yz = cbuffer.data(hd_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_z_yzzzz_zz = cbuffer.data(hd_geom_11_off + 1127 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzzzz_xx, g_z_z_yzzzz_xy, g_z_z_yzzzz_xz, g_z_z_yzzzz_yy, g_z_z_yzzzz_yz, g_z_z_yzzzz_zz, g_z_z_zzzz_xx, g_z_z_zzzz_xxy, g_z_z_zzzz_xy, g_z_z_zzzz_xyy, g_z_z_zzzz_xyz, g_z_z_zzzz_xz, g_z_z_zzzz_yy, g_z_z_zzzz_yyy, g_z_z_zzzz_yyz, g_z_z_zzzz_yz, g_z_z_zzzz_yzz, g_z_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzzzz_xx[k] = -g_z_z_zzzz_xx[k] * ab_y + g_z_z_zzzz_xxy[k];

                g_z_z_yzzzz_xy[k] = -g_z_z_zzzz_xy[k] * ab_y + g_z_z_zzzz_xyy[k];

                g_z_z_yzzzz_xz[k] = -g_z_z_zzzz_xz[k] * ab_y + g_z_z_zzzz_xyz[k];

                g_z_z_yzzzz_yy[k] = -g_z_z_zzzz_yy[k] * ab_y + g_z_z_zzzz_yyy[k];

                g_z_z_yzzzz_yz[k] = -g_z_z_zzzz_yz[k] * ab_y + g_z_z_zzzz_yyz[k];

                g_z_z_yzzzz_zz[k] = -g_z_z_zzzz_zz[k] * ab_y + g_z_z_zzzz_yzz[k];
            }

            /// Set up 1128-1134 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzzzz_xx = cbuffer.data(hd_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_z_zzzzz_xy = cbuffer.data(hd_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_z_zzzzz_xz = cbuffer.data(hd_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_z_zzzzz_yy = cbuffer.data(hd_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_z_zzzzz_yz = cbuffer.data(hd_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_z_zzzzz_zz = cbuffer.data(hd_geom_11_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xx, g_0_z_zzzz_xy, g_0_z_zzzz_xz, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_zz, g_z_0_zzzz_xx, g_z_0_zzzz_xy, g_z_0_zzzz_xz, g_z_0_zzzz_yy, g_z_0_zzzz_yz, g_z_0_zzzz_zz, g_z_z_zzzz_xx, g_z_z_zzzz_xxz, g_z_z_zzzz_xy, g_z_z_zzzz_xyz, g_z_z_zzzz_xz, g_z_z_zzzz_xzz, g_z_z_zzzz_yy, g_z_z_zzzz_yyz, g_z_z_zzzz_yz, g_z_z_zzzz_yzz, g_z_z_zzzz_zz, g_z_z_zzzz_zzz, g_z_z_zzzzz_xx, g_z_z_zzzzz_xy, g_z_z_zzzzz_xz, g_z_z_zzzzz_yy, g_z_z_zzzzz_yz, g_z_z_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzzzz_xx[k] = -g_0_z_zzzz_xx[k] + g_z_0_zzzz_xx[k] - g_z_z_zzzz_xx[k] * ab_z + g_z_z_zzzz_xxz[k];

                g_z_z_zzzzz_xy[k] = -g_0_z_zzzz_xy[k] + g_z_0_zzzz_xy[k] - g_z_z_zzzz_xy[k] * ab_z + g_z_z_zzzz_xyz[k];

                g_z_z_zzzzz_xz[k] = -g_0_z_zzzz_xz[k] + g_z_0_zzzz_xz[k] - g_z_z_zzzz_xz[k] * ab_z + g_z_z_zzzz_xzz[k];

                g_z_z_zzzzz_yy[k] = -g_0_z_zzzz_yy[k] + g_z_0_zzzz_yy[k] - g_z_z_zzzz_yy[k] * ab_z + g_z_z_zzzz_yyz[k];

                g_z_z_zzzzz_yz[k] = -g_0_z_zzzz_yz[k] + g_z_0_zzzz_yz[k] - g_z_z_zzzz_yz[k] * ab_z + g_z_z_zzzz_yzz[k];

                g_z_z_zzzzz_zz[k] = -g_0_z_zzzz_zz[k] + g_z_0_zzzz_zz[k] - g_z_z_zzzz_zz[k] * ab_z + g_z_z_zzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

