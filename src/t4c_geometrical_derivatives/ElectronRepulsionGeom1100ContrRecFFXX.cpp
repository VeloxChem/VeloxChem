#include "ElectronRepulsionGeom1100ContrRecFFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_ffxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_ffxx,
                                            const size_t idx_geom_01_dfxx,
                                            const size_t idx_geom_10_dfxx,
                                            const size_t idx_geom_11_dfxx,
                                            const size_t idx_geom_11_dgxx,
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
            /// Set up components of auxilary buffer : DFSS

            const auto df_geom_01_off = idx_geom_01_dfxx + i * dcomps + j;

            auto g_0_x_xx_xxx = cbuffer.data(df_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxy = cbuffer.data(df_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxz = cbuffer.data(df_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xyy = cbuffer.data(df_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xyz = cbuffer.data(df_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xzz = cbuffer.data(df_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_yyy = cbuffer.data(df_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_yyz = cbuffer.data(df_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_yzz = cbuffer.data(df_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_zzz = cbuffer.data(df_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xy_xxx = cbuffer.data(df_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xy_xxy = cbuffer.data(df_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xy_xxz = cbuffer.data(df_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xy_xyy = cbuffer.data(df_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xy_xyz = cbuffer.data(df_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xy_xzz = cbuffer.data(df_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xy_yyy = cbuffer.data(df_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xy_yyz = cbuffer.data(df_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xy_yzz = cbuffer.data(df_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xy_zzz = cbuffer.data(df_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xz_xxx = cbuffer.data(df_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xz_xxy = cbuffer.data(df_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xz_xxz = cbuffer.data(df_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xz_xyy = cbuffer.data(df_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xz_xyz = cbuffer.data(df_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xz_xzz = cbuffer.data(df_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xz_yyy = cbuffer.data(df_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xz_yyz = cbuffer.data(df_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xz_yzz = cbuffer.data(df_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xz_zzz = cbuffer.data(df_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_yy_xxx = cbuffer.data(df_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_yy_xxy = cbuffer.data(df_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_yy_xxz = cbuffer.data(df_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_yy_xyy = cbuffer.data(df_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_yy_xyz = cbuffer.data(df_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_yy_xzz = cbuffer.data(df_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_yy_yyy = cbuffer.data(df_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_yy_yyz = cbuffer.data(df_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_yy_yzz = cbuffer.data(df_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_yy_zzz = cbuffer.data(df_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_yz_xxx = cbuffer.data(df_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_yz_xxy = cbuffer.data(df_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_yz_xxz = cbuffer.data(df_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_yz_xyy = cbuffer.data(df_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_yz_xyz = cbuffer.data(df_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yz_xzz = cbuffer.data(df_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yz_yyy = cbuffer.data(df_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yz_yyz = cbuffer.data(df_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yz_yzz = cbuffer.data(df_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yz_zzz = cbuffer.data(df_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_zz_xxx = cbuffer.data(df_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_zz_xxy = cbuffer.data(df_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_zz_xxz = cbuffer.data(df_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_zz_xyy = cbuffer.data(df_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_zz_xyz = cbuffer.data(df_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_zz_xzz = cbuffer.data(df_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_zz_yyy = cbuffer.data(df_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_zz_yyz = cbuffer.data(df_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_zz_yzz = cbuffer.data(df_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_zz_zzz = cbuffer.data(df_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_xx_xxx = cbuffer.data(df_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_xx_xxy = cbuffer.data(df_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_xx_xxz = cbuffer.data(df_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_xx_xyy = cbuffer.data(df_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xx_xyz = cbuffer.data(df_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xx_xzz = cbuffer.data(df_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_xx_yyy = cbuffer.data(df_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xx_yyz = cbuffer.data(df_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xx_yzz = cbuffer.data(df_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_xx_zzz = cbuffer.data(df_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xy_xxx = cbuffer.data(df_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xy_xxy = cbuffer.data(df_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_xy_xxz = cbuffer.data(df_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xy_xyy = cbuffer.data(df_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xy_xyz = cbuffer.data(df_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_xy_xzz = cbuffer.data(df_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_xy_yyy = cbuffer.data(df_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_xy_yyz = cbuffer.data(df_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_xy_yzz = cbuffer.data(df_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_xy_zzz = cbuffer.data(df_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_xz_xxx = cbuffer.data(df_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_xz_xxy = cbuffer.data(df_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_xz_xxz = cbuffer.data(df_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_xz_xyy = cbuffer.data(df_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_xz_xyz = cbuffer.data(df_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xz_xzz = cbuffer.data(df_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xz_yyy = cbuffer.data(df_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_xz_yyz = cbuffer.data(df_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xz_yzz = cbuffer.data(df_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xz_zzz = cbuffer.data(df_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_yy_xxx = cbuffer.data(df_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_yy_xxy = cbuffer.data(df_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_yy_xxz = cbuffer.data(df_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_yy_xyy = cbuffer.data(df_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_yy_xyz = cbuffer.data(df_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_yy_xzz = cbuffer.data(df_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_yy_yyy = cbuffer.data(df_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_yy_yyz = cbuffer.data(df_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_yy_yzz = cbuffer.data(df_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_yy_zzz = cbuffer.data(df_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_yz_xxx = cbuffer.data(df_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_yz_xxy = cbuffer.data(df_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_yz_xxz = cbuffer.data(df_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_yz_xyy = cbuffer.data(df_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_yz_xyz = cbuffer.data(df_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_yz_xzz = cbuffer.data(df_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_yz_yyy = cbuffer.data(df_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_yz_yyz = cbuffer.data(df_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_yz_yzz = cbuffer.data(df_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_yz_zzz = cbuffer.data(df_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_zz_xxx = cbuffer.data(df_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_zz_xxy = cbuffer.data(df_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_zz_xxz = cbuffer.data(df_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_zz_xyy = cbuffer.data(df_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_zz_xyz = cbuffer.data(df_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_zz_xzz = cbuffer.data(df_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_zz_yyy = cbuffer.data(df_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_zz_yyz = cbuffer.data(df_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_zz_yzz = cbuffer.data(df_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_zz_zzz = cbuffer.data(df_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_z_xx_xxx = cbuffer.data(df_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_xx_xxy = cbuffer.data(df_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_xx_xxz = cbuffer.data(df_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_z_xx_xyy = cbuffer.data(df_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_xx_xyz = cbuffer.data(df_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_xx_xzz = cbuffer.data(df_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_xx_yyy = cbuffer.data(df_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_xx_yyz = cbuffer.data(df_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_xx_yzz = cbuffer.data(df_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_xx_zzz = cbuffer.data(df_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_xy_xxx = cbuffer.data(df_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_xy_xxy = cbuffer.data(df_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_xy_xxz = cbuffer.data(df_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_xy_xyy = cbuffer.data(df_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_xy_xyz = cbuffer.data(df_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_z_xy_xzz = cbuffer.data(df_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_xy_yyy = cbuffer.data(df_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_xy_yyz = cbuffer.data(df_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_z_xy_yzz = cbuffer.data(df_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_xy_zzz = cbuffer.data(df_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_xz_xxx = cbuffer.data(df_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_z_xz_xxy = cbuffer.data(df_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_xz_xxz = cbuffer.data(df_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_xz_xyy = cbuffer.data(df_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_z_xz_xyz = cbuffer.data(df_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_xz_xzz = cbuffer.data(df_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_xz_yyy = cbuffer.data(df_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_z_xz_yyz = cbuffer.data(df_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_xz_yzz = cbuffer.data(df_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_xz_zzz = cbuffer.data(df_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_z_yy_xxx = cbuffer.data(df_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_yy_xxy = cbuffer.data(df_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_yy_xxz = cbuffer.data(df_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_z_yy_xyy = cbuffer.data(df_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_yy_xyz = cbuffer.data(df_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_yy_xzz = cbuffer.data(df_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_z_yy_yyy = cbuffer.data(df_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_yy_yyz = cbuffer.data(df_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_yy_yzz = cbuffer.data(df_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_z_yy_zzz = cbuffer.data(df_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_yz_xxx = cbuffer.data(df_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_yz_xxy = cbuffer.data(df_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_z_yz_xxz = cbuffer.data(df_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_yz_xyy = cbuffer.data(df_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_yz_xyz = cbuffer.data(df_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_z_yz_xzz = cbuffer.data(df_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_yz_yyy = cbuffer.data(df_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_yz_yyz = cbuffer.data(df_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_yz_yzz = cbuffer.data(df_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_yz_zzz = cbuffer.data(df_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_zz_xxx = cbuffer.data(df_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_zz_xxy = cbuffer.data(df_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_zz_xxz = cbuffer.data(df_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_zz_xyy = cbuffer.data(df_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_zz_xyz = cbuffer.data(df_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_zz_xzz = cbuffer.data(df_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_zz_yyy = cbuffer.data(df_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_zz_yyz = cbuffer.data(df_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_zz_yzz = cbuffer.data(df_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_zz_zzz = cbuffer.data(df_geom_01_off + 179 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DFSS

            const auto df_geom_10_off = idx_geom_10_dfxx + i * dcomps + j;

            auto g_x_0_xx_xxx = cbuffer.data(df_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxy = cbuffer.data(df_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxz = cbuffer.data(df_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xyy = cbuffer.data(df_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xyz = cbuffer.data(df_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xzz = cbuffer.data(df_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_yyy = cbuffer.data(df_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_yyz = cbuffer.data(df_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_yzz = cbuffer.data(df_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_zzz = cbuffer.data(df_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xy_xxx = cbuffer.data(df_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xy_xxy = cbuffer.data(df_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xy_xxz = cbuffer.data(df_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xy_xyy = cbuffer.data(df_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xy_xyz = cbuffer.data(df_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xy_xzz = cbuffer.data(df_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xy_yyy = cbuffer.data(df_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xy_yyz = cbuffer.data(df_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xy_yzz = cbuffer.data(df_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xy_zzz = cbuffer.data(df_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xz_xxx = cbuffer.data(df_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xz_xxy = cbuffer.data(df_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xz_xxz = cbuffer.data(df_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xz_xyy = cbuffer.data(df_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xz_xyz = cbuffer.data(df_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xz_xzz = cbuffer.data(df_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xz_yyy = cbuffer.data(df_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xz_yyz = cbuffer.data(df_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xz_yzz = cbuffer.data(df_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xz_zzz = cbuffer.data(df_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_yy_xxx = cbuffer.data(df_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_yy_xxy = cbuffer.data(df_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_yy_xxz = cbuffer.data(df_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_yy_xyy = cbuffer.data(df_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_yy_xyz = cbuffer.data(df_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_yy_xzz = cbuffer.data(df_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_yy_yyy = cbuffer.data(df_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yy_yyz = cbuffer.data(df_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yy_yzz = cbuffer.data(df_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_yy_zzz = cbuffer.data(df_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yz_xxx = cbuffer.data(df_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yz_xxy = cbuffer.data(df_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_yz_xxz = cbuffer.data(df_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_yz_xyy = cbuffer.data(df_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_yz_xyz = cbuffer.data(df_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yz_xzz = cbuffer.data(df_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yz_yyy = cbuffer.data(df_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yz_yyz = cbuffer.data(df_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yz_yzz = cbuffer.data(df_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yz_zzz = cbuffer.data(df_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_zz_xxx = cbuffer.data(df_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_zz_xxy = cbuffer.data(df_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_zz_xxz = cbuffer.data(df_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_zz_xyy = cbuffer.data(df_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_zz_xyz = cbuffer.data(df_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_zz_xzz = cbuffer.data(df_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_zz_yyy = cbuffer.data(df_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_zz_yyz = cbuffer.data(df_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_zz_yzz = cbuffer.data(df_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_zz_zzz = cbuffer.data(df_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_xx_xxx = cbuffer.data(df_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xx_xxy = cbuffer.data(df_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xx_xxz = cbuffer.data(df_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xx_xyy = cbuffer.data(df_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xx_xyz = cbuffer.data(df_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xx_xzz = cbuffer.data(df_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xx_yyy = cbuffer.data(df_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xx_yyz = cbuffer.data(df_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xx_yzz = cbuffer.data(df_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xx_zzz = cbuffer.data(df_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xy_xxx = cbuffer.data(df_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xy_xxy = cbuffer.data(df_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xy_xxz = cbuffer.data(df_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xy_xyy = cbuffer.data(df_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xy_xyz = cbuffer.data(df_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xy_xzz = cbuffer.data(df_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xy_yyy = cbuffer.data(df_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xy_yyz = cbuffer.data(df_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_xy_yzz = cbuffer.data(df_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xy_zzz = cbuffer.data(df_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xz_xxx = cbuffer.data(df_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xz_xxy = cbuffer.data(df_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xz_xxz = cbuffer.data(df_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xz_xyy = cbuffer.data(df_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_xz_xyz = cbuffer.data(df_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xz_xzz = cbuffer.data(df_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xz_yyy = cbuffer.data(df_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xz_yyz = cbuffer.data(df_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xz_yzz = cbuffer.data(df_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xz_zzz = cbuffer.data(df_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_yy_xxx = cbuffer.data(df_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_yy_xxy = cbuffer.data(df_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_yy_xxz = cbuffer.data(df_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_yy_xyy = cbuffer.data(df_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_yy_xyz = cbuffer.data(df_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_yy_xzz = cbuffer.data(df_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_yy_yyy = cbuffer.data(df_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_yy_yyz = cbuffer.data(df_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_yy_yzz = cbuffer.data(df_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_yy_zzz = cbuffer.data(df_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_yz_xxx = cbuffer.data(df_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_yz_xxy = cbuffer.data(df_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_yz_xxz = cbuffer.data(df_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_yz_xyy = cbuffer.data(df_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_yz_xyz = cbuffer.data(df_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_yz_xzz = cbuffer.data(df_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_yz_yyy = cbuffer.data(df_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_yz_yyz = cbuffer.data(df_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yz_yzz = cbuffer.data(df_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yz_zzz = cbuffer.data(df_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_zz_xxx = cbuffer.data(df_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_zz_xxy = cbuffer.data(df_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_zz_xxz = cbuffer.data(df_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_zz_xyy = cbuffer.data(df_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_zz_xyz = cbuffer.data(df_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_zz_xzz = cbuffer.data(df_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_zz_yyy = cbuffer.data(df_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_zz_yyz = cbuffer.data(df_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_zz_yzz = cbuffer.data(df_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_zz_zzz = cbuffer.data(df_geom_10_off + 119 * ccomps * dcomps);

            auto g_z_0_xx_xxx = cbuffer.data(df_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_xx_xxy = cbuffer.data(df_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_xx_xxz = cbuffer.data(df_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_xx_xyy = cbuffer.data(df_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_xx_xyz = cbuffer.data(df_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_xx_xzz = cbuffer.data(df_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_xx_yyy = cbuffer.data(df_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xx_yyz = cbuffer.data(df_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xx_yzz = cbuffer.data(df_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xx_zzz = cbuffer.data(df_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_xy_xxx = cbuffer.data(df_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xy_xxy = cbuffer.data(df_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_xy_xxz = cbuffer.data(df_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xy_xyy = cbuffer.data(df_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xy_xyz = cbuffer.data(df_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xy_xzz = cbuffer.data(df_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xy_yyy = cbuffer.data(df_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xy_yyz = cbuffer.data(df_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_xy_yzz = cbuffer.data(df_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xy_zzz = cbuffer.data(df_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xz_xxx = cbuffer.data(df_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xz_xxy = cbuffer.data(df_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xz_xxz = cbuffer.data(df_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xz_xyy = cbuffer.data(df_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_xz_xyz = cbuffer.data(df_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xz_xzz = cbuffer.data(df_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xz_yyy = cbuffer.data(df_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xz_yyz = cbuffer.data(df_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xz_yzz = cbuffer.data(df_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xz_zzz = cbuffer.data(df_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_yy_xxx = cbuffer.data(df_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_yy_xxy = cbuffer.data(df_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_yy_xxz = cbuffer.data(df_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_yy_xyy = cbuffer.data(df_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_yy_xyz = cbuffer.data(df_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_yy_xzz = cbuffer.data(df_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_yy_yyy = cbuffer.data(df_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_yy_yyz = cbuffer.data(df_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_yy_yzz = cbuffer.data(df_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_yy_zzz = cbuffer.data(df_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_yz_xxx = cbuffer.data(df_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_yz_xxy = cbuffer.data(df_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_yz_xxz = cbuffer.data(df_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_yz_xyy = cbuffer.data(df_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_yz_xyz = cbuffer.data(df_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_yz_xzz = cbuffer.data(df_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_yz_yyy = cbuffer.data(df_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_yz_yyz = cbuffer.data(df_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_yz_yzz = cbuffer.data(df_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_yz_zzz = cbuffer.data(df_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_zz_xxx = cbuffer.data(df_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_zz_xxy = cbuffer.data(df_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_zz_xxz = cbuffer.data(df_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_zz_xyy = cbuffer.data(df_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_zz_xyz = cbuffer.data(df_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_zz_xzz = cbuffer.data(df_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_zz_yyy = cbuffer.data(df_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_zz_yyz = cbuffer.data(df_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_zz_yzz = cbuffer.data(df_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_zz_zzz = cbuffer.data(df_geom_10_off + 179 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DFSS

            const auto df_geom_11_off = idx_geom_11_dfxx + i * dcomps + j;

            auto g_x_x_xx_xxx = cbuffer.data(df_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxy = cbuffer.data(df_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxz = cbuffer.data(df_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xyy = cbuffer.data(df_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xyz = cbuffer.data(df_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xzz = cbuffer.data(df_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_yyy = cbuffer.data(df_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_yyz = cbuffer.data(df_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_yzz = cbuffer.data(df_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_zzz = cbuffer.data(df_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xy_xxx = cbuffer.data(df_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xy_xxy = cbuffer.data(df_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xy_xxz = cbuffer.data(df_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xy_xyy = cbuffer.data(df_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xy_xyz = cbuffer.data(df_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xy_xzz = cbuffer.data(df_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xy_yyy = cbuffer.data(df_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xy_yyz = cbuffer.data(df_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xy_yzz = cbuffer.data(df_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xy_zzz = cbuffer.data(df_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xz_xxx = cbuffer.data(df_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xz_xxy = cbuffer.data(df_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xz_xxz = cbuffer.data(df_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xz_xyy = cbuffer.data(df_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xz_xyz = cbuffer.data(df_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xz_xzz = cbuffer.data(df_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xz_yyy = cbuffer.data(df_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xz_yyz = cbuffer.data(df_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xz_yzz = cbuffer.data(df_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xz_zzz = cbuffer.data(df_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_yy_xxx = cbuffer.data(df_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_yy_xxy = cbuffer.data(df_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_yy_xxz = cbuffer.data(df_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_yy_xyy = cbuffer.data(df_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_yy_xyz = cbuffer.data(df_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_yy_xzz = cbuffer.data(df_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_yy_yyy = cbuffer.data(df_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_yy_yyz = cbuffer.data(df_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_yy_yzz = cbuffer.data(df_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_yy_zzz = cbuffer.data(df_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_yz_xxx = cbuffer.data(df_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_yz_xxy = cbuffer.data(df_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_yz_xxz = cbuffer.data(df_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_yz_xyy = cbuffer.data(df_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_yz_xyz = cbuffer.data(df_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_yz_xzz = cbuffer.data(df_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_yz_yyy = cbuffer.data(df_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_yz_yyz = cbuffer.data(df_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_yz_yzz = cbuffer.data(df_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_yz_zzz = cbuffer.data(df_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_zz_xxx = cbuffer.data(df_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_zz_xxy = cbuffer.data(df_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_zz_xxz = cbuffer.data(df_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_zz_xyy = cbuffer.data(df_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_zz_xyz = cbuffer.data(df_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_zz_xzz = cbuffer.data(df_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_zz_yyy = cbuffer.data(df_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_zz_yyz = cbuffer.data(df_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_zz_yzz = cbuffer.data(df_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_zz_zzz = cbuffer.data(df_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_xx_xxx = cbuffer.data(df_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_xx_xxy = cbuffer.data(df_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_xx_xxz = cbuffer.data(df_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_xx_xyy = cbuffer.data(df_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_xx_xyz = cbuffer.data(df_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_xx_xzz = cbuffer.data(df_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_xx_yyy = cbuffer.data(df_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_xx_yyz = cbuffer.data(df_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_xx_yzz = cbuffer.data(df_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_xx_zzz = cbuffer.data(df_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_xy_xxx = cbuffer.data(df_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_xy_xxy = cbuffer.data(df_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_xy_xxz = cbuffer.data(df_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_xy_xyy = cbuffer.data(df_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_xy_xyz = cbuffer.data(df_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_xy_xzz = cbuffer.data(df_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_xy_yyy = cbuffer.data(df_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_xy_yyz = cbuffer.data(df_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_xy_yzz = cbuffer.data(df_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_xy_zzz = cbuffer.data(df_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_xz_xxx = cbuffer.data(df_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_xz_xxy = cbuffer.data(df_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_xz_xxz = cbuffer.data(df_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_xz_xyy = cbuffer.data(df_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_xz_xyz = cbuffer.data(df_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_xz_xzz = cbuffer.data(df_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_xz_yyy = cbuffer.data(df_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_xz_yyz = cbuffer.data(df_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_xz_yzz = cbuffer.data(df_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_xz_zzz = cbuffer.data(df_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_yy_xxx = cbuffer.data(df_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_yy_xxy = cbuffer.data(df_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_yy_xxz = cbuffer.data(df_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_yy_xyy = cbuffer.data(df_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_yy_xyz = cbuffer.data(df_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_yy_xzz = cbuffer.data(df_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_yy_yyy = cbuffer.data(df_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_yy_yyz = cbuffer.data(df_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_yy_yzz = cbuffer.data(df_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_yy_zzz = cbuffer.data(df_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_yz_xxx = cbuffer.data(df_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_yz_xxy = cbuffer.data(df_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_yz_xxz = cbuffer.data(df_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_yz_xyy = cbuffer.data(df_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_yz_xyz = cbuffer.data(df_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_yz_xzz = cbuffer.data(df_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_yz_yyy = cbuffer.data(df_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_yz_yyz = cbuffer.data(df_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_yz_yzz = cbuffer.data(df_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_yz_zzz = cbuffer.data(df_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_zz_xxx = cbuffer.data(df_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_zz_xxy = cbuffer.data(df_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_zz_xxz = cbuffer.data(df_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_zz_xyy = cbuffer.data(df_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_zz_xyz = cbuffer.data(df_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_zz_xzz = cbuffer.data(df_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_zz_yyy = cbuffer.data(df_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_zz_yyz = cbuffer.data(df_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_zz_yzz = cbuffer.data(df_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_zz_zzz = cbuffer.data(df_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_z_xx_xxx = cbuffer.data(df_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_xx_xxy = cbuffer.data(df_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_xx_xxz = cbuffer.data(df_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_xx_xyy = cbuffer.data(df_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_xx_xyz = cbuffer.data(df_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_xx_xzz = cbuffer.data(df_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_xx_yyy = cbuffer.data(df_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_xx_yyz = cbuffer.data(df_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_xx_yzz = cbuffer.data(df_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_xx_zzz = cbuffer.data(df_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_xy_xxx = cbuffer.data(df_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_xy_xxy = cbuffer.data(df_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_xy_xxz = cbuffer.data(df_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_xy_xyy = cbuffer.data(df_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_xy_xyz = cbuffer.data(df_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_xy_xzz = cbuffer.data(df_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_xy_yyy = cbuffer.data(df_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_xy_yyz = cbuffer.data(df_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_xy_yzz = cbuffer.data(df_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_xy_zzz = cbuffer.data(df_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_z_xz_xxx = cbuffer.data(df_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_xz_xxy = cbuffer.data(df_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_xz_xxz = cbuffer.data(df_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_xz_xyy = cbuffer.data(df_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_xz_xyz = cbuffer.data(df_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_xz_xzz = cbuffer.data(df_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_xz_yyy = cbuffer.data(df_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_z_xz_yyz = cbuffer.data(df_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_xz_yzz = cbuffer.data(df_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_xz_zzz = cbuffer.data(df_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_z_yy_xxx = cbuffer.data(df_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_yy_xxy = cbuffer.data(df_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_yy_xxz = cbuffer.data(df_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_yy_xyy = cbuffer.data(df_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_yy_xyz = cbuffer.data(df_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_yy_xzz = cbuffer.data(df_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_yy_yyy = cbuffer.data(df_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_yy_yyz = cbuffer.data(df_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_yy_yzz = cbuffer.data(df_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_yy_zzz = cbuffer.data(df_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_z_yz_xxx = cbuffer.data(df_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_yz_xxy = cbuffer.data(df_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_yz_xxz = cbuffer.data(df_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_yz_xyy = cbuffer.data(df_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_yz_xyz = cbuffer.data(df_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_yz_xzz = cbuffer.data(df_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_yz_yyy = cbuffer.data(df_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_yz_yyz = cbuffer.data(df_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_z_yz_yzz = cbuffer.data(df_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_yz_zzz = cbuffer.data(df_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_zz_xxx = cbuffer.data(df_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_zz_xxy = cbuffer.data(df_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_zz_xxz = cbuffer.data(df_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_zz_xyy = cbuffer.data(df_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_zz_xyz = cbuffer.data(df_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_zz_xzz = cbuffer.data(df_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_zz_yyy = cbuffer.data(df_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_zz_yyz = cbuffer.data(df_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_zz_yzz = cbuffer.data(df_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_zz_zzz = cbuffer.data(df_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_x_xx_xxx = cbuffer.data(df_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_x_xx_xxy = cbuffer.data(df_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_x_xx_xxz = cbuffer.data(df_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_x_xx_xyy = cbuffer.data(df_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_x_xx_xyz = cbuffer.data(df_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_x_xx_xzz = cbuffer.data(df_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_x_xx_yyy = cbuffer.data(df_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_x_xx_yyz = cbuffer.data(df_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_x_xx_yzz = cbuffer.data(df_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_x_xx_zzz = cbuffer.data(df_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_x_xy_xxx = cbuffer.data(df_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_xy_xxy = cbuffer.data(df_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_xy_xxz = cbuffer.data(df_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_xy_xyy = cbuffer.data(df_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_xy_xyz = cbuffer.data(df_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_xy_xzz = cbuffer.data(df_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_xy_yyy = cbuffer.data(df_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_xy_yyz = cbuffer.data(df_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_xy_yzz = cbuffer.data(df_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_xy_zzz = cbuffer.data(df_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_x_xz_xxx = cbuffer.data(df_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_xz_xxy = cbuffer.data(df_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_xz_xxz = cbuffer.data(df_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_xz_xyy = cbuffer.data(df_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_xz_xyz = cbuffer.data(df_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_xz_xzz = cbuffer.data(df_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_xz_yyy = cbuffer.data(df_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_xz_yyz = cbuffer.data(df_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_xz_yzz = cbuffer.data(df_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_xz_zzz = cbuffer.data(df_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_x_yy_xxx = cbuffer.data(df_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_yy_xxy = cbuffer.data(df_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_yy_xxz = cbuffer.data(df_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_yy_xyy = cbuffer.data(df_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_yy_xyz = cbuffer.data(df_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_yy_xzz = cbuffer.data(df_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_yy_yyy = cbuffer.data(df_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_yy_yyz = cbuffer.data(df_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_yy_yzz = cbuffer.data(df_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_yy_zzz = cbuffer.data(df_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_x_yz_xxx = cbuffer.data(df_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_yz_xxy = cbuffer.data(df_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_yz_xxz = cbuffer.data(df_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_yz_xyy = cbuffer.data(df_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_yz_xyz = cbuffer.data(df_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_yz_xzz = cbuffer.data(df_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_yz_yyy = cbuffer.data(df_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_yz_yyz = cbuffer.data(df_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_yz_yzz = cbuffer.data(df_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_yz_zzz = cbuffer.data(df_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_x_zz_xxx = cbuffer.data(df_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_x_zz_xxy = cbuffer.data(df_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_zz_xxz = cbuffer.data(df_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_zz_xyy = cbuffer.data(df_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_zz_xyz = cbuffer.data(df_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_zz_xzz = cbuffer.data(df_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_zz_yyy = cbuffer.data(df_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_zz_yyz = cbuffer.data(df_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_zz_yzz = cbuffer.data(df_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_zz_zzz = cbuffer.data(df_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_y_xx_xxx = cbuffer.data(df_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_y_xx_xxy = cbuffer.data(df_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_y_xx_xxz = cbuffer.data(df_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_y_xx_xyy = cbuffer.data(df_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_y_xx_xyz = cbuffer.data(df_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_y_xx_xzz = cbuffer.data(df_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_y_xx_yyy = cbuffer.data(df_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_y_xx_yyz = cbuffer.data(df_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_y_xx_yzz = cbuffer.data(df_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_y_xx_zzz = cbuffer.data(df_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_y_xy_xxx = cbuffer.data(df_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_y_xy_xxy = cbuffer.data(df_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_y_xy_xxz = cbuffer.data(df_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_xy_xyy = cbuffer.data(df_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_xy_xyz = cbuffer.data(df_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_xy_xzz = cbuffer.data(df_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_xy_yyy = cbuffer.data(df_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_xy_yyz = cbuffer.data(df_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_xy_yzz = cbuffer.data(df_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_xy_zzz = cbuffer.data(df_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_y_xz_xxx = cbuffer.data(df_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_xz_xxy = cbuffer.data(df_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_xz_xxz = cbuffer.data(df_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_xz_xyy = cbuffer.data(df_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_xz_xyz = cbuffer.data(df_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_xz_xzz = cbuffer.data(df_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_xz_yyy = cbuffer.data(df_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_xz_yyz = cbuffer.data(df_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_xz_yzz = cbuffer.data(df_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_xz_zzz = cbuffer.data(df_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_y_yy_xxx = cbuffer.data(df_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_yy_xxy = cbuffer.data(df_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_yy_xxz = cbuffer.data(df_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_y_yy_xyy = cbuffer.data(df_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_yy_xyz = cbuffer.data(df_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_yy_xzz = cbuffer.data(df_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_yy_yyy = cbuffer.data(df_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_yy_yyz = cbuffer.data(df_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_yy_yzz = cbuffer.data(df_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_yy_zzz = cbuffer.data(df_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_y_yz_xxx = cbuffer.data(df_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_yz_xxy = cbuffer.data(df_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_yz_xxz = cbuffer.data(df_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_yz_xyy = cbuffer.data(df_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_yz_xyz = cbuffer.data(df_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_yz_xzz = cbuffer.data(df_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_yz_yyy = cbuffer.data(df_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_yz_yyz = cbuffer.data(df_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_yz_yzz = cbuffer.data(df_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_yz_zzz = cbuffer.data(df_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_y_zz_xxx = cbuffer.data(df_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_zz_xxy = cbuffer.data(df_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_zz_xxz = cbuffer.data(df_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_zz_xyy = cbuffer.data(df_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_y_zz_xyz = cbuffer.data(df_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_zz_xzz = cbuffer.data(df_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_zz_yyy = cbuffer.data(df_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_zz_yyz = cbuffer.data(df_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_zz_yzz = cbuffer.data(df_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_zz_zzz = cbuffer.data(df_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_z_xx_xxx = cbuffer.data(df_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_z_xx_xxy = cbuffer.data(df_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_z_xx_xxz = cbuffer.data(df_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_z_xx_xyy = cbuffer.data(df_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_z_xx_xyz = cbuffer.data(df_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_z_xx_xzz = cbuffer.data(df_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_z_xx_yyy = cbuffer.data(df_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_z_xx_yyz = cbuffer.data(df_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_z_xx_yzz = cbuffer.data(df_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_z_xx_zzz = cbuffer.data(df_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_z_xy_xxx = cbuffer.data(df_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_z_xy_xxy = cbuffer.data(df_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_z_xy_xxz = cbuffer.data(df_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_z_xy_xyy = cbuffer.data(df_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_z_xy_xyz = cbuffer.data(df_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_z_xy_xzz = cbuffer.data(df_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_xy_yyy = cbuffer.data(df_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_xy_yyz = cbuffer.data(df_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_xy_yzz = cbuffer.data(df_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_xy_zzz = cbuffer.data(df_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_z_xz_xxx = cbuffer.data(df_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_xz_xxy = cbuffer.data(df_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_xz_xxz = cbuffer.data(df_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_xz_xyy = cbuffer.data(df_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_xz_xyz = cbuffer.data(df_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_xz_xzz = cbuffer.data(df_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_xz_yyy = cbuffer.data(df_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_xz_yyz = cbuffer.data(df_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_xz_yzz = cbuffer.data(df_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_xz_zzz = cbuffer.data(df_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_z_yy_xxx = cbuffer.data(df_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_yy_xxy = cbuffer.data(df_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_yy_xxz = cbuffer.data(df_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_yy_xyy = cbuffer.data(df_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_yy_xyz = cbuffer.data(df_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_yy_xzz = cbuffer.data(df_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_z_yy_yyy = cbuffer.data(df_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_yy_yyz = cbuffer.data(df_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_yy_yzz = cbuffer.data(df_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_yy_zzz = cbuffer.data(df_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_z_yz_xxx = cbuffer.data(df_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_yz_xxy = cbuffer.data(df_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_yz_xxz = cbuffer.data(df_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_yz_xyy = cbuffer.data(df_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_yz_xyz = cbuffer.data(df_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_yz_xzz = cbuffer.data(df_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_yz_yyy = cbuffer.data(df_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_yz_yyz = cbuffer.data(df_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_yz_yzz = cbuffer.data(df_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_yz_zzz = cbuffer.data(df_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_z_zz_xxx = cbuffer.data(df_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_zz_xxy = cbuffer.data(df_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_zz_xxz = cbuffer.data(df_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_zz_xyy = cbuffer.data(df_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_zz_xyz = cbuffer.data(df_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_zz_xzz = cbuffer.data(df_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_zz_yyy = cbuffer.data(df_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_z_zz_yyz = cbuffer.data(df_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_zz_yzz = cbuffer.data(df_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_zz_zzz = cbuffer.data(df_geom_11_off + 359 * ccomps * dcomps);

            auto g_z_x_xx_xxx = cbuffer.data(df_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_x_xx_xxy = cbuffer.data(df_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_x_xx_xxz = cbuffer.data(df_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_x_xx_xyy = cbuffer.data(df_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_x_xx_xyz = cbuffer.data(df_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_x_xx_xzz = cbuffer.data(df_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_x_xx_yyy = cbuffer.data(df_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_x_xx_yyz = cbuffer.data(df_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_x_xx_yzz = cbuffer.data(df_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_x_xx_zzz = cbuffer.data(df_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_x_xy_xxx = cbuffer.data(df_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_x_xy_xxy = cbuffer.data(df_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_x_xy_xxz = cbuffer.data(df_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_x_xy_xyy = cbuffer.data(df_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_x_xy_xyz = cbuffer.data(df_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_x_xy_xzz = cbuffer.data(df_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_x_xy_yyy = cbuffer.data(df_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_x_xy_yyz = cbuffer.data(df_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_x_xy_yzz = cbuffer.data(df_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_xy_zzz = cbuffer.data(df_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_x_xz_xxx = cbuffer.data(df_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_xz_xxy = cbuffer.data(df_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_xz_xxz = cbuffer.data(df_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_xz_xyy = cbuffer.data(df_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_xz_xyz = cbuffer.data(df_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_xz_xzz = cbuffer.data(df_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_xz_yyy = cbuffer.data(df_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_xz_yyz = cbuffer.data(df_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_xz_yzz = cbuffer.data(df_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_xz_zzz = cbuffer.data(df_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_x_yy_xxx = cbuffer.data(df_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_yy_xxy = cbuffer.data(df_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_yy_xxz = cbuffer.data(df_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_yy_xyy = cbuffer.data(df_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_yy_xyz = cbuffer.data(df_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_yy_xzz = cbuffer.data(df_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_yy_yyy = cbuffer.data(df_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_yy_yyz = cbuffer.data(df_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_yy_yzz = cbuffer.data(df_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_x_yy_zzz = cbuffer.data(df_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_x_yz_xxx = cbuffer.data(df_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_yz_xxy = cbuffer.data(df_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_yz_xxz = cbuffer.data(df_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_yz_xyy = cbuffer.data(df_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_yz_xyz = cbuffer.data(df_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_yz_xzz = cbuffer.data(df_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_yz_yyy = cbuffer.data(df_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_yz_yyz = cbuffer.data(df_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_yz_yzz = cbuffer.data(df_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_yz_zzz = cbuffer.data(df_geom_11_off + 409 * ccomps * dcomps);

            auto g_z_x_zz_xxx = cbuffer.data(df_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_zz_xxy = cbuffer.data(df_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_zz_xxz = cbuffer.data(df_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_zz_xyy = cbuffer.data(df_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_zz_xyz = cbuffer.data(df_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_zz_xzz = cbuffer.data(df_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_zz_yyy = cbuffer.data(df_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_zz_yyz = cbuffer.data(df_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_zz_yzz = cbuffer.data(df_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_zz_zzz = cbuffer.data(df_geom_11_off + 419 * ccomps * dcomps);

            auto g_z_y_xx_xxx = cbuffer.data(df_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_y_xx_xxy = cbuffer.data(df_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_y_xx_xxz = cbuffer.data(df_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_y_xx_xyy = cbuffer.data(df_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_y_xx_xyz = cbuffer.data(df_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_y_xx_xzz = cbuffer.data(df_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_y_xx_yyy = cbuffer.data(df_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_y_xx_yyz = cbuffer.data(df_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_y_xx_yzz = cbuffer.data(df_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_y_xx_zzz = cbuffer.data(df_geom_11_off + 429 * ccomps * dcomps);

            auto g_z_y_xy_xxx = cbuffer.data(df_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_y_xy_xxy = cbuffer.data(df_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_y_xy_xxz = cbuffer.data(df_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_y_xy_xyy = cbuffer.data(df_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_y_xy_xyz = cbuffer.data(df_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_y_xy_xzz = cbuffer.data(df_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_y_xy_yyy = cbuffer.data(df_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_y_xy_yyz = cbuffer.data(df_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_y_xy_yzz = cbuffer.data(df_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_y_xy_zzz = cbuffer.data(df_geom_11_off + 439 * ccomps * dcomps);

            auto g_z_y_xz_xxx = cbuffer.data(df_geom_11_off + 440 * ccomps * dcomps);

            auto g_z_y_xz_xxy = cbuffer.data(df_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_xz_xxz = cbuffer.data(df_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_xz_xyy = cbuffer.data(df_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_xz_xyz = cbuffer.data(df_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_xz_xzz = cbuffer.data(df_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_xz_yyy = cbuffer.data(df_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_xz_yyz = cbuffer.data(df_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_xz_yzz = cbuffer.data(df_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_xz_zzz = cbuffer.data(df_geom_11_off + 449 * ccomps * dcomps);

            auto g_z_y_yy_xxx = cbuffer.data(df_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_yy_xxy = cbuffer.data(df_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_yy_xxz = cbuffer.data(df_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_yy_xyy = cbuffer.data(df_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_yy_xyz = cbuffer.data(df_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_yy_xzz = cbuffer.data(df_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_yy_yyy = cbuffer.data(df_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_yy_yyz = cbuffer.data(df_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_yy_yzz = cbuffer.data(df_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_yy_zzz = cbuffer.data(df_geom_11_off + 459 * ccomps * dcomps);

            auto g_z_y_yz_xxx = cbuffer.data(df_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_yz_xxy = cbuffer.data(df_geom_11_off + 461 * ccomps * dcomps);

            auto g_z_y_yz_xxz = cbuffer.data(df_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_yz_xyy = cbuffer.data(df_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_yz_xyz = cbuffer.data(df_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_yz_xzz = cbuffer.data(df_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_yz_yyy = cbuffer.data(df_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_yz_yyz = cbuffer.data(df_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_yz_yzz = cbuffer.data(df_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_yz_zzz = cbuffer.data(df_geom_11_off + 469 * ccomps * dcomps);

            auto g_z_y_zz_xxx = cbuffer.data(df_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_zz_xxy = cbuffer.data(df_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_zz_xxz = cbuffer.data(df_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_zz_xyy = cbuffer.data(df_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_zz_xyz = cbuffer.data(df_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_zz_xzz = cbuffer.data(df_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_zz_yyy = cbuffer.data(df_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_zz_yyz = cbuffer.data(df_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_zz_yzz = cbuffer.data(df_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_zz_zzz = cbuffer.data(df_geom_11_off + 479 * ccomps * dcomps);

            auto g_z_z_xx_xxx = cbuffer.data(df_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_z_xx_xxy = cbuffer.data(df_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_z_xx_xxz = cbuffer.data(df_geom_11_off + 482 * ccomps * dcomps);

            auto g_z_z_xx_xyy = cbuffer.data(df_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_z_xx_xyz = cbuffer.data(df_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_z_xx_xzz = cbuffer.data(df_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_z_xx_yyy = cbuffer.data(df_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_z_xx_yyz = cbuffer.data(df_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_z_xx_yzz = cbuffer.data(df_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_z_xx_zzz = cbuffer.data(df_geom_11_off + 489 * ccomps * dcomps);

            auto g_z_z_xy_xxx = cbuffer.data(df_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_z_xy_xxy = cbuffer.data(df_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_z_xy_xxz = cbuffer.data(df_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_z_xy_xyy = cbuffer.data(df_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_z_xy_xyz = cbuffer.data(df_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_z_xy_xzz = cbuffer.data(df_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_z_xy_yyy = cbuffer.data(df_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_z_xy_yyz = cbuffer.data(df_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_z_xy_yzz = cbuffer.data(df_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_z_xy_zzz = cbuffer.data(df_geom_11_off + 499 * ccomps * dcomps);

            auto g_z_z_xz_xxx = cbuffer.data(df_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_z_xz_xxy = cbuffer.data(df_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_z_xz_xxz = cbuffer.data(df_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_z_xz_xyy = cbuffer.data(df_geom_11_off + 503 * ccomps * dcomps);

            auto g_z_z_xz_xyz = cbuffer.data(df_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_xz_xzz = cbuffer.data(df_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_xz_yyy = cbuffer.data(df_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_xz_yyz = cbuffer.data(df_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_xz_yzz = cbuffer.data(df_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_xz_zzz = cbuffer.data(df_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_z_yy_xxx = cbuffer.data(df_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_yy_xxy = cbuffer.data(df_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_yy_xxz = cbuffer.data(df_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_yy_xyy = cbuffer.data(df_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_yy_xyz = cbuffer.data(df_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_yy_xzz = cbuffer.data(df_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_yy_yyy = cbuffer.data(df_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_yy_yyz = cbuffer.data(df_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_yy_yzz = cbuffer.data(df_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_yy_zzz = cbuffer.data(df_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_z_yz_xxx = cbuffer.data(df_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_yz_xxy = cbuffer.data(df_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_yz_xxz = cbuffer.data(df_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_yz_xyy = cbuffer.data(df_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_yz_xyz = cbuffer.data(df_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_z_yz_xzz = cbuffer.data(df_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_yz_yyy = cbuffer.data(df_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_yz_yyz = cbuffer.data(df_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_yz_yzz = cbuffer.data(df_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_yz_zzz = cbuffer.data(df_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_z_zz_xxx = cbuffer.data(df_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_zz_xxy = cbuffer.data(df_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_zz_xxz = cbuffer.data(df_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_zz_xyy = cbuffer.data(df_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_zz_xyz = cbuffer.data(df_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_zz_xzz = cbuffer.data(df_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_zz_yyy = cbuffer.data(df_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_zz_yyz = cbuffer.data(df_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_zz_yzz = cbuffer.data(df_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_zz_zzz = cbuffer.data(df_geom_11_off + 539 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_ffxx

            const auto ff_geom_11_off = idx_geom_11_ffxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxx_xxx = cbuffer.data(ff_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxx_xxy = cbuffer.data(ff_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxx_xxz = cbuffer.data(ff_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxx_xyy = cbuffer.data(ff_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxx_xyz = cbuffer.data(ff_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxx_xzz = cbuffer.data(ff_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxx_yyy = cbuffer.data(ff_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxx_yyz = cbuffer.data(ff_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxx_yzz = cbuffer.data(ff_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxx_zzz = cbuffer.data(ff_geom_11_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_xxx, g_0_x_xx_xxy, g_0_x_xx_xxz, g_0_x_xx_xyy, g_0_x_xx_xyz, g_0_x_xx_xzz, g_0_x_xx_yyy, g_0_x_xx_yyz, g_0_x_xx_yzz, g_0_x_xx_zzz, g_x_0_xx_xxx, g_x_0_xx_xxy, g_x_0_xx_xxz, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xzz, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yzz, g_x_0_xx_zzz, g_x_x_xx_xxx, g_x_x_xx_xxxx, g_x_x_xx_xxxy, g_x_x_xx_xxxz, g_x_x_xx_xxy, g_x_x_xx_xxyy, g_x_x_xx_xxyz, g_x_x_xx_xxz, g_x_x_xx_xxzz, g_x_x_xx_xyy, g_x_x_xx_xyyy, g_x_x_xx_xyyz, g_x_x_xx_xyz, g_x_x_xx_xyzz, g_x_x_xx_xzz, g_x_x_xx_xzzz, g_x_x_xx_yyy, g_x_x_xx_yyz, g_x_x_xx_yzz, g_x_x_xx_zzz, g_x_x_xxx_xxx, g_x_x_xxx_xxy, g_x_x_xxx_xxz, g_x_x_xxx_xyy, g_x_x_xxx_xyz, g_x_x_xxx_xzz, g_x_x_xxx_yyy, g_x_x_xxx_yyz, g_x_x_xxx_yzz, g_x_x_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxx_xxx[k] = -g_0_x_xx_xxx[k] + g_x_0_xx_xxx[k] - g_x_x_xx_xxx[k] * ab_x + g_x_x_xx_xxxx[k];

                g_x_x_xxx_xxy[k] = -g_0_x_xx_xxy[k] + g_x_0_xx_xxy[k] - g_x_x_xx_xxy[k] * ab_x + g_x_x_xx_xxxy[k];

                g_x_x_xxx_xxz[k] = -g_0_x_xx_xxz[k] + g_x_0_xx_xxz[k] - g_x_x_xx_xxz[k] * ab_x + g_x_x_xx_xxxz[k];

                g_x_x_xxx_xyy[k] = -g_0_x_xx_xyy[k] + g_x_0_xx_xyy[k] - g_x_x_xx_xyy[k] * ab_x + g_x_x_xx_xxyy[k];

                g_x_x_xxx_xyz[k] = -g_0_x_xx_xyz[k] + g_x_0_xx_xyz[k] - g_x_x_xx_xyz[k] * ab_x + g_x_x_xx_xxyz[k];

                g_x_x_xxx_xzz[k] = -g_0_x_xx_xzz[k] + g_x_0_xx_xzz[k] - g_x_x_xx_xzz[k] * ab_x + g_x_x_xx_xxzz[k];

                g_x_x_xxx_yyy[k] = -g_0_x_xx_yyy[k] + g_x_0_xx_yyy[k] - g_x_x_xx_yyy[k] * ab_x + g_x_x_xx_xyyy[k];

                g_x_x_xxx_yyz[k] = -g_0_x_xx_yyz[k] + g_x_0_xx_yyz[k] - g_x_x_xx_yyz[k] * ab_x + g_x_x_xx_xyyz[k];

                g_x_x_xxx_yzz[k] = -g_0_x_xx_yzz[k] + g_x_0_xx_yzz[k] - g_x_x_xx_yzz[k] * ab_x + g_x_x_xx_xyzz[k];

                g_x_x_xxx_zzz[k] = -g_0_x_xx_zzz[k] + g_x_0_xx_zzz[k] - g_x_x_xx_zzz[k] * ab_x + g_x_x_xx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxy_xxx = cbuffer.data(ff_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxy_xxy = cbuffer.data(ff_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxy_xxz = cbuffer.data(ff_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxy_xyy = cbuffer.data(ff_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxy_xyz = cbuffer.data(ff_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxy_xzz = cbuffer.data(ff_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxy_yyy = cbuffer.data(ff_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxy_yyz = cbuffer.data(ff_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xxy_yzz = cbuffer.data(ff_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxy_zzz = cbuffer.data(ff_geom_11_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_xxx, g_x_x_xx_xxxy, g_x_x_xx_xxy, g_x_x_xx_xxyy, g_x_x_xx_xxyz, g_x_x_xx_xxz, g_x_x_xx_xyy, g_x_x_xx_xyyy, g_x_x_xx_xyyz, g_x_x_xx_xyz, g_x_x_xx_xyzz, g_x_x_xx_xzz, g_x_x_xx_yyy, g_x_x_xx_yyyy, g_x_x_xx_yyyz, g_x_x_xx_yyz, g_x_x_xx_yyzz, g_x_x_xx_yzz, g_x_x_xx_yzzz, g_x_x_xx_zzz, g_x_x_xxy_xxx, g_x_x_xxy_xxy, g_x_x_xxy_xxz, g_x_x_xxy_xyy, g_x_x_xxy_xyz, g_x_x_xxy_xzz, g_x_x_xxy_yyy, g_x_x_xxy_yyz, g_x_x_xxy_yzz, g_x_x_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxy_xxx[k] = -g_x_x_xx_xxx[k] * ab_y + g_x_x_xx_xxxy[k];

                g_x_x_xxy_xxy[k] = -g_x_x_xx_xxy[k] * ab_y + g_x_x_xx_xxyy[k];

                g_x_x_xxy_xxz[k] = -g_x_x_xx_xxz[k] * ab_y + g_x_x_xx_xxyz[k];

                g_x_x_xxy_xyy[k] = -g_x_x_xx_xyy[k] * ab_y + g_x_x_xx_xyyy[k];

                g_x_x_xxy_xyz[k] = -g_x_x_xx_xyz[k] * ab_y + g_x_x_xx_xyyz[k];

                g_x_x_xxy_xzz[k] = -g_x_x_xx_xzz[k] * ab_y + g_x_x_xx_xyzz[k];

                g_x_x_xxy_yyy[k] = -g_x_x_xx_yyy[k] * ab_y + g_x_x_xx_yyyy[k];

                g_x_x_xxy_yyz[k] = -g_x_x_xx_yyz[k] * ab_y + g_x_x_xx_yyyz[k];

                g_x_x_xxy_yzz[k] = -g_x_x_xx_yzz[k] * ab_y + g_x_x_xx_yyzz[k];

                g_x_x_xxy_zzz[k] = -g_x_x_xx_zzz[k] * ab_y + g_x_x_xx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxz_xxx = cbuffer.data(ff_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxz_xxy = cbuffer.data(ff_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxz_xxz = cbuffer.data(ff_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxz_xyy = cbuffer.data(ff_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xxz_xyz = cbuffer.data(ff_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxz_xzz = cbuffer.data(ff_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxz_yyy = cbuffer.data(ff_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxz_yyz = cbuffer.data(ff_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxz_yzz = cbuffer.data(ff_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxz_zzz = cbuffer.data(ff_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xx_xxx, g_x_x_xx_xxxz, g_x_x_xx_xxy, g_x_x_xx_xxyz, g_x_x_xx_xxz, g_x_x_xx_xxzz, g_x_x_xx_xyy, g_x_x_xx_xyyz, g_x_x_xx_xyz, g_x_x_xx_xyzz, g_x_x_xx_xzz, g_x_x_xx_xzzz, g_x_x_xx_yyy, g_x_x_xx_yyyz, g_x_x_xx_yyz, g_x_x_xx_yyzz, g_x_x_xx_yzz, g_x_x_xx_yzzz, g_x_x_xx_zzz, g_x_x_xx_zzzz, g_x_x_xxz_xxx, g_x_x_xxz_xxy, g_x_x_xxz_xxz, g_x_x_xxz_xyy, g_x_x_xxz_xyz, g_x_x_xxz_xzz, g_x_x_xxz_yyy, g_x_x_xxz_yyz, g_x_x_xxz_yzz, g_x_x_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxz_xxx[k] = -g_x_x_xx_xxx[k] * ab_z + g_x_x_xx_xxxz[k];

                g_x_x_xxz_xxy[k] = -g_x_x_xx_xxy[k] * ab_z + g_x_x_xx_xxyz[k];

                g_x_x_xxz_xxz[k] = -g_x_x_xx_xxz[k] * ab_z + g_x_x_xx_xxzz[k];

                g_x_x_xxz_xyy[k] = -g_x_x_xx_xyy[k] * ab_z + g_x_x_xx_xyyz[k];

                g_x_x_xxz_xyz[k] = -g_x_x_xx_xyz[k] * ab_z + g_x_x_xx_xyzz[k];

                g_x_x_xxz_xzz[k] = -g_x_x_xx_xzz[k] * ab_z + g_x_x_xx_xzzz[k];

                g_x_x_xxz_yyy[k] = -g_x_x_xx_yyy[k] * ab_z + g_x_x_xx_yyyz[k];

                g_x_x_xxz_yyz[k] = -g_x_x_xx_yyz[k] * ab_z + g_x_x_xx_yyzz[k];

                g_x_x_xxz_yzz[k] = -g_x_x_xx_yzz[k] * ab_z + g_x_x_xx_yzzz[k];

                g_x_x_xxz_zzz[k] = -g_x_x_xx_zzz[k] * ab_z + g_x_x_xx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyy_xxx = cbuffer.data(ff_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xyy_xxy = cbuffer.data(ff_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xyy_xxz = cbuffer.data(ff_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xyy_xyy = cbuffer.data(ff_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xyy_xyz = cbuffer.data(ff_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xyy_xzz = cbuffer.data(ff_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xyy_yyy = cbuffer.data(ff_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xyy_yyz = cbuffer.data(ff_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xyy_yzz = cbuffer.data(ff_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xyy_zzz = cbuffer.data(ff_geom_11_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xy_xxx, g_x_x_xy_xxxy, g_x_x_xy_xxy, g_x_x_xy_xxyy, g_x_x_xy_xxyz, g_x_x_xy_xxz, g_x_x_xy_xyy, g_x_x_xy_xyyy, g_x_x_xy_xyyz, g_x_x_xy_xyz, g_x_x_xy_xyzz, g_x_x_xy_xzz, g_x_x_xy_yyy, g_x_x_xy_yyyy, g_x_x_xy_yyyz, g_x_x_xy_yyz, g_x_x_xy_yyzz, g_x_x_xy_yzz, g_x_x_xy_yzzz, g_x_x_xy_zzz, g_x_x_xyy_xxx, g_x_x_xyy_xxy, g_x_x_xyy_xxz, g_x_x_xyy_xyy, g_x_x_xyy_xyz, g_x_x_xyy_xzz, g_x_x_xyy_yyy, g_x_x_xyy_yyz, g_x_x_xyy_yzz, g_x_x_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyy_xxx[k] = -g_x_x_xy_xxx[k] * ab_y + g_x_x_xy_xxxy[k];

                g_x_x_xyy_xxy[k] = -g_x_x_xy_xxy[k] * ab_y + g_x_x_xy_xxyy[k];

                g_x_x_xyy_xxz[k] = -g_x_x_xy_xxz[k] * ab_y + g_x_x_xy_xxyz[k];

                g_x_x_xyy_xyy[k] = -g_x_x_xy_xyy[k] * ab_y + g_x_x_xy_xyyy[k];

                g_x_x_xyy_xyz[k] = -g_x_x_xy_xyz[k] * ab_y + g_x_x_xy_xyyz[k];

                g_x_x_xyy_xzz[k] = -g_x_x_xy_xzz[k] * ab_y + g_x_x_xy_xyzz[k];

                g_x_x_xyy_yyy[k] = -g_x_x_xy_yyy[k] * ab_y + g_x_x_xy_yyyy[k];

                g_x_x_xyy_yyz[k] = -g_x_x_xy_yyz[k] * ab_y + g_x_x_xy_yyyz[k];

                g_x_x_xyy_yzz[k] = -g_x_x_xy_yzz[k] * ab_y + g_x_x_xy_yyzz[k];

                g_x_x_xyy_zzz[k] = -g_x_x_xy_zzz[k] * ab_y + g_x_x_xy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyz_xxx = cbuffer.data(ff_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xyz_xxy = cbuffer.data(ff_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xyz_xxz = cbuffer.data(ff_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xyz_xyy = cbuffer.data(ff_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xyz_xyz = cbuffer.data(ff_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xyz_xzz = cbuffer.data(ff_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xyz_yyy = cbuffer.data(ff_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xyz_yyz = cbuffer.data(ff_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xyz_yzz = cbuffer.data(ff_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xyz_zzz = cbuffer.data(ff_geom_11_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyz_xxx, g_x_x_xyz_xxy, g_x_x_xyz_xxz, g_x_x_xyz_xyy, g_x_x_xyz_xyz, g_x_x_xyz_xzz, g_x_x_xyz_yyy, g_x_x_xyz_yyz, g_x_x_xyz_yzz, g_x_x_xyz_zzz, g_x_x_xz_xxx, g_x_x_xz_xxxy, g_x_x_xz_xxy, g_x_x_xz_xxyy, g_x_x_xz_xxyz, g_x_x_xz_xxz, g_x_x_xz_xyy, g_x_x_xz_xyyy, g_x_x_xz_xyyz, g_x_x_xz_xyz, g_x_x_xz_xyzz, g_x_x_xz_xzz, g_x_x_xz_yyy, g_x_x_xz_yyyy, g_x_x_xz_yyyz, g_x_x_xz_yyz, g_x_x_xz_yyzz, g_x_x_xz_yzz, g_x_x_xz_yzzz, g_x_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyz_xxx[k] = -g_x_x_xz_xxx[k] * ab_y + g_x_x_xz_xxxy[k];

                g_x_x_xyz_xxy[k] = -g_x_x_xz_xxy[k] * ab_y + g_x_x_xz_xxyy[k];

                g_x_x_xyz_xxz[k] = -g_x_x_xz_xxz[k] * ab_y + g_x_x_xz_xxyz[k];

                g_x_x_xyz_xyy[k] = -g_x_x_xz_xyy[k] * ab_y + g_x_x_xz_xyyy[k];

                g_x_x_xyz_xyz[k] = -g_x_x_xz_xyz[k] * ab_y + g_x_x_xz_xyyz[k];

                g_x_x_xyz_xzz[k] = -g_x_x_xz_xzz[k] * ab_y + g_x_x_xz_xyzz[k];

                g_x_x_xyz_yyy[k] = -g_x_x_xz_yyy[k] * ab_y + g_x_x_xz_yyyy[k];

                g_x_x_xyz_yyz[k] = -g_x_x_xz_yyz[k] * ab_y + g_x_x_xz_yyyz[k];

                g_x_x_xyz_yzz[k] = -g_x_x_xz_yzz[k] * ab_y + g_x_x_xz_yyzz[k];

                g_x_x_xyz_zzz[k] = -g_x_x_xz_zzz[k] * ab_y + g_x_x_xz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzz_xxx = cbuffer.data(ff_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xzz_xxy = cbuffer.data(ff_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xzz_xxz = cbuffer.data(ff_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xzz_xyy = cbuffer.data(ff_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xzz_xyz = cbuffer.data(ff_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xzz_xzz = cbuffer.data(ff_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xzz_yyy = cbuffer.data(ff_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xzz_yyz = cbuffer.data(ff_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xzz_yzz = cbuffer.data(ff_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xzz_zzz = cbuffer.data(ff_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xz_xxx, g_x_x_xz_xxxz, g_x_x_xz_xxy, g_x_x_xz_xxyz, g_x_x_xz_xxz, g_x_x_xz_xxzz, g_x_x_xz_xyy, g_x_x_xz_xyyz, g_x_x_xz_xyz, g_x_x_xz_xyzz, g_x_x_xz_xzz, g_x_x_xz_xzzz, g_x_x_xz_yyy, g_x_x_xz_yyyz, g_x_x_xz_yyz, g_x_x_xz_yyzz, g_x_x_xz_yzz, g_x_x_xz_yzzz, g_x_x_xz_zzz, g_x_x_xz_zzzz, g_x_x_xzz_xxx, g_x_x_xzz_xxy, g_x_x_xzz_xxz, g_x_x_xzz_xyy, g_x_x_xzz_xyz, g_x_x_xzz_xzz, g_x_x_xzz_yyy, g_x_x_xzz_yyz, g_x_x_xzz_yzz, g_x_x_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzz_xxx[k] = -g_x_x_xz_xxx[k] * ab_z + g_x_x_xz_xxxz[k];

                g_x_x_xzz_xxy[k] = -g_x_x_xz_xxy[k] * ab_z + g_x_x_xz_xxyz[k];

                g_x_x_xzz_xxz[k] = -g_x_x_xz_xxz[k] * ab_z + g_x_x_xz_xxzz[k];

                g_x_x_xzz_xyy[k] = -g_x_x_xz_xyy[k] * ab_z + g_x_x_xz_xyyz[k];

                g_x_x_xzz_xyz[k] = -g_x_x_xz_xyz[k] * ab_z + g_x_x_xz_xyzz[k];

                g_x_x_xzz_xzz[k] = -g_x_x_xz_xzz[k] * ab_z + g_x_x_xz_xzzz[k];

                g_x_x_xzz_yyy[k] = -g_x_x_xz_yyy[k] * ab_z + g_x_x_xz_yyyz[k];

                g_x_x_xzz_yyz[k] = -g_x_x_xz_yyz[k] * ab_z + g_x_x_xz_yyzz[k];

                g_x_x_xzz_yzz[k] = -g_x_x_xz_yzz[k] * ab_z + g_x_x_xz_yzzz[k];

                g_x_x_xzz_zzz[k] = -g_x_x_xz_zzz[k] * ab_z + g_x_x_xz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyy_xxx = cbuffer.data(ff_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_yyy_xxy = cbuffer.data(ff_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_yyy_xxz = cbuffer.data(ff_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yyy_xyy = cbuffer.data(ff_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yyy_xyz = cbuffer.data(ff_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yyy_xzz = cbuffer.data(ff_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yyy_yyy = cbuffer.data(ff_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yyy_yyz = cbuffer.data(ff_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yyy_yzz = cbuffer.data(ff_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yyy_zzz = cbuffer.data(ff_geom_11_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yy_xxx, g_x_x_yy_xxxy, g_x_x_yy_xxy, g_x_x_yy_xxyy, g_x_x_yy_xxyz, g_x_x_yy_xxz, g_x_x_yy_xyy, g_x_x_yy_xyyy, g_x_x_yy_xyyz, g_x_x_yy_xyz, g_x_x_yy_xyzz, g_x_x_yy_xzz, g_x_x_yy_yyy, g_x_x_yy_yyyy, g_x_x_yy_yyyz, g_x_x_yy_yyz, g_x_x_yy_yyzz, g_x_x_yy_yzz, g_x_x_yy_yzzz, g_x_x_yy_zzz, g_x_x_yyy_xxx, g_x_x_yyy_xxy, g_x_x_yyy_xxz, g_x_x_yyy_xyy, g_x_x_yyy_xyz, g_x_x_yyy_xzz, g_x_x_yyy_yyy, g_x_x_yyy_yyz, g_x_x_yyy_yzz, g_x_x_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyy_xxx[k] = -g_x_x_yy_xxx[k] * ab_y + g_x_x_yy_xxxy[k];

                g_x_x_yyy_xxy[k] = -g_x_x_yy_xxy[k] * ab_y + g_x_x_yy_xxyy[k];

                g_x_x_yyy_xxz[k] = -g_x_x_yy_xxz[k] * ab_y + g_x_x_yy_xxyz[k];

                g_x_x_yyy_xyy[k] = -g_x_x_yy_xyy[k] * ab_y + g_x_x_yy_xyyy[k];

                g_x_x_yyy_xyz[k] = -g_x_x_yy_xyz[k] * ab_y + g_x_x_yy_xyyz[k];

                g_x_x_yyy_xzz[k] = -g_x_x_yy_xzz[k] * ab_y + g_x_x_yy_xyzz[k];

                g_x_x_yyy_yyy[k] = -g_x_x_yy_yyy[k] * ab_y + g_x_x_yy_yyyy[k];

                g_x_x_yyy_yyz[k] = -g_x_x_yy_yyz[k] * ab_y + g_x_x_yy_yyyz[k];

                g_x_x_yyy_yzz[k] = -g_x_x_yy_yzz[k] * ab_y + g_x_x_yy_yyzz[k];

                g_x_x_yyy_zzz[k] = -g_x_x_yy_zzz[k] * ab_y + g_x_x_yy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyz_xxx = cbuffer.data(ff_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yyz_xxy = cbuffer.data(ff_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yyz_xxz = cbuffer.data(ff_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yyz_xyy = cbuffer.data(ff_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yyz_xyz = cbuffer.data(ff_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_yyz_xzz = cbuffer.data(ff_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_yyz_yyy = cbuffer.data(ff_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_yyz_yyz = cbuffer.data(ff_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_yyz_yzz = cbuffer.data(ff_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_yyz_zzz = cbuffer.data(ff_geom_11_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyz_xxx, g_x_x_yyz_xxy, g_x_x_yyz_xxz, g_x_x_yyz_xyy, g_x_x_yyz_xyz, g_x_x_yyz_xzz, g_x_x_yyz_yyy, g_x_x_yyz_yyz, g_x_x_yyz_yzz, g_x_x_yyz_zzz, g_x_x_yz_xxx, g_x_x_yz_xxxy, g_x_x_yz_xxy, g_x_x_yz_xxyy, g_x_x_yz_xxyz, g_x_x_yz_xxz, g_x_x_yz_xyy, g_x_x_yz_xyyy, g_x_x_yz_xyyz, g_x_x_yz_xyz, g_x_x_yz_xyzz, g_x_x_yz_xzz, g_x_x_yz_yyy, g_x_x_yz_yyyy, g_x_x_yz_yyyz, g_x_x_yz_yyz, g_x_x_yz_yyzz, g_x_x_yz_yzz, g_x_x_yz_yzzz, g_x_x_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyz_xxx[k] = -g_x_x_yz_xxx[k] * ab_y + g_x_x_yz_xxxy[k];

                g_x_x_yyz_xxy[k] = -g_x_x_yz_xxy[k] * ab_y + g_x_x_yz_xxyy[k];

                g_x_x_yyz_xxz[k] = -g_x_x_yz_xxz[k] * ab_y + g_x_x_yz_xxyz[k];

                g_x_x_yyz_xyy[k] = -g_x_x_yz_xyy[k] * ab_y + g_x_x_yz_xyyy[k];

                g_x_x_yyz_xyz[k] = -g_x_x_yz_xyz[k] * ab_y + g_x_x_yz_xyyz[k];

                g_x_x_yyz_xzz[k] = -g_x_x_yz_xzz[k] * ab_y + g_x_x_yz_xyzz[k];

                g_x_x_yyz_yyy[k] = -g_x_x_yz_yyy[k] * ab_y + g_x_x_yz_yyyy[k];

                g_x_x_yyz_yyz[k] = -g_x_x_yz_yyz[k] * ab_y + g_x_x_yz_yyyz[k];

                g_x_x_yyz_yzz[k] = -g_x_x_yz_yzz[k] * ab_y + g_x_x_yz_yyzz[k];

                g_x_x_yyz_zzz[k] = -g_x_x_yz_zzz[k] * ab_y + g_x_x_yz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzz_xxx = cbuffer.data(ff_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_yzz_xxy = cbuffer.data(ff_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_yzz_xxz = cbuffer.data(ff_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_yzz_xyy = cbuffer.data(ff_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_yzz_xyz = cbuffer.data(ff_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_yzz_xzz = cbuffer.data(ff_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_yzz_yyy = cbuffer.data(ff_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_yzz_yyz = cbuffer.data(ff_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_yzz_yzz = cbuffer.data(ff_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_yzz_zzz = cbuffer.data(ff_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzz_xxx, g_x_x_yzz_xxy, g_x_x_yzz_xxz, g_x_x_yzz_xyy, g_x_x_yzz_xyz, g_x_x_yzz_xzz, g_x_x_yzz_yyy, g_x_x_yzz_yyz, g_x_x_yzz_yzz, g_x_x_yzz_zzz, g_x_x_zz_xxx, g_x_x_zz_xxxy, g_x_x_zz_xxy, g_x_x_zz_xxyy, g_x_x_zz_xxyz, g_x_x_zz_xxz, g_x_x_zz_xyy, g_x_x_zz_xyyy, g_x_x_zz_xyyz, g_x_x_zz_xyz, g_x_x_zz_xyzz, g_x_x_zz_xzz, g_x_x_zz_yyy, g_x_x_zz_yyyy, g_x_x_zz_yyyz, g_x_x_zz_yyz, g_x_x_zz_yyzz, g_x_x_zz_yzz, g_x_x_zz_yzzz, g_x_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzz_xxx[k] = -g_x_x_zz_xxx[k] * ab_y + g_x_x_zz_xxxy[k];

                g_x_x_yzz_xxy[k] = -g_x_x_zz_xxy[k] * ab_y + g_x_x_zz_xxyy[k];

                g_x_x_yzz_xxz[k] = -g_x_x_zz_xxz[k] * ab_y + g_x_x_zz_xxyz[k];

                g_x_x_yzz_xyy[k] = -g_x_x_zz_xyy[k] * ab_y + g_x_x_zz_xyyy[k];

                g_x_x_yzz_xyz[k] = -g_x_x_zz_xyz[k] * ab_y + g_x_x_zz_xyyz[k];

                g_x_x_yzz_xzz[k] = -g_x_x_zz_xzz[k] * ab_y + g_x_x_zz_xyzz[k];

                g_x_x_yzz_yyy[k] = -g_x_x_zz_yyy[k] * ab_y + g_x_x_zz_yyyy[k];

                g_x_x_yzz_yyz[k] = -g_x_x_zz_yyz[k] * ab_y + g_x_x_zz_yyyz[k];

                g_x_x_yzz_yzz[k] = -g_x_x_zz_yzz[k] * ab_y + g_x_x_zz_yyzz[k];

                g_x_x_yzz_zzz[k] = -g_x_x_zz_zzz[k] * ab_y + g_x_x_zz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzz_xxx = cbuffer.data(ff_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_zzz_xxy = cbuffer.data(ff_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_zzz_xxz = cbuffer.data(ff_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_zzz_xyy = cbuffer.data(ff_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_zzz_xyz = cbuffer.data(ff_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_zzz_xzz = cbuffer.data(ff_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_zzz_yyy = cbuffer.data(ff_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_zzz_yyz = cbuffer.data(ff_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_zzz_yzz = cbuffer.data(ff_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_zzz_zzz = cbuffer.data(ff_geom_11_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zz_xxx, g_x_x_zz_xxxz, g_x_x_zz_xxy, g_x_x_zz_xxyz, g_x_x_zz_xxz, g_x_x_zz_xxzz, g_x_x_zz_xyy, g_x_x_zz_xyyz, g_x_x_zz_xyz, g_x_x_zz_xyzz, g_x_x_zz_xzz, g_x_x_zz_xzzz, g_x_x_zz_yyy, g_x_x_zz_yyyz, g_x_x_zz_yyz, g_x_x_zz_yyzz, g_x_x_zz_yzz, g_x_x_zz_yzzz, g_x_x_zz_zzz, g_x_x_zz_zzzz, g_x_x_zzz_xxx, g_x_x_zzz_xxy, g_x_x_zzz_xxz, g_x_x_zzz_xyy, g_x_x_zzz_xyz, g_x_x_zzz_xzz, g_x_x_zzz_yyy, g_x_x_zzz_yyz, g_x_x_zzz_yzz, g_x_x_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzz_xxx[k] = -g_x_x_zz_xxx[k] * ab_z + g_x_x_zz_xxxz[k];

                g_x_x_zzz_xxy[k] = -g_x_x_zz_xxy[k] * ab_z + g_x_x_zz_xxyz[k];

                g_x_x_zzz_xxz[k] = -g_x_x_zz_xxz[k] * ab_z + g_x_x_zz_xxzz[k];

                g_x_x_zzz_xyy[k] = -g_x_x_zz_xyy[k] * ab_z + g_x_x_zz_xyyz[k];

                g_x_x_zzz_xyz[k] = -g_x_x_zz_xyz[k] * ab_z + g_x_x_zz_xyzz[k];

                g_x_x_zzz_xzz[k] = -g_x_x_zz_xzz[k] * ab_z + g_x_x_zz_xzzz[k];

                g_x_x_zzz_yyy[k] = -g_x_x_zz_yyy[k] * ab_z + g_x_x_zz_yyyz[k];

                g_x_x_zzz_yyz[k] = -g_x_x_zz_yyz[k] * ab_z + g_x_x_zz_yyzz[k];

                g_x_x_zzz_yzz[k] = -g_x_x_zz_yzz[k] * ab_z + g_x_x_zz_yzzz[k];

                g_x_x_zzz_zzz[k] = -g_x_x_zz_zzz[k] * ab_z + g_x_x_zz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxx_xxx = cbuffer.data(ff_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_xxx_xxy = cbuffer.data(ff_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_xxx_xxz = cbuffer.data(ff_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_xxx_xyy = cbuffer.data(ff_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_xxx_xyz = cbuffer.data(ff_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_xxx_xzz = cbuffer.data(ff_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_xxx_yyy = cbuffer.data(ff_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_xxx_yyz = cbuffer.data(ff_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_xxx_yzz = cbuffer.data(ff_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_xxx_zzz = cbuffer.data(ff_geom_11_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_xxx, g_0_y_xx_xxy, g_0_y_xx_xxz, g_0_y_xx_xyy, g_0_y_xx_xyz, g_0_y_xx_xzz, g_0_y_xx_yyy, g_0_y_xx_yyz, g_0_y_xx_yzz, g_0_y_xx_zzz, g_x_y_xx_xxx, g_x_y_xx_xxxx, g_x_y_xx_xxxy, g_x_y_xx_xxxz, g_x_y_xx_xxy, g_x_y_xx_xxyy, g_x_y_xx_xxyz, g_x_y_xx_xxz, g_x_y_xx_xxzz, g_x_y_xx_xyy, g_x_y_xx_xyyy, g_x_y_xx_xyyz, g_x_y_xx_xyz, g_x_y_xx_xyzz, g_x_y_xx_xzz, g_x_y_xx_xzzz, g_x_y_xx_yyy, g_x_y_xx_yyz, g_x_y_xx_yzz, g_x_y_xx_zzz, g_x_y_xxx_xxx, g_x_y_xxx_xxy, g_x_y_xxx_xxz, g_x_y_xxx_xyy, g_x_y_xxx_xyz, g_x_y_xxx_xzz, g_x_y_xxx_yyy, g_x_y_xxx_yyz, g_x_y_xxx_yzz, g_x_y_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxx_xxx[k] = -g_0_y_xx_xxx[k] - g_x_y_xx_xxx[k] * ab_x + g_x_y_xx_xxxx[k];

                g_x_y_xxx_xxy[k] = -g_0_y_xx_xxy[k] - g_x_y_xx_xxy[k] * ab_x + g_x_y_xx_xxxy[k];

                g_x_y_xxx_xxz[k] = -g_0_y_xx_xxz[k] - g_x_y_xx_xxz[k] * ab_x + g_x_y_xx_xxxz[k];

                g_x_y_xxx_xyy[k] = -g_0_y_xx_xyy[k] - g_x_y_xx_xyy[k] * ab_x + g_x_y_xx_xxyy[k];

                g_x_y_xxx_xyz[k] = -g_0_y_xx_xyz[k] - g_x_y_xx_xyz[k] * ab_x + g_x_y_xx_xxyz[k];

                g_x_y_xxx_xzz[k] = -g_0_y_xx_xzz[k] - g_x_y_xx_xzz[k] * ab_x + g_x_y_xx_xxzz[k];

                g_x_y_xxx_yyy[k] = -g_0_y_xx_yyy[k] - g_x_y_xx_yyy[k] * ab_x + g_x_y_xx_xyyy[k];

                g_x_y_xxx_yyz[k] = -g_0_y_xx_yyz[k] - g_x_y_xx_yyz[k] * ab_x + g_x_y_xx_xyyz[k];

                g_x_y_xxx_yzz[k] = -g_0_y_xx_yzz[k] - g_x_y_xx_yzz[k] * ab_x + g_x_y_xx_xyzz[k];

                g_x_y_xxx_zzz[k] = -g_0_y_xx_zzz[k] - g_x_y_xx_zzz[k] * ab_x + g_x_y_xx_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxy_xxx = cbuffer.data(ff_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_xxy_xxy = cbuffer.data(ff_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_xxy_xxz = cbuffer.data(ff_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_xxy_xyy = cbuffer.data(ff_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_xxy_xyz = cbuffer.data(ff_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_xxy_xzz = cbuffer.data(ff_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_xxy_yyy = cbuffer.data(ff_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_xxy_yyz = cbuffer.data(ff_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_xxy_yzz = cbuffer.data(ff_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_xxy_zzz = cbuffer.data(ff_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_xxx, g_0_y_xy_xxy, g_0_y_xy_xxz, g_0_y_xy_xyy, g_0_y_xy_xyz, g_0_y_xy_xzz, g_0_y_xy_yyy, g_0_y_xy_yyz, g_0_y_xy_yzz, g_0_y_xy_zzz, g_x_y_xxy_xxx, g_x_y_xxy_xxy, g_x_y_xxy_xxz, g_x_y_xxy_xyy, g_x_y_xxy_xyz, g_x_y_xxy_xzz, g_x_y_xxy_yyy, g_x_y_xxy_yyz, g_x_y_xxy_yzz, g_x_y_xxy_zzz, g_x_y_xy_xxx, g_x_y_xy_xxxx, g_x_y_xy_xxxy, g_x_y_xy_xxxz, g_x_y_xy_xxy, g_x_y_xy_xxyy, g_x_y_xy_xxyz, g_x_y_xy_xxz, g_x_y_xy_xxzz, g_x_y_xy_xyy, g_x_y_xy_xyyy, g_x_y_xy_xyyz, g_x_y_xy_xyz, g_x_y_xy_xyzz, g_x_y_xy_xzz, g_x_y_xy_xzzz, g_x_y_xy_yyy, g_x_y_xy_yyz, g_x_y_xy_yzz, g_x_y_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxy_xxx[k] = -g_0_y_xy_xxx[k] - g_x_y_xy_xxx[k] * ab_x + g_x_y_xy_xxxx[k];

                g_x_y_xxy_xxy[k] = -g_0_y_xy_xxy[k] - g_x_y_xy_xxy[k] * ab_x + g_x_y_xy_xxxy[k];

                g_x_y_xxy_xxz[k] = -g_0_y_xy_xxz[k] - g_x_y_xy_xxz[k] * ab_x + g_x_y_xy_xxxz[k];

                g_x_y_xxy_xyy[k] = -g_0_y_xy_xyy[k] - g_x_y_xy_xyy[k] * ab_x + g_x_y_xy_xxyy[k];

                g_x_y_xxy_xyz[k] = -g_0_y_xy_xyz[k] - g_x_y_xy_xyz[k] * ab_x + g_x_y_xy_xxyz[k];

                g_x_y_xxy_xzz[k] = -g_0_y_xy_xzz[k] - g_x_y_xy_xzz[k] * ab_x + g_x_y_xy_xxzz[k];

                g_x_y_xxy_yyy[k] = -g_0_y_xy_yyy[k] - g_x_y_xy_yyy[k] * ab_x + g_x_y_xy_xyyy[k];

                g_x_y_xxy_yyz[k] = -g_0_y_xy_yyz[k] - g_x_y_xy_yyz[k] * ab_x + g_x_y_xy_xyyz[k];

                g_x_y_xxy_yzz[k] = -g_0_y_xy_yzz[k] - g_x_y_xy_yzz[k] * ab_x + g_x_y_xy_xyzz[k];

                g_x_y_xxy_zzz[k] = -g_0_y_xy_zzz[k] - g_x_y_xy_zzz[k] * ab_x + g_x_y_xy_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxz_xxx = cbuffer.data(ff_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_xxz_xxy = cbuffer.data(ff_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_xxz_xxz = cbuffer.data(ff_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_xxz_xyy = cbuffer.data(ff_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_xxz_xyz = cbuffer.data(ff_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_xxz_xzz = cbuffer.data(ff_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_xxz_yyy = cbuffer.data(ff_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xxz_yyz = cbuffer.data(ff_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xxz_yzz = cbuffer.data(ff_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xxz_zzz = cbuffer.data(ff_geom_11_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xx_xxx, g_x_y_xx_xxxz, g_x_y_xx_xxy, g_x_y_xx_xxyz, g_x_y_xx_xxz, g_x_y_xx_xxzz, g_x_y_xx_xyy, g_x_y_xx_xyyz, g_x_y_xx_xyz, g_x_y_xx_xyzz, g_x_y_xx_xzz, g_x_y_xx_xzzz, g_x_y_xx_yyy, g_x_y_xx_yyyz, g_x_y_xx_yyz, g_x_y_xx_yyzz, g_x_y_xx_yzz, g_x_y_xx_yzzz, g_x_y_xx_zzz, g_x_y_xx_zzzz, g_x_y_xxz_xxx, g_x_y_xxz_xxy, g_x_y_xxz_xxz, g_x_y_xxz_xyy, g_x_y_xxz_xyz, g_x_y_xxz_xzz, g_x_y_xxz_yyy, g_x_y_xxz_yyz, g_x_y_xxz_yzz, g_x_y_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxz_xxx[k] = -g_x_y_xx_xxx[k] * ab_z + g_x_y_xx_xxxz[k];

                g_x_y_xxz_xxy[k] = -g_x_y_xx_xxy[k] * ab_z + g_x_y_xx_xxyz[k];

                g_x_y_xxz_xxz[k] = -g_x_y_xx_xxz[k] * ab_z + g_x_y_xx_xxzz[k];

                g_x_y_xxz_xyy[k] = -g_x_y_xx_xyy[k] * ab_z + g_x_y_xx_xyyz[k];

                g_x_y_xxz_xyz[k] = -g_x_y_xx_xyz[k] * ab_z + g_x_y_xx_xyzz[k];

                g_x_y_xxz_xzz[k] = -g_x_y_xx_xzz[k] * ab_z + g_x_y_xx_xzzz[k];

                g_x_y_xxz_yyy[k] = -g_x_y_xx_yyy[k] * ab_z + g_x_y_xx_yyyz[k];

                g_x_y_xxz_yyz[k] = -g_x_y_xx_yyz[k] * ab_z + g_x_y_xx_yyzz[k];

                g_x_y_xxz_yzz[k] = -g_x_y_xx_yzz[k] * ab_z + g_x_y_xx_yzzz[k];

                g_x_y_xxz_zzz[k] = -g_x_y_xx_zzz[k] * ab_z + g_x_y_xx_zzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyy_xxx = cbuffer.data(ff_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xyy_xxy = cbuffer.data(ff_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xyy_xxz = cbuffer.data(ff_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xyy_xyy = cbuffer.data(ff_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xyy_xyz = cbuffer.data(ff_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xyy_xzz = cbuffer.data(ff_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xyy_yyy = cbuffer.data(ff_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xyy_yyz = cbuffer.data(ff_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_xyy_yzz = cbuffer.data(ff_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xyy_zzz = cbuffer.data(ff_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxx, g_0_y_yy_xxy, g_0_y_yy_xxz, g_0_y_yy_xyy, g_0_y_yy_xyz, g_0_y_yy_xzz, g_0_y_yy_yyy, g_0_y_yy_yyz, g_0_y_yy_yzz, g_0_y_yy_zzz, g_x_y_xyy_xxx, g_x_y_xyy_xxy, g_x_y_xyy_xxz, g_x_y_xyy_xyy, g_x_y_xyy_xyz, g_x_y_xyy_xzz, g_x_y_xyy_yyy, g_x_y_xyy_yyz, g_x_y_xyy_yzz, g_x_y_xyy_zzz, g_x_y_yy_xxx, g_x_y_yy_xxxx, g_x_y_yy_xxxy, g_x_y_yy_xxxz, g_x_y_yy_xxy, g_x_y_yy_xxyy, g_x_y_yy_xxyz, g_x_y_yy_xxz, g_x_y_yy_xxzz, g_x_y_yy_xyy, g_x_y_yy_xyyy, g_x_y_yy_xyyz, g_x_y_yy_xyz, g_x_y_yy_xyzz, g_x_y_yy_xzz, g_x_y_yy_xzzz, g_x_y_yy_yyy, g_x_y_yy_yyz, g_x_y_yy_yzz, g_x_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyy_xxx[k] = -g_0_y_yy_xxx[k] - g_x_y_yy_xxx[k] * ab_x + g_x_y_yy_xxxx[k];

                g_x_y_xyy_xxy[k] = -g_0_y_yy_xxy[k] - g_x_y_yy_xxy[k] * ab_x + g_x_y_yy_xxxy[k];

                g_x_y_xyy_xxz[k] = -g_0_y_yy_xxz[k] - g_x_y_yy_xxz[k] * ab_x + g_x_y_yy_xxxz[k];

                g_x_y_xyy_xyy[k] = -g_0_y_yy_xyy[k] - g_x_y_yy_xyy[k] * ab_x + g_x_y_yy_xxyy[k];

                g_x_y_xyy_xyz[k] = -g_0_y_yy_xyz[k] - g_x_y_yy_xyz[k] * ab_x + g_x_y_yy_xxyz[k];

                g_x_y_xyy_xzz[k] = -g_0_y_yy_xzz[k] - g_x_y_yy_xzz[k] * ab_x + g_x_y_yy_xxzz[k];

                g_x_y_xyy_yyy[k] = -g_0_y_yy_yyy[k] - g_x_y_yy_yyy[k] * ab_x + g_x_y_yy_xyyy[k];

                g_x_y_xyy_yyz[k] = -g_0_y_yy_yyz[k] - g_x_y_yy_yyz[k] * ab_x + g_x_y_yy_xyyz[k];

                g_x_y_xyy_yzz[k] = -g_0_y_yy_yzz[k] - g_x_y_yy_yzz[k] * ab_x + g_x_y_yy_xyzz[k];

                g_x_y_xyy_zzz[k] = -g_0_y_yy_zzz[k] - g_x_y_yy_zzz[k] * ab_x + g_x_y_yy_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyz_xxx = cbuffer.data(ff_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xyz_xxy = cbuffer.data(ff_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xyz_xxz = cbuffer.data(ff_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xyz_xyy = cbuffer.data(ff_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_xyz_xyz = cbuffer.data(ff_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xyz_xzz = cbuffer.data(ff_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xyz_yyy = cbuffer.data(ff_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_xyz_yyz = cbuffer.data(ff_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xyz_yzz = cbuffer.data(ff_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xyz_zzz = cbuffer.data(ff_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xy_xxx, g_x_y_xy_xxxz, g_x_y_xy_xxy, g_x_y_xy_xxyz, g_x_y_xy_xxz, g_x_y_xy_xxzz, g_x_y_xy_xyy, g_x_y_xy_xyyz, g_x_y_xy_xyz, g_x_y_xy_xyzz, g_x_y_xy_xzz, g_x_y_xy_xzzz, g_x_y_xy_yyy, g_x_y_xy_yyyz, g_x_y_xy_yyz, g_x_y_xy_yyzz, g_x_y_xy_yzz, g_x_y_xy_yzzz, g_x_y_xy_zzz, g_x_y_xy_zzzz, g_x_y_xyz_xxx, g_x_y_xyz_xxy, g_x_y_xyz_xxz, g_x_y_xyz_xyy, g_x_y_xyz_xyz, g_x_y_xyz_xzz, g_x_y_xyz_yyy, g_x_y_xyz_yyz, g_x_y_xyz_yzz, g_x_y_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyz_xxx[k] = -g_x_y_xy_xxx[k] * ab_z + g_x_y_xy_xxxz[k];

                g_x_y_xyz_xxy[k] = -g_x_y_xy_xxy[k] * ab_z + g_x_y_xy_xxyz[k];

                g_x_y_xyz_xxz[k] = -g_x_y_xy_xxz[k] * ab_z + g_x_y_xy_xxzz[k];

                g_x_y_xyz_xyy[k] = -g_x_y_xy_xyy[k] * ab_z + g_x_y_xy_xyyz[k];

                g_x_y_xyz_xyz[k] = -g_x_y_xy_xyz[k] * ab_z + g_x_y_xy_xyzz[k];

                g_x_y_xyz_xzz[k] = -g_x_y_xy_xzz[k] * ab_z + g_x_y_xy_xzzz[k];

                g_x_y_xyz_yyy[k] = -g_x_y_xy_yyy[k] * ab_z + g_x_y_xy_yyyz[k];

                g_x_y_xyz_yyz[k] = -g_x_y_xy_yyz[k] * ab_z + g_x_y_xy_yyzz[k];

                g_x_y_xyz_yzz[k] = -g_x_y_xy_yzz[k] * ab_z + g_x_y_xy_yzzz[k];

                g_x_y_xyz_zzz[k] = -g_x_y_xy_zzz[k] * ab_z + g_x_y_xy_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzz_xxx = cbuffer.data(ff_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xzz_xxy = cbuffer.data(ff_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xzz_xxz = cbuffer.data(ff_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xzz_xyy = cbuffer.data(ff_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xzz_xyz = cbuffer.data(ff_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xzz_xzz = cbuffer.data(ff_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_xzz_yyy = cbuffer.data(ff_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xzz_yyz = cbuffer.data(ff_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xzz_yzz = cbuffer.data(ff_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xzz_zzz = cbuffer.data(ff_geom_11_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xz_xxx, g_x_y_xz_xxxz, g_x_y_xz_xxy, g_x_y_xz_xxyz, g_x_y_xz_xxz, g_x_y_xz_xxzz, g_x_y_xz_xyy, g_x_y_xz_xyyz, g_x_y_xz_xyz, g_x_y_xz_xyzz, g_x_y_xz_xzz, g_x_y_xz_xzzz, g_x_y_xz_yyy, g_x_y_xz_yyyz, g_x_y_xz_yyz, g_x_y_xz_yyzz, g_x_y_xz_yzz, g_x_y_xz_yzzz, g_x_y_xz_zzz, g_x_y_xz_zzzz, g_x_y_xzz_xxx, g_x_y_xzz_xxy, g_x_y_xzz_xxz, g_x_y_xzz_xyy, g_x_y_xzz_xyz, g_x_y_xzz_xzz, g_x_y_xzz_yyy, g_x_y_xzz_yyz, g_x_y_xzz_yzz, g_x_y_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzz_xxx[k] = -g_x_y_xz_xxx[k] * ab_z + g_x_y_xz_xxxz[k];

                g_x_y_xzz_xxy[k] = -g_x_y_xz_xxy[k] * ab_z + g_x_y_xz_xxyz[k];

                g_x_y_xzz_xxz[k] = -g_x_y_xz_xxz[k] * ab_z + g_x_y_xz_xxzz[k];

                g_x_y_xzz_xyy[k] = -g_x_y_xz_xyy[k] * ab_z + g_x_y_xz_xyyz[k];

                g_x_y_xzz_xyz[k] = -g_x_y_xz_xyz[k] * ab_z + g_x_y_xz_xyzz[k];

                g_x_y_xzz_xzz[k] = -g_x_y_xz_xzz[k] * ab_z + g_x_y_xz_xzzz[k];

                g_x_y_xzz_yyy[k] = -g_x_y_xz_yyy[k] * ab_z + g_x_y_xz_yyyz[k];

                g_x_y_xzz_yyz[k] = -g_x_y_xz_yyz[k] * ab_z + g_x_y_xz_yyzz[k];

                g_x_y_xzz_yzz[k] = -g_x_y_xz_yzz[k] * ab_z + g_x_y_xz_yzzz[k];

                g_x_y_xzz_zzz[k] = -g_x_y_xz_zzz[k] * ab_z + g_x_y_xz_zzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyy_xxx = cbuffer.data(ff_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_yyy_xxy = cbuffer.data(ff_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_yyy_xxz = cbuffer.data(ff_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_yyy_xyy = cbuffer.data(ff_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_yyy_xyz = cbuffer.data(ff_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_yyy_xzz = cbuffer.data(ff_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_yyy_yyy = cbuffer.data(ff_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_yyy_yyz = cbuffer.data(ff_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_yyy_yzz = cbuffer.data(ff_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_yyy_zzz = cbuffer.data(ff_geom_11_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xxx, g_x_0_yy_xxy, g_x_0_yy_xxz, g_x_0_yy_xyy, g_x_0_yy_xyz, g_x_0_yy_xzz, g_x_0_yy_yyy, g_x_0_yy_yyz, g_x_0_yy_yzz, g_x_0_yy_zzz, g_x_y_yy_xxx, g_x_y_yy_xxxy, g_x_y_yy_xxy, g_x_y_yy_xxyy, g_x_y_yy_xxyz, g_x_y_yy_xxz, g_x_y_yy_xyy, g_x_y_yy_xyyy, g_x_y_yy_xyyz, g_x_y_yy_xyz, g_x_y_yy_xyzz, g_x_y_yy_xzz, g_x_y_yy_yyy, g_x_y_yy_yyyy, g_x_y_yy_yyyz, g_x_y_yy_yyz, g_x_y_yy_yyzz, g_x_y_yy_yzz, g_x_y_yy_yzzz, g_x_y_yy_zzz, g_x_y_yyy_xxx, g_x_y_yyy_xxy, g_x_y_yyy_xxz, g_x_y_yyy_xyy, g_x_y_yyy_xyz, g_x_y_yyy_xzz, g_x_y_yyy_yyy, g_x_y_yyy_yyz, g_x_y_yyy_yzz, g_x_y_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyy_xxx[k] = g_x_0_yy_xxx[k] - g_x_y_yy_xxx[k] * ab_y + g_x_y_yy_xxxy[k];

                g_x_y_yyy_xxy[k] = g_x_0_yy_xxy[k] - g_x_y_yy_xxy[k] * ab_y + g_x_y_yy_xxyy[k];

                g_x_y_yyy_xxz[k] = g_x_0_yy_xxz[k] - g_x_y_yy_xxz[k] * ab_y + g_x_y_yy_xxyz[k];

                g_x_y_yyy_xyy[k] = g_x_0_yy_xyy[k] - g_x_y_yy_xyy[k] * ab_y + g_x_y_yy_xyyy[k];

                g_x_y_yyy_xyz[k] = g_x_0_yy_xyz[k] - g_x_y_yy_xyz[k] * ab_y + g_x_y_yy_xyyz[k];

                g_x_y_yyy_xzz[k] = g_x_0_yy_xzz[k] - g_x_y_yy_xzz[k] * ab_y + g_x_y_yy_xyzz[k];

                g_x_y_yyy_yyy[k] = g_x_0_yy_yyy[k] - g_x_y_yy_yyy[k] * ab_y + g_x_y_yy_yyyy[k];

                g_x_y_yyy_yyz[k] = g_x_0_yy_yyz[k] - g_x_y_yy_yyz[k] * ab_y + g_x_y_yy_yyyz[k];

                g_x_y_yyy_yzz[k] = g_x_0_yy_yzz[k] - g_x_y_yy_yzz[k] * ab_y + g_x_y_yy_yyzz[k];

                g_x_y_yyy_zzz[k] = g_x_0_yy_zzz[k] - g_x_y_yy_zzz[k] * ab_y + g_x_y_yy_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyz_xxx = cbuffer.data(ff_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_yyz_xxy = cbuffer.data(ff_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_yyz_xxz = cbuffer.data(ff_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_yyz_xyy = cbuffer.data(ff_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_yyz_xyz = cbuffer.data(ff_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_yyz_xzz = cbuffer.data(ff_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_yyz_yyy = cbuffer.data(ff_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_yyz_yyz = cbuffer.data(ff_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_yyz_yzz = cbuffer.data(ff_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_yyz_zzz = cbuffer.data(ff_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yy_xxx, g_x_y_yy_xxxz, g_x_y_yy_xxy, g_x_y_yy_xxyz, g_x_y_yy_xxz, g_x_y_yy_xxzz, g_x_y_yy_xyy, g_x_y_yy_xyyz, g_x_y_yy_xyz, g_x_y_yy_xyzz, g_x_y_yy_xzz, g_x_y_yy_xzzz, g_x_y_yy_yyy, g_x_y_yy_yyyz, g_x_y_yy_yyz, g_x_y_yy_yyzz, g_x_y_yy_yzz, g_x_y_yy_yzzz, g_x_y_yy_zzz, g_x_y_yy_zzzz, g_x_y_yyz_xxx, g_x_y_yyz_xxy, g_x_y_yyz_xxz, g_x_y_yyz_xyy, g_x_y_yyz_xyz, g_x_y_yyz_xzz, g_x_y_yyz_yyy, g_x_y_yyz_yyz, g_x_y_yyz_yzz, g_x_y_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyz_xxx[k] = -g_x_y_yy_xxx[k] * ab_z + g_x_y_yy_xxxz[k];

                g_x_y_yyz_xxy[k] = -g_x_y_yy_xxy[k] * ab_z + g_x_y_yy_xxyz[k];

                g_x_y_yyz_xxz[k] = -g_x_y_yy_xxz[k] * ab_z + g_x_y_yy_xxzz[k];

                g_x_y_yyz_xyy[k] = -g_x_y_yy_xyy[k] * ab_z + g_x_y_yy_xyyz[k];

                g_x_y_yyz_xyz[k] = -g_x_y_yy_xyz[k] * ab_z + g_x_y_yy_xyzz[k];

                g_x_y_yyz_xzz[k] = -g_x_y_yy_xzz[k] * ab_z + g_x_y_yy_xzzz[k];

                g_x_y_yyz_yyy[k] = -g_x_y_yy_yyy[k] * ab_z + g_x_y_yy_yyyz[k];

                g_x_y_yyz_yyz[k] = -g_x_y_yy_yyz[k] * ab_z + g_x_y_yy_yyzz[k];

                g_x_y_yyz_yzz[k] = -g_x_y_yy_yzz[k] * ab_z + g_x_y_yy_yzzz[k];

                g_x_y_yyz_zzz[k] = -g_x_y_yy_zzz[k] * ab_z + g_x_y_yy_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzz_xxx = cbuffer.data(ff_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_yzz_xxy = cbuffer.data(ff_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_yzz_xxz = cbuffer.data(ff_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_yzz_xyy = cbuffer.data(ff_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_yzz_xyz = cbuffer.data(ff_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_yzz_xzz = cbuffer.data(ff_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_yzz_yyy = cbuffer.data(ff_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_yzz_yyz = cbuffer.data(ff_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_yzz_yzz = cbuffer.data(ff_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_yzz_zzz = cbuffer.data(ff_geom_11_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yz_xxx, g_x_y_yz_xxxz, g_x_y_yz_xxy, g_x_y_yz_xxyz, g_x_y_yz_xxz, g_x_y_yz_xxzz, g_x_y_yz_xyy, g_x_y_yz_xyyz, g_x_y_yz_xyz, g_x_y_yz_xyzz, g_x_y_yz_xzz, g_x_y_yz_xzzz, g_x_y_yz_yyy, g_x_y_yz_yyyz, g_x_y_yz_yyz, g_x_y_yz_yyzz, g_x_y_yz_yzz, g_x_y_yz_yzzz, g_x_y_yz_zzz, g_x_y_yz_zzzz, g_x_y_yzz_xxx, g_x_y_yzz_xxy, g_x_y_yzz_xxz, g_x_y_yzz_xyy, g_x_y_yzz_xyz, g_x_y_yzz_xzz, g_x_y_yzz_yyy, g_x_y_yzz_yyz, g_x_y_yzz_yzz, g_x_y_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzz_xxx[k] = -g_x_y_yz_xxx[k] * ab_z + g_x_y_yz_xxxz[k];

                g_x_y_yzz_xxy[k] = -g_x_y_yz_xxy[k] * ab_z + g_x_y_yz_xxyz[k];

                g_x_y_yzz_xxz[k] = -g_x_y_yz_xxz[k] * ab_z + g_x_y_yz_xxzz[k];

                g_x_y_yzz_xyy[k] = -g_x_y_yz_xyy[k] * ab_z + g_x_y_yz_xyyz[k];

                g_x_y_yzz_xyz[k] = -g_x_y_yz_xyz[k] * ab_z + g_x_y_yz_xyzz[k];

                g_x_y_yzz_xzz[k] = -g_x_y_yz_xzz[k] * ab_z + g_x_y_yz_xzzz[k];

                g_x_y_yzz_yyy[k] = -g_x_y_yz_yyy[k] * ab_z + g_x_y_yz_yyyz[k];

                g_x_y_yzz_yyz[k] = -g_x_y_yz_yyz[k] * ab_z + g_x_y_yz_yyzz[k];

                g_x_y_yzz_yzz[k] = -g_x_y_yz_yzz[k] * ab_z + g_x_y_yz_yzzz[k];

                g_x_y_yzz_zzz[k] = -g_x_y_yz_zzz[k] * ab_z + g_x_y_yz_zzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzz_xxx = cbuffer.data(ff_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_zzz_xxy = cbuffer.data(ff_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_zzz_xxz = cbuffer.data(ff_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_zzz_xyy = cbuffer.data(ff_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_zzz_xyz = cbuffer.data(ff_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_zzz_xzz = cbuffer.data(ff_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_zzz_yyy = cbuffer.data(ff_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_zzz_yyz = cbuffer.data(ff_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_zzz_yzz = cbuffer.data(ff_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_zzz_zzz = cbuffer.data(ff_geom_11_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zz_xxx, g_x_y_zz_xxxz, g_x_y_zz_xxy, g_x_y_zz_xxyz, g_x_y_zz_xxz, g_x_y_zz_xxzz, g_x_y_zz_xyy, g_x_y_zz_xyyz, g_x_y_zz_xyz, g_x_y_zz_xyzz, g_x_y_zz_xzz, g_x_y_zz_xzzz, g_x_y_zz_yyy, g_x_y_zz_yyyz, g_x_y_zz_yyz, g_x_y_zz_yyzz, g_x_y_zz_yzz, g_x_y_zz_yzzz, g_x_y_zz_zzz, g_x_y_zz_zzzz, g_x_y_zzz_xxx, g_x_y_zzz_xxy, g_x_y_zzz_xxz, g_x_y_zzz_xyy, g_x_y_zzz_xyz, g_x_y_zzz_xzz, g_x_y_zzz_yyy, g_x_y_zzz_yyz, g_x_y_zzz_yzz, g_x_y_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzz_xxx[k] = -g_x_y_zz_xxx[k] * ab_z + g_x_y_zz_xxxz[k];

                g_x_y_zzz_xxy[k] = -g_x_y_zz_xxy[k] * ab_z + g_x_y_zz_xxyz[k];

                g_x_y_zzz_xxz[k] = -g_x_y_zz_xxz[k] * ab_z + g_x_y_zz_xxzz[k];

                g_x_y_zzz_xyy[k] = -g_x_y_zz_xyy[k] * ab_z + g_x_y_zz_xyyz[k];

                g_x_y_zzz_xyz[k] = -g_x_y_zz_xyz[k] * ab_z + g_x_y_zz_xyzz[k];

                g_x_y_zzz_xzz[k] = -g_x_y_zz_xzz[k] * ab_z + g_x_y_zz_xzzz[k];

                g_x_y_zzz_yyy[k] = -g_x_y_zz_yyy[k] * ab_z + g_x_y_zz_yyyz[k];

                g_x_y_zzz_yyz[k] = -g_x_y_zz_yyz[k] * ab_z + g_x_y_zz_yyzz[k];

                g_x_y_zzz_yzz[k] = -g_x_y_zz_yzz[k] * ab_z + g_x_y_zz_yzzz[k];

                g_x_y_zzz_zzz[k] = -g_x_y_zz_zzz[k] * ab_z + g_x_y_zz_zzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxx_xxx = cbuffer.data(ff_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_xxx_xxy = cbuffer.data(ff_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_xxx_xxz = cbuffer.data(ff_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_xxx_xyy = cbuffer.data(ff_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_z_xxx_xyz = cbuffer.data(ff_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_xxx_xzz = cbuffer.data(ff_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_xxx_yyy = cbuffer.data(ff_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_xxx_yyz = cbuffer.data(ff_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_xxx_yzz = cbuffer.data(ff_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_xxx_zzz = cbuffer.data(ff_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_xxx, g_0_z_xx_xxy, g_0_z_xx_xxz, g_0_z_xx_xyy, g_0_z_xx_xyz, g_0_z_xx_xzz, g_0_z_xx_yyy, g_0_z_xx_yyz, g_0_z_xx_yzz, g_0_z_xx_zzz, g_x_z_xx_xxx, g_x_z_xx_xxxx, g_x_z_xx_xxxy, g_x_z_xx_xxxz, g_x_z_xx_xxy, g_x_z_xx_xxyy, g_x_z_xx_xxyz, g_x_z_xx_xxz, g_x_z_xx_xxzz, g_x_z_xx_xyy, g_x_z_xx_xyyy, g_x_z_xx_xyyz, g_x_z_xx_xyz, g_x_z_xx_xyzz, g_x_z_xx_xzz, g_x_z_xx_xzzz, g_x_z_xx_yyy, g_x_z_xx_yyz, g_x_z_xx_yzz, g_x_z_xx_zzz, g_x_z_xxx_xxx, g_x_z_xxx_xxy, g_x_z_xxx_xxz, g_x_z_xxx_xyy, g_x_z_xxx_xyz, g_x_z_xxx_xzz, g_x_z_xxx_yyy, g_x_z_xxx_yyz, g_x_z_xxx_yzz, g_x_z_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxx_xxx[k] = -g_0_z_xx_xxx[k] - g_x_z_xx_xxx[k] * ab_x + g_x_z_xx_xxxx[k];

                g_x_z_xxx_xxy[k] = -g_0_z_xx_xxy[k] - g_x_z_xx_xxy[k] * ab_x + g_x_z_xx_xxxy[k];

                g_x_z_xxx_xxz[k] = -g_0_z_xx_xxz[k] - g_x_z_xx_xxz[k] * ab_x + g_x_z_xx_xxxz[k];

                g_x_z_xxx_xyy[k] = -g_0_z_xx_xyy[k] - g_x_z_xx_xyy[k] * ab_x + g_x_z_xx_xxyy[k];

                g_x_z_xxx_xyz[k] = -g_0_z_xx_xyz[k] - g_x_z_xx_xyz[k] * ab_x + g_x_z_xx_xxyz[k];

                g_x_z_xxx_xzz[k] = -g_0_z_xx_xzz[k] - g_x_z_xx_xzz[k] * ab_x + g_x_z_xx_xxzz[k];

                g_x_z_xxx_yyy[k] = -g_0_z_xx_yyy[k] - g_x_z_xx_yyy[k] * ab_x + g_x_z_xx_xyyy[k];

                g_x_z_xxx_yyz[k] = -g_0_z_xx_yyz[k] - g_x_z_xx_yyz[k] * ab_x + g_x_z_xx_xyyz[k];

                g_x_z_xxx_yzz[k] = -g_0_z_xx_yzz[k] - g_x_z_xx_yzz[k] * ab_x + g_x_z_xx_xyzz[k];

                g_x_z_xxx_zzz[k] = -g_0_z_xx_zzz[k] - g_x_z_xx_zzz[k] * ab_x + g_x_z_xx_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxy_xxx = cbuffer.data(ff_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_xxy_xxy = cbuffer.data(ff_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_xxy_xxz = cbuffer.data(ff_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_xxy_xyy = cbuffer.data(ff_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_xxy_xyz = cbuffer.data(ff_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_xxy_xzz = cbuffer.data(ff_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_xxy_yyy = cbuffer.data(ff_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_xxy_yyz = cbuffer.data(ff_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_xxy_yzz = cbuffer.data(ff_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_xxy_zzz = cbuffer.data(ff_geom_11_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xx_xxx, g_x_z_xx_xxxy, g_x_z_xx_xxy, g_x_z_xx_xxyy, g_x_z_xx_xxyz, g_x_z_xx_xxz, g_x_z_xx_xyy, g_x_z_xx_xyyy, g_x_z_xx_xyyz, g_x_z_xx_xyz, g_x_z_xx_xyzz, g_x_z_xx_xzz, g_x_z_xx_yyy, g_x_z_xx_yyyy, g_x_z_xx_yyyz, g_x_z_xx_yyz, g_x_z_xx_yyzz, g_x_z_xx_yzz, g_x_z_xx_yzzz, g_x_z_xx_zzz, g_x_z_xxy_xxx, g_x_z_xxy_xxy, g_x_z_xxy_xxz, g_x_z_xxy_xyy, g_x_z_xxy_xyz, g_x_z_xxy_xzz, g_x_z_xxy_yyy, g_x_z_xxy_yyz, g_x_z_xxy_yzz, g_x_z_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxy_xxx[k] = -g_x_z_xx_xxx[k] * ab_y + g_x_z_xx_xxxy[k];

                g_x_z_xxy_xxy[k] = -g_x_z_xx_xxy[k] * ab_y + g_x_z_xx_xxyy[k];

                g_x_z_xxy_xxz[k] = -g_x_z_xx_xxz[k] * ab_y + g_x_z_xx_xxyz[k];

                g_x_z_xxy_xyy[k] = -g_x_z_xx_xyy[k] * ab_y + g_x_z_xx_xyyy[k];

                g_x_z_xxy_xyz[k] = -g_x_z_xx_xyz[k] * ab_y + g_x_z_xx_xyyz[k];

                g_x_z_xxy_xzz[k] = -g_x_z_xx_xzz[k] * ab_y + g_x_z_xx_xyzz[k];

                g_x_z_xxy_yyy[k] = -g_x_z_xx_yyy[k] * ab_y + g_x_z_xx_yyyy[k];

                g_x_z_xxy_yyz[k] = -g_x_z_xx_yyz[k] * ab_y + g_x_z_xx_yyyz[k];

                g_x_z_xxy_yzz[k] = -g_x_z_xx_yzz[k] * ab_y + g_x_z_xx_yyzz[k];

                g_x_z_xxy_zzz[k] = -g_x_z_xx_zzz[k] * ab_y + g_x_z_xx_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxz_xxx = cbuffer.data(ff_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_xxz_xxy = cbuffer.data(ff_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_xxz_xxz = cbuffer.data(ff_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_xxz_xyy = cbuffer.data(ff_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_xxz_xyz = cbuffer.data(ff_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_xxz_xzz = cbuffer.data(ff_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_xxz_yyy = cbuffer.data(ff_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_xxz_yyz = cbuffer.data(ff_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_xxz_yzz = cbuffer.data(ff_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_xxz_zzz = cbuffer.data(ff_geom_11_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_xxx, g_0_z_xz_xxy, g_0_z_xz_xxz, g_0_z_xz_xyy, g_0_z_xz_xyz, g_0_z_xz_xzz, g_0_z_xz_yyy, g_0_z_xz_yyz, g_0_z_xz_yzz, g_0_z_xz_zzz, g_x_z_xxz_xxx, g_x_z_xxz_xxy, g_x_z_xxz_xxz, g_x_z_xxz_xyy, g_x_z_xxz_xyz, g_x_z_xxz_xzz, g_x_z_xxz_yyy, g_x_z_xxz_yyz, g_x_z_xxz_yzz, g_x_z_xxz_zzz, g_x_z_xz_xxx, g_x_z_xz_xxxx, g_x_z_xz_xxxy, g_x_z_xz_xxxz, g_x_z_xz_xxy, g_x_z_xz_xxyy, g_x_z_xz_xxyz, g_x_z_xz_xxz, g_x_z_xz_xxzz, g_x_z_xz_xyy, g_x_z_xz_xyyy, g_x_z_xz_xyyz, g_x_z_xz_xyz, g_x_z_xz_xyzz, g_x_z_xz_xzz, g_x_z_xz_xzzz, g_x_z_xz_yyy, g_x_z_xz_yyz, g_x_z_xz_yzz, g_x_z_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxz_xxx[k] = -g_0_z_xz_xxx[k] - g_x_z_xz_xxx[k] * ab_x + g_x_z_xz_xxxx[k];

                g_x_z_xxz_xxy[k] = -g_0_z_xz_xxy[k] - g_x_z_xz_xxy[k] * ab_x + g_x_z_xz_xxxy[k];

                g_x_z_xxz_xxz[k] = -g_0_z_xz_xxz[k] - g_x_z_xz_xxz[k] * ab_x + g_x_z_xz_xxxz[k];

                g_x_z_xxz_xyy[k] = -g_0_z_xz_xyy[k] - g_x_z_xz_xyy[k] * ab_x + g_x_z_xz_xxyy[k];

                g_x_z_xxz_xyz[k] = -g_0_z_xz_xyz[k] - g_x_z_xz_xyz[k] * ab_x + g_x_z_xz_xxyz[k];

                g_x_z_xxz_xzz[k] = -g_0_z_xz_xzz[k] - g_x_z_xz_xzz[k] * ab_x + g_x_z_xz_xxzz[k];

                g_x_z_xxz_yyy[k] = -g_0_z_xz_yyy[k] - g_x_z_xz_yyy[k] * ab_x + g_x_z_xz_xyyy[k];

                g_x_z_xxz_yyz[k] = -g_0_z_xz_yyz[k] - g_x_z_xz_yyz[k] * ab_x + g_x_z_xz_xyyz[k];

                g_x_z_xxz_yzz[k] = -g_0_z_xz_yzz[k] - g_x_z_xz_yzz[k] * ab_x + g_x_z_xz_xyzz[k];

                g_x_z_xxz_zzz[k] = -g_0_z_xz_zzz[k] - g_x_z_xz_zzz[k] * ab_x + g_x_z_xz_xzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyy_xxx = cbuffer.data(ff_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_xyy_xxy = cbuffer.data(ff_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_xyy_xxz = cbuffer.data(ff_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_xyy_xyy = cbuffer.data(ff_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_xyy_xyz = cbuffer.data(ff_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_xyy_xzz = cbuffer.data(ff_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_xyy_yyy = cbuffer.data(ff_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_xyy_yyz = cbuffer.data(ff_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_xyy_yzz = cbuffer.data(ff_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_xyy_zzz = cbuffer.data(ff_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xy_xxx, g_x_z_xy_xxxy, g_x_z_xy_xxy, g_x_z_xy_xxyy, g_x_z_xy_xxyz, g_x_z_xy_xxz, g_x_z_xy_xyy, g_x_z_xy_xyyy, g_x_z_xy_xyyz, g_x_z_xy_xyz, g_x_z_xy_xyzz, g_x_z_xy_xzz, g_x_z_xy_yyy, g_x_z_xy_yyyy, g_x_z_xy_yyyz, g_x_z_xy_yyz, g_x_z_xy_yyzz, g_x_z_xy_yzz, g_x_z_xy_yzzz, g_x_z_xy_zzz, g_x_z_xyy_xxx, g_x_z_xyy_xxy, g_x_z_xyy_xxz, g_x_z_xyy_xyy, g_x_z_xyy_xyz, g_x_z_xyy_xzz, g_x_z_xyy_yyy, g_x_z_xyy_yyz, g_x_z_xyy_yzz, g_x_z_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyy_xxx[k] = -g_x_z_xy_xxx[k] * ab_y + g_x_z_xy_xxxy[k];

                g_x_z_xyy_xxy[k] = -g_x_z_xy_xxy[k] * ab_y + g_x_z_xy_xxyy[k];

                g_x_z_xyy_xxz[k] = -g_x_z_xy_xxz[k] * ab_y + g_x_z_xy_xxyz[k];

                g_x_z_xyy_xyy[k] = -g_x_z_xy_xyy[k] * ab_y + g_x_z_xy_xyyy[k];

                g_x_z_xyy_xyz[k] = -g_x_z_xy_xyz[k] * ab_y + g_x_z_xy_xyyz[k];

                g_x_z_xyy_xzz[k] = -g_x_z_xy_xzz[k] * ab_y + g_x_z_xy_xyzz[k];

                g_x_z_xyy_yyy[k] = -g_x_z_xy_yyy[k] * ab_y + g_x_z_xy_yyyy[k];

                g_x_z_xyy_yyz[k] = -g_x_z_xy_yyz[k] * ab_y + g_x_z_xy_yyyz[k];

                g_x_z_xyy_yzz[k] = -g_x_z_xy_yzz[k] * ab_y + g_x_z_xy_yyzz[k];

                g_x_z_xyy_zzz[k] = -g_x_z_xy_zzz[k] * ab_y + g_x_z_xy_yzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyz_xxx = cbuffer.data(ff_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_xyz_xxy = cbuffer.data(ff_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_xyz_xxz = cbuffer.data(ff_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_xyz_xyy = cbuffer.data(ff_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_xyz_xyz = cbuffer.data(ff_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_xyz_xzz = cbuffer.data(ff_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_xyz_yyy = cbuffer.data(ff_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_xyz_yyz = cbuffer.data(ff_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_xyz_yzz = cbuffer.data(ff_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_xyz_zzz = cbuffer.data(ff_geom_11_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyz_xxx, g_x_z_xyz_xxy, g_x_z_xyz_xxz, g_x_z_xyz_xyy, g_x_z_xyz_xyz, g_x_z_xyz_xzz, g_x_z_xyz_yyy, g_x_z_xyz_yyz, g_x_z_xyz_yzz, g_x_z_xyz_zzz, g_x_z_xz_xxx, g_x_z_xz_xxxy, g_x_z_xz_xxy, g_x_z_xz_xxyy, g_x_z_xz_xxyz, g_x_z_xz_xxz, g_x_z_xz_xyy, g_x_z_xz_xyyy, g_x_z_xz_xyyz, g_x_z_xz_xyz, g_x_z_xz_xyzz, g_x_z_xz_xzz, g_x_z_xz_yyy, g_x_z_xz_yyyy, g_x_z_xz_yyyz, g_x_z_xz_yyz, g_x_z_xz_yyzz, g_x_z_xz_yzz, g_x_z_xz_yzzz, g_x_z_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyz_xxx[k] = -g_x_z_xz_xxx[k] * ab_y + g_x_z_xz_xxxy[k];

                g_x_z_xyz_xxy[k] = -g_x_z_xz_xxy[k] * ab_y + g_x_z_xz_xxyy[k];

                g_x_z_xyz_xxz[k] = -g_x_z_xz_xxz[k] * ab_y + g_x_z_xz_xxyz[k];

                g_x_z_xyz_xyy[k] = -g_x_z_xz_xyy[k] * ab_y + g_x_z_xz_xyyy[k];

                g_x_z_xyz_xyz[k] = -g_x_z_xz_xyz[k] * ab_y + g_x_z_xz_xyyz[k];

                g_x_z_xyz_xzz[k] = -g_x_z_xz_xzz[k] * ab_y + g_x_z_xz_xyzz[k];

                g_x_z_xyz_yyy[k] = -g_x_z_xz_yyy[k] * ab_y + g_x_z_xz_yyyy[k];

                g_x_z_xyz_yyz[k] = -g_x_z_xz_yyz[k] * ab_y + g_x_z_xz_yyyz[k];

                g_x_z_xyz_yzz[k] = -g_x_z_xz_yzz[k] * ab_y + g_x_z_xz_yyzz[k];

                g_x_z_xyz_zzz[k] = -g_x_z_xz_zzz[k] * ab_y + g_x_z_xz_yzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzz_xxx = cbuffer.data(ff_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_xzz_xxy = cbuffer.data(ff_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_xzz_xxz = cbuffer.data(ff_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_xzz_xyy = cbuffer.data(ff_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_xzz_xyz = cbuffer.data(ff_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_xzz_xzz = cbuffer.data(ff_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_xzz_yyy = cbuffer.data(ff_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_xzz_yyz = cbuffer.data(ff_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_xzz_yzz = cbuffer.data(ff_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_xzz_zzz = cbuffer.data(ff_geom_11_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxx, g_0_z_zz_xxy, g_0_z_zz_xxz, g_0_z_zz_xyy, g_0_z_zz_xyz, g_0_z_zz_xzz, g_0_z_zz_yyy, g_0_z_zz_yyz, g_0_z_zz_yzz, g_0_z_zz_zzz, g_x_z_xzz_xxx, g_x_z_xzz_xxy, g_x_z_xzz_xxz, g_x_z_xzz_xyy, g_x_z_xzz_xyz, g_x_z_xzz_xzz, g_x_z_xzz_yyy, g_x_z_xzz_yyz, g_x_z_xzz_yzz, g_x_z_xzz_zzz, g_x_z_zz_xxx, g_x_z_zz_xxxx, g_x_z_zz_xxxy, g_x_z_zz_xxxz, g_x_z_zz_xxy, g_x_z_zz_xxyy, g_x_z_zz_xxyz, g_x_z_zz_xxz, g_x_z_zz_xxzz, g_x_z_zz_xyy, g_x_z_zz_xyyy, g_x_z_zz_xyyz, g_x_z_zz_xyz, g_x_z_zz_xyzz, g_x_z_zz_xzz, g_x_z_zz_xzzz, g_x_z_zz_yyy, g_x_z_zz_yyz, g_x_z_zz_yzz, g_x_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzz_xxx[k] = -g_0_z_zz_xxx[k] - g_x_z_zz_xxx[k] * ab_x + g_x_z_zz_xxxx[k];

                g_x_z_xzz_xxy[k] = -g_0_z_zz_xxy[k] - g_x_z_zz_xxy[k] * ab_x + g_x_z_zz_xxxy[k];

                g_x_z_xzz_xxz[k] = -g_0_z_zz_xxz[k] - g_x_z_zz_xxz[k] * ab_x + g_x_z_zz_xxxz[k];

                g_x_z_xzz_xyy[k] = -g_0_z_zz_xyy[k] - g_x_z_zz_xyy[k] * ab_x + g_x_z_zz_xxyy[k];

                g_x_z_xzz_xyz[k] = -g_0_z_zz_xyz[k] - g_x_z_zz_xyz[k] * ab_x + g_x_z_zz_xxyz[k];

                g_x_z_xzz_xzz[k] = -g_0_z_zz_xzz[k] - g_x_z_zz_xzz[k] * ab_x + g_x_z_zz_xxzz[k];

                g_x_z_xzz_yyy[k] = -g_0_z_zz_yyy[k] - g_x_z_zz_yyy[k] * ab_x + g_x_z_zz_xyyy[k];

                g_x_z_xzz_yyz[k] = -g_0_z_zz_yyz[k] - g_x_z_zz_yyz[k] * ab_x + g_x_z_zz_xyyz[k];

                g_x_z_xzz_yzz[k] = -g_0_z_zz_yzz[k] - g_x_z_zz_yzz[k] * ab_x + g_x_z_zz_xyzz[k];

                g_x_z_xzz_zzz[k] = -g_0_z_zz_zzz[k] - g_x_z_zz_zzz[k] * ab_x + g_x_z_zz_xzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyy_xxx = cbuffer.data(ff_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_yyy_xxy = cbuffer.data(ff_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_yyy_xxz = cbuffer.data(ff_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_yyy_xyy = cbuffer.data(ff_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_yyy_xyz = cbuffer.data(ff_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_yyy_xzz = cbuffer.data(ff_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_yyy_yyy = cbuffer.data(ff_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_yyy_yyz = cbuffer.data(ff_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_yyy_yzz = cbuffer.data(ff_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_yyy_zzz = cbuffer.data(ff_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yy_xxx, g_x_z_yy_xxxy, g_x_z_yy_xxy, g_x_z_yy_xxyy, g_x_z_yy_xxyz, g_x_z_yy_xxz, g_x_z_yy_xyy, g_x_z_yy_xyyy, g_x_z_yy_xyyz, g_x_z_yy_xyz, g_x_z_yy_xyzz, g_x_z_yy_xzz, g_x_z_yy_yyy, g_x_z_yy_yyyy, g_x_z_yy_yyyz, g_x_z_yy_yyz, g_x_z_yy_yyzz, g_x_z_yy_yzz, g_x_z_yy_yzzz, g_x_z_yy_zzz, g_x_z_yyy_xxx, g_x_z_yyy_xxy, g_x_z_yyy_xxz, g_x_z_yyy_xyy, g_x_z_yyy_xyz, g_x_z_yyy_xzz, g_x_z_yyy_yyy, g_x_z_yyy_yyz, g_x_z_yyy_yzz, g_x_z_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyy_xxx[k] = -g_x_z_yy_xxx[k] * ab_y + g_x_z_yy_xxxy[k];

                g_x_z_yyy_xxy[k] = -g_x_z_yy_xxy[k] * ab_y + g_x_z_yy_xxyy[k];

                g_x_z_yyy_xxz[k] = -g_x_z_yy_xxz[k] * ab_y + g_x_z_yy_xxyz[k];

                g_x_z_yyy_xyy[k] = -g_x_z_yy_xyy[k] * ab_y + g_x_z_yy_xyyy[k];

                g_x_z_yyy_xyz[k] = -g_x_z_yy_xyz[k] * ab_y + g_x_z_yy_xyyz[k];

                g_x_z_yyy_xzz[k] = -g_x_z_yy_xzz[k] * ab_y + g_x_z_yy_xyzz[k];

                g_x_z_yyy_yyy[k] = -g_x_z_yy_yyy[k] * ab_y + g_x_z_yy_yyyy[k];

                g_x_z_yyy_yyz[k] = -g_x_z_yy_yyz[k] * ab_y + g_x_z_yy_yyyz[k];

                g_x_z_yyy_yzz[k] = -g_x_z_yy_yzz[k] * ab_y + g_x_z_yy_yyzz[k];

                g_x_z_yyy_zzz[k] = -g_x_z_yy_zzz[k] * ab_y + g_x_z_yy_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyz_xxx = cbuffer.data(ff_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_z_yyz_xxy = cbuffer.data(ff_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_z_yyz_xxz = cbuffer.data(ff_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_z_yyz_xyy = cbuffer.data(ff_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_z_yyz_xyz = cbuffer.data(ff_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_z_yyz_xzz = cbuffer.data(ff_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_z_yyz_yyy = cbuffer.data(ff_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_z_yyz_yyz = cbuffer.data(ff_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_z_yyz_yzz = cbuffer.data(ff_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_z_yyz_zzz = cbuffer.data(ff_geom_11_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyz_xxx, g_x_z_yyz_xxy, g_x_z_yyz_xxz, g_x_z_yyz_xyy, g_x_z_yyz_xyz, g_x_z_yyz_xzz, g_x_z_yyz_yyy, g_x_z_yyz_yyz, g_x_z_yyz_yzz, g_x_z_yyz_zzz, g_x_z_yz_xxx, g_x_z_yz_xxxy, g_x_z_yz_xxy, g_x_z_yz_xxyy, g_x_z_yz_xxyz, g_x_z_yz_xxz, g_x_z_yz_xyy, g_x_z_yz_xyyy, g_x_z_yz_xyyz, g_x_z_yz_xyz, g_x_z_yz_xyzz, g_x_z_yz_xzz, g_x_z_yz_yyy, g_x_z_yz_yyyy, g_x_z_yz_yyyz, g_x_z_yz_yyz, g_x_z_yz_yyzz, g_x_z_yz_yzz, g_x_z_yz_yzzz, g_x_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyz_xxx[k] = -g_x_z_yz_xxx[k] * ab_y + g_x_z_yz_xxxy[k];

                g_x_z_yyz_xxy[k] = -g_x_z_yz_xxy[k] * ab_y + g_x_z_yz_xxyy[k];

                g_x_z_yyz_xxz[k] = -g_x_z_yz_xxz[k] * ab_y + g_x_z_yz_xxyz[k];

                g_x_z_yyz_xyy[k] = -g_x_z_yz_xyy[k] * ab_y + g_x_z_yz_xyyy[k];

                g_x_z_yyz_xyz[k] = -g_x_z_yz_xyz[k] * ab_y + g_x_z_yz_xyyz[k];

                g_x_z_yyz_xzz[k] = -g_x_z_yz_xzz[k] * ab_y + g_x_z_yz_xyzz[k];

                g_x_z_yyz_yyy[k] = -g_x_z_yz_yyy[k] * ab_y + g_x_z_yz_yyyy[k];

                g_x_z_yyz_yyz[k] = -g_x_z_yz_yyz[k] * ab_y + g_x_z_yz_yyyz[k];

                g_x_z_yyz_yzz[k] = -g_x_z_yz_yzz[k] * ab_y + g_x_z_yz_yyzz[k];

                g_x_z_yyz_zzz[k] = -g_x_z_yz_zzz[k] * ab_y + g_x_z_yz_yzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzz_xxx = cbuffer.data(ff_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_z_yzz_xxy = cbuffer.data(ff_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_z_yzz_xxz = cbuffer.data(ff_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_z_yzz_xyy = cbuffer.data(ff_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_z_yzz_xyz = cbuffer.data(ff_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_z_yzz_xzz = cbuffer.data(ff_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_z_yzz_yyy = cbuffer.data(ff_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_z_yzz_yyz = cbuffer.data(ff_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_z_yzz_yzz = cbuffer.data(ff_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_z_yzz_zzz = cbuffer.data(ff_geom_11_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzz_xxx, g_x_z_yzz_xxy, g_x_z_yzz_xxz, g_x_z_yzz_xyy, g_x_z_yzz_xyz, g_x_z_yzz_xzz, g_x_z_yzz_yyy, g_x_z_yzz_yyz, g_x_z_yzz_yzz, g_x_z_yzz_zzz, g_x_z_zz_xxx, g_x_z_zz_xxxy, g_x_z_zz_xxy, g_x_z_zz_xxyy, g_x_z_zz_xxyz, g_x_z_zz_xxz, g_x_z_zz_xyy, g_x_z_zz_xyyy, g_x_z_zz_xyyz, g_x_z_zz_xyz, g_x_z_zz_xyzz, g_x_z_zz_xzz, g_x_z_zz_yyy, g_x_z_zz_yyyy, g_x_z_zz_yyyz, g_x_z_zz_yyz, g_x_z_zz_yyzz, g_x_z_zz_yzz, g_x_z_zz_yzzz, g_x_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzz_xxx[k] = -g_x_z_zz_xxx[k] * ab_y + g_x_z_zz_xxxy[k];

                g_x_z_yzz_xxy[k] = -g_x_z_zz_xxy[k] * ab_y + g_x_z_zz_xxyy[k];

                g_x_z_yzz_xxz[k] = -g_x_z_zz_xxz[k] * ab_y + g_x_z_zz_xxyz[k];

                g_x_z_yzz_xyy[k] = -g_x_z_zz_xyy[k] * ab_y + g_x_z_zz_xyyy[k];

                g_x_z_yzz_xyz[k] = -g_x_z_zz_xyz[k] * ab_y + g_x_z_zz_xyyz[k];

                g_x_z_yzz_xzz[k] = -g_x_z_zz_xzz[k] * ab_y + g_x_z_zz_xyzz[k];

                g_x_z_yzz_yyy[k] = -g_x_z_zz_yyy[k] * ab_y + g_x_z_zz_yyyy[k];

                g_x_z_yzz_yyz[k] = -g_x_z_zz_yyz[k] * ab_y + g_x_z_zz_yyyz[k];

                g_x_z_yzz_yzz[k] = -g_x_z_zz_yzz[k] * ab_y + g_x_z_zz_yyzz[k];

                g_x_z_yzz_zzz[k] = -g_x_z_zz_zzz[k] * ab_y + g_x_z_zz_yzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzz_xxx = cbuffer.data(ff_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_z_zzz_xxy = cbuffer.data(ff_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_z_zzz_xxz = cbuffer.data(ff_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_z_zzz_xyy = cbuffer.data(ff_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_z_zzz_xyz = cbuffer.data(ff_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_z_zzz_xzz = cbuffer.data(ff_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_z_zzz_yyy = cbuffer.data(ff_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_z_zzz_yyz = cbuffer.data(ff_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_z_zzz_yzz = cbuffer.data(ff_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_z_zzz_zzz = cbuffer.data(ff_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xxx, g_x_0_zz_xxy, g_x_0_zz_xxz, g_x_0_zz_xyy, g_x_0_zz_xyz, g_x_0_zz_xzz, g_x_0_zz_yyy, g_x_0_zz_yyz, g_x_0_zz_yzz, g_x_0_zz_zzz, g_x_z_zz_xxx, g_x_z_zz_xxxz, g_x_z_zz_xxy, g_x_z_zz_xxyz, g_x_z_zz_xxz, g_x_z_zz_xxzz, g_x_z_zz_xyy, g_x_z_zz_xyyz, g_x_z_zz_xyz, g_x_z_zz_xyzz, g_x_z_zz_xzz, g_x_z_zz_xzzz, g_x_z_zz_yyy, g_x_z_zz_yyyz, g_x_z_zz_yyz, g_x_z_zz_yyzz, g_x_z_zz_yzz, g_x_z_zz_yzzz, g_x_z_zz_zzz, g_x_z_zz_zzzz, g_x_z_zzz_xxx, g_x_z_zzz_xxy, g_x_z_zzz_xxz, g_x_z_zzz_xyy, g_x_z_zzz_xyz, g_x_z_zzz_xzz, g_x_z_zzz_yyy, g_x_z_zzz_yyz, g_x_z_zzz_yzz, g_x_z_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzz_xxx[k] = g_x_0_zz_xxx[k] - g_x_z_zz_xxx[k] * ab_z + g_x_z_zz_xxxz[k];

                g_x_z_zzz_xxy[k] = g_x_0_zz_xxy[k] - g_x_z_zz_xxy[k] * ab_z + g_x_z_zz_xxyz[k];

                g_x_z_zzz_xxz[k] = g_x_0_zz_xxz[k] - g_x_z_zz_xxz[k] * ab_z + g_x_z_zz_xxzz[k];

                g_x_z_zzz_xyy[k] = g_x_0_zz_xyy[k] - g_x_z_zz_xyy[k] * ab_z + g_x_z_zz_xyyz[k];

                g_x_z_zzz_xyz[k] = g_x_0_zz_xyz[k] - g_x_z_zz_xyz[k] * ab_z + g_x_z_zz_xyzz[k];

                g_x_z_zzz_xzz[k] = g_x_0_zz_xzz[k] - g_x_z_zz_xzz[k] * ab_z + g_x_z_zz_xzzz[k];

                g_x_z_zzz_yyy[k] = g_x_0_zz_yyy[k] - g_x_z_zz_yyy[k] * ab_z + g_x_z_zz_yyyz[k];

                g_x_z_zzz_yyz[k] = g_x_0_zz_yyz[k] - g_x_z_zz_yyz[k] * ab_z + g_x_z_zz_yyzz[k];

                g_x_z_zzz_yzz[k] = g_x_0_zz_yzz[k] - g_x_z_zz_yzz[k] * ab_z + g_x_z_zz_yzzz[k];

                g_x_z_zzz_zzz[k] = g_x_0_zz_zzz[k] - g_x_z_zz_zzz[k] * ab_z + g_x_z_zz_zzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxx_xxx = cbuffer.data(ff_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_xxx_xxy = cbuffer.data(ff_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_xxx_xxz = cbuffer.data(ff_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_xxx_xyy = cbuffer.data(ff_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_xxx_xyz = cbuffer.data(ff_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_xxx_xzz = cbuffer.data(ff_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_x_xxx_yyy = cbuffer.data(ff_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_xxx_yyz = cbuffer.data(ff_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_x_xxx_yzz = cbuffer.data(ff_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_xxx_zzz = cbuffer.data(ff_geom_11_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_xxx, g_y_0_xx_xxy, g_y_0_xx_xxz, g_y_0_xx_xyy, g_y_0_xx_xyz, g_y_0_xx_xzz, g_y_0_xx_yyy, g_y_0_xx_yyz, g_y_0_xx_yzz, g_y_0_xx_zzz, g_y_x_xx_xxx, g_y_x_xx_xxxx, g_y_x_xx_xxxy, g_y_x_xx_xxxz, g_y_x_xx_xxy, g_y_x_xx_xxyy, g_y_x_xx_xxyz, g_y_x_xx_xxz, g_y_x_xx_xxzz, g_y_x_xx_xyy, g_y_x_xx_xyyy, g_y_x_xx_xyyz, g_y_x_xx_xyz, g_y_x_xx_xyzz, g_y_x_xx_xzz, g_y_x_xx_xzzz, g_y_x_xx_yyy, g_y_x_xx_yyz, g_y_x_xx_yzz, g_y_x_xx_zzz, g_y_x_xxx_xxx, g_y_x_xxx_xxy, g_y_x_xxx_xxz, g_y_x_xxx_xyy, g_y_x_xxx_xyz, g_y_x_xxx_xzz, g_y_x_xxx_yyy, g_y_x_xxx_yyz, g_y_x_xxx_yzz, g_y_x_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxx_xxx[k] = g_y_0_xx_xxx[k] - g_y_x_xx_xxx[k] * ab_x + g_y_x_xx_xxxx[k];

                g_y_x_xxx_xxy[k] = g_y_0_xx_xxy[k] - g_y_x_xx_xxy[k] * ab_x + g_y_x_xx_xxxy[k];

                g_y_x_xxx_xxz[k] = g_y_0_xx_xxz[k] - g_y_x_xx_xxz[k] * ab_x + g_y_x_xx_xxxz[k];

                g_y_x_xxx_xyy[k] = g_y_0_xx_xyy[k] - g_y_x_xx_xyy[k] * ab_x + g_y_x_xx_xxyy[k];

                g_y_x_xxx_xyz[k] = g_y_0_xx_xyz[k] - g_y_x_xx_xyz[k] * ab_x + g_y_x_xx_xxyz[k];

                g_y_x_xxx_xzz[k] = g_y_0_xx_xzz[k] - g_y_x_xx_xzz[k] * ab_x + g_y_x_xx_xxzz[k];

                g_y_x_xxx_yyy[k] = g_y_0_xx_yyy[k] - g_y_x_xx_yyy[k] * ab_x + g_y_x_xx_xyyy[k];

                g_y_x_xxx_yyz[k] = g_y_0_xx_yyz[k] - g_y_x_xx_yyz[k] * ab_x + g_y_x_xx_xyyz[k];

                g_y_x_xxx_yzz[k] = g_y_0_xx_yzz[k] - g_y_x_xx_yzz[k] * ab_x + g_y_x_xx_xyzz[k];

                g_y_x_xxx_zzz[k] = g_y_0_xx_zzz[k] - g_y_x_xx_zzz[k] * ab_x + g_y_x_xx_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxy_xxx = cbuffer.data(ff_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_xxy_xxy = cbuffer.data(ff_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_x_xxy_xxz = cbuffer.data(ff_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_xxy_xyy = cbuffer.data(ff_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_xxy_xyz = cbuffer.data(ff_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_x_xxy_xzz = cbuffer.data(ff_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_xxy_yyy = cbuffer.data(ff_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_xxy_yyz = cbuffer.data(ff_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_x_xxy_yzz = cbuffer.data(ff_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_xxy_zzz = cbuffer.data(ff_geom_11_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xxx, g_y_0_xy_xxy, g_y_0_xy_xxz, g_y_0_xy_xyy, g_y_0_xy_xyz, g_y_0_xy_xzz, g_y_0_xy_yyy, g_y_0_xy_yyz, g_y_0_xy_yzz, g_y_0_xy_zzz, g_y_x_xxy_xxx, g_y_x_xxy_xxy, g_y_x_xxy_xxz, g_y_x_xxy_xyy, g_y_x_xxy_xyz, g_y_x_xxy_xzz, g_y_x_xxy_yyy, g_y_x_xxy_yyz, g_y_x_xxy_yzz, g_y_x_xxy_zzz, g_y_x_xy_xxx, g_y_x_xy_xxxx, g_y_x_xy_xxxy, g_y_x_xy_xxxz, g_y_x_xy_xxy, g_y_x_xy_xxyy, g_y_x_xy_xxyz, g_y_x_xy_xxz, g_y_x_xy_xxzz, g_y_x_xy_xyy, g_y_x_xy_xyyy, g_y_x_xy_xyyz, g_y_x_xy_xyz, g_y_x_xy_xyzz, g_y_x_xy_xzz, g_y_x_xy_xzzz, g_y_x_xy_yyy, g_y_x_xy_yyz, g_y_x_xy_yzz, g_y_x_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxy_xxx[k] = g_y_0_xy_xxx[k] - g_y_x_xy_xxx[k] * ab_x + g_y_x_xy_xxxx[k];

                g_y_x_xxy_xxy[k] = g_y_0_xy_xxy[k] - g_y_x_xy_xxy[k] * ab_x + g_y_x_xy_xxxy[k];

                g_y_x_xxy_xxz[k] = g_y_0_xy_xxz[k] - g_y_x_xy_xxz[k] * ab_x + g_y_x_xy_xxxz[k];

                g_y_x_xxy_xyy[k] = g_y_0_xy_xyy[k] - g_y_x_xy_xyy[k] * ab_x + g_y_x_xy_xxyy[k];

                g_y_x_xxy_xyz[k] = g_y_0_xy_xyz[k] - g_y_x_xy_xyz[k] * ab_x + g_y_x_xy_xxyz[k];

                g_y_x_xxy_xzz[k] = g_y_0_xy_xzz[k] - g_y_x_xy_xzz[k] * ab_x + g_y_x_xy_xxzz[k];

                g_y_x_xxy_yyy[k] = g_y_0_xy_yyy[k] - g_y_x_xy_yyy[k] * ab_x + g_y_x_xy_xyyy[k];

                g_y_x_xxy_yyz[k] = g_y_0_xy_yyz[k] - g_y_x_xy_yyz[k] * ab_x + g_y_x_xy_xyyz[k];

                g_y_x_xxy_yzz[k] = g_y_0_xy_yzz[k] - g_y_x_xy_yzz[k] * ab_x + g_y_x_xy_xyzz[k];

                g_y_x_xxy_zzz[k] = g_y_0_xy_zzz[k] - g_y_x_xy_zzz[k] * ab_x + g_y_x_xy_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxz_xxx = cbuffer.data(ff_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_xxz_xxy = cbuffer.data(ff_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_xxz_xxz = cbuffer.data(ff_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_xxz_xyy = cbuffer.data(ff_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_xxz_xyz = cbuffer.data(ff_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_xxz_xzz = cbuffer.data(ff_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_xxz_yyy = cbuffer.data(ff_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_xxz_yyz = cbuffer.data(ff_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_xxz_yzz = cbuffer.data(ff_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_xxz_zzz = cbuffer.data(ff_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xx_xxx, g_y_x_xx_xxxz, g_y_x_xx_xxy, g_y_x_xx_xxyz, g_y_x_xx_xxz, g_y_x_xx_xxzz, g_y_x_xx_xyy, g_y_x_xx_xyyz, g_y_x_xx_xyz, g_y_x_xx_xyzz, g_y_x_xx_xzz, g_y_x_xx_xzzz, g_y_x_xx_yyy, g_y_x_xx_yyyz, g_y_x_xx_yyz, g_y_x_xx_yyzz, g_y_x_xx_yzz, g_y_x_xx_yzzz, g_y_x_xx_zzz, g_y_x_xx_zzzz, g_y_x_xxz_xxx, g_y_x_xxz_xxy, g_y_x_xxz_xxz, g_y_x_xxz_xyy, g_y_x_xxz_xyz, g_y_x_xxz_xzz, g_y_x_xxz_yyy, g_y_x_xxz_yyz, g_y_x_xxz_yzz, g_y_x_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxz_xxx[k] = -g_y_x_xx_xxx[k] * ab_z + g_y_x_xx_xxxz[k];

                g_y_x_xxz_xxy[k] = -g_y_x_xx_xxy[k] * ab_z + g_y_x_xx_xxyz[k];

                g_y_x_xxz_xxz[k] = -g_y_x_xx_xxz[k] * ab_z + g_y_x_xx_xxzz[k];

                g_y_x_xxz_xyy[k] = -g_y_x_xx_xyy[k] * ab_z + g_y_x_xx_xyyz[k];

                g_y_x_xxz_xyz[k] = -g_y_x_xx_xyz[k] * ab_z + g_y_x_xx_xyzz[k];

                g_y_x_xxz_xzz[k] = -g_y_x_xx_xzz[k] * ab_z + g_y_x_xx_xzzz[k];

                g_y_x_xxz_yyy[k] = -g_y_x_xx_yyy[k] * ab_z + g_y_x_xx_yyyz[k];

                g_y_x_xxz_yyz[k] = -g_y_x_xx_yyz[k] * ab_z + g_y_x_xx_yyzz[k];

                g_y_x_xxz_yzz[k] = -g_y_x_xx_yzz[k] * ab_z + g_y_x_xx_yzzz[k];

                g_y_x_xxz_zzz[k] = -g_y_x_xx_zzz[k] * ab_z + g_y_x_xx_zzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyy_xxx = cbuffer.data(ff_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_xyy_xxy = cbuffer.data(ff_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_xyy_xxz = cbuffer.data(ff_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_xyy_xyy = cbuffer.data(ff_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_xyy_xyz = cbuffer.data(ff_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_xyy_xzz = cbuffer.data(ff_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_x_xyy_yyy = cbuffer.data(ff_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_xyy_yyz = cbuffer.data(ff_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_xyy_yzz = cbuffer.data(ff_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_xyy_zzz = cbuffer.data(ff_geom_11_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xzz, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yzz, g_y_0_yy_zzz, g_y_x_xyy_xxx, g_y_x_xyy_xxy, g_y_x_xyy_xxz, g_y_x_xyy_xyy, g_y_x_xyy_xyz, g_y_x_xyy_xzz, g_y_x_xyy_yyy, g_y_x_xyy_yyz, g_y_x_xyy_yzz, g_y_x_xyy_zzz, g_y_x_yy_xxx, g_y_x_yy_xxxx, g_y_x_yy_xxxy, g_y_x_yy_xxxz, g_y_x_yy_xxy, g_y_x_yy_xxyy, g_y_x_yy_xxyz, g_y_x_yy_xxz, g_y_x_yy_xxzz, g_y_x_yy_xyy, g_y_x_yy_xyyy, g_y_x_yy_xyyz, g_y_x_yy_xyz, g_y_x_yy_xyzz, g_y_x_yy_xzz, g_y_x_yy_xzzz, g_y_x_yy_yyy, g_y_x_yy_yyz, g_y_x_yy_yzz, g_y_x_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyy_xxx[k] = g_y_0_yy_xxx[k] - g_y_x_yy_xxx[k] * ab_x + g_y_x_yy_xxxx[k];

                g_y_x_xyy_xxy[k] = g_y_0_yy_xxy[k] - g_y_x_yy_xxy[k] * ab_x + g_y_x_yy_xxxy[k];

                g_y_x_xyy_xxz[k] = g_y_0_yy_xxz[k] - g_y_x_yy_xxz[k] * ab_x + g_y_x_yy_xxxz[k];

                g_y_x_xyy_xyy[k] = g_y_0_yy_xyy[k] - g_y_x_yy_xyy[k] * ab_x + g_y_x_yy_xxyy[k];

                g_y_x_xyy_xyz[k] = g_y_0_yy_xyz[k] - g_y_x_yy_xyz[k] * ab_x + g_y_x_yy_xxyz[k];

                g_y_x_xyy_xzz[k] = g_y_0_yy_xzz[k] - g_y_x_yy_xzz[k] * ab_x + g_y_x_yy_xxzz[k];

                g_y_x_xyy_yyy[k] = g_y_0_yy_yyy[k] - g_y_x_yy_yyy[k] * ab_x + g_y_x_yy_xyyy[k];

                g_y_x_xyy_yyz[k] = g_y_0_yy_yyz[k] - g_y_x_yy_yyz[k] * ab_x + g_y_x_yy_xyyz[k];

                g_y_x_xyy_yzz[k] = g_y_0_yy_yzz[k] - g_y_x_yy_yzz[k] * ab_x + g_y_x_yy_xyzz[k];

                g_y_x_xyy_zzz[k] = g_y_0_yy_zzz[k] - g_y_x_yy_zzz[k] * ab_x + g_y_x_yy_xzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyz_xxx = cbuffer.data(ff_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_xyz_xxy = cbuffer.data(ff_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_x_xyz_xxz = cbuffer.data(ff_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_xyz_xyy = cbuffer.data(ff_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_xyz_xyz = cbuffer.data(ff_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_x_xyz_xzz = cbuffer.data(ff_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_xyz_yyy = cbuffer.data(ff_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_xyz_yyz = cbuffer.data(ff_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_x_xyz_yzz = cbuffer.data(ff_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_xyz_zzz = cbuffer.data(ff_geom_11_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xy_xxx, g_y_x_xy_xxxz, g_y_x_xy_xxy, g_y_x_xy_xxyz, g_y_x_xy_xxz, g_y_x_xy_xxzz, g_y_x_xy_xyy, g_y_x_xy_xyyz, g_y_x_xy_xyz, g_y_x_xy_xyzz, g_y_x_xy_xzz, g_y_x_xy_xzzz, g_y_x_xy_yyy, g_y_x_xy_yyyz, g_y_x_xy_yyz, g_y_x_xy_yyzz, g_y_x_xy_yzz, g_y_x_xy_yzzz, g_y_x_xy_zzz, g_y_x_xy_zzzz, g_y_x_xyz_xxx, g_y_x_xyz_xxy, g_y_x_xyz_xxz, g_y_x_xyz_xyy, g_y_x_xyz_xyz, g_y_x_xyz_xzz, g_y_x_xyz_yyy, g_y_x_xyz_yyz, g_y_x_xyz_yzz, g_y_x_xyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyz_xxx[k] = -g_y_x_xy_xxx[k] * ab_z + g_y_x_xy_xxxz[k];

                g_y_x_xyz_xxy[k] = -g_y_x_xy_xxy[k] * ab_z + g_y_x_xy_xxyz[k];

                g_y_x_xyz_xxz[k] = -g_y_x_xy_xxz[k] * ab_z + g_y_x_xy_xxzz[k];

                g_y_x_xyz_xyy[k] = -g_y_x_xy_xyy[k] * ab_z + g_y_x_xy_xyyz[k];

                g_y_x_xyz_xyz[k] = -g_y_x_xy_xyz[k] * ab_z + g_y_x_xy_xyzz[k];

                g_y_x_xyz_xzz[k] = -g_y_x_xy_xzz[k] * ab_z + g_y_x_xy_xzzz[k];

                g_y_x_xyz_yyy[k] = -g_y_x_xy_yyy[k] * ab_z + g_y_x_xy_yyyz[k];

                g_y_x_xyz_yyz[k] = -g_y_x_xy_yyz[k] * ab_z + g_y_x_xy_yyzz[k];

                g_y_x_xyz_yzz[k] = -g_y_x_xy_yzz[k] * ab_z + g_y_x_xy_yzzz[k];

                g_y_x_xyz_zzz[k] = -g_y_x_xy_zzz[k] * ab_z + g_y_x_xy_zzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzz_xxx = cbuffer.data(ff_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_xzz_xxy = cbuffer.data(ff_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_xzz_xxz = cbuffer.data(ff_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_xzz_xyy = cbuffer.data(ff_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_x_xzz_xyz = cbuffer.data(ff_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_xzz_xzz = cbuffer.data(ff_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_xzz_yyy = cbuffer.data(ff_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_xzz_yyz = cbuffer.data(ff_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_xzz_yzz = cbuffer.data(ff_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_xzz_zzz = cbuffer.data(ff_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xz_xxx, g_y_x_xz_xxxz, g_y_x_xz_xxy, g_y_x_xz_xxyz, g_y_x_xz_xxz, g_y_x_xz_xxzz, g_y_x_xz_xyy, g_y_x_xz_xyyz, g_y_x_xz_xyz, g_y_x_xz_xyzz, g_y_x_xz_xzz, g_y_x_xz_xzzz, g_y_x_xz_yyy, g_y_x_xz_yyyz, g_y_x_xz_yyz, g_y_x_xz_yyzz, g_y_x_xz_yzz, g_y_x_xz_yzzz, g_y_x_xz_zzz, g_y_x_xz_zzzz, g_y_x_xzz_xxx, g_y_x_xzz_xxy, g_y_x_xzz_xxz, g_y_x_xzz_xyy, g_y_x_xzz_xyz, g_y_x_xzz_xzz, g_y_x_xzz_yyy, g_y_x_xzz_yyz, g_y_x_xzz_yzz, g_y_x_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzz_xxx[k] = -g_y_x_xz_xxx[k] * ab_z + g_y_x_xz_xxxz[k];

                g_y_x_xzz_xxy[k] = -g_y_x_xz_xxy[k] * ab_z + g_y_x_xz_xxyz[k];

                g_y_x_xzz_xxz[k] = -g_y_x_xz_xxz[k] * ab_z + g_y_x_xz_xxzz[k];

                g_y_x_xzz_xyy[k] = -g_y_x_xz_xyy[k] * ab_z + g_y_x_xz_xyyz[k];

                g_y_x_xzz_xyz[k] = -g_y_x_xz_xyz[k] * ab_z + g_y_x_xz_xyzz[k];

                g_y_x_xzz_xzz[k] = -g_y_x_xz_xzz[k] * ab_z + g_y_x_xz_xzzz[k];

                g_y_x_xzz_yyy[k] = -g_y_x_xz_yyy[k] * ab_z + g_y_x_xz_yyyz[k];

                g_y_x_xzz_yyz[k] = -g_y_x_xz_yyz[k] * ab_z + g_y_x_xz_yyzz[k];

                g_y_x_xzz_yzz[k] = -g_y_x_xz_yzz[k] * ab_z + g_y_x_xz_yzzz[k];

                g_y_x_xzz_zzz[k] = -g_y_x_xz_zzz[k] * ab_z + g_y_x_xz_zzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyy_xxx = cbuffer.data(ff_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_x_yyy_xxy = cbuffer.data(ff_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_x_yyy_xxz = cbuffer.data(ff_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_x_yyy_xyy = cbuffer.data(ff_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_x_yyy_xyz = cbuffer.data(ff_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_x_yyy_xzz = cbuffer.data(ff_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_x_yyy_yyy = cbuffer.data(ff_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_x_yyy_yyz = cbuffer.data(ff_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_x_yyy_yzz = cbuffer.data(ff_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_x_yyy_zzz = cbuffer.data(ff_geom_11_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_xxx, g_0_x_yy_xxy, g_0_x_yy_xxz, g_0_x_yy_xyy, g_0_x_yy_xyz, g_0_x_yy_xzz, g_0_x_yy_yyy, g_0_x_yy_yyz, g_0_x_yy_yzz, g_0_x_yy_zzz, g_y_x_yy_xxx, g_y_x_yy_xxxy, g_y_x_yy_xxy, g_y_x_yy_xxyy, g_y_x_yy_xxyz, g_y_x_yy_xxz, g_y_x_yy_xyy, g_y_x_yy_xyyy, g_y_x_yy_xyyz, g_y_x_yy_xyz, g_y_x_yy_xyzz, g_y_x_yy_xzz, g_y_x_yy_yyy, g_y_x_yy_yyyy, g_y_x_yy_yyyz, g_y_x_yy_yyz, g_y_x_yy_yyzz, g_y_x_yy_yzz, g_y_x_yy_yzzz, g_y_x_yy_zzz, g_y_x_yyy_xxx, g_y_x_yyy_xxy, g_y_x_yyy_xxz, g_y_x_yyy_xyy, g_y_x_yyy_xyz, g_y_x_yyy_xzz, g_y_x_yyy_yyy, g_y_x_yyy_yyz, g_y_x_yyy_yzz, g_y_x_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyy_xxx[k] = -g_0_x_yy_xxx[k] - g_y_x_yy_xxx[k] * ab_y + g_y_x_yy_xxxy[k];

                g_y_x_yyy_xxy[k] = -g_0_x_yy_xxy[k] - g_y_x_yy_xxy[k] * ab_y + g_y_x_yy_xxyy[k];

                g_y_x_yyy_xxz[k] = -g_0_x_yy_xxz[k] - g_y_x_yy_xxz[k] * ab_y + g_y_x_yy_xxyz[k];

                g_y_x_yyy_xyy[k] = -g_0_x_yy_xyy[k] - g_y_x_yy_xyy[k] * ab_y + g_y_x_yy_xyyy[k];

                g_y_x_yyy_xyz[k] = -g_0_x_yy_xyz[k] - g_y_x_yy_xyz[k] * ab_y + g_y_x_yy_xyyz[k];

                g_y_x_yyy_xzz[k] = -g_0_x_yy_xzz[k] - g_y_x_yy_xzz[k] * ab_y + g_y_x_yy_xyzz[k];

                g_y_x_yyy_yyy[k] = -g_0_x_yy_yyy[k] - g_y_x_yy_yyy[k] * ab_y + g_y_x_yy_yyyy[k];

                g_y_x_yyy_yyz[k] = -g_0_x_yy_yyz[k] - g_y_x_yy_yyz[k] * ab_y + g_y_x_yy_yyyz[k];

                g_y_x_yyy_yzz[k] = -g_0_x_yy_yzz[k] - g_y_x_yy_yzz[k] * ab_y + g_y_x_yy_yyzz[k];

                g_y_x_yyy_zzz[k] = -g_0_x_yy_zzz[k] - g_y_x_yy_zzz[k] * ab_y + g_y_x_yy_yzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyz_xxx = cbuffer.data(ff_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_x_yyz_xxy = cbuffer.data(ff_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_x_yyz_xxz = cbuffer.data(ff_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_x_yyz_xyy = cbuffer.data(ff_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_x_yyz_xyz = cbuffer.data(ff_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_x_yyz_xzz = cbuffer.data(ff_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_x_yyz_yyy = cbuffer.data(ff_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_x_yyz_yyz = cbuffer.data(ff_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_x_yyz_yzz = cbuffer.data(ff_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_x_yyz_zzz = cbuffer.data(ff_geom_11_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yy_xxx, g_y_x_yy_xxxz, g_y_x_yy_xxy, g_y_x_yy_xxyz, g_y_x_yy_xxz, g_y_x_yy_xxzz, g_y_x_yy_xyy, g_y_x_yy_xyyz, g_y_x_yy_xyz, g_y_x_yy_xyzz, g_y_x_yy_xzz, g_y_x_yy_xzzz, g_y_x_yy_yyy, g_y_x_yy_yyyz, g_y_x_yy_yyz, g_y_x_yy_yyzz, g_y_x_yy_yzz, g_y_x_yy_yzzz, g_y_x_yy_zzz, g_y_x_yy_zzzz, g_y_x_yyz_xxx, g_y_x_yyz_xxy, g_y_x_yyz_xxz, g_y_x_yyz_xyy, g_y_x_yyz_xyz, g_y_x_yyz_xzz, g_y_x_yyz_yyy, g_y_x_yyz_yyz, g_y_x_yyz_yzz, g_y_x_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyz_xxx[k] = -g_y_x_yy_xxx[k] * ab_z + g_y_x_yy_xxxz[k];

                g_y_x_yyz_xxy[k] = -g_y_x_yy_xxy[k] * ab_z + g_y_x_yy_xxyz[k];

                g_y_x_yyz_xxz[k] = -g_y_x_yy_xxz[k] * ab_z + g_y_x_yy_xxzz[k];

                g_y_x_yyz_xyy[k] = -g_y_x_yy_xyy[k] * ab_z + g_y_x_yy_xyyz[k];

                g_y_x_yyz_xyz[k] = -g_y_x_yy_xyz[k] * ab_z + g_y_x_yy_xyzz[k];

                g_y_x_yyz_xzz[k] = -g_y_x_yy_xzz[k] * ab_z + g_y_x_yy_xzzz[k];

                g_y_x_yyz_yyy[k] = -g_y_x_yy_yyy[k] * ab_z + g_y_x_yy_yyyz[k];

                g_y_x_yyz_yyz[k] = -g_y_x_yy_yyz[k] * ab_z + g_y_x_yy_yyzz[k];

                g_y_x_yyz_yzz[k] = -g_y_x_yy_yzz[k] * ab_z + g_y_x_yy_yzzz[k];

                g_y_x_yyz_zzz[k] = -g_y_x_yy_zzz[k] * ab_z + g_y_x_yy_zzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzz_xxx = cbuffer.data(ff_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_x_yzz_xxy = cbuffer.data(ff_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_x_yzz_xxz = cbuffer.data(ff_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_x_yzz_xyy = cbuffer.data(ff_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_x_yzz_xyz = cbuffer.data(ff_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_x_yzz_xzz = cbuffer.data(ff_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_x_yzz_yyy = cbuffer.data(ff_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_x_yzz_yyz = cbuffer.data(ff_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_x_yzz_yzz = cbuffer.data(ff_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_x_yzz_zzz = cbuffer.data(ff_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yz_xxx, g_y_x_yz_xxxz, g_y_x_yz_xxy, g_y_x_yz_xxyz, g_y_x_yz_xxz, g_y_x_yz_xxzz, g_y_x_yz_xyy, g_y_x_yz_xyyz, g_y_x_yz_xyz, g_y_x_yz_xyzz, g_y_x_yz_xzz, g_y_x_yz_xzzz, g_y_x_yz_yyy, g_y_x_yz_yyyz, g_y_x_yz_yyz, g_y_x_yz_yyzz, g_y_x_yz_yzz, g_y_x_yz_yzzz, g_y_x_yz_zzz, g_y_x_yz_zzzz, g_y_x_yzz_xxx, g_y_x_yzz_xxy, g_y_x_yzz_xxz, g_y_x_yzz_xyy, g_y_x_yzz_xyz, g_y_x_yzz_xzz, g_y_x_yzz_yyy, g_y_x_yzz_yyz, g_y_x_yzz_yzz, g_y_x_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzz_xxx[k] = -g_y_x_yz_xxx[k] * ab_z + g_y_x_yz_xxxz[k];

                g_y_x_yzz_xxy[k] = -g_y_x_yz_xxy[k] * ab_z + g_y_x_yz_xxyz[k];

                g_y_x_yzz_xxz[k] = -g_y_x_yz_xxz[k] * ab_z + g_y_x_yz_xxzz[k];

                g_y_x_yzz_xyy[k] = -g_y_x_yz_xyy[k] * ab_z + g_y_x_yz_xyyz[k];

                g_y_x_yzz_xyz[k] = -g_y_x_yz_xyz[k] * ab_z + g_y_x_yz_xyzz[k];

                g_y_x_yzz_xzz[k] = -g_y_x_yz_xzz[k] * ab_z + g_y_x_yz_xzzz[k];

                g_y_x_yzz_yyy[k] = -g_y_x_yz_yyy[k] * ab_z + g_y_x_yz_yyyz[k];

                g_y_x_yzz_yyz[k] = -g_y_x_yz_yyz[k] * ab_z + g_y_x_yz_yyzz[k];

                g_y_x_yzz_yzz[k] = -g_y_x_yz_yzz[k] * ab_z + g_y_x_yz_yzzz[k];

                g_y_x_yzz_zzz[k] = -g_y_x_yz_zzz[k] * ab_z + g_y_x_yz_zzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzz_xxx = cbuffer.data(ff_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_x_zzz_xxy = cbuffer.data(ff_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_x_zzz_xxz = cbuffer.data(ff_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_x_zzz_xyy = cbuffer.data(ff_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_x_zzz_xyz = cbuffer.data(ff_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_x_zzz_xzz = cbuffer.data(ff_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_x_zzz_yyy = cbuffer.data(ff_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_x_zzz_yyz = cbuffer.data(ff_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_x_zzz_yzz = cbuffer.data(ff_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_x_zzz_zzz = cbuffer.data(ff_geom_11_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zz_xxx, g_y_x_zz_xxxz, g_y_x_zz_xxy, g_y_x_zz_xxyz, g_y_x_zz_xxz, g_y_x_zz_xxzz, g_y_x_zz_xyy, g_y_x_zz_xyyz, g_y_x_zz_xyz, g_y_x_zz_xyzz, g_y_x_zz_xzz, g_y_x_zz_xzzz, g_y_x_zz_yyy, g_y_x_zz_yyyz, g_y_x_zz_yyz, g_y_x_zz_yyzz, g_y_x_zz_yzz, g_y_x_zz_yzzz, g_y_x_zz_zzz, g_y_x_zz_zzzz, g_y_x_zzz_xxx, g_y_x_zzz_xxy, g_y_x_zzz_xxz, g_y_x_zzz_xyy, g_y_x_zzz_xyz, g_y_x_zzz_xzz, g_y_x_zzz_yyy, g_y_x_zzz_yyz, g_y_x_zzz_yzz, g_y_x_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzz_xxx[k] = -g_y_x_zz_xxx[k] * ab_z + g_y_x_zz_xxxz[k];

                g_y_x_zzz_xxy[k] = -g_y_x_zz_xxy[k] * ab_z + g_y_x_zz_xxyz[k];

                g_y_x_zzz_xxz[k] = -g_y_x_zz_xxz[k] * ab_z + g_y_x_zz_xxzz[k];

                g_y_x_zzz_xyy[k] = -g_y_x_zz_xyy[k] * ab_z + g_y_x_zz_xyyz[k];

                g_y_x_zzz_xyz[k] = -g_y_x_zz_xyz[k] * ab_z + g_y_x_zz_xyzz[k];

                g_y_x_zzz_xzz[k] = -g_y_x_zz_xzz[k] * ab_z + g_y_x_zz_xzzz[k];

                g_y_x_zzz_yyy[k] = -g_y_x_zz_yyy[k] * ab_z + g_y_x_zz_yyyz[k];

                g_y_x_zzz_yyz[k] = -g_y_x_zz_yyz[k] * ab_z + g_y_x_zz_yyzz[k];

                g_y_x_zzz_yzz[k] = -g_y_x_zz_yzz[k] * ab_z + g_y_x_zz_yzzz[k];

                g_y_x_zzz_zzz[k] = -g_y_x_zz_zzz[k] * ab_z + g_y_x_zz_zzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxx_xxx = cbuffer.data(ff_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_xxx_xxy = cbuffer.data(ff_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_y_xxx_xxz = cbuffer.data(ff_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_xxx_xyy = cbuffer.data(ff_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_xxx_xyz = cbuffer.data(ff_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_y_xxx_xzz = cbuffer.data(ff_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_xxx_yyy = cbuffer.data(ff_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_xxx_yyz = cbuffer.data(ff_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_y_xxx_yzz = cbuffer.data(ff_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_xxx_zzz = cbuffer.data(ff_geom_11_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xx_xxx, g_y_y_xx_xxxx, g_y_y_xx_xxxy, g_y_y_xx_xxxz, g_y_y_xx_xxy, g_y_y_xx_xxyy, g_y_y_xx_xxyz, g_y_y_xx_xxz, g_y_y_xx_xxzz, g_y_y_xx_xyy, g_y_y_xx_xyyy, g_y_y_xx_xyyz, g_y_y_xx_xyz, g_y_y_xx_xyzz, g_y_y_xx_xzz, g_y_y_xx_xzzz, g_y_y_xx_yyy, g_y_y_xx_yyz, g_y_y_xx_yzz, g_y_y_xx_zzz, g_y_y_xxx_xxx, g_y_y_xxx_xxy, g_y_y_xxx_xxz, g_y_y_xxx_xyy, g_y_y_xxx_xyz, g_y_y_xxx_xzz, g_y_y_xxx_yyy, g_y_y_xxx_yyz, g_y_y_xxx_yzz, g_y_y_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxx_xxx[k] = -g_y_y_xx_xxx[k] * ab_x + g_y_y_xx_xxxx[k];

                g_y_y_xxx_xxy[k] = -g_y_y_xx_xxy[k] * ab_x + g_y_y_xx_xxxy[k];

                g_y_y_xxx_xxz[k] = -g_y_y_xx_xxz[k] * ab_x + g_y_y_xx_xxxz[k];

                g_y_y_xxx_xyy[k] = -g_y_y_xx_xyy[k] * ab_x + g_y_y_xx_xxyy[k];

                g_y_y_xxx_xyz[k] = -g_y_y_xx_xyz[k] * ab_x + g_y_y_xx_xxyz[k];

                g_y_y_xxx_xzz[k] = -g_y_y_xx_xzz[k] * ab_x + g_y_y_xx_xxzz[k];

                g_y_y_xxx_yyy[k] = -g_y_y_xx_yyy[k] * ab_x + g_y_y_xx_xyyy[k];

                g_y_y_xxx_yyz[k] = -g_y_y_xx_yyz[k] * ab_x + g_y_y_xx_xyyz[k];

                g_y_y_xxx_yzz[k] = -g_y_y_xx_yzz[k] * ab_x + g_y_y_xx_xyzz[k];

                g_y_y_xxx_zzz[k] = -g_y_y_xx_zzz[k] * ab_x + g_y_y_xx_xzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxy_xxx = cbuffer.data(ff_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_xxy_xxy = cbuffer.data(ff_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_xxy_xxz = cbuffer.data(ff_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_xxy_xyy = cbuffer.data(ff_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_y_xxy_xyz = cbuffer.data(ff_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_xxy_xzz = cbuffer.data(ff_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_xxy_yyy = cbuffer.data(ff_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_xxy_yyz = cbuffer.data(ff_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_xxy_yzz = cbuffer.data(ff_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_xxy_zzz = cbuffer.data(ff_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxy_xxx, g_y_y_xxy_xxy, g_y_y_xxy_xxz, g_y_y_xxy_xyy, g_y_y_xxy_xyz, g_y_y_xxy_xzz, g_y_y_xxy_yyy, g_y_y_xxy_yyz, g_y_y_xxy_yzz, g_y_y_xxy_zzz, g_y_y_xy_xxx, g_y_y_xy_xxxx, g_y_y_xy_xxxy, g_y_y_xy_xxxz, g_y_y_xy_xxy, g_y_y_xy_xxyy, g_y_y_xy_xxyz, g_y_y_xy_xxz, g_y_y_xy_xxzz, g_y_y_xy_xyy, g_y_y_xy_xyyy, g_y_y_xy_xyyz, g_y_y_xy_xyz, g_y_y_xy_xyzz, g_y_y_xy_xzz, g_y_y_xy_xzzz, g_y_y_xy_yyy, g_y_y_xy_yyz, g_y_y_xy_yzz, g_y_y_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxy_xxx[k] = -g_y_y_xy_xxx[k] * ab_x + g_y_y_xy_xxxx[k];

                g_y_y_xxy_xxy[k] = -g_y_y_xy_xxy[k] * ab_x + g_y_y_xy_xxxy[k];

                g_y_y_xxy_xxz[k] = -g_y_y_xy_xxz[k] * ab_x + g_y_y_xy_xxxz[k];

                g_y_y_xxy_xyy[k] = -g_y_y_xy_xyy[k] * ab_x + g_y_y_xy_xxyy[k];

                g_y_y_xxy_xyz[k] = -g_y_y_xy_xyz[k] * ab_x + g_y_y_xy_xxyz[k];

                g_y_y_xxy_xzz[k] = -g_y_y_xy_xzz[k] * ab_x + g_y_y_xy_xxzz[k];

                g_y_y_xxy_yyy[k] = -g_y_y_xy_yyy[k] * ab_x + g_y_y_xy_xyyy[k];

                g_y_y_xxy_yyz[k] = -g_y_y_xy_yyz[k] * ab_x + g_y_y_xy_xyyz[k];

                g_y_y_xxy_yzz[k] = -g_y_y_xy_yzz[k] * ab_x + g_y_y_xy_xyzz[k];

                g_y_y_xxy_zzz[k] = -g_y_y_xy_zzz[k] * ab_x + g_y_y_xy_xzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxz_xxx = cbuffer.data(ff_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_y_xxz_xxy = cbuffer.data(ff_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_y_xxz_xxz = cbuffer.data(ff_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_y_xxz_xyy = cbuffer.data(ff_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_y_xxz_xyz = cbuffer.data(ff_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_y_xxz_xzz = cbuffer.data(ff_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_y_xxz_yyy = cbuffer.data(ff_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_y_xxz_yyz = cbuffer.data(ff_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_y_xxz_yzz = cbuffer.data(ff_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_y_xxz_zzz = cbuffer.data(ff_geom_11_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxz_xxx, g_y_y_xxz_xxy, g_y_y_xxz_xxz, g_y_y_xxz_xyy, g_y_y_xxz_xyz, g_y_y_xxz_xzz, g_y_y_xxz_yyy, g_y_y_xxz_yyz, g_y_y_xxz_yzz, g_y_y_xxz_zzz, g_y_y_xz_xxx, g_y_y_xz_xxxx, g_y_y_xz_xxxy, g_y_y_xz_xxxz, g_y_y_xz_xxy, g_y_y_xz_xxyy, g_y_y_xz_xxyz, g_y_y_xz_xxz, g_y_y_xz_xxzz, g_y_y_xz_xyy, g_y_y_xz_xyyy, g_y_y_xz_xyyz, g_y_y_xz_xyz, g_y_y_xz_xyzz, g_y_y_xz_xzz, g_y_y_xz_xzzz, g_y_y_xz_yyy, g_y_y_xz_yyz, g_y_y_xz_yzz, g_y_y_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxz_xxx[k] = -g_y_y_xz_xxx[k] * ab_x + g_y_y_xz_xxxx[k];

                g_y_y_xxz_xxy[k] = -g_y_y_xz_xxy[k] * ab_x + g_y_y_xz_xxxy[k];

                g_y_y_xxz_xxz[k] = -g_y_y_xz_xxz[k] * ab_x + g_y_y_xz_xxxz[k];

                g_y_y_xxz_xyy[k] = -g_y_y_xz_xyy[k] * ab_x + g_y_y_xz_xxyy[k];

                g_y_y_xxz_xyz[k] = -g_y_y_xz_xyz[k] * ab_x + g_y_y_xz_xxyz[k];

                g_y_y_xxz_xzz[k] = -g_y_y_xz_xzz[k] * ab_x + g_y_y_xz_xxzz[k];

                g_y_y_xxz_yyy[k] = -g_y_y_xz_yyy[k] * ab_x + g_y_y_xz_xyyy[k];

                g_y_y_xxz_yyz[k] = -g_y_y_xz_yyz[k] * ab_x + g_y_y_xz_xyyz[k];

                g_y_y_xxz_yzz[k] = -g_y_y_xz_yzz[k] * ab_x + g_y_y_xz_xyzz[k];

                g_y_y_xxz_zzz[k] = -g_y_y_xz_zzz[k] * ab_x + g_y_y_xz_xzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyy_xxx = cbuffer.data(ff_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_y_xyy_xxy = cbuffer.data(ff_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_y_xyy_xxz = cbuffer.data(ff_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_xyy_xyy = cbuffer.data(ff_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_xyy_xyz = cbuffer.data(ff_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_y_xyy_xzz = cbuffer.data(ff_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_xyy_yyy = cbuffer.data(ff_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_xyy_yyz = cbuffer.data(ff_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_y_xyy_yzz = cbuffer.data(ff_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_xyy_zzz = cbuffer.data(ff_geom_11_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyy_xxx, g_y_y_xyy_xxy, g_y_y_xyy_xxz, g_y_y_xyy_xyy, g_y_y_xyy_xyz, g_y_y_xyy_xzz, g_y_y_xyy_yyy, g_y_y_xyy_yyz, g_y_y_xyy_yzz, g_y_y_xyy_zzz, g_y_y_yy_xxx, g_y_y_yy_xxxx, g_y_y_yy_xxxy, g_y_y_yy_xxxz, g_y_y_yy_xxy, g_y_y_yy_xxyy, g_y_y_yy_xxyz, g_y_y_yy_xxz, g_y_y_yy_xxzz, g_y_y_yy_xyy, g_y_y_yy_xyyy, g_y_y_yy_xyyz, g_y_y_yy_xyz, g_y_y_yy_xyzz, g_y_y_yy_xzz, g_y_y_yy_xzzz, g_y_y_yy_yyy, g_y_y_yy_yyz, g_y_y_yy_yzz, g_y_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyy_xxx[k] = -g_y_y_yy_xxx[k] * ab_x + g_y_y_yy_xxxx[k];

                g_y_y_xyy_xxy[k] = -g_y_y_yy_xxy[k] * ab_x + g_y_y_yy_xxxy[k];

                g_y_y_xyy_xxz[k] = -g_y_y_yy_xxz[k] * ab_x + g_y_y_yy_xxxz[k];

                g_y_y_xyy_xyy[k] = -g_y_y_yy_xyy[k] * ab_x + g_y_y_yy_xxyy[k];

                g_y_y_xyy_xyz[k] = -g_y_y_yy_xyz[k] * ab_x + g_y_y_yy_xxyz[k];

                g_y_y_xyy_xzz[k] = -g_y_y_yy_xzz[k] * ab_x + g_y_y_yy_xxzz[k];

                g_y_y_xyy_yyy[k] = -g_y_y_yy_yyy[k] * ab_x + g_y_y_yy_xyyy[k];

                g_y_y_xyy_yyz[k] = -g_y_y_yy_yyz[k] * ab_x + g_y_y_yy_xyyz[k];

                g_y_y_xyy_yzz[k] = -g_y_y_yy_yzz[k] * ab_x + g_y_y_yy_xyzz[k];

                g_y_y_xyy_zzz[k] = -g_y_y_yy_zzz[k] * ab_x + g_y_y_yy_xzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyz_xxx = cbuffer.data(ff_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_xyz_xxy = cbuffer.data(ff_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_xyz_xxz = cbuffer.data(ff_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_xyz_xyy = cbuffer.data(ff_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_y_xyz_xyz = cbuffer.data(ff_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_xyz_xzz = cbuffer.data(ff_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_xyz_yyy = cbuffer.data(ff_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_xyz_yyz = cbuffer.data(ff_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_xyz_yzz = cbuffer.data(ff_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_xyz_zzz = cbuffer.data(ff_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyz_xxx, g_y_y_xyz_xxy, g_y_y_xyz_xxz, g_y_y_xyz_xyy, g_y_y_xyz_xyz, g_y_y_xyz_xzz, g_y_y_xyz_yyy, g_y_y_xyz_yyz, g_y_y_xyz_yzz, g_y_y_xyz_zzz, g_y_y_yz_xxx, g_y_y_yz_xxxx, g_y_y_yz_xxxy, g_y_y_yz_xxxz, g_y_y_yz_xxy, g_y_y_yz_xxyy, g_y_y_yz_xxyz, g_y_y_yz_xxz, g_y_y_yz_xxzz, g_y_y_yz_xyy, g_y_y_yz_xyyy, g_y_y_yz_xyyz, g_y_y_yz_xyz, g_y_y_yz_xyzz, g_y_y_yz_xzz, g_y_y_yz_xzzz, g_y_y_yz_yyy, g_y_y_yz_yyz, g_y_y_yz_yzz, g_y_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyz_xxx[k] = -g_y_y_yz_xxx[k] * ab_x + g_y_y_yz_xxxx[k];

                g_y_y_xyz_xxy[k] = -g_y_y_yz_xxy[k] * ab_x + g_y_y_yz_xxxy[k];

                g_y_y_xyz_xxz[k] = -g_y_y_yz_xxz[k] * ab_x + g_y_y_yz_xxxz[k];

                g_y_y_xyz_xyy[k] = -g_y_y_yz_xyy[k] * ab_x + g_y_y_yz_xxyy[k];

                g_y_y_xyz_xyz[k] = -g_y_y_yz_xyz[k] * ab_x + g_y_y_yz_xxyz[k];

                g_y_y_xyz_xzz[k] = -g_y_y_yz_xzz[k] * ab_x + g_y_y_yz_xxzz[k];

                g_y_y_xyz_yyy[k] = -g_y_y_yz_yyy[k] * ab_x + g_y_y_yz_xyyy[k];

                g_y_y_xyz_yyz[k] = -g_y_y_yz_yyz[k] * ab_x + g_y_y_yz_xyyz[k];

                g_y_y_xyz_yzz[k] = -g_y_y_yz_yzz[k] * ab_x + g_y_y_yz_xyzz[k];

                g_y_y_xyz_zzz[k] = -g_y_y_yz_zzz[k] * ab_x + g_y_y_yz_xzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzz_xxx = cbuffer.data(ff_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_y_xzz_xxy = cbuffer.data(ff_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_y_xzz_xxz = cbuffer.data(ff_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_y_xzz_xyy = cbuffer.data(ff_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_y_xzz_xyz = cbuffer.data(ff_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_y_xzz_xzz = cbuffer.data(ff_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_y_xzz_yyy = cbuffer.data(ff_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_y_xzz_yyz = cbuffer.data(ff_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_y_xzz_yzz = cbuffer.data(ff_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_y_xzz_zzz = cbuffer.data(ff_geom_11_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzz_xxx, g_y_y_xzz_xxy, g_y_y_xzz_xxz, g_y_y_xzz_xyy, g_y_y_xzz_xyz, g_y_y_xzz_xzz, g_y_y_xzz_yyy, g_y_y_xzz_yyz, g_y_y_xzz_yzz, g_y_y_xzz_zzz, g_y_y_zz_xxx, g_y_y_zz_xxxx, g_y_y_zz_xxxy, g_y_y_zz_xxxz, g_y_y_zz_xxy, g_y_y_zz_xxyy, g_y_y_zz_xxyz, g_y_y_zz_xxz, g_y_y_zz_xxzz, g_y_y_zz_xyy, g_y_y_zz_xyyy, g_y_y_zz_xyyz, g_y_y_zz_xyz, g_y_y_zz_xyzz, g_y_y_zz_xzz, g_y_y_zz_xzzz, g_y_y_zz_yyy, g_y_y_zz_yyz, g_y_y_zz_yzz, g_y_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzz_xxx[k] = -g_y_y_zz_xxx[k] * ab_x + g_y_y_zz_xxxx[k];

                g_y_y_xzz_xxy[k] = -g_y_y_zz_xxy[k] * ab_x + g_y_y_zz_xxxy[k];

                g_y_y_xzz_xxz[k] = -g_y_y_zz_xxz[k] * ab_x + g_y_y_zz_xxxz[k];

                g_y_y_xzz_xyy[k] = -g_y_y_zz_xyy[k] * ab_x + g_y_y_zz_xxyy[k];

                g_y_y_xzz_xyz[k] = -g_y_y_zz_xyz[k] * ab_x + g_y_y_zz_xxyz[k];

                g_y_y_xzz_xzz[k] = -g_y_y_zz_xzz[k] * ab_x + g_y_y_zz_xxzz[k];

                g_y_y_xzz_yyy[k] = -g_y_y_zz_yyy[k] * ab_x + g_y_y_zz_xyyy[k];

                g_y_y_xzz_yyz[k] = -g_y_y_zz_yyz[k] * ab_x + g_y_y_zz_xyyz[k];

                g_y_y_xzz_yzz[k] = -g_y_y_zz_yzz[k] * ab_x + g_y_y_zz_xyzz[k];

                g_y_y_xzz_zzz[k] = -g_y_y_zz_zzz[k] * ab_x + g_y_y_zz_xzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyy_xxx = cbuffer.data(ff_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_y_yyy_xxy = cbuffer.data(ff_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_y_yyy_xxz = cbuffer.data(ff_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_y_yyy_xyy = cbuffer.data(ff_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_y_yyy_xyz = cbuffer.data(ff_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_y_yyy_xzz = cbuffer.data(ff_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_y_yyy_yyy = cbuffer.data(ff_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_y_yyy_yyz = cbuffer.data(ff_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_y_yyy_yzz = cbuffer.data(ff_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_y_yyy_zzz = cbuffer.data(ff_geom_11_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_xxx, g_0_y_yy_xxy, g_0_y_yy_xxz, g_0_y_yy_xyy, g_0_y_yy_xyz, g_0_y_yy_xzz, g_0_y_yy_yyy, g_0_y_yy_yyz, g_0_y_yy_yzz, g_0_y_yy_zzz, g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xzz, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yzz, g_y_0_yy_zzz, g_y_y_yy_xxx, g_y_y_yy_xxxy, g_y_y_yy_xxy, g_y_y_yy_xxyy, g_y_y_yy_xxyz, g_y_y_yy_xxz, g_y_y_yy_xyy, g_y_y_yy_xyyy, g_y_y_yy_xyyz, g_y_y_yy_xyz, g_y_y_yy_xyzz, g_y_y_yy_xzz, g_y_y_yy_yyy, g_y_y_yy_yyyy, g_y_y_yy_yyyz, g_y_y_yy_yyz, g_y_y_yy_yyzz, g_y_y_yy_yzz, g_y_y_yy_yzzz, g_y_y_yy_zzz, g_y_y_yyy_xxx, g_y_y_yyy_xxy, g_y_y_yyy_xxz, g_y_y_yyy_xyy, g_y_y_yyy_xyz, g_y_y_yyy_xzz, g_y_y_yyy_yyy, g_y_y_yyy_yyz, g_y_y_yyy_yzz, g_y_y_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyy_xxx[k] = -g_0_y_yy_xxx[k] + g_y_0_yy_xxx[k] - g_y_y_yy_xxx[k] * ab_y + g_y_y_yy_xxxy[k];

                g_y_y_yyy_xxy[k] = -g_0_y_yy_xxy[k] + g_y_0_yy_xxy[k] - g_y_y_yy_xxy[k] * ab_y + g_y_y_yy_xxyy[k];

                g_y_y_yyy_xxz[k] = -g_0_y_yy_xxz[k] + g_y_0_yy_xxz[k] - g_y_y_yy_xxz[k] * ab_y + g_y_y_yy_xxyz[k];

                g_y_y_yyy_xyy[k] = -g_0_y_yy_xyy[k] + g_y_0_yy_xyy[k] - g_y_y_yy_xyy[k] * ab_y + g_y_y_yy_xyyy[k];

                g_y_y_yyy_xyz[k] = -g_0_y_yy_xyz[k] + g_y_0_yy_xyz[k] - g_y_y_yy_xyz[k] * ab_y + g_y_y_yy_xyyz[k];

                g_y_y_yyy_xzz[k] = -g_0_y_yy_xzz[k] + g_y_0_yy_xzz[k] - g_y_y_yy_xzz[k] * ab_y + g_y_y_yy_xyzz[k];

                g_y_y_yyy_yyy[k] = -g_0_y_yy_yyy[k] + g_y_0_yy_yyy[k] - g_y_y_yy_yyy[k] * ab_y + g_y_y_yy_yyyy[k];

                g_y_y_yyy_yyz[k] = -g_0_y_yy_yyz[k] + g_y_0_yy_yyz[k] - g_y_y_yy_yyz[k] * ab_y + g_y_y_yy_yyyz[k];

                g_y_y_yyy_yzz[k] = -g_0_y_yy_yzz[k] + g_y_0_yy_yzz[k] - g_y_y_yy_yzz[k] * ab_y + g_y_y_yy_yyzz[k];

                g_y_y_yyy_zzz[k] = -g_0_y_yy_zzz[k] + g_y_0_yy_zzz[k] - g_y_y_yy_zzz[k] * ab_y + g_y_y_yy_yzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyz_xxx = cbuffer.data(ff_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_y_yyz_xxy = cbuffer.data(ff_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_y_yyz_xxz = cbuffer.data(ff_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_y_yyz_xyy = cbuffer.data(ff_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_y_yyz_xyz = cbuffer.data(ff_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_y_yyz_xzz = cbuffer.data(ff_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_y_yyz_yyy = cbuffer.data(ff_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_y_yyz_yyz = cbuffer.data(ff_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_y_yyz_yzz = cbuffer.data(ff_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_y_yyz_zzz = cbuffer.data(ff_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yy_xxx, g_y_y_yy_xxxz, g_y_y_yy_xxy, g_y_y_yy_xxyz, g_y_y_yy_xxz, g_y_y_yy_xxzz, g_y_y_yy_xyy, g_y_y_yy_xyyz, g_y_y_yy_xyz, g_y_y_yy_xyzz, g_y_y_yy_xzz, g_y_y_yy_xzzz, g_y_y_yy_yyy, g_y_y_yy_yyyz, g_y_y_yy_yyz, g_y_y_yy_yyzz, g_y_y_yy_yzz, g_y_y_yy_yzzz, g_y_y_yy_zzz, g_y_y_yy_zzzz, g_y_y_yyz_xxx, g_y_y_yyz_xxy, g_y_y_yyz_xxz, g_y_y_yyz_xyy, g_y_y_yyz_xyz, g_y_y_yyz_xzz, g_y_y_yyz_yyy, g_y_y_yyz_yyz, g_y_y_yyz_yzz, g_y_y_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyz_xxx[k] = -g_y_y_yy_xxx[k] * ab_z + g_y_y_yy_xxxz[k];

                g_y_y_yyz_xxy[k] = -g_y_y_yy_xxy[k] * ab_z + g_y_y_yy_xxyz[k];

                g_y_y_yyz_xxz[k] = -g_y_y_yy_xxz[k] * ab_z + g_y_y_yy_xxzz[k];

                g_y_y_yyz_xyy[k] = -g_y_y_yy_xyy[k] * ab_z + g_y_y_yy_xyyz[k];

                g_y_y_yyz_xyz[k] = -g_y_y_yy_xyz[k] * ab_z + g_y_y_yy_xyzz[k];

                g_y_y_yyz_xzz[k] = -g_y_y_yy_xzz[k] * ab_z + g_y_y_yy_xzzz[k];

                g_y_y_yyz_yyy[k] = -g_y_y_yy_yyy[k] * ab_z + g_y_y_yy_yyyz[k];

                g_y_y_yyz_yyz[k] = -g_y_y_yy_yyz[k] * ab_z + g_y_y_yy_yyzz[k];

                g_y_y_yyz_yzz[k] = -g_y_y_yy_yzz[k] * ab_z + g_y_y_yy_yzzz[k];

                g_y_y_yyz_zzz[k] = -g_y_y_yy_zzz[k] * ab_z + g_y_y_yy_zzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzz_xxx = cbuffer.data(ff_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_y_yzz_xxy = cbuffer.data(ff_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_y_yzz_xxz = cbuffer.data(ff_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_y_yzz_xyy = cbuffer.data(ff_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_y_yzz_xyz = cbuffer.data(ff_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_y_yzz_xzz = cbuffer.data(ff_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_y_yzz_yyy = cbuffer.data(ff_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_y_yzz_yyz = cbuffer.data(ff_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_y_yzz_yzz = cbuffer.data(ff_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_y_yzz_zzz = cbuffer.data(ff_geom_11_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yz_xxx, g_y_y_yz_xxxz, g_y_y_yz_xxy, g_y_y_yz_xxyz, g_y_y_yz_xxz, g_y_y_yz_xxzz, g_y_y_yz_xyy, g_y_y_yz_xyyz, g_y_y_yz_xyz, g_y_y_yz_xyzz, g_y_y_yz_xzz, g_y_y_yz_xzzz, g_y_y_yz_yyy, g_y_y_yz_yyyz, g_y_y_yz_yyz, g_y_y_yz_yyzz, g_y_y_yz_yzz, g_y_y_yz_yzzz, g_y_y_yz_zzz, g_y_y_yz_zzzz, g_y_y_yzz_xxx, g_y_y_yzz_xxy, g_y_y_yzz_xxz, g_y_y_yzz_xyy, g_y_y_yzz_xyz, g_y_y_yzz_xzz, g_y_y_yzz_yyy, g_y_y_yzz_yyz, g_y_y_yzz_yzz, g_y_y_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzz_xxx[k] = -g_y_y_yz_xxx[k] * ab_z + g_y_y_yz_xxxz[k];

                g_y_y_yzz_xxy[k] = -g_y_y_yz_xxy[k] * ab_z + g_y_y_yz_xxyz[k];

                g_y_y_yzz_xxz[k] = -g_y_y_yz_xxz[k] * ab_z + g_y_y_yz_xxzz[k];

                g_y_y_yzz_xyy[k] = -g_y_y_yz_xyy[k] * ab_z + g_y_y_yz_xyyz[k];

                g_y_y_yzz_xyz[k] = -g_y_y_yz_xyz[k] * ab_z + g_y_y_yz_xyzz[k];

                g_y_y_yzz_xzz[k] = -g_y_y_yz_xzz[k] * ab_z + g_y_y_yz_xzzz[k];

                g_y_y_yzz_yyy[k] = -g_y_y_yz_yyy[k] * ab_z + g_y_y_yz_yyyz[k];

                g_y_y_yzz_yyz[k] = -g_y_y_yz_yyz[k] * ab_z + g_y_y_yz_yyzz[k];

                g_y_y_yzz_yzz[k] = -g_y_y_yz_yzz[k] * ab_z + g_y_y_yz_yzzz[k];

                g_y_y_yzz_zzz[k] = -g_y_y_yz_zzz[k] * ab_z + g_y_y_yz_zzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzz_xxx = cbuffer.data(ff_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_y_zzz_xxy = cbuffer.data(ff_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_y_zzz_xxz = cbuffer.data(ff_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_y_zzz_xyy = cbuffer.data(ff_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_y_zzz_xyz = cbuffer.data(ff_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_y_zzz_xzz = cbuffer.data(ff_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_y_zzz_yyy = cbuffer.data(ff_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_y_zzz_yyz = cbuffer.data(ff_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_y_zzz_yzz = cbuffer.data(ff_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_y_zzz_zzz = cbuffer.data(ff_geom_11_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zz_xxx, g_y_y_zz_xxxz, g_y_y_zz_xxy, g_y_y_zz_xxyz, g_y_y_zz_xxz, g_y_y_zz_xxzz, g_y_y_zz_xyy, g_y_y_zz_xyyz, g_y_y_zz_xyz, g_y_y_zz_xyzz, g_y_y_zz_xzz, g_y_y_zz_xzzz, g_y_y_zz_yyy, g_y_y_zz_yyyz, g_y_y_zz_yyz, g_y_y_zz_yyzz, g_y_y_zz_yzz, g_y_y_zz_yzzz, g_y_y_zz_zzz, g_y_y_zz_zzzz, g_y_y_zzz_xxx, g_y_y_zzz_xxy, g_y_y_zzz_xxz, g_y_y_zzz_xyy, g_y_y_zzz_xyz, g_y_y_zzz_xzz, g_y_y_zzz_yyy, g_y_y_zzz_yyz, g_y_y_zzz_yzz, g_y_y_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzz_xxx[k] = -g_y_y_zz_xxx[k] * ab_z + g_y_y_zz_xxxz[k];

                g_y_y_zzz_xxy[k] = -g_y_y_zz_xxy[k] * ab_z + g_y_y_zz_xxyz[k];

                g_y_y_zzz_xxz[k] = -g_y_y_zz_xxz[k] * ab_z + g_y_y_zz_xxzz[k];

                g_y_y_zzz_xyy[k] = -g_y_y_zz_xyy[k] * ab_z + g_y_y_zz_xyyz[k];

                g_y_y_zzz_xyz[k] = -g_y_y_zz_xyz[k] * ab_z + g_y_y_zz_xyzz[k];

                g_y_y_zzz_xzz[k] = -g_y_y_zz_xzz[k] * ab_z + g_y_y_zz_xzzz[k];

                g_y_y_zzz_yyy[k] = -g_y_y_zz_yyy[k] * ab_z + g_y_y_zz_yyyz[k];

                g_y_y_zzz_yyz[k] = -g_y_y_zz_yyz[k] * ab_z + g_y_y_zz_yyzz[k];

                g_y_y_zzz_yzz[k] = -g_y_y_zz_yzz[k] * ab_z + g_y_y_zz_yzzz[k];

                g_y_y_zzz_zzz[k] = -g_y_y_zz_zzz[k] * ab_z + g_y_y_zz_zzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxx_xxx = cbuffer.data(ff_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_xxx_xxy = cbuffer.data(ff_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_xxx_xxz = cbuffer.data(ff_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_xxx_xyy = cbuffer.data(ff_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_z_xxx_xyz = cbuffer.data(ff_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_z_xxx_xzz = cbuffer.data(ff_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_z_xxx_yyy = cbuffer.data(ff_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_z_xxx_yyz = cbuffer.data(ff_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_z_xxx_yzz = cbuffer.data(ff_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_z_xxx_zzz = cbuffer.data(ff_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xx_xxx, g_y_z_xx_xxxx, g_y_z_xx_xxxy, g_y_z_xx_xxxz, g_y_z_xx_xxy, g_y_z_xx_xxyy, g_y_z_xx_xxyz, g_y_z_xx_xxz, g_y_z_xx_xxzz, g_y_z_xx_xyy, g_y_z_xx_xyyy, g_y_z_xx_xyyz, g_y_z_xx_xyz, g_y_z_xx_xyzz, g_y_z_xx_xzz, g_y_z_xx_xzzz, g_y_z_xx_yyy, g_y_z_xx_yyz, g_y_z_xx_yzz, g_y_z_xx_zzz, g_y_z_xxx_xxx, g_y_z_xxx_xxy, g_y_z_xxx_xxz, g_y_z_xxx_xyy, g_y_z_xxx_xyz, g_y_z_xxx_xzz, g_y_z_xxx_yyy, g_y_z_xxx_yyz, g_y_z_xxx_yzz, g_y_z_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxx_xxx[k] = -g_y_z_xx_xxx[k] * ab_x + g_y_z_xx_xxxx[k];

                g_y_z_xxx_xxy[k] = -g_y_z_xx_xxy[k] * ab_x + g_y_z_xx_xxxy[k];

                g_y_z_xxx_xxz[k] = -g_y_z_xx_xxz[k] * ab_x + g_y_z_xx_xxxz[k];

                g_y_z_xxx_xyy[k] = -g_y_z_xx_xyy[k] * ab_x + g_y_z_xx_xxyy[k];

                g_y_z_xxx_xyz[k] = -g_y_z_xx_xyz[k] * ab_x + g_y_z_xx_xxyz[k];

                g_y_z_xxx_xzz[k] = -g_y_z_xx_xzz[k] * ab_x + g_y_z_xx_xxzz[k];

                g_y_z_xxx_yyy[k] = -g_y_z_xx_yyy[k] * ab_x + g_y_z_xx_xyyy[k];

                g_y_z_xxx_yyz[k] = -g_y_z_xx_yyz[k] * ab_x + g_y_z_xx_xyyz[k];

                g_y_z_xxx_yzz[k] = -g_y_z_xx_yzz[k] * ab_x + g_y_z_xx_xyzz[k];

                g_y_z_xxx_zzz[k] = -g_y_z_xx_zzz[k] * ab_x + g_y_z_xx_xzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxy_xxx = cbuffer.data(ff_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_z_xxy_xxy = cbuffer.data(ff_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_z_xxy_xxz = cbuffer.data(ff_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_z_xxy_xyy = cbuffer.data(ff_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_z_xxy_xyz = cbuffer.data(ff_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_z_xxy_xzz = cbuffer.data(ff_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_z_xxy_yyy = cbuffer.data(ff_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_z_xxy_yyz = cbuffer.data(ff_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_z_xxy_yzz = cbuffer.data(ff_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_z_xxy_zzz = cbuffer.data(ff_geom_11_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxy_xxx, g_y_z_xxy_xxy, g_y_z_xxy_xxz, g_y_z_xxy_xyy, g_y_z_xxy_xyz, g_y_z_xxy_xzz, g_y_z_xxy_yyy, g_y_z_xxy_yyz, g_y_z_xxy_yzz, g_y_z_xxy_zzz, g_y_z_xy_xxx, g_y_z_xy_xxxx, g_y_z_xy_xxxy, g_y_z_xy_xxxz, g_y_z_xy_xxy, g_y_z_xy_xxyy, g_y_z_xy_xxyz, g_y_z_xy_xxz, g_y_z_xy_xxzz, g_y_z_xy_xyy, g_y_z_xy_xyyy, g_y_z_xy_xyyz, g_y_z_xy_xyz, g_y_z_xy_xyzz, g_y_z_xy_xzz, g_y_z_xy_xzzz, g_y_z_xy_yyy, g_y_z_xy_yyz, g_y_z_xy_yzz, g_y_z_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxy_xxx[k] = -g_y_z_xy_xxx[k] * ab_x + g_y_z_xy_xxxx[k];

                g_y_z_xxy_xxy[k] = -g_y_z_xy_xxy[k] * ab_x + g_y_z_xy_xxxy[k];

                g_y_z_xxy_xxz[k] = -g_y_z_xy_xxz[k] * ab_x + g_y_z_xy_xxxz[k];

                g_y_z_xxy_xyy[k] = -g_y_z_xy_xyy[k] * ab_x + g_y_z_xy_xxyy[k];

                g_y_z_xxy_xyz[k] = -g_y_z_xy_xyz[k] * ab_x + g_y_z_xy_xxyz[k];

                g_y_z_xxy_xzz[k] = -g_y_z_xy_xzz[k] * ab_x + g_y_z_xy_xxzz[k];

                g_y_z_xxy_yyy[k] = -g_y_z_xy_yyy[k] * ab_x + g_y_z_xy_xyyy[k];

                g_y_z_xxy_yyz[k] = -g_y_z_xy_yyz[k] * ab_x + g_y_z_xy_xyyz[k];

                g_y_z_xxy_yzz[k] = -g_y_z_xy_yzz[k] * ab_x + g_y_z_xy_xyzz[k];

                g_y_z_xxy_zzz[k] = -g_y_z_xy_zzz[k] * ab_x + g_y_z_xy_xzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxz_xxx = cbuffer.data(ff_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_z_xxz_xxy = cbuffer.data(ff_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_z_xxz_xxz = cbuffer.data(ff_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_z_xxz_xyy = cbuffer.data(ff_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_z_xxz_xyz = cbuffer.data(ff_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_z_xxz_xzz = cbuffer.data(ff_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_z_xxz_yyy = cbuffer.data(ff_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_z_xxz_yyz = cbuffer.data(ff_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_z_xxz_yzz = cbuffer.data(ff_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_z_xxz_zzz = cbuffer.data(ff_geom_11_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxz_xxx, g_y_z_xxz_xxy, g_y_z_xxz_xxz, g_y_z_xxz_xyy, g_y_z_xxz_xyz, g_y_z_xxz_xzz, g_y_z_xxz_yyy, g_y_z_xxz_yyz, g_y_z_xxz_yzz, g_y_z_xxz_zzz, g_y_z_xz_xxx, g_y_z_xz_xxxx, g_y_z_xz_xxxy, g_y_z_xz_xxxz, g_y_z_xz_xxy, g_y_z_xz_xxyy, g_y_z_xz_xxyz, g_y_z_xz_xxz, g_y_z_xz_xxzz, g_y_z_xz_xyy, g_y_z_xz_xyyy, g_y_z_xz_xyyz, g_y_z_xz_xyz, g_y_z_xz_xyzz, g_y_z_xz_xzz, g_y_z_xz_xzzz, g_y_z_xz_yyy, g_y_z_xz_yyz, g_y_z_xz_yzz, g_y_z_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxz_xxx[k] = -g_y_z_xz_xxx[k] * ab_x + g_y_z_xz_xxxx[k];

                g_y_z_xxz_xxy[k] = -g_y_z_xz_xxy[k] * ab_x + g_y_z_xz_xxxy[k];

                g_y_z_xxz_xxz[k] = -g_y_z_xz_xxz[k] * ab_x + g_y_z_xz_xxxz[k];

                g_y_z_xxz_xyy[k] = -g_y_z_xz_xyy[k] * ab_x + g_y_z_xz_xxyy[k];

                g_y_z_xxz_xyz[k] = -g_y_z_xz_xyz[k] * ab_x + g_y_z_xz_xxyz[k];

                g_y_z_xxz_xzz[k] = -g_y_z_xz_xzz[k] * ab_x + g_y_z_xz_xxzz[k];

                g_y_z_xxz_yyy[k] = -g_y_z_xz_yyy[k] * ab_x + g_y_z_xz_xyyy[k];

                g_y_z_xxz_yyz[k] = -g_y_z_xz_yyz[k] * ab_x + g_y_z_xz_xyyz[k];

                g_y_z_xxz_yzz[k] = -g_y_z_xz_yzz[k] * ab_x + g_y_z_xz_xyzz[k];

                g_y_z_xxz_zzz[k] = -g_y_z_xz_zzz[k] * ab_x + g_y_z_xz_xzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyy_xxx = cbuffer.data(ff_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_z_xyy_xxy = cbuffer.data(ff_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_z_xyy_xxz = cbuffer.data(ff_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_z_xyy_xyy = cbuffer.data(ff_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_z_xyy_xyz = cbuffer.data(ff_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_z_xyy_xzz = cbuffer.data(ff_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_z_xyy_yyy = cbuffer.data(ff_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_z_xyy_yyz = cbuffer.data(ff_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_z_xyy_yzz = cbuffer.data(ff_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_z_xyy_zzz = cbuffer.data(ff_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyy_xxx, g_y_z_xyy_xxy, g_y_z_xyy_xxz, g_y_z_xyy_xyy, g_y_z_xyy_xyz, g_y_z_xyy_xzz, g_y_z_xyy_yyy, g_y_z_xyy_yyz, g_y_z_xyy_yzz, g_y_z_xyy_zzz, g_y_z_yy_xxx, g_y_z_yy_xxxx, g_y_z_yy_xxxy, g_y_z_yy_xxxz, g_y_z_yy_xxy, g_y_z_yy_xxyy, g_y_z_yy_xxyz, g_y_z_yy_xxz, g_y_z_yy_xxzz, g_y_z_yy_xyy, g_y_z_yy_xyyy, g_y_z_yy_xyyz, g_y_z_yy_xyz, g_y_z_yy_xyzz, g_y_z_yy_xzz, g_y_z_yy_xzzz, g_y_z_yy_yyy, g_y_z_yy_yyz, g_y_z_yy_yzz, g_y_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyy_xxx[k] = -g_y_z_yy_xxx[k] * ab_x + g_y_z_yy_xxxx[k];

                g_y_z_xyy_xxy[k] = -g_y_z_yy_xxy[k] * ab_x + g_y_z_yy_xxxy[k];

                g_y_z_xyy_xxz[k] = -g_y_z_yy_xxz[k] * ab_x + g_y_z_yy_xxxz[k];

                g_y_z_xyy_xyy[k] = -g_y_z_yy_xyy[k] * ab_x + g_y_z_yy_xxyy[k];

                g_y_z_xyy_xyz[k] = -g_y_z_yy_xyz[k] * ab_x + g_y_z_yy_xxyz[k];

                g_y_z_xyy_xzz[k] = -g_y_z_yy_xzz[k] * ab_x + g_y_z_yy_xxzz[k];

                g_y_z_xyy_yyy[k] = -g_y_z_yy_yyy[k] * ab_x + g_y_z_yy_xyyy[k];

                g_y_z_xyy_yyz[k] = -g_y_z_yy_yyz[k] * ab_x + g_y_z_yy_xyyz[k];

                g_y_z_xyy_yzz[k] = -g_y_z_yy_yzz[k] * ab_x + g_y_z_yy_xyzz[k];

                g_y_z_xyy_zzz[k] = -g_y_z_yy_zzz[k] * ab_x + g_y_z_yy_xzzz[k];
            }

            /// Set up 540-550 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyz_xxx = cbuffer.data(ff_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_z_xyz_xxy = cbuffer.data(ff_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_z_xyz_xxz = cbuffer.data(ff_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_z_xyz_xyy = cbuffer.data(ff_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_z_xyz_xyz = cbuffer.data(ff_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_z_xyz_xzz = cbuffer.data(ff_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_z_xyz_yyy = cbuffer.data(ff_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_z_xyz_yyz = cbuffer.data(ff_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_z_xyz_yzz = cbuffer.data(ff_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_z_xyz_zzz = cbuffer.data(ff_geom_11_off + 549 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyz_xxx, g_y_z_xyz_xxy, g_y_z_xyz_xxz, g_y_z_xyz_xyy, g_y_z_xyz_xyz, g_y_z_xyz_xzz, g_y_z_xyz_yyy, g_y_z_xyz_yyz, g_y_z_xyz_yzz, g_y_z_xyz_zzz, g_y_z_yz_xxx, g_y_z_yz_xxxx, g_y_z_yz_xxxy, g_y_z_yz_xxxz, g_y_z_yz_xxy, g_y_z_yz_xxyy, g_y_z_yz_xxyz, g_y_z_yz_xxz, g_y_z_yz_xxzz, g_y_z_yz_xyy, g_y_z_yz_xyyy, g_y_z_yz_xyyz, g_y_z_yz_xyz, g_y_z_yz_xyzz, g_y_z_yz_xzz, g_y_z_yz_xzzz, g_y_z_yz_yyy, g_y_z_yz_yyz, g_y_z_yz_yzz, g_y_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyz_xxx[k] = -g_y_z_yz_xxx[k] * ab_x + g_y_z_yz_xxxx[k];

                g_y_z_xyz_xxy[k] = -g_y_z_yz_xxy[k] * ab_x + g_y_z_yz_xxxy[k];

                g_y_z_xyz_xxz[k] = -g_y_z_yz_xxz[k] * ab_x + g_y_z_yz_xxxz[k];

                g_y_z_xyz_xyy[k] = -g_y_z_yz_xyy[k] * ab_x + g_y_z_yz_xxyy[k];

                g_y_z_xyz_xyz[k] = -g_y_z_yz_xyz[k] * ab_x + g_y_z_yz_xxyz[k];

                g_y_z_xyz_xzz[k] = -g_y_z_yz_xzz[k] * ab_x + g_y_z_yz_xxzz[k];

                g_y_z_xyz_yyy[k] = -g_y_z_yz_yyy[k] * ab_x + g_y_z_yz_xyyy[k];

                g_y_z_xyz_yyz[k] = -g_y_z_yz_yyz[k] * ab_x + g_y_z_yz_xyyz[k];

                g_y_z_xyz_yzz[k] = -g_y_z_yz_yzz[k] * ab_x + g_y_z_yz_xyzz[k];

                g_y_z_xyz_zzz[k] = -g_y_z_yz_zzz[k] * ab_x + g_y_z_yz_xzzz[k];
            }

            /// Set up 550-560 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzz_xxx = cbuffer.data(ff_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_z_xzz_xxy = cbuffer.data(ff_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_z_xzz_xxz = cbuffer.data(ff_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_z_xzz_xyy = cbuffer.data(ff_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_z_xzz_xyz = cbuffer.data(ff_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_z_xzz_xzz = cbuffer.data(ff_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_z_xzz_yyy = cbuffer.data(ff_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_z_xzz_yyz = cbuffer.data(ff_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_z_xzz_yzz = cbuffer.data(ff_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_z_xzz_zzz = cbuffer.data(ff_geom_11_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzz_xxx, g_y_z_xzz_xxy, g_y_z_xzz_xxz, g_y_z_xzz_xyy, g_y_z_xzz_xyz, g_y_z_xzz_xzz, g_y_z_xzz_yyy, g_y_z_xzz_yyz, g_y_z_xzz_yzz, g_y_z_xzz_zzz, g_y_z_zz_xxx, g_y_z_zz_xxxx, g_y_z_zz_xxxy, g_y_z_zz_xxxz, g_y_z_zz_xxy, g_y_z_zz_xxyy, g_y_z_zz_xxyz, g_y_z_zz_xxz, g_y_z_zz_xxzz, g_y_z_zz_xyy, g_y_z_zz_xyyy, g_y_z_zz_xyyz, g_y_z_zz_xyz, g_y_z_zz_xyzz, g_y_z_zz_xzz, g_y_z_zz_xzzz, g_y_z_zz_yyy, g_y_z_zz_yyz, g_y_z_zz_yzz, g_y_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzz_xxx[k] = -g_y_z_zz_xxx[k] * ab_x + g_y_z_zz_xxxx[k];

                g_y_z_xzz_xxy[k] = -g_y_z_zz_xxy[k] * ab_x + g_y_z_zz_xxxy[k];

                g_y_z_xzz_xxz[k] = -g_y_z_zz_xxz[k] * ab_x + g_y_z_zz_xxxz[k];

                g_y_z_xzz_xyy[k] = -g_y_z_zz_xyy[k] * ab_x + g_y_z_zz_xxyy[k];

                g_y_z_xzz_xyz[k] = -g_y_z_zz_xyz[k] * ab_x + g_y_z_zz_xxyz[k];

                g_y_z_xzz_xzz[k] = -g_y_z_zz_xzz[k] * ab_x + g_y_z_zz_xxzz[k];

                g_y_z_xzz_yyy[k] = -g_y_z_zz_yyy[k] * ab_x + g_y_z_zz_xyyy[k];

                g_y_z_xzz_yyz[k] = -g_y_z_zz_yyz[k] * ab_x + g_y_z_zz_xyyz[k];

                g_y_z_xzz_yzz[k] = -g_y_z_zz_yzz[k] * ab_x + g_y_z_zz_xyzz[k];

                g_y_z_xzz_zzz[k] = -g_y_z_zz_zzz[k] * ab_x + g_y_z_zz_xzzz[k];
            }

            /// Set up 560-570 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyy_xxx = cbuffer.data(ff_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_z_yyy_xxy = cbuffer.data(ff_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_z_yyy_xxz = cbuffer.data(ff_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_z_yyy_xyy = cbuffer.data(ff_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_z_yyy_xyz = cbuffer.data(ff_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_z_yyy_xzz = cbuffer.data(ff_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_z_yyy_yyy = cbuffer.data(ff_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_z_yyy_yyz = cbuffer.data(ff_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_z_yyy_yzz = cbuffer.data(ff_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_z_yyy_zzz = cbuffer.data(ff_geom_11_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_xxx, g_0_z_yy_xxy, g_0_z_yy_xxz, g_0_z_yy_xyy, g_0_z_yy_xyz, g_0_z_yy_xzz, g_0_z_yy_yyy, g_0_z_yy_yyz, g_0_z_yy_yzz, g_0_z_yy_zzz, g_y_z_yy_xxx, g_y_z_yy_xxxy, g_y_z_yy_xxy, g_y_z_yy_xxyy, g_y_z_yy_xxyz, g_y_z_yy_xxz, g_y_z_yy_xyy, g_y_z_yy_xyyy, g_y_z_yy_xyyz, g_y_z_yy_xyz, g_y_z_yy_xyzz, g_y_z_yy_xzz, g_y_z_yy_yyy, g_y_z_yy_yyyy, g_y_z_yy_yyyz, g_y_z_yy_yyz, g_y_z_yy_yyzz, g_y_z_yy_yzz, g_y_z_yy_yzzz, g_y_z_yy_zzz, g_y_z_yyy_xxx, g_y_z_yyy_xxy, g_y_z_yyy_xxz, g_y_z_yyy_xyy, g_y_z_yyy_xyz, g_y_z_yyy_xzz, g_y_z_yyy_yyy, g_y_z_yyy_yyz, g_y_z_yyy_yzz, g_y_z_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyy_xxx[k] = -g_0_z_yy_xxx[k] - g_y_z_yy_xxx[k] * ab_y + g_y_z_yy_xxxy[k];

                g_y_z_yyy_xxy[k] = -g_0_z_yy_xxy[k] - g_y_z_yy_xxy[k] * ab_y + g_y_z_yy_xxyy[k];

                g_y_z_yyy_xxz[k] = -g_0_z_yy_xxz[k] - g_y_z_yy_xxz[k] * ab_y + g_y_z_yy_xxyz[k];

                g_y_z_yyy_xyy[k] = -g_0_z_yy_xyy[k] - g_y_z_yy_xyy[k] * ab_y + g_y_z_yy_xyyy[k];

                g_y_z_yyy_xyz[k] = -g_0_z_yy_xyz[k] - g_y_z_yy_xyz[k] * ab_y + g_y_z_yy_xyyz[k];

                g_y_z_yyy_xzz[k] = -g_0_z_yy_xzz[k] - g_y_z_yy_xzz[k] * ab_y + g_y_z_yy_xyzz[k];

                g_y_z_yyy_yyy[k] = -g_0_z_yy_yyy[k] - g_y_z_yy_yyy[k] * ab_y + g_y_z_yy_yyyy[k];

                g_y_z_yyy_yyz[k] = -g_0_z_yy_yyz[k] - g_y_z_yy_yyz[k] * ab_y + g_y_z_yy_yyyz[k];

                g_y_z_yyy_yzz[k] = -g_0_z_yy_yzz[k] - g_y_z_yy_yzz[k] * ab_y + g_y_z_yy_yyzz[k];

                g_y_z_yyy_zzz[k] = -g_0_z_yy_zzz[k] - g_y_z_yy_zzz[k] * ab_y + g_y_z_yy_yzzz[k];
            }

            /// Set up 570-580 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyz_xxx = cbuffer.data(ff_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_z_yyz_xxy = cbuffer.data(ff_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_z_yyz_xxz = cbuffer.data(ff_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_z_yyz_xyy = cbuffer.data(ff_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_z_yyz_xyz = cbuffer.data(ff_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_z_yyz_xzz = cbuffer.data(ff_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_z_yyz_yyy = cbuffer.data(ff_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_z_yyz_yyz = cbuffer.data(ff_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_z_yyz_yzz = cbuffer.data(ff_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_z_yyz_zzz = cbuffer.data(ff_geom_11_off + 579 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_xxx, g_0_z_yz_xxy, g_0_z_yz_xxz, g_0_z_yz_xyy, g_0_z_yz_xyz, g_0_z_yz_xzz, g_0_z_yz_yyy, g_0_z_yz_yyz, g_0_z_yz_yzz, g_0_z_yz_zzz, g_y_z_yyz_xxx, g_y_z_yyz_xxy, g_y_z_yyz_xxz, g_y_z_yyz_xyy, g_y_z_yyz_xyz, g_y_z_yyz_xzz, g_y_z_yyz_yyy, g_y_z_yyz_yyz, g_y_z_yyz_yzz, g_y_z_yyz_zzz, g_y_z_yz_xxx, g_y_z_yz_xxxy, g_y_z_yz_xxy, g_y_z_yz_xxyy, g_y_z_yz_xxyz, g_y_z_yz_xxz, g_y_z_yz_xyy, g_y_z_yz_xyyy, g_y_z_yz_xyyz, g_y_z_yz_xyz, g_y_z_yz_xyzz, g_y_z_yz_xzz, g_y_z_yz_yyy, g_y_z_yz_yyyy, g_y_z_yz_yyyz, g_y_z_yz_yyz, g_y_z_yz_yyzz, g_y_z_yz_yzz, g_y_z_yz_yzzz, g_y_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyz_xxx[k] = -g_0_z_yz_xxx[k] - g_y_z_yz_xxx[k] * ab_y + g_y_z_yz_xxxy[k];

                g_y_z_yyz_xxy[k] = -g_0_z_yz_xxy[k] - g_y_z_yz_xxy[k] * ab_y + g_y_z_yz_xxyy[k];

                g_y_z_yyz_xxz[k] = -g_0_z_yz_xxz[k] - g_y_z_yz_xxz[k] * ab_y + g_y_z_yz_xxyz[k];

                g_y_z_yyz_xyy[k] = -g_0_z_yz_xyy[k] - g_y_z_yz_xyy[k] * ab_y + g_y_z_yz_xyyy[k];

                g_y_z_yyz_xyz[k] = -g_0_z_yz_xyz[k] - g_y_z_yz_xyz[k] * ab_y + g_y_z_yz_xyyz[k];

                g_y_z_yyz_xzz[k] = -g_0_z_yz_xzz[k] - g_y_z_yz_xzz[k] * ab_y + g_y_z_yz_xyzz[k];

                g_y_z_yyz_yyy[k] = -g_0_z_yz_yyy[k] - g_y_z_yz_yyy[k] * ab_y + g_y_z_yz_yyyy[k];

                g_y_z_yyz_yyz[k] = -g_0_z_yz_yyz[k] - g_y_z_yz_yyz[k] * ab_y + g_y_z_yz_yyyz[k];

                g_y_z_yyz_yzz[k] = -g_0_z_yz_yzz[k] - g_y_z_yz_yzz[k] * ab_y + g_y_z_yz_yyzz[k];

                g_y_z_yyz_zzz[k] = -g_0_z_yz_zzz[k] - g_y_z_yz_zzz[k] * ab_y + g_y_z_yz_yzzz[k];
            }

            /// Set up 580-590 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzz_xxx = cbuffer.data(ff_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_z_yzz_xxy = cbuffer.data(ff_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_z_yzz_xxz = cbuffer.data(ff_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_z_yzz_xyy = cbuffer.data(ff_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_z_yzz_xyz = cbuffer.data(ff_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_z_yzz_xzz = cbuffer.data(ff_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_z_yzz_yyy = cbuffer.data(ff_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_z_yzz_yyz = cbuffer.data(ff_geom_11_off + 587 * ccomps * dcomps);

            auto g_y_z_yzz_yzz = cbuffer.data(ff_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_z_yzz_zzz = cbuffer.data(ff_geom_11_off + 589 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxx, g_0_z_zz_xxy, g_0_z_zz_xxz, g_0_z_zz_xyy, g_0_z_zz_xyz, g_0_z_zz_xzz, g_0_z_zz_yyy, g_0_z_zz_yyz, g_0_z_zz_yzz, g_0_z_zz_zzz, g_y_z_yzz_xxx, g_y_z_yzz_xxy, g_y_z_yzz_xxz, g_y_z_yzz_xyy, g_y_z_yzz_xyz, g_y_z_yzz_xzz, g_y_z_yzz_yyy, g_y_z_yzz_yyz, g_y_z_yzz_yzz, g_y_z_yzz_zzz, g_y_z_zz_xxx, g_y_z_zz_xxxy, g_y_z_zz_xxy, g_y_z_zz_xxyy, g_y_z_zz_xxyz, g_y_z_zz_xxz, g_y_z_zz_xyy, g_y_z_zz_xyyy, g_y_z_zz_xyyz, g_y_z_zz_xyz, g_y_z_zz_xyzz, g_y_z_zz_xzz, g_y_z_zz_yyy, g_y_z_zz_yyyy, g_y_z_zz_yyyz, g_y_z_zz_yyz, g_y_z_zz_yyzz, g_y_z_zz_yzz, g_y_z_zz_yzzz, g_y_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzz_xxx[k] = -g_0_z_zz_xxx[k] - g_y_z_zz_xxx[k] * ab_y + g_y_z_zz_xxxy[k];

                g_y_z_yzz_xxy[k] = -g_0_z_zz_xxy[k] - g_y_z_zz_xxy[k] * ab_y + g_y_z_zz_xxyy[k];

                g_y_z_yzz_xxz[k] = -g_0_z_zz_xxz[k] - g_y_z_zz_xxz[k] * ab_y + g_y_z_zz_xxyz[k];

                g_y_z_yzz_xyy[k] = -g_0_z_zz_xyy[k] - g_y_z_zz_xyy[k] * ab_y + g_y_z_zz_xyyy[k];

                g_y_z_yzz_xyz[k] = -g_0_z_zz_xyz[k] - g_y_z_zz_xyz[k] * ab_y + g_y_z_zz_xyyz[k];

                g_y_z_yzz_xzz[k] = -g_0_z_zz_xzz[k] - g_y_z_zz_xzz[k] * ab_y + g_y_z_zz_xyzz[k];

                g_y_z_yzz_yyy[k] = -g_0_z_zz_yyy[k] - g_y_z_zz_yyy[k] * ab_y + g_y_z_zz_yyyy[k];

                g_y_z_yzz_yyz[k] = -g_0_z_zz_yyz[k] - g_y_z_zz_yyz[k] * ab_y + g_y_z_zz_yyyz[k];

                g_y_z_yzz_yzz[k] = -g_0_z_zz_yzz[k] - g_y_z_zz_yzz[k] * ab_y + g_y_z_zz_yyzz[k];

                g_y_z_yzz_zzz[k] = -g_0_z_zz_zzz[k] - g_y_z_zz_zzz[k] * ab_y + g_y_z_zz_yzzz[k];
            }

            /// Set up 590-600 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzz_xxx = cbuffer.data(ff_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_z_zzz_xxy = cbuffer.data(ff_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_z_zzz_xxz = cbuffer.data(ff_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_z_zzz_xyy = cbuffer.data(ff_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_z_zzz_xyz = cbuffer.data(ff_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_z_zzz_xzz = cbuffer.data(ff_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_z_zzz_yyy = cbuffer.data(ff_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_z_zzz_yyz = cbuffer.data(ff_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_z_zzz_yzz = cbuffer.data(ff_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_z_zzz_zzz = cbuffer.data(ff_geom_11_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xxx, g_y_0_zz_xxy, g_y_0_zz_xxz, g_y_0_zz_xyy, g_y_0_zz_xyz, g_y_0_zz_xzz, g_y_0_zz_yyy, g_y_0_zz_yyz, g_y_0_zz_yzz, g_y_0_zz_zzz, g_y_z_zz_xxx, g_y_z_zz_xxxz, g_y_z_zz_xxy, g_y_z_zz_xxyz, g_y_z_zz_xxz, g_y_z_zz_xxzz, g_y_z_zz_xyy, g_y_z_zz_xyyz, g_y_z_zz_xyz, g_y_z_zz_xyzz, g_y_z_zz_xzz, g_y_z_zz_xzzz, g_y_z_zz_yyy, g_y_z_zz_yyyz, g_y_z_zz_yyz, g_y_z_zz_yyzz, g_y_z_zz_yzz, g_y_z_zz_yzzz, g_y_z_zz_zzz, g_y_z_zz_zzzz, g_y_z_zzz_xxx, g_y_z_zzz_xxy, g_y_z_zzz_xxz, g_y_z_zzz_xyy, g_y_z_zzz_xyz, g_y_z_zzz_xzz, g_y_z_zzz_yyy, g_y_z_zzz_yyz, g_y_z_zzz_yzz, g_y_z_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzz_xxx[k] = g_y_0_zz_xxx[k] - g_y_z_zz_xxx[k] * ab_z + g_y_z_zz_xxxz[k];

                g_y_z_zzz_xxy[k] = g_y_0_zz_xxy[k] - g_y_z_zz_xxy[k] * ab_z + g_y_z_zz_xxyz[k];

                g_y_z_zzz_xxz[k] = g_y_0_zz_xxz[k] - g_y_z_zz_xxz[k] * ab_z + g_y_z_zz_xxzz[k];

                g_y_z_zzz_xyy[k] = g_y_0_zz_xyy[k] - g_y_z_zz_xyy[k] * ab_z + g_y_z_zz_xyyz[k];

                g_y_z_zzz_xyz[k] = g_y_0_zz_xyz[k] - g_y_z_zz_xyz[k] * ab_z + g_y_z_zz_xyzz[k];

                g_y_z_zzz_xzz[k] = g_y_0_zz_xzz[k] - g_y_z_zz_xzz[k] * ab_z + g_y_z_zz_xzzz[k];

                g_y_z_zzz_yyy[k] = g_y_0_zz_yyy[k] - g_y_z_zz_yyy[k] * ab_z + g_y_z_zz_yyyz[k];

                g_y_z_zzz_yyz[k] = g_y_0_zz_yyz[k] - g_y_z_zz_yyz[k] * ab_z + g_y_z_zz_yyzz[k];

                g_y_z_zzz_yzz[k] = g_y_0_zz_yzz[k] - g_y_z_zz_yzz[k] * ab_z + g_y_z_zz_yzzz[k];

                g_y_z_zzz_zzz[k] = g_y_0_zz_zzz[k] - g_y_z_zz_zzz[k] * ab_z + g_y_z_zz_zzzz[k];
            }

            /// Set up 600-610 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxx_xxx = cbuffer.data(ff_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_x_xxx_xxy = cbuffer.data(ff_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_x_xxx_xxz = cbuffer.data(ff_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_x_xxx_xyy = cbuffer.data(ff_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_x_xxx_xyz = cbuffer.data(ff_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_x_xxx_xzz = cbuffer.data(ff_geom_11_off + 605 * ccomps * dcomps);

            auto g_z_x_xxx_yyy = cbuffer.data(ff_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_x_xxx_yyz = cbuffer.data(ff_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_x_xxx_yzz = cbuffer.data(ff_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_x_xxx_zzz = cbuffer.data(ff_geom_11_off + 609 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_xxx, g_z_0_xx_xxy, g_z_0_xx_xxz, g_z_0_xx_xyy, g_z_0_xx_xyz, g_z_0_xx_xzz, g_z_0_xx_yyy, g_z_0_xx_yyz, g_z_0_xx_yzz, g_z_0_xx_zzz, g_z_x_xx_xxx, g_z_x_xx_xxxx, g_z_x_xx_xxxy, g_z_x_xx_xxxz, g_z_x_xx_xxy, g_z_x_xx_xxyy, g_z_x_xx_xxyz, g_z_x_xx_xxz, g_z_x_xx_xxzz, g_z_x_xx_xyy, g_z_x_xx_xyyy, g_z_x_xx_xyyz, g_z_x_xx_xyz, g_z_x_xx_xyzz, g_z_x_xx_xzz, g_z_x_xx_xzzz, g_z_x_xx_yyy, g_z_x_xx_yyz, g_z_x_xx_yzz, g_z_x_xx_zzz, g_z_x_xxx_xxx, g_z_x_xxx_xxy, g_z_x_xxx_xxz, g_z_x_xxx_xyy, g_z_x_xxx_xyz, g_z_x_xxx_xzz, g_z_x_xxx_yyy, g_z_x_xxx_yyz, g_z_x_xxx_yzz, g_z_x_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxx_xxx[k] = g_z_0_xx_xxx[k] - g_z_x_xx_xxx[k] * ab_x + g_z_x_xx_xxxx[k];

                g_z_x_xxx_xxy[k] = g_z_0_xx_xxy[k] - g_z_x_xx_xxy[k] * ab_x + g_z_x_xx_xxxy[k];

                g_z_x_xxx_xxz[k] = g_z_0_xx_xxz[k] - g_z_x_xx_xxz[k] * ab_x + g_z_x_xx_xxxz[k];

                g_z_x_xxx_xyy[k] = g_z_0_xx_xyy[k] - g_z_x_xx_xyy[k] * ab_x + g_z_x_xx_xxyy[k];

                g_z_x_xxx_xyz[k] = g_z_0_xx_xyz[k] - g_z_x_xx_xyz[k] * ab_x + g_z_x_xx_xxyz[k];

                g_z_x_xxx_xzz[k] = g_z_0_xx_xzz[k] - g_z_x_xx_xzz[k] * ab_x + g_z_x_xx_xxzz[k];

                g_z_x_xxx_yyy[k] = g_z_0_xx_yyy[k] - g_z_x_xx_yyy[k] * ab_x + g_z_x_xx_xyyy[k];

                g_z_x_xxx_yyz[k] = g_z_0_xx_yyz[k] - g_z_x_xx_yyz[k] * ab_x + g_z_x_xx_xyyz[k];

                g_z_x_xxx_yzz[k] = g_z_0_xx_yzz[k] - g_z_x_xx_yzz[k] * ab_x + g_z_x_xx_xyzz[k];

                g_z_x_xxx_zzz[k] = g_z_0_xx_zzz[k] - g_z_x_xx_zzz[k] * ab_x + g_z_x_xx_xzzz[k];
            }

            /// Set up 610-620 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxy_xxx = cbuffer.data(ff_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_x_xxy_xxy = cbuffer.data(ff_geom_11_off + 611 * ccomps * dcomps);

            auto g_z_x_xxy_xxz = cbuffer.data(ff_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_x_xxy_xyy = cbuffer.data(ff_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_x_xxy_xyz = cbuffer.data(ff_geom_11_off + 614 * ccomps * dcomps);

            auto g_z_x_xxy_xzz = cbuffer.data(ff_geom_11_off + 615 * ccomps * dcomps);

            auto g_z_x_xxy_yyy = cbuffer.data(ff_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_x_xxy_yyz = cbuffer.data(ff_geom_11_off + 617 * ccomps * dcomps);

            auto g_z_x_xxy_yzz = cbuffer.data(ff_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_x_xxy_zzz = cbuffer.data(ff_geom_11_off + 619 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xx_xxx, g_z_x_xx_xxxy, g_z_x_xx_xxy, g_z_x_xx_xxyy, g_z_x_xx_xxyz, g_z_x_xx_xxz, g_z_x_xx_xyy, g_z_x_xx_xyyy, g_z_x_xx_xyyz, g_z_x_xx_xyz, g_z_x_xx_xyzz, g_z_x_xx_xzz, g_z_x_xx_yyy, g_z_x_xx_yyyy, g_z_x_xx_yyyz, g_z_x_xx_yyz, g_z_x_xx_yyzz, g_z_x_xx_yzz, g_z_x_xx_yzzz, g_z_x_xx_zzz, g_z_x_xxy_xxx, g_z_x_xxy_xxy, g_z_x_xxy_xxz, g_z_x_xxy_xyy, g_z_x_xxy_xyz, g_z_x_xxy_xzz, g_z_x_xxy_yyy, g_z_x_xxy_yyz, g_z_x_xxy_yzz, g_z_x_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxy_xxx[k] = -g_z_x_xx_xxx[k] * ab_y + g_z_x_xx_xxxy[k];

                g_z_x_xxy_xxy[k] = -g_z_x_xx_xxy[k] * ab_y + g_z_x_xx_xxyy[k];

                g_z_x_xxy_xxz[k] = -g_z_x_xx_xxz[k] * ab_y + g_z_x_xx_xxyz[k];

                g_z_x_xxy_xyy[k] = -g_z_x_xx_xyy[k] * ab_y + g_z_x_xx_xyyy[k];

                g_z_x_xxy_xyz[k] = -g_z_x_xx_xyz[k] * ab_y + g_z_x_xx_xyyz[k];

                g_z_x_xxy_xzz[k] = -g_z_x_xx_xzz[k] * ab_y + g_z_x_xx_xyzz[k];

                g_z_x_xxy_yyy[k] = -g_z_x_xx_yyy[k] * ab_y + g_z_x_xx_yyyy[k];

                g_z_x_xxy_yyz[k] = -g_z_x_xx_yyz[k] * ab_y + g_z_x_xx_yyyz[k];

                g_z_x_xxy_yzz[k] = -g_z_x_xx_yzz[k] * ab_y + g_z_x_xx_yyzz[k];

                g_z_x_xxy_zzz[k] = -g_z_x_xx_zzz[k] * ab_y + g_z_x_xx_yzzz[k];
            }

            /// Set up 620-630 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxz_xxx = cbuffer.data(ff_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_x_xxz_xxy = cbuffer.data(ff_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_x_xxz_xxz = cbuffer.data(ff_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_x_xxz_xyy = cbuffer.data(ff_geom_11_off + 623 * ccomps * dcomps);

            auto g_z_x_xxz_xyz = cbuffer.data(ff_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_x_xxz_xzz = cbuffer.data(ff_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_x_xxz_yyy = cbuffer.data(ff_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_x_xxz_yyz = cbuffer.data(ff_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_x_xxz_yzz = cbuffer.data(ff_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_x_xxz_zzz = cbuffer.data(ff_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xxx, g_z_0_xz_xxy, g_z_0_xz_xxz, g_z_0_xz_xyy, g_z_0_xz_xyz, g_z_0_xz_xzz, g_z_0_xz_yyy, g_z_0_xz_yyz, g_z_0_xz_yzz, g_z_0_xz_zzz, g_z_x_xxz_xxx, g_z_x_xxz_xxy, g_z_x_xxz_xxz, g_z_x_xxz_xyy, g_z_x_xxz_xyz, g_z_x_xxz_xzz, g_z_x_xxz_yyy, g_z_x_xxz_yyz, g_z_x_xxz_yzz, g_z_x_xxz_zzz, g_z_x_xz_xxx, g_z_x_xz_xxxx, g_z_x_xz_xxxy, g_z_x_xz_xxxz, g_z_x_xz_xxy, g_z_x_xz_xxyy, g_z_x_xz_xxyz, g_z_x_xz_xxz, g_z_x_xz_xxzz, g_z_x_xz_xyy, g_z_x_xz_xyyy, g_z_x_xz_xyyz, g_z_x_xz_xyz, g_z_x_xz_xyzz, g_z_x_xz_xzz, g_z_x_xz_xzzz, g_z_x_xz_yyy, g_z_x_xz_yyz, g_z_x_xz_yzz, g_z_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxz_xxx[k] = g_z_0_xz_xxx[k] - g_z_x_xz_xxx[k] * ab_x + g_z_x_xz_xxxx[k];

                g_z_x_xxz_xxy[k] = g_z_0_xz_xxy[k] - g_z_x_xz_xxy[k] * ab_x + g_z_x_xz_xxxy[k];

                g_z_x_xxz_xxz[k] = g_z_0_xz_xxz[k] - g_z_x_xz_xxz[k] * ab_x + g_z_x_xz_xxxz[k];

                g_z_x_xxz_xyy[k] = g_z_0_xz_xyy[k] - g_z_x_xz_xyy[k] * ab_x + g_z_x_xz_xxyy[k];

                g_z_x_xxz_xyz[k] = g_z_0_xz_xyz[k] - g_z_x_xz_xyz[k] * ab_x + g_z_x_xz_xxyz[k];

                g_z_x_xxz_xzz[k] = g_z_0_xz_xzz[k] - g_z_x_xz_xzz[k] * ab_x + g_z_x_xz_xxzz[k];

                g_z_x_xxz_yyy[k] = g_z_0_xz_yyy[k] - g_z_x_xz_yyy[k] * ab_x + g_z_x_xz_xyyy[k];

                g_z_x_xxz_yyz[k] = g_z_0_xz_yyz[k] - g_z_x_xz_yyz[k] * ab_x + g_z_x_xz_xyyz[k];

                g_z_x_xxz_yzz[k] = g_z_0_xz_yzz[k] - g_z_x_xz_yzz[k] * ab_x + g_z_x_xz_xyzz[k];

                g_z_x_xxz_zzz[k] = g_z_0_xz_zzz[k] - g_z_x_xz_zzz[k] * ab_x + g_z_x_xz_xzzz[k];
            }

            /// Set up 630-640 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyy_xxx = cbuffer.data(ff_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_x_xyy_xxy = cbuffer.data(ff_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_x_xyy_xxz = cbuffer.data(ff_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_x_xyy_xyy = cbuffer.data(ff_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_x_xyy_xyz = cbuffer.data(ff_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_x_xyy_xzz = cbuffer.data(ff_geom_11_off + 635 * ccomps * dcomps);

            auto g_z_x_xyy_yyy = cbuffer.data(ff_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_x_xyy_yyz = cbuffer.data(ff_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_x_xyy_yzz = cbuffer.data(ff_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_x_xyy_zzz = cbuffer.data(ff_geom_11_off + 639 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xy_xxx, g_z_x_xy_xxxy, g_z_x_xy_xxy, g_z_x_xy_xxyy, g_z_x_xy_xxyz, g_z_x_xy_xxz, g_z_x_xy_xyy, g_z_x_xy_xyyy, g_z_x_xy_xyyz, g_z_x_xy_xyz, g_z_x_xy_xyzz, g_z_x_xy_xzz, g_z_x_xy_yyy, g_z_x_xy_yyyy, g_z_x_xy_yyyz, g_z_x_xy_yyz, g_z_x_xy_yyzz, g_z_x_xy_yzz, g_z_x_xy_yzzz, g_z_x_xy_zzz, g_z_x_xyy_xxx, g_z_x_xyy_xxy, g_z_x_xyy_xxz, g_z_x_xyy_xyy, g_z_x_xyy_xyz, g_z_x_xyy_xzz, g_z_x_xyy_yyy, g_z_x_xyy_yyz, g_z_x_xyy_yzz, g_z_x_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyy_xxx[k] = -g_z_x_xy_xxx[k] * ab_y + g_z_x_xy_xxxy[k];

                g_z_x_xyy_xxy[k] = -g_z_x_xy_xxy[k] * ab_y + g_z_x_xy_xxyy[k];

                g_z_x_xyy_xxz[k] = -g_z_x_xy_xxz[k] * ab_y + g_z_x_xy_xxyz[k];

                g_z_x_xyy_xyy[k] = -g_z_x_xy_xyy[k] * ab_y + g_z_x_xy_xyyy[k];

                g_z_x_xyy_xyz[k] = -g_z_x_xy_xyz[k] * ab_y + g_z_x_xy_xyyz[k];

                g_z_x_xyy_xzz[k] = -g_z_x_xy_xzz[k] * ab_y + g_z_x_xy_xyzz[k];

                g_z_x_xyy_yyy[k] = -g_z_x_xy_yyy[k] * ab_y + g_z_x_xy_yyyy[k];

                g_z_x_xyy_yyz[k] = -g_z_x_xy_yyz[k] * ab_y + g_z_x_xy_yyyz[k];

                g_z_x_xyy_yzz[k] = -g_z_x_xy_yzz[k] * ab_y + g_z_x_xy_yyzz[k];

                g_z_x_xyy_zzz[k] = -g_z_x_xy_zzz[k] * ab_y + g_z_x_xy_yzzz[k];
            }

            /// Set up 640-650 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyz_xxx = cbuffer.data(ff_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_x_xyz_xxy = cbuffer.data(ff_geom_11_off + 641 * ccomps * dcomps);

            auto g_z_x_xyz_xxz = cbuffer.data(ff_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_x_xyz_xyy = cbuffer.data(ff_geom_11_off + 643 * ccomps * dcomps);

            auto g_z_x_xyz_xyz = cbuffer.data(ff_geom_11_off + 644 * ccomps * dcomps);

            auto g_z_x_xyz_xzz = cbuffer.data(ff_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_x_xyz_yyy = cbuffer.data(ff_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_x_xyz_yyz = cbuffer.data(ff_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_x_xyz_yzz = cbuffer.data(ff_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_x_xyz_zzz = cbuffer.data(ff_geom_11_off + 649 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyz_xxx, g_z_x_xyz_xxy, g_z_x_xyz_xxz, g_z_x_xyz_xyy, g_z_x_xyz_xyz, g_z_x_xyz_xzz, g_z_x_xyz_yyy, g_z_x_xyz_yyz, g_z_x_xyz_yzz, g_z_x_xyz_zzz, g_z_x_xz_xxx, g_z_x_xz_xxxy, g_z_x_xz_xxy, g_z_x_xz_xxyy, g_z_x_xz_xxyz, g_z_x_xz_xxz, g_z_x_xz_xyy, g_z_x_xz_xyyy, g_z_x_xz_xyyz, g_z_x_xz_xyz, g_z_x_xz_xyzz, g_z_x_xz_xzz, g_z_x_xz_yyy, g_z_x_xz_yyyy, g_z_x_xz_yyyz, g_z_x_xz_yyz, g_z_x_xz_yyzz, g_z_x_xz_yzz, g_z_x_xz_yzzz, g_z_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyz_xxx[k] = -g_z_x_xz_xxx[k] * ab_y + g_z_x_xz_xxxy[k];

                g_z_x_xyz_xxy[k] = -g_z_x_xz_xxy[k] * ab_y + g_z_x_xz_xxyy[k];

                g_z_x_xyz_xxz[k] = -g_z_x_xz_xxz[k] * ab_y + g_z_x_xz_xxyz[k];

                g_z_x_xyz_xyy[k] = -g_z_x_xz_xyy[k] * ab_y + g_z_x_xz_xyyy[k];

                g_z_x_xyz_xyz[k] = -g_z_x_xz_xyz[k] * ab_y + g_z_x_xz_xyyz[k];

                g_z_x_xyz_xzz[k] = -g_z_x_xz_xzz[k] * ab_y + g_z_x_xz_xyzz[k];

                g_z_x_xyz_yyy[k] = -g_z_x_xz_yyy[k] * ab_y + g_z_x_xz_yyyy[k];

                g_z_x_xyz_yyz[k] = -g_z_x_xz_yyz[k] * ab_y + g_z_x_xz_yyyz[k];

                g_z_x_xyz_yzz[k] = -g_z_x_xz_yzz[k] * ab_y + g_z_x_xz_yyzz[k];

                g_z_x_xyz_zzz[k] = -g_z_x_xz_zzz[k] * ab_y + g_z_x_xz_yzzz[k];
            }

            /// Set up 650-660 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzz_xxx = cbuffer.data(ff_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_x_xzz_xxy = cbuffer.data(ff_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_x_xzz_xxz = cbuffer.data(ff_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_x_xzz_xyy = cbuffer.data(ff_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_x_xzz_xyz = cbuffer.data(ff_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_x_xzz_xzz = cbuffer.data(ff_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_x_xzz_yyy = cbuffer.data(ff_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_x_xzz_yyz = cbuffer.data(ff_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_x_xzz_yzz = cbuffer.data(ff_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_x_xzz_zzz = cbuffer.data(ff_geom_11_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz, g_z_x_xzz_xxx, g_z_x_xzz_xxy, g_z_x_xzz_xxz, g_z_x_xzz_xyy, g_z_x_xzz_xyz, g_z_x_xzz_xzz, g_z_x_xzz_yyy, g_z_x_xzz_yyz, g_z_x_xzz_yzz, g_z_x_xzz_zzz, g_z_x_zz_xxx, g_z_x_zz_xxxx, g_z_x_zz_xxxy, g_z_x_zz_xxxz, g_z_x_zz_xxy, g_z_x_zz_xxyy, g_z_x_zz_xxyz, g_z_x_zz_xxz, g_z_x_zz_xxzz, g_z_x_zz_xyy, g_z_x_zz_xyyy, g_z_x_zz_xyyz, g_z_x_zz_xyz, g_z_x_zz_xyzz, g_z_x_zz_xzz, g_z_x_zz_xzzz, g_z_x_zz_yyy, g_z_x_zz_yyz, g_z_x_zz_yzz, g_z_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzz_xxx[k] = g_z_0_zz_xxx[k] - g_z_x_zz_xxx[k] * ab_x + g_z_x_zz_xxxx[k];

                g_z_x_xzz_xxy[k] = g_z_0_zz_xxy[k] - g_z_x_zz_xxy[k] * ab_x + g_z_x_zz_xxxy[k];

                g_z_x_xzz_xxz[k] = g_z_0_zz_xxz[k] - g_z_x_zz_xxz[k] * ab_x + g_z_x_zz_xxxz[k];

                g_z_x_xzz_xyy[k] = g_z_0_zz_xyy[k] - g_z_x_zz_xyy[k] * ab_x + g_z_x_zz_xxyy[k];

                g_z_x_xzz_xyz[k] = g_z_0_zz_xyz[k] - g_z_x_zz_xyz[k] * ab_x + g_z_x_zz_xxyz[k];

                g_z_x_xzz_xzz[k] = g_z_0_zz_xzz[k] - g_z_x_zz_xzz[k] * ab_x + g_z_x_zz_xxzz[k];

                g_z_x_xzz_yyy[k] = g_z_0_zz_yyy[k] - g_z_x_zz_yyy[k] * ab_x + g_z_x_zz_xyyy[k];

                g_z_x_xzz_yyz[k] = g_z_0_zz_yyz[k] - g_z_x_zz_yyz[k] * ab_x + g_z_x_zz_xyyz[k];

                g_z_x_xzz_yzz[k] = g_z_0_zz_yzz[k] - g_z_x_zz_yzz[k] * ab_x + g_z_x_zz_xyzz[k];

                g_z_x_xzz_zzz[k] = g_z_0_zz_zzz[k] - g_z_x_zz_zzz[k] * ab_x + g_z_x_zz_xzzz[k];
            }

            /// Set up 660-670 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyy_xxx = cbuffer.data(ff_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_x_yyy_xxy = cbuffer.data(ff_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_x_yyy_xxz = cbuffer.data(ff_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_x_yyy_xyy = cbuffer.data(ff_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_x_yyy_xyz = cbuffer.data(ff_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_x_yyy_xzz = cbuffer.data(ff_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_x_yyy_yyy = cbuffer.data(ff_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_x_yyy_yyz = cbuffer.data(ff_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_x_yyy_yzz = cbuffer.data(ff_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_x_yyy_zzz = cbuffer.data(ff_geom_11_off + 669 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yy_xxx, g_z_x_yy_xxxy, g_z_x_yy_xxy, g_z_x_yy_xxyy, g_z_x_yy_xxyz, g_z_x_yy_xxz, g_z_x_yy_xyy, g_z_x_yy_xyyy, g_z_x_yy_xyyz, g_z_x_yy_xyz, g_z_x_yy_xyzz, g_z_x_yy_xzz, g_z_x_yy_yyy, g_z_x_yy_yyyy, g_z_x_yy_yyyz, g_z_x_yy_yyz, g_z_x_yy_yyzz, g_z_x_yy_yzz, g_z_x_yy_yzzz, g_z_x_yy_zzz, g_z_x_yyy_xxx, g_z_x_yyy_xxy, g_z_x_yyy_xxz, g_z_x_yyy_xyy, g_z_x_yyy_xyz, g_z_x_yyy_xzz, g_z_x_yyy_yyy, g_z_x_yyy_yyz, g_z_x_yyy_yzz, g_z_x_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyy_xxx[k] = -g_z_x_yy_xxx[k] * ab_y + g_z_x_yy_xxxy[k];

                g_z_x_yyy_xxy[k] = -g_z_x_yy_xxy[k] * ab_y + g_z_x_yy_xxyy[k];

                g_z_x_yyy_xxz[k] = -g_z_x_yy_xxz[k] * ab_y + g_z_x_yy_xxyz[k];

                g_z_x_yyy_xyy[k] = -g_z_x_yy_xyy[k] * ab_y + g_z_x_yy_xyyy[k];

                g_z_x_yyy_xyz[k] = -g_z_x_yy_xyz[k] * ab_y + g_z_x_yy_xyyz[k];

                g_z_x_yyy_xzz[k] = -g_z_x_yy_xzz[k] * ab_y + g_z_x_yy_xyzz[k];

                g_z_x_yyy_yyy[k] = -g_z_x_yy_yyy[k] * ab_y + g_z_x_yy_yyyy[k];

                g_z_x_yyy_yyz[k] = -g_z_x_yy_yyz[k] * ab_y + g_z_x_yy_yyyz[k];

                g_z_x_yyy_yzz[k] = -g_z_x_yy_yzz[k] * ab_y + g_z_x_yy_yyzz[k];

                g_z_x_yyy_zzz[k] = -g_z_x_yy_zzz[k] * ab_y + g_z_x_yy_yzzz[k];
            }

            /// Set up 670-680 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyz_xxx = cbuffer.data(ff_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_x_yyz_xxy = cbuffer.data(ff_geom_11_off + 671 * ccomps * dcomps);

            auto g_z_x_yyz_xxz = cbuffer.data(ff_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_x_yyz_xyy = cbuffer.data(ff_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_x_yyz_xyz = cbuffer.data(ff_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_x_yyz_xzz = cbuffer.data(ff_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_x_yyz_yyy = cbuffer.data(ff_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_x_yyz_yyz = cbuffer.data(ff_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_x_yyz_yzz = cbuffer.data(ff_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_x_yyz_zzz = cbuffer.data(ff_geom_11_off + 679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyz_xxx, g_z_x_yyz_xxy, g_z_x_yyz_xxz, g_z_x_yyz_xyy, g_z_x_yyz_xyz, g_z_x_yyz_xzz, g_z_x_yyz_yyy, g_z_x_yyz_yyz, g_z_x_yyz_yzz, g_z_x_yyz_zzz, g_z_x_yz_xxx, g_z_x_yz_xxxy, g_z_x_yz_xxy, g_z_x_yz_xxyy, g_z_x_yz_xxyz, g_z_x_yz_xxz, g_z_x_yz_xyy, g_z_x_yz_xyyy, g_z_x_yz_xyyz, g_z_x_yz_xyz, g_z_x_yz_xyzz, g_z_x_yz_xzz, g_z_x_yz_yyy, g_z_x_yz_yyyy, g_z_x_yz_yyyz, g_z_x_yz_yyz, g_z_x_yz_yyzz, g_z_x_yz_yzz, g_z_x_yz_yzzz, g_z_x_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyz_xxx[k] = -g_z_x_yz_xxx[k] * ab_y + g_z_x_yz_xxxy[k];

                g_z_x_yyz_xxy[k] = -g_z_x_yz_xxy[k] * ab_y + g_z_x_yz_xxyy[k];

                g_z_x_yyz_xxz[k] = -g_z_x_yz_xxz[k] * ab_y + g_z_x_yz_xxyz[k];

                g_z_x_yyz_xyy[k] = -g_z_x_yz_xyy[k] * ab_y + g_z_x_yz_xyyy[k];

                g_z_x_yyz_xyz[k] = -g_z_x_yz_xyz[k] * ab_y + g_z_x_yz_xyyz[k];

                g_z_x_yyz_xzz[k] = -g_z_x_yz_xzz[k] * ab_y + g_z_x_yz_xyzz[k];

                g_z_x_yyz_yyy[k] = -g_z_x_yz_yyy[k] * ab_y + g_z_x_yz_yyyy[k];

                g_z_x_yyz_yyz[k] = -g_z_x_yz_yyz[k] * ab_y + g_z_x_yz_yyyz[k];

                g_z_x_yyz_yzz[k] = -g_z_x_yz_yzz[k] * ab_y + g_z_x_yz_yyzz[k];

                g_z_x_yyz_zzz[k] = -g_z_x_yz_zzz[k] * ab_y + g_z_x_yz_yzzz[k];
            }

            /// Set up 680-690 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzz_xxx = cbuffer.data(ff_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_x_yzz_xxy = cbuffer.data(ff_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_x_yzz_xxz = cbuffer.data(ff_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_x_yzz_xyy = cbuffer.data(ff_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_x_yzz_xyz = cbuffer.data(ff_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_x_yzz_xzz = cbuffer.data(ff_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_x_yzz_yyy = cbuffer.data(ff_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_x_yzz_yyz = cbuffer.data(ff_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_x_yzz_yzz = cbuffer.data(ff_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_x_yzz_zzz = cbuffer.data(ff_geom_11_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzz_xxx, g_z_x_yzz_xxy, g_z_x_yzz_xxz, g_z_x_yzz_xyy, g_z_x_yzz_xyz, g_z_x_yzz_xzz, g_z_x_yzz_yyy, g_z_x_yzz_yyz, g_z_x_yzz_yzz, g_z_x_yzz_zzz, g_z_x_zz_xxx, g_z_x_zz_xxxy, g_z_x_zz_xxy, g_z_x_zz_xxyy, g_z_x_zz_xxyz, g_z_x_zz_xxz, g_z_x_zz_xyy, g_z_x_zz_xyyy, g_z_x_zz_xyyz, g_z_x_zz_xyz, g_z_x_zz_xyzz, g_z_x_zz_xzz, g_z_x_zz_yyy, g_z_x_zz_yyyy, g_z_x_zz_yyyz, g_z_x_zz_yyz, g_z_x_zz_yyzz, g_z_x_zz_yzz, g_z_x_zz_yzzz, g_z_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzz_xxx[k] = -g_z_x_zz_xxx[k] * ab_y + g_z_x_zz_xxxy[k];

                g_z_x_yzz_xxy[k] = -g_z_x_zz_xxy[k] * ab_y + g_z_x_zz_xxyy[k];

                g_z_x_yzz_xxz[k] = -g_z_x_zz_xxz[k] * ab_y + g_z_x_zz_xxyz[k];

                g_z_x_yzz_xyy[k] = -g_z_x_zz_xyy[k] * ab_y + g_z_x_zz_xyyy[k];

                g_z_x_yzz_xyz[k] = -g_z_x_zz_xyz[k] * ab_y + g_z_x_zz_xyyz[k];

                g_z_x_yzz_xzz[k] = -g_z_x_zz_xzz[k] * ab_y + g_z_x_zz_xyzz[k];

                g_z_x_yzz_yyy[k] = -g_z_x_zz_yyy[k] * ab_y + g_z_x_zz_yyyy[k];

                g_z_x_yzz_yyz[k] = -g_z_x_zz_yyz[k] * ab_y + g_z_x_zz_yyyz[k];

                g_z_x_yzz_yzz[k] = -g_z_x_zz_yzz[k] * ab_y + g_z_x_zz_yyzz[k];

                g_z_x_yzz_zzz[k] = -g_z_x_zz_zzz[k] * ab_y + g_z_x_zz_yzzz[k];
            }

            /// Set up 690-700 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzz_xxx = cbuffer.data(ff_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_x_zzz_xxy = cbuffer.data(ff_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_x_zzz_xxz = cbuffer.data(ff_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_x_zzz_xyy = cbuffer.data(ff_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_x_zzz_xyz = cbuffer.data(ff_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_x_zzz_xzz = cbuffer.data(ff_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_x_zzz_yyy = cbuffer.data(ff_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_x_zzz_yyz = cbuffer.data(ff_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_x_zzz_yzz = cbuffer.data(ff_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_x_zzz_zzz = cbuffer.data(ff_geom_11_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_xxx, g_0_x_zz_xxy, g_0_x_zz_xxz, g_0_x_zz_xyy, g_0_x_zz_xyz, g_0_x_zz_xzz, g_0_x_zz_yyy, g_0_x_zz_yyz, g_0_x_zz_yzz, g_0_x_zz_zzz, g_z_x_zz_xxx, g_z_x_zz_xxxz, g_z_x_zz_xxy, g_z_x_zz_xxyz, g_z_x_zz_xxz, g_z_x_zz_xxzz, g_z_x_zz_xyy, g_z_x_zz_xyyz, g_z_x_zz_xyz, g_z_x_zz_xyzz, g_z_x_zz_xzz, g_z_x_zz_xzzz, g_z_x_zz_yyy, g_z_x_zz_yyyz, g_z_x_zz_yyz, g_z_x_zz_yyzz, g_z_x_zz_yzz, g_z_x_zz_yzzz, g_z_x_zz_zzz, g_z_x_zz_zzzz, g_z_x_zzz_xxx, g_z_x_zzz_xxy, g_z_x_zzz_xxz, g_z_x_zzz_xyy, g_z_x_zzz_xyz, g_z_x_zzz_xzz, g_z_x_zzz_yyy, g_z_x_zzz_yyz, g_z_x_zzz_yzz, g_z_x_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzz_xxx[k] = -g_0_x_zz_xxx[k] - g_z_x_zz_xxx[k] * ab_z + g_z_x_zz_xxxz[k];

                g_z_x_zzz_xxy[k] = -g_0_x_zz_xxy[k] - g_z_x_zz_xxy[k] * ab_z + g_z_x_zz_xxyz[k];

                g_z_x_zzz_xxz[k] = -g_0_x_zz_xxz[k] - g_z_x_zz_xxz[k] * ab_z + g_z_x_zz_xxzz[k];

                g_z_x_zzz_xyy[k] = -g_0_x_zz_xyy[k] - g_z_x_zz_xyy[k] * ab_z + g_z_x_zz_xyyz[k];

                g_z_x_zzz_xyz[k] = -g_0_x_zz_xyz[k] - g_z_x_zz_xyz[k] * ab_z + g_z_x_zz_xyzz[k];

                g_z_x_zzz_xzz[k] = -g_0_x_zz_xzz[k] - g_z_x_zz_xzz[k] * ab_z + g_z_x_zz_xzzz[k];

                g_z_x_zzz_yyy[k] = -g_0_x_zz_yyy[k] - g_z_x_zz_yyy[k] * ab_z + g_z_x_zz_yyyz[k];

                g_z_x_zzz_yyz[k] = -g_0_x_zz_yyz[k] - g_z_x_zz_yyz[k] * ab_z + g_z_x_zz_yyzz[k];

                g_z_x_zzz_yzz[k] = -g_0_x_zz_yzz[k] - g_z_x_zz_yzz[k] * ab_z + g_z_x_zz_yzzz[k];

                g_z_x_zzz_zzz[k] = -g_0_x_zz_zzz[k] - g_z_x_zz_zzz[k] * ab_z + g_z_x_zz_zzzz[k];
            }

            /// Set up 700-710 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxx_xxx = cbuffer.data(ff_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_y_xxx_xxy = cbuffer.data(ff_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_y_xxx_xxz = cbuffer.data(ff_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_y_xxx_xyy = cbuffer.data(ff_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_y_xxx_xyz = cbuffer.data(ff_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_y_xxx_xzz = cbuffer.data(ff_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_y_xxx_yyy = cbuffer.data(ff_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_y_xxx_yyz = cbuffer.data(ff_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_y_xxx_yzz = cbuffer.data(ff_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_y_xxx_zzz = cbuffer.data(ff_geom_11_off + 709 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xx_xxx, g_z_y_xx_xxxx, g_z_y_xx_xxxy, g_z_y_xx_xxxz, g_z_y_xx_xxy, g_z_y_xx_xxyy, g_z_y_xx_xxyz, g_z_y_xx_xxz, g_z_y_xx_xxzz, g_z_y_xx_xyy, g_z_y_xx_xyyy, g_z_y_xx_xyyz, g_z_y_xx_xyz, g_z_y_xx_xyzz, g_z_y_xx_xzz, g_z_y_xx_xzzz, g_z_y_xx_yyy, g_z_y_xx_yyz, g_z_y_xx_yzz, g_z_y_xx_zzz, g_z_y_xxx_xxx, g_z_y_xxx_xxy, g_z_y_xxx_xxz, g_z_y_xxx_xyy, g_z_y_xxx_xyz, g_z_y_xxx_xzz, g_z_y_xxx_yyy, g_z_y_xxx_yyz, g_z_y_xxx_yzz, g_z_y_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxx_xxx[k] = -g_z_y_xx_xxx[k] * ab_x + g_z_y_xx_xxxx[k];

                g_z_y_xxx_xxy[k] = -g_z_y_xx_xxy[k] * ab_x + g_z_y_xx_xxxy[k];

                g_z_y_xxx_xxz[k] = -g_z_y_xx_xxz[k] * ab_x + g_z_y_xx_xxxz[k];

                g_z_y_xxx_xyy[k] = -g_z_y_xx_xyy[k] * ab_x + g_z_y_xx_xxyy[k];

                g_z_y_xxx_xyz[k] = -g_z_y_xx_xyz[k] * ab_x + g_z_y_xx_xxyz[k];

                g_z_y_xxx_xzz[k] = -g_z_y_xx_xzz[k] * ab_x + g_z_y_xx_xxzz[k];

                g_z_y_xxx_yyy[k] = -g_z_y_xx_yyy[k] * ab_x + g_z_y_xx_xyyy[k];

                g_z_y_xxx_yyz[k] = -g_z_y_xx_yyz[k] * ab_x + g_z_y_xx_xyyz[k];

                g_z_y_xxx_yzz[k] = -g_z_y_xx_yzz[k] * ab_x + g_z_y_xx_xyzz[k];

                g_z_y_xxx_zzz[k] = -g_z_y_xx_zzz[k] * ab_x + g_z_y_xx_xzzz[k];
            }

            /// Set up 710-720 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxy_xxx = cbuffer.data(ff_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_y_xxy_xxy = cbuffer.data(ff_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_y_xxy_xxz = cbuffer.data(ff_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_y_xxy_xyy = cbuffer.data(ff_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_y_xxy_xyz = cbuffer.data(ff_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_y_xxy_xzz = cbuffer.data(ff_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_y_xxy_yyy = cbuffer.data(ff_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_y_xxy_yyz = cbuffer.data(ff_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_y_xxy_yzz = cbuffer.data(ff_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_y_xxy_zzz = cbuffer.data(ff_geom_11_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxy_xxx, g_z_y_xxy_xxy, g_z_y_xxy_xxz, g_z_y_xxy_xyy, g_z_y_xxy_xyz, g_z_y_xxy_xzz, g_z_y_xxy_yyy, g_z_y_xxy_yyz, g_z_y_xxy_yzz, g_z_y_xxy_zzz, g_z_y_xy_xxx, g_z_y_xy_xxxx, g_z_y_xy_xxxy, g_z_y_xy_xxxz, g_z_y_xy_xxy, g_z_y_xy_xxyy, g_z_y_xy_xxyz, g_z_y_xy_xxz, g_z_y_xy_xxzz, g_z_y_xy_xyy, g_z_y_xy_xyyy, g_z_y_xy_xyyz, g_z_y_xy_xyz, g_z_y_xy_xyzz, g_z_y_xy_xzz, g_z_y_xy_xzzz, g_z_y_xy_yyy, g_z_y_xy_yyz, g_z_y_xy_yzz, g_z_y_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxy_xxx[k] = -g_z_y_xy_xxx[k] * ab_x + g_z_y_xy_xxxx[k];

                g_z_y_xxy_xxy[k] = -g_z_y_xy_xxy[k] * ab_x + g_z_y_xy_xxxy[k];

                g_z_y_xxy_xxz[k] = -g_z_y_xy_xxz[k] * ab_x + g_z_y_xy_xxxz[k];

                g_z_y_xxy_xyy[k] = -g_z_y_xy_xyy[k] * ab_x + g_z_y_xy_xxyy[k];

                g_z_y_xxy_xyz[k] = -g_z_y_xy_xyz[k] * ab_x + g_z_y_xy_xxyz[k];

                g_z_y_xxy_xzz[k] = -g_z_y_xy_xzz[k] * ab_x + g_z_y_xy_xxzz[k];

                g_z_y_xxy_yyy[k] = -g_z_y_xy_yyy[k] * ab_x + g_z_y_xy_xyyy[k];

                g_z_y_xxy_yyz[k] = -g_z_y_xy_yyz[k] * ab_x + g_z_y_xy_xyyz[k];

                g_z_y_xxy_yzz[k] = -g_z_y_xy_yzz[k] * ab_x + g_z_y_xy_xyzz[k];

                g_z_y_xxy_zzz[k] = -g_z_y_xy_zzz[k] * ab_x + g_z_y_xy_xzzz[k];
            }

            /// Set up 720-730 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxz_xxx = cbuffer.data(ff_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_y_xxz_xxy = cbuffer.data(ff_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_y_xxz_xxz = cbuffer.data(ff_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_y_xxz_xyy = cbuffer.data(ff_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_y_xxz_xyz = cbuffer.data(ff_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_y_xxz_xzz = cbuffer.data(ff_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_y_xxz_yyy = cbuffer.data(ff_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_y_xxz_yyz = cbuffer.data(ff_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_y_xxz_yzz = cbuffer.data(ff_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_y_xxz_zzz = cbuffer.data(ff_geom_11_off + 729 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxz_xxx, g_z_y_xxz_xxy, g_z_y_xxz_xxz, g_z_y_xxz_xyy, g_z_y_xxz_xyz, g_z_y_xxz_xzz, g_z_y_xxz_yyy, g_z_y_xxz_yyz, g_z_y_xxz_yzz, g_z_y_xxz_zzz, g_z_y_xz_xxx, g_z_y_xz_xxxx, g_z_y_xz_xxxy, g_z_y_xz_xxxz, g_z_y_xz_xxy, g_z_y_xz_xxyy, g_z_y_xz_xxyz, g_z_y_xz_xxz, g_z_y_xz_xxzz, g_z_y_xz_xyy, g_z_y_xz_xyyy, g_z_y_xz_xyyz, g_z_y_xz_xyz, g_z_y_xz_xyzz, g_z_y_xz_xzz, g_z_y_xz_xzzz, g_z_y_xz_yyy, g_z_y_xz_yyz, g_z_y_xz_yzz, g_z_y_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxz_xxx[k] = -g_z_y_xz_xxx[k] * ab_x + g_z_y_xz_xxxx[k];

                g_z_y_xxz_xxy[k] = -g_z_y_xz_xxy[k] * ab_x + g_z_y_xz_xxxy[k];

                g_z_y_xxz_xxz[k] = -g_z_y_xz_xxz[k] * ab_x + g_z_y_xz_xxxz[k];

                g_z_y_xxz_xyy[k] = -g_z_y_xz_xyy[k] * ab_x + g_z_y_xz_xxyy[k];

                g_z_y_xxz_xyz[k] = -g_z_y_xz_xyz[k] * ab_x + g_z_y_xz_xxyz[k];

                g_z_y_xxz_xzz[k] = -g_z_y_xz_xzz[k] * ab_x + g_z_y_xz_xxzz[k];

                g_z_y_xxz_yyy[k] = -g_z_y_xz_yyy[k] * ab_x + g_z_y_xz_xyyy[k];

                g_z_y_xxz_yyz[k] = -g_z_y_xz_yyz[k] * ab_x + g_z_y_xz_xyyz[k];

                g_z_y_xxz_yzz[k] = -g_z_y_xz_yzz[k] * ab_x + g_z_y_xz_xyzz[k];

                g_z_y_xxz_zzz[k] = -g_z_y_xz_zzz[k] * ab_x + g_z_y_xz_xzzz[k];
            }

            /// Set up 730-740 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyy_xxx = cbuffer.data(ff_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_y_xyy_xxy = cbuffer.data(ff_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_y_xyy_xxz = cbuffer.data(ff_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_y_xyy_xyy = cbuffer.data(ff_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_y_xyy_xyz = cbuffer.data(ff_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_y_xyy_xzz = cbuffer.data(ff_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_y_xyy_yyy = cbuffer.data(ff_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_y_xyy_yyz = cbuffer.data(ff_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_y_xyy_yzz = cbuffer.data(ff_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_y_xyy_zzz = cbuffer.data(ff_geom_11_off + 739 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyy_xxx, g_z_y_xyy_xxy, g_z_y_xyy_xxz, g_z_y_xyy_xyy, g_z_y_xyy_xyz, g_z_y_xyy_xzz, g_z_y_xyy_yyy, g_z_y_xyy_yyz, g_z_y_xyy_yzz, g_z_y_xyy_zzz, g_z_y_yy_xxx, g_z_y_yy_xxxx, g_z_y_yy_xxxy, g_z_y_yy_xxxz, g_z_y_yy_xxy, g_z_y_yy_xxyy, g_z_y_yy_xxyz, g_z_y_yy_xxz, g_z_y_yy_xxzz, g_z_y_yy_xyy, g_z_y_yy_xyyy, g_z_y_yy_xyyz, g_z_y_yy_xyz, g_z_y_yy_xyzz, g_z_y_yy_xzz, g_z_y_yy_xzzz, g_z_y_yy_yyy, g_z_y_yy_yyz, g_z_y_yy_yzz, g_z_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyy_xxx[k] = -g_z_y_yy_xxx[k] * ab_x + g_z_y_yy_xxxx[k];

                g_z_y_xyy_xxy[k] = -g_z_y_yy_xxy[k] * ab_x + g_z_y_yy_xxxy[k];

                g_z_y_xyy_xxz[k] = -g_z_y_yy_xxz[k] * ab_x + g_z_y_yy_xxxz[k];

                g_z_y_xyy_xyy[k] = -g_z_y_yy_xyy[k] * ab_x + g_z_y_yy_xxyy[k];

                g_z_y_xyy_xyz[k] = -g_z_y_yy_xyz[k] * ab_x + g_z_y_yy_xxyz[k];

                g_z_y_xyy_xzz[k] = -g_z_y_yy_xzz[k] * ab_x + g_z_y_yy_xxzz[k];

                g_z_y_xyy_yyy[k] = -g_z_y_yy_yyy[k] * ab_x + g_z_y_yy_xyyy[k];

                g_z_y_xyy_yyz[k] = -g_z_y_yy_yyz[k] * ab_x + g_z_y_yy_xyyz[k];

                g_z_y_xyy_yzz[k] = -g_z_y_yy_yzz[k] * ab_x + g_z_y_yy_xyzz[k];

                g_z_y_xyy_zzz[k] = -g_z_y_yy_zzz[k] * ab_x + g_z_y_yy_xzzz[k];
            }

            /// Set up 740-750 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyz_xxx = cbuffer.data(ff_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_y_xyz_xxy = cbuffer.data(ff_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_y_xyz_xxz = cbuffer.data(ff_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_y_xyz_xyy = cbuffer.data(ff_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_y_xyz_xyz = cbuffer.data(ff_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_y_xyz_xzz = cbuffer.data(ff_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_y_xyz_yyy = cbuffer.data(ff_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_y_xyz_yyz = cbuffer.data(ff_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_y_xyz_yzz = cbuffer.data(ff_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_y_xyz_zzz = cbuffer.data(ff_geom_11_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyz_xxx, g_z_y_xyz_xxy, g_z_y_xyz_xxz, g_z_y_xyz_xyy, g_z_y_xyz_xyz, g_z_y_xyz_xzz, g_z_y_xyz_yyy, g_z_y_xyz_yyz, g_z_y_xyz_yzz, g_z_y_xyz_zzz, g_z_y_yz_xxx, g_z_y_yz_xxxx, g_z_y_yz_xxxy, g_z_y_yz_xxxz, g_z_y_yz_xxy, g_z_y_yz_xxyy, g_z_y_yz_xxyz, g_z_y_yz_xxz, g_z_y_yz_xxzz, g_z_y_yz_xyy, g_z_y_yz_xyyy, g_z_y_yz_xyyz, g_z_y_yz_xyz, g_z_y_yz_xyzz, g_z_y_yz_xzz, g_z_y_yz_xzzz, g_z_y_yz_yyy, g_z_y_yz_yyz, g_z_y_yz_yzz, g_z_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyz_xxx[k] = -g_z_y_yz_xxx[k] * ab_x + g_z_y_yz_xxxx[k];

                g_z_y_xyz_xxy[k] = -g_z_y_yz_xxy[k] * ab_x + g_z_y_yz_xxxy[k];

                g_z_y_xyz_xxz[k] = -g_z_y_yz_xxz[k] * ab_x + g_z_y_yz_xxxz[k];

                g_z_y_xyz_xyy[k] = -g_z_y_yz_xyy[k] * ab_x + g_z_y_yz_xxyy[k];

                g_z_y_xyz_xyz[k] = -g_z_y_yz_xyz[k] * ab_x + g_z_y_yz_xxyz[k];

                g_z_y_xyz_xzz[k] = -g_z_y_yz_xzz[k] * ab_x + g_z_y_yz_xxzz[k];

                g_z_y_xyz_yyy[k] = -g_z_y_yz_yyy[k] * ab_x + g_z_y_yz_xyyy[k];

                g_z_y_xyz_yyz[k] = -g_z_y_yz_yyz[k] * ab_x + g_z_y_yz_xyyz[k];

                g_z_y_xyz_yzz[k] = -g_z_y_yz_yzz[k] * ab_x + g_z_y_yz_xyzz[k];

                g_z_y_xyz_zzz[k] = -g_z_y_yz_zzz[k] * ab_x + g_z_y_yz_xzzz[k];
            }

            /// Set up 750-760 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzz_xxx = cbuffer.data(ff_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_y_xzz_xxy = cbuffer.data(ff_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_y_xzz_xxz = cbuffer.data(ff_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_y_xzz_xyy = cbuffer.data(ff_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_y_xzz_xyz = cbuffer.data(ff_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_y_xzz_xzz = cbuffer.data(ff_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_y_xzz_yyy = cbuffer.data(ff_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_y_xzz_yyz = cbuffer.data(ff_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_y_xzz_yzz = cbuffer.data(ff_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_y_xzz_zzz = cbuffer.data(ff_geom_11_off + 759 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzz_xxx, g_z_y_xzz_xxy, g_z_y_xzz_xxz, g_z_y_xzz_xyy, g_z_y_xzz_xyz, g_z_y_xzz_xzz, g_z_y_xzz_yyy, g_z_y_xzz_yyz, g_z_y_xzz_yzz, g_z_y_xzz_zzz, g_z_y_zz_xxx, g_z_y_zz_xxxx, g_z_y_zz_xxxy, g_z_y_zz_xxxz, g_z_y_zz_xxy, g_z_y_zz_xxyy, g_z_y_zz_xxyz, g_z_y_zz_xxz, g_z_y_zz_xxzz, g_z_y_zz_xyy, g_z_y_zz_xyyy, g_z_y_zz_xyyz, g_z_y_zz_xyz, g_z_y_zz_xyzz, g_z_y_zz_xzz, g_z_y_zz_xzzz, g_z_y_zz_yyy, g_z_y_zz_yyz, g_z_y_zz_yzz, g_z_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzz_xxx[k] = -g_z_y_zz_xxx[k] * ab_x + g_z_y_zz_xxxx[k];

                g_z_y_xzz_xxy[k] = -g_z_y_zz_xxy[k] * ab_x + g_z_y_zz_xxxy[k];

                g_z_y_xzz_xxz[k] = -g_z_y_zz_xxz[k] * ab_x + g_z_y_zz_xxxz[k];

                g_z_y_xzz_xyy[k] = -g_z_y_zz_xyy[k] * ab_x + g_z_y_zz_xxyy[k];

                g_z_y_xzz_xyz[k] = -g_z_y_zz_xyz[k] * ab_x + g_z_y_zz_xxyz[k];

                g_z_y_xzz_xzz[k] = -g_z_y_zz_xzz[k] * ab_x + g_z_y_zz_xxzz[k];

                g_z_y_xzz_yyy[k] = -g_z_y_zz_yyy[k] * ab_x + g_z_y_zz_xyyy[k];

                g_z_y_xzz_yyz[k] = -g_z_y_zz_yyz[k] * ab_x + g_z_y_zz_xyyz[k];

                g_z_y_xzz_yzz[k] = -g_z_y_zz_yzz[k] * ab_x + g_z_y_zz_xyzz[k];

                g_z_y_xzz_zzz[k] = -g_z_y_zz_zzz[k] * ab_x + g_z_y_zz_xzzz[k];
            }

            /// Set up 760-770 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyy_xxx = cbuffer.data(ff_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_y_yyy_xxy = cbuffer.data(ff_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_y_yyy_xxz = cbuffer.data(ff_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_y_yyy_xyy = cbuffer.data(ff_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_y_yyy_xyz = cbuffer.data(ff_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_y_yyy_xzz = cbuffer.data(ff_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_y_yyy_yyy = cbuffer.data(ff_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_y_yyy_yyz = cbuffer.data(ff_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_y_yyy_yzz = cbuffer.data(ff_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_y_yyy_zzz = cbuffer.data(ff_geom_11_off + 769 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_xxx, g_z_0_yy_xxy, g_z_0_yy_xxz, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xzz, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yzz, g_z_0_yy_zzz, g_z_y_yy_xxx, g_z_y_yy_xxxy, g_z_y_yy_xxy, g_z_y_yy_xxyy, g_z_y_yy_xxyz, g_z_y_yy_xxz, g_z_y_yy_xyy, g_z_y_yy_xyyy, g_z_y_yy_xyyz, g_z_y_yy_xyz, g_z_y_yy_xyzz, g_z_y_yy_xzz, g_z_y_yy_yyy, g_z_y_yy_yyyy, g_z_y_yy_yyyz, g_z_y_yy_yyz, g_z_y_yy_yyzz, g_z_y_yy_yzz, g_z_y_yy_yzzz, g_z_y_yy_zzz, g_z_y_yyy_xxx, g_z_y_yyy_xxy, g_z_y_yyy_xxz, g_z_y_yyy_xyy, g_z_y_yyy_xyz, g_z_y_yyy_xzz, g_z_y_yyy_yyy, g_z_y_yyy_yyz, g_z_y_yyy_yzz, g_z_y_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyy_xxx[k] = g_z_0_yy_xxx[k] - g_z_y_yy_xxx[k] * ab_y + g_z_y_yy_xxxy[k];

                g_z_y_yyy_xxy[k] = g_z_0_yy_xxy[k] - g_z_y_yy_xxy[k] * ab_y + g_z_y_yy_xxyy[k];

                g_z_y_yyy_xxz[k] = g_z_0_yy_xxz[k] - g_z_y_yy_xxz[k] * ab_y + g_z_y_yy_xxyz[k];

                g_z_y_yyy_xyy[k] = g_z_0_yy_xyy[k] - g_z_y_yy_xyy[k] * ab_y + g_z_y_yy_xyyy[k];

                g_z_y_yyy_xyz[k] = g_z_0_yy_xyz[k] - g_z_y_yy_xyz[k] * ab_y + g_z_y_yy_xyyz[k];

                g_z_y_yyy_xzz[k] = g_z_0_yy_xzz[k] - g_z_y_yy_xzz[k] * ab_y + g_z_y_yy_xyzz[k];

                g_z_y_yyy_yyy[k] = g_z_0_yy_yyy[k] - g_z_y_yy_yyy[k] * ab_y + g_z_y_yy_yyyy[k];

                g_z_y_yyy_yyz[k] = g_z_0_yy_yyz[k] - g_z_y_yy_yyz[k] * ab_y + g_z_y_yy_yyyz[k];

                g_z_y_yyy_yzz[k] = g_z_0_yy_yzz[k] - g_z_y_yy_yzz[k] * ab_y + g_z_y_yy_yyzz[k];

                g_z_y_yyy_zzz[k] = g_z_0_yy_zzz[k] - g_z_y_yy_zzz[k] * ab_y + g_z_y_yy_yzzz[k];
            }

            /// Set up 770-780 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyz_xxx = cbuffer.data(ff_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_y_yyz_xxy = cbuffer.data(ff_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_y_yyz_xxz = cbuffer.data(ff_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_y_yyz_xyy = cbuffer.data(ff_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_y_yyz_xyz = cbuffer.data(ff_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_y_yyz_xzz = cbuffer.data(ff_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_y_yyz_yyy = cbuffer.data(ff_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_y_yyz_yyz = cbuffer.data(ff_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_y_yyz_yzz = cbuffer.data(ff_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_y_yyz_zzz = cbuffer.data(ff_geom_11_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xxx, g_z_0_yz_xxy, g_z_0_yz_xxz, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xzz, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yzz, g_z_0_yz_zzz, g_z_y_yyz_xxx, g_z_y_yyz_xxy, g_z_y_yyz_xxz, g_z_y_yyz_xyy, g_z_y_yyz_xyz, g_z_y_yyz_xzz, g_z_y_yyz_yyy, g_z_y_yyz_yyz, g_z_y_yyz_yzz, g_z_y_yyz_zzz, g_z_y_yz_xxx, g_z_y_yz_xxxy, g_z_y_yz_xxy, g_z_y_yz_xxyy, g_z_y_yz_xxyz, g_z_y_yz_xxz, g_z_y_yz_xyy, g_z_y_yz_xyyy, g_z_y_yz_xyyz, g_z_y_yz_xyz, g_z_y_yz_xyzz, g_z_y_yz_xzz, g_z_y_yz_yyy, g_z_y_yz_yyyy, g_z_y_yz_yyyz, g_z_y_yz_yyz, g_z_y_yz_yyzz, g_z_y_yz_yzz, g_z_y_yz_yzzz, g_z_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyz_xxx[k] = g_z_0_yz_xxx[k] - g_z_y_yz_xxx[k] * ab_y + g_z_y_yz_xxxy[k];

                g_z_y_yyz_xxy[k] = g_z_0_yz_xxy[k] - g_z_y_yz_xxy[k] * ab_y + g_z_y_yz_xxyy[k];

                g_z_y_yyz_xxz[k] = g_z_0_yz_xxz[k] - g_z_y_yz_xxz[k] * ab_y + g_z_y_yz_xxyz[k];

                g_z_y_yyz_xyy[k] = g_z_0_yz_xyy[k] - g_z_y_yz_xyy[k] * ab_y + g_z_y_yz_xyyy[k];

                g_z_y_yyz_xyz[k] = g_z_0_yz_xyz[k] - g_z_y_yz_xyz[k] * ab_y + g_z_y_yz_xyyz[k];

                g_z_y_yyz_xzz[k] = g_z_0_yz_xzz[k] - g_z_y_yz_xzz[k] * ab_y + g_z_y_yz_xyzz[k];

                g_z_y_yyz_yyy[k] = g_z_0_yz_yyy[k] - g_z_y_yz_yyy[k] * ab_y + g_z_y_yz_yyyy[k];

                g_z_y_yyz_yyz[k] = g_z_0_yz_yyz[k] - g_z_y_yz_yyz[k] * ab_y + g_z_y_yz_yyyz[k];

                g_z_y_yyz_yzz[k] = g_z_0_yz_yzz[k] - g_z_y_yz_yzz[k] * ab_y + g_z_y_yz_yyzz[k];

                g_z_y_yyz_zzz[k] = g_z_0_yz_zzz[k] - g_z_y_yz_zzz[k] * ab_y + g_z_y_yz_yzzz[k];
            }

            /// Set up 780-790 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzz_xxx = cbuffer.data(ff_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_y_yzz_xxy = cbuffer.data(ff_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_y_yzz_xxz = cbuffer.data(ff_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_y_yzz_xyy = cbuffer.data(ff_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_y_yzz_xyz = cbuffer.data(ff_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_y_yzz_xzz = cbuffer.data(ff_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_y_yzz_yyy = cbuffer.data(ff_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_y_yzz_yyz = cbuffer.data(ff_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_y_yzz_yzz = cbuffer.data(ff_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_y_yzz_zzz = cbuffer.data(ff_geom_11_off + 789 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz, g_z_y_yzz_xxx, g_z_y_yzz_xxy, g_z_y_yzz_xxz, g_z_y_yzz_xyy, g_z_y_yzz_xyz, g_z_y_yzz_xzz, g_z_y_yzz_yyy, g_z_y_yzz_yyz, g_z_y_yzz_yzz, g_z_y_yzz_zzz, g_z_y_zz_xxx, g_z_y_zz_xxxy, g_z_y_zz_xxy, g_z_y_zz_xxyy, g_z_y_zz_xxyz, g_z_y_zz_xxz, g_z_y_zz_xyy, g_z_y_zz_xyyy, g_z_y_zz_xyyz, g_z_y_zz_xyz, g_z_y_zz_xyzz, g_z_y_zz_xzz, g_z_y_zz_yyy, g_z_y_zz_yyyy, g_z_y_zz_yyyz, g_z_y_zz_yyz, g_z_y_zz_yyzz, g_z_y_zz_yzz, g_z_y_zz_yzzz, g_z_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzz_xxx[k] = g_z_0_zz_xxx[k] - g_z_y_zz_xxx[k] * ab_y + g_z_y_zz_xxxy[k];

                g_z_y_yzz_xxy[k] = g_z_0_zz_xxy[k] - g_z_y_zz_xxy[k] * ab_y + g_z_y_zz_xxyy[k];

                g_z_y_yzz_xxz[k] = g_z_0_zz_xxz[k] - g_z_y_zz_xxz[k] * ab_y + g_z_y_zz_xxyz[k];

                g_z_y_yzz_xyy[k] = g_z_0_zz_xyy[k] - g_z_y_zz_xyy[k] * ab_y + g_z_y_zz_xyyy[k];

                g_z_y_yzz_xyz[k] = g_z_0_zz_xyz[k] - g_z_y_zz_xyz[k] * ab_y + g_z_y_zz_xyyz[k];

                g_z_y_yzz_xzz[k] = g_z_0_zz_xzz[k] - g_z_y_zz_xzz[k] * ab_y + g_z_y_zz_xyzz[k];

                g_z_y_yzz_yyy[k] = g_z_0_zz_yyy[k] - g_z_y_zz_yyy[k] * ab_y + g_z_y_zz_yyyy[k];

                g_z_y_yzz_yyz[k] = g_z_0_zz_yyz[k] - g_z_y_zz_yyz[k] * ab_y + g_z_y_zz_yyyz[k];

                g_z_y_yzz_yzz[k] = g_z_0_zz_yzz[k] - g_z_y_zz_yzz[k] * ab_y + g_z_y_zz_yyzz[k];

                g_z_y_yzz_zzz[k] = g_z_0_zz_zzz[k] - g_z_y_zz_zzz[k] * ab_y + g_z_y_zz_yzzz[k];
            }

            /// Set up 790-800 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzz_xxx = cbuffer.data(ff_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_y_zzz_xxy = cbuffer.data(ff_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_y_zzz_xxz = cbuffer.data(ff_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_y_zzz_xyy = cbuffer.data(ff_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_y_zzz_xyz = cbuffer.data(ff_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_y_zzz_xzz = cbuffer.data(ff_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_y_zzz_yyy = cbuffer.data(ff_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_y_zzz_yyz = cbuffer.data(ff_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_y_zzz_yzz = cbuffer.data(ff_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_y_zzz_zzz = cbuffer.data(ff_geom_11_off + 799 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_xxx, g_0_y_zz_xxy, g_0_y_zz_xxz, g_0_y_zz_xyy, g_0_y_zz_xyz, g_0_y_zz_xzz, g_0_y_zz_yyy, g_0_y_zz_yyz, g_0_y_zz_yzz, g_0_y_zz_zzz, g_z_y_zz_xxx, g_z_y_zz_xxxz, g_z_y_zz_xxy, g_z_y_zz_xxyz, g_z_y_zz_xxz, g_z_y_zz_xxzz, g_z_y_zz_xyy, g_z_y_zz_xyyz, g_z_y_zz_xyz, g_z_y_zz_xyzz, g_z_y_zz_xzz, g_z_y_zz_xzzz, g_z_y_zz_yyy, g_z_y_zz_yyyz, g_z_y_zz_yyz, g_z_y_zz_yyzz, g_z_y_zz_yzz, g_z_y_zz_yzzz, g_z_y_zz_zzz, g_z_y_zz_zzzz, g_z_y_zzz_xxx, g_z_y_zzz_xxy, g_z_y_zzz_xxz, g_z_y_zzz_xyy, g_z_y_zzz_xyz, g_z_y_zzz_xzz, g_z_y_zzz_yyy, g_z_y_zzz_yyz, g_z_y_zzz_yzz, g_z_y_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzz_xxx[k] = -g_0_y_zz_xxx[k] - g_z_y_zz_xxx[k] * ab_z + g_z_y_zz_xxxz[k];

                g_z_y_zzz_xxy[k] = -g_0_y_zz_xxy[k] - g_z_y_zz_xxy[k] * ab_z + g_z_y_zz_xxyz[k];

                g_z_y_zzz_xxz[k] = -g_0_y_zz_xxz[k] - g_z_y_zz_xxz[k] * ab_z + g_z_y_zz_xxzz[k];

                g_z_y_zzz_xyy[k] = -g_0_y_zz_xyy[k] - g_z_y_zz_xyy[k] * ab_z + g_z_y_zz_xyyz[k];

                g_z_y_zzz_xyz[k] = -g_0_y_zz_xyz[k] - g_z_y_zz_xyz[k] * ab_z + g_z_y_zz_xyzz[k];

                g_z_y_zzz_xzz[k] = -g_0_y_zz_xzz[k] - g_z_y_zz_xzz[k] * ab_z + g_z_y_zz_xzzz[k];

                g_z_y_zzz_yyy[k] = -g_0_y_zz_yyy[k] - g_z_y_zz_yyy[k] * ab_z + g_z_y_zz_yyyz[k];

                g_z_y_zzz_yyz[k] = -g_0_y_zz_yyz[k] - g_z_y_zz_yyz[k] * ab_z + g_z_y_zz_yyzz[k];

                g_z_y_zzz_yzz[k] = -g_0_y_zz_yzz[k] - g_z_y_zz_yzz[k] * ab_z + g_z_y_zz_yzzz[k];

                g_z_y_zzz_zzz[k] = -g_0_y_zz_zzz[k] - g_z_y_zz_zzz[k] * ab_z + g_z_y_zz_zzzz[k];
            }

            /// Set up 800-810 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxx_xxx = cbuffer.data(ff_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_z_xxx_xxy = cbuffer.data(ff_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_z_xxx_xxz = cbuffer.data(ff_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_z_xxx_xyy = cbuffer.data(ff_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_z_xxx_xyz = cbuffer.data(ff_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_z_xxx_xzz = cbuffer.data(ff_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_z_xxx_yyy = cbuffer.data(ff_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_z_xxx_yyz = cbuffer.data(ff_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_z_xxx_yzz = cbuffer.data(ff_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_z_xxx_zzz = cbuffer.data(ff_geom_11_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xx_xxx, g_z_z_xx_xxxx, g_z_z_xx_xxxy, g_z_z_xx_xxxz, g_z_z_xx_xxy, g_z_z_xx_xxyy, g_z_z_xx_xxyz, g_z_z_xx_xxz, g_z_z_xx_xxzz, g_z_z_xx_xyy, g_z_z_xx_xyyy, g_z_z_xx_xyyz, g_z_z_xx_xyz, g_z_z_xx_xyzz, g_z_z_xx_xzz, g_z_z_xx_xzzz, g_z_z_xx_yyy, g_z_z_xx_yyz, g_z_z_xx_yzz, g_z_z_xx_zzz, g_z_z_xxx_xxx, g_z_z_xxx_xxy, g_z_z_xxx_xxz, g_z_z_xxx_xyy, g_z_z_xxx_xyz, g_z_z_xxx_xzz, g_z_z_xxx_yyy, g_z_z_xxx_yyz, g_z_z_xxx_yzz, g_z_z_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxx_xxx[k] = -g_z_z_xx_xxx[k] * ab_x + g_z_z_xx_xxxx[k];

                g_z_z_xxx_xxy[k] = -g_z_z_xx_xxy[k] * ab_x + g_z_z_xx_xxxy[k];

                g_z_z_xxx_xxz[k] = -g_z_z_xx_xxz[k] * ab_x + g_z_z_xx_xxxz[k];

                g_z_z_xxx_xyy[k] = -g_z_z_xx_xyy[k] * ab_x + g_z_z_xx_xxyy[k];

                g_z_z_xxx_xyz[k] = -g_z_z_xx_xyz[k] * ab_x + g_z_z_xx_xxyz[k];

                g_z_z_xxx_xzz[k] = -g_z_z_xx_xzz[k] * ab_x + g_z_z_xx_xxzz[k];

                g_z_z_xxx_yyy[k] = -g_z_z_xx_yyy[k] * ab_x + g_z_z_xx_xyyy[k];

                g_z_z_xxx_yyz[k] = -g_z_z_xx_yyz[k] * ab_x + g_z_z_xx_xyyz[k];

                g_z_z_xxx_yzz[k] = -g_z_z_xx_yzz[k] * ab_x + g_z_z_xx_xyzz[k];

                g_z_z_xxx_zzz[k] = -g_z_z_xx_zzz[k] * ab_x + g_z_z_xx_xzzz[k];
            }

            /// Set up 810-820 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxy_xxx = cbuffer.data(ff_geom_11_off + 810 * ccomps * dcomps);

            auto g_z_z_xxy_xxy = cbuffer.data(ff_geom_11_off + 811 * ccomps * dcomps);

            auto g_z_z_xxy_xxz = cbuffer.data(ff_geom_11_off + 812 * ccomps * dcomps);

            auto g_z_z_xxy_xyy = cbuffer.data(ff_geom_11_off + 813 * ccomps * dcomps);

            auto g_z_z_xxy_xyz = cbuffer.data(ff_geom_11_off + 814 * ccomps * dcomps);

            auto g_z_z_xxy_xzz = cbuffer.data(ff_geom_11_off + 815 * ccomps * dcomps);

            auto g_z_z_xxy_yyy = cbuffer.data(ff_geom_11_off + 816 * ccomps * dcomps);

            auto g_z_z_xxy_yyz = cbuffer.data(ff_geom_11_off + 817 * ccomps * dcomps);

            auto g_z_z_xxy_yzz = cbuffer.data(ff_geom_11_off + 818 * ccomps * dcomps);

            auto g_z_z_xxy_zzz = cbuffer.data(ff_geom_11_off + 819 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxy_xxx, g_z_z_xxy_xxy, g_z_z_xxy_xxz, g_z_z_xxy_xyy, g_z_z_xxy_xyz, g_z_z_xxy_xzz, g_z_z_xxy_yyy, g_z_z_xxy_yyz, g_z_z_xxy_yzz, g_z_z_xxy_zzz, g_z_z_xy_xxx, g_z_z_xy_xxxx, g_z_z_xy_xxxy, g_z_z_xy_xxxz, g_z_z_xy_xxy, g_z_z_xy_xxyy, g_z_z_xy_xxyz, g_z_z_xy_xxz, g_z_z_xy_xxzz, g_z_z_xy_xyy, g_z_z_xy_xyyy, g_z_z_xy_xyyz, g_z_z_xy_xyz, g_z_z_xy_xyzz, g_z_z_xy_xzz, g_z_z_xy_xzzz, g_z_z_xy_yyy, g_z_z_xy_yyz, g_z_z_xy_yzz, g_z_z_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxy_xxx[k] = -g_z_z_xy_xxx[k] * ab_x + g_z_z_xy_xxxx[k];

                g_z_z_xxy_xxy[k] = -g_z_z_xy_xxy[k] * ab_x + g_z_z_xy_xxxy[k];

                g_z_z_xxy_xxz[k] = -g_z_z_xy_xxz[k] * ab_x + g_z_z_xy_xxxz[k];

                g_z_z_xxy_xyy[k] = -g_z_z_xy_xyy[k] * ab_x + g_z_z_xy_xxyy[k];

                g_z_z_xxy_xyz[k] = -g_z_z_xy_xyz[k] * ab_x + g_z_z_xy_xxyz[k];

                g_z_z_xxy_xzz[k] = -g_z_z_xy_xzz[k] * ab_x + g_z_z_xy_xxzz[k];

                g_z_z_xxy_yyy[k] = -g_z_z_xy_yyy[k] * ab_x + g_z_z_xy_xyyy[k];

                g_z_z_xxy_yyz[k] = -g_z_z_xy_yyz[k] * ab_x + g_z_z_xy_xyyz[k];

                g_z_z_xxy_yzz[k] = -g_z_z_xy_yzz[k] * ab_x + g_z_z_xy_xyzz[k];

                g_z_z_xxy_zzz[k] = -g_z_z_xy_zzz[k] * ab_x + g_z_z_xy_xzzz[k];
            }

            /// Set up 820-830 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxz_xxx = cbuffer.data(ff_geom_11_off + 820 * ccomps * dcomps);

            auto g_z_z_xxz_xxy = cbuffer.data(ff_geom_11_off + 821 * ccomps * dcomps);

            auto g_z_z_xxz_xxz = cbuffer.data(ff_geom_11_off + 822 * ccomps * dcomps);

            auto g_z_z_xxz_xyy = cbuffer.data(ff_geom_11_off + 823 * ccomps * dcomps);

            auto g_z_z_xxz_xyz = cbuffer.data(ff_geom_11_off + 824 * ccomps * dcomps);

            auto g_z_z_xxz_xzz = cbuffer.data(ff_geom_11_off + 825 * ccomps * dcomps);

            auto g_z_z_xxz_yyy = cbuffer.data(ff_geom_11_off + 826 * ccomps * dcomps);

            auto g_z_z_xxz_yyz = cbuffer.data(ff_geom_11_off + 827 * ccomps * dcomps);

            auto g_z_z_xxz_yzz = cbuffer.data(ff_geom_11_off + 828 * ccomps * dcomps);

            auto g_z_z_xxz_zzz = cbuffer.data(ff_geom_11_off + 829 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxz_xxx, g_z_z_xxz_xxy, g_z_z_xxz_xxz, g_z_z_xxz_xyy, g_z_z_xxz_xyz, g_z_z_xxz_xzz, g_z_z_xxz_yyy, g_z_z_xxz_yyz, g_z_z_xxz_yzz, g_z_z_xxz_zzz, g_z_z_xz_xxx, g_z_z_xz_xxxx, g_z_z_xz_xxxy, g_z_z_xz_xxxz, g_z_z_xz_xxy, g_z_z_xz_xxyy, g_z_z_xz_xxyz, g_z_z_xz_xxz, g_z_z_xz_xxzz, g_z_z_xz_xyy, g_z_z_xz_xyyy, g_z_z_xz_xyyz, g_z_z_xz_xyz, g_z_z_xz_xyzz, g_z_z_xz_xzz, g_z_z_xz_xzzz, g_z_z_xz_yyy, g_z_z_xz_yyz, g_z_z_xz_yzz, g_z_z_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxz_xxx[k] = -g_z_z_xz_xxx[k] * ab_x + g_z_z_xz_xxxx[k];

                g_z_z_xxz_xxy[k] = -g_z_z_xz_xxy[k] * ab_x + g_z_z_xz_xxxy[k];

                g_z_z_xxz_xxz[k] = -g_z_z_xz_xxz[k] * ab_x + g_z_z_xz_xxxz[k];

                g_z_z_xxz_xyy[k] = -g_z_z_xz_xyy[k] * ab_x + g_z_z_xz_xxyy[k];

                g_z_z_xxz_xyz[k] = -g_z_z_xz_xyz[k] * ab_x + g_z_z_xz_xxyz[k];

                g_z_z_xxz_xzz[k] = -g_z_z_xz_xzz[k] * ab_x + g_z_z_xz_xxzz[k];

                g_z_z_xxz_yyy[k] = -g_z_z_xz_yyy[k] * ab_x + g_z_z_xz_xyyy[k];

                g_z_z_xxz_yyz[k] = -g_z_z_xz_yyz[k] * ab_x + g_z_z_xz_xyyz[k];

                g_z_z_xxz_yzz[k] = -g_z_z_xz_yzz[k] * ab_x + g_z_z_xz_xyzz[k];

                g_z_z_xxz_zzz[k] = -g_z_z_xz_zzz[k] * ab_x + g_z_z_xz_xzzz[k];
            }

            /// Set up 830-840 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyy_xxx = cbuffer.data(ff_geom_11_off + 830 * ccomps * dcomps);

            auto g_z_z_xyy_xxy = cbuffer.data(ff_geom_11_off + 831 * ccomps * dcomps);

            auto g_z_z_xyy_xxz = cbuffer.data(ff_geom_11_off + 832 * ccomps * dcomps);

            auto g_z_z_xyy_xyy = cbuffer.data(ff_geom_11_off + 833 * ccomps * dcomps);

            auto g_z_z_xyy_xyz = cbuffer.data(ff_geom_11_off + 834 * ccomps * dcomps);

            auto g_z_z_xyy_xzz = cbuffer.data(ff_geom_11_off + 835 * ccomps * dcomps);

            auto g_z_z_xyy_yyy = cbuffer.data(ff_geom_11_off + 836 * ccomps * dcomps);

            auto g_z_z_xyy_yyz = cbuffer.data(ff_geom_11_off + 837 * ccomps * dcomps);

            auto g_z_z_xyy_yzz = cbuffer.data(ff_geom_11_off + 838 * ccomps * dcomps);

            auto g_z_z_xyy_zzz = cbuffer.data(ff_geom_11_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyy_xxx, g_z_z_xyy_xxy, g_z_z_xyy_xxz, g_z_z_xyy_xyy, g_z_z_xyy_xyz, g_z_z_xyy_xzz, g_z_z_xyy_yyy, g_z_z_xyy_yyz, g_z_z_xyy_yzz, g_z_z_xyy_zzz, g_z_z_yy_xxx, g_z_z_yy_xxxx, g_z_z_yy_xxxy, g_z_z_yy_xxxz, g_z_z_yy_xxy, g_z_z_yy_xxyy, g_z_z_yy_xxyz, g_z_z_yy_xxz, g_z_z_yy_xxzz, g_z_z_yy_xyy, g_z_z_yy_xyyy, g_z_z_yy_xyyz, g_z_z_yy_xyz, g_z_z_yy_xyzz, g_z_z_yy_xzz, g_z_z_yy_xzzz, g_z_z_yy_yyy, g_z_z_yy_yyz, g_z_z_yy_yzz, g_z_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyy_xxx[k] = -g_z_z_yy_xxx[k] * ab_x + g_z_z_yy_xxxx[k];

                g_z_z_xyy_xxy[k] = -g_z_z_yy_xxy[k] * ab_x + g_z_z_yy_xxxy[k];

                g_z_z_xyy_xxz[k] = -g_z_z_yy_xxz[k] * ab_x + g_z_z_yy_xxxz[k];

                g_z_z_xyy_xyy[k] = -g_z_z_yy_xyy[k] * ab_x + g_z_z_yy_xxyy[k];

                g_z_z_xyy_xyz[k] = -g_z_z_yy_xyz[k] * ab_x + g_z_z_yy_xxyz[k];

                g_z_z_xyy_xzz[k] = -g_z_z_yy_xzz[k] * ab_x + g_z_z_yy_xxzz[k];

                g_z_z_xyy_yyy[k] = -g_z_z_yy_yyy[k] * ab_x + g_z_z_yy_xyyy[k];

                g_z_z_xyy_yyz[k] = -g_z_z_yy_yyz[k] * ab_x + g_z_z_yy_xyyz[k];

                g_z_z_xyy_yzz[k] = -g_z_z_yy_yzz[k] * ab_x + g_z_z_yy_xyzz[k];

                g_z_z_xyy_zzz[k] = -g_z_z_yy_zzz[k] * ab_x + g_z_z_yy_xzzz[k];
            }

            /// Set up 840-850 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyz_xxx = cbuffer.data(ff_geom_11_off + 840 * ccomps * dcomps);

            auto g_z_z_xyz_xxy = cbuffer.data(ff_geom_11_off + 841 * ccomps * dcomps);

            auto g_z_z_xyz_xxz = cbuffer.data(ff_geom_11_off + 842 * ccomps * dcomps);

            auto g_z_z_xyz_xyy = cbuffer.data(ff_geom_11_off + 843 * ccomps * dcomps);

            auto g_z_z_xyz_xyz = cbuffer.data(ff_geom_11_off + 844 * ccomps * dcomps);

            auto g_z_z_xyz_xzz = cbuffer.data(ff_geom_11_off + 845 * ccomps * dcomps);

            auto g_z_z_xyz_yyy = cbuffer.data(ff_geom_11_off + 846 * ccomps * dcomps);

            auto g_z_z_xyz_yyz = cbuffer.data(ff_geom_11_off + 847 * ccomps * dcomps);

            auto g_z_z_xyz_yzz = cbuffer.data(ff_geom_11_off + 848 * ccomps * dcomps);

            auto g_z_z_xyz_zzz = cbuffer.data(ff_geom_11_off + 849 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyz_xxx, g_z_z_xyz_xxy, g_z_z_xyz_xxz, g_z_z_xyz_xyy, g_z_z_xyz_xyz, g_z_z_xyz_xzz, g_z_z_xyz_yyy, g_z_z_xyz_yyz, g_z_z_xyz_yzz, g_z_z_xyz_zzz, g_z_z_yz_xxx, g_z_z_yz_xxxx, g_z_z_yz_xxxy, g_z_z_yz_xxxz, g_z_z_yz_xxy, g_z_z_yz_xxyy, g_z_z_yz_xxyz, g_z_z_yz_xxz, g_z_z_yz_xxzz, g_z_z_yz_xyy, g_z_z_yz_xyyy, g_z_z_yz_xyyz, g_z_z_yz_xyz, g_z_z_yz_xyzz, g_z_z_yz_xzz, g_z_z_yz_xzzz, g_z_z_yz_yyy, g_z_z_yz_yyz, g_z_z_yz_yzz, g_z_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyz_xxx[k] = -g_z_z_yz_xxx[k] * ab_x + g_z_z_yz_xxxx[k];

                g_z_z_xyz_xxy[k] = -g_z_z_yz_xxy[k] * ab_x + g_z_z_yz_xxxy[k];

                g_z_z_xyz_xxz[k] = -g_z_z_yz_xxz[k] * ab_x + g_z_z_yz_xxxz[k];

                g_z_z_xyz_xyy[k] = -g_z_z_yz_xyy[k] * ab_x + g_z_z_yz_xxyy[k];

                g_z_z_xyz_xyz[k] = -g_z_z_yz_xyz[k] * ab_x + g_z_z_yz_xxyz[k];

                g_z_z_xyz_xzz[k] = -g_z_z_yz_xzz[k] * ab_x + g_z_z_yz_xxzz[k];

                g_z_z_xyz_yyy[k] = -g_z_z_yz_yyy[k] * ab_x + g_z_z_yz_xyyy[k];

                g_z_z_xyz_yyz[k] = -g_z_z_yz_yyz[k] * ab_x + g_z_z_yz_xyyz[k];

                g_z_z_xyz_yzz[k] = -g_z_z_yz_yzz[k] * ab_x + g_z_z_yz_xyzz[k];

                g_z_z_xyz_zzz[k] = -g_z_z_yz_zzz[k] * ab_x + g_z_z_yz_xzzz[k];
            }

            /// Set up 850-860 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzz_xxx = cbuffer.data(ff_geom_11_off + 850 * ccomps * dcomps);

            auto g_z_z_xzz_xxy = cbuffer.data(ff_geom_11_off + 851 * ccomps * dcomps);

            auto g_z_z_xzz_xxz = cbuffer.data(ff_geom_11_off + 852 * ccomps * dcomps);

            auto g_z_z_xzz_xyy = cbuffer.data(ff_geom_11_off + 853 * ccomps * dcomps);

            auto g_z_z_xzz_xyz = cbuffer.data(ff_geom_11_off + 854 * ccomps * dcomps);

            auto g_z_z_xzz_xzz = cbuffer.data(ff_geom_11_off + 855 * ccomps * dcomps);

            auto g_z_z_xzz_yyy = cbuffer.data(ff_geom_11_off + 856 * ccomps * dcomps);

            auto g_z_z_xzz_yyz = cbuffer.data(ff_geom_11_off + 857 * ccomps * dcomps);

            auto g_z_z_xzz_yzz = cbuffer.data(ff_geom_11_off + 858 * ccomps * dcomps);

            auto g_z_z_xzz_zzz = cbuffer.data(ff_geom_11_off + 859 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzz_xxx, g_z_z_xzz_xxy, g_z_z_xzz_xxz, g_z_z_xzz_xyy, g_z_z_xzz_xyz, g_z_z_xzz_xzz, g_z_z_xzz_yyy, g_z_z_xzz_yyz, g_z_z_xzz_yzz, g_z_z_xzz_zzz, g_z_z_zz_xxx, g_z_z_zz_xxxx, g_z_z_zz_xxxy, g_z_z_zz_xxxz, g_z_z_zz_xxy, g_z_z_zz_xxyy, g_z_z_zz_xxyz, g_z_z_zz_xxz, g_z_z_zz_xxzz, g_z_z_zz_xyy, g_z_z_zz_xyyy, g_z_z_zz_xyyz, g_z_z_zz_xyz, g_z_z_zz_xyzz, g_z_z_zz_xzz, g_z_z_zz_xzzz, g_z_z_zz_yyy, g_z_z_zz_yyz, g_z_z_zz_yzz, g_z_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzz_xxx[k] = -g_z_z_zz_xxx[k] * ab_x + g_z_z_zz_xxxx[k];

                g_z_z_xzz_xxy[k] = -g_z_z_zz_xxy[k] * ab_x + g_z_z_zz_xxxy[k];

                g_z_z_xzz_xxz[k] = -g_z_z_zz_xxz[k] * ab_x + g_z_z_zz_xxxz[k];

                g_z_z_xzz_xyy[k] = -g_z_z_zz_xyy[k] * ab_x + g_z_z_zz_xxyy[k];

                g_z_z_xzz_xyz[k] = -g_z_z_zz_xyz[k] * ab_x + g_z_z_zz_xxyz[k];

                g_z_z_xzz_xzz[k] = -g_z_z_zz_xzz[k] * ab_x + g_z_z_zz_xxzz[k];

                g_z_z_xzz_yyy[k] = -g_z_z_zz_yyy[k] * ab_x + g_z_z_zz_xyyy[k];

                g_z_z_xzz_yyz[k] = -g_z_z_zz_yyz[k] * ab_x + g_z_z_zz_xyyz[k];

                g_z_z_xzz_yzz[k] = -g_z_z_zz_yzz[k] * ab_x + g_z_z_zz_xyzz[k];

                g_z_z_xzz_zzz[k] = -g_z_z_zz_zzz[k] * ab_x + g_z_z_zz_xzzz[k];
            }

            /// Set up 860-870 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyy_xxx = cbuffer.data(ff_geom_11_off + 860 * ccomps * dcomps);

            auto g_z_z_yyy_xxy = cbuffer.data(ff_geom_11_off + 861 * ccomps * dcomps);

            auto g_z_z_yyy_xxz = cbuffer.data(ff_geom_11_off + 862 * ccomps * dcomps);

            auto g_z_z_yyy_xyy = cbuffer.data(ff_geom_11_off + 863 * ccomps * dcomps);

            auto g_z_z_yyy_xyz = cbuffer.data(ff_geom_11_off + 864 * ccomps * dcomps);

            auto g_z_z_yyy_xzz = cbuffer.data(ff_geom_11_off + 865 * ccomps * dcomps);

            auto g_z_z_yyy_yyy = cbuffer.data(ff_geom_11_off + 866 * ccomps * dcomps);

            auto g_z_z_yyy_yyz = cbuffer.data(ff_geom_11_off + 867 * ccomps * dcomps);

            auto g_z_z_yyy_yzz = cbuffer.data(ff_geom_11_off + 868 * ccomps * dcomps);

            auto g_z_z_yyy_zzz = cbuffer.data(ff_geom_11_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yy_xxx, g_z_z_yy_xxxy, g_z_z_yy_xxy, g_z_z_yy_xxyy, g_z_z_yy_xxyz, g_z_z_yy_xxz, g_z_z_yy_xyy, g_z_z_yy_xyyy, g_z_z_yy_xyyz, g_z_z_yy_xyz, g_z_z_yy_xyzz, g_z_z_yy_xzz, g_z_z_yy_yyy, g_z_z_yy_yyyy, g_z_z_yy_yyyz, g_z_z_yy_yyz, g_z_z_yy_yyzz, g_z_z_yy_yzz, g_z_z_yy_yzzz, g_z_z_yy_zzz, g_z_z_yyy_xxx, g_z_z_yyy_xxy, g_z_z_yyy_xxz, g_z_z_yyy_xyy, g_z_z_yyy_xyz, g_z_z_yyy_xzz, g_z_z_yyy_yyy, g_z_z_yyy_yyz, g_z_z_yyy_yzz, g_z_z_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyy_xxx[k] = -g_z_z_yy_xxx[k] * ab_y + g_z_z_yy_xxxy[k];

                g_z_z_yyy_xxy[k] = -g_z_z_yy_xxy[k] * ab_y + g_z_z_yy_xxyy[k];

                g_z_z_yyy_xxz[k] = -g_z_z_yy_xxz[k] * ab_y + g_z_z_yy_xxyz[k];

                g_z_z_yyy_xyy[k] = -g_z_z_yy_xyy[k] * ab_y + g_z_z_yy_xyyy[k];

                g_z_z_yyy_xyz[k] = -g_z_z_yy_xyz[k] * ab_y + g_z_z_yy_xyyz[k];

                g_z_z_yyy_xzz[k] = -g_z_z_yy_xzz[k] * ab_y + g_z_z_yy_xyzz[k];

                g_z_z_yyy_yyy[k] = -g_z_z_yy_yyy[k] * ab_y + g_z_z_yy_yyyy[k];

                g_z_z_yyy_yyz[k] = -g_z_z_yy_yyz[k] * ab_y + g_z_z_yy_yyyz[k];

                g_z_z_yyy_yzz[k] = -g_z_z_yy_yzz[k] * ab_y + g_z_z_yy_yyzz[k];

                g_z_z_yyy_zzz[k] = -g_z_z_yy_zzz[k] * ab_y + g_z_z_yy_yzzz[k];
            }

            /// Set up 870-880 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyz_xxx = cbuffer.data(ff_geom_11_off + 870 * ccomps * dcomps);

            auto g_z_z_yyz_xxy = cbuffer.data(ff_geom_11_off + 871 * ccomps * dcomps);

            auto g_z_z_yyz_xxz = cbuffer.data(ff_geom_11_off + 872 * ccomps * dcomps);

            auto g_z_z_yyz_xyy = cbuffer.data(ff_geom_11_off + 873 * ccomps * dcomps);

            auto g_z_z_yyz_xyz = cbuffer.data(ff_geom_11_off + 874 * ccomps * dcomps);

            auto g_z_z_yyz_xzz = cbuffer.data(ff_geom_11_off + 875 * ccomps * dcomps);

            auto g_z_z_yyz_yyy = cbuffer.data(ff_geom_11_off + 876 * ccomps * dcomps);

            auto g_z_z_yyz_yyz = cbuffer.data(ff_geom_11_off + 877 * ccomps * dcomps);

            auto g_z_z_yyz_yzz = cbuffer.data(ff_geom_11_off + 878 * ccomps * dcomps);

            auto g_z_z_yyz_zzz = cbuffer.data(ff_geom_11_off + 879 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyz_xxx, g_z_z_yyz_xxy, g_z_z_yyz_xxz, g_z_z_yyz_xyy, g_z_z_yyz_xyz, g_z_z_yyz_xzz, g_z_z_yyz_yyy, g_z_z_yyz_yyz, g_z_z_yyz_yzz, g_z_z_yyz_zzz, g_z_z_yz_xxx, g_z_z_yz_xxxy, g_z_z_yz_xxy, g_z_z_yz_xxyy, g_z_z_yz_xxyz, g_z_z_yz_xxz, g_z_z_yz_xyy, g_z_z_yz_xyyy, g_z_z_yz_xyyz, g_z_z_yz_xyz, g_z_z_yz_xyzz, g_z_z_yz_xzz, g_z_z_yz_yyy, g_z_z_yz_yyyy, g_z_z_yz_yyyz, g_z_z_yz_yyz, g_z_z_yz_yyzz, g_z_z_yz_yzz, g_z_z_yz_yzzz, g_z_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyz_xxx[k] = -g_z_z_yz_xxx[k] * ab_y + g_z_z_yz_xxxy[k];

                g_z_z_yyz_xxy[k] = -g_z_z_yz_xxy[k] * ab_y + g_z_z_yz_xxyy[k];

                g_z_z_yyz_xxz[k] = -g_z_z_yz_xxz[k] * ab_y + g_z_z_yz_xxyz[k];

                g_z_z_yyz_xyy[k] = -g_z_z_yz_xyy[k] * ab_y + g_z_z_yz_xyyy[k];

                g_z_z_yyz_xyz[k] = -g_z_z_yz_xyz[k] * ab_y + g_z_z_yz_xyyz[k];

                g_z_z_yyz_xzz[k] = -g_z_z_yz_xzz[k] * ab_y + g_z_z_yz_xyzz[k];

                g_z_z_yyz_yyy[k] = -g_z_z_yz_yyy[k] * ab_y + g_z_z_yz_yyyy[k];

                g_z_z_yyz_yyz[k] = -g_z_z_yz_yyz[k] * ab_y + g_z_z_yz_yyyz[k];

                g_z_z_yyz_yzz[k] = -g_z_z_yz_yzz[k] * ab_y + g_z_z_yz_yyzz[k];

                g_z_z_yyz_zzz[k] = -g_z_z_yz_zzz[k] * ab_y + g_z_z_yz_yzzz[k];
            }

            /// Set up 880-890 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzz_xxx = cbuffer.data(ff_geom_11_off + 880 * ccomps * dcomps);

            auto g_z_z_yzz_xxy = cbuffer.data(ff_geom_11_off + 881 * ccomps * dcomps);

            auto g_z_z_yzz_xxz = cbuffer.data(ff_geom_11_off + 882 * ccomps * dcomps);

            auto g_z_z_yzz_xyy = cbuffer.data(ff_geom_11_off + 883 * ccomps * dcomps);

            auto g_z_z_yzz_xyz = cbuffer.data(ff_geom_11_off + 884 * ccomps * dcomps);

            auto g_z_z_yzz_xzz = cbuffer.data(ff_geom_11_off + 885 * ccomps * dcomps);

            auto g_z_z_yzz_yyy = cbuffer.data(ff_geom_11_off + 886 * ccomps * dcomps);

            auto g_z_z_yzz_yyz = cbuffer.data(ff_geom_11_off + 887 * ccomps * dcomps);

            auto g_z_z_yzz_yzz = cbuffer.data(ff_geom_11_off + 888 * ccomps * dcomps);

            auto g_z_z_yzz_zzz = cbuffer.data(ff_geom_11_off + 889 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzz_xxx, g_z_z_yzz_xxy, g_z_z_yzz_xxz, g_z_z_yzz_xyy, g_z_z_yzz_xyz, g_z_z_yzz_xzz, g_z_z_yzz_yyy, g_z_z_yzz_yyz, g_z_z_yzz_yzz, g_z_z_yzz_zzz, g_z_z_zz_xxx, g_z_z_zz_xxxy, g_z_z_zz_xxy, g_z_z_zz_xxyy, g_z_z_zz_xxyz, g_z_z_zz_xxz, g_z_z_zz_xyy, g_z_z_zz_xyyy, g_z_z_zz_xyyz, g_z_z_zz_xyz, g_z_z_zz_xyzz, g_z_z_zz_xzz, g_z_z_zz_yyy, g_z_z_zz_yyyy, g_z_z_zz_yyyz, g_z_z_zz_yyz, g_z_z_zz_yyzz, g_z_z_zz_yzz, g_z_z_zz_yzzz, g_z_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzz_xxx[k] = -g_z_z_zz_xxx[k] * ab_y + g_z_z_zz_xxxy[k];

                g_z_z_yzz_xxy[k] = -g_z_z_zz_xxy[k] * ab_y + g_z_z_zz_xxyy[k];

                g_z_z_yzz_xxz[k] = -g_z_z_zz_xxz[k] * ab_y + g_z_z_zz_xxyz[k];

                g_z_z_yzz_xyy[k] = -g_z_z_zz_xyy[k] * ab_y + g_z_z_zz_xyyy[k];

                g_z_z_yzz_xyz[k] = -g_z_z_zz_xyz[k] * ab_y + g_z_z_zz_xyyz[k];

                g_z_z_yzz_xzz[k] = -g_z_z_zz_xzz[k] * ab_y + g_z_z_zz_xyzz[k];

                g_z_z_yzz_yyy[k] = -g_z_z_zz_yyy[k] * ab_y + g_z_z_zz_yyyy[k];

                g_z_z_yzz_yyz[k] = -g_z_z_zz_yyz[k] * ab_y + g_z_z_zz_yyyz[k];

                g_z_z_yzz_yzz[k] = -g_z_z_zz_yzz[k] * ab_y + g_z_z_zz_yyzz[k];

                g_z_z_yzz_zzz[k] = -g_z_z_zz_zzz[k] * ab_y + g_z_z_zz_yzzz[k];
            }

            /// Set up 890-900 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzz_xxx = cbuffer.data(ff_geom_11_off + 890 * ccomps * dcomps);

            auto g_z_z_zzz_xxy = cbuffer.data(ff_geom_11_off + 891 * ccomps * dcomps);

            auto g_z_z_zzz_xxz = cbuffer.data(ff_geom_11_off + 892 * ccomps * dcomps);

            auto g_z_z_zzz_xyy = cbuffer.data(ff_geom_11_off + 893 * ccomps * dcomps);

            auto g_z_z_zzz_xyz = cbuffer.data(ff_geom_11_off + 894 * ccomps * dcomps);

            auto g_z_z_zzz_xzz = cbuffer.data(ff_geom_11_off + 895 * ccomps * dcomps);

            auto g_z_z_zzz_yyy = cbuffer.data(ff_geom_11_off + 896 * ccomps * dcomps);

            auto g_z_z_zzz_yyz = cbuffer.data(ff_geom_11_off + 897 * ccomps * dcomps);

            auto g_z_z_zzz_yzz = cbuffer.data(ff_geom_11_off + 898 * ccomps * dcomps);

            auto g_z_z_zzz_zzz = cbuffer.data(ff_geom_11_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_xxx, g_0_z_zz_xxy, g_0_z_zz_xxz, g_0_z_zz_xyy, g_0_z_zz_xyz, g_0_z_zz_xzz, g_0_z_zz_yyy, g_0_z_zz_yyz, g_0_z_zz_yzz, g_0_z_zz_zzz, g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz, g_z_z_zz_xxx, g_z_z_zz_xxxz, g_z_z_zz_xxy, g_z_z_zz_xxyz, g_z_z_zz_xxz, g_z_z_zz_xxzz, g_z_z_zz_xyy, g_z_z_zz_xyyz, g_z_z_zz_xyz, g_z_z_zz_xyzz, g_z_z_zz_xzz, g_z_z_zz_xzzz, g_z_z_zz_yyy, g_z_z_zz_yyyz, g_z_z_zz_yyz, g_z_z_zz_yyzz, g_z_z_zz_yzz, g_z_z_zz_yzzz, g_z_z_zz_zzz, g_z_z_zz_zzzz, g_z_z_zzz_xxx, g_z_z_zzz_xxy, g_z_z_zzz_xxz, g_z_z_zzz_xyy, g_z_z_zzz_xyz, g_z_z_zzz_xzz, g_z_z_zzz_yyy, g_z_z_zzz_yyz, g_z_z_zzz_yzz, g_z_z_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzz_xxx[k] = -g_0_z_zz_xxx[k] + g_z_0_zz_xxx[k] - g_z_z_zz_xxx[k] * ab_z + g_z_z_zz_xxxz[k];

                g_z_z_zzz_xxy[k] = -g_0_z_zz_xxy[k] + g_z_0_zz_xxy[k] - g_z_z_zz_xxy[k] * ab_z + g_z_z_zz_xxyz[k];

                g_z_z_zzz_xxz[k] = -g_0_z_zz_xxz[k] + g_z_0_zz_xxz[k] - g_z_z_zz_xxz[k] * ab_z + g_z_z_zz_xxzz[k];

                g_z_z_zzz_xyy[k] = -g_0_z_zz_xyy[k] + g_z_0_zz_xyy[k] - g_z_z_zz_xyy[k] * ab_z + g_z_z_zz_xyyz[k];

                g_z_z_zzz_xyz[k] = -g_0_z_zz_xyz[k] + g_z_0_zz_xyz[k] - g_z_z_zz_xyz[k] * ab_z + g_z_z_zz_xyzz[k];

                g_z_z_zzz_xzz[k] = -g_0_z_zz_xzz[k] + g_z_0_zz_xzz[k] - g_z_z_zz_xzz[k] * ab_z + g_z_z_zz_xzzz[k];

                g_z_z_zzz_yyy[k] = -g_0_z_zz_yyy[k] + g_z_0_zz_yyy[k] - g_z_z_zz_yyy[k] * ab_z + g_z_z_zz_yyyz[k];

                g_z_z_zzz_yyz[k] = -g_0_z_zz_yyz[k] + g_z_0_zz_yyz[k] - g_z_z_zz_yyz[k] * ab_z + g_z_z_zz_yyzz[k];

                g_z_z_zzz_yzz[k] = -g_0_z_zz_yzz[k] + g_z_0_zz_yzz[k] - g_z_z_zz_yzz[k] * ab_z + g_z_z_zz_yzzz[k];

                g_z_z_zzz_zzz[k] = -g_0_z_zz_zzz[k] + g_z_0_zz_zzz[k] - g_z_z_zz_zzz[k] * ab_z + g_z_z_zz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

