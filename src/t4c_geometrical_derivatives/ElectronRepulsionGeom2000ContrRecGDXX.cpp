#include "ElectronRepulsionGeom2000ContrRecGDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_gdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_gdxx,
                                            const size_t idx_geom_10_fdxx,
                                            const size_t idx_geom_20_fdxx,
                                            const size_t idx_geom_20_ffxx,
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
            /// Set up components of auxilary buffer : FDSS

            const auto fd_geom_10_off = idx_geom_10_fdxx + i * dcomps + j;

            auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxy_xx = cbuffer.data(fd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxy_xz = cbuffer.data(fd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxy_zz = cbuffer.data(fd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxz_xx = cbuffer.data(fd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xyy_xx = cbuffer.data(fd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xyy_xz = cbuffer.data(fd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xyy_zz = cbuffer.data(fd_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xyz_xx = cbuffer.data(fd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xyz_xz = cbuffer.data(fd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xyz_zz = cbuffer.data(fd_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xzz_xx = cbuffer.data(fd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_yyy_xx = cbuffer.data(fd_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yyy_xz = cbuffer.data(fd_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yyy_zz = cbuffer.data(fd_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_yyz_xx = cbuffer.data(fd_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_yyz_xz = cbuffer.data(fd_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yyz_zz = cbuffer.data(fd_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yzz_xx = cbuffer.data(fd_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yzz_xz = cbuffer.data(fd_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yzz_zz = cbuffer.data(fd_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_zzz_xx = cbuffer.data(fd_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xxx_yy = cbuffer.data(fd_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xxx_yz = cbuffer.data(fd_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xxx_zz = cbuffer.data(fd_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xxy_yy = cbuffer.data(fd_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xxy_yz = cbuffer.data(fd_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xxy_zz = cbuffer.data(fd_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xxz_yy = cbuffer.data(fd_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xxz_yz = cbuffer.data(fd_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xxz_zz = cbuffer.data(fd_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xyy_yy = cbuffer.data(fd_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xyy_yz = cbuffer.data(fd_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xyy_zz = cbuffer.data(fd_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xyz_yy = cbuffer.data(fd_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xyz_yz = cbuffer.data(fd_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xyz_zz = cbuffer.data(fd_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xzz_yy = cbuffer.data(fd_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xzz_yz = cbuffer.data(fd_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xzz_zz = cbuffer.data(fd_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_yyz_yy = cbuffer.data(fd_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_yzz_yy = cbuffer.data(fd_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_zzz_yy = cbuffer.data(fd_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 119 * ccomps * dcomps);

            auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_xxx_yy = cbuffer.data(fd_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_xxx_yz = cbuffer.data(fd_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_xxx_zz = cbuffer.data(fd_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xxy_yy = cbuffer.data(fd_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_xxy_yz = cbuffer.data(fd_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xxy_zz = cbuffer.data(fd_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xxz_yy = cbuffer.data(fd_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xxz_yz = cbuffer.data(fd_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xxz_zz = cbuffer.data(fd_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xyy_yy = cbuffer.data(fd_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xyy_yz = cbuffer.data(fd_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xyy_zz = cbuffer.data(fd_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xyz_yy = cbuffer.data(fd_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xyz_yz = cbuffer.data(fd_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xyz_zz = cbuffer.data(fd_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_xzz_yy = cbuffer.data(fd_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_xzz_yz = cbuffer.data(fd_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_xzz_zz = cbuffer.data(fd_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_yyy_zz = cbuffer.data(fd_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_yyz_zz = cbuffer.data(fd_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_yzz_zz = cbuffer.data(fd_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 179 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FDSS

            const auto fd_geom_20_off = idx_geom_20_fdxx + i * dcomps + j;

            auto g_xx_0_xxx_xx = cbuffer.data(fd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_xy = cbuffer.data(fd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_xz = cbuffer.data(fd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxx_yy = cbuffer.data(fd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxx_yz = cbuffer.data(fd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxx_zz = cbuffer.data(fd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxy_xx = cbuffer.data(fd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxy_xy = cbuffer.data(fd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxy_xz = cbuffer.data(fd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxy_yy = cbuffer.data(fd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxy_yz = cbuffer.data(fd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxy_zz = cbuffer.data(fd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxz_xx = cbuffer.data(fd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxz_xy = cbuffer.data(fd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxz_xz = cbuffer.data(fd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxz_yy = cbuffer.data(fd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxz_yz = cbuffer.data(fd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxz_zz = cbuffer.data(fd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xyy_xx = cbuffer.data(fd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xyy_xy = cbuffer.data(fd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xyy_xz = cbuffer.data(fd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xyy_yy = cbuffer.data(fd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xyy_yz = cbuffer.data(fd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xyy_zz = cbuffer.data(fd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xyz_xx = cbuffer.data(fd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xyz_xy = cbuffer.data(fd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xyz_xz = cbuffer.data(fd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xyz_yy = cbuffer.data(fd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xyz_yz = cbuffer.data(fd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xyz_zz = cbuffer.data(fd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xzz_xx = cbuffer.data(fd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xzz_xy = cbuffer.data(fd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xzz_xz = cbuffer.data(fd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xzz_yy = cbuffer.data(fd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xzz_yz = cbuffer.data(fd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xzz_zz = cbuffer.data(fd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_yyy_xx = cbuffer.data(fd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_yyy_xy = cbuffer.data(fd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_yyy_xz = cbuffer.data(fd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_yyy_yy = cbuffer.data(fd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_yyy_yz = cbuffer.data(fd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_yyy_zz = cbuffer.data(fd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_yyz_xx = cbuffer.data(fd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_yyz_xy = cbuffer.data(fd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_yyz_xz = cbuffer.data(fd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_yyz_yy = cbuffer.data(fd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yyz_yz = cbuffer.data(fd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yyz_zz = cbuffer.data(fd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yzz_xx = cbuffer.data(fd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yzz_xy = cbuffer.data(fd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yzz_xz = cbuffer.data(fd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_yzz_yy = cbuffer.data(fd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yzz_yz = cbuffer.data(fd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yzz_zz = cbuffer.data(fd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_zzz_xx = cbuffer.data(fd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_zzz_xy = cbuffer.data(fd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_zzz_xz = cbuffer.data(fd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_zzz_yy = cbuffer.data(fd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_zzz_yz = cbuffer.data(fd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_zzz_zz = cbuffer.data(fd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_xxx_xx = cbuffer.data(fd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_xxx_xy = cbuffer.data(fd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_xxx_xz = cbuffer.data(fd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_xxx_yy = cbuffer.data(fd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xxx_yz = cbuffer.data(fd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xxx_zz = cbuffer.data(fd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_xxy_xx = cbuffer.data(fd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xxy_xy = cbuffer.data(fd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xxy_xz = cbuffer.data(fd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_xxy_yy = cbuffer.data(fd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xxy_yz = cbuffer.data(fd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xxy_zz = cbuffer.data(fd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_xxz_xx = cbuffer.data(fd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xxz_xy = cbuffer.data(fd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xxz_xz = cbuffer.data(fd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_xxz_yy = cbuffer.data(fd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_xxz_yz = cbuffer.data(fd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_xxz_zz = cbuffer.data(fd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_xyy_xx = cbuffer.data(fd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_xyy_xy = cbuffer.data(fd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_xyy_xz = cbuffer.data(fd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_xyy_yy = cbuffer.data(fd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_xyy_yz = cbuffer.data(fd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_xyy_zz = cbuffer.data(fd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_xyz_xx = cbuffer.data(fd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xyz_xy = cbuffer.data(fd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xyz_xz = cbuffer.data(fd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_xyz_yy = cbuffer.data(fd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xyz_yz = cbuffer.data(fd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xyz_zz = cbuffer.data(fd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xzz_xx = cbuffer.data(fd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xzz_xy = cbuffer.data(fd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xzz_xz = cbuffer.data(fd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xzz_yy = cbuffer.data(fd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xzz_yz = cbuffer.data(fd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xzz_zz = cbuffer.data(fd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_yyy_xx = cbuffer.data(fd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_yyy_xy = cbuffer.data(fd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_yyy_xz = cbuffer.data(fd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_yyy_yy = cbuffer.data(fd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_yyy_yz = cbuffer.data(fd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_yyy_zz = cbuffer.data(fd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_yyz_xx = cbuffer.data(fd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_yyz_xy = cbuffer.data(fd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_yyz_xz = cbuffer.data(fd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_yyz_yy = cbuffer.data(fd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_yyz_yz = cbuffer.data(fd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_yyz_zz = cbuffer.data(fd_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_yzz_xx = cbuffer.data(fd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_yzz_xy = cbuffer.data(fd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_yzz_xz = cbuffer.data(fd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_yzz_yy = cbuffer.data(fd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_yzz_yz = cbuffer.data(fd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_yzz_zz = cbuffer.data(fd_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_zzz_xx = cbuffer.data(fd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_zzz_xy = cbuffer.data(fd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_zzz_xz = cbuffer.data(fd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_zzz_yy = cbuffer.data(fd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_zzz_yz = cbuffer.data(fd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_zzz_zz = cbuffer.data(fd_geom_20_off + 119 * ccomps * dcomps);

            auto g_xz_0_xxx_xx = cbuffer.data(fd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_xxx_xy = cbuffer.data(fd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_xxx_xz = cbuffer.data(fd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xz_0_xxx_yy = cbuffer.data(fd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_xxx_yz = cbuffer.data(fd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_xxx_zz = cbuffer.data(fd_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_xxy_xx = cbuffer.data(fd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_xxy_xy = cbuffer.data(fd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_xxy_xz = cbuffer.data(fd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_xxy_yy = cbuffer.data(fd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_xxy_yz = cbuffer.data(fd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_xxy_zz = cbuffer.data(fd_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_xxz_xx = cbuffer.data(fd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_xxz_xy = cbuffer.data(fd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_xxz_xz = cbuffer.data(fd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_xxz_yy = cbuffer.data(fd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_xxz_yz = cbuffer.data(fd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_xxz_zz = cbuffer.data(fd_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_xyy_xx = cbuffer.data(fd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_xyy_xy = cbuffer.data(fd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_xyy_xz = cbuffer.data(fd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_xyy_yy = cbuffer.data(fd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_xyy_yz = cbuffer.data(fd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_xyy_zz = cbuffer.data(fd_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_xyz_xx = cbuffer.data(fd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_xyz_xy = cbuffer.data(fd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_xyz_xz = cbuffer.data(fd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_xyz_yy = cbuffer.data(fd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_xyz_yz = cbuffer.data(fd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_xyz_zz = cbuffer.data(fd_geom_20_off + 149 * ccomps * dcomps);

            auto g_xz_0_xzz_xx = cbuffer.data(fd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_xzz_xy = cbuffer.data(fd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_xzz_xz = cbuffer.data(fd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_xzz_yy = cbuffer.data(fd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_xzz_yz = cbuffer.data(fd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_xzz_zz = cbuffer.data(fd_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_yyy_xx = cbuffer.data(fd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_yyy_xy = cbuffer.data(fd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_yyy_xz = cbuffer.data(fd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_yyy_yy = cbuffer.data(fd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_yyy_yz = cbuffer.data(fd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_yyy_zz = cbuffer.data(fd_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_yyz_xx = cbuffer.data(fd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_yyz_xy = cbuffer.data(fd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_yyz_xz = cbuffer.data(fd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_yyz_yy = cbuffer.data(fd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_yyz_yz = cbuffer.data(fd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_yyz_zz = cbuffer.data(fd_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_yzz_xx = cbuffer.data(fd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_yzz_xy = cbuffer.data(fd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_yzz_xz = cbuffer.data(fd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_yzz_yy = cbuffer.data(fd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_yzz_yz = cbuffer.data(fd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_yzz_zz = cbuffer.data(fd_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_zzz_xx = cbuffer.data(fd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_zzz_xy = cbuffer.data(fd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_zzz_xz = cbuffer.data(fd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_zzz_yy = cbuffer.data(fd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_zzz_yz = cbuffer.data(fd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_zzz_zz = cbuffer.data(fd_geom_20_off + 179 * ccomps * dcomps);

            auto g_yy_0_xxx_xx = cbuffer.data(fd_geom_20_off + 180 * ccomps * dcomps);

            auto g_yy_0_xxx_xy = cbuffer.data(fd_geom_20_off + 181 * ccomps * dcomps);

            auto g_yy_0_xxx_xz = cbuffer.data(fd_geom_20_off + 182 * ccomps * dcomps);

            auto g_yy_0_xxx_yy = cbuffer.data(fd_geom_20_off + 183 * ccomps * dcomps);

            auto g_yy_0_xxx_yz = cbuffer.data(fd_geom_20_off + 184 * ccomps * dcomps);

            auto g_yy_0_xxx_zz = cbuffer.data(fd_geom_20_off + 185 * ccomps * dcomps);

            auto g_yy_0_xxy_xx = cbuffer.data(fd_geom_20_off + 186 * ccomps * dcomps);

            auto g_yy_0_xxy_xy = cbuffer.data(fd_geom_20_off + 187 * ccomps * dcomps);

            auto g_yy_0_xxy_xz = cbuffer.data(fd_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_xxy_yy = cbuffer.data(fd_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_xxy_yz = cbuffer.data(fd_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_xxy_zz = cbuffer.data(fd_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_xxz_xx = cbuffer.data(fd_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_xxz_xy = cbuffer.data(fd_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_xxz_xz = cbuffer.data(fd_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_xxz_yy = cbuffer.data(fd_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_xxz_yz = cbuffer.data(fd_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_xxz_zz = cbuffer.data(fd_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_xyy_xx = cbuffer.data(fd_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_xyy_xy = cbuffer.data(fd_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_xyy_xz = cbuffer.data(fd_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_xyy_yy = cbuffer.data(fd_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_xyy_yz = cbuffer.data(fd_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_xyy_zz = cbuffer.data(fd_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_xyz_xx = cbuffer.data(fd_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_xyz_xy = cbuffer.data(fd_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_xyz_xz = cbuffer.data(fd_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_xyz_yy = cbuffer.data(fd_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_xyz_yz = cbuffer.data(fd_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_xyz_zz = cbuffer.data(fd_geom_20_off + 209 * ccomps * dcomps);

            auto g_yy_0_xzz_xx = cbuffer.data(fd_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_xzz_xy = cbuffer.data(fd_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_xzz_xz = cbuffer.data(fd_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_xzz_yy = cbuffer.data(fd_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_xzz_yz = cbuffer.data(fd_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_xzz_zz = cbuffer.data(fd_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_yyy_xx = cbuffer.data(fd_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_yyy_xy = cbuffer.data(fd_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_yyy_xz = cbuffer.data(fd_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_yyy_yy = cbuffer.data(fd_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_yyy_yz = cbuffer.data(fd_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_yyy_zz = cbuffer.data(fd_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_yyz_xx = cbuffer.data(fd_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_yyz_xy = cbuffer.data(fd_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_yyz_xz = cbuffer.data(fd_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_yyz_yy = cbuffer.data(fd_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_yyz_yz = cbuffer.data(fd_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_yyz_zz = cbuffer.data(fd_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_yzz_xx = cbuffer.data(fd_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_yzz_xy = cbuffer.data(fd_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_yzz_xz = cbuffer.data(fd_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_yzz_yy = cbuffer.data(fd_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_yzz_yz = cbuffer.data(fd_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_yzz_zz = cbuffer.data(fd_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_zzz_xx = cbuffer.data(fd_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_zzz_xy = cbuffer.data(fd_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_zzz_xz = cbuffer.data(fd_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_zzz_yy = cbuffer.data(fd_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_zzz_yz = cbuffer.data(fd_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_zzz_zz = cbuffer.data(fd_geom_20_off + 239 * ccomps * dcomps);

            auto g_yz_0_xxx_xx = cbuffer.data(fd_geom_20_off + 240 * ccomps * dcomps);

            auto g_yz_0_xxx_xy = cbuffer.data(fd_geom_20_off + 241 * ccomps * dcomps);

            auto g_yz_0_xxx_xz = cbuffer.data(fd_geom_20_off + 242 * ccomps * dcomps);

            auto g_yz_0_xxx_yy = cbuffer.data(fd_geom_20_off + 243 * ccomps * dcomps);

            auto g_yz_0_xxx_yz = cbuffer.data(fd_geom_20_off + 244 * ccomps * dcomps);

            auto g_yz_0_xxx_zz = cbuffer.data(fd_geom_20_off + 245 * ccomps * dcomps);

            auto g_yz_0_xxy_xx = cbuffer.data(fd_geom_20_off + 246 * ccomps * dcomps);

            auto g_yz_0_xxy_xy = cbuffer.data(fd_geom_20_off + 247 * ccomps * dcomps);

            auto g_yz_0_xxy_xz = cbuffer.data(fd_geom_20_off + 248 * ccomps * dcomps);

            auto g_yz_0_xxy_yy = cbuffer.data(fd_geom_20_off + 249 * ccomps * dcomps);

            auto g_yz_0_xxy_yz = cbuffer.data(fd_geom_20_off + 250 * ccomps * dcomps);

            auto g_yz_0_xxy_zz = cbuffer.data(fd_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_xxz_xx = cbuffer.data(fd_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_xxz_xy = cbuffer.data(fd_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_xxz_xz = cbuffer.data(fd_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_xxz_yy = cbuffer.data(fd_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_xxz_yz = cbuffer.data(fd_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_xxz_zz = cbuffer.data(fd_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_xyy_xx = cbuffer.data(fd_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_xyy_xy = cbuffer.data(fd_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_xyy_xz = cbuffer.data(fd_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_xyy_yy = cbuffer.data(fd_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_xyy_yz = cbuffer.data(fd_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_xyy_zz = cbuffer.data(fd_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_xyz_xx = cbuffer.data(fd_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_xyz_xy = cbuffer.data(fd_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_xyz_xz = cbuffer.data(fd_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_xyz_yy = cbuffer.data(fd_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_xyz_yz = cbuffer.data(fd_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_xyz_zz = cbuffer.data(fd_geom_20_off + 269 * ccomps * dcomps);

            auto g_yz_0_xzz_xx = cbuffer.data(fd_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_xzz_xy = cbuffer.data(fd_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_xzz_xz = cbuffer.data(fd_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_xzz_yy = cbuffer.data(fd_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_xzz_yz = cbuffer.data(fd_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_xzz_zz = cbuffer.data(fd_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_yyy_xx = cbuffer.data(fd_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_yyy_xy = cbuffer.data(fd_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_yyy_xz = cbuffer.data(fd_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_yyy_yy = cbuffer.data(fd_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_yyy_yz = cbuffer.data(fd_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_yyy_zz = cbuffer.data(fd_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_yyz_xx = cbuffer.data(fd_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_yyz_xy = cbuffer.data(fd_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_yyz_xz = cbuffer.data(fd_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_yyz_yy = cbuffer.data(fd_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_yyz_yz = cbuffer.data(fd_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_yyz_zz = cbuffer.data(fd_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_yzz_xx = cbuffer.data(fd_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_yzz_xy = cbuffer.data(fd_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_yzz_xz = cbuffer.data(fd_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_yzz_yy = cbuffer.data(fd_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_yzz_yz = cbuffer.data(fd_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_yzz_zz = cbuffer.data(fd_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_zzz_xx = cbuffer.data(fd_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_zzz_xy = cbuffer.data(fd_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_zzz_xz = cbuffer.data(fd_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_zzz_yy = cbuffer.data(fd_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_zzz_yz = cbuffer.data(fd_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_zzz_zz = cbuffer.data(fd_geom_20_off + 299 * ccomps * dcomps);

            auto g_zz_0_xxx_xx = cbuffer.data(fd_geom_20_off + 300 * ccomps * dcomps);

            auto g_zz_0_xxx_xy = cbuffer.data(fd_geom_20_off + 301 * ccomps * dcomps);

            auto g_zz_0_xxx_xz = cbuffer.data(fd_geom_20_off + 302 * ccomps * dcomps);

            auto g_zz_0_xxx_yy = cbuffer.data(fd_geom_20_off + 303 * ccomps * dcomps);

            auto g_zz_0_xxx_yz = cbuffer.data(fd_geom_20_off + 304 * ccomps * dcomps);

            auto g_zz_0_xxx_zz = cbuffer.data(fd_geom_20_off + 305 * ccomps * dcomps);

            auto g_zz_0_xxy_xx = cbuffer.data(fd_geom_20_off + 306 * ccomps * dcomps);

            auto g_zz_0_xxy_xy = cbuffer.data(fd_geom_20_off + 307 * ccomps * dcomps);

            auto g_zz_0_xxy_xz = cbuffer.data(fd_geom_20_off + 308 * ccomps * dcomps);

            auto g_zz_0_xxy_yy = cbuffer.data(fd_geom_20_off + 309 * ccomps * dcomps);

            auto g_zz_0_xxy_yz = cbuffer.data(fd_geom_20_off + 310 * ccomps * dcomps);

            auto g_zz_0_xxy_zz = cbuffer.data(fd_geom_20_off + 311 * ccomps * dcomps);

            auto g_zz_0_xxz_xx = cbuffer.data(fd_geom_20_off + 312 * ccomps * dcomps);

            auto g_zz_0_xxz_xy = cbuffer.data(fd_geom_20_off + 313 * ccomps * dcomps);

            auto g_zz_0_xxz_xz = cbuffer.data(fd_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_xxz_yy = cbuffer.data(fd_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_xxz_yz = cbuffer.data(fd_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_xxz_zz = cbuffer.data(fd_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_xyy_xx = cbuffer.data(fd_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_xyy_xy = cbuffer.data(fd_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_xyy_xz = cbuffer.data(fd_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_xyy_yy = cbuffer.data(fd_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_xyy_yz = cbuffer.data(fd_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_xyy_zz = cbuffer.data(fd_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_xyz_xx = cbuffer.data(fd_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_xyz_xy = cbuffer.data(fd_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_xyz_xz = cbuffer.data(fd_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_xyz_yy = cbuffer.data(fd_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_xyz_yz = cbuffer.data(fd_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_xyz_zz = cbuffer.data(fd_geom_20_off + 329 * ccomps * dcomps);

            auto g_zz_0_xzz_xx = cbuffer.data(fd_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_xzz_xy = cbuffer.data(fd_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_xzz_xz = cbuffer.data(fd_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_xzz_yy = cbuffer.data(fd_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_xzz_yz = cbuffer.data(fd_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_xzz_zz = cbuffer.data(fd_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_yyy_xx = cbuffer.data(fd_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_yyy_xy = cbuffer.data(fd_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_yyy_xz = cbuffer.data(fd_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_yyy_yy = cbuffer.data(fd_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_yyy_yz = cbuffer.data(fd_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_yyy_zz = cbuffer.data(fd_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_yyz_xx = cbuffer.data(fd_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_yyz_xy = cbuffer.data(fd_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_yyz_xz = cbuffer.data(fd_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_yyz_yy = cbuffer.data(fd_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_yyz_yz = cbuffer.data(fd_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_yyz_zz = cbuffer.data(fd_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_yzz_xx = cbuffer.data(fd_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_yzz_xy = cbuffer.data(fd_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_yzz_xz = cbuffer.data(fd_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_yzz_yy = cbuffer.data(fd_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_yzz_yz = cbuffer.data(fd_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_yzz_zz = cbuffer.data(fd_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_zzz_xx = cbuffer.data(fd_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_zzz_xy = cbuffer.data(fd_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_zzz_xz = cbuffer.data(fd_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_zzz_yy = cbuffer.data(fd_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_zzz_yz = cbuffer.data(fd_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_zzz_zz = cbuffer.data(fd_geom_20_off + 359 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_20_off = idx_geom_20_ffxx + i * dcomps + j;

            auto g_xx_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 359 * ccomps * dcomps);

            auto g_yy_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 360 * ccomps * dcomps);

            auto g_yy_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 361 * ccomps * dcomps);

            auto g_yy_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 362 * ccomps * dcomps);

            auto g_yy_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 363 * ccomps * dcomps);

            auto g_yy_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 364 * ccomps * dcomps);

            auto g_yy_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 365 * ccomps * dcomps);

            auto g_yy_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 366 * ccomps * dcomps);

            auto g_yy_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 367 * ccomps * dcomps);

            auto g_yy_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 368 * ccomps * dcomps);

            auto g_yy_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 369 * ccomps * dcomps);

            auto g_yy_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 370 * ccomps * dcomps);

            auto g_yy_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 371 * ccomps * dcomps);

            auto g_yy_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 372 * ccomps * dcomps);

            auto g_yy_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 373 * ccomps * dcomps);

            auto g_yy_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 374 * ccomps * dcomps);

            auto g_yy_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 375 * ccomps * dcomps);

            auto g_yy_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 376 * ccomps * dcomps);

            auto g_yy_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 419 * ccomps * dcomps);

            auto g_yz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 425 * ccomps * dcomps);

            auto g_yz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 449 * ccomps * dcomps);

            auto g_yz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 450 * ccomps * dcomps);

            auto g_yz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 451 * ccomps * dcomps);

            auto g_yz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 452 * ccomps * dcomps);

            auto g_yz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 453 * ccomps * dcomps);

            auto g_yz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 454 * ccomps * dcomps);

            auto g_yz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 455 * ccomps * dcomps);

            auto g_yz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 456 * ccomps * dcomps);

            auto g_yz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 457 * ccomps * dcomps);

            auto g_yz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 458 * ccomps * dcomps);

            auto g_yz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 459 * ccomps * dcomps);

            auto g_yz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 460 * ccomps * dcomps);

            auto g_yz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 461 * ccomps * dcomps);

            auto g_yz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 462 * ccomps * dcomps);

            auto g_yz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 463 * ccomps * dcomps);

            auto g_yz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 464 * ccomps * dcomps);

            auto g_yz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 465 * ccomps * dcomps);

            auto g_yz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 466 * ccomps * dcomps);

            auto g_yz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 467 * ccomps * dcomps);

            auto g_yz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 468 * ccomps * dcomps);

            auto g_yz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 469 * ccomps * dcomps);

            auto g_yz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 470 * ccomps * dcomps);

            auto g_yz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 471 * ccomps * dcomps);

            auto g_yz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 472 * ccomps * dcomps);

            auto g_yz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 473 * ccomps * dcomps);

            auto g_yz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 474 * ccomps * dcomps);

            auto g_yz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 475 * ccomps * dcomps);

            auto g_yz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 476 * ccomps * dcomps);

            auto g_yz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 477 * ccomps * dcomps);

            auto g_yz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 478 * ccomps * dcomps);

            auto g_yz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 479 * ccomps * dcomps);

            auto g_yz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 480 * ccomps * dcomps);

            auto g_yz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 481 * ccomps * dcomps);

            auto g_yz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 482 * ccomps * dcomps);

            auto g_yz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 483 * ccomps * dcomps);

            auto g_yz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 484 * ccomps * dcomps);

            auto g_yz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 485 * ccomps * dcomps);

            auto g_yz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 486 * ccomps * dcomps);

            auto g_yz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 487 * ccomps * dcomps);

            auto g_yz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 488 * ccomps * dcomps);

            auto g_yz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 489 * ccomps * dcomps);

            auto g_yz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 490 * ccomps * dcomps);

            auto g_yz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 491 * ccomps * dcomps);

            auto g_yz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 492 * ccomps * dcomps);

            auto g_yz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 493 * ccomps * dcomps);

            auto g_yz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 494 * ccomps * dcomps);

            auto g_yz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 495 * ccomps * dcomps);

            auto g_yz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 496 * ccomps * dcomps);

            auto g_yz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 497 * ccomps * dcomps);

            auto g_yz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 498 * ccomps * dcomps);

            auto g_yz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_xxx_xxx = cbuffer.data(ff_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_xxx_xxy = cbuffer.data(ff_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_xxx_xxz = cbuffer.data(ff_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_xxx_xyy = cbuffer.data(ff_geom_20_off + 503 * ccomps * dcomps);

            auto g_zz_0_xxx_xyz = cbuffer.data(ff_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_xxx_xzz = cbuffer.data(ff_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_xxx_yyy = cbuffer.data(ff_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_xxx_yyz = cbuffer.data(ff_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_xxx_yzz = cbuffer.data(ff_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_xxx_zzz = cbuffer.data(ff_geom_20_off + 509 * ccomps * dcomps);

            auto g_zz_0_xxy_xxx = cbuffer.data(ff_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_xxy_xxy = cbuffer.data(ff_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_xxy_xxz = cbuffer.data(ff_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_xxy_xyy = cbuffer.data(ff_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_xxy_xyz = cbuffer.data(ff_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_xxy_xzz = cbuffer.data(ff_geom_20_off + 515 * ccomps * dcomps);

            auto g_zz_0_xxy_yyy = cbuffer.data(ff_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_xxy_yyz = cbuffer.data(ff_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_xxy_yzz = cbuffer.data(ff_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_xxy_zzz = cbuffer.data(ff_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_xxz_xxx = cbuffer.data(ff_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_xxz_xxy = cbuffer.data(ff_geom_20_off + 521 * ccomps * dcomps);

            auto g_zz_0_xxz_xxz = cbuffer.data(ff_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_xxz_xyy = cbuffer.data(ff_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_xxz_xyz = cbuffer.data(ff_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_xxz_xzz = cbuffer.data(ff_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_xxz_yyy = cbuffer.data(ff_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_xxz_yyz = cbuffer.data(ff_geom_20_off + 527 * ccomps * dcomps);

            auto g_zz_0_xxz_yzz = cbuffer.data(ff_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_xxz_zzz = cbuffer.data(ff_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_xyy_xxx = cbuffer.data(ff_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_xyy_xxy = cbuffer.data(ff_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_xyy_xxz = cbuffer.data(ff_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_xyy_xyy = cbuffer.data(ff_geom_20_off + 533 * ccomps * dcomps);

            auto g_zz_0_xyy_xyz = cbuffer.data(ff_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_xyy_xzz = cbuffer.data(ff_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_xyy_yyy = cbuffer.data(ff_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_xyy_yyz = cbuffer.data(ff_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_xyy_yzz = cbuffer.data(ff_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_xyy_zzz = cbuffer.data(ff_geom_20_off + 539 * ccomps * dcomps);

            auto g_zz_0_xyz_xxx = cbuffer.data(ff_geom_20_off + 540 * ccomps * dcomps);

            auto g_zz_0_xyz_xxy = cbuffer.data(ff_geom_20_off + 541 * ccomps * dcomps);

            auto g_zz_0_xyz_xxz = cbuffer.data(ff_geom_20_off + 542 * ccomps * dcomps);

            auto g_zz_0_xyz_xyy = cbuffer.data(ff_geom_20_off + 543 * ccomps * dcomps);

            auto g_zz_0_xyz_xyz = cbuffer.data(ff_geom_20_off + 544 * ccomps * dcomps);

            auto g_zz_0_xyz_xzz = cbuffer.data(ff_geom_20_off + 545 * ccomps * dcomps);

            auto g_zz_0_xyz_yyy = cbuffer.data(ff_geom_20_off + 546 * ccomps * dcomps);

            auto g_zz_0_xyz_yyz = cbuffer.data(ff_geom_20_off + 547 * ccomps * dcomps);

            auto g_zz_0_xyz_yzz = cbuffer.data(ff_geom_20_off + 548 * ccomps * dcomps);

            auto g_zz_0_xyz_zzz = cbuffer.data(ff_geom_20_off + 549 * ccomps * dcomps);

            auto g_zz_0_xzz_xxx = cbuffer.data(ff_geom_20_off + 550 * ccomps * dcomps);

            auto g_zz_0_xzz_xxy = cbuffer.data(ff_geom_20_off + 551 * ccomps * dcomps);

            auto g_zz_0_xzz_xxz = cbuffer.data(ff_geom_20_off + 552 * ccomps * dcomps);

            auto g_zz_0_xzz_xyy = cbuffer.data(ff_geom_20_off + 553 * ccomps * dcomps);

            auto g_zz_0_xzz_xyz = cbuffer.data(ff_geom_20_off + 554 * ccomps * dcomps);

            auto g_zz_0_xzz_xzz = cbuffer.data(ff_geom_20_off + 555 * ccomps * dcomps);

            auto g_zz_0_xzz_yyy = cbuffer.data(ff_geom_20_off + 556 * ccomps * dcomps);

            auto g_zz_0_xzz_yyz = cbuffer.data(ff_geom_20_off + 557 * ccomps * dcomps);

            auto g_zz_0_xzz_yzz = cbuffer.data(ff_geom_20_off + 558 * ccomps * dcomps);

            auto g_zz_0_xzz_zzz = cbuffer.data(ff_geom_20_off + 559 * ccomps * dcomps);

            auto g_zz_0_yyy_xxx = cbuffer.data(ff_geom_20_off + 560 * ccomps * dcomps);

            auto g_zz_0_yyy_xxy = cbuffer.data(ff_geom_20_off + 561 * ccomps * dcomps);

            auto g_zz_0_yyy_xxz = cbuffer.data(ff_geom_20_off + 562 * ccomps * dcomps);

            auto g_zz_0_yyy_xyy = cbuffer.data(ff_geom_20_off + 563 * ccomps * dcomps);

            auto g_zz_0_yyy_xyz = cbuffer.data(ff_geom_20_off + 564 * ccomps * dcomps);

            auto g_zz_0_yyy_xzz = cbuffer.data(ff_geom_20_off + 565 * ccomps * dcomps);

            auto g_zz_0_yyy_yyy = cbuffer.data(ff_geom_20_off + 566 * ccomps * dcomps);

            auto g_zz_0_yyy_yyz = cbuffer.data(ff_geom_20_off + 567 * ccomps * dcomps);

            auto g_zz_0_yyy_yzz = cbuffer.data(ff_geom_20_off + 568 * ccomps * dcomps);

            auto g_zz_0_yyy_zzz = cbuffer.data(ff_geom_20_off + 569 * ccomps * dcomps);

            auto g_zz_0_yyz_xxx = cbuffer.data(ff_geom_20_off + 570 * ccomps * dcomps);

            auto g_zz_0_yyz_xxy = cbuffer.data(ff_geom_20_off + 571 * ccomps * dcomps);

            auto g_zz_0_yyz_xxz = cbuffer.data(ff_geom_20_off + 572 * ccomps * dcomps);

            auto g_zz_0_yyz_xyy = cbuffer.data(ff_geom_20_off + 573 * ccomps * dcomps);

            auto g_zz_0_yyz_xyz = cbuffer.data(ff_geom_20_off + 574 * ccomps * dcomps);

            auto g_zz_0_yyz_xzz = cbuffer.data(ff_geom_20_off + 575 * ccomps * dcomps);

            auto g_zz_0_yyz_yyy = cbuffer.data(ff_geom_20_off + 576 * ccomps * dcomps);

            auto g_zz_0_yyz_yyz = cbuffer.data(ff_geom_20_off + 577 * ccomps * dcomps);

            auto g_zz_0_yyz_yzz = cbuffer.data(ff_geom_20_off + 578 * ccomps * dcomps);

            auto g_zz_0_yyz_zzz = cbuffer.data(ff_geom_20_off + 579 * ccomps * dcomps);

            auto g_zz_0_yzz_xxx = cbuffer.data(ff_geom_20_off + 580 * ccomps * dcomps);

            auto g_zz_0_yzz_xxy = cbuffer.data(ff_geom_20_off + 581 * ccomps * dcomps);

            auto g_zz_0_yzz_xxz = cbuffer.data(ff_geom_20_off + 582 * ccomps * dcomps);

            auto g_zz_0_yzz_xyy = cbuffer.data(ff_geom_20_off + 583 * ccomps * dcomps);

            auto g_zz_0_yzz_xyz = cbuffer.data(ff_geom_20_off + 584 * ccomps * dcomps);

            auto g_zz_0_yzz_xzz = cbuffer.data(ff_geom_20_off + 585 * ccomps * dcomps);

            auto g_zz_0_yzz_yyy = cbuffer.data(ff_geom_20_off + 586 * ccomps * dcomps);

            auto g_zz_0_yzz_yyz = cbuffer.data(ff_geom_20_off + 587 * ccomps * dcomps);

            auto g_zz_0_yzz_yzz = cbuffer.data(ff_geom_20_off + 588 * ccomps * dcomps);

            auto g_zz_0_yzz_zzz = cbuffer.data(ff_geom_20_off + 589 * ccomps * dcomps);

            auto g_zz_0_zzz_xxx = cbuffer.data(ff_geom_20_off + 590 * ccomps * dcomps);

            auto g_zz_0_zzz_xxy = cbuffer.data(ff_geom_20_off + 591 * ccomps * dcomps);

            auto g_zz_0_zzz_xxz = cbuffer.data(ff_geom_20_off + 592 * ccomps * dcomps);

            auto g_zz_0_zzz_xyy = cbuffer.data(ff_geom_20_off + 593 * ccomps * dcomps);

            auto g_zz_0_zzz_xyz = cbuffer.data(ff_geom_20_off + 594 * ccomps * dcomps);

            auto g_zz_0_zzz_xzz = cbuffer.data(ff_geom_20_off + 595 * ccomps * dcomps);

            auto g_zz_0_zzz_yyy = cbuffer.data(ff_geom_20_off + 596 * ccomps * dcomps);

            auto g_zz_0_zzz_yyz = cbuffer.data(ff_geom_20_off + 597 * ccomps * dcomps);

            auto g_zz_0_zzz_yzz = cbuffer.data(ff_geom_20_off + 598 * ccomps * dcomps);

            auto g_zz_0_zzz_zzz = cbuffer.data(ff_geom_20_off + 599 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gdxx

            const auto gd_geom_20_off = idx_geom_20_gdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, g_xx_0_xxx_xx, g_xx_0_xxx_xxx, g_xx_0_xxx_xxy, g_xx_0_xxx_xxz, g_xx_0_xxx_xy, g_xx_0_xxx_xyy, g_xx_0_xxx_xyz, g_xx_0_xxx_xz, g_xx_0_xxx_xzz, g_xx_0_xxx_yy, g_xx_0_xxx_yz, g_xx_0_xxx_zz, g_xx_0_xxxx_xx, g_xx_0_xxxx_xy, g_xx_0_xxxx_xz, g_xx_0_xxxx_yy, g_xx_0_xxxx_yz, g_xx_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxx_xx[k] = -2.0 * g_x_0_xxx_xx[k] - g_xx_0_xxx_xx[k] * ab_x + g_xx_0_xxx_xxx[k];

                g_xx_0_xxxx_xy[k] = -2.0 * g_x_0_xxx_xy[k] - g_xx_0_xxx_xy[k] * ab_x + g_xx_0_xxx_xxy[k];

                g_xx_0_xxxx_xz[k] = -2.0 * g_x_0_xxx_xz[k] - g_xx_0_xxx_xz[k] * ab_x + g_xx_0_xxx_xxz[k];

                g_xx_0_xxxx_yy[k] = -2.0 * g_x_0_xxx_yy[k] - g_xx_0_xxx_yy[k] * ab_x + g_xx_0_xxx_xyy[k];

                g_xx_0_xxxx_yz[k] = -2.0 * g_x_0_xxx_yz[k] - g_xx_0_xxx_yz[k] * ab_x + g_xx_0_xxx_xyz[k];

                g_xx_0_xxxx_zz[k] = -2.0 * g_x_0_xxx_zz[k] - g_xx_0_xxx_zz[k] * ab_x + g_xx_0_xxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_xx, g_xx_0_xxx_xxy, g_xx_0_xxx_xy, g_xx_0_xxx_xyy, g_xx_0_xxx_xyz, g_xx_0_xxx_xz, g_xx_0_xxx_yy, g_xx_0_xxx_yyy, g_xx_0_xxx_yyz, g_xx_0_xxx_yz, g_xx_0_xxx_yzz, g_xx_0_xxx_zz, g_xx_0_xxxy_xx, g_xx_0_xxxy_xy, g_xx_0_xxxy_xz, g_xx_0_xxxy_yy, g_xx_0_xxxy_yz, g_xx_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxy_xx[k] = -g_xx_0_xxx_xx[k] * ab_y + g_xx_0_xxx_xxy[k];

                g_xx_0_xxxy_xy[k] = -g_xx_0_xxx_xy[k] * ab_y + g_xx_0_xxx_xyy[k];

                g_xx_0_xxxy_xz[k] = -g_xx_0_xxx_xz[k] * ab_y + g_xx_0_xxx_xyz[k];

                g_xx_0_xxxy_yy[k] = -g_xx_0_xxx_yy[k] * ab_y + g_xx_0_xxx_yyy[k];

                g_xx_0_xxxy_yz[k] = -g_xx_0_xxx_yz[k] * ab_y + g_xx_0_xxx_yyz[k];

                g_xx_0_xxxy_zz[k] = -g_xx_0_xxx_zz[k] * ab_y + g_xx_0_xxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_xx, g_xx_0_xxx_xxz, g_xx_0_xxx_xy, g_xx_0_xxx_xyz, g_xx_0_xxx_xz, g_xx_0_xxx_xzz, g_xx_0_xxx_yy, g_xx_0_xxx_yyz, g_xx_0_xxx_yz, g_xx_0_xxx_yzz, g_xx_0_xxx_zz, g_xx_0_xxx_zzz, g_xx_0_xxxz_xx, g_xx_0_xxxz_xy, g_xx_0_xxxz_xz, g_xx_0_xxxz_yy, g_xx_0_xxxz_yz, g_xx_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxz_xx[k] = -g_xx_0_xxx_xx[k] * ab_z + g_xx_0_xxx_xxz[k];

                g_xx_0_xxxz_xy[k] = -g_xx_0_xxx_xy[k] * ab_z + g_xx_0_xxx_xyz[k];

                g_xx_0_xxxz_xz[k] = -g_xx_0_xxx_xz[k] * ab_z + g_xx_0_xxx_xzz[k];

                g_xx_0_xxxz_yy[k] = -g_xx_0_xxx_yy[k] * ab_z + g_xx_0_xxx_yyz[k];

                g_xx_0_xxxz_yz[k] = -g_xx_0_xxx_yz[k] * ab_z + g_xx_0_xxx_yzz[k];

                g_xx_0_xxxz_zz[k] = -g_xx_0_xxx_zz[k] * ab_z + g_xx_0_xxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxy_xx, g_xx_0_xxy_xxy, g_xx_0_xxy_xy, g_xx_0_xxy_xyy, g_xx_0_xxy_xyz, g_xx_0_xxy_xz, g_xx_0_xxy_yy, g_xx_0_xxy_yyy, g_xx_0_xxy_yyz, g_xx_0_xxy_yz, g_xx_0_xxy_yzz, g_xx_0_xxy_zz, g_xx_0_xxyy_xx, g_xx_0_xxyy_xy, g_xx_0_xxyy_xz, g_xx_0_xxyy_yy, g_xx_0_xxyy_yz, g_xx_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyy_xx[k] = -g_xx_0_xxy_xx[k] * ab_y + g_xx_0_xxy_xxy[k];

                g_xx_0_xxyy_xy[k] = -g_xx_0_xxy_xy[k] * ab_y + g_xx_0_xxy_xyy[k];

                g_xx_0_xxyy_xz[k] = -g_xx_0_xxy_xz[k] * ab_y + g_xx_0_xxy_xyz[k];

                g_xx_0_xxyy_yy[k] = -g_xx_0_xxy_yy[k] * ab_y + g_xx_0_xxy_yyy[k];

                g_xx_0_xxyy_yz[k] = -g_xx_0_xxy_yz[k] * ab_y + g_xx_0_xxy_yyz[k];

                g_xx_0_xxyy_zz[k] = -g_xx_0_xxy_zz[k] * ab_y + g_xx_0_xxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyz_xx, g_xx_0_xxyz_xy, g_xx_0_xxyz_xz, g_xx_0_xxyz_yy, g_xx_0_xxyz_yz, g_xx_0_xxyz_zz, g_xx_0_xxz_xx, g_xx_0_xxz_xxy, g_xx_0_xxz_xy, g_xx_0_xxz_xyy, g_xx_0_xxz_xyz, g_xx_0_xxz_xz, g_xx_0_xxz_yy, g_xx_0_xxz_yyy, g_xx_0_xxz_yyz, g_xx_0_xxz_yz, g_xx_0_xxz_yzz, g_xx_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyz_xx[k] = -g_xx_0_xxz_xx[k] * ab_y + g_xx_0_xxz_xxy[k];

                g_xx_0_xxyz_xy[k] = -g_xx_0_xxz_xy[k] * ab_y + g_xx_0_xxz_xyy[k];

                g_xx_0_xxyz_xz[k] = -g_xx_0_xxz_xz[k] * ab_y + g_xx_0_xxz_xyz[k];

                g_xx_0_xxyz_yy[k] = -g_xx_0_xxz_yy[k] * ab_y + g_xx_0_xxz_yyy[k];

                g_xx_0_xxyz_yz[k] = -g_xx_0_xxz_yz[k] * ab_y + g_xx_0_xxz_yyz[k];

                g_xx_0_xxyz_zz[k] = -g_xx_0_xxz_zz[k] * ab_y + g_xx_0_xxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxz_xx, g_xx_0_xxz_xxz, g_xx_0_xxz_xy, g_xx_0_xxz_xyz, g_xx_0_xxz_xz, g_xx_0_xxz_xzz, g_xx_0_xxz_yy, g_xx_0_xxz_yyz, g_xx_0_xxz_yz, g_xx_0_xxz_yzz, g_xx_0_xxz_zz, g_xx_0_xxz_zzz, g_xx_0_xxzz_xx, g_xx_0_xxzz_xy, g_xx_0_xxzz_xz, g_xx_0_xxzz_yy, g_xx_0_xxzz_yz, g_xx_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzz_xx[k] = -g_xx_0_xxz_xx[k] * ab_z + g_xx_0_xxz_xxz[k];

                g_xx_0_xxzz_xy[k] = -g_xx_0_xxz_xy[k] * ab_z + g_xx_0_xxz_xyz[k];

                g_xx_0_xxzz_xz[k] = -g_xx_0_xxz_xz[k] * ab_z + g_xx_0_xxz_xzz[k];

                g_xx_0_xxzz_yy[k] = -g_xx_0_xxz_yy[k] * ab_z + g_xx_0_xxz_yyz[k];

                g_xx_0_xxzz_yz[k] = -g_xx_0_xxz_yz[k] * ab_z + g_xx_0_xxz_yzz[k];

                g_xx_0_xxzz_zz[k] = -g_xx_0_xxz_zz[k] * ab_z + g_xx_0_xxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyy_xx, g_xx_0_xyy_xxy, g_xx_0_xyy_xy, g_xx_0_xyy_xyy, g_xx_0_xyy_xyz, g_xx_0_xyy_xz, g_xx_0_xyy_yy, g_xx_0_xyy_yyy, g_xx_0_xyy_yyz, g_xx_0_xyy_yz, g_xx_0_xyy_yzz, g_xx_0_xyy_zz, g_xx_0_xyyy_xx, g_xx_0_xyyy_xy, g_xx_0_xyyy_xz, g_xx_0_xyyy_yy, g_xx_0_xyyy_yz, g_xx_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyy_xx[k] = -g_xx_0_xyy_xx[k] * ab_y + g_xx_0_xyy_xxy[k];

                g_xx_0_xyyy_xy[k] = -g_xx_0_xyy_xy[k] * ab_y + g_xx_0_xyy_xyy[k];

                g_xx_0_xyyy_xz[k] = -g_xx_0_xyy_xz[k] * ab_y + g_xx_0_xyy_xyz[k];

                g_xx_0_xyyy_yy[k] = -g_xx_0_xyy_yy[k] * ab_y + g_xx_0_xyy_yyy[k];

                g_xx_0_xyyy_yz[k] = -g_xx_0_xyy_yz[k] * ab_y + g_xx_0_xyy_yyz[k];

                g_xx_0_xyyy_zz[k] = -g_xx_0_xyy_zz[k] * ab_y + g_xx_0_xyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyz_xx, g_xx_0_xyyz_xy, g_xx_0_xyyz_xz, g_xx_0_xyyz_yy, g_xx_0_xyyz_yz, g_xx_0_xyyz_zz, g_xx_0_xyz_xx, g_xx_0_xyz_xxy, g_xx_0_xyz_xy, g_xx_0_xyz_xyy, g_xx_0_xyz_xyz, g_xx_0_xyz_xz, g_xx_0_xyz_yy, g_xx_0_xyz_yyy, g_xx_0_xyz_yyz, g_xx_0_xyz_yz, g_xx_0_xyz_yzz, g_xx_0_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyz_xx[k] = -g_xx_0_xyz_xx[k] * ab_y + g_xx_0_xyz_xxy[k];

                g_xx_0_xyyz_xy[k] = -g_xx_0_xyz_xy[k] * ab_y + g_xx_0_xyz_xyy[k];

                g_xx_0_xyyz_xz[k] = -g_xx_0_xyz_xz[k] * ab_y + g_xx_0_xyz_xyz[k];

                g_xx_0_xyyz_yy[k] = -g_xx_0_xyz_yy[k] * ab_y + g_xx_0_xyz_yyy[k];

                g_xx_0_xyyz_yz[k] = -g_xx_0_xyz_yz[k] * ab_y + g_xx_0_xyz_yyz[k];

                g_xx_0_xyyz_zz[k] = -g_xx_0_xyz_zz[k] * ab_y + g_xx_0_xyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzz_xx, g_xx_0_xyzz_xy, g_xx_0_xyzz_xz, g_xx_0_xyzz_yy, g_xx_0_xyzz_yz, g_xx_0_xyzz_zz, g_xx_0_xzz_xx, g_xx_0_xzz_xxy, g_xx_0_xzz_xy, g_xx_0_xzz_xyy, g_xx_0_xzz_xyz, g_xx_0_xzz_xz, g_xx_0_xzz_yy, g_xx_0_xzz_yyy, g_xx_0_xzz_yyz, g_xx_0_xzz_yz, g_xx_0_xzz_yzz, g_xx_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzz_xx[k] = -g_xx_0_xzz_xx[k] * ab_y + g_xx_0_xzz_xxy[k];

                g_xx_0_xyzz_xy[k] = -g_xx_0_xzz_xy[k] * ab_y + g_xx_0_xzz_xyy[k];

                g_xx_0_xyzz_xz[k] = -g_xx_0_xzz_xz[k] * ab_y + g_xx_0_xzz_xyz[k];

                g_xx_0_xyzz_yy[k] = -g_xx_0_xzz_yy[k] * ab_y + g_xx_0_xzz_yyy[k];

                g_xx_0_xyzz_yz[k] = -g_xx_0_xzz_yz[k] * ab_y + g_xx_0_xzz_yyz[k];

                g_xx_0_xyzz_zz[k] = -g_xx_0_xzz_zz[k] * ab_y + g_xx_0_xzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzz_xx, g_xx_0_xzz_xxz, g_xx_0_xzz_xy, g_xx_0_xzz_xyz, g_xx_0_xzz_xz, g_xx_0_xzz_xzz, g_xx_0_xzz_yy, g_xx_0_xzz_yyz, g_xx_0_xzz_yz, g_xx_0_xzz_yzz, g_xx_0_xzz_zz, g_xx_0_xzz_zzz, g_xx_0_xzzz_xx, g_xx_0_xzzz_xy, g_xx_0_xzzz_xz, g_xx_0_xzzz_yy, g_xx_0_xzzz_yz, g_xx_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzz_xx[k] = -g_xx_0_xzz_xx[k] * ab_z + g_xx_0_xzz_xxz[k];

                g_xx_0_xzzz_xy[k] = -g_xx_0_xzz_xy[k] * ab_z + g_xx_0_xzz_xyz[k];

                g_xx_0_xzzz_xz[k] = -g_xx_0_xzz_xz[k] * ab_z + g_xx_0_xzz_xzz[k];

                g_xx_0_xzzz_yy[k] = -g_xx_0_xzz_yy[k] * ab_z + g_xx_0_xzz_yyz[k];

                g_xx_0_xzzz_yz[k] = -g_xx_0_xzz_yz[k] * ab_z + g_xx_0_xzz_yzz[k];

                g_xx_0_xzzz_zz[k] = -g_xx_0_xzz_zz[k] * ab_z + g_xx_0_xzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyy_xx, g_xx_0_yyy_xxy, g_xx_0_yyy_xy, g_xx_0_yyy_xyy, g_xx_0_yyy_xyz, g_xx_0_yyy_xz, g_xx_0_yyy_yy, g_xx_0_yyy_yyy, g_xx_0_yyy_yyz, g_xx_0_yyy_yz, g_xx_0_yyy_yzz, g_xx_0_yyy_zz, g_xx_0_yyyy_xx, g_xx_0_yyyy_xy, g_xx_0_yyyy_xz, g_xx_0_yyyy_yy, g_xx_0_yyyy_yz, g_xx_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyy_xx[k] = -g_xx_0_yyy_xx[k] * ab_y + g_xx_0_yyy_xxy[k];

                g_xx_0_yyyy_xy[k] = -g_xx_0_yyy_xy[k] * ab_y + g_xx_0_yyy_xyy[k];

                g_xx_0_yyyy_xz[k] = -g_xx_0_yyy_xz[k] * ab_y + g_xx_0_yyy_xyz[k];

                g_xx_0_yyyy_yy[k] = -g_xx_0_yyy_yy[k] * ab_y + g_xx_0_yyy_yyy[k];

                g_xx_0_yyyy_yz[k] = -g_xx_0_yyy_yz[k] * ab_y + g_xx_0_yyy_yyz[k];

                g_xx_0_yyyy_zz[k] = -g_xx_0_yyy_zz[k] * ab_y + g_xx_0_yyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyz_xx, g_xx_0_yyyz_xy, g_xx_0_yyyz_xz, g_xx_0_yyyz_yy, g_xx_0_yyyz_yz, g_xx_0_yyyz_zz, g_xx_0_yyz_xx, g_xx_0_yyz_xxy, g_xx_0_yyz_xy, g_xx_0_yyz_xyy, g_xx_0_yyz_xyz, g_xx_0_yyz_xz, g_xx_0_yyz_yy, g_xx_0_yyz_yyy, g_xx_0_yyz_yyz, g_xx_0_yyz_yz, g_xx_0_yyz_yzz, g_xx_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyz_xx[k] = -g_xx_0_yyz_xx[k] * ab_y + g_xx_0_yyz_xxy[k];

                g_xx_0_yyyz_xy[k] = -g_xx_0_yyz_xy[k] * ab_y + g_xx_0_yyz_xyy[k];

                g_xx_0_yyyz_xz[k] = -g_xx_0_yyz_xz[k] * ab_y + g_xx_0_yyz_xyz[k];

                g_xx_0_yyyz_yy[k] = -g_xx_0_yyz_yy[k] * ab_y + g_xx_0_yyz_yyy[k];

                g_xx_0_yyyz_yz[k] = -g_xx_0_yyz_yz[k] * ab_y + g_xx_0_yyz_yyz[k];

                g_xx_0_yyyz_zz[k] = -g_xx_0_yyz_zz[k] * ab_y + g_xx_0_yyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzz_xx, g_xx_0_yyzz_xy, g_xx_0_yyzz_xz, g_xx_0_yyzz_yy, g_xx_0_yyzz_yz, g_xx_0_yyzz_zz, g_xx_0_yzz_xx, g_xx_0_yzz_xxy, g_xx_0_yzz_xy, g_xx_0_yzz_xyy, g_xx_0_yzz_xyz, g_xx_0_yzz_xz, g_xx_0_yzz_yy, g_xx_0_yzz_yyy, g_xx_0_yzz_yyz, g_xx_0_yzz_yz, g_xx_0_yzz_yzz, g_xx_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzz_xx[k] = -g_xx_0_yzz_xx[k] * ab_y + g_xx_0_yzz_xxy[k];

                g_xx_0_yyzz_xy[k] = -g_xx_0_yzz_xy[k] * ab_y + g_xx_0_yzz_xyy[k];

                g_xx_0_yyzz_xz[k] = -g_xx_0_yzz_xz[k] * ab_y + g_xx_0_yzz_xyz[k];

                g_xx_0_yyzz_yy[k] = -g_xx_0_yzz_yy[k] * ab_y + g_xx_0_yzz_yyy[k];

                g_xx_0_yyzz_yz[k] = -g_xx_0_yzz_yz[k] * ab_y + g_xx_0_yzz_yyz[k];

                g_xx_0_yyzz_zz[k] = -g_xx_0_yzz_zz[k] * ab_y + g_xx_0_yzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzz_xx, g_xx_0_yzzz_xy, g_xx_0_yzzz_xz, g_xx_0_yzzz_yy, g_xx_0_yzzz_yz, g_xx_0_yzzz_zz, g_xx_0_zzz_xx, g_xx_0_zzz_xxy, g_xx_0_zzz_xy, g_xx_0_zzz_xyy, g_xx_0_zzz_xyz, g_xx_0_zzz_xz, g_xx_0_zzz_yy, g_xx_0_zzz_yyy, g_xx_0_zzz_yyz, g_xx_0_zzz_yz, g_xx_0_zzz_yzz, g_xx_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzz_xx[k] = -g_xx_0_zzz_xx[k] * ab_y + g_xx_0_zzz_xxy[k];

                g_xx_0_yzzz_xy[k] = -g_xx_0_zzz_xy[k] * ab_y + g_xx_0_zzz_xyy[k];

                g_xx_0_yzzz_xz[k] = -g_xx_0_zzz_xz[k] * ab_y + g_xx_0_zzz_xyz[k];

                g_xx_0_yzzz_yy[k] = -g_xx_0_zzz_yy[k] * ab_y + g_xx_0_zzz_yyy[k];

                g_xx_0_yzzz_yz[k] = -g_xx_0_zzz_yz[k] * ab_y + g_xx_0_zzz_yyz[k];

                g_xx_0_yzzz_zz[k] = -g_xx_0_zzz_zz[k] * ab_y + g_xx_0_zzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzz_xx, g_xx_0_zzz_xxz, g_xx_0_zzz_xy, g_xx_0_zzz_xyz, g_xx_0_zzz_xz, g_xx_0_zzz_xzz, g_xx_0_zzz_yy, g_xx_0_zzz_yyz, g_xx_0_zzz_yz, g_xx_0_zzz_yzz, g_xx_0_zzz_zz, g_xx_0_zzz_zzz, g_xx_0_zzzz_xx, g_xx_0_zzzz_xy, g_xx_0_zzzz_xz, g_xx_0_zzzz_yy, g_xx_0_zzzz_yz, g_xx_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzz_xx[k] = -g_xx_0_zzz_xx[k] * ab_z + g_xx_0_zzz_xxz[k];

                g_xx_0_zzzz_xy[k] = -g_xx_0_zzz_xy[k] * ab_z + g_xx_0_zzz_xyz[k];

                g_xx_0_zzzz_xz[k] = -g_xx_0_zzz_xz[k] * ab_z + g_xx_0_zzz_xzz[k];

                g_xx_0_zzzz_yy[k] = -g_xx_0_zzz_yy[k] * ab_z + g_xx_0_zzz_yyz[k];

                g_xx_0_zzzz_yz[k] = -g_xx_0_zzz_yz[k] * ab_z + g_xx_0_zzz_yzz[k];

                g_xx_0_zzzz_zz[k] = -g_xx_0_zzz_zz[k] * ab_z + g_xx_0_zzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_xx, g_xy_0_xxx_xxx, g_xy_0_xxx_xxy, g_xy_0_xxx_xxz, g_xy_0_xxx_xy, g_xy_0_xxx_xyy, g_xy_0_xxx_xyz, g_xy_0_xxx_xz, g_xy_0_xxx_xzz, g_xy_0_xxx_yy, g_xy_0_xxx_yz, g_xy_0_xxx_zz, g_xy_0_xxxx_xx, g_xy_0_xxxx_xy, g_xy_0_xxxx_xz, g_xy_0_xxxx_yy, g_xy_0_xxxx_yz, g_xy_0_xxxx_zz, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxx_xx[k] = -g_y_0_xxx_xx[k] - g_xy_0_xxx_xx[k] * ab_x + g_xy_0_xxx_xxx[k];

                g_xy_0_xxxx_xy[k] = -g_y_0_xxx_xy[k] - g_xy_0_xxx_xy[k] * ab_x + g_xy_0_xxx_xxy[k];

                g_xy_0_xxxx_xz[k] = -g_y_0_xxx_xz[k] - g_xy_0_xxx_xz[k] * ab_x + g_xy_0_xxx_xxz[k];

                g_xy_0_xxxx_yy[k] = -g_y_0_xxx_yy[k] - g_xy_0_xxx_yy[k] * ab_x + g_xy_0_xxx_xyy[k];

                g_xy_0_xxxx_yz[k] = -g_y_0_xxx_yz[k] - g_xy_0_xxx_yz[k] * ab_x + g_xy_0_xxx_xyz[k];

                g_xy_0_xxxx_zz[k] = -g_y_0_xxx_zz[k] - g_xy_0_xxx_zz[k] * ab_x + g_xy_0_xxx_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_xx, g_xy_0_xxxy_xy, g_xy_0_xxxy_xz, g_xy_0_xxxy_yy, g_xy_0_xxxy_yz, g_xy_0_xxxy_zz, g_xy_0_xxy_xx, g_xy_0_xxy_xxx, g_xy_0_xxy_xxy, g_xy_0_xxy_xxz, g_xy_0_xxy_xy, g_xy_0_xxy_xyy, g_xy_0_xxy_xyz, g_xy_0_xxy_xz, g_xy_0_xxy_xzz, g_xy_0_xxy_yy, g_xy_0_xxy_yz, g_xy_0_xxy_zz, g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxy_xx[k] = -g_y_0_xxy_xx[k] - g_xy_0_xxy_xx[k] * ab_x + g_xy_0_xxy_xxx[k];

                g_xy_0_xxxy_xy[k] = -g_y_0_xxy_xy[k] - g_xy_0_xxy_xy[k] * ab_x + g_xy_0_xxy_xxy[k];

                g_xy_0_xxxy_xz[k] = -g_y_0_xxy_xz[k] - g_xy_0_xxy_xz[k] * ab_x + g_xy_0_xxy_xxz[k];

                g_xy_0_xxxy_yy[k] = -g_y_0_xxy_yy[k] - g_xy_0_xxy_yy[k] * ab_x + g_xy_0_xxy_xyy[k];

                g_xy_0_xxxy_yz[k] = -g_y_0_xxy_yz[k] - g_xy_0_xxy_yz[k] * ab_x + g_xy_0_xxy_xyz[k];

                g_xy_0_xxxy_zz[k] = -g_y_0_xxy_zz[k] - g_xy_0_xxy_zz[k] * ab_x + g_xy_0_xxy_xzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_xx, g_xy_0_xxx_xxz, g_xy_0_xxx_xy, g_xy_0_xxx_xyz, g_xy_0_xxx_xz, g_xy_0_xxx_xzz, g_xy_0_xxx_yy, g_xy_0_xxx_yyz, g_xy_0_xxx_yz, g_xy_0_xxx_yzz, g_xy_0_xxx_zz, g_xy_0_xxx_zzz, g_xy_0_xxxz_xx, g_xy_0_xxxz_xy, g_xy_0_xxxz_xz, g_xy_0_xxxz_yy, g_xy_0_xxxz_yz, g_xy_0_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxz_xx[k] = -g_xy_0_xxx_xx[k] * ab_z + g_xy_0_xxx_xxz[k];

                g_xy_0_xxxz_xy[k] = -g_xy_0_xxx_xy[k] * ab_z + g_xy_0_xxx_xyz[k];

                g_xy_0_xxxz_xz[k] = -g_xy_0_xxx_xz[k] * ab_z + g_xy_0_xxx_xzz[k];

                g_xy_0_xxxz_yy[k] = -g_xy_0_xxx_yy[k] * ab_z + g_xy_0_xxx_yyz[k];

                g_xy_0_xxxz_yz[k] = -g_xy_0_xxx_yz[k] * ab_z + g_xy_0_xxx_yzz[k];

                g_xy_0_xxxz_zz[k] = -g_xy_0_xxx_zz[k] * ab_z + g_xy_0_xxx_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_xx, g_xy_0_xxyy_xy, g_xy_0_xxyy_xz, g_xy_0_xxyy_yy, g_xy_0_xxyy_yz, g_xy_0_xxyy_zz, g_xy_0_xyy_xx, g_xy_0_xyy_xxx, g_xy_0_xyy_xxy, g_xy_0_xyy_xxz, g_xy_0_xyy_xy, g_xy_0_xyy_xyy, g_xy_0_xyy_xyz, g_xy_0_xyy_xz, g_xy_0_xyy_xzz, g_xy_0_xyy_yy, g_xy_0_xyy_yz, g_xy_0_xyy_zz, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyy_xx[k] = -g_y_0_xyy_xx[k] - g_xy_0_xyy_xx[k] * ab_x + g_xy_0_xyy_xxx[k];

                g_xy_0_xxyy_xy[k] = -g_y_0_xyy_xy[k] - g_xy_0_xyy_xy[k] * ab_x + g_xy_0_xyy_xxy[k];

                g_xy_0_xxyy_xz[k] = -g_y_0_xyy_xz[k] - g_xy_0_xyy_xz[k] * ab_x + g_xy_0_xyy_xxz[k];

                g_xy_0_xxyy_yy[k] = -g_y_0_xyy_yy[k] - g_xy_0_xyy_yy[k] * ab_x + g_xy_0_xyy_xyy[k];

                g_xy_0_xxyy_yz[k] = -g_y_0_xyy_yz[k] - g_xy_0_xyy_yz[k] * ab_x + g_xy_0_xyy_xyz[k];

                g_xy_0_xxyy_zz[k] = -g_y_0_xyy_zz[k] - g_xy_0_xyy_zz[k] * ab_x + g_xy_0_xyy_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_xx, g_xy_0_xxy_xxz, g_xy_0_xxy_xy, g_xy_0_xxy_xyz, g_xy_0_xxy_xz, g_xy_0_xxy_xzz, g_xy_0_xxy_yy, g_xy_0_xxy_yyz, g_xy_0_xxy_yz, g_xy_0_xxy_yzz, g_xy_0_xxy_zz, g_xy_0_xxy_zzz, g_xy_0_xxyz_xx, g_xy_0_xxyz_xy, g_xy_0_xxyz_xz, g_xy_0_xxyz_yy, g_xy_0_xxyz_yz, g_xy_0_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyz_xx[k] = -g_xy_0_xxy_xx[k] * ab_z + g_xy_0_xxy_xxz[k];

                g_xy_0_xxyz_xy[k] = -g_xy_0_xxy_xy[k] * ab_z + g_xy_0_xxy_xyz[k];

                g_xy_0_xxyz_xz[k] = -g_xy_0_xxy_xz[k] * ab_z + g_xy_0_xxy_xzz[k];

                g_xy_0_xxyz_yy[k] = -g_xy_0_xxy_yy[k] * ab_z + g_xy_0_xxy_yyz[k];

                g_xy_0_xxyz_yz[k] = -g_xy_0_xxy_yz[k] * ab_z + g_xy_0_xxy_yzz[k];

                g_xy_0_xxyz_zz[k] = -g_xy_0_xxy_zz[k] * ab_z + g_xy_0_xxy_zzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxz_xx, g_xy_0_xxz_xxz, g_xy_0_xxz_xy, g_xy_0_xxz_xyz, g_xy_0_xxz_xz, g_xy_0_xxz_xzz, g_xy_0_xxz_yy, g_xy_0_xxz_yyz, g_xy_0_xxz_yz, g_xy_0_xxz_yzz, g_xy_0_xxz_zz, g_xy_0_xxz_zzz, g_xy_0_xxzz_xx, g_xy_0_xxzz_xy, g_xy_0_xxzz_xz, g_xy_0_xxzz_yy, g_xy_0_xxzz_yz, g_xy_0_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzz_xx[k] = -g_xy_0_xxz_xx[k] * ab_z + g_xy_0_xxz_xxz[k];

                g_xy_0_xxzz_xy[k] = -g_xy_0_xxz_xy[k] * ab_z + g_xy_0_xxz_xyz[k];

                g_xy_0_xxzz_xz[k] = -g_xy_0_xxz_xz[k] * ab_z + g_xy_0_xxz_xzz[k];

                g_xy_0_xxzz_yy[k] = -g_xy_0_xxz_yy[k] * ab_z + g_xy_0_xxz_yyz[k];

                g_xy_0_xxzz_yz[k] = -g_xy_0_xxz_yz[k] * ab_z + g_xy_0_xxz_yzz[k];

                g_xy_0_xxzz_zz[k] = -g_xy_0_xxz_zz[k] * ab_z + g_xy_0_xxz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_xx, g_xy_0_xyyy_xy, g_xy_0_xyyy_xz, g_xy_0_xyyy_yy, g_xy_0_xyyy_yz, g_xy_0_xyyy_zz, g_xy_0_yyy_xx, g_xy_0_yyy_xxx, g_xy_0_yyy_xxy, g_xy_0_yyy_xxz, g_xy_0_yyy_xy, g_xy_0_yyy_xyy, g_xy_0_yyy_xyz, g_xy_0_yyy_xz, g_xy_0_yyy_xzz, g_xy_0_yyy_yy, g_xy_0_yyy_yz, g_xy_0_yyy_zz, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyy_xx[k] = -g_y_0_yyy_xx[k] - g_xy_0_yyy_xx[k] * ab_x + g_xy_0_yyy_xxx[k];

                g_xy_0_xyyy_xy[k] = -g_y_0_yyy_xy[k] - g_xy_0_yyy_xy[k] * ab_x + g_xy_0_yyy_xxy[k];

                g_xy_0_xyyy_xz[k] = -g_y_0_yyy_xz[k] - g_xy_0_yyy_xz[k] * ab_x + g_xy_0_yyy_xxz[k];

                g_xy_0_xyyy_yy[k] = -g_y_0_yyy_yy[k] - g_xy_0_yyy_yy[k] * ab_x + g_xy_0_yyy_xyy[k];

                g_xy_0_xyyy_yz[k] = -g_y_0_yyy_yz[k] - g_xy_0_yyy_yz[k] * ab_x + g_xy_0_yyy_xyz[k];

                g_xy_0_xyyy_zz[k] = -g_y_0_yyy_zz[k] - g_xy_0_yyy_zz[k] * ab_x + g_xy_0_yyy_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_xx, g_xy_0_xyy_xxz, g_xy_0_xyy_xy, g_xy_0_xyy_xyz, g_xy_0_xyy_xz, g_xy_0_xyy_xzz, g_xy_0_xyy_yy, g_xy_0_xyy_yyz, g_xy_0_xyy_yz, g_xy_0_xyy_yzz, g_xy_0_xyy_zz, g_xy_0_xyy_zzz, g_xy_0_xyyz_xx, g_xy_0_xyyz_xy, g_xy_0_xyyz_xz, g_xy_0_xyyz_yy, g_xy_0_xyyz_yz, g_xy_0_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyz_xx[k] = -g_xy_0_xyy_xx[k] * ab_z + g_xy_0_xyy_xxz[k];

                g_xy_0_xyyz_xy[k] = -g_xy_0_xyy_xy[k] * ab_z + g_xy_0_xyy_xyz[k];

                g_xy_0_xyyz_xz[k] = -g_xy_0_xyy_xz[k] * ab_z + g_xy_0_xyy_xzz[k];

                g_xy_0_xyyz_yy[k] = -g_xy_0_xyy_yy[k] * ab_z + g_xy_0_xyy_yyz[k];

                g_xy_0_xyyz_yz[k] = -g_xy_0_xyy_yz[k] * ab_z + g_xy_0_xyy_yzz[k];

                g_xy_0_xyyz_zz[k] = -g_xy_0_xyy_zz[k] * ab_z + g_xy_0_xyy_zzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyz_xx, g_xy_0_xyz_xxz, g_xy_0_xyz_xy, g_xy_0_xyz_xyz, g_xy_0_xyz_xz, g_xy_0_xyz_xzz, g_xy_0_xyz_yy, g_xy_0_xyz_yyz, g_xy_0_xyz_yz, g_xy_0_xyz_yzz, g_xy_0_xyz_zz, g_xy_0_xyz_zzz, g_xy_0_xyzz_xx, g_xy_0_xyzz_xy, g_xy_0_xyzz_xz, g_xy_0_xyzz_yy, g_xy_0_xyzz_yz, g_xy_0_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzz_xx[k] = -g_xy_0_xyz_xx[k] * ab_z + g_xy_0_xyz_xxz[k];

                g_xy_0_xyzz_xy[k] = -g_xy_0_xyz_xy[k] * ab_z + g_xy_0_xyz_xyz[k];

                g_xy_0_xyzz_xz[k] = -g_xy_0_xyz_xz[k] * ab_z + g_xy_0_xyz_xzz[k];

                g_xy_0_xyzz_yy[k] = -g_xy_0_xyz_yy[k] * ab_z + g_xy_0_xyz_yyz[k];

                g_xy_0_xyzz_yz[k] = -g_xy_0_xyz_yz[k] * ab_z + g_xy_0_xyz_yzz[k];

                g_xy_0_xyzz_zz[k] = -g_xy_0_xyz_zz[k] * ab_z + g_xy_0_xyz_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzz_xx, g_xy_0_xzz_xxz, g_xy_0_xzz_xy, g_xy_0_xzz_xyz, g_xy_0_xzz_xz, g_xy_0_xzz_xzz, g_xy_0_xzz_yy, g_xy_0_xzz_yyz, g_xy_0_xzz_yz, g_xy_0_xzz_yzz, g_xy_0_xzz_zz, g_xy_0_xzz_zzz, g_xy_0_xzzz_xx, g_xy_0_xzzz_xy, g_xy_0_xzzz_xz, g_xy_0_xzzz_yy, g_xy_0_xzzz_yz, g_xy_0_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzz_xx[k] = -g_xy_0_xzz_xx[k] * ab_z + g_xy_0_xzz_xxz[k];

                g_xy_0_xzzz_xy[k] = -g_xy_0_xzz_xy[k] * ab_z + g_xy_0_xzz_xyz[k];

                g_xy_0_xzzz_xz[k] = -g_xy_0_xzz_xz[k] * ab_z + g_xy_0_xzz_xzz[k];

                g_xy_0_xzzz_yy[k] = -g_xy_0_xzz_yy[k] * ab_z + g_xy_0_xzz_yyz[k];

                g_xy_0_xzzz_yz[k] = -g_xy_0_xzz_yz[k] * ab_z + g_xy_0_xzz_yzz[k];

                g_xy_0_xzzz_zz[k] = -g_xy_0_xzz_zz[k] * ab_z + g_xy_0_xzz_zzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_xx, g_x_0_yyy_xy, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_zz, g_xy_0_yyy_xx, g_xy_0_yyy_xxy, g_xy_0_yyy_xy, g_xy_0_yyy_xyy, g_xy_0_yyy_xyz, g_xy_0_yyy_xz, g_xy_0_yyy_yy, g_xy_0_yyy_yyy, g_xy_0_yyy_yyz, g_xy_0_yyy_yz, g_xy_0_yyy_yzz, g_xy_0_yyy_zz, g_xy_0_yyyy_xx, g_xy_0_yyyy_xy, g_xy_0_yyyy_xz, g_xy_0_yyyy_yy, g_xy_0_yyyy_yz, g_xy_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyy_xx[k] = -g_x_0_yyy_xx[k] - g_xy_0_yyy_xx[k] * ab_y + g_xy_0_yyy_xxy[k];

                g_xy_0_yyyy_xy[k] = -g_x_0_yyy_xy[k] - g_xy_0_yyy_xy[k] * ab_y + g_xy_0_yyy_xyy[k];

                g_xy_0_yyyy_xz[k] = -g_x_0_yyy_xz[k] - g_xy_0_yyy_xz[k] * ab_y + g_xy_0_yyy_xyz[k];

                g_xy_0_yyyy_yy[k] = -g_x_0_yyy_yy[k] - g_xy_0_yyy_yy[k] * ab_y + g_xy_0_yyy_yyy[k];

                g_xy_0_yyyy_yz[k] = -g_x_0_yyy_yz[k] - g_xy_0_yyy_yz[k] * ab_y + g_xy_0_yyy_yyz[k];

                g_xy_0_yyyy_zz[k] = -g_x_0_yyy_zz[k] - g_xy_0_yyy_zz[k] * ab_y + g_xy_0_yyy_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyy_xx, g_xy_0_yyy_xxz, g_xy_0_yyy_xy, g_xy_0_yyy_xyz, g_xy_0_yyy_xz, g_xy_0_yyy_xzz, g_xy_0_yyy_yy, g_xy_0_yyy_yyz, g_xy_0_yyy_yz, g_xy_0_yyy_yzz, g_xy_0_yyy_zz, g_xy_0_yyy_zzz, g_xy_0_yyyz_xx, g_xy_0_yyyz_xy, g_xy_0_yyyz_xz, g_xy_0_yyyz_yy, g_xy_0_yyyz_yz, g_xy_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyz_xx[k] = -g_xy_0_yyy_xx[k] * ab_z + g_xy_0_yyy_xxz[k];

                g_xy_0_yyyz_xy[k] = -g_xy_0_yyy_xy[k] * ab_z + g_xy_0_yyy_xyz[k];

                g_xy_0_yyyz_xz[k] = -g_xy_0_yyy_xz[k] * ab_z + g_xy_0_yyy_xzz[k];

                g_xy_0_yyyz_yy[k] = -g_xy_0_yyy_yy[k] * ab_z + g_xy_0_yyy_yyz[k];

                g_xy_0_yyyz_yz[k] = -g_xy_0_yyy_yz[k] * ab_z + g_xy_0_yyy_yzz[k];

                g_xy_0_yyyz_zz[k] = -g_xy_0_yyy_zz[k] * ab_z + g_xy_0_yyy_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyz_xx, g_xy_0_yyz_xxz, g_xy_0_yyz_xy, g_xy_0_yyz_xyz, g_xy_0_yyz_xz, g_xy_0_yyz_xzz, g_xy_0_yyz_yy, g_xy_0_yyz_yyz, g_xy_0_yyz_yz, g_xy_0_yyz_yzz, g_xy_0_yyz_zz, g_xy_0_yyz_zzz, g_xy_0_yyzz_xx, g_xy_0_yyzz_xy, g_xy_0_yyzz_xz, g_xy_0_yyzz_yy, g_xy_0_yyzz_yz, g_xy_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzz_xx[k] = -g_xy_0_yyz_xx[k] * ab_z + g_xy_0_yyz_xxz[k];

                g_xy_0_yyzz_xy[k] = -g_xy_0_yyz_xy[k] * ab_z + g_xy_0_yyz_xyz[k];

                g_xy_0_yyzz_xz[k] = -g_xy_0_yyz_xz[k] * ab_z + g_xy_0_yyz_xzz[k];

                g_xy_0_yyzz_yy[k] = -g_xy_0_yyz_yy[k] * ab_z + g_xy_0_yyz_yyz[k];

                g_xy_0_yyzz_yz[k] = -g_xy_0_yyz_yz[k] * ab_z + g_xy_0_yyz_yzz[k];

                g_xy_0_yyzz_zz[k] = -g_xy_0_yyz_zz[k] * ab_z + g_xy_0_yyz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzz_xx, g_xy_0_yzz_xxz, g_xy_0_yzz_xy, g_xy_0_yzz_xyz, g_xy_0_yzz_xz, g_xy_0_yzz_xzz, g_xy_0_yzz_yy, g_xy_0_yzz_yyz, g_xy_0_yzz_yz, g_xy_0_yzz_yzz, g_xy_0_yzz_zz, g_xy_0_yzz_zzz, g_xy_0_yzzz_xx, g_xy_0_yzzz_xy, g_xy_0_yzzz_xz, g_xy_0_yzzz_yy, g_xy_0_yzzz_yz, g_xy_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzz_xx[k] = -g_xy_0_yzz_xx[k] * ab_z + g_xy_0_yzz_xxz[k];

                g_xy_0_yzzz_xy[k] = -g_xy_0_yzz_xy[k] * ab_z + g_xy_0_yzz_xyz[k];

                g_xy_0_yzzz_xz[k] = -g_xy_0_yzz_xz[k] * ab_z + g_xy_0_yzz_xzz[k];

                g_xy_0_yzzz_yy[k] = -g_xy_0_yzz_yy[k] * ab_z + g_xy_0_yzz_yyz[k];

                g_xy_0_yzzz_yz[k] = -g_xy_0_yzz_yz[k] * ab_z + g_xy_0_yzz_yzz[k];

                g_xy_0_yzzz_zz[k] = -g_xy_0_yzz_zz[k] * ab_z + g_xy_0_yzz_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzz_xx, g_xy_0_zzz_xxz, g_xy_0_zzz_xy, g_xy_0_zzz_xyz, g_xy_0_zzz_xz, g_xy_0_zzz_xzz, g_xy_0_zzz_yy, g_xy_0_zzz_yyz, g_xy_0_zzz_yz, g_xy_0_zzz_yzz, g_xy_0_zzz_zz, g_xy_0_zzz_zzz, g_xy_0_zzzz_xx, g_xy_0_zzzz_xy, g_xy_0_zzzz_xz, g_xy_0_zzzz_yy, g_xy_0_zzzz_yz, g_xy_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzz_xx[k] = -g_xy_0_zzz_xx[k] * ab_z + g_xy_0_zzz_xxz[k];

                g_xy_0_zzzz_xy[k] = -g_xy_0_zzz_xy[k] * ab_z + g_xy_0_zzz_xyz[k];

                g_xy_0_zzzz_xz[k] = -g_xy_0_zzz_xz[k] * ab_z + g_xy_0_zzz_xzz[k];

                g_xy_0_zzzz_yy[k] = -g_xy_0_zzz_yy[k] * ab_z + g_xy_0_zzz_yyz[k];

                g_xy_0_zzzz_yz[k] = -g_xy_0_zzz_yz[k] * ab_z + g_xy_0_zzz_yzz[k];

                g_xy_0_zzzz_zz[k] = -g_xy_0_zzz_zz[k] * ab_z + g_xy_0_zzz_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_xx, g_xz_0_xxx_xxx, g_xz_0_xxx_xxy, g_xz_0_xxx_xxz, g_xz_0_xxx_xy, g_xz_0_xxx_xyy, g_xz_0_xxx_xyz, g_xz_0_xxx_xz, g_xz_0_xxx_xzz, g_xz_0_xxx_yy, g_xz_0_xxx_yz, g_xz_0_xxx_zz, g_xz_0_xxxx_xx, g_xz_0_xxxx_xy, g_xz_0_xxxx_xz, g_xz_0_xxxx_yy, g_xz_0_xxxx_yz, g_xz_0_xxxx_zz, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxx_xx[k] = -g_z_0_xxx_xx[k] - g_xz_0_xxx_xx[k] * ab_x + g_xz_0_xxx_xxx[k];

                g_xz_0_xxxx_xy[k] = -g_z_0_xxx_xy[k] - g_xz_0_xxx_xy[k] * ab_x + g_xz_0_xxx_xxy[k];

                g_xz_0_xxxx_xz[k] = -g_z_0_xxx_xz[k] - g_xz_0_xxx_xz[k] * ab_x + g_xz_0_xxx_xxz[k];

                g_xz_0_xxxx_yy[k] = -g_z_0_xxx_yy[k] - g_xz_0_xxx_yy[k] * ab_x + g_xz_0_xxx_xyy[k];

                g_xz_0_xxxx_yz[k] = -g_z_0_xxx_yz[k] - g_xz_0_xxx_yz[k] * ab_x + g_xz_0_xxx_xyz[k];

                g_xz_0_xxxx_zz[k] = -g_z_0_xxx_zz[k] - g_xz_0_xxx_zz[k] * ab_x + g_xz_0_xxx_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_xx, g_xz_0_xxx_xxy, g_xz_0_xxx_xy, g_xz_0_xxx_xyy, g_xz_0_xxx_xyz, g_xz_0_xxx_xz, g_xz_0_xxx_yy, g_xz_0_xxx_yyy, g_xz_0_xxx_yyz, g_xz_0_xxx_yz, g_xz_0_xxx_yzz, g_xz_0_xxx_zz, g_xz_0_xxxy_xx, g_xz_0_xxxy_xy, g_xz_0_xxxy_xz, g_xz_0_xxxy_yy, g_xz_0_xxxy_yz, g_xz_0_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxy_xx[k] = -g_xz_0_xxx_xx[k] * ab_y + g_xz_0_xxx_xxy[k];

                g_xz_0_xxxy_xy[k] = -g_xz_0_xxx_xy[k] * ab_y + g_xz_0_xxx_xyy[k];

                g_xz_0_xxxy_xz[k] = -g_xz_0_xxx_xz[k] * ab_y + g_xz_0_xxx_xyz[k];

                g_xz_0_xxxy_yy[k] = -g_xz_0_xxx_yy[k] * ab_y + g_xz_0_xxx_yyy[k];

                g_xz_0_xxxy_yz[k] = -g_xz_0_xxx_yz[k] * ab_y + g_xz_0_xxx_yyz[k];

                g_xz_0_xxxy_zz[k] = -g_xz_0_xxx_zz[k] * ab_y + g_xz_0_xxx_yzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxz_xx, g_xz_0_xxxz_xy, g_xz_0_xxxz_xz, g_xz_0_xxxz_yy, g_xz_0_xxxz_yz, g_xz_0_xxxz_zz, g_xz_0_xxz_xx, g_xz_0_xxz_xxx, g_xz_0_xxz_xxy, g_xz_0_xxz_xxz, g_xz_0_xxz_xy, g_xz_0_xxz_xyy, g_xz_0_xxz_xyz, g_xz_0_xxz_xz, g_xz_0_xxz_xzz, g_xz_0_xxz_yy, g_xz_0_xxz_yz, g_xz_0_xxz_zz, g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxz_xx[k] = -g_z_0_xxz_xx[k] - g_xz_0_xxz_xx[k] * ab_x + g_xz_0_xxz_xxx[k];

                g_xz_0_xxxz_xy[k] = -g_z_0_xxz_xy[k] - g_xz_0_xxz_xy[k] * ab_x + g_xz_0_xxz_xxy[k];

                g_xz_0_xxxz_xz[k] = -g_z_0_xxz_xz[k] - g_xz_0_xxz_xz[k] * ab_x + g_xz_0_xxz_xxz[k];

                g_xz_0_xxxz_yy[k] = -g_z_0_xxz_yy[k] - g_xz_0_xxz_yy[k] * ab_x + g_xz_0_xxz_xyy[k];

                g_xz_0_xxxz_yz[k] = -g_z_0_xxz_yz[k] - g_xz_0_xxz_yz[k] * ab_x + g_xz_0_xxz_xyz[k];

                g_xz_0_xxxz_zz[k] = -g_z_0_xxz_zz[k] - g_xz_0_xxz_zz[k] * ab_x + g_xz_0_xxz_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxy_xx, g_xz_0_xxy_xxy, g_xz_0_xxy_xy, g_xz_0_xxy_xyy, g_xz_0_xxy_xyz, g_xz_0_xxy_xz, g_xz_0_xxy_yy, g_xz_0_xxy_yyy, g_xz_0_xxy_yyz, g_xz_0_xxy_yz, g_xz_0_xxy_yzz, g_xz_0_xxy_zz, g_xz_0_xxyy_xx, g_xz_0_xxyy_xy, g_xz_0_xxyy_xz, g_xz_0_xxyy_yy, g_xz_0_xxyy_yz, g_xz_0_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyy_xx[k] = -g_xz_0_xxy_xx[k] * ab_y + g_xz_0_xxy_xxy[k];

                g_xz_0_xxyy_xy[k] = -g_xz_0_xxy_xy[k] * ab_y + g_xz_0_xxy_xyy[k];

                g_xz_0_xxyy_xz[k] = -g_xz_0_xxy_xz[k] * ab_y + g_xz_0_xxy_xyz[k];

                g_xz_0_xxyy_yy[k] = -g_xz_0_xxy_yy[k] * ab_y + g_xz_0_xxy_yyy[k];

                g_xz_0_xxyy_yz[k] = -g_xz_0_xxy_yz[k] * ab_y + g_xz_0_xxy_yyz[k];

                g_xz_0_xxyy_zz[k] = -g_xz_0_xxy_zz[k] * ab_y + g_xz_0_xxy_yzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyz_xx, g_xz_0_xxyz_xy, g_xz_0_xxyz_xz, g_xz_0_xxyz_yy, g_xz_0_xxyz_yz, g_xz_0_xxyz_zz, g_xz_0_xxz_xx, g_xz_0_xxz_xxy, g_xz_0_xxz_xy, g_xz_0_xxz_xyy, g_xz_0_xxz_xyz, g_xz_0_xxz_xz, g_xz_0_xxz_yy, g_xz_0_xxz_yyy, g_xz_0_xxz_yyz, g_xz_0_xxz_yz, g_xz_0_xxz_yzz, g_xz_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyz_xx[k] = -g_xz_0_xxz_xx[k] * ab_y + g_xz_0_xxz_xxy[k];

                g_xz_0_xxyz_xy[k] = -g_xz_0_xxz_xy[k] * ab_y + g_xz_0_xxz_xyy[k];

                g_xz_0_xxyz_xz[k] = -g_xz_0_xxz_xz[k] * ab_y + g_xz_0_xxz_xyz[k];

                g_xz_0_xxyz_yy[k] = -g_xz_0_xxz_yy[k] * ab_y + g_xz_0_xxz_yyy[k];

                g_xz_0_xxyz_yz[k] = -g_xz_0_xxz_yz[k] * ab_y + g_xz_0_xxz_yyz[k];

                g_xz_0_xxyz_zz[k] = -g_xz_0_xxz_zz[k] * ab_y + g_xz_0_xxz_yzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzz_xx, g_xz_0_xxzz_xy, g_xz_0_xxzz_xz, g_xz_0_xxzz_yy, g_xz_0_xxzz_yz, g_xz_0_xxzz_zz, g_xz_0_xzz_xx, g_xz_0_xzz_xxx, g_xz_0_xzz_xxy, g_xz_0_xzz_xxz, g_xz_0_xzz_xy, g_xz_0_xzz_xyy, g_xz_0_xzz_xyz, g_xz_0_xzz_xz, g_xz_0_xzz_xzz, g_xz_0_xzz_yy, g_xz_0_xzz_yz, g_xz_0_xzz_zz, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzz_xx[k] = -g_z_0_xzz_xx[k] - g_xz_0_xzz_xx[k] * ab_x + g_xz_0_xzz_xxx[k];

                g_xz_0_xxzz_xy[k] = -g_z_0_xzz_xy[k] - g_xz_0_xzz_xy[k] * ab_x + g_xz_0_xzz_xxy[k];

                g_xz_0_xxzz_xz[k] = -g_z_0_xzz_xz[k] - g_xz_0_xzz_xz[k] * ab_x + g_xz_0_xzz_xxz[k];

                g_xz_0_xxzz_yy[k] = -g_z_0_xzz_yy[k] - g_xz_0_xzz_yy[k] * ab_x + g_xz_0_xzz_xyy[k];

                g_xz_0_xxzz_yz[k] = -g_z_0_xzz_yz[k] - g_xz_0_xzz_yz[k] * ab_x + g_xz_0_xzz_xyz[k];

                g_xz_0_xxzz_zz[k] = -g_z_0_xzz_zz[k] - g_xz_0_xzz_zz[k] * ab_x + g_xz_0_xzz_xzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyy_xx, g_xz_0_xyy_xxy, g_xz_0_xyy_xy, g_xz_0_xyy_xyy, g_xz_0_xyy_xyz, g_xz_0_xyy_xz, g_xz_0_xyy_yy, g_xz_0_xyy_yyy, g_xz_0_xyy_yyz, g_xz_0_xyy_yz, g_xz_0_xyy_yzz, g_xz_0_xyy_zz, g_xz_0_xyyy_xx, g_xz_0_xyyy_xy, g_xz_0_xyyy_xz, g_xz_0_xyyy_yy, g_xz_0_xyyy_yz, g_xz_0_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyy_xx[k] = -g_xz_0_xyy_xx[k] * ab_y + g_xz_0_xyy_xxy[k];

                g_xz_0_xyyy_xy[k] = -g_xz_0_xyy_xy[k] * ab_y + g_xz_0_xyy_xyy[k];

                g_xz_0_xyyy_xz[k] = -g_xz_0_xyy_xz[k] * ab_y + g_xz_0_xyy_xyz[k];

                g_xz_0_xyyy_yy[k] = -g_xz_0_xyy_yy[k] * ab_y + g_xz_0_xyy_yyy[k];

                g_xz_0_xyyy_yz[k] = -g_xz_0_xyy_yz[k] * ab_y + g_xz_0_xyy_yyz[k];

                g_xz_0_xyyy_zz[k] = -g_xz_0_xyy_zz[k] * ab_y + g_xz_0_xyy_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyz_xx, g_xz_0_xyyz_xy, g_xz_0_xyyz_xz, g_xz_0_xyyz_yy, g_xz_0_xyyz_yz, g_xz_0_xyyz_zz, g_xz_0_xyz_xx, g_xz_0_xyz_xxy, g_xz_0_xyz_xy, g_xz_0_xyz_xyy, g_xz_0_xyz_xyz, g_xz_0_xyz_xz, g_xz_0_xyz_yy, g_xz_0_xyz_yyy, g_xz_0_xyz_yyz, g_xz_0_xyz_yz, g_xz_0_xyz_yzz, g_xz_0_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyz_xx[k] = -g_xz_0_xyz_xx[k] * ab_y + g_xz_0_xyz_xxy[k];

                g_xz_0_xyyz_xy[k] = -g_xz_0_xyz_xy[k] * ab_y + g_xz_0_xyz_xyy[k];

                g_xz_0_xyyz_xz[k] = -g_xz_0_xyz_xz[k] * ab_y + g_xz_0_xyz_xyz[k];

                g_xz_0_xyyz_yy[k] = -g_xz_0_xyz_yy[k] * ab_y + g_xz_0_xyz_yyy[k];

                g_xz_0_xyyz_yz[k] = -g_xz_0_xyz_yz[k] * ab_y + g_xz_0_xyz_yyz[k];

                g_xz_0_xyyz_zz[k] = -g_xz_0_xyz_zz[k] * ab_y + g_xz_0_xyz_yzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzz_xx, g_xz_0_xyzz_xy, g_xz_0_xyzz_xz, g_xz_0_xyzz_yy, g_xz_0_xyzz_yz, g_xz_0_xyzz_zz, g_xz_0_xzz_xx, g_xz_0_xzz_xxy, g_xz_0_xzz_xy, g_xz_0_xzz_xyy, g_xz_0_xzz_xyz, g_xz_0_xzz_xz, g_xz_0_xzz_yy, g_xz_0_xzz_yyy, g_xz_0_xzz_yyz, g_xz_0_xzz_yz, g_xz_0_xzz_yzz, g_xz_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzz_xx[k] = -g_xz_0_xzz_xx[k] * ab_y + g_xz_0_xzz_xxy[k];

                g_xz_0_xyzz_xy[k] = -g_xz_0_xzz_xy[k] * ab_y + g_xz_0_xzz_xyy[k];

                g_xz_0_xyzz_xz[k] = -g_xz_0_xzz_xz[k] * ab_y + g_xz_0_xzz_xyz[k];

                g_xz_0_xyzz_yy[k] = -g_xz_0_xzz_yy[k] * ab_y + g_xz_0_xzz_yyy[k];

                g_xz_0_xyzz_yz[k] = -g_xz_0_xzz_yz[k] * ab_y + g_xz_0_xzz_yyz[k];

                g_xz_0_xyzz_zz[k] = -g_xz_0_xzz_zz[k] * ab_y + g_xz_0_xzz_yzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzz_xx, g_xz_0_xzzz_xy, g_xz_0_xzzz_xz, g_xz_0_xzzz_yy, g_xz_0_xzzz_yz, g_xz_0_xzzz_zz, g_xz_0_zzz_xx, g_xz_0_zzz_xxx, g_xz_0_zzz_xxy, g_xz_0_zzz_xxz, g_xz_0_zzz_xy, g_xz_0_zzz_xyy, g_xz_0_zzz_xyz, g_xz_0_zzz_xz, g_xz_0_zzz_xzz, g_xz_0_zzz_yy, g_xz_0_zzz_yz, g_xz_0_zzz_zz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzz_xx[k] = -g_z_0_zzz_xx[k] - g_xz_0_zzz_xx[k] * ab_x + g_xz_0_zzz_xxx[k];

                g_xz_0_xzzz_xy[k] = -g_z_0_zzz_xy[k] - g_xz_0_zzz_xy[k] * ab_x + g_xz_0_zzz_xxy[k];

                g_xz_0_xzzz_xz[k] = -g_z_0_zzz_xz[k] - g_xz_0_zzz_xz[k] * ab_x + g_xz_0_zzz_xxz[k];

                g_xz_0_xzzz_yy[k] = -g_z_0_zzz_yy[k] - g_xz_0_zzz_yy[k] * ab_x + g_xz_0_zzz_xyy[k];

                g_xz_0_xzzz_yz[k] = -g_z_0_zzz_yz[k] - g_xz_0_zzz_yz[k] * ab_x + g_xz_0_zzz_xyz[k];

                g_xz_0_xzzz_zz[k] = -g_z_0_zzz_zz[k] - g_xz_0_zzz_zz[k] * ab_x + g_xz_0_zzz_xzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyy_xx, g_xz_0_yyy_xxy, g_xz_0_yyy_xy, g_xz_0_yyy_xyy, g_xz_0_yyy_xyz, g_xz_0_yyy_xz, g_xz_0_yyy_yy, g_xz_0_yyy_yyy, g_xz_0_yyy_yyz, g_xz_0_yyy_yz, g_xz_0_yyy_yzz, g_xz_0_yyy_zz, g_xz_0_yyyy_xx, g_xz_0_yyyy_xy, g_xz_0_yyyy_xz, g_xz_0_yyyy_yy, g_xz_0_yyyy_yz, g_xz_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyy_xx[k] = -g_xz_0_yyy_xx[k] * ab_y + g_xz_0_yyy_xxy[k];

                g_xz_0_yyyy_xy[k] = -g_xz_0_yyy_xy[k] * ab_y + g_xz_0_yyy_xyy[k];

                g_xz_0_yyyy_xz[k] = -g_xz_0_yyy_xz[k] * ab_y + g_xz_0_yyy_xyz[k];

                g_xz_0_yyyy_yy[k] = -g_xz_0_yyy_yy[k] * ab_y + g_xz_0_yyy_yyy[k];

                g_xz_0_yyyy_yz[k] = -g_xz_0_yyy_yz[k] * ab_y + g_xz_0_yyy_yyz[k];

                g_xz_0_yyyy_zz[k] = -g_xz_0_yyy_zz[k] * ab_y + g_xz_0_yyy_yzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyz_xx, g_xz_0_yyyz_xy, g_xz_0_yyyz_xz, g_xz_0_yyyz_yy, g_xz_0_yyyz_yz, g_xz_0_yyyz_zz, g_xz_0_yyz_xx, g_xz_0_yyz_xxy, g_xz_0_yyz_xy, g_xz_0_yyz_xyy, g_xz_0_yyz_xyz, g_xz_0_yyz_xz, g_xz_0_yyz_yy, g_xz_0_yyz_yyy, g_xz_0_yyz_yyz, g_xz_0_yyz_yz, g_xz_0_yyz_yzz, g_xz_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyz_xx[k] = -g_xz_0_yyz_xx[k] * ab_y + g_xz_0_yyz_xxy[k];

                g_xz_0_yyyz_xy[k] = -g_xz_0_yyz_xy[k] * ab_y + g_xz_0_yyz_xyy[k];

                g_xz_0_yyyz_xz[k] = -g_xz_0_yyz_xz[k] * ab_y + g_xz_0_yyz_xyz[k];

                g_xz_0_yyyz_yy[k] = -g_xz_0_yyz_yy[k] * ab_y + g_xz_0_yyz_yyy[k];

                g_xz_0_yyyz_yz[k] = -g_xz_0_yyz_yz[k] * ab_y + g_xz_0_yyz_yyz[k];

                g_xz_0_yyyz_zz[k] = -g_xz_0_yyz_zz[k] * ab_y + g_xz_0_yyz_yzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzz_xx, g_xz_0_yyzz_xy, g_xz_0_yyzz_xz, g_xz_0_yyzz_yy, g_xz_0_yyzz_yz, g_xz_0_yyzz_zz, g_xz_0_yzz_xx, g_xz_0_yzz_xxy, g_xz_0_yzz_xy, g_xz_0_yzz_xyy, g_xz_0_yzz_xyz, g_xz_0_yzz_xz, g_xz_0_yzz_yy, g_xz_0_yzz_yyy, g_xz_0_yzz_yyz, g_xz_0_yzz_yz, g_xz_0_yzz_yzz, g_xz_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzz_xx[k] = -g_xz_0_yzz_xx[k] * ab_y + g_xz_0_yzz_xxy[k];

                g_xz_0_yyzz_xy[k] = -g_xz_0_yzz_xy[k] * ab_y + g_xz_0_yzz_xyy[k];

                g_xz_0_yyzz_xz[k] = -g_xz_0_yzz_xz[k] * ab_y + g_xz_0_yzz_xyz[k];

                g_xz_0_yyzz_yy[k] = -g_xz_0_yzz_yy[k] * ab_y + g_xz_0_yzz_yyy[k];

                g_xz_0_yyzz_yz[k] = -g_xz_0_yzz_yz[k] * ab_y + g_xz_0_yzz_yyz[k];

                g_xz_0_yyzz_zz[k] = -g_xz_0_yzz_zz[k] * ab_y + g_xz_0_yzz_yzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzz_xx, g_xz_0_yzzz_xy, g_xz_0_yzzz_xz, g_xz_0_yzzz_yy, g_xz_0_yzzz_yz, g_xz_0_yzzz_zz, g_xz_0_zzz_xx, g_xz_0_zzz_xxy, g_xz_0_zzz_xy, g_xz_0_zzz_xyy, g_xz_0_zzz_xyz, g_xz_0_zzz_xz, g_xz_0_zzz_yy, g_xz_0_zzz_yyy, g_xz_0_zzz_yyz, g_xz_0_zzz_yz, g_xz_0_zzz_yzz, g_xz_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzz_xx[k] = -g_xz_0_zzz_xx[k] * ab_y + g_xz_0_zzz_xxy[k];

                g_xz_0_yzzz_xy[k] = -g_xz_0_zzz_xy[k] * ab_y + g_xz_0_zzz_xyy[k];

                g_xz_0_yzzz_xz[k] = -g_xz_0_zzz_xz[k] * ab_y + g_xz_0_zzz_xyz[k];

                g_xz_0_yzzz_yy[k] = -g_xz_0_zzz_yy[k] * ab_y + g_xz_0_zzz_yyy[k];

                g_xz_0_yzzz_yz[k] = -g_xz_0_zzz_yz[k] * ab_y + g_xz_0_zzz_yyz[k];

                g_xz_0_yzzz_zz[k] = -g_xz_0_zzz_zz[k] * ab_y + g_xz_0_zzz_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_xx, g_x_0_zzz_xy, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_zz, g_xz_0_zzz_xx, g_xz_0_zzz_xxz, g_xz_0_zzz_xy, g_xz_0_zzz_xyz, g_xz_0_zzz_xz, g_xz_0_zzz_xzz, g_xz_0_zzz_yy, g_xz_0_zzz_yyz, g_xz_0_zzz_yz, g_xz_0_zzz_yzz, g_xz_0_zzz_zz, g_xz_0_zzz_zzz, g_xz_0_zzzz_xx, g_xz_0_zzzz_xy, g_xz_0_zzzz_xz, g_xz_0_zzzz_yy, g_xz_0_zzzz_yz, g_xz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzz_xx[k] = -g_x_0_zzz_xx[k] - g_xz_0_zzz_xx[k] * ab_z + g_xz_0_zzz_xxz[k];

                g_xz_0_zzzz_xy[k] = -g_x_0_zzz_xy[k] - g_xz_0_zzz_xy[k] * ab_z + g_xz_0_zzz_xyz[k];

                g_xz_0_zzzz_xz[k] = -g_x_0_zzz_xz[k] - g_xz_0_zzz_xz[k] * ab_z + g_xz_0_zzz_xzz[k];

                g_xz_0_zzzz_yy[k] = -g_x_0_zzz_yy[k] - g_xz_0_zzz_yy[k] * ab_z + g_xz_0_zzz_yyz[k];

                g_xz_0_zzzz_yz[k] = -g_x_0_zzz_yz[k] - g_xz_0_zzz_yz[k] * ab_z + g_xz_0_zzz_yzz[k];

                g_xz_0_zzzz_zz[k] = -g_x_0_zzz_zz[k] - g_xz_0_zzz_zz[k] * ab_z + g_xz_0_zzz_zzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxx_xx, g_yy_0_xxx_xxx, g_yy_0_xxx_xxy, g_yy_0_xxx_xxz, g_yy_0_xxx_xy, g_yy_0_xxx_xyy, g_yy_0_xxx_xyz, g_yy_0_xxx_xz, g_yy_0_xxx_xzz, g_yy_0_xxx_yy, g_yy_0_xxx_yz, g_yy_0_xxx_zz, g_yy_0_xxxx_xx, g_yy_0_xxxx_xy, g_yy_0_xxxx_xz, g_yy_0_xxxx_yy, g_yy_0_xxxx_yz, g_yy_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxx_xx[k] = -g_yy_0_xxx_xx[k] * ab_x + g_yy_0_xxx_xxx[k];

                g_yy_0_xxxx_xy[k] = -g_yy_0_xxx_xy[k] * ab_x + g_yy_0_xxx_xxy[k];

                g_yy_0_xxxx_xz[k] = -g_yy_0_xxx_xz[k] * ab_x + g_yy_0_xxx_xxz[k];

                g_yy_0_xxxx_yy[k] = -g_yy_0_xxx_yy[k] * ab_x + g_yy_0_xxx_xyy[k];

                g_yy_0_xxxx_yz[k] = -g_yy_0_xxx_yz[k] * ab_x + g_yy_0_xxx_xyz[k];

                g_yy_0_xxxx_zz[k] = -g_yy_0_xxx_zz[k] * ab_x + g_yy_0_xxx_xzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxy_xx, g_yy_0_xxxy_xy, g_yy_0_xxxy_xz, g_yy_0_xxxy_yy, g_yy_0_xxxy_yz, g_yy_0_xxxy_zz, g_yy_0_xxy_xx, g_yy_0_xxy_xxx, g_yy_0_xxy_xxy, g_yy_0_xxy_xxz, g_yy_0_xxy_xy, g_yy_0_xxy_xyy, g_yy_0_xxy_xyz, g_yy_0_xxy_xz, g_yy_0_xxy_xzz, g_yy_0_xxy_yy, g_yy_0_xxy_yz, g_yy_0_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxy_xx[k] = -g_yy_0_xxy_xx[k] * ab_x + g_yy_0_xxy_xxx[k];

                g_yy_0_xxxy_xy[k] = -g_yy_0_xxy_xy[k] * ab_x + g_yy_0_xxy_xxy[k];

                g_yy_0_xxxy_xz[k] = -g_yy_0_xxy_xz[k] * ab_x + g_yy_0_xxy_xxz[k];

                g_yy_0_xxxy_yy[k] = -g_yy_0_xxy_yy[k] * ab_x + g_yy_0_xxy_xyy[k];

                g_yy_0_xxxy_yz[k] = -g_yy_0_xxy_yz[k] * ab_x + g_yy_0_xxy_xyz[k];

                g_yy_0_xxxy_zz[k] = -g_yy_0_xxy_zz[k] * ab_x + g_yy_0_xxy_xzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 284 * ccomps * dcomps);

            auto g_yy_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxz_xx, g_yy_0_xxxz_xy, g_yy_0_xxxz_xz, g_yy_0_xxxz_yy, g_yy_0_xxxz_yz, g_yy_0_xxxz_zz, g_yy_0_xxz_xx, g_yy_0_xxz_xxx, g_yy_0_xxz_xxy, g_yy_0_xxz_xxz, g_yy_0_xxz_xy, g_yy_0_xxz_xyy, g_yy_0_xxz_xyz, g_yy_0_xxz_xz, g_yy_0_xxz_xzz, g_yy_0_xxz_yy, g_yy_0_xxz_yz, g_yy_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxz_xx[k] = -g_yy_0_xxz_xx[k] * ab_x + g_yy_0_xxz_xxx[k];

                g_yy_0_xxxz_xy[k] = -g_yy_0_xxz_xy[k] * ab_x + g_yy_0_xxz_xxy[k];

                g_yy_0_xxxz_xz[k] = -g_yy_0_xxz_xz[k] * ab_x + g_yy_0_xxz_xxz[k];

                g_yy_0_xxxz_yy[k] = -g_yy_0_xxz_yy[k] * ab_x + g_yy_0_xxz_xyy[k];

                g_yy_0_xxxz_yz[k] = -g_yy_0_xxz_yz[k] * ab_x + g_yy_0_xxz_xyz[k];

                g_yy_0_xxxz_zz[k] = -g_yy_0_xxz_zz[k] * ab_x + g_yy_0_xxz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyy_xx, g_yy_0_xxyy_xy, g_yy_0_xxyy_xz, g_yy_0_xxyy_yy, g_yy_0_xxyy_yz, g_yy_0_xxyy_zz, g_yy_0_xyy_xx, g_yy_0_xyy_xxx, g_yy_0_xyy_xxy, g_yy_0_xyy_xxz, g_yy_0_xyy_xy, g_yy_0_xyy_xyy, g_yy_0_xyy_xyz, g_yy_0_xyy_xz, g_yy_0_xyy_xzz, g_yy_0_xyy_yy, g_yy_0_xyy_yz, g_yy_0_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyy_xx[k] = -g_yy_0_xyy_xx[k] * ab_x + g_yy_0_xyy_xxx[k];

                g_yy_0_xxyy_xy[k] = -g_yy_0_xyy_xy[k] * ab_x + g_yy_0_xyy_xxy[k];

                g_yy_0_xxyy_xz[k] = -g_yy_0_xyy_xz[k] * ab_x + g_yy_0_xyy_xxz[k];

                g_yy_0_xxyy_yy[k] = -g_yy_0_xyy_yy[k] * ab_x + g_yy_0_xyy_xyy[k];

                g_yy_0_xxyy_yz[k] = -g_yy_0_xyy_yz[k] * ab_x + g_yy_0_xyy_xyz[k];

                g_yy_0_xxyy_zz[k] = -g_yy_0_xyy_zz[k] * ab_x + g_yy_0_xyy_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyz_xx, g_yy_0_xxyz_xy, g_yy_0_xxyz_xz, g_yy_0_xxyz_yy, g_yy_0_xxyz_yz, g_yy_0_xxyz_zz, g_yy_0_xyz_xx, g_yy_0_xyz_xxx, g_yy_0_xyz_xxy, g_yy_0_xyz_xxz, g_yy_0_xyz_xy, g_yy_0_xyz_xyy, g_yy_0_xyz_xyz, g_yy_0_xyz_xz, g_yy_0_xyz_xzz, g_yy_0_xyz_yy, g_yy_0_xyz_yz, g_yy_0_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyz_xx[k] = -g_yy_0_xyz_xx[k] * ab_x + g_yy_0_xyz_xxx[k];

                g_yy_0_xxyz_xy[k] = -g_yy_0_xyz_xy[k] * ab_x + g_yy_0_xyz_xxy[k];

                g_yy_0_xxyz_xz[k] = -g_yy_0_xyz_xz[k] * ab_x + g_yy_0_xyz_xxz[k];

                g_yy_0_xxyz_yy[k] = -g_yy_0_xyz_yy[k] * ab_x + g_yy_0_xyz_xyy[k];

                g_yy_0_xxyz_yz[k] = -g_yy_0_xyz_yz[k] * ab_x + g_yy_0_xyz_xyz[k];

                g_yy_0_xxyz_zz[k] = -g_yy_0_xyz_zz[k] * ab_x + g_yy_0_xyz_xzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzz_xx, g_yy_0_xxzz_xy, g_yy_0_xxzz_xz, g_yy_0_xxzz_yy, g_yy_0_xxzz_yz, g_yy_0_xxzz_zz, g_yy_0_xzz_xx, g_yy_0_xzz_xxx, g_yy_0_xzz_xxy, g_yy_0_xzz_xxz, g_yy_0_xzz_xy, g_yy_0_xzz_xyy, g_yy_0_xzz_xyz, g_yy_0_xzz_xz, g_yy_0_xzz_xzz, g_yy_0_xzz_yy, g_yy_0_xzz_yz, g_yy_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzz_xx[k] = -g_yy_0_xzz_xx[k] * ab_x + g_yy_0_xzz_xxx[k];

                g_yy_0_xxzz_xy[k] = -g_yy_0_xzz_xy[k] * ab_x + g_yy_0_xzz_xxy[k];

                g_yy_0_xxzz_xz[k] = -g_yy_0_xzz_xz[k] * ab_x + g_yy_0_xzz_xxz[k];

                g_yy_0_xxzz_yy[k] = -g_yy_0_xzz_yy[k] * ab_x + g_yy_0_xzz_xyy[k];

                g_yy_0_xxzz_yz[k] = -g_yy_0_xzz_yz[k] * ab_x + g_yy_0_xzz_xyz[k];

                g_yy_0_xxzz_zz[k] = -g_yy_0_xzz_zz[k] * ab_x + g_yy_0_xzz_xzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyy_xx, g_yy_0_xyyy_xy, g_yy_0_xyyy_xz, g_yy_0_xyyy_yy, g_yy_0_xyyy_yz, g_yy_0_xyyy_zz, g_yy_0_yyy_xx, g_yy_0_yyy_xxx, g_yy_0_yyy_xxy, g_yy_0_yyy_xxz, g_yy_0_yyy_xy, g_yy_0_yyy_xyy, g_yy_0_yyy_xyz, g_yy_0_yyy_xz, g_yy_0_yyy_xzz, g_yy_0_yyy_yy, g_yy_0_yyy_yz, g_yy_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyy_xx[k] = -g_yy_0_yyy_xx[k] * ab_x + g_yy_0_yyy_xxx[k];

                g_yy_0_xyyy_xy[k] = -g_yy_0_yyy_xy[k] * ab_x + g_yy_0_yyy_xxy[k];

                g_yy_0_xyyy_xz[k] = -g_yy_0_yyy_xz[k] * ab_x + g_yy_0_yyy_xxz[k];

                g_yy_0_xyyy_yy[k] = -g_yy_0_yyy_yy[k] * ab_x + g_yy_0_yyy_xyy[k];

                g_yy_0_xyyy_yz[k] = -g_yy_0_yyy_yz[k] * ab_x + g_yy_0_yyy_xyz[k];

                g_yy_0_xyyy_zz[k] = -g_yy_0_yyy_zz[k] * ab_x + g_yy_0_yyy_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyz_xx, g_yy_0_xyyz_xy, g_yy_0_xyyz_xz, g_yy_0_xyyz_yy, g_yy_0_xyyz_yz, g_yy_0_xyyz_zz, g_yy_0_yyz_xx, g_yy_0_yyz_xxx, g_yy_0_yyz_xxy, g_yy_0_yyz_xxz, g_yy_0_yyz_xy, g_yy_0_yyz_xyy, g_yy_0_yyz_xyz, g_yy_0_yyz_xz, g_yy_0_yyz_xzz, g_yy_0_yyz_yy, g_yy_0_yyz_yz, g_yy_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyz_xx[k] = -g_yy_0_yyz_xx[k] * ab_x + g_yy_0_yyz_xxx[k];

                g_yy_0_xyyz_xy[k] = -g_yy_0_yyz_xy[k] * ab_x + g_yy_0_yyz_xxy[k];

                g_yy_0_xyyz_xz[k] = -g_yy_0_yyz_xz[k] * ab_x + g_yy_0_yyz_xxz[k];

                g_yy_0_xyyz_yy[k] = -g_yy_0_yyz_yy[k] * ab_x + g_yy_0_yyz_xyy[k];

                g_yy_0_xyyz_yz[k] = -g_yy_0_yyz_yz[k] * ab_x + g_yy_0_yyz_xyz[k];

                g_yy_0_xyyz_zz[k] = -g_yy_0_yyz_zz[k] * ab_x + g_yy_0_yyz_xzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzz_xx, g_yy_0_xyzz_xy, g_yy_0_xyzz_xz, g_yy_0_xyzz_yy, g_yy_0_xyzz_yz, g_yy_0_xyzz_zz, g_yy_0_yzz_xx, g_yy_0_yzz_xxx, g_yy_0_yzz_xxy, g_yy_0_yzz_xxz, g_yy_0_yzz_xy, g_yy_0_yzz_xyy, g_yy_0_yzz_xyz, g_yy_0_yzz_xz, g_yy_0_yzz_xzz, g_yy_0_yzz_yy, g_yy_0_yzz_yz, g_yy_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzz_xx[k] = -g_yy_0_yzz_xx[k] * ab_x + g_yy_0_yzz_xxx[k];

                g_yy_0_xyzz_xy[k] = -g_yy_0_yzz_xy[k] * ab_x + g_yy_0_yzz_xxy[k];

                g_yy_0_xyzz_xz[k] = -g_yy_0_yzz_xz[k] * ab_x + g_yy_0_yzz_xxz[k];

                g_yy_0_xyzz_yy[k] = -g_yy_0_yzz_yy[k] * ab_x + g_yy_0_yzz_xyy[k];

                g_yy_0_xyzz_yz[k] = -g_yy_0_yzz_yz[k] * ab_x + g_yy_0_yzz_xyz[k];

                g_yy_0_xyzz_zz[k] = -g_yy_0_yzz_zz[k] * ab_x + g_yy_0_yzz_xzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzz_xx, g_yy_0_xzzz_xy, g_yy_0_xzzz_xz, g_yy_0_xzzz_yy, g_yy_0_xzzz_yz, g_yy_0_xzzz_zz, g_yy_0_zzz_xx, g_yy_0_zzz_xxx, g_yy_0_zzz_xxy, g_yy_0_zzz_xxz, g_yy_0_zzz_xy, g_yy_0_zzz_xyy, g_yy_0_zzz_xyz, g_yy_0_zzz_xz, g_yy_0_zzz_xzz, g_yy_0_zzz_yy, g_yy_0_zzz_yz, g_yy_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzz_xx[k] = -g_yy_0_zzz_xx[k] * ab_x + g_yy_0_zzz_xxx[k];

                g_yy_0_xzzz_xy[k] = -g_yy_0_zzz_xy[k] * ab_x + g_yy_0_zzz_xxy[k];

                g_yy_0_xzzz_xz[k] = -g_yy_0_zzz_xz[k] * ab_x + g_yy_0_zzz_xxz[k];

                g_yy_0_xzzz_yy[k] = -g_yy_0_zzz_yy[k] * ab_x + g_yy_0_zzz_xyy[k];

                g_yy_0_xzzz_yz[k] = -g_yy_0_zzz_yz[k] * ab_x + g_yy_0_zzz_xyz[k];

                g_yy_0_xzzz_zz[k] = -g_yy_0_zzz_zz[k] * ab_x + g_yy_0_zzz_xzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, g_yy_0_yyy_xx, g_yy_0_yyy_xxy, g_yy_0_yyy_xy, g_yy_0_yyy_xyy, g_yy_0_yyy_xyz, g_yy_0_yyy_xz, g_yy_0_yyy_yy, g_yy_0_yyy_yyy, g_yy_0_yyy_yyz, g_yy_0_yyy_yz, g_yy_0_yyy_yzz, g_yy_0_yyy_zz, g_yy_0_yyyy_xx, g_yy_0_yyyy_xy, g_yy_0_yyyy_xz, g_yy_0_yyyy_yy, g_yy_0_yyyy_yz, g_yy_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyy_xx[k] = -2.0 * g_y_0_yyy_xx[k] - g_yy_0_yyy_xx[k] * ab_y + g_yy_0_yyy_xxy[k];

                g_yy_0_yyyy_xy[k] = -2.0 * g_y_0_yyy_xy[k] - g_yy_0_yyy_xy[k] * ab_y + g_yy_0_yyy_xyy[k];

                g_yy_0_yyyy_xz[k] = -2.0 * g_y_0_yyy_xz[k] - g_yy_0_yyy_xz[k] * ab_y + g_yy_0_yyy_xyz[k];

                g_yy_0_yyyy_yy[k] = -2.0 * g_y_0_yyy_yy[k] - g_yy_0_yyy_yy[k] * ab_y + g_yy_0_yyy_yyy[k];

                g_yy_0_yyyy_yz[k] = -2.0 * g_y_0_yyy_yz[k] - g_yy_0_yyy_yz[k] * ab_y + g_yy_0_yyy_yyz[k];

                g_yy_0_yyyy_zz[k] = -2.0 * g_y_0_yyy_zz[k] - g_yy_0_yyy_zz[k] * ab_y + g_yy_0_yyy_yzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyy_xx, g_yy_0_yyy_xxz, g_yy_0_yyy_xy, g_yy_0_yyy_xyz, g_yy_0_yyy_xz, g_yy_0_yyy_xzz, g_yy_0_yyy_yy, g_yy_0_yyy_yyz, g_yy_0_yyy_yz, g_yy_0_yyy_yzz, g_yy_0_yyy_zz, g_yy_0_yyy_zzz, g_yy_0_yyyz_xx, g_yy_0_yyyz_xy, g_yy_0_yyyz_xz, g_yy_0_yyyz_yy, g_yy_0_yyyz_yz, g_yy_0_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyz_xx[k] = -g_yy_0_yyy_xx[k] * ab_z + g_yy_0_yyy_xxz[k];

                g_yy_0_yyyz_xy[k] = -g_yy_0_yyy_xy[k] * ab_z + g_yy_0_yyy_xyz[k];

                g_yy_0_yyyz_xz[k] = -g_yy_0_yyy_xz[k] * ab_z + g_yy_0_yyy_xzz[k];

                g_yy_0_yyyz_yy[k] = -g_yy_0_yyy_yy[k] * ab_z + g_yy_0_yyy_yyz[k];

                g_yy_0_yyyz_yz[k] = -g_yy_0_yyy_yz[k] * ab_z + g_yy_0_yyy_yzz[k];

                g_yy_0_yyyz_zz[k] = -g_yy_0_yyy_zz[k] * ab_z + g_yy_0_yyy_zzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyz_xx, g_yy_0_yyz_xxz, g_yy_0_yyz_xy, g_yy_0_yyz_xyz, g_yy_0_yyz_xz, g_yy_0_yyz_xzz, g_yy_0_yyz_yy, g_yy_0_yyz_yyz, g_yy_0_yyz_yz, g_yy_0_yyz_yzz, g_yy_0_yyz_zz, g_yy_0_yyz_zzz, g_yy_0_yyzz_xx, g_yy_0_yyzz_xy, g_yy_0_yyzz_xz, g_yy_0_yyzz_yy, g_yy_0_yyzz_yz, g_yy_0_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzz_xx[k] = -g_yy_0_yyz_xx[k] * ab_z + g_yy_0_yyz_xxz[k];

                g_yy_0_yyzz_xy[k] = -g_yy_0_yyz_xy[k] * ab_z + g_yy_0_yyz_xyz[k];

                g_yy_0_yyzz_xz[k] = -g_yy_0_yyz_xz[k] * ab_z + g_yy_0_yyz_xzz[k];

                g_yy_0_yyzz_yy[k] = -g_yy_0_yyz_yy[k] * ab_z + g_yy_0_yyz_yyz[k];

                g_yy_0_yyzz_yz[k] = -g_yy_0_yyz_yz[k] * ab_z + g_yy_0_yyz_yzz[k];

                g_yy_0_yyzz_zz[k] = -g_yy_0_yyz_zz[k] * ab_z + g_yy_0_yyz_zzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzz_xx, g_yy_0_yzz_xxz, g_yy_0_yzz_xy, g_yy_0_yzz_xyz, g_yy_0_yzz_xz, g_yy_0_yzz_xzz, g_yy_0_yzz_yy, g_yy_0_yzz_yyz, g_yy_0_yzz_yz, g_yy_0_yzz_yzz, g_yy_0_yzz_zz, g_yy_0_yzz_zzz, g_yy_0_yzzz_xx, g_yy_0_yzzz_xy, g_yy_0_yzzz_xz, g_yy_0_yzzz_yy, g_yy_0_yzzz_yz, g_yy_0_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzz_xx[k] = -g_yy_0_yzz_xx[k] * ab_z + g_yy_0_yzz_xxz[k];

                g_yy_0_yzzz_xy[k] = -g_yy_0_yzz_xy[k] * ab_z + g_yy_0_yzz_xyz[k];

                g_yy_0_yzzz_xz[k] = -g_yy_0_yzz_xz[k] * ab_z + g_yy_0_yzz_xzz[k];

                g_yy_0_yzzz_yy[k] = -g_yy_0_yzz_yy[k] * ab_z + g_yy_0_yzz_yyz[k];

                g_yy_0_yzzz_yz[k] = -g_yy_0_yzz_yz[k] * ab_z + g_yy_0_yzz_yzz[k];

                g_yy_0_yzzz_zz[k] = -g_yy_0_yzz_zz[k] * ab_z + g_yy_0_yzz_zzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzz_xx, g_yy_0_zzz_xxz, g_yy_0_zzz_xy, g_yy_0_zzz_xyz, g_yy_0_zzz_xz, g_yy_0_zzz_xzz, g_yy_0_zzz_yy, g_yy_0_zzz_yyz, g_yy_0_zzz_yz, g_yy_0_zzz_yzz, g_yy_0_zzz_zz, g_yy_0_zzz_zzz, g_yy_0_zzzz_xx, g_yy_0_zzzz_xy, g_yy_0_zzzz_xz, g_yy_0_zzzz_yy, g_yy_0_zzzz_yz, g_yy_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzz_xx[k] = -g_yy_0_zzz_xx[k] * ab_z + g_yy_0_zzz_xxz[k];

                g_yy_0_zzzz_xy[k] = -g_yy_0_zzz_xy[k] * ab_z + g_yy_0_zzz_xyz[k];

                g_yy_0_zzzz_xz[k] = -g_yy_0_zzz_xz[k] * ab_z + g_yy_0_zzz_xzz[k];

                g_yy_0_zzzz_yy[k] = -g_yy_0_zzz_yy[k] * ab_z + g_yy_0_zzz_yyz[k];

                g_yy_0_zzzz_yz[k] = -g_yy_0_zzz_yz[k] * ab_z + g_yy_0_zzz_yzz[k];

                g_yy_0_zzzz_zz[k] = -g_yy_0_zzz_zz[k] * ab_z + g_yy_0_zzz_zzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxx_xx, g_yz_0_xxx_xxx, g_yz_0_xxx_xxy, g_yz_0_xxx_xxz, g_yz_0_xxx_xy, g_yz_0_xxx_xyy, g_yz_0_xxx_xyz, g_yz_0_xxx_xz, g_yz_0_xxx_xzz, g_yz_0_xxx_yy, g_yz_0_xxx_yz, g_yz_0_xxx_zz, g_yz_0_xxxx_xx, g_yz_0_xxxx_xy, g_yz_0_xxxx_xz, g_yz_0_xxxx_yy, g_yz_0_xxxx_yz, g_yz_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxx_xx[k] = -g_yz_0_xxx_xx[k] * ab_x + g_yz_0_xxx_xxx[k];

                g_yz_0_xxxx_xy[k] = -g_yz_0_xxx_xy[k] * ab_x + g_yz_0_xxx_xxy[k];

                g_yz_0_xxxx_xz[k] = -g_yz_0_xxx_xz[k] * ab_x + g_yz_0_xxx_xxz[k];

                g_yz_0_xxxx_yy[k] = -g_yz_0_xxx_yy[k] * ab_x + g_yz_0_xxx_xyy[k];

                g_yz_0_xxxx_yz[k] = -g_yz_0_xxx_yz[k] * ab_x + g_yz_0_xxx_xyz[k];

                g_yz_0_xxxx_zz[k] = -g_yz_0_xxx_zz[k] * ab_x + g_yz_0_xxx_xzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxy_xx, g_yz_0_xxxy_xy, g_yz_0_xxxy_xz, g_yz_0_xxxy_yy, g_yz_0_xxxy_yz, g_yz_0_xxxy_zz, g_yz_0_xxy_xx, g_yz_0_xxy_xxx, g_yz_0_xxy_xxy, g_yz_0_xxy_xxz, g_yz_0_xxy_xy, g_yz_0_xxy_xyy, g_yz_0_xxy_xyz, g_yz_0_xxy_xz, g_yz_0_xxy_xzz, g_yz_0_xxy_yy, g_yz_0_xxy_yz, g_yz_0_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxy_xx[k] = -g_yz_0_xxy_xx[k] * ab_x + g_yz_0_xxy_xxx[k];

                g_yz_0_xxxy_xy[k] = -g_yz_0_xxy_xy[k] * ab_x + g_yz_0_xxy_xxy[k];

                g_yz_0_xxxy_xz[k] = -g_yz_0_xxy_xz[k] * ab_x + g_yz_0_xxy_xxz[k];

                g_yz_0_xxxy_yy[k] = -g_yz_0_xxy_yy[k] * ab_x + g_yz_0_xxy_xyy[k];

                g_yz_0_xxxy_yz[k] = -g_yz_0_xxy_yz[k] * ab_x + g_yz_0_xxy_xyz[k];

                g_yz_0_xxxy_zz[k] = -g_yz_0_xxy_zz[k] * ab_x + g_yz_0_xxy_xzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 374 * ccomps * dcomps);

            auto g_yz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxz_xx, g_yz_0_xxxz_xy, g_yz_0_xxxz_xz, g_yz_0_xxxz_yy, g_yz_0_xxxz_yz, g_yz_0_xxxz_zz, g_yz_0_xxz_xx, g_yz_0_xxz_xxx, g_yz_0_xxz_xxy, g_yz_0_xxz_xxz, g_yz_0_xxz_xy, g_yz_0_xxz_xyy, g_yz_0_xxz_xyz, g_yz_0_xxz_xz, g_yz_0_xxz_xzz, g_yz_0_xxz_yy, g_yz_0_xxz_yz, g_yz_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxz_xx[k] = -g_yz_0_xxz_xx[k] * ab_x + g_yz_0_xxz_xxx[k];

                g_yz_0_xxxz_xy[k] = -g_yz_0_xxz_xy[k] * ab_x + g_yz_0_xxz_xxy[k];

                g_yz_0_xxxz_xz[k] = -g_yz_0_xxz_xz[k] * ab_x + g_yz_0_xxz_xxz[k];

                g_yz_0_xxxz_yy[k] = -g_yz_0_xxz_yy[k] * ab_x + g_yz_0_xxz_xyy[k];

                g_yz_0_xxxz_yz[k] = -g_yz_0_xxz_yz[k] * ab_x + g_yz_0_xxz_xyz[k];

                g_yz_0_xxxz_zz[k] = -g_yz_0_xxz_zz[k] * ab_x + g_yz_0_xxz_xzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyy_xx, g_yz_0_xxyy_xy, g_yz_0_xxyy_xz, g_yz_0_xxyy_yy, g_yz_0_xxyy_yz, g_yz_0_xxyy_zz, g_yz_0_xyy_xx, g_yz_0_xyy_xxx, g_yz_0_xyy_xxy, g_yz_0_xyy_xxz, g_yz_0_xyy_xy, g_yz_0_xyy_xyy, g_yz_0_xyy_xyz, g_yz_0_xyy_xz, g_yz_0_xyy_xzz, g_yz_0_xyy_yy, g_yz_0_xyy_yz, g_yz_0_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyy_xx[k] = -g_yz_0_xyy_xx[k] * ab_x + g_yz_0_xyy_xxx[k];

                g_yz_0_xxyy_xy[k] = -g_yz_0_xyy_xy[k] * ab_x + g_yz_0_xyy_xxy[k];

                g_yz_0_xxyy_xz[k] = -g_yz_0_xyy_xz[k] * ab_x + g_yz_0_xyy_xxz[k];

                g_yz_0_xxyy_yy[k] = -g_yz_0_xyy_yy[k] * ab_x + g_yz_0_xyy_xyy[k];

                g_yz_0_xxyy_yz[k] = -g_yz_0_xyy_yz[k] * ab_x + g_yz_0_xyy_xyz[k];

                g_yz_0_xxyy_zz[k] = -g_yz_0_xyy_zz[k] * ab_x + g_yz_0_xyy_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyz_xx, g_yz_0_xxyz_xy, g_yz_0_xxyz_xz, g_yz_0_xxyz_yy, g_yz_0_xxyz_yz, g_yz_0_xxyz_zz, g_yz_0_xyz_xx, g_yz_0_xyz_xxx, g_yz_0_xyz_xxy, g_yz_0_xyz_xxz, g_yz_0_xyz_xy, g_yz_0_xyz_xyy, g_yz_0_xyz_xyz, g_yz_0_xyz_xz, g_yz_0_xyz_xzz, g_yz_0_xyz_yy, g_yz_0_xyz_yz, g_yz_0_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyz_xx[k] = -g_yz_0_xyz_xx[k] * ab_x + g_yz_0_xyz_xxx[k];

                g_yz_0_xxyz_xy[k] = -g_yz_0_xyz_xy[k] * ab_x + g_yz_0_xyz_xxy[k];

                g_yz_0_xxyz_xz[k] = -g_yz_0_xyz_xz[k] * ab_x + g_yz_0_xyz_xxz[k];

                g_yz_0_xxyz_yy[k] = -g_yz_0_xyz_yy[k] * ab_x + g_yz_0_xyz_xyy[k];

                g_yz_0_xxyz_yz[k] = -g_yz_0_xyz_yz[k] * ab_x + g_yz_0_xyz_xyz[k];

                g_yz_0_xxyz_zz[k] = -g_yz_0_xyz_zz[k] * ab_x + g_yz_0_xyz_xzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzz_xx, g_yz_0_xxzz_xy, g_yz_0_xxzz_xz, g_yz_0_xxzz_yy, g_yz_0_xxzz_yz, g_yz_0_xxzz_zz, g_yz_0_xzz_xx, g_yz_0_xzz_xxx, g_yz_0_xzz_xxy, g_yz_0_xzz_xxz, g_yz_0_xzz_xy, g_yz_0_xzz_xyy, g_yz_0_xzz_xyz, g_yz_0_xzz_xz, g_yz_0_xzz_xzz, g_yz_0_xzz_yy, g_yz_0_xzz_yz, g_yz_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzz_xx[k] = -g_yz_0_xzz_xx[k] * ab_x + g_yz_0_xzz_xxx[k];

                g_yz_0_xxzz_xy[k] = -g_yz_0_xzz_xy[k] * ab_x + g_yz_0_xzz_xxy[k];

                g_yz_0_xxzz_xz[k] = -g_yz_0_xzz_xz[k] * ab_x + g_yz_0_xzz_xxz[k];

                g_yz_0_xxzz_yy[k] = -g_yz_0_xzz_yy[k] * ab_x + g_yz_0_xzz_xyy[k];

                g_yz_0_xxzz_yz[k] = -g_yz_0_xzz_yz[k] * ab_x + g_yz_0_xzz_xyz[k];

                g_yz_0_xxzz_zz[k] = -g_yz_0_xzz_zz[k] * ab_x + g_yz_0_xzz_xzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyy_xx, g_yz_0_xyyy_xy, g_yz_0_xyyy_xz, g_yz_0_xyyy_yy, g_yz_0_xyyy_yz, g_yz_0_xyyy_zz, g_yz_0_yyy_xx, g_yz_0_yyy_xxx, g_yz_0_yyy_xxy, g_yz_0_yyy_xxz, g_yz_0_yyy_xy, g_yz_0_yyy_xyy, g_yz_0_yyy_xyz, g_yz_0_yyy_xz, g_yz_0_yyy_xzz, g_yz_0_yyy_yy, g_yz_0_yyy_yz, g_yz_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyy_xx[k] = -g_yz_0_yyy_xx[k] * ab_x + g_yz_0_yyy_xxx[k];

                g_yz_0_xyyy_xy[k] = -g_yz_0_yyy_xy[k] * ab_x + g_yz_0_yyy_xxy[k];

                g_yz_0_xyyy_xz[k] = -g_yz_0_yyy_xz[k] * ab_x + g_yz_0_yyy_xxz[k];

                g_yz_0_xyyy_yy[k] = -g_yz_0_yyy_yy[k] * ab_x + g_yz_0_yyy_xyy[k];

                g_yz_0_xyyy_yz[k] = -g_yz_0_yyy_yz[k] * ab_x + g_yz_0_yyy_xyz[k];

                g_yz_0_xyyy_zz[k] = -g_yz_0_yyy_zz[k] * ab_x + g_yz_0_yyy_xzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyz_xx, g_yz_0_xyyz_xy, g_yz_0_xyyz_xz, g_yz_0_xyyz_yy, g_yz_0_xyyz_yz, g_yz_0_xyyz_zz, g_yz_0_yyz_xx, g_yz_0_yyz_xxx, g_yz_0_yyz_xxy, g_yz_0_yyz_xxz, g_yz_0_yyz_xy, g_yz_0_yyz_xyy, g_yz_0_yyz_xyz, g_yz_0_yyz_xz, g_yz_0_yyz_xzz, g_yz_0_yyz_yy, g_yz_0_yyz_yz, g_yz_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyz_xx[k] = -g_yz_0_yyz_xx[k] * ab_x + g_yz_0_yyz_xxx[k];

                g_yz_0_xyyz_xy[k] = -g_yz_0_yyz_xy[k] * ab_x + g_yz_0_yyz_xxy[k];

                g_yz_0_xyyz_xz[k] = -g_yz_0_yyz_xz[k] * ab_x + g_yz_0_yyz_xxz[k];

                g_yz_0_xyyz_yy[k] = -g_yz_0_yyz_yy[k] * ab_x + g_yz_0_yyz_xyy[k];

                g_yz_0_xyyz_yz[k] = -g_yz_0_yyz_yz[k] * ab_x + g_yz_0_yyz_xyz[k];

                g_yz_0_xyyz_zz[k] = -g_yz_0_yyz_zz[k] * ab_x + g_yz_0_yyz_xzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzz_xx, g_yz_0_xyzz_xy, g_yz_0_xyzz_xz, g_yz_0_xyzz_yy, g_yz_0_xyzz_yz, g_yz_0_xyzz_zz, g_yz_0_yzz_xx, g_yz_0_yzz_xxx, g_yz_0_yzz_xxy, g_yz_0_yzz_xxz, g_yz_0_yzz_xy, g_yz_0_yzz_xyy, g_yz_0_yzz_xyz, g_yz_0_yzz_xz, g_yz_0_yzz_xzz, g_yz_0_yzz_yy, g_yz_0_yzz_yz, g_yz_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzz_xx[k] = -g_yz_0_yzz_xx[k] * ab_x + g_yz_0_yzz_xxx[k];

                g_yz_0_xyzz_xy[k] = -g_yz_0_yzz_xy[k] * ab_x + g_yz_0_yzz_xxy[k];

                g_yz_0_xyzz_xz[k] = -g_yz_0_yzz_xz[k] * ab_x + g_yz_0_yzz_xxz[k];

                g_yz_0_xyzz_yy[k] = -g_yz_0_yzz_yy[k] * ab_x + g_yz_0_yzz_xyy[k];

                g_yz_0_xyzz_yz[k] = -g_yz_0_yzz_yz[k] * ab_x + g_yz_0_yzz_xyz[k];

                g_yz_0_xyzz_zz[k] = -g_yz_0_yzz_zz[k] * ab_x + g_yz_0_yzz_xzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzz_xx, g_yz_0_xzzz_xy, g_yz_0_xzzz_xz, g_yz_0_xzzz_yy, g_yz_0_xzzz_yz, g_yz_0_xzzz_zz, g_yz_0_zzz_xx, g_yz_0_zzz_xxx, g_yz_0_zzz_xxy, g_yz_0_zzz_xxz, g_yz_0_zzz_xy, g_yz_0_zzz_xyy, g_yz_0_zzz_xyz, g_yz_0_zzz_xz, g_yz_0_zzz_xzz, g_yz_0_zzz_yy, g_yz_0_zzz_yz, g_yz_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzz_xx[k] = -g_yz_0_zzz_xx[k] * ab_x + g_yz_0_zzz_xxx[k];

                g_yz_0_xzzz_xy[k] = -g_yz_0_zzz_xy[k] * ab_x + g_yz_0_zzz_xxy[k];

                g_yz_0_xzzz_xz[k] = -g_yz_0_zzz_xz[k] * ab_x + g_yz_0_zzz_xxz[k];

                g_yz_0_xzzz_yy[k] = -g_yz_0_zzz_yy[k] * ab_x + g_yz_0_zzz_xyy[k];

                g_yz_0_xzzz_yz[k] = -g_yz_0_zzz_yz[k] * ab_x + g_yz_0_zzz_xyz[k];

                g_yz_0_xzzz_zz[k] = -g_yz_0_zzz_zz[k] * ab_x + g_yz_0_zzz_xzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyy_xx, g_yz_0_yyy_xxy, g_yz_0_yyy_xy, g_yz_0_yyy_xyy, g_yz_0_yyy_xyz, g_yz_0_yyy_xz, g_yz_0_yyy_yy, g_yz_0_yyy_yyy, g_yz_0_yyy_yyz, g_yz_0_yyy_yz, g_yz_0_yyy_yzz, g_yz_0_yyy_zz, g_yz_0_yyyy_xx, g_yz_0_yyyy_xy, g_yz_0_yyyy_xz, g_yz_0_yyyy_yy, g_yz_0_yyyy_yz, g_yz_0_yyyy_zz, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyy_xx[k] = -g_z_0_yyy_xx[k] - g_yz_0_yyy_xx[k] * ab_y + g_yz_0_yyy_xxy[k];

                g_yz_0_yyyy_xy[k] = -g_z_0_yyy_xy[k] - g_yz_0_yyy_xy[k] * ab_y + g_yz_0_yyy_xyy[k];

                g_yz_0_yyyy_xz[k] = -g_z_0_yyy_xz[k] - g_yz_0_yyy_xz[k] * ab_y + g_yz_0_yyy_xyz[k];

                g_yz_0_yyyy_yy[k] = -g_z_0_yyy_yy[k] - g_yz_0_yyy_yy[k] * ab_y + g_yz_0_yyy_yyy[k];

                g_yz_0_yyyy_yz[k] = -g_z_0_yyy_yz[k] - g_yz_0_yyy_yz[k] * ab_y + g_yz_0_yyy_yyz[k];

                g_yz_0_yyyy_zz[k] = -g_z_0_yyy_zz[k] - g_yz_0_yyy_zz[k] * ab_y + g_yz_0_yyy_yzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyz_xx, g_yz_0_yyyz_xy, g_yz_0_yyyz_xz, g_yz_0_yyyz_yy, g_yz_0_yyyz_yz, g_yz_0_yyyz_zz, g_yz_0_yyz_xx, g_yz_0_yyz_xxy, g_yz_0_yyz_xy, g_yz_0_yyz_xyy, g_yz_0_yyz_xyz, g_yz_0_yyz_xz, g_yz_0_yyz_yy, g_yz_0_yyz_yyy, g_yz_0_yyz_yyz, g_yz_0_yyz_yz, g_yz_0_yyz_yzz, g_yz_0_yyz_zz, g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyz_xx[k] = -g_z_0_yyz_xx[k] - g_yz_0_yyz_xx[k] * ab_y + g_yz_0_yyz_xxy[k];

                g_yz_0_yyyz_xy[k] = -g_z_0_yyz_xy[k] - g_yz_0_yyz_xy[k] * ab_y + g_yz_0_yyz_xyy[k];

                g_yz_0_yyyz_xz[k] = -g_z_0_yyz_xz[k] - g_yz_0_yyz_xz[k] * ab_y + g_yz_0_yyz_xyz[k];

                g_yz_0_yyyz_yy[k] = -g_z_0_yyz_yy[k] - g_yz_0_yyz_yy[k] * ab_y + g_yz_0_yyz_yyy[k];

                g_yz_0_yyyz_yz[k] = -g_z_0_yyz_yz[k] - g_yz_0_yyz_yz[k] * ab_y + g_yz_0_yyz_yyz[k];

                g_yz_0_yyyz_zz[k] = -g_z_0_yyz_zz[k] - g_yz_0_yyz_zz[k] * ab_y + g_yz_0_yyz_yzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzz_xx, g_yz_0_yyzz_xy, g_yz_0_yyzz_xz, g_yz_0_yyzz_yy, g_yz_0_yyzz_yz, g_yz_0_yyzz_zz, g_yz_0_yzz_xx, g_yz_0_yzz_xxy, g_yz_0_yzz_xy, g_yz_0_yzz_xyy, g_yz_0_yzz_xyz, g_yz_0_yzz_xz, g_yz_0_yzz_yy, g_yz_0_yzz_yyy, g_yz_0_yzz_yyz, g_yz_0_yzz_yz, g_yz_0_yzz_yzz, g_yz_0_yzz_zz, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzz_xx[k] = -g_z_0_yzz_xx[k] - g_yz_0_yzz_xx[k] * ab_y + g_yz_0_yzz_xxy[k];

                g_yz_0_yyzz_xy[k] = -g_z_0_yzz_xy[k] - g_yz_0_yzz_xy[k] * ab_y + g_yz_0_yzz_xyy[k];

                g_yz_0_yyzz_xz[k] = -g_z_0_yzz_xz[k] - g_yz_0_yzz_xz[k] * ab_y + g_yz_0_yzz_xyz[k];

                g_yz_0_yyzz_yy[k] = -g_z_0_yzz_yy[k] - g_yz_0_yzz_yy[k] * ab_y + g_yz_0_yzz_yyy[k];

                g_yz_0_yyzz_yz[k] = -g_z_0_yzz_yz[k] - g_yz_0_yzz_yz[k] * ab_y + g_yz_0_yzz_yyz[k];

                g_yz_0_yyzz_zz[k] = -g_z_0_yzz_zz[k] - g_yz_0_yzz_zz[k] * ab_y + g_yz_0_yzz_yzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzz_xx, g_yz_0_yzzz_xy, g_yz_0_yzzz_xz, g_yz_0_yzzz_yy, g_yz_0_yzzz_yz, g_yz_0_yzzz_zz, g_yz_0_zzz_xx, g_yz_0_zzz_xxy, g_yz_0_zzz_xy, g_yz_0_zzz_xyy, g_yz_0_zzz_xyz, g_yz_0_zzz_xz, g_yz_0_zzz_yy, g_yz_0_zzz_yyy, g_yz_0_zzz_yyz, g_yz_0_zzz_yz, g_yz_0_zzz_yzz, g_yz_0_zzz_zz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzz_xx[k] = -g_z_0_zzz_xx[k] - g_yz_0_zzz_xx[k] * ab_y + g_yz_0_zzz_xxy[k];

                g_yz_0_yzzz_xy[k] = -g_z_0_zzz_xy[k] - g_yz_0_zzz_xy[k] * ab_y + g_yz_0_zzz_xyy[k];

                g_yz_0_yzzz_xz[k] = -g_z_0_zzz_xz[k] - g_yz_0_zzz_xz[k] * ab_y + g_yz_0_zzz_xyz[k];

                g_yz_0_yzzz_yy[k] = -g_z_0_zzz_yy[k] - g_yz_0_zzz_yy[k] * ab_y + g_yz_0_zzz_yyy[k];

                g_yz_0_yzzz_yz[k] = -g_z_0_zzz_yz[k] - g_yz_0_zzz_yz[k] * ab_y + g_yz_0_zzz_yyz[k];

                g_yz_0_yzzz_zz[k] = -g_z_0_zzz_zz[k] - g_yz_0_zzz_zz[k] * ab_y + g_yz_0_zzz_yzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz, g_yz_0_zzz_xx, g_yz_0_zzz_xxz, g_yz_0_zzz_xy, g_yz_0_zzz_xyz, g_yz_0_zzz_xz, g_yz_0_zzz_xzz, g_yz_0_zzz_yy, g_yz_0_zzz_yyz, g_yz_0_zzz_yz, g_yz_0_zzz_yzz, g_yz_0_zzz_zz, g_yz_0_zzz_zzz, g_yz_0_zzzz_xx, g_yz_0_zzzz_xy, g_yz_0_zzzz_xz, g_yz_0_zzzz_yy, g_yz_0_zzzz_yz, g_yz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzz_xx[k] = -g_y_0_zzz_xx[k] - g_yz_0_zzz_xx[k] * ab_z + g_yz_0_zzz_xxz[k];

                g_yz_0_zzzz_xy[k] = -g_y_0_zzz_xy[k] - g_yz_0_zzz_xy[k] * ab_z + g_yz_0_zzz_xyz[k];

                g_yz_0_zzzz_xz[k] = -g_y_0_zzz_xz[k] - g_yz_0_zzz_xz[k] * ab_z + g_yz_0_zzz_xzz[k];

                g_yz_0_zzzz_yy[k] = -g_y_0_zzz_yy[k] - g_yz_0_zzz_yy[k] * ab_z + g_yz_0_zzz_yyz[k];

                g_yz_0_zzzz_yz[k] = -g_y_0_zzz_yz[k] - g_yz_0_zzz_yz[k] * ab_z + g_yz_0_zzz_yzz[k];

                g_yz_0_zzzz_zz[k] = -g_y_0_zzz_zz[k] - g_yz_0_zzz_zz[k] * ab_z + g_yz_0_zzz_zzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxx_xx = cbuffer.data(gd_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xxxx_xy = cbuffer.data(gd_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xxxx_xz = cbuffer.data(gd_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_xxxx_yy = cbuffer.data(gd_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xxxx_yz = cbuffer.data(gd_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xxxx_zz = cbuffer.data(gd_geom_20_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxx_xx, g_zz_0_xxx_xxx, g_zz_0_xxx_xxy, g_zz_0_xxx_xxz, g_zz_0_xxx_xy, g_zz_0_xxx_xyy, g_zz_0_xxx_xyz, g_zz_0_xxx_xz, g_zz_0_xxx_xzz, g_zz_0_xxx_yy, g_zz_0_xxx_yz, g_zz_0_xxx_zz, g_zz_0_xxxx_xx, g_zz_0_xxxx_xy, g_zz_0_xxxx_xz, g_zz_0_xxxx_yy, g_zz_0_xxxx_yz, g_zz_0_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxx_xx[k] = -g_zz_0_xxx_xx[k] * ab_x + g_zz_0_xxx_xxx[k];

                g_zz_0_xxxx_xy[k] = -g_zz_0_xxx_xy[k] * ab_x + g_zz_0_xxx_xxy[k];

                g_zz_0_xxxx_xz[k] = -g_zz_0_xxx_xz[k] * ab_x + g_zz_0_xxx_xxz[k];

                g_zz_0_xxxx_yy[k] = -g_zz_0_xxx_yy[k] * ab_x + g_zz_0_xxx_xyy[k];

                g_zz_0_xxxx_yz[k] = -g_zz_0_xxx_yz[k] * ab_x + g_zz_0_xxx_xyz[k];

                g_zz_0_xxxx_zz[k] = -g_zz_0_xxx_zz[k] * ab_x + g_zz_0_xxx_xzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxy_xx = cbuffer.data(gd_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xxxy_xy = cbuffer.data(gd_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xxxy_xz = cbuffer.data(gd_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_xxxy_yy = cbuffer.data(gd_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xxxy_yz = cbuffer.data(gd_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xxxy_zz = cbuffer.data(gd_geom_20_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxy_xx, g_zz_0_xxxy_xy, g_zz_0_xxxy_xz, g_zz_0_xxxy_yy, g_zz_0_xxxy_yz, g_zz_0_xxxy_zz, g_zz_0_xxy_xx, g_zz_0_xxy_xxx, g_zz_0_xxy_xxy, g_zz_0_xxy_xxz, g_zz_0_xxy_xy, g_zz_0_xxy_xyy, g_zz_0_xxy_xyz, g_zz_0_xxy_xz, g_zz_0_xxy_xzz, g_zz_0_xxy_yy, g_zz_0_xxy_yz, g_zz_0_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxy_xx[k] = -g_zz_0_xxy_xx[k] * ab_x + g_zz_0_xxy_xxx[k];

                g_zz_0_xxxy_xy[k] = -g_zz_0_xxy_xy[k] * ab_x + g_zz_0_xxy_xxy[k];

                g_zz_0_xxxy_xz[k] = -g_zz_0_xxy_xz[k] * ab_x + g_zz_0_xxy_xxz[k];

                g_zz_0_xxxy_yy[k] = -g_zz_0_xxy_yy[k] * ab_x + g_zz_0_xxy_xyy[k];

                g_zz_0_xxxy_yz[k] = -g_zz_0_xxy_yz[k] * ab_x + g_zz_0_xxy_xyz[k];

                g_zz_0_xxxy_zz[k] = -g_zz_0_xxy_zz[k] * ab_x + g_zz_0_xxy_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxz_xx = cbuffer.data(gd_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xxxz_xy = cbuffer.data(gd_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xxxz_xz = cbuffer.data(gd_geom_20_off + 464 * ccomps * dcomps);

            auto g_zz_0_xxxz_yy = cbuffer.data(gd_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xxxz_yz = cbuffer.data(gd_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xxxz_zz = cbuffer.data(gd_geom_20_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxz_xx, g_zz_0_xxxz_xy, g_zz_0_xxxz_xz, g_zz_0_xxxz_yy, g_zz_0_xxxz_yz, g_zz_0_xxxz_zz, g_zz_0_xxz_xx, g_zz_0_xxz_xxx, g_zz_0_xxz_xxy, g_zz_0_xxz_xxz, g_zz_0_xxz_xy, g_zz_0_xxz_xyy, g_zz_0_xxz_xyz, g_zz_0_xxz_xz, g_zz_0_xxz_xzz, g_zz_0_xxz_yy, g_zz_0_xxz_yz, g_zz_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxz_xx[k] = -g_zz_0_xxz_xx[k] * ab_x + g_zz_0_xxz_xxx[k];

                g_zz_0_xxxz_xy[k] = -g_zz_0_xxz_xy[k] * ab_x + g_zz_0_xxz_xxy[k];

                g_zz_0_xxxz_xz[k] = -g_zz_0_xxz_xz[k] * ab_x + g_zz_0_xxz_xxz[k];

                g_zz_0_xxxz_yy[k] = -g_zz_0_xxz_yy[k] * ab_x + g_zz_0_xxz_xyy[k];

                g_zz_0_xxxz_yz[k] = -g_zz_0_xxz_yz[k] * ab_x + g_zz_0_xxz_xyz[k];

                g_zz_0_xxxz_zz[k] = -g_zz_0_xxz_zz[k] * ab_x + g_zz_0_xxz_xzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyy_xx = cbuffer.data(gd_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xxyy_xy = cbuffer.data(gd_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xxyy_xz = cbuffer.data(gd_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_xxyy_yy = cbuffer.data(gd_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xxyy_yz = cbuffer.data(gd_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xxyy_zz = cbuffer.data(gd_geom_20_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyy_xx, g_zz_0_xxyy_xy, g_zz_0_xxyy_xz, g_zz_0_xxyy_yy, g_zz_0_xxyy_yz, g_zz_0_xxyy_zz, g_zz_0_xyy_xx, g_zz_0_xyy_xxx, g_zz_0_xyy_xxy, g_zz_0_xyy_xxz, g_zz_0_xyy_xy, g_zz_0_xyy_xyy, g_zz_0_xyy_xyz, g_zz_0_xyy_xz, g_zz_0_xyy_xzz, g_zz_0_xyy_yy, g_zz_0_xyy_yz, g_zz_0_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyy_xx[k] = -g_zz_0_xyy_xx[k] * ab_x + g_zz_0_xyy_xxx[k];

                g_zz_0_xxyy_xy[k] = -g_zz_0_xyy_xy[k] * ab_x + g_zz_0_xyy_xxy[k];

                g_zz_0_xxyy_xz[k] = -g_zz_0_xyy_xz[k] * ab_x + g_zz_0_xyy_xxz[k];

                g_zz_0_xxyy_yy[k] = -g_zz_0_xyy_yy[k] * ab_x + g_zz_0_xyy_xyy[k];

                g_zz_0_xxyy_yz[k] = -g_zz_0_xyy_yz[k] * ab_x + g_zz_0_xyy_xyz[k];

                g_zz_0_xxyy_zz[k] = -g_zz_0_xyy_zz[k] * ab_x + g_zz_0_xyy_xzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyz_xx = cbuffer.data(gd_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xxyz_xy = cbuffer.data(gd_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xxyz_xz = cbuffer.data(gd_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_xxyz_yy = cbuffer.data(gd_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xxyz_yz = cbuffer.data(gd_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xxyz_zz = cbuffer.data(gd_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyz_xx, g_zz_0_xxyz_xy, g_zz_0_xxyz_xz, g_zz_0_xxyz_yy, g_zz_0_xxyz_yz, g_zz_0_xxyz_zz, g_zz_0_xyz_xx, g_zz_0_xyz_xxx, g_zz_0_xyz_xxy, g_zz_0_xyz_xxz, g_zz_0_xyz_xy, g_zz_0_xyz_xyy, g_zz_0_xyz_xyz, g_zz_0_xyz_xz, g_zz_0_xyz_xzz, g_zz_0_xyz_yy, g_zz_0_xyz_yz, g_zz_0_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyz_xx[k] = -g_zz_0_xyz_xx[k] * ab_x + g_zz_0_xyz_xxx[k];

                g_zz_0_xxyz_xy[k] = -g_zz_0_xyz_xy[k] * ab_x + g_zz_0_xyz_xxy[k];

                g_zz_0_xxyz_xz[k] = -g_zz_0_xyz_xz[k] * ab_x + g_zz_0_xyz_xxz[k];

                g_zz_0_xxyz_yy[k] = -g_zz_0_xyz_yy[k] * ab_x + g_zz_0_xyz_xyy[k];

                g_zz_0_xxyz_yz[k] = -g_zz_0_xyz_yz[k] * ab_x + g_zz_0_xyz_xyz[k];

                g_zz_0_xxyz_zz[k] = -g_zz_0_xyz_zz[k] * ab_x + g_zz_0_xyz_xzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzz_xx = cbuffer.data(gd_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xxzz_xy = cbuffer.data(gd_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xxzz_xz = cbuffer.data(gd_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_xxzz_yy = cbuffer.data(gd_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_xxzz_yz = cbuffer.data(gd_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_xxzz_zz = cbuffer.data(gd_geom_20_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzz_xx, g_zz_0_xxzz_xy, g_zz_0_xxzz_xz, g_zz_0_xxzz_yy, g_zz_0_xxzz_yz, g_zz_0_xxzz_zz, g_zz_0_xzz_xx, g_zz_0_xzz_xxx, g_zz_0_xzz_xxy, g_zz_0_xzz_xxz, g_zz_0_xzz_xy, g_zz_0_xzz_xyy, g_zz_0_xzz_xyz, g_zz_0_xzz_xz, g_zz_0_xzz_xzz, g_zz_0_xzz_yy, g_zz_0_xzz_yz, g_zz_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzz_xx[k] = -g_zz_0_xzz_xx[k] * ab_x + g_zz_0_xzz_xxx[k];

                g_zz_0_xxzz_xy[k] = -g_zz_0_xzz_xy[k] * ab_x + g_zz_0_xzz_xxy[k];

                g_zz_0_xxzz_xz[k] = -g_zz_0_xzz_xz[k] * ab_x + g_zz_0_xzz_xxz[k];

                g_zz_0_xxzz_yy[k] = -g_zz_0_xzz_yy[k] * ab_x + g_zz_0_xzz_xyy[k];

                g_zz_0_xxzz_yz[k] = -g_zz_0_xzz_yz[k] * ab_x + g_zz_0_xzz_xyz[k];

                g_zz_0_xxzz_zz[k] = -g_zz_0_xzz_zz[k] * ab_x + g_zz_0_xzz_xzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyy_xx = cbuffer.data(gd_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_xyyy_xy = cbuffer.data(gd_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_xyyy_xz = cbuffer.data(gd_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_xyyy_yy = cbuffer.data(gd_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_xyyy_yz = cbuffer.data(gd_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_xyyy_zz = cbuffer.data(gd_geom_20_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyy_xx, g_zz_0_xyyy_xy, g_zz_0_xyyy_xz, g_zz_0_xyyy_yy, g_zz_0_xyyy_yz, g_zz_0_xyyy_zz, g_zz_0_yyy_xx, g_zz_0_yyy_xxx, g_zz_0_yyy_xxy, g_zz_0_yyy_xxz, g_zz_0_yyy_xy, g_zz_0_yyy_xyy, g_zz_0_yyy_xyz, g_zz_0_yyy_xz, g_zz_0_yyy_xzz, g_zz_0_yyy_yy, g_zz_0_yyy_yz, g_zz_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyy_xx[k] = -g_zz_0_yyy_xx[k] * ab_x + g_zz_0_yyy_xxx[k];

                g_zz_0_xyyy_xy[k] = -g_zz_0_yyy_xy[k] * ab_x + g_zz_0_yyy_xxy[k];

                g_zz_0_xyyy_xz[k] = -g_zz_0_yyy_xz[k] * ab_x + g_zz_0_yyy_xxz[k];

                g_zz_0_xyyy_yy[k] = -g_zz_0_yyy_yy[k] * ab_x + g_zz_0_yyy_xyy[k];

                g_zz_0_xyyy_yz[k] = -g_zz_0_yyy_yz[k] * ab_x + g_zz_0_yyy_xyz[k];

                g_zz_0_xyyy_zz[k] = -g_zz_0_yyy_zz[k] * ab_x + g_zz_0_yyy_xzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyz_xx = cbuffer.data(gd_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_xyyz_xy = cbuffer.data(gd_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_xyyz_xz = cbuffer.data(gd_geom_20_off + 494 * ccomps * dcomps);

            auto g_zz_0_xyyz_yy = cbuffer.data(gd_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_xyyz_yz = cbuffer.data(gd_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_xyyz_zz = cbuffer.data(gd_geom_20_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyz_xx, g_zz_0_xyyz_xy, g_zz_0_xyyz_xz, g_zz_0_xyyz_yy, g_zz_0_xyyz_yz, g_zz_0_xyyz_zz, g_zz_0_yyz_xx, g_zz_0_yyz_xxx, g_zz_0_yyz_xxy, g_zz_0_yyz_xxz, g_zz_0_yyz_xy, g_zz_0_yyz_xyy, g_zz_0_yyz_xyz, g_zz_0_yyz_xz, g_zz_0_yyz_xzz, g_zz_0_yyz_yy, g_zz_0_yyz_yz, g_zz_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyz_xx[k] = -g_zz_0_yyz_xx[k] * ab_x + g_zz_0_yyz_xxx[k];

                g_zz_0_xyyz_xy[k] = -g_zz_0_yyz_xy[k] * ab_x + g_zz_0_yyz_xxy[k];

                g_zz_0_xyyz_xz[k] = -g_zz_0_yyz_xz[k] * ab_x + g_zz_0_yyz_xxz[k];

                g_zz_0_xyyz_yy[k] = -g_zz_0_yyz_yy[k] * ab_x + g_zz_0_yyz_xyy[k];

                g_zz_0_xyyz_yz[k] = -g_zz_0_yyz_yz[k] * ab_x + g_zz_0_yyz_xyz[k];

                g_zz_0_xyyz_zz[k] = -g_zz_0_yyz_zz[k] * ab_x + g_zz_0_yyz_xzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzz_xx = cbuffer.data(gd_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_xyzz_xy = cbuffer.data(gd_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_xyzz_xz = cbuffer.data(gd_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_xyzz_yy = cbuffer.data(gd_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_xyzz_yz = cbuffer.data(gd_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_xyzz_zz = cbuffer.data(gd_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzz_xx, g_zz_0_xyzz_xy, g_zz_0_xyzz_xz, g_zz_0_xyzz_yy, g_zz_0_xyzz_yz, g_zz_0_xyzz_zz, g_zz_0_yzz_xx, g_zz_0_yzz_xxx, g_zz_0_yzz_xxy, g_zz_0_yzz_xxz, g_zz_0_yzz_xy, g_zz_0_yzz_xyy, g_zz_0_yzz_xyz, g_zz_0_yzz_xz, g_zz_0_yzz_xzz, g_zz_0_yzz_yy, g_zz_0_yzz_yz, g_zz_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzz_xx[k] = -g_zz_0_yzz_xx[k] * ab_x + g_zz_0_yzz_xxx[k];

                g_zz_0_xyzz_xy[k] = -g_zz_0_yzz_xy[k] * ab_x + g_zz_0_yzz_xxy[k];

                g_zz_0_xyzz_xz[k] = -g_zz_0_yzz_xz[k] * ab_x + g_zz_0_yzz_xxz[k];

                g_zz_0_xyzz_yy[k] = -g_zz_0_yzz_yy[k] * ab_x + g_zz_0_yzz_xyy[k];

                g_zz_0_xyzz_yz[k] = -g_zz_0_yzz_yz[k] * ab_x + g_zz_0_yzz_xyz[k];

                g_zz_0_xyzz_zz[k] = -g_zz_0_yzz_zz[k] * ab_x + g_zz_0_yzz_xzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzz_xx = cbuffer.data(gd_geom_20_off + 504 * ccomps * dcomps);

            auto g_zz_0_xzzz_xy = cbuffer.data(gd_geom_20_off + 505 * ccomps * dcomps);

            auto g_zz_0_xzzz_xz = cbuffer.data(gd_geom_20_off + 506 * ccomps * dcomps);

            auto g_zz_0_xzzz_yy = cbuffer.data(gd_geom_20_off + 507 * ccomps * dcomps);

            auto g_zz_0_xzzz_yz = cbuffer.data(gd_geom_20_off + 508 * ccomps * dcomps);

            auto g_zz_0_xzzz_zz = cbuffer.data(gd_geom_20_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzz_xx, g_zz_0_xzzz_xy, g_zz_0_xzzz_xz, g_zz_0_xzzz_yy, g_zz_0_xzzz_yz, g_zz_0_xzzz_zz, g_zz_0_zzz_xx, g_zz_0_zzz_xxx, g_zz_0_zzz_xxy, g_zz_0_zzz_xxz, g_zz_0_zzz_xy, g_zz_0_zzz_xyy, g_zz_0_zzz_xyz, g_zz_0_zzz_xz, g_zz_0_zzz_xzz, g_zz_0_zzz_yy, g_zz_0_zzz_yz, g_zz_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzz_xx[k] = -g_zz_0_zzz_xx[k] * ab_x + g_zz_0_zzz_xxx[k];

                g_zz_0_xzzz_xy[k] = -g_zz_0_zzz_xy[k] * ab_x + g_zz_0_zzz_xxy[k];

                g_zz_0_xzzz_xz[k] = -g_zz_0_zzz_xz[k] * ab_x + g_zz_0_zzz_xxz[k];

                g_zz_0_xzzz_yy[k] = -g_zz_0_zzz_yy[k] * ab_x + g_zz_0_zzz_xyy[k];

                g_zz_0_xzzz_yz[k] = -g_zz_0_zzz_yz[k] * ab_x + g_zz_0_zzz_xyz[k];

                g_zz_0_xzzz_zz[k] = -g_zz_0_zzz_zz[k] * ab_x + g_zz_0_zzz_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyy_xx = cbuffer.data(gd_geom_20_off + 510 * ccomps * dcomps);

            auto g_zz_0_yyyy_xy = cbuffer.data(gd_geom_20_off + 511 * ccomps * dcomps);

            auto g_zz_0_yyyy_xz = cbuffer.data(gd_geom_20_off + 512 * ccomps * dcomps);

            auto g_zz_0_yyyy_yy = cbuffer.data(gd_geom_20_off + 513 * ccomps * dcomps);

            auto g_zz_0_yyyy_yz = cbuffer.data(gd_geom_20_off + 514 * ccomps * dcomps);

            auto g_zz_0_yyyy_zz = cbuffer.data(gd_geom_20_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyy_xx, g_zz_0_yyy_xxy, g_zz_0_yyy_xy, g_zz_0_yyy_xyy, g_zz_0_yyy_xyz, g_zz_0_yyy_xz, g_zz_0_yyy_yy, g_zz_0_yyy_yyy, g_zz_0_yyy_yyz, g_zz_0_yyy_yz, g_zz_0_yyy_yzz, g_zz_0_yyy_zz, g_zz_0_yyyy_xx, g_zz_0_yyyy_xy, g_zz_0_yyyy_xz, g_zz_0_yyyy_yy, g_zz_0_yyyy_yz, g_zz_0_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyy_xx[k] = -g_zz_0_yyy_xx[k] * ab_y + g_zz_0_yyy_xxy[k];

                g_zz_0_yyyy_xy[k] = -g_zz_0_yyy_xy[k] * ab_y + g_zz_0_yyy_xyy[k];

                g_zz_0_yyyy_xz[k] = -g_zz_0_yyy_xz[k] * ab_y + g_zz_0_yyy_xyz[k];

                g_zz_0_yyyy_yy[k] = -g_zz_0_yyy_yy[k] * ab_y + g_zz_0_yyy_yyy[k];

                g_zz_0_yyyy_yz[k] = -g_zz_0_yyy_yz[k] * ab_y + g_zz_0_yyy_yyz[k];

                g_zz_0_yyyy_zz[k] = -g_zz_0_yyy_zz[k] * ab_y + g_zz_0_yyy_yzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyz_xx = cbuffer.data(gd_geom_20_off + 516 * ccomps * dcomps);

            auto g_zz_0_yyyz_xy = cbuffer.data(gd_geom_20_off + 517 * ccomps * dcomps);

            auto g_zz_0_yyyz_xz = cbuffer.data(gd_geom_20_off + 518 * ccomps * dcomps);

            auto g_zz_0_yyyz_yy = cbuffer.data(gd_geom_20_off + 519 * ccomps * dcomps);

            auto g_zz_0_yyyz_yz = cbuffer.data(gd_geom_20_off + 520 * ccomps * dcomps);

            auto g_zz_0_yyyz_zz = cbuffer.data(gd_geom_20_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyz_xx, g_zz_0_yyyz_xy, g_zz_0_yyyz_xz, g_zz_0_yyyz_yy, g_zz_0_yyyz_yz, g_zz_0_yyyz_zz, g_zz_0_yyz_xx, g_zz_0_yyz_xxy, g_zz_0_yyz_xy, g_zz_0_yyz_xyy, g_zz_0_yyz_xyz, g_zz_0_yyz_xz, g_zz_0_yyz_yy, g_zz_0_yyz_yyy, g_zz_0_yyz_yyz, g_zz_0_yyz_yz, g_zz_0_yyz_yzz, g_zz_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyz_xx[k] = -g_zz_0_yyz_xx[k] * ab_y + g_zz_0_yyz_xxy[k];

                g_zz_0_yyyz_xy[k] = -g_zz_0_yyz_xy[k] * ab_y + g_zz_0_yyz_xyy[k];

                g_zz_0_yyyz_xz[k] = -g_zz_0_yyz_xz[k] * ab_y + g_zz_0_yyz_xyz[k];

                g_zz_0_yyyz_yy[k] = -g_zz_0_yyz_yy[k] * ab_y + g_zz_0_yyz_yyy[k];

                g_zz_0_yyyz_yz[k] = -g_zz_0_yyz_yz[k] * ab_y + g_zz_0_yyz_yyz[k];

                g_zz_0_yyyz_zz[k] = -g_zz_0_yyz_zz[k] * ab_y + g_zz_0_yyz_yzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzz_xx = cbuffer.data(gd_geom_20_off + 522 * ccomps * dcomps);

            auto g_zz_0_yyzz_xy = cbuffer.data(gd_geom_20_off + 523 * ccomps * dcomps);

            auto g_zz_0_yyzz_xz = cbuffer.data(gd_geom_20_off + 524 * ccomps * dcomps);

            auto g_zz_0_yyzz_yy = cbuffer.data(gd_geom_20_off + 525 * ccomps * dcomps);

            auto g_zz_0_yyzz_yz = cbuffer.data(gd_geom_20_off + 526 * ccomps * dcomps);

            auto g_zz_0_yyzz_zz = cbuffer.data(gd_geom_20_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzz_xx, g_zz_0_yyzz_xy, g_zz_0_yyzz_xz, g_zz_0_yyzz_yy, g_zz_0_yyzz_yz, g_zz_0_yyzz_zz, g_zz_0_yzz_xx, g_zz_0_yzz_xxy, g_zz_0_yzz_xy, g_zz_0_yzz_xyy, g_zz_0_yzz_xyz, g_zz_0_yzz_xz, g_zz_0_yzz_yy, g_zz_0_yzz_yyy, g_zz_0_yzz_yyz, g_zz_0_yzz_yz, g_zz_0_yzz_yzz, g_zz_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzz_xx[k] = -g_zz_0_yzz_xx[k] * ab_y + g_zz_0_yzz_xxy[k];

                g_zz_0_yyzz_xy[k] = -g_zz_0_yzz_xy[k] * ab_y + g_zz_0_yzz_xyy[k];

                g_zz_0_yyzz_xz[k] = -g_zz_0_yzz_xz[k] * ab_y + g_zz_0_yzz_xyz[k];

                g_zz_0_yyzz_yy[k] = -g_zz_0_yzz_yy[k] * ab_y + g_zz_0_yzz_yyy[k];

                g_zz_0_yyzz_yz[k] = -g_zz_0_yzz_yz[k] * ab_y + g_zz_0_yzz_yyz[k];

                g_zz_0_yyzz_zz[k] = -g_zz_0_yzz_zz[k] * ab_y + g_zz_0_yzz_yzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzz_xx = cbuffer.data(gd_geom_20_off + 528 * ccomps * dcomps);

            auto g_zz_0_yzzz_xy = cbuffer.data(gd_geom_20_off + 529 * ccomps * dcomps);

            auto g_zz_0_yzzz_xz = cbuffer.data(gd_geom_20_off + 530 * ccomps * dcomps);

            auto g_zz_0_yzzz_yy = cbuffer.data(gd_geom_20_off + 531 * ccomps * dcomps);

            auto g_zz_0_yzzz_yz = cbuffer.data(gd_geom_20_off + 532 * ccomps * dcomps);

            auto g_zz_0_yzzz_zz = cbuffer.data(gd_geom_20_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzz_xx, g_zz_0_yzzz_xy, g_zz_0_yzzz_xz, g_zz_0_yzzz_yy, g_zz_0_yzzz_yz, g_zz_0_yzzz_zz, g_zz_0_zzz_xx, g_zz_0_zzz_xxy, g_zz_0_zzz_xy, g_zz_0_zzz_xyy, g_zz_0_zzz_xyz, g_zz_0_zzz_xz, g_zz_0_zzz_yy, g_zz_0_zzz_yyy, g_zz_0_zzz_yyz, g_zz_0_zzz_yz, g_zz_0_zzz_yzz, g_zz_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzz_xx[k] = -g_zz_0_zzz_xx[k] * ab_y + g_zz_0_zzz_xxy[k];

                g_zz_0_yzzz_xy[k] = -g_zz_0_zzz_xy[k] * ab_y + g_zz_0_zzz_xyy[k];

                g_zz_0_yzzz_xz[k] = -g_zz_0_zzz_xz[k] * ab_y + g_zz_0_zzz_xyz[k];

                g_zz_0_yzzz_yy[k] = -g_zz_0_zzz_yy[k] * ab_y + g_zz_0_zzz_yyy[k];

                g_zz_0_yzzz_yz[k] = -g_zz_0_zzz_yz[k] * ab_y + g_zz_0_zzz_yyz[k];

                g_zz_0_yzzz_zz[k] = -g_zz_0_zzz_zz[k] * ab_y + g_zz_0_zzz_yzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzz_xx = cbuffer.data(gd_geom_20_off + 534 * ccomps * dcomps);

            auto g_zz_0_zzzz_xy = cbuffer.data(gd_geom_20_off + 535 * ccomps * dcomps);

            auto g_zz_0_zzzz_xz = cbuffer.data(gd_geom_20_off + 536 * ccomps * dcomps);

            auto g_zz_0_zzzz_yy = cbuffer.data(gd_geom_20_off + 537 * ccomps * dcomps);

            auto g_zz_0_zzzz_yz = cbuffer.data(gd_geom_20_off + 538 * ccomps * dcomps);

            auto g_zz_0_zzzz_zz = cbuffer.data(gd_geom_20_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_zz_0_zzz_xx, g_zz_0_zzz_xxz, g_zz_0_zzz_xy, g_zz_0_zzz_xyz, g_zz_0_zzz_xz, g_zz_0_zzz_xzz, g_zz_0_zzz_yy, g_zz_0_zzz_yyz, g_zz_0_zzz_yz, g_zz_0_zzz_yzz, g_zz_0_zzz_zz, g_zz_0_zzz_zzz, g_zz_0_zzzz_xx, g_zz_0_zzzz_xy, g_zz_0_zzzz_xz, g_zz_0_zzzz_yy, g_zz_0_zzzz_yz, g_zz_0_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzz_xx[k] = -2.0 * g_z_0_zzz_xx[k] - g_zz_0_zzz_xx[k] * ab_z + g_zz_0_zzz_xxz[k];

                g_zz_0_zzzz_xy[k] = -2.0 * g_z_0_zzz_xy[k] - g_zz_0_zzz_xy[k] * ab_z + g_zz_0_zzz_xyz[k];

                g_zz_0_zzzz_xz[k] = -2.0 * g_z_0_zzz_xz[k] - g_zz_0_zzz_xz[k] * ab_z + g_zz_0_zzz_xzz[k];

                g_zz_0_zzzz_yy[k] = -2.0 * g_z_0_zzz_yy[k] - g_zz_0_zzz_yy[k] * ab_z + g_zz_0_zzz_yyz[k];

                g_zz_0_zzzz_yz[k] = -2.0 * g_z_0_zzz_yz[k] - g_zz_0_zzz_yz[k] * ab_z + g_zz_0_zzz_yzz[k];

                g_zz_0_zzzz_zz[k] = -2.0 * g_z_0_zzz_zz[k] - g_zz_0_zzz_zz[k] * ab_z + g_zz_0_zzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

